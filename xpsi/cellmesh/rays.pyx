#cython: cdivision=True
#cython: boundscheck=False
#cython: nonecheck=False
#cython: wraparound=False

import numpy as np
cimport numpy as np
cimport cython
from cython.parallel cimport *
from libc.math cimport M_PI, sqrt, cos, asin, acos, log, atan, NAN, pow
from libc.stdlib cimport malloc, free
from libc.stdio cimport printf
from libcpp cimport bool
import xpsi

from ..tools.core cimport are_equal

from GSL cimport (gsl_function,
                  gsl_integration_workspace,
                  gsl_integration_workspace_alloc,
                  gsl_integration_workspace_free,
                  gsl_integration_cquad_workspace,
                  gsl_integration_qag,
                  gsl_interp_accel,
                  gsl_interp_accel_alloc,
                  gsl_interp_accel_free,
                  gsl_interp_accel_reset,
                  gsl_interp_eval,
                  gsl_interp_eval_deriv,
                  gsl_interp_alloc,
                  gsl_interp,
                  gsl_interp_steffen,
                  gsl_interp_init,
                  gsl_interp_free,
                  gsl_isnan,
                  GSL_INTEG_GAUSS61,
                  GSL_EFAILED,
                  gsl_set_error_handler,
                  gsl_set_error_handler_off,
                  gsl_error_handler_t)

ctypedef gsl_interp interp
ctypedef gsl_interp_accel accel
ctypedef gsl_integration_workspace gsl_work
ctypedef gsl_integration_cquad_workspace gsl_cq_work

cdef double _pi = M_PI
cdef double _hlfpi = M_PI / 2.0
cdef double c = 2.99792458e8

cdef int ERROR = 1

# Ray integrals in the Schwarzschild metric for the Oblate Schwarzschild +
# Doppler approximation (Morsink et al. 2007, ApJ 663, 1244; hereafter M07).
# For each surface radius R and emission angle alpha (from the radial
# direction) we compute the light-bending angle psi (M07 Eqs. 19-20) and the
# arrival-time delay relative to a radial photon from the same R (M07
# Eq. 38). The integrands below are those of Salmi et al. (2018, A&A 618,
# A161; hereafter S18), Appendix A, which substitutes y = R/r and then
# x = sqrt(1 - y), so that they stay finite at the turning point. Their
# constant prefactors are applied in rayIntegrator, lengths are in units of
# r_s, and u/(1 - u) in S18 (u = r_s/R) appears as 1/(R/r_s - 1). As in
# S18 Eq. (9), each lag integral is relative to a radial photon emitted from
# the radius at its lower integration limit; the shift to a common reference
# radius (S18 Eq. 10) is applied in HotRegion.__calibrate_lag.

# Map the primary-image deflection psi in [0, pi] to the deflection of image
# order `order` (photons bent by more than pi; order 0 is the primary image).
cdef double eval_image_deflection(int order, double psi) noexcept nogil:
    if order % 2 != 0:
        return <double>(order + 1) * _pi + pow(-1.0, <double>order) * psi
    else:
        return <double>(order) * _pi + pow(-1.0, <double>order) * psi

cdef double b_phsph_over_r_s = 3.0 * sqrt(3.0) / 2.0

# inDef and inLag are outDef and outLag with the emission radius R replaced by
# the turning radius r_c, where alpha = pi/2 (cos alpha = 0), i.e. with
# y = r_c/r (r_c is p, the distance of closest approach, in S18 Eq. 4).
# Integrating from x = 0 gives psi_p(r_c, pi/2) and Delta t_p(r_c, pi/2),
# used for alpha > pi/2 (S18 Eqs. 3 and 12); integrating from
# x = sqrt(1 - r_c/R) covers R -> infinity, i.e. the same as outDef/outLag
# for the same impact parameter, but without their steep rise near x = 0 as
# alpha -> pi/2. This second use is not in S18, which uses Eqs. (A.2) and
# (A.7) at R for all alpha <= pi/2. The r_c form is used whenever
# b > b_photon-sphere, i.e. whenever the geodesic has a turning point r_c
# (reached by the photon only if alpha > pi/2).

# Deflection integrand: that of S18 Eq. (A.2) with cos alpha = 0, which
# reduces to 1/sqrt(q), with q = 2 - x^2 - u(1 - x^2)^2/(1 - u) and
# u = r_s/r_c (S18 Eq. A.3 with sin alpha = 1).
cdef double inDef(double x, void *params) noexcept nogil:

    cdef double r_c_over_r_s = (<double*> params)[0]
    return 1.0 / sqrt(2.0 - x * x - (1.0 - x * x) * (1.0 - x * x) / (r_c_over_r_s - 1.0))

# Lag integrand: that of S18 Eq. (A.7) with cos alpha = 0, which reduces to
# 1/(q x + sqrt(q)), written as X^2/(x + X) with X = 1/sqrt(q), and with
# q = 2 - x^2 - u(1 - x^2)^2/(1 - u) and u = r_s/r_c (S18 Eq. A.3 with
# sin alpha = 1).
cdef double inLag(double x, void *params) noexcept nogil:

    cdef double r_c_over_r_s = (<double*> params)[0]

    cdef double X = 1.0 / sqrt(2.0 - x * x - (1.0 - x * x) * (1.0 - x * x) / (r_c_over_r_s - 1.0))

    return X * X / (x + X)

# Deflection integrand: S18 Eqs. (A.2)-(A.3), x/sqrt(cos^2 alpha + x^2 q),
# with sas = sin^2 alpha and q = sas * [2 - x^2 - u(1 - x^2)^2/(1 - u)].
cdef double outDef(double x, void *params) noexcept nogil:

    cdef double sas = (<double*> params)[0]
    cdef double R_over_r_s = (<double*> params)[1]

    return x / sqrt(1.0 - sas + x * x * sas * (2.0 - x * x - (1.0 - x * x) * (1.0 - x * x) / (R_over_r_s - 1.0)))

# Lag integrand: S18 Eq. (A.7), x/(f + f^2) with f = sqrt(cos^2 alpha + x^2 q),
# sas = sin^2 alpha and q = sas * [2 - x^2 - u(1 - x^2)^2/(1 - u)].
cdef double outLag(double x, void *params) noexcept nogil:

    cdef double sas = (<double*> params)[0]
    cdef double R_over_r_s = (<double*> params)[1]

    cdef double f = sqrt(1.0 - sas + x * x * sas * (2.0 - x * x - (1.0 - x * x) * (1.0 - x * x) / (R_over_r_s - 1.0)))

    return x / (f + f*f)

cdef void rayIntegrator(size_t thread,
                             gsl_work *w,
                             double cos_alpha,
                             double r_s,
                             double r_s_over_R,
                             double *rayParams) noexcept nogil:

    cdef int status
    cdef double outParams[2]
    cdef double inParams[1]
    cdef gsl_function f_outDef, f_outLag, f_inDef, f_inLag
    cdef:
        double DeltaDef[2]
        double DeltaDefError
        double DeltaLag[2]
        double DeltaLagError, DeltaLagCor
        double r_c_over_r_s, b_over_r_s, w_R, sin_alpha, sin_alpha_sq, alpha
        double deltaDeflection, deltaLag

    # Returns, per thread, rayParams = [psi, error, lag in seconds, error].
    sin_alpha_sq = 1.0 - cos_alpha * cos_alpha
    sin_alpha = sqrt(sin_alpha_sq)
    alpha = acos(cos_alpha)

    # Impact parameter b = R sin(alpha) / sqrt(1 - r_s/R) (M07 Eq. 19).
    b_over_r_s = sin_alpha / (r_s_over_R * sqrt(1.0 - r_s_over_R))

    # Three cases:
    # 1. b <= b_photon-sphere (3*sqrt(3)/2 r_s): outgoing photon with no
    #    turning point; integrate M07 Eqs. 20 and 38 from R to infinity.
    # 2. b > b_photon-sphere, alpha <= pi/2: outgoing photon. The integrals
    #    from R to infinity are written in terms of the turning radius r_c
    #    (< R), which keeps them well behaved as alpha -> pi/2.
    # 3. alpha > pi/2: initially ingoing photon, possible on an oblate
    #    surface (M07 Sect. 3.1). It moves inwards to r_c and then out, so
    #    psi = 2 * (r_c -> infinity) - (R -> infinity), i.e.
    #    2 Delta-psi + psi(b, R) as in M07 Eq. 44; likewise for the lag
    #    (M07 Eq. 45); equivalently S18 Eqs. (3) and (12). These are
    #    evaluated numerically, not with the series approximations of M07
    #    Eqs. 46-47.
    # r_c is the largest root of r_c^3 - b^2 r_c + b^2 r_s = 0, i.e.
    # r_c = b sqrt(1 - r_s/r_c) (M07 Eq. 41), from the trigonometric solution
    # (S18 Eq. 4, in an equivalent form).
    if (b_over_r_s <= b_phsph_over_r_s):

        outParams[0] = sin_alpha_sq
        outParams[1] = 1.0 / r_s_over_R

        f_outDef.function = &outDef
        f_outLag.function = &outLag

        f_outDef.params = &outParams
        f_outLag.params = &outParams

        status = gsl_integration_qag(&f_outDef, 0.0, 1.0,
                                     0, 1.0E-12, 100, GSL_INTEG_GAUSS61, w,
                                     &rayParams[4 * thread],
                                     &rayParams[4 * thread + 1])

        if (status == GSL_EFAILED):
            rayParams[4 * thread] = NAN

        if (status != GSL_EFAILED):
            status = gsl_integration_qag(&f_outLag, 0.0, 1.0,
                                         0, 1.0E-12, 100, GSL_INTEG_GAUSS61, w,
                                         &rayParams[4 * thread + 2],
                                         &rayParams[4 * thread + 3])

            if (status == GSL_EFAILED):
                rayParams[4 * thread] = NAN

        if (status != GSL_EFAILED):
            # S18 Eq. (A.2) prefactor 2 sin(alpha)/sqrt(1 - u) = 2 b/R.
            rayParams[4 * thread] *= 2.0 * b_over_r_s * r_s_over_R
            # S18 Eq. (A.7) prefactor 2 R sin^2(alpha)/(1 - u) = 2 b^2/R
            # (its 1/c is applied later in this function).
            rayParams[4 * thread + 2] *= 2.0 * b_over_r_s * b_over_r_s * r_s_over_R * r_s

    elif (alpha <= _pi / 2.0):

        f_outDef.function = &inDef
        f_outLag.function = &inLag

        if are_equal(sin_alpha,1.0):
            r_c_over_r_s = 1.0 / r_s_over_R
            w_R = 0.0
        else:
            r_c_over_r_s = 2.0 * b_over_r_s * cos((atan(sqrt(4.0 * b_over_r_s * b_over_r_s / 27.0 - 1.0)) - _pi) / 3.0) / sqrt(3.0)
            w_R = sqrt(1.0 - r_c_over_r_s * r_s_over_R)

        inParams[0] = r_c_over_r_s

        f_outDef.params = &inParams
        f_outLag.params = &inParams

        if gsl_isnan(w_R):
            w_R = 0.0

        status = gsl_integration_qag(&f_outDef, w_R, 1.0,
                                     0, 1.0E-12, 100, GSL_INTEG_GAUSS61, w,
                                     &rayParams[4 * thread],
                                     &rayParams[4 * thread + 1])

        if (status == GSL_EFAILED):
            rayParams[4 * thread] = NAN

        if (status != GSL_EFAILED):
            status = gsl_integration_qag(&f_outLag, w_R, 1.0,
                                         0, 1.0E-12, 100, GSL_INTEG_GAUSS61, w,
                                         &rayParams[4 * thread + 2],
                                         &rayParams[4 * thread + 3])

            if (status == GSL_EFAILED):
                rayParams[4 * thread] = NAN

        if (status != GSL_EFAILED):
            # Prefactors of S18 Eqs. (A.2) and (A.7) with R -> r_c and
            # sin(alpha) = 1: 2/sqrt(1 - r_s/r_c) = 2 b/r_c and
            # 2 r_c/(1 - r_s/r_c) = 2 b^2/r_c (the 1/c of Eq. (A.7) is
            # applied later in this function).
            rayParams[4 * thread] *= 2.0 * b_over_r_s / r_c_over_r_s

            rayParams[4 * thread + 2] *= 2.0 * b_over_r_s * b_over_r_s * r_s / r_c_over_r_s

    elif (alpha > _pi / 2.0):

        f_inDef.function = &inDef
        f_inLag.function = &inLag

        r_c_over_r_s = 2.0 * b_over_r_s * cos((atan(sqrt(4.0 * b_over_r_s * b_over_r_s / 27.0 - 1.0)) - _pi) / 3.0) / sqrt(3.0)

        inParams[0] = r_c_over_r_s

        f_inDef.params = &inParams
        f_inLag.params = &inParams

        w_R = sqrt(1.0 - r_c_over_r_s * r_s_over_R)

        if gsl_isnan(w_R):
            w_R = 0.0

        status = gsl_integration_qag(&f_inDef, 0.0, 1.0,
                                     0, 1.0E-12, 100, GSL_INTEG_GAUSS61, w,
                                     &DeltaDef[1], &DeltaDefError)

        if (status == GSL_EFAILED):
            rayParams[4 * thread] = NAN

        if (status != GSL_EFAILED):
            status = gsl_integration_qag(&f_inDef, w_R, 1.0,
                                         0, 1.0E-12, 100, GSL_INTEG_GAUSS61, w,
                                         &DeltaDef[0], &DeltaDefError)

            if (status == GSL_EFAILED):
                rayParams[4 * thread] = NAN

        if (status != GSL_EFAILED):
            status = gsl_integration_qag(&f_inLag, 0.0, 1.0,
                                         0, 1.0E-12, 100, GSL_INTEG_GAUSS61, w,
                                         &DeltaLag[1], &DeltaLagError)

            if (status == GSL_EFAILED):
                rayParams[4 * thread] = NAN

        if (status != GSL_EFAILED):
            status = gsl_integration_qag(&f_inLag, w_R, 1.0,
                                         0, 1.0E-12, 100, GSL_INTEG_GAUSS61, w,
                                         &DeltaLag[0], &DeltaLagError)

            if (status == GSL_EFAILED):
                rayParams[4 * thread] = NAN

        if (status != GSL_EFAILED):
            # r_c-form prefactors 2 b/r_c and 2 b^2/r_c, as above.
            rayParams[4 * thread] = 2.0 * b_over_r_s * (2.0 * DeltaDef[1] - DeltaDef[0]) / r_c_over_r_s

            rayParams[4 * thread + 2] = 2.0 * b_over_r_s * b_over_r_s * r_s * (2.0 * DeltaLag[1] - DeltaLag[0]) / r_c_over_r_s

            # Radial light-travel time between r_c and R (in units of
            # r_s/c), added twice so that the lag stays relative to a radial
            # photon emitted from R. This is 2 delta t(r_c, R) of S18
            # Eq. (10), i.e. with r_ref = R.
            DeltaLagCor = 1.0/r_s_over_R - r_c_over_r_s + log((1.0/r_s_over_R - 1.0) / (r_c_over_r_s - 1.0))

            rayParams[4 * thread + 2] += 2.0 * r_s * DeltaLagCor

    if (status != GSL_EFAILED):
        rayParams[4 * thread + 2] /= c

def compute_rays(size_t N_T,
                 size_t N_u,
                 double r_s,
                 double[::1] r_s_over_R,
                 double[::1] maxAlpha,
                 size_t N_R):
    """
    Compute ray properties in the Schwarzschild spacetime.

    Used by the Oblate Schwarzschild + Doppler approximation (Morsink et al.
    2007): rays are traced in the Schwarzschild metric from each surface
    radius :math:`R(\\theta)` of the oblate star. See
    :ref:`spacetime_models`.

    :param int: Number of OpenMP threads.
    :param int: Number of :math:`u\mathrel{:=}r_{s}/R` coordinate values
    :param float: Schwarzschild radius :math:`r_{s}` in metres.
    :param obj: 1D :class:`numpy.ndarray` :math:`u`
    :param obj: 1D :class:`numpy.ndarray` of the *maximum* emission angle
                :math:`\\alpha` (radians) for each :math:`u`; this exceeds
                :math:`\\pi/2` on an oblate surface, where initially ingoing
                photons can be visible.
    :param int: Number of rays per :math:`u` value

    :return: Termination flag
    :return: Deflection angle as a 2D :class:`numpy.ndarray`
    :return: Cosine of local ray angle as a 2D :class:`numpy.ndarray`
    :return: Ray lag in seconds, relative to a radial ray from the same
             :math:`R` (Morsink et al. 2007, Eq. 38), as a 2D
             :class:`numpy.ndarray`
    :return: Maximum deflection angle as a 1D :class:`numpy.ndarray`

    """

    cdef:
        signed int ii
        size_t i, j, k, T
        gsl_work **w = <gsl_work**> malloc(N_T * sizeof(gsl_work*))
        double alpha, cos_alpha_inc, delta, extreme
        int terminate = 0
        int *terminate_thread = <int*> malloc(N_T * sizeof(int))
        cdef double _epsilon = 1.0e-8

    for i in range(N_T):
        w[i] = gsl_integration_workspace_alloc(100)
        terminate_thread[i] = 0

    cdef double *rayParams = <double*> malloc(N_T * 4 * sizeof(double))
    cdef size_t *interp_index = <size_t*> malloc(N_T * sizeof(size_t))
    cdef size_t *interp_counter = <size_t*> malloc(N_T * sizeof(size_t))

    for i in range(N_T):
        interp_counter[i] = 0

    cdef:
        double[:,::1] deflection = np.empty((N_u, N_R), dtype = np.double)
        double[:,::1] cos_alpha = np.empty((N_u, N_R), dtype = np.double)
        double[:,::1] lag = np.empty((N_u, N_R), dtype = np.double)
        double[::1] maxDeflection = np.empty(N_u, dtype = np.double)

    # set_gsl_error_handler outside multi-threading
    cdef gsl_error_handler_t *handler = gsl_set_error_handler_off()

    for ii in prange(<signed int>N_u,
                     nogil = True,
                     schedule = 'static',
                     num_threads = N_T,
                     chunksize = 1):

        T = threadid()
        i = <size_t> ii

        if r_s_over_R[i] < 2.0/3.0:
            extreme = _pi - asin(sqrt(1.0 - r_s_over_R[i]) * b_phsph_over_r_s * r_s_over_R[i])
        elif are_equal(r_s_over_R[i], 2.0/3.0):
            extreme = _pi/2.0
        elif 2.0/3.0 < r_s_over_R[i] < 1.0:
            extreme = asin(sqrt(1.0 - r_s_over_R[i]) * b_phsph_over_r_s * r_s_over_R[i])

        if maxAlpha[i] >= (1.0 - _epsilon) * extreme:
            maxAlpha[i] = (1.0 - _epsilon) * extreme

        cos_alpha_inc = (1.0 - cos(maxAlpha[i])) / (<double>N_R - 1.0)

        for j in range(N_R):

            cos_alpha[i,j] = 1.0 - (<double>j) * cos_alpha_inc

            if j == 0:
                cos_alpha[i,j] = 1.0
                deflection[i,j] = 0.0
                lag[i,j] = 0.0
            else:
                if terminate_thread[T] == 0:
                    rayIntegrator(T, w[T], cos_alpha[i,j], r_s, r_s_over_R[i], rayParams)
                    if j == 0 and gsl_isnan(rayParams[4*T]) == 1:
                        terminate_thread[T] = ERROR
                    elif (j == N_R - 1 and
                          (gsl_isnan(rayParams[4*T]) == 1 or
                           rayParams[4*T] <= deflection[i,j - 1])):

                        terminate_thread[T] = ERROR
                    elif (interp_counter[T] > 0 and
                          gsl_isnan(rayParams[4*T]) == 0 and
                          (0.0 <= rayParams[4*T] < 100.0 * _pi) and
                          rayParams[4*T] > deflection[i,interp_index[T]]):

                        deflection[i,j] = rayParams[4*T]
                        lag[i,j] = rayParams[4*T + 2]

                        delta = <double>interp_counter[T] + 1
                        for k in range(1, interp_counter[T]+1):
                            deflection[i,interp_index[T] + k] = deflection[i,j] * <double>k
                            deflection[i,interp_index[T] + k] += deflection[i,interp_index[T]] * (delta - <double>k)
                            deflection[i,interp_index[T] + k] /= delta

                            lag[i,interp_index[T] + k] = lag[i,j] * <double>k
                            lag[i,interp_index[T] + k] += lag[i,interp_index[T]] * (delta - <double>k)
                            lag[i,interp_index[T] + k] /= delta

                        interp_counter[T] = 0 # reset counter
                    elif (j > 0 and
                          (gsl_isnan(rayParams[4*T]) == 1 or
                           rayParams[4*T] < 0.0 or
                           rayParams[4*T] > 100.0 * _pi)):

                        if interp_counter[T] == 0:
                            interp_index[T] = j - 1
                        interp_counter[T] += 1
                    else:
                        deflection[i,j] = rayParams[4*T]
                        lag[i,j] = rayParams[4*T + 2]

        maxDeflection[i] = deflection[i, N_R - 1]

    gsl_set_error_handler(handler)

    for i in range(N_T):
        gsl_integration_workspace_free(w[i])
        if terminate_thread[i] == ERROR:
            terminate = ERROR

    free(w)
    free(rayParams)
    free(interp_index)
    free(interp_counter)
    free(terminate_thread)

    return (terminate,
            np.asarray(deflection, dtype = np.double, order = 'C'),
            np.asarray(cos_alpha, dtype = np.double, order = 'C'),
            np.asarray(lag, dtype = np.double, order = 'C'),
            np.asarray(maxDeflection, dtype = np.double, order = 'C'))
