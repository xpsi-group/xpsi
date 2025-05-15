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

cdef double eval_image_deflection(int order, double psi) noexcept nogil:
    if order % 2 != 0:
        return <double>(order + 1) * _pi + pow(-1.0, <double>order) * psi
    else:
        return <double>(order) * _pi + pow(-1.0, <double>order) * psi

cdef double b_phsph_over_r_s = 3.0 * sqrt(3.0) / 2.0

cdef double inDef(double x, void *params) noexcept nogil:

    cdef double r_c_over_r_s = (<double*> params)[0]
    return 1.0 / sqrt(2.0 - x * x - (1.0 - x * x) * (1.0 - x * x) / (r_c_over_r_s - 1.0))

cdef double inLag(double x, void *params) noexcept nogil:

    cdef double r_c_over_r_s = (<double*> params)[0]

    cdef double X = 1.0 / sqrt(2.0 - x * x - (1.0 - x * x) * (1.0 - x * x) / (r_c_over_r_s - 1.0))

    return X * X / (x + X)

cdef double outDef(double x, void *params) noexcept nogil:

    cdef double sas = (<double*> params)[0]
    cdef double R_over_r_s = (<double*> params)[1]

    return x / sqrt(1.0 - sas + x * x * sas * (2.0 - x * x - (1.0 - x * x) * (1.0 - x * x) / (R_over_r_s - 1.0)))

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

    sin_alpha_sq = 1.0 - cos_alpha * cos_alpha
    sin_alpha = sqrt(sin_alpha_sq)
    alpha = acos(cos_alpha)

    b_over_r_s = sin_alpha / (r_s_over_R * sqrt(1.0 - r_s_over_R))

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
            rayParams[4 * thread] *= 2.0 * b_over_r_s * r_s_over_R
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
            rayParams[4 * thread] = 2.0 * b_over_r_s * (2.0 * DeltaDef[1] - DeltaDef[0]) / r_c_over_r_s

            rayParams[4 * thread + 2] = 2.0 * b_over_r_s * b_over_r_s * r_s * (2.0 * DeltaLag[1] - DeltaLag[0]) / r_c_over_r_s

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

    :param int: Number of OpenMP threads.
    :param int: Number of :math:`u\mathrel{:=}r_{s}/R` coordinate values
    :param obj: 1D :class:`numpy.ndarray` :math:`u`
    :param obj: 1D :class:`numpy.ndarray` of *maximum* :math:`\cos\\alpha`
    :param int: Number of rays per :math:`u` value

    :return: Termination flag
    :return: Deflection angle as a 2D :class:`numpy.ndarray`
    :return: Cosine of local ray angle as a 2D :class:`numpy.ndarray`
    :return: Ray lag relative to radial ray as a 2D :class:`numpy.ndarray`
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
