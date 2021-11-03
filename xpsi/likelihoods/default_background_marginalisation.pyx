#cython: cdivision=True
#cython: boundscheck=False
#cython: nonecheck=False
#cython: wraparound=False
#cython: embedsignature=True

import numpy as np
cimport numpy as np
from libc.stdlib cimport malloc, free
from libc.math cimport exp, pow, log, sqrt, fabs, floor

ctypedef np.uint8_t uint8

from GSL cimport (gsl_interp,
                   gsl_interp_alloc,
                   gsl_interp_init,
                   gsl_interp_free,
                   gsl_interp_eval,
                   gsl_interp_eval_integ,
                   gsl_interp_accel,
                   gsl_interp_accel_alloc,
                   gsl_interp_accel_free,
                   gsl_interp_accel_reset,
                   gsl_function,
                   gsl_integration_cquad_workspace,
                   gsl_integration_cquad_workspace_alloc,
                   gsl_integration_cquad_workspace_free,
                   gsl_integration_cquad)

ctypedef gsl_interp_accel accel

cdef extern from "gsl/gsl_sf_gamma.h":

    double gsl_sf_lnfact(const unsigned int n)

from ..tools cimport _get_phase_interpolant, gsl_interp_type

def precomputation(int[:,::1] data):
    """ Compute negative of sum of log-factorials of data count numbers.

    Use this function to perform precomputation before repeatedly calling
    :func:`~.eval_marginal_likelihood`.

    :param int[:,::1] data:
        Phase-channel resolved count numbers (integers). Phase increases with
        column number.

    :returns:
        A 1D :class:`numpy.ndarray` information, one element per channel. Each
        element is the negative of the sum (over phase intervals) of
        log-factorials of data count numbers.

    """

    cdef:
        size_t i, j
        double[::1] precomp = np.zeros(data.shape[0], dtype = np.double)

    for i in range(data.shape[0]):
        for j in range(data.shape[1]):
            precomp[i] += gsl_sf_lnfact(<unsigned int>(data[i,j]))

        precomp[i] *= -1.0

    return np.asarray(precomp, dtype = np.double, order = 'C')

ctypedef struct args:
    size_t n
    size_t i
    double SCALE
    double *precomp
    double *data
    double *star
    double B
    double interval
    double T_exp
    double std
    double A

cdef double marginal_integrand(double B, void *params) nogil:

    cdef int j
    cdef double c, x = 0.0
    cdef args *a  = <args*> params

    for j in range(a.n):
        c = a.SCALE * (a.star[j] + B)
        if c == 0.0 and a.data[j] == 0.0:
            x += 1.0
        else:
            x += a.data[j] * log(c) - c

    #with gil:
    #    if x - a.A > 0.0:
    #        print("a.i, x-a.A, a.B, off, B = %i, %.8f, %.8f, %.8f, %.8f" %
    #              (a.i, x - a.A, a.B, (a.B - B)/a.interval, B))

    return exp(x - a.A)

cdef double delta(double B, void *params) nogil:

    cdef int j
    cdef double x = 0.0, y = 0.0
    cdef args *a  = <args*> params

    for j in range(a.n):
        y += a.data[j] / (a.star[j] + B)
        x += a.data[j] / pow(a.star[j] + B, 2.0)

    y = 2.0 * a.T_exp - 2.0 * y
    x *= 2.0

    a.std = sqrt(2.0 / x)

    return -1.0 * y / x

def eval_marginal_likelihood(double exposure_time,
                             double[::1] phases,
                             double[:,::1] counts,
                             components,
                             component_phases,
                             phase_shifts,
                             double[::1] neg_sum_ln_data_factorial,
                             double[:,::1] support,
                             size_t workspace_intervals,
                             double epsabs,
                             double epsrel,
                             double epsilon,
                             double sigmas,
                             double llzero,
                             allow_negative = False,
                             background = None):
    """ Evaluate the Poisson likelihood.

    The count rate is integrated over phase intervals.

    A Newton iteration procedure is implemented channel-by-channel to
    approximate the conditional-ML background count-rate variable given
    a fiducial estimate. Marginalisation over each background variable (in each
    channel) is then executed numerically between bounds centred on the ML
    background estimate. The bounds are based on a Gaussian expansion of the
    conditional likelihood function, and are :obj:`sigmas` standard deviations
    above and below the ML estimate. The marginalisation is with respect to a
    flat bounded or unbounded prior density function of each background
    variable.

    :param double exposure_time:
        Exposure time in seconds by which to scale the expected count rate
        in each phase interval.

    :param double[::1] phases:
        A :class:`numpy.ndarray` of phase interval edges in cycles.

    :param tuple components:
        Component signals, each a C-contiguous :class:`numpy.ndarray` of
        signal count rates where phase increases with column number.

    :param tuple component_phases:
        For each component, a C-contiguous :class:`numpy.ndarray` of phases
        in cycles at which the model :obj:`signal` is evaluated on
        the interval ``[0,1]``.

    :param array-like phase_shifts:
        Phase shifts in cycles, such as on the interval ``[-0.5,0.5]``, for
        the component signals.

    :param double[::1] neg_sum_ln_data_factorial:
        The precomputed output of :func:`~.precomputation` given the data count
        numbers.

    :param double[:,::1] support:
        The prior support of the background *count-rate* variables. The first
        column contains lower-bounds as a function of channel number. Lower-
        bounds must be greater than or equal to zero. The second column contains
        the upper-bounds as a function of channel number. Upper-bounds for a
        proper prior must be greater than zero and greater than the
        corresponding lower-bound but finite. If the upper-bound is less than
        zero the prior is semi-unbounded and thus improper. Must be a
        C-contiguous :class:`numpy.ndarray`.

    :param size_t workspace_intervals:
        The size of the workspace to allocate for marginalisation via numerical
        quadrature using the GSL ``cquad`` integration routine.

    :param double epsabs:
        The absolute tolerance for marginalisation via numerical quadrature
        using the GSL ``cquad`` integration routine.

    :param double epsrel:
        The relative tolerance for marginalisation via numerical quadrature
        using the GSL ``cquad`` integration routine.

    :param double epsilon:
        The fraction of the standard deviation of the approximating Gaussian
        function to tolerate when a new step of the quadratic maximisation
        scheme is proposed. This fraction should be some adequately small
        number. The standard deviation is recalculated with each iteration
        as the position (in each background variable) evolves.

    :param double sigmas:
        The number of approximating Gaussian standard deviations to expand
        the integration domain about the point estimated to maximise the
        conditional likelihood function in each channel. This number should
        probably be at least five but no more than ten.

    :param double llzero:
        The log-likelihood that a MultiNest process treats as zero and thus
        ignores points with smaller log-likelihood. This number will be *very*
        negative. This is useful because if the likelihood is predicted, in
        advance of full executation of this function, to be incredibly small,
        computation can be avoided, returning a number slightly above this
        zero-threshold.

    :param obj allow_negative:
        A boolean or an array of booleans, one per component, declaring whether
        to allow negative phase interpolant integrals. If the interpolant is
        not a Steffen spline, then the interpolant of a non-negative function
        can be negative due to oscillations. For the default Akima Periodic
        spline from GSL, such oscillations should manifest as small relative
        to those present in cubic splines, for instance, because it is
        designed to handle a rapidly changing second-order derivative.

    :param obj background:
        If not ``None``, then a C-contiguous :class:`numpy.ndarray` of
        background count rates where phase interval increases with column
        number. Useful for phase-dependent backgrounds, or a phase-independent
        background if the channel-by-channel background variable prior support
        is restricted.

    :returns:
        A tuple ``(double, 2D ndarray, 1D ndarray)``. The first element is
        the logarithm of the marginal likelihood. The second element is the
        expected count numbers in joint phase-channel intervals from the star
        (the target source). The last element is the vector of background
        count numbers that are estimated to maximise the conditional
        likelihood function, one per channel.

    """

    cdef:
        size_t i, j, p
        double LOGLIKE = 0.0, av_STAR, av_DATA
        double pa, pb, c

        double[:,::1] STAR = np.zeros((components[0].shape[0], phases.shape[0] - 1),
                                      dtype = np.double)
        double[::1] MCL_BACKGROUND = np.zeros(components[0].shape[0], dtype = np.double)
        double[::1] MCL_BACKGROUND_GIVEN_SUPPORT = np.zeros(components[0].shape[0], dtype = np.double)

        double[:,::1] _background

        double n = <double>(phases.shape[0] - 1)
        double SCALE = exposure_time / n
        double _val

    if background is not None:
        _background = background

    cdef args a
    a.n = <size_t>n
    a.precomp = &(neg_sum_ln_data_factorial[0])
    a.SCALE = SCALE

    cdef double result, abserr, upper, lower, std_est
    cdef size_t nevals
    cdef gsl_function f

    f.function = &marginal_integrand

    cdef const gsl_interp_type *_interpolant

    _interpolant = _get_phase_interpolant()

    cdef double *phases_ptr = NULL
    cdef double *pulse_ptr = NULL

    cdef double B, B_for_integrand, dB, B_min, min_counts, limit
    cdef int counter

    cdef size_t num_components = len(components)
    cdef double[:,::1] pulse
    cdef double[::1] pulse_phase_set
    cdef double phase_shift
    cdef uint8[::1] _allow_negative = np.zeros(num_components, dtype=np.uint8)

    if isinstance(allow_negative, bool):
        for i in range(num_components):
            _allow_negative[i] = <uint8>allow_negative
    else:
        try:
            len(allow_negative)
        except TypeError:
            raise TypeError('An iterable is required to specify component-by-'
                            'component positivity.')
        else:
            if len(allow_negative) != num_components:
                raise ValueError('Number of allow_negative declarations does '
                                 'not match the number of components..')

            for i in range(num_components):
                _allow_negative[i] = allow_negative[i]

    cdef gsl_interp **interp = <gsl_interp**> malloc(num_components * sizeof(gsl_interp*))
    cdef accel **acc =  <accel**> malloc(num_components * sizeof(accel*))
    cdef gsl_integration_cquad_workspace *w = gsl_integration_cquad_workspace_alloc(workspace_intervals)

    for p in range(num_components):
        pulse_phase_set = component_phases[p]
        interp[p] = gsl_interp_alloc(_interpolant,
                                     pulse_phase_set.shape[0])
        acc[p] = gsl_interp_accel_alloc()
        gsl_interp_accel_reset(acc[p])

    cdef gsl_interp *inter_ptr = NULL
    cdef accel *acc_ptr = NULL

    for i in range(STAR.shape[0]):
        for p in range(num_components):
            pulse = components[p]
            pulse_phase_set = component_phases[p]
            phase_shift = phase_shifts[p]

            interp_ptr = interp[p]
            acc_ptr = acc[p]
            phases_ptr = &(pulse_phase_set[0])
            pulse_ptr = &(pulse[i,0])

            gsl_interp_init(interp_ptr, phases_ptr, pulse_ptr,
                            pulse_phase_set.shape[0])

            for j in range(phases.shape[0] - 1):
                pa = phases[j] + phase_shift
                pb = phases[j+1] + phase_shift

                if pb - pa == 1.0:
                    pa = 0.0
                    pb = 1.0
                else:
                    pa -= floor(pa)
                    pb -= floor(pb)

                if pa < pb:
                    _val = gsl_interp_eval_integ(interp_ptr,
                                                       phases_ptr,
                                                       pulse_ptr,
                                                       pa, pb,
                                                       acc_ptr)
                    if _val > 0.0 or _allow_negative[p] == 1:
                        STAR[i,j] += _val
                else:
                    _val = gsl_interp_eval_integ(interp_ptr,
                                                       phases_ptr,
                                                       pulse_ptr,
                                                       pa, 1.0,
                                                       acc_ptr)
                    if _val > 0.0 or _allow_negative[p] == 1:
                        STAR[i,j] += _val

                    _val = gsl_interp_eval_integ(interp_ptr,
                                                       phases_ptr,
                                                       pulse_ptr,
                                                       0.0, pb,
                                                       acc_ptr)
                    if _val > 0.0 or _allow_negative[p] == 1:
                        STAR[i,j] += _val

        for j in range(phases.shape[0] - 1): # interpolant safety procedure
            if STAR[i,j] < 0.0:
                STAR[i,j] = 0.0

        av_DATA = 0.0; av_STAR = 0.0

        for j in range(phases.shape[0] - 1):
            STAR[i,j] *= n
            if background is not None:
                STAR[i,j] += _background[i,j]

            av_STAR += STAR[i,j]
            av_DATA += counts[i,j]

        limit = av_STAR * SCALE - 20.0 * sqrt(av_STAR * SCALE) - av_DATA

        if limit > 0.0:
            LOGLIKE = llzero * (0.1 + 0.9 * np.random.rand(1))
            break

        av_STAR /= n
        av_DATA /= exposure_time

        a.i = i
        a.data = &(counts[i,0])
        a.star = &(STAR[i,0])
        a.T_exp = exposure_time

        if av_DATA == 0.0 and av_STAR == 0.0:
            # zero counts in channel and zero hot region signal

            lower = 0.0
            if lower < support[i,0]:
                lower = support[i,0]
            upper = 10.0 / exposure_time
            if upper > support[i,1]:
                upper = support[i,1]

            B = 0.0

            LOGLIKE += log( ( exp(-1.0 * lower * exposure_time) - exp(-1.0 * upper * exposure_time)) / exposure_time )
        else:
            B_min = 0.0
            B = av_DATA - av_STAR

            if B <= B_min:
                min_counts = -1.0
                for j in range(0, a.n):
                    if STAR[i,j] == 0.0:
                        min_counts = -2.0
                        break
                if min_counts == -2.0:
                    for j in range(0, a.n):
                        if (min_counts == -2.0 and counts[i,j] > 0.0) or (0.0 < counts[i,j] < min_counts):
                            min_counts = counts[i,j]
                if min_counts != -1.0:
                    B = 0.01 * min_counts / SCALE
                    B_min = 0.1 * B
                    #print("i, B = %i, %.16f" % (a.i,B))
                else:
                    B = B_min

            dB = delta(B, &a)

            #print("B_0 = %.16e" % B)
            counter = 0
            while fabs(dB) > epsilon*a.std and counter < 2:
                B += dB
                #print("i, B = %i, %.8f" % (a.i,B))
                if B < B_min:
                    if B_min > 0.0:
                        counter += 1
                    else:
                        counter = 2
                    B = B_min

                dB = delta(B, &a)

            std_est = 0.0
            for j in range(a.n):
                #print("i, j, STAR[i,j] + B = %i, %i, %.16e, %.16e" % (a.i,j,STAR[i,j], B))
                std_est += counts[i,j] / pow(STAR[i,j] + B, 2.0)

            std_est = sqrt(1.0 / std_est)
            #print("i, B, std = %i, %.16e, %.16e" % (a.i,B,std_est))

            a.B = B
            lower = B - sigmas * std_est

            if lower < B_min:
                lower = B_min

            upper = B + sigmas * std_est

            #a.interval = upper - lower

            B_for_integrand = B

            if lower < support[i,0]:
                lower = support[i,0]
                if upper < support[i,0] and support[i,1] > 0.0:
                    upper = support[i,1]
                    B_for_integrand = support[i,0]
                elif upper < support[i,0]:
                    upper = support[i,0] + sigmas * std_est
                    B_for_integrand = support[i,0]

            if upper > support[i,1] and support[i,1] > 0.0:
                upper = support[i,1]
                if lower > support[i,1]:
                    lower = support[i,0]
                    B_for_integrand = support[i,1]

            f.params = &a

            a.A = 0.0
            for j in range(a.n):
                c = a.SCALE * (a.star[j] + B_for_integrand)
                if c == 0.0 and a.data[j] == 0.0:
                    a.A += 1.0
                else:
                    a.A += a.data[j] * log(c) - c

            gsl_integration_cquad(&f, lower, upper,
                                  epsabs, epsrel,
                                  w, &result,
                                  &abserr, &nevals)
            if result > 0.0:
                LOGLIKE += log(result) + a.A + a.precomp[a.i]
            else:
                LOGLIKE = llzero * (0.1 + 0.9 * np.random.rand(1))
                break

        MCL_BACKGROUND[i] = B * exposure_time

        if B < support[i,0]:
            B = support[i,0]
        elif B > support[i,1] and support[i,1] > 0.0:
            B = support[i,1]

        for j in range(a.n):
            STAR[i,j] = a.SCALE * (a.star[j] + B)

        MCL_BACKGROUND_GIVEN_SUPPORT[i] = B * exposure_time

    for p in range(num_components):
        gsl_interp_accel_free(acc[p])
        gsl_interp_free(interp[p])

    free(acc)
    free(interp)

    gsl_integration_cquad_workspace_free(w)

    return (LOGLIKE,
            np.asarray(STAR, dtype=np.double, order='C'),
            np.asarray(MCL_BACKGROUND, dtype=np.double, order='C'),
            np.asarray(MCL_BACKGROUND_GIVEN_SUPPORT, dtype=np.double, order='C'))
