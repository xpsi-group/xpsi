#cython: cdivision=True
#cython: boundscheck=False
#cython: nonecheck=False
#cython: wraparound=False
#cython: embedsignature=True

import numpy as np
cimport numpy as np
from libc.math cimport pow, log, floor
from libc.stdlib cimport malloc, free

from GSL cimport (gsl_interp,
                   gsl_interp_alloc,
                   gsl_interp_init,
                   gsl_interp_free,
                   gsl_interp_eval,
                   gsl_interp_eval_integ,
                   gsl_interp_accel,
                   gsl_interp_accel_alloc,
                   gsl_interp_accel_free,
                   gsl_interp_accel_reset)

ctypedef gsl_interp_accel accel

ctypedef np.uint8_t uint8

from ..tools.core cimport _get_phase_interpolant, gsl_interp_type

def poisson_likelihood_given_background(double exposure_time,
                                        double[::1] phases,
                                        double[:,::1] counts,
                                        components,
                                        component_phases,
                                        phase_shifts,
                                        double[:,::1] background,
                                        allow_negative = False):
    """ Evaluate the Poisson likelihood.

    The count rate is integrated over phase intervals.

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

    :param double[:,::1] background:
        A C-contiguous :class:`numpy.ndarray` of background expected
        *counts*, whose shape matches the number of channels in each element
        of :obj:`components` and the number of phase intervals constructed
        from :obj:`phases`.

    :param obj allow_negative:
        A boolean or an array of booleans, one per component, declaring whether
        to allow negative phase interpolant integrals. If the interpolant is
        not a Steffen spline, then the interpolant of a non-negative function
        can be negative due to oscillations. For the default Akima Periodic
        spline from GSL, such oscillations should manifest as small relative
        to those present in cubic splines, for instance, because it is
        designed to handle a rapidly changing second-order derivative.

    :returns:
        A tuple ``(double, 2D ndarray)``. The first element is
        the logarithm of the marginal likelihood. The second element is the
        expected count numbers in joint phase-channel intervals from the star
        (the target source).

    """

    cdef:
        size_t i, j, p, num_components = len(components)
        double a, b

        double[:,::1] STAR = np.zeros((components[0].shape[0], phases.shape[0]-1),
                                       dtype = np.double)

    cdef const gsl_interp_type *_interpolant

    _interpolant = _get_phase_interpolant()

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
            if <size_t>len(allow_negative) != num_components:
                raise ValueError('Number of allow_negative declarations does '
                                 'not match the number of components..')

            for i in range(num_components):
                _allow_negative[i] = allow_negative[i]

    cdef double *phases_ptr = NULL
    cdef double *signal_ptr = NULL

    cdef double[:,::1] signal
    cdef double[::1] signal_phase_set
    cdef double phase_shift, _val

    cdef gsl_interp **interp = <gsl_interp**> malloc(num_components * sizeof(gsl_interp*))
    cdef accel **acc =  <accel**> malloc(num_components * sizeof(accel*))

    for p in range(num_components):
        signal_phase_set = component_phases[p]
        interp[p] = gsl_interp_alloc(_interpolant, signal_phase_set.shape[0])
        acc[p] = gsl_interp_accel_alloc()
        gsl_interp_accel_reset(acc[p])

    cdef gsl_interp *inter_ptr = NULL
    cdef accel *acc_ptr = NULL

    for i in range(<size_t>STAR.shape[0]):
        for p in range(num_components):
            signal = components[p]
            signal_phase_set = component_phases[p]
            phase_shift = phase_shifts[p]

            interp_ptr = interp[p]
            acc_ptr = acc[p]
            phases_ptr = &(signal_phase_set[0])
            signal_ptr = &(signal[i,0])

            gsl_interp_init(interp_ptr, phases_ptr, signal_ptr,
                            signal_phase_set.shape[0])

            for j in range(<size_t>phases.shape[0] - 1):
                a = phases[j] + phase_shift
                b = phases[j+1] + phase_shift

                if b - a == 1.0:
                    a = 0.0
                    b = 1.0
                else:
                    a -= floor(a)
                    b -= floor(b)

                if a < b:
                    _val = gsl_interp_eval_integ(interp_ptr,
                                                       phases_ptr,
                                                       signal_ptr,
                                                       a, b,
                                                       acc_ptr)
                    if _val > 0.0 or _allow_negative[p] == 1:
                        STAR[i,j] += _val
                else:
                    _val =  gsl_interp_eval_integ(interp_ptr,
                                                       phases_ptr,
                                                       signal_ptr,
                                                       a, 1.0,
                                                       acc_ptr)
                    if _val > 0.0 or _allow_negative[p] == 1:
                        STAR[i,j] += _val

                    _val = gsl_interp_eval_integ(interp_ptr,
                                                       phases_ptr,
                                                       signal_ptr,
                                                       0.0, b,
                                                       acc_ptr)
                    if _val > 0.0 or _allow_negative[p] == 1:
                        STAR[i,j] += _val

        for j in range(<size_t>phases.shape[0] - 1): # interpolant safety procedure
            if STAR[i,j] < 0.0:
                STAR[i,j] = 0.0

    for p in range(num_components):
        gsl_interp_accel_free(acc[p])
        gsl_interp_free(interp[p])

    free(acc)
    free(interp)

    cdef:
        double LOGLIKE = 0.0, EXPEC = 0.0
        double n = <double>(phases.shape[0] - 1)

    for i in range(<size_t>STAR.shape[0]):
        for j in range(<size_t>STAR.shape[1]):
            EXPEC = (STAR[i,j] + background[i,j]/n) * exposure_time

            # Ensuring that the log likelihood doesn't crash
            # even when model counts is zero in a given phase bin.
            if(EXPEC > 0.0):
                LOGLIKE -= EXPEC
                LOGLIKE += counts[i,j] * log(EXPEC)
    
            elif(counts[i,j] == 0 and EXPEC == 0):
                pass
            
            # Penalizing situations where the data counts are not zero 
            # but the model counts are.
            else:
                LOGLIKE += -1.0e90 * (0.1 + 0.9*np.random.rand(1))

            STAR[i,j] += background[i,j]/n
            STAR[i,j] *= exposure_time

    return (LOGLIKE, np.asarray(STAR, order='C', dtype=np.double))
