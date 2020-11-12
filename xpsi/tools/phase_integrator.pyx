#cython: cdivision=True
#cython: boundscheck=False
#cython: nonecheck=False
#cython: wraparound=False
#cython: embedsignature=True

import numpy as np
from libc.math cimport floor

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

from . cimport _get_phase_interpolant, gsl_interp_type

def phase_integrator(double exposure_time,
                     double[::1] phases,
                     double[:,::1] signal,
                     double[::1] signal_phases,
                     double phase_shift,
                     bint allow_negative = 0):
    """ Integrate a signal over phase intervals.

    :param double exposure_time:
        Total exposure time in seconds. The exposure time scales the integrated
        signal.

    :param double[::1] phases:
        A :class:`numpy.ndarray` of phase interval edges in rotational cycles.

    :param double[:,::1] signal:
        A C-contiguous :class:`numpy.ndarray` of signal count rates. Phase
        increases with column number.

    :param double[::1] signal_phases:
        A C-contiguous :class:`numpy.ndarray` of phases in cycles at which
        the model :obj:`signal` is evaluated on the interval ``[0,1]``.

    :param double phase_shift:
        A phase shift in cycles such as on the interval ``[-0.5,0.5]``.

    :param obj allow_negative:
        A boolean declaring whether to allow negative phase interpolant
        integrals. If the interpolant is not a Steffen spline, then the
        interpolant of a non-negative function can be negative due to
        oscillations. For the default Akima Periodic spline from GSL, such
        oscillations should manifest as small relative to those present in
        cubic splines, for instance, because it is designed to handle a rapidly
        changing second-order derivative.

    :returns:
        A 2D :class:`numpy.ndarray` of the phase-shifted signal integrated
        over phase intervals and scaled by the exposure time. Phase interval
        number increases with column number.

    """
    cdef const gsl_interp_type *_interpolant

    _interpolant = _get_phase_interpolant()

    cdef:
        size_t i, j
        double a, b
        double _val

        gsl_interp *interp = gsl_interp_alloc(_interpolant, signal_phases.shape[0])
        gsl_interp_accel *acc = gsl_interp_accel_alloc()

        double[:,::1] _signal = np.zeros((signal.shape[0], phases.shape[0] - 1),
                                       dtype = np.double)

        double *phase_ptr
        double *signal_ptr
    for i in range(signal.shape[0]):
        gsl_interp_accel_reset(acc)
        phase_ptr = &(signal_phases[0])
        signal_ptr = &(signal[i,0])
        gsl_interp_init(interp, phase_ptr, signal_ptr, signal_phases.shape[0])

        for j in range(phases.shape[0] - 1):
            a = phases[j] + phase_shift
            b = phases[j+1] + phase_shift

            if b - a == 1.0: # for safety
                a = 0.0
                b = 1.0
            else:
                a -= floor(a)
                b -= floor(b)

            if a < b:
                _val = gsl_interp_eval_integ(interp, phase_ptr,
                                                     signal_ptr, a, b, acc)

                if _val > 0.0 or allow_negative:
                    _signal[i,j] = _val
            else:
                _val = gsl_interp_eval_integ(interp, phase_ptr,
                                                     signal_ptr, a, 1.0, acc)
                if _val > 0.0 or allow_negative:
                    _signal[i,j] = _val

                _val = gsl_interp_eval_integ(interp, phase_ptr,
                                                      signal_ptr, 0.0, b, acc)

                if _val > 0.0 or allow_negative:
                    _signal[i,j] += _val

            _signal[i,j] *= exposure_time

    gsl_interp_free(interp)
    gsl_interp_accel_free(acc)

    return np.asarray(_signal, dtype = np.double, order = 'C')
