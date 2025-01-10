#cython: cdivision=True
#cython: boundscheck=False
#cython: nonecheck=False
#cython: wraparound=False
#cython: embedsignature=True

from __future__ import division, print_function
import numpy as np
cimport numpy as np

from libc.math cimport floor

from GSL cimport (gsl_interp_eval,
                   gsl_interp_alloc,
                   gsl_interp_accel,
                   gsl_interp_accel_alloc,
                   gsl_interp,
                   gsl_interp_init,
                   gsl_interp_free,
                   gsl_interp_accel_free,
                   gsl_interp_accel_reset)

from .core cimport _get_phase_interpolant, gsl_interp_type

def phase_interpolator(double[::1] new_phases,
                       double[::1] phases,
                       double[:,::1] signal,
                       double phase_shift,
                       bint allow_negative = 0):
    """ Interpolate a signal in phase.

    :param double[::1] new_phases:
        A :class:`numpy.ndarray` of phases in rotational cycles at which to
        interpolate.

    :param double[::1] signal_phases:
        A C-contiguous :class:`numpy.ndarray` of phases in cycles at which
        the model :obj:`signal` is evaluated on the interval ``[0,1]``.

    :param double[:,::1] signal:
        A C-contiguous :class:`numpy.ndarray` of signal count rates. Phase
        increases with column number.

    :param double phase_shift:
        A phase shift in cycles, such as on the interval ``[-0.5,0.5]``.

    :param obj allow_negative:
        A boolean declaring whether to allow negative phase interpolant
        integrals. If the interpolant is not a Steffen spline, then the
        interpolant of a non-negative function can be negative due to
        oscillations. For the default Akima Periodic spline from GSL, such
        oscillations should manifest as small relative to those present in
        cubic splines, for instance, because it is designed to handle a rapidly
        changing second-order derivative.

    :returns:
        A 2D :class:`numpy.ndarray` of the phase-shifted signal interpolated
        at the new set of phases. Phase increases with column number.

    """
    cdef const gsl_interp_type *_interpolant

    _interpolant = _get_phase_interpolant()

    cdef:
        size_t i, j
        double PHASE
        double[:,::1] new_signal = np.zeros((signal.shape[0],
                                            new_phases.shape[0]),
                                            dtype = np.double)

        gsl_interp_accel *accel = gsl_interp_accel_alloc()
        gsl_interp *interp = gsl_interp_alloc(_interpolant,
                                              phases.shape[0])

    cdef double *phase_ptr
    cdef double *signal_ptr
    cdef double _val

    for i in range(<size_t>signal.shape[0]):
        gsl_interp_accel_reset(accel)
        phase_ptr = &phases[0]
        signal_ptr = &signal[i][0]
        gsl_interp_init(interp, phase_ptr, signal_ptr, phases.shape[0])

        for j in range(<size_t>new_phases.shape[0]):
            PHASE = new_phases[j] + phase_shift
            PHASE -= floor(PHASE)

            _val =  gsl_interp_eval(interp, phase_ptr, signal_ptr,
                                             PHASE, accel)
            if _val > 0.0 or allow_negative:
                new_signal[i,j] = _val

    gsl_interp_free(interp)
    gsl_interp_accel_free(accel)

    return np.asarray(new_signal, dtype = np.double, order = 'C')
