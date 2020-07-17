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
                   gsl_interp_steffen,
                   gsl_interp_init,
                   gsl_interp_free,
                   gsl_interp_accel_free,
                   gsl_interp_accel_reset)

def phase_interpolator(double[::1] new_phases,
                       double[::1] phases,
                       double[:,::1] signal,
                       double phase_shift):
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

    :returns:
        A 2D :class:`numpy.ndarray` of the phase-shifted signal interpolated
        at the new set of phases. Phase increases with column number.

    """

    cdef:
        size_t i, j
        double PHASE
        double[:,::1] new_signal = np.zeros((signal.shape[0],
                                            new_phases.shape[0]),
                                            dtype = np.double)

        gsl_interp_accel *accel = gsl_interp_accel_alloc()
        gsl_interp *interp = gsl_interp_alloc(gsl_interp_steffen, phases.shape[0])

    cdef double *phase_ptr
    cdef double *signal_ptr

    for i in range(signal.shape[0]):
        gsl_interp_accel_reset(accel)
        phase_ptr = &phases[0]
        signal_ptr = &signal[i][0]
        gsl_interp_init(interp, phase_ptr, signal_ptr, phases.shape[0])

        for j in range(new_phases.shape[0]):
            PHASE = new_phases[j] + phase_shift
            PHASE -= floor(PHASE)

            new_signal[i,j] = gsl_interp_eval(interp, phase_ptr, signal_ptr,
                                             PHASE, accel)

    gsl_interp_free(interp)
    gsl_interp_accel_free(accel)

    return np.asarray(new_signal, dtype = np.double, order = 'C')
