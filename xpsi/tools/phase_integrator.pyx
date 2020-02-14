#cython: cdivision=True
#cython: boundscheck=False
#cython: nonecheck=False
#cython: wraparound=False

import numpy as np

from GSL cimport (gsl_interp,
                   gsl_interp_alloc,
                   gsl_interp_init,
                   gsl_interp_free,
                   gsl_interp_eval,
                   gsl_interp_eval_integ,
                   gsl_interp_steffen,
                   gsl_interp_accel,
                   gsl_interp_accel_alloc,
                   gsl_interp_accel_free,
                   gsl_interp_accel_reset)

def phase_integrator(double T_exp,
                     double[::1] phases,
                     double[:,::1] pulse,
                     double[::1] pulse_phases,
                     double phase_shift):
    """ Integrate pulsation over phase intervals.

    :param float T_exp: Total exposure time in seconds.
    :param phases: A C-contiguous :class:`numpy.ndarray` of phase interval edges.
    :param pulse: A C-contiguous :class:`numpy.ndarray` of pulse count rates.
    :param pulse_phases: A C-contiguous :class:`numpy.ndarray` of phases at which
                         the model :obj:`pulse` is evaluated on the interval
                         ``[0,1]``.
    :param double phase_shift: A phase shift on the interval ``[-0.5,0.5]``.

    """

    cdef:
        size_t i, j
        double a, b

        gsl_interp *interp = gsl_interp_alloc(gsl_interp_steffen, pulse_phases.shape[0])
        gsl_interp_accel *acc = gsl_interp_accel_alloc()

        double[:,::1] PULSE = np.zeros((pulse.shape[0], phases.shape[0] - 1),
                                       dtype = np.double)

        double *phase_ptr
        double *pulse_ptr

    for i in range(pulse.shape[0]):
        gsl_interp_accel_reset(acc)
        phase_ptr = &(pulse_phases[0])
        pulse_ptr = &(pulse[i,0])
        gsl_interp_init(interp, phase_ptr, pulse_ptr, pulse_phases.shape[0])

        for j in range(phases.shape[0] - 1):
            a = phases[j] + phase_shift
            b = phases[j+1] + phase_shift

            if a > 1.0:
                while a > 1.0:
                    a -= 1.0
            elif a < 0.0:
                while a < 0.0:
                    a += 1.0

            if b > 1.0:
                while b > 1.0:
                    b -= 1.0
            elif b < 0.0:
                while b < 0.0:
                    b += 1.0

            if a < b:
                PULSE[i,j] = gsl_interp_eval_integ(interp, phase_ptr,
                                                   pulse_ptr, a, b, acc)
            else:
                PULSE[i,j] = gsl_interp_eval_integ(interp, phase_ptr,
                                                   pulse_ptr, a, 1.0, acc)
                PULSE[i,j] += gsl_interp_eval_integ(interp, phase_ptr,
                                                    pulse_ptr, 0.0, b, acc)

            PULSE[i,j] *= T_exp

    gsl_interp_free(interp)
    gsl_interp_accel_free(acc)

    return np.asarray(PULSE, dtype = np.double, order = 'C')
