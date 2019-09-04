#cython: cdivision=True
#cython: boundscheck=False
#cython: nonecheck=False
#cython: wraparound=False

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

def interpolate_pulse(double[::1] new_phases,
                      double[::1] phases,
                      double[:,::1] pulse,
                      double shift):

    cdef:
        size_t i, j
        double PHASE
        double[:,::1] new_pulse = np.zeros((pulse.shape[0],
                                            new_phases.shape[0]),
                                           dtype = np.double)

        gsl_interp_accel *accel = gsl_interp_accel_alloc()
        gsl_interp *interp = gsl_interp_alloc(gsl_interp_steffen, phases.shape[0])

    cdef double *phase_ptr
    cdef double *pulse_ptr

    for i in range(pulse.shape[0]):
        gsl_interp_accel_reset(accel)
        phase_ptr = &phases[0]
        pulse_ptr = &pulse[i][0]
        gsl_interp_init(interp, phase_ptr, pulse_ptr, phases.shape[0])

        for j in range(new_phases.shape[0]):
            PHASE = new_phases[j] - floor(new_phases[j]) + shift

            if PHASE > 1.0:
                PHASE -= 1.0
            elif PHASE < 0.0:
                PHASE += 1.0

            new_pulse[i,j] = gsl_interp_eval(interp, phase_ptr, pulse_ptr,
                                             PHASE, accel)

    gsl_interp_free(interp)
    gsl_interp_accel_free(accel)

    return np.asarray(new_pulse, dtype = np.double, order = 'C')
