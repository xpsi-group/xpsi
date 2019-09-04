#cython: cdivision=True
#cython: boundscheck=False
#cython: nonecheck=False
#cython: wraparound=False

from __future__ import division, print_function

import numpy as np
from cython.parallel cimport *
from libc.math cimport pow, log, floor
from libc.stdio cimport printf
from libc.stdlib cimport calloc, malloc, free

from GSL cimport (gsl_spline,
                   gsl_spline_alloc,
                   gsl_spline_init,
                   gsl_spline_free,
                   gsl_spline_eval,
                   gsl_spline_eval_integ,
                   gsl_interp_steffen,
                   gsl_interp_accel,
                   gsl_interp_accel_alloc,
                   gsl_interp_accel_free,
                   gsl_interp_accel_reset,
                   gsl_rng,
                   gsl_rng_alloc,
                   gsl_rng_type,
                   gsl_rng_env_setup,
                   gsl_rng_free,
                   gsl_rng_default,
                   gsl_ran_poisson)

ctypedef gsl_interp_accel accel

def eval_loglike(double T_exp,
                 double[::1] phases,
                 double[:,::1] counts,
                 pulses,
                 double[::1] pulse_phases,
                 phase_shifts,
                 double[:,::1] background):
    """ Evaluate the Poisson likelihood. """

    cdef:
        size_t i, j, p, num_pulses = len(pulses)

        double LOGLIKE = 0.0, EXPEC = 0.0
        double PHASE, a, b, phase_shift

        double n = <double>(phases.shape[0] - 1)

        double[:,::1] pulse

        gsl_spline *spline = gsl_spline_alloc(gsl_interp_steffen, pulse_phases.shape[0])
        accel *acc = gsl_interp_accel_alloc()

        double[:,::1] PULSE = np.zeros((pulses[0].shape[0], phases.shape[0]-1),
                                       dtype = np.double)

    for i in range(PULSE.shape[0]):
        for p in range(num_pulses):
            pulse = pulses[p]
            phase_shift = phase_shifts[p]

            gsl_interp_accel_reset(acc)
            gsl_spline_init(spline, &(pulse_phases[0]), &(pulse[i,0]), pulse_phases.shape[0])

            for j in range(phases.shape[0] - 1):
                a = phases[j] + phase_shift
                b = phases[j+1] + phase_shift

                if a > 1.0:
                    a -= 1.0
                elif a < 0.0:
                    a += 1.0

                if b > 1.0:
                    b -= 1.0
                elif b < 0.0:
                    b += 1.0

                if a < b:
                    PULSE[i,j] += gsl_spline_eval_integ(spline, a, b, acc)
                else:
                    PULSE[i,j] += gsl_spline_eval_integ(spline, a, 1.0, acc)
                    PULSE[i,j] += gsl_spline_eval_integ(spline, 0.0, b, acc)

    gsl_spline_free(spline)
    gsl_interp_accel_free(acc)

    for i in range(PULSE.shape[0]):
        for j in range(PULSE.shape[1]):
            EXPEC = (PULSE[i,j] + background[i,j]/n) * T_exp
            LOGLIKE -= EXPEC
            LOGLIKE += counts[i,j] * log(EXPEC)

            PULSE[i,j] += background[i,j]/n
            PULSE[i,j] *= T_exp

    return (LOGLIKE, np.asarray(PULSE, order='C', dtype=np.double))







