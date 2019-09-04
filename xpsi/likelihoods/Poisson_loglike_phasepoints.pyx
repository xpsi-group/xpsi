#cython: cdivision=True
#cython: boundscheck=False
#cython: nonecheck=False
#cython: wraparound=False

from __future__ import print_function

import numpy as np
from libc.math cimport pow, log

from GSL cimport (gsl_spline,
                   gsl_spline_alloc,
                   gsl_spline_init,
                   gsl_spline_free,
                   gsl_spline_eval,
                   gsl_interp_steffen,
                   gsl_interp_accel,
                   gsl_interp_accel_alloc,
                   gsl_interp_accel_free,
                   gsl_interp_accel_reset,
                   gsl_isnan,
                   gsl_isinf)

ctypedef gsl_interp_accel accel

def eval_loglike(double exposure_time,
                 double[::1] phases,
                 double[:,::1] counts,
                 double[:,::1] pulse,
                 double[::1] pulse_phases,
                 double phase_shift):
    """ Evaluate the Poisson likelihood. """

    cdef:
        size_t i, j
        double LOGLIKE = 0.0

        gsl_spline *spline = gsl_spline_alloc(gsl_interp_steffen, pulse_phases.shape[0])
        accel *acc = gsl_interp_accel_alloc()

        double[:,::1] STAR = np.zeros((pulse.shape[0], phases.shape[0]),
                                      dtype = np.double)

    cdef double PHASE, SCALE = exposure_time / <double>(phases.shape[0])

    cdef double total_star = 0.0

    for i in range(pulse.shape[0]):

        gsl_interp_accel_reset(acc)
        gsl_spline_init(spline, &(pulse_phases[0]), &(pulse[i,0]), pulse_phases.shape[0])

        for j in range(phases.shape[0]):
            PHASE = phases[j] + phase_shift

            if PHASE > 1.0:
                PHASE = PHASE - 1.0
            elif PHASE < 0.0:
                PHASE = 1.0 + PHASE

            STAR[i,j] = gsl_spline_eval(spline, PHASE, acc)

            total_star += STAR[i,j]

    gsl_interp_accel_free(acc)
    gsl_spline_free(spline)

    cdef:
        double B_prime, sum_data, sum_star, BppS, s1, s2
        double n = <double>(phases.shape[0])

    for i in range(STAR.shape[0]):
        sum_data = 0.0; sum_star = 0.0
        for j in range(STAR.shape[1]):
            sum_data += counts[i,j]
            sum_star += STAR[i,j]

        B_prime = sum_data / exposure_time - sum_star / n
        print('B_prime = %.6f' % B_prime)

        if B_prime < 0.0:
            B_prime = 0.0

        s1 = 0.0; s2 = 0.0
        for j in range(STAR.shape[1]):
            BppS = B_prime + STAR[i,j]
            LOGLIKE += counts[i,j] * log(SCALE * BppS)
            LOGLIKE -= SCALE * BppS

            s1 += counts[i,j]/BppS
            s2 += counts[i,j]/pow(BppS,2.0)

        LOGLIKE += 0.5 * pow(exposure_time - s1, 2.0) / s2
        LOGLIKE -= 0.5 * log(s2)

    print('10000/total_star = %.6f' % (10000.0/total_star))

    return LOGLIKE

