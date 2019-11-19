#cython: cdivision=True
#cython: boundscheck=False
#cython: nonecheck=False
#cython: wraparound=False

from __future__ import print_function

import numpy as np
from libc.math cimport exp, pow, log, sqrt, fabs

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

def synthesise_exposure(double[::1] phases,
                        double exposure_time,
                        double expected_background_counts,
                        pulses,
                        double[::1] pulse_phases,
                        double[:,::1] background,
                        phase_shifts):
    """ Synthesise from Poisson generative model given an exposure time.

    :param phases:
        A C-contiguous :class:`numpy.ndarray` of phase interval edges.
    :param float exposure_time:
        Exposure time in seconds.
    :param pulse:
        A C-contiguous :class:`numpy.ndarray` of pulse count rates.
    :param pulse_phases:
        A C-contiguous :class:`numpy.ndarray` of phases at which
        the model :obj:`pulse` is evaluated on the interval ``[0,1]``.
    :param background:
        A C-contiguous :class:`numpy.ndarray` of background expected
        counts, whose shape matches :obj:`counts`.
    :param array-like phase_shift:
        Phase shifts in cycles.

    """

    cdef:
        size_t i, j, p, num_pulses = len(pulses)
        double BACKGROUND, a, b

        gsl_spline *spline = gsl_spline_alloc(gsl_interp_steffen, pulse_phases.shape[0])
        accel *acc = gsl_interp_accel_alloc()

        double[:,::1] PULSE = np.zeros((pulses[0].shape[0], phases.shape[0]-1),
                                      dtype = np.double)
        double[:,::1] SYNTHETIC = np.zeros((pulses[0].shape[0], phases.shape[0]-1),
                                           dtype = np.double)

    cdef double PHASE, SCALE_BACKGROUND
    cdef double[:,::1] pulse
    cdef double phase_shift

    BACKGROUND = 0.0

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

        for j in range(phases.shape[0] - 1):
            BACKGROUND += background[i,j]

    gsl_spline_free(spline)
    gsl_interp_accel_free(acc)

    SCALE_BACKGROUND = expected_background_counts / BACKGROUND

    cdef:
        const gsl_rng_type *T
        gsl_rng *r

    gsl_rng_env_setup()

    T = gsl_rng_default
    r = gsl_rng_alloc(T)

    for i in range(PULSE.shape[0]):
        for j in range(PULSE.shape[1]):
            PULSE[i,j] *= exposure_time
            background[i,j] *= SCALE_BACKGROUND

            PULSE[i,j] += background[i,j]

            SYNTHETIC[i,j] = gsl_ran_poisson(r, PULSE[i,j])

    gsl_rng_free(r)

    return (np.asarray(PULSE, order='C', dtype=np.double),
            np.asarray(SYNTHETIC, order='C', dtype=np.double),
            SCALE_BACKGROUND)
            #np.asarray(background, order='C', dtype=np.double))

def synthesise(double[::1] phases,
               double expected_star_counts,
               double expected_background_counts,
               pulses,
               double[::1] pulse_phases,
               double[:,::1] background,
               phase_shifts):
    """ Synthesise from Poisson generative model.

    :param phases:
        A C-contiguous :class:`numpy.ndarray` of phase interval edges.
    :param float expected_star_counts:
        Total number of expected counts from the star to require.
    :param float expected_background_stars:
        Total number of expected background counts to require.
    :param tuple pulses:
        A tuple of C-contiguous :class:`numpy.ndarray` of pulse count rates.
    :param pulse_phases:
        A C-contiguous :class:`numpy.ndarray` of phases at which
        the model :obj:`pulse` is evaluated on the interval ``[0,1]``.
    :param background:
        A C-contiguous :class:`numpy.ndarray` of background expected
        counts, whose shape matches :obj:`counts`.
    :param array-like phase_shifts: Phase shifts in cycles.

    """

    cdef:
        size_t i, j, p, num_pulses = len(pulses)
        double STAR, BACKGROUND, a, b

        gsl_spline *spline = gsl_spline_alloc(gsl_interp_steffen, pulse_phases.shape[0])
        accel *acc = gsl_interp_accel_alloc()

        double[:,::1] PULSE = np.zeros((pulses[0].shape[0], phases.shape[0]-1),
                                      dtype = np.double)
        double[:,::1] SYNTHETIC = np.zeros((pulses[0].shape[0], phases.shape[0]-1),
                                           dtype = np.double)

    cdef double PHASE, SCALE_STAR, SCALE_BACKGROUND
    cdef double[:,::1] pulse
    cdef double phase_shift

    STAR = 0.0
    BACKGROUND = 0.0

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

        for j in range(phases.shape[0] - 1):
            STAR += PULSE[i,j]
            BACKGROUND += background[i,j]

    gsl_spline_free(spline)
    gsl_interp_accel_free(acc)

    SCALE_STAR = expected_star_counts / STAR
    SCALE_BACKGROUND = expected_background_counts / BACKGROUND

    print("Exposure time: %.6f [s]" % SCALE_STAR)
    print("Background normalisation: %.8e" % (SCALE_BACKGROUND/SCALE_STAR))

    cdef:
        const gsl_rng_type *T
        gsl_rng *r

    gsl_rng_env_setup()

    T = gsl_rng_default
    r = gsl_rng_alloc(T)

    for i in range(pulse.shape[0]):
        for j in range(phases.shape[0] - 1):
            PULSE[i,j] *= SCALE_STAR
            background[i,j] *= SCALE_BACKGROUND

            PULSE[i,j] += background[i,j]

            SYNTHETIC[i,j] = gsl_ran_poisson(r, PULSE[i,j])

    gsl_rng_free(r)

    return (np.asarray(PULSE, order='C', dtype=np.double),
            np.asarray(SYNTHETIC, order='C', dtype=np.double))

