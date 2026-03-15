#cython: cdivision=True
#cython: boundscheck=False
#cython: nonecheck=False
#cython: wraparound=False
#cython: embedsignature=True

from __future__ import print_function

import numpy as np
cimport numpy as np
from libc.stdlib cimport malloc, free
from libc.math cimport exp, pow, log, sqrt, fabs, floor
import time

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
                   gsl_rng,
                   gsl_rng_set,
                   gsl_rng_alloc,
                   gsl_rng_type,
                   gsl_rng_env_setup,
                   gsl_rng_free,
                   gsl_rng_default,
                   gsl_ran_poisson)

ctypedef gsl_interp_accel accel

ctypedef np.uint8_t uint8

from .compute_expected_counts cimport compute_expected_counts

def synthesise_exposure(double exposure_time,
                        double[::1] phases,
                        components,
                        component_phases,
                        phase_shifts,
                        double expected_background_counts,
                        double[:,::1] background,
                        allow_negative = False,
                        gsl_seed=None):
    """ Synthesise Poissonian count numbers given an exposure time.

    :param double exposure_time:
        Exposure time in seconds by which to scale the expected count rate

    :param double[::1] phases:
        A :class:`numpy.ndarray` of phase interval edges in cycles.

    :param tuple components:
        Component signals, each a C-contiguous :class:`numpy.ndarray` of
        signal count rates where phase increases with column number.

    :param tuple component_phases:
        For each component, a C-contiguous :class:`numpy.ndarray` of phases
        in cycles at which the model :obj:`signal` is evaluated on
        the interval ``[0,1]``. Typically this is more finely spaced than the
        resulting synthesised data.

    :param array-like phase_shift:
        Phase shifts in cycles, such as on the interval ``[-0.5,0.5]``, for
        the component signals.

    :param double expected_background_counts:
        The total expected number of background counts to set the background
        normalisation (given the exposure time).

    :param double[:,::1] background:
        A C-contiguous :class:`numpy.ndarray` of background expected
        *counts*, whose shape matches the number of channels in each element
        of :obj:`components` and the number of phase intervals (i.e. one element
        fewer than phase interval edges) constructed from :obj:`phases`.

    :param int gsl_seed:
        Seed number for adding random noise to the data. If not specified,
        seed is based on the clock time.

    :returns:
        A tuple ``(2D ndarray, 2D ndarray, double)``. The first element is
        the expected count numbers in joint phase-channel intervals. Note those
        are intervals, not edges. The second element is a stochastic realisation
        of those count numbers. The last element is the required normalisation of
        the background.
    """

    # Prepare variables
    cdef:
        double BACKGROUND, SCALE_BACKGROUND
        double[:,::1] EXPEC, SYNTHETIC, background_count_rate
        size_t i,j

    # Scale background
    for i in range(<size_t>components[0].shape[0] ):
        for j in range(<size_t>phases.shape[0] - 1):
            BACKGROUND += background[i,j]

    if BACKGROUND == 0.0: # allow zero background
        SCALE_BACKGROUND = 0.0
    else:
        SCALE_BACKGROUND = expected_background_counts / BACKGROUND

    background_count_rate = np.zeros((components[0].shape[0], phases.shape[0]-1), dtype = np.double)
    for i in range(<size_t>components[0].shape[0] ):
        for j in range(<size_t>phases.shape[0] - 1):
            background_count_rate[i,j] = background[i,j] * SCALE_BACKGROUND / exposure_time

    # Compute the expected counts from star + background
    EXPEC = compute_expected_counts(exposure_time,
                                    phases,
                                    components,
                                    component_phases,
                                    phase_shifts,
                                    background_count_rate,
                                    allow_negative)

    # Compute the Poisson realization
    cdef:
        const gsl_rng_type *T
        gsl_rng *r

    # Setup random number generator
    gsl_rng_env_setup()
    T = gsl_rng_default
    r = gsl_rng_alloc(T)
    if gsl_seed is not None:
        gsl_rng_set(r, gsl_seed)
    else:
        gsl_rng_set(r, time.time())

    # Do a Poisson realization     
    SYNTHETIC = np.zeros((components[0].shape[0], phases.shape[0]-1), dtype = np.double)
    for i in range(<size_t>EXPEC.shape[0]):
        for j in range(<size_t>EXPEC.shape[1]):
            SYNTHETIC[i,j] = gsl_ran_poisson(r, EXPEC[i,j])

    # Free memory
    gsl_rng_free(r)

    # Return result  
    return (np.asarray(EXPEC, order='C', dtype=np.double),
            np.asarray(SYNTHETIC, order='C', dtype=np.double),
            SCALE_BACKGROUND)


def synthesise_given_total_count_number(double[::1] phases,
                                        double expected_star_counts,
                                        components,
                                        component_phases,
                                        phase_shifts,
                                        double expected_background_counts,
                                        double[:,::1] background,
                                        allow_negative = False,
                                        gsl_seed=None):
    """ Synthesise Poissonian count numbers given expected target source counts.

    :param double[::1] phases:
        A :class:`numpy.ndarray` of phase interval edges in cycles.

    :param float expected_star_counts:
        Total number of expected counts from the star (the target source) to
        require.

    :param float expected_background_stars:
        Total number of expected background counts to require.

    :param tuple components:
        Component signals, each a C-contiguous :class:`numpy.ndarray` of
        signal count rates where phase increases with column number.

    :param tuple component_phases:
        For each component, a C-contiguous :class:`numpy.ndarray` of phases
        in cycles at which the model :obj:`signal` is evaluated on
        the interval ``[0,1]``. Typically this is more finely spaced than the
        resulting synthesised data.

    :param array-like phase_shift:
        Phase shifts in cycles, such as on the interval ``[-0.5,0.5]``, for
        the component signals.

    :param double expected_background_counts:
        The total expected number of background counts to set the background
        normalisation (given the exposure time).

    :param double[:,::1] background:
        A C-contiguous :class:`numpy.ndarray` of background expected
        *counts*, whose shape matches the number of channels in each element
        of :obj:`components` and the number of phase intervals (i.e. one element
        fewer than phase interval edges) constructed from :obj:`phases`.

    :param obj allow_negative:
        A boolean or an array of booleans, one per component, declaring whether
        to allow negative phase interpolant integrals. If the interpolant is
        not a Steffen spline, then the interpolant of a non-negative function
        can be negative due to oscillations. For the default Akima Periodic
        spline from GSL, such oscillations should manifest as small relative
        to those present in cubic splines, for instance, because it is
        designed to handle a rapidly changing second-order derivative.

    :param int gsl_seed:
        Seed number for adding random noise to the data. If not specified,
        seed is based on the clock time.

    :returns:
        A tuple ``(2D ndarray, 2D ndarray, double, double)``. The first element
        is the expected count numbers in joint phase-channel intervals. Note those
        are intervals, not edges. The second element is a stochastic realisation
        of those count numbers. The third element is the required exposure time.
        The last element is the required normalisation of the background.
    """

    # Prepare variables
    cdef:
        double BACKGROUND, SCALE_BACKGROUND, SCALE_STAR, STAR_total
        double[:,::1] EXPEC, SYNTHETIC
        size_t i,j

    # Compute the expected counts from star (with 0 background)
    STAR = compute_expected_counts(1e5,
                                    phases,
                                    components,
                                    component_phases,
                                    phase_shifts,
                                    np.zeros((components[0].shape[0], phases.shape[0]-1), dtype = np.double),
                                    allow_negative)

    # Get the total star and background count rates for normalization
    for i in range(<size_t>components[0].shape[0] ):
        for j in range(<size_t>phases.shape[0] - 1):
            STAR_total += STAR[i,j]
            BACKGROUND += background[i,j]
    SCALE_STAR = expected_star_counts / STAR_total

    if BACKGROUND == 0.0: # allow zero background
        SCALE_BACKGROUND = 0.0
    else:
        SCALE_BACKGROUND = expected_background_counts / BACKGROUND

    # Rescale and add up
    EXPEC = np.zeros((components[0].shape[0], phases.shape[0]-1), dtype = np.double)
    for i in range(<size_t>components[0].shape[0] ):
        for j in range(<size_t>phases.shape[0] - 1):
            STAR[i,j] = STAR[i,j] * SCALE_STAR
            EXPEC[i,j] = STAR[i,j] + background[i,j] * SCALE_BACKGROUND

    # Compute the Poisson realization
    cdef:
        const gsl_rng_type *T
        gsl_rng *r

    # Setup random number generator
    gsl_rng_env_setup()
    T = gsl_rng_default
    r = gsl_rng_alloc(T)
    if gsl_seed is not None:
        gsl_rng_set(r, gsl_seed)
    else:
        gsl_rng_set(r, time.time())

    # Do a Poisson realization     
    SYNTHETIC = np.zeros((components[0].shape[0], phases.shape[0]-1), dtype = np.double)
    for i in range(<size_t>EXPEC.shape[0]):
        for j in range(<size_t>EXPEC.shape[1]):
            SYNTHETIC[i,j] = gsl_ran_poisson(r, EXPEC[i,j])

    gsl_rng_free(r)

    return (np.asarray(EXPEC, order='C', dtype=np.double),
            np.asarray(SYNTHETIC, order='C', dtype=np.double),
            SCALE_STAR,
            SCALE_BACKGROUND/SCALE_STAR)