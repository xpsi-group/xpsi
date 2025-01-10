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

from .core cimport _get_phase_interpolant, gsl_interp_type

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
    cdef const gsl_interp_type *_interpolant

    _interpolant = _get_phase_interpolant()

    cdef:
        size_t i, j, p, num_components = len(components)
        double BACKGROUND, a, b

        double[:,::1] STAR = np.zeros((components[0].shape[0], phases.shape[0]-1),
                                       dtype = np.double)
        double[:,::1] SYNTHETIC = np.zeros((components[0].shape[0], phases.shape[0]-1),
                                           dtype = np.double)

    cdef double *phases_ptr = NULL
    cdef double *signal_ptr = NULL

    cdef double[:,::1] signal
    cdef double[::1] signal_phase_set
    cdef double phase_shift
    cdef double _val
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

    cdef gsl_interp **interp = <gsl_interp**> malloc(num_components * sizeof(gsl_interp*))
    cdef accel **acc =  <accel**> malloc(num_components * sizeof(accel*))

    for p in range(num_components):
        signal_phase_set = component_phases[p]
        interp[p] = gsl_interp_alloc(_interpolant, signal_phase_set.shape[0])
        acc[p] = gsl_interp_accel_alloc()
        gsl_interp_accel_reset(acc[p])

    cdef gsl_interp *inter_ptr = NULL
    cdef accel *acc_ptr = NULL

    cdef double SCALE_BACKGROUND
    BACKGROUND = 0.0

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
                    _val = gsl_interp_eval_integ(interp_ptr,
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


        for j in range(<size_t>phases.shape[0] - 1):
            BACKGROUND += background[i,j]

    for p in range(num_components):
        gsl_interp_accel_free(acc[p])
        gsl_interp_free(interp[p])

    free(acc)
    free(interp)

    if BACKGROUND == 0.0: # allow zero background
        SCALE_BACKGROUND = 0.0
    else:
        SCALE_BACKGROUND = expected_background_counts / BACKGROUND

    cdef:
        const gsl_rng_type *T
        gsl_rng *r

    gsl_rng_env_setup()

    T = gsl_rng_default
    r = gsl_rng_alloc(T)

    if gsl_seed is not None:
        gsl_rng_set(r, gsl_seed);
    else:
        gsl_rng_set(r, time.time());

    for i in range(<size_t>STAR.shape[0]):
        for j in range(<size_t>STAR.shape[1]):
            STAR[i,j] *= exposure_time

            STAR[i,j] += background[i,j] * SCALE_BACKGROUND

            SYNTHETIC[i,j] = gsl_ran_poisson(r, STAR[i,j])

    gsl_rng_free(r)

    return (np.asarray(STAR, order='C', dtype=np.double),
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
    cdef const gsl_interp_type *_interpolant

    _interpolant = _get_phase_interpolant()

    cdef:
        size_t i, j, p, num_components = len(components)
        double STAR, BACKGROUND, a, b

        double[:,::1] _signal = np.zeros((components[0].shape[0], phases.shape[0]-1),
                                         dtype = np.double)
        double[:,::1] SYNTHETIC = np.zeros((components[0].shape[0], phases.shape[0]-1),
                                           dtype = np.double)

    cdef double *phases_ptr = NULL
    cdef double *signal_ptr = NULL

    cdef double[:,::1] signal
    cdef double[::1] signal_phase_set
    cdef double phase_shift
    cdef double _val
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

    cdef gsl_interp **interp = <gsl_interp**> malloc(num_components * sizeof(gsl_interp*))
    cdef accel **acc =  <accel**> malloc(num_components * sizeof(accel*))

    for p in range(num_components):
        signal_phase_set = component_phases[p]
        interp[p] = gsl_interp_alloc(_interpolant, signal_phase_set.shape[0])
        acc[p] = gsl_interp_accel_alloc()
        gsl_interp_accel_reset(acc[p])

    cdef gsl_interp *inter_ptr = NULL
    cdef accel *acc_ptr = NULL

    cdef double SCALE_STAR, SCALE_BACKGROUND
    STAR = 0.0
    BACKGROUND = 0.0

    for i in range(<size_t>_signal.shape[0]):
        for p in range(<size_t>num_components):
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

                a -= floor(a)
                b -= floor(b)

                if a < b:
                    _val = gsl_interp_eval_integ(interp_ptr,
                                                          phases_ptr,
                                                          signal_ptr,
                                                          a, b,
                                                          acc_ptr)
                    if _val > 0.0 or _allow_negative[p] == 1:
                        _signal[i,j] += _val
                else:
                    _val = gsl_interp_eval_integ(interp_ptr,
                                                          phases_ptr,
                                                          signal_ptr,
                                                          a, 1.0,
                                                          acc_ptr)
                    if _val > 0.0 or _allow_negative[p] == 1:
                        _signal[i,j] += _val

                    _val = gsl_interp_eval_integ(interp_ptr,
                                                          phases_ptr,
                                                          signal_ptr,
                                                          0.0, b,
                                                          acc_ptr)
                    if _val > 0.0 or _allow_negative[p] == 1:
                        _signal[i,j] += _val

        for j in range(<size_t>phases.shape[0] - 1): # interpolant safety procedure
            if _signal[i,j] < 0.0:
                _signal[i,j] = 0.0

        for j in range(<size_t>phases.shape[0] - 1):
            STAR += _signal[i,j]
            BACKGROUND += background[i,j]

    for p in range(num_components):
        gsl_interp_accel_free(acc[p])
        gsl_interp_free(interp[p])

    free(acc)
    free(interp)

    SCALE_STAR = expected_star_counts / STAR

    if BACKGROUND == 0.0: # allow zero background
        SCALE_BACKGROUND = 0.0
    else:
        SCALE_BACKGROUND = expected_background_counts / BACKGROUND

    print('Exposure time: %.6f [s]' % SCALE_STAR)
    print('Background normalisation: %.8e' % (SCALE_BACKGROUND/SCALE_STAR))

    cdef:
        const gsl_rng_type *T
        gsl_rng *r

    gsl_rng_env_setup()

    T = gsl_rng_default
    r = gsl_rng_alloc(T)

    if gsl_seed is not None:
        gsl_rng_set(r, gsl_seed);
    else:
        gsl_rng_set(r, time.time());

    for i in range(<size_t>_signal.shape[0]):
        for j in range(<size_t>phases.shape[0] - 1):
            _signal[i,j] *= SCALE_STAR
            _signal[i,j] += background[i,j] * SCALE_BACKGROUND

            SYNTHETIC[i,j] = gsl_ran_poisson(r, _signal[i,j])

    gsl_rng_free(r)

    return (np.asarray(_signal, order='C', dtype=np.double),
            np.asarray(SYNTHETIC, order='C', dtype=np.double),
            SCALE_STAR,
            SCALE_BACKGROUND/SCALE_STAR)
