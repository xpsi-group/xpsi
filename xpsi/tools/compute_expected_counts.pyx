#cython: cdivision=True
#cython: boundscheck=False
#cython: nonecheck=False
#cython: wraparound=False
#cython: embedsignature=True

import numpy as np
cimport numpy as np
from libc.math cimport floor
from libc.stdlib cimport malloc, free

from GSL cimport (gsl_interp,
                   gsl_interp_alloc,
                   gsl_interp_init,
                   gsl_interp_free,
                   gsl_interp_eval_integ,
                   gsl_interp_accel,
                   gsl_interp_accel_alloc,
                   gsl_interp_accel_free,
                   gsl_interp_accel_reset)

ctypedef gsl_interp_accel accel

from ..tools.core cimport _get_phase_interpolant, gsl_interp_type, are_equal

cdef uint8[::1] build_allow_negative_array( object allow_negative , 
                                            size_t num_components):
    """ Build the array of allowed negative array to be compliant with different input types
    
    :param object allow_negative:
        A boolean or an array of booleans, one per component, declaring whether
        to allow negative phase interpolant integrals. If the interpolant is
        not a Steffen spline, then the interpolant of a non-negative function
        can be negative due to oscillations. For the default Akima Periodic
        spline from GSL, such oscillations should manifest as small relative
        to those present in cubic splines, for instance, because it is
        designed to handle a rapidly changing second-order derivative.
    
    :param size_t num_components:
        Integer of the number of different components being included in the star model.
    """
    
    cdef uint8[::1] _allow_negative = np.zeros(num_components, dtype=np.uint8)
    cdef size_t i

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
            if <size_t> len(allow_negative) != num_components:
                raise ValueError('Number of allow_negative declarations does '
                                 'not match the number of components..')

            for i in range(num_components):
                _allow_negative[i] = allow_negative[i]

    return _allow_negative


cdef void compute_expected_star_count_rate_single_channel(double[::1] phases,
                                         size_t num_components,
                                         object components,
                                         object component_phases,
                                         double[::1] phase_shifts,
                                         gsl_interp **interp,
                                         gsl_interp_accel **acc,
                                         double *STAR,
                                         int channel,
                                         uint8[::1] allow_negative):

    """ Compute the expected signal from the star in cts/s for a single channel.
    
    :param double[::1] phases:
        A :class:`numpy.ndarray` of phase interval edges in cycles.

    :param size_t num_components:
        Number of components to loop over.
    
    :param tuple components:
        Component signals, each a C-contiguous :class:`numpy.ndarray` of
        signal count rates where phase increases with column number.

    :param tuple component_phases:
        For each component, a C-contiguous :class:`numpy.ndarray` of phases
        in cycles at which the model :obj:`signal` is evaluated on
        the interval ``[0,1]``.

    :param double[::1] phase_shifts:
        Phase shifts in cycles, such as on the interval ``[-0.5,0.5]``, for
        the component signals.

    :param gsl_interp **interp:
        Pointer to the gsl interpolator to use.

    :param gsl_interp_accel **acc:
        Pointer to the accelerator of the gsl interpolator to use.
 
    :param double *STAR:
        Pointer to the STAR array to update with the expected star count rate.
        Points to the right channel so only needs 1D input on phase.

    :param int channel:
        Channel number to compute the count rate for.

    :param uint8[::1] allow_negative:
        An array of unint (0 or 1), one per component, declaring whether
        to allow negative phase interpolant integrals. 
     """

    # Define variables
    cdef:
        double[:,::1] pulse
        double[::1] pulse_phase_set
        double phase_shift
        double *phases_ptr = NULL
        double *pulse_ptr = NULL
        gsl_interp *inter_ptr = NULL
        accel *acc_ptr = NULL
        double pa, pb, _val
        size_t p, j

    # Loop over the different components and retrieve useful information
    for p in range(num_components):
        pulse = components[p]
        pulse_phase_set = component_phases[p]
        phase_shift = phase_shifts[p]

        # Do some integration if there is more than one pulse phase 
        if pulse.shape[1] > 1:

            # Prepare pointers for the integrator
            interp_ptr = interp[p]
            acc_ptr = acc[p]
            phases_ptr = &(pulse_phase_set[0])
            pulse_ptr = &(pulse[channel,0])

            gsl_interp_init(interp_ptr, phases_ptr, pulse_ptr,
                            pulse_phase_set.shape[0])

            # Loop through the instrument phase bins to perform integration 
            for j in range(<size_t> (phases.shape[0] - 1)):

                # Prepare edges of integration by shifting the phase
                pa = phases[j] + phase_shift
                pb = phases[j+1] + phase_shift

                if are_equal(pb - pa, 1.0):
                    pa = 0.0
                    pb = 1.0
                else:
                    pa -= floor(pa)
                    pb -= floor(pb)

                # Do the integration
                if pa < pb:

                    # Integrate the pulse 
                    _val = gsl_interp_eval_integ(interp_ptr,
                                                phases_ptr,
                                                pulse_ptr,
                                                pa, pb,
                                                acc_ptr)
                    if _val > 0.0 or allow_negative[p] == 1:
                        STAR[j] += _val

                # Boundary case when reaching the periodic phase=1 boundary
                else:
                    _val = gsl_interp_eval_integ(interp_ptr,
                                                    phases_ptr,
                                                    pulse_ptr,
                                                    pa, 1.0,
                                                    acc_ptr)
                    if _val > 0.0 or allow_negative[p] == 1:
                        STAR[j] += _val

                    _val = gsl_interp_eval_integ(interp_ptr,
                                                    phases_ptr,
                                                    pulse_ptr,
                                                    0.0, pb,
                                                    acc_ptr)
                    if _val > 0.0 or allow_negative[p] == 1:
                        STAR[j] += _val
        
        # Otherwise, just store (e.g. for Elsewhere)
        else:
            STAR[0] = pulse[channel,0] 

    # Safety to clip negative values to 0
    for j in range(<size_t> (phases.shape[0] - 1)): # interpolant safety procedure
        if STAR[j] < 0.0:
            STAR[j] = 0.0


cdef double[:,::1] compute_expected_counts(double exposure_time,
                            double[::1] phases,
                            object components,
                            object component_phases,
                            double[::1] phase_shifts,
                            double[:,::1] background,
                            object allow_negative):

    """ Compute the expected counts given various components and a background count rate.
    
    :param double exposure_time:
        Exposure time in seconds by which to scale the expected count rate
        in each phase interval.
    
    :param double[::1] phases:
        A :class:`numpy.ndarray` of phase interval edges in cycles.

    :param tuple components:
        Component signals, each a C-contiguous :class:`numpy.ndarray` of
        signal count rates where phase increases with column number.

    :param tuple component_phases:
        For each component, a C-contiguous :class:`numpy.ndarray` of phases
        in cycles at which the model :obj:`signal` is evaluated on
        the interval ``[0,1]``.

    :param array-like phase_shifts:
        Phase shifts in cycles, such as on the interval ``[-0.5,0.5]``, for
        the component signals.

    :param obj background:
        If not ``None``, then a C-contiguous :class:`numpy.ndarray` of
        background count rates where phase interval increases with column
        number. Useful for phase-dependent backgrounds, or a phase-independent
        background if the channel-by-channel background variable prior support
        is restricted.

    :param obj allow_negative:
        A boolean or an array of booleans, one per component, declaring whether
        to allow negative phase interpolant integrals. If the interpolant is
        not a Steffen spline, then the interpolant of a non-negative function
        can be negative due to oscillations. For the default Akima Periodic
        spline from GSL, such oscillations should manifest as small relative
        to those present in cubic splines, for instance, because it is
        designed to handle a rapidly changing second-order derivative.

    :return double 2Dndarray:
        The expected count numbers in joint phase-channel intervals from the star
        (the target source) added with background.

    """ 

    # Define the expected count array and number of component and phases
    cdef:
        size_t num_components = len(components)
        size_t i, j, p
        double n = <double>(phases.shape[0] - 1)
        double[:,::1] STAR = np.zeros((components[0].shape[0], phases.shape[0]-1),
                                       dtype = np.double)

    # Set the interpolant to use for the integration
    cdef const gsl_interp_type *_interpolant
    _interpolant = _get_phase_interpolant()

    # Define the allow_negative when using multiple components
    cdef uint8[::1] _allow_negative = build_allow_negative_array( allow_negative , num_components )

    # Allocate the space for interpolator and accelerator
    cdef gsl_interp **interp = <gsl_interp**> malloc(num_components * sizeof(gsl_interp*))
    cdef accel **acc =  <accel**> malloc(num_components * sizeof(accel*))

    for p in range(num_components):
        pulse_phase_set = component_phases[p]
        if components[p].shape[1] > 1:
            interp[p] = gsl_interp_alloc(_interpolant,
                                        pulse_phase_set.shape[0])
            acc[p] = gsl_interp_accel_alloc()
            gsl_interp_accel_reset(acc[p])

    cdef gsl_interp *inter_ptr = NULL
    cdef accel *acc_ptr = NULL

    # Finally we will perform the integration
    # Loop over the energy channels
    for i in range(<size_t> STAR.shape[0]):

        # Compute the expected count rate
        compute_expected_star_count_rate_single_channel(phases,
                                         num_components,
                                         components,
                                         component_phases,
                                         phase_shifts,
                                         interp,
                                         acc,
                                         &(STAR[i,0]),
                                         i,
                                         _allow_negative)

        # Add the background and transform to counts
        for j in range(<size_t> STAR.shape[1]):
            
            # Add background
            STAR[i,j] += background[i,j]/n

            # From count rate to counts
            STAR[i,j] *= exposure_time

    # Clean up the work space to avoid memory leakage
    for p in range(num_components):
        gsl_interp_accel_free(acc[p])
        gsl_interp_free(interp[p])
    free(acc)
    free(interp)

    # Return the computed expected counts
    return STAR
