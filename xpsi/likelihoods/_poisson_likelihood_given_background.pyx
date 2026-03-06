#cython: cdivision=True
#cython: boundscheck=False
#cython: nonecheck=False
#cython: wraparound=False
#cython: embedsignature=True

import numpy as np
cimport numpy as np
from libc.math cimport log

from ..tools.compute_expected_counts cimport compute_expected_counts


def poisson_likelihood_given_background(double exposure_time,
                                        double[::1] phases,
                                        double[:,::1] counts,
                                        components,
                                        component_phases,
                                        double[::1] phase_shifts,
                                        double[:,::1] background,
                                        double[::1] neg_sum_ln_data_factorial=None,
                                        allow_negative = False):
    """ Evaluate the Poisson likelihood and compute the expected star + background counts.

    The count rate is integrated over phase intervals.

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

    :param double[:,::1] background:
        A C-contiguous :class:`numpy.ndarray` of background expected
        *counts*, whose shape matches the number of channels in each element
        of :obj:`components` and the number of phase intervals constructed
        from :obj:`phases`.

    :param double[::1] neg_sum_ln_data_factorial:
        The precomputed output of :func:`~.precomputation` (found in
        default_background_marginalisation.pyx) given the data count numbers.

    :param obj allow_negative:
        A boolean or an array of booleans, one per component, declaring whether
        to allow negative phase interpolant integrals. If the interpolant is
        not a Steffen spline, then the interpolant of a non-negative function
        can be negative due to oscillations. For the default Akima Periodic
        spline from GSL, such oscillations should manifest as small relative
        to those present in cubic splines, for instance, because it is
        designed to handle a rapidly changing second-order derivative.

    :returns:
        A tuple ``(double, 2D ndarray)``. The first element is
        the logarithm of the marginal likelihood. The second element is the
        expected count numbers in joint phase-channel intervals from the star
        (the target source).

    """

    # Prepare variables
    cdef:
        double LOGLIKE = 0.0
        double[:,::1] EXPEC
        size_t i,j

    # Compute the expected counts from star + background
    EXPEC = compute_expected_counts(exposure_time,
                                    phases,
                                    components,
                                    component_phases,
                                    phase_shifts,
                                    background,
                                    allow_negative)

    # Compute the poisson likelihood
    # First loop over channels
    for i in range(<size_t>EXPEC.shape[0]):

        # Add the constant term if any
        if neg_sum_ln_data_factorial is not None:
            LOGLIKE += neg_sum_ln_data_factorial[i]

        # Loop over phases
        for j in range(<size_t>EXPEC.shape[1]):

            # Ensuring that the log likelihood doesn't crash
            # even when model counts is zero in a given phase bin.
            if(EXPEC[i,j] > 0.0):
                LOGLIKE -= EXPEC[i,j]
                LOGLIKE += counts[i,j] * log(EXPEC[i,j])
            elif(counts[i,j] == 0 and EXPEC[i,j] == 0):
                pass
            
            # Penalizing situations where the data counts are not zero 
            # but the model counts are.
            else:
                LOGLIKE += -1.0e90 * (0.1 + 0.9*np.random.rand(1))


    return (LOGLIKE, np.asarray(EXPEC, order='C', dtype=np.double))