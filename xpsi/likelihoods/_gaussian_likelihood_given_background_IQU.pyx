#cython: cdivision=True
#cython: boundscheck=False
#cython: nonecheck=False
#cython: wraparound=False
#cython: embedsignature=True

from __future__ import division

import numpy as np
cimport numpy as np
from libc.math cimport pow, log, pi

from ..tools.compute_expected_counts cimport compute_expected_counts

def gaussian_likelihood_given_background(double exposure_time,
                                        double[::1] phases,
                                        double[:,::1] counts,
                                        double[:,::1] errors,
                                        components,
                                        component_phases,
                                        double[::1] phase_shifts,
                                        double[:,::1] background,
                                        allow_negative = False):
    """ Evaluate the Gaussian likelihood and compute the expected star + background counts.

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
        double LOGLIKE = 0.0, sigma_tot2 = 1.0, norm
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

    # Compute the gaussian likelihood
    for i in range(<size_t> EXPEC.shape[0]):
        for j in range(<size_t> EXPEC.shape[1]):
            sigma_tot2 = pow(errors[i,j],2.0)
            norm = 0.5 * log(2.0*pi*sigma_tot2)
            LOGLIKE -= (EXPEC[i,j]-counts[i,j])**2/(2.0*sigma_tot2)+norm
    return (LOGLIKE, np.asarray(EXPEC, order='C', dtype=np.double))