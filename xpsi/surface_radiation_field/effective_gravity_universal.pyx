#cython: cdivision=True
#cython: boundscheck=False
#cython: nonecheck=False
#cython: wraparound=False

from libc.math cimport sqrt, log10, fabs

cdef double c = 2.99792458e8

cdef double effectiveGravity(double mu,
                             double R_eq,
                             double x,
                             double epsilon) noexcept nogil:
    """
    Calculate the effective surface gravity log-likelihood value based on the approximation from
    AlGendy & Morsink (2014) (see Eq. 21 and Table 5 for coefficient values).

    This function computes the logarithmic surface gravity (log10(g / g0)) for a rotating neutron star,
    incorporating the effects of spin and compactness. The result is expressed in units compatible with
    tabulated atmosphere models.

    :param double mu: 
        The cosine of the colatitude angle (theta). This parameter defines the angular
        position on the star's surface, where mu is between -1 and 1.

    :param double R_eq: 
        The equatorial radius of the star.

    :param double x: 
        The compactness parameter, defined as (G * M) / (R_eq * c**2) (see Eq. 1 or Eq. 9 from 
        Morsink et al. (2007)).

    :param double epsilon: 
        The dimensionless spin parameter, defined as (omega**2 * R_eq**3) / (G * M) (see Eq. 2 or 
        Eq. 10 from Morsink et al. (2007)).

    :return: 
        The effective surface gravities in log10(g*g0). The result is scaled to centimeters (conversion
        factor 2.0 is applied) to match the format used in tabulated atmosphere models.
    :rtype: double
    """

    cdef:
        double g_0 = x * c * c / (R_eq * sqrt(1.0 - 2.0 * x))
        double esq = epsilon
        double esqsq = epsilon * epsilon
        double c_e = -0.791 + 0.776 * x
        double c_p = 1.138 - 1.431 * x
        double d_e = (-1.315 + 2.431 * x) * esq * x
        double d_p = (0.653 - 2.864 * x) * esq * x
        double d_60 = (13.47 - 27.13 * x) * esq * x
        double f_e = -1.172 * x * esqsq
        double f_p = 0.975 * x * esqsq
        double g = 1.0

    g += (c_e + d_e + f_e) * esq * (1.0 - mu * mu)
    g += (c_p + d_p + f_p - d_60) * esq * mu * mu
    g += d_60 * esq * fabs(mu)

    return log10(g * g_0) + 2.0
