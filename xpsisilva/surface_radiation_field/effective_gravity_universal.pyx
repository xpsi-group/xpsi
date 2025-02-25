#cython: cdivision=True
#cython: boundscheck=False
#cython: nonecheck=False
#cython: wraparound=False

from libc.math cimport sqrt, log10, pow 

cdef double c = 2.99792458e8

cdef double effectiveGravity(double mu, 
                             double R_eq, 
                             double x, 
                             double epsilon) nogil:  
    
    """
    Calculate the effective surface gravity loglikelihood value based on the slow-elliptical approximation 
    from Silva et al. (2021) (see equation 17). This approximation is obtained using stars with epsilon <= 0.25.   

    :param mu: Colatitude (cos(theta))
    :type mu: double

    :param R_eq: Equatorial radius of the star 
    :type R_eq: double

    :param x: Compactness = (G*M)/(R_eq*c**2) defined as zeta in AlGendy & Morsink (2007) and as 
        kappa in Silva et al. (2021) 
    :type x: double 

    :param epsilon: Dimensionless spin parameter = (omega**2*R_eq**3)/(G*M) defined as sigma in Silva et al. (2021).
    :type epsilon: double

    :return: surface gravities in log10(g*g0) in units of centimeters (note that 2.0 is to convert from meters to 
        centimeters). This is used in the format used in the tabulated atmosphere tables. 
    :rtype: double 
    """

    cdef:
        double g_0 = x * c * c / (R_eq * sqrt(1.0 - 2.0 * x))
        double esq = epsilon 
        double esqsq = epsilon * epsilon 

        # slow-elliptical fit (epsilon<=0.25)
        double e = 1.089 * sqrt(esq) + 0.168 * esq - 0.685 * esq * x - 0.802 * esqsq
        double a2 = -1.013 - 0.312 *  esq + 0.930 * esq * x - 1.596 * esqsq 
        double a4 = 0.016 + 0.301 *  esq - 1.261 * esq * x + 2.728 * esqsq 

        # fast-elliptical fit (epsilon>0.25, corresponds to min ~700-800 Hz) 
        # double e = 0.251 + 0.935 * esq + 0.709 * x + 0.030 * esq * x - 0.472 * esqsq - 2.427 * x * x
        # double a2 = -1.265 + 0.220 * esq + 2.651 * x + 1.010 * esq * x - 1.815 * esqsq - 7.657 * x * x
        # double a4 = 0.556 - 1.465 * esq - 4.260 * x - 2.327 * esq * x + 4.921 * esqsq + 12.98 * x * x

        double g = 1 + a2 * mu**2 + a4 * pow(mu, 4) - (1 + a2 + a4) * pow(mu, 6)

    return log10(g * g_0) + 2.0
