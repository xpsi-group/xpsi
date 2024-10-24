#cython: cdivision=True
#cython: boundscheck=False
#cython: nonecheck=False
#cython: wraparound=False

from libc.math cimport sqrt, log10, pow # fabs

cdef double c = 2.99792458e8

cdef double effectiveGravity(double mu, # cos(colatitude: cos(theta)
                             double R_eq, # equatorial radius 
                             double x, # compactness: zeta (in AlGendy) = kappa (in Silva) = (G*M)/(R_eq*c**2)
                             double epsilon) nogil: # dimensionless spin parameter: sigma (in Silva) = (omega**2*R_eq**3)/(G*M)

    cdef:
        double g_0 = x * c * c / (R_eq * sqrt(1.0 - 2.0 * x))
        double esq = epsilon # why do this? 
        double esqsq = epsilon * epsilon # why seperately define this? 

        # slow-elliptical fit (epsilon<=0.25)
        double e = 1.089 * sqrt(esq) + 0.168 * esq - 0.685 * esq * x - 0.802 * esqsq
        double a2 = -1.013 - 0.312 *  esq + 0.930 * esq * x - 1.596 * esqsq 
        double a4 = 0.016 + 0.301 *  esq - 1.261 * esq * x + 2.728 * esqsq 

        # fast-elliptical fit (epsilon>0.25, corresponds to min ~700-800 Hz) 
        # double e = 0.251 + 0.935 * esq + 0.709 * x + 0.030 * esq * x - 0.472 * esqsq - 2.427 * x * x
        # double a2 = -1.265 + 0.220 * esq + 2.651 * x + 1.010 * esq * x - 1.815 * esqsq - 7.657 * x * x
        # double a4 = 0.556 - 1.465 * esq - 4.260 * x - 2.327 * esq * x + 4.921 * esqsq + 12.98 * x * x

        double g = 1 + a2 * mu**2 + a4 * pow(mu, 4) - (1 + a2 + a4) * pow(mu, 6)

    #     double c_e = -0.791 + 0.776 * x
    #     double c_p = 1.138 - 1.431 * x
    #     double d_e = (-1.315 + 2.431 * x) * esq * x
    #     double d_p = (0.653 - 2.864 * x) * esq * x
    #     double d_60 = (13.47 - 27.13 * x) * esq * x
    #     double f_e = -1.172 * x * esqsq
    #     double f_p = 0.975 * x * esqsq
    #     double g = 1.0

    # g += (c_e + d_e + f_e) * esq * (1.0 - mu * mu)
    # g += (c_p + d_p + f_p - d_60) * esq * mu * mu
    # g += d_60 * esq * fabs(mu)

    return log10(g * g_0) + 2.0
