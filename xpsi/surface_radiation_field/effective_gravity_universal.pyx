#cython: cdivision=True
#cython: boundscheck=False
#cython: nonecheck=False
#cython: wraparound=False

from libc.math cimport sqrt, log10, fabs

cdef double c = 2.99792458e8

cdef double effectiveGravity(double mu,
                             double R_eq,
                             double x,
                             double epsilon) nogil:

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
