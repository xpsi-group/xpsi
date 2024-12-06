#cython: cdivision=True
#cython: boundscheck=False
#cython: nonecheck=False
#cython: wraparound=False

from libc.math cimport log

cdef extern from "math.h":
    long double logl(long double x) noexcept nogil

# Quadrupole moment correction F(r) functions
cdef double F1(double r, double r_s) noexcept nogil:

    return -5.0*(2.0*r - r_s) * (r_s*r_s + 6.0*r_s*r - 6.0*r*r) / (16.0*r_s*r*(r - r_s)) - 15.0*r*(r - r_s)*log(r/(r - r_s)) / (4.0*r_s*r_s)

cdef double F2(double r, double r_s) noexcept nogil:

    return 5.0*(r_s*r_s - 3.0*r_s*r - 6.0*r*r) / (8.0*r_s*r) + 15.0*(2.0*r*r - r_s*r_s)*log(r/(r - r_s)) / (8.0*r_s*r_s)

cdef double dF1dr(double r, double r_s) noexcept nogil:

    cdef double temp = (-5.0/(16.0*r_s*r*(r - r_s)))*(4.0*(-9.0*r*r + 9.0*r_s*r - r_s*r_s) - (2.0*r - r_s)*(2.0*r - r_s)*(r_s*r_s + 6.0*r_s*r - 6.0*r*r) / (r*(r - r_s)))

    return temp - (15.0/(4.0*r_s*r_s))*(log(r/(r - r_s))*(2.0*r - r_s) - r_s)

cdef double dF2dr(double r, double r_s) noexcept nogil:

    return (-5.0/(8.0*r_s))*(6.0 + r_s*r_s / (r*r)) + (15.0/(8.0*r_s*r_s))*(4.0*r*log(r/(r - r_s)) - r_s*(2.0*r*r - r_s*r_s) / (r*(r - r_s)))

# Quadrupole moment correction F(r) functions for more quadrupole precision
cdef long double F1_l(long double r, long double r_s) noexcept nogil:

    return -5.0*(2.0*r - r_s) * (r_s*r_s + 6.0*r_s*r - 6.0*r*r) / (16.0*r_s*r*(r - r_s)) - 15.0*r*(r - r_s)*logl(r/(r - r_s)) / (4.0*r_s*r_s)

cdef long double F2_l(long double r, long double r_s) noexcept nogil:

    return 5.0*(r_s*r_s - 3.0*r_s*r - 6.0*r*r) / (8.0*r_s*r) + 15.0*(2.0*r*r - r_s*r_s)*logl(r/(r - r_s)) / (8.0*r_s*r_s)

cdef long double dF1dr_l(long double r, long double r_s) noexcept nogil:

    cdef long double temp = (-5.0/(16.0*r_s*r*(r - r_s)))*(4.0*(-9.0*r*r + 9*r_s*r - r_s*r_s) - (2.0*r - r_s)*(2.0*r - r_s)*(r_s*r_s + 6.0*r_s*r - 6.0*r*r) / (r*(r - r_s)))

    return temp - (15.0/(4.0*r_s*r_s))*(logl(r/(r - r_s))*(2.0*r - r_s) - r_s)

cdef long double dF2dr_l(long double r, long double r_s) noexcept nogil:

    return (-5.0/(8.0*r_s))*(6.0 + r_s*r_s / (r*r)) + (15.0/(8.0*r_s*r_s))*(4.0*r*logl(r/(r - r_s)) - r_s*(2.0*r*r - r_s*r_s) / (r*(r - r_s)))

# quasi-Kerr contravariant metric tensor components g_uv = g_uv^K + \kappa h_uv
cdef double g_11(double r,
                 double theta,
                 double r_s,
                 double a,
                 double Sigma,
                 double Delta,
                 double kappa,
                 double sin_theta,
                 double func_theta,
                 double F1_r,
                 double func_r) noexcept nogil:

    cdef double temp = r*r + a*a + r_s * r * a * a * sin_theta * sin_theta / Sigma

    temp *= -1.0 / Delta

    return temp + kappa * func_theta * F1_r / func_r

cdef double g_14(double r,
                 double theta,
                 double r_s,
                 double a,
                 double Sigma,
                 double Delta) noexcept nogil:

    return -r_s*a*r / (Sigma*Delta)

cdef double g_22(double r,
                 double theta,
                 double Sigma,
                 double Delta,
                 double kappa,
                 double func_theta,
                 double F1_r,
                 double func_r) noexcept nogil:

    return Delta / Sigma + kappa*F1_r*func_theta*func_r

cdef double g_33(double r,
                 double theta,
                 double Sigma,
                 double kappa,
                 double func_theta,
                 double F2_r) noexcept nogil:

    return 1.0 / Sigma - kappa * F2_r * func_theta / (r * r)

cdef double g_44(double r,
                 double theta,
                 double a,
                 double Sigma,
                 double Delta,
                 double kappa,
                 double sin_theta,
                 double func_theta,
                 double F2_r) noexcept nogil:

    cdef double temp = (Delta - a * a * sin_theta * sin_theta) / (Sigma * Delta)

    temp -= kappa * F2_r * func_theta / (r * r)

    return temp / (sin_theta * sin_theta)

# quasi-Kerr reduced (t, \phi) contravariant metric tensor determinant
cdef double det_g(double r,
                  double theta,
                  double r_s,
                  double a,
                  double Sigma,
                  double Delta,
                  double kappa,
                  double func_theta,
                  double F1_r,
                  double F2_r,
                  double func_r,
                  double sin_theta) noexcept nogil:

    cdef:
        double G_11 = g_11(r, theta, r_s, a, Sigma, Delta, kappa, sin_theta, func_theta, F1_r, func_r)
        double G_14 = g_14(r, theta, r_s, a, Sigma, Delta)
        double G_44 = g_44(r, theta, a, Sigma, Delta, kappa, sin_theta, func_theta, F2_r)

    return G_11*G_44 - G_14*G_14

# quasi-Kerr contravariant metric tensor components g_uv = g_uv^K + \kappa h_uv
cdef long double g_11_l(long double r,
                        long double theta,
                        long double r_s,
                        long double a,
                        long double Sigma,
                        long double Delta,
                        long double kappa,
                        long double sin_theta,
                        long double func_theta,
                        long double F1_r,
                        long double func_r) noexcept nogil:

    return -1.0*(r*r + a*a + r_s*r*a*a*sin_theta*sin_theta / Sigma) / Delta + kappa*func_theta*F1_r / func_r

cdef long double g_14_l(long double r,
                        long double theta,
                        long double r_s,
                        long double a,
                        long double Sigma,
                        long double Delta) noexcept nogil:

    return -r_s*a*r / (Sigma*Delta)

cdef long double g_22_l(long double r,
                        long double theta,
                        long double Sigma,
                        long double Delta,
                        long double kappa,
                        long double func_theta,
                        long double F1_r,
                        long double func_r) noexcept nogil:

    return Delta / Sigma + kappa*F1_r*func_theta*func_r

cdef long double g_33_l(long double r,
                        long double theta,
                        long double Sigma,
                        long double kappa,
                        long double func_theta,
                        long double F2_r) noexcept nogil:

    return 1.0 / Sigma - kappa*F2_r*func_theta / (r*r)

cdef long double g_44_l(long double r,
                        long double theta,
                        long double a,
                        long double Sigma,
                        long double Delta,
                        long double kappa,
                        long double sin_theta,
                        long double func_theta,
                        long double F2_r) noexcept nogil:

    cdef long double temp = Delta - a * a * sin_theta * sin_theta

    temp /= Sigma * Delta * sin_theta * sin_theta

    return  temp - kappa * F2_r * func_theta / (r * r * sin_theta * sin_theta)

# quasi-Kerr reduced (t, \phi) contravariant metric tensor determinant 
cdef long double det_g_l(long double r,
                         long double theta,
                         long double r_s,
                         long double a,
                         long double Sigma,
                         long double Delta,
                         long double kappa,
                         long double func_theta,
                         long double F1_r,
                         long double F2_r,
                         long double func_r,
                         long double sin_theta) noexcept nogil:

    cdef:
        long double G_11 = g_11_l(r, theta, r_s, a, Sigma, Delta, kappa, sin_theta, func_theta, F1_r, func_r)
        long double G_14 = g_14_l(r, theta, r_s, a, Sigma, Delta)
        long double G_44 = g_44_l(r, theta, a, Sigma, Delta, kappa, sin_theta, func_theta, F2_r)

    return G_11 * G_44 - G_14 * G_14
