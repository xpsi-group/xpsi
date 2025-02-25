#cython: cdivision=True
#cython: boundscheck=False
#cython: nonecheck=False
#cython: wraparound=False

from libc.math cimport *
from libc.stdio cimport printf
from libc.math cimport sqrt, sin, cos, acos, pow
from xpsisilva.pixelmesh.METRIC_qK cimport *
from GSL cimport gsl_isnan

cdef double _pi = M_PI
cdef double c = 2.99792458e8
cdef int SUCCESS = 1

cdef double RADIUS(double x, double epsilon, double zeta) nogil:

    return 1.0 + epsilon * (-0.788 + 1.030 * zeta) * x * x

cdef double LINTERP(double t,
                    double t1,
                    double t2,
                    double x1,
                    double x2) nogil:

    return x1*(t2 - t) / (t2 - t1) + x2*(t - t1) / (t2 - t1)

cdef int BISECT(const double *const y_p,
                const double *const y,
                double *const y_new,
                double affine_p,
                double affine,
                double a,
                double R_eq,
                double epsilon,
                double zeta) nogil:

    cdef:
        int bisect = 1
        size_t i
        double sin_theta, r_sph, cos_theta_sph, frac_diff
        double affine_bisect, affine_1 = affine_p, affine_2 = affine

    while bisect == 1:
        affine_bisect = 0.5 * (affine_1 + affine_2)

        for i in range(6):
            y_new[i] = LINTERP(affine_bisect, affine_p, affine, y_p[i], y[i])

        sin_theta = sin(y_new[4])
        r_sph = sqrt(y_new[2]*y_new[2] + a*a*sin_theta*sin_theta)
        cos_theta_sph = cos(y_new[4]) / sqrt(1.0 + a*a*sin_theta*sin_theta / (y_new[2]*y_new[2]))
        frac_diff = 1.0 - RADIUS(cos_theta_sph, epsilon, zeta) * R_eq / r_sph

        if fabs(frac_diff) < 1.0e-8:
            bisect = 0
        else:
            if frac_diff > 0.0:
                affine_1 = affine_bisect
            else:
                affine_2 = affine_bisect

    return SUCCESS

cdef double ZABB(const _GEOM *const GEOM,
                 const double *const y,
                 double b,
                 double *const Z,
                 double *const ABB) nogil:

    cdef:
        double R = y[2] # BL
        double sin_theta = sin(y[4]) # BL
        double cos_theta = cos(y[4]) # BL
        double u = sqrt(y[2]*y[2] + GEOM.asq * sin_theta * sin_theta) # SPH
        double v = acos(cos_theta / sqrt(1.0 + GEOM.asq * sin_theta * sin_theta / (y[2]*y[2]))) # SPH
        double sin_v = sin(v)
        double cos_v = cos(v)

    cdef:
        double func_theta = 1.0 - 3.0 * cos_theta * cos_theta
        double func_r = 1.0 - GEOM.r_s/y[2]
        double F1_r = F1(y[2], GEOM.r_s)
        double F2_r = F2(y[2], GEOM.r_s)
        double dF1dr_r = dF1dr(y[2], GEOM.r_s)
        double dF2dr_r = dF2dr(y[2], GEOM.r_s)
        double Sigma = y[2]*y[2] + GEOM.asq * cos_theta * cos_theta
        double Delta  = y[2]*y[2] - GEOM.r_s * y[2] + GEOM.asq

    cdef double det_metric = det_g(y[2],
                                   y[4],
                                   GEOM.r_s,
                                   GEOM.a,
                                   Sigma,
                                   Delta,
                                   GEOM.kappa,
                                   func_theta,
                                   F1_r,
                                   F2_r,
                                   func_r,
                                   sin_theta)

    cdef:
        double G_11 = g_44(y[2],
                           y[4],
                           GEOM.a,
                           Sigma,
                           Delta,
                           GEOM.kappa,
                           sin_theta,
                           func_theta,
                           F2_r) / det_metric

        double G_14 = -g_14(y[2],
                            y[4],
                            GEOM.r_s,
                            GEOM.a,
                            Sigma,
                            Delta) / det_metric

        double G_22 = 1.0 / g_22(y[2], y[4], Sigma, Delta, GEOM.kappa, func_theta, F1_r, func_r)

        double G_33 = 1.0 / g_33(y[2], y[4], Sigma, GEOM.kappa, func_theta, F2_r)

        double G_44 = g_11(y[2],
                           y[4],
                           GEOM.r_s,
                           GEOM.a,
                           Sigma,
                           Delta,
                           GEOM.kappa,
                           sin_theta,
                           func_theta,
                           F1_r,
                           func_r) / det_metric

    #cdef double R_SPH = GEOM.R_eq * RADIUS(cos_theta_SPH, GEOM.epsilon, GEOM.zeta)
    cdef double R_prime = -2.0 * GEOM.R_eq * GEOM.epsilon * (-0.788 + 1.030 * GEOM.zeta) * cos_v * sin_v

    # Outgoing null geodesic has impact parameter -b, where b is that of the
    # null geodesic of the incoming photon
    cdef double k_t, k_phi

    k_t = (G_44 + b*G_14) / (-G_11*G_44 + G_14*G_14)
    k_phi = (b*G_11 + G_14) / (G_11*G_44 - G_14*G_14)

    #if k_t < 0.0:
    #    printf("\nWarning: -ve k_t = %.8e", k_t)

    Z[0] = G_11 * k_t
    Z[0] += G_44 * GEOM.omega * k_phi / c
    Z[0] += G_14 * GEOM.omega * k_t / c
    Z[0] += G_14 * k_phi

    Z[0] /= -1.0 * sqrt(-G_11 - 2.0*G_14*GEOM.omega / c - G_44*GEOM.omega*GEOM.omega / (c*c))
    #Z[0] = (1.0 - b*GEOM.omega/c) / sqrt(-G_11 - 2.0*G_14*GEOM.omega / c - G_44*GEOM.omega*GEOM.omega / (c*c))

    cdef double drdu, drdv, dtdu, dtdv, g_uu, g_vv, f

    drdu = u / R + GEOM.asq * u * cos_v * cos_v / pow(R, 3.0)
    drdu /= 1.0 + GEOM.asq * u * u * cos_v * cos_v / pow(R, 4.0)

    drdv = -GEOM.asq * u * u * cos_v * sin_v / pow(R, 3.0)
    drdv /= 1.0 + GEOM.asq * u * u * cos_v * cos_v / pow(R, 4.0)

    dtdu = cos_v / R - u * cos_v * drdu / pow(R, 2.0)
    dtdu /= -sin_theta

    dtdv = -u * sin_v / R - u * cos_v * drdv / pow(R, 2.0)
    dtdv /= -sin_theta

    g_uu = drdu * drdu * G_22 + dtdu * dtdu * G_33
    g_vv = dtdv * dtdv * G_33 + drdv * drdv * G_22

    f = sqrt(g_uu / g_vv) * R_prime

    ABB[0] = G_22 * y[3] * (drdu / sqrt(g_uu) - drdv * f / sqrt(g_vv))
    ABB[0] += G_33 * y[5] * (dtdu / sqrt(g_uu) - dtdv * f / sqrt(g_vv))
    ABB[0] /= -1.0 * Z[0] * sqrt(1.0 + f * f)
    # the factor of \pi is required because the traced photon wavevector is
    # ingoing at the NS surface

    #ABB[0] = -1.0 * sqrt(G_22)*(r_dot - theta_dot*R_prime) / (Z[0]*sqrt(1.0 + R_prime*R_prime / (R*R)))

    #cdef double dR_dtheta, dt_dtsph, temp

    #dr_dtheta = 2.0 * R * R * R_SPH * R_prime
    #dr_dtheta += 2.0 * GEOM.asq * R_SPH * R_prime * cos_theta_SPH * cos_theta_SPH
    #dr_dtheta -= 2.0 * GEOM.asq * R_SPH * R_SPH * sin_theta_SPH * cos_theta_SPH

    #dr_dtheta /= 4.0 * pow(R, 3.0) - 2.0 * R * (R_SPH * R_SPH - GEOM.asq)

    #dt_dtsph = 2.0 * R_SPH * R_prime * cos_theta_SPH * cos_theta_SPH
    #dt_dtsph -= 2.0 * R_SPH * R_SPH * sin_theta_SPH * cos_theta_SPH
    #dt_dtsph -= 2.0 * R_SPH * R_prime * cos_theta * cos_theta

    #temp = 2.0 * GEOM.asq * pow(sin_theta, 3.0) * cos_theta
    #temp -= 2.0 * R_SPH * R_SPH * sin_theta * cos_theta
    #temp -= 2.0 * GEOM.asq * sin_theta * pow(cos_theta, 3.0)

    #dt_dtsph /= temp

    #ABB[0] = -1.0 * sqrt(G_22) * (y[3] - y[5] * dr_dtheta)
    #ABB[0] /= Z[0] * sqrt(1.0 + G_22 * dr_dtheta * dr_dtheta / G_33)

    #if gsl_isnan(acos(sqrt(G_22)*(r_dot - theta_dot*R_prime) / (Z[0]*sqrt(1.0 + R_prime*R_prime / (R*R))))) == 1:
    #    printf("\narccos argument: %.8e", sqrt(G_22)*(r_dot - theta_dot*R_prime) / (Z[0]*sqrt(1.0 + R_prime*R_prime / (R*R))))

    if ABB[0] < 0.0 or ABB[0] > 1.0:
        printf("\nWarning: cosine aberrated surface zenith angle = %.8e...", ABB[0])
        #printf("\nderiv: %.8e", dr_dtheta)
        #printf("\nratio: %.8e", G_22 * R * R/G_33)

        if ABB[0] > 1.0:
            printf("\nWarning: forcing the aberrated ray to be normal.")
            ABB[0] = 1.0

        if ABB[0] < 0.0:
            printf("\nWarning: forcing the aberrated ray to be tangential.")
            ABB[0] = 0.0

    #if Z[0] < 0.0:
    #    printf("\nWarning: -ve redshift %.8e", Z[0])
    #    printf("\nWarning: forcing the redshift to zero.")
    #    Z[0] = 0.0

    if gsl_isnan(G_22) == 1:
        printf("\nG_22 is NaN.")

    if gsl_isnan(R_prime) == 1:
        printf("\nR_prime is NaN.")

    if gsl_isnan(R) == 1:
        printf("\nR is NaN.")

    if gsl_isnan(Z[0]) == 1:
        printf("\nZ is NaN.")
