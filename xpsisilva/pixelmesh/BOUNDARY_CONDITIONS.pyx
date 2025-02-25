#cython: cdivision=True
#cython: boundscheck=False
#cython: nonecheck=False
#cython: wraparound=False

from libc.math cimport M_PI
from libc.stdio cimport printf
from GSL cimport GSL_SUCCESS

cdef extern from "math.h" nogil:
    long double sinl(long double x)
    long double asinl(long double x)
    long double cosl(long double x)
    long double acosl(long double x)
    long double tanl(long double x)
    long double atanl(long double x)
    long double sqrtl(long double x)

from xpsisilva.pixelmesh.METRIC_qK cimport *

cdef double _pi = M_PI

cdef void COMPUTE_BCs(_RAY *const RAY,
                      const _GEOM *const GEOM) nogil:

    cdef:
        long double d = <long double>GEOM.d
        long double inclination = <long double>GEOM.inclination
        long double r_s = <long double>GEOM.r_s
        long double a = <long double>GEOM.a
        long double kappa = <long double>GEOM.kappa
        long double X_IP = <long double>RAY.X_IP
        long double Y_IP = <long double>RAY.Y_IP

        long double sin_inclination = sinl(inclination)
        long double cos_inclination = cosl(inclination)
        long double xsq = X_IP * X_IP
        long double ysq = Y_IP * Y_IP
        long double dsq = d * d
        long double dsin = d * sin_inclination
        long double dcos = d * cos_inclination
        long double Xsin = X_IP * sin_inclination
        long double Ycos = Y_IP * cos_inclination
        long double Ysin = Y_IP * sin_inclination
        long double asq = a * a
        long double Sigma_l, Delta_l, sin_theta_l, cos_theta_l
        long double func_theta_l, func_r_l, F1_r_l, F2_r_l, dF1dr_r_l, dF2dr_r_l
        long double Y[6]

    # Compute boundary conditions at image plane
    Y[0] = 0.0

    Y[1] = atanl(X_IP / (dsin - Ycos))

    Y[2] = sqrtl(dsq + xsq + ysq)
    Y[3] = -d / sqrtl(dsq + xsq + ysq)
    Y[4] = acosl((dcos + Ysin) / sqrtl(dsq + xsq + ysq))

    Y[5] = (cos_inclination - (dsq*cos_inclination + Y_IP*dsin) / (dsq + xsq + ysq)) / sqrtl(xsq + (dsin - Ycos)*(dsin - Ycos))

    sin_theta_l = sinl(Y[4])
    cos_theta_l = cosl(Y[4])
    Sigma_l = Y[2] * Y[2] + asq * cos_theta_l * cos_theta_l
    Delta_l  = Y[2] * Y[2] - r_s * Y[2] + asq
    func_theta_l = 1.0 - 3.0 * cos_theta_l * cos_theta_l
    func_r_l = 1.0 - r_s / Y[2]
    F1_r_l = F1_l(Y[2], r_s)
    F2_r_l = F2_l(Y[2], r_s)
    dF1dr_r_l = dF1dr_l(Y[2], r_s)
    dF2dr_r_l = dF2dr_l(Y[2], r_s)

    cdef long double det_metric_l = det_g_l(Y[2],
                                            Y[4],
                                            r_s,
                                            a,
                                            Sigma_l,
                                            Delta_l,
                                            kappa,
                                            func_theta_l,
                                            F1_r_l,
                                            F2_r_l,
                                            func_r_l,
                                            sin_theta_l)

    # Covariant metric components
    cdef:
        long double G_11_l = g_44_l(Y[2],
                                    Y[4],
                                    a,
                                    Sigma_l,
                                    Delta_l,
                                    kappa,
                                    sin_theta_l,
                                    func_theta_l,
                                    F2_r_l)

        long double G_14_l = -g_14_l(Y[2],
                                     Y[4],
                                     r_s,
                                     a,
                                     Sigma_l,
                                     Delta_l)
        long double G_22_l = 1.0
        long double G_33_l = 1.0

        long double G_44_l = g_11_l(Y[2],
                                    Y[4],
                                    r_s,
                                    a,
                                    Sigma_l,
                                    Delta_l,
                                    kappa,
                                    sin_theta_l,
                                    func_theta_l,
                                    F1_r_l,
                                    func_r_l)

    G_11_l /= det_metric_l
    G_14_l /= det_metric_l
    G_44_l /= det_metric_l

    G_22_l /= g_22_l(Y[2],
                     Y[4],
                     Sigma_l,
                     Delta_l,
                     kappa,
                     func_theta_l,
                     F1_r_l,
                     func_r_l)

    G_33_l /= g_33_l(Y[2],
                     Y[4],
                     Sigma_l,
                     kappa,
                     func_theta_l,
                     F2_r_l)

    # Require that the photon wavevector is null
    cdef long double k_phi_l = -1.0 * Xsin / ((dsin - Ycos)*(dsin - Ycos) + xsq)
    cdef long double temp_l = G_14_l * k_phi_l / G_11_l
    cdef long double k_t_l = -temp_l + sqrtl(temp_l*temp_l - (G_22_l*Y[3]*Y[3] + G_33_l*Y[5]*Y[5] + G_44_l*k_phi_l*k_phi_l) / G_11_l)

    #printf("+ve root: %.8e", <double>k_t_l)
    #printf("-ve root: %.8e", <double>k_t_l)

    #printf("constant: %.8e", <double>(G_11_l * k_t_l + G_14_l * k_phi_l))

    RAY.IMPACT = <double>(G_44_l * k_phi_l + G_14_l * k_t_l)
    RAY.XI = <double>((G_22_l*Y[3]*Y[3] + G_33_l*Y[5]*Y[5] + G_44_l*k_phi_l*k_phi_l + 2.0*G_14_l*k_t_l*k_phi_l) / (G_11_l*k_t_l*k_t_l))
    RAY.XI += 1.0
    RAY.STATE[0] = <double>Y[0]
    RAY.STATE[1] = <double>Y[1]
    RAY.STATE[2] = <double>Y[2]
    RAY.STATE[3] = <double>Y[3]
    RAY.STATE[4] = <double>Y[4]
    RAY.STATE[5] = <double>Y[5]
