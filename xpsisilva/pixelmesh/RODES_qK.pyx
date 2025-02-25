#cython: cdivision=True
#cython: boundscheck=False
#cython: nonecheck=False
#cython: wraparound=False

from libc.math cimport M_PI, sin, cos
from libc.stdio cimport printf
from GSL cimport GSL_SUCCESS
from xpsisilva.pixelmesh.METRIC_qK cimport *

cdef double _pi = M_PI

# Partial derivatives of the metric tensor with respect to coordinates r and \theta
cdef double dg_11dr(double r,
                    double theta,
                    double r_s,
                    double a,
                    double Sigma,
                    double Delta,
                    double kappa,
                    double func_theta,
                    double F1_r,
                    double dF1dr_r,
                    double func_r,
                    double sin_theta) nogil:

    cdef double temp = (2.0*r - r_s)*(r*r + a*a + r_s*r*a*a*sin_theta*sin_theta / Sigma) / (Delta*Delta)

    temp -= (2.0*r + r_s*a*a*sin_theta*sin_theta / Sigma - 2.0*r_s*r*r*a*a*sin_theta*sin_theta / (Sigma*Sigma)) / Delta # Reduce #flops

    return temp + kappa*func_theta*(dF1dr_r - F1_r*r_s / (r*r*func_r)) / func_r

cdef double dg_11dtheta(double r,
                        double theta,
                        double r_s,
                        double a,
                        double Sigma,
                        double Delta,
                        double kappa,
                        double F1_r,
                        double func_r,
                        double sin_theta,
                        double cos_theta) nogil:

    return 2.0*cos_theta*sin_theta*(3.0*kappa*F1_r / func_r - r_s*a*a*r*(1.0 + a*a*sin_theta*sin_theta / Sigma) / (Delta*Sigma))

cdef double dg_22dr(double r,
                    double theta,
                    double r_s,
                    double Sigma,
                    double Delta,
                    double kappa,
                    double func_theta,
                    double F1_r,
                    double dF1dr_r,
                    double func_r) nogil:

    return (2.0*r*(1.0 - Delta / Sigma) - r_s) / Sigma + kappa*func_theta*(func_r*dF1dr_r + r_s*F1_r / (r*r))

cdef double dg_22dtheta(double r,
                        double theta,
                        double a,
                        double Sigma,
                        double Delta,
                        double kappa,
                        double F1_r,
                        double func_r,
                        double sin_theta,
                        double cos_theta) nogil:

    return 2.0*cos_theta*sin_theta*(Delta*a*a / (Sigma*Sigma) + 3.0*kappa*func_r*F1_r)

cdef double dg_33dr(double r,
                    double theta,
                    double Sigma,
                    double kappa,
                    double func_theta,
                    double F2_r,
                    double dF2dr_r) nogil:

    return -2.0*r / (Sigma*Sigma) - kappa*func_theta*(dF2dr_r - 2.0*F2_r / r) / (r*r)

cdef double dg_33dtheta(double r,
                        double theta,
                        double a,
                        double Sigma,
                        double kappa,
                        double F2_r,
                        double sin_theta,
                        double cos_theta) nogil:

    return 2.0*sin_theta*cos_theta*(a*a / (Sigma*Sigma) - 3.0*kappa*F2_r / (r*r))

cdef double dg_44dr(double r,
                    double theta,
                    double r_s,
                    double a,
                    double Sigma,
                    double Delta,
                    double kappa,
                    double sin_theta,
                    double func_theta,
                    double F2_r,
                    double dF2dr_r) nogil:

    cdef double temp = 2.0*r - r_s - (Delta - a*a*sin_theta*sin_theta)*((2.0*r - r_s)*Sigma + 2.0*r*Delta) / (Delta*Sigma)

    temp /= (Delta*Sigma*sin_theta*sin_theta)

    return temp - kappa*func_theta*(dF2dr_r - 2.0*F2_r / r) / (r*r*sin_theta*sin_theta)

cdef double dg_44dtheta(double r,
                        double theta,
                        double a,
                        double Sigma,
                        double Delta,
                        double kappa,
                        double F2_r,
                        double sin_theta,
                        double cos_theta,
                        double func_theta) nogil:

    cdef double temp = (-a*a - (Delta - a*a*sin_theta*sin_theta)*(Sigma - a*a*sin_theta*sin_theta) / (Sigma*sin_theta*sin_theta)) / (Delta*Sigma*sin_theta*sin_theta)

    temp += kappa*F2_r*(func_theta / (sin_theta*sin_theta) - 3.0) / (r*r*sin_theta*sin_theta)

    return temp * 2.0*sin_theta*cos_theta

cdef double dg_14dr(double r,
                    double theta,
                    double r_s,
                    double a,
                    double Sigma,
                    double Delta) nogil:

    return r_s*a*(r*((2.0*r - r_s) / Delta + 2.0*r / Sigma) - 1.0) / (Delta*Sigma)

cdef double dg_14dtheta(double r,
                        double theta,
                        double r_s,
                        double a,
                        double Sigma,
                        double Delta,
                        double sin_theta,
                        double cos_theta) nogil:

    return -2.0*r_s*r*a*a*a*cos_theta*sin_theta / (Delta*Sigma*Sigma)

# quasi-Kerr connection coefficients (Christoffel symbols of the second kind)
cdef double Gamma_211(double *m,
                      double *dr,
                      double *dt,
                      double det) nogil:

    return dr[4] / det - m[4]*(m[0]*dr[4] + m[4]*dr[0] - 2.0*m[1]*dr[1]) / (det*det)

cdef double Gamma_222(double *m,
                      double *dr,
                      double *dt,
                      double det) nogil:

    return -dr[2] / (m[2] * m[2])

cdef double Gamma_233(double * m, double * dr, double * dt, double det) nogil:

    return -dr[3] / (m[3]*m[3])

cdef double Gamma_244(double *m,
                      double *dr,
                      double *dt,
                      double det) nogil:

    return dr[0] / det - m[0]*(m[0]*dr[4] + m[4]*dr[0] - 2.0*m[1]*dr[1]) / (det*det)

cdef double Gamma_241(double *m,
                      double *dr,
                      double *dt,
                      double det) nogil:

    return -dr[1] / det + m[1]*(m[0]*dr[4] + m[4]*dr[0] - 2.0*m[1]*dr[1]) / (det*det)

cdef double Gamma_232(double *m,
                      double *dr,
                      double *dt,
                      double det) nogil:

    return -dt[2] / (m[2]*m[2])

cdef double Gamma_311(double *m,
                      double *dr,
                      double *dt,
                      double det) nogil:

    return dt[4] / det - m[4]*(m[0]*dt[4] + m[4]*dt[0] - 2.0*m[1]*dt[1]) / (det*det)

cdef double Gamma_322(double *m,
                      double *dr,
                      double *dt,
                      double det) nogil:

    return -dt[2] / (m[2] * m[2])

cdef double Gamma_333(double *m,
                      double *dr,
                      double *dt,double det) nogil:

    return -dt[3] / (m[3] * m[3])

cdef double Gamma_344(double *m,
                      double *dr,
                      double *dt,
                      double det) nogil:

    return dt[0] / det - m[0]*(m[0]*dt[4] + m[4]*dt[0] - 2.0*m[1]*dt[1]) / (det*det)

cdef double Gamma_341(double *m,
                      double *dr,
                      double *dt,
                      double det) nogil:

    return -dt[1] / det + m[1]*(m[0]*dt[4] + m[4]*dt[0] - 2.0*m[1]*dt[1]) / (det*det)

cdef double Gamma_332(double *m,
                      double *dr,
                      double *dt,
                      double det) nogil:

    return -dr[3] / (m[3] * m[3])

# Six-dimensional coupled second-order non-linear system of differential equations
# Struct for spacetime parameters useful in this module?
cdef int RODES(double t,
                  const double *const y,
                  double *const dydl,
                  void *const params) nogil:

    # y[] := y[t, phi, r, r', theta, theta']

    cdef:
        double r_s = (<double*> params)[0]
        double a = (<double*> params)[1]
        double b = (<double*> params)[2]
        double kappa = (<double*> params)[3]
        double sin_theta = sin(y[4])
        double cos_theta = cos(y[4])
        double Sigma = y[2]*y[2] + a*a*cos_theta*cos_theta
        double Delta  = y[2]*y[2] - r_s*y[2] + a*a
        double func_theta = 1.0 - 3.0*cos_theta*cos_theta
        double func_r = 1.0 - r_s/y[2]
        double F1_r = F1(y[2], r_s)
        double F2_r = F2(y[2], r_s)
        double dF1dr_r = dF1dr(y[2], r_s)
        double dF2dr_r = dF2dr(y[2], r_s)
        double X_IP = (<double*> params)[4]
        double Y_IP = (<double*> params)[5]
        double inclination = (<double*> params)[6]

    cdef double metric[5]
    metric[0] = g_11(y[2], y[4], r_s, a, Sigma, Delta, kappa, sin_theta, func_theta, F1_r, func_r)
    metric[1] = g_14(y[2], y[4], r_s, a, Sigma, Delta)
    metric[2] = g_22(y[2], y[4], Sigma, Delta, kappa, func_theta, F1_r, func_r)
    metric[3] = g_33(y[2], y[4], Sigma, kappa, func_theta, F2_r)
    metric[4] = g_44(y[2], y[4], a, Sigma, Delta, kappa, sin_theta, func_theta, F2_r)

    cdef double deriv_r[5]
    deriv_r[0] = dg_11dr(y[2], y[4], r_s, a, Sigma, Delta, kappa, func_theta, F1_r, dF1dr_r, func_r, sin_theta)
    deriv_r[1] = dg_14dr(y[2], y[4], r_s, a, Sigma, Delta)
    deriv_r[2] = dg_22dr(y[2], y[4], r_s, Sigma, Delta, kappa, func_theta, F1_r, dF1dr_r, func_r) 
    deriv_r[3] = dg_33dr(y[2], y[4], Sigma, kappa, func_theta, F2_r, dF2dr_r) 
    deriv_r[4] = dg_44dr(y[2], y[4], r_s, a, Sigma, Delta, kappa, sin_theta, func_theta, F2_r, dF2dr_r)

    cdef double deriv_theta[5]
    deriv_theta[0] = dg_11dtheta(y[2], y[4], r_s, a, Sigma, Delta, kappa, F1_r, func_r, sin_theta, cos_theta)
    deriv_theta[1] = dg_14dtheta(y[2], y[4], r_s, a, Sigma, Delta, sin_theta, cos_theta) 
    deriv_theta[2] = dg_22dtheta(y[2], y[4], a, Sigma, Delta, kappa, F1_r, func_r, sin_theta, cos_theta) 
    deriv_theta[3] = dg_33dtheta(y[2], y[4], a, Sigma, kappa, F2_r, sin_theta, cos_theta)
    deriv_theta[4] = dg_44dtheta(y[2], y[4], a, Sigma, Delta, kappa, F2_r, sin_theta, cos_theta, func_theta)

    cdef double det_metric = det_g(y[2], y[4], r_s, a, Sigma, Delta, kappa, func_theta, F1_r, F2_r, func_r, sin_theta) ##!

    #printf("\n%.8f", det_metric)

    cdef double Gamma_r[6]
    Gamma_r[0] = -0.5 * metric[2] * Gamma_211(metric, deriv_r, deriv_theta, det_metric)
    Gamma_r[1] = 0.5 * metric[2] * Gamma_222(metric, deriv_r, deriv_theta, det_metric)
    Gamma_r[2] = -0.5 * metric[2] * Gamma_233(metric, deriv_r, deriv_theta, det_metric)
    Gamma_r[3] = -0.5 * metric[2] * Gamma_244(metric, deriv_r, deriv_theta, det_metric)
    Gamma_r[4] = -0.5 * metric[2] * Gamma_241(metric, deriv_r, deriv_theta, det_metric)
    Gamma_r[5] = 0.5 * metric[2] * Gamma_232(metric, deriv_r, deriv_theta, det_metric)

    #printf("\nGamma_r: %.8f, %.8f, %.8f, %.8f, %.8f, %.8f",
           #Gamma_r[0],
           #Gamma_r[1],
           #Gamma_r[2],
           #Gamma_r[3],
           #Gamma_r[4],
           #Gamma_r[5])

    cdef double Gamma_theta[6]
    Gamma_theta[0] = -0.5 * metric[3] * Gamma_311(metric, deriv_r, deriv_theta, det_metric)
    Gamma_theta[1] = -0.5 * metric[3] * Gamma_322(metric, deriv_r, deriv_theta, det_metric)
    Gamma_theta[2] = 0.5 * metric[3] * Gamma_333(metric, deriv_r, deriv_theta, det_metric)
    Gamma_theta[3] = -0.5 * metric[3] * Gamma_344(metric, deriv_r, deriv_theta, det_metric)
    Gamma_theta[4] = -0.5 * metric[3] * Gamma_341(metric, deriv_r, deriv_theta, det_metric)
    Gamma_theta[5] = 0.5 * metric[3] * Gamma_332(metric, deriv_r, deriv_theta, det_metric)

    #printf("\nGamma_theta: %.8f, %.8f, %.8f, %.8f, %.8f, %.8f",
           #Gamma_theta[0],
           #Gamma_theta[1],
           #Gamma_theta[2],
           #Gamma_theta[3],
           #Gamma_theta[4],
           #Gamma_r[5])

    # Covariant metric components
    cdef:
        double G_11 = metric[4] / det_metric
        double G_14 = -metric[1] / det_metric
        double G_44 = metric[0] / det_metric

    dydl[0] = (G_44 + b*G_14) / (G_11*G_44 - G_14*G_14) # Reduce #flops here

    # Had only the second condition, thus phase-shifting for finite omega parameter...
    if a == 0.0 and X_IP == 0.0:
        dydl[1] = 0.0
    else:
        dydl[1] = -1.0 * (b*G_11 + G_14) / (G_11*G_44 - G_14*G_14)

    dydl[2] = y[3]
    dydl[3] = (-Gamma_r[0]*dydl[0]*dydl[0]
               - Gamma_r[1]*y[3]*y[3]
               - Gamma_r[2]*y[5]*y[5]
               - Gamma_r[3]*dydl[1]*dydl[1]
               - 2.0*Gamma_r[4]*dydl[0]*dydl[1]
               - 2.0*Gamma_r[5]*y[3]*y[5])
    dydl[4] = y[5]

    # For the Schwarzschild metric, if rays are equatorial, or radial, one or
    # both of the angles are zero, leading to slow evolution because the
    # relative error is small but the derivative(s) are not set to precisely
    # zero.
    if a == 0.0 and X_IP == 0.0 and Y_IP == 0.0:
        dydl[5] = 0.0
    elif a == 0.0 and Y_IP == 0.0 and inclination == _pi/2.0:
        dydl[5] = 0.0
    else:
        dydl[5] = (-Gamma_theta[0]*dydl[0]*dydl[0]
                   - Gamma_theta[1]*y[3]*y[3]
                   - Gamma_theta[2]*y[5]*y[5]
                   - Gamma_theta[3]*dydl[1]*dydl[1]
                   - 2.0*Gamma_theta[4]*dydl[0]*dydl[1]
                   - 2.0*Gamma_theta[5]*y[3]*y[5])

    #printf("\ndydl: %.8f, %.8f, %.8f, %.8f, %.8f, %.8f",
           #dydl[0], dydl[1], dydl[2], dydl[3], dydl[4], dydl[5])

    return GSL_SUCCESS
