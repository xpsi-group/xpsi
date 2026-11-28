#cython: cdivision=True
#cython: boundscheck=False
#cython: nonecheck=False
#cython: wraparound=False

## Right-hand side of the photon geodesic equations in the quasi-Kerr metric,
## for the GSL ODE integrator (called from RK_IP2S_tracer.pyx).
##
## Expressions are linked to equations in:
##   Psaltis & Johannsen 2012, ApJ 745, 1 (PJ12) - the ray-tracing scheme
##   Bauböck et al. 2012, ApJ 753, 175 (B12) - the neutron-star application
##   Glampedakis & Babak 2006, CQG 23, 4167 (GB06) - the quasi-Kerr metric
## Conventions (units, index labels 1-4 = t, r, theta, phi, the helper
## quantities Sigma, Delta, func_theta, func_r) are as in METRIC_qK.pyx.
##
## Structure of this module:
##   1. dg_XXdr, dg_XXdtheta: partial derivatives of the CONTRAVARIANT metric
##      components g^uv (METRIC_qK.pyx) with respect to r and theta.
##      The metric is stationary and axisymmetric, so d/dt = d/dphi = 0.
##   2. Gamma_XYZ: derivatives of the COVARIANT metric components g_uv,
##      built from item 1 (see the note above Gamma_211).
##   3. RODES: assembles the Christoffel symbols and returns dy/dlambda
##      for PJ12 Eqs. (9)-(12).
##
## All derivatives and Christoffel symbols in this module have been checked
## against independent symbolic/numerical differentiation of the metric.
##
from libc.math cimport M_PI, sin, cos
from libc.stdio cimport printf
from GSL cimport GSL_SUCCESS
from xpsi.pixelmesh.METRIC_qK cimport *

from ..tools.core cimport are_equal

cdef double _pi = M_PI

## Partial derivatives of the CONTRAVARIANT metric components g^uv
## with respect to r and theta (labels 1-4 = t, r, th, ph).
## Kerr part: d/dr of -(r^2 + a^2 + r_s r a^2 sin^2/Sigma) / Delta.
## Quadrupole part: d/dr of kappa (1-3cos^2) F1 / (1 - r_s/r), using
## d(1 - r_s/r)/dr = r_s/r^2 and dF1/dr from METRIC_qK.dF1dr.
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
                    double sin_theta) noexcept nogil:

    cdef double temp = (2.0*r - r_s)*(r*r + a*a + r_s*r*a*a*sin_theta*sin_theta / Sigma) / (Delta*Delta)

    temp -= (2.0*r + r_s*a*a*sin_theta*sin_theta / Sigma - 2.0*r_s*r*r*a*a*sin_theta*sin_theta / (Sigma*Sigma)) / Delta # Reduce #flops

    return temp + kappa*func_theta*(dF1dr_r - F1_r*r_s / (r*r*func_r)) / func_r

## d g^tt / d theta. Uses d(1-3cos^2)/dtheta = 6 cos sin and
## dSigma/dtheta = -2 a^2 cos sin; overall factor 2 cos sin pulled out.
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
                        double cos_theta) noexcept nogil:

    return 2.0*cos_theta*sin_theta*(3.0*kappa*F1_r / func_r - r_s*a*a*r*(1.0 + a*a*sin_theta*sin_theta / Sigma) / (Delta*Sigma))

## d g^rr / dr: Kerr part d(Delta/Sigma)/dr; quadrupole part
## d/dr of kappa (1-3cos^2) F1 (1 - r_s/r).
cdef double dg_22dr(double r,
                    double theta,
                    double r_s,
                    double Sigma,
                    double Delta,
                    double kappa,
                    double func_theta,
                    double F1_r,
                    double dF1dr_r,
                    double func_r) noexcept nogil:

    return (2.0*r*(1.0 - Delta / Sigma) - r_s) / Sigma + kappa*func_theta*(func_r*dF1dr_r + r_s*F1_r / (r*r))

## d g^rr / d theta (Kerr part via Sigma; quadrupole part via 1-3cos^2).
cdef double dg_22dtheta(double r,
                        double theta,
                        double a,
                        double Sigma,
                        double Delta,
                        double kappa,
                        double F1_r,
                        double func_r,
                        double sin_theta,
                        double cos_theta) noexcept nogil:

    return 2.0*cos_theta*sin_theta*(Delta*a*a / (Sigma*Sigma) + 3.0*kappa*func_r*F1_r)

## d g^thth / dr: Kerr part d(1/Sigma)/dr = -2r/Sigma^2; quadrupole part
## d/dr of -kappa (1-3cos^2) F2 / r^2.
cdef double dg_33dr(double r,
                    double theta,
                    double Sigma,
                    double kappa,
                    double func_theta,
                    double F2_r,
                    double dF2dr_r) noexcept nogil:

    return -2.0*r / (Sigma*Sigma) - kappa*func_theta*(dF2dr_r - 2.0*F2_r / r) / (r*r)

## d g^thth / d theta.
cdef double dg_33dtheta(double r,
                        double theta,
                        double a,
                        double Sigma,
                        double kappa,
                        double F2_r,
                        double sin_theta,
                        double cos_theta) noexcept nogil:

    return 2.0*sin_theta*cos_theta*(a*a / (Sigma*Sigma) - 3.0*kappa*F2_r / (r*r))

## d g^phph / dr: Kerr part d/dr of (Delta - a^2 sin^2)/(Sigma Delta sin^2);
## quadrupole part d/dr of -kappa (1-3cos^2) F2 / (r^2 sin^2).
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
                    double dF2dr_r) noexcept nogil:

    cdef double temp = 2.0*r - r_s - (Delta - a*a*sin_theta*sin_theta)*((2.0*r - r_s)*Sigma + 2.0*r*Delta) / (Delta*Sigma)

    temp /= (Delta*Sigma*sin_theta*sin_theta)

    return temp - kappa*func_theta*(dF2dr_r - 2.0*F2_r / r) / (r*r*sin_theta*sin_theta)

## d g^phph / d theta, with the overall factor 2 sin cos pulled out.
## The quadrupole term combines d(1-3cos^2)/dtheta and d(1/sin^2)/dtheta.
cdef double dg_44dtheta(double r,
                        double theta,
                        double a,
                        double Sigma,
                        double Delta,
                        double kappa,
                        double F2_r,
                        double sin_theta,
                        double cos_theta,
                        double func_theta) noexcept nogil:

    cdef double temp = (-a*a - (Delta - a*a*sin_theta*sin_theta)*(Sigma - a*a*sin_theta*sin_theta) / (Sigma*sin_theta*sin_theta)) / (Delta*Sigma*sin_theta*sin_theta)

    temp += kappa*F2_r*(func_theta / (sin_theta*sin_theta) - 3.0) / (r*r*sin_theta*sin_theta)

    return temp * 2.0*sin_theta*cos_theta

## d g^tph / dr. Kerr only: h^tph = 0 (GB06 Eq. 2.16, PJ12 Eq. 5), so
## the quadrupole correction does not enter the frame-dragging term.
cdef double dg_14dr(double r,
                    double theta,
                    double r_s,
                    double a,
                    double Sigma,
                    double Delta) noexcept nogil:

    return r_s*a*(r*((2.0*r - r_s) / Delta + 2.0*r / Sigma) - 1.0) / (Delta*Sigma)

## d g^tph / d theta (Kerr only, as above).
cdef double dg_14dtheta(double r,
                        double theta,
                        double r_s,
                        double a,
                        double Sigma,
                        double Delta,
                        double sin_theta,
                        double cos_theta) noexcept nogil:

    return -2.0*r_s*r*a*a*a*cos_theta*sin_theta / (Delta*Sigma*Sigma)

## Derivatives of the covariant quasi-Kerr metric components, used to build
## the Christoffel symbols in the geodesic equations.
##
## Naming: the Gamma_XYZ functions do not return Christoffel symbols.
## Gamma_XYZ returns the one covariant-metric derivative
## that the Christoffel symbol Gamma^X_YZ depends on (labels 1-4 = t, r, th, ph).
## For most this is d_X g_YZ, e.g. Gamma_211 = d g_tt / dr. The two mixed
## r-theta ones differ: Gamma_232 = d g_rr / d theta (for Gamma^r_rth) and
## Gamma_332 = d g_thth / dr (for Gamma^th_rth).
## RODES (below) multiplies each by +/- (1/2) g^rr or g^thth.
##
## Arguments:
##   m[0..4]  = contravariant g^tt, g^tph, g^rr, g^thth, g^phph
##   dr[0..4] = their r-derivatives (dg_XXdr), same ordering
##   dt[0..4] = their theta-derivatives (dg_XXdtheta), same ordering
##   det      = m[0] m[4] - m[1]^2 (METRIC_qK.det_g)
##
## The (t, phi) block is inverted as in METRIC_qK.det_g:
##   g_tt = g^phph/det,  g_tph = -g^tph/det,  g_phph = g^tt/det,
## with d(det) = m[0] d m[4] + m[4] d m[0] - 2 m[1] d m[1] (the bracket
## that appears in several functions below). The (r, theta) block is
## diagonal, so g_rr = 1/g^rr and g_thth = 1/g^thth, giving
## d g_rr = -d g^rr / (g^rr)^2, and similarly for theta.
##
## Gamma_211 = d g_tt / dr
cdef double Gamma_211(double *m,
                      double *dr,
                      double *dt,
                      double det) noexcept nogil:

    return dr[4] / det - m[4]*(m[0]*dr[4] + m[4]*dr[0] - 2.0*m[1]*dr[1]) / (det*det)

## Gamma_222 = d g_rr / dr
cdef double Gamma_222(double *m,
                      double *dr,
                      double *dt,
                      double det) noexcept nogil:

    return -dr[2] / (m[2] * m[2])

## Gamma_233 = d g_thth / dr
cdef double Gamma_233(double * m, double * dr, double * dt, double det) noexcept nogil:

    return -dr[3] / (m[3]*m[3])

## Gamma_244 = d g_phph / dr
cdef double Gamma_244(double *m,
                      double *dr,
                      double *dt,
                      double det) noexcept nogil:

    return dr[0] / det - m[0]*(m[0]*dr[4] + m[4]*dr[0] - 2.0*m[1]*dr[1]) / (det*det)

## Gamma_241 = d g_phi t / dr
cdef double Gamma_241(double *m,
                      double *dr,
                      double *dt,
                      double det) noexcept nogil:

    return -dr[1] / det + m[1]*(m[0]*dr[4] + m[4]*dr[0] - 2.0*m[1]*dr[1]) / (det*det)

## Gamma_232 = d g_rr / d theta (for Gamma^r_rth)
cdef double Gamma_232(double *m,
                      double *dr,
                      double *dt,
                      double det) noexcept nogil:

    return -dt[2] / (m[2]*m[2])

## Gamma_311 = d g_tt / d theta
cdef double Gamma_311(double *m,
                      double *dr,
                      double *dt,
                      double det) noexcept nogil:

    return dt[4] / det - m[4]*(m[0]*dt[4] + m[4]*dt[0] - 2.0*m[1]*dt[1]) / (det*det)

## Gamma_322 = d g_rr / d theta (for Gamma^th_rr; same value as Gamma_232)
cdef double Gamma_322(double *m,
                      double *dr,
                      double *dt,
                      double det) noexcept nogil:

    return -dt[2] / (m[2] * m[2])

## Gamma_333 = d g_thth / d theta
cdef double Gamma_333(double *m,
                      double *dr,
                      double *dt,double det) noexcept nogil:

    return -dt[3] / (m[3] * m[3])

## Gamma_344 = d g_phph / d theta
cdef double Gamma_344(double *m,
                      double *dr,
                      double *dt,
                      double det) noexcept nogil:

    return dt[0] / det - m[0]*(m[0]*dt[4] + m[4]*dt[0] - 2.0*m[1]*dt[1]) / (det*det)

## Gamma_341 = d g_phi t / d theta
cdef double Gamma_341(double *m,
                      double *dr,
                      double *dt,
                      double det) noexcept nogil:

    return -dt[1] / det + m[1]*(m[0]*dt[4] + m[4]*dt[0] - 2.0*m[1]*dt[1]) / (det*det)

## Gamma_332 = d g_thth / dr (for Gamma^th_rth; same value as Gamma_233)
cdef double Gamma_332(double *m,
                      double *dr,
                      double *dt,
                      double det) noexcept nogil:

    return -dr[3] / (m[3] * m[3])

## RODES: right-hand side dy/dlambda of the photon equations of motion,
## in the signature required by GSL's gsl_odeiv2_system.
##
## State vector y[0..5] = (t, phi, r, dr/dlambda, theta, dtheta/dlambda).
## t and phi obey first-order equations from the conserved energy and
## angular momentum (PJ12 Eqs. 7-10); r and theta obey the second-order
## geodesic equations (PJ12 Eqs. 11-12), written as two first-order pairs.
##
## The ray is traced BACKWARDS from the image plane to the star, so the
## affine parameter runs opposite to the photon's motion. This flips the
## sign of dt/dlambda and dphi/dlambda relative to PJ12 Eqs. (9)-(10);
## the initial dr/dlambda and dtheta/dlambda are flipped to match in
## BOUNDARY_CONDITIONS.pyx.
##
## params[0..6] = (r_s, a, b, kappa, X_IP, Y_IP, inclination), where
## b = L/E is the photon's conserved impact parameter (PJ12, below Eq. 10),
## computed in BOUNDARY_CONDITIONS.pyx (RAY.IMPACT). X_IP, Y_IP and the
## inclination are used only to detect the special cases below.
cdef int RODES(double t,
                  const double *const y,
                  double *const dydl,
                  void *const params) noexcept nogil:

    ## y[] := y[t, phi, r, r', theta, theta'], with ' = d/dlambda

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

    ## Contravariant metric components g^uv (METRIC_qK.pyx),
    ## ordered (g^tt, g^tph, g^rr, g^thth, g^phph)
    cdef double metric[5]
    metric[0] = g_11(y[2], y[4], r_s, a, Sigma, Delta, kappa, sin_theta, func_theta, F1_r, func_r)
    metric[1] = g_14(y[2], y[4], r_s, a, Sigma, Delta)
    metric[2] = g_22(y[2], y[4], Sigma, Delta, kappa, func_theta, F1_r, func_r)
    metric[3] = g_33(y[2], y[4], Sigma, kappa, func_theta, F2_r)
    metric[4] = g_44(y[2], y[4], a, Sigma, Delta, kappa, sin_theta, func_theta, F2_r)

    ## Their partial derivatives with respect to r, same ordering
    cdef double deriv_r[5]
    deriv_r[0] = dg_11dr(y[2], y[4], r_s, a, Sigma, Delta, kappa, func_theta, F1_r, dF1dr_r, func_r, sin_theta)
    deriv_r[1] = dg_14dr(y[2], y[4], r_s, a, Sigma, Delta)
    deriv_r[2] = dg_22dr(y[2], y[4], r_s, Sigma, Delta, kappa, func_theta, F1_r, dF1dr_r, func_r) 
    deriv_r[3] = dg_33dr(y[2], y[4], Sigma, kappa, func_theta, F2_r, dF2dr_r) 
    deriv_r[4] = dg_44dr(y[2], y[4], r_s, a, Sigma, Delta, kappa, sin_theta, func_theta, F2_r, dF2dr_r)

    ## Their partial derivatives with respect to theta, same ordering
    cdef double deriv_theta[5]
    deriv_theta[0] = dg_11dtheta(y[2], y[4], r_s, a, Sigma, Delta, kappa, F1_r, func_r, sin_theta, cos_theta)
    deriv_theta[1] = dg_14dtheta(y[2], y[4], r_s, a, Sigma, Delta, sin_theta, cos_theta) 
    deriv_theta[2] = dg_22dtheta(y[2], y[4], a, Sigma, Delta, kappa, F1_r, func_r, sin_theta, cos_theta) 
    deriv_theta[3] = dg_33dtheta(y[2], y[4], a, Sigma, kappa, F2_r, sin_theta, cos_theta)
    deriv_theta[4] = dg_44dtheta(y[2], y[4], a, Sigma, Delta, kappa, F2_r, sin_theta, cos_theta, func_theta)

    ## Determinant of the contravariant (t, phi) block, used to invert it
    cdef double det_metric = det_g(y[2], y[4], r_s, a, Sigma, Delta, kappa, func_theta, F1_r, F2_r, func_r, sin_theta) ##!

    #printf("\n%.8f", det_metric)

    ## Christoffel symbols Gamma^r_uv, appearing in PJ12 Eq. (11).
    ## Gamma^r_uv = (1/2) g^rr (d_u g_rv + d_v g_ru - d_r g_uv), and g^rr is the
    ## only non-zero g^r. component, so each is +/- (1/2) g^rr times one
    ## covariant derivative from the Gamma_XYZ helpers:
    ##   Gamma_r[0] = Gamma^r_tt     = -(1/2) g^rr d_r g_tt
    ##   Gamma_r[1] = Gamma^r_rr     = +(1/2) g^rr d_r g_rr
    ##   Gamma_r[2] = Gamma^r_thth   = -(1/2) g^rr d_r g_thth
    ##   Gamma_r[3] = Gamma^r_phph   = -(1/2) g^rr d_r g_phph
    ##   Gamma_r[4] = Gamma^r_tph    = -(1/2) g^rr d_r g_tph
    ##   Gamma_r[5] = Gamma^r_rth    = +(1/2) g^rr d_th g_rr
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

    ## Christoffel symbols Gamma^theta_uv, appearing in PJ12 Eq. (12), same scheme:
    ##   Gamma_theta[0] = Gamma^th_tt   = -(1/2) g^thth d_th g_tt
    ##   Gamma_theta[1] = Gamma^th_rr   = -(1/2) g^thth d_th g_rr
    ##   Gamma_theta[2] = Gamma^th_thth = +(1/2) g^thth d_th g_thth
    ##   Gamma_theta[3] = Gamma^th_phph = -(1/2) g^thth d_th g_phph
    ##   Gamma_theta[4] = Gamma^th_tph  = -(1/2) g^thth d_th g_tph
    ##   Gamma_theta[5] = Gamma^th_rth  = +(1/2) g^thth d_r g_thth
    ## (The printf below labels its last entry Gamma_r[5]; harmless, it is commented out.)
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

    ## Covariant (t, phi) components from inverting the contravariant block:
    ## G_11 = g_tt, G_14 = g_tph, G_44 = g_phph.
    ## Note G_11*G_44 - G_14*G_14 = 1/det_metric.
    # Covariant metric components
    cdef:
        double G_11 = metric[4] / det_metric
        double G_14 = -metric[1] / det_metric
        double G_44 = metric[0] / det_metric

    ## dt/dlambda: PJ12 Eq. (9) with the overall sign flipped for backward tracing.
    dydl[0] = (G_44 + b*G_14) / (G_11*G_44 - G_14*G_14) # Reduce #flops here

    ## dphi/dlambda: PJ12 Eq. (10) with the sign flipped for backward tracing.
    ## The numerator here is b g_tt + g_tph, which is what solving PJ12
    ## Eqs. (7)-(8) gives. PJ12 Eq. (10) as printed has b g_tt + b g_tph, which
    ## appears to be a typo in the paper; the code is correct.
    ## Special case: for a = 0 and a ray in the x = 0 plane, the photon
    ## stays in a plane of constant phi, so dphi/dlambda is set exactly to 0.
    # Had only the second condition, thus phase-shifting for finite omega parameter...
    if are_equal(a, 0.0) and are_equal(X_IP, 0.0):
        dydl[1] = 0.0
    else:
        dydl[1] = -1.0 * (b*G_11 + G_14) / (G_11*G_44 - G_14*G_14)

    ## dr/dlambda (first-order form of PJ12 Eq. 11)
    dydl[2] = y[3]
    ## d^2r/dlambda^2: PJ12 Eq. (11). Terms in order: tt, rr, thth, phph, 2 ph t, 2 th r.
    dydl[3] = (-Gamma_r[0]*dydl[0]*dydl[0]
               - Gamma_r[1]*y[3]*y[3]
               - Gamma_r[2]*y[5]*y[5]
               - Gamma_r[3]*dydl[1]*dydl[1]
               - 2.0*Gamma_r[4]*dydl[0]*dydl[1]
               - 2.0*Gamma_r[5]*y[3]*y[5])
    ## dtheta/dlambda (first-order form of PJ12 Eq. 12)
    dydl[4] = y[5]

    # For the Schwarzschild metric, if rays are equatorial, or radial, one or
    # both of the angles are zero, leading to slow evolution because the
    # relative error is small but the derivative(s) are not set to precisely
    # zero.
    if are_equal(a, 0.0) and are_equal(X_IP, 0.0) and are_equal(Y_IP, 0.0):
        dydl[5] = 0.0
    elif are_equal(a, 0.0) and are_equal(Y_IP, 0.0) and are_equal(inclination, _pi/2.0):
        dydl[5] = 0.0
    else:
        ## d^2theta/dlambda^2: PJ12 Eq. (12), same term order as Eq. (11).
        dydl[5] = (-Gamma_theta[0]*dydl[0]*dydl[0]
                   - Gamma_theta[1]*y[3]*y[3]
                   - Gamma_theta[2]*y[5]*y[5]
                   - Gamma_theta[3]*dydl[1]*dydl[1]
                   - 2.0*Gamma_theta[4]*dydl[0]*dydl[1]
                   - 2.0*Gamma_theta[5]*y[3]*y[5])

    #printf("\ndydl: %.8f, %.8f, %.8f, %.8f, %.8f, %.8f",
           #dydl[0], dydl[1], dydl[2], dydl[3], dydl[4], dydl[5])

    return GSL_SUCCESS
