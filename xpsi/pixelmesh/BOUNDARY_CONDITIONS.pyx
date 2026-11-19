#cython: cdivision=True
#cython: boundscheck=False
#cython: nonecheck=False
#cython: wraparound=False

## Initial conditions for a photon leaving the image plane of a distant
## observer, for backward ray tracing in the quasi-Kerr metric.
##
## Expressions are linked to equations in:
##   Johannsen & Psaltis 2010, ApJ 718, 446, Appendix (JP10)
##   Psaltis & Johannsen 2012, ApJ 745, 1 (PJ12)
## Conventions (units, index labels 1-4 = t, r, theta, phi, the helper
## quantities Sigma, Delta, func_theta, func_r) are as in METRIC_qK.pyx.
##
## Geometry (JP10 Fig. 1 and Eq. A1; PJ12 Fig. 1): the image plane is a
## distance d from the star, at inclination i to the spin axis, centred
## on phi = 0. A pixel has image-plane coordinates (X_IP, Y_IP), which
## are (x', y') in JP10 and (alpha_0, beta_0) in PJ12. In the code the
## plane sits at d = 10^4 b_max (integrator.pyx), where b_max is the
## image-plane radius of the star (get_IP_radius.pyx).
##
## The calculation is done in long double because d/r_i -> 1: PJ12 (below
## Eq. 20) warn that the initial momentum is prone to round-off error.
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
    long double atan2l(long double y, long double x)
    long double sqrtl(long double x)

from xpsi.pixelmesh.METRIC_qK cimport *

cdef double _pi = M_PI

cdef void COMPUTE_BCs(_RAY *const RAY,
                      const _GEOM *const GEOM) noexcept nogil:

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

    ## Y[0..5] = (t, phi, r, dr/dlambda, theta, dtheta/dlambda), the state
    ## vector integrated in RODES_qK.pyx.
    ##
    ## t = 0 at the image plane (arbitrary origin; time delays are measured
    ## relative to a reference ray in globalRayMap.pyx).
    Y[0] = 0.0

    ## phi_i: JP10 Eq. (A3); PJ12 Eq. (17), tan(phi_i) = X / (d sin i - Y cos i).
    ## atan2 is used so that phi_i lands in the correct quadrant for any sign
    ## of the denominator (the pixel's x-coordinate, JP10 Eq. A1). For
    ## inclinations above ~b_max/d = 1e-4 rad the denominator is always
    ## positive and atan2 gives exactly the same result as atan(X/(...)).
    ## Near-pole-on inclinations are additionally clamped in integrator.pyx,
    ## because the central pixel would otherwise start on the spin axis.
    Y[1] = atan2l(X_IP, dsin - Ycos)

    ## r_i: JP10 Eq. (A3); PJ12 Eq. (15).
    Y[2] = sqrtl(dsq + xsq + ysq)
    ## dr/dlambda = -d/r_i: JP10 Eq. (A4) with k_0 = 1. The minus sign is because
    ## the photon is traced backwards, towards the star. PJ12 Eq. (18) gives
    ## +d/r_i for the forward-moving photon.
    Y[3] = -d / sqrtl(dsq + xsq + ysq)
    ## theta_i: JP10 Eq. (A3); PJ12 Eq. (16).
    Y[4] = acosl((dcos + Ysin) / sqrtl(dsq + xsq + ysq))

    ## dtheta/dlambda: JP10 Eq. (A4) with k_0 = 1, traced backwards (so minus
    ## PJ12 Eq. 19). The denominator is the PJ12 form
    ## [r_i^2 - (d cos i + Y sin i)^2]^(1/2) multiplied out, which equals
    ## [X^2 + (d sin i - Y cos i)^2]^(1/2) and avoids the cancellation PJ12
    ## warn about.
    Y[5] = (cos_inclination - (dsq*cos_inclination + Y_IP*dsin) / (dsq + xsq + ysq)) / sqrtl(xsq + (dsin - Ycos)*(dsin - Ycos))

    ## Metric helper quantities at the initial position (see METRIC_qK.pyx)
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

    ## Determinant of the contravariant (t, phi) block, used to invert it
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

    ## Covariant components from the contravariant ones (METRIC_qK.det_g):
    ##   G_11 = g_tt     = g^phph / det
    ##   G_14 = g_tph    = -g^tph / det
    ##   G_44 = g_phph   = g^tt / det
    ##   G_22 = g_rr     = 1 / g^rr
    ##   G_33 = g_thth   = 1 / g^thth
    ## (G_22 and G_33 are initialised to 1 and then divided below.)
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

    ## Remaining momentum components and the conserved impact parameter.
    ##
    ## SIGN CONVENTION. The state vector carries the BACKWARD-traced
    ## dr/dlambda and dtheta/dlambda (Y[3], Y[5] above). k_phi and k_t below
    ## instead belong to the FORWARD-moving photon (travelling towards the
    ## observer):
    ##   - k_phi = -X sin(i) / [(d sin i - Y cos i)^2 + X^2] is PJ12 Eq. (20),
    ##     i.e. minus JP10 Eq. (A4).
    ##   - k_t is the positive (future-directed) root of the null condition
    ##     PJ12 Eq. (13). Only squares of dr/dlambda and dtheta/dlambda enter
    ##     that condition, so mixing directions there is harmless.
    ## This is deliberate: the impact parameter b is computed from the forward
    ## photon, and RODES_qK.pyx then flips the signs of dt/dlambda and
    ## dphi/dlambda for backward tracing. At the image plane, RODES's
    ## dphi/dlambda equals JP10's k^phi, as it should.
    ##
    ## k_t: solving g_tt k_t^2 + 2 g_tph k_t k_phi + (g_rr k_r^2 + g_thth k_th^2
    ## + g_phph k_phi^2) = 0, with temp = g_tph k_phi / g_tt.
    # Require that the photon wavevector is null
    cdef long double k_phi_l = -1.0 * Xsin / ((dsin - Ycos)*(dsin - Ycos) + xsq)
    cdef long double temp_l = G_14_l * k_phi_l / G_11_l
    cdef long double k_t_l = -temp_l + sqrtl(temp_l*temp_l - (G_22_l*Y[3]*Y[3] + G_33_l*Y[5]*Y[5] + G_44_l*k_phi_l*k_phi_l) / G_11_l)

    #printf("+ve root: %.8e", <double>k_t_l)
    #printf("-ve root: %.8e", <double>k_t_l)

    #printf("constant: %.8e", <double>(G_11_l * k_t_l + G_14_l * k_phi_l))

    ## Impact parameter b = L/E (PJ12, below Eq. 10), where
    ## L = k_phi (covariant) = g_phph k^phi + g_tph k^t and E = -k_t (covariant).
    ## The code stores L and does not divide by E. That is safe because, with
    ## |k_spatial| = 1 at d = 10^4 b_max, E = 1 to about 1 part in 10^11.
    ## (Checked numerically for a = 0.15 r_s, kappa = 0.2 over the image plane.)
    ## RODES_qK.pyx and IS_NULL in RK_IP2S_tracer.pyx both assume E = 1.
    ##
    ## XI: the null-check ratio of PJ12 Eq. (14), which should equal -1.
    ## The code stores xi + 1, so a perfectly null vector gives RAY.XI = 0.
    ## RK_IP2S_tracer.pyx rejects the ray if |RAY.XI| > 1e-8 here.
    RAY.IMPACT = <double>(G_44_l * k_phi_l + G_14_l * k_t_l)
    RAY.XI = <double>((G_22_l*Y[3]*Y[3] + G_33_l*Y[5]*Y[5] + G_44_l*k_phi_l*k_phi_l + 2.0*G_14_l*k_t_l*k_phi_l) / (G_11_l*k_t_l*k_t_l))
    RAY.XI += 1.0
    ## Hand the initial state to the integrator (in double precision).
    RAY.STATE[0] = <double>Y[0]
    RAY.STATE[1] = <double>Y[1]
    RAY.STATE[2] = <double>Y[2]
    RAY.STATE[3] = <double>Y[3]
    RAY.STATE[4] = <double>Y[4]
    RAY.STATE[5] = <double>Y[5]
