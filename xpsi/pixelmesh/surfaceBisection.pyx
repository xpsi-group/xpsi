#cython: cdivision=True
#cython: boundscheck=False
#cython: nonecheck=False
#cython: wraparound=False

## Stellar surface, ray-surface intersection, redshift and emission angle.
##
## Expressions are linked to:
##   AlGendy & Morsink 2014, ApJ 791, 78 (AM14) - surface shape
##   Cadeau et al. 2007, ApJ 654, 458 (C07) - redshift (Eq. 16) and emission
##                                            angle on an oblate star (Eqs. 18-20)
##   Psaltis & Johannsen 2012, ApJ 745, 1 (PJ12) - redshift (Eq. 24)
##   Bauböck et al. 2012, ApJ 753, 175 (B12) - Eq. 13, the same angle
## Conventions as in METRIC_qK.pyx (SI lengths; t, r, affine parameter in m).
##
## COORDINATES. The metric is in Boyer-Lindquist (BL) coordinates (r, theta).
## The surface is defined in 'spherical' coordinates (u, v) built from BL via
## the Kerr-Schild embedding (also coordinateTransformation.pyx):
##   u = r_sph = sqrt(r^2 + a^2 sin^2 theta),   u cos(v) = r cos(theta).
## AM14 fit the shape in the circumferential radius of the Butterworth-Ipser
## metric, not in this u, so this is an approximation (see the note on
## RADIUS below).
from libc.math cimport *
from libc.stdio cimport printf
from libc.math cimport sqrt, sin, cos, acos, pow
from xpsi.pixelmesh.METRIC_qK cimport *
from GSL cimport gsl_isnan
from xpsi.surface_radiation_field.effective_gravity_universal cimport oblateness_func_o2

cdef double _pi = M_PI
cdef double c = 2.99792458e8
cdef int SUCCESS = 1

cdef int MAX_BISECT_ITERATIONS = 100

## Surface shape R(theta)/R_eq = 1 + o2 cos^2(theta): AM14 Eq. (20).
## o2 = Omega-bar^2 * oblateness_func_o2(zeta, obl_surfgrav_ind), where
## obl_surfgrav_ind is the index of the oblateness relation chosen in Star
## (Spacetime.STAR_OBLATENESS_SURFGRAV_CANONICAL). For 'universal'
## (index 0) this is AM14 Eq. (21) with o20 = -0.788 and o21 = 1.030
## (AM14 Table 1); 'sphere' (index 1) gives o2 = 0. The other relations are
## in effective_gravity_universal.pyx.
## Arguments: x = cos(theta), epsilon = Omega-bar^2 = Omega^2 R_eq^3/(GM),
## zeta = compactness GM/(R_eq c^2), obl_surfgrav_ind = Spacetime.obl_surfgrav_ind.
##
## MODELLING CHOICE: WHICH RADIUS. AM14's R(theta) and R_eq are
## circumferential radii in the Butterworth-Ipser metric, and R_eq is what
## X-PSI takes as input everywhere else (e.g. the effective gravity in
## surface_radiation_field, which reads GEOM.R_eq). The ray tracer works in
## Boyer-Lindquist coordinates of the quasi-Kerr metric and imposes the
## shape on u = r_sph (see header). Since 2026-11 the surface is
##   r_sph(v) = R_eq_sph * RADIUS(cos v),
## where R_eq_sph (GEOM.R_eq_sph, set by SURFACE_SCALE below) is chosen so
## that the CIRCUMFERENTIAL equatorial radius, sqrt(g_phph) at theta = pi/2,
## equals the input R_eq exactly. GEOM.R_eq itself is left as the physical
## radius. Previously R_eq was used directly as the r_sph scale. For
## 1.4 Msun, R_eq = 12 km, 600 Hz (a = 651 m, kappa = 0.40) that made the
## circumferential equatorial radius 12026.9 m: 0.22% (27 m) too large,
## about 3% of the oblateness R_eq - R_pole = 969 m.
## What this does NOT do: match AM14's shape exactly away from the equator.
## That would need AM14's (Butterworth-Ipser) colatitude mapped onto
## Boyer-Lindquist theta, which differ at O(Omega-bar^2). The residual
## mismatch is expected to be a fraction of the 0.22% equatorial correction.
## To restore the old behaviour, set GEOM.R_eq_sph = GEOM.R_eq in
## integrator.pyx.
## B12 also impose the shape in 'spherical coordinates' without saying
## which radius they mean.
cdef double RADIUS(double x, double epsilon, double zeta, int obl_surfgrav_ind) noexcept nogil:

    return 1.0 + epsilon * oblateness_func_o2(zeta, obl_surfgrav_ind) * x * x

## Linear interpolation of x between (t1, x1) and (t2, x2)
cdef double LINTERP(double t,
                    double t1,
                    double t2,
                    double x1,
                    double x2) noexcept nogil:

    return x1*(t2 - t) / (t2 - t1) + x2*(t - t1) / (t2 - t1)

## Find where the ray crosses the surface, between the last state outside
## (y_p, at affine_p) and the first state inside (y, at affine).
## As in B12 Section 3: interpolate the state linearly in the affine
## parameter between the two Runge-Kutta steps, then bisect until
## |1 - R(theta_sph)/r_sph| < 1e-8 (at most MAX_BISECT_ITERATIONS steps).
## The R_eq passed in is the surface scale GEOM.R_eq_sph (see RADIUS).
## Accuracy: the interpolation follows the chord of the step, not the curved
## geodesic, so the crossing point carries an error of order the step's
## sagitta. In tests at the default tolerance (1e-12) the hit point agreed
## with a 1e-14 run to ~6e-7 R_eq in radius and ~2e-5 rad in angle, which
## includes this effect.
cdef int BISECT(const double *const y_p,
                const double *const y,
                double *const y_new,
                double affine_p,
                double affine,
                double a,
                double R_eq,
                double epsilon,
                double zeta,
                int obl_surfgrav_ind) noexcept nogil:

    cdef:
        int bisect = 1
        int iterations = 0
        size_t i
        double sin_theta, r_sph, cos_theta_sph, frac_diff
        double affine_bisect, affine_1 = affine_p, affine_2 = affine

    while bisect == 1:
        iterations += 1
        affine_bisect = 0.5 * (affine_1 + affine_2)

        for i in range(6):
            y_new[i] = LINTERP(affine_bisect, affine_p, affine, y_p[i], y[i])

        sin_theta = sin(y_new[4])
        r_sph = sqrt(y_new[2]*y_new[2] + a*a*sin_theta*sin_theta)
        cos_theta_sph = cos(y_new[4]) / sqrt(1.0 + a*a*sin_theta*sin_theta / (y_new[2]*y_new[2]))
        ## > 0 outside the surface, < 0 inside
        frac_diff = 1.0 - RADIUS(cos_theta_sph, epsilon, zeta, obl_surfgrav_ind) * R_eq / r_sph

        if fabs(frac_diff) < 1.0e-8:
            bisect = 0
        elif iterations >= MAX_BISECT_ITERATIONS:
            ## Safety cap (added 2026-11). The crossing is always bracketed, so
            ## bisection normally converges in ~30 iterations; reaching the cap
            ## means the tolerance cannot be met in double precision, and the
            ## current midpoint is kept.
            printf("\nWarning: surface bisection stopped after %d iterations; |frac_diff| = %.3e", iterations, fabs(frac_diff))
            bisect = 0
        else:
            if frac_diff > 0.0:
                affine_1 = affine_bisect
            else:
                affine_2 = affine_bisect

    return SUCCESS

cdef double SURFACE_SCALE(const _GEOM *const GEOM) noexcept nogil:

    ## Return R_eq_sph, the equatorial radius of the surface in r_sph, such
    ## that the circumferential equatorial radius sqrt(g_phph(r, pi/2))
    ## equals the physical input GEOM.R_eq (see the MODELLING CHOICE note on
    ## RADIUS). At the equator r_sph = sqrt(r^2 + a^2) with r the BL radius.
    ## Fixed-point iteration R_sph <- R_sph * R_eq / R_circ(R_sph); R_circ is
    ## very nearly proportional to R_sph, so this converges in a few steps.
    cdef:
        double R_sph = GEOM.R_eq
        double r, Sigma, Delta, F1_r, F2_r, func_r, det_metric, R_circ
        int k

    for k in range(50):
        r = sqrt(R_sph * R_sph - GEOM.asq)
        Sigma = r * r
        Delta = r * r - GEOM.r_s * r + GEOM.asq
        func_r = 1.0 - GEOM.r_s / r
        F1_r = F1(r, GEOM.r_s)
        F2_r = F2(r, GEOM.r_s)
        ## theta = pi/2: sin = 1, cos = 0, func_theta = 1
        det_metric = det_g(r, _pi/2.0, GEOM.r_s, GEOM.a, Sigma, Delta,
                           GEOM.kappa, 1.0, F1_r, F2_r, func_r, 1.0)
        ## covariant g_phph = g^tt / det (METRIC_qK.det_g)
        R_circ = sqrt(g_11(r, _pi/2.0, GEOM.r_s, GEOM.a, Sigma, Delta,
                           GEOM.kappa, 1.0, 1.0, F1_r, func_r) / det_metric)
        if fabs(R_circ / GEOM.R_eq - 1.0) < 1.0e-14:
            break
        R_sph *= GEOM.R_eq / R_circ

    return R_sph

## Redshift factor Z = 1 + z = E_emitted / E_observed, and the cosine of the
## emission angle ABB, measured in the comoving frame of the surface from
## the local surface normal, for a ray that has hit the star at state y.
## b is the ray's impact parameter (RAY.IMPACT, with E = 1).
## Results go into Z[0] and ABB[0].
cdef void ZABB(const _GEOM *const GEOM,
               const double *const y,
               double b,
               double *const Z,
               double *const ABB, int obl_surfgrav_ind) noexcept nogil:

    cdef:
        ## Position in BL (R, theta) and in spherical (u, v) coordinates
        double R = y[2] # BL
        double sin_theta = sin(y[4]) # BL
        double cos_theta = cos(y[4]) # BL
        double u = sqrt(y[2]*y[2] + GEOM.asq * sin_theta * sin_theta) # SPH
        double v = acos(cos_theta / sqrt(1.0 + GEOM.asq * sin_theta * sin_theta / (y[2]*y[2]))) # SPH
        double sin_v = sin(v)
        double cos_v = cos(v)

    cdef:
        ## Metric helper quantities (METRIC_qK.pyx)
        double func_theta = 1.0 - 3.0 * cos_theta * cos_theta
        double func_r = 1.0 - GEOM.r_s/y[2]
        double F1_r = F1(y[2], GEOM.r_s)
        double F2_r = F2(y[2], GEOM.r_s)
        double dF1dr_r = dF1dr(y[2], GEOM.r_s)
        double dF2dr_r = dF2dr(y[2], GEOM.r_s)
        double Sigma = y[2]*y[2] + GEOM.asq * cos_theta * cos_theta
        double Delta  = y[2]*y[2] - GEOM.r_s * y[2] + GEOM.asq

    ## Covariant components: G_11 = g_tt, G_14 = g_tph, G_44 = g_phph, G_22 = g_rr, G_33 = g_thth
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
    ## R'(v) = dR_surface/dv for r_sph(v) = R_eq_sph (1 + o2 cos^2 v)
    cdef double R_prime = -2.0 * GEOM.R_eq_sph * GEOM.epsilon * oblateness_func_o2(GEOM.zeta, obl_surfgrav_ind) * cos_v * sin_v

    ## (The original comment below is misleading: the code uses +b. k_t and
    ## k_phi are the FORWARD photon's contravariant components with E = 1,
    ## from PJ12 Eqs. (9)-(10), as in IS_NULL in RK_IP2S_tracer.pyx.)
    # Outgoing null geodesic has impact parameter -b, where b is that of the
    # null geodesic of the incoming photon
    cdef double k_t, k_phi

    k_t = (G_44 + b*G_14) / (-G_11*G_44 + G_14*G_14)
    k_phi = (b*G_11 + G_14) / (G_11*G_44 - G_14*G_14)

    #if k_t < 0.0:
    #    printf("\nWarning: -ve k_t = %.8e", k_t)

    ## Redshift. The surface rotates rigidly with angular velocity Omega, so its
    ## 4-velocity is u = gamma (1, 0, 0, Omega/c) with
    ## gamma = 1 / sqrt(-g_tt - 2 g_tph Omega/c - g_phph Omega^2/c^2).
    ## The numerator below is k_t + (Omega/c) k_phi (covariant) = -E + (Omega/c) L
    ## = -(1 - b Omega/c) for E = 1, so
    ##   Z = -k.u / E_obs = (1 - b Omega/c) / sqrt(...),
    ## the commented-out line. This is C07 Eq. (16) written for this metric, and
    ## the energy ratio of PJ12 Eq. (24) for a static observer at infinity.
    ## integrator.pyx then uses
    ## I_obs(E) = I_emit(E Z) / Z^3 (invariance of I/E^3).
    Z[0] = G_11 * k_t
    Z[0] += G_44 * GEOM.omega * k_phi / c
    Z[0] += G_14 * GEOM.omega * k_t / c
    Z[0] += G_14 * k_phi

    Z[0] /= -1.0 * sqrt(-G_11 - 2.0*G_14*GEOM.omega / c - G_44*GEOM.omega*GEOM.omega / (c*c))
    #Z[0] = (1.0 - b*GEOM.omega/c) / sqrt(-G_11 - 2.0*G_14*GEOM.omega / c - G_44*GEOM.omega*GEOM.omega / (c*c))

    ## Emission angle. The surface is u = R(v), so its normal covector in (u, v)
    ## is n_i = (1, -R'). The unit normal vector is n^i = g^ij n_j / |n|, using
    ## the full inverse of the (u, v) metric (see below).
    ## n has only (r, theta) components and u only (t, phi) components, and the
    ## metric has no (t,r), (t,theta), (phi,r), (phi,theta) terms, so n.u = 0
    ## exactly and the comoving-frame angle is cos(alpha) = k.n / (-k.u) with no
    ## separate aberration step. This is C07 Eq. (19): the angle between the
    ## photon and the normal after projecting both into the corotating frame
    ## (the projection does nothing to n, since n.u = 0). C07 Eq. (20) and B12
    ## Eq. (13) are explicit forms of the same thing; B12 Eq. (13) as printed
    ## omits the metric factor inside its square root, which enters only at
    ## O(Omega-bar^4).
    ##
    ## Jacobian of BL (R, theta) with respect to (u, v), from
    ## R^4 - R^2 (u^2 - a^2) - a^2 u^2 cos^2 v = 0 and R cos(theta) = u cos(v).
    ## (Checked against finite differences to 1e-9.)
    cdef double drdu, drdv, dtdu, dtdv, g_uu, g_vv, g_uv, det_uv, n_u, n_v, norm

    drdu = u / R + GEOM.asq * u * cos_v * cos_v / pow(R, 3.0)
    drdu /= 1.0 + GEOM.asq * u * u * cos_v * cos_v / pow(R, 4.0)

    drdv = -GEOM.asq * u * u * cos_v * sin_v / pow(R, 3.0)
    drdv /= 1.0 + GEOM.asq * u * u * cos_v * cos_v / pow(R, 4.0)

    dtdu = cos_v / R - u * cos_v * drdu / pow(R, 2.0)
    dtdu /= -sin_theta

    dtdv = -u * sin_v / R - u * cos_v * drdv / pow(R, 2.0)
    dtdv /= -sin_theta

    ## (u, v) metric components pulled back from the BL (r, theta) block,
    ## including the cross term g_uv, which is non-zero when a != 0. (Before
    ## 2026-11 g_uv was dropped and the (u, v) metric treated as diagonal; for
    ## the 600 Hz example above that changed cos(alpha) by up to ~7e-4.)
    g_uu = drdu * drdu * G_22 + dtdu * dtdu * G_33
    g_vv = dtdv * dtdv * G_33 + drdv * drdv * G_22
    g_uv = drdu * drdv * G_22 + dtdu * dtdv * G_33
    det_uv = g_uu * g_vv - g_uv * g_uv

    ## Contravariant normal n^i = g^ij n_j with n_j = (1, -R'), where
    ## g^uu = g_vv/det, g^uv = -g_uv/det, g^vv = g_uu/det; then normalise
    ## with |n|^2 = n^i n_i = n^u - R' n^v. Reduces to the old diagonal form
    ## when g_uv = 0.
    n_u = (g_vv + g_uv * R_prime) / det_uv
    n_v = (-g_uv - g_uu * R_prime) / det_uv
    norm = sqrt(n_u - R_prime * n_v)
    n_u /= norm
    n_v /= norm

    ## k.n = k_r n^r + k_theta n^theta, with k_r = g_rr dr/dlambda,
    ## k_theta = g_thth dtheta/dlambda, and n^r = (dr/du) n^u + (dr/dv) n^v,
    ## n^theta likewise, through the Jacobian. The state's dr/dlambda,
    ## dtheta/dlambda are for the
    ## backward-traced (ingoing) ray, hence the extra minus sign (the original
    ## comment's 'factor of pi': the angle is pi - alpha before the flip).
    ABB[0] = G_22 * y[3] * (drdu * n_u + drdv * n_v)
    ABB[0] += G_33 * y[5] * (dtdu * n_u + dtdv * n_v)
    ABB[0] /= -1.0 * Z[0]
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

    ## Guard against round-off at the limb (cos(alpha) -> 0) and at normal
    ## emission (cos(alpha) -> 1).
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
