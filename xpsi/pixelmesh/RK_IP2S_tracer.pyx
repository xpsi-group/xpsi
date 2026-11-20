#cython: cdivision=True
#cython: boundscheck=False
#cython: nonecheck=False
#cython: wraparound=False

## Backward ray tracer: integrates one photon from the distant image plane
## until it hits the stellar surface, misses the star, or fails.
##
## Expressions and methods are linked to:
##   Psaltis & Johannsen 2012, ApJ 745, 1 (PJ12) - scheme and null check
##   Chan, Psaltis & Ozel 2013, ApJ 777, 13 (C13) - handling of the polar axis
##
## Flow for each ray (RK below):
##   1. reset_RAY: initial conditions (BOUNDARY_CONDITIONS.pyx) and a check
##      that the initial wavevector is null.
##   2. Adaptive Runge-Kutta-Cash-Karp steps of the equations in RODES_qK.pyx.
##   3. After each successful step: reflect theta back into [0, pi] if the
##      ray crossed the axis; then check, in order, max steps, null
##      condition, escape, and intersection with the stellar surface.
##   4. If GSL cannot take a step (in practice near the spin axis), fall
##      back to a forward-Euler step across the axis.
##
## On exit, STATE holds the ray at the surface (or a failure flag):
##   -102  maximum number of steps exceeded
##   -103  ray escaped: it missed the star
##   -104  initial wavevector not null (set in reset_RAY)
##   -105  ray could not be carried across the spin axis (see the fallback)
## Downstream (globalRayMap.pyx) any ray with STATE[2] > 0 is treated as
## having hit the surface, so every failure must set a negative flag.
##
## Units: SI lengths throughout (r, t and the affine parameter in metres;
## t is c times coordinate time). See METRIC_qK.pyx for conventions.
cimport cython
from cython.parallel cimport *
from libc.math cimport M_PI, fabs, sqrt, sin, cos, acos, pow
from libc.stdio cimport printf
from libc.stdlib cimport malloc, free
from GSL cimport GSL_SUCCESS

from xpsi.pixelmesh.surfaceBisection cimport RADIUS, BISECT
from xpsi.pixelmesh.METRIC_qK cimport *
from xpsi.pixelmesh.RODES_qK cimport RODES
from xpsi.pixelmesh.BOUNDARY_CONDITIONS cimport COMPUTE_BCs

from GSL cimport (gsl_strerror,
                  gsl_isnan,
                  gsl_odeiv2_system,
                  gsl_odeiv2_step,
                  gsl_odeiv2_control,
                  gsl_odeiv2_evolve,
                  gsl_odeiv2_step_free,
                  gsl_odeiv2_control_free,
                  gsl_odeiv2_evolve_free,
                  gsl_odeiv2_step_reset,
                  gsl_odeiv2_evolve_reset,
                  gsl_odeiv2_step_alloc,
                  gsl_odeiv2_step_rkck,
                  gsl_odeiv2_control_y_new,
                  gsl_odeiv2_control_standard_new,
                  gsl_odeiv2_evolve_alloc,
                  gsl_odeiv2_control_init,
                  gsl_odeiv2_evolve_apply,
                  gsl_odeiv2_evolve_apply_fixed_step)

cdef double _pi = M_PI
cdef double _2pi = 2.0 * M_PI
cdef double c = 2.99792458e8
cdef double G = 6.6730831e-11
cdef int SUCCESS = 0
cdef int ERROR = 1
## Maximum number of polar-axis fallbacks per ray (see RK). A ray that
## crosses the axis normally needs one; repeated GSL failures mean the
## fallback is not helping, and the ray is terminated.
cdef int MAX_SINGULARITIES = 10

cdef inline double _log_rate(double dy, double y) noexcept nogil:
    if y != 0.0:
        return fabs(dy / y)
    return 0.0

cdef _RAY* alloc_RAY(double epsabs_ray,
                     double epsrel_ray,
                     double INIT_STEP,
                     size_t MAXSTEPS) noexcept nogil:

    cdef _RAY *RAY = <_RAY*> malloc(sizeof(_RAY))

    RAY.INIT_STEP = INIT_STEP
    RAY.MAXSTEPS = MAXSTEPS

    ## Embedded Runge-Kutta-Cash-Karp (4,5) stepper with adaptive step size.
    ## PJ12 (Section 3) used a fixed-order 4th-order RK with their own step
    ## rule (Eq. 21); here step control is left to GSL instead.
    RAY.S = gsl_odeiv2_step_alloc(gsl_odeiv2_step_rkck, 6)
    ## Standard GSL error control with a_y = a_dydt = 1: the allowed local
    ## error in component i is epsabs + epsrel * (|y_i| + h |dy_i/dlambda|).
    ## The state mixes metres (t, r) and radians (phi, theta), so in practice
    ## epsrel sets the accuracy of t and r, and epsabs that of the angles.
    RAY.C = gsl_odeiv2_control_standard_new(epsabs_ray,
                                            epsrel_ray,
                                            1.0,
                                            1.0)
    RAY.E = gsl_odeiv2_evolve_alloc(6)

    ## Right-hand side of PJ12 Eqs. (9)-(12): see RODES_qK.pyx
    RAY.SYS.function = RODES
    RAY.SYS.jacobian = NULL
    RAY.SYS.dimension = 6
    ## Upper limit on the affine parameter (metres). The image plane is at
    ## d = 10^4 b_max ~ 10^8 m, so this is never reached in practice.
    RAY.MAX_AFFINE = 1.0e15

    return RAY

cdef void reset_RAY(_RAY *const RAY,
                    const _GEOM *const GEOM) noexcept nogil:

    RAY.NUMSTEPS = 0
    RAY.EVOLVE = 1
    RAY.AFFINE = 0.0
    RAY.STEP = RAY.INIT_STEP
    gsl_odeiv2_step_reset(RAY.S)
    gsl_odeiv2_evolve_reset(RAY.E)
    RAY.NUM_SINGULARITIES = 0

    # compute ray boundary conditions at image plane
    ## Initial state, impact parameter b (RAY.IMPACT) and RAY.XI = xi + 1,
    ## where xi is the null-check ratio of PJ12 Eq. (14); see BOUNDARY_CONDITIONS.pyx.
    COMPUTE_BCs(RAY, GEOM)

    ## Reject the ray if the initial wavevector is not null (flag -104).
    if fabs(RAY.XI) > 1.0e-8:
        RAY.EVOLVE = 0
        RAY.STATE[0] = -104.0
        RAY.STATE[1] = -104.0
        RAY.STATE[2] = -104.0
        RAY.STATE[3] = -104.0
        RAY.STATE[4] = -104.0
        RAY.STATE[5] = -104.0
        printf("\nError: Tangent 4-vector at domain boundary is non-null.")

    ## Parameters passed to RODES (see RODES_qK.pyx)
    RAY.PARAMS[0] = GEOM.r_s
    RAY.PARAMS[1] = GEOM.a
    RAY.PARAMS[2] = RAY.IMPACT
    RAY.PARAMS[3] = GEOM.kappa
    RAY.PARAMS[4] = RAY.X_IP
    RAY.PARAMS[5] = RAY.Y_IP
    RAY.PARAMS[6] = GEOM.inclination

    RAY.SYS.params = &RAY.PARAMS

cdef void free_RAY(_RAY *const RAY) noexcept nogil:

    gsl_odeiv2_step_free(RAY.S)
    gsl_odeiv2_control_free(RAY.C)
    gsl_odeiv2_evolve_free(RAY.E)

    free(<_RAY*> RAY)

cdef int RK(_RAY *const RAY,
            const _GEOM *const GEOM, int obl_surfgrav_ind) noexcept nogil:

    reset_RAY(RAY, GEOM)

    cdef:
        int status
        double r_sph, cos_theta_sph
        double rsq, sin_theta_sq
        ## Only test for the surface once the ray is within 10 R_eq (saves the
        ## coordinate transformation on the long approach from the image plane).
        double PROXIMITY_4_COORD_TRANSFORM = 10.0 * GEOM.R_eq
        ## A ray that came within 10 R_eq and then moves back out past it is taken
        ## to have missed the star. Every pixel has impact parameter <= b_max
        ## ~ 1.3 R_eq, so every ray does get inside this radius.
        double BYPASS_THRESHOLD = 10.0 * GEOM.R_eq
        ## Image-plane coordinates in units of b_max (used only in diagnostics)
        double X_IP = RAY.X_IP / GEOM.b_max
        double Y_IP = RAY.Y_IP / GEOM.b_max
        double temp

    while RAY.EVOLVE == 1:
        ## Save the state so that a step can be undone.
        RAY.PREVIOUS_AFFINE = RAY.AFFINE
        RAY.PREVIOUS[0] = RAY.STATE[0]
        RAY.PREVIOUS[1] = RAY.STATE[1]
        RAY.PREVIOUS[2] = RAY.STATE[2]
        RAY.PREVIOUS[3] = RAY.STATE[3]
        RAY.PREVIOUS[4] = RAY.STATE[4]
        RAY.PREVIOUS[5] = RAY.STATE[5]
        RAY.PREVIOUS_XI = RAY.XI
        RAY.PREVIOUS_STEP = RAY.STEP

        ## One adaptive step. RODES always returns GSL_SUCCESS, so a failure here
        ## means GSL could not find an acceptable step. In practice this happens
        ## near the spin axis, where dphi/dlambda ~ 1/sin^2(theta) grows without
        ## bound; it is handled in the 'status != GSL_SUCCESS' branch below.
        status = gsl_odeiv2_evolve_apply(RAY.E,
                                         RAY.C,
                                         RAY.S,
                                         &RAY.SYS,
                                         &RAY.AFFINE,
                                         RAY.MAX_AFFINE,
                                         &RAY.STEP,
                                         RAY.STATE)


        ## Polar-axis crossing: theta has left [0, pi]. Map it back with
        ## theta -> -theta (or 2 pi - theta), reverse dtheta/dlambda, and
        ## phi -> phi + pi, as in C13 Eq. (3). The phi shift is skipped if phi
        ## already jumped by more than pi/2 during the step (the integrator
        ## followed the rapid swing in phi near the axis itself).
        ## Note: dtheta/dlambda is forced positive after crossing theta = 0 and
        ## negative after crossing theta = pi, which equals reversing it whenever
        ## the ray was heading towards the axis.
        if RAY.STATE[4] <= 0.0:
            RAY.STATE[4] *= -1.0
            if RAY.STATE[5] < 0.0:
                RAY.STATE[5] *= -1.0
            if fabs(RAY.PREVIOUS[1] - RAY.STATE[1]) < _pi/2.0:
                RAY.STATE[1] += _pi
        elif RAY.STATE[4] >= _pi:
            RAY.STATE[4] = _2pi - RAY.STATE[4]
            if RAY.STATE[5] > 0.0:
                RAY.STATE[5] *= -1.0
            if fabs(RAY.PREVIOUS[1] - RAY.STATE[1]) < _pi/2.0:
                RAY.STATE[1] += _pi

        if status == GSL_SUCCESS:
            RAY.NUMSTEPS += 1

            ## Too many steps: give up on this ray (flag -102).
            if <size_t> RAY.NUMSTEPS >= RAY.MAXSTEPS:
                printf("\n\nSteps exceeded... %d", RAY.NUMSTEPS)
                printf("\nx: %.8e; y: %.8e", X_IP, Y_IP)
                printf("\ny[]: %.8e, %.8e, %.8e, %.8e, %.8e, %.8e",
                        RAY.STATE[0]/c,
                        RAY.STATE[1]/_2pi,
                        RAY.STATE[2]/GEOM.R_eq,
                        RAY.STATE[3]/GEOM.R_eq,
                        RAY.STATE[4]/_pi,
                        RAY.STATE[5]/_pi)
                RAY.EVOLVE = 0
                RAY.STATE[0] = -102.0
                RAY.STATE[1] = -102.0
                RAY.STATE[2] = -102.0
                RAY.STATE[3] = -102.0
                RAY.STATE[4] = -102.0
                RAY.STATE[5] = -102.0

            ## Null check (PJ12 Eq. 13; C13 Section 2). If the wavevector has drifted
            ## from null, undo the step and retry with a step 10 times smaller.
            ## NUMSTEPS still counts the rejected step, so this cannot loop forever.
            ## The check is switched off for the rest of the ray once it has crossed
            ## the axis via the forward-Euler fallback (NUM_SINGULARITIES >= 1).
            elif (RAY.NUM_SINGULARITIES < 1 and IS_NULL(RAY.STATE,
                                                           &RAY.XI,
                                                           &RAY.IMPACT,
                                                           GEOM.r_s,
                                                           GEOM.kappa,
                                                           GEOM.a,
                                                           GEOM.asq) == 0):
                # Periodically check that tangent 4-vector is null
                printf("\n\nTangent 4-vector has become non-null due to precision loss.")
                printf("\nx: %.8e; y: %.8e", X_IP, Y_IP)
                printf("\ny[]: %.8e, %.8e, %.8e, %.8e, %.8e, %.8e",
                        RAY.STATE[0]/c,
                        RAY.STATE[1]/_2pi,
                        RAY.STATE[2]/GEOM.R_eq,
                        RAY.STATE[3]/GEOM.R_eq,
                        RAY.STATE[4]/_pi,
                        RAY.STATE[5]/_pi)

                RAY.AFFINE = RAY.PREVIOUS_AFFINE
                RAY.STATE[0] = RAY.PREVIOUS[0]
                RAY.STATE[1] = RAY.PREVIOUS[1]
                RAY.STATE[2] = RAY.PREVIOUS[2]
                RAY.STATE[3] = RAY.PREVIOUS[3]
                RAY.STATE[4] = RAY.PREVIOUS[4]
                RAY.STATE[5] = RAY.PREVIOUS[5]
                RAY.XI       = RAY.PREVIOUS_XI

                RAY.STEP *= 0.1

                gsl_odeiv2_step_reset(RAY.S)
                gsl_odeiv2_evolve_reset(RAY.E)

            ## Ray has turned round and is heading back out: it missed the star (flag -103).
            elif RAY.PREVIOUS[2] < BYPASS_THRESHOLD and RAY.STATE[2] >= BYPASS_THRESHOLD:
                # Photon is likely on a scattering trajectory
                #printf("\nBYPASS")

                #printf("\nx: %.8e; y: %.8e", X_IP, Y_IP)
                RAY.EVOLVE = 0

                if RAY.NUM_SINGULARITIES > 0:
                    printf("\n\nRay encountered polar singularity and scattered...: %d", RAY.NUM_SINGULARITIES)
                    printf("\nx: %.8e; y: %.8e", X_IP, Y_IP)
                    printf("\ny[]: %.8e, %.8e, %.8e, %.8e, %.8e, %.8e",
                            RAY.STATE[0]/c,
                            RAY.STATE[1]/_2pi,
                            RAY.STATE[2]/GEOM.R_eq,
                            RAY.STATE[3]/GEOM.R_eq,
                            RAY.STATE[4]/_pi,
                            RAY.STATE[5]/_pi)

                RAY.STATE[0] = -103.0
                RAY.STATE[1] = -103.0
                RAY.STATE[2] = -103.0
                RAY.STATE[3] = -103.0
                RAY.STATE[4] = -103.0
                RAY.STATE[5] = -103.0

            ## Surface test. The surface shape R(theta) (AlGendy & Morsink 2014, see
            ## surfaceBisection.RADIUS) is applied in spherical coordinates built from
            ## Boyer-Lindquist (r, theta) via the Kerr-Schild embedding:
            ##   r_sph = sqrt(r^2 + a^2 sin^2 theta),
            ##   cos(theta_sph) = cos(theta) / sqrt(1 + a^2 sin^2(theta) / r^2)
            ## (same transformation as coordinateTransformation.pyx).
            ## The surface is r_sph = R_eq_sph * RADIUS(cos theta_sph), with
            ## R_eq_sph set so that the circumferential equatorial radius is
            ## R_eq (modelling choice documented at surfaceBisection.RADIUS).
            ## If the ray is now inside the surface, BISECT finds the crossing point
            ## between the previous and current states and the ray stops there.
            ## Caveat: only the end of each step is tested, so a step whose chord
            ## enters and leaves the star (a grazing ray at the limb) is not caught.
            ## This shifts the limb by at most ~h^2/(8R) for a step of length h.
            elif RAY.STATE[2] < PROXIMITY_4_COORD_TRANSFORM:
                sin_theta_sq = pow(sin(RAY.STATE[4]), 2.0)
                rsq = RAY.STATE[2] * RAY.STATE[2]
                r_sph = sqrt(rsq + GEOM.asq * sin_theta_sq)
                cos_theta_sph = cos(RAY.STATE[4]) / sqrt(1.0 + GEOM.asq * sin_theta_sq / rsq)

                #printf("\nx: %.8e; y: %.8e", X_IP, Y_IP)

                if r_sph / GEOM.R_eq_sph <= RADIUS(cos_theta_sph, GEOM.epsilon, GEOM.zeta, obl_surfgrav_ind):
                    RAY.EVOLVE = 0
                    status = BISECT(RAY.PREVIOUS,
                                    RAY.STATE,
                                    RAY.BISECT_STATE,
                                    RAY.PREVIOUS_AFFINE,
                                    RAY.AFFINE,
                                    GEOM.a,
                                    GEOM.R_eq_sph,
                                    GEOM.epsilon,
                                    GEOM.zeta,
                                    obl_surfgrav_ind)
                    RAY.STATE[0] = RAY.BISECT_STATE[0]
                    RAY.STATE[1] = RAY.BISECT_STATE[1]
                    RAY.STATE[2] = RAY.BISECT_STATE[2]
                    RAY.STATE[3] = RAY.BISECT_STATE[3]
                    RAY.STATE[4] = RAY.BISECT_STATE[4]
                    RAY.STATE[5] = RAY.BISECT_STATE[5]

                    if RAY.NUM_SINGULARITIES > 0:
                        printf("\n\nRay encountered polar singularity...: %d", RAY.NUM_SINGULARITIES)
                        printf("\nx: %.8e; y: %.8e", X_IP, Y_IP)
                        printf("\ny[]: %.8e, %.8e, %.8e, %.8e, %.8e, %.8e",
                                RAY.STATE[0]/c,
                                RAY.STATE[1]/_2pi,
                                RAY.STATE[2]/GEOM.R_eq,
                                RAY.STATE[3]/GEOM.R_eq,
                                RAY.STATE[4]/_pi,
                                RAY.STATE[5]/_pi)

        ## Fallback when GSL cannot step (near the spin axis). Following the idea
        ## of C13 (Section 2, Fig. 2), take a single first-order forward-Euler
        ## step from the last good state so that the ray jumps over the axis,
        ## then hand back to the adaptive Runge-Kutta stepper.
        ##
        ## Differences from C13: C13 take an Euler step of Delta_lambda/9 when
        ## |xi + 1| > ~1e-3, where Delta_lambda is their step rule (C13 Eq. 2):
        ## Delta / (|d ln r/dlambda| + |dtheta/dlambda| + |dphi/dlambda|) with
        ## Delta ~ 1/32. Here the first trial step is 1 / sum |dy_i/dlambda / y_i|
        ## over the five non-phi state variables (phi is excluded because
        ## dphi/dlambda diverges at the axis and is set by hand below). It is
        ## deliberately large, so that the step can carry the ray across the
        ## axis, and is then cut by factors of 10 until |norm| <= 1e-3, the same
        ## tolerance C13 use.
        ##
        ## At most MAX_SINGULARITIES fallbacks are allowed per ray; after that
        ## the ray is terminated with flag -105 (final 'else' branch).
        elif status != GSL_SUCCESS and RAY.NUM_SINGULARITIES < MAX_SINGULARITIES:
            # Integrating towards polar coordinate singularity...
            #printf("\n\n----------------------------------->>> Evolving new ray...")
            printf("\nPolar singularity encountered...")

            RAY.NUM_SINGULARITIES += 1
            RAY.NUM_SINGULARITY_STEPS = 0

            RODES(RAY.AFFINE, RAY.STATE, RAY.dSTATEdAFFINE, RAY.PARAMS)

            # Chan et al. (2013) use forward Euler...
            ## Trial step from the absolute logarithmic rates of change of t, r,
            ## dr/dlambda, theta and dtheta/dlambda (absolute values as in C13 Eq. 2,
            ## so the terms cannot cancel or give a negative step). _log_rate skips
            ## any component that is exactly zero. If every term is zero, fall back
            ## to the last adaptive step size.
            temp =  _log_rate(RAY.dSTATEdAFFINE[0], RAY.STATE[0])
            #temp += _log_rate(RAY.dSTATEdAFFINE[1], RAY.STATE[1])
            temp += _log_rate(RAY.dSTATEdAFFINE[2], RAY.STATE[2])
            temp += _log_rate(RAY.dSTATEdAFFINE[3], RAY.STATE[3])
            temp += _log_rate(RAY.dSTATEdAFFINE[4], RAY.STATE[4])
            temp += _log_rate(RAY.dSTATEdAFFINE[5], RAY.STATE[5])

            if temp > 0.0:
                RAY.SINGULARITY_AFFINE_STEP_SIZE = 1.0 / temp
            ## Fallback ran out of steps (NUM_SINGULARITY_STEPS == MAXSTEPS):
            ## terminate the ray and flag it (-105) so it is not shaded as a surface hit.
            else:
                RAY.SINGULARITY_AFFINE_STEP_SIZE = RAY.PREVIOUS_STEP

            while 1:
                RAY.NUM_SINGULARITY_STEPS += 1
                #printf("\nSingularity steps: %d", RAY.NUM_SINGULARITY_STEPS)

                if <size_t>RAY.NUM_SINGULARITY_STEPS == RAY.MAXSTEPS:
                    RAY.EVOLVE = 0
                    break

                RAY.AFFINE = RAY.PREVIOUS_AFFINE
                RAY.STATE[0] = RAY.PREVIOUS[0]
                RAY.STATE[1] = RAY.PREVIOUS[1]
                RAY.STATE[2] = RAY.PREVIOUS[2]
                RAY.STATE[3] = RAY.PREVIOUS[3]
                RAY.STATE[4] = RAY.PREVIOUS[4]
                RAY.STATE[5] = RAY.PREVIOUS[5]

                ## Forward-Euler step for everything except phi
                RAY.STATE[0] += RAY.SINGULARITY_AFFINE_STEP_SIZE * RAY.dSTATEdAFFINE[0]
                #RAY.STATE[1] += RAY.SINGULARITY_AFFINE_STEP_SIZE * RAY.dSTATEdAFFINE[1]
                RAY.STATE[2] += RAY.SINGULARITY_AFFINE_STEP_SIZE * RAY.dSTATEdAFFINE[2]
                RAY.STATE[3] += RAY.SINGULARITY_AFFINE_STEP_SIZE * RAY.dSTATEdAFFINE[3]
                RAY.STATE[4] += RAY.SINGULARITY_AFFINE_STEP_SIZE * RAY.dSTATEdAFFINE[4]
                RAY.STATE[5] += RAY.SINGULARITY_AFFINE_STEP_SIZE * RAY.dSTATEdAFFINE[5]

                ## Map theta back into [0, pi] as in the main loop (C13 Eq. 3)
                if RAY.STATE[4] <= 0.0:
                    RAY.STATE[4] *= -1.0
                    if RAY.STATE[5] < 0.0:
                        RAY.STATE[5] *= -1.0
                elif RAY.STATE[4] >= _pi:
                    RAY.STATE[4] = _2pi - RAY.STATE[4]
                    if RAY.STATE[5] > 0.0:
                        RAY.STATE[5] *= -1.0

                # only pass pole first time the RK evolution failed
                # this does not account for two or more interesctions of a
                # ray with poles at affine values with finite separation
                ## Crossing the axis shifts phi by pi (C13 Eq. 3). As the original comment
                ## says, this is only applied the first time the fallback is used.
                if RAY.NUM_SINGULARITIES < 2:
                    RAY.STATE[1] = RAY.PREVIOUS[1] + _pi

                ## Accept the Euler step if |g_uv k^u k^v| <= 1e-3 (RAY.XI is set to the
                ## norm by IS_NULL; the test uses its absolute value because the norm can
                ## have either sign). Otherwise shrink the step by 10 and retry from the
                ## saved state. The 1e-3 tolerance matches C13 Section 2.
                IS_NULL(RAY.STATE,
                        &RAY.XI,
                        &RAY.IMPACT,
                        GEOM.r_s,
                        GEOM.kappa,
                        GEOM.a,
                        GEOM.asq)

                if fabs(RAY.XI) > 1.0e-3:
                    RAY.SINGULARITY_AFFINE_STEP_SIZE *= 0.1
                else:
                    RAY.AFFINE += RAY.SINGULARITY_AFFINE_STEP_SIZE
                    break

            ## Diagnostics printed after a successful axis crossing; then resume Runge-Kutta stepping.
            if RAY.EVOLVE == 1:
                gsl_odeiv2_step_reset(RAY.S)
                gsl_odeiv2_evolve_reset(RAY.E)

                printf("\n\nRay state upon encountering polar singularity:")
                printf("\nx: %.8e; y: %.8e", X_IP, Y_IP)
                printf("\ny[]: %.8e, %.8e, %.8e, %.8e, %.8e, %.8e",
                        RAY.PREVIOUS[0]/c,
                        RAY.PREVIOUS[1]/_2pi,
                        RAY.PREVIOUS[2]/GEOM.R_eq,
                        RAY.PREVIOUS[3]/GEOM.R_eq,
                        RAY.PREVIOUS[4]/_pi,
                        RAY.PREVIOUS[5]/_pi)

                RODES(RAY.AFFINE, RAY.STATE, RAY.dSTATEdAFFINE, RAY.PARAMS)

                printf("\n\nTangent 4-vector and derivatives upon encountering polar singularity:")
                printf("\ny[]: %.8e, %.8e, %.8e, %.8e, %.8e, %.8e",
                        RAY.dSTATEdAFFINE[0]/c,
                        RAY.dSTATEdAFFINE[1]/_2pi,
                        RAY.dSTATEdAFFINE[2]/GEOM.R_eq,
                        RAY.dSTATEdAFFINE[3]/GEOM.R_eq,
                        RAY.dSTATEdAFFINE[4]/_pi,
                        RAY.dSTATEdAFFINE[5]/_pi)

                RAY.NUMSTEPS += 1
                #RAY.STEP = RAY.INIT_STEP
                printf("\n\nPolar singularity handled.")
                printf("\nRay state upon bypassing polar singularity:")
                printf("\nx: %.8e; y: %.8e", X_IP, Y_IP)
                printf("\ny[]: %.8e, %.8e, %.8e, %.8e, %.8e, %.8e",
                    RAY.STATE[0]/c,
                    RAY.STATE[1]/_2pi,
                    RAY.STATE[2]/GEOM.R_eq,
                    RAY.STATE[3]/GEOM.R_eq,
                    RAY.STATE[4]/_pi,
                    RAY.STATE[5]/_pi)

                RODES(RAY.AFFINE, RAY.STATE, RAY.dSTATEdAFFINE, RAY.PARAMS)

                printf("\n\nTangent 4-vector and derivatives upon bypassing polar singularity:")
                printf("\ny[]: %.8e, %.8e, %.8e, %.8e, %.8e, %.8e",
                        RAY.dSTATEdAFFINE[0]/c,
                        RAY.dSTATEdAFFINE[1]/_2pi,
                        RAY.dSTATEdAFFINE[2]/GEOM.R_eq,
                        RAY.dSTATEdAFFINE[3]/GEOM.R_eq,
                        RAY.dSTATEdAFFINE[4]/_pi,
                        RAY.dSTATEdAFFINE[5]/_pi)
            else:
                printf("\nCannot handle polar singularity... terminating ray.")
                RAY.STATE[0] = -105.0
                RAY.STATE[1] = -105.0
                RAY.STATE[2] = -105.0
                RAY.STATE[3] = -105.0
                RAY.STATE[4] = -105.0
                RAY.STATE[5] = -105.0
        ## GSL failed again after MAX_SINGULARITIES fallbacks: terminate the ray
        ## and flag it (-105) so it is not shaded as a surface hit. (Before
        ## 2026-11 this branch was unreachable, because the fallback condition did
        ## not include the NUM_SINGULARITIES limit.)
        else:
            printf("\n\nMultiple singularities encountered... terminating ray.")
            printf("\n\nRay state upon re-encountering polar singularity:")
            printf("\nx: %.8e; y: %.8e", X_IP, Y_IP)
            printf("\ny[]: %.8e, %.8e, %.8e, %.8e, %.8e, %.8e",
                    RAY.STATE[0]/c,
                    RAY.STATE[1]/_2pi,
                    RAY.STATE[2]/GEOM.R_eq,
                    RAY.STATE[3]/GEOM.R_eq,
                    RAY.STATE[4]/_pi,
                    RAY.STATE[5]/_pi)
            RODES(RAY.AFFINE, RAY.STATE, RAY.dSTATEdAFFINE, RAY.PARAMS)
            printf("\n\nTangent 4-vector and derivatives upon re-encountering polar singularity:")
            printf("\ny[]: %.8e, %.8e, %.8e, %.8e, %.8e, %.8e",
                    RAY.dSTATEdAFFINE[0]/c,
                    RAY.dSTATEdAFFINE[1]/_2pi,
                    RAY.dSTATEdAFFINE[2]/GEOM.R_eq,
                    RAY.dSTATEdAFFINE[3]/GEOM.R_eq,
                    RAY.dSTATEdAFFINE[4]/_pi,
                    RAY.dSTATEdAFFINE[5]/_pi)
            RAY.EVOLVE = 0
            RAY.STATE[0] = -105.0
            RAY.STATE[1] = -105.0
            RAY.STATE[2] = -105.0
            RAY.STATE[3] = -105.0
            RAY.STATE[4] = -105.0
            RAY.STATE[5] = -105.0

    return status

## Null check during integration (PJ12 Eq. 13).
##
## k^t and k^phi are rebuilt from the impact parameter b using PJ12
## Eqs. (9)-(10) (with the corrected numerator b g_tt + g_tph), i.e. for
## the forward photon with energy E = 1. That normalisation matches the
## initial conditions, where E = 1 to ~1e-11 (BOUNDARY_CONDITIONS.pyx).
## dr/dlambda and dtheta/dlambda come from the state; only their squares
## enter, so the backward-tracing sign does not matter.
##
## xi[0] is set to the NORM g_uv k^u k^v itself, which should be 0. This
## is a different quantity from RAY.XI set in BOUNDARY_CONDITIONS.pyx
## (xi + 1 with xi the ratio of PJ12 Eq. 14): xi + 1 = norm / (g_tt (k^t)^2).
## Since g_tt < 0 that factor is negative and O(1), so the two have
## opposite signs and similar size. Both vanish for a null vector.
## Returns 1 if |norm| <= 1e-6 (null to tolerance), else 0.
##
## Covariant components as in METRIC_qK.det_g: G_11 = g_tt, G_14 = g_tph,
## G_44 = g_phph, G_22 = g_rr, G_33 = g_thth. (G_14 is not in the cdef
## list below; Cython infers it as a double.)
cdef int IS_NULL(const double *const y,
                 double *const xi,
                 double *const b,
                 double r_s,
                 double kappa,
                 double a,
                 double asq) noexcept nogil:

    cdef:
        double Sigma, Delta, sin_theta, cos_theta
        double func_theta, func_r, F1_r, F2_r, dF1dr_r, dF2dr_r
        double det_metric, G_11, G_22, G_33, G_44, k_phi, k_t

    sin_theta = sin(y[4])
    cos_theta = cos(y[4])
    Sigma = y[2] * y[2] + asq * cos_theta * cos_theta
    Delta = y[2] * y[2] - r_s * y[2] + asq
    func_theta = 1.0 - 3.0 * cos_theta * cos_theta
    func_r = 1.0 - r_s / y[2]
    F1_r = F1(y[2], r_s)
    F2_r = F2(y[2], r_s)
    dF1dr_r = dF1dr(y[2], r_s)
    dF2dr_r = dF2dr(y[2], r_s)

    det_metric = det_g(y[2], y[4], r_s, a, Sigma, Delta, kappa, func_theta, F1_r, F2_r, func_r, sin_theta)

    G_11 = g_44(y[2], y[4], a, Sigma, Delta, kappa, sin_theta, func_theta, F2_r) / det_metric
    G_14 = -g_14(y[2], y[4], r_s, a, Sigma, Delta) / det_metric
    G_22 = 1.0 / g_22(y[2], y[4], Sigma, Delta, kappa, func_theta, F1_r, func_r)
    G_33 = 1.0 / g_33(y[2], y[4], Sigma, kappa, func_theta, F2_r)
    G_44 = g_11(y[2], y[4], r_s, a, Sigma, Delta, kappa, sin_theta, func_theta, F1_r, func_r) / det_metric

    ## k^t and k^phi from PJ12 Eqs. (9) and (10) with E = 1 (forward photon)
    # signs of both can be inverted without changing norm
    k_t = -1.0 * (G_44 + b[0]*G_14) / (G_11*G_44 - G_14*G_14)
    k_phi = (b[0]*G_11 + G_14) / (G_11*G_44 - G_14*G_14)

    xi[0] = (G_11 * k_t * k_t
             + G_22 * y[3] * y[3]
             + G_33 * y[5] * y[5]
             + G_44 * k_phi * k_phi
             + 2.0 * G_14 * k_t * k_phi)

    if fabs(xi[0]) > 1.0e-6:
        #printf("\n\nTangent 4-vector has become non-null due to precision loss: %.8e", fabs(xi[0]))
        return 0
    else:
        #printf("\nxi: %.8e", xi[0])
        return 1
