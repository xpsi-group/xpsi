#cython: cdivision=True
#cython: boundscheck=False
#cython: nonecheck=False
#cython: wraparound=False

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

cdef _RAY* alloc_RAY(double epsabs_ray,
                     double epsrel_ray,
                     double INIT_STEP,
                     size_t MAXSTEPS) nogil:

    cdef _RAY *RAY = <_RAY*> malloc(sizeof(_RAY))

    RAY.INIT_STEP = INIT_STEP
    RAY.MAXSTEPS = MAXSTEPS

    RAY.S = gsl_odeiv2_step_alloc(gsl_odeiv2_step_rkck, 6)
    RAY.C = gsl_odeiv2_control_standard_new(epsabs_ray,
                                            epsrel_ray,
                                            1.0,
                                            1.0)
    RAY.E = gsl_odeiv2_evolve_alloc(6)

    RAY.SYS.function = RODES
    RAY.SYS.jacobian = NULL
    RAY.SYS.dimension = 6
    RAY.MAX_AFFINE = 1.0e15

    return RAY

cdef void reset_RAY(_RAY *const RAY,
                    const _GEOM *const GEOM) nogil:

    RAY.NUMSTEPS = 0
    RAY.EVOLVE = 1
    RAY.AFFINE = 0.0
    RAY.STEP = RAY.INIT_STEP
    gsl_odeiv2_step_reset(RAY.S)
    gsl_odeiv2_evolve_reset(RAY.E)
    RAY.NUM_SINGULARITIES = 0

    # compute ray boundary conditions at image plane
    COMPUTE_BCs(RAY, GEOM)

    if fabs(RAY.XI) > 1.0e-8:
        RAY.EVOLVE = 0
        RAY.STATE[0] = -104.0
        RAY.STATE[1] = -104.0
        RAY.STATE[2] = -104.0
        RAY.STATE[3] = -104.0
        RAY.STATE[4] = -104.0
        RAY.STATE[5] = -104.0
        printf("\nError: Tangent 4-vector at domain boundary is non-null.")

    RAY.PARAMS[0] = GEOM.r_s
    RAY.PARAMS[1] = GEOM.a
    RAY.PARAMS[2] = RAY.IMPACT
    RAY.PARAMS[3] = GEOM.kappa
    RAY.PARAMS[4] = RAY.X_IP
    RAY.PARAMS[5] = RAY.Y_IP
    RAY.PARAMS[6] = GEOM.inclination

    RAY.SYS.params = &RAY.PARAMS

cdef void free_RAY(_RAY *const RAY) nogil:

    gsl_odeiv2_step_free(RAY.S)
    gsl_odeiv2_control_free(RAY.C)
    gsl_odeiv2_evolve_free(RAY.E)

    free(<_RAY*> RAY)

cdef int RK(_RAY *const RAY,
            const _GEOM *const GEOM) nogil:

    reset_RAY(RAY, GEOM)

    cdef:
        int status
        double r_sph, cos_theta_sph
        double rsq, sin_theta_sq
        double PROXIMITY_4_COORD_TRANSFORM = 10.0 * GEOM.R_eq
        double BYPASS_THRESHOLD = 10.0 * GEOM.R_eq
        double X_IP = RAY.X_IP / GEOM.b_max
        double Y_IP = RAY.Y_IP / GEOM.b_max
        double temp

    while RAY.EVOLVE == 1:
        RAY.PREVIOUS_AFFINE = RAY.AFFINE
        RAY.PREVIOUS[0] = RAY.STATE[0]
        RAY.PREVIOUS[1] = RAY.STATE[1]
        RAY.PREVIOUS[2] = RAY.STATE[2]
        RAY.PREVIOUS[3] = RAY.STATE[3]
        RAY.PREVIOUS[4] = RAY.STATE[4]
        RAY.PREVIOUS[5] = RAY.STATE[5]
        RAY.PREVIOUS_XI = RAY.XI
        RAY.PREVIOUS_STEP = RAY.STEP

        status = gsl_odeiv2_evolve_apply(RAY.E,
                                         RAY.C,
                                         RAY.S,
                                         &RAY.SYS,
                                         &RAY.AFFINE,
                                         RAY.MAX_AFFINE,
                                         &RAY.STEP,
                                         RAY.STATE)


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

            if <unsigned int> RAY.NUMSTEPS >= RAY.MAXSTEPS:
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

            elif RAY.STATE[2] < PROXIMITY_4_COORD_TRANSFORM:
                sin_theta_sq = pow(sin(RAY.STATE[4]), 2.0)
                rsq = RAY.STATE[2] * RAY.STATE[2]
                r_sph = sqrt(rsq + GEOM.asq * sin_theta_sq)
                cos_theta_sph = cos(RAY.STATE[4]) / sqrt(1.0 + GEOM.asq * sin_theta_sq / rsq)

                #printf("\nx: %.8e; y: %.8e", X_IP, Y_IP)

                if r_sph / GEOM.R_eq <= RADIUS(cos_theta_sph, GEOM.epsilon, GEOM.zeta):
                    RAY.EVOLVE = 0
                    status = BISECT(RAY.PREVIOUS,
                                    RAY.STATE,
                                    RAY.BISECT_STATE,
                                    RAY.PREVIOUS_AFFINE,
                                    RAY.AFFINE,
                                    GEOM.a,
                                    GEOM.R_eq,
                                    GEOM.epsilon,
                                    GEOM.zeta)
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

        elif status != GSL_SUCCESS: #and RAY.NUM_SINGULARITIES < 1000:
            # Integrating towards polar coordinate singularity...
            #printf("\n\n----------------------------------->>> Evolving new ray...")
            printf("\nPolar singularity encountered...")

            RAY.NUM_SINGULARITIES += 1
            RAY.NUM_SINGULARITY_STEPS = 0

            RODES(RAY.AFFINE, RAY.STATE, RAY.dSTATEdAFFINE, RAY.PARAMS)

            # Chan et al. (2013) use forward Euler...
            temp =  RAY.dSTATEdAFFINE[0] / RAY.STATE[0]
            #temp += RAY.dSTATEdAFFINE[1] / RAY.STATE[1]
            temp += RAY.dSTATEdAFFINE[2] / RAY.STATE[2]
            temp += RAY.dSTATEdAFFINE[3] / RAY.STATE[3]
            temp += RAY.dSTATEdAFFINE[4] / RAY.STATE[4]
            temp += RAY.dSTATEdAFFINE[5] / RAY.STATE[5]

            RAY.SINGULARITY_AFFINE_STEP_SIZE = 1.0 / temp

            while 1:
                RAY.NUM_SINGULARITY_STEPS += 1
                #printf("\nSingularity steps: %d", RAY.NUM_SINGULARITY_STEPS)

                if <unsigned int> RAY.NUM_SINGULARITY_STEPS == RAY.MAXSTEPS:
                    RAY.EVOLVE = 0
                    break

                RAY.AFFINE = RAY.PREVIOUS_AFFINE
                RAY.STATE[0] = RAY.PREVIOUS[0]
                RAY.STATE[1] = RAY.PREVIOUS[1]
                RAY.STATE[2] = RAY.PREVIOUS[2]
                RAY.STATE[3] = RAY.PREVIOUS[3]
                RAY.STATE[4] = RAY.PREVIOUS[4]
                RAY.STATE[5] = RAY.PREVIOUS[5]

                RAY.STATE[0] += RAY.SINGULARITY_AFFINE_STEP_SIZE * RAY.dSTATEdAFFINE[0]
                #RAY.STATE[1] += RAY.SINGULARITY_AFFINE_STEP_SIZE * RAY.dSTATEdAFFINE[1]
                RAY.STATE[2] += RAY.SINGULARITY_AFFINE_STEP_SIZE * RAY.dSTATEdAFFINE[2]
                RAY.STATE[3] += RAY.SINGULARITY_AFFINE_STEP_SIZE * RAY.dSTATEdAFFINE[3]
                RAY.STATE[4] += RAY.SINGULARITY_AFFINE_STEP_SIZE * RAY.dSTATEdAFFINE[4]
                RAY.STATE[5] += RAY.SINGULARITY_AFFINE_STEP_SIZE * RAY.dSTATEdAFFINE[5]

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
                if RAY.NUM_SINGULARITIES < 2:
                    RAY.STATE[1] = RAY.PREVIOUS[1] + _pi

                IS_NULL(RAY.STATE,
                        &RAY.XI,
                        &RAY.IMPACT,
                        GEOM.r_s,
                        GEOM.kappa,
                        GEOM.a,
                        GEOM.asq)

                if (RAY.XI - RAY.PREVIOUS_XI) / RAY.PREVIOUS_XI > 1.0e-3:
                    RAY.SINGULARITY_AFFINE_STEP_SIZE *= 0.1
                elif RAY.XI > 1.0e-3:
                    RAY.SINGULARITY_AFFINE_STEP_SIZE *= 0.1
                else:
                    RAY.AFFINE += RAY.SINGULARITY_AFFINE_STEP_SIZE
                    break

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

    return status

cdef int IS_NULL(const double *const y,
                 double *const xi,
                 double *const b,
                 double r_s,
                 double kappa,
                 double a,
                 double asq) nogil:

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
