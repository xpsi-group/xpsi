#cython: cdivision=True
#cython: boundscheck=False
#cython: nonecheck=False
#cython: wraparound=False

## The 'global ray map': trace a grid of rays over the image plane and
## record, for each pixel, where its ray meets the star and with what
## redshift, emission angle and time delay. integrator.pyx then turns these
## into intensities.
##
## Steps (compute_globalRayMap):
##   1. Screen: centre and semi-axes of the image (get_IP_radius.pyx).
##   2. Grid: an elliptical polar grid of N x N rays on that screen (init_RAY_MAP).
##   3. Centre ray: traced separately; it is the reference for time delays.
##   4. All rays, in parallel over rings: surface point (r, theta, phi),
##      redshift Z and cos(emission angle) (surfaceBisection.ZABB).
##   5. Time delays converted to phase lags relative to the reference ray.
## Afterwards integrator.pyx converts (r, theta) from Boyer-Lindquist to the
## spherical coordinates of the surface (coordinateTransformation.pyx).
##
## Units: image-plane coordinates in metres; t in the ray state is c times
## coordinate time (metres), converted to seconds here by dividing by c.
from cython.parallel cimport *
from libc.math cimport M_PI, sin, cos, log, pow, fabs
from libc.stdio cimport printf
from libc.stdlib cimport calloc, malloc, free

from GSL cimport (gsl_strerror,
                  gsl_spline_alloc,
                  gsl_interp_accel,
                  gsl_interp_accel_alloc,
                  gsl_spline,
                  gsl_spline2d,
                  gsl_interp_linear,
                  gsl_spline_init,
                  gsl_spline_free,
                  gsl_interp_accel_free,
                  gsl_interp_accel_reset,
                  gsl_isnan)

from xpsi.pixelmesh.RK_IP2S_tracer cimport RK
from xpsi.pixelmesh.surfaceBisection cimport ZABB
from xpsi.pixelmesh.get_IP_radius cimport compute_imagePlane_radius

from ..tools.core cimport are_equal

cdef int SUCCESS = 0

cdef double _2pi = 2.0 * M_PI
cdef double c = 2.99792458e8

## Arrays per pixel (N*N entries, INDEX = i*N + j for ring i and angle j):
##   X_MESH, Y_MESH  image-plane coordinates (m)
##   X_MODDED        X_MESH with |X| < 1e-4 m set to exactly 0, so that the
##                   X = 0 special cases in RODES_qK.pyx are triggered
##   RADIAL, THETA, PHI  surface point of the ray (BL, later spherical)
##   LAG             time delay, then phase lag (see below)
##   Z, ABB          redshift factor 1+z and cos(emission angle);
##                   -1 and -2 for rays that miss the star
##   XI              per-pixel accuracy diagnostic: the ray's last null-check
##                   value (see below)
## ORIGIN[0..5] holds the same quantities for the centre ray:
## (lag, phi, r, theta, Z, ABB).
cdef RAY_MAP* alloc_RAY_MAP(size_t NGR) noexcept nogil:

    cdef RAY_MAP *MAP = <RAY_MAP*> malloc(sizeof(RAY_MAP))

    MAP.numRays = NGR

    MAP.r_MESH = <double*> malloc((NGR + 1) * sizeof(double))
    MAP.t_MESH = <double*> malloc((NGR + 1) * sizeof(double))

    MAP.X_MESH = <double*> malloc(NGR * NGR * sizeof(double))
    MAP.X_MODDED = <double*> malloc(NGR * NGR * sizeof(double))
    MAP.Y_MESH = <double*> malloc(NGR * NGR * sizeof(double))
    MAP.RADIAL = <double*> malloc(NGR * NGR * sizeof(double))
    MAP.THETA = <double*> malloc(NGR * NGR * sizeof(double))
    MAP.PHI = <double*> malloc(NGR * NGR * sizeof(double))
    MAP.LAG = <double*> malloc(NGR * NGR * sizeof(double))
    MAP.ABB = <double*> malloc(NGR * NGR * sizeof(double))
    MAP.Z = <double*> malloc(NGR * NGR * sizeof(double))
    MAP.XI = <double*> malloc(NGR * NGR * sizeof(double))

    return MAP

cdef void free_RAY_MAP(RAY_MAP *const MAP) noexcept nogil:

    free(MAP.r_MESH)
    free(MAP.t_MESH)

    free(MAP.X_MESH)
    free(MAP.Y_MESH)
    free(MAP.X_MODDED) # for coordinates forced to zero
    free(MAP.RADIAL)
    free(MAP.THETA)
    free(MAP.PHI)
    free(MAP.LAG)
    free(MAP.ABB)
    free(MAP.Z)
    free(MAP.XI)

    free(<RAY_MAP*> MAP)

## Elliptical polar grid: ring i = 1..N at fractional radius
##   r_i = (i/N)^p,  p = radialIncrementExponent,
## and N azimuths t_j = 2 pi [ (N mod 4 - 2)/(4N) + j/N ], j = 0..N-1;
## pixel (i, j) is at
##   X = ORIGIN_X + SEMI_MAJOR r_i cos t_j,   Y = ORIGIN_Y + SEMI_MINOR r_i sin t_j.
## With the default p = 1/2 (PhotospherePlotter.image) r_i^2 is linear in
## i, so every annulus between successive rings has the same area. The
## centre (r = 0) is the separately traced origin ray; the outer ring r_N = 1
## lies on the screen edge, outside the star by the 0.5% margin, so it
## always misses.
## The azimuthal offset keeps t = pi/2 off the grid for every N. Which of
## the other axis directions fall on the grid depends on N mod 4 (checked
## for N = 4..40): N = 0 mod 4 (e.g. the default N = 100): none of
## 0, pi/2, pi, 3 pi/2; N = 2 mod 4: t = 0 and pi; N odd: t = 3 pi/2.
cdef void init_RAY_MAP(RAY_MAP *const MAP,
                       double radialIncrementExponent) noexcept nogil:

    cdef:
        size_t i, j
        double temp = 0.0
        double r_increment = 1.0 / pow(<double>(MAP.numRays), radialIncrementExponent)
        double t_increment = 1.0 / <double>(MAP.numRays)
        size_t INDEX = 0

    MAP.r_MESH[0] = 0.0
    # initialise so that no value of t is pi/2
    MAP.t_MESH[0] = (<double>(MAP.numRays%4) - 2.0) / (4.0 * <double>(MAP.numRays))

    for i in range(1, MAP.numRays + 1):
        temp += 1.0
        MAP.r_MESH[i] = pow(temp, radialIncrementExponent) * r_increment
        #printf("r: %.8e", MAP.r_MESH[i])
        MAP.t_MESH[i] = MAP.t_MESH[i - 1] + t_increment

    for i in range(0, MAP.numRays + 1):
        MAP.t_MESH[i] *= _2pi

    for i in range(MAP.numRays):
        for j in range(MAP.numRays):
            MAP.X_MESH[INDEX] = MAP.SEMI_MAJOR * MAP.r_MESH[i + 1] * cos(MAP.t_MESH[j])
            MAP.X_MESH[INDEX] += MAP.ORIGIN_X
            if fabs(MAP.X_MESH[INDEX]) < 1.0e-4:
                MAP.X_MODDED[INDEX] = 0.0
            else:
                MAP.X_MODDED[INDEX] = MAP.X_MESH[INDEX]
            MAP.Y_MESH[INDEX] = MAP.SEMI_MINOR * MAP.r_MESH[i + 1] * sin(MAP.t_MESH[j])
            MAP.Y_MESH[INDEX] += MAP.ORIGIN_Y
            INDEX += 1

cdef int compute_globalRayMap(size_t numThreads,
                              double radialIncrementExponent,
                              const _GEOM *const GEOM,
                              RAY_MAP *const MAP,
                              _RAY **RAYS,
                              int force_circular, int obl_surfgrav_ind) noexcept nogil:

    cdef:
        signed int ii
        size_t i, j, T
        int status
        size_t N = MAP.numRays
        size_t ROOT, INDEX

    compute_imagePlane_radius(GEOM, RAYS[0], MAP, force_circular, obl_surfgrav_ind)

    init_RAY_MAP(MAP, radialIncrementExponent)

    ## Centre ray: the reference for time delays. NOTE: ZABB is called on it
    ## without checking that it hit the star; if it failed, ORIGIN holds
    ## failure flags.
    RAYS[0].X_IP = MAP.ORIGIN_X
    RAYS[0].Y_IP = MAP.ORIGIN_Y
    # Compute centre ray; (X,Y)=(0,0) shifted
    status = RK(RAYS[0], GEOM, obl_surfgrav_ind)

    MAP.ORIGIN[0] = -1.0 * RAYS[0].STATE[0] / c
    MAP.ORIGIN[1] = RAYS[0].STATE[1]
    MAP.ORIGIN[2] = RAYS[0].STATE[2]
    MAP.ORIGIN[3] = RAYS[0].STATE[4]

    ZABB(GEOM,
         RAYS[0].STATE,
         RAYS[0].IMPACT,
         MAP.ORIGIN + 4,
         MAP.ORIGIN + 5,
         obl_surfgrav_ind)

    #printf("\nOrigin (Z, ABB): (%.8e, %.8e)", MAP.ORIGIN[4], MAP.ORIGIN[5])

    ## Reference travel time. In Schwarzschild (a = kappa = 0, i.e. a
    ## non-rotating, spherical star) use the exact coordinate time of a radial
    ## null ray from r = d to r = R_eq:
    ##   t = [ (d - R_eq) + r_s ln((d - r_s)/(R_eq - r_s)) ] / c.
    ## Otherwise use the numerically traced centre ray, whose integration error
    ## then largely cancels against that of the other rays.
    ## (Check: for a 12 km, 1.4 Msun star the numerical centre ray
    ## differs from this formula by 6.8 m (2.3e-8 s) at tolerance 1e-12 and
    ## 1.4 m at 1e-13, i.e. ~1e-4 rad of phase at 600 Hz: negligible.)
    if are_equal(GEOM.a, 0.0) and are_equal(GEOM.kappa, 0.0):
        # calibration for Schwarzschild equatorial b=0 ray
        MAP.refRayTime = GEOM.d - GEOM.R_eq
        MAP.refRayTime += GEOM.r_s * log((GEOM.d - GEOM.r_s) / (GEOM.R_eq - GEOM.r_s))
        MAP.refRayTime /= c
    else: # more general lag calibration
        MAP.refRayTime = MAP.ORIGIN[0]

    ## Trace all rays, rings distributed over threads (each thread has its own
    ## _RAY workspace RAYS[T]; every pixel writes only its own INDEX, so there
    ## are no races). Rays that miss or fail (STATE[2] <= 0) get Z = -1,
    ## ABB = -2 and are skipped by the integrator.
    for ii in prange(<signed int>N,
                     nogil = True,
                     schedule = 'static',
                     num_threads = numThreads,
                     chunksize = 1):
        i = <size_t> ii
        T = threadid()
        ROOT = i * N
        if i%100 == 0:
            printf("\nThread %d is tracing annulus #%d of rays.", <int>T, <int>i)

        for j in range(N):
            INDEX = ROOT + j

            RAYS[T].X_IP = MAP.X_MODDED[INDEX]
            RAYS[T].Y_IP = MAP.Y_MESH[INDEX]

            status = RK(RAYS[T], GEOM, obl_surfgrav_ind)

            ## The ray is traced backwards, so t decreases from 0 at the image plane;
            ## -t/c is the photon's travel time (s) from the surface to the image plane.
            MAP.LAG[INDEX] = -1.0 * RAYS[T].STATE[0] / c
            MAP.PHI[INDEX] = RAYS[T].STATE[1]
            MAP.RADIAL[INDEX] = RAYS[T].STATE[2]
            MAP.THETA[INDEX] = RAYS[T].STATE[4]
            ## Accuracy diagnostic: the last value of the null
            ## check for this ray. Normally the norm g_uv k^u k^v from IS_NULL
            ## (RK_IP2S_tracer.pyx) at the ray's last checked step, which should
            ## be ~0 (the tracer rejects steps with |norm| > 1e-6). It is NOT
            ## updated after a polar-axis fallback, and for a ray rejected at
            ## the image plane it is xi + 1 from BOUNDARY_CONDITIONS.pyx.
            ## Values far from 0 flag pixels where the integration struggled.
            MAP.XI[INDEX] = RAYS[T].XI

            if RAYS[T].STATE[2] > 0.0:
                ZABB(GEOM,
                     RAYS[T].STATE,
                     RAYS[T].IMPACT,
                     MAP.Z + INDEX,
                     MAP.ABB + INDEX, obl_surfgrav_ind)
            else:
                (MAP.Z + INDEX)[0] = -1.0
                (MAP.ABB + INDEX)[0] = -2.0

    #cdef size_t i_min, j_min

    # for i in range(MAP.numRays):
    #     for j in range(MAP.numRays):
    #         if lag[i,j] < refRayTime and lag[i,j] > 0.0:
    #             refRayTime = lag[i,j]
    #             i_min = i
    #             j_min = j

    ## Convert travel times to phase lags relative to the reference ray:
    ##   lag = (travel time - reference travel time) * mode_omega,
    ## where mode_omega = 2 pi * mode_frequency (normally the spin frequency;
    ## a different frequency can be used for an oscillation mode). A pixel
    ## whose photons take longer to arrive is seen at a later rotational phase.
    MAP.ORIGIN[0] -= MAP.refRayTime
    MAP.ORIGIN[0] *= GEOM.mode_omega # use mode frequency instead

    INDEX = 0
    for i in range(N):
        for j in range(N):
            MAP.LAG[INDEX] -= MAP.refRayTime
            MAP.LAG[INDEX] *= GEOM.mode_omega # use mode frequency instead
            INDEX += 1

    return SUCCESS
