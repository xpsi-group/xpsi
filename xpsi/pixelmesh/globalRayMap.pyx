#cython: cdivision=True
#cython: boundscheck=False
#cython: nonecheck=False
#cython: wraparound=False

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

cdef int SUCCESS = 0

cdef double _2pi = 2.0 * M_PI
cdef double c = 2.99792458e8

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
                              int force_circular) noexcept nogil:

    cdef:
        signed int ii
        size_t i, j, T
        int status
        size_t N = MAP.numRays
        size_t ROOT, INDEX

    compute_imagePlane_radius(GEOM, RAYS[0], MAP, force_circular)

    init_RAY_MAP(MAP, radialIncrementExponent)

    RAYS[0].X_IP = MAP.ORIGIN_X
    RAYS[0].Y_IP = MAP.ORIGIN_Y
    # Compute centre ray; (X,Y)=(0,0) shifted
    status = RK(RAYS[0], GEOM)

    MAP.ORIGIN[0] = -1.0 * RAYS[0].STATE[0] / c
    MAP.ORIGIN[1] = RAYS[0].STATE[1]
    MAP.ORIGIN[2] = RAYS[0].STATE[2]
    MAP.ORIGIN[3] = RAYS[0].STATE[4]

    ZABB(GEOM,
         RAYS[0].STATE,
         RAYS[0].IMPACT,
         MAP.ORIGIN + 4,
         MAP.ORIGIN + 5)

    #printf("\nOrigin (Z, ABB): (%.8e, %.8e)", MAP.ORIGIN[4], MAP.ORIGIN[5])

    if GEOM.a == 0.0 and GEOM.kappa == 0.0:
        # calibration for Schwarzschild equatorial b=0 ray
        MAP.refRayTime = GEOM.d - GEOM.R_eq
        MAP.refRayTime += GEOM.r_s * log((GEOM.d - GEOM.r_s) / (GEOM.R_eq - GEOM.r_s))
        MAP.refRayTime /= c
    else: # more general lag calibration
        MAP.refRayTime = MAP.ORIGIN[0]

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

            status = RK(RAYS[T], GEOM)

            MAP.LAG[INDEX] = -1.0 * RAYS[T].STATE[0] / c
            MAP.PHI[INDEX] = RAYS[T].STATE[1]
            MAP.RADIAL[INDEX] = RAYS[T].STATE[2]
            MAP.THETA[INDEX] = RAYS[T].STATE[4]

            if RAYS[T].STATE[2] > 0.0:
                ZABB(GEOM,
                     RAYS[T].STATE,
                     RAYS[T].IMPACT,
                     MAP.Z + INDEX,
                     MAP.ABB + INDEX)
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

    MAP.ORIGIN[0] -= MAP.refRayTime
    MAP.ORIGIN[0] *= GEOM.mode_omega # use mode frequency instead

    INDEX = 0
    for i in range(N):
        for j in range(N):
            MAP.LAG[INDEX] -= MAP.refRayTime
            MAP.LAG[INDEX] *= GEOM.mode_omega # use mode frequency instead
            INDEX += 1

    return SUCCESS
