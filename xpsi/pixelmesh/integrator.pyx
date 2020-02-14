#cython: cdivision=True
#cython: boundscheck=False
#cython: nonecheck=False
#cython: wraparound=False

from __future__ import division, print_function
import numpy as np
from cython.parallel cimport *
from libc.math cimport M_PI, log10, pow, sqrt
from libc.stdlib cimport calloc, malloc, free
from libc.stdio cimport printf, setbuf, stdout
from GSL cimport gsl_isnan

from xpsi.pixelmesh.globalRayMap cimport (RAY_MAP,
                                          alloc_RAY_MAP,
                                          compute_globalRayMap,
                                          free_RAY_MAP)

from xpsi.pixelmesh.coordinateTransformation cimport BOYERLINDQUIST_2_SPHERICAL

from xpsi.surface_radiation_field.hot_radiation_field cimport (init_hotRadField,
                                                       eval_hotRadField,
                                                       eval_hotRadField_norm,
                                                       free_hotRadField,
                                                       hotRadField_PRELOAD)

from xpsi.surface_radiation_field.local_variables cimport (HIT_or_MISS,
                                                           init_local_variables,
                                                           free_local_variables,
                                                           eval_local_variables)

from xpsi.pixelmesh.RK_IP2S_tracer cimport _RAY, alloc_RAY, free_RAY
from xpsi.pixelmesh.geometricConfiguration cimport _GEOM

cdef int SUCCESS = 0
cdef double keV = 1.60217662e-16
cdef double _2pi = 2.0 * M_PI

def integrate(size_t numThreads,
              double r_s,
              double R_eq,
              double omega,
              double mode_frequency,
              double zeta,
              double epsilon,
              double a,
              double kappa,
              double distance,
              double inclination,
              size_t numGlobalRays,
              double epsabs_ray,
              double epsrel_ray,
              size_t max_steps,
              double init_step,
              double radialIncrementExponent,
              double[::1] radiation_field_global_variables,
              double[::1] energies,
              double[::1] phases,
              int cache_intensities,
              atmosphere):

    """
    Compute image-plane-to-star X-ray light-curves when the photospheric
    radiation field is global in nature and or if adaptive focussing is not
    appropriate.

    For instances where a hot, localised region of plasma dominates the flux
    at higher X-ray energies, it might be advisable to use two integration
    routines, one which tracks the spot and a second global integration routine
    which neglects the spot, but because it exhibits a small angular extent,
    the error in the soft X-rays due to not resolving the boundary well is
    small relative to the total flux from the star.

    This routine is OpenMP-parallelised.

    """
    setbuf(stdout, NULL)
    printf("Imaging the star and computing light-curves...")

    cdef:
        signed int kk
        size_t i, j, k, p
        size_t T, THREAD_INDEX, ROOT, INDEX
        size_t N_T = numThreads
        size_t N_E = energies.shape[0]
        size_t N_P = phases.shape[0]
        size_t NGR = numGlobalRays
        size_t NGR_SQ = NGR * NGR + 1
        _GEOM GEOM

        double[::1] X = np.zeros(NGR_SQ, dtype = np.double)
        double[::1] Y = np.zeros(NGR_SQ, dtype = np.double)
        double[::1] THETA = np.zeros(NGR_SQ, dtype = np.double)
        double[::1] PHI = np.zeros(NGR_SQ, dtype = np.double)
        double[::1] cos_zenith = np.zeros(NGR_SQ, dtype = np.double)
        double[::1] LAG = np.zeros(NGR_SQ, dtype = np.double)
        double[::1] Z = np.zeros(NGR_SQ, dtype = np.double)
        double[:,:,::1] IMAGE

    if cache_intensities == 1:
        IMAGE = -1.0 * np.zeros((N_P, N_E, NGR_SQ), dtype = np.double)

    GEOM.r_s = r_s
    GEOM.R_eq = R_eq
    GEOM.omega = omega
    GEOM.mode_omega = mode_frequency * _2pi
    GEOM.epsilon = epsilon
    GEOM.zeta = zeta
    GEOM.a = a
    GEOM.asq = a * a
    GEOM.kappa = kappa
    GEOM.inclination = inclination

    if R_eq > 1.5 * r_s: # compare to Schwarzschild photon sphere
        GEOM.b_max = R_eq / sqrt(1.0 - r_s / R_eq)
    else: # limiting value
        GEOM.b_max = sqrt(3) * 1.5 * r_s

    GEOM.d = 10000.0 * GEOM.b_max

    distance -= GEOM.d

    cdef _RAY **RAYS = <_RAY**> malloc(N_T * sizeof(_RAY*))
    for T in range(N_T):
        RAYS[T] = alloc_RAY(epsabs_ray,
                            epsrel_ray,
                            init_step,
                            max_steps)
    #----------------------------------------------------------------------->>>
    # >>> Compute sparse set of rays distributed over global image plane.
    #----------------------------------------------------------------------->>>
    cdef RAY_MAP *MAP = alloc_RAY_MAP(NGR)
    # The last argument forces image plane to be circular, allowing for a
    # trivial mapping between mesh coordinates and a plotting library.
    compute_globalRayMap(N_T,
                         radialIncrementExponent,
                         &GEOM,
                         MAP,
                         RAYS,
                         0)

    printf("\n\nGlobal ray map computed.")

    BOYERLINDQUIST_2_SPHERICAL(MAP, &GEOM)

    printf("\nCoordinates transformed from Boyer-Lindquist to spherical.")

    # now store transformed spherical coordinates for return

    # Deal with origin and NGR + 1
    X[0] = MAP.ORIGIN_X / GEOM.b_max
    Y[0] = MAP.ORIGIN_Y / GEOM.b_max
    LAG[0] = MAP.ORIGIN[0]
    THETA[0] = MAP.ORIGIN[3]
    PHI[0] = MAP.ORIGIN[1]
    Z[0] = MAP.ORIGIN[4]
    cos_zenith[0] = MAP.ORIGIN[5]

    INDEX = 1
    while INDEX < NGR_SQ:
        X[INDEX] = MAP.X_MESH[INDEX - 1] / GEOM.b_max
        Y[INDEX] = MAP.Y_MESH[INDEX - 1] / GEOM.b_max
        LAG[INDEX] = MAP.LAG[INDEX - 1]
        THETA[INDEX] = MAP.THETA[INDEX - 1]
        PHI[INDEX] = MAP.PHI[INDEX - 1]
        Z[INDEX] = MAP.Z[INDEX - 1]
        cos_zenith[INDEX] = MAP.ABB[INDEX - 1]
        INDEX += 1

    cdef:
        size_t N = NGR
        size_t NSQ = NGR * NGR
        int ORIGIN_HIT
        int *HIT = <int*> calloc(N_T * NSQ, sizeof(int))
        double area, r1, r2

        unsigned long NUM_HITS, NUM_MISSES
    #----------------------------------------------------------------------->>>
    # >>> Integrate over approximate incident radiation field.
    #----------------------------------------------------------------------->>>

    # Initialise the source radiation field
    cdef void *local_variables = init_local_variables(N_T)
    cdef hotRadField_PRELOAD *src_preload = NULL
    cdef double[::1] cast
    cdef double[::1] intensity
    cdef void *data = NULL
    cdef void *ext_data = NULL
    cdef size_t num_args
    if atmosphere:
        args = atmosphere
        num_args = len(args)
        src_preload = <hotRadField_PRELOAD*> malloc(sizeof(hotRadField_PRELOAD))
        src_preload.params = <double**> malloc(sizeof(double*) * (num_args - 1))
        src_preload.S = <size_t*> malloc(sizeof(size_t) * (num_args - 2))
        for i in range(num_args - 1):
            cast = args[i]
            src_preload.params[i] = &cast[0]
            if i < num_args - 2:
                cast = args[i+1]
                src_preload.S[i] = cast.shape[0]
                if i < num_args - 3:
                    for j in range(i+2, num_args - 1):
                        cast = args[j]
                        src_preload.S[i] *= cast.shape[0]
        intensity = args[i+1]
        src_preload.I = &intensity[0]
        data = init_hotRadField(N_T, src_preload)
    else:
        data = init_hotRadField(N_T, NULL)

    cdef double[:,::1] integrated_flux = np.zeros((N_P, N_E), dtype = np.double)

    cdef double E_prime, I_E

    cdef double Delta_t = _2pi / NGR

    printf("\nCommencing imaging...")

    for kk in prange(N_P,
                     nogil = True,
                     schedule = 'dynamic',
                     num_threads = N_T,
                     chunksize = 1):

        k = <size_t> kk
        T = threadid()
        THREAD_INDEX = T * NSQ

        if MAP.ORIGIN[2] > 0.0:
            if HIT_or_MISS(MAP.ORIGIN[3],
                           MAP.ORIGIN[1],
                           phases[k] - MAP.ORIGIN[0],
                           &(radiation_field_global_variables[0])) == 1:
                NUM_HITS = 1
                NUM_MISSES = 0
                ORIGIN_HIT = 1
            else:
                NUM_HITS = 0
                NUM_MISSES = 1
                ORIGIN_HIT = 0
        else:
            NUM_HITS = 0
            NUM_MISSES = 0
            ORIGIN_HIT = 0

        for i in range(NGR):
            ROOT = i * N
            for j in range(NGR):
                INDEX = ROOT + j
                if MAP.RADIAL[INDEX] > 0.0:
                    if HIT_or_MISS(MAP.THETA[INDEX],
                                   MAP.PHI[INDEX],
                                   phases[k] - MAP.LAG[INDEX],
                                   &(radiation_field_global_variables[0])) == 1:

                        NUM_HITS = NUM_HITS + 1
                        HIT[THREAD_INDEX + INDEX] = 1
                    else:
                        NUM_MISSES = NUM_MISSES + 1
                        HIT[THREAD_INDEX + INDEX] = 0
                else:
                    NUM_MISSES = NUM_MISSES + 1
                    HIT[THREAD_INDEX + INDEX] = 0

        if NUM_HITS >= 4:
            for p in range(N_E): # possible to reorder loops, energy innermost?
                for i in range(NGR):
                    ROOT = i * N

                    r1 = 0.5 * (MAP.r_MESH[i+1] + MAP.r_MESH[i])
                    if i < NGR - 1:
                        r2 = 0.5 * (MAP.r_MESH[i+2] + MAP.r_MESH[i+1])
                    else:
                        r2 = MAP.r_MESH[i + 1] + (MAP.r_MESH[i + 1] - r1)

                    area = 0.5 * (pow(r2, 2.0) - pow(r1, 2.0))

                    for j in range(NGR):
                        INDEX = ROOT + j
                        if HIT[THREAD_INDEX + INDEX] == 1:
                            E_prime = energies[p] * MAP.Z[INDEX]
                            #printf("\n(i,j): %d, %d", <int>i, <int>j)
                            # The following function calls are split so that
                            # the module the latter calls is shared by the
                            # PixelMesh and CellMesh routines
                            # Move this first evaluation to outside the loop
                            # over energies since it is achromatic?
                            eval_local_variables(MAP.THETA[INDEX],
                                                 MAP.PHI[INDEX],
                                                 phases[k] - MAP.LAG[INDEX],
                                                 &GEOM,
                                                 &(radiation_field_global_variables[0]),
                                                 T,
                                                 local_variables)

                            I_E = eval_hotRadField(T,
                                                   E_prime,
                                                   MAP.ABB[INDEX],
                                                   (<double**>local_variables)[T],
                                                   data)

                            I_E = I_E * pow(MAP.Z[INDEX], -3.0)

                            integrated_flux[k,p] += I_E * area * Delta_t

                            if cache_intensities == 1:
                                IMAGE[k, p, INDEX + 1] = I_E / energies[p]

                if ORIGIN_HIT == 1:
                    E_prime = energies[p] * MAP.ORIGIN[4]

                    eval_local_variables(MAP.ORIGIN[3],
                                         MAP.ORIGIN[1],
                                         phases[k] - MAP.ORIGIN[0],
                                         &GEOM,
                                         &(radiation_field_global_variables[0]),
                                         T,
                                         local_variables)

                    I_E = eval_hotRadField(T,
                                           E_prime,
                                           MAP.ORIGIN[5],
                                           (<double**>local_variables)[T],
                                           data)

                    I_E = I_E * pow(MAP.ORIGIN[4], -3.0)

                    r2 = 0.5 * MAP.r_MESH[1]

                    integrated_flux[k,p] += I_E * pow(r2, 2.0) * M_PI

                    if cache_intensities == 1:
                        IMAGE[k, p, 0] = I_E / energies[p]

                integrated_flux[k,p] *= MAP.SEMI_MINOR * MAP.SEMI_MAJOR
                integrated_flux[k,p] *= eval_hotRadField_norm() / (distance * distance * energies[p] * keV)
        else:
            for p in range(N_E):
                integrated_flux[k,p] = 0.0

        printf("\n ")
        printf("\n------------------------------------------------------------->>>")
        printf("\nPhase index: %d", <int>k)
        printf("\nFraction of hot rays: %.6f", <double>(<double>NUM_HITS/<double>(NUM_MISSES + NUM_HITS)))
        if NUM_HITS < 4:
            printf("\nRadiating region poorly resolved... integrator not called.")
        else:
            printf("\nIntegrator called.")
        printf("\n------------------------------------------------------------->>>")
        printf("\n ")

    #----------------------------------------------------------------------->>>
    # >>> Free memory.
    #----------------------------------------------------------------------->>>
    for T in range(N_T):
        free_RAY(RAYS[T])
    free(RAYS)
    free_RAY_MAP(MAP)

    free_hotRadField(N_T, data)
    free_local_variables(N_T, local_variables)

    if cache_intensities == 1:
        _IMAGE = np.asarray(IMAGE, dtype = np.double, order = 'C')
    else:
        _IMAGE = None

    return (SUCCESS,
            np.asarray(integrated_flux, dtype = np.double, order = 'C'),
            np.asarray(X, dtype = np.double, order = 'C'),
            np.asarray(Y, dtype = np.double, order = 'C'),
            np.asarray(THETA, dtype = np.double, order = 'C'),
            np.asarray(PHI, dtype = np.double, order = 'C'),
            np.asarray(LAG, dtype = np.double, order = 'C'),
            np.asarray(Z, dtype = np.double, order = 'C'),
            np.asarray(cos_zenith, dtype = np.double, order = 'C'),
            _IMAGE)

