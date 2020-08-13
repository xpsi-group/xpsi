#cython: cdivision=True
#cython: boundscheck=False
#cython: nonecheck=False
#cython: wraparound=False
#cython: embedsignature=True

""" Image a star on a distant observer's sky. """

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

from xpsi.surface_radiation_field.preload cimport (_preloaded,
                                                   init_preload,
                                                   free_preload)

from xpsi.surface_radiation_field.hot cimport (init_hot,
                                               eval_hot,
                                               eval_hot_norm,
                                               free_hot)

from xpsi.surface_radiation_field.local_variables cimport (storage,
                                                           HIT_or_MISS,
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
              reuse_ray_map = None,
              global_to_local_file = None,
              atmosphere = None):

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

        double[::1] X
        double[::1] Y
        double[::1] THETA
        double[::1] PHI
        double[::1] RADIAL
        double[::1] cos_zenith
        double[::1] LAG
        double[::1] Z
        double[::1] r_MESH
        double[:,:,::1] IMAGE
        double SEMI_MAJOR, SEMI_MINOR

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

    cdef _RAY **RAYS = NULL
    cdef RAY_MAP *MAP = NULL

    if reuse_ray_map is None:
        X = np.zeros(NGR_SQ, dtype = np.double)
        Y = np.zeros(NGR_SQ, dtype = np.double)
        THETA = np.zeros(NGR_SQ, dtype = np.double)
        PHI = np.zeros(NGR_SQ, dtype = np.double)
        RADIAL = np.zeros(NGR_SQ, dtype = np.double)
        cos_zenith = np.zeros(NGR_SQ, dtype = np.double)
        LAG = np.zeros(NGR_SQ, dtype = np.double)
        Z = np.zeros(NGR_SQ, dtype = np.double)
        r_MESH = np.zeros(NGR, dtype = np.double)

        RAYS = <_RAY**> malloc(N_T * sizeof(_RAY*))
        for T in range(N_T):
            RAYS[T] = alloc_RAY(epsabs_ray,
                                epsrel_ray,
                                init_step,
                                max_steps)
        #-------------------------------------------------------------------->>>
        # >>> Compute sparse set of rays distributed over global image plane.
        #-------------------------------------------------------------------->>>
        MAP = alloc_RAY_MAP(NGR)
        # The last argument if one can force the image plane to be circular,
        # allowing for a trivial mapping between mesh coordinates and a
        # plotting library, but we use elliptical image plane with matplotlib.
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
        RADIAL[0] = MAP.ORIGIN[2]
        Z[0] = MAP.ORIGIN[4]
        cos_zenith[0] = MAP.ORIGIN[5]

        SEMI_MAJOR = MAP.SEMI_MAJOR
        SEMI_MINOR = MAP.SEMI_MINOR

        INDEX = 1
        while INDEX < NGR_SQ:
            X[INDEX] = MAP.X_MESH[INDEX - 1] / GEOM.b_max
            Y[INDEX] = MAP.Y_MESH[INDEX - 1] / GEOM.b_max
            LAG[INDEX] = MAP.LAG[INDEX - 1]
            THETA[INDEX] = MAP.THETA[INDEX - 1]
            PHI[INDEX] = MAP.PHI[INDEX - 1]
            RADIAL[INDEX] = MAP.RADIAL[INDEX - 1]
            Z[INDEX] = MAP.Z[INDEX - 1]
            cos_zenith[INDEX] = MAP.ABB[INDEX - 1]
            if INDEX <= NGR:
                r_MESH[INDEX - 1] = MAP.r_MESH[INDEX - 1]
            INDEX += 1

        # free memory
        for T in range(N_T):
            free_RAY(RAYS[T])
        free(RAYS)
        free_RAY_MAP(MAP)

    else: # reuse ray map
        X = reuse_ray_map[0]
        Y = reuse_ray_map[1]
        THETA = reuse_ray_map[2]
        PHI = reuse_ray_map[3]
        RADIAL = reuse_ray_map[4]
        LAG = reuse_ray_map[5]
        Z = reuse_ray_map[6]
        cos_zenith = reuse_ray_map[7]
        r_MESH = reuse_ray_map[8]
        SEMI_MAJOR = reuse_ray_map[9]
        SEMI_MINOR = reuse_ray_map[10]

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
    cdef const unsigned char[:] _filepath
    cdef const char *filepath = NULL
    if global_to_local_file is not None:
        if isinstance(global_to_local_file, unicode):
            global_to_local_file = (<unicode>global_to_local_file).encode('utf8')
        _filepath = global_to_local_file
        filepath = <const char*>&(_filepath[0])

    # Initialise the source radiation field
    cdef storage *local_vars_buf = init_local_variables(N_T, filepath)
    cdef _preloaded *preloaded = NULL
    cdef void *data = NULL
    if atmosphere:
        preloaded = init_preload(atmosphere)
        data = init_hot(N_T, preloaded)
    else:
        data = init_hot(N_T, NULL)

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

        if RADIAL[0] > 0.0:
            if HIT_or_MISS(THETA[0],
                           PHI[0],
                           phases[k] - LAG[0],
                           &(radiation_field_global_variables[0]),
                           local_vars_buf) == 1:
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
                INDEX = ROOT + j + 1
                if RADIAL[INDEX] > 0.0:
                    if HIT_or_MISS(THETA[INDEX],
                                   PHI[INDEX],
                                   phases[k] - LAG[INDEX],
                                   &(radiation_field_global_variables[0]),
                                   local_vars_buf) == 1:

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

                    r1 = 0.5 * (r_MESH[i+1] + r_MESH[i])
                    if i < NGR - 1:
                        r2 = 0.5 * (r_MESH[i+2] + r_MESH[i+1])
                    else:
                        r2 = r_MESH[i + 1] + (r_MESH[i + 1] - r1)

                    area = 0.5 * (pow(r2, 2.0) - pow(r1, 2.0))

                    for j in range(NGR):
                        INDEX = ROOT + j + 1
                        if HIT[THREAD_INDEX + INDEX] == 1:
                            E_prime = energies[p] * Z[INDEX]
                            #printf("\n(i,j): %d, %d", <int>i, <int>j)
                            # The following function calls are split so that
                            # the module the latter calls is shared by the
                            # PixelMesh and CellMesh routines
                            # Move this first evaluation to outside the loop
                            # over energies since it is achromatic?
                            eval_local_variables(THETA[INDEX],
                                                 PHI[INDEX],
                                                 phases[k] - LAG[INDEX],
                                                 &GEOM,
                                                 &(radiation_field_global_variables[0]),
                                                 local_vars_buf,
                                                 T)

                            I_E = eval_hot(T,
                                           E_prime,
                                           ABB[INDEX],
                                           local_vars_buf.local_variables[T],
                                           data)

                            I_E = I_E * pow(Z[INDEX], -3.0)

                            integrated_flux[k,p] += I_E * area * Delta_t

                            if cache_intensities == 1:
                                IMAGE[k, p, INDEX + 1] = I_E / energies[p]

                if ORIGIN_HIT == 1:
                    E_prime = energies[p] * Z[0]

                    eval_local_variables(THETA[0],
                                         PHI[0],
                                         phases[k] - LAG[0],
                                         &GEOM,
                                         &(radiation_field_global_variables[0]),
                                         local_vars_buf,
                                         T)

                    I_E = eval_hot(T,
                                   E_prime,
                                   cos_zenith[0],
                                   local_vars_buf.local_variables[T],
                                   data)

                    I_E = I_E * pow(Z[0], -3.0)

                    r2 = 0.5 * r_MESH[1]

                    integrated_flux[k,p] += I_E * pow(r2, 2.0) * M_PI

                    if cache_intensities == 1:
                        IMAGE[k, p, 0] = I_E / energies[p]

                integrated_flux[k,p] *= SEMI_MINOR * SEMI_MAJOR
                integrated_flux[k,p] *= eval_hot_norm() / (distance * distance * energies[p] * keV)
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
    if atmosphere:
        free_preload(preloaded)

    free_hot(N_T, data)
    free_local_variables(N_T, local_vars_buf)

    if cache_intensities == 1:
        _IMAGE = np.asarray(IMAGE, dtype = np.double, order = 'C')
    else:
        _IMAGE = None

    if reuse_ray_map is None:
        return (SUCCESS,
                np.asarray(integrated_flux, dtype = np.double, order = 'C'),
                np.asarray(X, dtype = np.double, order = 'C'),
                np.asarray(Y, dtype = np.double, order = 'C'),
                np.asarray(THETA, dtype = np.double, order = 'C'),
                np.asarray(PHI, dtype = np.double, order = 'C'),
                np.asarray(RADIAL, dtype = np.double, order = 'C'),
                np.asarray(LAG, dtype = np.double, order = 'C'),
                np.asarray(Z, dtype = np.double, order = 'C'),
                np.asarray(cos_zenith, dtype = np.double, order = 'C'),
                np.asarray(r_MESH = np.double, order= 'C')
                SEMI_MAJOR,
                SEMI_MINOR,
                _IMAGE)
    else:
        return (SUCCESS,
                np.asarray(integrated_flux, dtype = np.double, order = 'C'),
                _IMAGE)
