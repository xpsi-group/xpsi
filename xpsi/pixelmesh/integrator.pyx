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
from xpsi.pixelmesh.surfaceBisection cimport SURFACE_SCALE

from xpsi.surface_radiation_field.preload cimport (_preloaded,
                                                   init_preload,
                                                   free_preload)

from xpsi.surface_radiation_field.hot_wrapper cimport (init_hot,
                                               eval_hot_I,
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
              int[::1] cache_energy_indices,
              int[::1] cache_phase_indices,
              int single_precision_cache,
              int obl_surfgrav_ind,
              reuse_ray_map = None,
              global_to_local_file = None,
              atmosphere = None,
              atmosphere_extension = 1):

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
        size_t i, j, k, p, search
        size_t T, THREAD_INDEX, ROOT, INDEX
        size_t N_T = numThreads
        size_t N_E = energies.shape[0]
        size_t N_CE = cache_energy_indices.shape[0]
        size_t N_CP = cache_phase_indices.shape[0]
        size_t cache_counter_E
        size_t N_P = phases.shape[0]
        size_t NGR
        size_t NGR_SQ
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
        float[:,:,::1] IMAGE
        double[:,:,::1] IMAGE_DUB
        double SEMI_MAJOR, SEMI_MINOR

    GEOM.r_s = r_s
    GEOM.R_eq = R_eq
    GEOM.omega = omega
    GEOM.mode_omega = mode_frequency * _2pi
    GEOM.epsilon = epsilon
    GEOM.zeta = zeta
    GEOM.a = a
    GEOM.asq = a * a
    GEOM.kappa = kappa
    ## Pole-on guard for the image-plane ray tracer. At inclination 0 or pi
    ## the central pixel starts on the spin axis, where sin(theta) = 0 and
    ## the 1/sin^2(theta) terms in the metric (METRIC_qK.pyx) diverge.
    ## Clamping to 1e-3 rad (~0.06 deg) from the axis avoids this and is
    ## visually indistinguishable from an exactly pole-on view.
    if inclination < 1.0e-3:
        printf("\nWarning: inclination clamped to 1e-3 rad for imaging (pole-on).")
        inclination = 1.0e-3
    elif inclination > M_PI - 1.0e-3:
        printf("\nWarning: inclination clamped to pi - 1e-3 rad for imaging (pole-on).")
        inclination = M_PI - 1.0e-3

    GEOM.inclination = inclination
    GEOM.obl_surfgrav_ind = obl_surfgrav_ind

    ## Surface scale for the ray tracer (added 2026-11). GEOM.R_eq stays the
    ## physical equatorial circumferential radius (it is also used for the
    ## effective gravity in surface_radiation_field). GEOM.R_eq_sph is the
    ## equatorial radius of the surface in the spherical coordinate r_sph
    ## built from Boyer-Lindquist r, chosen so that the circumferential
    ## equatorial radius equals R_eq. See the MODELLING CHOICE note at
    ## surfaceBisection.RADIUS. To restore the pre-2026-11 behaviour, use
    ## GEOM.R_eq_sph = GEOM.R_eq instead.
    GEOM.R_eq_sph = SURFACE_SCALE(&GEOM)

    if R_eq > 1.5 * r_s: # compare to Schwarzschild photon sphere
        GEOM.b_max = R_eq / sqrt(1.0 - r_s / R_eq)
    else: # limiting value
        GEOM.b_max = sqrt(3) * 1.5 * r_s

    ## Image plane at d = 10^4 b_max from the star (BOUNDARY_CONDITIONS.pyx).
    GEOM.d = 10000.0 * GEOM.b_max

    ## 'distance' is the observer's distance (m). Image-plane coordinates are
    ## impact parameters, so the solid angle of a pixel is (area)/distance^2
    ## (PJ12 Eq. 25); subtracting d changes the result by ~d/distance ~ 1e-11
    ## and is immaterial.
    distance -= GEOM.d

    cdef _RAY **RAYS = NULL
    cdef RAY_MAP *MAP = NULL

    ## Either trace a new ray map (globalRayMap.pyx) or reuse one returned by a
    ## previous call. Returned arrays have N*N + 1 entries: index 0 is the
    ## centre ray, index 1 + i*N + j is ring i, azimuth j. X, Y are in units
    ## of b_max; THETA, PHI, RADIAL are the surface point in the spherical
    ## coordinates of the surface; LAG is the phase lag (rad); Z = 1 + z; and
    ## cos_zenith is the cosine of the comoving-frame emission angle.
    if reuse_ray_map is None:
        NGR = numGlobalRays
        NGR_SQ = NGR * NGR + 1

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
                             0,
                             obl_surfgrav_ind)

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

        NGR = r_MESH.shape[0]
        NGR_SQ = NGR * NGR + 1

    cdef:
        size_t N = NGR
        size_t NSQ = NGR * NGR
        int ORIGIN_HIT
        ## HIT: per-thread flags, 1 if pixel INDEX (1..N*N) sees radiation at the
        ## current phase. Pixels run 1..N*N, so each thread needs N*N + 1 slots.
        ## (Before 2026-11 this was N_T * N*N, so the last pixel of the last thread
        ## wrote one element past the end; it was also never freed.)
        int *HIT = <int*> calloc(N_T * (NSQ + 1), sizeof(int))
        double area, r1, r2
        double *r_FULL = <double*> malloc((NGR + 1) * sizeof(double))

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
    ## Surface radiation field: which pixels radiate (HIT_or_MISS), their local
    ## variables such as temperature (eval_local_variables), and the specific
    ## intensity (eval_hot_I); shared with the cellmesh pulse-profile code.
    cdef storage *local_vars_buf = init_local_variables(N_T, filepath)
    cdef _preloaded *preloaded = NULL
    cdef void *data = NULL
    if atmosphere:
        preloaded = init_preload(atmosphere)
        data = init_hot(N_T, preloaded, atmosphere_extension)
    else:
        data = init_hot(N_T, NULL, atmosphere_extension)

    cdef double[:,::1] integrated_flux = np.zeros((N_P, N_E), dtype = np.double)

    cdef double E_prime, I_E

    ## Azimuthal width of each pixel
    cdef double Delta_t = _2pi / NGR

    ## Ring radii r_0 = 0, r_1, ..., r_N = 1 in units of the screen semi-axes.
    ## r_MESH (returned to Python, and passed back in via reuse_ray_map) holds
    ## only r_0 .. r_(N-1): the outer ring r_N = 1 by construction
    ## (globalRayMap.init_RAY_MAP). The area weights below need r_(i+2) for
    ## ring i, i.e. up to r_N, so the full set is rebuilt here.
    ## (Before 2026-11 the loop read r_MESH[N], one element past the end of
    ## the array; see the note at the area calculation.)
    for i in range(NGR):
        r_FULL[i] = r_MESH[i]
    r_FULL[NGR] = 1.0

    if cache_intensities == 1:
        if single_precision_cache == 0:
            ## BUG FIXED 2026-11: the double-precision cache was allocated with
            ## dtype np.single, so the assignment to the double memoryview raised
            ## 'Buffer dtype mismatch' and single_precision_cache = 0 always failed.
            ## (-1.0 * zeros is simply zeros: pixels with no radiation stay 0.)
            IMAGE_DUB = -1.0 * np.zeros((N_CP, N_CE, NGR_SQ), dtype = np.double)
        else:
            IMAGE = -1.0 * np.zeros((N_CP, N_CE, NGR_SQ), dtype = np.single)

    printf("\nCommencing imaging...")

    ## Loop over observer rotational phases, in parallel. For each phase:
    ##   1. Flag the pixels whose surface point radiates at the EMISSION phase
    ##      phases[k] - LAG (photons from pixels with a longer travel time left
    ##      the star at an earlier rotational phase).
    ##   2. If at least 4 pixels radiate, sum the observed specific intensity
    ##      times pixel area over the image, for each energy. Otherwise the flux
    ##      is set to 0 (and a warning printed): a small hot region can vanish
    ##      if the grid is too coarse to resolve it.
    for kk in prange(N_P,
                     nogil = True,
                     schedule = 'dynamic',
                     num_threads = N_T,
                     chunksize = 1):

        k = <size_t> kk
        T = threadid()
        THREAD_INDEX = T * (NSQ + 1)

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
            cache_counter_E = 0
            for p in range(N_E): # possible to reorder loops, energy innermost?
                for i in range(NGR):
                    ROOT = i * N

                    ## Area weight of ring i (pixel radius r_(i+1)): the annulus between the
                    ## midpoints to the neighbouring rings, times the azimuthal width Delta_t:
                    ##   area = 0.5 (r2^2 - r1^2) Delta_t   (in units of SEMI_MAJOR * SEMI_MINOR).
                    ## The outermost ring's outer edge is placed symmetrically.
                    ## BUG FIXED 2026-11: this used r_MESH[i+2] and r_MESH[i+1] directly, but
                    ## r_MESH has only N entries, so ring N-2 read r_MESH[N], one element past
                    ## the end of the array (uninitialised memory). That ring lies just inside
                    ## the limb (radius 0.995 with the defaults), so whenever the limb radiates
                    ## its area weight was garbage. For a uniformly hot, non-rotating star the
                    ## total flux came out at 0.25 of the analytic value; it now agrees with
                    ## the analytic value to 0.5%.
                    ##
                    ## LIMB ACCURACY: a pixel counts with its full area if its central ray
                    ## hits the star, so the area just outside the limb is over-counted by up
                    ## to one ring width. With the defaults (N = 100, p = 1/2) the limb falls
                    ## just outside ring N-1, the worst case, and a uniformly bright disc comes
                    ## out 0.5% too bright; the error scales with the ring spacing at the limb.
                    r1 = 0.5 * (r_FULL[i+1] + r_FULL[i])
                    if i < NGR - 1:
                        r2 = 0.5 * (r_FULL[i+2] + r_FULL[i+1])
                    else:
                        r2 = r_FULL[i + 1] + (r_FULL[i + 1] - r1)

                    area = 0.5 * (pow(r2, 2.0) - pow(r1, 2.0))

                    for j in range(NGR):
                        INDEX = ROOT + j + 1
                        if HIT[THREAD_INDEX + INDEX] == 1:
                            ## Observed energy E corresponds to emitted energy E' = E (1 + z) = E Z.
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

                            I_E = eval_hot_I(T,
                                           E_prime,
                                           cos_zenith[INDEX],
                                           local_vars_buf.local_variables[T],
                                           data,
                                           0)

                            ## Observed specific intensity: I_obs(E) = I_emit(E Z) / Z^3, from
                            ## invariance of I/E^3 (PJ12, Section 4).
                            I_E = I_E * pow(Z[INDEX], -3.0)

                            ## Flux = sum of I_obs times solid angle; normalised to physical units after
                            ## the loop over rings (semi-axes, atmosphere normalisation, distance).
                            integrated_flux[k,p] += I_E * area * Delta_t

                            ## Sky-map cache: I_obs / E (photon specific intensity, unnormalised) for
                            ## the requested phases and energies, stored at the pixel's own INDEX so
                            ## that it lines up with X[INDEX], Y[INDEX] (PhotospherePlotter plots them
                            ## together). BUG FIXED 2026-11: it was stored at INDEX + 1, which shifted
                            ## every sky map by one pixel in azimuth (3.6 deg about the screen centre
                            ## for N = 100), wrapping from the end of each ring onto the next.
                            if cache_intensities == 1:
                                for search in range(N_CP):
                                    if k == <size_t> cache_phase_indices[search]:
                                        if cache_counter_E < N_CE and p == <size_t> cache_energy_indices[cache_counter_E]:
                                            if single_precision_cache == 0:
                                                IMAGE_DUB[cache_phase_indices[search], cache_counter_E, INDEX] = I_E / energies[p]
                                            else:
                                                IMAGE[cache_phase_indices[search], cache_counter_E, INDEX] = I_E / energies[p]

                ## Centre ray: a disc of radius half the first ring radius.
                if ORIGIN_HIT == 1:
                    E_prime = energies[p] * Z[0]

                    eval_local_variables(THETA[0],
                                         PHI[0],
                                         phases[k] - LAG[0],
                                         &GEOM,
                                         &(radiation_field_global_variables[0]),
                                         local_vars_buf,
                                         T)

                    I_E = eval_hot_I(T,
                                   E_prime,
                                   cos_zenith[0],
                                   local_vars_buf.local_variables[T],
                                   data,
                                   0)

                    I_E = I_E * pow(Z[0], -3.0)

                    r2 = 0.5 * r_FULL[1]

                    integrated_flux[k,p] += I_E * pow(r2, 2.0) * M_PI

                    if cache_intensities == 1:
                        for search in range(N_CP):
                            if k == <size_t> cache_phase_indices[search]:
                                if cache_counter_E < N_CE and p == <size_t> cache_energy_indices[cache_counter_E]:
                                    if single_precision_cache == 0:
                                        IMAGE_DUB[cache_phase_indices[search], cache_counter_E, 0] = I_E / energies[p]
                                    else:
                                        IMAGE[cache_phase_indices[search], cache_counter_E, 0] = I_E / energies[p]

                if cache_intensities == 1:
                    for search in range(N_CP):
                        if k == <size_t> cache_phase_indices[search]:
                            if cache_counter_E < N_CE and p == <size_t> cache_energy_indices[cache_counter_E]:
                                cache_counter_E = cache_counter_E + 1

                ## Convert to physical units: image-plane area in m^2 (semi-axes in m),
                ## intensity normalisation of the atmosphere (eval_hot_norm), division by
                ## distance^2 (solid angle) and by E (photon flux).
                integrated_flux[k,p] *= SEMI_MINOR * SEMI_MAJOR
                integrated_flux[k,p] *= eval_hot_norm() / (distance * distance * energies[p] * keV)
        else:
            for p in range(N_E):
                integrated_flux[k,p] = 0.0

        # printf("\n ")
        # printf("\n------------------------------------------------------------->>>")
        # printf("\nPhase index: %d", <int>k)
        # printf("\nFraction of hot rays: %.6f", <double>(<double>NUM_HITS/<double>(NUM_MISSES + NUM_HITS)))
        if NUM_HITS < 4:
            printf("\nRadiating region poorly resolved... integrator not called.")
        # else:
        #     printf("\nIntegrator called.")
        # printf("\n------------------------------------------------------------->>>")
        # printf("\n ")

    #----------------------------------------------------------------------->>>
    # >>> Free memory.
    #----------------------------------------------------------------------->>>
    if atmosphere:
        free_preload(preloaded)

    free_hot(N_T, data)
    free_local_variables(N_T, local_vars_buf)
    free(HIT)
    free(r_FULL)

    if cache_intensities == 1:
        if single_precision_cache == 0:
            _IMAGE = np.asarray(IMAGE_DUB, dtype = np.double, order = 'C')
        else:
            _IMAGE = np.asarray(IMAGE, dtype = np.single, order = 'C')
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
                np.asarray(r_MESH, dtype = np.double, order= 'C'),
                SEMI_MAJOR,
                SEMI_MINOR,
                _IMAGE)
    else:
        return (SUCCESS,
                np.asarray(integrated_flux, dtype = np.double, order = 'C'),
                _IMAGE)
