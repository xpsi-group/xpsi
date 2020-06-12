#cython: cdivision=True
#cython: boundscheck=False
#cython: nonecheck=False
#cython: wraparound=False
#cython: embedsignature=True

""" Integrate over the image of a star on a distant observer's sky. """

from __future__ import division
import numpy as np
cimport numpy as np
from cython.parallel cimport *
from libc.math cimport M_PI, sqrt, sin, cos, acos, log10, pow, exp
from libc.stdlib cimport malloc, free
from libc.stdio cimport printf

cdef double _pi = M_PI
cdef double _2pi = 2.0 * M_PI
cdef double c = 2.99792458e8
cdef double keV = 1.60217662e-16
cdef double erg = 1.0e-7
cdef double G = 6.6730831e-11
cdef double _h_keV = 4.135667662e-18
cdef double k_B = 1.38064852e-23
cdef double _h = 6.62607004e-34
cdef double SB = 5.6704e-5 # cgs
cdef double Planck_dist_const = 5.040366110812353e22
cdef double k_B_over_keV = k_B / keV

cdef int SUCCESS = 0
cdef int ERROR = 1

cdef int VERBOSE = 1
cdef int QUIET = 0

from xpsi.surface_radiation_field.preload cimport (_preloaded,
                                                   init_preload,
                                                   free_preload)

from xpsi.surface_radiation_field.hot cimport (init_hot,
                                               eval_hot,
                                               eval_hot_norm,
                                               free_hot)

from xpsi.surface_radiation_field.elsewhere cimport (init_elsewhere,
                                                     free_elsewhere,
                                                     eval_elsewhere,
                                                     eval_elsewhere_norm)

cdef void INVIS(size_t k,
                size_t N_L,
                int *const InvisFlag,
                double *const InvisStep,
                size_t *const InvisPhase,
                double *const _PHASE,
                double *const _GEOM,
                double *const _Z,
                double *const _ABB) nogil:

    cdef size_t m

    _GEOM[k] = 0.0
    _ABB[k] = 0.0

    if InvisPhase[0] == 0:
        _Z[k] = 0.0
    else:
        _Z[k] = _Z[InvisPhase[0] - 1]

    if k == N_L - 1 and InvisFlag[0] == 0:
        _PHASE[k] = _PHASE[k - 1] + InvisStep[0]
    elif InvisFlag[0] == 0:
        InvisFlag[0] = 1
        InvisPhase[0] = k
    elif k == N_L - 1 and InvisFlag[0] == 1:
        if InvisStep[0] == 0.0:
            InvisFlag[0] = 2
        else:
            for m in range(InvisPhase[0], k + 1):
                _PHASE[m] = _PHASE[m - 1] + InvisStep[0]

def integrate(size_t numThreads,
              double R,
              double omega,
              double r_s,
              double inclination,
              double[:,::1] cellArea,
              double[::1] radialCoords_of_parallels,
              double[::1] r_s_over_r,
              double[:,::1] theta,
              double[:,::1] phi,
              double[:,:,::1] srcCellParams,
              int[:,::1] CELL_RADIATES,
              correction_srcCellParams,
              int numRays,
              double[:,::1] deflection,
              double[:,::1] cos_alphaMatrix,
              double[:,::1] lag,
              double[::1] maxDeflection,
              double[::1] cos_gammaArray,
              double[::1] energies,
              double[::1] leaves,
              double[::1] phases,
              hot_atmosphere,
              elsewhere_atmosphere):

    #----------------------------------------------------------------------->>>
    # >>> General memory allocation.
    # >>>
    #----------------------------------------------------------------------->>>
    cdef:
        signed int ii
        size_t i, j, k, m, n, p # Array indexing
        size_t T, twoT # Used globally to represent thread index
        size_t N_T = numThreads # shared
        size_t N_R = numRays # shared
        size_t N_E = energies.shape[0] # shared
        size_t N_L = leaves.shape[0] # shared
        size_t N_P = phases.shape[0] # shared
        double sin_i = sin(inclination) # shared
        double cos_i = cos(inclination) # shared
        double Grav_z # Gravitational redshift; thread private (hereafter TP)
        double radius # Radial coordinate of pixel; TP
        double psi # Colatitude relative to star-observer direction; TP
        double cos_psi, sin_psi
        double deriv # $\frac{d\cos\alpha}{d\cos\phi}$; TP
        double beta # Surface velocity in the local NRF; TP
        double cos_alpha # Emission angle w.r.t outward radial direction in CRF; TP
        double sin_alpha
        double cos_delta # Unit sphere trigonometric identity; TP
        double cos_gamma # Surface normal tilt w.r.t outward radial direction; TP
        double cos_xi # Emission angle relative to the NRF surface velocity; TP
        double eta # Doppler boost factor in the NRF; TP
        double mu # NRF and then CRF emission angle w.r.t surface normal; TP
        double E_prime # Photon energy in the CRF, given an energy at infinity; TP
        double I_E # Radiant and or spectral intensity in the CRF; TP
        double __PHASE, __PHASE_plusShift, __GEOM, __Z, __ABB # TP
        double phi_shift # TP
        double superlum # TP
        double cos_gamma_sq, sin_gamma_sq
        double cos_theta_ij, sin_theta_ij
        double sqrt_cos_gamma_sq
        double theta_ij_over_pi
        double beta_sq
        double Lorentz
        double correction_I_E

        double[:,:,::1] privateFlux = np.zeros((N_T, N_P, N_E), dtype = np.double)
        double[:,::1] flux = np.zeros((N_E, N_P), dtype = np.double)

        int *terminate = <int*> malloc(N_T * sizeof(int))

        int *InvisFlag = <int*> malloc(N_T * sizeof(int))
        double *InvisStep = <double*> malloc(N_T * sizeof(double))
        size_t *InvisPhase = <size_t*> malloc(N_T * 2 * sizeof(size_t))

        accel **accel_alpha = <accel**> malloc(N_T * sizeof(accel*))
        accel **accel_lag = <accel**> malloc(N_T * sizeof(accel*))
        interp **interp_alpha = <interp**> malloc(N_T * sizeof(interp*))
        interp **interp_lag = <interp**> malloc(N_T * sizeof(interp*))

        interp *interp_alpha_store
        interp *interp_lag_store

        # Geometric quantity memory allocation
        double **_GEOM = <double**> malloc(N_T * sizeof(double*))
        double **_PHASE = <double**> malloc(N_T * sizeof(double*))
        double **_Z = <double**> malloc(N_T * sizeof(double*))
        double **_ABB = <double**> malloc(N_T * sizeof(double*))

        # Splines for geometric quantites
        accel **accel_GEOM = <accel**> malloc(N_T * sizeof(accel*))
        interp **interp_GEOM = <interp**> malloc(N_T * sizeof(interp*))
        accel **accel_Z = <accel**> malloc(N_T * sizeof(accel*))
        interp **interp_Z = <interp**> malloc(N_T * sizeof(interp*))
        accel **accel_ABB = <accel**> malloc(N_T * sizeof(accel*))
        interp **interp_ABB = <interp**> malloc(N_T * sizeof(interp*))

        double *defl_ptr
        double *alpha_ptr
        double *lag_ptr
        double *GEOM_ptr
        double *Z_ptr
        double *ABB_ptr
        double *phase_ptr

    for T in range(N_T):
        terminate[T] = 0
        accel_alpha[T] = gsl_interp_accel_alloc()
        interp_alpha[T] = gsl_interp_alloc(gsl_interp_steffen, N_R)
        accel_lag[T] = gsl_interp_accel_alloc()
        interp_lag[T] = gsl_interp_alloc(gsl_interp_steffen, N_R)

        _GEOM[T] = <double*> malloc(N_L * sizeof(double))
        _PHASE[T] = <double*> malloc(N_L * sizeof(double))
        _Z[T] = <double*> malloc(N_L * sizeof(double))
        _ABB[T] = <double*> malloc(N_L * sizeof(double))

        accel_Z[T] = gsl_interp_accel_alloc()
        interp_Z[T] = gsl_interp_alloc(gsl_interp_steffen, N_L)
        accel_ABB[T] = gsl_interp_accel_alloc()
        interp_ABB[T] = gsl_interp_alloc(gsl_interp_steffen, N_L)
        accel_GEOM[T] = gsl_interp_accel_alloc()
        interp_GEOM[T] = gsl_interp_alloc(gsl_interp_steffen, N_L)

    cdef double[:,::1] cos_deflection = np.zeros((deflection.shape[0],
                                                  deflection.shape[1]),
                                                 dtype = np.double)

    for i in range(deflection.shape[0]):
        for j in range(deflection.shape[1]):
            cos_deflection[i,j] = cos(deflection[i,j])

    # initialise the source radiation field
    cdef _preloaded *hot_preloaded = NULL
    cdef _preloaded *ext_preloaded = NULL
    cdef void *hot_data = NULL
    cdef void *ext_data = NULL

    if hot_atmosphere:
        hot_preloaded = init_preload(hot_atmosphere)
        hot_data = init_hot(N_T, hot_preloaded)
    else:
        hot_data = init_hot(N_T, NULL)

    cdef double[:,:,::1] correction
    cdef int perform_correction

    if correction_srcCellParams is not None:
        perform_correction = 1
        correction = correction_srcCellParams
    else:
        perform_correction = 0

    if perform_correction == 1:
        if elsewhere_atmosphere:
            ext_preloaded = init_preload(elsewhere_atmosphere)
            ext_data = init_elsewhere(N_T, ext_preloaded)
        else:
            ext_data = init_elsewhere(N_T, NULL)

    #----------------------------------------------------------------------->>>
    # >>> Integrate.
    # >>>
    #----------------------------------------------------------------------->>>
    for ii in prange(<signed int>cellArea.shape[0],
                     nogil = True,
                     schedule = 'static',
                     num_threads = N_T,
                     chunksize = 1):

        # Thread index
        T = threadid(); twoT = 2*T
        i = <size_t> ii

        radius = radialCoords_of_parallels[i]
        Grav_z = sqrt(1.0 - r_s_over_r[i])
        cos_gamma = cos_gammaArray[i]
        cos_gamma_sq = cos_gamma * cos_gamma
        sin_gamma_sq = sqrt(1.0 - cos_gamma_sq)

        gsl_interp_accel_reset(accel_alpha[T])
        gsl_interp_accel_reset(accel_lag[T])

        j = 0
        while deflection[i,j] > _pi:
            j = j + 1

        if j != 0:
            interp_alpha_store = interp_alpha[T]
            interp_alpha[T] = gsl_interp_alloc(gsl_interp_steffen, N_R - j)
            interp_lag_store = interp_lag[T]
            interp_lag[T] = gsl_interp_alloc(gsl_interp_steffen, N_R - j)

            defl_ptr = &(cos_deflection[i,j])
            alpha_ptr = &(cos_alphaMatrix[i,j])
            gsl_interp_init(interp_alpha[T], defl_ptr, alpha_ptr,  N_R - j)
            lag_ptr = &(lag[i,j])
            gsl_interp_init(interp_lag[T], defl_ptr, lag_ptr, N_R - j)
        else:
            interp_alpha_store = NULL
            interp_lag_store = NULL

            defl_ptr = &(cos_deflection[i,0])
            alpha_ptr = &(cos_alphaMatrix[i,0])
            gsl_interp_init(interp_alpha[T], defl_ptr, alpha_ptr, N_R)
            lag_ptr = &(lag[i,0])
            gsl_interp_init(interp_lag[T], defl_ptr, lag_ptr, N_R)

        j = 0
        # use this to decide whether or not to compute parallel:
        # Does the local vicinity of the parallel contain radiating material?
        while j < cellArea.shape[1]:
            if CELL_RADIATES[i,j] == 1:
                break
            j = j + 1

        #if j == cellArea.shape[1]:
        #    printf('%i', <int>i)

        InvisFlag[T] = 2
        InvisPhase[twoT] = 0
        InvisPhase[twoT + 1] = 0
        InvisStep[T] = 0.0

        cos_theta_ij = cos(theta[i,j])
        sin_theta_ij = sin(theta[i,j])
        theta_ij_over_pi = theta[i,j] / _pi
        beta = radius * omega * sin_theta_ij / (c * Grav_z)
        beta_sq = beta * beta
        Lorentz = sqrt(1.0 - beta_sq)

        correction_I_E = 0.0

        for k in range(N_L):
            cos_psi = cos_i * cos_theta_ij + sin_i * sin_theta_ij * cos(leaves[k])
            psi = acos(cos_psi)
            sin_psi = sin(psi)

            if psi <= maxDeflection[i]:
                if (cos_psi < interp_alpha[T].xmin or cos_psi > interp_alpha[T].xmax):
                    printf("cos_psi: %.16e\n", cos_psi)
                    printf("min: %.16e\n", interp_alpha[T].xmin)
                    printf("max: %.16e\n", interp_alpha[T].xmax)
                    terminate[T] = 1
                    break
                else:
                    cos_alpha = gsl_interp_eval(interp_alpha[T], defl_ptr, alpha_ptr, cos_psi, accel_alpha[T])

                sin_alpha = sqrt(1.0 - cos_alpha * cos_alpha)
                mu = cos_alpha * cos_gamma

                if sin_psi != 0.0:
                    cos_delta = (cos_i - cos_theta_ij * cos_psi) / (sin_theta_ij * sin_psi)
                    if theta_ij_over_pi < 0.5:
                        mu = mu + sin_alpha * sin_gamma_sq * cos_delta
                    else:
                        mu = mu - sin_alpha * sin_gamma_sq * cos_delta

                if mu > 0.0:
                    if sin_psi != 0.0:
                        cos_xi = sin_alpha * sin_i * sin(leaves[k]) / sin_psi
                        superlum = (1.0 + beta * cos_xi)
                        eta = Lorentz / superlum
                    else:
                        cos_xi = 0.0
                        superlum = 1.0
                        eta = Lorentz

                    # Beloborodov (2002): deriv = 1.0 - r_s_over_r[i]
                    deriv = gsl_interp_eval_deriv(interp_alpha[T], defl_ptr, alpha_ptr, cos_psi, accel_alpha[T])

                    _Z[T][k] = eta * Grav_z
                    _ABB[T][k] = mu * eta
                    _GEOM[T][k] = mu * deriv * Grav_z * eta * eta * eta / superlum

                    if (cos_psi < interp_lag[T].xmin or cos_psi > interp_lag[T].xmax):
                        printf("lag: %.16e\n", cos_psi)
                        printf("min: %.16e\n", interp_lag[T].xmin)
                        printf("max: %.16e\n", interp_lag[T].xmax)
                        terminate[T] = 1
                        break
                    else:
                        _PHASE[T][k] = leaves[k] + gsl_interp_eval(interp_lag[T], defl_ptr, lag_ptr, cos_psi, accel_lag[T])
                        #printf("\nphase, cos_psi, i, k, leaf: %.8e, %.8e, %i, %i, %.8e", _PHASE[T][k], cos_psi, <int>i, <int>k, leaves[k])

                    # Check whether cell was visible at previous rotation step
                    if InvisFlag[T] == 1 or (k != 0 and InvisFlag[T] == 2):
                        InvisPhase[twoT + 1] = k
                        if InvisPhase[twoT] == 0:
                            _PHASE[T][0] = leaves[0]
                            InvisStep[T] = (_PHASE[T][k] - _PHASE[T][0]) / (<double>k)
                            for m in range(1, k):
                                _PHASE[T][m] = _PHASE[T][m - 1] + InvisStep[T]
                        else:
                            InvisStep[T] = ((_PHASE[T][k] - _PHASE[T][InvisPhase[twoT] - 1])
                                                    / (<double>k - <double>InvisPhase[twoT] + 1.0))
                            for m in range(InvisPhase[twoT], k):
                                _PHASE[T][m] = _PHASE[T][m - 1] + InvisStep[T]

                    # Reset visibility flag
                    InvisFlag[T] = 0

                else:
                    #printf("\ncos_psi, i, k, leaf: %.8e, %i, %i, %.8e", cos_psi, <int>i, <int>k, leaves[k])
                    INVIS(k,
                          N_L,
                          InvisFlag + T,
                          InvisStep + T,
                          InvisPhase + twoT,
                          _PHASE[T],
                          _GEOM[T],
                          _Z[T],
                          _ABB[T])
            else:
                #printf("\ncos_psi, i, k, leaf: %.8e, %i, %i, %.8e", cos_psi, <int>i, <int>k, leaves[k])
                INVIS(k,
                      N_L,
                      InvisFlag + T,
                      InvisStep + T,
                      InvisPhase + twoT,
                      _PHASE[T],
                      _GEOM[T],
                      _Z[T],
                      _ABB[T])

        if interp_alpha_store != NULL:
            gsl_interp_free(interp_alpha[T])
            gsl_interp_free(interp_lag[T])
            interp_alpha[T] = interp_alpha_store
            interp_lag[T] = interp_lag_store

        if terminate[T] == 1:
            break
        elif terminate[T] == 0 and InvisFlag[T] != 2:
            for n in range(1, N_L):
                if _PHASE[T][n] <= _PHASE[T][n-1]:
                    terminate[T] = 1
                    break
            if terminate[T] == 0:
                # Initialise geometric interps
                phase_ptr = _PHASE[T]
                GEOM_ptr = _GEOM[T]
                Z_ptr = _Z[T]
                ABB_ptr = _ABB[T]

                gsl_interp_accel_reset(accel_Z[T])
                gsl_interp_init(interp_Z[T], phase_ptr, Z_ptr, N_L)
                gsl_interp_accel_reset(accel_ABB[T])
                gsl_interp_init(interp_ABB[T], phase_ptr, ABB_ptr, N_L)
                gsl_interp_accel_reset(accel_GEOM[T])
                gsl_interp_init(interp_GEOM[T], phase_ptr, GEOM_ptr, N_L)

                while j < cellArea.shape[1] and terminate[T] == 0:
                    if CELL_RADIATES[i,j] == 1:
                        phi_shift = phi[i,j]
                        for k in range(N_P):
                            __PHASE = phases[k]
                            __PHASE_plusShift = __PHASE + phi_shift
                            if __PHASE_plusShift > _PHASE[T][N_L - 1]:
                                while __PHASE_plusShift > _PHASE[T][N_L - 1]:
                                    __PHASE_plusShift = __PHASE_plusShift - _2pi
                            elif __PHASE_plusShift < _PHASE[T][0]:
                                while __PHASE_plusShift < _PHASE[T][0]:
                                    __PHASE_plusShift = __PHASE_plusShift + _2pi

                            __GEOM = gsl_interp_eval(interp_GEOM[T], phase_ptr, GEOM_ptr, __PHASE_plusShift, accel_GEOM[T])

                            if __GEOM > 0.0:
                                __Z = gsl_interp_eval(interp_Z[T], phase_ptr, Z_ptr, __PHASE_plusShift, accel_Z[T])
                                __ABB = gsl_interp_eval(interp_ABB[T], phase_ptr, ABB_ptr, __PHASE_plusShift, accel_ABB[T])

                                for p in range(N_E):
                                    E_prime = energies[p] / __Z
                                    I_E = eval_hot(T,
                                                   E_prime,
                                                   __ABB,
                                                   &(srcCellParams[i,j,0]),
                                                   hot_data)

                                    I_E = I_E * eval_hot_norm()

                                    if perform_correction == 1:
                                        correction_I_E = eval_elsewhere(T,
                                                               E_prime,
                                                               __ABB,
                                                               &(correction[i,j,0]),
                                                               ext_data)

                                        correction_I_E = correction_I_E * eval_elsewhere_norm()

                                    privateFlux[T,k,p] += cellArea[i,j] * (I_E - correction_I_E) * __GEOM
                            if terminate[T] == 1:
                                break

                    j = j + 1

    for i in range(N_E):
        for T in range(N_T):
            for k in range(N_P):
                flux[i,k] += privateFlux[T,k,i]

    for p in range(N_E):
        for k in range(N_P):
            flux[p,k] /= (energies[p] * keV)

    for T in range(N_T):
        gsl_interp_free(interp_alpha[T])
        gsl_interp_accel_free(accel_alpha[T])
        gsl_interp_free(interp_lag[T])
        gsl_interp_accel_free(accel_lag[T])

        free(_GEOM[T])
        free(_PHASE[T])
        free(_Z[T])
        free(_ABB[T])

        gsl_interp_free(interp_GEOM[T])
        gsl_interp_accel_free(accel_GEOM[T])
        gsl_interp_free(interp_Z[T])
        gsl_interp_accel_free(accel_Z[T])
        gsl_interp_free(interp_ABB[T])
        gsl_interp_accel_free(accel_ABB[T])

    free(interp_alpha)
    free(accel_alpha)
    free(interp_lag)
    free(accel_lag)

    free(_GEOM)
    free(_PHASE)
    free(_Z)
    free(_ABB)

    free(interp_GEOM)
    free(accel_GEOM)
    free(interp_Z)
    free(accel_Z)
    free(interp_ABB)
    free(accel_ABB)

    free(InvisFlag)
    free(InvisStep)
    free(InvisPhase)

    if hot_atmosphere:
        free_preload(hot_preloaded)

    free_hot(N_T, hot_data)

    if perform_correction == 1:
        if elsewhere_atmosphere:
            free_preload(ext_preloaded)

        free_elsewhere(N_T, ext_data)

    for T in range(N_T):
        if terminate[T] == 1:
            free(terminate)
            return (ERROR, None)

    #printf("\nIntegration completed.")
    return (SUCCESS, np.asarray(flux, dtype = np.double, order = 'C'))
