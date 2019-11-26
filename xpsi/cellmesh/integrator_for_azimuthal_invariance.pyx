#cython: cdivision=True
#cython: boundscheck=False
#cython: nonecheck=False
#cython: wraparound=False

from __future__ import division, print_function
import numpy as np
cimport numpy as np
from cython.parallel cimport *
from libc.math cimport M_PI, sqrt, sin, cos, acos, log10, pow, exp
from libc.stdlib cimport malloc, free
from libc.stdio cimport printf, setbuf, stdout

from xpsi.cellmesh.integrator cimport (gsl_interp_eval,
                                       gsl_interp_eval_deriv,
                                       gsl_interp_alloc,
                                       gsl_interp_accel,
                                       gsl_interp_accel_alloc,
                                       gsl_interp,
                                       gsl_interp_steffen,
                                       gsl_interp_init,
                                       gsl_interp_free,
                                       gsl_interp_accel_free,
                                       gsl_interp_accel_reset,
                                       gsl_isnan,
                                       gsl_isinf)

ctypedef gsl_interp_accel accel
ctypedef gsl_interp interp

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

cdef int SUCCESS = 0
cdef int ERROR = 1

cdef int VERBOSE = 1
cdef int QUIET = 0

from xpsi.surface_radiation_field.hot_radiation_field cimport (init_hotRadField,
                                                       free_hotRadField,
                                                       eval_hotRadField,
                                                       eval_hotRadField_norm,
                                                       hotRadField_PRELOAD)

from xpsi.surface_radiation_field.elsewhere_radiation_field cimport (init_elsewhereRadField,
                                                    free_elsewhereRadField,
                                                    eval_elsewhereRadField,
                                                    eval_elsewhereRadField_norm,
                                                    elsewhereRadField_PRELOAD)

cdef void INVIS(size_t k,
                size_t N_L,
                size_t N_E,
                const size_t *const BLOCK,
                int *const InvisFlag,
                double *const InvisStep,
                size_t *const InvisPhase,
                double *const PHASE,
                double *const PROFILE) nogil:
    # Revisit this... since one starts aligned to observer, some conditional
    # statements may not be required...
    cdef size_t p, m

    for p in range(N_E):
        (PROFILE + BLOCK[p] + k)[0] = 0.0

    if k == N_L - 1 and InvisFlag[0] == 0:
        PHASE[k] = PHASE[k - 1] + InvisStep[0]
    elif InvisFlag[0] == 0:
        InvisFlag[0] = 1
        InvisPhase[0] = k
    elif k == N_L - 1 and InvisFlag[0] == 1:
        if InvisStep[0] == 0.0:
            InvisFlag[0] = 2
        else:
            for m in range(InvisPhase[0], k + 1):
                PHASE[m] = PHASE[m - 1] + InvisStep[0]

#----------------------------------------------------------------------->>>
# >>> Integrate over the celestial sphere of distant observer.
# >>>
#----------------------------------------------------------------------->>>
def integrate_radField(size_t numThreads,
                         double M,
                         double R,
                         double omega,
                         double r_s,
                         double zeta,
                         double epsilon,
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
                         spot_atmosphere,
                         elsewhere_atmosphere):

    #----------------------------------------------------------------------->>>
    # >>> General memory allocation.
    # >>>
    #----------------------------------------------------------------------->>>
    cdef:
        signed int ii
        size_t i, j, J, k, m, n, p # Array indexing
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
        double _PHASE, _PHASE_plusShift, _GEOM, _Z, _ABB # TP
        double phi_shift # TP
        double superlum # TP
        double cos_gamma_sq, sin_gamma_sq
        double cos_theta_i, sin_theta_i
        double sqrt_cos_gamma_sq
        double theta_i_over_pi
        double beta_sq
        double Lorentz
        double correction_I_E

        double[:,:,::1] privateFlux = np.zeros((N_T, N_E, N_P), dtype = np.double)
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

        # Intensity spline interpolation
        double **PHASE = <double**> malloc(N_T * sizeof(double*))
        double **PROFILE = <double**> malloc(N_T * sizeof(double*))

        accel **accel_PROFILE = <accel**> malloc(N_T * sizeof(accel*))
        interp **interp_PROFILE = <interp**> malloc(N_T * sizeof(interp*))

        size_t *BLOCK = <size_t*> malloc(N_E * sizeof(size_t))

        double *defl_ptr
        double *alpha_ptr
        double *lag_ptr
        double *profile_ptr
        double *phase_ptr

    for T in range(N_T):
        terminate[T] = 0
        accel_alpha[T] = gsl_interp_accel_alloc()
        interp_alpha[T] = gsl_interp_alloc(gsl_interp_steffen, N_R)
        accel_lag[T] = gsl_interp_accel_alloc()
        interp_lag[T] = gsl_interp_alloc(gsl_interp_steffen, N_R)

        PHASE[T] = <double*> malloc(N_L * sizeof(double))
        PROFILE[T] = <double*> malloc(N_E * N_L * sizeof(double))
        accel_PROFILE[T] = gsl_interp_accel_alloc()
        interp_PROFILE[T] = gsl_interp_alloc(gsl_interp_steffen, N_L)

    for p in range(N_E):
        BLOCK[p] = p * N_L

    cdef double[:,::1] cos_deflection = np.zeros((deflection.shape[0],
                                                  deflection.shape[1]),
                                                 dtype = np.double)

    for i in range(deflection.shape[0]):
        for j in range(deflection.shape[1]):
            cos_deflection[i,j] = cos(deflection[i,j])

    # Initialise the source radiation field
    cdef hotRadField_PRELOAD *src_preload = NULL
    cdef elsewhereRadField_PRELOAD *ext_preload = NULL
    cdef double[::1] cast
    cdef double[::1] intensity
    cdef void *data = NULL
    cdef void *ext_data = NULL
    cdef size_t num_args
    if spot_atmosphere:
        args = spot_atmosphere
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

    cdef double[:,:,::1] correction
    cdef int perform_correction

    if correction_srcCellParams is not None:
        perform_correction = 1
        correction = correction_srcCellParams
    else:
        perform_correction = 0

    if perform_correction == 1:
        if elsewhere_atmosphere:
            args = elsewhere_atmosphere
            num_args = len(args)
            ext_preload = <elsewhereRadField_PRELOAD*> malloc(sizeof(elsewhereRadField_PRELOAD))
            ext_preload.params = <double**> malloc(sizeof(double*) * (num_args - 1))
            ext_preload.S = <size_t*> malloc(sizeof(size_t) * (num_args - 2))
            for i in range(num_args - 1):
                cast = args[i]
                ext_preload.params[i] = &cast[0]
                if i < num_args - 2:
                    cast = args[i+1]
                    ext_preload.S[i] = cast.shape[0]
                    if i < num_args - 3:
                        for j in range(i+2, num_args - 1):
                            cast = args[j]
                            ext_preload.S[i] *= cast.shape[0]
            intensity = args[i+1]
            ext_preload.I = &intensity[0]
            ext_data = init_elsewhereRadField(N_T, ext_preload)
        else:
            ext_data = init_elsewhereRadField(N_T, NULL)
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

        j = 0
        # use this to decide whether or not to compute parallel:
        # Does the local vicinity of the parallel contain radiating material?
        while j < cellArea.shape[1]:
            if CELL_RADIATES[i,j] == 1:
                J = j
                break
            j = j + 1

        if j == cellArea.shape[1]:
            continue

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

        InvisFlag[T] = 2
        InvisPhase[twoT] = 0
        InvisPhase[twoT + 1] = 0
        InvisStep[T] = 0.0

        cos_theta_i = cos(theta[i,0])
        sin_theta_i = sin(theta[i,0])
        theta_i_over_pi = theta[i,0] / _pi
        beta = radius * omega * sin_theta_i / (c * Grav_z)
        beta_sq = beta * beta
        Lorentz = sqrt(1.0 - beta_sq)

        correction_I_E = 0.0

        for k in range(N_L):
            cos_psi = cos_i * cos_theta_i + sin_i * sin_theta_i * cos(leaves[k])
            psi = acos(cos_psi)
            sin_psi = sin(psi)

            if psi < maxDeflection[i]:
                if (cos_psi < interp_alpha[T].xmin or cos_psi > interp_alpha[T].xmax):
                    #printf("cos_psi: %.16e\n", cos_psi)
                    #printf("min: %.16e\n", interp_alpha[T].xmin)
                    #printf("max: %.16e\n", interp_alpha[T].xmax)
                    terminate[T] = 1
                    break
                else:
                    cos_alpha = gsl_interp_eval(interp_alpha[T], defl_ptr, alpha_ptr, cos_psi, accel_alpha[T])
                sin_alpha = sqrt(1.0 - cos_alpha * cos_alpha)
                mu = cos_alpha * cos_gamma

                if sin_psi != 0.0:
                    cos_delta = (cos_i - cos_theta_i * cos_psi) / (sin_theta_i * sin_psi)
                    if theta_i_over_pi < 0.5:
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

                    _Z = eta * Grav_z
                    _ABB = mu * eta
                    _GEOM = mu * deriv * Grav_z * eta * eta * eta / superlum

                    if (cos_psi < interp_lag[T].xmin or cos_psi > interp_lag[T].xmax):
                        #printf("lag: %.16e\n", cos_psi)
                        #printf("min: %.16e\n", interp_lag[T].xmin)
                        #printf("max: %.16e\n", interp_lag[T].xmax)
                        terminate[T] = 1
                        break
                    else:
                        PHASE[T][k] = leaves[k] + gsl_interp_eval(interp_lag[T], defl_ptr, lag_ptr, cos_psi, accel_lag[T])

                    for p in range(N_E):
                        E_prime = energies[p] / _Z

                        I_E = eval_hotRadField(T,
                                               E_prime,
                                               _ABB,
                                               &(srcCellParams[i,J,0]),
                                               data)

                        if perform_correction == 1:
                            correction_I_E = eval_elsewhereRadField(T,
                                                              E_prime,
                                                              _ABB,
                                                              &(correction[i,J,0]),
                                                              ext_data)
                            correction_I_E = correction_I_E * eval_elsewhereRadField_norm()

                        (PROFILE[T] + BLOCK[p] + k)[0] = (I_E * eval_hotRadField_norm() - correction_I_E) * _GEOM

                    # Check whether cell was visible at previous rotation step
                    if InvisFlag[T] == 1 or (k != 0 and InvisFlag[T] == 2):
                        InvisPhase[twoT + 1] = k
                        if InvisPhase[twoT] == 0:
                            PHASE[T][0] = leaves[0]
                            InvisStep[T] = (PHASE[T][k] - PHASE[T][0]) / (<double>k)
                            for m in range(1, k):
                                PHASE[T][m] = PHASE[T][m - 1] + InvisStep[T]
                        else:
                            InvisStep[T] = ((PHASE[T][k] - PHASE[T][InvisPhase[twoT] - 1])
                                                    / (<double>k - <double>InvisPhase[twoT] + 1.0))
                            for m in range(InvisPhase[twoT], k):
                                PHASE[T][m] = PHASE[T][m - 1] + InvisStep[T]

                    # Reset visibility flag
                    InvisFlag[T] = 0

                else:
                    INVIS(k,
                          N_L,
                          N_E,
                          BLOCK,
                          InvisFlag + T,
                          InvisStep + T,
                          InvisPhase + twoT,
                          PHASE[T],
                          PROFILE[T])
            else:
                INVIS(k,
                      N_L,
                      N_E,
                      BLOCK,
                      InvisFlag + T,
                      InvisStep + T,
                      InvisPhase + twoT,
                      PHASE[T],
                      PROFILE[T])

        if interp_alpha_store != NULL:
            gsl_interp_free(interp_alpha[T])
            gsl_interp_free(interp_lag[T])
            interp_alpha[T] = interp_alpha_store
            interp_lag[T] = interp_lag_store

        #with gil:
        #    print("i = %i checkpoint a" % i)
        if terminate[T] == 1:
           break
        elif terminate[T] == 0 and InvisFlag[T] != 2:
            for n in range(1, N_L):
                if PHASE[T][n] <= PHASE[T][n-1]:
                    terminate[T] = 1
                    break
            if terminate[T] == 0:
                phase_ptr = PHASE[T]
                for p in range(N_E):
                    gsl_interp_accel_reset(accel_PROFILE[T])
                    profile_ptr = PROFILE[T] + BLOCK[p]
                    gsl_interp_init(interp_PROFILE[T], phase_ptr, profile_ptr, N_L)

                    j = 0
                    while j < cellArea.shape[1] and terminate[T] == 0:
                        if CELL_RADIATES[i,j] == 1:
                            phi_shift = phi[i,j]
                            for k in range(N_P):
                                _PHASE = phases[k]
                                _PHASE_plusShift = _PHASE + phi_shift
                                if _PHASE_plusShift > PHASE[T][N_L - 1]:
                                    while _PHASE_plusShift > PHASE[T][N_L - 1]:
                                        _PHASE_plusShift = _PHASE_plusShift - _2pi
                                    if (_PHASE_plusShift  < interp_PROFILE[T].xmin or _PHASE_plusShift > interp_PROFILE[T].xmax):
                                        printf("profile_1: %.16e\n", _PHASE_plusShift)
                                        printf("min: %.16e\n", interp_PROFILE[T].xmin)
                                        printf("max: %.16e\n", interp_PROFILE[T].xmax)
                                        terminate[T] = 1
                                        break
                                    else:
                                        privateFlux[T,p,k] += cellArea[i,j] * gsl_interp_eval(interp_PROFILE[T], phase_ptr, profile_ptr, _PHASE_plusShift, accel_PROFILE[T])
                                elif _PHASE_plusShift < PHASE[T][0]:
                                    while _PHASE_plusShift < PHASE[T][0]:
                                        _PHASE_plusShift = _PHASE_plusShift + _2pi
                                    if (_PHASE_plusShift  < interp_PROFILE[T].xmin or _PHASE_plusShift > interp_PROFILE[T].xmax):
                                        printf("profile_2: %.16e\n", _PHASE_plusShift)
                                        printf("min: %.16e\n", interp_PROFILE[T].xmin)
                                        printf("max: %.16e\n", interp_PROFILE[T].xmax)
                                        terminate[T] = 1
                                        break
                                    else:
                                        privateFlux[T,p,k] += cellArea[i,j] * gsl_interp_eval(interp_PROFILE[T], phase_ptr, profile_ptr, _PHASE_plusShift, accel_PROFILE[T])
                                else:
                                    if (_PHASE_plusShift < interp_PROFILE[T].xmin or _PHASE_plusShift > interp_PROFILE[T].xmax):
                                        printf("profile_3: %.16e\n", _PHASE_plusShift)
                                        printf("min: %.16e\n", interp_PROFILE[T].xmin)
                                        printf("max: %.16e\n", interp_PROFILE[T].xmax)
                                        terminate[T] = 1
                                        break
                                    else:
                                        privateFlux[T,p,k] += cellArea[i,j] * gsl_interp_eval(interp_PROFILE[T], phase_ptr, profile_ptr, _PHASE_plusShift, accel_PROFILE[T])
                        j = j + 1
                    if terminate[T] == 1:
                        break
       # with gil:
       #     print("i = %i checkpoint b" % i)

    for i in range(N_E):
        for T in range(N_T):
            for k in range(N_P):
                flux[i,k] += privateFlux[T,i,k]

    for p in range(N_E):
        for k in range(N_P):
            flux[p,k] /= (energies[p] * keV)

    for T in range(N_T):
        gsl_interp_free(interp_alpha[T])
        gsl_interp_accel_free(accel_alpha[T])
        gsl_interp_free(interp_lag[T])
        gsl_interp_accel_free(accel_lag[T])

        free(PHASE[T])
        free(PROFILE[T])

        gsl_interp_free(interp_PROFILE[T])
        gsl_interp_accel_free(accel_PROFILE[T])

    free(interp_alpha)
    free(accel_alpha)
    free(interp_lag)
    free(accel_lag)

    free(PHASE)
    free(PROFILE)

    free(interp_PROFILE)
    free(accel_PROFILE)

    free(InvisFlag)
    free(InvisStep)
    free(InvisPhase)

    free(BLOCK)

    if spot_atmosphere:
        free(src_preload.params)
        free(src_preload.S)
        free(src_preload)

    free_hotRadField(N_T, data)

    if perform_correction == 1:
        if elsewhere_atmosphere:
            free(ext_preload.params)
            free(ext_preload.S)
            free(ext_preload)

        free_elsewhereRadField(N_T, ext_data)

    for T in range(N_T):
        if terminate[T] == 1:
            free(terminate)
            return (ERROR, None)

    free(terminate)

    return (SUCCESS, np.asarray(flux, dtype = np.double, order = 'C'))

