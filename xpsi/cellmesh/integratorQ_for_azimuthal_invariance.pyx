#cython: cdivision=True
#cython: boundscheck=False
#cython: nonecheck=False
#cython: wraparound=False
#cython: embedsignature=True

""" Integrate over the image of a star on a distant observer's sky. """

from __future__ import division, print_function
import numpy as np
cimport numpy as np
from cython.parallel cimport *
from libc.math cimport M_PI, sqrt, sin, cos, acos, log10, pow, exp, fabs, ceil, log, atan2
from libc.stdlib cimport malloc, free
from libc.stdio cimport printf
import xpsi

cdef double _pi = M_PI
cdef double _hlfpi = M_PI / 2.0
cdef double _2pi = 2.0 * M_PI
cdef double keV = xpsi.global_imports._keV
cdef double c = xpsi.global_imports._c

cdef int SUCCESS = 0
cdef int ERROR = 1

cdef int VERBOSE = 1
cdef int QUIET = 0

from xpsi.cellmesh.integrator cimport (gsl_interp_eval,
                                       gsl_interp_eval_deriv,
                                       gsl_interp_alloc,
                                       gsl_interp_accel,
                                       gsl_interp_accel_alloc,
                                       gsl_interp_steffen,
                                       gsl_interp,
                                       gsl_interp_init,
                                       gsl_interp_free,
                                       gsl_interp_accel_free,
                                       gsl_interp_accel_reset,
                                       gsl_isnan,
                                       gsl_isinf)

ctypedef gsl_interp_accel accel
ctypedef gsl_interp interp

from .rays cimport eval_image_deflection, invert, link_rayXpanda

from xpsi.surface_radiation_field.preload cimport (_preloaded,
                                                   init_preload,
                                                   free_preload)

from xpsi.surface_radiation_field.hot cimport (init_hot,
                                               eval_hot,
                                               eval_hot_norm,
                                               eval_hot_PD,
                                               eval_hot_bbeam,
                                               free_hot)

from xpsi.surface_radiation_field.elsewhere cimport (init_elsewhere,
                                                     free_elsewhere,
                                                     eval_elsewhere,
                                                     eval_elsewhere_norm)

from ..tools cimport _get_phase_interpolant, gsl_interp_type

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
              double[:,::1] cos_alpha,
              double[:,::1] lag,
              double[::1] maxDeflection,
              double[::1] cos_gammaArray,
              double[::1] energies,
              double[::1] leaves,
              double[::1] phases,
              hot_atmosphere,
              elsewhere_atmosphere,
              image_order_limit = None):

    # check for rayXpanda explicitly in case of some linker issue
    cdef double rayXpanda_defl_lim
    cdef bint _use_rayXpanda
    try:
        xpsi.cellmesh.__deactivate_rayXpanda__
    except AttributeError:
        _try_rayXpanda = True
    else:
        if xpsi.cellmesh.__deactivate_rayXpanda__:
            _try_rayXpanda = False
            _use_rayXpanda = 0
        else:
            _try_rayXpanda = True
    finally:
        if _try_rayXpanda:
            link_rayXpanda(&_use_rayXpanda, &rayXpanda_defl_lim)

    cdef const gsl_interp_type *_interpolant

    _interpolant = _get_phase_interpolant()

    #----------------------------------------------------------------------->>>
    # >>> General memory allocation.
    # >>>
    #----------------------------------------------------------------------->>>
    cdef:
        signed int ii
        size_t i, j, J, k, ks, _kdx, m, p # Array indexing
        size_t T, twoT # Used globally to represent thread index
        size_t N_T = numThreads # shared
        size_t N_R = numRays # shared
        size_t N_E = energies.shape[0] # shared
        size_t N_L = leaves.shape[0] # shared
        size_t leaf_lim
        size_t N_P = phases.shape[0] # shared
        double sin_i = sin(inclination) # shared
        double cos_i = cos(inclination) # shared
        double Grav_z # Gravitational redshift; thread private (hereafter TP)
        double radius # Radial coordinate of pixel; TP
        double psi # Colatitude relative to star-observer direction; TP
        double cos_psi, sin_psi, _i, _cos_i, _sin_i
        double deriv # $\frac{d\cos\alpha}{d\cos\phi}$; TP
        double beta # Surface velocity in the local NRF; TP
        double _cos_alpha # Emission angle w.r.t outward radial direction in CRF; TP
        double sin_alpha
        double cos_delta # Unit sphere trigonometric identity; TP
        double cos_gamma # Surface normal tilt w.r.t outward radial direction; TP
        double cos_xi # Emission angle relative to the NRF surface velocity; TP
        double eta # Doppler boost factor in the NRF; TP
        double mu # NRF and then CRF emission angle w.r.t surface normal; TP
        double E_prime # Photon energy in the CRF, given an energy at infinity; TP
        double I_E # Radiant and or spectral intensity in the CRF; TP
        double PD # Polarization degree in the CRF; TP
        double Q_obs # Stokes Q in the observer frame; TP
        double U_obs # Stokes Q in the observer frame; TP
        double sin_chi_0, cos_chi_0, chi_0, chi_1, chi_prime, chi
        double sin_chi_1, cos_chi_1, sin_chi_prime, cos_chi_prime
        double sin_lambda, cos_lambda, cos_eps
        double sin_2chi # sine of 2*PA; TP
        double cos_2chi # cosine of 2*PA; TP
        double _PHASE, _PHASE_plusShift, _GEOM, _Z, _ABB # TP
        double phi_shift # TP
        double superlum # TP
        double cos_gamma_sq, sin_gamma
        double cos_theta_i, sin_theta_i
        double theta_i_over_pi
        double beta_sq
        double Lorentz
        double correction_Q
        int I, image_order, _IO
        double _phase_lag
        double _specific_flux
        size_t _InvisPhase

        double[:,:,::1] privateFlux = np.zeros((N_T, N_E, N_P), dtype = np.double)
        double[:,::1] flux = np.zeros((N_E, N_P), dtype = np.double)

        int *terminate = <int*> malloc(N_T * sizeof(int))

        int *InvisFlag = <int*> malloc(N_T * sizeof(int))
        double *InvisStep = <double*> malloc(N_T * sizeof(double))

        accel **accel_alpha = <accel**> malloc(N_T * sizeof(accel*))
        accel **accel_alpha_alt = NULL
        accel **accel_lag = <accel**> malloc(N_T * sizeof(accel*))
        interp **interp_alpha = <interp**> malloc(N_T * sizeof(interp*))
        interp **interp_alpha_alt = NULL
        interp **interp_lag = <interp**> malloc(N_T * sizeof(interp*))

        # Intensity spline interpolation
        double **PHASE = <double**> malloc(N_T * sizeof(double*))
        double **PROFILE = <double**> malloc(N_T * sizeof(double*))

        accel **accel_PROFILE = <accel**> malloc(N_T * sizeof(accel*))
        interp **interp_PROFILE = <interp**> malloc(N_T * sizeof(interp*))

        size_t *BLOCK = <size_t*> malloc(N_E * sizeof(size_t))

        double *defl_ptr
        double *alpha_ptr
        double *defl_alt_ptr
        double *alpha_alt_ptr
        double *lag_ptr
        double *profile_ptr
        double *phase_ptr

    if not _use_rayXpanda:
        accel_alpha_alt = <accel**> malloc(N_T * sizeof(accel*))
        interp_alpha_alt = <interp**> malloc(N_T * sizeof(interp*))

    for T in range(N_T):
        terminate[T] = 0
        accel_alpha[T] = gsl_interp_accel_alloc()
        interp_alpha[T] = gsl_interp_alloc(gsl_interp_steffen, N_R)
        if not _use_rayXpanda:
            accel_alpha_alt[T] = gsl_interp_accel_alloc()
        accel_lag[T] = gsl_interp_accel_alloc()
        interp_lag[T] = gsl_interp_alloc(gsl_interp_steffen, N_R)

        PHASE[T] = <double*> malloc(N_L * sizeof(double))
        PROFILE[T] = <double*> malloc(N_E * N_L * sizeof(double))
        accel_PROFILE[T] = gsl_interp_accel_alloc()
        interp_PROFILE[T] = gsl_interp_alloc(_interpolant, N_L)

    for p in range(N_E):
        BLOCK[p] = p * N_L

    cdef double[:,::1] cos_alpha_alt
    cdef double[:,::1] cos_deflection
    if not _use_rayXpanda:
        cos_deflection = np.zeros((deflection.shape[0],
                                   deflection.shape[1]),
                                   dtype = np.double)
        for i in range(deflection.shape[0]):
            for j in range(deflection.shape[1]):
                cos_deflection[i,j] = cos(deflection[i, N_R - j - 1])

        cos_alpha_alt = np.zeros((cos_alpha.shape[0],
                                  cos_alpha.shape[1]),
                                  dtype = np.double)

        for i in range(cos_alpha.shape[0]):
            for j in range(cos_alpha.shape[1]):
                cos_alpha_alt[i,j] = cos_alpha[i, N_R - j - 1]

    if image_order_limit is not None:
        image_order = image_order_limit
    else:
        image_order = 0

    if N_L%2 == 0:
        leaf_lim = N_L / 2
    else:
        leaf_lim = (N_L + 1)/2

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

        gsl_interp_accel_reset(accel_alpha[T])
        gsl_interp_accel_reset(accel_lag[T])

        if not _use_rayXpanda:
            j = 0
            while deflection[i,j] <= _hlfpi:
                j = j + 1

            defl_alt_ptr = &(cos_deflection[i, N_R - j - 1])
            alpha_alt_ptr = &(cos_alpha_alt[i, N_R - j - 1])
            interp_alpha_alt[T] = gsl_interp_alloc(gsl_interp_steffen, j + 1)
            gsl_interp_init(interp_alpha_alt[T], defl_alt_ptr, alpha_alt_ptr, j + 1)
            gsl_interp_accel_reset(accel_alpha_alt[T])

        defl_ptr = &(deflection[i,0])
        alpha_ptr = &(cos_alpha[i,0])
        gsl_interp_init(interp_alpha[T], defl_ptr, alpha_ptr, N_R)

        lag_ptr = &(lag[i,0])
        gsl_interp_init(interp_lag[T], defl_ptr, lag_ptr, N_R)

        radius = radialCoords_of_parallels[i]
        Grav_z = sqrt(1.0 - r_s_over_r[i])
        cos_gamma = cos_gammaArray[i]
        cos_gamma_sq = cos_gamma * cos_gamma
        sin_gamma = sqrt(1.0 - cos_gamma_sq)

        cos_theta_i = cos(theta[i,0])
        sin_theta_i = sin(theta[i,0])
        theta_i_over_pi = theta[i,0] / _pi
        beta = radius * omega * sin_theta_i / (c * Grav_z)
        beta_sq = beta * beta
        Lorentz = sqrt(1.0 - beta_sq) #inverse of Lorentz factor
        _cos_alpha = -1.0 # lastprivate
        deriv = -1.0 # lastprivate

        if image_order == 0: # infer maximum possible image order
            _IO = <int>ceil(maxDeflection[i] / _pi)
            # explanation: image_order = 1 means primary image only
        else:
            _IO = image_order
        for I in range(_IO): # loop over images
            InvisFlag[T] = 2 # initialise image order as not visible
            correction_Q = 0.0

            for k in range(leaf_lim):
                #printf("leaves = %.6e\n", leaves[k])
                #printf("k = %d\n", k)
                cos_psi = cos_i * cos_theta_i + sin_i * sin_theta_i * cos(leaves[k])
                psi = eval_image_deflection(I, acos(cos_psi))
                sin_psi = sin(psi)

                if psi != 0.0 and sin_psi == 0.0: # singularity at poles
                    # hack bypass by slight change of viewing angle
                    #printf("sinpsi = %.6e\n", sin_psi)
                    #printf("psi = %.6e\n", psi)
                    if cos_i >= 0.0:
                        _i = inclination + inclination * 1.0e-6 # arbitrary small
                        _cos_i = cos(_i)
                        _sin_i = sin(_i)
                    else:
                        _i = inclination - inclination * 1.0e-6
                        _cos_i = cos(_i)
                        _sin_i = sin(_i)

                    cos_psi =  _cos_i * cos_theta_i + _sin_i * sin_theta_i * cos(leaves[k])
                    psi = eval_image_deflection(I, acos(cos_psi))
                    sin_psi = sin(psi)
                #if sin_psi <= 1.0e-2: #Examine more this ...
                #    printf("sinpsi = %.6e\n", sin_psi)
                #    printf("psi = %.6e\n", psi)

                if psi <= maxDeflection[i]:
                    if (psi < interp_alpha[T].xmin or psi > interp_alpha[T].xmax):
                        # some crude diagnostic output
                        printf("Interpolation error: deflection = %.16e\n", psi)
                        printf("Out of bounds: min = %.16e\n", interp_alpha[T].xmin)
                        printf("Out of bounds: max = %.16e\n", interp_alpha[T].xmax)
                        terminate[T] = 1
                        break # out of phase loop
                    else:
                        if _use_rayXpanda and psi <= rayXpanda_defl_lim:
                            invert(cos_psi, r_s_over_r[i], &_cos_alpha, &deriv)
                            deriv = deriv * (1.0 - r_s_over_r[i])
                        elif not _use_rayXpanda and psi <= _hlfpi and cos_psi >= interp_alpha_alt[T].xmin:
                            _cos_alpha = gsl_interp_eval(interp_alpha_alt[T], defl_alt_ptr, alpha_alt_ptr, cos_psi, accel_alpha_alt[T])
                        else:
                            _cos_alpha = gsl_interp_eval(interp_alpha[T], defl_ptr, alpha_ptr, psi, accel_alpha[T])

                    sin_alpha = sqrt(1.0 - _cos_alpha * _cos_alpha)
                    mu = _cos_alpha * cos_gamma

                    # for spherical stars mu is defined, but for tilted local
                    # surface, there is not one unique value for mu because
                    # Einstein ring(s)
                    if psi != 0.0:
                        cos_delta = (cos_i - cos_theta_i * cos_psi) / (sin_theta_i * sin_psi)
                        if theta_i_over_pi < 0.5:
                            mu = mu + sin_alpha * sin_gamma * cos_delta
                        else:
                            mu = mu - sin_alpha * sin_gamma * cos_delta

                    if mu > 0.0:
                        if _use_rayXpanda and psi <= rayXpanda_defl_lim:
                            pass
                        elif not _use_rayXpanda and psi <= _hlfpi and cos_psi >= interp_alpha_alt[T].xmin:
                            deriv = gsl_interp_eval_deriv(interp_alpha_alt[T], defl_alt_ptr, alpha_alt_ptr, cos_psi, accel_alpha_alt[T])
                        else:
                            deriv = gsl_interp_eval_deriv(interp_alpha[T], defl_ptr, alpha_ptr, psi, accel_alpha[T])
                            deriv = exp( log(fabs(deriv)) - log(fabs(sin_psi)) ) # singularity hack above

                        if (psi < interp_lag[T].xmin or psi > interp_lag[T].xmax):
                            printf("Interpolation error: deflection = %.16e\n", psi)
                            printf("Out of bounds: min = %.16e\n", interp_lag[T].xmin)
                            printf("Out of bounds: max = %.16e\n", interp_lag[T].xmax)
                            terminate[T] = 1
                            break # out of phase loop
                        else:
                            _phase_lag = gsl_interp_eval(interp_lag[T], defl_ptr, lag_ptr, psi, accel_lag[T])

                        for ks in range(2):
                            if (0 < k < leaf_lim - 1
                                    or (k == 0 and ks == 0)
                                    or (k == leaf_lim - 1 and N_L%2 == 1 and ks == 0)
                                    or (k == leaf_lim - 1 and N_L%2 == 0)):

                                if ks == 0:
                                    _kdx = k
                                else:
                                    _kdx = N_L - 1 - k # switch due to symmetry

                                # phase asymmetric now
                                if psi != 0.0:
                                    cos_xi = sin_alpha * sin_i * sin(leaves[_kdx]) / sin_psi
                                    superlum = (1.0 + beta * cos_xi)
                                    eta = Lorentz / superlum
                                else:
                                    cos_xi = 0.0
                                    superlum = 1.0
                                    eta = Lorentz

                                _Z = eta * Grav_z
                                _ABB = mu * eta
                                _GEOM = mu * fabs(deriv) * Grav_z * eta * eta * eta / superlum

                                PHASE[T][_kdx] = leaves[_kdx] + _phase_lag

                                #sin_2chi, cos_2chi = eval_PA()
                                #Or calculaing the PA here in place:
                                sin_chi_0 = - sin_theta_i*sin(leaves[_kdx]) 
                                cos_chi_0 = sin_i*cos_theta_i - sin_theta_i*cos_i*cos(leaves[_kdx])
                                chi_0 = atan2(sin_chi_0,cos_chi_0)

                                #Notes: mu = cos_sigma , Lorentz = 1/Gamma, mu0=eta*mu, cos_xi defined with no minus sign
                                #TBD: how to get sinalpha/sinpsi as sin_alpha_over_sin_psi (when sinpsi -> 0)?
                                sin_chi_1 = sin_gamma*sin_i*sin(leaves[_kdx])*sin_alpha/sin_psi #times sin alpha sin sigma
                                cos_chi_1 = cos_gamma - _cos_alpha*mu  #times sin alpha sin sigma 
                                chi_1 = atan2(sin_chi_1,cos_chi_1)

                                sin_lambda = sin_theta_i*cos_gamma - sin_gamma*cos_theta_i
                                cos_lambda = cos_theta_i*cos_gamma + sin_theta_i*sin_gamma
                                cos_eps = (sin_alpha/sin_psi)*(cos_i*sin_lambda - sin_i*cos_lambda*cos(leaves[_kdx]) + cos_psi*sin_gamma) - _cos_alpha*sin_gamma

                                sin_chi_prime = cos_eps*eta*mu*beta/Lorentz
                                cos_chi_prime = (1. - mu**2 /(1. + beta*cos_xi))
                                chi_prime = atan2(sin_chi_prime,cos_chi_prime)

                                chi = chi_0+chi_1+chi_prime

                                #printf("leaves[k] = %.6e ",leaves[k])
                                #printf("chi_0 = %.6e ",chi_0)
                                #printf("chi_1 = %.6e ",chi_1)
                                #printf("chi_prime = %.6e ",chi_prime)
                                #printf("PA_tot = %.6e\n",chi)
                                cos_2chi = cos(2*chi)
                                sin_2chi = sin(2*chi)

                                # specific intensities
                                for p in range(N_E):
                                    E_prime = energies[p] / _Z

                                    I_E = eval_hot_bbeam(T,
                                                   E_prime,
                                                   _ABB,
                                                   &(srcCellParams[i,J,0]),
                                                   hot_data)

                                    PD = eval_hot_PD(T,
                                                   E_prime,
                                                   _ABB,
                                                   &(srcCellParams[i,J,0]),
                                                   hot_data)


                                    Q_obs = PD*I_E*cos_2chi
                                    #U_obs = PD*I_E*sin_2chi
                                    #printf("Q_obs = %.6e\n",Q_obs)

                                    if perform_correction == 1:
                                        correction_Q = eval_elsewhere(T,
                                                                        E_prime,
                                                                        _ABB,
                                                                        &(correction[i,J,0]),
                                                                        ext_data)
                                        correction_Q = correction_Q * eval_elsewhere_norm()

                                    (PROFILE[T] + BLOCK[p] + _kdx)[0] = (Q_obs * eval_hot_norm() - correction_Q) * _GEOM
                                    #printf("Q_obs = %.6e\n",(Q_obs * eval_hot_norm() - correction_Q) * _GEOM)

                        if k == 0: # if initially visible at first/last phase steps
                            # periodic
                            PHASE[T][N_L - 1] = PHASE[T][0] + _2pi
                            for p in range(N_E):
                                (PROFILE[T] + BLOCK[p] + N_L - 1)[0] = (PROFILE[T] + BLOCK[p])[0]
                        elif k > 0 and InvisFlag[T] == 2: # initially not visible
                            # calculate the appropriate phase increment for
                            # phase steps through non-visible fraction of cycle
                            InvisStep[T] = leaves[k] / <double>k

                            # increment phase from start to end of non-visible
                            # interval
                            # first up to the periodic boundary
                            for m in range(N_L - k, N_L):
                                PHASE[T][m] = PHASE[T][m - 1] + InvisStep[T]

                            # set the specific intensities to zero
                            for p in range(N_E):
                                for m in range(N_L - k, N_L):
                                    (PROFILE[T] + BLOCK[p] + m)[0] = 0.0

                            # handle the duplicate points at the periodic
                            # boundary which are needed for interpolation
                            PHASE[T][0] = PHASE[T][N_L - 1] - _2pi

                            # now after the periodic boundary up to the step
                            # where image becomes visible
                            for m in range(1, k):
                                PHASE[T][m] = PHASE[T][m - 1] + InvisStep[T]

                            # set the reminaing specific intensities to zero
                            for p in range(N_E):
                                (PROFILE[T] + BLOCK[p])[0] = 0.0
                                for m in range(1, k):
                                    (PROFILE[T] + BLOCK[p] + m)[0] = 0.0
                        elif InvisFlag[T] == 1: # handle linearly spaced phases
                            InvisStep[T] = PHASE[T][k] - PHASE[T][_InvisPhase - 1]
                            InvisStep[T] = InvisStep[T] / <double>(k - _InvisPhase + 1)

                            for m in range(_InvisPhase, k):
                                PHASE[T][m] = PHASE[T][m - 1] + InvisStep[T]

                            InvisStep[T] = PHASE[T][N_L - _InvisPhase] - PHASE[T][N_L - 1 - k]
                            InvisStep[T] = InvisStep[T] / <double>(k - _InvisPhase + 1)

                            for m in range(N_L - k, N_L - _InvisPhase):
                                PHASE[T][m] = PHASE[T][m - 1] + InvisStep[T]

                        # reset visibility flag
                        InvisFlag[T] = 0

                    else:
                        # check whether cell was visible at previous phase step
                        if InvisFlag[T] == 0:
                            # if image was visible, calculate the appropriate
                            # phase step for the fraction of the cycle when
                            # image is not visible
                            InvisStep[T] = PHASE[T][N_L - k] - PHASE[T][k - 1]
                            InvisStep[T] = InvisStep[T] / <double>(N_L - 2*k + 1)

                            # step in phase between the phases at which image
                            # is visible
                            for m in range(k, N_L - k):
                                PHASE[T][m] = PHASE[T][m - 1] + InvisStep[T]

                            # set the specific intensities to zero when image
                            # is not visible
                            for p in range(N_E):
                                for m in range(k, N_L - k):
                                    (PROFILE[T] + BLOCK[p] + m)[0] = 0.0

                            InvisFlag[T] = 1 # declare not visible
                            _InvisPhase = k
                else:
                    if InvisFlag[T] == 0:
                        InvisStep[T] = PHASE[T][N_L - k] - PHASE[T][k - 1]
                        InvisStep[T] = InvisStep[T] / <double>(N_L - 2*k + 1)

                        for m in range(k, N_L - k):
                            PHASE[T][m] = PHASE[T][m - 1] + InvisStep[T]

                        for p in range(N_E):
                            for m in range(k, N_L - k):
                                (PROFILE[T] + BLOCK[p] + m)[0] = 0.0

                        InvisFlag[T] = 1
                        _InvisPhase = k

            if terminate[T] == 1:
                break # out of image loop
            elif InvisFlag[T] == 2: # no visibility detected
                break # ignore higher order images, assume no visiblity
            else: # proceed to sum over images
                for m in range(1, N_L):
                    if PHASE[T][m] <= PHASE[T][m - 1]:
                        printf("Interpolation error: phases are not strictly increasing.\n")
                        printf('%.8e -> %.8e\n', PHASE[T][m - 1], PHASE[T][m])
                        terminate[T] = 1
                        break # out of phase loop
                if terminate[T] == 1:
                    break # out of image loop
                else:
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
                                    elif _PHASE_plusShift < PHASE[T][0]:
                                        while _PHASE_plusShift < PHASE[T][0]:
                                            _PHASE_plusShift = _PHASE_plusShift + _2pi

                                    if (_PHASE_plusShift < interp_PROFILE[T].xmin or _PHASE_plusShift > interp_PROFILE[T].xmax):
                                        printf("Interpolation error: phase = %.16e\n", _PHASE_plusShift)
                                        printf("Out of bounds: min = %.16e\n", interp_PROFILE[T].xmin)
                                        printf("Out of bounds: max = %.16e\n", interp_PROFILE[T].xmax)
                                        terminate[T] = 1
                                        break # out of phase loop

                                    _specific_flux = gsl_interp_eval(interp_PROFILE[T], phase_ptr, profile_ptr, _PHASE_plusShift, accel_PROFILE[T])
                                    #if _specific_flux > 0.0 or perform_correction == 1:
                                    if True or perform_correction == 1:
                                        privateFlux[T,p,k] += cellArea[i,j] * _specific_flux

                            j = j + 1
                        if terminate[T] == 1:
                            break # out of energy loop
            if terminate[T] == 1:
                break # out of image loop
        if not _use_rayXpanda:
            gsl_interp_free(interp_alpha_alt[T])
        if terminate[T] == 1:
           break # out of colatitude loop

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
        if not _use_rayXpanda:
            gsl_interp_accel_free(accel_alpha_alt[T])
        gsl_interp_free(interp_lag[T])
        gsl_interp_accel_free(accel_lag[T])

        free(PHASE[T])
        free(PROFILE[T])

        gsl_interp_free(interp_PROFILE[T])
        gsl_interp_accel_free(accel_PROFILE[T])

    free(interp_alpha)
    free(accel_alpha)
    if not _use_rayXpanda:
        free(interp_alpha_alt)
        free(accel_alpha_alt)
    free(interp_lag)
    free(accel_lag)

    free(PHASE)
    free(PROFILE)

    free(interp_PROFILE)
    free(accel_PROFILE)

    free(InvisFlag)
    free(InvisStep)

    free(BLOCK)

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

    free(terminate)

    #printf("FluxQ = %.6e\n", flux)

    return (SUCCESS, np.asarray(flux, dtype = np.double, order = 'C'))


#sin_chi_0= - sin_theta*sin_phi # times sin psi
#cos_chi_0=sin_i*cos_theta - sin_theta*cos_i*cos_phi # times sin psi 
#chi_0=arctan2(sin_chi_0,cos_chi_0)

#sin_chi_1=sin_gamma*sin_i*sin_phi*sin_alpha_over_sin_psi #times sin alpha sin sigma 
#cos_chi_1=cos_gamma - cos_alpha*cos_sigma  #times sin alpha sin sigma 
#chi_1=arctan2(sin_chi_1,cos_chi_1)

#sin_lambda=sin_theta*cos_gamma - sin_gamma*cos_theta
#cos_lambda=cos_theta*cos_gamma + sin_theta*sin_gamma
#cos_eps = sin_alpha_over_sin_psi*(cos_i*sin_lambda - sin_i*cos_lambda*cos_phi + cos_psi*sin_gamma) - cos_alpha*sin_gamma

#sin_chi_prime=cos_eps*mu0*Gamma*beta#*the_thing#*delta**3#*delta
#cos_chi_prime=(1. - cos_sigma**2 /(1. - beta*cos_xi))#*the_thing
#chi_prime=arctan2(sin_chi_prime,cos_chi_prime)   
#chi=chi_0 + chi_1 + chi_prime
#  print(chi,'\t',chi_0/chi,'\t',chi_1/chi ,'\t',  chi_prime/chi )

#sin_2chi=sin(2*chi)
#cos_2chi=cos(2*chi)


     
