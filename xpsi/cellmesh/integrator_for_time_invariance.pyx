#cython: cdivision=True
#cython: boundscheck=False
#cython: nonecheck=False
#cython: wraparound=False

from __future__ import division, print_function
import numpy as np
cimport numpy as np
from cython.parallel cimport *
from libc.math cimport M_PI, sqrt, sin, cos, acos
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

from xpsi.surface_radiation_field.elsewhere_radiation_field cimport (init_elsewhereRadField,
                                                                eval_elsewhereRadField,
                                                                eval_elsewhereRadField_norm,
                                                                free_elsewhereRadField,
                                                                elsewhereRadField_PRELOAD)

#----------------------------------------------------------------------->>>
# >>> Integrate over the celestial sphere of distant observer.
# >>> 
#----------------------------------------------------------------------->>>
def integrate_radField(size_t numThreads,
                         double R,
                         double omega,
                         double r_s,
                         double inclination,
                         size_t sqrt_numPix,
                         double cellArea,
                         double[::1] radialCoords_of_parallels,
                         double[::1] r_s_over_r,
                         double[:,::1] theta,
                         double[:,::1] phi,
                         double[:,:,::1] srcCellParams,
                         int numRays,
                         double[:,::1] deflection,
                         double[:,::1] cos_alphaMatrix,
                         double[::1] maxDeflection,
                         double[::1] cos_gammaArray,
                         double[::1] energies,
                         *args):

    #----------------------------------------------------------------------->>>
    # >>> General memory allocation.
    # >>>
    #----------------------------------------------------------------------->>>
    cdef:
        signed int ii
        size_t i, j, e # Array indexing
        size_t T # Used globally to represent thread index
        size_t N_T = numThreads # shared
        size_t N_R = numRays # shared
        size_t N_E = energies.shape[0] # shared
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
        double _GEOM, _Z, _ABB # TP
        double superlum # TP
        double cos_gamma_sq, sin_gamma_sq
        double cos_theta_ij, sin_theta_ij
        double sqrt_cos_gamma_sq
        double theta_ij_over_pi
        double beta_sq
        double Lorentz

        double[:,::1] privateFlux = np.zeros((N_T, N_E), dtype = np.double)
        double[::1] flux = np.zeros(N_E, dtype = np.double)

        int *terminate = <int*> malloc(N_T * sizeof(int))

        accel** accel_alpha = <accel**> malloc(N_T * sizeof(accel*))
        interp** interp_alpha = <interp**> malloc(N_T * sizeof(interp*))

        interp *interp_alpha_store

        double *defl_ptr
        double *alpha_ptr

    for T in range(N_T):
        terminate[T] = 0
        accel_alpha[T] = gsl_interp_accel_alloc()
        interp_alpha[T] = gsl_interp_alloc(gsl_interp_steffen, N_R)

    cdef double[:,::1] cos_deflection = np.zeros((deflection.shape[0],
                                                  deflection.shape[1]),
                                                 dtype = np.double)

    for i in range(deflection.shape[0]):
        for j in range(deflection.shape[1]):
            cos_deflection[i,j] = cos(deflection[i,j])

    # Initialise the source radiation field
    cdef elsewhereRadField_PRELOAD *preload = NULL
    cdef double[::1] cast
    cdef double[::1] intensity
    cdef void *data = NULL
    cdef size_t num_args
    if args:
        num_args = len(args)
        preload = <elsewhereRadField_PRELOAD*> malloc(sizeof(elsewhereRadField_PRELOAD))
        preload.params = <double**> malloc(sizeof(double*) * (num_args - 1))
        preload.S = <size_t*> malloc(sizeof(size_t) * (num_args - 2))
        for i in range(num_args - 1):
            cast = args[i]
            preload.params[i] = &cast[0]
            if i < num_args - 2:
                cast = args[i+1]
                preload.S[i] = cast.shape[0]
                if i < num_args - 3:
                    for j in range(i+2, num_args - 1):
                        cast = args[j]
                        preload.S[i] *= cast.shape[0]
        intensity = args[i+1]
        preload.I = &intensity[0]
        data = init_elsewhereRadField(N_T, preload)
    else:
        data = init_elsewhereRadField(N_T, NULL)

    #----------------------------------------------------------------------->>>
    # >>> Integrate.
    # >>>
    #----------------------------------------------------------------------->>>
    for ii in prange(<signed int>sqrt_numPix,
                     nogil = True,
                     schedule = 'static',
                     num_threads = N_T,
                     chunksize = 1):

        T = threadid()
        i = <size_t> ii

        radius = radialCoords_of_parallels[i]
        Grav_z = sqrt(1.0 - r_s_over_r[i])
        cos_gamma = cos_gammaArray[i]
        cos_gamma_sq = cos_gamma * cos_gamma
        sin_gamma_sq = sqrt(1.0 - cos_gamma_sq)

        gsl_interp_accel_reset(accel_alpha[T])

        j = 0
        while deflection[i,j] > _pi:
            j = j + 1

        if j != 0:
            interp_alpha_store = interp_alpha[T]
            interp_alpha[T] = gsl_interp_alloc(gsl_interp_steffen, N_R - j)

            defl_ptr = &(cos_deflection[i,j])
            alpha_ptr = &(cos_alphaMatrix[i,j])
            gsl_interp_init(interp_alpha[T], defl_ptr, alpha_ptr, N_R - j)
        else:
            interp_alpha_store = NULL

            defl_ptr = &(cos_deflection[i,0])
            alpha_ptr = &(cos_alphaMatrix[i,0])
            gsl_interp_init(interp_alpha[T], defl_ptr, alpha_ptr, N_R)

        cos_theta_ij = cos(theta[i,0])
        sin_theta_ij = sin(theta[i,0])
        theta_ij_over_pi = theta[i,0] / _pi
        beta = radius * omega * sin_theta_ij / (c * Grav_z)
        beta_sq = beta * beta
        Lorentz = sqrt(1.0 - beta_sq)

        for j in range(sqrt_numPix):
            cos_psi = cos_i * cos_theta_ij + sin_i * sin_theta_ij * cos(phi[i,j])
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
                        cos_xi = sin_alpha * sin_i * sin(phi[i,j]) / sin_psi
                        superlum = (1.0 + beta * cos_xi)
                        eta = Lorentz / superlum
                    else:
                        cos_xi = 0.0
                        superlum = 1.0
                        eta = Lorentz

                    deriv = gsl_interp_eval_deriv(interp_alpha[T], defl_ptr, alpha_ptr, cos_psi, accel_alpha[T])

                    _Z = eta * Grav_z
                    _ABB = mu * eta
                    _GEOM = mu * deriv * Grav_z * eta * eta * eta / superlum

                    for e in range(N_E):
                        E_prime = energies[e] / _Z

                        I_E = eval_elsewhereRadField(T,
                                               E_prime,
                                               _ABB,
                                               &(srcCellParams[i,j,0]),
                                               data)

                        privateFlux[T,e] += I_E * _GEOM

        if interp_alpha_store != NULL:
            gsl_interp_free(interp_alpha[T])
            interp_alpha[T] = interp_alpha_store

        if terminate[T] == 1:
            break

    for T in range(N_T):
        gsl_interp_free(interp_alpha[T])
        gsl_interp_accel_free(accel_alpha[T])
        for e in range(N_E):
            flux[e] += privateFlux[T,e]

    for e in range(N_E):
        flux[e] *= cellArea * eval_elsewhereRadField_norm() / (energies[e] * keV)

    free(interp_alpha)
    free(accel_alpha)

    if args:
        free(preload.params)
        free(preload.S)
        free(preload)

    free_elsewhereRadField(N_T, data)

    for T in range(N_T):
        if terminate[T] == 1:
            free(terminate)
            return (ERROR, None)

    free(terminate)

    return (SUCCESS, np.asarray(flux, dtype = np.double, order = 'C'))

