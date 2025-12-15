#cython: cdivision=True
#cython: boundscheck=False
#cython: nonecheck=False
#cython: wraparound=False
#cython: embedsignature=True

""" Integrate over the image of a star on a distant observer's sky. """

import numpy as np
cimport numpy as np
from cython.parallel cimport *
from libc.math cimport M_PI, sqrt, sin, cos, acos, fabs, pow, ceil, log, exp
from libc.stdlib cimport malloc, free
from libc.stdio cimport printf, setbuf, stdout
import xpsi

cdef double _pi = M_PI
cdef double _hlfpi = M_PI / 2.0
cdef double _2pi = 2.0 * M_PI
cdef double keV = xpsi.global_imports._keV
cdef double c = xpsi.global_imports._c

cdef int SUCCESS = 0
cdef int ERROR = 1

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

from .rays cimport eval_image_deflection
from ..tools.core cimport are_equal

from xpsi.surface_radiation_field.preload cimport (_preloaded,
                                                   init_preload,
                                                   free_preload)

from xpsi.surface_radiation_field.elsewhere_wrapper cimport (init_elsewhere,
                                                     free_elsewhere,
                                                     eval_elsewhere,
                                                     eval_elsewhere_norm)

#----------------------------------------------------------------------->>>
# >>> Integrate over the celestial sphere of distant observer.
# >>> 
#----------------------------------------------------------------------->>>
def integrate(size_t numThreads,
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
              double[:,::1] cos_alpha,
              double[::1] maxDeflection,
              double[::1] cos_gammaArray,
              double[::1] energies,
              atmosphere,
              atm_ext,
              image_order_limit = None,
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
        double psi, _psi # Colatitude relative to star-observer direction; TP
        double cos_psi, sin_psi, _cos_psi, _i, _cos_i, _sin_i
        double deriv # $\frac{d\cos\alpha}{d\cos\psi}$; TP
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
        double _GEOM, _Z, _ABB # TP
        double superlum # TP
        double cos_gamma_sq, sin_gamma
        double cos_theta_i, sin_theta_i
        double theta_i_over_pi
        double beta_sq
        double Lorentz
        int I, image_order, _IO

        double[:,::1] privateFlux = np.zeros((N_T, N_E), dtype = np.double)
        double[::1] flux = np.zeros(N_E, dtype = np.double)

        int *terminate = <int*> malloc(N_T * sizeof(int))

        accel **accel_alpha = <accel**> malloc(N_T * sizeof(accel*))
        interp **interp_alpha = <interp**> malloc(N_T * sizeof(interp*))
        accel **accel_alpha_alt = NULL
        interp **interp_alpha_alt = NULL

        double *defl_ptr
        double *alpha_ptr
        double *defl_alt_ptr
        double *alpha_alt_ptr

    accel_alpha_alt = <accel**> malloc(N_T * sizeof(accel*))
    interp_alpha_alt = <interp**> malloc(N_T * sizeof(interp*))

    for T in range(N_T):
        terminate[T] = 0
        accel_alpha[T] = gsl_interp_accel_alloc()
        interp_alpha[T] = gsl_interp_alloc(gsl_interp_steffen, N_R)

        accel_alpha_alt[T] = gsl_interp_accel_alloc()

    cdef double[:,::1] cos_alpha_alt
    cdef double[:,::1] cos_deflection

    cos_deflection = np.zeros((deflection.shape[0],
                               deflection.shape[1]),
                               dtype = np.double)
    for i in range(<size_t>deflection.shape[0]):
        for j in range(<size_t>deflection.shape[1]):
            cos_deflection[i,j] = cos(deflection[i, N_R - j - 1])

    cos_alpha_alt = np.zeros((cos_alpha.shape[0],
                              cos_alpha.shape[1]),
                              dtype = np.double)

    for i in range(<size_t>cos_alpha.shape[0]):
        for j in range(<size_t>cos_alpha.shape[1]):
            cos_alpha_alt[i,j] = cos_alpha[i, N_R - j - 1]

    if image_order_limit is not None:
        image_order = image_order_limit
    else:
        image_order = 0

    # Initialise the source radiation field
    cdef _preloaded *preloaded = NULL
    cdef void *data = NULL

    if atmosphere:
        preloaded = init_preload(atmosphere)
        data = init_elsewhere(N_T, preloaded, atm_ext)
    else:
        data = init_elsewhere(N_T, NULL, atm_ext)

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

        gsl_interp_accel_reset(accel_alpha[T])

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
        Lorentz = sqrt(1.0 - beta_sq)
        _cos_alpha = -1.0 # lastprivate
        deriv = -1.0 # lastprivate

        for j in range(sqrt_numPix):
            _cos_psi = cos_i * cos_theta_i + sin_i * sin_theta_i * cos(phi[i,j])
            _psi = acos(_cos_psi)
            if image_order == 0: # infer maximum possible image order
                _IO = <int>ceil(maxDeflection[i] / _pi)
                # explanation: image_order = 1 means primary image only
            else:
                _IO = image_order
            for I in range(_IO):
                cos_psi = _cos_psi
                psi = eval_image_deflection(I, _psi)
                sin_psi = sin(psi)
                if not are_equal(psi, 0.0) and are_equal(sin_psi, 0.0): # singularity at poles
                    # hack bypass by slight change of viewing angle
                    if cos_i >= 0.0:
                        _i = inclination + inclination * 1.0e-6 # arbitrary small
                        _cos_i = cos(_i)
                        _sin_i = sin(_i)
                    else:
                        _i = inclination - inclination * 1.0e-6
                        _cos_i = cos(_i)
                        _sin_i = sin(_i)

                    cos_psi = _cos_i * cos_theta_i + _sin_i * sin_theta_i * cos(phi[i,j])
                    psi = eval_image_deflection(I, acos(cos_psi))
                    sin_psi = sin(psi)

                if psi > maxDeflection[i]:
                    break # higher order images not visible
                else:
                    if (psi < interp_alpha[T].xmin or psi > interp_alpha[T].xmax):
                        printf("psi: %.16e\n", psi)
                        printf("min: %.16e\n", interp_alpha[T].xmin)
                        printf("max: %.16e\n", interp_alpha[T].xmax)
                        terminate[T] = 1
                        break # out of image loop
                    else:
                        if psi <= _hlfpi and cos_psi >= interp_alpha_alt[T].xmin:
                            _cos_alpha = gsl_interp_eval(interp_alpha_alt[T], defl_alt_ptr, alpha_alt_ptr, cos_psi, accel_alpha_alt[T])
                        else:
                            _cos_alpha = gsl_interp_eval(interp_alpha[T], defl_ptr, alpha_ptr, psi, accel_alpha[T])

                    sin_alpha = sqrt(1.0 - _cos_alpha * _cos_alpha)
                    mu = _cos_alpha * cos_gamma

                    if not are_equal(psi, 0.0):
                        cos_delta = (cos_i - cos_theta_i * cos_psi) / (sin_theta_i * sin_psi)
                        if theta_i_over_pi < 0.5:
                            mu = mu + sin_alpha * sin_gamma * cos_delta
                        else:
                            mu = mu - sin_alpha * sin_gamma * cos_delta

                    if mu > 0.0:
                        if not are_equal(sin_psi, 0.0):
                            cos_xi = sin_alpha * sin_i * sin(phi[i,j]) / sin_psi
                            superlum = (1.0 + beta * cos_xi)
                            eta = Lorentz / superlum
                        else:
                            cos_xi = 0.0
                            superlum = 1.0
                            eta = Lorentz

                        if psi <= _hlfpi and cos_psi >= interp_alpha_alt[T].xmin:
                            deriv = gsl_interp_eval_deriv(interp_alpha_alt[T], defl_alt_ptr, alpha_alt_ptr, cos_psi, accel_alpha_alt[T])
                        else:
                            deriv = gsl_interp_eval_deriv(interp_alpha[T], defl_ptr, alpha_ptr, psi, accel_alpha[T])
                            deriv = exp( log(fabs(deriv)) - log(fabs(sin_psi)) ) # singularity hack above

                        _Z = eta * Grav_z
                        _ABB = mu * eta
                        _GEOM = mu * fabs(deriv) * Grav_z * eta * eta * eta #/ superlum

                        for e in range(N_E):
                            E_prime = energies[e] / _Z

                            I_E = eval_elsewhere(T,
                                                 E_prime,
                                                 _ABB,
                                                 &(srcCellParams[i,j,0]),
                                                 data,
                                                 0)

                            privateFlux[T,e] += I_E * _GEOM
            if terminate[T] == 1:
                break # out of azimuth loop
        gsl_interp_free(interp_alpha_alt[T])
        if terminate[T] == 1:
            break # out of colatitude loop

    for T in range(N_T):
        gsl_interp_free(interp_alpha[T])
        gsl_interp_accel_free(accel_alpha[T])
        gsl_interp_accel_free(accel_alpha_alt[T])
        for e in range(N_E):
            flux[e] += privateFlux[T,e]

    for e in range(N_E):
        flux[e] *= cellArea * eval_elsewhere_norm() / (energies[e] * keV)

    free(interp_alpha)
    free(accel_alpha)
    free(interp_alpha_alt)
    free(accel_alpha_alt)

    if atmosphere:
        free_preload(preloaded)

    free_elsewhere(N_T, data)

    for T in range(N_T):
        if terminate[T] == 1:
            free(terminate)
            return (ERROR, None)

    free(terminate)

    return (SUCCESS, np.asarray(flux, dtype = np.double, order = 'C'))

