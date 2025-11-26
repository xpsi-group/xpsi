#cython: cdivision=True
#cython: boundscheck=False
#cython: nonecheck=False
#cython: wraparound=False

"""
A linear atmosphere model, used in Molkov et al. 2024.
Returns specific intensity, but units were never given even in that paper,
so only consider normalised pulse profiles.

Functions:
    - init_hot_user: does nothing but must exist for placeholder purposes.
    - free_hot_user: does nothing but must exist for placeholder purposes.
    - eval_hot_user_I: evaluate the specific intensity.
    - eval_hot_user_Q: evaluate the specific intensity.
    - eval_hot_user_norm: normalize.

"""


from libc.math cimport exp, pow
from libc.stdio cimport printf

from xpsi.global_imports import _keV, _k_B

cdef int SUCCESS = 0

cdef double erg = 1.0e-7
cdef double Planck_dist_const = 5.040366110812353e22

cdef double k_B = _k_B
cdef double keV = _keV
cdef double k_B_over_keV = k_B / keV

cdef void* init_hot_user(size_t numThreads, 
                       const _preloaded *const preloaded) noexcept nogil:
    """
    This function always returns NULL because no memory needs to be 
    initialised, but must exist for placeholder purposes.
    """

    if preloaded != NULL :
        printf("WARNING: Numerical atmosphere data were preloaded, \
               even though those are not used by this atmosphere extension.\n") 
    return NULL

cdef int free_hot_user(size_t numThreads, void *const data) noexcept nogil:
    """
    This function always returns SUCCESS = 0 because no memory needs to be 
    freed, but must exist for placeholder purposes.
    """

    return SUCCESS

cdef double eval_hot_user_I(size_t THREAD,
                     double E,
                     double mu,
                     const double *const VEC,
                     void *const data) noexcept nogil:
    """
    Evaluate the specific intensity of hot regions with angular emission pattern.
    I = I_E * (1.0 + h * mu)  

    Arguments:
        THREAD (size_t): Thread ID used for parallel execution.
        E (double): Photon energy in electron rest energy for AMXPs. [electron rest mass energy]
        mu (double): Cosine of the ray zenith angle (angle to surface normal).
        VEC (const double *const): Pointer to variable h [unitless]
        data (void *const): This will be pointing to NULL.

    Returns:
        The output intensity in units J/cm^2/s/keV/steradian. 
    """

   

    # Angular anisotropy parameter h (stored in VEC[0])
    #cdef double h = -0.63
    h = VEC[0]
   

    # Apply angular emission pattern
    return (1.0 + h * mu)

cdef double eval_hot_user_Q(size_t THREAD,
                     double E,
                     double mu,
                     const double *const VEC,
                     void *const data) noexcept nogil:
    """
    Evaluate the specific intensity of hot regions with angular emission pattern.
    I = I_E * (1.0 + h * mu)  

    Arguments:
        THREAD (size_t): Thread ID used for parallel execution.
        E (double): Photon energy in electron rest energy for AMXPs. [electron rest mass energy]
        mu (double): Cosine of the ray zenith angle (angle to surface normal).
        VEC (const double *const): Pointer to variable h [unitless]
        data (void *const): This will be pointing to NULL.

    Returns:
        The output intensity in units J/cm^2/s/keV/steradian. 
    """

   

    # Angular anisotropy parameter h (stored in VEC[0])
    #cdef double h = -0.63
    h = VEC[0]
   

    # Apply angular emission pattern
    return (1.0 + h * mu)


cdef double eval_hot_norm_user() noexcept nogil:
    return erg * Planck_dist_const