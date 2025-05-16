#cython: cdivision=True
#cython: boundscheck=False
#cython: nonecheck=False
#cython: wraparound=False

"""
The default atmosphere module, corresponding to an analytical blackbody. 
Returns specific intensity, which after normalisation should be 
J/cm^2/s/keV/steradian.

Functions:
    - init_hot_BB: does nothing but must exist for placeholder purposes.
    - free_hot_BB: does nothing but must exist for placeholder purposes.
    - eval_hot_BB: evaluate the specific intensty.
    - eval_hot_BB_norm: normalize.

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

cdef void* init_hot_BB(size_t numThreads, 
                       const _preloaded *const preloaded) noexcept nogil:
    """
    This function always returns NULL because no memory needs to be 
    initialised, but must exist for placeholder purposes.
    """

    if preloaded != NULL :
        printf("WARNING: Numerical atmosphere data were preloaded, \
               even though those are not used by this atmosphere extension.\n") 
    return NULL

cdef int free_hot_BB(size_t numThreads, void *const data) noexcept nogil:
    """
    This function always returns SUCCESS = 0 because no memory needs to be 
    freed, but must exist for placeholder purposes.
    """

    return SUCCESS

cdef double eval_hot_BB(size_t THREAD,
                     double E,
                     double mu,
                     const double *const VEC,
                     void *const data) noexcept nogil:
    """
    Evaluate the specific intensity of hot regions based on photon E and temp.
    mu and other parameters in VEC don't play role here.
    
    Arguments:
        THREAD (size_t): Thread ID used for parallel execution.
        E (double): Photon energy in electron rest energy for AMXPs.
        mu (double): Cosine of the ray zenith angle (angle to surface normal).
        VEC (const double *const): Pointer to variables (e.g., temperature, 
                                                         effective gravity).
        data (void *const): This will be pointing to NULL.
        
    Returns:
        double: Calculated intensity value to be normalised.
        
    Attributes:
        noexcept nogil: Indicates that the function does not throw exceptions 
                       and can be executed without the Global Interpreter Lock 
                       (GIL).
    """

    cdef double temp = k_B_over_keV * pow(10.0, VEC[0])

    return E * E * E / ( exp(E / temp) - 1.0 )

cdef double eval_hot_norm_BB() noexcept nogil:
    """
    Source radiation field normalisation which is independent of the parameters
    of the parametrised model -- i.e. cell properties, energy, and angle.
    Writing the normalisation here reduces the number of operations required
    during integration. The units of the specific intensity need to be 
    J/cm^2/s/keV/steradian.
    
    Attributes:
        noexcept nogil: Indicates that the function does not throw exceptions 
                       and can be executed without the Global Interpreter Lock
                       (GIL).
    """

    return erg * Planck_dist_const
