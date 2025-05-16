#cython: cdivision=True
#cython: boundscheck=False
#cython: nonecheck=False
#cython: wraparound=False

"""
The BB atmosphere module, but modified for Bursts. Returns specific intensity, 
which after normalisation should be J/cm^2/s/keV/steradian.

Functions:
    - init_hot_BB_burst: does nothing but must exist for placeholder purposes.
    - free_hot_BB_burst: does nothing but must exists for placeholder purposes.
    - eval_hot_BB_burst_I: evaluate the specific intensty.
    - eval_hot_BB_burst_Q: evaluate stokes Q.
    - eval_hot_BB_norm: normalize.

"""


from libc.math cimport exp, pow

from xpsi.global_imports import _keV, _k_B

cdef int SUCCESS = 0

cdef double erg = 1.0e-7
cdef double Planck_dist_const = 5.040366110812353e22

cdef double k_B = _k_B
cdef double keV = _keV
cdef double k_B_over_keV = k_B / keV

cdef void* init_hot_BB_burst(size_t numThreads, const _preloaded *const preloaded) noexcept nogil:
    """
    This function always returns NULL because no memory needs to be 
    initialised, but must exist for placeholder purposes.
    """
    
    return NULL

cdef int free_hot_BB_burst(size_t numThreads, void *const data) noexcept nogil:
    """
    This function always returns SUCCESS = 0 because no memory needs to be 
    freed, but must exist for placeholder purposes.
    """

    return SUCCESS

cdef double eval_hot_norm_BB_burst() noexcept nogil:
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
    
    
cdef double eval_hot_BB_burst_I(size_t THREAD,
                     double E,
                     double mu,
                     const double *const VEC,
                     void *const data) noexcept nogil:
    """
    Evaluate the specific intensity of hot regions based on photon E and temp
    and mu. Other parameters in VEC don't play role here.
    
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

    return E * E * E / ( exp(E / temp) - 1.0 )*(0.421+0.868*mu)    
    
    
cdef double eval_hot_BB_burst_Q(size_t THREAD,
                     double E,
                     double mu,
                     const double *const VEC,
                     void *const data) noexcept nogil:
    """
    Computes stokes Q. Missing: why are these numbers correct?
    
    """

    cdef double I_E
    I_E = eval_hot_BB_burst_I(THREAD,E,mu,VEC,data)
    cdef double PD = 0.1171*(mu - 1.)/(1. + 3.582*mu)
    return PD*I_E

    
    


