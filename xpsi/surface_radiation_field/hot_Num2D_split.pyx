#cython: cdivision=True
#cython: boundscheck=False
#cython: nonecheck=False
#cython: wraparound=False

"""
This module does 2D interpolation. It also accompanies hot_Num5D_split, 
intended for interpolation for AMXPs, that's why it's called 'split' at the 
moment. Please look at hot_num5D_split or hot_num4D for detailed comments.


"""

from libc.math cimport M_PI, sqrt, sin, cos, acos, log10, pow, exp, fabs
from libc.stdio cimport printf, fopen, fclose, fread, FILE
from libc.stdlib cimport malloc, free

from xpsi.global_imports import _keV, _k_B

cdef int SUCCESS = 0
cdef int ERROR = 1

cdef double erg = 1.0e-7
cdef double k_B = _k_B
cdef double keV = _keV
cdef double k_B_over_keV = k_B / keV
cdef int VERBOSE = 0

ctypedef struct ACCELERATE:
    size_t **BN
    double **node_vals	
    double **SPACE
    double **DIFF
    double **INTENSITY_CACHE
    double **VEC_CACHE

ctypedef struct DATA:
    const _preloaded *p 	
    ACCELERATE acc		

cdef void* init_hot_2D(size_t numThreads, 
                       const _preloaded *const preloaded) noexcept nogil:
    """
    Initialize data for hot region intensity evaluation in a multi-threaded 
    environment.
    
    User-Modifiable Function:
        The user is free to wrap thread-safe and non-parallel external C 
        routines from an external library. The function body does not 
        necessarily need to be written in Cython, allowing flexibility 
        in integrating external code.
    
    Memory Management:
        This function must match the free management routine `free_hot()` to 
        properly free dynamically allocated memory. Managing memory correctly 
        is the user's responsibility. The function should return NULL if 
        dynamic memory is not required for the model.
    
    Arguments:
        numThreads (size_t): Number of threads for parallel execution.
        preloaded (const _preloaded *const): Preloaded data structure 
        containing model parameters.
    
    Returns:
        void*: Pointer to the dynamically allocated DATA structure.
    
    Attributes:
        noexcept nogil: Indicates that the function does not throw exceptions 
                       and can be executed without the Global Interpreter 
                       Lock (GIL).
    """
    
    if preloaded == NULL :
        printf("ERROR: The numerical atmosphere data were not preloaded, \
               which are required by this extension.\n")

    
    cdef DATA *D = <DATA*> malloc(sizeof(DATA))
    D.p = preloaded

    # Notice here fewer blocks
    D.p.BLOCKS[0] = 4    


    cdef size_t T, i, j

    D.acc.BN = <size_t**> malloc(numThreads * sizeof(size_t*))
    D.acc.node_vals = <double**> malloc(numThreads * sizeof(double*))
    D.acc.SPACE = <double**> malloc(numThreads * sizeof(double*))
    D.acc.DIFF = <double**> malloc(numThreads * sizeof(double*))
    D.acc.INTENSITY_CACHE = <double**> malloc(numThreads * sizeof(double*))
    D.acc.VEC_CACHE = <double**> malloc(numThreads * sizeof(double*))

    for T in range(numThreads):						
        D.acc.BN[T] = <size_t*> malloc(D.p.ndims * sizeof(size_t)) 		
        D.acc.node_vals[T] = <double*> malloc(2 * D.p.ndims * sizeof(double))	
        D.acc.SPACE[T] = <double*> malloc(4 * D.p.ndims * sizeof(double)) 	
        D.acc.DIFF[T] = <double*> malloc(4 * D.p.ndims * sizeof(double))	
        # Notice here 16 points for the hypercube with 2D interpolation.
        D.acc.INTENSITY_CACHE[T] = <double*> malloc(16 * sizeof(double)) 
        D.acc.VEC_CACHE[T] = <double*> malloc(D.p.ndims * sizeof(double))	

        for i in range(D.p.ndims):						
            D.acc.BN[T][i] = 0							
            D.acc.VEC_CACHE[T][i] = D.p.params[i][1]			
            D.acc.node_vals[T][2*i] = D.p.params[i][1]			
            D.acc.node_vals[T][2*i + 1] = D.p.params[i][2]

            j = 4*i

            D.acc.SPACE[T][j] = 1.0 / (D.p.params[i][0] - D.p.params[i][1])
            D.acc.SPACE[T][j] /= D.p.params[i][0] - D.p.params[i][2]
            D.acc.SPACE[T][j] /= D.p.params[i][0] - D.p.params[i][3]

            D.acc.SPACE[T][j + 1] = 1.0 / (D.p.params[i][1] - D.p.params[i][0])
            D.acc.SPACE[T][j + 1] /= D.p.params[i][1] - D.p.params[i][2]
            D.acc.SPACE[T][j + 1] /= D.p.params[i][1] - D.p.params[i][3]

            D.acc.SPACE[T][j + 2] = 1.0 / (D.p.params[i][2] - D.p.params[i][0])
            D.acc.SPACE[T][j + 2] /= D.p.params[i][2] - D.p.params[i][1]
            D.acc.SPACE[T][j + 2] /= D.p.params[i][2] - D.p.params[i][3]

            D.acc.SPACE[T][j + 3] = 1.0 / (D.p.params[i][3] - D.p.params[i][0])
            D.acc.SPACE[T][j + 3] /= D.p.params[i][3] - D.p.params[i][1]
            D.acc.SPACE[T][j + 3] /= D.p.params[i][3] - D.p.params[i][2]

            D.acc.DIFF[T][j] = D.acc.VEC_CACHE[T][i] - D.p.params[i][1]
            D.acc.DIFF[T][j] *= D.acc.VEC_CACHE[T][i] - D.p.params[i][2]
            D.acc.DIFF[T][j] *= D.acc.VEC_CACHE[T][i] - D.p.params[i][3]

            D.acc.DIFF[T][j + 1] = D.acc.VEC_CACHE[T][i] - D.p.params[i][0]
            D.acc.DIFF[T][j + 1] *= D.acc.VEC_CACHE[T][i] - D.p.params[i][2]
            D.acc.DIFF[T][j + 1] *= D.acc.VEC_CACHE[T][i] - D.p.params[i][3]

            D.acc.DIFF[T][j + 2] = D.acc.VEC_CACHE[T][i] - D.p.params[i][0]
            D.acc.DIFF[T][j + 2] *= D.acc.VEC_CACHE[T][i] - D.p.params[i][1]
            D.acc.DIFF[T][j + 2] *= D.acc.VEC_CACHE[T][i] - D.p.params[i][3]

            D.acc.DIFF[T][j + 3] = D.acc.VEC_CACHE[T][i] - D.p.params[i][0]
            D.acc.DIFF[T][j + 3] *= D.acc.VEC_CACHE[T][i] - D.p.params[i][1]
            D.acc.DIFF[T][j + 3] *= D.acc.VEC_CACHE[T][i] - D.p.params[i][2]

    cdef double *address = NULL

    # Notice here the nested for loops have only 2 layers for 2D interpolation.
    for T in range(numThreads): 
        for i in range(4):
            for j in range(4):
                address = D.p.intensity + (D.acc.BN[T][0] + i) * D.p.S[0]
                address += D.acc.BN[T][1] + j
                D.acc.INTENSITY_CACHE[T][i * D.p.BLOCKS[0] + j] = address[0]
    return <void*> D


cdef int free_hot_2D(size_t numThreads, void *const data) noexcept nogil:
    """
    Free dynamically allocated memory associated with hot region data.
    
    This function must match the initialization routine `init_hot()` in terms
    of freeing dynamically allocated memory. It is entirely the user's 
    responsibility to manage memory correctly.
    
    Arguments:
        numThreads (size_t): Number of threads used for parallel execution.
        data (void *const): Pointer to dynamically allocated memory. 
            - The void pointer must be appropriately cast before freeing the 
            memory.
            - Use `free(<void*> data)` only if no memory was dynamically 
            allocated in the function `init_hot()`, as `data` is expected to be
            NULL in this case.
    
    Attributes:
        noexcept nogil: Indicates that the function does not throw exceptions 
                       and can be executed without the Global Interpreter Lock
                       (GIL).
    """

    cdef DATA *D = <DATA*> data
    cdef size_t T

    for T in range(numThreads):
        free(D.acc.BN[T])
        free(D.acc.node_vals[T])
        free(D.acc.SPACE[T])
        free(D.acc.DIFF[T])
        free(D.acc.INTENSITY_CACHE[T])
        free(D.acc.VEC_CACHE[T])

    free(D.acc.BN)
    free(D.acc.node_vals)
    free(D.acc.SPACE)
    free(D.acc.DIFF)
    free(D.acc.INTENSITY_CACHE)
    free(D.acc.VEC_CACHE)
    free(D)

    return SUCCESS


cdef double eval_hot_2D(size_t THREAD,
                     double E,
                     double mu,
                     void *const data) noexcept nogil:
    
    """
    Evaluate the intensity of hot regions based on given parameters.
    
    Cubic polynomial interpolation:
    This function implements cubic polynomial interpolation to improve 
    acceleration properties. Specifically, it avoids recomputing numerical 
    weights or re-reading intensities when not necessary, optimizing 
    performance.
    
    Arguments:
        THREAD (size_t): Thread ID used for parallel execution.
        E (double): Photon energy in units provided by the integrator. 
            For AMXPs, the integrator provides electron rest energy, which is
            the same as in the atmosphere table.
        mu (double): Cosine of the ray zenith angle (angle to surface normal).
            The function must appropriately cast the void pointer for use.
        data (void *const): Numerical model data required for intensity 
            evaluation.    
        
    Attributes:
        noexcept nogil: Indicates that the function does not throw exceptions 
                       and can be executed without the Global Interpreter Lock
                       (GIL).
    """

    cdef DATA *D = <DATA*> data

    cdef: 
        size_t i = 0, ii
        double I = 0.0, temp
        double *node_vals = D.acc.node_vals[THREAD] 
        size_t *BN = D.acc.BN[THREAD]
        double *SPACE = D.acc.SPACE[THREAD]
        double *DIFF = D.acc.DIFF[THREAD]
        double *I_CACHE = D.acc.INTENSITY_CACHE[THREAD]
        double *V_CACHE = D.acc.VEC_CACHE[THREAD]
        double vec[2] 
        int update_baseNode[2]  
        int CACHE = 0

    vec[0] = mu
    vec[1] = E

    while i < D.p.ndims: 					
        update_baseNode[i] = 0					
        if vec[i] < node_vals[2*i] and BN[i] != 0:		
            update_baseNode[i] = 1				
            while vec[i] < D.p.params[i][BN[i] + 1]:	
                if BN[i] > 0:					
                    BN[i] -= 1				
                elif vec[i] <= D.p.params[i][0]:		
                    vec[i] = D.p.params[i][0]			
                    break
                elif BN[i] == 0:				
                    break

            node_vals[2*i] = D.p.params[i][BN[i] + 1] 		
            node_vals[2*i + 1] = D.p.params[i][BN[i] + 2]

        elif vec[i] > node_vals[2*i + 1] and BN[i] != D.p.N[i] - 4: 	
            update_baseNode[i] = 1 				
            while vec[i] > D.p.params[i][BN[i] + 2]:	
                if BN[i] < D.p.N[i] - 4:			
                    BN[i] += 1					
                elif vec[i] >= D.p.params[i][D.p.N[i] - 1]:
                    vec[i] = D.p.params[i][D.p.N[i] - 1]	
                    break
                elif BN[i] == D.p.N[i] - 4:			
                    break

            node_vals[2*i] = D.p.params[i][BN[i] + 1]	
            node_vals[2*i + 1] = D.p.params[i][BN[i] + 2]

        if V_CACHE[i] != vec[i] or update_baseNode[i] == 1:
            ii = 4*i
            DIFF[ii] = vec[i] - D.p.params[i][BN[i] + 1]	
            DIFF[ii] *= vec[i] - D.p.params[i][BN[i] + 2]
            DIFF[ii] *= vec[i] - D.p.params[i][BN[i] + 3]

            DIFF[ii + 1] = vec[i] - D.p.params[i][BN[i]]
            DIFF[ii + 1] *= vec[i] - D.p.params[i][BN[i] + 2]
            DIFF[ii + 1] *= vec[i] - D.p.params[i][BN[i] + 3]

            DIFF[ii + 2] = vec[i] - D.p.params[i][BN[i]]
            DIFF[ii + 2] *= vec[i] - D.p.params[i][BN[i] + 1]
            DIFF[ii + 2] *= vec[i] - D.p.params[i][BN[i] + 3]

            DIFF[ii + 3] = vec[i] - D.p.params[i][BN[i]]
            DIFF[ii + 3] *= vec[i] - D.p.params[i][BN[i] + 1]
            DIFF[ii + 3] *= vec[i] - D.p.params[i][BN[i] + 2]

            V_CACHE[i] = vec[i]				

        if update_baseNode[i] == 1:			

            CACHE = 1						
            SPACE[ii] = 1.0 / (D.p.params[i][BN[i]] - D.p.params[i][BN[i] + 1])
            SPACE[ii] /= D.p.params[i][BN[i]] - D.p.params[i][BN[i] + 2]
            SPACE[ii] /= D.p.params[i][BN[i]] - D.p.params[i][BN[i] + 3]

            SPACE[ii + 1] = 1.0 / (D.p.params[i][BN[i] + 1] - 
                                   D.p.params[i][BN[i]])
            SPACE[ii + 1] /= (D.p.params[i][BN[i] + 1] - 
                              D.p.params[i][BN[i] + 2])
            SPACE[ii + 1] /= (D.p.params[i][BN[i] + 1] - 
                              D.p.params[i][BN[i] + 3])

            SPACE[ii + 2] = 1.0 / (D.p.params[i][BN[i] + 2] - 
                                   D.p.params[i][BN[i]])
            SPACE[ii + 2] /= (D.p.params[i][BN[i] + 2] - 
                              D.p.params[i][BN[i] + 1])
            SPACE[ii + 2] /= (D.p.params[i][BN[i] + 2] - 
                              D.p.params[i][BN[i] + 3])

            SPACE[ii + 3] = 1.0 / (D.p.params[i][BN[i] + 3] - 
                                   D.p.params[i][BN[i]])
            SPACE[ii + 3] /= (D.p.params[i][BN[i] + 3] - 
                              D.p.params[i][BN[i] + 1])
            SPACE[ii + 3] /= (D.p.params[i][BN[i] + 3] - 
                              D.p.params[i][BN[i] + 2])

        i += 1							

    cdef size_t j, INDEX, II
    cdef double *address = NULL
        
    
    for i in range(4):
        II = i * D.p.BLOCKS[0]
        for j in range(4):
            address = D.p.intensity + (BN[0] + i) * D.p.S[0]
            address += BN[1] + j 			

            temp = DIFF[i] * DIFF[4 + j] 
            temp *= SPACE[i] * SPACE[4 + j] 
            INDEX = II + j 		
            
            if CACHE == 1:				
                I_CACHE[INDEX] = address[0]	

            I += temp * I_CACHE[INDEX]

    return I

cdef double eval_hot_2D_I(size_t THREAD,
                     double E,
                     double mu,
                     void *const data) noexcept nogil:
    """
    Evaluate the intensity for hot regions using given parameters.
    
    This function calculates the intensity based on the provided photon energy, 
    zenith angle, and variable parameters. It wraps the core intensity 
    evaluation function `eval_hot` and ensures that the returned intensity is 
    non-negative.
    
    Arguments:
        THREAD (size_t): Thread ID used for parallel execution.
        E (double): Photon energy in electron rest energy for AMXPs.
        mu (double): Cosine of the ray zenith angle (angle to surface normal).
        data (void *const): Numerical model data required for intensity 
        evaluation.
    
    Returns:
        double: Calculated intensity value, or 0.0 if the evaluated intensity 
        is negative.
    
    Attributes:
        noexcept nogil: Indicates that the function does not throw exceptions 
                       and can be executed without the Global Interpreter Lock 
                       (GIL).
    """

    cdef double I = eval_hot_2D(THREAD,E,mu,data)

    if I < 0.0:
        return 0.0

    return I


cdef double eval_hot_2D_Q(size_t THREAD,
                     double E,
                     double mu,
                     void *const data) noexcept nogil:
    """
    Same as eval_hot_2D_I but returns the stokes Q, which doesn't have to
    be above 0. 
    """

    return eval_hot_2D(THREAD,E,mu,data)

cdef double eval_hot_2D_norm() noexcept nogil:
    """
    Same multiplicative normalisation as in the hot_Num5D_split interpolator.
    """
    return erg / 4.135667662e-18

