#cython: cdivision=True
#cython: boundscheck=False
#cython: nonecheck=False
#cython: wraparound=False

"""
This module implements 4 dimensional cubic polynomial Lagrangian interpolation.
Cubic means there are 4 data points in each dimension for which a polynomial 
will be built (https://en.wikipedia.org/wiki/Lagrange_polynomial). We have been
using this specifically for RMPs, which have a 4D atmosphere table.

Functions:
    - init_hot_Num4D: initialize the memory-space that will be used by the 
    variables in the interpolation.
    - free_hot_Num4D: free the memory-space
    - eval_hot_Num4D: computes the 4D interpolatin
    - eval_hot_Num4D_I: wrapper for eval_hot
    - eval_hot_Num4D_Q: wrapper for eval_hot
    - eval_hot_norm_Num4D: any multiplicative normalisation of the output of 
    eval_hot.
    
"""

from libc.math cimport M_PI, sqrt, sin, cos, acos, log10, pow, exp, fabs
from libc.stdio cimport printf, fopen, fclose, fread, FILE
from GSL cimport gsl_isnan, gsl_isinf
from libc.stdlib cimport malloc, free

from xpsi.global_imports import _keV, _k_B, _h_keV

cdef int SUCCESS = 0
cdef int ERROR = 1

cdef double erg = 1.0e-7
cdef double k_B = _k_B
cdef double keV = _keV
cdef double h_keV = _h_keV
cdef double k_B_over_keV = k_B / keV
cdef int VERBOSE = 0

ctypedef struct ACCELERATE:
    size_t **BN		# BN is a pointer to a pointer to a size_t object. 
    # So that allows to create a dynamic 2D array. The array BN will point to 
    # data points of a 4*n interpolation hypercube, where n is the amount of 
    # parameters in the interpolation. BN is later used as an index for node 
    # vals array.
    double **node_vals		# For each BN value, there are two double objects 
    # node_vals. The second one is the next value in the contiguous atmosphere 
    # array.
    double **SPACE		# SPACE are the denominators of the Lagrange 
    # polynomials.
    double **DIFF		# DIFF are the numerators of the Lagrange polynomials.
    double **INTENSITY_CACHE	# Caches intensity values for reuse.
    double **VEC_CACHE		# Cache VEC values for reuse.

# Modify this struct if useful for the user-defined source radiation field.
# Note that the members of DATA will be shared by all threads and are
# statically allocated, whereas the members of ACCELERATE will point to
# dynamically allocated memory, not shared by threads.

ctypedef struct DATA:
    const _preloaded *p 	# _preloaded is a struct defined in preload.pxd
    ACCELERATE acc		# Make ACCELERATE object


cdef void* init_hot_Num4D(size_t numThreads, 
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

    cdef DATA *D = <DATA*> malloc(sizeof(DATA))	 # Define DATA object
    D.p = preloaded  # Store preloaded information from function call in DATA 
    # object. See also preload.pyx.
    
    # These BLOCKS are related to the number of interpolation points needed in 
    # a hypercube.
    D.p.BLOCKS[0] = 64
    D.p.BLOCKS[1] = 16
    D.p.BLOCKS[2] = 4

    cdef size_t T, i, j, k, l

    # Prepare memory slots for ACCELERATE object
    D.acc.BN = <size_t**> malloc(numThreads * sizeof(size_t*))
    D.acc.node_vals = <double**> malloc(numThreads * sizeof(double*))
    D.acc.SPACE = <double**> malloc(numThreads * sizeof(double*))
    D.acc.DIFF = <double**> malloc(numThreads * sizeof(double*))
    D.acc.INTENSITY_CACHE = <double**> malloc(numThreads * sizeof(double*))
    D.acc.VEC_CACHE = <double**> malloc(numThreads * sizeof(double*))

    # Allocate memory for each thread, initializing data structures used for 
    # interpolation and caching
    for T in range(numThreads):
        # Allocate memory for base nodes (BN), node values, space differences 
        # (SPACE), and interpolation differences (DIFF)
        D.acc.BN[T] = <size_t*> malloc(D.p.ndims * sizeof(size_t))
        D.acc.node_vals[T] = <double*> malloc(2 * D.p.ndims * sizeof(double))
        D.acc.SPACE[T] = <double*> malloc(4 * D.p.ndims * sizeof(double))  
        # Cubic interpolation: 4 nodes per dimension
        D.acc.DIFF[T] = <double*> malloc(4 * D.p.ndims * sizeof(double))   
        # Cubic interpolation: 4 difference values per dimension
        D.acc.INTENSITY_CACHE[T] = <double*> malloc(256 * sizeof(double)) 
        # Full hypercube: 4^5 = 1024 intensities
        D.acc.VEC_CACHE[T] = <double*> malloc(D.p.ndims * sizeof(double))
        
        # Initialize base nodes and cache for each dimension
        for i in range(D.p.ndims):
             D.acc.BN[T][i] = 0  # Set initial base nodes to zero
             D.acc.VEC_CACHE[T][i] = D.p.params[i][1]  # Initialize with 
             # atmosphere value
             D.acc.node_vals[T][2 * i] = D.p.params[i][1]   # Store current and
             # next atmosphere values as node values
             D.acc.node_vals[T][2 * i + 1] = D.p.params[i][2]
     
             # Calculate the reciprocal of the difference for Lagrangian 
             # interpolation
             j = 4 * i
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
    # For the full interpolation hypercube, store all intensities in an array
    # with the right shape, so all values are lookupable later by knowing the 
    # i,j,k,l address.
    for T in range(numThreads):
        for i in range(4):
            for j in range(4):
                for k in range(4):
                    for l in range(4):
                        address = D.p.intensity + (D.acc.BN[T][0] + i) * D.p.S[0]
                        address += (D.acc.BN[T][1] + j) * D.p.S[1]
                        address += (D.acc.BN[T][2] + k) * D.p.S[2]
                        address += D.acc.BN[T][3] + l
                        D.acc.INTENSITY_CACHE[T][i * D.p.BLOCKS[0] + j * D.p.BLOCKS[1] + k * D.p.BLOCKS[2] + l] = address[0]

    # Cast for generalised usage in integration routines
    return <void*> D


cdef int free_hot_Num4D(size_t numThreads, void *const data) noexcept nogil:
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


cdef double eval_hot_Num4D(size_t THREAD,
                     double E,
                     double mu,
                     const double *const VEC,
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
       E (double): Photon energy in units provided by the integrator. For RMPs
       it is keV. Some conversion is done to match the units in the atmosphere
       table.
       mu (double): Cosine of the ray zenith angle (angle to surface normal).
       VEC (const double *const): Pointer to variables (e.g., temperature in 
        log10 of K, effective gravity log10 of cm/s^2).
       data (void *const): Numerical model data required for intensity 
       evaluation.
           The function must appropriately cast the void pointer for use.
           
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
        double vec[4]
        double E_eff = k_B_over_keV * pow(10.0, VEC[0])
        int update_baseNode[4]
        int CACHE = 0

    vec[0] = VEC[0]
    vec[1] = VEC[1]
    vec[2] = mu
    vec[3] = log10(E / E_eff)

    # Loop through each dimension and update the base node if necessary
    while i < D.p.ndims:
        update_baseNode[i] = 0  # Initialize: no change to base node
    
        # Check if the input value (vec) is smaller than the base node value 
        # (node_vals[2*i]) and the current base node (BN) is not the first (BN 
        # is set to zero in init_hot)
        if vec[i] < node_vals[2 * i] and BN[i] != 0:
            update_baseNode[i] = 1  # Mark the base node for change
    
            # Adjust the base node while the input value remains below the next 
            # base node value
            while vec[i] < D.p.params[i][BN[i] + 1]:
                if BN[i] > 0:  # If not at the first base node, decrement
                    BN[i] -= 1
                elif vec[i] <= D.p.params[i][0]:  # If vec is smaller than the 
                # first parameter value
                    vec[i] = D.p.params[i][0]  # Clamp to the first value
                    break
                elif BN[i] == 0:  # If already at the first base node, no 
                # further adjustment needed
                    break

            # Update node values based on the updated base node (BN). The +1 
            # and +2 ensure that the node values correctly surround the input 
            # value.
            node_vals[2 * i] = D.p.params[i][BN[i] + 1]     
            node_vals[2 * i + 1] = D.p.params[i][BN[i] + 2] 
        
        # Check if the input value (vec) exceeds the upper node value 
        # (node_vals[2 * i + 1]) and the current base node (BN) is not the last
        # one that allows for an update.
        elif vec[i] > node_vals[2 * i + 1] and BN[i] != D.p.N[i] - 4:
            update_baseNode[i] = 1  # Mark the base node for change
        
            # Move the base node up while the input value remains above the 
            # next base node value.
            while vec[i] > D.p.params[i][BN[i] + 2]:
                if BN[i] < D.p.N[i] - 4:  # If not at the last allowed base 
                # node, increment BN.
                    BN[i] += 1
                elif vec[i] >= D.p.params[i][D.p.N[i] - 1]:  # If vec is at the
                # upper limit, clamp it.
                    vec[i] = D.p.params[i][D.p.N[i] - 1]
                    break
                elif BN[i] == D.p.N[i] - 4:  # If at the last base node, stop
                # adjusting.
                    break
        
            # Update node values based on the new base node position.
            node_vals[2 * i] = D.p.params[i][BN[i] + 1]
            node_vals[2 * i + 1] = D.p.params[i][BN[i] + 2]
            
        # Here we avoid extra work with the caching system. We only compute 
        # this if the query value is not the same as the cached value, or more
        # obviously if the base node was changed.
        if V_CACHE[i] != vec[i] or update_baseNode[i] == 1:
            ii = 4*i
            # Go through the work of fetching the numerators. You will need the
            # query value for that.
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


            V_CACHE[i] = vec[i]	 # Store this input value in the cache for next
            # time so that work can be skipped.


        # For the denominators you have to redo the work only if the basenode 
        # was changed. vec[i] is not present in the denominators.
        if update_baseNode[i] == 1:	
            CACHE = 1  # If the basenode was changed, this is a cache flag 
            # indicating that it is the case, to be used later.
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

        i += 1	# For each dimension (while loop)
    cdef size_t j, k, l, INDEX, II, JJ, KK
    cdef double *address = NULL
    
    # Combinatorics over nodes of hypercube; weight cgs intensities
    for i in range(4):
        II = i * D.p.BLOCKS[0]
        for j in range(4):
            JJ = j * D.p.BLOCKS[1]
            for k in range(4):
                KK = k * D.p.BLOCKS[2]
                for l in range(4):
                    address = D.p.intensity + (BN[0] + i) * D.p.S[0]
                    address += (BN[1] + j) * D.p.S[1]
                    address += (BN[2] + k) * D.p.S[2]
                    address += BN[3] + l

                    temp = DIFF[i] * DIFF[4 + j] * DIFF[8 + k] * DIFF[12 + l]
                    temp *= (SPACE[i] * SPACE[4 + j] * SPACE[8 + k] * 
                             SPACE[12 + l])
                    INDEX = II + JJ + KK + l
                    if CACHE == 1:
                        I_CACHE[INDEX] = address[0]
                    I += temp * I_CACHE[INDEX]

    if I < 0.0:
        return 0.0
    return I * pow(10.0, 3.0 * vec[0]) # back to units required by integrator, 
    # except normalisation

cdef double eval_hot_norm_Num4D() noexcept nogil:
    """
    Calculate the normalization factor for the source radiation field. The
    normalisation step could be different depending on what is in the
    atmosphere data.
    
    Purpose:
        This function computes a normalization factor for the specific 
        intensity of the source radiation field.
    
    Returns:
        double: Normalization factor for the specific intensity in units of 
                J/cm²/s/keV/steradian.
    
    Attributes:
        noexcept nogil: Indicates that the function does not throw exceptions 
                       and can be executed without the Global Interpreter Lock
                       (GIL).
    """
    return erg / h_keV

