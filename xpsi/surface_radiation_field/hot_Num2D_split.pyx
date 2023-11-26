#cython: cdivision=True
#cython: boundscheck=False
#cython: nonecheck=False
#cython: wraparound=False

from libc.math cimport M_PI, sqrt, sin, cos, acos, log10, pow, exp, fabs
from libc.stdio cimport printf, fopen, fclose, fread, FILE
#from GSL cimport gsl_isnan, gsl_isinf
from libc.stdlib cimport malloc, free

from xpsi.global_imports import _keV, _k_B

cdef int SUCCESS = 0
cdef int ERROR = 1

cdef double erg = 1.0e-7
cdef double k_B = _k_B
cdef double keV = _keV
cdef double k_B_over_keV = k_B / keV
cdef int VERBOSE = 0

# We are doing cubic polynomial Lagrangian interpolation. There are 4 data points in each dimension for which a polynomial will be built. https://en.wikipedia.org/wiki/Lagrange_polynomial

ctypedef struct ACCELERATE:
    size_t **BN		# BN is a pointer to a pointer to a size_t object. So that allows to create a dynamic 2D array. The array BN will point to data points of a 4*n interpolation hypercube, where n is the amount of parameters in the interpolation. BN is later used as an index for node vals array.
    double **node_vals		# For each BN value, there are two double objects node_vals. The second one is the next value in the contiguous atmosphere array.
    double **SPACE		# SPACE are the denominators of the Lagrange polynomials.
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

#----------------------------------------------------------------------->>>
# >>> User modifiable functions.
# >>> Note that the user is entirely free to wrap thread-safe and
# ... non-parallel external C routines from an external library.
# >>> Thus the bodies of the following need not be written explicitly in
# ... the Cython language.
#----------------------------------------------------------------------->>>
cdef void* init_hot_2D(size_t numThreads, const _preloaded *const preloaded) nogil:
    # This function must match the free management routine free_hot()
    # in terms of freeing dynamically allocated memory. This is entirely
    # the user's responsibility to manage.
    # Return NULL if dynamic memory is not required for the model

    #printf("inside init_hot()")
    cdef DATA *D = <DATA*> malloc(sizeof(DATA))	# Define DATA object
    D.p = preloaded 					# Store preloaded information from function call in DATA object. See also preload.pyx.
    
    # (1) These BLOCKS appear to be related to the number of interpolation
    # points needed in a "hypercube". However, I would expect this to be 256 
    # for 4 dimensional interpolation already..

    # D.p.BLOCKS[0] = 64
    # D.p.BLOCKS[1] = 16
    # D.p.BLOCKS[2] = 4

    # By analogy, remove two factors of four.

    D.p.BLOCKS[0] = 4    


    cdef size_t T, i, j

    # Prepare memory slots for ACCELERATE object
    D.acc.BN = <size_t**> malloc(numThreads * sizeof(size_t*))
    D.acc.node_vals = <double**> malloc(numThreads * sizeof(double*))
    D.acc.SPACE = <double**> malloc(numThreads * sizeof(double*))
    D.acc.DIFF = <double**> malloc(numThreads * sizeof(double*))
    D.acc.INTENSITY_CACHE = <double**> malloc(numThreads * sizeof(double*))
    D.acc.VEC_CACHE = <double**> malloc(numThreads * sizeof(double*))

    for T in range(numThreads):						# Same for each thread
        D.acc.BN[T] = <size_t*> malloc(D.p.ndims * sizeof(size_t)) 		# See ACCELERATE definition above.
        D.acc.node_vals[T] = <double*> malloc(2 * D.p.ndims * sizeof(double))	# See ACCELERATE definition above.
        D.acc.SPACE[T] = <double*> malloc(4 * D.p.ndims * sizeof(double)) 	# k=3 (cubic), so we have k+1=4 nodes per dimension, and 4 SPACEs
        D.acc.DIFF[T] = <double*> malloc(4 * D.p.ndims * sizeof(double))	# k=3 (cubic), so we have k+1=4 nodes per dimension, and 4 DIFFs
        #D.acc.INTENSITY_CACHE[T] = <double*> malloc(256 * sizeof(double)) 	# (2) N.B. This has been changed compared to 4D
        D.acc.INTENSITY_CACHE[T] = <double*> malloc(16 * sizeof(double))	# There are 4^2 (nodes^dimensions) intensities (the full hypercube)  
        D.acc.VEC_CACHE[T] = <double*> malloc(D.p.ndims * sizeof(double))	# See ACCELERATE definition above.

        for i in range(D.p.ndims):						# For each dimension
            D.acc.BN[T][i] = 0							# Initially store zeros in the basenodes of the hypercube
            D.acc.VEC_CACHE[T][i] = D.p.params[i][1]				# Store atmosphere value in the VEC_CACHE
            D.acc.node_vals[T][2*i] = D.p.params[i][1]				# Store atmosphere value and the next one as the current node_vals
            D.acc.node_vals[T][2*i + 1] = D.p.params[i][2]

            j = 4*i								# So we can do j, j+1, j+2, j+3

            D.acc.SPACE[T][j] = 1.0 / (D.p.params[i][0] - D.p.params[i][1])	# See also definition of Lagrangian interpolation.
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
    
        # printf("diagnostics for initialization\n")
        # for i in range(D.p.ndims):
        #     printf("i=%d, ", <int>i)
        #     printf("D.p.params[i][0]: %.2e\n", D.p.params[i][0])
        #     printf("VEC_CACHE[T][i]: %.2e\n", D.acc.VEC_CACHE[T][i])
        

    cdef double *address = NULL
    # Here we produce the intensity cache.
    # printf("\ncommencing cache intensity")
    
    # (3) For every dimension, we have a D.acc.BN[T][i], so I add a for-loop to 
    # this. It pains me to have so many forloops.
    
    # for T in range(numThreads):
    #     for i in range(4):
    #         for j in range(4):
    #             for k in range(4):
    #                 for l in range(4):
    #                     address = D.p.I + (D.acc.BN[T][0] + i) * D.p.S[0]
    #                     address += (D.acc.BN[T][1] + j) * D.p.S[1]
    #                     address += (D.acc.BN[T][2] + k) * D.p.S[2]
    #                     address += D.acc.BN[T][3] + l
    #                     D.acc.INTENSITY_CACHE[T][i * D.p.BLOCKS[0] + j * D.p.BLOCKS[1] + k * D.p.BLOCKS[2] + l] = address[0]

    for T in range(numThreads): #For the full interpolation hypercube, store all intensities in an array with the right shape, so all values are lookupable later by knowing the i,j,k,l,m address.
        for i in range(4):
            for j in range(4):
                address = D.p.I + (D.acc.BN[T][0] + i) * D.p.S[0] 
                address += D.acc.BN[T][1] + j
                D.acc.INTENSITY_CACHE[T][i * D.p.BLOCKS[0] + j] = address[0]


    # Cast for generalised usage in integration routines
    return <void*> D


cdef int free_hot_2D(size_t numThreads, void *const data) nogil:
    # This function must match the initialisation routine init_hot()
    # in terms of freeing dynamically allocated memory. This is entirely
    # the user's responsibility to manage.
    # The void pointer must be appropriately cast before memory is freed --
    # only the user can know this at compile time.
    # Just use free(<void*> data) iff no memory was dynamically
    # allocated in the function:
    #   init_hot()
    # because data is expected to be NULL in this case

    # printf("inside free_hot()")

    cdef DATA *D = <DATA*> data

    cdef size_t T

    for T in range(numThreads):
        # printf("freeing thread specific memory")
        free(D.acc.BN[T])
        free(D.acc.node_vals[T])
        free(D.acc.SPACE[T])
        free(D.acc.DIFF[T])
        free(D.acc.INTENSITY_CACHE[T])
        free(D.acc.VEC_CACHE[T])

    # printf("freeing D.acc...")
    free(D.acc.BN)
    free(D.acc.node_vals)
    free(D.acc.SPACE)
    free(D.acc.DIFF)
    free(D.acc.INTENSITY_CACHE)
    free(D.acc.VEC_CACHE)

    # printf("freeing D...")
    free(D)

    return SUCCESS

#----------------------------------------------------------------------->>>
# >>> Cubic polynomial interpolation.
# >>> Improve acceleration properties... i.e. do not recompute numerical
# ... weights or re-read intensities if not necessary.
#----------------------------------------------------------------------->>>
    
cdef double eval_hot_2D(size_t THREAD,
                     double E,
                     double mu,
                     void *const data) nogil:
    
    # Arguments:
    # E = photon energy in keV
    # mu = cosine of ray zenith angle (i.e., angle to surface normal)
    # VEC = variables such as temperature, effective gravity, ...
    # data = numerical model data required for intensity evaluation
    # This function must cast the void pointer appropriately for use.
    cdef DATA *D = <DATA*> data

    cdef: 
        size_t i = 0, ii
        double I = 0.0, temp
        double *node_vals = D.acc.node_vals[THREAD] # unpacking information from D
        size_t *BN = D.acc.BN[THREAD]
        double *SPACE = D.acc.SPACE[THREAD]
        double *DIFF = D.acc.DIFF[THREAD]
        double *I_CACHE = D.acc.INTENSITY_CACHE[THREAD]
        double *V_CACHE = D.acc.VEC_CACHE[THREAD]
        double vec[2] # (4) should be = ndims
        # double E_eff = k_B_over_keV * pow(10.0, VEC[0])
        # double E_eff = k_B_over_keV * pow(10.0, Temperature)
        int update_baseNode[2]  # (5) should be = ndims
        int CACHE = 0

    # cdef double te, tbb, tau # I have three parameters in VEC
    # te = VEC[0]
    # tbb = VEC[1]
    # tau = VEC[2]
    
    cdef double evere = 0.5109989e6 # electron volts in elecron rest energy

    # The input value of the parameter (vec) to be interpolated. Note this is the order of *._hot_atmosphere
    # vec[0] = te
    # vec[1] = tbb
    # vec[2] = tau
    vec[0] = mu
    vec[1] = E#*1e3/evere # conversion from keV to electron rest energy
 
    # printf("Bobrikova atmosphere interpolator")
    
    # printf("diagnostics 0:\n")
    # printf("E: %.2e, ", E)
    # printf("Temperature: %.2e, ", Temperature)
    # printf("k_B_over_keV: %.2e, ", k_B_over_keV)
    # printf("E_eff: %.2e, ", E_eff)
    # printf("vec[4]: %.2e\n", vec[4])
    

    # printf("\nvec[0]: %.8e, ", vec[0])
    # printf("vec[1]: %.8e, ", vec[1])
    # printf("vec[2]: %.8e, ", vec[2])
    # printf("vec[3]: %.8e, ", vec[3])
    # printf("vec[4]: %.8e, ", vec[4])
    
    
    # printf("\neval_hot() called")
    # printf("\nvec[0]: %f", vec[0])
    # printf("\nvec[1]: %f", vec[1])
    # printf('ndims')
    # printf('\nD.p.ndims 2D: %ld', D.p.ndims)


    while i < D.p.ndims: 					# For each dimension 
        # if parallel == 31:
        # printf("\nDimension: %d", <int>i)
        update_baseNode[i] = 0					# Flag to change base node
        if vec[i] < node_vals[2*i] and BN[i] != 0:		# If the input value of the parameter (vec) is smaller than the base node value (we are valued below the hypercube/basenode interval), and BN is not the first (it is zet zero in init_hot)
            # if parallel == 31:
            # printf("\nExecute block 1: %d", <int>i)
            update_baseNode[i] = 1				# then the base node should be changed.
            while vec[i] < D.p.params[i][BN[i] + 1]:		# while the input value remains smaller than the value at the next basenode
                # if parallel == 31:
                #     printf("\n!")
                #     printf("\nvec i: %.8e", vec[i])
                #     printf("\nBase node: %d", <int>BN[i])
                if BN[i] > 0:					# and if the basenode is not the first
                    BN[i] -= 1					# decrement by 1
                elif vec[i] <= D.p.params[i][0]:		# or, if it IS not bigger than zero, and if vec is smaller than the first value
                    vec[i] = D.p.params[i][0]			# make it equal to the first value
                    break
                elif BN[i] == 0:				# or, if it IS not bigger than zero, and vec is bigger than the first value, then it should be zero.
                    break

            node_vals[2*i] = D.p.params[i][BN[i] + 1] 		# now update node values with updated basenode values. Foot note: Not sure why it is +1 and +2. It must be correctly surrounding the input value. Oh perhaps this is so that the first two are below and the second two are above?
            node_vals[2*i + 1] = D.p.params[i][BN[i] + 2]

            # if parallel == 31:
            # printf("\nEnd Block 1: %d", <int>i)

        elif vec[i] > node_vals[2*i + 1] and BN[i] != D.p.N[i] - 4: 	# Opposite case of above. If the input value is larger than the second node_vals and BN is not the last
            # if parallel == 31:
            # printf("\nExecute block 2: %d", <int>i)
            update_baseNode[i] = 1 				# The base node should be changed.
            while vec[i] > D.p.params[i][BN[i] + 2]:		# While the input value remains larger than the value at the third base node
                if BN[i] < D.p.N[i] - 4:			# if the base node is smaller than the last - 4
                    BN[i] += 1					# then it may be incremented.
                elif vec[i] >= D.p.params[i][D.p.N[i] - 1]:	# or, if the base node IS NOT smaller than the last - 4, and IF it is larger than one before the last value (foot note: i don't think it could be larger than the last value because the while statement already checks with 2 values above. So here we are at the limit, but not over it.)
                    vec[i] = D.p.params[i][D.p.N[i] - 1]	# then it should be set to that value.
                    break
                elif BN[i] == D.p.N[i] - 4:			# Or, if it IS NOT smaller than the last - 4 and vec is not too big, then it should be exactly the last base node.
                    break

            node_vals[2*i] = D.p.params[i][BN[i] + 1]		# This is exactly the same, update the node_vals from the base nodes.
            node_vals[2*i + 1] = D.p.params[i][BN[i] + 2]

            # if parallel == 31:
            # printf("\nEnd Block 2: %d", <int>i)

        # if parallel == 31:
        # printf("\nTry block 3: %d", <int>i)

        if V_CACHE[i] != vec[i] or update_baseNode[i] == 1:	# If the input value is not already the same as the cached value (so e.g. if the the change in the vec was super duper small and data is spread out), or more obviously if the base node was changed.
            # if parallel == 31:
            # printf("\nExecute block 3: %d", <int>i)
            ii = 4*i
            DIFF[ii] = vec[i] - D.p.params[i][BN[i] + 1]	# Go through the work of fetching the numerators. You will need the input value for that.
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

            # printf("\nupdating V_CACHE")


            V_CACHE[i] = vec[i]				# Store this input value in the cache for next time so that work can be skipped.

            # if parallel == 31:
            # printf("\nEnd block 3: %d", <int>i)

        # if parallel == 31:
        #     printf("\nTry block 4: %d", <int>i)

        if update_baseNode[i] == 1:				# For the denominators you have to redo the work if the basenode was changed.
            # if parallel == 31:
            # printf("\nExecute block 4: %d", <int>i)
            # printf("i=%d, ", <int>i)
            # printf("D.p.params[i][BN[i]]: %.2e\n", D.p.params[i][BN[i]])
            CACHE = 1						# If the basenode was changed, this is a cache flag indicating that we've done so.
            SPACE[ii] = 1.0 / (D.p.params[i][BN[i]] - D.p.params[i][BN[i] + 1])
            SPACE[ii] /= D.p.params[i][BN[i]] - D.p.params[i][BN[i] + 2]
            SPACE[ii] /= D.p.params[i][BN[i]] - D.p.params[i][BN[i] + 3]

            SPACE[ii + 1] = 1.0 / (D.p.params[i][BN[i] + 1] - D.p.params[i][BN[i]])
            SPACE[ii + 1] /= D.p.params[i][BN[i] + 1] - D.p.params[i][BN[i] + 2]
            SPACE[ii + 1] /= D.p.params[i][BN[i] + 1] - D.p.params[i][BN[i] + 3]

            SPACE[ii + 2] = 1.0 / (D.p.params[i][BN[i] + 2] - D.p.params[i][BN[i]])
            SPACE[ii + 2] /= D.p.params[i][BN[i] + 2] - D.p.params[i][BN[i] + 1]
            SPACE[ii + 2] /= D.p.params[i][BN[i] + 2] - D.p.params[i][BN[i] + 3]

            SPACE[ii + 3] = 1.0 / (D.p.params[i][BN[i] + 3] - D.p.params[i][BN[i]])
            SPACE[ii + 3] /= D.p.params[i][BN[i] + 3] - D.p.params[i][BN[i] + 1]
            SPACE[ii + 3] /= D.p.params[i][BN[i] + 3] - D.p.params[i][BN[i] + 2]

            # if parallel == 31:
            # printf("\nEnd block 4: %d", <int>i)

        # printf("\ncomputing DIFFs and SPACEs\n")
        # printf("DIFF[ii]: %.2e, ", DIFF[ii])
        # printf("DIFF[ii+1]: %.2e, ", DIFF[ii+1])
        # printf("DIFF[ii+2]: %.2e, ", DIFF[ii+2])
        # printf("DIFF[ii+3]: %.2e\n", DIFF[ii+3])
        
        
        # printf("SPACE[ii]: %.2e, ", SPACE[ii])
        # printf("SPACE[ii+1]: %.2e, ", SPACE[ii+1])
        # printf("SPACE[ii+2]: %.2e, ", SPACE[ii+2])
        # printf("SPACE[ii+3]: %.2e\n", SPACE[ii+3])

        i += 1							# For each dimension (while loop)

    # printf("Diagnostics: 2\n")
    # for i in range(D.p.ndims):
    #     printf("i=%d, ", <int>i)
    #     printf("vec[i]: %.2e, ", vec[i])
    #     printf("V_CACHE[i]: %.2e\n, ", V_CACHE[i])
    #     printf("D.p.params[i][BN[i]]: %.2e, ", D.p.params[i][BN[i]])
    #     printf("D.p.params[i][BN[i]+1]: %.2e, ", D.p.params[i][BN[i]+1])
    #     printf("D.p.params[i][BN[i]+2]: %.2e, ", D.p.params[i][BN[i]+2])
    #     printf("D.p.params[i][BN[i]+3]: %.2e\n", D.p.params[i][BN[i]+3])

    cdef size_t j, INDEX, II
    cdef double *address = NULL

    # (6) Here again, I need to iterate over an additional dimension.
    
    # Combinatorics over nodes of hypercube; weight cgs intensities
    # for i in range(4):
    #     II = i * D.p.BLOCKS[0]
    #     for j in range(4):
    #         JJ = j * D.p.BLOCKS[1]
    #         for k in range(4):
    #             KK = k * D.p.BLOCKS[2]
    #             for l in range(4):
    #                 address = D.p.I + (BN[0] + i) * D.p.S[0]
    #                 address += (BN[1] + j) * D.p.S[1]
    #                 address += (BN[2] + k) * D.p.S[2]
    #                 address += BN[3] + l

    #                 temp = DIFF[i] * DIFF[4 + j] * DIFF[8 + k] * DIFF[12 + l]
    #                 temp *= SPACE[i] * SPACE[4 + j] * SPACE[8 + k] * SPACE[12 + l]
    #                 INDEX = II + JJ + KK + l
    #                 if CACHE == 1:
    #                     I_CACHE[INDEX] = address[0]
    #                 I += temp * I_CACHE[INDEX]
    
    # printf("diagnostics for interpolation\n")
    # printf("D.p.S[0]: %d, ", <int>D.p.S[0])
    # printf("BN[0]: %d\n", <int>BN[0])
    # printf("D.p.S[1]: %d, ", <int>D.p.S[1])
    # printf("BN[1]: %d\n", <int>BN[1])
    # printf("D.p.S[2]: %d, ", <int>D.p.S[2])
    # printf("BN[2]: %d\n", <int>BN[2])
    # printf("D.p.S[3]: %d, ", <int>D.p.S[3])
    # printf("BN[3]: %d\n", <int>BN[3])
    # printf("BN[4]: %d\n", <int>BN[4])
        
    # Combinatorics over nodes of hypercube; weight cgs intensities
    for i in range(4):
        II = i * D.p.BLOCKS[0]
        for j in range(4):
            address = D.p.I + (BN[0] + i) * D.p.S[0] 
            address += BN[1] + j 			# fecthing the memory address such that we can grab the intensity from the data

            temp = DIFF[i] * DIFF[4 + j] # set up Lagrange polynomial numerators.
            temp *= SPACE[i] * SPACE[4 + j] # set up Lagrange polynomial denominators.
            INDEX = II + j 		# Corresponding index in the congiguous intensity array we've built in the init hot.
            # if temp == 0.0: 
            #     printf("\n temp is zero!")
            #     printf(" INDEX: %lu",INDEX)
            #     if DIFF[12+l] == 0.0: printf('DIFF[12+l] is zero')

            
            if CACHE == 1:				# Cache flag that indicates base node was changed.
                I_CACHE[INDEX] = address[0]	# So the intensity value of this value is presumably new and value at this address can be saved in the cache. The only problem I have with this is that if we ever change and then return to this exact value, work will not have been saved.

            I += temp * I_CACHE[INDEX]		# Last step lagrange interpolation to recover intensity value.
                        
                        #printf('i=%d,j=%d,k=%d,l=%d,m=%d, ', <int>i, <int>j, <int>k, <int>l, <int>m)
                        #printf('address = %d, ', <int>(address-D.p.I))   
                        #printf('I_CACHE[INDEX] = %d, ', <int>I_CACHE[INDEX])
                        #printf('temp = %0.2e, ', temp)                         
                        #printf('dI = %0.2e\n', temp * I_CACHE[INDEX])

    #if gsl_isnan(I) == 1:
        #printf("\nIntensity: NaN; Index [%d,%d,%d,%d] ",
                #<int>BN[0], <int>BN[1], <int>BN[2], <int>BN[3])

    #printf("\nBase-nodes [%d,%d,%d,%d] ",
                #<int>BN[0], <int>BN[1], <int>BN[2], <int>BN[3])

    # printf("\nI:  %.8e", I)
    # printf("\nvec[0]: %f", vec[0])
    # printf("\nvec[1]: %f", vec[1])

    if I < 0.0:
        return 0.0
    # Vec here should be the temperature!
    # printf(" I_out: %.8e, ", I * pow(10.0, 3.0 * vec[1]))
    # printf(" I_out: %f, ", I)
    #return I * pow(10.0, 3.0 * Temperature)
    return I

cdef double eval_hot_2D_norm() nogil:
    # Source radiation field normalisation which is independent of the
    # parameters of the parametrised model -- i.e. cell properties, energy,
    # and angle.
    # Writing the normalisation here reduces the number of operations required
    # during integration.
    # The units of the specific intensity need to be J/cm^2/s/keV/steradian.

    return erg / 4.135667662e-18


# cdef double get_param(int param_idx, void *const data) nogil:
#     cdef DATA *D = <DATA*> data
    
#     cdef int i
#     cdef double[:] param
    
#     for i in range(D.p.N[param_idx]):
#          param[i] = D.p.params[param_idx][i]
    
#     return param
    
# cdef double eval_hot_faster(size_t THREAD,
#                      double E,
#                      double mu,
#                      const double *const VEC,
#                      void *const data) nogil:
    
#     # Arguments:
#     # E = photon energy in keV
#     # mu = cosine of ray zenith angle (i.e., angle to surface normal)
#     # VEC = variables such as temperature, effective gravity, ...
#     # data = numerical model data required for intensity evaluation
#     # This function must cast the void pointer appropriately for use.
#     cdef DATA *D = <DATA*> data

#     cdef:
#         size_t i = 0, ii
#         double I = 0.0, temp
#         double *node_vals = D.acc.node_vals[THREAD]
#         size_t *BN = D.acc.BN[THREAD]
#         double *SPACE = D.acc.SPACE[THREAD]
#         double *DIFF = D.acc.DIFF[THREAD]
#         double *I_CACHE = D.acc.INTENSITY_CACHE[THREAD]
#         double *V_CACHE = D.acc.VEC_CACHE[THREAD]
#         double vec[5] # should be = ndims
#         # double E_eff = k_B_over_keV * pow(10.0, VEC[0])
#         # double E_eff = k_B_over_keV * pow(10.0, Temperature)
#         int update_baseNode[5]  # should be = ndims
#         int CACHE 
# cdef double eval_ho= 0

#     cdef double te, tbb, tau
#     te = VEC[0]
#     tbb = VEC[1]
#     tau = VEC[2]
    
#     cdef double evere = 0.5109989e6 # electron volts in elecron rest energy

#     # take into account the order of *._hot_atmosphere
#     vec[0] = te
#     vec[1] = tbb
#     vec[2] = tau
#     vec[3] = mu #0.5  # mu
#     vec[4] = E*1e3/evere #1.01088  # E
#     # printf("Bobrikova atmosphere interpolator")
    
#     # printf("diagnostics 0:\n")
#     # printf("E: %.2e, ", E)
#     # printf("Temperature: %.2e, ", Temperature)
#     # printf("k_B_over_keV: %.2e, ", k_B_over_keV)
#     # printf("E_eff: %.2e, ", E_eff)
#     # printf("vec[4]: %.2e\n", vec[4])
    

#     # printf("\nvec[0]: %.8e, ", vec[0])
#     # printf("vec[1]: %.8e, ", vec[1])
#     # printf("vec[2]: %.8e, ", vec[2])
#     # printf("vec[3]: %.8e, ", vec[3])
#     # printf("vec[4]: %.8e, ", vec[4])
    
    
#     #printf("\neval_hot() called")
#     #printf("\nVEC[0]: %f", VEC[0])
#     #printf("\nVEC[1]: %f", VEC[1])

#     while i < D.p.ndims:
#         # if parallel == 31:
#         # printf("\nDimension: %d", <int>i)
#         update_baseNode[i] = 0
#         if vec[i] < node_vals[2*i] and BN[i] != 0:
#             # if parallel == 31:
#             # printf("\nExecute block 1: %d", <int>i)
#             update_baseNode[i] = 1
#             while vec[i] < D.p.params[i][BN[i] + 1]:
#                 # if parallel == 31:
#                 #     printf("\n!")
#                 #     printf("\nvec i: %.8e", vec[i])
#                 #     printf("\nBase node: %d", <int>BN[i])
#                 if BN[i] > 0:
#                     BN[i] -= 1
#                 elif vec[i] <= D.p.params[i][0]:
#                     vec[i] = D.p.params[i][0]
#                     break
#                 elif BN[i] == 0:
#                     break

#             node_vals[2*i] = D.p.params[i][BN[i] + 1]
#             node_vals[2*i + 1] = D.p.params[i][BN[i] + 2]

#             # if parallel == 31:
#             # printf("\nEnd Block 1: %d", <int>i)

#         elif vec[i] > node_vals[2*i + 1] and BN[i] != D.p.N[i] - 4: # I believe this has to do with the cubic interpolation points, so this remains 4
#             # if parallel == 31:
#             # printf("\nExecute block 2: %d", <int>i)
#             update_baseNode[i] = 1
#             while vec[i] > D.p.params[i][BN[i] + 2]:
#                 if BN[i] < D.p.N[i] - 4:
#                     BN[i] += 1
#                 elif vec[i] >= D.p.params[i][D.p.N[i] - 1]:
#                     vec[i] = D.p.params[i][D.p.N[i] - 1]
#                     break
#                 elif BN[i] == D.p.N[i] - 4:
#                     break

#             node_vals[2*i] = D.p.params[i][BN[i] + 1]
#             node_vals[2*i + 1] = D.p.params[i][BN[i] + 2]

#             # if parallel == 31:
#             # printf("\nEnd Block 2: %d", <int>i)

#         # if parallel == 31:
#         # printf("\nTry block 3: %d", <int>i)

#         if V_CACHE[i] != vec[i] or update_baseNode[i] == 1:
#             # if parallel == 31:
#             # printf("\nExecute block 3: %d", <int>i)
#             ii = 4*i
#             DIFF[ii] = vec[i] - D.p.params[i][BN[i] + 1]
#             DIFF[ii] *= vec[i] - D.p.params[i][BN[i] + 2]
#             DIFF[ii] *= vec[i] - D.p.params[i][BN[i] + 3]

#             DIFF[ii + 1] = vec[i] - D.p.params[i][BN[i]]
#             DIFF[ii + 1] *= vec[i] - D.p.params[i][BN[i] + 2]
#             DIFF[ii + 1] *= vec[i] - D.p.params[i][BN[i] + 3]

#             DIFF[ii + 2] = vec[i] - D.p.params[i][BN[i]]
#             DIFF[ii + 2] *= vec[i] - D.p.params[i][BN[i] + 1]
#             DIFF[ii + 2] *= vec[i] - D.p.params[i][BN[i] + 3]

#             DIFF[ii + 3] = vec[i] - D.p.params[i][BN[i]]
#             DIFF[ii + 3] *= vec[i] - D.p.params[i][BN[i] + 1]
#             DIFF[ii + 3] *= vec[i] - D.p.params[i][BN[i] + 2]

#             # printf("\nupdating V_CACHE")


#             V_CACHE[i] = vec[i]

#             # if parallel == 31:
#             # printf("\nEnd block 3: %d", <int>i)

#         # if parallel == 31:
#         #     printf("\nTry block 4: %d", <int>i)

#         if update_baseNode[i] == 1:
#             # if parallel == 31:
#             # printf("\nExecute block 4: %d", <int>i)
#             # printf("i=%d, ", <int>i)
#             # printf("D.p.params[i][BN[i]]: %.2e\n", D.p.params[i][BN[i]])
#             CACHE = 1
#             SPACE[ii] = 1.0 / (D.p.params[i][BN[i]] - D.p.params[i][BN[i] + 1])
#             SPACE[ii] /= D.p.params[i][BN[i]] - D.p.params[i][BN[i] + 2]
#             SPACE[ii] /= D.p.params[i][BN[i]] - D.p.params[i][BN[i] + 3]

#             SPACE[ii + 1] = 1.0 / (D.p.params[i][BN[i] + 1] - D.p.params[i][BN[i]])
#             SPACE[ii + 1] /= D.p.params[i][BN[i] + 1] - D.p.params[i][BN[i] + 2]
#             SPACE[ii + 1] /= D.p.params[i][BN[i] + 1] - D.p.params[i][BN[i] + 3]

#             SPACE[ii + 2] = 1.0 / (D.p.params[i][BN[i] + 2] - D.p.params[i][BN[i]])
#             SPACE[ii + 2] /= D.p.params[i][BN[i] + 2] - D.p.params[i][BN[i] + 1]
#             SPACE[ii + 2] /= D.p.params[i][BN[i] + 2] - D.p.params[i][BN[i] + 3]

#             SPACE[ii + 3] = 1.0 / (D.p.params[i][BN[i] + 3] - D.p.params[i][BN[i]])
#             SPACE[ii + 3] /= D.p.params[i][BN[i] + 3] - D.p.params[i][BN[i] + 1]
#             SPACE[ii + 3] /= D.p.params[i][BN[i] + 3] - D.p.params[i][BN[i] + 2]

#             # if parallel == 31:
#             # printf("\nEnd block 4: %d", <int>i)

#         # printf("\ncomputing DIFFs and SPACEs\n")
#         # printf("DIFF[ii]: %.2e, ", DIFF[ii])
#         # printf("DIFF[ii+1]: %.2e, ", DIFF[ii+1])
#         # printf("DIFF[ii+2]: %.2e, ", DIFF[ii+2])
#         # printf("DIFF[ii+3]: %.2e\n", DIFF[ii+3])
        
        
#         # printf("SPACE[ii]: %.2e, ", SPACE[ii])
#         # printf("SPACE[ii+1]: %.2e, ", SPACE[ii+1])
#         # printf("SPACE[ii+2]: %.2e, ", SPACE[ii+2])
#         # printf("SPACE[ii+3]: %.2e\n", SPACE[ii+3])

#         i += 1

#     cdef size_t j, k, l, m, INDEX, II, JJ, KK, LL
#     cdef double *address = NULL
#     # cdef double *address = <double*>malloc(iteration_size * sizeof(double))


#     cdef int iterator
#     cdef int iteration_size = 1024
#     cdef double* temp2 = <double*>malloc(iteration_size * sizeof(double))

#     for iterator in range(iteration_size):
#         # printf("i:%lu, ",i)
#         m = iterator % 4
#         l = iterator / 4 % 4
#         k = iterator / 16 % 4
#         j = iterator / 64 % 4
#         i = iterator / 256 % 4
        
        
#         address = D.p.I + (BN[0] + i) * D.p.S[0] + (BN[1] + j) * D.p.S[1] + (BN[2] + k) * D.p.S[2] + (BN[3] + l) * D.p.S[3] + BN[4] + m
 
#     # for iterator in range(iteration_size):
#     #     # printf("i:%lu, ",i)
#     #     m = iterator % 4
#     #     l = iterator / 4 % 4
#     #     k = iterator / 16 % 4
#     #     j = iterator / 64 % 4
#     #     i = iterator / 256 % 4
    
#         temp2[iterator] = DIFF[i] * DIFF[4 + j] * DIFF[8 + k] * DIFF[12 + l] * DIFF[16 + m] * SPACE[i] * SPACE[4 + j] * SPACE[8 + k] * SPACE[12 + l] * SPACE[16 + m]
        
        
        
#     # for iterator in range(iteration_size):
#     #     # printf("i:%lu, ",i)
#     #     m = iterator % 4
#     #     l = iterator / 4 % 4
#     #     k = iterator / 16 % 4
#     #     j = iterator / 64 % 4
#     #     i = iterator / 256 % 4

#         if CACHE == 1:
#             # printf('replace %f ',I_CACHE[INDEX])
#             # printf('with %f\n', address[0])
#             I_CACHE[iterator] = address[0]
#         # printf('%f\n', I2)
#         # if temp2[iterator] != 0.0:
#         I += temp2[iterator] *  I_CACHE[iterator]
#         # I2 += temp *  I_CACHE[iterator]
#         # printf('I2: %f\n',I2)
#     free(temp2)
#     # free(address)

#     if I < 0.0:
#          return 0.0

#     return I

