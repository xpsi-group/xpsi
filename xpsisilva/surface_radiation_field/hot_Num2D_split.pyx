#cython: cdivision=True
#cython: boundscheck=False
#cython: nonecheck=False
#cython: wraparound=False

from libc.math cimport M_PI, sqrt, sin, cos, acos, log10, pow, exp, fabs
from libc.stdio cimport printf, fopen, fclose, fread, FILE
#from GSL cimport gsl_isnan, gsl_isinf
from libc.stdlib cimport malloc, free

from xpsisilva.global_imports import _keV, _k_B

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

    cdef double *address = NULL

    for T in range(numThreads): #For the full interpolation hypercube, store all intensities in an array with the right shape, so all values are lookupable later by knowing the i,j address.
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
        double vec[2] # should be = ndims
        int update_baseNode[2]  # should be = ndims
        int CACHE = 0

    # The input value of the parameter (vec) to be interpolated. Note this is the order of *._hot_atmosphere
    vec[0] = mu
    vec[1] = E

    while i < D.p.ndims: 					# For each dimension 
        update_baseNode[i] = 0					# Flag to change base node
        if vec[i] < node_vals[2*i] and BN[i] != 0:		# If the input value of the parameter (vec) is smaller than the base node value (we are valued below the hypercube/basenode interval), and BN is not the first (it is zet zero in init_hot)
            update_baseNode[i] = 1				# then the base node should be changed.
            while vec[i] < D.p.params[i][BN[i] + 1]:		# while the input value remains smaller than the value at the next basenode
                if BN[i] > 0:					# and if the basenode is not the first
                    BN[i] -= 1					# decrement by 1
                elif vec[i] <= D.p.params[i][0]:		# or, if it IS not bigger than zero, and if vec is smaller than the first value
                    vec[i] = D.p.params[i][0]			# make it equal to the first value
                    break
                elif BN[i] == 0:				# or, if it IS not bigger than zero, and vec is bigger than the first value, then it should be zero.
                    break

            node_vals[2*i] = D.p.params[i][BN[i] + 1] 		# now update node values with updated basenode values. Footnote: Not sure why it is +1 and +2. It must be correctly surrounding the input value. Oh perhaps this is so that the first two are below and the second two are above?
            node_vals[2*i + 1] = D.p.params[i][BN[i] + 2]

        elif vec[i] > node_vals[2*i + 1] and BN[i] != D.p.N[i] - 4: 	# Opposite case of above. If the input value is larger than the second node_vals and BN is not the last
            update_baseNode[i] = 1 				# The base node should be changed.
            while vec[i] > D.p.params[i][BN[i] + 2]:		# While the input value remains larger than the value at the third base node
                if BN[i] < D.p.N[i] - 4:			# if the base node is smaller than the last - 4
                    BN[i] += 1					# then it may be incremented.
                elif vec[i] >= D.p.params[i][D.p.N[i] - 1]:	# or, if the base node IS NOT smaller than the last - 4, and IF it is larger than one before the last value (footnote: i don't think it could be larger than the last value because the while statement already checks with 2 values above. So here we are at the limit, but not over it.)
                    vec[i] = D.p.params[i][D.p.N[i] - 1]	# then it should be set to that value.
                    break
                elif BN[i] == D.p.N[i] - 4:			# Or, if it IS NOT smaller than the last - 4 and vec is not too big, then it should be exactly the last base node.
                    break

            node_vals[2*i] = D.p.params[i][BN[i] + 1]		# This is exactly the same, update the node_vals from the base nodes.
            node_vals[2*i + 1] = D.p.params[i][BN[i] + 2]

        if V_CACHE[i] != vec[i] or update_baseNode[i] == 1:	# If the input value is not already the same as the cached value (so e.g. if the the change in the vec was super super small and data is spread out), or more obviously if the base node was changed.
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

            V_CACHE[i] = vec[i]				# Store this input value in the cache for next time so that work can be skipped.

        if update_baseNode[i] == 1:				# For the denominators you have to redo the work if the basenode was changed.

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

        i += 1							# For each dimension (while loop)

    cdef size_t j, INDEX, II
    cdef double *address = NULL
        
    # Combinatorics over nodes of hypercube; weight cgs intensities
    for i in range(4):
        II = i * D.p.BLOCKS[0]
        for j in range(4):
            address = D.p.I + (BN[0] + i) * D.p.S[0] 
            address += BN[1] + j 			# fecthing the memory address such that we can grab the intensity from the data

            temp = DIFF[i] * DIFF[4 + j] # set up Lagrange polynomial numerators.
            temp *= SPACE[i] * SPACE[4 + j] # set up Lagrange polynomial denominators.
            INDEX = II + j 		# Corresponding index in the congiguous intensity array we've built in the init hot.
            
            if CACHE == 1:				# Cache flag that indicates base node was changed.
                I_CACHE[INDEX] = address[0]	# So the intensity value of this value is presumably new and value at this address can be saved in the cache. The only problem I have with this is that if we ever change and then return to this exact value, work will not have been saved.

            I += temp * I_CACHE[INDEX]		# Last step lagrange interpolation to recover intensity value.

    return I

cdef double eval_hot_2D_I(size_t THREAD,
                     double E,
                     double mu,
                     void *const data) nogil:
    # Arguments:
    # E = photon energy in keV
    # mu = cosine of ray zenith angle (i.e., angle to surface normal)
    # VEC = variables such as temperature, effective gravity, ...
    # data = numerical model data required for intensity evaluation

    cdef double I = eval_hot_2D(THREAD,E,mu,data)

    if I < 0.0:
        return 0.0

    return I


cdef double eval_hot_2D_Q(size_t THREAD,
                     double E,
                     double mu,
                     void *const data) nogil:
    # Arguments:
    # E = photon energy in keV
    # mu = cosine of ray zenith angle (i.e., angle to surface normal)
    # VEC = variables such as temperature, effective gravity, ...
    # data = numerical model data required for intensity evaluation

    return eval_hot_2D(THREAD,E,mu,data)

cdef double eval_hot_2D_norm() nogil:
    # Source radiation field normalisation which is independent of the
    # parameters of the parametrised model -- i.e. cell properties, energy,
    # and angle.
    # Writing the normalisation here reduces the number of operations required
    # during integration.
    # The units of the specific intensity need to be J/cm^2/s/keV/steradian.

    return erg / 4.135667662e-18

