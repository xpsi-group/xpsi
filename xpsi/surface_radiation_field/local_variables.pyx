#cython: cdivision=True
#cython: boundscheck=False
#cython: nonecheck=False
#cython: wraparound=False

from libc.math cimport sin, cos, acos, log, exp, pow, sqrt
from libc.stdlib cimport malloc, free
from libc.stdio cimport printf

from xpsi.tools.effective_gravity_universal cimport effectiveGravity

cdef int SUCCESS = 0
cdef int ERROR = 1
cdef int HIT = 1
cdef int MISS = 0

cdef int HIT_or_MISS(double theta,
                     double phi,
                     double HYPERSLICE,
                     const double *const global_variables) nogil:
    # Use this function to determine whether an arbitrary input ray transports
    # a FINITE quantity of radiation
    # HYPERSLICE is a phase parameter with units of radians that the
    # time evolution of the photosphere as it undergoes bulk rotation

    return HIT

#----------------------------------------------------------------------->>>
# >>> User modifiable functions.
# >>> Note that the user is entirely free to wrap thread-safe and
# ... non-parallel external C routines from an external library.
# >>> Thus the bodies of the following need not be written explicitly in
# ... the Cython language.
#----------------------------------------------------------------------->>>
cdef void* init_local_variables(size_t numThreads) nogil:
    # Return NULL if dynamic memory is not required for the model

    cdef double **local_variables = <double**> malloc(numThreads * sizeof(double*))
    cdef size_t T

    for T in range(numThreads):
        local_variables[T]= <double*> malloc(2 * sizeof(double))

    # Explicit cast not required, but explicit is clearer
    return <void*> local_variables

cdef int free_local_variables(size_t numThreads,
                              void *const local_variables) nogil:
    # Just use free(<void*> local_variables) iff no memory was dynamically
    # allocated in the function:
    #   init_local_variables()
    # because local_variables is expected to be NULL in this case

    cdef double** local_vars = <double**> local_variables
    cdef size_t T

    for T in range(numThreads):
        free(local_vars[T])

    free(local_vars)

    return SUCCESS

cdef int eval_local_variables(double theta,
                                      double phi,
                                      double HYPERSLICE,
                                      const _GEOM *const GEOM,
                                      const double *const global_variables,
                                      size_t THREAD,
                                      void *const local_variables) nogil:

    cdef:
        double *local_vars = (<double**> local_variables)[THREAD]

    local_vars[0] = global_variables[0]

    local_vars[1] = effectiveGravity(cos(theta),
                                          GEOM.R_eq,
                                          GEOM.zeta,
                                          GEOM.epsilon)

    return SUCCESS
