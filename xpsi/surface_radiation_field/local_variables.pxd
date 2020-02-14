from xpsi.pixelmesh.geometricConfiguration cimport _GEOM

cdef int eval_local_variables(double theta,
                              double phi,
                              double HYPERSLICE,
                              const _GEOM *const GEOM,
                              const double *const global_variables,
                              size_t THREAD,
                              void *const local_variables) nogil

cdef void* init_local_variables(size_t numThreads) nogil

cdef int free_local_variables(size_t numThreads,
                              void *const local_variables) nogil

cdef int HIT_or_MISS(double theta,
                     double phi,
                     double HYPERSLICE,
                     const double *const global_variables) nogil
