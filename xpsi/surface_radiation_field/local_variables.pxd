from xpsi.pixelmesh.geometricConfiguration cimport _GEOM

ctypedef struct storage:
    double **local_variables
    void *ptr

cdef int eval_local_variables(double theta,
                              double phi,
                              double HYPERSLICE,
                              const _GEOM *const GEOM,
                              const double *const global_variables,
                              const storage *const buf,
                              size_t THREAD) noexcept nogil

cdef storage* init_local_variables(size_t numTHREADS,
                                   const char *const filepath) noexcept nogil

cdef int free_local_variables(size_t numTHREADS, storage *const buf) noexcept nogil

cdef int HIT_or_MISS(double theta,
                     double phi,
                     double HYPERSLICE,
                     const double *const global_variables,
                     const storage *const buf) noexcept nogil
