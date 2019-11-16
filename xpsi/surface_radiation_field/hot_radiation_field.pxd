cdef double eval_hotRadField(size_t THREAD,
                             double E,
                             double mu,
                             const double *const VEC,
                             void *const data) nogil

cdef double eval_hotRadField_norm() nogil

cdef void* init_hotRadField(size_t numThreads,
                            const hotRadField_PRELOAD *const preload) nogil

cdef int free_hotRadField(size_t numThreads, void *const data) nogil

ctypedef struct hotRadField_PRELOAD:
    double **params
    size_t *S
    double *I
