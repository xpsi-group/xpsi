cdef double eval_srcRadField(size_t THREAD,
                             double E,
                             double mu,
                             const double *const VEC,
                             void *const data) nogil

cdef double eval_srcRadField_norm() nogil

cdef void* init_srcRadField(size_t numThreads,
                            const srcRadField_PRELOAD *const preload) nogil

cdef int free_srcRadField(size_t numThreads, void *const data) nogil

ctypedef struct srcRadField_PRELOAD:
    double **params
    size_t *S
    double *I
