cdef double eval_extRadField(size_t THREAD,
                             double E,
                             double mu,
                             const double *const VEC,
                             void *const data) nogil

cdef double eval_extRadField_norm() nogil

cdef void* init_extRadField(size_t numThreads,
                            const extRadField_PRELOAD *const preload) nogil

cdef int free_extRadField(size_t numThreads, void *const data) nogil

ctypedef struct extRadField_PRELOAD:
    double **params
    size_t *S
    double *I
