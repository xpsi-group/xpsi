cdef double eval_elsewhereRadField(size_t THREAD,
                                   double E,
                                   double mu,
                                   const double *const VEC,
                                   void *const data) nogil

cdef double eval_elsewhereRadField_norm() nogil

cdef void* init_elsewhereRadField(size_t numThreads,
                                  const elsewhereRadField_PRELOAD *const preload) nogil

cdef int free_elsewhereRadField(size_t numThreads, void *const data) nogil

ctypedef struct elsewhereRadField_PRELOAD:
    double **params
    size_t *S
    double *I
