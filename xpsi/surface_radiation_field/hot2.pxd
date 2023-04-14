from .preload cimport _preloaded

cdef double eval_hot2(size_t THREAD,
                     double E,
                     double mu,
                     const double *const VEC,
                     void *const data) nogil

cdef double eval_hot_norm2() nogil

cdef void* init_hot2(size_t numThreads, const _preloaded *const preloaded) nogil

cdef int free_hot2(size_t numThreads, void *const data) nogil
