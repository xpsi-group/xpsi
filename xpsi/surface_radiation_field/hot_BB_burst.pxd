from .preload cimport _preloaded

cdef double eval_hot_BB_burst_I(size_t THREAD,
                     double E,
                     double mu,
                     const double *const VEC,
                     void *const data) nogil

cdef double eval_hot_BB_burst_Q(size_t THREAD,
                     double E,
                     double mu,
                     const double *const VEC,
                     void *const data) nogil

cdef double eval_hot_norm_BB_burst() nogil

cdef void* init_hot_BB_burst(size_t numThreads, const _preloaded *const preloaded) nogil

cdef int free_hot_BB_burst(size_t numThreads, void *const data) nogil
