from .preload cimport _preloaded

cdef double eval_hot_Num4D(size_t THREAD,
                     double E,
                     double mu,
                     const double *const VEC,
                     void *const data) noexcept nogil

cdef double eval_hot_norm_Num4D() noexcept nogil

cdef void* init_hot_Num4D(size_t numThreads, const _preloaded *const preloaded) noexcept nogil

cdef int free_hot_Num4D(size_t numThreads, void *const data) noexcept nogil
