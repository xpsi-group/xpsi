from .preload cimport _preloaded

cdef double eval_elsewhere_user(size_t THREAD,
                           double E,
                           double mu,
                           const double *const VEC,
                           void *const data) noexcept nogil

cdef double eval_elsewhere_norm_user() noexcept nogil

cdef void* init_elsewhere_user(size_t numThreads,
                          const _preloaded *const preloaded) noexcept nogil

cdef int free_elsewhere_user(size_t numThreads, void *const data) noexcept nogil
