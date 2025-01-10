from .preload cimport _preloaded

cdef double eval_hot_I(size_t THREAD,
                     double E,
                     double mu,
                     const double *const VEC,
                     void *const data,
                     size_t beam_opt) noexcept nogil

cdef double eval_hot_Q(size_t THREAD,
                     double E,
                     double mu,
                     const double *const VEC,
                     void *const data,
                     size_t beam_opt) noexcept nogil

cdef double eval_hot_norm() noexcept nogil

cdef void* init_hot(size_t numThreads, const _preloaded *const preloaded, size_t atm_ext) noexcept nogil

cdef int free_hot(size_t numThreads, void *const data) noexcept nogil
