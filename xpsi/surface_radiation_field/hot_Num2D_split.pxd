from .preload cimport _preloaded

cdef double eval_hot_2D(size_t THREAD,
                     double E,
                     double mu,
                     void *const data) noexcept nogil

cdef double eval_hot_2D_I(size_t THREAD,
                     double E,
                     double mu,
                     void *const data) noexcept nogil

cdef double eval_hot_2D_Q(size_t THREAD,
                     double E,
                     double mu,
                     void *const data) noexcept nogil

cdef double eval_hot_2D_norm() noexcept nogil

cdef void* init_hot_2D(size_t numThreads, const _preloaded *const preloaded) noexcept nogil

cdef int free_hot_2D(size_t numThreads, void *const data) noexcept nogil
