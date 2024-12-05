from .preload cimport _preloaded

cdef double eval_hot(size_t THREAD,
                     double E,
                     double mu,
                     const double *const VEC,
                     void *const data) noexcept nogil

cdef double eval_hot_Num5D_I(size_t THREAD,
                     double E,
                     double mu,
                     const double *const VEC,
                     void *const data) noexcept nogil

cdef double eval_hot_Num5D_Q(size_t THREAD,
                     double E,
                     double mu,
                     const double *const VEC,
                     void *const data) noexcept nogil

cdef double eval_hot_norm_Num5D() noexcept nogil

cdef void* init_hot_Num5D(size_t numThreads, const _preloaded *const preloaded) noexcept nogil

cdef int free_hot_Num5D(size_t numThreads, void *const data) noexcept nogil

cdef double* produce_2D_data(size_t THREAD, const double *const VEC, void *const data) noexcept nogil

cdef object make_atmosphere_2D(double *I_data, void *const data)
