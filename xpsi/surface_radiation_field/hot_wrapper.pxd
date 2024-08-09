from .preload cimport _preloaded

cdef double eval_hot_I(size_t THREAD,
                     double E,
                     double mu,
                     const double *const VEC,
                     void *const data,
                     size_t beam_opt) nogil

cdef double eval_hot_Q(size_t THREAD,
                     double E,
                     double mu,
                     const double *const VEC,
                     void *const data,
                     size_t beam_opt) nogil

cdef double eval_hot_norm() nogil

cdef void* init_hot(size_t numThreads, const _preloaded *const preloaded, size_t atm_ext) nogil

cdef int free_hot(size_t numThreads, void *const data) nogil

cdef double* produce_2D_data(size_t THREAD, const double *const VEC, void *const data) nogil

cdef object make_atmosphere_2D(double *I_data, void *const data)