from .preload cimport _preloaded

cdef double eval_hot_Num4D(size_t THREAD,
                     double E,
                     double mu,
                     const double *const VEC,
                     void *const data) nogil

cdef double eval_hot_norm_Num4D() nogil

cdef void* init_hot_Num4D(size_t numThreads, const _preloaded *const preloaded) nogil

cdef int free_hot_Num4D(size_t numThreads, void *const data) nogil

cdef double* produce_2D_data_Num4D(size_t THREAD, const double *const VEC, void *const data) nogil

cdef object make_atmosphere_2D_Num4D(double *I_data, void *const data)