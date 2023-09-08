from .preload cimport _preloaded

cdef double eval_hot_Num2D(size_t THREAD,
                     double E,
                     double mu,
                     const double *const VEC,
                     void *const data) nogil
                     
cdef double eval_hot_Num2D(size_t THREAD,
                     double E,
                     double mu,
                     const double *const VEC,
                     void *const data) nogil

cdef double eval_hot_Num2D_I(size_t THREAD,
                     double E,
                     double mu,
                     const double *const VEC,
                     void *const data) nogil

cdef double eval_hot_Num2D_Q(size_t THREAD,
                     double E,
                     double mu,
                     const double *const VEC,
                     void *const data) nogil

cdef double eval_hot_norm_Num2D() nogil

cdef void* init_hot_Num2D(size_t numThreads, const _preloaded *const preloaded) nogil

cdef int free_hot_Num2D(size_t numThreads, void *const data) nogil
