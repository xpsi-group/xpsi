cdef size_t atmos_extension = 1

from libc.stdio cimport printf

#Blackbody
from xpsi.surface_radiation_field.hot1 cimport (init_hot1,
                                                     free_hot1,
                                                     eval_hot1,
                                                     eval_hot_norm1)
#Numerical+beaming
from xpsi.surface_radiation_field.hot2 cimport (init_hot2,
                                                     free_hot2,
                                                     eval_hot2,
                                                     eval_hot_norm2)

                                                     
#----------------------------------------------------------------------->>>
cdef void* init_hot(size_t numThreads, const _preloaded *const preloaded, size_t atm_ext) nogil:
    global atmos_extension
    atmos_extension=atm_ext
    if atmos_extension == 1:
        return init_hot1(numThreads, preloaded)
    else:
        return init_hot2(numThreads, preloaded)

cdef int free_hot(size_t numThreads, void *const data) nogil:
    if atmos_extension == 1:
        return free_hot1(numThreads, data)
    else:
        return free_hot2(numThreads, data)

cdef double eval_hot(size_t THREAD,
                     double E,
                     double mu,
                     const double *const VEC,
                     void *const data) nogil:
    if atmos_extension == 1:
        return eval_hot1(THREAD,E,mu,VEC,data)
    else: 
        return eval_hot2(THREAD,E,mu,VEC,data)

cdef double eval_hot_norm() nogil:
    if atmos_extension == 1:
        return eval_hot_norm1()
    else:
        return eval_hot_norm2()
