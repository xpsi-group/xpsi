cdef size_t atmos_extension_elsewhere = 1

from libc.stdio cimport printf

#Blackbody
from xpsi.surface_radiation_field.hot_BB cimport (init_hot_BB,
                                                     free_hot_BB,
                                                     eval_hot_BB,
                                                     eval_hot_norm_BB)
#4D-Numerical
from xpsi.surface_radiation_field.hot_Num4D cimport (init_hot_Num4D,
                                                     free_hot_Num4D,
                                                     eval_hot_Num4D,
                                                     eval_hot_norm_Num4D)

                                                     
#----------------------------------------------------------------------->>>
cdef void* init_elsewhere(size_t numThreads, const _preloaded *const preloaded, size_t atm_ext) nogil:
    global atmos_extension_elsewhere
    atmos_extension_elsewhere=atm_ext
    if atmos_extension_elsewhere == 1:
        return init_hot_BB(numThreads, preloaded)
    else:
        return init_hot_Num4D(numThreads, preloaded)

cdef int free_elsewhere(size_t numThreads, void *const data) nogil:
    if atmos_extension_elsewhere == 1:
        return free_hot_BB(numThreads, data)
    else:
        return free_hot_Num4D(numThreads, data)

cdef double eval_elsewhere(size_t THREAD,
                     double E,
                     double mu,
                     const double *const VEC,
                     void *const data) nogil:

    if atmos_extension_elsewhere == 1:
        I_hot = eval_hot_BB(THREAD,E,mu,VEC,data)
    else:
        I_hot = eval_hot_Num4D(THREAD,E,mu,VEC,data)
    return I_hot

cdef double eval_elsewhere_norm() nogil:
    if atmos_extension_elsewhere == 1:
        return eval_hot_norm_BB()
    else:
        return eval_hot_norm_Num4D()
