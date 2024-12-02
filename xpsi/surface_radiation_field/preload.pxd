from libc.stdlib cimport malloc, free

ctypedef struct _preloaded:
    size_t ndims
    double *intensity

    double **params
    size_t *S
    size_t *N
    size_t *BLOCKS

cdef _preloaded* init_preload(object atmosphere)

cdef int free_preload(_preloaded *const preloaded) nogil
