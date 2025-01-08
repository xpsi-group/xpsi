#cython: cdivision=True
#cython: boundscheck=False
#cython: nonecheck=False
#cython: wraparound=False

cdef _preloaded* init_preload(object atmosphere):
    cdef size_t i, j
    cdef double[::1] cast
    cdef double[::1] intensity

    preloaded = <_preloaded*> malloc(sizeof(_preloaded))

    preloaded.ndims = <size_t>(len(atmosphere) - 1)

    preloaded.params = <double**> malloc(sizeof(double*) * preloaded.ndims)
    preloaded.N = <size_t*> malloc(sizeof(size_t) * (preloaded.ndims))
    preloaded.S = <size_t*> malloc(sizeof(size_t) * (preloaded.ndims - 1))
    preloaded.BLOCKS = <size_t*> malloc(sizeof(size_t) * (preloaded.ndims - 1))

    for i in range(preloaded.ndims):
        cast = atmosphere[i]
        preloaded.N[i] = cast.shape[0]
        preloaded.params[i] = &cast[0]
        if i < preloaded.ndims - 1:
            cast = atmosphere[i+1]
            preloaded.S[i] = cast.shape[0]
            if i < preloaded.ndims - 2:
                for j in range(i+2, preloaded.ndims):
                    cast = atmosphere[j]
                    preloaded.S[i] *= cast.shape[0]

    intensity = atmosphere[i+1]
    preloaded.intensity = &intensity[0]

    return preloaded

cdef int free_preload(_preloaded *const preloaded) noexcept nogil:
    free(preloaded.params)
    free(preloaded.N)
    free(preloaded.S)
    free(preloaded.BLOCKS)
    free(preloaded)

    return 0
