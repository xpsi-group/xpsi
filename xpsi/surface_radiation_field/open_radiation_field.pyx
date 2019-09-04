#cython: cdivision=True
#cython: boundscheck=False
#cython: nonecheck=False
#cython: wraparound=False

from libc.math cimport exp, pow

cdef int SUCCESS = 0

cdef double erg = 1.0e-7
cdef double Planck_dist_const = 5.040366110812353e22
cdef double k_B = 1.38064852e-23
cdef double keV = 1.60217662e-16
cdef double k_B_over_keV = k_B / keV

#----------------------------------------------------------------------->>>
# >>> User modifiable functions.
# >>> Note that the user is entirely free to wrap thread-safe and
# ... non-parallel external C routines from an external library.
# >>> Thus the bodies of the following need not be written explicitly in
# ... the Cython language.
#----------------------------------------------------------------------->>>
cdef void* init_extRadField(size_t numThreads, const extRadField_PRELOAD *const preload) nogil:
    # This function must match the free management routine free_extRadField()
    # in terms of freeing dynamically allocated memory. This is entirely
    # the user's responsibility to manage.

    # Return NULL if dynamic memory is not required for the model.
    return NULL

cdef int free_extRadField(size_t numThreads, void *const data) nogil:
    # This function must match the initialisation routine init_extRadField()
    # in terms of freeing dynamically allocated memory. This is entirely
    # the user's responsibility to manage.
    # The void pointer must be appropriately cast before memory is freed --
    # only the user can know this at compile time.
    # Just use free(<void*> data) iff no memory was dynamically
    # allocated in the function:
    #   init_extRadField()
    # because data is expected to be NULL in this case

    #printf("\nNo data to be freed.")

    return SUCCESS

cdef double eval_extRadField(size_t THREAD,
                             double E,
                             double mu,
                             const double *const VEC,
                             void *const data) nogil:

    cdef double E_eff = k_B_over_keV * pow(10.0, VEC[0])

    return E * E * E / ( exp(E / E_eff) - 1.0 )

cdef double eval_extRadField_norm() nogil:
    # Source radiation field normalisation which is independent of the
    # parameters of the parametrised model -- i.e. cell properties, energy,
    # and angle.
    # Writing the normalisation here reduces the number of operations required
    # during integration.

    return erg * Planck_dist_const
