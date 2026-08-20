.. _surface_radiation_field:

Surface radiation field
=======================

.. automodule:: xpsi.surface_radiation_field

Wrappers
^^^^^^^^

.. autofunction:: xpsi.surface_radiation_field.intensity

.. autofunction:: xpsi.surface_radiation_field.intensity_from_globals

.. autofunction:: xpsi.surface_radiation_field.effective_gravity

.. _numerical_atmosphere:

A numerical atmosphere
^^^^^^^^^^^^^^^^^^^^^^

The following is the extension module ``hot_Num4D.pyx`` for surface radiation field evaluation
that may be found on path ``surface_radiation_field/hot_Num4D.pyx``.
This extension works with numerical atmosphere preloading to execute cubic polynomial interpolation
in four dimensions.

.. code-block:: cython

    #cython: cdivision=True
    #cython: boundscheck=False
    #cython: nonecheck=False
    #cython: wraparound=False

    from libc.math cimport M_PI, sqrt, sin, cos, acos, log10, pow, exp, fabs
    from libc.stdio cimport printf, fopen, fclose, fread, FILE
    from GSL cimport gsl_isnan, gsl_isinf
    from libc.stdlib cimport malloc, free

    from xpsi.global_imports import _keV, _k_B, _h_keV

    cdef int SUCCESS = 0
    cdef int ERROR = 1

    cdef double erg = 1.0e-7
    cdef double k_B = _k_B
    cdef double keV = _keV
    cdef double h_keV = _h_keV
    cdef double k_B_over_keV = k_B / keV
    cdef int VERBOSE = 0

    ctypedef struct ACCELERATE:
        size_t **BN                # base node for interpolation
        double **node_vals
        double **SPACE
        double **DIFF
        double **INTENSITY_CACHE
        double **VEC_CACHE

    # Modify this struct if useful for the user-defined source radiation field.
    # Note that the members of DATA will be shared by all threads and are
    # statically allocated, whereas the members of ACCELERATE will point to
    # dynamically allocated memory, not shared by threads.

    ctypedef struct DATA:
        const _preloaded *p
        ACCELERATE acc

    #----------------------------------------------------------------------->>>
    # >>> User modifiable functions.
    # >>> Note that the user is entirely free to wrap thread-safe and
    # ... non-parallel external C routines from an external library.
    # >>> Thus the bodies of the following need not be written explicitly in
    # ... the Cython language.
    #----------------------------------------------------------------------->>>
    cdef void* init_hot_Num4D(size_t numThreads, const _preloaded *const preloaded) nogil:
        # This function must match the free management routine free_hot()
        # in terms of freeing dynamically allocated memory. This is entirely
        # the user's responsibility to manage.
        # Return NULL if dynamic memory is not required for the model

        if preloaded == NULL :
            printf("ERROR: The numerical atmosphere data were not preloaded, which are required by this extension.\n")

        cdef DATA *D = <DATA*> malloc(sizeof(DATA))
        D.p = preloaded

        D.p.BLOCKS[0] = 64
        D.p.BLOCKS[1] = 16
        D.p.BLOCKS[2] = 4

        cdef size_t T, i, j, k, l

        D.acc.BN = <size_t**> malloc(numThreads * sizeof(size_t*))
        D.acc.node_vals = <double**> malloc(numThreads * sizeof(double*))
        D.acc.SPACE = <double**> malloc(numThreads * sizeof(double*))
        D.acc.DIFF = <double**> malloc(numThreads * sizeof(double*))
        D.acc.INTENSITY_CACHE = <double**> malloc(numThreads * sizeof(double*))
        D.acc.VEC_CACHE = <double**> malloc(numThreads * sizeof(double*))

        for T in range(numThreads):
            D.acc.BN[T] = <size_t*> malloc(D.p.ndims * sizeof(size_t))
            D.acc.node_vals[T] = <double*> malloc(2 * D.p.ndims * sizeof(double))
            D.acc.SPACE[T] = <double*> malloc(4 * D.p.ndims * sizeof(double))
            D.acc.DIFF[T] = <double*> malloc(4 * D.p.ndims * sizeof(double))
            D.acc.INTENSITY_CACHE[T] = <double*> malloc(256 * sizeof(double))
            D.acc.VEC_CACHE[T] = <double*> malloc(D.p.ndims * sizeof(double))
            for i in range(D.p.ndims):
                D.acc.BN[T][i] = 0
                D.acc.VEC_CACHE[T][i] = D.p.params[i][1]
                D.acc.node_vals[T][2*i] = D.p.params[i][1]
                D.acc.node_vals[T][2*i + 1] = D.p.params[i][2]

                j = 4*i

                D.acc.SPACE[T][j] = 1.0 / (D.p.params[i][0] - D.p.params[i][1])
                D.acc.SPACE[T][j] /= D.p.params[i][0] - D.p.params[i][2]
                D.acc.SPACE[T][j] /= D.p.params[i][0] - D.p.params[i][3]

                D.acc.SPACE[T][j + 1] = 1.0 / (D.p.params[i][1] - D.p.params[i][0])
                D.acc.SPACE[T][j + 1] /= D.p.params[i][1] - D.p.params[i][2]
                D.acc.SPACE[T][j + 1] /= D.p.params[i][1] - D.p.params[i][3]

                D.acc.SPACE[T][j + 2] = 1.0 / (D.p.params[i][2] - D.p.params[i][0])
                D.acc.SPACE[T][j + 2] /= D.p.params[i][2] - D.p.params[i][1]
                D.acc.SPACE[T][j + 2] /= D.p.params[i][2] - D.p.params[i][3]

                D.acc.SPACE[T][j + 3] = 1.0 / (D.p.params[i][3] - D.p.params[i][0])
                D.acc.SPACE[T][j + 3] /= D.p.params[i][3] - D.p.params[i][1]
                D.acc.SPACE[T][j + 3] /= D.p.params[i][3] - D.p.params[i][2]

                D.acc.DIFF[T][j] = D.acc.VEC_CACHE[T][i] - D.p.params[i][1]
                D.acc.DIFF[T][j] *= D.acc.VEC_CACHE[T][i] - D.p.params[i][2]
                D.acc.DIFF[T][j] *= D.acc.VEC_CACHE[T][i] - D.p.params[i][3]

                D.acc.DIFF[T][j + 1] = D.acc.VEC_CACHE[T][i] - D.p.params[i][0]
                D.acc.DIFF[T][j + 1] *= D.acc.VEC_CACHE[T][i] - D.p.params[i][2]
                D.acc.DIFF[T][j + 1] *= D.acc.VEC_CACHE[T][i] - D.p.params[i][3]

                D.acc.DIFF[T][j + 2] = D.acc.VEC_CACHE[T][i] - D.p.params[i][0]
                D.acc.DIFF[T][j + 2] *= D.acc.VEC_CACHE[T][i] - D.p.params[i][1]
                D.acc.DIFF[T][j + 2] *= D.acc.VEC_CACHE[T][i] - D.p.params[i][3]

                D.acc.DIFF[T][j + 3] = D.acc.VEC_CACHE[T][i] - D.p.params[i][0]
                D.acc.DIFF[T][j + 3] *= D.acc.VEC_CACHE[T][i] - D.p.params[i][1]
                D.acc.DIFF[T][j + 3] *= D.acc.VEC_CACHE[T][i] - D.p.params[i][2]

        cdef double *address = NULL
        # Cache intensity
        for T in range(numThreads):
            for i in range(4):
                for j in range(4):
                    for k in range(4):
                        for l in range(4):
                            address = D.p.I + (D.acc.BN[T][0] + i) * D.p.S[0]
                            address += (D.acc.BN[T][1] + j) * D.p.S[1]
                            address += (D.acc.BN[T][2] + k) * D.p.S[2]
                            address += D.acc.BN[T][3] + l
                            D.acc.INTENSITY_CACHE[T][i * D.p.BLOCKS[0] + j * D.p.BLOCKS[1] + k * D.p.BLOCKS[2] + l] = address[0]

        # Cast for generalised usage in integration routines
        return <void*> D


    cdef int free_hot_Num4D(size_t numThreads, void *const data) nogil:
        # This function must match the initialisation routine init_hot()
        # in terms of freeing dynamically allocated memory. This is entirely
        # the user's responsibility to manage.
        # The void pointer must be appropriately cast before memory is freed --
        # only the user can know this at compile time.
        # Just use free(<void*> data) iff no memory was dynamically
        # allocated in the function:
        #   init_hot()
        # because data is expected to be NULL in this case

        cdef DATA *D = <DATA*> data

        cdef size_t T

        for T in range(numThreads):
            free(D.acc.BN[T])
            free(D.acc.node_vals[T])
            free(D.acc.SPACE[T])
            free(D.acc.DIFF[T])
            free(D.acc.INTENSITY_CACHE[T])
            free(D.acc.VEC_CACHE[T])

        free(D.acc.BN)
        free(D.acc.node_vals)
        free(D.acc.SPACE)
        free(D.acc.DIFF)
        free(D.acc.INTENSITY_CACHE)
        free(D.acc.VEC_CACHE)

        free(D)

        return SUCCESS

    #----------------------------------------------------------------------->>>
    # >>> Cubic polynomial interpolation.
    # >>> Improve acceleration properties... i.e. do not recompute numerical
    # ... weights or re-read intensities if not necessary.
    #----------------------------------------------------------------------->>>
    cdef double eval_hot_Num4D(size_t THREAD,
                         double E,
                         double mu,
                         const double *const VEC,
                         void *const data) nogil:
        # Arguments:
        # E = photon energy in keV
        # mu = cosine of ray zenith angle (i.e., angle to surface normal)
        # VEC = variables such as temperature, effective gravity, ...
        # data = numerical model data required for intensity evaluation

        # This function must cast the void pointer appropriately for use.
        cdef DATA *D = <DATA*> data

        cdef:
            size_t i = 0, ii
            double I = 0.0, temp
            double *node_vals = D.acc.node_vals[THREAD]
            size_t *BN = D.acc.BN[THREAD]
            double *SPACE = D.acc.SPACE[THREAD]
            double *DIFF = D.acc.DIFF[THREAD]
            double *I_CACHE = D.acc.INTENSITY_CACHE[THREAD]
            double *V_CACHE = D.acc.VEC_CACHE[THREAD]
            double vec[4]
            double E_eff = k_B_over_keV * pow(10.0, VEC[0])
            int update_baseNode[4]
            int CACHE = 0

        vec[0] = VEC[0]
        vec[1] = VEC[1]
        vec[2] = mu
        vec[3] = log10(E / E_eff)

        while i < D.p.ndims:
            # if parallel == 31:
            #     printf("\nDimension: %d", <int>i)
            update_baseNode[i] = 0
            if vec[i] < node_vals[2*i] and BN[i] != 0:
                # if parallel == 31:
                #     printf("\nExecute block 1: %d", <int>i)
                update_baseNode[i] = 1
                while vec[i] < D.p.params[i][BN[i] + 1]:
                    # if parallel == 31:
                    #     printf("\n!")
                    #     printf("\nvec i: %.8e", vec[i])
                    #     printf("\nBase node: %d", <int>BN[i])
                    if BN[i] > 0:
                        BN[i] -= 1
                    elif vec[i] <= D.p.params[i][0]:
                        vec[i] = D.p.params[i][0]
                        break
                    elif BN[i] == 0:
                        break

                node_vals[2*i] = D.p.params[i][BN[i] + 1]
                node_vals[2*i + 1] = D.p.params[i][BN[i] + 2]

                # if parallel == 31:
                #     printf("\nEnd Block 1: %d", <int>i)

            elif vec[i] > node_vals[2*i + 1] and BN[i] != D.p.N[i] - 4:
                # if parallel == 31:
                #     printf("\nExecute block 2: %d", <int>i)
                update_baseNode[i] = 1
                while vec[i] > D.p.params[i][BN[i] + 2]:
                    if BN[i] < D.p.N[i] - 4:
                        BN[i] += 1
                    elif vec[i] >= D.p.params[i][D.p.N[i] - 1]:
                        vec[i] = D.p.params[i][D.p.N[i] - 1]
                        break
                    elif BN[i] == D.p.N[i] - 4:
                        break

                node_vals[2*i] = D.p.params[i][BN[i] + 1]
                node_vals[2*i + 1] = D.p.params[i][BN[i] + 2]

                # if parallel == 31:
                #     printf("\nEnd Block 2: %d", <int>i)

            # if parallel == 31:
            #     printf("\nTry block 3: %d", <int>i)

            if V_CACHE[i] != vec[i] or update_baseNode[i] == 1:
                # if parallel == 31:
                #     printf("\nExecute block 3: %d", <int>i)
                ii = 4*i
                DIFF[ii] = vec[i] - D.p.params[i][BN[i] + 1]
                DIFF[ii] *= vec[i] - D.p.params[i][BN[i] + 2]
                DIFF[ii] *= vec[i] - D.p.params[i][BN[i] + 3]

                DIFF[ii + 1] = vec[i] - D.p.params[i][BN[i]]
                DIFF[ii + 1] *= vec[i] - D.p.params[i][BN[i] + 2]
                DIFF[ii + 1] *= vec[i] - D.p.params[i][BN[i] + 3]

                DIFF[ii + 2] = vec[i] - D.p.params[i][BN[i]]
                DIFF[ii + 2] *= vec[i] - D.p.params[i][BN[i] + 1]
                DIFF[ii + 2] *= vec[i] - D.p.params[i][BN[i] + 3]

                DIFF[ii + 3] = vec[i] - D.p.params[i][BN[i]]
                DIFF[ii + 3] *= vec[i] - D.p.params[i][BN[i] + 1]
                DIFF[ii + 3] *= vec[i] - D.p.params[i][BN[i] + 2]

                V_CACHE[i] = vec[i]

                # if parallel == 31:
                #     printf("\nEnd block 3: %d", <int>i)

            # if parallel == 31:
            #     printf("\nTry block 4: %d", <int>i)

            if update_baseNode[i] == 1:
                # if parallel == 31:
                #     printf("\nExecute block 4: %d", <int>i)
                CACHE = 1
                SPACE[ii] = 1.0 / (D.p.params[i][BN[i]] - D.p.params[i][BN[i] + 1])
                SPACE[ii] /= D.p.params[i][BN[i]] - D.p.params[i][BN[i] + 2]
                SPACE[ii] /= D.p.params[i][BN[i]] - D.p.params[i][BN[i] + 3]

                SPACE[ii + 1] = 1.0 / (D.p.params[i][BN[i] + 1] - D.p.params[i][BN[i]])
                SPACE[ii + 1] /= D.p.params[i][BN[i] + 1] - D.p.params[i][BN[i] + 2]
                SPACE[ii + 1] /= D.p.params[i][BN[i] + 1] - D.p.params[i][BN[i] + 3]

                SPACE[ii + 2] = 1.0 / (D.p.params[i][BN[i] + 2] - D.p.params[i][BN[i]])
                SPACE[ii + 2] /= D.p.params[i][BN[i] + 2] - D.p.params[i][BN[i] + 1]
                SPACE[ii + 2] /= D.p.params[i][BN[i] + 2] - D.p.params[i][BN[i] + 3]

                SPACE[ii + 3] = 1.0 / (D.p.params[i][BN[i] + 3] - D.p.params[i][BN[i]])
                SPACE[ii + 3] /= D.p.params[i][BN[i] + 3] - D.p.params[i][BN[i] + 1]
                SPACE[ii + 3] /= D.p.params[i][BN[i] + 3] - D.p.params[i][BN[i] + 2]

                # if parallel == 31:
                #     printf("\nEnd block 4: %d", <int>i)

            i += 1

        cdef size_t j, k, l, INDEX, II, JJ, KK
        cdef double *address = NULL
        # Combinatorics over nodes of hypercube; weight cgs intensities
        for i in range(4):
            II = i * D.p.BLOCKS[0]
            for j in range(4):
                JJ = j * D.p.BLOCKS[1]
                for k in range(4):
                    KK = k * D.p.BLOCKS[2]
                    for l in range(4):
                        address = D.p.I + (BN[0] + i) * D.p.S[0]
                        address += (BN[1] + j) * D.p.S[1]
                        address += (BN[2] + k) * D.p.S[2]
                        address += BN[3] + l

                        temp = DIFF[i] * DIFF[4 + j] * DIFF[8 + k] * DIFF[12 + l]
                        temp *= SPACE[i] * SPACE[4 + j] * SPACE[8 + k] * SPACE[12 + l]
                        INDEX = II + JJ + KK + l
                        if CACHE == 1:
                            I_CACHE[INDEX] = address[0]
                        I += temp * I_CACHE[INDEX]

        #if gsl_isnan(I) == 1:
            #printf("\nIntensity: NaN; Index [%d,%d,%d,%d] ",
                    #<int>BN[0], <int>BN[1], <int>BN[2], <int>BN[3])

        #printf("\nBase-nodes [%d,%d,%d,%d] ",
                    #<int>BN[0], <int>BN[1], <int>BN[2], <int>BN[3])

        if I < 0.0:
            return 0.0
        return I * pow(10.0, 3.0 * vec[0])

    cdef double eval_hot_norm_Num4D() nogil:
        # Source radiation field normalisation which is independent of the
        # parameters of the parametrised model -- i.e. cell properties, energy,
        # and angle.
        # Writing the normalisation here reduces the number of operations required
        # during integration.
        # The units of the specific intensity need to be J/cm^2/s/keV/steradian.

        return erg / h_keV

