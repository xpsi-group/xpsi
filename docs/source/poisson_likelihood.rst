.. module:: XPSI

.. _poisson_likelihood:

Poisson likelihood
==================

.. code-block:: cython

    #cython: cdivision=True
    #cython: boundscheck=False
    #cython: nonecheck=False
    #cython: wraparound=False

    from __future__ import division, print_function

    import numpy as np
    from cython.parallel cimport *
    from libc.math cimport pow, log, floor
    from libc.stdio cimport printf
    from libc.stdlib cimport calloc, malloc, free

    from GSL cimport (gsl_strerror,
                      gsl_spline,
                      gsl_spline_alloc, 
                      gsl_spline_init,
                      gsl_spline_free,
                      gsl_spline_eval,
                      gsl_interp_steffen,
                      gsl_interp_accel,
                      gsl_interp_accel_alloc,
                      gsl_interp_accel_free,
                      gsl_interp_accel_reset,
                      gsl_isnan,
                      gsl_isinf)

    ctypedef gsl_interp_accel accel

    def eval_loglike(size_t N_Ts,
                     double[::1] phase_edges,
                     double[:,::1] counts,
                     double[:,::1] source,
                     double[:,::1] background,
                     double[::1] source_phases,
                     double phase_shift):

        """ Evaluate the Poisson likelihood. """

        cdef:
            signed int ii
            size_t i, j, T
            double LOGLIKE

            gsl_spline **spline = <gsl_spline**> malloc(sizeof(gsl_spline*) * N_Ts)
            accel **acc = <accel**> malloc(N_Ts * sizeof(accel*))

            double[:,::1] STAR = np.zeros((source.shape[0], phase_edges.shape[0]),
                                          dtype = np.double)

        for T in range(N_Ts):
            acc[T] = gsl_interp_accel_alloc()
            spline[T] = gsl_spline_alloc(gsl_interp_steffen, source_phases.shape[0])

        cdef double PHASE, TOTAL_STAR = 0.0, SCALE_STAR, TOTAL_BG = 0.0, SCALE_BG

        for ii in prange(<signed int>source.shape[0],
                         nogil = True,
                         schedule = 'static',
                         num_threads = N_Ts,
                         chunksize = 1):
            i = <size_t> ii
            T = threadid()

            gsl_interp_accel_reset(acc[T])
            gsl_spline_init(spline[T], &(source_phases[0]), &(source[i,0]), source_phases.shape[0])

            for j in range(phase_edges.shape[0]):
                PHASE = phase_edges[j] + phase_shift

                if PHASE > 1.0:
                    PHASE = PHASE - 1.0
                elif PHASE < 0.0:
                    PHASE = 1.0 + PHASE

                STAR[i,j] = gsl_spline_eval(spline[T], PHASE, acc[T])

        for T in range(N_Ts):
            gsl_interp_accel_free(acc[T])
            gsl_spline_free(spline[T])

        free(spline)
        free(acc)

        for i in range(STAR.shape[0]):
            for j in range(STAR.shape[1]):
                TOTAL_STAR += STAR[i,j]
                TOTAL_BG += background[i,j]

        if TOTAL_STAR > 0.0:
            SCALE_STAR = 1.0e4 / TOTAL_STAR

        SCALE_BG = 1.0e4 / TOTAL_BG

        LOGLIKE = 0.0

        cdef double EXPEC = 0.0
        for i in range(STAR.shape[0]):
            for j in range(STAR.shape[1]):
                EXPEC = background[i,j] * SCALE_BG
                if TOTAL_STAR > 0.0:
                    EXPEC += STAR[i,j] * SCALE_STAR

                LOGLIKE -= EXPEC
                LOGLIKE += counts[i,j] * log(EXPEC)

        return LOGLIKE