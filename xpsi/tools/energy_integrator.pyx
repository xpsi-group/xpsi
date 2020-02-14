#cython: cdivision=True
#cython: boundscheck=False
#cython: nonecheck=False
#cython: wraparound=False

import numpy as np
from cython.parallel cimport *
from libc.stdlib cimport malloc, free
from libc.math cimport pow, log

from GSL cimport (gsl_interp,
                   gsl_interp_steffen,
                   gsl_interp_alloc,
                   gsl_interp_init,
                   gsl_interp_free,
                   gsl_interp_eval,
                   gsl_interp_eval_integ,
                   gsl_interp_accel,
                   gsl_interp_accel_alloc,
                   gsl_interp_accel_free,
                   gsl_interp_accel_reset,
                   gsl_isnan,
                   gsl_isinf,
                   GSL_SUCCESS)

ctypedef gsl_interp_accel accel

def energy_integrator(size_t N_Ts,
                       double[:,::1] integrand,
                       double[::1] energies,
                       double a,
                       double b):

    cdef:
        signed int ii
        size_t i, j, T
        double *cpy

        double **flux = <double**> malloc(sizeof(double*) * N_Ts)
        gsl_interp **interp = <gsl_interp**> malloc(sizeof(gsl_interp*) * N_Ts)
        accel **acc = <accel**> malloc(N_Ts * sizeof(accel*))

        double[::1] integrated = np.zeros(integrand.shape[0], dtype = np.double)

    for T in range(N_Ts):
        acc[T] = gsl_interp_accel_alloc()
        interp[T] = gsl_interp_alloc(gsl_interp_steffen, energies.shape[0])
        flux[T] = <double*> malloc(sizeof(double) * energies.shape[0])

    for ii in prange(<signed int>integrand.shape[0],
                     nogil = True,
                     schedule = 'static',
                     num_threads = N_Ts,
                     chunksize = 1):
        i = <size_t> ii
        T = threadid()
        cpy = flux[T]

        for j in range(energies.shape[0]):
            cpy[j] = pow(10.0, energies[j]) * integrand[i,j] * log(10.0)

        gsl_interp_accel_reset(acc[T])
        gsl_interp_init(interp[T], &(energies[0]), cpy, energies.shape[0])

        integrated[i] = gsl_interp_eval_integ(interp[T],
                                              &(energies[0]),
                                              cpy,
                                              a,
                                              b,
                                              acc[T])
    for T in range(N_Ts):
        gsl_interp_accel_free(acc[T])
        gsl_interp_free(interp[T])
        free(flux[T])

    free(flux)
    free(interp)
    free(acc)

    return np.asarray(integrated, dtype = np.double, order = 'C')
