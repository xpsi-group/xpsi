#cython: cdivision=True
#cython: boundscheck=False
#cython: nonecheck=False
#cython: wraparound=False

import numpy as np
from cython.parallel cimport *
from libc.stdlib cimport malloc, free
from libc.math cimport log10, pow

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

def energy_interpolator(size_t N_Ts,
                        double[:,::1] FLUX,
                        double[::1] energies,
                        double[::1] new_energies):

    cdef:
        signed int ii
        int mode
        size_t i, j, T
        double *cpy

        double **flux = <double**> malloc(sizeof(double*) * N_Ts)
        gsl_interp **interp = <gsl_interp**> malloc(sizeof(gsl_interp*) * N_Ts)
        accel **acc = <accel**> malloc(N_Ts * sizeof(accel*))

        double[:,::1] interp_flux = np.zeros((FLUX.shape[1],
                                              new_energies.shape[0]),
                                              dtype = np.double)

    for T in range(N_Ts):
        acc[T] = gsl_interp_accel_alloc()
        interp[T] = gsl_interp_alloc(gsl_interp_steffen, energies.shape[0])
        flux[T] = <double*> malloc(sizeof(double) * energies.shape[0])

    for ii in prange(<signed int>FLUX.shape[1],
                     nogil = True,
                     schedule = 'static',
                     num_threads = N_Ts,
                     chunksize = 1):
        i = <size_t> ii
        T = threadid()
        cpy = flux[T]

        mode = 1
        for j in range(energies.shape[0]):
            if FLUX[j,i] <= 0.0:
                mode = 0

        for j in range(energies.shape[0]):
            if mode == 1:
                cpy[j] = log10(FLUX[j,i])
            else:
                cpy[j] = FLUX[j,i]

        gsl_interp_accel_reset(acc[T])
        gsl_interp_init(interp[T], &(energies[0]), cpy, energies.shape[0])

        for j in range(new_energies.shape[0]):
            interp_flux[i,j] = gsl_interp_eval(interp[T],
                                               &(energies[0]),
                                               cpy,
                                               new_energies[j],
                                               acc[T])

            if mode == 1:
                interp_flux[i,j] = pow(10.0, interp_flux[i,j])
            else:
                interp_flux[i,j] = interp_flux[i,j]

    for T in range(N_Ts):
        gsl_interp_accel_free(acc[T])
        gsl_interp_free(interp[T])
        free(flux[T])

    free(flux)
    free(interp)
    free(acc)

    return np.asarray(interp_flux.T, dtype = np.double, order = 'C')
