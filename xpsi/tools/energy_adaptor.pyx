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

def energy_adaptor(double[::1] signal,
                   double[::1] energies,
                   size_t num_energies):

    cdef:
        size_t i

        gsl_interp *interp = gsl_interp_alloc(gsl_interp_steffen, energies.shape[0])
        accel *acc = gsl_interp_accel_alloc()

        double norm

        double[::1] masses = np.zeros(energies.shape[0], dtype = np.double)
        double[::1] target_masses = np.linspace(0.0, 1.0, num_energies,
                                                dtype = np.double)
        double[::1] adapted_energies = np.zeros(num_energies,
                                                dtype = np.double)

    gsl_interp_accel_reset(acc)
    gsl_interp_init(interp, &(energies[0]), &(signal[0]), energies.shape[0])

    norm = gsl_interp_eval_integ(interp,
                                 &(energies[0]),
                                 &(signal[0]),
                                 energies[0],
                                 energies[energies.shape[0] - 1],
                                 acc)

    for i in range(energies.shape[0]):
        masses[i] = gsl_interp_eval_integ(interp,
                                          &(energies[0]),
                                          &(signal[0]),
                                          energies[0],
                                          energies[i],
                                          acc)

        masses[i] /= norm

    for i in range(masses.shape[0] - 1):
        if masses[i] >= masses[i+1]:
            gsl_interp_accel_free(acc)
            gsl_interp_free(interp)
            return np.linspace(energies[0],
                               energies[energies.shape[0] - 1],
                               num_energies - 1,
                               dtype = np.double)

    gsl_interp_accel_reset(acc)
    gsl_interp_init(interp, &(masses[0]), &(energies[0]), masses.shape[0])

    for i in range(target_masses.shape[0]):
        adapted_energies[i] = gsl_interp_eval(interp,
                                              &(masses[0]),
                                              &(energies[0]),
                                              target_masses[i],
                                              acc)

    gsl_interp_accel_free(acc)
    gsl_interp_free(interp)

    return np.asarray(adapted_energies, dtype = np.double, order = 'C')
