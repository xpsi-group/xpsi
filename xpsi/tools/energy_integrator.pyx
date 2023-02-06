#cython: cdivision=True
#cython: boundscheck=False
#cython: nonecheck=False
#cython: wraparound=False
#cython: embedsignature=True

import numpy as np
from cython.parallel cimport *
from libc.stdlib cimport malloc, free
from libc.math cimport pow, log

from GSL cimport (gsl_interp,
                   gsl_interp_alloc,
                   gsl_interp_init,
                   gsl_interp_free,
                   gsl_interp_eval,
                   gsl_interp_eval_integ,
                   gsl_interp_accel,
                   gsl_interp_accel_alloc,
                   gsl_interp_accel_free,
                   gsl_interp_accel_reset)

ctypedef gsl_interp_accel accel

from . cimport _get_phase_interpolant, gsl_interp_type

def energy_integrator(size_t N_Ts,
                       double[:,::1] signal,
                       double[::1] energies,
                       double[::1] energy_edges):
    """ Integrate a signal over energy intervals.

    :param size_t N_Ts:
        Number of OpenMP threads to spawn.

    :param double[:,::1] signal:
        A C-contiguous :class:`numpy.ndarray` of an energy-resolved (specific
        flux) signal. Energy increases with row number.

    :param double[::1] energies:
        A :class:`numpy.ndarray` of the logarithms (base 10) of the energies.

    :param double[::1] energy_edges:
        A :class:`numpy.ndarray` of the logarithm (base 10) of the energy
        interval edges.

    :returns:
        A 2D :class:`numpy.ndarray` of the signal integrated over energy
        intervals. Energy interval number increases with row number.

    """
    cdef const gsl_interp_type *_interpolant

    _interpolant = _get_phase_interpolant()

    cdef:
        signed int ii
        size_t i, T #j, T
        unsigned int j
        double *cpy
        double upper_energy
        double max_energy = energies[energies.shape[0] - 1]

        double **_signal = <double**> malloc(sizeof(double*) * N_Ts)
        gsl_interp **interp = <gsl_interp**> malloc(sizeof(gsl_interp*) * N_Ts)
        accel **acc = <accel**> malloc(N_Ts * sizeof(accel*))

        double[:,::1] binned_signal = np.zeros((signal.shape[1],
                                               energy_edges.shape[0] - 1),
                                               dtype = np.double)

    for T in range(N_Ts):
        acc[T] = gsl_interp_accel_alloc()
        interp[T] = gsl_interp_alloc(_interpolant, energies.shape[0])
        _signal[T] = <double*> malloc(sizeof(double) * energies.shape[0])

    for ii in prange(<signed int>signal.shape[1],
                     nogil = True,
                     schedule = 'static',
                     num_threads = N_Ts,
                     chunksize = 1):
        i = <size_t> ii
        T = threadid()
        cpy = _signal[T]

        for j in range(energies.shape[0]):
            cpy[j] = pow(10.0, energies[j]) * signal[j,i] * log(10.0)

        gsl_interp_accel_reset(acc[T])
        gsl_interp_init(interp[T], &(energies[0]), cpy, energies.shape[0])

        for j in range(energy_edges.shape[0] - 1):
            if energy_edges[j + 1] > max_energy:
                upper_energy = max_energy
            else:
                upper_energy = energy_edges[j + 1]
            binned_signal[i,j] = gsl_interp_eval_integ(interp[T],
                                                       &(energies[0]),
                                                       cpy,
                                                       energy_edges[j],
                                                       upper_energy,
                                                       acc[T])
            if energy_edges[j + 1] > max_energy:
                break

    for T in range(N_Ts):
        gsl_interp_accel_free(acc[T])
        gsl_interp_free(interp[T])
        free(_signal[T])

    free(_signal)
    free(interp)
    free(acc)

    return np.asarray(binned_signal.T, dtype = np.double, order = 'C')
