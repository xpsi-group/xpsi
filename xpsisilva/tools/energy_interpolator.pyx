#cython: cdivision=True
#cython: boundscheck=False
#cython: nonecheck=False
#cython: wraparound=False
#cython: embedsignature=True

import numpy as np
from cython.parallel cimport *
from libc.stdlib cimport malloc, free
from libc.math cimport log10, pow

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

from . cimport _get_energy_interpolant, gsl_interp_type

def energy_interpolator(size_t N_Ts,
                        double[:,::1] signal,
                        double[::1] energies,
                        double[::1] new_energies):
    """ Interpolate a signal in energy.

    :param size_t N_Ts:
        Number of OpenMP threads to spawn.

    :param double[:,::1] signal:
        A C-contiguous :class:`numpy.ndarray` of an energy-resolved (specific
        _signal) signal. Energy increases with row number.

    :param double[::1] energies:
        A :class:`numpy.ndarray` of the logarithms (base 10) of the energies at
        which the :obj:`signal` is resolved.

    :param double[::1] new_energies:
        A :class:`numpy.ndarray` of the logarithm (base 10) of the energies at
        which to interpolate.

    :returns:
        A 2D :class:`numpy.ndarray` of the signal interpolated at the new set
        of energies. Energy increases with row number.

    """
    cdef const gsl_interp_type *_interpolant

    _interpolant = _get_energy_interpolant()

    cdef:
        signed int ii
        int mode
        size_t i, j, T
        double *cpy
        double max_energy = energies[energies.shape[0] - 1]

        double **_signal = <double**> malloc(sizeof(double*) * N_Ts)
        gsl_interp **interp = <gsl_interp**> malloc(sizeof(gsl_interp*) * N_Ts)
        accel **acc = <accel**> malloc(N_Ts * sizeof(accel*))

        double[:,::1] interp_signal = np.zeros((signal.shape[1],
                                               new_energies.shape[0]),
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

        mode = 1
        for j in range(<size_t>energies.shape[0]):
            if signal[j,i] <= 0.0:
                mode = 0

        for j in range(<size_t>energies.shape[0]):
            if mode == 1:
                cpy[j] = log10(signal[j,i])
            else:
                cpy[j] = signal[j,i]

        gsl_interp_accel_reset(acc[T])
        gsl_interp_init(interp[T], &(energies[0]), cpy, energies.shape[0])

        for j in range(<size_t>new_energies.shape[0]):

            if new_energies[j] > max_energy: # extrapolate by setting to zero
                continue

            interp_signal[i,j] = gsl_interp_eval(interp[T],
                                               &(energies[0]),
                                               cpy,
                                               new_energies[j],
                                               acc[T])

            if mode == 1:
                interp_signal[i,j] = pow(10.0, interp_signal[i,j])
            else:
                interp_signal[i,j] = interp_signal[i,j]

    for T in range(N_Ts):
        gsl_interp_accel_free(acc[T])
        gsl_interp_free(interp[T])
        free(_signal[T])

    free(_signal)
    free(interp)
    free(acc)

    return np.asarray(interp_signal.T, dtype = np.double, order = 'C')
