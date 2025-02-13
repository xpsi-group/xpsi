.. _likelihoods:

Likelihoods
===========

.. automodule:: xpsi.likelihoods

Extensions
^^^^^^^^^^

.. autofunction:: xpsi.likelihoods.precomputation

.. autofunction:: xpsi.likelihoods.eval_marginal_likelihood

.. autofunction:: xpsi.likelihoods.poisson_likelihood_given_background

Example
^^^^^^^

The following is an example of a custom extension module to evaluate a
simple likelihood function based on Poisson sampling distribution of count
data.

.. code-block:: cython

    #cython: cdivision=True
    #cython: boundscheck=False
    #cython: nonecheck=False
    #cython: wraparound=False
    #cython: embedsignature=True

    from __future__ import division

    import numpy as np
    from libc.math cimport pow, log, floor
    from libc.stdlib cimport malloc, free

    from GSL cimport (gsl_interp,
                       gsl_interp_alloc,
                       gsl_interp_init,
                       gsl_interp_free,
                       gsl_interp_eval,
                       gsl_interp_eval_integ,
                       gsl_interp_steffen,
                       gsl_interp_accel,
                       gsl_interp_accel_alloc,
                       gsl_interp_accel_free,
                       gsl_interp_accel_reset)

    ctypedef gsl_interp_accel accel

    def poisson_likelihood_given_background(double exposure_time,
                                            double[::1] phases,
                                            double[:,::1] counts,
                                            components,
                                            component_phases,
                                            phase_shifts,
                                            double[:,::1] background):
        """ Evaluate the Poisson likelihood.

        The count rate is integrated over phase intervals.

        :param double exposure_time:
            Exposure time in seconds by which to scale the expected count rate
            in each phase interval.

        :param phases:
            A C-contiguous :class:`numpy.ndarray` of phase interval edges in cycles.

        :param tuple components:
            Component signals, each a C-contiguous :class:`numpy.ndarray` of
            signal count rates where phase increases with column number.

        :param tuple component_phases:
            For each component, a C-contiguous :class:`numpy.ndarray` of phases
            in cycles at which the model :obj:`signal` is evaluated on
            the interval ``[0,1]``.

        :param array-like phase_shifts:
            Phase shifts in cycles, such as on the interval ``[-0.5,0.5]``, for
            the component signals.

        :param background:
            A C-contiguous :class:`numpy.ndarray` of background expected
            *counts*, whose shape matches the number of channels in each element
            of :obj:`components` and the number of phase intervals constructed
            from :obj:`phases`.

        """

        cdef:
            size_t i, j, p, num_components = len(components)
            double a, b

            double[:,::1] STAR = np.zeros((components[0].shape[0], phases.shape[0]-1),
                                           dtype = np.double)

        cdef double *phases_ptr = NULL
        cdef double *signal_ptr = NULL

        cdef double[:,::1] signal
        cdef double[::1] signal_phase_set
        cdef double phase_shift

        cdef gsl_interp **interp = <gsl_interp**> malloc(num_components * sizeof(gsl_interp*))
        cdef accel **acc =  <accel**> malloc(num_components * sizeof(accel*))

        for p in range(num_components):
            signal_phase_set = component_phases[p]
            interp[p] = gsl_interp_alloc(gsl_interp_steffen, signal_phase_set.shape[0])
            acc[p] = gsl_interp_accel_alloc()
            gsl_interp_accel_reset(acc[p])

        cdef gsl_interp *inter_ptr = NULL
        cdef accel *acc_ptr = NULL

        for i in range(STAR.shape[0]):
            for p in range(num_components):
                signal = components[p]
                signal_phase_set = component_phases[p]
                phase_shift = phase_shifts[p]

                interp_ptr = interp[p]
                acc_ptr = acc[p]
                phases_ptr = &(signal_phase_set[0])
                signal_ptr = &(signal[i,0])

                gsl_interp_init(interp_ptr, phases_ptr, signal_ptr,
                                signal_phase_set.shape[0])

                for j in range(phases.shape[0] - 1):
                    a = phases[j] + phase_shift
                    b = phases[j+1] + phase_shift

                    a -= floor(a)
                    b -= floor(b)

                    if a < b:
                        STAR[i,j] += gsl_interp_eval_integ(interp_ptr,
                                                           phases_ptr,
                                                           signal_ptr,
                                                           a, b,
                                                           acc_ptr)
                    else:
                        STAR[i,j] += gsl_interp_eval_integ(interp_ptr,
                                                           phases_ptr,
                                                           signal_ptr,
                                                           a, 1.0,
                                                           acc_ptr)

                        STAR[i,j] += gsl_interp_eval_integ(interp_ptr,
                                                           phases_ptr,
                                                           signal_ptr,
                                                           0.0, b,
                                                           acc_ptr)

        for p in range(num_components):
            gsl_interp_accel_free(acc[p])
            gsl_interp_free(interp[p])

        free(acc)
        free(interp)

        cdef:
            double LOGLIKE = 0.0, EXPEC = 0.0
            double n = <double>(phases.shape[0] - 1)

        for i in range(STAR.shape[0]):
            for j in range(STAR.shape[1]):
                EXPEC = (STAR[i,j] + background[i,j]/n) * exposure_time
                LOGLIKE -= EXPEC
                LOGLIKE += counts[i,j] * log(EXPEC)

                STAR[i,j] += background[i,j]/n
                STAR[i,j] *= exposure_time

        return (LOGLIKE, np.asarray(STAR, order='C', dtype=np.double))
