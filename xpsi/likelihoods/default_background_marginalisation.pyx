#cython: cdivision=True
#cython: boundscheck=False
#cython: nonecheck=False
#cython: wraparound=False

from __future__ import print_function

import numpy as np
from libc.math cimport exp, pow, log, sqrt, fabs

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
                   gsl_interp_accel_reset,
                   gsl_isnan,
                   gsl_isinf,
                   gsl_function,
                   gsl_integration_cquad_workspace,
                   gsl_integration_cquad_workspace_alloc,
                   gsl_integration_cquad_workspace_free,
                   gsl_integration_cquad)

ctypedef gsl_interp_accel accel

cdef extern from "gsl/gsl_sf_gamma.h":

    double gsl_sf_lnfact(const unsigned int n)

def precomputation(int[:,::1] data):
    """ Compute negative of sum of log-factorials, channel-by-channel.

    Also precompute maximum data counts in each channel.

    """

    cdef:
        size_t i, j
        double[::1] precomp = np.zeros(data.shape[0], dtype = np.double)

    for i in range(data.shape[0]):
        for j in range(data.shape[1]):
            precomp[i] += gsl_sf_lnfact(<unsigned int>(data[i,j]))

        precomp[i] *= -1.0

    return np.asarray(precomp, dtype = np.double, order = 'C')

ctypedef struct args:
    size_t n
    size_t i
    double SCALE
    double *precomp
    double *data
    double *star
    double B
    double interval
    double T_exp
    double std
    double A

cdef double marginal_integrand(double B, void *params) nogil:

    cdef int j
    cdef double c, x = 0.0
    cdef args *a  = <args*> params

    for j in range(a.n):
        c = a.SCALE * (a.star[j] + B)
        x += a.data[j] * log(c) - c

    #with gil:
    #    if x - a.A > 0.0:
    #        print("a.i, x-a.A, a.B, off, B = %i, %.8f, %.8f, %.8f, %.8f" %
    #              (a.i, x - a.A, a.B, (a.B - B)/a.interval, B))

    return exp(x - a.A)

cdef double delta(double B, void *params) nogil:

    cdef int j
    cdef double x = 0.0, y = 0.0
    cdef args *a  = <args*> params

    for j in range(a.n):
        y += a.data[j] / (a.star[j] + B)
        x += a.data[j] / pow(a.star[j] + B, 2.0)

    y = 2.0 * a.T_exp - 2.0 * y
    x *= 2.0

    a.std = sqrt(2.0 / x)

    return -1.0 * y / x

def eval_loglike_phaseIntervals_maximise(double exposure_time,
                                         double[::1] phases,
                                         double[:,::1] counts,
                                         pulses,
                                         double[::1] pulse_phases,
                                         phase_shifts,
                                         double[::1] neg_sum_ln_data_factorial,
                                         size_t workspace_intervals,
                                         double epsabs,
                                         double epsrel,
                                         double epsilon,
                                         double sigmas,
                                         double llzero):

    """ Evaluate the Poisson likelihood.

    The count rate is integrated over phase intervals.

    A Newton iteration procedure is implemented to approximate the ML point
    given a fiducial background count rate in each channel.

    """

    cdef:
        size_t i, j, p
        double LOGLIKE = 0.0, av_STAR, av_DATA
        double pa, pb, c

        gsl_interp *interp = gsl_interp_alloc(gsl_interp_steffen, pulse_phases.shape[0])
        accel *acc = gsl_interp_accel_alloc()

        gsl_integration_cquad_workspace *w = gsl_integration_cquad_workspace_alloc(workspace_intervals)

        double[:,::1] STAR = np.zeros((pulses[0].shape[0], phases.shape[0] - 1),
                                      dtype = np.double)
        double[::1] MCL_BACKGROUND = np.zeros(pulses[0].shape[0], dtype = np.double)

        double n = <double>(phases.shape[0] - 1)
        double SCALE = exposure_time / n

    cdef args a
    a.n = <size_t>n
    a.precomp = &(neg_sum_ln_data_factorial[0])
    a.SCALE = SCALE

    cdef double result, abserr, upper, lower, std_est
    cdef size_t nevals
    cdef gsl_function f

    f.function = &marginal_integrand

    cdef double *phases_ptr
    cdef double *pulse_ptr

    cdef double dB, B_min, min_counts, limit
    cdef int counter
    cdef size_t num_pulses = len(pulses)
    cdef double[:,::1] pulse
    cdef double phase_shift

    for i in range(STAR.shape[0]):
        for p in range(num_pulses):
            pulse = pulses[p]
            phase_shift = phase_shifts[p]

            gsl_interp_accel_reset(acc)
            phases_ptr = &(pulse_phases[0])
            pulse_ptr = &(pulse[i,0])
            gsl_interp_init(interp, phases_ptr, pulse_ptr, pulse_phases.shape[0])

            for j in range(phases.shape[0] - 1):
                pa = phases[j] + phase_shift
                pb = phases[j+1] + phase_shift

                if pa > 1.0:
                    pa -= 1.0
                elif pa < 0.0:
                    pa += 1.0

                if pb > 1.0:
                    pb -= 1.0
                elif pb < 0.0:
                    pb += 1.0

                if pa < pb:
                    STAR[i,j] += gsl_interp_eval_integ(interp, phases_ptr,
                                                      pulse_ptr, pa, pb, acc)
                else:
                    STAR[i,j] += gsl_interp_eval_integ(interp, phases_ptr,
                                                      pulse_ptr, pa, 1.0, acc)
                    STAR[i,j] += gsl_interp_eval_integ(interp, phases_ptr,
                                                       pulse_ptr,  0.0, pb, acc)

        av_DATA = 0.0; av_STAR = 0.0

        for j in range(phases.shape[0] - 1):
            STAR[i,j] *= n

            av_STAR += STAR[i,j]
            av_DATA += counts[i,j]

        limit = av_STAR * SCALE - 20.0 * sqrt(av_STAR * SCALE) - av_DATA

        if limit > 0.0:
            LOGLIKE = llzero * (0.1 + 0.9 * np.random.rand(1))
            break

        av_STAR /= n
        av_DATA /= exposure_time

        a.i = i
        a.data = &(counts[i,0])
        a.star = &(STAR[i,0])
        a.T_exp = exposure_time

        B_min = 0.0
        B = av_DATA - av_STAR
        if B <= B_min:
            min_counts = -1.0
            for j in range(0, a.n):
                if STAR[i,j] == 0.0:
                    min_counts = -2.0
                    break
            if min_counts == -2.0:
                for j in range(0, a.n):
                    if min_counts == -2.0 or 0.0 < counts[i,j] < min_counts:
                        min_counts = counts[i,j]
            if min_counts != -1.0:
                B = 0.01 * min_counts / SCALE
                B_min = 0.1 * B
                #print("i, B = %i, %.16f" % (a.i,B))
            else:
                B = B_min

        dB = delta(B, &a)

        #print("B_0 = %.16e" % B)
        counter = 0
        while fabs(dB) > epsilon*a.std and counter < 2:
            B += dB
            #print("i, B = %i, %.8f" % (a.i,B))
            if B < B_min:
                if B_min > 0.0:
                    counter += 1
                else:
                    counter = 2
                B = B_min

            dB = delta(B, &a)

        std_est = 0.0
        for j in range(a.n):
            #print("i, j, STAR[i,j] + B = %i, %i, %.16e, %.16e" % (a.i,j,STAR[i,j], B))
            std_est += counts[i,j] / pow(STAR[i,j] + B, 2.0)

        std_est = sqrt(1.0 / std_est)
        #print("i, B, std = %i, %.16e, %.16e" % (a.i,B,std_est))

        a.B = B
        lower = B - sigmas * std_est

        if lower < B_min:
            lower = B_min

        upper = B + sigmas * std_est

        a.interval = upper - lower

        f.params = &a

        a.A = 0.0
        for j in range(a.n):
            c = a.SCALE * (a.star[j] + B)
            a.A += a.data[j] * log(c) - c

        gsl_integration_cquad(&f, lower, upper,
                              epsabs, epsrel,
                              w, &result,
                              &abserr, &nevals)
        if result > 0.0:
            LOGLIKE += log(result) + a.A + a.precomp[a.i]
        else:
            LOGLIKE += llzero

        for j in range(a.n):
            STAR[i,j] = a.SCALE * (a.star[j] + B)

        MCL_BACKGROUND[i] = B * exposure_time

    gsl_interp_accel_free(acc)
    gsl_interp_free(interp)

    gsl_integration_cquad_workspace_free(w)

    return (LOGLIKE,
            np.asarray(STAR, dtype=np.double, order='C'),
            np.asarray(MCL_BACKGROUND, dtype=np.double, order='C'))
