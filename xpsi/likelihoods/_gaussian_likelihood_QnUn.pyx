#cython: cdivision=True
#cython: boundscheck=False
#cython: nonecheck=False
#cython: wraparound=False
#cython: embedsignature=True

from __future__ import division

import numpy as np
cimport numpy as np
from libc.math cimport pow, log, floor, fabs, pi
from libc.stdlib cimport malloc, free
from libc.stdio cimport printf

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

ctypedef np.uint8_t uint8

from ..tools.core cimport _get_phase_interpolant, gsl_interp_type

def gaussian_likelihood_QnUn(double[::1] phases,
                                        double[:,::1] counts,
                                        double[:,::1] errors,
                                        double[:] star):
    """ Evaluate the Gaussian likelihood.

    The count rate is integrated over phase intervals.

    :param double[::1] phases:
        A :class:`numpy.ndarray` of phase interval edges in cycles.

    :param  double[::1] star:
        Model normalized Stokes signal converted to data phases
        
    :returns:
        A tuple ``(double, 2D ndarray)``. The first element is
        the logarithm of the marginal likelihood. The second element is the
        expected count numbers in joint phase-channel intervals from the star
        (the target source).

    """
    cdef:
        double LOGLIKE = 0.0, EXPEC = 0.0, sigma_tot2 = 1.0
        double n = <double>(phases.shape[0] - 1)


    for i in range(counts.shape[0]):
        for j in range(counts.shape[1]):
            sigma_tot2 = pow(errors[i,j],2.0)
            norm = 0.5 * log(2.0*pi*sigma_tot2)
            LOGLIKE -= ((star[j])-counts[i,j])**2/(2.0*sigma_tot2)-norm
            #printf("\nLOGLIKE, star, counts, err: %f, %f, %f, %f\n", ((star[j])-counts[i,j])**2/(2.0*sigma_tot2)-norm, star[j],counts[i,j],errors[i,j])            


    return (LOGLIKE, np.asarray(star, order='C', dtype=np.double))
