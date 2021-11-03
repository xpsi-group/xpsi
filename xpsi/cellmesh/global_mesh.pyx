#cython: cdivision=True
#cython: boundscheck=False
#cython: nonecheck=False
#cython: wraparound=False

import numpy as np
cimport numpy as np
cimport cython
from cython.parallel cimport *
from libc.math cimport M_PI, sqrt, sin, cos, tan, asin, acos, atan, fabs
from libc.stdlib cimport malloc, free
from libc.stdio cimport printf

from xpsi.cellmesh.mesh_tools cimport *

cdef double _pi = M_PI
cdef double _2pi = 2.0 * M_PI
cdef double _hpi = 0.5 * M_PI
cdef double c = 2.99792458e8

def construct_closed_cellMesh(size_t numThreads,
                              size_t sqrt_numCell,
                              size_t numCell,
                              double mass,
                              double r_s,
                              double R_eq,
                              double zeta,
                              double epsilon):
    """
    Construct a closed photospheric cell mesh.

    """

    cdef size_t numCellsCompute

    numCellsCompute = sqrt_numCell / 2

    cdef:
        signed int ii
        size_t i, thread
        double cellArea, mu
        double f, radius, r_s_over_r
        double colat_Limits[2]
        double[::1] areaInterpArray = np.zeros(1000, dtype = np.double)
        double[::1] maxEmissionAngle = np.zeros(sqrt_numCell, dtype = np.double)
        double[::1] cellRadialCoord = np.zeros(sqrt_numCell, dtype = np.double)
        double[::1] cos_gamma = np.zeros(sqrt_numCell, dtype = np.double)
        double[::1] effGrav = np.zeros(sqrt_numCell, dtype = np.double)
        double[::1] colatInterpArray = np.linspace(_hpi, 0.0, 1000)
        double[::1] parallels = np.zeros(numCellsCompute - 1, dtype = np.double)
        double[::1] cellColatitudes = np.zeros(sqrt_numCell, dtype = np.double)
        gsl_interp_accel *accel = gsl_interp_accel_alloc()
        gsl_interp *interp = gsl_interp_alloc(gsl_interp_steffen, 1000)

    cdef gsl_cq_work **w = <gsl_cq_work**> malloc(numThreads * sizeof(gsl_cq_work*))
    for k in range(numThreads):
        w[k] = gsl_integration_cquad_workspace_alloc(100)

    cellArea = _2pi * integrateArea(0.0, _hpi, R_eq, epsilon, zeta, 0, w[0]) / (<double> numCell)

    cellArea *= 2.0

    eta = cellArea * (<double> sqrt_numCell) / _2pi

    for ii in prange(1000,
                     nogil = True,
                     schedule = 'static',
                     num_threads = numThreads,
                     chunksize = 1):
        i = <size_t> ii
        thread = threadid()
        areaInterpArray[i] = integrateArea(colatInterpArray[i],
                                           _hpi,
                                           R_eq,
                                           epsilon,
                                           zeta,
                                           0,
                                           w[thread])

    # Compute cell boundary parallels
    cdef double i_eta = (<double> numCellsCompute) * eta
    cdef double *area_ptr = &areaInterpArray[0]
    cdef double *colat_ptr = &colatInterpArray[0]
    gsl_interp_init(interp, area_ptr, colat_ptr, 1000)

    for i in range(numCellsCompute - 1):
        i_eta -= eta
        parallels[i] = gsl_interp_eval(interp, area_ptr, colat_ptr,
                                       i_eta, accel)

    gsl_interp_free(interp)
    gsl_interp_accel_free(accel)

    for ii in prange(<signed int>numCellsCompute, nogil = True, schedule = 'static',
                    num_threads = numThreads, chunksize = 1):
        i = <size_t> ii
        thread = threadid()
        if (i == numCellsCompute - 1):
            cellColatitudes[i] = integrateArea(parallels[i - 1], _hpi,
                                                R_eq, epsilon, zeta, 1,
                                                w[thread]) / eta
        elif (i == 0):
            cellColatitudes[i] = integrateArea(0.0, parallels[i],
                                                R_eq, epsilon, zeta, 1,
                                                w[thread]) / eta
        else:
            cellColatitudes[i] = integrateArea(parallels[i - 1], parallels[i],
                                                R_eq, epsilon, zeta, 1,
                                                w[thread]) / eta

        mu = cos(cellColatitudes[i])
        radius = radiusNormalised(mu, epsilon, zeta)
        cellRadialCoord[i] = radius * R_eq
        #r_s_over_r = r_s / cellRadialCoord[i]
        effGrav[i] = effectiveGravity(mu, R_eq, zeta, epsilon)
        f = f_theta(mu, radius, epsilon, zeta)
        cos_gamma[i] = 1.0 / sqrt(1.0 + f * f)
        maxEmissionAngle[i] = _hpi + acos(cos_gamma[i])

    cdef double delta_phi

    delta_phi = _2pi / (<double>sqrt_numCell)
    phi = np.linspace(_2pi - 0.5 * delta_phi, 0.5 * delta_phi, sqrt_numCell)

    for i in range(numCellsCompute, sqrt_numCell):
        cellColatitudes[i] = _pi - cellColatitudes[sqrt_numCell - i - 1]
        cellRadialCoord[i] = cellRadialCoord[sqrt_numCell - i - 1]
        cos_gamma[i] = cos_gamma[sqrt_numCell - i - 1]
        effGrav[i] = effGrav[sqrt_numCell - i - 1]
        maxEmissionAngle[i] = maxEmissionAngle[sqrt_numCell - i - 1]

    phi, theta = np.meshgrid(phi,
                             np.asarray(cellColatitudes, dtype = np.double, order='C'))

    for k in range(numThreads):
        gsl_integration_cquad_workspace_free(w[k])
    free(w)

    return (theta,
            phi,
            np.asarray(cellRadialCoord, dtype = np.double, order='C'),
            cellArea,
            np.asarray(maxEmissionAngle, dtype = np.double, order='C'),
            np.asarray(cos_gamma, dtype = np.double, order='C'),
            np.asarray(effGrav, dtype = np.double, order='C'))
