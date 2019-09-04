#cython: cdivision=True
#cython: boundscheck=False
#cython: nonecheck=False
#cython: wraparound=False

from __future__ import division, print_function

import numpy as np
cimport numpy as np
cimport cython
from cython.parallel cimport *
from libc.math cimport M_PI, sqrt, sin, cos, tan, asin, acos, atan, fabs
from libc.stdlib cimport malloc, free
from libc.stdio cimport printf

from xpsi.cellmesh.mesh_tools cimport *

cdef double _2pi = 2.0 * M_PI
cdef double _hpi = 0.5 * M_PI
cdef double c = 2.99792458e8

cdef double M_2PI = 2.0 * M_PI

def construct_spot_cellMesh(size_t numThreads,
                            size_t numCell,
                            size_t sqrt_numCell,
                            double mass,
                            double r_s,
                            double R_eq,
                            double zeta,
                            double epsilon,
                            double cedeRadius,
                            double cedeColatitude,
                            double superRadius,
                            double superColatitude,
                            double superAzimuth):
    """
    Construct a photospheric cell mesh enclosing a spot.

    """
    cdef:
        signed int ii
        size_t i, j, thread
        double cellArea, boundary_phi, mu, delta_phi
        double f, radius, r_s_over_r
        double colat_lims[2]

        gsl_interp_accel *accel = gsl_interp_accel_alloc()
        gsl_interp *interp = gsl_interp_alloc(gsl_interp_steffen, 1000)

    colat_lims[0] = cedeColatitude - cedeRadius
    colat_lims[1] = cedeColatitude + cedeRadius

    cdef gsl_cq_work **w = <gsl_cq_work**> malloc(numThreads * sizeof(gsl_cq_work*))
    for thread in range(numThreads):
        w[thread] = gsl_integration_cquad_workspace_alloc(100)

    # Right-spherical triangle
    boundary_phi = asin(sin(cedeRadius) / sin(cedeColatitude))

    cellArea = 2.0 * boundary_phi

    cellArea *= integrateArea(colat_lims[0],
                              colat_lims[1],
                              R_eq,
                              epsilon,
                              zeta,
                              0,
                              w[0])

    delta_phi = 2.0 * boundary_phi / (<double>sqrt_numCell)
    phi = np.linspace(-boundary_phi + 0.5 * delta_phi,
                       boundary_phi - 0.5 * delta_phi,
                       sqrt_numCell)

    cdef:
        double[:,::1] cellAreas = np.zeros((sqrt_numCell, sqrt_numCell), dtype = np.double)
        double[::1] areaInterpArray = np.zeros(1000, dtype = np.double)
        double[::1] maxEmissionAngle = np.zeros(sqrt_numCell, dtype = np.double)
        double[::1] cellRadialCoord = np.zeros(sqrt_numCell, dtype = np.double)
        double[::1] cos_gamma = np.zeros(sqrt_numCell, dtype = np.double)
        double[::1] effGrav = np.zeros(sqrt_numCell, dtype = np.double)
        double[::1] colatInterpArray
        double[::1] parallels = np.zeros(sqrt_numCell - 1, dtype = np.double)
        double[::1] cellColatitudes = np.zeros(sqrt_numCell, dtype = np.double)

    colatInterpArray = np.linspace(colat_lims[1],
                                   colat_lims[0],
                                   1000)

    cellArea /= (<double> numCell)
    eta = cellArea / delta_phi

    for ii in prange(1000,
                     nogil = True,
                     schedule = 'static',
                     num_threads = numThreads,
                     chunksize = 1):
        i = <size_t> ii
        thread = threadid()
        areaInterpArray[i] = integrateArea(colatInterpArray[i],
                                           colat_lims[1],
                                           R_eq,
                                           epsilon,
                                           zeta,
                                           0,
                                           w[thread])

    # Compute cell boundary parallels
    cdef double i_eta = (<double>sqrt_numCell) * eta
    cdef double *area_ptr = &areaInterpArray[0]
    cdef double *colat_ptr = &colatInterpArray[0]
    gsl_interp_init(interp, area_ptr, colat_ptr, 1000)

    for i in range(sqrt_numCell - 1):
        i_eta -= eta
        parallels[i] = gsl_interp_eval(interp, area_ptr, colat_ptr,
                                       i_eta, accel)

    gsl_interp_free(interp)
    gsl_interp_accel_free(accel)

    # cell vertices
    cdef double lower, upper, leftmost, rightmost
    cdef int p1, p2, p3, p4, p5
    cdef size_t j_max

    # superseding region boundary points where tangent to iso-coordinate curve
    cdef double super_1[2]
    cdef double super_2[2]
    cdef double super_3[2]
    cdef double super_4[2]
    cdef int s1, s2, s3, s4

    if superRadius > 0.0:
        if superColatitude - superRadius < 0.0:
            super_1[0] = superRadius - superColatitude
            super_1[1] = superAzimuth + M_PI
            if super_1[1] > M_PI:
                super_1[1] = super_1[1] - M_2PI

            super_2[0] = superColatitude + superRadius
            super_2[1] = superAzimuth

            # skip point checks if encompasses pole
            super_3[0] = 0.0
            super_3[1] = 0.0
            super_4[0] = 0.0
            super_4[1] = 0.0
        elif superColatitude + superRadius > M_PI:
            super_1[0] = superColatitude - superRadius
            super_1[1] = superAzimuth

            super_2[0] = M_2PI - superColatitude - superRadius
            super_2[1] = superAzimuth + M_PI
            if super_2[1] > M_PI:
                super_2[1] = super_2[1] - M_2PI

            # skip point checks if encompasses pole
            super_3[0] = 0.0
            super_3[1] = 0.0
            super_4[0] = 0.0
            super_4[1] = 0.0
        else:
            super_1[0] = superColatitude - superRadius
            super_1[1] = superAzimuth

            # right-spherical triangle identity; principal solution is correct
            super_2[0] = acos(cos(superColatitude)/cos(superRadius))
            super_2[1] = superAzimuth + asin(sin(superRadius) / sin(superColatitude))
            if super_2[1] > M_PI:
                super_2[1] = super_2[1] - M_2PI

            super_3[0] = superColatitude + superRadius
            super_3[1] = superAzimuth

            # symmetry
            super_4[0] = super_2[0]
            super_4[1] = superAzimuth - asin(sin(superRadius) / sin(superColatitude))
            if super_4[1] < -M_PI:
                super_4[1] = super_4[1] + M_2PI

    # exploit symmetry if no superseding region offset
    if superColatitude == cedeColatitude and superAzimuth == 0.0:
        j_max = <size_t>(sqrt_numCell/2)
    else:
        j_max = sqrt_numCell

    # Compute cell centre colatitudes
    for ii in prange(<signed int>sqrt_numCell,
                     nogil = True,
                     schedule = 'static',
                     num_threads = numThreads,
                     chunksize = 1):
        thread = threadid()
        i = <size_t> ii
        if i == sqrt_numCell - 1:
            lower = parallels[i-1]
            upper = colat_lims[1]
        elif i == 0:
            lower = colat_lims[0]
            upper = parallels[i]
        else:
            lower = parallels[i - 1]
            upper = parallels[i]

        cellColatitudes[i] = integrateArea(lower,
                                           upper,
                                           R_eq,
                                           epsilon,
                                           zeta,
                                           1,
                                           w[thread]) / eta

        mu = cos(cellColatitudes[i])
        radius = radiusNormalised(mu, epsilon, zeta)
        cellRadialCoord[i] = radius * R_eq
        #r_s_over_r = r_s / cellRadialCoord[i]
        effGrav[i] = effectiveGravity(mu, R_eq, zeta, epsilon)
        f = f_theta(mu, radius, epsilon, zeta)
        cos_gamma[i] = 1.0 / sqrt(1.0 + f * f)
        maxEmissionAngle[i] = _hpi + acos(cos_gamma[i])

        leftmost = -boundary_phi
        for j in range(j_max):
            rightmost = leftmost + delta_phi
            p1 = (eval_psi(lower, leftmost, cedeColatitude) <= cedeRadius)
            p2 = (eval_psi(upper, leftmost, cedeColatitude) <= cedeRadius)
            p3 = (eval_psi(lower, rightmost, cedeColatitude) <= cedeRadius)
            p4 = (eval_psi(upper, rightmost, cedeColatitude) <= cedeRadius)
            p5 = (eval_psi(cellColatitudes[i], leftmost+0.5*delta_phi, cedeColatitude) <= cedeRadius)

            if p1 == 1:
                p1 = p1 * (superRadius <= eval_psi(lower, leftmost - superAzimuth, superColatitude))
            else:
                p1 = 2
            if p2 == 1:
                p2 = p2 * (superRadius <= eval_psi(upper, leftmost - superAzimuth, superColatitude))
            else:
                p2 = 2
            if p3 == 1:
                p3 = p3 * (superRadius <= eval_psi(lower, rightmost - superAzimuth, superColatitude))
            else:
                p3 = 2
            if p4 == 1:
                p4 = p4 * (superRadius <= eval_psi(upper, rightmost - superAzimuth, superColatitude))
            else:
                p4 = 2
            if p5 == 1:
                p5 = p5 * (superRadius <= eval_psi(cellColatitudes[i], leftmost+0.5*delta_phi - superAzimuth, superColatitude))
            else:
                p5 = 2

            #with gil:
                #print("i=%i, j=%i, p1=%i, p2=%i, p3=%i, p4=%i, p5=%i" % (i,j,p1,p2,p3,p4,p5))

            if not (p1 == p2 == p3 == p4 == p5):
                cellAreas[i,j] = integrateCell(lower,
                                               upper,
                                               leftmost,
                                               rightmost,
                                               R_eq,
                                               epsilon,
                                               zeta,
                                               cedeColatitude,
                                               cedeRadius,
                                               superRadius,
                                               superAzimuth,
                                               superColatitude,
                                               w[thread])
                #if cellAreas[i,j] == 0.0:
                #    with gil:
                #        print('i=%d, j=%d, cellAreas[i,j]=%.8e' % (i,j, cellAreas[i,j]))
            elif superRadius > 0.0:
                if (lower <= super_1[0] <= upper) and (leftmost <= super_1[1] <= rightmost):
                    s1 = 1
                else:
                    s1 = 0

                if (lower <= super_2[0] <= upper) and (leftmost <= super_2[1] <= rightmost):
                    s2 = 1
                else:
                    s2 = 0

                if (lower <= super_3[0] <= upper) and (leftmost <= super_3[1] <= rightmost):
                    s3 = 1
                else:
                    s3 = 0

                if (lower <= super_4[0] <= upper) and (leftmost <= super_4[1] <= rightmost):
                    s4 = 1
                else:
                    s4 = 0

                if s1==1 or s2==1 or s3==1 or s4==1:
                    cellAreas[i,j] = integrateCell(lower,
                                                   upper,
                                                   leftmost,
                                                   rightmost,
                                                   R_eq,
                                                   epsilon,
                                                   zeta,
                                                   cedeColatitude,
                                                   cedeRadius,
                                                   superRadius,
                                                   superAzimuth,
                                                   superColatitude,
                                                   w[thread])

                elif p5 == 1:
                    cellAreas[i,j] = cellArea
            elif p5 == 1:
                cellAreas[i,j] = cellArea

            if superColatitude == cedeColatitude and superAzimuth == 0.0:
                cellAreas[i,sqrt_numCell - 1 - j] = cellAreas[i,j]

            leftmost = rightmost

    phi, theta = np.meshgrid(phi,
                             np.asarray(cellColatitudes, dtype = np.double, order='C'))

    #cdef double total, superRegion
    ## testing
    #total = integrateSpot(super_colat_lims[0],
    #                         supecolat_lims[1],
    #                         -boundary_phi,
    #                         boundary_phi,
    #                         R_eq,
    #                         epsilon,
    #                         zeta,
    #                         colatitude,
    #                         cedeRadius,
    #                         w1[0],
    #                         w2[0])

    #print("\n(Total) %.16e" % total)

    #superRegion = integrateSpot(colat_lims[0],
    #                         colat_lims[1],
    #                         -boundary_phi,
    #                         boundary_phi,
    #                         R_eq,
    #                         epsilon,
    #                         zeta,
    #                         colatitude,
    #                         superRadius,
    #                         w1[0],
    #                         w2[0])

    #print("(SuperRegion) %.16e" % superRegion)

    #spotArea = integrateCell(colat_lims[0],
    #                         colat_lims[1],
    #                         -boundary_phi,
    #                         boundary_phi,
    #                         R_eq,
    #                         epsilon,
    #                         zeta,
    #                         colatitude,
    #                         superRadius,
    #                         cedeRadius,
    #                         w1[0],
    #                         w2[0])

    #print("(CedeRegion direct) %.16e" % spotArea)
    #print("(Ratio) %.16e" % ((total - superRegion)/spotArea))

    #cdef double totalCellArea
    #totalCellArea = 0.0
    #for i in range(sqrt_numCell):
    #   for j in range(sqrt_numCell):
    #       totalCellArea += cellAreas[i,j]
    #print("(Spot) spotArea / totalCellArea: %.16e" % ((total - superRegion) / totalCellArea))
    ## end testing

    # testing
    #spotArea = integrateSpot(colat_lims[0],
    #                         colat_lims[1],
    #                         R_eq,
    #                         epsilon,
    #                         zeta,
    #                         cedeColatitude,
    #                         cedeRadius,
    #                         superRadius,
    #                         superColatitude,
    #                         superAzimuth,
    #                         0,
    #                         0.0,
    #                         w[0])

    #cdef double totalCellArea
    #totalCellArea = 0.0
    #for i in range(sqrt_numCell):
    #   for j in range(sqrt_numCell):
    #       totalCellArea += cellAreas[i,j]
    #print("(Spot) spotArea / totalCellArea: %.16e" % (spotArea / totalCellArea))
    # end testing


    for thread in range(numThreads):
        gsl_integration_cquad_workspace_free(w[thread])

    free(w)

    return (theta,
            phi,
            np.asarray(cellRadialCoord, dtype = np.double, order='C'),
            np.asarray(cellAreas, dtype = np.double, order='C'),
            np.asarray(maxEmissionAngle, dtype = np.double, order='C'),
            np.asarray(cos_gamma, dtype = np.double, order='C'),
            np.asarray(effGrav, dtype = np.double, order='C'))
