#cython: cdivision=True
#cython: boundscheck=False
#cython: nonecheck=False
#cython: wraparound=False

from libc.math cimport sqrt, sin, cos, acos, pow

cdef void BOYERLINDQUIST_2_SPHERICAL(RAY_MAP *const MAP,
                                     const _GEOM *const GEOM) nogil:
    cdef:
        size_t i, j
        double sin_theta_sq
        size_t N = MAP.numRays
        size_t INDEX = 0
        double rsq

    for i in range(N):
        for j in range(N):
            if MAP.THETA[INDEX] >= 0.0:

                rsq = pow(MAP.RADIAL[INDEX], 2.0)
                sin_theta_sq = pow(sin(MAP.THETA[INDEX]), 2.0)

                MAP.THETA[INDEX] = acos(cos(MAP.THETA[INDEX]) / sqrt(1.0 + GEOM.asq * sin_theta_sq / rsq))
                MAP.RADIAL[INDEX] = sqrt(rsq + GEOM.asq * sin_theta_sq)
            INDEX += 1

    if MAP.ORIGIN[3] >= 0.0:

        rsq = MAP.ORIGIN[2] * MAP.ORIGIN[2]
        sin_theta_sq = pow(sin(MAP.ORIGIN[3]), 2.0)

        MAP.ORIGIN[3] = acos(cos(MAP.ORIGIN[3]) / sqrt(1.0 + GEOM.asq * sin_theta_sq / rsq))
        MAP.ORIGIN[2] = sqrt(rsq + GEOM.asq * sin_theta_sq)

#cdef void BOYERLINDQUIST_2_SPHERICAL_FOCUS(FOCUS_RAY_MAP *const F,
#                                           const _GEOM *const GEOM) nogil:
#    cdef:
#        size_t i, j
#        double sin_theta_sq
#        size_t N = F.numFocusRays
#        size_t INDEX = 0
#        double rsq
#
#    for i in range(N):
#        for j in range(N):
#            if F.F_RADIAL[INDEX] >= 0.0:
#
#                rsq = pow(F.F_RADIAL[INDEX], 2.0)
#                sin_theta_sq = pow(sin(F.F_THETA[INDEX]), 2.0)
#
#                F.F_THETA[INDEX] = acos(cos(F.F_THETA[INDEX]) / sqrt(1.0 + GEOM.asq * sin_theta_sq / rsq))
#                F.F_RADIAL[INDEX] = sqrt(rsq + GEOM.asq * sin_theta_sq)
#            INDEX += 1
#
#    if F.ORIGIN[3] >= 0.0:
#        rsq = F.ORIGIN[2] * F.ORIGIN[2]
#        sin_theta_sq = pow(sin(F.ORIGIN[3]), 2.0)
#
#        F.ORIGIN[3] = acos(cos(F.ORIGIN[3]) / sqrt(1.0 + GEOM.asq * sin_theta_sq / rsq))
#        F.ORIGIN[2] = sqrt(rsq + GEOM.asq * sin_theta_sq)
#
#cdef void BOYERLINDQUIST_2_SPHERICAL_DENSE(FOCUS_RAY_MAP *const F,
#                                           const _GEOM *const GEOM) nogil:
#    cdef:
#        size_t i, j
#        double sin_theta_sq
#        size_t N = F.numDenseRays
#        size_t INDEX = 0
#        double rsq
#
#    for i in range(N):
#        for j in range(N):
#            if F.RADIAL[INDEX] >= 0.0:
#
#                rsq = pow(F.RADIAL[INDEX], 2.0)
#                sin_theta_sq = pow(sin(F.THETA[INDEX]), 2.0)
#
#                F.THETA[INDEX] = acos(cos(F.THETA[INDEX]) / sqrt(1.0 + GEOM.asq * sin_theta_sq / rsq))
#                F.RADIAL[INDEX] = sqrt(rsq + GEOM.asq * sin_theta_sq)
#            INDEX += 1
#
#    if F.ORIGIN[3] >= 0.0:
#        rsq = F.ORIGIN[2] * F.ORIGIN[2]
#        sin_theta_sq = pow(sin(F.ORIGIN[3]), 2.0)
#
#        F.ORIGIN[3] = acos(cos(F.ORIGIN[3]) / sqrt(1.0 + GEOM.asq * sin_theta_sq / rsq))
#        F.ORIGIN[2] = sqrt(rsq + GEOM.asq * sin_theta_sq)






