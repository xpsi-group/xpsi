from xpsi.pixelmesh.RK_IP2S_tracer cimport _RAY
from .geometricConfiguration cimport _GEOM

ctypedef struct RAY_MAP:
    size_t numRays
    double *r_MESH
    double *t_MESH
    double *X_MESH
    double *X_MODDED
    double *Y_MESH
    double *RADIAL
    double *THETA
    double *PHI
    double *LAG
    double *ABB
    double *Z
    double *XI
    double ORIGIN[6]
    double refRayTime
    double SEMI_MAJOR
    double SEMI_MINOR
    double ORIGIN_X
    double ORIGIN_Y

cdef RAY_MAP* alloc_RAY_MAP(size_t NGR) nogil

cdef void free_RAY_MAP(RAY_MAP *const MAP) nogil

cdef void init_RAY_MAP(RAY_MAP *const MAP,
                       double radialIncrementExponent) nogil

cdef int compute_globalRayMap(size_t numThreads,
                              double radialIncrementExponent,
                              const _GEOM *const GEOM,
                              RAY_MAP *const MAP,
                              _RAY **RAYS,
                              int force_circular) nogil
