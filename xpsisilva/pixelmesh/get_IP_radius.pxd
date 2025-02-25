from xpsisilva.pixelmesh.RK_IP2S_tracer cimport _RAY
from .geometricConfiguration cimport _GEOM
from xpsisilva.pixelmesh.globalRayMap cimport RAY_MAP

cdef double compute_imagePlane_radius(const _GEOM *const GEOM,
                                      _RAY *const RAY,
                                      RAY_MAP *const MAP,
				      int force_circular) nogil
