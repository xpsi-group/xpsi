from .geometricConfiguration cimport _GEOM
from xpsisilva.pixelmesh.RK_IP2S_tracer cimport _RAY

cdef void COMPUTE_BCs(_RAY *const RAY,
                      const _GEOM *const GEOM) nogil
