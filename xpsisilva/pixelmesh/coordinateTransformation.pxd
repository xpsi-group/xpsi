from .globalRayMap cimport RAY_MAP
from .geometricConfiguration cimport _GEOM
#from focusRayMap cimport FOCUS_RAY_MAP

cdef void BOYERLINDQUIST_2_SPHERICAL(RAY_MAP *const MAP,
                                     const _GEOM *const GEOM) nogil

#cdef void BOYERLINDQUIST_2_SPHERICAL_FOCUS(FOCUS_RAY_MAP *const F,
#                                           const _GEOM *const GEOM) nogil
#
#cdef void BOYERLINDQUIST_2_SPHERICAL_DENSE(FOCUS_RAY_MAP *const F,
#                                           const _GEOM *const GEOM) nogil
