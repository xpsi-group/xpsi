from .geometricConfiguration cimport _GEOM

cdef double RADIUS(double x, double epsilon, double zeta) noexcept nogil
    
cdef int BISECT(const double *const y_p,
                const double *const y,
                double *const y_new,
                double affine_p,
                double affine,
                double a,
                double R_eq,
                double epsilon,
                double zeta) noexcept nogil

cdef double ZABB(const _GEOM *const GEOM,
                 const double *const y,
                 double b,
                 double *const Z,
                 double *const ABB) noexcept nogil
