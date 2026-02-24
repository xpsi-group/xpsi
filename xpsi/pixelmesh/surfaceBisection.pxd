from .geometricConfiguration cimport _GEOM

cdef double RADIUS(double x, double epsilon, double zeta, int obl_surfgrav_ind) noexcept nogil

cdef int BISECT(const double *const y_p,
                const double *const y,
                double *const y_new,
                double affine_p,
                double affine,
                double a,
                double R_eq,
                double epsilon,
                double zeta,
                int obl_surfgrav_ind) noexcept nogil

cdef double ZABB(const _GEOM *const GEOM,
                 const double *const y,
                 double b,
                 double *const Z,
                 double *const ABB, int obl_surfgrav_ind) noexcept nogil
