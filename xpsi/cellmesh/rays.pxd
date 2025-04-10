cdef extern from "rayXpanda/inversion.h":
    void __pyx_f_9rayXpanda_9inversion_c_invert(double, double, double *, double *) noexcept nogil

cdef extern from "rayXpanda/deflection.h":
    void __pyx_f_9rayXpanda_10deflection_c_deflect(double, double, double *, double *) noexcept nogil

cdef void invert(double a, double b, double *c, double *d) noexcept nogil
cdef void deflect(double a, double b, double *c, double *d) noexcept nogil

cdef double eval_image_deflection(int order, double psi) noexcept nogil

cdef void link_rayXpanda(bint *use_rayXpanda, double *rayXpanda_defl_lim) except *
