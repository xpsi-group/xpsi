# Gravity function pointer type
ctypedef double (*grav_fn)(double mu,
                           double R_eq,
                           double x,
                           double epsilon) noexcept nogil

# Selector (called once with GIL)
cdef grav_fn select_gravity_fn(int model_id)

# Compatibility wrapper (used elsewhere in X-PSI)
cdef double effectiveGravity(double mu,
                             double R_eq,
                             double x,
                             double epsilon,
                             int model_id) noexcept nogil



cdef  double oblateness_func_o2(double x, int obl_surfgrav_ind) noexcept nogil
cdef  double dimless_moment_of_inertia_i(double x, int obl_surfgrav_ind) noexcept nogil
cdef  double beta1_coeff(int obl_surfgrav_ind) noexcept nogil

cpdef double py_dimless_moment_of_inertia_i(double x, int obl_surfgrav_ind)
cpdef double py_oblateness_func_o2(double x, int obl_surfgrav_ind)
cpdef double py_beta1_coeff(int obl_surfgrav_ind)
