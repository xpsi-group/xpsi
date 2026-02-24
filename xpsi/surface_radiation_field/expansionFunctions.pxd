cdef  double oblateness_func_o2(double x, int obl_surfgrav_ind) noexcept nogil
cdef  double dimless_moment_of_inertia_i(double x, int obl_surfgrav_ind) noexcept nogil
cdef  double beta1_coeff(int obl_surfgrav_ind) noexcept nogil

cpdef double py_dimless_moment_of_inertia_i(double x, int obl_surfgrav_ind)
cpdef double py_oblateness_func_o2(double x, int obl_surfgrav_ind)
cpdef double py_beta1_coeff(int obl_surfgrav_ind)
