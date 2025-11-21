from GSL cimport (gsl_interp_type,
                  gsl_interp_akima_periodic,
                  gsl_interp_akima,
                  gsl_interp_steffen,
                  gsl_interp_linear,
                  gsl_interp_cspline_periodic,
                  gsl_interp_cspline)

cdef const gsl_interp_type* _get_phase_interpolant() except *

cdef const gsl_interp_type* _get_energy_interpolant() except *

cdef bint are_equal(double x, double y, double epsilon=*) noexcept nogil
