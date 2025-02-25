from GSL cimport (gsl_interp_type,
                  gsl_interp_akima_periodic,
                  gsl_interp_akima,
                  gsl_interp_steffen,
                  gsl_interp_cspline_periodic,
                  gsl_interp_cspline)

cdef gsl_interp_type* _get_phase_interpolant() except *

cdef gsl_interp_type* _get_energy_interpolant() except *
