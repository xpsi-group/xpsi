from GSL cimport (gsl_function,
                   gsl_integration_cquad,
                   gsl_integration_cquad_workspace,
                   gsl_integration_cquad_workspace_alloc,
                   gsl_integration_cquad_workspace_free,
                   gsl_interp_eval,
                   gsl_integration_cquad_workspace,
                   gsl_interp_alloc,
                   gsl_interp_accel,
                   gsl_interp_accel_alloc,
                   gsl_interp,
                   gsl_interp_linear,
                   gsl_interp_steffen,
                   gsl_interp_init,
                   gsl_interp_free,
                   gsl_interp_accel_free,
                   gsl_interp_accel_reset,
                   gsl_isnan)

ctypedef gsl_integration_cquad_workspace gsl_cq_work

from xpsi.surface_radiation_field.effective_gravity_universal cimport effectiveGravity

cdef double eval_psi(double theta, double phi, double THETA) nogil
cdef double eval_phi(double theta, double THETA, double psi) nogil
cdef double eval_cedeAzi(double theta, double phi, double psi, double THETA) nogil

cdef double radiusNormalised(double mu,
                             double epsilon,
                             double zeta) nogil
cdef double f_theta(double mu,
                    double radiusNormed,
                    double epsilon,
                    double zeta) nogil

cdef double integrateArea(double lower_lim,
                          double upper_lim,
                          double R_eq,
                          double epsilon,
                          double zeta,
                          int av,
                          gsl_integration_cquad_workspace *w) nogil

cdef double integrateCell(double theta_a,
                          double theta_b,
                          double phi_a,
                          double phi_b,
                          double R_eq,
                          double epsilon,
                          double zeta,
                          double colat,
                          double cedeRadius,
                          double superRadius,
                          double superCentreAzi,
                          double superCentreColat,
                          gsl_cq_work *w) nogil

cdef double integrateSpot(double theta_a,
                          double theta_b,
                          double R_eq,
                          double epsilon,
                          double zeta,
                          double colat,
                          double cedeRadius,
                          double superRadius,
                          double superCentreAzi,
                          double superCentreColat,
                          int Lorentz,
                          double Omega,
                          gsl_cq_work *w) nogil
