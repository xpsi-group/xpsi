from GSL cimport (gsl_interp_eval,
                   gsl_interp_eval_deriv,
                   gsl_interp_alloc,
                   gsl_interp_accel,
                   gsl_interp_accel_alloc,
                   gsl_interp_steffen,
                   gsl_interp_linear,
                   gsl_interp,
                   gsl_interp_init,
                   gsl_interp_free,
                   gsl_interp_accel_free,
                   gsl_interp_accel_reset,
                   gsl_isnan,
                   gsl_isinf)

ctypedef gsl_interp_accel accel
ctypedef gsl_interp interp
