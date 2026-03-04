cimport numpy as np

from GSL cimport (gsl_interp,
                   gsl_interp_accel)
ctypedef np.uint8_t uint8

cdef uint8[::1] build_allow_negative_array( allow_negative , 
                                            size_t num_components)

cdef void compute_expected_star_count_rate_single_channel(double[::1] phases,
                                         size_t num_components,
                                         object components,
                                         object component_phases,
                                         double[::1] phase_shifts,
                                         gsl_interp **interp,
                                         gsl_interp_accel **acc,
                                         double *STAR,
                                         int channel,
                                         uint8[::1] allow_negative)

cdef double[:,::1] compute_expected_counts(double exposure_time,
                            double[::1] phases,
                            object components,
                            object component_phases,
                            double[::1] phase_shifts,
                            double[:,::1] background,
                            object allow_negative)