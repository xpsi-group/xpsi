cdef double compute_pol_ang(
              double leaves_kdx,
              double sin_psi,
              double cos_psi,
              double sin_alpha,
              double cos_alpha,
              double sin_theta_i,
              double cos_theta_i,
              double sin_i,
              double cos_i,
              double sin_gamma,
              double cos_gamma,
              double Grav_z,
              double mu,
              double eta,
              double beta,
              double Lorentz,
              double cos_xi) noexcept nogil
              
cdef int disk_block(
              double R_in,
              double cos_i,
              double cos_psi,
              double cos_theta_i,
              double r_s_over_r,
              double radius,
              double sin_alpha,
              double theta_i_over_pi
              ) noexcept nogil             
              
