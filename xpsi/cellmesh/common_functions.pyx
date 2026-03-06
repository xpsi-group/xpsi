from libc.math cimport M_PI, sqrt, sin, cos, acos, log10, pow, exp, fabs, ceil, log, atan2
from ..tools.core cimport are_equal
from libc.stdio cimport printf

cdef double _2pi = 2.0 * M_PI

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
              double cos_xi) noexcept nogil:
    """
    Calculation of polarization angle using the polarized Oblate Schwarcshild approximation as in Loktev et al. (2020).
    """
     
    cdef:         
        double sin_chi_0, cos_chi_0, chi_0, chi_1, chi_prime, chi
        double sin_chi_1, cos_chi_1, sin_chi_prime, cos_chi_prime
        double sin_lambda, cos_lambda, cos_eps
        double sin_alpha_over_sin_psi        

    sin_chi_0 = - sin_theta_i*sin(leaves_kdx) 
    cos_chi_0 = sin_i*cos_theta_i - sin_theta_i*cos_i*cos(leaves_kdx)
    chi_0 = atan2(sin_chi_0,cos_chi_0)

    if not are_equal(sin_psi, 0.0):
        sin_alpha_over_sin_psi = sin_alpha/sin_psi
    else: #using small-angle limit of the Beloborodov (2002) approximation
        sin_alpha_over_sin_psi = Grav_z

    #Notes: mu = cos_sigma , Lorentz = 1/Gamma, mu0=eta*mu, cos_xi defined with no minus sign
    sin_chi_1 = sin_gamma*sin_i*sin(leaves_kdx)*sin_alpha_over_sin_psi #times sin alpha sin sigma
    cos_chi_1 = cos_gamma - cos_alpha*mu  #times sin alpha sin sigma 
    chi_1 = atan2(sin_chi_1,cos_chi_1)

    sin_lambda = sin_theta_i*cos_gamma - sin_gamma*cos_theta_i
    cos_lambda = cos_theta_i*cos_gamma + sin_theta_i*sin_gamma
    cos_eps = sin_alpha_over_sin_psi*(cos_i*sin_lambda - sin_i*cos_lambda*cos(leaves_kdx) + cos_psi*sin_gamma) - cos_alpha*sin_gamma

    sin_chi_prime = cos_eps*eta*mu*beta/Lorentz
    cos_chi_prime = (1. - mu**2 /(1. + beta*cos_xi))
    chi_prime = atan2(sin_chi_prime,cos_chi_prime)

    chi = chi_0+chi_1+chi_prime

    #printf("leaves[_kdx] = %.6e ",leaves_kdx/_2pi)
    #printf("chi_0 = %.6e ",chi_0)
    #printf("chi_1 = %.6e ",chi_1)
    #printf("chi_prime = %.6e ",chi_prime)
    #printf("PA_tot = %.6e\n",chi)

    return chi

# Python wrapper for testing
cpdef double compute_pol_ang_py(
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
              double cos_xi):
    """
    Python-callable wrapper for compute_pol_ang.
    Useful for unit testing.
    """
    return compute_pol_ang(
        leaves_kdx,
        sin_psi,
        cos_psi,
        sin_alpha,
        cos_alpha,
        sin_theta_i,
        cos_theta_i,
        sin_i,
        cos_i,
        sin_gamma,
        cos_gamma,
        Grav_z,
        mu,
        eta,
        beta,
        Lorentz,
        cos_xi
    )

cdef int disk_block(
              double R_in,
              double cos_i,
              double cos_psi,
              double cos_theta_i,
              double r_s_over_r_i,
              double radius,
              double sin_alpha,
              double theta_i_over_pi
              ) noexcept nogil:
    """
    Checking whether an accretion disk blocks certain rays based on Ibragimov & Poutanen (2009). Returns 1 if the ray is not blocked.
    """
     
    cdef:         
        double cos_psi_d, sin_psi_d      # geometric quantities
        double impact_b, r_psi_d         # impact parameter and radial coordinate
        double r_s_i                     # schwarzschild radius at the cell location  

    cos_psi_d = (cos_i * cos_psi - cos_theta_i) / sqrt(cos_i * cos_i + cos_theta_i * cos_theta_i - 2 * cos_i * cos_theta_i * cos_psi) #Ibragimov & Poutanen (2009), Equation (C2)
    sin_psi_d = sqrt(1 - cos_psi_d * cos_psi_d)
    r_s_i = r_s_over_r_i*radius
    impact_b = radius * sin_alpha / sqrt(1 - r_s_over_r_i) # impact parameter
    r_psi_d = sqrt((r_s_i * r_s_i * (1 - cos_psi_d) * (1 - cos_psi_d)) / (4 * (1 + cos_psi_d) * (1 + cos_psi_d)) +  ((impact_b * impact_b) / (sin_psi_d * sin_psi_d))) - (r_s_i * (1 - cos_psi_d)) / (2 * (1 + cos_psi_d)) ##Ibragimov & Poutanen (2009), Equation (B9)
    printf("r_psi_d = %.6e\n",r_psi_d)
    if theta_i_over_pi < 0.5 or (theta_i_over_pi > 0.5 and r_psi_d < R_in):
        return 1
    else: #theta > pi/2 and (theta < pi/2 or r_psi_d > R_in), don't calculate.
        return 0

# Python wrapper for testing
cpdef int disk_block_py(
              double R_in,
              double cos_i,
              double cos_psi,
              double cos_theta_i,
              double r_s_over_r_i,
              double radius,
              double sin_alpha,
              double theta_i_over_pi
              ):
    """
    Python-callable wrapper for disk_ block.
    Useful for unit testing.
    """
    return disk_block(
            R_in,
            cos_i,
            cos_psi,
            cos_theta_i,
            r_s_over_r_i,
            radius,
            sin_alpha,
            theta_i_over_pi
    )

