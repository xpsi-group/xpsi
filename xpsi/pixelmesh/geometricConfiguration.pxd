ctypedef struct _GEOM:
    double mass
    double r_s
    double R_eq       # equatorial CIRCUMFERENTIAL radius (physical input)
    double R_eq_sph   # pixelmesh only: equatorial radius of the surface in the
                      # spherical coordinate r_sph built from Boyer-Lindquist r,
                      # chosen so that the circumferential equatorial radius
                      # equals R_eq (set by surfaceBisection.SURFACE_SCALE)
    double omega
    double mode_omega
    double epsilon
    double zeta
    double a
    double asq
    double kappa
    double inclination
    double d
    double b_max
    int obl_surfgrav_ind
