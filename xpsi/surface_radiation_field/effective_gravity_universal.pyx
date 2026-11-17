#cython: cdivision=True
#cython: boundscheck=False
#cython: nonecheck=False
#cython: wraparound=False

from libc.math cimport sqrt, log10, fabs, NAN

cdef double c = 2.99792458e8

cdef inline double g0(double R_eq, double x) noexcept nogil:
    """
    Compute the relativistic surface gravity normalization g0.

    This corresponds to:
        g0 = (x c²) / [R_eq √(1 − 2x)]

    where:
      - x = GM / (R_eq c²) is the compactness
      - R_eq is the equatorial radius
      - c is the speed of light

    This is the Schwarzschild surface gravity for a non-rotating star,
    used as the baseline for all rotating corrections.

    Parameters
    ----------
    R_eq : double
        Equatorial radius of the star.
    x : double
        Compactness GM / (R_eq c²).

    Returns
    -------
    double
        The surface gravity normalization g₀.
    """
    return x * c * c / (R_eq * sqrt(1.0 - 2.0 * x))

cdef double grav_sphere(double mu, double R_eq, double x, double epsilon) noexcept nogil:
    """
    Effective surface gravity for a non-rotating (spherical) neutron star.

    This returns the logarithm of the surface gravity assuming a perfectly
    spherical star, i.e. ignoring rotational corrections and oblateness.

    The gravity is given by:
        log10(g0) + 2

    where:
        g0 = (x c²) / [R_eq √(1 − 2x)]

    is the Schwarzschild surface gravity normalization, and the +2 converts
    from SI units (m s⁻²) to cgs (cm s⁻²), as required by tabulated
    atmosphere models.
    Parameters
    ----------
    mu : double
        Cosine of the colatitude (cos θ). Unused here since the star is spherical,
        but kept for API consistency with rotating models.
    R_eq : double
        Equatorial radius of the star.
    x : double
        Compactness GM / (R_eq c²).
    epsilon : double
        Dimensionless spin parameter. Ignored in this model.

    Returns
    -------
    double
        log10 of the effective surface gravity in cgs units.
    """
    # For a spherical star, the surface gravity is constant over the surface,
    # so no dependence on mu or epsilon enters.

    return log10(g0(R_eq, x)) + 2.0

cdef double grav_algendy_morsink_2014(double mu, double R_eq, double x, double epsilon) noexcept nogil:
    """
    Calculate the effective surface gravity log-likelihood value based on the approximation from
    AlGendy & Morsink (2014) (see Eq. 21 and Table 5 for coefficient values), or assuming a spherical
    star.

    This function computes the logarithmic surface gravity (log10(g / g0)) for a rotating neutron star,
    incorporating the effects of spin and compactness. The result is expressed in units compatible with
    tabulated atmosphere models.

    :param double mu:
        The cosine of the colatitude angle (theta). This parameter defines the angular
        position on the star's surface, where mu is between -1 and 1.

    :param double R_eq:
        The equatorial radius of the star.

    :param double x:
        The compactness parameter, defined as (G * M) / (R_eq * c**2) (see Eq. 1 or Eq. 9 from
        Morsink et al. (2007)).

    :param double epsilon:
        The dimensionless spin parameter, defined as (omega**2 * R_eq**3) / (G * M) (see Eq. 2 or
        Eq. 10 from Morsink et al. (2007)).

    :return:
        The effective surface gravities in log10(g*g0). The result is scaled to centimeters (conversion
        factor 2.0 is applied) to match the format used in tabulated atmosphere models.
    :rtype: double

    """

    cdef:
        double g_0 = g0(R_eq, x)
        double esq = epsilon
        double esqsq = epsilon * epsilon

        double c_e = -0.791 + 0.776 * x
        double c_p = 1.138 - 1.431 * x
        double d_e = (-1.315 + 2.431 * x) * esq * x
        double d_p = (0.653 - 2.864 * x) * esq * x
        double d_60 = (13.47 - 27.13 * x) * esq * x
        double f_e = -1.172 * x * esqsq
        double f_p = 0.975 * x * esqsq
        double g = 1.0

    g += (c_e + d_e + f_e) * esq * (1.0 - mu * mu)
    g += (c_p + d_p + f_p - d_60) * esq * mu * mu
    g += d_60 * esq * fabs(mu)

    return log10(g * g_0) + 2.0


cdef double grav_j1231(double mu, double R_eq, double x, double epsilon) noexcept nogil:

    """
    Compute the effective surface gravity for PSR J1231-1411 using
     tailored approximation.

    The functional form follows AlGendy & Morsink (2014), Eq. 21, but
    uses best-fit coefficients specifically calibrated for PSR J1231-1411.
    These coefficients were obtained using XXX (ask Lami).

    The approximation incorporates rotational corrections up to sixth
    order in the dimensionless angular velocity and accounts for the
    dependence on stellar compactness and colatitude.

    Parameters
    ----------
    mu : double
        Cosine of the colatitude (cos theta), ranging from -1 to 1.
    R_eq : double
        Equatorial radius of the star in meters.
    x : double
        Stellar compactness, GM / (R_eq c²).
    epsilon : double
        Dimensionless squared angular velocity,
        Omega² R_eq³ / (GM).

    Returns
    -------
    double
        Base-10 logarithm of the effective surface gravity in cgs units
        (cm s⁻²), as required by tabulated atmosphere models.

    Notes
    -----
    The effective surface gravity is calculated as

        g_eff = g0 * [1 + rotational corrections],

    where g0 is the Schwarzschild surface gravity. The rotational
    corrections depend on compactness, colatitude, and powers of epsilon.
    The conversion from SI to cgs units introduces an additive factor
    of 2 in log10(g_eff).
    """
    cdef:
        double g_0   = g0(R_eq, x)

        # epsilon is Ω̄²
        double esq   = epsilon           # Ω̄²
        double esqsq = epsilon * epsilon # Ω̄⁴

        # ---- best-fit coefficients (hard-coded constants) ----
        double C0E  = -1.06853309
        double C1E  =  2.61536903
        double C0P  =  1.07372218
        double C1P  = -1.00497676

        double D1E  =  21.80711091
        double D2E  = -99.07086300
        double F1E  = -69.84926438

        double D1P  =  -6.11239894
        double D2P  =  10.26276504
        double F1P  =  40.79848048

        double D160 =  76.56620037
        double D260 = -380.74234012

        # ---- build coefficients ----
        double c_e  = C0E + C1E * x
        double c_p  = C0P + C1P * x

        # These are the Ω̄⁴ contributions
        double d_e  = (D1E + D2E * x) * esq * x
        double d_p  = (D1P + D2P * x) * esq * x
        double d_60 = (D160 + D260 * x) * esq * x

        # These are the Ω̄⁶ contributions
        double f_e  = F1E * x * esqsq
        double f_p  = F1P * x * esqsq

        double g = 1.0

    # mu = cos(theta)
    g += (c_e + d_e + f_e) * esq * (1.0 - mu * mu)          # sin²θ term
    g += (c_p + d_p + f_p - d_60) * esq * (mu * mu)         # cos²θ term
    g += d_60 * esq * fabs(mu)                              # cosθ term (symmetrized)

    return log10(g * g_0) + 2.0


cdef double grav_j1614(double mu, double R_eq, double x, double epsilon) noexcept nogil:
    """
    Compute the effective surface gravity for PSR J1614-2230 using
    tailored approximation.

    The surface-gravity prescription is given in:
    https://arxiv.org/pdf/2609.00172

    The approximation accounts for the dependence on stellar compactness,
    colatitude, and rotation through the dimensionless squared angular
    velocity.

    Parameters
    ----------
    mu : double
        Cosine of the colatitude (cos theta), ranging from -1 to 1.
    R_eq : double
        Equatorial radius of the star in meters.
    x : double
        Stellar compactness, GM / (R_eq c²).
    epsilon : double
        Dimensionless squared angular velocity,
        Omega² R_eq³ / (GM).

    Returns
    -------
    double
        Base-10 logarithm of the effective surface gravity in cgs units
        (cm s⁻²), as required by tabulated atmosphere models.

    Notes
    -----
    The effective surface gravity is calculated as

        g_eff = g0 * [0.9985 + epsilon * (-0.5454 + (1.504 + 0.2422*x)*mu² - 0.6679*x*mu⁴)],

    where g0 is the Schwarzschild surface gravity. The conversion from
    SI to cgs units introduces an additive factor of 2 in log10(g_eff).
    """

    cdef:
        double g_0 = g0(R_eq, x)
        double mu2 = mu * mu
        double g = 0.9985 + epsilon * (-0.5454 + (1.504 + 0.2422 * x) * mu2 - 0.6679 * x * mu2 * mu2)
    return log10(g * g_0) + 2.0


cdef grav_fn select_gravity_fn(int obl_surfgrav_ind):
    """
    Select gravity function pointer.

    Must be called with GIL (not inside nogil/prange).
    Raises if model_id is invalid.
    """
    if obl_surfgrav_ind == 0:
        return grav_algendy_morsink_2014
    elif obl_surfgrav_ind == 1:
        return grav_sphere
    elif obl_surfgrav_ind == 5:
        return grav_j1231
    elif obl_surfgrav_ind == 7:
        return grav_j1614

    else:
        raise ValueError("Unknown gravity model_id")


cdef inline double effectiveGravity(double mu,
                                    double R_eq,
                                    double x,
                                    double epsilon,
                                    int obl_surfgrav_ind) noexcept nogil:


    if obl_surfgrav_ind ==0 :
        return grav_algendy_morsink_2014(mu, R_eq, x, epsilon)
    elif obl_surfgrav_ind == 1:
        return grav_sphere(mu, R_eq, x, epsilon)
    elif obl_surfgrav_ind == 5:
        return grav_j1231(mu, R_eq, x, epsilon)
    elif obl_surfgrav_ind == 7:
        return grav_j1614(mu, R_eq, x, epsilon)
    else:
        raise ValueError("Unknown gravity model_id")




cdef double oblateness_func_o2(double x, int obl_surfgrav_ind) noexcept nogil:
    """
    Oblateness coefficient o2(x) entering the surface shape expansion:

        R(θ) = Re * [1 + o2(x) * cos^2(θ)]

    Parameters
    ----------
    x : double
        Compactness ratio x = M / Re.
    model : obl_surfgrav
        0 : AlGendy & Morsink (2014) slow-rotation universal fit
        5   : Empirical fit from LORENE + CUTTER

    Returns
    -------
    double
        o2(x) coefficient (dimensionless).
    """
    cdef double c1, c2, c3
    if obl_surfgrav_ind == 0:
        # Eq. (21): o2 = o20 + o21 x   (overall Ω̄^2 factor is applied elsewhere)
        return -0.788 + 1.030 * x
    elif obl_surfgrav_ind == 1:
        return 1.
    elif obl_surfgrav_ind == 5:
        # Computed by Yves (Kini et al. in prep)
        return -0.37562887294361336 - 4.605853993313363*x + 21.39216085992957*x*x - 24.252054038949396*x*x*x
    elif obl_surfgrav_ind == 7:
        # Computed by Lucien (https://arxiv.org/pdf/2609.00172)
        return 0.615 -2.0*x-0.0539/(x*x)+0.00079/(x*x*x*x)


cdef double dimless_moment_of_inertia_i(double x, int obl_surfgrav_ind) noexcept nogil:
    """
    Dimensionless moment of inertia i(x) defined by:

        I = i(x) * M * Re^2

    Parameters
    ----------
    x : double
        Compactness ratio x = M / Re.
    model : obl_surfgrav
        0 : AlGendy & Morsink (2014) slow-rotation universal fit
        5   : Empirical fit from LORENE + CUTTER

    Returns
    -------
    double
        i(x) (dimensionless).
    """
    cdef double a0, a1, a2

    if obl_surfgrav_ind == 0:
        # Eq. (15): i(x) = sqrt(x) * (a0 + a1 x + a2 x^2)
        a0 = 1.136
        a1 = -2.53
        a2 = 5.6
        return sqrt(x) * (a0 + a1 * x + a2 * x * x)

    elif obl_surfgrav_ind == 5:
        # I have not computed the MoI for this star
        return NAN



cdef double beta1_coeff(int obl_surfgrav_ind) noexcept nogil:
    """
    Leading-order coefficient beta1 in the slow-rotation expansion:

        beta = beta1 * x * Omega_bar^2

    This returns beta1 only (NOT beta). This is intentional: it prevents
    accidental double-counting or missing factors of x and Ω̄^2.

    Parameters
    ----------
    model : obl_surfgrav
        0 : AlGendy & Morsink (2014) slow-rotation universal fit
        5   : Empirical fit from LORENE + CUTTER

    Returns
    -------
    double
        beta1 (dimensionless).
    """
    if obl_surfgrav_ind == 0:
        return 0.4454
    elif obl_surfgrav_ind == 5:
        # Not computed for J1231
        return NAN

cpdef double py_dimless_moment_of_inertia_i(double x, int obl_surfgrav_ind):
    return dimless_moment_of_inertia_i(x, obl_surfgrav_ind)

cpdef double py_oblateness_func_o2(double x, int obl_surfgrav_ind):
    return oblateness_func_o2(x, obl_surfgrav_ind)

cpdef double py_beta1_coeff(int obl_surfgrav_ind):
    return beta1_coeff(obl_surfgrav_ind)
