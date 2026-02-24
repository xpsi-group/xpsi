# ns_universal_fits.pyx
# cython: language_level=3
# cython: boundscheck=False, wraparound=False, cdivision=True, nonecheck=False, initializedcheck=False

from libc.math cimport sqrt, NAN



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
    if obl_surfgrav_ind == 0:
        # Eq. (21): o2 = o20 + o21 x   (overall Ω̄^2 factor is applied elsewhere)
        return -0.788 + 1.030 * x

    elif obl_surfgrav_ind == 1:

        return 1.


    elif obl_surfgrav_ind == 0:
        # Empirical rational + quadratic fit (your existing coefficients)
        cdef double c1 = -0.7940771
        cdef double c2 =  0.019817498
        cdef double c3 =  1.7524167
        return c1 * x / (x + c2) + c3 * x * x


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
        # If you truly want a different normalization for LORENE fits:
        a0 = 1.176
        a1 = -2.53
        a2 = 5.6
        return sqrt(x) * (a0 + a1 * x + a2 * x * x)



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
        return 0.4554
