from xpsi.global_imports import *

from xpsi.utils import make_verbose

from xpsi.Parameter import Parameter
from xpsi.ParameterSubspace import ParameterSubspace

from xpsi.surface_radiation_field.effective_gravity_universal import (py_dimless_moment_of_inertia_i as dimless_moment_of_inertia_i,
                                                             py_oblateness_func_o2 as oblateness_func_o2,
                                                             py_beta1_coeff as beta1_coeff)


def _normalise_key(name):
    """Normalise user input for robust key matching."""
    return " ".join(name.strip().lower().split())


def _format_allowed(keys):
    return ", ".join(sorted(keys))


class Spacetime(ParameterSubspace):
    """ The ambient Schwarzschild spacetime and Earth coordinates.

    :param dict bounds:
        Tuples of hard bounds on parameters with keys matching the names
        in the initialiser body. A ``None`` bound will be interpreted as
        fixed variable, in which case a value must be supplied in the
        corresponding dictionary. If there is no entry for a default
        parameter, it will be assumed to be fixed, but an exception
        will be raised if there is no associated value.

    :param dict values:
        Values of fixed parameters or initial values for free parameters,
        with keys matching the names in the initialiser body. If there is
        no entry for a default parameter, no initial value will be specified,
        but an exception will be raised if there is also no bound specified.

    The prescription used for stellar oblateness and effective surface
    gravity is set with the ``name`` argument of :class:`~.Star.Star`, and
    can be read from :attr:`obl_surfgrav`. It is ``'universal'``
    (AlGendy & Morsink 2014) until a :class:`~.Star.Star` sets it.

    We define a property for parameters and combinations of parameters to
    shortcut access, given that this subspace is passed to other subspaces
    for model computation. We would like to access the values with fewer
    characters than ``self[<name>]``.

    """

    STAR_OBLATENESS_SURFGRAV_CANONICAL = {
        "universal": 0,
        "sphere": 1,
        "psr j0030+0451": 2,
        "psr j0740+6620": 3,
        "psr j0437-4715": 4,
        "psr j1231-1411": 5,
        "psr j0614-3329": 6,
        "psr j1614-2230": 7,}

    # Reserved for relations whose coefficients are not implemented yet.
    STAR_OBLATENESS_SURFGRAV_NOT_IMPLEMENTED = (
        "psr j0030+0451",
        "psr j0740+6620",
        "psr j0437-4715",
        "psr j0614-3329",)

    # Relations without a moment of inertia and quadrupole coefficient, so
    # the spin (a) and quadrupole (q) needed by the PhotospherePlotter are
    # not available. The spherical star is only meant for comparison
    # calculations with the Schwarzschild + Doppler approximation.
    STAR_OBLATENESS_SURFGRAV_NO_SPIN_QUADRUPOLE = (
        "sphere",
        "psr j1231-1411",
        "psr j1614-2230",)

    # Aliases -> canonical key
    STAR_OBLATENESS_SURFGRAV_ALIASES = {
        # universal/AGM-type synonyms
        "default": "universal",
        "agm_14": "universal",
        "algendy": "universal",
        "morsink": "universal",
        "universal relation": "universal",

        # shorthand pulsar names
        "j0030": "psr j0030+0451",
        "j0740": "psr j0740+6620",
        "j0437": "psr j0437-4715",
        "j1231": "psr j1231-1411",
        "j1614": "psr j1614-2230",}

    required_names = ['frequency',
                      'mass',
                      'radius',
                      'distance',
                      'cos_inclination']

    def __init__(self, bounds, values):

        # Default until a Star sets it (see Star.__init__).
        self._set_obl_surfgrav('universal')

        if not isinstance(bounds, dict) or not isinstance(values, dict):
             raise TypeError("Both bounds and values need to be dictionaries.")

        f = Parameter('frequency',
                      strict_bounds = (0.0, 800.0),
                      bounds = bounds.get('frequency', None),
                      doc = 'Spin frequency [Hz]',
                      symbol = r'$f$',
                      value = values.get('frequency', None))

        M = Parameter('mass',
                      strict_bounds = (0.001, 3.0),
                      bounds = bounds.get('mass', None),
                      doc = 'Gravitational mass [solar masses]',
                      symbol = r'$M$',
                      value = values.get('mass', None))

        R = Parameter('radius',
                      strict_bounds = (1.0, 20.0),
                      bounds = bounds.get('radius', None),
                      doc = 'Circumferential equatorial radius [km]',
                      symbol = r'$R_{\rm eq}$',
                      value = values.get('radius', None))

        D = Parameter('distance',
                      strict_bounds = (0.01, 30.0), # inside Milky Way
                      bounds = bounds.get('distance', None),
                      doc = 'Earth distance [kpc]',
                      symbol = r'$D$',
                      value = values.get('distance', None))

        cosi = Parameter('cos_inclination',
                         strict_bounds = (-1.0, 1.0),
                         bounds = bounds.get('cos_inclination', None),
                         doc = 'Cosine of Earth inclination to rotation axis',
                         symbol = r'$\cos(i)$',
                         value = values.get('cos_inclination', None))

        super(Spacetime, self).__init__(f, M, R, D, cosi)

    @classmethod
    def resolve_obl_surfgrav(cls, name):
        """Resolve user-provided name/alias to (canonical_key, backend_index)."""
        if not isinstance(name, str) or not name.strip():
            raise TypeError("obl_surfgrav must be a non-empty string.")

        key = _normalise_key(name)
        key = cls.STAR_OBLATENESS_SURFGRAV_ALIASES.get(key, key)

        allowed = _format_allowed(cls.available_obl_surfgrav())

        try:
            ind = cls.STAR_OBLATENESS_SURFGRAV_CANONICAL[key]
        except KeyError as e:
            raise ValueError(
                f"Unknown oblateness and surface gravity for '{name}'. "
                f"Allowed canonical keys: {allowed}") from e

        if key in cls.STAR_OBLATENESS_SURFGRAV_NOT_IMPLEMENTED:
            raise NotImplementedError(
                f"The oblateness and surface-gravity relation for '{key}' "
                f"is not implemented yet. Available: {allowed}")

        return key, int(ind)

    @classmethod
    def available_obl_surfgrav(cls):
        """ Get the names of the implemented prescriptions. """
        return tuple(sorted(k for k in cls.STAR_OBLATENESS_SURFGRAV_CANONICAL
                            if k not in cls.STAR_OBLATENESS_SURFGRAV_NOT_IMPLEMENTED))

    def _set_obl_surfgrav(self, name):
        """ Set the oblateness and surface-gravity prescription.

        Only meant to be called by :class:`~.Star.Star`, which is where users
        choose the prescription.

        """
        key, ind = self.resolve_obl_surfgrav(name)
        self._obl_surfgrav = key
        self._obl_surfgrav_ind = ind

    @property
    def obl_surfgrav(self):
        """ Get the oblateness and surface-gravity prescription (read-only).

        Set it with the ``name`` argument of :class:`~.Star.Star`.

        """
        return self._obl_surfgrav

    @property
    def obl_surfgrav_ind(self):
        """ Get the backend index of :attr:`obl_surfgrav` (read-only). """
        return self._obl_surfgrav_ind

    @property
    def M(self):
        """ Get the (rotationally deformed) gravitational mass in SI. """
        return self['mass'] * _GM * _csq / _G

    @property
    def r_g(self):
        """ Get the Schwarzschild gravitational radius in SI. """
        return self['mass'] * _GM

    @property
    def r_s(self):
        """ Get the Schwarzschild radius in SI. """
        return 2.0 * self.r_g

    @property
    def R(self):
        """ Get the circumferential equatorial radius in SI. """
        return self['radius'] * _km

    @property
    def R_r_s(self):
        """ Get the ratio of the equatorial radius to the Schwarzschild radius.
        """
        return self.R / self.r_s

    @property
    def f(self):
        """ Get the coordinate rotation frequency. """
        return self['frequency']

    @property
    def Omega(self):
        """ Get the coordinate *angular* rotation frequency. """
        return _2pi * self['frequency']

    @property
    def i(self):
        """ Get the inclination of the Earth to the rotational axis. """
        return _m.acos(self['cos_inclination'])

    @property
    def d(self):
        """ Get the distance of the Earth to the star in SI. """
        return self['distance'] * _kpc

    @property
    def d_sq(self):
        """ Get the squared distance of the Earth to the star in SI. """
        return self.d * self.d

    @property
    def zeta(self):
        """ Get the derived parameter ``zeta`` for universal relations.

        A dimensionless function of stellar properties.

        See Morsink et al. (2007), and AlGendy & Morsink (2014).

        """
        try:
            return self._zeta
        except AttributeError:
            return self.r_g / self.R

    @property
    def epsilon(self):
        """ Get the derived parameter ``epsilon`` for universal relations.

        A dimensionless function of stellar properties.

        See Morsink et al. (2007), and AlGendy & Morsink (2014).

        """
        try:
            return self._epsilon
        except AttributeError:
            return self.Omega**2.0 * self.R**3.0 / (_G * self.M)

    @property
    def a(self):
        """ Get the spin parameter, first order in spin.

        See AlGendy & Morsink (2014). a=J/(Mc)

        """
        try:
            return self._a
        except AttributeError:
            zeta = self.zeta

            # The returned i(x) already includes the factor sqrt(zeta).
            I_dimless = dimless_moment_of_inertia_i(zeta, self._obl_surfgrav_ind)

            a = self.R * self.R * self.Omega * I_dimless / _c

            return a

    @a.setter
    def a(self, a):
        """ Set the spin parameter. """
        self._a = a

    @a.deleter
    def a(self):
        """ Delete the spin parameter. """
        try:
            del self._a
        except AttributeError:
            pass # silently do nothing

    @property
    def q(self):
        """ Get the quadrupole deviation parameter of the quasi-Kerr metric,
        second order in spin.

        This is not the star's quadrupole moment. It is the parameter that
        measures how far the quadrupole differs from the Kerr value, called
        epsilon in Glampedakis & Babak (2006, Eq. 2.7) and in Psaltis &
        Johannsen (2012, Eq. 6). (It is unrelated to :attr:`epsilon`
        below.) The star's quadrupole moment is then::

            Q = -M (a^2 + q M^2)

        so ``q = 0`` gives the Kerr quadrupole.

        The value comes from the slow-rotation fits of AlGendy & Morsink
        (2014)::

            q = -q_AM + (4/3) beta - j^2

        where

        - ``q_AM = q2 * Omega_bar^2 / x^2`` is their dimensionless quadrupole
          (Eq. 17). It is not coordinate invariant, and it is negative for
          an oblate star.
        - ``beta = beta1 * Omega_bar^2 * x`` (Eq. 19).
        - ``j = a / r_g`` is the dimensionless spin.
        - ``x`` is :attr:`zeta` and ``Omega_bar^2`` is :attr:`epsilon`.

        ``q_AM - (4/3) beta`` is the coordinate-invariant (Geroch-Hansen)
        quadrupole, divided by M^3. Subtracting ``j^2`` removes the Kerr
        part, which the metric already contains.

        Note the sign of the beta term. AlGendy & Morsink (2014, Eq. 18)
        print ``q_inv = q_AM + (4/3) beta``, and Nattila & Pihajoki (2018,
        Eq. 6) repeat it. Their own definitions of q_AM and beta, together
        with Pappas & Apostolatos (2012, Eq. 11), give a minus sign. The
        exact Kerr limit (q_AM = -2j^2/3, beta = j^2/4) confirms it:
        q_AM - (4/3) beta = -j^2, the Kerr value, whereas the plus sign gives
        -j^2/3. Do not "correct" the sign below to match the paper.

        Three different "q" conventions appear in these papers:

        - AlGendy & Morsink's q_AM is negative for an oblate star.
        - Glampedakis & Babak's own q = -Q/M^3 = j^2 + epsilon
          (their Eqs. 2.5 and 2.7) is positive.
        - This property is the deviation epsilon, which is neither.

        """
        try:
            return self._q
        except AttributeError:
            # -q_AM: q2 = -0.11 (AlGendy & Morsink 2014, Table 1)
            temp = self.epsilon * 0.11 / (self.zeta * self.zeta)
            # -j^2, with j = a / r_g
            temp -= self.a * self.a / (self.r_g * self.r_g)
            # +(4/3) beta: minus the beta term of q_inv (see docstring)
            return temp + beta1_coeff(self._obl_surfgrav_ind)* 4.0 * self.epsilon * self.zeta / 3.0


    @q.setter
    def q(self, q):
        """ Set the mass quadrupole moment. """
        self._q = q

    @q.deleter
    def q(self):
        """ Delete the mass quadrupole moment. """
        try:
            del self._q
        except AttributeError:
            pass # silently do nothing

Spacetime._update_doc()
