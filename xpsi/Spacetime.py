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

    :param string star_shape:
        A string specifying the assumed shape of the star. Options 'AGM_14'
        (an oblate spheroid from Algendy & Morsink 2014) or 'sphere' are
        currently allowed.

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
        "psr j0614-3329": 6,}

    # Aliases -> canonical key
    STAR_OBLATENESS_SURFGRAV_ALIASES = {
        # universal/AGM-type synonyms
        "default": "universal",
        "agm_14": "universal",
        "algendy": "universal",
        "morsink": "universal",
        "morsink 2007": "universal",
        "universal relation": "universal",

        # shorthand pulsar names
        "j0030": "psr j0030+0451",
        "j0740": "psr j0740+6620",
        "j0437": "psr j0437-4715",
        "j1231": "psr j1231-1411",
        "j0614": "psr j0614-3329",}

    required_names = ['frequency',
                      'mass',
                      'radius',
                      'distance',
                      'cos_inclination']

    def __init__(self, bounds, values, obl_surfgrav=None):

        self._obl_surfgrav = None
        self._obl_surfgrav_ind = None

        if obl_surfgrav is not None:
            self.obl_surfgrav = obl_surfgrav


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
                      strict_bounds = (1.0, 40.0),
                      bounds = bounds.get('radius', None),
                      doc = 'Coordinate equatorial radius [km]',
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

        try:
            ind = cls.STAR_OBLATENESS_SURFGRAV_CANONICAL[key]
        except KeyError as e:
            allowed = _format_allowed(cls.STAR_OBLATENESS_SURFGRAV_CANONICAL.keys())
            raise ValueError(
                f"Unknown oblateness and surface gravity for '{name}'. "
                f"Allowed canonical keys: {allowed}") from e

        return key, int(ind)

    @classmethod
    def available_obl_surfgrav(cls):
        return tuple(sorted(cls.STAR_OBLATENESS_SURFGRAV_CANONICAL.keys()))

    @property
    def obl_surfgrav(self):
        return self._obl_surfgrav

    @obl_surfgrav.setter
    def obl_surfgrav(self, name):
        key, ind = self.resolve_obl_surfgrav(name)
        self._obl_surfgrav = key
        self._obl_surfgrav_ind = ind

    @property
    def obl_surfgrav_ind(self):
        if self._obl_surfgrav_ind is None:
            raise RuntimeError(
                "obl_surfgrav was not set. Set it via Star(name=...) or spacetime.obl_surfgrav=...")
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
        """ Get the coordinate equatorial radius in SI. """
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

        A dimensionless function of stellar properties. This is the  compactness

        See Morsink et al. (2007), and AlGendy & Morsink (2014).

        """
        try:
            return self._zeta
        except AttributeError:
            return self.r_g / self.R

    @property
    def epsilon(self):
        """ Get the derived parameter ``epsilon`` for universal relations.

        A dimensionless function of stellar properties. This is omega bar squared

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

            I_dimless = _m.sqrt(zeta) * dimless_moment_of_inertia_i(zeta, self._obl_surfgrav_ind)

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
        """ Get the dimensionless mass quadrupole, second order in spin.

        See AlGendy & Morsink (2014).

        """
        try:
            return self._q
        except AttributeError:
             temp = self.epsilon * 0.11 / (self.zeta * self.zeta)
             temp -= self.a * self.a / (self.r_g * self.r_g)
             return temp + beta1_coeff(self._obl_surfgrav_ind) * 4.0 * self.epsilon * self.zeta / 3.0

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
