from __future__ import division, print_function

from .global_imports import *
from . import global_imports

from . import make_verbose

from .Parameter import Parameter
from .ParameterSubspace import ParameterSubspace

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

    We define a property for parameters and combinations of parameters to
    shortcut access, given that this subspace is passed to other subspaces
    for model computation. We would like to access the values with fewer
    characters than ``self[<name>]``.

    """
    required_names = ['frequency',
                      'mass',
                      'radius',
                      'distance',
                      'cos_inclination']

    def __init__(self, bounds, values):

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
    @make_verbose('Configuring default bounds with fixed spin',
                  'Spacetime configured')
    def fixed_spin(cls, frequency):
        """
        :param float frequency:
            The fixed coordinate spin frequency in Hz.

        .. note::

            The degeneracy due to equatorially-reflection symmetric physics
            is eliminated here by declaring prior support for Earth inclination
            from northern rotation pole to equatorial plane.

        """

        bounds = dict(mass = (1.0, 3.0),
                      radius = (gravradius(1.0), 16.0),
                      distance = (0.05, 2.0),
                      cos_inclination = (0.0, 1.0))

        return cls(bounds, dict(frequency = frequency))

    @property
    def M(self):
        """ Get the (rotationally deformed) gravitational mass in SI. """
        return self['mass'] * _M_s

    @property
    def r_g(self):
        """ Get the Schwarzschild gravitational radius in SI. """
        return _G * self.M / _csq

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

        See AlGendy & Morsink (2014).

        """
        try:
            return self._a
        except AttributeError:
            zeta = self.zeta

            I_dimless = _m.sqrt(zeta) * (1.136 - 2.53 * zeta + 5.6 * zeta * zeta)

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
             return temp + 0.4554 * 4.0 * self.epsilon * self.zeta / 3.0

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
