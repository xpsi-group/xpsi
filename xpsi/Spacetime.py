from __future__ import division, print_function

from .global_imports import *
from . import global_imports

from abc import ABCMeta, abstractmethod
from .ParameterSubspace import ParameterSubspace, BoundsError

class Spacetime(ParameterSubspace):
    """ The ambient Schwarzschild spacetime. """

    @abstractmethod
    def __init__(self, num_params, bounds):
        """
        :param int num_params: Number of free spacetime parameters.

        :param list bounds: Tuples of hard bounds on parameters.

        """
        super(Spacetime, self).__init__(num_params, bounds)

        if self._num_params < 4:
            raise BoundsError('A spacetime requires at least three parameters.')

        for mass in self._bounds[1]:
            if not 0.8 <= mass <= 3.0:
                raise BoundsError('Invalid Schwarzschild mass bound.')

        for radius in self._bounds[2]:
            if not 0.0 <= radius <= 30.0:
                raise BoundsError('Invalid coordinate equatorial radius bound.')

        r_g = _G * self._bounds[1][0] * _M_s / _csq

        if self._bounds[2][0] * _km <= 3.0 * r_g:
            raise BoundsError('Lower radius bound is at or within the'
                              'photon sphere for the lower mass bound.')

        for inclination in self._bounds[3]:
            if not 0.0 < inclination < _pi:
                raise BoundsError('Invalid inclination bound.')

    def update(self, d, M, R, i):
        """
        :param M: The (rotationally deformed) gravitational mass (solar masses).
        :param R: The coordinate equatorial radius (km).
        :param i: Inclination of the Earth to the rotational axis (radians).

        :raises TypeError: If the number of arguments passed based on the
                           ``num_params`` property does not match the
                           number of arguments in the function signature.

        :raises AttributeError: If the coordinate angular frequency is not
                                set in the initialiser of a custom subclass, or
                                in a custom :meth:`update` method of a subclass.

        """
        self._M = M * _M_s
        self._r_g = _G * self._M / _csq
        self._r_s = 2.0 * self._r_g

        self._R = R * _km
        self._R_r_s = self._R / self._r_s

        # observer parameters
        self._i = i
        self.d = d

        self._zeta = self._r_g / self._R
        self._epsilon = self._Omega**2.0 * self._R**3.0 / (_G * self._M)

    @property
    def num_params(self):
        """ Get the number of spacetime parameters. """
        return self._num_params

    @property
    def M(self):
        """ Get the (rotationally deformed) gravitational mass. """
        return self._M

    @property
    def r_g(self):
        """ Get the Schwarzschild gravitational radius. """
        return self._r_g

    @property
    def r_s(self):
        """ Get the Schwarzschild radius. """
        return self._r_s

    @property
    def R(self):
        """ Get the coordinate equatorial radius. """
        return self._R

    @property
    def R_r_s(self):
        """ Get the ratio of the radius to the Schwarzschild radius. """
        return self._R_r_s

    @property
    def S(self):
        """ Get the coordinate rotation frequency.

        :raises AttributeError: If the coordinate angular frequency is not
                                set in the initialiser of a custom subclass, or
                                in a custom :meth:`update` method of a subclass.

        """
        return self._S

    @property
    def Omega(self):
        """Get the coordinate *angular* rotation frequency.

        :raises AttributeError: If the coordinate angular frequency is not
                                set in the initialiser of a custom subclass, or
                                in a custom :meth:`update` method of a subclass.

        """
        return self._Omega

    @property
    def i(self):
        """ Get the inclination of the Earth to the rotational axis. """
        return self._i

    @property
    def d(self):
        """ Get the distance of the Earth to the rotational axis. """
        return self._d

    @d.setter
    def d(self, value):
        self._d = value * _kpc

    @property
    def d_sq(self):
        return self._d * self._d

    @property
    def zeta(self):
        """ Get the derived parameter ``zeta`` for universality relation. """
        return self._zeta

    @property
    def epsilon(self):
        """ Get the derived parameter ``epsilon`` for universality relation. """
        return self._epsilon




