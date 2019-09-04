from __future__ import division, print_function

from .global_imports import *
from . import global_imports

from .Spot import Spot

from .ParameterSubspace import ParameterSubspace

class PulseError(xpsiError):
    """ Raised if a numerical problems encountered during integration. """

class TwoSpots(ParameterSubspace):
    """ Two photospheric spots, where the spots are objects.

    """

    def __init__(self, spots):
        """
        :param tuple num_params: Number of parameters for spot model, including
                                 the primary spot centre colatitude and spot
                                 angular radius. The tuple must have two
                                 elements, specifying the number of parameters
                                 per spot.

        :param list bounds: Hard parameter bounds for the instance of
                            :class:`.ParameterSubspace.ParameterSubspace`.

        :param bool antipodal_symmetry: Apply antipodal reflection symmetry?
                            This means that all parameters describing the
                            secondary spot are derived from the primary spot:
                            the secondary is antipodal and given by an
                            equatorial reflection of the primary spot,
                            followed by a rotation by :math:`\pi` radians
                            about the rotation axis. The mirroring is thus
                            with respect to a 2-plane through the coordinate
                            origin which is perpendicular to the line
                            through the origin and the primary spot centre.

        :param kwargs: Keyword arguments passed to :class:`.Spot.Spot` class.

        .. note:: The parameter bounds of the secondary spot must satisfy:

            * the colatitude is restricted to :math:`\Theta\in(0.0,\pi)`
            * the angular radius is restricted to :math:`\zeta\in(0.0,\pi/2)`

        """
        self.objects = spots

        self._num_primary_params = self._objects[0].num_params
        self._num_secondary_params = self._objects[1].num_params

        super(TwoSpots, self).__init__(self._num_primary_params\
                                       + self._num_secondary_params,
                                       self._objects[0].bounds\
                                       + self._objects[1].bounds)

        self.phases_in_cycles = self._objects[0].phases_in_cycles
        self.fast_phases_in_cycles = self._objects[0].fast_phases_in_cycles

    @property
    def objects(self):
        return self._objects

    @objects.setter
    def objects(self, objs):
        if isinstance(objs, tuple) or isinstance(objs, list):
            for obj in objs:
                if not isinstance(obj, Spot):
                    raise ValueError('Invalid type for spot object.')
        else:
            raise ValueError('Spot container must be iterable.')

        self._objects = objs

    @property
    def do_fast(self):
        for obj in self._objects:
            if obj.do_fast:
                return True
        return False

    @property
    def fast_mode(self):
        for obj in self._objects:
            if obj.fast_mode:
                return True
        return False

    @fast_mode.setter
    def fast_mode(self, activate):
        for obj in self._objects:
            obj.fast_mode = activate

    def print_settings(self):
        """ Print numerical settings. """
        for obj in self._objects:
            obj.print_settings()

    def embed(self, spacetime, p, fast_total_counts, threads, *args):
        """ Embed the spots. """

        if fast_total_counts is not None:
            fast_primary_total_counts = fast_total_counts[0]
            fast_secondary_total_counts = fast_total_counts[1]
        else:
            fast_primary_total_counts = None
            fast_secondary_total_counts = None

        pp = p[:self._num_primary_params]
        ps = p[self._num_primary_params:]

        self._objects[0].embed(spacetime, pp,
                               fast_primary_total_counts,
                               threads, *args)
        self._objects[1].embed(spacetime, ps,
                               fast_secondary_total_counts,
                               threads, *args)

    def integrate(self, st, energies, threads,
                  spot_atmosphere, elsewhere_atmosphere):
        """ Integrate over the photospheric radiation field.

        Calls the CellMesh integrator, with or without exploitation of
        azimuthal invariance of the spot radiation field.

        :param st: Instance of :class:`~.Spacetime.Spacetime`.

        :param energies: A one-dimensional :class:`numpy.ndarray` of energies
                         in keV.

        :param int threads: Number of ``OpenMP`` threads for pulse
                            integration.

        """
        if isinstance(energies, tuple):
            primary_energies = energies[0]
            secondary_energies = energies[1]
        else:
            primary_energies = secondary_energies = energies

        primary = self._objects[0].integrate(st, primary_energies,
                                             threads,
                                             spot_atmosphere,
                                             elsewhere_atmosphere)
        secondary = self._objects[1].integrate(st, secondary_energies,
                                               threads,
                                               spot_atmosphere,
                                               elsewhere_atmosphere)

        return (primary, secondary)
