from __future__ import division, print_function

from .global_imports import *
from . import global_imports

from .HotRegion import HotRegion

from .ParameterSubspace import ParameterSubspace

class PulseError(xpsiError):
    """ Raised if a numerical problems encountered during integration. """

class HotRegions(ParameterSubspace):
    """ Two photospheric hot regions, where the hot regions are objects.

    This class could be extended in principle to operate with multiple hot
    regions, but the computational expense scales ~linearly with number.
    Applications thus far have used two distinct hot regions with equal
    and unequal complexity.

    :param tuple hotregions:
            Two-element container of :class:`.HotRegion.HotRegion` instances.

    """

    def __init__(self, hotregions):
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
                if not isinstance(obj, HotRegion):
                    raise ValueError('Invalid type for hot-region object.')
        else:
            raise ValueError('Hot-region container must be iterable.')

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
        """ Embed the hot regions. """

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
                  hot_atmosphere, elsewhere_atmosphere):
        """ Integrate over the photospheric radiation field.

        Calls the CellMesh integrator, with or without exploitation of
        azimuthal invariance of the radiation field.

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
                                             hot_atmosphere,
                                             elsewhere_atmosphere)
        secondary = self._objects[1].integrate(st, secondary_energies,
                                               threads,
                                               hot_atmosphere,
                                               elsewhere_atmosphere)

        return (primary, secondary)
