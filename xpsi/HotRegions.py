
from xpsi.ParameterSubspace import ParameterSubspace

from xpsi.HotRegion import HotRegion

class HotRegions(ParameterSubspace):
    """ Two or more photospheric hot regions, where the hot regions are objects.

    The hot regions can be equal or different in parametrisation complexity.

    In principle you can construct multiple hot regions, but the computational
    expense scales ~linearly with number. Applications thus far have used two
    distinct hot regions with equal and unequal complexities.

    :param tuple hotregions:
        Iterable container of :class:`.HotRegion.HotRegion` instances.

    """

    def __init__(self, hotregions):

        self.objects = hotregions
        super(HotRegions, self).__init__(hotregions)

    @property
    def objects(self):
        """ Get an iterable of hot region objects. """
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
    def phases_in_cycles(self):
        return [obj.phases_in_cycles for obj in self.objects]

    @property
    def fast_phases_in_cycles(self):
        """  """
        if self.do_fast:
            return [obj.fast_phases_in_cycles for obj in self.objects]

    @property
    def do_fast(self):
        """ Should fast-mode be invoked for some subset of hot regions. """
        for obj in self._objects:
            if obj.do_fast:
                return True
        return False

    @property
    def fast_mode(self):
        """ Is fast-mode currently activated? """
        for obj in self._objects:
            if obj.fast_mode:
                return True
        return False

    @fast_mode.setter
    def fast_mode(self, activate):
        """ Activate or deactivate fast-mode. """
        for obj in self._objects:
            obj.fast_mode = activate

    def print_settings(self):
        """ Print numerical settings. """
        for obj in self._objects:
            obj.print_settings()

    def __len__(self):
        """ Get the number of hot region objects referenced.

        Overwrites the base class magic method.

        """
        return len(self._objects)

    def embed(self, spacetime, photosphere, fast_total_counts, threads, *args):
        """ Embed the hot regions. """

        if fast_total_counts is None:
            fast_total_counts = [None] * len(self)

        for obj, fast in zip(self._objects, fast_total_counts):
            obj.embed(spacetime,
                      photosphere,
                      fast,
                      threads, *args)

    def integrate(self, st, energies, threads,
                  hot_atmosphere, elsewhere_atmosphere):
        """ Integrate over the photospheric radiation field.

        Calls the CellMesh integrator, with or without exploitation of
        azimuthal invariance of the radiation field.

        :param st: Instance of :class:`~.Spacetime.Spacetime`.

        :param energies:
            A one-dimensional :class:`numpy.ndarray` of energies in keV.

        :param int threads:
            Number of ``OpenMP`` threads for pulse integration.

        """
        if not isinstance(energies, tuple):
            energies = [energies] * len(self)

        signals = []
        for obj, E in zip(self._objects, energies):
            signals.append(obj.integrate(st, E,
                                         threads,
                                         hot_atmosphere,
                                         elsewhere_atmosphere))

        return tuple(signals)
