from __future__ import division, print_function

from .global_imports import *
from . import global_imports

from .ParameterSubspace import ParameterSubspace

from .Spacetime import Spacetime
from .Photosphere import Photosphere

class Star(ParameterSubspace):
    """ Instances of :class:`Star` represent model stars.

    Each model star is abstractly constructed from an ambient spacetime,
    :class:`Spacetime`, and a collection of photosphere objects which are each
    embedded in that ambient spacetime for disjoint intervals of coordinate
    time.

    :param obj spacetime: An instance of :class:`~.Spacetime.Spacetime`.

    :param list photospheres: Each element must be an instance of
                              :class:`~.Photosphere`.

    """
    def __init__(self, spacetime, photospheres):
        if not isinstance(spacetime, Spacetime):
            raise TypeError('Invalid type for ambient spacetime object.')

        self._spacetime = spacetime

        try:
            for photosphere in photospheres:
                assert isinstance(photosphere, Photosphere)
        except AssertionError:
            try:
                assert isinstance(photospheres, Photosphere)
            except AssertionError:
                raise TypeError('Invalid type for a photosphere object.')
            else:
                self._photospheres = [photospheres]
        else:
            self._photospheres = photospheres

        # checks passed, so store reference between objects
        for photosphere in self._photospheres:
            photosphere.spacetime = self._spacetime

        # by default turn off fast_mode in HotRegion objects
        # the likelihood object will activate this if it is needed
        self.activate_fast_mode(False)

        super(Star, self).__init__(spacetime, photospheres)

    @property
    def spacetime(self):
        """ Get the ambient spacetime object. """
        return self._spacetime

    @property
    def photospheres(self):
        """ Get the list of photosphere objects. """
        return self._photospheres

    def activate_fast_mode(self, activate):
        try:
            for photosphere in self._photospheres:
                photosphere.hot.fast_mode = activate
        except AttributeError:
            pass # no hot regions to worry about

    def update(self, fast_counts=None, threads=1):
        """ Update the star.

        :param tuple fast_counts:
            Total count numbers from two or more temperature components of
            a hot region. Each *tuple* element is itself a *tuple* over hot
            regions, and if a hot region has two temperature components
            that were integrated over during the fast mode, the total counts
            per component are elements of another *tuple*. This nested
            containerisation is handled automatically.

        :param int threads:
            Number of ``OpenMP`` threads to spawn for embedding
            photosphere objects in the ambient spacetime.

        """

        if fast_counts is None:
            fast_counts = tuple([None]*len(self._photospheres))

        # Iteratively embed each photosphere (that needs to be updated)
        # in the ambient spacetime
        for photosphere, fast_count in zip(self._photospheres, fast_counts):
            if photosphere.needs_update or self._spacetime.needs_update:
                photosphere.embed(fast_count, threads)
