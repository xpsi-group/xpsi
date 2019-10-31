from __future__ import division, print_function

from .global_imports import *
from . import global_imports

from .Spacetime import Spacetime
from .Photosphere import Photosphere

class Star(object):
    """ Instances of :class:`Star` represent model stars.

    Each model star is abstractly constructed from an ambient spacetime,
    :class:`Spacetime`, and a collection of photosphere objects which are each
    embedded in that ambient spacetime for disjoint intervals of coordinate
    time.

    """
    def __init__(self, spacetime, photospheres):
        """
        :param spacetime: An instance of :class:`~.Spacetime.Spacetime`.

        :param list photospheres: Each element must be an instance of
                                  :class:`~.Photosphere`.

        """
        try:
            assert isinstance(spacetime, Spacetime)
        except AssertionError:
            raise TypeError('Invalid type for ambient spacetime object.')
        else:
            self._spacetime = spacetime

        try:
            for photosphere in photospheres:
                assert isinstance(photosphere, Photosphere)
        except AssertionError:
            raise TypeError('Invalid type for a photosphere object.')
        except TypeError:
            try:
                assert isinstance(photospheres, Photosphere)
            except AssertionError:
                raise TypeError('Invalid type for a photosphere object.')
            else:
                self._photospheres = [photospheres]
        else:
            self._photospheres = photospheres

        self._eval_num_params()

    def _eval_num_params(self):
        """ Evaluate the number of *slow* parameters defining the star. """
        self._num_params = self._spacetime.num_params

        for photosphere in self._photospheres:
            self._num_params += photosphere.total_params

    @property
    def num_params(self):
        """ Get the number of (slow) parameters describing the star. """
        return self._num_params

    @property
    def spacetime(self):
        """ Get the ambient spacetime object. """
        return self._spacetime

    @property
    def photospheres(self):
        """ Get the list of photosphere objects. """
        return self._photospheres

    def activate_fast_mode(self, activate):
        for photosphere in self._photospheres:
            photosphere.hot.fast_mode = activate

    def update(self, p, fast_counts=None, threads=1):
        """ Update the star.

        :param list p: Model parameters relevant to ambient spacetime solution
                       and photosphere objects.

        :param int threads: Number of ``OpenMP`` threads to spawn for embedding
                            photosphere objects in the ambient spacetime.


        Parameter vector:

        * ``p[:i]`` = spacetime parameters
        * ``p[i:]`` = photospheric radiation field parameters

        where

        .. code-block:: python

            i = self._spacetime.num_params

        """

        # Update the ambient spacetime
        i = self._spacetime.num_params
        self._spacetime.update(*p[:i])

        # Iteratively embed each photosphere in the ambient spacetime
        for photosphere, fast_count in zip(self._photospheres, fast_counts):
            photosphere.embed(self._spacetime,
                              p[i:i + photosphere.total_params],
                              fast_count,
                              threads)
            i += photosphere.total_params





