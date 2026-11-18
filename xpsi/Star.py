
from xpsi.ParameterSubspace import ParameterSubspace

from xpsi.Spacetime import Spacetime
from xpsi.Photosphere import Photosphere

class Star(ParameterSubspace):
    """ Instances of :class:`Star` represent model stars.

    Each model star is abstractly constructed from an ambient spacetime,
    :class:`Spacetime`, and a collection of photosphere objects which are each
    embedded in that ambient spacetime for disjoint intervals of coordinate
    time.

    :param obj spacetime: An instance of :class:`~.Spacetime.Spacetime`.

    :param list photospheres: Each element must be an instance of
                              :class:`~.Photosphere`.

    :param str name:
        The prescription for the stellar oblateness and effective surface
        gravity. Options are ``'universal'`` (AlGendy & Morsink 2014; the
        default, also used when ``None``), ``'sphere'``, and the
        pulsar-specific relations ``'psr j1231-1411'`` and
        ``'psr j1614-2230'`` (see also the aliases in
        :attr:`~.Spacetime.Spacetime.STAR_OBLATENESS_SURFGRAV_ALIASES`).
        The pulsar-specific relations cannot yet be used with
        :class:`~.PhotospherePlotter.PhotospherePlotter`.
        The selected prescription is stored in the spacetime and can be read
        from :attr:`~.Spacetime.Spacetime.obl_surfgrav`.

    """
    def __init__(self, spacetime, photospheres, name='universal'):
        if not isinstance(spacetime, Spacetime):
            raise TypeError('Invalid type for ambient spacetime object.')

        self._spacetime = spacetime

        self._spacetime._set_obl_surfgrav('universal' if name is None else name)

        try:
            assert isinstance(photospheres, list)
        except AssertionError:
            try:
                assert isinstance(photospheres, Photosphere)
            except AssertionError:
                raise TypeError('Invalid type for a photosphere object')
            else:
                self._photospheres = [photospheres]
        else:
            try:
                for photosphere in photospheres:
                    assert isinstance(photosphere, Photosphere)
            except AssertionError:
                raise TypeError('Invalid type for a photosphere object')
            else:
                self._photospheres = photospheres

        # checks passed, so store reference between objects
        for photosphere in self._photospheres:
            photosphere.spacetime = self._spacetime

        super(Star, self).__init__(spacetime, photospheres)

    @property
    def spacetime(self):
        """ Get the ambient spacetime object. """
        return self._spacetime

    @property
    def photospheres(self):
        """ Get the list of photosphere objects. """
        return self._photospheres

    @photospheres.setter
    def photospheres(self, obj):
        self._photospheres = [obj]

    def update(self, threads=1, force_update=False):
        """ Update the star.

        :param int threads:
            Number of ``OpenMP`` threads to spawn for embedding
            photosphere objects in the ambient spacetime.

        :param bool force_update:
            Setting to force update even if both the photosphere and spacetime say they do not need updating.

        """

        # Iteratively embed each photosphere (that needs to be updated)
        # in the ambient spacetime
        for photosphere in self._photospheres:
            if photosphere.needs_update or self._spacetime.needs_update or force_update:
                photosphere.embed(threads)
