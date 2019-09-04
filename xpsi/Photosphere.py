from __future__ import division, print_function

from .global_imports import *
from . import global_imports

from .Spot import Spot
from .Elsewhere import Elsewhere

from .ParameterSubspace import ParameterSubspace, BoundsError

class Photosphere(ParameterSubspace):
    """
    A photosphere embedded in an ambient spacetime for use with the CellMesh
    algorithm.

    The exterior domain of the photosphere is everywhere vacuum.

    """
    def __init__(self, num_params, bounds, tag, spot, elsewhere=None):
        """
        :param num_params: Number of parameters for the *global* photosphere,
                           which are not defined by the :class:`~.Spot.Spot`
                           or :class:`~.Elsewhere.Elsewhere` parameter
                           subspaces. At present this is restricted to the
                           coordinate rotation frequency if not *locked*
                           to the stellar rotation frequency.

        :param list bounds: Hard parameter bounds for the instance of
                            :class:`.ParameterSubspace.ParameterSubspace`.

        :param str tag: An identification tag to enforce intended pairing with
                        a :class:`~.Pulse.Pulse` object.

        :param spot: An instance of :class:`~.Spot.Spot` (or a derived class).

        :param elsewhere: An instance of :class:`~.Elsewhere.Elsewhere` (or a
                          derived class).

        """
        super(Photosphere, self).__init__(num_params, bounds)

        if self._num_params > 1:
            raise ValueError('Photosphere parameter subspace dimensionality '
                             'can be a maximum of one.')
        elif self._num_params == 1:
            for frequency in self._bounds[0]:
                if not 0.0 < frequency <= 1000.0:
                    raise BoundsError('Invalid spot rotation frequency bound.')

        try:
            self._tag = str(tag)
        except TypeError:
            raise TypeError('Incompatible type for identification tag.')

        if elsewhere is not None:
            try:
                assert isinstance(elsewhere, Elsewhere)
            except AssertionError:
                raise TypeError('Invalid type for an elsewhere object.')
            else:
                self._elsewhere = elsewhere
        else:
            self._elsewhere = None

        if not isinstance(spot, Spot): # including derived classes
            if hasattr(spot, 'objects'):
                for obj in getattr(spot, 'objects'):
                    if not isinstance(obj, Spot):
                        raise TypeError('Invalid type for a spot object')
            else:
                raise TypeError('Invalid type for a spot object.')

        self._spot = spot

        self._total_params = self._num_params + self._spot.num_params
        if self._elsewhere is not None:
            self._total_params += self._elsewhere.num_params

        self._spot_atmosphere = ()
        self._elsewhere_atmosphere = ()

    @property
    def total_params(self):
        """ Get the total number of photospheric radiation field parameters. """
        return self._total_params

    @property
    def tag(self):
        """Get the ID tag of the photosphere object for pulse pairing. """
        return self._tag

    @property
    def spot_atmosphere(self):
        """ Get the numerical atmosphere buffers for spots if used.

        To preload a numerical atmosphere into a buffer, subclass and
        overwrite the setter. The underscore attribute set by the setter
        must be an :math:`n`-tuple whose :math:`n^{th}` element is an
        :math:`(n-1)`-dimensional array flattened into a one-dimensional
        :class:`numpy.ndarray`. The first :math:`n-1`
        elements of the :math:`n`-tuple must each be an ordered one-dimensional
        :class:`numpy.ndarray` of parameter values for the purpose of
        multi-dimensional interpolation in the :math:`n^{th}` buffer. The
        first :math:`n-1` elements must be ordered to match the index
        arithmetic applied to the :math:`n^{th}` buffer. An example would be
        ``self._spot_atmosphere = (logT, logg, mu, logE, buf)``, where:
        ``logT`` is a logarithm of local comoving effective temperature;
        ``logg`` is a logarithm of effective surface gravity;
        ``mu`` is the cosine of the angle from the local surface normal;
        ``logE`` is a logarithm of the photon energy; and
        ``buf`` is a one-dimensional buffer of intensities of size given by
        the product of sizes of the first :math:`n-1` tuple elements.

        It is highly recommended that buffer preloading is used, instead
        of loading from disk in the customisable radiation field extension
        module, to avoid reading from disk for every spot and signal
        (likelihood) evaluation. This can be a non-negligible waste of compute
        resources. By preloading in Python, the memory is allocated and
        references to that memory are not in general deleted until a sampling
        script exits and the kernel stops. The likelihood callback accesses
        the same memory upon each call without I/O.

        """
        return self._spot_atmosphere

    @spot_atmosphere.setter
    def spot_atmosphere(self, path):
        """ Implement if required. """
        pass

    @property
    def elsewhere_atmosphere(self):
        """ Get the numerical atmosphere buffers for elsewhere if used.

        To preload a numerical atmosphere into a buffer, subclass and
        overwrite the setter. The underscore attribute set by the setter
        must be an :math:`n`-tuple whose :math:`n^{th}` element is an
        :math:`(n-1)`-dimensional array flattened into a one-dimensional
        :class:`numpy.ndarray`. The first :math:`n-1`
        elements of the :math:`n`-tuple must each be an ordered one-dimensional
        :class:`numpy.ndarray` of parameter values for the purpose of
        multi-dimensional interpolation in the :math:`n^{th}` buffer. The
        first :math:`n-1` elements must be ordered to match the index
        arithmetic applied to the :math:`n^{th}` buffer. An example would be
        ``self._spot_atmosphere = (logT, logg, mu, logE, buf)``, where:
        ``logT`` is a logarithm of local comoving effective temperature;
        ``logg`` is a logarithm of effective surface gravity;
        ``mu`` is the cosine of the angle from the local surface normal;
        ``logE`` is a logarithm of the photon energy; and
        ``buf`` is a one-dimensional buffer of intensities of size given by
        the product of sizes of the first :math:`n-1` tuple elements.

        It is highly recommended that buffer preloading is used, instead
        of loading from disk in the customisable radiation field extension
        module, to avoid reading from disk for every spot and signal
        (likelihood) evaluation. This can be a non-negligible waste of compute
        resources. By preloading in Python, the memory is allocated and
        references to that memory are not in general deleted until a sampling
        script exits and the kernel stops. The likelihood callback accesses
        the same memory upon each call without I/O.

        """
        return self._elsewhere_atmosphere

    @elsewhere_atmosphere.setter
    def elsewhere_atmosphere(self, path):
        """ Implement if required. """
        pass

    @property
    def spot(self):
        """ Get the instance of :class:`~.Spot.Spot`. """
        return self._spot

    @property
    def elsewhere(self):
        """ Get the instance of :class:`~.Elsewhere.Elsewhere`. """
        return self._elsewhere

    def embed(self, spacetime, p, fast_total_counts, threads):
        """ Embed the photosphere in an ambient Schwarzschild spacetime.

        In other words, generate a discrete representation of the photospheric
        radiation field and the null mapping from the photosphere to infinity,
        for use in flux integrators called by distant observers.

        :param obj spacetime: Instance of :class:`~.Spacetime.Spacetime`.

        :param list p: Parameters of photospheric radiation field.

        Parameter vector:

        * ``p[i]`` = coordinate mode/oscillation frequency, *if not locked*
        * ``p[i:j]`` = parameters for photospheric spot(s)
        * ``p[j:]`` = parameters for radiation field elsewhere

        where ``i = self._num_params`` and ``j = 1 + self._spot.num_params``.

        """
        i = self._num_params

        if self._elsewhere is not None:
            j = i + self._spot.num_params
            self._spot.embed(spacetime, p[i:j],
                             fast_total_counts,
                             threads,
                             self._elsewhere.eval_srcRadFieldParamVectors,
                             p[j:])
            self._elsewhere.embed(spacetime, p[j:], threads)
        else:
            self._spot.embed(spacetime, p[i:],
                             fast_total_counts,
                             threads)

        self._spacetime = spacetime # store a reference to the spacetime object

    def integrate(self, energies, threads):
        """ Integrate over the photospheric radiation field.

        :param energies: A one-dimensional :class:`numpy.ndarray` of energies
                         in keV.

        :param int threads: Number of ``OpenMP`` threads to spawn for pulse
                            integration.

        """
        if self._elsewhere is not None:
            time_invariant = self._elsewhere.integrate(self._spacetime,
                                                       energies,
                                                       threads,
                                                       *self._elsewhere_atmosphere)

        self._pulse = self._spot.integrate(self._spacetime,
                                           energies,
                                           threads,
                                           self._spot_atmosphere,
                                           self._elsewhere_atmosphere)

        # add time-invariant component to first time-dependent component
        if self._elsewhere is not None:
            for i in range(self._pulse[0].shape[1]):
                self._pulse[0][:,i] += time_invariant

    @property
    def pulse(self):
        """ Get the stored pulse.

        :return: A :class:`numpy.ndarray` of size ``m x n``, where ``m`` is
                 the number of energies, and ``n`` is the number of phases.
                 Units are photon/s/keV; the distance is a fast parameter so
                 the area units are not yet factored in.

        """
        return self._pulse





