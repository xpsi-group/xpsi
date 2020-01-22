from __future__ import division, print_function

from .global_imports import *
from . import global_imports

from .Spacetime import Spacetime
from .HotRegion import HotRegion
from .Elsewhere import Elsewhere

from .Parameter import Parameter
from .ParameterSubspace import ParameterSubspace

class Photosphere(ParameterSubspace):
    """ A photosphere embedded in an ambient Schwarzschild spacetime.

    :param obj hot:
        An instance of :class:`~.HotRegion.HotRegion` (or a
        derived class). This objects represents the hot
        regions of the surface that in most use-cases will be
        assumed to contain radiating material that is hotter
        than that *elsewhere*.

    :param obj elsewhere:
        An instance of :class:`~.Elsewhere.Elsewhere` (or a derived class).

    :param dict bounds:
        Bounds are supplied for instantiation of a frequency parameter.
        The parameter name ``'mode_frequency'`` must be a key in the
        dictionary unless the parameter is *fixed* or *derived*. If a bound
        is ``None`` that bound is set equal to a strict hard-coded bound.
        If ``None``, lock the coordinate rotation frequency of a mode of
        asymmetry in the photosphere to a fixed frequency, e.g., the stellar
        rotation frequency. If bounds are passed, the frequency is interpreted
        as a free parameter.

    :param dict values:
        Either the fixed value of the mode frequency, a callable if the
        frequency is *derived*, or a value upon initialisation if the
        frequency is free. The dictionary must have a key with name
        ``'mode_frequency'`` if it is *fixed* or *derived*.
        If the asymmetry is locked to the stellar spin, then you need to pass
        the spin frequency. If fixed but different to the spin frequency, this
        value needs to be passed instead. In the hot region base class this
        mode frequency is applied to normalise the ray lags instead of the
        stellar rotation frequency.

    .. note::

        In basic modelling patterns the frequency is the spin frequency,
        and thus you only need to explicitly pass the spin as ``value`` whilst
        leaving ``bounds`` to default. If the spin frequency happens to be a
        free parameter (perhaps with informative prior information), then
        pass a callable instead that can be used to get the spin frequency
        dynamically when the derived mode frequency variable is called for.

    """
    required_names = ['mode_frequency']

    def __init__(self,
                 hot = None, elsewhere = None,
                 bounds = {}, values = {},
                 **kwargs):

        if elsewhere is not None:
            if not isinstance(elsewhere, Elsewhere):
                raise TypeError('Invalid type for an elsewhere object.')
            else:
                self._elsewhere = elsewhere
        else:
            self._elsewhere = None
            if hot is None:
                raise ValueError('The photosphere must radiate.')

                                          # including derived classes
        if hot is not None and hot is not isinstance(hot, HotRegion):
            if hasattr(hot, 'objects'):
                for obj in getattr(hot, 'objects'):
                    if not isinstance(obj, HotRegion):
                        raise TypeError('Invalid object for the hot '
                                        'region(s).')
            else:
                raise TypeError('Invalid object for the hot region(s).')

        self._hot = hot

        self._hot_atmosphere = ()
        self._elsewhere_atmosphere = ()

        doc = """
        Coordinate frequency of the mode of radiative asymmetry in the
        photosphere that is assumed to generate the pulsed signal [Hz].
        """
        mode_frequency = Parameter('mode_frequency',
                                   strict_bounds = (0.0, 2000.0),
                                   bounds = bounds.get('mode_frequency', None),
                                   doc = doc,
                                   symbol = r'$f_{\rm mode}$',
                                   value = values.get('mode_frequency', None))

        super(Photosphere, self).__init__(mode_frequency,
                                          hot, elsewhere, **kwargs)

    @property
    def hot_atmosphere(self):
        """ Get the numerical atmosphere buffers for hot regions if used.

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
        ``self._hot_atmosphere = (logT, logg, mu, logE, buf)``, where:
        ``logT`` is a logarithm of local comoving effective temperature;
        ``logg`` is a logarithm of effective surface gravity;
        ``mu`` is the cosine of the angle from the local surface normal;
        ``logE`` is a logarithm of the photon energy; and
        ``buf`` is a one-dimensional buffer of intensities of size given by
        the product of sizes of the first :math:`n-1` tuple elements.

        It is highly recommended that buffer preloading is used, instead
        of loading from disk in the customisable radiation field extension
        module, to avoid reading from disk for every signal
        (likelihood) evaluation. This can be a non-negligible waste of compute
        resources. By preloading in Python, the memory is allocated and
        references to that memory are not in general deleted until a sampling
        script exits and the kernel stops. The likelihood callback accesses
        the same memory upon each call without I/O.

        """
        return self._hot_atmosphere

    @hot_atmosphere.setter
    def hot_atmosphere(self, path):
        """ Implement if required. """
        raise NotImplementedError('Implement setter if required.')

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
        ``self._hot_atmosphere = (logT, logg, mu, logE, buf)``, where:
        ``logT`` is a logarithm of local comoving effective temperature;
        ``logg`` is a logarithm of effective surface gravity;
        ``mu`` is the cosine of the angle from the local surface normal;
        ``logE`` is a logarithm of the photon energy; and
        ``buf`` is a one-dimensional buffer of intensities of size given by
        the product of sizes of the first :math:`n-1` tuple elements.

        It is highly recommended that buffer preloading is used, instead
        of loading from disk in the customisable radiation field extension
        module, to avoid reading from disk for every signal
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
        raise NotImplementedError('Implement setter if required.')

    @property
    def hot(self):
        """ Get the instance of :class:`~.HotRegion.HotRegion`. """
        return self._hot

    @property
    def elsewhere(self):
        """ Get the instance of :class:`~.Elsewhere.Elsewhere`. """
        return self._elsewhere

    @property
    def spacetime(self):
        """ Return instance of :class:`~.Spacetime.Spacetime`. """
        return self._spacetime

    @spacetime.setter
    def spacetime(self, obj):
        if not isinstance(obj, Spacetime):
            raise TypeError('Invalid type for spacetime object.')
        # otherwise store a reference to the spacetime object
        self._spacetime = obj

    def embed(self, fast_total_counts, threads):
        """ Embed the photosphere in an ambient Schwarzschild spacetime.

        In other words, generate a discrete representation of the photospheric
        radiation field and the null mapping from the photosphere to infinity,
        for use in flux integrators called by distant observers.

        """
        if self._elsewhere is not None:
            self._elsewhere.embed(self._spacetime, threads)

            if self._hot is not None:
                self._hot.embed(self._spacetime,
                                self,
                                fast_total_counts,
                                threads,
                                self._elsewhere._compute_cellParamVecs)
        else:
            self._hot.embed(self._spacetime,
                            self,
                            fast_total_counts,
                            threads)

    def integrate(self, energies, threads):
        """ Integrate over the photospheric radiation field.

        :param energies:
            A one-dimensional :class:`numpy.ndarray` of energies in keV.

        :param int threads:
            Number of ``OpenMP`` threads to spawn for pulse integration.

        """
        if self._elsewhere is not None:
            if isinstance(energies, tuple): # resolve energy container type
                if not isinstance(energies[0], tuple):
                    _energies = energies[0]
                else:
                    _energies = energies[0][0]
            else:
                _energies = energies
            self._time_invariant = self._elsewhere.integrate(self._spacetime,
                                                   _energies,
                                                   threads,
                                                   *self._elsewhere_atmosphere)

        if self._hot is not None:
            self._pulse = self._hot.integrate(self._spacetime,
                                              energies,
                                              threads,
                                              self._hot_atmosphere,
                                              self._elsewhere_atmosphere)

            if not isinstance(self._pulse[0], tuple):
                self._pulse = (self._pulse,)

            # add time-invariant component to first time-dependent component
            if self._elsewhere is not None:
                for i in range(self._pulse[0][0].shape[1]):
                    self._pulse[0][0][:,i] += self._time_invariant

    @property
    def pulse(self):
        """ Get the stored pulse.

        :return: *ndarray[m,n]*, where :math:`m` is
                 the number of energies, and :math:`n` is the number of phases.
                 Units are photon/s/keV; the distance is a fast parameter so
                 the areal units are not yet factored in.

        """
        return self._pulse

    @property
    def time_invariant(self):
        """ Get the stored time-invariant signal.


        :returns:
            *ndarray[n]* containing the time-invariant signal at
            :math:`n` energies. Units are photon/s/keV. The distance is a
            fast parameter so the areal units are not yet factored in.

        """
        return self._time_invariant

Photosphere._update_doc()
