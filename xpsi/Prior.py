__all__ = ["Prior"]

from xpsi import _comm, _size, _rank # of MPI.COMM_WORLD

from xpsi.global_imports import *

from xpsi.utils import make_verbose

from abc import ABCMeta, abstractmethod
from xpsi.ParameterSubspace import ParameterSubspace

class Prior(ParameterSubspace, metaclass=ABCMeta):
    """ The joint prior distribution of parameters (including hyperparameters).

    Methods to both evaluate the distribution (required by MCMC) and
    draw uniformly from the distribution (required for nested sampling, and
    by default is used to initialise an ensemble of MCMC chains).

    The distribution must be integrable (proper) and is thus
    *usually* bounded (compactly supported on a space).

    In the ``__call__()`` method the parameter values are checked against the
    hard parameter bounds of the prior support.

    .. note::

        If you wish to check bounds manually, implement the ``__init__()`` and
        ``__call__()`` methods, and do not access the default code (e.g., via
        ``super()``) or only use the default code for some parameters.

    :param obj parameters:
        An optional instance of :class:`~.ParameterSubspace.ParameterSubspace`.

    :param obj hyperparameters:
        Positional arguments that are hyperparameters (parameters of the prior
        distribution), or iterables of hyperparameters.

    """
    __derived_names__ = None

    __draws_from_support__ = 5

    def __init__(self, parameters = None, *hyperparameters):
        """
        You might want to overwrite this initialiser to do some custom
        setup for priors whose handling is more involved.

        If you overwrite this initialiser, you do not need to take parameters
        as an argument if you prefer not to. When the instance is passed to
        a likelihood object, the likelihood object will use itself to store
        a reference to the parameters to ensure interoperability.

        """
        if parameters is not None:
            self.parameters = parameters

        super(Prior, self).__init__(*hyperparameters)

        for hyperparameter in self:
            if not hyperparameter.is_hyperparameter:
                raise ValueError('Parameter is not declared as a '
                                 'hyperparameter but is a parameter of the '
                                 'prior distribution.')

    def __len__(self):
        """ Number of parameter + hyperparameter dimensions. """
        return len(self._parameters) # redirect, overwriting base class

    @property
    def parameters(self):
        """ Get parameter subspace object. """
        return self._parameters

    @parameters.setter
    def parameters(self, obj):
        if not isinstance(obj, ParameterSubspace):
            raise TypeError('A ParameterSubspace object is required.')

        self._parameters = obj

    @abstractmethod
    def __call__(self, p = None):
        """ Evaluate distribution at :obj:`p` and store it as a property.

        :param list p:
            Vector of model parameter values, but typically unused. If you
            use it, handle it in a custom implementation of this method.

        """
        for param in self.parameters:
            if param.bounds[0] is not None:
                if not param.bounds[0] <= param.value:
                    return -_np.inf
            if param.bounds[1] is not None:
                if not param.value <= param.bounds[1]:
                    return -_np.inf

        return 0.0

    @abstractmethod
    def inverse_sample(self, hypercube = None):
        """ Draw sample uniformly from the distribution via inverse sampling.

        By default, implements a flat density between bounds. If ``None`` is
        in a tuple of bounds, ``None`` is assigned to the corresponding
        coordinate and the user must handle in a custom subclass.

        :param iterable hypercube:
            A pseudorandom point in an n-dimensional hypercube.

        :returns: A parameter *list*.

        .. note::

            If you call this base method via the ``super()`` built-in, the
            current parameter values will be cached when new values are
            assigned. If you then assign to a parameter again, the current
            value will be automatically cached, thus overwriting the cache
            established in the body of this present method. If you want to
            call this present method and then assign again to a parameter,
            you can restore the cached value so that it is pushed to the
            cache when you reassign.

        """
        if hypercube is None:
            hypercube = _np.random.rand(len(self))

        p = [None] * len(self)

        for i, param in enumerate(self.parameters):
            if None not in param.bounds:
                b = param.bounds
                param.value = b[0] + (b[1] - b[0]) * hypercube[i]
            else:
                raise ValueError('Compact support required.')

            p[i] = param.value

        return p

    def inverse_sample_and_transform(self, hypercube = None):
        """ Inverse sample and then transform.

        This method is useful for drawing from the prior when overplotting
        posterior and prior density functions.

        """

        p = self.transform(self.inverse_sample(hypercube))

        return p

    @staticmethod
    def transform(p, **kwargs):
        """ A transformation for post-processing.

        A subclass can implement this attribute as an instance method: it
        does not need to be a static method.

        :param list p:
            Parameter vector.
        :kwargs:
            A key sometimes passed is ``old_API``, which flags whether a
            transformation needs to account for a parameter vector written to
            file by an older software version, which might be different in
            due to transformations of parameters defined in the current
            software version.

        .. note::

            As an example, the `Riley et al. 2019 <https://ui.adsabs.harvard
            .edu/abs/2019ApJ...887L..21R/abstract>`_ samples are in inclination
            :math:`i` instead of :math:`\cos(i)` which is the current inclination
            parameter in the API. Therefore the transformation needed depends
            on the source of the parameter vector. If the vector is from the
            original sample files, then it needs to be transformed to have
            the same parameter definitions as the current API. However, when
            drawing samples from the prior in the current API, no such
            transformation needs to be performed because these are the
            definitions we need to match. Refer to the dummy example code
            in the method body.

        :returns: Transformed vector ``p`` where ``len(p) > len(self)``.
        :rtype: *list*

        """
        # it is suggested that you copy like this when you overwrite
        # to make a mutable container that can be appended to and returned
        p = list(p)

        # example transformation of parameters to match API definitions if
        # this method is bound
        if kwargs.get(old_API, False):
            idx = self.parameters.index('cos_inclination')
            p[idx] = math.cos(p[idx])

        # used ordered names and values if this method is bound
        ref = dict(list(zip(self.parameters.names, p)))

        raise NotImplementedError('Define a transformation.')

        return p

    @property
    def derived_names(self):
        if self.__derived_names__ is None:
            try:
                self.transform(self.inverse_sample())
            except NotImplementedError:
                pass
            else:
                raise AttributeError('A transformation has been implemented '
                                     'in class %s, but no names have been '
                                     'declared for the derived parameters.'
                                     % type(self).__name__)

        return self.__derived_names__

    def index(self, name):
        return len(self) + self.derived_names.index(name)

    @make_verbose('Drawing samples from the joint prior','Samples drawn')
    def draw(self, ndraws, transform=False):
        """ Draw samples uniformly from the prior via inverse sampling.

        :param int ndraws: Number of draws.
        :return: (samples, acceptance fraction)
        :rtype: (*ndarray[ndraws, len(self)]*, *float*)

        """

        h = _np.random.rand(int(ndraws), len(self))

        if transform:
            try:
                p = self.inverse_sample_and_transform(h[0, :])
            except AttributeError:
                yield 'Cannot transform... falling back to default...'
                p = self.inverse_sample(h[0, :])
        else:
            p = self.inverse_sample(h[0, :])

        try:
            samples = _np.zeros((int(ndraws), len(p)),
                                dtype = _np.double)
        except TypeError:
            yield 'An error occurred when handling the output of a custom '\
                  'method'

        finite_counter = counter = index = 0
        while finite_counter < ndraws:
            try:
                if transform:
                    try:
                        p = self.inverse_sample_and_transform(h[index, :])
                    except NotImplementedError:
                        p = self.inverse_sample(h[index, :])
                else:
                    p = self.inverse_sample(h[index, :])
            except IndexError: # use estimate to draw more from hypercube
                if finite_counter > 0:
                    redraw = float(counter) * ndraws / finite_counter
                else:
                    redraw = ndraws
                redraw -= finite_counter
                h = _np.random.rand(int(redraw)+1, len(self))
                index = 0
                if transform:
                    try:
                        p = self.inverse_sample_and_transform(h[index, :])
                    except NotImplementedError:
                        p = self.inverse_sample(h[index, :])
                else:
                    p = self.inverse_sample(h[index, :])

            if _np.isfinite(self.__call__(p)):
                samples[finite_counter,:] = p
                finite_counter += 1
            counter += 1
            index += 1

        yield samples, float(finite_counter) / counter

    @make_verbose('Estimating fractional hypervolume of the unit hypercube '
                  'with finite prior density:',
                  'Fractional hypervolume estimated')
    def estimate_hypercube_frac(self, ndraws=5):
        """
        Estimate using Monte Carlo integration the fractional hypervolume
        within a unit hypercube at which prior density is finite.

        :param optional[int] ndraws:
            Base-10 logarithm of number of draws from the prior to require
            at which the density is finite.

        """
        if _rank == 0:
            ndraws = 10**ndraws
            yield ('Requiring %.E draws from the prior support '
                   'for Monte Carlo estimation' % ndraws)

            self._unit_hypercube_frac = self.draw(ndraws)[1]
        else:
            self._unit_hypercube_frac = None

        if _size > 1:
            self._unit_hypercube_frac = _comm.bcast(self._unit_hypercube_frac,
                                                    root=0)

        yield ('The support occupies an estimated '
               '%.1f%% of the hypervolume within the unit hypercube'
               % (self._unit_hypercube_frac*100.0))

        yield self._unit_hypercube_frac

    @property
    def unit_hypercube_frac(self):
        """ Get the fractional hypervolume with finite prior density. """

        try:
            return self._unit_hypercube_frac
        except AttributeError:
            try:
                self.estimate_hypercube_frac(self.__draws_from_support__)
            except AttributeError:
                print('Cannot locate method for estimating fraction.')
            else:
                return self._unit_hypercube_frac
