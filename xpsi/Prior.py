from __future__ import division, print_function

__all__ = ["Prior", "PriorError", "BoundsError"]

from . import _comm, _size, _rank # of MPI.COMM_WORLD

from .global_imports import *
from . import global_imports

from . import make_verbose

from abc import ABCMeta, abstractmethod

class Prior(object):
    """ The joint prior distribution.

    Methods to both evaluate the distribution (required by MCMC) and
    draw uniformly from the distribution (required for nested sampling, and
    by default is used to initialise an ensemble of MCMC chains).

    The distribution must be integrable (proper) and is thus
    *usually* bounded (compactly supported on a space).

    :param list bounds: The set of 2-tuples of hard parameter bounds. The
                        list has length :math:`d`, the dimensionality of
                        the space :math:`\mathbb{R}^{d}`. Used to rapidly
                        evaluate for finiteness.

    .. note:: If you wish to check bounds manually, implement
              the ``__init__()`` and ``__call__()`` methods, and do not
              access the default code (e.g., via ``super()``).

    """
    __metaclass__ = ABCMeta

    @abstractmethod
    def __init__(self, bounds):
        self._bounds = bounds

    @property
    def ndims(self):
        """ Get the dimensionality of the parameter space. """
        return len(self._bounds)

    @abstractmethod
    def __call__(self, p):
        """ Evaluate distribution at :obj:`p` and store it as a property.

        :param list p: Vector of model parameter values.

        """
        for i, b in enumerate(self._bounds):
            if None not in b:
                if not b[0] <= p[i] <= b[1]:
                    return -_np.inf

        return 0.0

    @abstractmethod
    def inverse_sample(self, hypercube):
        """ Draw sample uniformly from the distribution via inverse sampling.

        By default, implements a flat density between bounds. If ``None`` is
        in a tuple of bounds, ``None`` is assigned to the corresponding
        coordinate and the user must handle in a custom subclass.

        :param iterable hypercube:
            A pseudorandom point in an n-dimensional hypercube.

        :return: A parameter *list*.

        """
        p = [0.0] * len(hypercube)

        for i, b in enumerate(self._bounds):
            if None not in b:
                p[i] = b[0] + (b[1] - b[0]) * hypercube[i]
            else:
                p[i] = None

        return p

    def inverse_sample_and_transform(self, hypercube):
        """ Inverse sample and then transform.

        This method is useful for drawing from the prior when overplotting
        posterior and prior density functions.

        """

        p = self.transform(self.inverse_sample(hypercube))

        return p

    @staticmethod
    def transform(p):
        """ A transformation for post-processing. """

        if not isinstance(p, list):
            p = list(p)

        raise NotImplementedError('Define a transformation.')

    @make_verbose('Drawing samples from the joint prior','Samples drawn')
    def draw(self, ndraws, transform=False):
        """ Draw samples uniformly from the prior via inverse sampling.

        :param ndraws: Number of draws.
        :type ndraws: int
        :return: (samples, acceptance fraction)
        :rtype: (*ndarray[ndraws, self.ndims]*, *float*)

        """

        h = _np.random.rand(int(ndraws), self.ndims)

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
                                dtype=_np.double)
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
                h = _np.random.rand(int(redraw)+1, self.ndims)
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
                self.estimate_hypercube_frac()
            except AttributeError:
                print('Cannot locate method for estimating fraction.')
            else:
                return self._unit_hypercube_frac
