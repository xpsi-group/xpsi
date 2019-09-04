from __future__ import division, print_function

from .global_imports import *
from . import global_imports

from abc import ABCMeta, abstractmethod

class BoundsError(Exception):
    """
    Raised if an implementation of the abstract initialiser of
    :class:`~.ParameterSubspace.ParameterSubspace` passes an invalid bounds
    object to the abstract initialiser (e.g., via ``super()``).

    Can also be imported into another module and raised if there is an error
    on one of the bounds.

    """

class ParameterSubspace(object):
    """ Abstract parameter subspace. """

    __metaclass__ = ABCMeta

    @abstractmethod
    def __init__(self, num_params, bounds):
        """
        :param int num_params: The number of parameters in the subspace of
                               :math:`\mathbb{R}^{d}`.

        :param list bounds: One 2-tuple of hard bounds per parameter. Can 
                            be unbounded *in principle*, but read the
                            documentation for the :class:`~.Prior.Prior` class
                            first.

        """
        try:
            self._num_params = int(num_params)
        except TypeError:
            raise TypeError('The number of parameters must be an integer.')

        try:
            assert isinstance(bounds, list)
            assert len(bounds) == self._num_params
            for b in bounds:
                assert len(b) == 2
                assert b[0] < b[1]
        except AssertionError:
            raise BoundsError('The bounds must be passed as a list of lists or '
                              'a list of tuples, where each list or tuple has '
                              'a length of two, and a first element less than '
                              'the second element.')
        else:
            self._bounds = bounds

    @property
    def num_params(self):
        """ Get the number of dimensions of the parameter subspace. """
        return self._num_params

    @property
    def bounds(self):
        """ Get a list of tuples of hard parameter bounds. """
        return self._bounds






