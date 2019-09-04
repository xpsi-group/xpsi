from __future__ import division, print_function

from .global_imports import *
from . import global_imports

from abc import ABCMeta, abstractproperty

class ParameterSpace(object):
    """ Abstract global model parameter space.

    The parameter space is :math:`\mathbb{R}^{d}`, where :math:`d` is the
    total number of model parameters.

    """

    __metaclass__ = ABCMeta

    @abstractproperty
    def num_params(self):
        """ Get the number of dimensions of the parameter space. """

    @abstractproperty
    def bounds(self):
        """ Get a list of tuples of hard parameter bounds. """
