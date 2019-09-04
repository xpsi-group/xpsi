from __future__ import division, print_function

from .global_imports import *
import global_imports

from . import make_verbose, _verbose

from .Likelihood import Likelihood
from .Prior import Prior

try:
    import pymultinest
except ImportError:
    print('Check your PyMultiNest installation.')
    raise
else:
    if _verbose:
        print('Imported PyMultiNest.')

class MultiNestIntegrator(object):
    """ Extended MultiNest wrapper.

    """
    def __init__(self,
                 ndims,
                 likelihood,
                 prior):
        """
        :param int ndims: Number of parameters.

        :param likelihood: An instance of :class:`~.Likelihood.Likelihood`.

        :param prior: An instance of :class:`~.Prior.Prior`.

        :param kwargs: A dictionary of MultiNest runtime settings, to be
                       handled via the PyMultiNest wrapper.

        """

        try:
            self._ndims = int(ndims)
        except TypeError:
            raise TypeError('Dimensionality must be an integer.')

        try:
            assert isinstance(likelihood, Likelihood)
        except AttributeError:
            raise TypeError('Invalid type for likelihood object.')
        else:
            self._likelihood = likelihood

        try:
            assert isinstance(prior, Prior)
        except AttributeError:
            raise TypeError('Invalid type for prior object.')
        else:
            self._prior = prior

    @make_verbose('Commencing integration', 'Integration completed')
    def __call__(self, **kwargs):
        """ Integrate. """

        try:
            kwargs['sampling_efficiency'] /= self._prior.unit_hypercube_frac
        except KeyError:
            kwargs['sampling_efficiency'] = 1.0/self._prior.unit_hypercube_frac

        yield 'Sampling efficiency set to: %.4f' % kwargs['sampling_efficiency']

        # ignore the pymultinest output object
        # it has never been used in the author's workflows but change this
        # it is useful to you
        _ = pymultinest.solve(self._likelihood,
                              self._prior.inverse_sample,
                              self._ndims,
                              log_zero = self._likelihood.llzero,
                              **kwargs)
