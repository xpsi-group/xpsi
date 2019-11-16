from __future__ import division, print_function

__all__ = ["Posterior"]

from .global_imports import *
from . import global_imports

from .Likelihood import Likelihood
from .Pulse import LikelihoodError

from .Prior import Prior

class PriorError(ValueError):
    """ Raised if there is a problem with the value of the log-prior. """

class Posterior(object):
    """ The (joint) posterior distribution.

    A callable instance is required by `emcee <http://dfm.io/emcee/current/>`_
    (but is not required for nested sampling).

    :param likelihood: An instance of :class:`~.Likelihood.Likelihood`.

    :param prior: An instance of :class:`~.Prior.Prior`.

    """
    def __init__(self,
                 likelihood,
                 prior,
                 **kwargs):

        if not isinstance(likelihood, Likelihood):
            raise TypeError('Invalid type for likelihood object.')
        else:
            self._likelihood = likelihood

        if not isinstance(prior, Prior):
            raise TypeError('Invalid type for prior object.')
        else:
            self._prior = prior

    @property
    def likelihood(self):
        """ Get the likelihood object. """
        return self._likelihood

    @property
    def prior(self):
        """ Get the prior object. """
        return self._prior

    def __call__(self, p):
        """ Special method to make callable instances to feed to ``emcee``.

        :return: The logarithm of the posterior density (up to a
                 normalising constant), evaluated at :obj:`p`.
        """
        lp = self._prior(p)

        if _np.isnan(lp):
            raise PriorError('Log-prior is not a number.')

        if _np.isfinite(lp):
            try:
                ll = self._likelihood(p)
            except LikelihoodError:
                print('Parameter vector: ', p)
                raise
        else:
            ll = 0.0

        return ll + lp
