__all__ = ["Posterior"]

from xpsi.global_imports import *

from xpsi.Likelihood import Likelihood
from xpsi.Signal import LikelihoodError

from xpsi.Prior import Prior
from xpsi.ParameterSubspace import ParameterSubspace
from xpsi.Parameter import StrictBoundsError

class PriorError(ValueError):
    """ Raised if there is a problem with the value of the log-prior. """

class Posterior(object):
    """ The (joint) posterior distribution.

    A callable instance is required by `emcee <http://dfm.io/emcee/current/>`_
    (but is not required for nested sampling).

    :param likelihood:
        An instance of :class:`~.Likelihood.Likelihood`.

    :param prior:
        An instance of :class:`~.Prior.Prior`.

    """
    def __init__(self,
                 likelihood,
                 prior,
                 **kwargs):

        if not isinstance(likelihood, Likelihood):
            raise TypeError('Invalid type for likelihood object.')

        self._likelihood = likelihood
        # no inverse sampling but best to update before prior check
        # and then revert if not in prior support
        self._likelihood.externally_updated = True

        # we do not need to the likelihood object to call the prior for
        # this type of sampling because it is handled in the present class
        # instead; for nested sampling we want to call the prior from the
        # likelihood object in order to check inclusion in the prior support
        del self._likelihood.prior

        if not isinstance(prior, Prior):
            raise TypeError('Invalid type for prior object.')

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

        :returns:
            The logarithm of the posterior density (up to a normalising
            constant), evaluated at :obj:`p`.

        """
        q = self._likelihood.vector # make our own cache here

        try:
            ParameterSubspace.__call__(self._likelihood, p)
        except StrictBoundsError:
            ParameterSubspace.__call__(self._likelihood, q)
            return -_np.inf

        lp = self._prior(p)

        if _np.isnan(lp):
            raise PriorError('Log-prior is not a number.')

        if _np.isfinite(lp):
            try:
                ll = self._likelihood()
            except LikelihoodError:
                print('Parameter vector: ', p)
                raise
        else:
            ParameterSubspace.__call__(self._likelihood, q)
            ll = 0.0

        return ll + lp
