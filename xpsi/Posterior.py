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

        # single-point cache for the split log-likelihood/log-prior evaluation
        # used by parallel-tempered sampling (see :meth:`_evaluate_split`)
        self._split_point = None
        self._split_values = (0.0, -_np.inf)

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

    def _evaluate_split(self, p):
        """ Evaluate the log-likelihood and log-prior together at :obj:`p`.

        Parallel-tempered sampling (e.g. `ptemcee`_) requires the prior and
        likelihood as distinct callables rather than the combined posterior in
        :meth:`__call__`, and it evaluates a point by calling :meth:`logprior`
        immediately followed by :meth:`loglikelihood` for the *same* ``p``.

        The parameter vector must be pushed into the likelihood object only
        *once* per point: the X-PSI likelihood detects a change by comparing a
        parameter's value with its cached previous value, and pushing the same
        vector a second time would make them equal again, clearing the
        update flag and leaving the likelihood stale (returning the previously
        computed value). This method therefore does the single update, computes
        both quantities, and caches them so the paired call reuses the result
        instead of updating the vector again.

        .. _ptemcee: https://github.com/willvousden/ptemcee

        :returns: The tuple ``(loglikelihood, logprior)``.

        """
        key = tuple(p)
        if key == self._split_point:
            return self._split_values

        q = self._likelihood.vector # make our own cache here

        try:
            ParameterSubspace.__call__(self._likelihood, p)
        except StrictBoundsError:
            ParameterSubspace.__call__(self._likelihood, q)
            self._split_point = key
            self._split_values = (0.0, -_np.inf)
            return self._split_values

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

        self._split_point = key
        self._split_values = (ll, lp)
        return self._split_values

    def logprior(self, p):
        """ The logarithm of the prior density, evaluated at :obj:`p`.

        Companion to :meth:`loglikelihood` for parallel-tempered sampling; see
        :meth:`_evaluate_split`.

        :returns:
            The logarithm of the prior density (up to a normalising constant),
            or minus infinity if :obj:`p` lies outside the strict parameter
            bounds.

        """
        return self._evaluate_split(p)[1]

    def loglikelihood(self, p):
        """ The logarithm of the likelihood function, evaluated at :obj:`p`.

        Companion to :meth:`logprior` for parallel-tempered sampling; see
        :meth:`_evaluate_split`.

        :returns:
            The logarithm of the likelihood. ``ptemcee`` only requires it where
            the log-prior is finite.

        """
        return self._evaluate_split(p)[0]
