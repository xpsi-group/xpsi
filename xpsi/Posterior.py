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

    A callable instance is required by `emcee <http://dfm.io/emcee/current/>`_.

    """
    def __init__(self,
                 likelihood,
                 prior,
                 source_pulse_blobs = False,
                 **kwargs):
        """
        :param likelihood: An instance of :class:`~.Likelihood.Likelihood`.

        :param prior: An instance of :class:`~.Prior.Prior`.

        :param bool source_pulse_blobs: Store ``emcee`` blobs with the
                                        incident source pulses.

        """
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

        try:
            assert isinstance(source_pulse_blobs, bool)
        except AssertionError:
            self._source_pulse_blobs = False
        else:
            self._source_pulse_blobs = source_pulse_blobs

            if source_pulse_blobs:
                self._blobs = [None] * len(self._likelihood.pulses)

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
            raise PriorError('Log-prior is ``NaN``.')

        if _np.isfinite(lp):
            try:
                ll = self._likelihood(p)
            except LikelihoodError:
                print('Parameter vector: ', p)
                raise
        else:
            ll = 0.0

        if self._source_pulse_blobs:
            for i, pulse in enumerate(self._likelihood.pulses):
                self._blobs[i] = pulse.derived
            return ll + lp, self._blobs
        else:
            return ll + lp








