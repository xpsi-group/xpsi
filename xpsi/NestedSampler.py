
from xpsi import _verbose
from xpsi.utils import make_verbose

from xpsi.Likelihood import Likelihood
from xpsi.Prior import Prior

try:
    import pymultinest
except ImportError:
    print('Check your PyMultiNest installation.')
    raise
else:
    if _verbose:
        print('Imported PyMultiNest.')

class NestedSampler(object):
    """ Extended MultiNest wrapper.

    :param int ndims: Number of parameters.

    :param likelihood: An instance of :class:`~.Likelihood.Likelihood`.

    :param prior: An instance of :class:`~.Prior.Prior`.

    :param kwargs: A dictionary of MultiNest runtime settings, to be
                   handled via the PyMultiNest wrapper.

    """
    def __init__(self,
                 ndims,
                 likelihood,
                 prior):

        try:
            self._ndims = int(ndims)
        except TypeError:
            raise TypeError('Dimensionality must be an integer.')

        if not isinstance(likelihood, Likelihood):
            raise TypeError('Invalid type for likelihood object.')
        else:
            self._likelihood = likelihood

        if not isinstance(prior, Prior):
            raise TypeError('Invalid type for prior object.')
        else:
            self._prior = prior

    @make_verbose('Commencing integration', 'Integration completed')
    def __call__(self, **kwargs):
        """ Integrate.

        :param kwargs:
            Keyword arguments passed to :func:`pymultinest.solve`.

        """

        LHS_seed = LHS_seed if LHS_seed is not None else None
        kwargs.setdefault('sampling_efficiency', 0.8)
        kwargs['sampling_efficiency'] /= self._prior.unit_hypercube_frac(LHS_seed=LHS_seed)

        yield 'Sampling efficiency set to: %.4f.'%kwargs['sampling_efficiency']

        # ignore the pymultinest output object
        # it has never been used in the author's workflows but change this
        # it is useful to you
        _ = pymultinest.solve(self._likelihood,
                              self._prior.inverse_sample,
                              self._ndims,
                              log_zero = self._likelihood.llzero,
                              **kwargs)
