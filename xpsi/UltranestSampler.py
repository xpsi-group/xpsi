from xpsi import _verbose
from xpsi.Likelihood import Likelihood
from xpsi.Prior import Prior

try:
    import ultranest
except ImportError:
    print('Check your UltraNest installation.')
    raise
else:
    if _verbose:
        print('Imported UltraNest.')

class UltranestSampler():
    """ Initiate Ultranest sampler (from https://johannesbuchner.github.io/UltraNest/ultranest.html)

    :param likelihood: An instance of :class:`~.Likelihood.Likelihood`.

    :param prior: An instance of :class:`~.Prior.Prior`.

    :param sampler_params: Keyword arguments passed instance of :class:`~.ultranest.ReactiveNestedSampler`.

    """

    def __init__(self, 
                 likelihood,
                 prior,
                 sampler_params):

        if not isinstance(likelihood, Likelihood):
            raise TypeError('Invalid type for likelihood object.')
        self._likelihood = likelihood

        if not isinstance(prior, Prior):
            raise TypeError('Invalid type for prior object.')
        self._prior = prior

        # initialise sampler 
        self._sampler = ultranest.ReactiveNestedSampler(param_names=self._likelihood.params, loglike=self._likelihood, transform=self._prior.inverse_sample, **sampler_params)

    def __call__(self, runtime_params):
        """ Start the sampling.
        
        :param runtime_params: Keyword arguments passed passed to :func:`run`.
        """

        # run sampler with given runtime params
        _ = self._sampler.run(**runtime_params)