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

class UltranestSampler(ultranest.ReactiveNestedSampler):
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
        
        self._param_names = self._likelihood.names

        # initialise sampler 
        super().__init__(param_names=self._param_names, 
                         loglike=self.my_likelihood, 
                         transform=self._prior.inverse_sample, 
                         **sampler_params)

    def __call__(self, runtime_params):
        """ Start the sampling. --> Say what the output is 
        
        :param runtime_params: Keyword arguments passed passed to :func:`run`.
        
        """

        # run sampler with given runtime params
        self.run(**runtime_params)

    def my_likelihood(self, params):
        """Create a non-xpsi likelihood object that ultranest understands. """

        arg1, *args = params

        # calculate the log-likelihood
        ultranest_likelihood = self._likelihood(p=[arg1, *args], reinitialise=True)

        return ultranest_likelihood
