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

class UltranestSampler(object):

    def __init__(self, 
                 likelihood,
                 prior):

        if not isinstance(likelihood, Likelihood):
            raise TypeError('Invalid type for likelihood object.')
        else:
            self._likelihood = likelihood

        if not isinstance(prior, Prior):
            raise TypeError('Invalid type for prior object.')
        else:
            self._prior = prior

        # initialise sampler 
        self._sampler = ultranest.ReactiveNestedSampler(param_names=self._likelihood.params, loglike=self._likelihood, transform=self._prior.inverse_sample)

        pass

    def __call__(self, *args, **kwargs):
        # start sampling using region sampler MLFriends as default
        # result = sampler.run()

        pass