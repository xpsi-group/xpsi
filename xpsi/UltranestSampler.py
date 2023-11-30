from xpsi import _verbose
from xpsi.Likelihood import Likelihood
from xpsi.Prior import Prior

import numpy as np
import os

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
                 sampler_params,
                 step):

        if not isinstance(likelihood, Likelihood):
            raise TypeError('Invalid type for likelihood object.')
        self._likelihood = likelihood

        if not isinstance(prior, Prior):
            raise TypeError('Invalid type for prior object.')
        self._prior = prior
        
        self._param_names = self._likelihood.names

        #  avoid getting output from ultranest since different format is required for xpsi
        adjusted_sampler_params = sampler_params.copy()
        adjusted_sampler_params['log_dir'] = None

        # initialise (region) sampler
        super().__init__(param_names=self._param_names, 
                         loglike=self.my_likelihood, 
                         transform=self._prior.inverse_sample, 
                         **adjusted_sampler_params)
        
        # change region sampler to step sampler 
        if step: 
            nsteps = 2*len(self._param_names)
            self.stepsampler = ultranest.stepsampler.SliceSampler(
                nsteps=nsteps,
                generate_direction=ultranest.stepsampler.generate_mixture_random_direction,
                # adaptive_nsteps=False,
                # max_nsteps=400
                )

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
    
    def write_results(self, sampler_params, out_filename):

        # extract results
        data = np.array(self.results["weighted_samples"]["points"])
        weights = np.array(self.results["weighted_samples"]["weights"])
        logl = np.array(self.results["weighted_samples"]["logl"])

        print("sampler params 2", sampler_params)

        if 'log_dir' in sampler_params is not None:
            log_dir = sampler_params['log_dir']
        else:
            log_dir = os.mkdir("results/")
        
        file_path = os.path.join(log_dir, out_filename)

        output = np.column_stack((weights, -2*logl, data))

        np.savetxt(file_path, output, delimiter=' ')
