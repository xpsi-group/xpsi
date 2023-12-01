from xpsi import _verbose
from xpsi.Likelihood import Likelihood
from xpsi.Prior import Prior

import numpy as np
import os

try:
    import ultranest
    import ultranest.stepsampler
except ImportError:
    print('Check your UltraNest installation.')
    raise
else:
    if _verbose:
        print('Imported UltraNest.')

class UltranestSampler(ultranest.ReactiveNestedSampler):
    """ Wrapper for the Ultranest sampler (from https://johannesbuchner.github.io/UltraNest/ultranest.html)

    :param likelihood: An instance of :class:`~.Likelihood.Likelihood`.

    :param prior: An instance of :class:`~.Prior.Prior`.

    :param sampler_params: A dictionary of the keyword arguments passed to the instance of :class:`~.ultranest.ReactiveNestedSampler`.

    :param use_stepsampler: Boolean indicating if the stepsampler is used. 

    :param stepsampler_params: A dictionary of the keyword arguments passed to the stepsampler :`~.ultranest.stepsampler.SliceSampler`.

    """

    def __init__(self, 
                 likelihood,
                 prior,
                 sampler_params,
                 use_stepsampler,
                 stepsampler_params):

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
        if use_stepsampler: 
            # set default stepsampler parameters if not specified 
            stepsampler_params.setdefault('nsteps', 2*len(self._param_names))
            stepsampler_params.setdefault('generate_direction', ultranest.stepsampler.generate_mixture_random_direction)
            
            # initialise step sampler 
            self.stepsampler = ultranest.stepsampler.SliceSampler(**stepsampler_params)

    def __call__(self, runtime_params):
        """ Run the sampler until target convergence criteria are fulfilled with the 
        given runtime parameters. 
        
        :param runtime_params: A dictionary of the keyword arguments passed to :func:`run`.
        
        """

        # run sampler with given runtime params
        self.run(**runtime_params)

    def my_likelihood(self, params):
        """Calculate the loglikelihood value for a given set of parameter values.  
        
        :param params: List of parameter values. 

        :returns: Float with loglikelihood value for given set of parameter values. 

        """

        arg1, *args = params

        # calculate the log-likelihood
        ultranest_likelihood = self._likelihood(p=[arg1, *args], reinitialise=True)

        return ultranest_likelihood
    
    def write_results(self, sampler_params, out_filename):
        """ Get output in txt file with columns containing weights, -2*loglikelihood, 
        and parameters (this is the format required for post-processing within X-PSI). 

        :param sampler_params: A dictionary of the keyword arguments passed to the instance of :class:`~.ultranest.ReactiveNestedSampler`.

        :param out_filename: String of the output filename. 

        """
        # extract results
        data = np.array(self.results["weighted_samples"]["points"])
        weights = np.array(self.results["weighted_samples"]["weights"])
        logl = np.array(self.results["weighted_samples"]["logl"])

        output = np.column_stack((weights, -2*logl, data))

        # set default sampler parameters if not specified
        sampler_params.setdefault('log_dir', 'results/')
        
        # check if directory exists, otherwise create one
        log_dir = os.makedirs(sampler_params['log_dir'], exist_ok=True)
        file_path = os.path.join(sampler_params['log_dir'], out_filename)

        np.savetxt(file_path, output)
