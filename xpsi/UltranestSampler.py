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
    """ Wrapper for the UltraNest (https://johannesbuchner.github.io/UltraNest/) package (Buchner 2021). 

    :param likelihood: An instance of :class:`~.Likelihood.Likelihood`.

    :param prior: An instance of :class:`~.Prior.Prior`.

    :param sampler_params: A dictionary of the keyword arguments passed to the instance of :class:`ultranest.ReactiveNestedSampler`.

    :param use_stepsampler: Boolean indicating if the stepsampler is used. In this case the :class:`ultranest.stepsampler.SliceSampler` is used. 

    :param stepsampler_params: A dictionary of the keyword arguments passed to the stepsampler :class:`ultranest.stepsampler.SliceSampler`.

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

        # make default folder for output if not specified
        sampler_params.setdefault('log_dir', 'output')
        os.makedirs(sampler_params['log_dir'], exist_ok=True)

        # initialise (region) sampler
        super().__init__(param_names=self._param_names, 
                         loglike=self.my_likelihood, 
                         transform=self._prior.inverse_sample, 
                         **sampler_params)
                
        # change region sampler to step sampler 
        if use_stepsampler: 
            # set default stepsampler parameters if not specified 
            stepsampler_params.setdefault('nsteps', 2*len(self._param_names))
            stepsampler_params.setdefault('generate_direction', ultranest.stepsampler.generate_mixture_random_direction)
            stepsampler_params.setdefault('check_nsteps','move-distance')
            stepsampler_params.setdefault('log', open(f"{sampler_params['log_dir']}stepsampler_diagnostics.txt", "w") )

            # initialise step sampler 
            self.stepsampler = ultranest.stepsampler.SliceSampler(**stepsampler_params)

    def __call__(self, runtime_params):
        """ Run the sampler until target convergence criteria are fulfilled with the 
        given runtime parameters. 
        
        :param runtime_params: A dictionary of the keyword arguments passed to :func:`run`.
        
        """

        self.run(**runtime_params)

    def my_likelihood(self, params):
        """Calculate the loglikelihood value for a given set of parameter values.  
        
        :param params: List of parameter values. 

        :returns: Float with loglikelihood value for given set of parameter values. 

        """

        arg1, *args = params

        # calculate the log-likelihood
        ultranest_likelihood = self._likelihood(p=[arg1, *args], reinitialise=True)

        if isinstance(ultranest_likelihood, np.ndarray): # hackish way around: need to find better solution for this
            ultranest_likelihood = ultranest_likelihood[0]
        
        return ultranest_likelihood
    
    def write_results(self, sampler_params, out_filename):
        """ Get output txt file with columns containing weights, -2*loglikelihood, 
        and parameters, which is the format required for post-processing within X-PSI. 
        This is additional to output files UltraNest produces.  

        :param sampler_params: A dictionary of the keyword arguments passed to the instance of :class:`ultranest.ReactiveNestedSampler`.

        :param out_filename: String specifying the name of the output file. 

        """
        # extract results
        data = np.array(self.results["weighted_samples"]["points"])
        weights = np.array(self.results["weighted_samples"]["weights"])
        logl = np.array(self.results["weighted_samples"]["logl"])

        # save extra output file special for xpsi post-processing
        output = np.column_stack((weights, -2*logl, data))        
        file_path = os.path.join(sampler_params['log_dir'], out_filename)
        np.savetxt(f"{file_path}.txt", output)
