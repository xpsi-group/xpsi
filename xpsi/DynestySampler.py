from xpsi import _verbose
from xpsi.Likelihood import Likelihood
from xpsi.Prior import Prior

import numpy as np
import os

try:
    import dynesty
    from dynesty.pool import Pool
except ImportError:
    print('Check your Dynesty installation.')
    raise
else:
    if _verbose:
        print('Imported Dynesty.')

class DynestySampler(dynesty.DynamicNestedSampler):
    """  Wrapper for Dynesty (https://dynesty.readthedocs.io/en/stable/) package.

    :param likelihood: An instance of :class:`~.Likelihood.Likelihood`.

    :param prior: An instance of :class:`~.Prior.Prior`.

    :param sampler_params: A dictionary of the keyword arguments passed to the instance of :class:`dynesty.DynamicNestedSampler`.
    
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

        # initialise dynamic nested sampler 
        super().__init__(loglikelihood=self.my_likelihood, 
                         prior_transform=self._prior.inverse_sample, 
                         ndim=len(self._likelihood),
                         **sampler_params)
            
    def __call__(self, runtime_params):
        """ Run the sampler until target convergence criteria are fulfilled with the 
        given runtime parameters. 
        
        :param runtime_params: A dictionary of the keyword arguments passed to :func:`run`.
        
        """
    
        self.run_nested(**runtime_params)
        

    def my_likelihood(self, params):
        """Calculate the loglikelihood value for a given set of parameter values.  
        
        :param params: List of parameter values. 

        :returns: Float with loglikelihood value for given set of parameter values. 

        """

        arg1, *args = params

        # calculate the log-likelihood
        dynesty_likelihood = self._likelihood(p=[arg1, *args], reinitialise=True)
        
        return dynesty_likelihood
    
    def write_results(self, results, out_filename):
        """ Get output txt file with columns containing weights, -2*loglikelihood, 
        and parameters, which is the format required for post-processing within X-PSI. 

        :param sampler_params: A dictionary of the keyword arguments passed to the instance of :class:`ultranest.ReactiveNestedSampler`.

        :param out_filename: String specifying the name of the output file. 

        """
        # extract results
        data = results.samples
        weights = results.importance_weights()
        logl = results.logl

        # save extra output file special for xpsi post-processing
        output = np.column_stack((weights, -2*logl, data))        
        # file_path = os.path.join(sampler_params['log_dir'], out_filename)
        np.savetxt(out_filename, output)

        return 

