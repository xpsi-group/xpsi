from xpsi import _verbose
from xpsi.Likelihood import Likelihood
from xpsi.Prior import Prior
from xpsi import _size
from xpsi.global_imports import *

import numpy as np
import os

try:
    import dynesty
except ImportError:
    print('Check your Dynesty installation.')   # crtl+C to stop run before finished
    raise
else:
    if _verbose:
        print("Imported Dynesty!")

class DynestySampler(dynesty.DynamicNestedSampler):   # inherits dynamic nested sampler from Dynesty
    """ Wrapper for ...

    :param likelihood: An instance of :class:`~.Likelihood.Likelihood`.

    :param prior: An instance of :class:`~.Prior.Prior`.

    ~:param n_params: amount of prior parameters

    ~:param other_params: 

    """
    def __init__(self,
                 likelihood,
                 prior,
                 other_params,
                 MPI):   # do I also need possible extra params? if those wish to be set, library -> dictionary with live points or something?
        
        # raise errors if necessary
        if not isinstance(likelihood, Likelihood):  # both likelihood and prior need to be functions, does this still work then?
            raise TypeError('Invalid type for likelihood object.')
        self._likelihood = likelihood

        if not isinstance(prior, Prior):
            raise TypeError('Invalid type for prior object.')
        self._prior = prior

        
        # Could make some default options here for other_params if wanted

        # initialise sampler with super
        if MPI or _size > 1: # or _size
            try:
                from schwimmbad import MPIPool
            except ImportError:
                raise ImportError('Check your schwimmbad package installation.')

            # Initialize the MPI-based pool used for parallelisation.
            with MPIPool() as pool:
                if not pool.is_master():
                    # Wait for instructions from the master process.
                    pool.wait()
                    _sys.exit(0)

                #sampling
                super().__init__(loglikelihood=self.likelihood_function,
                                prior_transform=self._prior.inverse_sample,   # function translating a unit cube to the parameter space according to the prior. The input is a 1-d numpy array with length ndim, where each value is in the range [0, 1)
                                 ndim=len(self._likelihood.prior),   # similar as in main.py
                                pool=pool,
                                **other_params)  # unpacking dictionary of parameters so that they can be used as all seperate parameters of the function
 
        else:
            # Initialise emcee sampler
            # original code without parallelization 
            super().__init__(loglikelihood=self.likelihood_function,
                            prior_transform=self._prior.inverse_sample,   # function translating a unit cube to the parameter space according to the prior. The input is a 1-d numpy array with length ndim, where each value is in the range [0, 1)
                            ndim=len(self._likelihood.prior),   # similar as in main.py
                            **other_params)  # unpacking dictionary of parameters so that they can be used as all seperate parameters of the function
    
    def __call__(self, runtime_params):
        """ Run the sampler until target convergence criteria are fulfilled with the 
        given runtime parameters. 
        
        :param runtime_params: A dictionary of the keyword arguments passed to :func:`run_nested`.
        
        """

        self.run_nested(**runtime_params)

    def likelihood_function(self, params):
        """Calculate the loglikelihood value for a given set of parameter values.  
        
        :param params: 1D numpy array of length ...
 
        :returns: Float with loglikelihood value for given set of parameter values. 

        """

        arg1, *args = params   # arg 1 of params + rest

        # calculate the log-likelihood
        log_likelihood = self._likelihood(p=[arg1, *args], reinitialise=True)   #params now in list p so they can be used by Likelihood class of xpsi

        if isinstance(log_likelihood, np.ndarray): # hackish way around: need to find better solution for this
            log_likelihood = log_likelihood.item()   # handles possibility of unexpected return of numpy array of size 1 (turns it into float)
        
        return log_likelihood
    
    def write_results(self, log_dir, out_filename):
        """ Get output txt file with columns containing weights, -2*loglikelihood, 
        and parameters, which is the format required for post-processing within X-PSI. 
        This is additional to output files UltraNest produces.  

        :param other_params:

        :param out_filename: String specifying the name of the output file. 

        """
        # extract results
        data = np.array(self.results.samples)
        weights = np.array(self.results.importance_weights())
        logl = np.array(self.results.logl)

        #file_path = os.path.join(other_params['log_dir'], out_filename)
        output = np.column_stack((weights, -2*logl, data))
        os.makedirs(log_dir, exist_ok=True)
        file_path = os.path.join(log_dir, out_filename)
        np.savetxt(file_path, output)

        # output has colums of length of amount of samples that have been done
        # output has rows with a value for weights, a value for -2*logl and a value for every parameter that is being sampled for
        # so rows are length of # of parameters + 2