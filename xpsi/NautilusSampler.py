from xpsi import _verbose
from xpsi.Likelihood import Likelihood
from xpsi.Prior import Prior
from xpsi.utils import make_verbose

import numpy as np
import os

try:
    import nautilus

except ImportError:
    print('Check your nautilus-sampler installation. (pip install nautilus-sampler)')
    raise
else:
    if _verbose:
        print('Imported Nautilus.')


class NautilusSampler:
    """ Wrapper for the Nautilus (https://nautilus-sampler.readthedocs.io/) package (Lange 2023).

    :param likelihood: An instance of :class:`~.Likelihood.Likelihood`.

    :param prior: An instance of :class:`~.Prior.Prior`.

    :param sampler_build_kwargs: A dictionary of the keyword arguments passed to the instance of :class:`nautilus.Sampler`.
                                Some recommended parameters:
                                - "n_live" : 4000   Number of live points
                                - "n_networks" :    Number of networks used in the estimator. Default is 4.
                                - "n_batch" :       Number of likelihood evaluations that are performed at each step. 
                                            If likelihood evaluations are parallelized, should be multiple of the number of parallel processes. 
                                            If None, will be the smallest multiple of the pool size used for likelihood calls that is at least 100. 
                                            Recommended a few times the pool size. Default is None.
                                - "resume":True     Resume from previous run if filepath exists. Default is True
                                - "seed":0          Seed for random number generation used for reproducible results accross different runs. 
                                                    If None, results are not reproducible. Default is None.

    :param pool: Pool used for parallelization of likelihood calls and sampler calculations. 
                If None, no parallelization is performed. Default is None
                If an integer, the sampler will use a multiprocessing.Pool object with the specified number of processes.

    :param filepath: Path to the file where results are saved. If None, no results are written.
                    If not None, two files will be written: '.hdf5' file used to resume run, and '.txt' file with results for PostProcessing.                                                   

    """

    def __init__(
        self,
        likelihood : Likelihood,
        prior : Prior,
        sampler_build_kwargs,
        pool, 
        filepath=None,
    ):

        if not isinstance(likelihood, Likelihood):
            raise TypeError('Invalid type for likelihood object.')
        self._likelihood = likelihood

        if not isinstance(prior, Prior):
            raise TypeError('Invalid type for prior object.')
        self._prior = prior

        self._param_names = self._likelihood.names

        if pool is not None:
            print('Attention: NautilusSampler does not work yet on more than one node. Make sure that `pool` corresponds to the cpu-per-node.')

        self._sampler = nautilus.Sampler(
            self._prior.inverse_sample, # prior
            self.likelihood_wrapper,
            n_dim=len(self._param_names),
            pass_dict=False,
            pool=pool,
            filepath=filepath,
            **sampler_build_kwargs
        )

    @make_verbose('Commencing sampling', 'Sampling completed')
    def __call__(self, runtime_params):
        """ Run the sampler until target convergence criteria are fulfilled with the
        given runtime parameters.

        :param runtime_params: A dictionary of the keyword arguments passed to :func:`run`.
                        Some recommended parameters: 
                        - "discard_exploration" : True 
                             --Whether to discard points drawn in the exploration phase. 
                             This is required for a fully unbiased posterior and evidence estimate. 
                             Default is False.

                        - "verbose" : True
                             -- If True, print information about sampler progress. Default is False.

        """

        self._sampler.run(**runtime_params)

    def likelihood_wrapper(self, params):
        """Calculate the loglikelihood value for a given set of parameter values.
        Copy-pasted from UltraNestSampler interface.

        :param params: List of parameter values.

        :returns: Float with loglikelihood value for given set of parameter values.

        """

        arg1, *args = params

        # calculate the log-likelihood
        nautilus_likelihood = self._likelihood(p=[arg1, *args], reinitialise=True)

        if isinstance(nautilus_likelihood, np.ndarray):  # hackish way around: need to find better solution for this
            nautilus_likelihood = nautilus_likelihood[0]

        return nautilus_likelihood

    def write_results(self, out_filename):
        """ Get output txt file with columns containing weights, -2*loglikelihood,
        and parameters, which is the format required for post-processing within X-PSI.
        This is additional to output files UltraNest produces.

        :param out_filename: String specifying the name of the output file.

        """
        # extract results
        data, weights, logl = self._sampler.posterior()

        # save extra output file special for xpsi post-processing
        output = np.column_stack((np.exp(weights), -2 * logl, data))
        np.savetxt(f"{out_filename}.txt", output)
