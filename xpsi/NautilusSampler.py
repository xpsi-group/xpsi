from xpsi import _verbose
from xpsi.Likelihood import Likelihood
from xpsi.Prior import Prior

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

    """

    def __init__(
        self,
        likelihood : Likelihood,
        prior : Prior,
        sampler_build_kwargs,
        MPI=True,
    ):

        if not isinstance(likelihood, Likelihood):
            raise TypeError('Invalid type for likelihood object.')
        self._likelihood = likelihood

        if not isinstance(prior, Prior):
            raise TypeError('Invalid type for prior object.')
        self._prior = prior

        self._param_names = self._likelihood.names

        if MPI:
            from mpi4py.futures import MPIPoolExecutor

            self._sampler = nautilus.Sampler(
                self._prior.inverse_sample, # prior
                self.likelihood_wrapper,
                n_dim=len(self._param_names),
                pass_dict=False,
                pool = MPIPoolExecutor(),
                **sampler_build_kwargs
            )

        else:

            self._sampler = nautilus.Sampler(
                self._prior.inverse_sample, # prior
                self.likelihood_wrapper,
                n_dim=len(self._param_names),
                pass_dict=False,
                **sampler_build_kwargs
            )

    def __call__(self, runtime_params):
        """ Run the sampler until target convergence criteria are fulfilled with the
        given runtime parameters.

        :param runtime_params: A dictionary of the keyword arguments passed to :func:`run`.

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
        output = np.column_stack((weights, -2 * logl, data))
        np.savetxt(f"{out_filename}.txt", output)
