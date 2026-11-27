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


# Module-level likelihood and log-likelihood function. A module-level function
# is pickled by reference (module + name only), so with an MPI pool only the
# parameter vectors are sent to the workers, never the X-PSI likelihood object.
# Every MPI rank must call ``set_likelihood`` before the pool is created
# (``run_nautilus`` does this).
_likelihood = None

def set_likelihood(likelihood):
    """ Register the likelihood used by :func:`loglike` on this process. """
    global _likelihood
    _likelihood = likelihood

def loglike(params):
    """ Log-likelihood of the registered X-PSI likelihood for one point. """
    value = _likelihood(p=list(params), reinitialise=True)
    if isinstance(value, np.ndarray):
        value = value[0]
    return value


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
                If an integer, the sampler will use a multiprocessing.Pool object with the specified number of processes
                (one node only).
                For several nodes pass an MPI pool (e.g. ``schwimmbad.MPIPool``);
                ``run_nautilus(..., MPI=True)`` does this.

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

        if isinstance(pool, int):
            print('Attention: an integer `pool` uses multiprocessing and only works on '
                  'one node. Use run_nautilus(..., MPI=True) under mpiexec/srun for '
                  'several nodes.')

        set_likelihood(self._likelihood)

        self._sampler = nautilus.Sampler(
            self._prior.inverse_sample, # prior
            loglike,                    # module-level: cheap to pickle (MPI and multiprocessing)
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

    @make_verbose('Commencing sampling (with temporary output)', 'Sampling completed')
    def run_with_checkpoints(self, runtime_params, out_filename,
                             checkpoint_every=None, checkpoint_interval=None):
        """ Run the sampler in chunks and write a temporary posterior file
        (``<out_filename>_temp.txt``, same format as :meth:`write_results`)
        after each chunk, similar to the intermediate output of MultiNest or
        the ``checkpoint`` option of the PTMC sampler.

        The ``.hdf5`` file is still written by Nautilus itself after every
        batch and is what is used to resume. Splitting the run into chunks
        does not change the sampling, it only adds file writes.

        :param runtime_params: Dictionary of keyword arguments for :func:`run`.
            ``n_like_max`` and ``timeout`` are respected as overall limits.

        :param out_filename: Root of the output files.

        :param int checkpoint_every: Write the temporary file every
            ``checkpoint_every`` likelihood evaluations (rounded up to whole
            batches by Nautilus).

        :param float checkpoint_interval: Write the temporary file every
            ``checkpoint_interval`` seconds. Can be combined with
            ``checkpoint_every`` (whichever is reached first).

        :returns: ``True`` if the run converged, ``False`` if it stopped
            because of ``n_like_max`` or ``timeout``.
        """
        if checkpoint_every is None and checkpoint_interval is None:
            return self._sampler.run(**runtime_params)

        params = dict(runtime_params)
        n_like_max = params.pop('n_like_max', np.inf)
        timeout = params.pop('timeout', np.inf)

        import time
        t_start = time.time()
        temp_root = f"{out_filename}_temp"

        while True:
            target = n_like_max
            if checkpoint_every is not None:
                target = min(n_like_max, self._sampler.n_like + int(checkpoint_every))
            remaining = timeout - (time.time() - t_start)
            chunk_timeout = remaining
            if checkpoint_interval is not None:
                chunk_timeout = min(remaining, float(checkpoint_interval))

            success = self._sampler.run(n_like_max=target,
                                        timeout=chunk_timeout,
                                        **params)

            self.write_temp_results(temp_root)

            if success:
                return True
            if self._sampler.n_like >= n_like_max or time.time() - t_start >= timeout:
                return False

    def write_temp_results(self, temp_root):
        """ Write the current (unconverged) posterior to ``<temp_root>.txt``,
        overwriting the previous temporary file. Failures are reported but do
        not stop the run. """
        try:
            if self._sampler.n_like == 0:
                return
            data, weights, logl = self._sampler.posterior()
            output = np.column_stack((np.exp(weights), -2 * logl, data))
            tmp = f"{temp_root}.txt.part"
            with open(tmp, 'w') as f:  # write then rename: never a half-written file
                np.savetxt(f, output)
            os.replace(tmp, f"{temp_root}.txt")
            if _verbose:
                print(f'Temporary results written to {temp_root}.txt '
                      f'({self._sampler.n_like} likelihood calls, '
                      f'log Z = {self._sampler.log_z:.3f}).')
        except Exception as e:
            print(f'Could not write temporary Nautilus results: {e}')

    def close_pool(self):
        """ Close the multiprocessing pool(s) Nautilus created when ``pool``
        was given as an integer. Without this the pool is only cleaned up at
        interpreter shutdown, which prints a harmless but noisy
        ``OSError: [Errno 9] Bad file descriptor`` from ``Pool.__del__``. """
        from multiprocessing.pool import Pool
        seen = set()
        for npool in (getattr(self._sampler, 'pool_l', None),
                      getattr(self._sampler, 'pool_s', None)):
            inner = getattr(npool, 'pool', None)
            if isinstance(inner, Pool) and id(inner) not in seen:
                seen.add(id(inner))
                inner.close()
                inner.join()

    def likelihood_wrapper(self, params):
        """ Log-likelihood for a given set of parameter values.

        :param params: List of parameter values.

        :returns: Float with loglikelihood value for given set of parameter values.

        """
        set_likelihood(self._likelihood)
        return loglike(params)

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
