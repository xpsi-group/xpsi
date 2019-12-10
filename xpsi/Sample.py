from __future__ import division, print_function

__all__ = ["nested", "ensemble"]

from .global_imports import *
from . import global_imports

posterior = None

# Function defined on the module-level for function and module name pickling
# during point-to-point MPI comminications. The global variable is set in
# emcee() function before creation of an MPIPool.
def func(p):
    global posterior
    return posterior(p)

def ensemble(likelihood, prior, MPI = True, **kwargs):
    """ Initialise `emcee <http://dfm.io/emcee/current/>`_ and sample.

    :param likelihood: An instance of :class:`~.Likelihood.Likelihood`.

    :param prior: An instance of :class:`~.Prior.Prior`.

    :param bool MPI: Parallelise with MPI? If calling script not lauched with
                     an MPI directive, sampling will not commence because there
                     is only one process. Default is ``True`` since only in
                     testing is it justifiable to use a single process.


    :param kwargs: Passed to initialisers of appropriate classes.

       * boolean to resume, under keyword :obj:`resume`;
       * number of steps, under keyword :obj:`nsteps`;
       * number of walkers, under keyword :obj:`nwalkers`;
       * moments of initial walker multivariate Gaussian distribution, under
         keyword :obj:`walker_dist_moments` (can be ``None``);
       * root directory for output, under keyword :obj:`root_dir`;

    The above objects are used to instantiate :class:`~.Posterior.Posterior`.

    :return: An instance of :class:`emcee.backends.HDFBackend`.

    """
    # use the globally scoped variable
    global posterior
    from xpsi.Posterior import Posterior
    # Callable instance of the posterior
    posterior = Posterior(likelihood,
                          prior,
                          **kwargs)

    if MPI or global_imports._size > 1:
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

            from xpsi.EnsembleSampler import EnsembleSampler

            # Initialise emcee sampler
            sampler = EnsembleSampler(ndims = len(likelihood),
                                      posterior = func,
                                      prior = prior,
                                      pool = pool,
                                      **kwargs)

            # Commence emcee sampling process
            sampler()
    else:
        from xpsi.EnsembleSampler import EnsembleSampler

        # Initialise emcee sampler
        sampler = EnsembleSampler(ndims = len(likelihood),
                                  posterior = posterior,
                                  prior = prior,
                                  pool = None,
                                  **kwargs)

        # Commence emcee sampling process
        sampler()

    return sampler.backend

def nested(likelihood, prior, check_kwargs={}, **kwargs):
    """ Initialise MultiNest and integrate.

    :param dict likelihood:
        Keyword arguments required for instantiation of
        :class:`~.Likelihood.Likelihood`.

    :param dict prior:
        Keyword arguments required for instantiation of
        :class:`~.Prior.Prior`.

    :param dict check_kwargs:
        Keywords for likelihood function checker. Passed to checker before
        commencement of sampling.

    :param kwargs: Keyword arguments for PyMultiNest.

    """
    from xpsi.NestedSampler import NestedSampler

    if check_kwargs:
        likelihood.check(**check_kwargs)

    sampler = NestedSampler(len(likelihood),
                            likelihood,
                            prior)
    sampler(**kwargs)
