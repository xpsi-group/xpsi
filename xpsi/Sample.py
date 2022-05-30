__all__ = ["nested", "ensemble"]

from xpsi.global_imports import *
#from xpsi import global_imports
from xpsi import _comm, _rank, _size, make_verbose
from xpsi import Likelihood

posterior = None

# Function defined on the module-level for function and module name pickling
# during point-to-point MPI comminications. The global variable is set in
# emcee() function before creation of an MPIPool.
def func(p):
    global posterior
    return posterior(p)

def ensemble(likelihood, prior, MPI = True, **kwargs):
    """ Initialise `emcee <http://dfm.io/emcee/current/>`_ and sample.

    :param likelihood:
        An instance of :class:`~.Likelihood.Likelihood`.

    :param prior:
        An instance of :class:`~.Prior.Prior`.

    :param bool MPI:
        Parallelise with MPI? If calling script not lauched with an MPI
        directive, sampling will not commence because there is only one
        process. Default is ``True`` since only in testing is it justifiable to
        use a single process.

    :param kwargs:
        Passed to initialisers of appropriate classes:

            * boolean to resume, under keyword :obj:`resume`;
            * number of steps, under keyword :obj:`nsteps`;
            * number of walkers, under keyword :obj:`nwalkers`;
            * moments of initial walker multivariate Gaussian distribution,
              under keyword :obj:`walker_dist_moments` (can be ``None``);
            * root directory for output, under keyword :obj:`root_dir`.

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

    if MPI or _size > 1:
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
                                      pool = pool,
                                      _posterior = posterior,
                                      **kwargs)

            # Commence emcee sampling process
            sampler()
    else:
        from xpsi.EnsembleSampler import EnsembleSampler

        # Initialise emcee sampler
        sampler = EnsembleSampler(ndims = len(likelihood),
                                  posterior = posterior,
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

@make_verbose('Importance sampling commencing',
              'Importance sampling completed')
def importance(target, importance,
               sample_root, names,
               weight_threshold=1.0e-3,
               likelihood_change=False,
               prior_change=False,
               overwrite=True):
    """ Importance sample from target using approximative distribution.

    :param obj target:
        A :class:`~.Likelihood` object for the target posterior distribution.

    :param obj importance:
        A :class:`~.Likelihood` object for the importance posterior
        distribution. This distribution should be some sufficient approxmation
        to the target distribution for validity.

    :param str sample_root:
        Path to sample files, with the file root but without the file extension.

    :param list names:
        The parameter names for the samples on disk. These names must match
        the parameter names in the :class:`~.Likelihood` objects.

    :param float weight_threshold:
        Only use samples with weight greater than this fraction of the maximum
        weight amongst samples.

    :param bool likelihood_change:
        Declare whether the likelihood function has changed. If so, it is
        recommended that workload is distributed amongst a world of MPI
        processes on a cluster.

    :param bool prior_change:
         Declare whether the prior density function has changed. Generally,
         this reweighting factor is inexpensive to evaluate and can be
         executed on a desktop in minutes.

    :param bool overwrite:
        Overwrite an existing importance sample file on disk?

    """

    yield 'Cross-checking parameter names'
    if len(target) > len(importance):
        raise TypeError('Every target distribution parameter must '
                        'be an importance distribution parameter.')

    for name in target.names:
        if name not in names:
            raise NameError('Every target distribution parameter must be '
                            'an importance sample space dimension.')

    assert len(importance) == len(names), ('The importance samples must be '
                                           'draws from the the importance '
                                           'distribution')

    for name in importance.names:
        if name not in names:
            raise NameError('Every importance distribution parameter must be '
                            'an importance sample space dimension.')

    _target_state_to_restore = target.externally_updated
    target.externally_updated = False

    _importance_state_to_restore = importance.externally_updated
    importance.externally_updated = False

    yield 'Loading samples'
    if _rank == 0:
        samples = _np.loadtxt(sample_root + '.txt', dtype=_np.double)

        _max = _np.max(samples[:,0])
        _ref = importance_sampled = samples[samples[:,0]/_max >= weight_threshold,:]

        yield 'Importance sample fraction = %.3f%%.'%(100.0 * _ref.shape[0]/samples.shape[0])

    def _prior_reweight(theta):
        """ Helper function for to break conditional statement. """
        if not hasattr(target, 'prior') or not hasattr(importance, 'prior'):
            return 1.0 # silently terminate
        if not hasattr(target.prior, 'density') or not hasattr(importance.prior, 'density'):
            return 1.0 # silently terminate

        _target_theta = [theta[names.index(n)] for n in target.names]
        _importance_theta = [theta[names.index(n)] for n in importance.names]
        if not likelihood_change:
            super(Likelihood, target).__call__(_target_theta)
            super(Likelihood, importance).__call__(_importance_theta)

        _target = target.prior.density(_target_theta)
        _importance = importance.prior.density(_importance_theta)

        if not likelihood_change: # memoization: revert parameter value cache
            super(Likelihood, target).__call__(target.cached)
            super(Likelihood, target).__call__(target.vector)

            super(Likelihood, importance).__call__(importance.cached)
            super(Likelihood, importance).__call__(importance.vector)

        return _target/_importance

    yield 'Distributing workload amongst %i MPI processes'%_size
    if _rank == 0:
        # number of MPI scatterings
        _iterations = int(_m.ceil(_ref.shape[0]/_size))
    else:
        _iterations = None

    _iterations = _comm.bcast(_iterations, root=0)

    for i in range(_iterations):
        if _rank == 0:
            if i == _iterations - 1 and _ref.shape[0]%_size != 0:
                _remaining = _ref.shape[0]%_size
                theta = [list(vector) for vector in _ref[i*_size:i*_size+_remaining,2:]]
                theta += [None]*(_size - _remaining)
            else:
                theta = [list(vector) for vector in _ref[i*_size:(i+1)*_size,2:]]
        else:
            theta = None

        theta = _comm.scatter(theta, root=0)
        weight = None
        _target = None

        if theta is not None:
            if likelihood_change:
                _target = target([theta[names.index(n)] for n in target.names])
                _importance = importance([theta[names.index(n)] for n in importance.names])
                # objects return logarithm of likelihood
                weight = _np.exp(_target - _importance)

            if prior_change:
                if weight is not None:
                    weight *= _prior_reweight(theta)
                else:
                    weight = _prior_reweight(theta)

        weight = _comm.gather(weight, root=0)
        if likelihood_change:
            _target = _comm.gather(_target, root=0)

        if _rank == 0:
            if i == _iterations - 1 and _ref.shape[0]%_size != 0:
                _ref[i*_size:i*_size+_remaining,0] *= _np.array(weight[:_remaining])
                if likelihood_change:
                    _ref[i*_size:i*_size+_remaining,1] = -2.0*_np.array(_target[:_remaining])
            else:
                _ref[i*_size:(i+1)*_size,0] *= _np.array(weight)
                if likelihood_change:
                    _ref[i*_size:(i+1)*_size,1] = -2.0*_np.array(_target)

    yield 'Renormalising weights'
    if _rank == 0:
        # renormalise
        _ref[:,0] /= _np.sum(_ref[:,0])

        _file = sample_root + '__importance_sampled' + '.txt'
        if _os.path.isfile(_file):
            _write = True if overwrite else False
        else:
            _write = True

        if _write:
            yield 'Writing to disk'
            _np.savetxt(_file, _ref)

    yield 'Restoring likelihood-object parameter update protocol attributes'
    target.externally_updated = _target_state_to_restore
    importance.externally_updated = _importance_state_to_restore

    yield
