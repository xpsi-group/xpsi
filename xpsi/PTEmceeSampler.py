
from xpsi.global_imports import *
from xpsi import _verbose

try:
    import numpy as _numpy
    # ptemcee 1.0.0 still uses a couple of numpy names that are deprecated,
    # so restore them (when absent) before importing ptemcee. 
    if not hasattr(_numpy, 'float'):
        _numpy.float = float
    if not hasattr(_numpy, 'trapz') and hasattr(_numpy, 'trapezoid'):
        _numpy.trapz = _numpy.trapezoid
    import ptemcee
    from ptemcee import Sampler as _PTSampler
except ImportError:
    print('Check your ptemcee installation.')
    raise
else:
    if _verbose:
        print('Imported ptemcee version: %s' % ptemcee.__version__)

try:
    from tqdm.auto import tqdm as _tqdm # detects Jupyter notebook vs terminal
except ImportError:
    _tqdm = None

from xpsi.Parameter import StrictBoundsError
from xpsi.ParameterSubspace import ParameterSubspace
from xpsi.Posterior import PriorError

class PTEmceeSampler(_PTSampler):
    """ Derives from `ptemcee`_'s :class:`~ptemcee.Sampler`.

    Parallel-tempered affine-invariant ensemble sampling. Unlike
    :class:`~.EmceeSampler.EmceeSampler`, which samples the combined posterior,
    :mod:`ptemcee` requires the log-likelihood and log-prior as separate
    callables and operates an ensemble per inverse temperature. The cold chain
    (:math:`\\beta=1`) yields the posterior samples; the full temperature
    ladder additionally provides a thermodynamic-integration estimate of the
    evidence.

    .. _ptemcee: https://github.com/willvousden/ptemcee

    :param int ndims:
        Number of model parameters.

    :param callable loglikelihood:
        The log-likelihood function (for MPI pickling).

    :param callable logprior:
        The log-prior function (for MPI pickling).

    :param pool:
        Instance of :class:`schwimmbad.MPIPool`.

    :param bool resume:
        Resume a sampling run from file?

    :param int nwalkers:
        Number of ensemble walkers per temperature.

    :param int nsteps:
        Number of iterations to perform.

    :param str root_dir:
        Directory to be written to (and read from when resuming).

    :param int ntemps:
        Number of temperatures in the ladder. If ``None`` (and ``betas`` is
        not supplied), :mod:`ptemcee` chooses the number from ``ndims`` and
        ``Tmax``.

    :param float Tmax:
        Maximum temperature of the ladder. Passed to :mod:`ptemcee`.

    :param betas:
        Optional array of inverse temperatures, overriding ``ntemps``/``Tmax``.

    :param bool adapt:
        Adapt the temperature ladder during sampling? Default is ``True``.

    :param str root:
        File-root basename for output files, written as
        ``<root_dir>/<root>.txt``, ``<root_dir>/<root>.npy`` (the
        post-processable cold chain) and ``<root_dir>/<root>_ptemcee.npz``
        (the full temperature ladder, also used to resume).

    :param int checkpoint:
        Flush output to disk every ``checkpoint`` iterations. Default is 100.

    :param bool progress:
        Display a :mod:`tqdm` progress bar over the iterations? Default is
        ``True`` (silently ignored if :mod:`tqdm` is not installed).

    :param int burn_in:
        Number of initial iterations to discard when writing the
        post-processable cold chain. Default is 0.

    :param int thin:
        Thinning factor applied to the cold chain when writing the
        post-processable output. Default is 1.

    :param list walker_dist_moments:
        Optional set of tuples of means and standard deviations.
        One per parameter.

    :param obj _posterior:
        Reserved for internal usage when invoking MPI.

    """
    def __init__(self,
                 ndims: int,
                 loglikelihood,
                 logprior,
                 pool,
                 resume: bool,
                 nwalkers: int,
                 nsteps: int,
                 root_dir: str,
                 ntemps: int = None,
                 Tmax: float = None,
                 betas = None,
                 adapt: bool = True,
                 root: str = 'run',
                 checkpoint: int = 500,
                 burn_in: int = 1000,
                 thin: int = 5,
                 progress: bool = True,
                 _posterior = None,
                 walker_dist_moments = None,
                 **kwargs):

        if not callable(loglikelihood) or not callable(logprior):
            raise RuntimeError('Log-likelihood and log-prior must be callable.')

        if _posterior is None:
            raise RuntimeError('A posterior object is required.')

        self._posterior = _posterior
        self._likelihood = _posterior.likelihood
        self._prior = _posterior.prior # get in case needed for initialisation

        try:
            self._root_dir = str(root_dir)
            self._root = str(root)
        except TypeError:
            raise TypeError('Root directory and root must be strings.')

        try:
            self._ndims = int(ndims)
            self._nwalkers = int(nwalkers)
            self._nsteps = int(nsteps)
            self._checkpoint = int(checkpoint)
            self._burn_in = int(burn_in)
            self._thin = int(thin)
        except TypeError:
            raise TypeError('Numbers of dimensions, walkers, steps, checkpoint '
                            'interval, burn-in and thinning must all be '
                            'integers.')

        try:
            assert isinstance(resume, bool)
        except AssertionError:
            raise AssertionError('The ``resume`` argument is not boolean.')

        self._adapt = bool(adapt)
        self._progress = bool(progress) and _tqdm is not None
        self._resume_data = None

        if resume:
            self._resume_data = self._load_resume()
            betas = self._resume_data['betas']
        else:
            if not _os.path.isdir(self._root_dir):
                _os.mkdir(self._root_dir)
            else:
                self._archive_existing()

        super(PTEmceeSampler, self).__init__(self._nwalkers,
                                             self._ndims,
                                             loglikelihood,
                                             logprior,
                                             ntemps = ntemps,
                                             Tmax = Tmax,
                                             betas = betas,
                                             pool = pool)

        # the number of temperatures may have been chosen internally
        self._ntemps = len(self.betas)

        if resume:
            # continue from the last stored positions of every temperature
            self._p0 = self._resume_data['chain'][:, :, -1, :]
        else:
            try:
                self._p0 = self._nd_ball(walker_dist_moments)
            except (AssertionError, PriorError):
                print('Attempting to initialise the walker positions via '
                      'inverse sampling...')
                draws = self._prior.draw(self._ntemps * self._nwalkers)[0]
                self._p0 = draws.reshape(self._ntemps,
                                         self._nwalkers,
                                         self._ndims)
                for p in self._p0.reshape(-1, self._ndims):
                    if not _np.isfinite(self._prior(p)):
                        raise PriorError('Failed to initialise walkers.')
                print('Walker positions successfully initialised via inverse '
                      'sampling.')

    def _resume_path(self):
        """ Path to the full-ladder resume/diagnostics file. """
        return self._root_dir + self._root + '_ptemcee.npz'

    def _load_resume(self):
        """ Load the stored temperature ladder to resume from. """
        path = self._resume_path()
        if not _os.path.isfile(path):
            raise RuntimeError('Cannot resume: no ``%s`` file found.' % path)

        with _np.load(path) as _f:
            data = {key: _f[key] for key in _f.files}

        return data

    def _archive_existing(self):
        """ Move any existing output for this root to a subdirectory. """
        if not _os.path.isfile(self._resume_path()):
            return

        print('\nWarning: output for root ``%s`` exists '
              'in ``%s``.' % (self._root, self._root_dir))
        print('Attempting to move existing output to a subdirectory '
              'of ``%s``.' % self._root_dir)
        try: # to archive the existing output files
            if not _os.path.isdir(self._root_dir + "old_samples/"):
                _os.mkdir(self._root_dir + "old_samples/")
            from datetime import datetime
            obj = datetime.now()
            temp = '__datetime__%i.%i.%i__%i.%i.%i' % (obj.day,
                                                       obj.month,
                                                       obj.year,
                                                       obj.hour,
                                                       obj.minute,
                                                       obj.second)

            for ext in ('_ptemcee.npz', '.txt', '.npy'):
                src = self._root_dir + self._root + ext
                if _os.path.isfile(src):
                    _os.rename(src, (self._root_dir +
                                     ("old_samples/%s%s%s" %
                                      (self._root, temp, ext))))
        except Exception:
            print('Failed: existing output will be overwritten.')
        else:
            print('Output archived in subdirectory ``%s``.' %
                  (self._root_dir + 'old_samples'))

    def _nd_ball(self, moments):
        """ Draw points in parameter space from a n-dimensional Gaussian.

        :param list moments: Tuples of means and standard deviations.
                             One per parameter.

        :return: An array of shape ``(ntemps, nwalkers, ndims)``.

        """
        if moments is not None:
            try:
                if isinstance(moments, zip):
                    moments = list(moments)
                assert isinstance(moments, list)
                assert len(moments) == self._ndims
                for m in moments:
                    assert len(m) == 2
            except AssertionError:
                print('\nInvalid specification of initial walker position '
                      'distribution.')
                raise
        else:
            raise PriorError

        moments = list(map(list, zip(*moments)))

        def helper(p):
            """ Check inclusion in prior support. """

            try:
                ParameterSubspace.__call__(self._likelihood, p)
            except StrictBoundsError:
                return False

            lp = self._prior(p) # p object not used explicitly, but for clarity

            if _np.isnan(lp):
                raise PriorError('Failed to initialise walkers.')

            if not _np.isfinite(lp):
                return False

            return True

        p0 = _np.empty((self._ntemps, self._nwalkers, self._ndims),
                       dtype = _np.double)

        q = self._likelihood.vector # make our own cache here

        for t in range(self._ntemps):
            for i in range(self._nwalkers):

                included = False
                counter = 0

                try:
                    while not included and counter < 1000:
                        p = moments[0] + moments[1] * _np.random.normal(size=self._ndims)
                        included = helper(p)
                        counter += 1

                    if not included:
                        raise PriorError('Failed to initialise walkers.')

                    ParameterSubspace.__call__(self._likelihood, q)
                except PriorError:
                    ParameterSubspace.__call__(self._likelihood, q)
                    raise

                p0[t,i,:] = p

        return p0

    def __call__(self):
        """ Hammer the nail. """
        print('\nCommencing parallel-tempered posterior sampling.')

        # drive the run through ptemcee's per-iteration generator so we can
        # both update the progress bar and flush output at each checkpoint
        iterator = self.sample(self._p0,
                               iterations = self._nsteps,
                               adapt = self._adapt)

        if self._progress:
            iterator = _tqdm(iterator,
                             total = self._nsteps,
                             desc = 'Sampling',
                             dynamic_ncols = True)

        for i, _ in enumerate(iterator):
            if (i + 1) % self._checkpoint == 0 or (i + 1) == self._nsteps:
                self._flush()

        try:
            logZ, dlogZ = self.log_evidence_estimate()
            print('\nThermodynamic-integration log-evidence estimate: '
                  '%.6g +/- %.6g' % (logZ, dlogZ))
        except Exception:
            print('\nCould not compute a log-evidence estimate.')

        print('\nSampling complete.')

    def _assemble(self):
        """ Return the full-ladder arrays, concatenated with resumed samples. """
        chain = self.chain
        loglike = self.loglikelihood
        logprob = self.logprobability

        if self._resume_data is not None:
            chain = _np.concatenate((self._resume_data['chain'], chain), axis=2)
            loglike = _np.concatenate((self._resume_data['loglikelihood'],
                                       loglike), axis=2)
            logprob = _np.concatenate((self._resume_data['logprobability'],
                                       logprob), axis=2)

        return chain, loglike, logprob

    def _build_posterior_array(self, chain, loglike):
        """ Build the post-processable cold-chain array.

        Columns: sample weight, minus-two-log-likelihood, then the parameter
        vector -- the layout read by :mod:`xpsi.PostProcessing` (and by
        MultiNest-style getdist output).

        """
        cold = chain[0, :, self._burn_in::self._thin, :]
        cold_logl = loglike[0, :, self._burn_in::self._thin]

        nwalkers, nsteps, ndims = cold.shape
        flat = cold.reshape(nwalkers * nsteps, ndims)
        flat_logl = cold_logl.reshape(nwalkers * nsteps)

        out = _np.zeros((flat.shape[0], ndims + 2), dtype = _np.double)
        out[:, 0] = 1.0                 # equal-weight MCMC samples
        out[:, 1] = -2.0 * flat_logl    # minus two log-likelihood
        out[:, 2:] = flat

        return out

    @staticmethod
    def _atomic_write(path, writer):
        """ Write via a temporary file and atomically replace the target.

        The ``writer`` receives an open binary file handle (not a path) so
        :func:`numpy.save` does not append a spurious ``.npy`` to the temporary
        name and the replacement is seen by readers as a single atomic update.

        """
        tmp = path + '.tmp'
        with open(tmp, 'wb') as f:
            writer(f)
        _os.replace(tmp, path)

    def _flush(self):
        """ Write the current state to disk (checkpoint). """
        chain, loglike, logprob = self._assemble()

        # full temperature ladder, for resuming and diagnostics
        self._atomic_write(self._resume_path(),
            lambda f: _np.savez(f,
                                chain = chain,
                                loglikelihood = loglike,
                                logprobability = logprob,
                                betas = self.betas,
                                acceptance_fraction = self.acceptance_fraction,
                                tswap_acceptance_fraction =
                                    self.tswap_acceptance_fraction))

        # post-processable cold chain
        out = self._build_posterior_array(chain, loglike)
        self._atomic_write(self._root_dir + self._root + '.txt',
                           lambda f: _np.savetxt(f, out))
        self._atomic_write(self._root_dir + self._root + '.npy',
                           lambda f: _np.save(f, out))

    @property
    def posterior(self):
        """ Get the posterior object. """
        return self._posterior

    def get_chain(self):
        """ Get the sampled chain, shape ``(ntemps, nwalkers, nsteps, ndims)``. """
        try:
            return self.chain
        except AttributeError:
            raise AttributeError('ptemcee has not been run.')
