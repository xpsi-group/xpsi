
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
        ``<root_dir>/<root>.txt`` (the post-processable cold chain) and
        ``<root_dir>/<root>_fullchain.npz``
        (the full temperature ladder, also used to resume).

    :param int checkpoint:
        Flush output to disk every ``checkpoint`` iterations. At each
        checkpoint the cold-chain posterior files hold all samples so far
        (after burn-in and thinning), and ``<root>_current_walkers.npz``/``.txt``
        are replaced with where every walker is at that iteration, at all
        temperatures (``ntemps * nwalkers`` rows in the ``.txt``).

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
        One per parameter. Used when ``init_mode='gaussian'``.

    :param str init_mode:
        How to initialise the walker positions (ignored when resuming):

            * ``'prior'``: random draws from the prior (inverse sampling);
            * ``'gaussian'``: a Gaussian ball defined by
              ``walker_dist_moments``;
            * ``'samples'``: draws from existing posterior samples given by
              ``init_samples`` (e.g. a MultiNest run);
            * ``None`` (default): ``'samples'`` if ``init_samples`` is given,
              else ``'gaussian'`` if ``walker_dist_moments`` is given, else
              ``'prior'``.

    :param init_samples:
        Posterior samples to initialise from, for ``init_mode='samples'``.
        Either a path to a text (``.txt``/``.dat``), ``.npy`` or walkers
        ``.npz`` file, or an array. Accepted layouts:

            * ``(nsamples, ndims + 2)``: X-PSI/MultiNest format, i.e. columns
              of sample weight, minus-two-log-likelihood, then the physical
              parameter vector (e.g. ``<root>.txt`` from a MultiNest run);
            * ``(nsamples, ndims)``: equally weighted parameter vectors (e.g.
              ``post_equal_weights.dat`` with its last column removed);
            * ``(ntemps, nwalkers, ndims)``: exact walker positions, e.g. a
              ``<root>_initial_walkers.npz`` or ``<root>_current_walkers.npz``
              file saved by an earlier run, used as is.

        Parameter columns must be in the order of the likelihood parameters.
        Samples are drawn in proportion to their weights, without replacement
        where enough samples with non-zero weight exist.

    :param float init_scatter:
        For ``init_mode='samples'``: scatter each drawn point by a Gaussian
        with standard deviation ``init_scatter`` times the (weighted) posterior
        standard deviation of each parameter, redrawing points that leave the
        prior support. Default is ``0.0`` (start exactly on the samples).

    :param bool save_init:
        Save the initial walker positions to ``<root_dir>/<root>_initial_walkers.npz``
        (array ``walkers`` of shape ``(ntemps, nwalkers, ndims)``, can be passed
        back as ``init_samples``) and a readable
        ``<root_dir>/<root>_initial_walkers.txt``?
        Default is ``False``.

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
                 init_mode: str = None,
                 init_samples = None,
                 init_scatter: float = 0.0,
                 save_init: bool = False,
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
            if init_mode is None:
                if init_samples is not None:
                    init_mode = 'samples'
                elif walker_dist_moments is not None:
                    init_mode = 'gaussian'
                else:
                    init_mode = 'prior'

            if init_mode == 'samples':
                if init_samples is None:
                    raise ValueError('``init_mode=\'samples\'`` requires '
                                     '``init_samples``.')
                self._p0 = self._from_samples(init_samples,
                                              float(init_scatter))
            elif init_mode == 'gaussian':
                if walker_dist_moments is None:
                    raise ValueError('``init_mode=\'gaussian\'`` requires '
                                     '``walker_dist_moments``.')
                try:
                    self._p0 = self._nd_ball(walker_dist_moments)
                except (AssertionError, PriorError):
                    self._p0 = self._from_prior()
            elif init_mode == 'prior':
                self._p0 = self._from_prior()
            else:
                raise ValueError('Unknown ``init_mode`` ``%s``: use '
                                 '``\'prior\'``, ``\'gaussian\'`` or '
                                 '``\'samples\'``.' % init_mode)

            if save_init:
                self._save_initial()

    def _from_prior(self):
        """ Draw the initial walker positions from the prior. """
        print('Attempting to initialise the walker positions via '
              'inverse sampling...')
        draws = self._prior.draw(self._ntemps * self._nwalkers)[0]
        p0 = draws.reshape(self._ntemps, self._nwalkers, self._ndims)
        for p in p0.reshape(-1, self._ndims):
            if not _np.isfinite(self._prior(p)):
                raise PriorError('Failed to initialise walkers.')
        print('Walker positions successfully initialised via inverse '
              'sampling.')
        return p0

    def _in_support(self, p):
        """ Check inclusion of ``p`` in the prior support.

        The parameter vector of the likelihood object is restored afterwards.

        """
        q = self._likelihood.vector # make our own cache here
        try:
            try:
                ParameterSubspace.__call__(self._likelihood, p)
            except StrictBoundsError:
                return False

            lp = self._prior(p)

            if _np.isnan(lp):
                raise PriorError('Failed to initialise walkers.')

            return bool(_np.isfinite(lp))
        finally:
            ParameterSubspace.__call__(self._likelihood, q)

    def _from_samples(self, samples, scatter):
        """ Draw the initial walker positions from existing posterior samples.

        :param samples: Path or array; see the class docstring.
        :param float scatter: Scatter in units of the posterior std.

        :return: An array of shape ``(ntemps, nwalkers, ndims)``.

        """
        source = 'array'
        if isinstance(samples, str):
            source = samples
            if not _os.path.isfile(samples):
                raise IOError('Initial samples file ``%s`` not found.'
                              % samples)
            if samples.endswith('.npz'):
                with _np.load(samples) as data:
                    if 'walkers' not in data.files:
                        raise ValueError('``%s`` holds no ``walkers`` array: '
                                         'pass a ``_initial_walkers.npz`` or '
                                         '``_current_walkers.npz`` file.'
                                         % samples)
                    samples = data['walkers']
            elif samples.endswith('.npy'):
                samples = _np.load(samples)
            else:
                samples = _np.loadtxt(samples)

        samples = _np.asarray(samples, dtype = _np.double)
        shape = (self._ntemps, self._nwalkers, self._ndims)

        # exact walker positions, e.g. saved by ``save_init``
        if samples.ndim == 3:
            if samples.shape != shape:
                raise ValueError('Initial walker array has shape %s but %s '
                                 'is required.' % (samples.shape, shape))
            for p in samples.reshape(-1, self._ndims):
                if not self._in_support(p):
                    raise PriorError('An initial walker position lies outside '
                                     'the prior support.')
            print('Walker positions initialised from ``%s``.' % source)
            return samples.copy()

        if samples.ndim == 1: # a single row
            samples = samples.reshape(1, -1)

        if samples.shape[1] == self._ndims + 2:
            weights = samples[:, 0].copy()
            points = samples[:, 2:]
        elif samples.shape[1] == self._ndims:
            weights = _np.ones(samples.shape[0])
            points = samples
        else:
            raise ValueError('Initial samples have %i columns but %i (weight, '
                             '-2logL, parameters) or %i (parameters) are '
                             'required.' % (samples.shape[1],
                                            self._ndims + 2, self._ndims))

        good = _np.isfinite(weights) & (weights > 0.0) \
                    & _np.all(_np.isfinite(points), axis = 1)
        weights = weights[good]
        points = points[good]

        # discard samples outside the prior support of the present model
        inside = _np.array([self._in_support(p) for p in points], dtype=bool)
        if not _np.all(inside):
            print('Warning: discarding %i initial samples outside the prior '
                  'support.' % _np.sum(~inside))
        weights = weights[inside]
        points = points[inside]

        if points.shape[0] == 0:
            raise PriorError('No usable initial samples in ``%s``.' % source)

        weights /= _np.sum(weights)
        n = self._ntemps * self._nwalkers

        mean = _np.sum(weights[:, None] * points, axis = 0)
        std = _np.sqrt(_np.sum(weights[:, None] * (points - mean)**2, axis=0))

        replace = points.shape[0] < n
        if replace:
            print('Warning: only %i usable initial samples for %i walkers; '
                  'drawing with replacement.' % (points.shape[0], n))
            if scatter <= 0.0:
                scatter = 1.0e-3
                print('Warning: setting ``init_scatter=%g`` to avoid '
                      'duplicate walkers.' % scatter)

        idx = _np.random.choice(points.shape[0], size = n,
                                replace = replace, p = weights)
        p0 = points[idx].copy()

        if scatter > 0.0:
            for i in range(n):
                for _ in range(1000):
                    p = points[idx[i]] + scatter * std \
                            * _np.random.normal(size = self._ndims)
                    if self._in_support(p):
                        p0[i] = p
                        break
                else: # keep the unscattered sample, which is in support
                    print('Warning: could not scatter initial walker %i '
                          'within the prior support.' % i)

        print('Walker positions initialised from %i posterior samples in '
              '``%s`` (scatter: %g posterior std).'
              % (points.shape[0], source, scatter))

        return p0.reshape(shape)

    def _walkers_path(self, which, ext):
        """ Path to saved walker positions (``initial`` or ``current``). """
        return self._root_dir + self._root + '_%s_walkers' % which + ext

    def _init_path(self, ext):
        """ Path to the saved initial walker positions. """
        return self._walkers_path('initial', ext)

    def _save_walkers(self, which, p, betas, label, iteration):
        """ Save walker positions of shape ``(ntemps, nwalkers, ndims)``.

        Writes ``<root>_<which>_walkers.npz`` (arrays ``walkers``, ``betas``
        and ``iteration``; usable as ``init_samples``) and a readable
        ``<root>_<which>_walkers.txt`` with one row per walker: temperature
        index (0 is the cold chain), walker index, then the parameters.

        """
        self._atomic_write(self._walkers_path(which, '.npz'),
                           lambda f: _np.savez(f, walkers = p,
                                               betas = _np.asarray(betas),
                                               iteration = iteration))

        ntemps, nwalkers, ndims = p.shape
        t, w = _np.meshgrid(_np.arange(ntemps), _np.arange(nwalkers),
                            indexing = 'ij')
        table = _np.column_stack((t.ravel(), w.ravel(), p.reshape(-1, ndims)))
        try:
            names = ' '.join(self._likelihood.names)
        except Exception:
            names = ' '.join('p%i' % i for i in range(ndims))
        header = ('%s, shape (ntemps=%i, nwalkers=%i, ndims=%i); betas: %s'
                  '\ntemp_index walker_index %s'
                  % (label, ntemps, nwalkers, ndims,
                     ' '.join('%.6g' % b for b in betas), names))
        fmt = ['%d', '%d'] + ['%.16e'] * ndims
        self._atomic_write(self._walkers_path(which, '.txt'),
            lambda f: _np.savetxt(f, table, fmt = fmt, header = header))

    def _save_initial(self):
        """ Save the initial walker positions. """
        self._save_walkers('initial', self._p0, self.betas,
                           'Initial walker positions', 0)
        print('Initial walker positions saved to ``%s`` and ``%s``.'
              % (self._init_path('.npz'), self._init_path('.txt')))

    def _resume_path(self):
        """ Path to the full-ladder resume/diagnostics file. """
        return self._root_dir + self._root + '_fullchain.npz'

    # names of the full-chain file written by earlier versions
    _legacy_fullchain = ('_ptemcee.npz', '_.npz')

    def _load_resume(self):
        """ Load the stored temperature ladder to resume from. """
        path = self._resume_path()
        if not _os.path.isfile(path):
            for ext in self._legacy_fullchain: # runs of earlier versions
                old = self._root_dir + self._root + ext
                if _os.path.isfile(old):
                    print('Resuming from ``%s`` (earlier file name).' % old)
                    path = old
                    break
            else:
                raise RuntimeError('Cannot resume: no ``%s`` file found.'
                                   % path)

        with _np.load(path) as _f:
            data = {key: _f[key] for key in _f.files}

        return data

    def _archive_existing(self):
        """ Move any existing output for this root to a subdirectory. """
        fullchain = ((self._resume_path(),) +
                     tuple(self._root_dir + self._root + ext
                           for ext in self._legacy_fullchain))
        if not any(_os.path.isfile(f) for f in fullchain):
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

            # .npy entries and the earlier full-chain names: files written by
            # earlier versions; archived so that, e.g., a stale <root>.npy is
            # not picked up by post-processing
            for ext in (('_fullchain.npz',) + self._legacy_fullchain +
                        ('.txt', '.npy',
                        '_initial_walkers.npz', '_initial_walkers.txt',
                        '_current_walkers.npz', '_current_walkers.txt',
                        '_initial_walkers.npy', '_current_walkers.npy')):
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
        """ Return the full-ladder arrays, concatenated with resumed samples.

        ptemcee preallocates its arrays with zeros for all ``nsteps`` at the
        start of a run, so they are cut to the iterations completed so far;
        otherwise a checkpoint written mid-run contains zero-filled samples.

        """
        n = self.time # iterations completed by this sampler instance
        chain = self.chain[:, :, :n, :]
        loglike = self.loglikelihood[:, :, :n]
        logprob = self.logprobability[:, :, :n]

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
        :func:`numpy.savez` does not append a spurious ``.npz`` to the temporary
        name and the replacement is seen by readers as a single atomic update.

        """
        tmp = path + '.tmp'
        with open(tmp, 'wb') as f:
            writer(f)
        _os.replace(tmp, path)

    def _flush(self):
        """ Write the current state to disk (checkpoint). """
        chain, loglike, logprob = self._assemble()

        # where every walker is now (all temperatures), replaced each time
        if chain.shape[2] > 0:
            self._save_walkers('current', chain[:, :, -1, :], self.betas,
                               'Walker positions after iteration %i'
                               % chain.shape[2], chain.shape[2])

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
