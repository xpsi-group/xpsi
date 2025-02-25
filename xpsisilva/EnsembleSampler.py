
from xpsisilva.global_imports import *
from xpsisilva import _verbose

try:
    import emcee
    from emcee import EnsembleSampler as _EnsembleSampler
    from emcee.backends import HDFBackend
except ImportError:
    print('Check your emcee installation.')
    raise
else:
    if _verbose:
        print('Imported emcee version: %s' % emcee.__version__)

from xpsisilva.Parameter import StrictBoundsError
from xpsisilva.ParameterSubspace import ParameterSubspace
from xpsisilva.Posterior import PriorError

class EnsembleSampler(_EnsembleSampler):
    """ Derives from `emcee`_'s :class:`~emcee.EnsembleSampler`.

    .. _emcee: http://emcee.readthedocs.io/en/latest/

    :param callable posterior:
        The posterior function for MPI pickling.

    :param int ndims:
        Number of model parameters.

    :param bool resume:
        Resume a sampling run from file?

    :param int nwalkers:
        Number of ensemble walkers.

    :param nsteps:
        Number of iterations to perform.

    :param pool:
        Instance of :class:`schwimmbad.MPIPool`.

    :param str root_dir:
        Path to a directory containing a ``samples.h5`` to resume from;
        also the directory to be written to.

    :param list walker_dist_moments:
        Optional set of tuples of means and standard deviations.
        One per parameter.

    :param obj _posterior:
        Reserved for internal usage when invoking MPI.

    """
    def __init__(self,
                 ndims,
                 posterior,
                 pool,
                 resume,
                 nwalkers,
                 nsteps,
                 root_dir,
                 _posterior=None,
                 walker_dist_moments = None,
                 **kwargs):

        if _posterior is None:
            _posterior = posterior

        if not callable(posterior):
            raise RuntimeError('Posterior object is not callable.')

        self._posterior = posterior

        self._likelihood = _posterior.likelihood
        # get the prior object in case needed for initialisation
        self._prior = _posterior.prior

        try:
            self._root_dir = str(root_dir)
        except TypeError:
            raise TypeError('Root directory must be a string.')

        try:
            self._ndims = int(ndims)
            self._nwalkers = int(nwalkers)
            self._nsteps = int(nsteps)
        except TypeError:
            raise TypeError('Numbers of dimensions, walkers, and steps must '
                            'all be integers.')

        try:
            assert isinstance(resume, bool)
        except AssertionError:
            raise AssertionError('The ``resume`` argument is not boolean.')

        if resume:
            backend = HDFBackend(self._root_dir + "samples.h5")
            self._p0 = None
        else:
            if not _os.path.isdir(self._root_dir):
                _os.mkdir(self._root_dir)
            elif _os.path.isfile(self._root_dir + "samples.h5"):
                print('\nWarning: a ``samples.h5`` file exists '
                      'in ``%s``.' % self._root_dir)
                print('Attempting to move ``samples.h5`` to a subdirectory '
                      'of ``%s``.' % self._root_dir)
                try: # to archive the existing ``samples.h5`` file
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

                    _os.rename(self._root_dir + "samples.h5",
                               (self._root_dir +
                                ("old_samples/samples%s.h5" % temp)))
                except Exception:
                    print('Failed: file ``samples.h5`` will be overwritten.')
                else:
                    print('File archived in subdirectory ``%s``.' %
                          (self._root_dir + 'old_samples'))

            backend = HDFBackend(self._root_dir + "samples.h5")

            # Initialise walkers
            try:
                self._p0 = self._nd_ball(walker_dist_moments)
            except (AssertionError, PriorError):
                print('Attempting to initialise the walker positions via '
                      'inverse sampling...')
                self._p0 = self._prior.draw(nwalkers)[0]
                for i in range(nwalkers):
                    if not _np.isfinite(self._prior(self._p0[i,:])):
                        raise PriorError('Failed to initialise walkers.')
                print('Walker positions successfully initialised via inverse '
                      'sampling.')

        super(EnsembleSampler, self).__init__(nwalkers = nwalkers,
                                              ndim = ndims,
                                              log_prob_fn = posterior,
                                              pool = pool,
                                              backend = backend)

    def _nd_ball(self, moments):
        """ Draw points in parameter space from a n-dimensional Gaussian.

        :param list moments: Tuples of means and standard deviations.
                             One per parameter.

        :return:

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

        p0 = _np.empty((self._nwalkers, self._ndims), dtype = _np.double)

        q = self._likelihood.vector # make our own cache here

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

            p0[i,:] = p

        return p0

    def __call__(self):
        """ Hammer the nail. """
        print('\nCommencing posterior sampling.')

        self.run_mcmc(initial_state = self._p0,
                      nsteps = self._nsteps,
                      thin_by = 1,
                      store = True,
                      progress = True)

        print('\nSampling complete.')

    @property
    def posterior(self):
        """ Get the posterior object. """
        return self._posterior

    def get_backend(self):
        """ Get the :class:`emcee.backends.HDFBackend` instance. """
        try:
            return self.backend
        except AttributeError:
            raise AttributeError('emcee has not been run.')

