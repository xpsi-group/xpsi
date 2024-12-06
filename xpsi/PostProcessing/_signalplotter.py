from .. import _warning

try:
    import h5py
except ImportError:
    _warning('Cannot import h5py for signal caching.')
    h5py = None

try:
    from tqdm.auto import trange as _range # detect if Jupyter notebook
except ImportError:
    def _range(n, *args, **kwargs): # ignore extra (kw)args
        return range(n)

from ._global_imports import *

from .. import Likelihood
from ._postprocessor import PostProcessor
from ._cache import _Cache
from ._signalplot import SignalPlot

class SignalPlotter(PostProcessor):
    """ Plot conditional posterior distributions of thematic X-ray signals.

    Methods to plot the data and model for posterior checking.

    Plots are generated for each posterior selected using the associated
    likelihood object.

    For a given model, there may be multiple :class:`~xpsi.Signal.Signal`
    instances per likelihood object. If this is the case, you need to
    reduce the model down to only the objects needed for functioning of
    the :class:`~xpsi.Signal.Signal` object to be handled. To do this, simply
    remove supply to your likelihood object the references to this minimal
    set of objects. This minimises computation time and ensures explicit
    declaration of the signals to be plotted.

    .. note::

        If a model has multiple instruments, then the energy set for signal
        integration is calculated based on waveband coverage union.
        If instruments are omitted from the likelihood object in order
        to execute posterior signal plotting, the number of energies
        that span the waveband of the remaining instrument waveband should
        be set to match the number in the full model (with all instruments)
        if the likelihood factor for the remaining instrument is to be
        exactly the same.

    """
    @fix_random_seed
    @make_verbose('Plotting signals for posterior checking',
                  'Plotted signals for posterior checking')
    def plot(self,
             plots,
             IDs=None,
             combine=False,
             combine_all=False,
             only_combined=False,
             force_combine=True,
             nsamples=200,
             cache=True,
             force_cache=False,
             cache_dir='./',
             read_only=False,
             archive=True):
        """ Compute and plot signals *a posteriori*.

        :param dict plots:
            Dictionary of lists of plot objects, where each dictionary key
            must match a posterior ID.

        :param OrderedDict IDs:
            Keys must be string identifiers of :class:`Runs` instances.
            Each dictionary element must be a list of string identifiers,
            each matching objects collected in :class:`Runs` instance
            corresponding to the key. Defaults to ``None``, meaning attempt to
            use as many runs as possible subject to plotting restrictions.

        :param bool cache:
            Cache intermediary model objects to accelerate post-processing?

        :param bool force_cache:
            Force caching irrespective of checks of existing cache files.
            Useful if some part of the model is tweaked and the cache file with
            the same name and sample set is not manually moved from the
            designated directory..

        :param int nsamples:
            Number of samples to use. Equally-weighted samples are generated,
            thus introducing a additional Monte Carlo noise which is ignored.

        :param int num_phases:
            Number of phases to interpolate at on the interval [0,2] cycles.

        :param str filename:
            Filename of cache.

        :param str cache_dir:
            Directory to write cache to.

        :param bool read_only:
            Do not write to cache file?

        :param bool archive:
            If not read-only, then archive an existing cache file found at the
            same path?

        """
        self.set_subset(IDs, combine, combine_all,
                        force_combine, only_combined)

        for posterior in self.subset:
            yield 'Handling posterior %s' % posterior.ID
            run = posterior.subset_to_plot[0]

            try:
                likelihood = posterior.likelihood
            except AttributeError as e:
                print('Supply a likelihood object to proceed.')
                raise e

            state_to_restore = likelihood.externally_updated
            likelihood.externally_updated = False

            try:
                _plots = plots[posterior.ID]
            except (TypeError, KeyError) as e:
                print('Invalid plot object specification.')
                raise e
            else:
                try:
                    iter(_plots)
                except TypeError:
                    if not isinstance(_plots, SignalPlot):
                        raise TypeError('Invalid plot object type.')
                    _plots.posterior = run.prepend_ID
                    _plots = [_plots]
                else:
                    for plot in _plots:
                        if not isinstance(plot, SignalPlot):
                            raise TypeError('Invalid plot object type.')
                        plot.posterior = run.prepend_ID

            # assimilate all declared targets
            _caching_targets = []
            for plot in _plots:
                try:
                    _caching_targets += plot.__caching_targets__
                except (TypeError, AttributeError) as e:
                    print('Invalid specification of caching targets.')
                    raise e

            # eliminate duplicates (some plot types can share caching targets)
            caching_targets = []
            for target in _caching_targets:
                if target not in caching_targets:
                    caching_targets.append(target)

            # caching targets are assimilated from plot objects
            likelihood.signal.caching_targets = caching_targets

            self._driver(run,
                         likelihood,
                         nsamples,
                         cache,
                         force_cache,
                         cache_dir,
                         read_only,
                         archive,
                         _plots)

            likelihood.externally_updated = state_to_restore

            yield 'Handled posterior %s.' % posterior.ID

        self._plots = plots

        # in case user needs a handle (e.g., if plot objects created
        # via a classmethod), could return handle here
        yield

    @property
    def plots(self):
        """ Get the dictionary of plot objects last processed. """
        return self._plots

    @staticmethod
    def _draw_equally_weighted(samples, nsamples, num_params):
        """ Get a set of equally weighted samples from a weighted set.

        ..note::

            Additional Monte Carlo noise from trimming sample set ignored.

        """

        assert nsamples < samples.shape[0], 'Number of samples for plotting \
                                             cannot exceed number of nested \
                                             samples.'

        weights = samples[:,0].copy()
        if _np.abs(1.0 - _np.sum(weights)) > 0.01:
            print('Warning: 1 - (sum of weights) = %.8e'%(1.0-_np.sum(weights)))
            print('Weights renormalized to sum to unity.')

        weights /= _np.sum(weights)

        indices = _np.random.choice(samples.shape[0], nsamples,
                                    replace=True, p=weights)

        thetas = samples[indices, 2:2+num_params]

        return thetas

    def _driver(self,
                run,
                likelihood,
                nsamples,
                cache,
                force_cache,
                cache_dir,
                read_only,
                archive,
                plots):
        """ Execute plotting loop given samples. """

        # If has a parameter to plot, use it. Else, draw samples !
        if hasattr(plots[0],'parameters_vector'):
            try:
                print( plots[0].parameters_vector )
                _ = plots[0].parameters_vector[0]
                assert len(likelihood) == len(plots[0].parameters_vector) , 'Length of the parameter vector does not match likelihood one'
            except (TypeError) as e:
                raise e('When provided, theta to plot must be given as an iterable and have matching length with likelihood' )
            else:
                thetas = _np.array( [plots[0].parameters_vector] )
        else:
            thetas = self._draw_equally_weighted(run.samples, nsamples,
                                                len(likelihood))


        signal = likelihood.signal # should only be one signal object available
        for p in plots:
            p.signal = likelihood.signal

        names = likelihood.names

        if cache and h5py is not None:
            try:
                s = signal.prefix
            except AttributeError:
                s = ''
            filename = run.prepend_ID.replace(' ', '_')
            filename += '__signal' + ('_' + s + '__' if s else '__')
            filename += 'cached__'
            filename += '__'.join(signal.caching_target_names)

            cache = _Cache(filename,
                           cache_dir,
                           read_only,
                           archive)
            if cache.do_caching(thetas, force_cache):
                # skips body if can simply read cache
                for i in _range(thetas.shape[0], desc='Signal caching loop'):
                    likelihood([thetas[i,run.get_index(n)] for n in names])
                    cache.cache(signal.caching_targets)
        elif cache and h5py is None:
            raise ImportError('You need to install h5py to use caching.')

        def update(theta):
            # order the parameter vector appropriately
            vector = [theta[run.get_index(n)] for n in names]

            if cache: # use the cache if available
                # set the parameter values in case needed by plot objects
                # that call methods of the signal object, or methods of objects
                # that the signal object encapsulates refernences to, an
                # example being the interstellar attenuation applied for
                # the spectrum plot type
                super(Likelihood, likelihood).__call__(vector)

                # now restore the signals objects that were cached
                cached = next(cache)
                for key, value in cached.items():
                    try:
                        delattr(signal, key)
                    except AttributeError:
                        pass

                    if len(value.shape) == 3:
                        for i in range(value.shape[0]):
                            setattr(signal, key, value[i,...])
                    else:
                        setattr(signal, key, value)
            else: # otherwise resort to likelihood evaluations
                likelihood(vector)

        def wrapper(plot_obj, delete_me=None, index=None):
            """ Wrap a plot obj's next method into a cache-enabled callback. """
            if cache:
                cache.reset_iterator()

            if delete_me is not None:
                try:
                    iter(delete_me)
                except TypeError:
                    delete_me = [delete_me]

                for attr in delete_me:
                    try:
                        delattr(plot_obj, attr)
                    except AttributeError:
                        pass

            # ignore x because it is already known by a plot object in
            # this more object oriented approach to calling fgivenx:
            if index is not None:
                def callback(x, theta):
                    update(theta)
                    return next(plot_obj)[index]
                return callback # for fgivenx
            else:
                def callback(x, theta):
                    update(theta)
                    return next(plot_obj)
                return callback

        for plot in plots:
            with plot: # acquire plot as context manager
                plot.execute(thetas, wrapper) # iterate over samples

                truths = [run.truth_vector[run.get_index(n)] for n in likelihood.names]
                if None not in truths:
                    likelihood(truths)
