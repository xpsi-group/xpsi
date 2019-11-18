from __future__ import division, print_function

__all__ = ["Run",
           "NSBackend",
           "Runs",
           "PostProcessor",
           "publication_rc_settings",
           "slide_rc_settings"]

from .global_imports import *
from . import global_imports
from . import make_verbose, verbose, fragile, _verbose
from . import __version__

import six
from functools import wraps

import copy as _copy

from scipy.special import logsumexp

from matplotlib import pyplot as plt
from matplotlib import rc, rcParams
from matplotlib.ticker import MultipleLocator, MaxNLocator, AutoMinorLocator,\
                                AutoLocator, ScalarFormatter, LogLocator,\
                                NullFormatter

def _get_default_locator(prune):
    return MaxNLocator(nbins=4, min_n_ticks=3, prune=prune)

def _get_default_formatter():
    default = ScalarFormatter(useOffset=False)
    default.set_powerlimits(lims = (-2.0,3.0))
    return default

from matplotlib import gridspec
from matplotlib import cm
import matplotlib.patches as mpatches

from .Likelihood import Likelihood
from .Prior import Prior
from . import Pulse
from .tools.phase_interpolator import interpolate_pulse as interp
from .tools.phase_integrator import phase_integrator
from .tools.energy_interpolator import energy_interpolator

try:
    import h5py
except ImportError:
    if _verbose:
        print('Warning: Cannot import h5py for caching.')

try:
    import getdist
except ImportError:
    if _verbose:
        print('Warning: Cannot import GetDist.')
else:
    from getdist import loadMCSamples, plots
    # the following disables getdist.chains.print_load_details
    getdist.chains.print_load_details = False
    if _verbose:
        print('Imported GetDist version: %s' % getdist.__version__)

try:
    import nestcheck
except ImportError:
    if _verbose:
        print('Warning: Cannot import nestcheck.')
else:
    if _verbose:
        print('Imported nestcheck version: %s' % nestcheck.__version__)
    from nestcheck.data_processing import (process_multinest_run,
                                           process_polychord_run)
    from nestcheck.ns_run_utils import combine_ns_runs, get_logw, get_w_rel
    from nestcheck.write_polychord_output import write_run_output
    from nestcheck.plots import (bs_param_dists,
                                 param_logx_diagram,
                                 weighted_1d_gaussian_kde)
    from nestcheck.error_analysis import run_ci_bootstrap
    from nestcheck.estimators import param_cred, logz

    try:
        from nestcheck.plots import getdist_kde
    except ImportError:
        if _verbose:
            print('Warning: Using native nestcheck KDE instead of GetDist KDE.')

try:
    import fgivenx
except ImportError:
    if _verbose:
        print('Warning: Cannot import fgivenx.')

try:
    from emcee.backends import HDFBackend
except ImportError:
    if _verbose:
        print('Warning: Cannot import emcee.')

def publication_rc_settings():
    """ Serif font with TeX, for papers. """
    rc('font',**{'family':'serif',
                 'serif':['Computer Modern Roman']})
    rc('text', usetex=True)

def slide_rc_settings():
    """ Sans-serif font without Tex, for slides. """
    rc('font',**{'family':'sans-serif',
                 'sans-serif':['Computer Modern Sans serif']})
    rc('text', usetex=False)

random_seed = None
def fix_random_seed(func):
    @wraps(func)
    def wrapper(*args, **kwargs):
        global random_seed
        state = _np.random.get_state()
        _np.random.seed(random_seed)
        _ = func(*args, **kwargs)
        _np.random.set_state(state)
        return _
    return wrapper

class AmbiguityError(Exception):
    """ Thrown if ambiguous sample-set IDs are given when plotting. """

class ParameterError(Exception):
    """ Thrown if inconsistent parameter names are specified for plotting. """

class Run(object):
    """ Base class for a sample object, containing basic information.

    :param str ID: For identification of the set of the samples if multiple
                   sets are jointly plotted.

    :param str filepath: Path to samples file.

    :param list names: An ordered list of ``str`` parameter names. The
                       ordering must match the parameter vector ordering
                       defined in sample backend objects.

    :param dict bounds: A dictionary of one-dimensional *hard* parameter
                        bounds. See :class:`getdist.mcsamples.MCSamples`;
                        the keys must match the :obj:`names` list. For
                        the purpose of density estimation plots these
                        bounds can be viewing bounds.

    :param list labels: An ordered list of (LaTeX compatible) ``str``
                        literal parameter labels.

    :param str implementation: Sampling software applied to generate the
                               samples. Known options are
                               ``['multinest', 'polychord', 'emcee']``.

    :param dict truths: Optional dictionary of parameter truths, if *known*;
                        if *unknown* leave as ``None``. The keys
                        must be names which match those in
                        :obj:`names`.

    """

    def __init__(self, filepath, ID, names, bounds, labels=None,
                 implementation=None, truths=None):

        self.ID = ID
        self.samples = filepath
        self.names = names
        self.bounds = bounds

        if labels is None:
            self._labels = dict(zip(self.names, self.names))
        else:
            self.labels = labels

        if implementation is not None:
            self.implementation = implementation
        else:
            self._implementation = None

        if truths is not None and None not in truths.itervalues():
            self.truths = truths
        else:
            self._truths = dict(zip(self.names, [None] * len(self.names)))

    @property
    def ID(self):
        """ Get the identification ``str`` of the sample set. """
        return self._ID

    @ID.setter
    def ID(self, obj):
        """ Set the identification ``str`` of the sample set. """

        if not isinstance(obj, str):
            raise TypeError('Invalid sample set ID specification. '
                            'See the docstring.')
        self._ID = obj

    @property
    def names(self):
        """ Get the parameter names. """

        return self._names

    @names.setter
    def names(self, obj):
        """ Set the parameter names. """

        try:
            for name in obj:
                if not isinstance(name, str):
                    raise TypeError
        except TypeError:
            print('Invalid parameter name specification. See the docstring.')
            raise

        self._names =  obj

    @property
    def bounds(self):
        """ Get the parameter bounds. """

        return self._bounds

    @bounds.setter
    def bounds(self, obj):
        """ Set the parameter bounds. """

        try:
            if len(obj) != len(self.names):
                raise TypeError
            for key in obj:
                if not isinstance(obj[key], (list, tuple)):
                    raise TypeError
                for x in obj[key]:
                    if x is not None: float(x)
        except TypeError:
            print('Invalid parameter bound specification. See the docstring.')
            raise

        self._bounds =  obj

    @property
    def labels(self):
        """ Get the parameter labels. """

        return self._labels

    @labels.setter
    def labels(self, obj):
        """ Set the parameter labels. """

        try:
            if len(obj) != len(self.names):
                raise TypeError
            for key in obj:
                if key not in self.names:
                    raise TypeError
                if not isinstance(obj[key], str):
                    raise TypeError
        except TypeError:
            print('Invalid parameter name specification. See the docstring.')
            raise

        self._labels =  obj

    @property
    def implementation(self):
        """ Sampling software applied to generate the samples of the run. """

        return self._implementation

    @implementation.setter
    def implementation(self, obj):
        """ Set the name of the sampling software. """

        if obj in ['multinest', 'polychord', 'emcee']:
            self._implementation = obj
        else:
            print('Warning: Unrecognised software was used to generate the '
                  'samples of run with ID %s. The functionality of this '
                  'module may be incompatible with the sample information '
                  % self.ID)

    @property
    def truths(self):
        """ Get the parameter truths as a dictionary. """
        return self._truths

    @truths.setter
    def truths(self, obj):
        """ Set the parameter truths. """

        try:
            if len(obj) != len(self.names):
                raise TypeError
            for key in obj:
                float(obj[key])
                if key not in self.names:
                    raise TypeError
        except TypeError:
            print('Invalid parameter truth specification. See the docstring.')
            raise

        self._truths = obj

    @property
    def truth_vector(self):
        """ Get the parameter truths as a list. """
        v = [None] * len(self.names)
        for i, name in enumerate(self.names):
            v[i] = self.truths[name]
        return v

    @property
    def samples(self):
        """ Get a copy of the samples array. """

        return self._samples.copy()

    @samples.setter
    def samples(self, filepath):
        if _os.path.isfile(filepath):
            self._samples = _np.loadtxt(filepath)
        else:
            raise ValueError('File %s does not exist.' % filepath)

    @property
    def lines(self):
        """ Get the dictionary of line arguments for :mod:`getdist`.

        :param dict lines: :mod:`getdist`-compatible dictionary of parameters
                           specifying the properties of the smooth lines
                           representing one-dimensional marginal distributions
                           of parameters.
        """
        return self._lines

    @lines.setter
    def lines(self, obj):
        if not isinstance(obj, dict):
            raise TypeError('Invalid line argument specification. '
                            'See the docstring.')
        self._lines = obj

    @property
    def contours(self):
        """ Get the dictionary of contour arguments for :mod:`getdist`.

        :param dict contours: :mod:`getdist`-compatible dictionary of parameters
                              specifying the properties of the contours
                              representing credible regions of the
                              two-dimensional marginal distributions of
                              parameters.
        """
        return self._contours

    @contours.setter
    def contours(self, obj):

        if not isinstance(obj, dict):
            raise TypeError('Invalid contour argument specification. '
                            'See the docstring.')
        self._contours = obj

class NSBackend(Run):
    """
    Container for nested samples generated by a single run, and backends
    for analysis of the run.

    :param dict kde_settings:
        Settings for instantiation of :class:`getdist.mcsamples.MCSamples`.

    The other keyword arguments are generic properties passed to the parent
    class, such as the identification (ID) string of the run.

    """

    def __init__(self, root, base_dir, kde_settings, use_nestcheck,
                 transform=None, **kwargs):
        filerootpath =_os.path.join(base_dir, root)
        _filerootpath = filerootpath

        if transform is not None:
            samples = _np.loadtxt(filerootpath+'.txt')
            ndims = samples.shape[1] - 2
            temp = transform(samples[0,2:])
            ntransform = len(temp) - ndims

            if not _os.path.isfile(filerootpath+'_transformed.txt'):
                transformed = _np.zeros((samples.shape[0],
                                         samples.shape[1] + ntransform))
                transformed[:,:2] = samples[:,:2]
                for i in range(samples.shape[0]):
                    transformed[i,2:] = transform(samples[i,2:])
                _np.savetxt(filerootpath+'_transformed.txt', transformed)

            filerootpath += '_transformed'
            root += '_transformed'

        super(NSBackend, self).__init__(filepath=filerootpath+'.txt',**kwargs)

        # getdist backend
        self._gd_bcknd = getdist.mcsamples.MCSamples(root=filerootpath,
                             settings=kde_settings,
                             sampler='nested',
                             names=self.names,
                             ranges=self.bounds,
                             labels=[self.labels[name] for name in self.names])
        self._gd_bcknd.readChains(getdist.chains.chainFiles(filerootpath))

        self._kde_settings = kde_settings

        self.use_nestcheck = use_nestcheck

        if self.use_nestcheck: # nestcheck backend
            if transform is not None:
                for ext in ['dead-birth.txt', 'phys_live-birth.txt']:
                    if not _os.path.isfile(filerootpath + ext):
                        samples = _np.loadtxt(_filerootpath + ext)
                        transformed = _np.zeros((samples.shape[0],
                                                 samples.shape[1] + ntransform))
                        transformed[:,ndims+ntransform:] = samples[:,ndims:]
                        for i in range(samples.shape[0]):
                            transformed[i,:ndims+ntransform] =\
                                                transform(samples[i,:ndims])

                        _np.savetxt(filerootpath + ext, transformed)

                if not _os.path.isfile(filerootpath + '.stats'):
                    if _os.path.isfile(_filerootpath + '.stats'):
                        try:
                            from shutil import copyfile as _copyfile
                        except ImportError:
                            pass
                        else:
                            _copyfile(_filerootpath + '.stats',
                                      filerootpath + '.stats')
            try:
                kwargs['implementation']
            except KeyError:
                print('Root %r sampling implementation not specified... '
                      'assuming MultiNest for nestcheck...')
                self._nc_bcknd = process_multinest_run(root,
                                                       base_dir=base_dir)
            else:
                if kwargs['implementation'] == 'multinest':
                    self._nc_bcknd = process_multinest_run(root,
                                                           base_dir=base_dir)
                elif kwargs['implementation'] == 'polychord':
                    self._nc_bcknd = process_polychord_run(root,
                                                           base_dir=base_dir)
                else:
                    raise ValueError('Cannot process with nestcheck.')

    @property
    def getdist_backend(self):
        """ Get the :class:`getdist.mcsamples.MCSamples` instance. """
        return self._gd_bcknd

    @property
    def kde_settings(self):
        """ Get the input :mod:`getdist` KDE settings dictionary. """
        return self._kde_settings

    @property
    def nestcheck_backend(self):
        """ Get the :mod:`nestcheck` backend for the nested samples. """
        return self._nc_bcknd

    @property
    def margeStats(self):
        """ Return the marginal statistics using :mod:`getdist`. """
        return self._mcsamples.getMargeStats()

class Params(object):
    """ Information about parameters shared by runs for plotting. """

    def __init__(self, names):
        self._names = names

    @property
    def names(self):
        return self._names

    @property
    def idxs(self):
        return self._idxs

    @property
    def bounds(self):
        return self._bounds

    @property
    def labels(self):
        return self._labels

    @property
    def truths(self):
        return self._truths

class Runs(object):
    """ Container for nested sampling run objects. """

    def __init__(self, runs):

        try:
            iter(runs)
        except TypeError:
            if isinstance(runs, Run):
                self._runs = [runs]
            else:
                raise TypeError('Run objects must be instances of the ``Run`` '
                                'class.')
        else:
            IDs = []
            for run in runs:
                if not isinstance(run, Run):
                    raise TypeError('Run objects must be instances of the '
                                    '``Run`` class.')
                if run.ID not in IDs:
                    IDs.append(run.ID)
                else:
                    raise AmbiguityError('Use distinct IDs for distinct sets '
                                         'of samples.')
            self._runs = runs

    @classmethod
    def load_all(cls, roots, IDs, use_nestcheck, *args, **kwargs):
        """ Construct a :class:`Run` instance by loading distinct runs with
            equivalent arguments. """

        runs = []
        for root, ID, use_nestcheck in zip(roots, IDs, use_nestcheck):
            runs.append(NSBackend(root, *args,
                                  ID=ID,
                                  use_nestcheck=use_nestcheck,
                                  **kwargs))

        return cls(runs)

    def set_subset(self, IDs=None, combine=False, combine_all=False,
                   force_combine=False):
        """ Set a current list of :class:`Run` instances."""
        if IDs is None:
            self._subset = self._runs
        else:
            self._subset = []
            for ID in IDs:
                for run in self._runs:
                    if run.ID == ID:
                        self._subset.append(run)

        if force_combine: # create new run object
            self._combine(combine_all)
        elif combine:
            if not hasattr(self, '_combined'):
                self._combine(combine_all)
        else:
            self._combined = None

    def _combine(self, combine_all):
        if combine_all:
            str_IDs = '_'.join(str(ID) for ID in self.get_attr('ID',
                                                    current=False,
                                                    nestcheck_compatible=True))
        else:
            str_IDs = '_'.join(str(ID) for ID in self.get_attr('ID',
                                                    current=True,
                                                    nestcheck_compatible=True))

        for run in self._subset:
            if run.use_nestcheck:
                break

        base_dir = run.nestcheck_backend['output']['base_dir']
        file_root = 'combined_IDs_%s_' % str_IDs
        if not _os.path.isfile(_os.path.join(base_dir, file_root+'.txt')):
            if combine_all:
                run = combine_ns_runs(self.get_attr('nestcheck_backend',
                                                    current=False,
                                                    nestcheck_compatible=True))
            else:
                run = combine_ns_runs(self.get_attr('nestcheck_backend',
                                                    current=True,
                                                    nestcheck_compatible=True))

            run['output']['base_dir'] = base_dir
            run['output']['file_root'] = file_root

            write_run_output(run,
                             write_dead = True,
                             write_stats = True,
                             posteriors = True,
                             stats_means_errs = True,
                             n_simulate = 1000)

        kwargs = {'kde_settings': _.kde_settings,
                  'ID': 'combined',
                  'implementation': 'polychord',
                  'names': _.names,
                  'bounds': _.bounds,
                  'labels': _.labels,
                  'truths': _.truths}

        self._combined = NSBackend(file_root,
                                   base_dir = base_dir,
                                   use_nestcheck = True,
                                   **kwargs)
    @property
    def subset(self):
        """ Get the current subset of runs and parameters for plotting. """

        return self._subset

    def __getitem__(self, ID):
        """ Get a :class:`Run` instance using the associated ID. """
        for run in self._runs:
            if ID == run.ID:
                return run

        raise KeyError('No run with ID matching key.')

    @property
    def combined(self):
        """ Try to get a combined run in the form of a :mod:`nestcheck`
            backend. """
        try:
            return self._combined
        except AttributeError:
            print('No combined run available. Set the run-subset and combine.')
            raise
        finally:
            if self._combined is None:
                raise ValueError('You must combine runs.')

    def get_attr(self, attribute, current=True, nestcheck_compatible=False):
        """ Get a list of attributes of the :class:`Run` instances stored as
            the current subset. """

        if nestcheck_compatible:
            return [getattr(run, attribute) for run in \
                 (self._subset if current else self._runs) if run.use_nestcheck]
        else:
            return [getattr(run, attribute) for run in \
                                     (self._subset if current else self._runs)]

    @property
    def params(self):
        """ Get the current parameter information. """
        return self._params

    def set_params(self, names):
        """ Set current parameters for plotting which are shared. """

        self._params = Params(names)
        self._params._idxs = [None] * len(names)
        for i, param in enumerate(names):
            for run in self._subset:
                for k, name in enumerate(run.names):
                    if param == name:
                        if self._params._idxs[i] is None:
                            self._params._idxs[i] = k
                        elif self._params._idxs[i] != k:
                            self._params = None
                            raise ValueError('Parameter %s index not uniform '
                                             'across runs.' %
                                             (param,
                                              run.ID,
                                              self._subset[0].ID))
                        break
                    elif k == len(run.names) - 1:
                        self._params = None
                        raise ValueError('No parameter name matching %s in run '
                                         'with ID %s.' % (param, run.ID))

        try:
            self._set_params('bounds')
            self._set_params('labels')
            self._set_params('truths')
        except Exception:
            self._params = None
            raise

    def _set_params(self, attrs):
        """ Check consistency of a parameter across runs. """

        _attrs = '_' + attrs

        setattr(self._params, _attrs, [None] * len(self._params._names))
        try:
            for i, param in enumerate(self._params._names):
                for run in self._subset[1:]:
                    if (getattr(self._subset[0], attrs)[param] \
                            != getattr(run, attrs)[param]):
                        raise ValueError('Inconsistent %s for parameter'
                                         ' %s between runs %s and %s.' %
                                         (attrs, param,
                                          self._subset[0].ID, run.ID))
                getattr(self._params, _attrs)[i] = \
                        getattr(self._subset[0], attrs)[param]
        except (AttributeError, KeyError):
            print('Parameter %s not specified correctly.' % attrs)
            raise

class PostProcessor(object):
    """ Post-process samples for inference and posterior checking.

    Calculate inferences, usually in the form of approximate integrals over
    the posterior distribution. Also visualise posterior information and
    derived quantities.

    Perform posterior predictive checking.

    :param runs: An instance of :class:`Runs` or an instance of :class:`Run`.

    :param likelihood: The instance of :class:`~.Likelihood` used for
                       sampling, or a clone of that instance. Defaults to
                       ``None``. If multiple likelihood functions are
                       associated with the runs due to discreteness in the
                       model space, one can pass a dictionary with keys
                       matching the IDs of the runs.

    """
    def __init__(self, runs = None, likelihood = None):

        if runs is not None:
            self.runs = runs

        if likelihood is not None:
            self._likelihood = likelihood

    @classmethod
    def load_runs(cls,
                  *args,
                  **kwargs):
        """ Construct :class:`Runs` instance and then :class:`PostProcessor`
        instance.

        :param `*args`: Passed directly to :meth:`Runs.load_all`.

        :param `**kwargs`:
            Passed directly to :meth:`Runs.load_all`. Also a ``'likelihood'``
            argument can be handled for ``__init__``.

        """
        likelihood = kwargs.pop('likelihood')
        return cls(Runs.load_all(*args, **kwargs), likelihood)

    @property
    def runs(self):
        """ Get the nested sampling runs. """
        return self._runs

    @runs.setter
    def runs(self, obj):
        """ Set the nested sampling runs object. """
        if isinstance(obj, Runs):
            self._runs = obj
        else:
            if isinstance(obj, Run):
                self._runs = Runs(obj)
            else:
                raise TypeError('An instance of :class:`Runs` must be '
                                'supplied.')
    @property
    def likelihood(self):
        """ Get the likelihood instance. """

        try:
            return self._likelihood
        except AttributeError:
            return None

    @likelihood.setter
    def likelihood(self, obj): # make dictionary possible
        """ Set the likelihood object. """
        if isinstance(obj, Likelihood):
            self._likelihood = obj
        else:
             raise TypeError('The likelihood object needs to derive from '
                             'xpsi.Likelihood.Likelihood.')

    @make_verbose('Curating set of runs for posterior plotting',
                  'Run set curated')
    def _set_runs_to_plot(self, IDs, combine, combine_all):
        """ Helper function to get and notify which runs will be plotted. """

        if IDs is None:
            IDs = self._runs.get_attr('ID', current=False)

        if len(IDs) >= 5 and not combine:
            self._runs.set_subset(IDs[:5], combine, combine_all)
            print('Warning: Only the first five positional runs will be '
                  'plotted individually, with IDs'
                  + ', '.join(str(ID) for ID in self._runs.get_attr('ID')))

        elif len(IDs) == 5 and combine:
            self._runs.set_subset(IDs[:4], combine, combine_all)
            print('Warning: Only the first four positional runs will be '
                  'plotted individually, with IDs '
                  + ', '.join(str(ID) for ID in self._runs.get_attr('ID')))
        else:
            self._runs.set_subset(IDs, combine, combine_all)

    def _set_line_and_contour_args(self, lw=1.0, alpha=1.0, **kwargs):
        """ Match the :mod:`nestcheck` color scheme.

        Always assigns reds to a combined run if it is found to exist.

        """

        nestcheck_colors = ['darkred', 'darkblue', 'darkgrey', 'darkgreen',
                            'darkorange']

        try:
            self._runs.combined.lines = {'lw': kwargs.get('lw_1d', lw),
                                          'color': 'darkred',
                                          'alpha': alpha}
            self._runs.combined.contours = {'lw': lw,
                                            'color': 'darkred',
                                            'alpha': alpha}
        except ValueError:
            i=0
        else:
            i=1

        for run, color in zip(self._runs.subset,
                              nestcheck_colors[i:i+len(self._runs.subset)]):
            run.lines =  {'lw': kwargs.get('lw_1d', lw),
                          'color': color,
                          'alpha': alpha}
            run.contours = {'lw': lw, 'color': color, 'alpha': alpha}

    @fix_random_seed
    @make_verbose('Executing posterior density estimation',
                  'Posterior density estimation complete')
    def plot_posteriorDensity(self, params,
                               run_IDs=None,
                               combine=False,
                               combine_all=False,
                               only_combined=False,
                               bootstrap_estimators=True,
                               bootstrap_density=False,
                               separate_plots=False,
                               write=False,
                               root_filename='',
                               directory='./',
                               ext='.pdf',
                               dpi=300,
                               maxdots=2000,
                               **kwargs):
        """
        Generate posterior density plots with :mod:`getdist` and :mod:`nestcheck` for
        nested sampling runs.

        Up to five runs can be plotted natively via :mod:`nestcheck`; beyond
        such a number the plots are generally display too much information and
        clarity is lost.

        :param list params:
            List of parameter strings for plotting. Must match
            identifier strings which are attributes of :class:`Run`
            instances, and those parameters must occupy the same
            indices across the :class:`Run` instances (a limitation,
            but the :mod:`nestcheck` package assumes a set of realisations
            of the same nested sampling process).

        :param list run_IDs:
            An iterable containing string identifiers of instances
            of :class:`Run` instances. These identifiers must
            match objects collected in the :class:`Runs` instance
            stored as an attribute of an instance of the
            :class:`PostProcessing` class. Defaults to ``None``,
            meaning attempt to use all runs in the attribute
            instance of :class:`Runs`.

        :param bool combine:
            Additionally combine the runs into a single run for overplotting?

        :param bool combine_all:
            Combine all runs in the :class:`Runs` instance or only those
            for which IDs are provided? Ignored if ``combine`` is ``False``.

        :param bool only_combine:
            Only plot the combined run? Ignored if ``combine`` is ``False``.

        :param bool bootstrap:
            Use :mod:`nestcheck` and :mod:`fgivenx` to bootstrap the runs for
            posterior density error estimation?

        :param bool separate_plots:
            Generate a lower-triangle plot with :mod:`getdist`, and a separate
            error plot with :mod:`nestcheck` (with :mod:`fgivenx` and
            :mod:`getdist`). If ``False`` (default), the diagonal panels of the
            lower-triangle plot are modified by adding the nestcheck output.
            Ignored if ``bootstrap`` is ``False``.

        :param bool write: Export the figure?

        :param str root_filename:
            Root filename to prepend to automatically generated name. Can be,
            e.g., a model and/or data set identifier.

        :param str directory: If ``None`` defaults to current directory.

        :param str ext: File extension for writing. E.g., ``'.png'``.

        :param int dpi: Dots-per-square-inch settings for exporting plots.

        :param kwargs:
            Keyword arguments for the :meth:`_plot_density_with_error` and
            :meth:`_plot_triangle` methods. Keyword arguments for line
            properties (width and alpha) for :mod:`getdist` contours and density
            distributions. If ``bootstrap and not separate_plots`` then
            the density distribution linewidth is set to zero if not
            explicitly specified with kward ``lw_1d``.

        """
        self._set_runs_to_plot(run_IDs, combine, combine_all)
        self._runs.set_params(params)

        if bootstrap_density and not separate_plots:
            if 'lw_1d' not in kwargs: kwargs['lw_1d'] = 0.0
        self._set_line_and_contour_args(**kwargs)

        self._plotter = self._plot_triangle(only_combined,
                                      bootstrap_estimators,
                                      **kwargs)

        if bootstrap_density and separate_plots:
            figs = self._plot_density_with_error(only_combined=only_combined,
                                                 **kwargs)
        elif bootstrap_density:
            figs = self._plot_density_with_error(plotter=self._plotter,
                                                 only_combined=only_combined,
                                                 **kwargs)
        if write:
            root_filename = (root_filename + '__' if root_filename else '') + \
                'posteriorDensity__runs_' + \
                '_'.join(str(ID).replace(' ', '') for
                         ID in self._runs.get_attr('ID')) + '__'


            _dpi = dpi
            if maxdots > 0:
                ndots = dpi * len(self._runs.params.names)
                ndots *= self._plotter.settings.subplot_size_inch
                if ndots > maxdots:
                    dpi = int(maxdots * _dpi / ndots)

            self._plotter.export(fname=root_filename+'triangle'+ext,
                                   adir=directory, dpi=dpi)
            try:
                figs[1].savefig(_os.path.join(directory,
                                              root_filename+'fthetas_1d.pdf'),
                                dpi=_dpi,
                                bbox_inches='tight')
            except IndexError:
                if separate_plots:
                    fname = root_filename + 'params_1d.pdf'
                else:
                    fname = root_filename + 'fthetas_1d.pdf'
                figs[0].savefig(_os.path.join(directory, fname),
                                dpi=_dpi,
                                bbox_inches='tight')
            except (TypeError, NameError):
                pass

        return self._plotter

    @make_verbose('Simulating nested sampling realisations for '
                  'posterior density error visualisation',
                  'Simulated nested sampling realisations and '
                  'plotted posterior density error estimates')
    def _plot_density_with_error(self,
                                 only_combined,
                                 plotter = None,
                                 fthetas = None,
                                 kde_func = None,
                                 kde_settings = None,
                                 **kwargs):
        """
        :param bool only_combined:
            Only plot the combination of all runs?

        :param plotter:
            A :attr:`getdist.GetDistPlotter` instance if the :mod:`nestcheck`
            output is to be displayed on a lower-triangle plot. Assumes the
            parameter order for the ``GetDistPlotter`` is equal to the order
            of parameters in the ``Runs`` instance.

        :param list fthetas:
            Iterable containing functions of parameter vector for which
            density is to be plotted with :mod:`nestcheck` via :mod:`fgivenx`.
            The parameters themselves are handled automatically with ``lambda``
            functions. Additional functions are always plotted using the
            native :mod:`nestcheck` matplotlib figure; the parameter densities
            are be added to a :mod:`getdist` lower-triangle plot is supplied.

        :param func kde_func:
            Function for KDE compatible with :mod:`nestcheck.plots`. Must
            be *weighted* KDE (Higson et al. 2018). If ``None``, uses
            :mod:`getdist` if available, or the native KDE function otherwise. If
            using :mod:`getdist`, the KDE settings are automatically retrieved from
            the first run and applied to :mod:`nestcheck` and :mod:`fgivenx`
            *for all runs*.

        :param kwargs:
            Keyword arguments for :func:`nestcheck.plots.bs_param_dists`.

        TODO:
        ----

            * lims based on credible interval estimate for efficiency?

        """
        try:
            for run in self._runs.subset:
                if not isinstance(run, NSBackend):
                    raise TypeError('Nested sampling backends are required.')
        except AttributeError:
            print('Nested sampling runs are required.')
            raise
        else:
            if only_combined:
                nestcheck_bcknds = self._runs.combined.nestcheck_backend
            else:
                # ensure combined run plotted on top with nestcheck
                try:
                    nestcheck_bcknds = [self._runs.combined.nestcheck_backend]
                except ValueError:
                    nestcheck_bcknds = self._runs.get_attr('nestcheck_backend',
                                                           nestcheck_compatible=True)
                else:
                    nestcheck_bcknds += self._runs.get_attr('nestcheck_backend',
                                                            nestcheck_compatible=True)

        nx = kwargs.pop('nx', 200); ny = kwargs.pop('ny', nx)
        scale_ymax = kwargs.pop('scale_ymax', 1.1)
        n_simulate = kwargs.pop('n_simulate', 200)

        params, idxs = self._runs.params.names, self._runs.params.idxs
        labels = self._runs.params.labels
        bounds = self._runs.params.bounds

        lims = [plotter.subplots[i,i].get_xlim() for i in range(len(params))]
        for l, b in zip(lims, bounds):
            l = list(l)
            l[0] = (l[0] if l[0] > b[0] else b[0])
            l[1] = (l[1] if l[1] < b[1] else b[1])

        if kde_func is None:
            try:
                kde_func = getdist_kde
            except NameError:
                kde_func = weighted_1d_gaussian_kde
                kde_kwargs = None
            else:
                try:
                    kde_settings = [self._runs.combined.kde_settings]
                except ValueError:
                    kde_settings = self._runs.get_attr('kde_settings')
                else:
                    kde_settings += self._runs.get_attr('kde_settings')
                kde_kwargs = {'settings': kde_settings,
                              'ranges': bounds,
                              'normalize': kwargs.pop('normalize', False)}

        lines = kwargs.pop('lines', False)
        parallel = kwargs.pop('parallel', True)
        rasterize_contours = kwargs.pop('rasterize_contours', True)
        tqdm_kwargs = kwargs.pop('tqdm_kwargs', None)

        figsize = kwargs.pop('figsize', _np.array([6.0 * len(params),
                                                   3.0 * len(params)]))

        figs = []
        with verbose(plotter is not None,
                     'Adding density error information to triangle plot',
                     'Added density error information'):
            fig = nestcheck.plots.bs_param_dists(nestcheck_bcknds,
                fthetas=[(lambda y: (lambda theta: theta[:,y]))(i) for i in idxs],
                kde_func=kde_func,
                kde_kwargs=kde_kwargs,
                ftheta_lims=lims,
                nx=nx,
                ny=ny,
                scale_ymax=scale_ymax,
                n_simulate=n_simulate,
                simulate_weights=True,
                getdist_plotter=plotter,
                figsize=figsize,
                lines=lines,
                parallel=parallel,
                rasterize_contours=rasterize_contours,
                labels=labels,
                no_means=True,
                tqdm_kwargs=tqdm_kwargs)

        if fig: figs.append(fig)

        if fthetas:
            try:
                iter(fthetas)
            except TypeError:
                fthetas = [fthetas]
            ftheta_lims = kwargs.pop('ftheta_lims', None)
            if not ftheta_lims:
                raise ValueError('Supply ftheta limits.')
            else:
                kde_kwargs['ranges'] = ftheta_lims
            figsize *= float(len(fthetas))/len(params)
            if 'ftheta_labels' in kwargs:
                kwargs = {'labels': kwargs['ftheta_labels']}
            else:
                kwargs = {}

            fig = nestcheck.plots.bs_param_dists(nestcheck_bcknds,
                                        fthetas=fthetas,
                                        kde_func=kde_func,
                                        kde_kwargs=kde_kwargs,
                                        ftheta_lims=ftheta_lims,
                                        nx=nx,
                                        ny=ny,
                                        scale_ymax=scale_ymax,
                                        n_simulate=n_simulate,
                                        simulate_weights=True,
                                        figsize=figsize,
                                        lines=lines,
                                        parallel=parallel,
                                        rasterize_contours=rasterize_contours,
                                        **kwargs)
            figs.append(fig)

        return figs if figs else None

    @make_verbose('Constructing lower-triangle posterior density '
                  'plot via Gaussian KDE:',
                  'Constructed lower-triangle posterior density plot')
    def _plot_triangle(self,
                       only_combined,
                       bootstrap,
                       prior_density = True,
                       KL_divergence = True,
                       KL_base = 'bits',
                       ndraws = int(1e6),
                       param_plot_lims = None,
                       crosshairs = False,
                       filled = False,
                       legend_loc = 'lower left',
                       legend_corner_coords = (0.75,0.75),
                       legend_frameon = False,
                       scale_attrs = {},
                       normalize = True,
                       veneer = False,
                       no_zero = True,
                       no_ylabel = False,
                       label_right = True,
                       no_ytick = False,
                       credible_interval_1d = True,
                       ID_for_1D_estimators = None,
                       annotate_credible_interval = True,
                       annotate_xy=(0.025,0.915),
                       sixtyeight = True,
                       ninety = False,
                       compute_all_intervals=True,
                       **kwargs):
        """ Call :meth:`getdist.plots.GetDistPlotter.triangle_plot`.

        :param bool only_combined:
            Plot only the combined run?

        :param bool prior_density:
            If ``True`` tries to draw samples from the joint prior and plot
            marginal 1D prior densit functions. Silently fails if attempt
            fails due to missing likelihood and prior callables.

        :param bool KL_divergence:
            If `True` and `prior_density` is `True`, estimate and annotate
            credible interval for Kullback-Leibler divergence for each
            parameter in triangle plot.

        :param str KL_base:
            Base for Kullback-Leibler divergence. Options are {'bits', 'nats'}.

        :param int ndraws:
            Number of samples drawn from the joint prior. Ignored if
            ``prior_density is False`` attempt to plot density fails.

        :param dict param_plot_lims:
            Dictionary of viewing ranges for plotting. Keys must be parameter
            names.

        :param bool crosshairs: Display parameter truth crosshairs?

        :param bool filled: Specify whether the contours are to be filled.

        :param str legend_loc: Specify the legend location. Defaults to
                               ``upper right`` because the plot is a lower-
                               triangle array.

        :param tuple legend_corner_coords:
            Modifies meaning of ``legend_loc`` to be the coordinates of the
            point on the legend box specified by ``legend_loc``. Pass ``None``
            if not applicable. Defaults to place legend in empty upper region
            of lower-triangle plot.

        :param dict scale_attrs:
            Scale :class:`getdist.plots.GetDistPlotterSettings` attributes
            from the automatic values. E.g., ``{'legend_fontsize': 2.0}``.
            Use string values to set the key attribute to the value attribute.
            Caution: do not rely on ordering of pairs in the dictionary.

        :param bool normalize:
            Normalize density distribution in on-diagonal 1D panels?

        :param bool no_zero:
            Hide axes zeros within on-diagonal 1D panels?

        :param bool no_ylabel:
            Hide *probability density* label on diagonal 1D marginal panels?

        :param bool label_right:
            Display *probability density* label on diagonal 1D marginal plots?

        :param bool no_ytick:
            Hide y-axis ticks on diagonal 1D marginal plots?

        :param bool credible_interval_1d:
            Estimate and plot 1D marginal credible intervals? The upper
            and lower quantiles of the interval are estimated via bootstrapping
            with :mod:`nestcheck`, each bounding quantile plotted as a (narrow)
            band bounded by the same quantiles with respect to the bootstrap
            realisations. The interval between the two bands is generally much
            wider and is shaded lighter.

        :param str ID_for_1D_estimators:
            To avoid displaying too much information and resultant plot
            confusion, only one credible interval is displayed per
            on-diagonal plot. If a combined run is found to exist, the
            credible intervals are computed for that combined run; otherwise
            a run ID string needs to be passed using this argument.

        :param bool annotate_credible_interval:
            Annotate each on-diagonal panel with numeric credible interval
            as median +/- distance to symmetric quantiles. Each quantile,
            including the median, is estimated via bootstrapping with
            :mod:`nestcheck`, and the median of each quantile from the bootstrap
            realisations is used for the reported numeric credible interval.

        :param tuple annotate_xy:
            Coordinates as axis fractions for annotation of credible intervals.

        :param bool sixtyeight:
            Should the credible interval, which is symmetric in quantiles about
            the mean, be approximately the 1-\sigma credible interval thus
            containing ~68% of the posterior mass? If ``False`` the interval
            computed is the approximate 2-\sigma interval containing ~95% of
            the posterior mass.

        :param kwargs:

            * additional keyword arguments passed to
              :meth:`getdist.GetDistPlotter.triangle_plot`
            * settings for :mod:`getdist` posterior lower-triangle plotting, applied
              to a :class:`getdist.plots.GetDistPlotSettings` instance

        .. note::
            Using `subplot_size` keyword argument (specify in inches) invokes
            automated labelling fontsizes and tick sizes. If `width_inch` is
            used instead 

        """
        try:
            for run in self._runs.subset:
                if not isinstance(run, NSBackend):
                    raise TypeError('Nested sampling backends are required.')
        except AttributeError:
            print('Nested sampling runs are required.')
            raise
        else:
            if only_combined:
                getdist_bcknds = self._runs.combined.getdist_backend
                legend_labels = None
                line_args = self._runs.combined.lines
                contour_args = self._runs.combined.contours
            else:
                # ensure combined run plotted on top with nestcheck
                try:
                    getdist_bcknds = [self._runs.combined.getdist_backend]
                except ValueError:
                    getdist_bcknds = self._runs.get_attr('getdist_backend')

                    legend_labels = self._runs.get_attr('ID')
                    line_args = self._runs.get_attr('lines')
                    contour_args = self._runs.get_attr('contours')

                    getdist_bcknds.reverse()
                    legend_labels.reverse()
                    line_args.reverse()
                    contour_args.reverse()

                else:
                    getdist_bcknds += self._runs.get_attr('getdist_backend')

                    legend_labels = [self._runs.combined.ID]
                    legend_labels += self._runs.get_attr('ID')
                    line_args = [self._runs.combined.lines]
                    line_args += self._runs.get_attr('lines')
                    contour_args = [self._runs.combined.contours]
                    contour_args += self._runs.get_attr('contours')

                    getdist_bcknds.reverse()
                    legend_labels.reverse()
                    line_args.reverse()
                    contour_args.reverse()

        if param_plot_lims:
            prune = kwargs.get('tick_prune')

        # try to set matching :class:`getdist.plots.GetDistPlotSettings` attrs
        plotter = plots.getSubplotPlotter(kwargs.pop('subplot_size', 2),
                                          kwargs.pop('width_inch', None))
        setattr(plotter.settings, 'progress', True)
        setattr(plotter.settings, 'norm_prob_label', 'Probability density')
        setattr(plotter.settings, 'prob_y_ticks', True)

        for key in kwargs.copy():
            if hasattr(plotter.settings, key):
                setattr(plotter.settings, key, kwargs[key])
                del kwargs[key]

        for key, value in scale_attrs.iteritems():
            if hasattr(plotter.settings, key):
                if isinstance(value, float) or isinstance(value, int):
                    setattr(plotter.settings, key,
                            getattr(plotter.settings, key) * value)
                elif isinstance(value, six.string_types):
                    if hasattr(plotter.settings, value):
                        setattr(plotter.settings, key,
                                getattr(plotter.settings, value))

        if isinstance(normalize, bool):
            diag1d_kwargs = {'normalized': normalize}
        if isinstance(no_zero, bool):
            diag1d_kwargs['no_zero'] = no_zero
        if isinstance(no_ylabel, bool):
            diag1d_kwargs['no_ylabel'] = no_ylabel
        if isinstance(label_right, bool):
            diag1d_kwargs['label_right'] = label_right
        if isinstance(no_ytick, bool):
            diag1d_kwargs['no_ytick'] = no_ytick

        plotter.triangle_plot(getdist_bcknds,
                               legend_labels = legend_labels,
                               params = self._runs.params.names,
                               filled = filled,
                               legend_loc = legend_loc,
                               line_args = line_args,
                               contour_args = contour_args,
                               diag1d_kwargs = diag1d_kwargs,
                               **kwargs)
        try:
            if legend_corner_coords:
                plotter.legend.set_bbox_to_anchor(legend_corner_coords)
        except AttributeError:
            pass
        else:
            plotter.legend.set_frame_on(legend_frameon)

        # add custom parameter plotting limits and updated autolocation
        with fragile(verbose(param_plot_lims,
                             'Applying bespoke parameter viewing intervals',
                             'Viewing intervals applied')) as condition:
            if not condition: fragile.Break

            params = self._runs.params
            for param, l in param_plot_lims.items():
                j = params.names.index(param)
                for i in range(j,len(params.names)):
                    ax = plotter.subplots[i,j]
                    ax.set_xlim(l)
                    ax.xaxis.set_major_locator(_get_default_locator(prune))
                    ax.xaxis.set_minor_locator(AutoMinorLocator())
                for i in range(j):
                    ax = plotter.subplots[j,i]
                    ax.set_ylim(l)
                    ax.yaxis.set_major_locator(_get_default_locator(prune))
                    ax.yaxis.set_minor_locator(AutoMinorLocator())

        if prior_density:
            self._add_prior_density(plotter, ndraws, normalize,
                                    KL_divergence, KL_base,
                                    ID=ID_for_1D_estimators,
                                    bootstrap=bootstrap,
                                    n_simulate = kwargs.get('n_simulate'))
        if veneer:
            self._veneer_spines_ticks(plotter, **kwargs)
        if crosshairs:
            self._add_crosshairs(plotter, self._runs.params.truths)
        if credible_interval_1d: # include nestcheck estimator bootstrap error
            self._add_credible_interval(plotter,
                                        bootstrap=bootstrap,
                                        n_simulate=kwargs.get('n_simulate'),
                                        ID=ID_for_1D_estimators,
                                        annotate=annotate_credible_interval,
                                        annotate_xy=annotate_xy,
                                        sixtyeight=sixtyeight,
                                        ninety=ninety,
                                        compute_all_intervals=compute_all_intervals)

        self._plotter = plotter
        return plotter

    @make_verbose('Adding 1D marginal prior density functions',
                  'Added 1D marginal prior density functions')
    def _add_prior_density(self, plotter, ndraws, normalize,
                           KL_divergence, KL_base,
                           ID, bootstrap, n_simulate):
        """ Crudely estimate the prior density.

        Kullback-Leibler divergence estimated in bits for a combined run or
        the same run for which the credible intervals are calculated.

        .. todo::

            * Apply GetDist KDE to 1D marginal distributions manually so that
              density at boundaries is accurate. SciPy probably will not
              handle boundary correction and density is often appreciable
              at boundaries for weakly informative prior density functions.
            * Plot multiple priors if there are multiple likelihood functions?
            * Add prior-to-posterior KL-divergence functionality for 1D
              marginal distributions using Guassian KDE for prior density
              function.

        """
        try:
            l = self._likelihood[self._runs.subset[0].ID]
        except TypeError:
            l = self._likelihood

        if l is None:
            return
        elif not hasattr(l, 'prior'):
            return
        elif not callable(l.prior):
            return
        elif not hasattr(l.prior, 'inverse_sample'):
            return
        elif not callable(l.prior.inverse_sample):
            return

        h = _np.random.rand(int(ndraws), int(l.num_params))

        try:
            p = l.prior.inverse_sample_and_transform(h[0, :])
        except AttributeError:
            p = l.prior.inverse_sample(h[0, :])

        try:
            samples = _np.zeros((int(ndraws), len(p)),
                                dtype=_np.double)
        except TypeError:
            return

        finite_counter = counter = index = 0
        while finite_counter < ndraws:
            try:
                try:
                    p = l.prior.inverse_sample_and_transform(h[index, :])
                except AttributeError:
                    p = l.prior.inverse_sample(h[index, :])
            except IndexError: # use estimate to draw more from hypercube
                redraw = float(counter) * ndraws / finite_counter
                redraw -= finite_counter
                h = _np.random.rand(int(redraw)+1, int(l.num_params))
                index = 0
                try:
                    p = l.prior.inverse_sample_and_transform(h[index, :])
                except AttributeError:
                    p = l.prior.inverse_sample(h[index, :])

            if _np.isfinite(l.prior(p)):
                samples[finite_counter,:] = p
                finite_counter += 1

            counter += 1
            index += 1

        try:
            run = self._runs.combined
        except (AttributeError, ValueError):
            run = self._runs[ID]
        color, lw = (run.contours[key] for key in ('color', 'lw'))

        quantiles = [None] * 3

        with verbose(KL_divergence,
                     'Estimating 1D marginal KL-divergences in %s' % KL_base,
                     'Estimated 1D marginal KL-divergences') as condition:
            for i, ax in enumerate([plotter.subplots[i,i] \
                                for i in range(plotter.subplots.shape[0])]):

                name = self._runs.params.names[i]
                bounds = {name: self._runs.params.bounds[i]}
                settings = {'fine_bins': 1024,
                            'smooth_scale_1D': 0.3,
                            'boundary_correction_order': 1,
                            'mult_bias_correction_order': 1} # adopt from posterior settings or take custom input?

                bcknd = getdist.mcsamples.MCSamples(sampler='uncorrelated',
                                samples=samples[:,self._runs.params.idxs[i]],
                                weights=None,
                                names=[name],
                                ranges=bounds,
                                settings=settings)

                if normalize:
                    bcknd.get1DDensity(name).normalize(by='integral',
                                                       in_place=True)

                x = _np.linspace(ax.xaxis.get_view_interval()[0],
                                 ax.xaxis.get_view_interval()[1],
                                 1000)

                ax.plot(x, bcknd.get1DDensity(name).Prob(x),
                        ls='-.', color=color, lw=lw)

                if not condition: continue

                # a prototype Kullback-Leibler divergence callback
                # information in bits
                def KL(ns_run, logw):
                    x = ns_run['theta'][:,self._runs.params.idxs[i]]
                    w_rel = _np.exp(logw - logw.max())
                    where = w_rel > run.kde_settings.get('min_weight_ratio',
                                                         1.0e-30)
                    prior = bcknd.get1DDensity(name).Prob(x[where])
                    posterior = getdist_kde(x[where], x, w_rel,
                                        ranges=[self._runs.params.bounds[i]],
                                        idx=0,
                                        normalize=normalize,
                                        settings=run.kde_settings)
                    # Due to spline interpolation, very small densities can be
                    # negative, so manually give a small postive value which
                    # does not affect KL integral approximation
                    posterior[posterior<=0.0] = posterior[posterior>0.0].min()

                    KL = _np.sum(w_rel[where] \
                                   * (_np.log(posterior) - _np.log(prior))) \
                                   /_np.sum(w_rel[where])

                    if KL_base == 'bits':
                        return KL / _np.log(2.0)
                    elif KL_base == 'nats':
                        return KL
                    else:
                        raise ValueError('Invalid base for KL-divergence.')

                if bootstrap:
                    for j, cred_int in enumerate([0.025, 0.5, 0.975]):
                        quantiles[j] = run_ci_bootstrap(run.nestcheck_backend,
                                                     estimator_list=[KL],
                                                     cred_int=cred_int,
                                                     n_simulate=n_simulate,
                                                     simulate_weights=True,
                                                     flip_skew=True)
                    # KL in bits
                    interval = r'$D_{\mathrm{KL}}=%.2f_{-%.2f}^{+%.2f}$' \
                                                  % (quantiles[1],
                                                     quantiles[1] - quantiles[0],
                                                     quantiles[2] - quantiles[1])

                    yield ('%s KL-divergence = %.4f/-%.4f/+%.4f'
                            % (name,
                               quantiles[1],
                               quantiles[1] - quantiles[0],
                               quantiles[2] - quantiles[1]))

                    if not rcParams['text.usetex']:
                        fontsize = plotter.settings.lab_fontsize - 1
                    else:
                        fontsize = plotter.settings.lab_fontsize

                    ax.set_title(interval, color=color,
                                 fontsize=fontsize)
                else:
                    where = run.samples[:,0] > 0.0
                    _ = run.samples[where,2:]
                    __ = _np.log(run.samples[where,0])

                    ns_run = {'theta': _}
                    divergence = KL(ns_run, __)

                    divergence = (r'$D_{\mathrm{KL}}=%.2f$' % divergence)

                    yield ('%s %s' % (name, divergence))

                    if not rcParams['text.usetex']:
                        fontsize = plotter.settings.lab_fontsize - 1
                    else:
                        fontsize = plotter.settings.lab_fontsize

                    ax.set_title(divergence, color=color,
                                 fontsize=fontsize)

        yield None

    def KL_divergence(self, base='bits', bootstrap=False,
                      ID=None, n_simulate=200):
        """ Kullback-Leibler divergence integral jointly for all parameters. """
        if ID is None:
            run = self._runs.combined
        else:
            run = self._runs[ID]

        def estimator(ns_run, logw):
            w_rel = _np.exp(logw - logw.max())
            KL = _np.sum(w_rel * ns_run['logl']) / _np.sum(w_rel)
            KL -= logsumexp(logw)

            if base == 'bits':
                return KL / _np.log(2.0)
            elif base == 'nats':
                return KL
            else:
                raise ValueError('Invalid base for KL-divergence.')

        if bootstrap:
            quantiles = [None] * 3
            for j, cred_int in enumerate([0.025, 0.5, 0.975]):
                quantiles[j] = run_ci_bootstrap(run.nestcheck_backend,
                                                 estimator_list=[estimator],
                                                 cred_int=cred_int,
                                                 n_simulate=n_simulate,
                                                 simulate_weights=True,
                                                 flip_skew=True)
            return quantiles
        else:
            divergence = estimator(run.nestcheck_backend,
                                   get_logw(run.nestcheck_backend))

        return divergence

    @make_verbose('Adding 1D marginal credible intervals',
                  'Added 1D marginal credible intervals')
    def _add_credible_interval(self, plotter, bootstrap, n_simulate,
                               ID, annotate, annotate_xy, sixtyeight,
                               ninety, compute_all_intervals):
        """
        Estimate 1-\sigma credible interval in one-dimension on a
        combined run, or if such a run does not exist, on the run with
        the specified ID.

        Calls :func:`nestcheck.estimators.param_cred` for one-tailed weighted
        estimate; two such estimates give a credible interval which is
        symmetric in quantiles with respect to the median. Also calls
        :func:`nestcheck.error_analysis.run_ci_bootstrap` for
        credible interval on quantiles.

        :param bool sixtyeight:
            Plot the 68% credible interval about the median in 1D plots? If
            ``False`` plots 95% credible interval about the median -- i.e.,
            symmetric quantiles about the median.
        """
        diag = [plotter.subplots[i,i] for i in range(plotter.subplots.shape[0])]

        # get the color of the run if ID is used, else red for combined
        try:
            run = self._runs.combined
            color = self._runs.combined.lines['color']
        except (AttributeError, ValueError):
            run = self._runs[ID]
            color = self._runs[ID].lines['color']

        # estimator requires closure to be changable
        def get_estimator(quantile, param_ind):
            def estimator(*args, **kwargs):
                return param_cred(*args,
                                  probability=quantile,
                                  param_ind=param_ind,
                                  **kwargs)
            return estimator

        quantiles = [0.159, 0.5, 0.841] if sixtyeight else ([0.05,0.5,0.95] if ninety else [0.025, 0.5, 0.975])

        def format_CI(name, cred, summary):
            if len(cred.shape) > 1:
                msg = ('%s CI_{%i\%%} = %.4f/-%.4f/+%.4f'
                        % (name,
                           summary,
                           cred[1,1],
                           cred[1,1] - cred[0,1],
                           cred[2,1] - cred[1,1]))
            else:
                msg = ('%s CI_{%i\%%} = %.4f/-%.4f/+%.4f'
                        % (name,
                           summary,
                           cred[1],
                           cred[1] - cred[0],
                           cred[2] - cred[1]))


            return msg

        if bootstrap:
            for i, (ax, ind) in enumerate(zip(diag, self._runs.params.idxs)):
                def calculate_intervals(quantiles):
                    cred = _np.zeros((len(quantiles), len(quantiles)), dtype=_np.double)
                    for j, p in enumerate(quantiles):
                        for k, q in enumerate(quantiles):
                            cred[j,k] = run_ci_bootstrap(run.nestcheck_backend,
                                             estimator_list=[get_estimator(p, ind)],
                                             cred_int=q,
                                             n_simulate=n_simulate,
                                             simulate_weights=True,
                                             flip_skew=True)[0]
                    return cred

                cred = calculate_intervals(quantiles)
                zorder = max([_.zorder for _ in ax.get_children()]) + 1

                ax.axvspan(cred[0,0], cred[0,2], alpha=0.5,
                           facecolor=color,
                           edgecolor=color,
                           linewidth=0,
                           rasterized=True,
                           zorder=zorder)

                ax.axvspan(cred[2,0], cred[2,2], alpha=0.5,
                           facecolor=color,
                           edgecolor=color,
                           linewidth=0,
                           rasterized=True,
                           zorder=zorder)

                ax.axvspan(cred[0,2], cred[2,0], alpha=0.25,
                           facecolor=color,
                           edgecolor=color,
                           linewidth=0,
                           rasterized=True,
                           zorder=zorder)

                if annotate:
                    stats = r'CI$_{%d\%%}=%.2f_{-%.2f}^{+%.2f}$' \
                                % (68 if sixtyeight else (90 if ninety else 95),
                                   cred[1,1],
                                   cred[1,1] - cred[0,1],
                                   cred[2,1] - cred[1,1])
                    title = ax.get_title()
                    if title:
                        title = stats.center(30) + '\n' + title.center(30)
                    else:
                        title = stats

                    if not rcParams['text.usetex']:
                        fontsize = plotter.settings.lab_fontsize - 1
                    else:
                        fontsize = plotter.settings.lab_fontsize

                    ax.set_title(title, color=color,
                                 fontsize=fontsize)

                if compute_all_intervals:
                    yield format_CI(self._runs.params.names[i],
                                    cred,
                                    68 if sixtyeight else (90 if ninety else 95))

                    if sixtyeight:
                        yield format_CI(self._runs.params.names[i],
                                        calculate_intervals([0.05, 0.5, 0.95]),
                                        90)
                        yield format_CI(self._runs.params.names[i],
                                        calculate_intervals([0.025, 0.5, 0.975]),
                                        95)
                    elif ninety:
                        yield format_CI(self._runs.params.names[i],
                                        calculate_intervals([0.159, 0.5, 0.841]),
                                        68)
                        yield format_CI(self._runs.params.names[i],
                                        calculate_intervals([0.025, 0.5, 0.975]),
                                        95)
                    else:
                        yield format_CI(self._runs.params.names[i],
                                        calculate_intervals([0.159, 0.5, 0.841]),
                                        68)
                        yield format_CI(self._runs.params.names[i],
                                        calculate_intervals([0.05, 0.5, 0.95]),
                                        90)

        else:
            for i, (ax, ind) in enumerate(zip(diag, self._runs.params.idxs)):
                def calculate_intervals(quantiles):
                    cred = _np.zeros(len(quantiles), dtype=_np.double)
                    for j, p in enumerate(quantiles):
                        where = run.samples[:,0] > 0.0
                        _ = run.samples[where,2:]
                        __ = _np.log(run.samples[where,0])
                        cred[j] = get_estimator(p, ind)({'theta': _}, __)

                    return cred

                cred = calculate_intervals(quantiles)
                zorder = max([_.zorder for _ in ax.get_children()]) + 1

                ax.axvspan(cred[0], cred[2], alpha=0.25,
                           facecolor=color,
                           edgecolor=color,
                           linewidth=0,
                           rasterized=True,
                           zorder=zorder)

                if annotate:
                    stats = r'CI$_{%d\%%}=%.2f_{-%.2f}^{+%.2f}$' \
                                                % (68 if sixtyeight else 95,
                                                   cred[1],
                                                   cred[1] - cred[0],
                                                   cred[2] - cred[1])

                    title = ax.get_title()
                    if title:
                        title = stats + r'  $|$  ' + title
                    else:
                        title = stats

                    if not rcParams['text.usetex']:
                        fontsize = plotter.settings.lab_fontsize - 1
                    else:
                        fontsize = plotter.settings.lab_fontsize

                    ax.set_title(title, color=color,
                                 fontsize=fontsize)

                if compute_all_intervals:
                    yield format_CI(self._runs.params.names[i],
                                    cred,
                                    68 if sixtyeight else (90 if ninety else 95))

                    if sixtyeight:
                        yield format_CI(self._runs.params.names[i],
                                        calculate_intervals([0.05, 0.5, 0.95]),
                                        90)
                        yield format_CI(self._runs.params.names[i],
                                        calculate_intervals([0.025, 0.5, 0.975]),
                                        95)
                    elif ninety:
                        yield format_CI(self._runs.params.names[i],
                                        calculate_intervals([0.159, 0.5, 0.841]),
                                        68)
                        yield format_CI(self._runs.params.names[i],
                                        calculate_intervals([0.025, 0.5, 0.975]),
                                        95)
                    else:
                        yield format_CI(self._runs.params.names[i],
                                        calculate_intervals([0.159, 0.5, 0.841]),
                                        68)
                        yield format_CI(self._runs.params.names[i],
                                        calculate_intervals([0.05, 0.5, 0.95]),
                                        90)

        yield None

    @staticmethod
    @make_verbose('Adding parameter truth crosshairs',
                  'Added crosshairs')
    def _add_crosshairs(plotter, truths):
        """ Add parameter crosshairs to triangle plot. """
        spine = next(plotter.subplots[0,0].spines.itervalues())
        lw = spine.get_linewidth()
        for i, truth in enumerate(truths):
            if truth is not None:
                for ax in plotter.subplots[:,i]:
                    if ax is not None:
                        ax.axvline(truth, color='black', ls='-', lw=lw)
                if i > 0:
                    for ax in plotter.subplots[i,:i]:
                        if ax is not None:
                            ax.axhline(truth, color='black', ls='-', lw=lw)

    @staticmethod
    @make_verbose('Veneering spines and axis ticks',
                  'Veneered')
    def _veneer_spines_ticks(plotter, lengthen=2.0, embolden=2.0,
                             **kwargs):
        """ Embolden spines, and embolden and lengthen ticks. """

        ax = plotter.subplots[0,0]

        major_length = ax.xaxis.majorTicks[0].tick1line.get_markersize()
        major_length *= lengthen
        minor_length = ax.xaxis.minorTicks[0].tick1line.get_markersize()
        minor_length *= lengthen
        lw = ax.xaxis.majorTicks[0].tick1line.get_markeredgewidth()
        lw *= embolden

        for i in range(plotter.subplots.shape[0]):
            for j in range(i+1):
                ax = plotter.subplots[i,j]
                ax.tick_params(which='major', colors='black', length=major_length)
                ax.tick_params(which='minor', colors='black', length=minor_length)
                ax.xaxis.set_tick_params(which='both', width=lw)
                ax.yaxis.set_tick_params(which='both', width=lw)
                for spine in ax.spines.itervalues():
                    spine.set_linewidth(lw)

    @fix_random_seed
    @make_verbose('Plotting data and model for posterior checking',
                  'Plotted data and model for posterior checking')
    def plot_pulse_and_spectrum(self,
                                run_IDs,
                                combine,
                                combine_all,
                                nsamples,
                                num_phases=1000,
                                num_energies=1000,
                                cache=True,
                                write=True,
                                fig_dir = './',
                                root_filename = '',
                                extension = '.pdf',
                                dpi=300,
                                use_fgivenx=False,
                                random_seed=0,
                                absorbed_spectrum=False,
                                **kwargs):
        """
        Plot data and model for posterior checking.

        If a combined run is available, only this run is used. Otherwise
        a plot is generated for each run using the corresponding likelihood
        object. There may in principle be multiple :class:`~.Pulse.Pulse`
        instances per likelihood object, but only the first instance is used.

        Compute and plot posterior-averaged channel-summed count-rate pulse.
        The figure contains four panels which share phase as an x-axis:

            * the first (topmost) panel displays the data in joint channel-phase
              intervals
            * the second panel displays the posterior mean expectation in joint
              channel-phase intervals
            * the third panel displays the standardised residuals in joint
              channel-phase intervals
            * the last (bottommost) panel displays the channel-summed pulse as
              a function of phase for a subset of samples, optionally using
              :mod:`fgivenx`; the injected pulse is plotted by default if the
              injection (truth) vector is known.

        :param bool cache:
            Cache intermediary model objects to accelerate post-processing?

        :param int nsamples:
            Number of samples to use. Equally-weighted samples are generated,
            thus introducing a additional Monte Carlo noise which is ignored.

        :param int num_phases:
            Number of phases to interpolate at on the interval [0,2] cycles.

        :param bool use_fgivenx:
            Use :mod:`fgivenx` to plot iso-probability credible interval
            contours in target-source count rate as a function of phase.

        :param int random_seed:
            Seed :mod:`numpy.random` for determinism.

        :param kwargs:
            Keyword arguments for the following:

                * :meth:`_SignalPlot.add_expected_source`

        .. todo::

            * Use tqdm for progress bar.
            * Bootstrap and weight simulation via nestcheck for additional
              variation?
            * Return plot objects for the user handling?

        """

        self._set_runs_to_plot(run_IDs, combine, combine_all)

        if combine:
            try:
                l = self._likelihood[run_IDs[0]]
            except TypeError:
                l = self._likelihood

            self._signalPlot = _SignalPlot(l.pulses[0],
                                           int(num_phases),
                                           int(num_energies),
                                           use_fgivenx, **kwargs)

            self._pulse_and_spectrum_driver(self._runs.combined.ID,
                                            l, self._runs.combined.samples,
                                            nsamples, self._runs.combined.truth_vector,
                                            self._signalPlot, use_fgivenx, cache,
                                            absorbed_spectrum,
                                            **kwargs)

            root_filename = (root_filename + '__' if root_filename else '') + \
                            str(self._runs.combined.ID).replace(' ', '')

            if write:
                self._signalPlot.savefig(fig_dir, root_filename, extension,
                                         dpi=dpi, bbox_inches='tight')

        else:
            for run in self._runs.subset:
                try:
                    l = self._likelihood[run.ID]
                except TypeError:
                    l = self._likelihood

                self._signalPlot = _SignalPlot(l.pulses[0],
                                               int(num_phases),
                                               int(num_energies),
                                               use_fgivenx,
                                               **kwargs)

                self._pulse_and_spectrum_driver(run.ID, l,
                                                run.samples, nsamples,
                                                run.truth_vector,
                                                self._signalPlot,
                                                use_fgivenx, cache,
                                                absorbed_spectrum,
                                                **kwargs)

                if write:
                    self._signalPlot.savefig(fig_dir, run.ID, extension,
                                             dpi=dpi, bbox_inches='tight')

    @property
    def signalPlot(self):
        return self._signalPlot

    @signalPlot.setter
    def signalPlot(self, obj):
        if isinstance(obj, _SignalPlot):
            self._signalPlot = obj

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

    def _pulse_and_spectrum_driver(self, ID, likelihood, samples, nsamples,
                                   truths, signalplot,
                                   use_fgivenx, cache,
                                   absorbed_spectrum, **kwargs):
        """ Execute plotting loop given samples.

        .. todo::
            Use tqdm for progress bar whilst caching.

        """

        thetas = self._draw_equally_weighted(samples, nsamples,
                                             likelihood.num_params)

        pulse = signalplot.pulse

        if cache and h5py is not None:
            filename = ID + '__pulse_' + pulse.tag
            cache = _Cache(filename, **kwargs)
            if cache.do_caching(thetas): # skips body if can simply read cache
                for i in range(thetas.shape[0]):
                    likelihood(thetas[i,:])
                    cache.cache(pulse.caching_targets)

        elif cache and h5py is None:
            raise ImportError('You need to install h5py to use caching.')

        def update(theta):
            if cache: # use the cache if available
                cached = next(cache)
                #unique_keys = []
                for key, value in cached.iteritems():
                    try:
                        delattr(pulse, key)
                    except AttributeError:
                        pass

                    if len(value.shape) == 3:
                        for i in range(value.shape[0]):
                            setattr(pulse, key, value[i,...])
                    else:
                        setattr(pulse, key, value)
            else: # otherwise resort to likelihood evaluations
                likelihood(theta)

        if use_fgivenx:
            def wrap(index):
                if cache:
                    cache.reset_iterator()

                try:
                    del signalplot.model_sum
                except AttributeError:
                    pass
                def callback(x, theta):
                    # Ignore x because it is already known by a plot object in
                    # this more object oriented approach to calling fgivenx
                    update(theta)
                    return signalplot.update_expectations(absorbed_spectrum, **kwargs)[index]
                return callback

            signalplot.add_source_contours(wrap(0), thetas, **kwargs)
            signalplot.add_folded_contours(wrap(1), thetas, **kwargs)
            if not absorbed_spectrum:
                signalplot.add_spectrum_contours(wrap(2), thetas, **kwargs)
            else:
                signalplot.add_spectrum_contours(wrap(3), thetas, **kwargs)

        else:
            for i in range(thetas.shape[0]):
                update(thetas[i,:])
                signalplot.update_expectations(**kwargs)

        if None not in truths:
            likelihood(truths[:likelihood.num_params])

        signalplot.add_data()
        signalplot.add_expected_counts()
        signalplot.add_expected_source(absorbed_spectrum=absorbed_spectrum,
                                       **kwargs)

class _Cache(object):
    """ Cache numerical model objects computed during likelihood evaluation.


    """

    def __init__(self, filename, base_dir='./',
                 read_only=False, archive=True, **kwargs):

        if isinstance(filename, six.string_types):
            if filename[-3:] != '.h5':
                self._filename = filename + '.h5'
            else:
                self._filename = filename

        self._base_dir = base_dir
        self._path = _os.path.join(self._base_dir, self._filename)
        self._read_only = read_only
        self._archive_if_incompatible = archive

    def __enter__(self):
        return self

    def __exit__(self, exc, exc_value, traceback):
        if exc:
            print('Encountered problem whilst caching:')

    def _open(self, mode='r'):
        """ Get the :mod:`h5py` context manager. """
        if self._read_only and mode != 'r':
            raise RuntimeError('The cache is in read-only mode.')
        return h5py.File(self._path, mode)

    def cache(self, data):
        """ Cache the computational data. """

        with self._open('r+') as f:
            g = f['data']
            for key, value in data.iteritems():
                if isinstance(value, tuple) or isinstance(value, list):
                    if key not in g.keys():
                        shape = [f.attrs['n'], len(value)]
                        shape += [s for s in value[0].shape]
                        g.create_dataset(key, shape=shape, dtype='float64')

                    for j, v in enumerate(value):
                        g[key][self.i,j,...] = v
                else:
                    if key not in g.keys():
                        shape = [f.attrs['n']] + [s for s in value.shape]
                        g.create_dataset(key, shape=shape, dtype='float64')

                    g[key][self.i,...] = value

        self.i += 1

    def reset_iterator(self):
        """ Reset the counter for the cache iterator. """
        self.i = 0

    def __iter__(self):
        self.reset_iterator()
        return self

    def __next__(self):
        """ Read from the cache. """

        cached = {}

        with self._open('r') as f:
            g = f['data']
            for key in g.keys():
                cached[key] = g[key][self.i,...]

        self.i += 1

        return cached

    def next(self):
        """ Python 2.x compatibility. """
        return self.__next__()

    @make_verbose('Checking whether an existing cache can be read:',
                  'Cache state determined')
    def do_caching(self, samples):
        """ Check whether a new cache is required or whether an exising
            cache can be read without additional computation.

        :return: Boolean indicating whether to read (``False``) or write.

        """

        try: # try reading file and checking keys
            with self._open('r') as f:
                if 'thetas' not in f.keys():
                    self._new(samples)
                    return True
        except IOError: # create new cache file
            self._new(samples)
            return True
        else: # can be read, so check if samples array are matching
            if self._changed(samples):
                self._new(samples)
                return True
            else:
                return False

    @make_verbose('Creating new cache file', 'Cache file created')
    def _new(self, samples):
        """ Prepare a new cache file. """

        if not _os.path.isdir(self._base_dir):
            _os.mkdir(self._base_dir)

        if self._archive_if_incompatible:
            try:
                with self._open('r'):
                    pass
            except IOError:
                self._initialise(samples)
            else:
                self._archive()
                self._initialise(samples)
        else:
            self._initialise(samples)

    @make_verbose('Initialising cache file', 'Cache file initialised')
    def _initialise(self, samples):
        """ Initialise the cache. """

        with self._open('w') as f:
            f.attrs['version'] = __version__
            f.attrs['n'] = samples.shape[0]
            f.create_dataset('thetas', data=samples)
            f.create_group('/data')

        self.reset_iterator()

    def _changed(self, samples):
        """ Check whether sample set has changed. """
        with self._open('r') as f:
            if f.attrs['version'] != __version__:
                return True
            if not _np.array_equal(f['thetas'], samples):
                return True
        return False

    @make_verbose('Attempting to archive existing cache file in '
                  'a subdirectory')
    def _archive(self):
        """ Archive an existing cache file. """

        # to archive the existing cache file
        archive_dir = _os.path.join(self._base_dir, 'archive')

        try:
            if not _os.path.isdir(archive_dir):
                _os.mkdir(archive_dir)
        except OSError:
            yield ('Archiving failed... cache file %s will be '
                   'overwritten.' % self._filename)
            yield
        else:
            yield 'Targeting subdirectory: %s.' % archive_dir

        try:
            from datetime import datetime
        except ImportError:
            yield ('Archiving failed... cache file %s will be '
                   'overwritten.' % self._filename)
            yield
        else:
            name_archived = self._filename[:-3] + '__archive__'
            name_archived += 'xpsi_version_%s__' % __version__
            obj = datetime.now()
            name_archived += 'datetime__%i.%i.%i__%i.%i.%i' % (obj.day,
                                                               obj.month,
                                                               obj.year,
                                                               obj.hour,
                                                               obj.minute,
                                                               obj.second)

            try:
                _os.rename(self._filename,
                           _os.path.join(archive_dir, name_archived + '.h5'))
            except OSError:
                yield ('Archiving failed... cache file %s will be '
                       'overwritten.' % self._filename)
            else:
                yield ('Exisiting cache file archived in '
                       'subdirectory %s.' % archive_dir)

        yield None

class _SignalPlot(object):
    """ A signal plot.

    Transform a single :class:`.Pulse.Pulse` instance into crude posterior
    predictive form.

    Compute and plot posterior-averaged channel count-rate spectrum.

    :param obj pulse: An instance of :class:`~.Pulse.Pulse`.

    :param bool rasterized:
        Beware if ``False``white polygon edges may appear in figures.

    :param tuple figsize:
        The size, ``(width, height)``, in inches per panel. Each panel
        extends over two phase cycles.

    .. todo::
        Set figure size and label fontsizes automatically.

    """

    @make_verbose('Instantiating a pulse plotter for posterior checking',
                  'Pulse plotter instantiated')
    def __init__(self, pulse, num_phases=1000, num_energies=1000,
                 use_fgivenx=False, rasterized=True,
                 figsize=(8,4), **kwargs):

        if figsize[1] >= figsize[0]:
            yield ('Warning: pulse plots are generally better rendered such '
                   'that panel width is greater than panel height. Use the '
                   '``figsize`` argument to specify figure size per panel.')

        fscale = kwargs.pop('fscale', 1.0)
        fontsize = fscale * 3.25 * min(figsize[0], figsize[1])
        rc('font',**{'size': fontsize})

        self.pulse = pulse

        # figure for data, posterior expected model counts, and residuals
        self._fig_data = plt.figure(figsize = (figsize[0], figsize[1] * 3))
        self._gs_data = gridspec.GridSpec(3, 2,
                                          height_ratios=[1]*3,
                                          width_ratios=[50,1])
        self._gs_data.update(wspace=0.025, hspace=0.125*fscale)

        #temp refs for readability:
        fig = self._fig_data
        gs = self._gs_data

        self._ax_data = fig.add_subplot(gs[0,0])
        self._ax_data_cb = fig.add_subplot(gs[0,1])

        self._ax_model = fig.add_subplot(gs[1,0])
        self._ax_model_cb = fig.add_subplot(gs[1,1])

        self._ax_resid = fig.add_subplot(gs[2,0])
        self._ax_resid_cb = fig.add_subplot(gs[2,1])

        # figure for source signals, without interstellar effects, and folded
        # but no background
        self._fig_source = plt.figure(figsize = (figsize[0], figsize[1] * 4))
        self._gs_source = gridspec.GridSpec(4, 2,
                                            height_ratios=[1]*4,
                                            width_ratios=[50,1])
        self._gs_source.update(wspace=0.025, hspace=0.125*fscale)

        # temp refs for readability:
        fig = self._fig_source
        gs = self._gs_source

        self._ax_source = fig.add_subplot(gs[0,0])
        self._ax_source_cb = fig.add_subplot(gs[0,1])

        self._ax_source_1d = fig.add_subplot(gs[1,0])
        if use_fgivenx:
            self._ax_source_1d_cb = fig.add_subplot(gs[1,1])

        self._ax_folded = fig.add_subplot(gs[2,0])
        self._ax_folded_cb = fig.add_subplot(gs[2,1])

        self._ax_folded_1d = fig.add_subplot(gs[3,0])
        if use_fgivenx:
            self._ax_folded_1d_cb = fig.add_subplot(gs[3,1])

        # figure for source signals, without interstellar effects, and folded
        # and background
        self._fig_spectrum = plt.figure(figsize = (figsize[0], figsize[1] * 3))
        self._gs_spectrum = gridspec.GridSpec(2, 1,
                                             height_ratios=[1,2])
        self._gs_spectrum.update(hspace=0.175*fscale)

        # temp refs for readability
        fig = self._fig_spectrum
        gs = self._gs_spectrum

        gs_top = gridspec.GridSpecFromSubplotSpec(1, 2,
                                                  subplot_spec=gs[0,0],
                                                  wspace=0.025,
                                                  width_ratios=[50,1])
        self._ax_spectrum_incident = fig.add_subplot(gs_top[0,0])
        if use_fgivenx:
            self._ax_spectrum_incident_cb = fig.add_subplot(gs_top[0,1])

        gs_bottom = gridspec.GridSpecFromSubplotSpec(2, 2,
                                                     subplot_spec=gs[1,0],
                                                     wspace=0.025,
                                                     hspace=0.125*fscale,
                                                     height_ratios=[1]*2,
                                                     width_ratios=[50,1])

        self._ax_spectrum_folded = fig.add_subplot(gs_bottom[0,0])
        self._ax_spectrum_folded_cb = fig.add_subplot(gs_bottom[0,1])

        self._ax_spectrum_folded_1d = fig.add_subplot(gs_bottom[1,0])

        # focus only on axes objects now
        # axes objects are split over figures

        axes = [self._ax_data,
                self._ax_model,
                self._ax_resid,
                self._ax_source,
                self._ax_source_1d,
                self._ax_folded,
                self._ax_folded_1d]

        saxes = [self._ax_spectrum_incident,
                 self._ax_spectrum_folded,
                 self._ax_spectrum_folded_1d]

        cbs = [self._ax_data_cb,
               self._ax_model_cb,
               self._ax_resid_cb,
               self._ax_source_cb,
               self._ax_folded_cb,
               self._ax_spectrum_folded_cb]

        try:
            cbs.append(self._ax_source_1d_cb)
            cbs.append(self._ax_folded_1d_cb)
            cbs.append(self._ax_spectrum_incident_cb)
        except AttributeError:
            pass

        # these are globally applicable settings
        for ax in axes + saxes + cbs:
            ax.tick_params(axis='both', which='major',
                           colors='black', length=8, direction='out')
            ax.tick_params(axis='both', which='minor',
                           colors='black', length=4, direction='out')
            ax.xaxis.set_tick_params(which='both', width=1.0)
            ax.yaxis.set_tick_params(which='both', width=1.0)
            for spine in ax.spines.itervalues():
                spine.set_linewidth(1.0)

        # properties for phase each axis
        for ax in axes:
            ax.xaxis.set_major_locator(MultipleLocator(0.2))
            ax.xaxis.set_minor_locator(MultipleLocator(0.05))
            ax.set_xlim([0.0,2.0])
            if ax not in (self._ax_resid, self._ax_folded_1d):
                ax.tick_params(axis='x', labelbottom=False)
            else:
                ax.set_xlabel('$\phi$ [cycles]')

        # set more specific properties
        for ax in (self._ax_data, self._ax_model, self._ax_resid,
                   self._ax_folded):
            ax.set_ylabel('Channel (PI)')

        self._ax_source.set_ylabel(r'$E$ [keV]')
        self._ax_source_1d.set_ylabel(r'photons/cm$^{2}$/s')
        self._ax_folded_1d.set_ylabel('counts/s')

        self._ax_spectrum_incident.set_xlabel(r'$E$ [keV]')
        self._ax_spectrum_incident.set_ylabel(r'photons/keV/cm$^{2}$/s')

        self._ax_spectrum_incident.set_xlim([self.pulse.energy_edges[0],
                                             self.pulse.energy_edges[-1]])
        self._ax_spectrum_incident.set_xscale('log')
        self._ax_spectrum_incident.set_yscale('log')

        self._ax_spectrum_folded.set_ylabel(r'$\phi$ [cycles]')
        self._ax_spectrum_folded.yaxis.set_major_locator(MultipleLocator(0.5))
        self._ax_spectrum_folded.yaxis.set_minor_locator(MultipleLocator(0.1))
        self._ax_spectrum_folded.set_ylim([0.0,2.0])

        self._ax_spectrum_folded_1d.set_xlabel('Channel (PI)')
        self._ax_spectrum_folded_1d.set_ylabel('counts/s')
        self._ax_spectrum_folded_1d.set_yscale('log')

        for ax in saxes[1:]:
            ax.set_xlim([self.pulse.instrument.channels[0],
                         self.pulse.instrument.channels[-1]])
            ax.set_xscale('log')

        self._ax_spectrum_folded.tick_params(axis='x', which='both',
                                             labelbottom=False)

        self._phases = _np.linspace(0.0, 2.0, int(num_phases))
        self._energies = self.pulse.logspace_energies_hires

        #self._energies = _np.logspace(_np.log10(self.pulse.logspace_energies[0]),
        #                              _np.log10(self.pulse.logspace_energies[-1]),
        #                              int(num_energies),
        #                              base=10.0)

        self._use_fgivenx = use_fgivenx
        self._rasterized = rasterized

        yield

    @property
    def pulse(self):
        return self._pulse

    @pulse.setter
    def pulse(self, obj):
        if not isinstance(obj, Pulse):
            raise TypeError('Invalid type for Pulse object.')
        else:
            self._pulse = obj

    def _set_vminmax(self):
        """ Compute minimum and maximum for data and model colorbars. """
        self._vmin = min(_np.min(self.expected_counts/2.0),
                         _np.min(self._pulse.data.counts/2.0))
        self._vmax = max(_np.max(self.expected_counts/2.0),
                         _np.max(self._pulse.data.counts/2.0))

    @make_verbose('Adding count data to topmost panel split '
                  'over two phase cycles',
                  'Data added')
    def add_data(self):
        """ Display data in topmost panel. """

        try:
            self._vmin
        except AttributeError:
            self._set_vminmax()

        data = self._ax_data.pcolormesh(self._pulse.data.phases,
                                        self._pulse.instrument.channels,
                                        self._pulse.data.counts/2.0,
                                        cmap = cm.get_cmap('inferno'),
                                        vmin = self._vmin,
                                        vmax = self._vmax,
                                        linewidth = 0,
                                        rasterized = self._rasterized)
        data.set_edgecolor('face')

        data = self._ax_data.pcolormesh(self._pulse.data.phases + 1.0,
                                        self._pulse.instrument.channels,
                                        self._pulse.data.counts/2.0,
                                        cmap = cm.get_cmap('inferno'),
                                        vmin = self._vmin,
                                        vmax = self._vmax,
                                        linewidth = 0,
                                        rasterized = self._rasterized)
        data.set_edgecolor('face')

        self._ax_data.set_ylim([self._pulse.instrument.channels[0],
                                self._pulse.instrument.channels[-1]])
        self._ax_data.set_yscale('log')

        self._data_cb = plt.colorbar(data, cax=self._ax_data_cb,
                                     ticks=_get_default_locator(None),
                                     format=_get_default_formatter())
        self._data_cb.ax.set_frame_on(True)
        self._data_cb.ax.yaxis.set_minor_locator(AutoMinorLocator())
        self._data_cb.set_label(label=r'Counts', labelpad=15)

    def update_expectations(self, absorbed_spectrum=True, **kwargs):
        """ Update expected quantities. """

        try:
            self._model_sum
        except AttributeError:
            self._model_sum = self._pulse.expected_counts.copy()

            self._source_sums = [None]
            self._source_sums *= len(self._pulse.raw_signals_energy_intervals)

            source = None
            for i, (p, s) in enumerate(zip(self._pulse.raw_signals_energy_intervals,
                                           self._pulse.shift)):
                temp = interp(self._phases,
                              self._pulse.phases,
                              p, s)

                if source is None:
                    source = temp.copy()
                else:
                    source += temp

                self._source_sums[i] = temp

            self._folded_sums = [None] * len(self._pulse.pulse)
            folded = None
            for i, (p, s) in enumerate(zip(self._pulse.pulse,
                                           self._pulse.shift)):
                temp = interp(self._phases,
                              self._pulse.phases,
                              p, s)

                if folded is None:
                    folded = temp.copy()
                else:
                    folded += temp

                self._folded_sums[i] = temp

            self._spectrum_incident_sums = [None] * len(self._pulse.raw_signals)
            incident = None
            for i, p in enumerate(self._pulse.raw_signals):
                temp = phase_integrator(1.0,
                                        _np.array([0.0,1.0]),
                                        p,
                                        self._pulse.phases,
                                        0.0)

                #temp = energy_interpolator(1,
                #                           temp,
                #                           _np.log10(self._pulse.logspace_energies),
                #                           _np.log10(self._energies))
                temp = temp.reshape(-1)

                if incident is None:
                    incident = temp.copy()
                else:
                    incident += temp

                self._spectrum_incident_sums[i] = temp

            if absorbed_spectrum:
                self._absorbed_spectrum_incident_sums = [None] * len(self._pulse.absorbed_raw_signals)
                incident_absorbed = None
                for i, p in enumerate(self._pulse.absorbed_raw_signals):
                    temp = phase_integrator(1.0,
                                            _np.array([0.0,1.0]),
                                            p,
                                            self._pulse.phases,
                                            0.0)

                    #temp = energy_interpolator(1,
                    #                           temp,
                    #                           _np.log10(self._pulse.logspace_energies),
                    #                           _np.log10(self._energies))
                    temp = temp.reshape(-1)

                    if incident_absorbed is None:
                        incident_absorbed = temp.copy()
                    else:
                        incident_absorbed += temp

                    self._absorbed_spectrum_incident_sums[i] = temp

            spectrum_folded = None
            for i, (p, s) in enumerate(zip(self._pulse.pulse,
                                           self._pulse.shift)):
                temp = interp(self._phases,
                              self._pulse.phases,
                              p, s)

                if spectrum_folded is None:
                    spectrum_folded = temp.copy()
                else:
                    spectrum_folded += temp

            for i in range(spectrum_folded.shape[1]):
                spectrum_folded[:,i] += self._pulse.background_signal/\
                                                self._pulse.data.exposure_time
            self._spectrum_folded_sum = spectrum_folded.T

            self._spectrum_folded_1d_sums = [None] * len(self._pulse.pulse)
            for i, p in enumerate(self._pulse.pulse):
                temp = phase_integrator(1.0,
                                        _np.array([0.0,1.0]),
                                        p,
                                        self._pulse.phases,
                                        0.0)
                self._spectrum_folded_1d_sums[i] = temp.reshape(-1)

            self._counter = 1
        else:
            self._model_sum += self._pulse.expected_counts
            source = None
            for i, (p, s) in enumerate(zip(self._pulse.raw_signals_energy_intervals,
                                           self._pulse.shift)):
                temp = interp(self._phases,
                              self._pulse.phases,
                              p, s)

                if source is None:
                    source = temp.copy()
                else:
                    source += temp

                self._source_sums[i] += temp

            folded = None
            for i, (p, s) in enumerate(zip(self._pulse.pulse,
                                           self._pulse.shift)):
                temp = interp(self._phases,
                              self._pulse.phases,
                              p, s)

                if folded is None:
                    folded = temp.copy()
                else:
                    folded += temp

                self._folded_sums[i] += temp

            incident = None
            for i, p in enumerate(self._pulse.raw_signals):
                temp = phase_integrator(1.0,
                                        _np.array([0.0,1.0]),
                                        p,
                                        self._pulse.phases,
                                        0.0)

                #temp = energy_interpolator(1,
                #                           temp,
                #                           _np.log10(self._pulse.logspace_energies),
                #                           _np.log10(self._energies))
                temp = temp.reshape(-1)

                if incident is None:
                    incident = temp.copy()
                else:
                    incident += temp

                self._spectrum_incident_sums[i] += temp

            if absorbed_spectrum:
                incident_absorbed = None
                for i, p in enumerate(self._pulse.absorbed_raw_signals):
                    temp = phase_integrator(1.0,
                                            _np.array([0.0,1.0]),
                                            p,
                                            self._pulse.phases,
                                            0.0)

                    #temp = energy_interpolator(1,
                    #                           temp,
                    #                           _np.log10(self._pulse.logspace_energies),
                    #                           _np.log10(self._energies))
                    temp = temp.reshape(-1)

                    if incident_absorbed is None:
                        incident_absorbed = temp.copy()
                    else:
                        incident_absorbed += temp

                    self._absorbed_spectrum_incident_sums[i] += temp

            spectrum_folded = None
            for i, (p, s) in enumerate(zip(self._pulse.pulse,
                                           self._pulse.shift)):
                temp = interp(self._phases,
                              self._pulse.phases,
                              p, s)

                if spectrum_folded is None:
                    spectrum_folded = temp.copy()
                else:
                    spectrum_folded += temp

            for i in range(spectrum_folded.shape[1]):
                spectrum_folded[:,i] += self._pulse.background_signal/\
                                                self._pulse.data.exposure_time
            self._spectrum_folded_sum += spectrum_folded.T

            for i, p in enumerate(self._pulse.pulse):
                temp = phase_integrator(1.0,
                                        _np.array([0.0,1.0]),
                                        p,
                                        self._pulse.phases,
                                        0.0)
                self._spectrum_folded_1d_sums[i] += temp.reshape(-1)

            self._counter += 1

        if not self._use_fgivenx:
            self._add_pulse(self._ax_source_1d, source, **kwargs)
            self._add_pulse(self._ax_folded_1d, folded, **kwargs)

        try:
            incident_absorbed
        except NameError:
            incident_absorbed = None

        if self._use_fgivenx:
            return (_np.sum(source, axis=0),
                    _np.sum(folded, axis=0),
                    incident,
                    incident_absorbed)
        else:
            return None

    @property
    def model_sum(self):
        """ Get the current posterior sum of the count numbers. """
        return self._model_sum

    @model_sum.deleter
    def model_sum(self):
        del self._model_sum

    @property
    def expected_counts(self):
        """ Get the estimated posterior expectation of the count numbers. """
        return self._model_sum / self._counter

    @make_verbose('Adding posterior expectation of phase-channel count signal'
                  'to second panel',
                  'Posterior expected signal added')
    def add_expected_counts(self):
        """ Display posterior expectation of model in second panel. """

        try:
            self._vmin
        except AttributeError:
            self._set_vminmax()

        model = self._ax_model.pcolormesh(self._pulse.data.phases,
                                          self._pulse.instrument.channels,
                                          self.expected_counts/2.0,
                                          cmap = cm.get_cmap('inferno'),
                                          vmin = self._vmin,
                                          vmax = self._vmax,
                                          linewidth = 0,
                                          rasterized = self._rasterized)
        model.set_edgecolor('face')

        model = self._ax_model.pcolormesh(self._pulse.data.phases + 1.0,
                                          self._pulse.instrument.channels,
                                          self.expected_counts/2.0,
                                          cmap = cm.get_cmap('inferno'),
                                          vmin = self._vmin,
                                          vmax = self._vmax,
                                          linewidth = 0,
                                          rasterized = self._rasterized)
        model.set_edgecolor('face')

        self._ax_model.set_ylim([self._pulse.instrument.channels[0],
                                self._pulse.instrument.channels[-1]])
        self._ax_model.set_yscale('log')

        self._model_cb = plt.colorbar(model, cax=self._ax_model_cb,
                                      ticks=_get_default_locator(None),
                                      format=_get_default_formatter())
        self._model_cb.ax.set_frame_on(True)
        self._model_cb.ax.yaxis.set_minor_locator(AutoMinorLocator())

        self._model_cb.set_label(label=r'Counts', labelpad=15)

        self._add_resid()

    @make_verbose('Adding residuals between data and '
                  'posterior expected signal to third panel',
                  'Residuals added')
    def _add_resid(self):
        """ Display the residuals in the third panel. """

        self._residuals = self.expected_counts - self._pulse.data.counts
        self._residuals /= _np.sqrt(self.expected_counts)

        resid = self._ax_resid.pcolormesh(self._pulse.data.phases,
                                          self._pulse.instrument.channels, # change to channels for data only assuming contiguity
                                          self._residuals,
                                          cmap = cm.get_cmap('PuOr'),
                                          vmin=-_np.max(_np.abs(self._residuals)),
                                          vmax=_np.max(_np.abs(self._residuals)),
                                          linewidth = 0,
                                          rasterized = self._rasterized)
        resid.set_edgecolor('face')

        resid = self._ax_resid.pcolormesh(self._pulse.data.phases + 1.0,
                                          self._pulse.instrument.channels,
                                          self._residuals,
                                          cmap = cm.get_cmap('PuOr'),
                                          vmin=-_np.max(_np.abs(self._residuals)),
                                          vmax=_np.max(_np.abs(self._residuals)),
                                          linewidth = 0,
                                          rasterized = self._rasterized)
        resid.set_edgecolor('face')

        self._ax_resid.set_ylim([self._pulse.instrument.channels[0],
                                 self._pulse.instrument.channels[-1]])
        self._ax_resid.set_yscale('log')

        self._resid_cb = plt.colorbar(resid, cax = self._ax_resid_cb,
                                      ticks=AutoLocator())
        self._resid_cb.ax.set_frame_on(True)
        self._resid_cb.ax.yaxis.set_minor_locator(AutoMinorLocator())

        self._resid_cb.set_label(label=r'$(c_{ik}-d_{ik})/\sqrt{c_{ik}}$',
                                 labelpad=15)

    @property
    def expected_source(self):
        """ Get the expectations of the source components. """
        return [component/self._counter for component in self._source_sums]

    @property
    def expected_folded(self):
        """ Get the expectations of the folded source components. """
        return [component/self._counter for component in self._folded_sums]

    @property
    def expected_incident(self):
        """ Get the expectations of the incident source component spectra. """
        return [component/self._counter for component \
                                             in self._spectrum_incident_sums]

    @property
    def expected_incident_absorbed(self):
        """ Get the expectations of the incident source component spectra. """
        return [component/self._counter for component \
                                    in self._absorbed_spectrum_incident_sums]

    @property
    def expected_spectrum_folded(self):
        """ Get the estimated posterior expectation of the count numbers. """
        return self._spectrum_folded_sum / self._counter

    @property
    def expected_spectrum_folded_1d(self):
        """ Get the expectations of the incident source component spectra. """
        return [component/self._counter for component \
                                             in self._spectrum_folded_1d_sums]

    @make_verbose('Adding posterior expectation of source signals',
                  'Added posterior expected source signals')
    def add_expected_source(self,
                            expectation_line_kwargs=None,
                            plot_truth=True,
                            truth_line_kwargs=None,
                            show_components=False,
                            total_ls='-',
                            component_ls='-',
                            absorbed_spectrum=False,
                            **kwargs):
        """ Display the posterior expectation of the source signal in the
            fourth panel.

        :param bool show_components:
            If the :class:`~.Pulse.Pulse` instance has multiple components,
            display the posterior expectations of those components as a
            function of phase.

        """
        if plot_truth and truth_line_kwargs is None:
            truth_line_kwargs = \
                dict(color=('b' if self._use_fgivenx else 'darkgreen'),
                     ls=component_ls,
                     lw=1.0,
                     alpha=1.0)

        if expectation_line_kwargs is None:
            expectation_line_kwargs = dict(color='k',
                                           ls=component_ls,
                                           lw=1.0,
                                           alpha=1.0)

        if show_components and len(self._source_sums) > 1:
            if plot_truth:
                for p, s in zip(self._pulse.raw_signals_energy_intervals,
                                self._pulse.shift):
                    temp = interp(self._phases,
                                  self._pulse.phases,
                                  p, s)
                    self._add_pulse(self._ax_source_1d, temp,
                                    **truth_line_kwargs)

            for component in self.expected_source:
                self._add_pulse(self._ax_source_1d, component,
                                **expectation_line_kwargs)

        if plot_truth:
            total = None
            for p, s in zip(self._pulse.raw_signals_energy_intervals,
                            self._pulse.shift):
                temp = interp(self._phases,
                                self._pulse.phases,
                                p, s)

                if total is None:
                    total = temp.copy()
                else:
                    total += temp

            truth_line_kwargs['ls'] = total_ls
            self._add_pulse(self._ax_source_1d, total, **truth_line_kwargs)

        total = None
        for component in self.expected_source:
            if total is None:
                total = component
            else:
                total += component

        expectation_line_kwargs['ls'] = total_ls
        self._add_pulse(self._ax_source_1d, total, **expectation_line_kwargs)

        source = self._ax_source.pcolormesh(self._phases,
                                            self._pulse.energy_edges,
                                            total,
                                            cmap = cm.get_cmap('inferno'),
                                            linewidth = 0,
                                            rasterized = self._rasterized)

        source.set_edgecolor('face')
        self._ax_source.set_ylim([self._pulse.energy_edges[0],
                                  self._pulse.energy_edges[-1]])
        self._ax_source.set_yscale('log')

        self._source_cb = plt.colorbar(source, cax=self._ax_source_cb,
                                       ticks=_get_default_locator(None),
                                       format=_get_default_formatter())
        self._source_cb.ax.set_frame_on(True)
        self._source_cb.ax.yaxis.set_minor_locator(AutoMinorLocator())

        self._source_cb.set_label(label=r'photons/cm$^{2}$/s',
                                  labelpad=15)

        self._ax_source_1d.yaxis.set_major_locator(_get_default_locator(None))
        self._ax_source_1d.yaxis.set_major_formatter(_get_default_formatter())
        self._ax_source_1d.yaxis.set_minor_locator(AutoMinorLocator())

        # add folded source signal
        if plot_truth and truth_line_kwargs:
            truth_line_kwargs['ls'] = component_ls

        expectation_line_kwargs['ls'] = component_ls

        if show_components and len(self._folded_sums) > 1:
            if plot_truth:
                for p, s in zip(self._pulse.pulse, self._pulse.shift):
                    temp = interp(self._phases,
                                  self._pulse.phases,
                                  p, s)
                    self._add_pulse(self._ax_folded_1d, temp,
                                    **truth_line_kwargs)

            for component in self.expected_folded:
                self._add_pulse(self._ax_folded_1d, component,
                                **expectation_line_kwargs)

        if plot_truth:
            total = None
            for p, s in zip(self._pulse.pulse, self._pulse.shift):
                temp = interp(self._phases,
                                self._pulse.phases,
                                p, s)

                if total is None:
                    total = temp.copy()
                else:
                    total += temp

            truth_line_kwargs['ls'] = total_ls
            self._add_pulse(self._ax_folded_1d, total, **truth_line_kwargs)

        total = None
        for component in self.expected_folded:
            if total is None:
                total = component
            else:
                total += component

        expectation_line_kwargs['ls'] = total_ls
        self._add_pulse(self._ax_folded_1d, total, **expectation_line_kwargs)

        folded = self._ax_folded.pcolormesh(self._phases,
                                            self._pulse.instrument.channels,
                                            total,
                                            cmap = cm.get_cmap('inferno'),
                                            linewidth = 0,
                                            rasterized = self._rasterized)

        folded.set_edgecolor('face')
        self._ax_folded.set_ylim([self._pulse.instrument.channels[0],
                                  self._pulse.instrument.channels[-1]])
        self._ax_folded.set_yscale('log')

        self._folded_cb = plt.colorbar(folded, cax=self._ax_folded_cb,
                                       ticks=_get_default_locator(None),
                                       format=_get_default_formatter())
        self._folded_cb.ax.set_frame_on(True)
        self._folded_cb.ax.yaxis.set_minor_locator(AutoMinorLocator())

        self._folded_cb.set_label(label=r'counts/s', labelpad=15)

        self._ax_folded_1d.yaxis.set_major_locator(_get_default_locator(None))
        self._ax_folded_1d.yaxis.set_major_formatter(_get_default_formatter())

        self._ax_folded_1d.yaxis.set_minor_locator(AutoMinorLocator())

        # add spectra
        if plot_truth and truth_line_kwargs:
            truth_line_kwargs['ls'] = component_ls
            truth_line_kwargs['color'] = 'darkblue'

        expectation_line_kwargs['ls'] = component_ls
        expectation_line_kwargs['color'] = 'black'

        view_y_bottom = self._ax_spectrum_incident.yaxis.get_view_interval()[0]

        if show_components and len(self._spectrum_incident_sums) > 1:
            if plot_truth:
                for p in self._pulse.raw_signals:
                    temp = phase_integrator(1.0,
                                            _np.array([0.0,1.0]),
                                            p,
                                            self._pulse.phases,
                                            0.0)

                    #temp = energy_interpolator(1,
                    #                           temp,
                    #                           _np.log10(self._pulse.logspace_energies),
                    #                           _np.log10(self._energies))
                    temp = temp.reshape(-1)

                    self._add_incident_spectrum(self._ax_spectrum_incident,
                                                temp,
                                                **truth_line_kwargs)

            for component in self.expected_incident:
                self._add_incident_spectrum(self._ax_spectrum_incident,
                                            component,
                                            **expectation_line_kwargs)

            if absorbed_spectrum:
                for component in self.expected_incident_absorbed:
                    self._add_incident_spectrum(self._ax_spectrum_incident,
                                                component,
                                                **expectation_line_kwargs)

            for component in self.expected_spectrum_folded_1d:
                self._add_folded_spectrum(self._ax_spectrum_folded_1d,
                                          component,
                                          **expectation_line_kwargs)

        if plot_truth:
            total = None
            for p in self._pulse.raw_signals:
                temp = phase_integrator(1.0,
                                        _np.array([0.0,1.0]),
                                        p,
                                        self._pulse.phases,
                                        0.0)

                #temp = energy_interpolator(1,
                #                           temp,
                #                           _np.log10(self._pulse.logspace_energies),
                #                           _np.log10(self._energies))
                temp = temp.reshape(-1)

                if total is None:
                    total = temp.copy()
                else:
                    total += temp

            truth_line_kwargs['ls'] = total_ls
            self._add_incident_spectrum(self._ax_spectrum_incident,
                                        total,
                                        **truth_line_kwargs)

        total = None
        for component in self.expected_incident:
            if total is None:
                total = component
            else:
                total += component

        expectation_line_kwargs['ls'] = total_ls
        self._add_incident_spectrum(self._ax_spectrum_incident,
                                    total,
                                    **expectation_line_kwargs)
        if absorbed_spectrum:
            total = None
            for component in self.expected_incident_absorbed:
                if total is None:
                    total = component
                else:
                    total += component

            expectation_line_kwargs['ls'] = total_ls
            self._add_incident_spectrum(self._ax_spectrum_incident,
                                        total,
                                        **expectation_line_kwargs)

        self._ax_spectrum_incident.set_ylim(bottom=view_y_bottom)

        folded = self._ax_spectrum_folded.pcolormesh(\
                                            self._pulse.instrument.channels,
                                            self._phases,
                                            self.expected_spectrum_folded,
                                            cmap = cm.get_cmap('inferno'),
                                            linewidth = 0,
                                            rasterized = self._rasterized)

        folded.set_edgecolor('face')
        self._ax_spectrum_folded.set_xlim([self._pulse.instrument.channels[0],
                                           self._pulse.instrument.channels[-1]])

        self._spectrum_folded_cb = plt.colorbar(folded,
                                           cax=self._ax_spectrum_folded_cb,
                                           ticks=_get_default_locator(None),
                                           format=_get_default_formatter())
        self._spectrum_folded_cb.ax.set_frame_on(True)
        self._spectrum_folded_cb.ax.yaxis.set_minor_locator(AutoMinorLocator())

        self._spectrum_folded_cb.set_label(label=r'counts/s', labelpad=15)

        self._add_folded_spectrum(self._ax_spectrum_folded_1d,
          _np.sum(self.expected_counts, axis=1)/self._pulse.data.exposure_time,
          **expectation_line_kwargs)

        total = None
        for component in self.expected_spectrum_folded_1d:
            if total is None:
                total = component
            else:
                total += component

        self._add_folded_spectrum(self._ax_spectrum_folded_1d, total,
                                  **expectation_line_kwargs)

        for ax in [self._ax_spectrum_folded_1d, self._ax_spectrum_incident]:
            locmaj = LogLocator(base=10.0, numticks=100)
            ax.yaxis.set_major_locator(locmaj)

            locmin = LogLocator(base=10.0, subs=_np.arange(2,10)*0.1,
                                numticks=100)
            ax.yaxis.set_minor_locator(locmin)
            ax.yaxis.set_minor_formatter(NullFormatter())

    def _add_pulse(self, ax, signal, **kwargs):
        """ Add signal line as a function of phase. """
        if not kwargs:
            kwargs.update(dict(color='k', ls='-', lw=0.05, alpha=1.0))

        ax.plot(self._phases, _np.sum(signal, axis=0), **kwargs)

    def _add_incident_spectrum(self, ax, spectrum, **kwargs):
        """ Add incident spectrum line as a function of energy. """
        if not kwargs:
            kwargs.update(dict(color='k', ls='-', lw=0.05, alpha=1.0))

        ax.plot(self._energies, spectrum, **kwargs)

    def _add_folded_spectrum(self, ax, spectrum, **kwargs):
        """ Add folded spectrum line as a function of channel. """
        if not kwargs:
            kwargs.update(dict(color='k', ls='-', lw=0.05, alpha=1.0))

        ax.step(self._pulse.instrument.channels, spectrum,
                  where='mid', **kwargs)

    @make_verbose('Adding credible intervals on source photon flux '
                  'signal as function of phase',
                  'Credible intervals added')
    def add_source_contours(self, callback, thetas, **kwargs):
        """ Add contours to 1D source photon flux signal axes objects. """
        self._add_contours(callback, thetas, self._phases,
                           self._ax_source_1d, self._ax_source_1d_cb,
                           **kwargs)
        label = r'$\pi(\mathrm{photons/cm}^{2}\mathrm{/s};\phi)$'
        self._ax_source_1d_cb.set_ylabel(label)

    @make_verbose('Adding credible intervals on source count rate '
                  'signal as function of phase',
                  'Credible intervals added')
    def add_folded_contours(self, callback, thetas, **kwargs):
        """ Add contours to 1D source photon flux signal axes objects. """
        self._add_contours(callback, thetas, self._phases,
                           self._ax_folded_1d, self._ax_folded_1d_cb,
                           **kwargs)
        self._ax_folded_1d_cb.set_ylabel(r'$\pi(\mathrm{counts/s};\phi)$')

    @make_verbose('Adding credible intervals on source photon specific flux '
                  'spectrum',
                  'Credible intervals added')
    def add_spectrum_contours(self, callback, thetas, **kwargs):
        """ Add contours to 1D source photon flux signal axes objects. """
        self._add_contours(callback, thetas, self._energies,
                           self._ax_spectrum_incident,
                           self._ax_spectrum_incident_cb,
                           logspace_y=True,
                           cmap='RdPu_r',
                           **kwargs)
        label = r'$\pi(\mathrm{photons/keV/cm}^{2}\mathrm{/s};E)$'
        self._ax_spectrum_incident_cb.set_ylabel(label)

    @staticmethod
    def _add_contours(callback, thetas, x, ax, ax_cb,
                      logspace_y=False,
                      tqdm_kwargs={'disable': False},
                      parallel=False,
                      **kwargs):
        """ Plot contours with :mod:`fgivenx`.

        In order to increase control, we bypass direct use of some functions
        from :mod:`fgivenx.drivers`.

        Avoid fgivenx caching until it can also be implemented for other
        derived quantities which are need outside of fgivenx.

        .. todo::

            Determine whether parallelistion is straighforwardly possible
            within :mod:`fgivenx` given that the callback is complicated
            object. Currently ``OpenMP`` threading is enabled by the
            likelihood object instead.

        """

        fsamps = fgivenx.drivers.compute_samples(f=callback,
                                                 x=x,
                                                 samples=thetas,
                                                 cache=None,
                                                 parallel=parallel,
                                                 tqdm_kwargs=tqdm_kwargs)

        ymin = fsamps[~_np.isnan(fsamps)].min(axis=None)
        ymin *= kwargs.get('scale_ymin', 0.9)
        ymax = fsamps[~_np.isnan(fsamps)].max(axis=None)
        ymax *= kwargs.get('scale_ymax', 1.1)
        if not logspace_y:
            y = _np.linspace(ymin, ymax, kwargs.get('ny', 100))
        else:
            y = _np.logspace(_np.log10(ymin),
                             _np.log10(ymax),
                             kwargs.get('ny', 100),
                             base=10.0)

        pmf = fgivenx.mass.compute_pmf(fsamps=fsamps, y=y, parallel=parallel,
                                       cache=None, tqdm_kwargs=tqdm_kwargs)

        cb = fgivenx.plot.plot(x=x, y=y, z=pmf,
                             ax=ax,
                             colors=cm.get_cmap(kwargs.get('cmap', 'RdPu_r')),
                             lines=kwargs.get('lines_on', False),
                             rasterize_contours=True,
                             smooth=False)

        cb = plt.colorbar(cb, cax=ax_cb, ticks=[1,2,3])
        cb.ax.set_frame_on(True)
        cb.ax.yaxis.set_minor_locator(AutoMinorLocator())
        if kwargs.get('add_label', True):
            cb.ax.set_yticklabels([r'$1\sigma$', r'$2\sigma$',
                                            r'$3\sigma$'])

    @make_verbose('Writing to disk', 'Written')
    def savefig(self, base_dir, file_root, ext, **kwargs):
        """ Write figure to file. """

        filename = _os.path.join(base_dir,
                                 file_root + '__signalplot_residuals' + ext)
        self._fig_data.savefig(filename, **kwargs)
        filename = _os.path.join(base_dir,
                                 file_root + '__signalplot_source' + ext)
        self._fig_source.savefig(filename, **kwargs)
        filename = _os.path.join(base_dir,
                                 file_root + '__signalplot_spectrum' + ext)
        self._fig_spectrum.savefig(filename, **kwargs)

class EMCBackend(Run):
    """ Prepare `emcee`_ samples for use with :mod:`getdist`.

    .. _emcee: http://dfm.io/emcee/current/

    .. todo::
        Update to ensure compatibility, although nested sampling is
        recommended and has most features.

    :param backend: A :class:`emcee.backends.HDFBackend` instance.

    :param list names: An ordered list of ``str`` parameter names.
                       The ordering must match the parameter vector 
                       ordering defined in :obj:`backend`.

    :param int burn: If ``None``, defaults to the minumum of zero.

    :param int thin: If ``None``, defaults to the minumum of one.

    """
    def __init__(self, ID, backend, names, bounds, labels, burn, thin,
                 lines, contours, walker_discard=None,
                 settings = {}, truths = None):
        SampleContainer.__init__(self, ID, names, bounds, labels,
                                  lines, contours, truths)

        self.backend = backend

        self._mcsamples = MCSamples(samples = self._extract_samples(burn, thin, walker_discard),
                                    names = self.names,
                                    ranges = self.bounds,
                                    labels = self.labels,
                                    settings = settings)

    def _extract_samples(self, burn, thin, walker_discard):
        """ Extract samples from a :class:`~emcee.backends.HDFBackend` instance.

        :param int burn: Number of iterations to discard from the start of the
                         sampling process. Can be ``'last'``, in which case
                         only the last ensemble state is used; in this case
                         :obj:`thin` is ignored.

        :param int thin: Number of iterations to thin by to reduce
                         autocorrelation between samples such that they are
                         approximately i.i.d.

        :param walker_discard: A :class:`numpy.ndarray` of walker numbers to
                               discard from the flattened set of samples. This
                               is useful if a subset of walkers get stuck in
                               parameter space due to complex structures and
                               an inefficient proposal distribution for
                               global mapping when such structure exists.
                               If ``None`` no walkers are discarded.

        :return: An ``s x d`` :class:`numpy.ndarray` of samples, where ``n`` is
                 the number of samples and ``d`` is the number of parameters.

        """
        try:
            samples = self.backend.get_chain(discard = int(burn),
                                             thin = int(thin),
                                             flat = False)
        except TypeError:
            if burn is None: burn = 0
            if thin is None: thin = 1

            samples = self.backend.get_chain(discard = burn,
                                             thin = thin,
                                             flat = False)
        if walker_discard is not None:
            try:
                walkers = _np.arange(samples.shape[1])
                walkers = _np.delete(walkers, walker_discard)
                samples = samples[:,walkers,:]
            except IndexError:
                print('Invalid indexing instructions for discarding walkers.')
                raise

        s = list(samples.shape[1:])
        s[0] = _np.prod(samples.shape[:2])

        return samples.reshape(s)

    @property
    def backend(self):
        """ Get the :class:`emcee.backends.HDFBackend` instance. """
        return self._backend

    @backend.setter
    def backend(self, obj):
        """ Set the sample storage backend. """

        try:
            from emcee.backends import HDFBackend
        except ImportError:
            raise ImportError('Check your ``emcee`` installation.')
        else:
            if not isinstance(obj, HDFBackend):
                raise TypeError('Incompatible backend for samples.')
        self._backend = obj

    @property
    def reader(self):
        """ Get the reader. """
        return self._reader

    @reader.setter
    def reader(self, obj):
        """ Set the reader. """
        if isinstance(obj, HDFBackend):
            self._reader = obj
        else:
            raise TypeError('Invalid backend object.')



