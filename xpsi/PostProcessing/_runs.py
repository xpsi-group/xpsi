from .. import Likelihood

from ._global_imports import *

try:
    from nestcheck.ns_run_utils import combine_ns_runs
    from nestcheck.write_polychord_output import write_run_output
except ImportError:
    _warning('Cannot use nestcheck to combine runs.')

from ._metadata import Metadata
from ._run import Run
from ._backends import NestedBackend

class ParameterError(xpsiError):
    """ Thrown if inconsistent parameter names are specified for plotting. """

class Runs(Metadata):
    """ Container for nested sampling runs on a shared posterior.

    :param likelihood:
        The instance of :class:`~.Likelihood` used for sampling, or a clone of
        that instance. Defaults to ``None``. If multiple likelihood functions
        are associated with the runs due to discreteness in the model space,
        one can pass a dictionary with keys matching the IDs of the runs.

    """

    def __init__(self, runs, likelihood, *args, **kwargs):

        super(Runs, self).__init__(*args, **kwargs)

        try:
            iter(runs)
        except TypeError:
            if isinstance(runs, Run):
                self._runs = [runs]
            else:
                raise TypeError('Run objects must be instances of the ``Run`` '
                                'class.')
        else:
            if len(runs) < 1:
                raise TypeError('No objects of type ``Run`` supplied.')
            run_IDs = []
            for run in runs:
                if not isinstance(run, Run):
                    raise TypeError('Run objects must be instances of the '
                                    '``Run`` class.')
                if run.ID not in run_IDs:
                    run_IDs.append(run.ID)
                else:
                    raise AmbiguityError('Use distinct IDs for distinct sets '
                                         'of samples.')
            for run in runs:
                run.parent_ID = self.ID

            self._runs = runs

        if likelihood is not None:
            self.likelihood = likelihood

    @property
    def runs(self):
        """ Get a copy of the list of runs. """
        return list(self._runs)

    @property
    def likelihood(self):
        """ Get the likelihood instance. """

        try:
            return self._likelihood
        except AttributeError:
            return None

    @likelihood.setter
    def likelihood(self, obj):
        """ Set the likelihood object. """
        if isinstance(obj, Likelihood):
            self._likelihood = obj
        else:
             raise TypeError('The likelihood object needs to derive from '
                             'xpsi.Likelihood.Likelihood.')

        # weak check that all likelihood free parameters have associated
        # name in sample set; doesn't check that all names in self match
        # a parameter referenced by the likelihood object because of
        # derived parameters calculated and given names in the sample sets
        # that are not free parameters referenced by the likelihood object
        for name in obj.names:
            if name not in self.names:
                raise ParameterError('The likelihood object must have free '
                                     'parameters that share a name with a '
                                     'parameter in the sample set.')

    @classmethod
    def load_runs(cls, ID, run_IDs, roots, base_dirs, use_nestcheck,
                  likelihood=None, **kwargs):
        """ Construct a :class:`~.Runs` instance by loading distinct runs.

        The kwargs will be shared by nested sampling runs. The args must be
        lists that will be zipped to instantiate a set of run backends.

        """
        # if there is a transform method available, try to wrap it
        # so that error due to a mismatch in parameter order is bypassed
        if likelihood is not None:
            try:
                transform = likelihood.prior.transform
            except AttributeError: # quietly assume no transformation desired
                _transform = None
                _overwrite = False
            else:
                names = kwargs.get('names')
                def _transform(q, **kwargs):
                    p = [q[names.index(name)] for name in likelihood.names]
                    p = transform(p, **kwargs)
                    _q = [p[likelihood.names.index(name)] for name in names[:len(likelihood)]]
                    return _np.concatenate((_q, p[len(likelihood):]))
                _overwrite = kwargs.pop('overwrite_transformed', False)
        else:
            _transform = None
            _overwrite = False

        runs = []
        for root, run_ID, base_dir, check in zip(roots, run_IDs,
                                                 base_dirs,
                                                 use_nestcheck):
            runs.append(NestedBackend(root, base_dir,
                                      ID=run_ID,
                                      use_nestcheck=check,
                                      transform=_transform,
                                      overwrite_transformed=_overwrite,
                                      **kwargs))
        
        return cls(runs, likelihood, ID, **kwargs)

    def set_subset(self, IDs=None, combine=False, combine_all=False,
                   force_combine=False, only_combined=False,
                   only_principal=False, overwrite=False):
        """ Set a current list of :class:`~.Run` instances."""

        if IDs is None:
            self._subset = self._runs
        else:
            self._subset = [self[ID] for ID in IDs]

        if combine and force_combine: # create new run object
            self._combine(combine_all, overwrite)
        elif combine:
            if getattr(self, '_combined', None) is None:
                self._combine(combine_all, overwrite)
        else:
            self._combined = None

        self._only_combined = only_combined
        self._only_principal = only_principal

    def _combine(self, combine_all, overwrite):
        """ Helper method. """

        IDs = self.get_attr('ID',
                            current = False if combine_all else True,
                            nestcheck_compatible = True)

        if len(IDs) <= 1:
            self._combined = None
            return

        str_IDs = '_'.join(str(ID) for ID in IDs)

        for run in self._subset:
            if run.use_nestcheck:
                break # found nestcheck-compatible run

        # write combined run to the same directory
        base_dir = run.nestcheck_backend['output']['base_dir']
        file_root = self.ID + '_combined_IDs_%s_' % str_IDs
        _exists = _os.path.isfile(_os.path.join(base_dir, file_root+'.txt'))
        if not _exists or overwrite:
            run = combine_ns_runs(self.get_attr('nestcheck_backend',
                                    current = False if combine_all else True,
                                    nestcheck_compatible=True))

            run['output']['base_dir'] = base_dir
            run['output']['file_root'] = file_root

            try:
                #use MultiNest initial likelihood (logl_init) instead of the default PolyChord
                write_run_output(run,
                                 write_dead = True,
                                 write_stats = True,
                                 posteriors = True,
                                 stats_means_errs = True,
                                 n_simulate = 1000,
                                 logl_init = -0.179769313486231571E+309)
            except TypeError as e:
                if str(e) == "Unexpected **kwargs: {'logl_init': -1.7976931348623157e+308}":
                    raise TypeError("The used nestcheck version does not support combining "
                    "MultiNest runs in X-PSI. To use this feature, nestcheck version newer "
                    "than in this commit: "
                    "https://github.com/ejhigson/nestcheck/commit/513ef962ef7b0d66377686f9fe0a9e354dad48b3 "
                    "should be used (see installation instructions for installing it from github).")
                else:
                    raise

        kwargs = {'kde_settings': self.kde_settings,
                  'ID': 'combined',
                  'implementation': 'polychord', # match nestcheck output format
                  'names': self.names,
                  'bounds': self.bounds,
                  'labels': self.labels,
                  'truths': self.truths,
                  'precisions': self.precisions}

        self._combined = NestedBackend(file_root,
                                       base_dir = base_dir,
                                       use_nestcheck = True,
                                       **kwargs)

        self._combined.parent_ID = self.ID

    @property
    def subset(self):
        """ Get the current subset of runs for plotting. """
        return self._subset

    @property
    def subset_to_plot(self):
        """ Get the current subset of runs (+ combined) for plotting. """

        if self._combined is not None:
            if self._only_combined:
                return [self._combined]
            else:
                return [self._combined] + self._subset
        elif self._only_principal:
            return [self._subset[0]]
        else:
            return self._subset

    def __getitem__(self, ID):
        """ Get a :class:`~.Run` instance using the associated ID. """
        if isinstance(ID, _six.string_types):
            for run in self._runs:
                if ID == run.ID:
                    return run
        elif isinstance(ID, int):
            return self._runs[ID]

        # if get this far there was lookup error
        raise KeyError('No run with ID matching request.')

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
        """ Get a list of attributes of the :class:`~.Run` instances stored as
            the current subset. """

        if nestcheck_compatible:
            return [getattr(run, attribute) for run in \
                 (self._subset if current else self._runs) if run.use_nestcheck]
        else:
            return [getattr(run, attribute) for run in \
                                     (self._subset if current else self._runs)]
