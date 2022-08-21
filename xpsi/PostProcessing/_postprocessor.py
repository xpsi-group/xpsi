from __future__ import division, print_function

from ._global_imports import *

from ._backends import NestedBackend
from ._runs import Runs, ParameterError
from ._run import Run # for Sphinx docstring cross-references

class Params(object):
    """ Information about parameters shared by runs for plotting. """

    def __init__(self, names):
        self._names = names

    @property
    def names(self):
        return self._names

    @property
    def labels(self):
        return self._labels

    def __len__(self):
        return len(self._names)


class PostProcessor(object):
    """ Post-process samples for inference and posterior checking.

    Base class that functions as a container for nested sampling run objects
    representing some set of posterior distributions that share some subset of
    parameters of interest that are to be plotted. The parameters can be
    shared in the sense that their relationship is the identity relation, or
    multiple physical objects in a population context each have an instance
    of a parameter type. An example of the former is the mass of a star, for
    which there are two or more distinct posterior distributions, due, e.g.:

        * to sequential inference conditional on independent data sets, where
          each posterior is an updated of the previous in the sequence;
        * inference conditional on independent data sets, where the different
          analyses invoked same [different] prior models and one simply wants
          to overlay the posteriors e.g., to gauge the synergy of different
          experiments [and argue the priors are weakly informative/
          diffuse in the context of the likelihood functions];
        * to analysis of a single data set but the effective existence of
          discrete hyperparameter resulting in posteriors conditional on
          the discrete label over population-level prior distributions that
          has not been marginalized over (i.e., the problem of choosing an
          appropriate the hierarchical prior model);
        * to competing models, where models can differ in the likelihood
          function, prior density, or both.

    An example of the latter is two more (neutron) stars, each with a mass;
    these masses collectively encode interesting population-level information.

    Subclasses will calculate inferences, usually in the form of approximate
    integrals over the posterior distribution. They will also provide tools
    for visualisation of posterior information and derived quantities, which
    is useful for basic posterior predictive checking.

    :param iterable posteriors:
        An iterable of instances of :class:`~.Runs`.

    """
    def __init__(self, posteriors):
        self.posteriors = posteriors

        self.val_cred=[]
	self.samples={}
    @property
    def posteriors(self):
        """ Get the posteriors. """
        return self._posteriors

    @posteriors.setter
    def posteriors(self, obj):
        """ Set the posteriors attribute. """

        if isinstance(obj, (list, tuple)):
            if not obj:
                raise TypeError('No objects of type ``Runs`` supplied.')
            IDs = []
            for o in obj:
                if not isinstance(o, Runs):
                    raise TypeError('Objects must be instances of the '
                                    '``Runs`` class.')
                if o.ID not in IDs:
                    IDs.append(o.ID)
                else:
                    raise AmbiguityError('Use distinct IDs for distinct '
                                         'posteriors (run sets).')
            self._posteriors = obj
        elif isinstance(obj, Runs):
                self._posteriors = [obj]
        else:
            raise TypeError('Supply an instance or instances of the '
                            '``Runs`` class.')

    @make_verbose('Curating set of runs for posterior plotting',
                  'Run set curated')
    def set_subset(self, IDs=None,
                   combine=False, combine_all=False,
                   force_combine=True, only_combined=False,
                   overwrite=False):
        """ Set a current list of :class:`~.Runs` instances.

        Helper function to get and notify which runs will be plotted.

        :param OrderedDict IDs:
            A dictionary of lists, where keys match :class:`~.Runs` instance
            IDs, and list elements match :class:`~.Run` instance IDs. If no
            key matches a :class:`~.Runs` instance ID, it is assumed that
            all associated :class:`~.Run` instances are to be selected as the
            current subset for plotting applications.

        :param bool combine:
            Additionally combine the runs into a single run for overplotting?
            The overriding setting if there is more than one underlying
            posterior to be plotted, is to attempt to combine runs on each
            posterior if multiple such runs are available, in order to
            reduce information density. If there is a single underlying
            posterior, the user-specified value is respected.

        :param bool combine_all:
            Combine all runs in each :class:`Runs` instance or only those
            for which IDs are provided? Ignored if ``combine`` is ``False``.

        :param bool force_combine:
            Force recombination of elligible run sets, even if a
            combined run is already cached?

        :param bool only_combined:
            Only plot the combined run? Only heeded if a single posterior
            is selected for plotting, and in that case is ignored if
            ``combine`` is ``False``.

        :param bool overwrite:
            Overwrite combined-sample files on disk with the same filename?

        """

        if IDs is None:
            self._subset = self._posteriors
            IDs = {}
        else:
            if isinstance(IDs, OrderedDict):
                self._subset = [self[ID] for ID in IDs]
            else:
                raise TypeError('IDs must be supplied in a '
                                'collections.OrderedDict container.')

        if len(self._subset) > 1:
            if len(self._subset) > 5:
                print('Warning: Only the first five positional posteriors '
                      'will be plotted, with IDs'
                      + ', '.join(str(p.ID) for p in self._subset[:5]))

            for posterior in self._subset[:5]:
                posterior.set_subset(IDs.get(posterior.ID, None),
                                     combine = True,
                                     combine_all = combine_all,
                                     force_combine = True,
                                     only_combined = True,
                                     only_principal = True,
                                     overwrite = overwrite)
        else:
            posterior = self._subset[0]
            ids = IDs.get(posterior.ID, None)

            if only_combined:
                pass

            elif len(ids) > 5 and not combine:
                print('Warning: Only the first five positional runs '
                      'will be plotted individually, with IDs'
                      + ', '.join(ids[:5]))

                ids = ids[:5]

            elif len(ids) == 5 and combine:
                print('Warning: Only the first four positional runs will be '
                      'plotted individually, with IDs '
                      + ', '.join(ids[:4]))

                ids = ids[:4]

            posterior.set_subset(ids,
                                 combine = combine,
                                 combine_all = combine_all,
                                 force_combine = force_combine,
                                 only_combined = only_combined,
                                 only_principal = False,
                                 overwrite = overwrite)

    @property
    def subset(self):
        """ Get the current subset of posteriors for plotting. """
        return self._subset

    @property
    def subset_to_plot(self):
        """ Get the current subset of runs for plotting. """

        if len(self._subset) > 1:
            runs = []
            for posterior in self._subset:
                runs += posterior.subset_to_plot

            # reorder so nestcheck-compatible last
            reordered_runs = []
            temp_runs = []
            reordered_posteriors = []
            temp_posteriors = []
            for i, run in enumerate(runs):
                if run.use_nestcheck:
                    reordered_runs.append(run)
                    reordered_posteriors.append(self._subset[i])
                else:
                    temp_runs.append(run)
                    temp_posteriors.append(self._subset[i])
            reordered_runs += temp_runs
            reordered_posteriors += temp_posteriors
            self._subset = reordered_posteriors
        else:
            posterior = self._subset[0]
            runs = posterior.subset_to_plot

            # reorder so nestcheck-compatible first
            reordered_runs = []
            temp_runs = []
            for run in runs:
                if run.use_nestcheck:
                    reordered_runs.append(run)
                else:
                    temp_runs.append(run)
            reordered_runs += temp_runs

        return reordered_runs

    def __getitem__(self, ID):
        """ Get a :class:`~.Runs` instance using the associated ID. """

        def search(ID):
            if isinstance(ID, _six.string_types):
                for posterior in self._posteriors:
                    if ID == posterior.ID:
                        return posterior
            raise KeyError('No posterior with ID matching key.')

        if isinstance(ID, _six.string_types):
            return search(ID)
        elif isinstance(ID, tuple):
            if len(ID) == 2:
                return search(ID[0])[ID[1]]
            else:
                raise TypeError('Invalid run ID specification.')

    @property
    def params(self):
        """ Get the current parameter information. """
        return self._params

    def set_params(self, names):
        """ Set current parameters for plotting, which must be shared. """

        self._params = Params(names)
        for name in names:
            for posterior in self._subset:
                if name not in posterior.names:
                    self._params = None
                    raise ParameterError('No parameter name matching %s in run '
                                         'with ID %s.' % (name, posterior.ID))

        try:
            self._check_params('labels')
        except Exception:
            self._params = None
            raise

    def _check_params(self, attrs):
        """ Check consistency of parameter information across posteriors. """

        _attrs = '_' + attrs

        setattr(self._params, _attrs, [None] * len(self._params._names))
        try:
            for i, param in enumerate(self._params._names):
                for posterior in self._subset[1:]:
                    if (getattr(self._subset[0], attrs)[param] \
                            != getattr(posterior, attrs)[param]):
                        raise ValueError('Inconsistent %s for parameter'
                                         ' %s between posteriors %s and %s.' %
                                         (attrs, param,
                                          self._subset[0].ID, posterior.ID))
                getattr(self._params, _attrs)[i] = \
                        getattr(self._subset[0], attrs)[param]
        except (AttributeError, KeyError):
            print('Parameter %s not specified correctly.' % attrs)
            raise

    def get_attr(self, attribute):
        """ Get a list of attributes of the :class:`~.Runs` instances stored as
            the current subset. """

        return [getattr(run, attribute) for run in self.subset_to_plot]

    def _filter_nestcheck_compatible(self):
        """ Return only runs to plot that are compatible with nestcheck. """
        try:
            for run in self.subset_to_plot:
                if not isinstance(run, NestedBackend):
                    raise TypeError('Nested sampling backends are required.')
        except AttributeError:
            print('Nested sampling runs are required.')
            raise
        else:
            runs = [] # space to cache refs to only nestcheck-compatible runs
            nestcheck_bcknds = []
            # we only want the runs that are nestcheck-compatible
            for run in self.subset_to_plot:
                if run.use_nestcheck:
                    nestcheck_bcknds.append(run.nestcheck_backend)
                    runs.append(run)

        return nestcheck_bcknds, runs
