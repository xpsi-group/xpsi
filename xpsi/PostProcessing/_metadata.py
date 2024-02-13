from ._global_imports import *

class Metadata(object):
    """ Base class to record basic information about sampling runs.

    :param str ID:
        For identification of the object.

    :param list names:
        An ordered list of ``str`` parameter names. The ordering must match
        the parameter vector ordering defined in sample backend objects.

    :param dict bounds:
        A dictionary of one-dimensional *hard* parameter bounds. See
        :class:`getdist.mcsamples.MCSamples`; the keys must match the
        :obj:`names` list. For the purpose of density estimation plots these
        bounds can be viewing bounds.

    :param list labels:
        An ordered list of (LaTeX compatible) ``str`` literal parameter labels.

    :param str implementation:
        Sampling software applied to generate the samples. Known options are
        ``['multinest', 'polychord', 'emcee']``.

    :param dict kde_settings:
        Settings for instantiation of :class:`getdist.mcsamples.MCSamples`.

    :param dict truths:
        Optional dictionary of parameter truths, if *known*; if *unknown*
        leave as ``None``. The keys must be names which match those in
        :obj:`names`.

    """

    def __init__(self, ID, names, bounds=None, labels=None,
                 implementation=None, kde_settings=None, truths=None):

        self.ID = ID
        self.names = names

        if bounds is not None:
            self.bounds = bounds

        if labels is None:
            self._labels = dict(zip(self.names, self.names))
        else:
            self.labels = labels

        if implementation is not None:
            self.implementation = implementation
        else:
            self._implementation = None

        if kde_settings is not None:
            self.kde_settings = kde_settings

        if truths is not None:
            self.truths = truths
        else:
            self._truths = dict(list(zip(self.names, [None] * len(self.names))))

    @property
    def ID(self):
        """ Get the identification ``str`` of the sample set. """
        return self._ID

    @ID.setter
    def ID(self, obj):
        """ Set the identification ``str`` of the sample set. """

        if not isinstance(obj, _six.string_types):
            raise TypeError('Invalid sample set ID specification. '
                            'See the docstring.')
        self._ID = obj

    @property
    def prepend_ID(self):
        """ Prepend the run ID with the parent ID."""
        try:
            return self._parent_ID + '.' + self._ID
        except AttributeError:
            return self.ID

    @property
    def parent_ID(self):
        """ The ID of the parent container with which a run is associated. """
        return self._parent_ID

    @parent_ID.setter
    def parent_ID(self, ID):
        if not isinstance(ID, _six.string_types):
            raise TypeError('Invalid parent ID specification.')
        self._parent_ID = ID

    @property
    def names(self):
        """ Get the parameter names. """

        return self._names

    @names.setter
    def names(self, obj):
        """ Set the parameter names. """

        try:
            for name in obj:
                if not isinstance(name, _six.string_types):
                    raise TypeError
        except TypeError:
            print('Invalid parameter name specification. See the docstring.')
            raise

        self._names =  obj

    def get_index(self, name):
        try:
            return self._names.index(name)
        except ValueError:
            raise KeyError('No parameter with matching name.')

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
    def kde_settings(self):
        """ Get the input :mod:`getdist` KDE settings dictionary. """
        return self._kde_settings

    @kde_settings.setter
    def kde_settings(self, obj):
        if not isinstance(obj, dict):
            raise TypeError('KDE settings for %s must be specified in a '
                            'dictionary.' % self.prepend_ID)
        self._kde_settings = obj

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

