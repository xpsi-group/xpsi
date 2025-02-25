from .. import __version__

from ._global_imports import *

try:
    import h5py
except ImportError:
    print('Install h5py to enable signal caching.')
    raise

class _Cache(object):
    """ Cache numerical model objects computed during likelihood evaluation.

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

    def __init__(self, filename, cache_dir='./',
                 read_only=False, archive=True):

        if isinstance(filename, _six.string_types):
            if filename[-3:] != '.h5':
                self._filename = filename + '.h5'
            else:
                self._filename = filename

        self._cache_dir = cache_dir
        self._path = _os.path.join(self._cache_dir, self._filename)
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
            for key, value in data.items():
                if isinstance(value, tuple) or isinstance(value, list):
                    if key not in list(g.keys()):
                        shape = [f.attrs['n'], len(value)]
                        shape += [s for s in value[0].shape]
                        g.create_dataset(key, shape=shape, dtype='float64')

                    for j, v in enumerate(value):
                        g[key][self.i,j,...] = v
                else:
                    if key not in list(g.keys()):
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

    @make_verbose('Checking whether an existing cache can be read:',
                  'Cache state determined')
    def do_caching(self, samples, force=False):
        """ Check whether a new cache is required or whether an exising
            cache can be read without additional computation.

        :return: Boolean indicating whether to read (``False``) or write.

        """
        if force:
            self._new(samples)
            return True

        try: # try reading file and checking keys
            with self._open('r') as f:
                if 'thetas' not in list(f.keys()):
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

        if not _os.path.isdir(self._cache_dir):
            _os.mkdir(self._cache_dir)

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
        """ Check whether software version or sample set has changed. """
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
        archive_dir = _os.path.join(self._cache_dir, 'archive')

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
