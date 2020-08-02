from __future__ import division, print_function

from ._global_imports import *

from ._metadata import Metadata

class Run(Metadata):
    """ Base class for a sample container. """

    def __init__(self, filepath, **kwargs):

        self.ID = kwargs.pop('ID')

        if kwargs:
            super(Run, self).__init__(self.ID, **kwargs)

        self.samples = filepath

    @property
    def samples(self):
        """ Get a copy of the samples array. """

        return self._samples.copy()

    @samples.setter
    def samples(self, filepath):
        if _os.path.isfile(filepath):
            # we should be able to load samples straightforwardly into
            # a NumPy array for the most basic type of sample lookup
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

        :param dict contours:
            :mod:`getdist`-compatible dictionary of parameters specifying the
            properties of the contours representing credible regions of the
            two-dimensional marginal distributions of parameters.

        """
        return self._contours

    @contours.setter
    def contours(self, obj):

        if not isinstance(obj, dict):
            raise TypeError('Invalid contour argument specification. '
                            'See the docstring.')
        self._contours = obj
