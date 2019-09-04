from __future__ import print_function, division

import numpy as np
import math

import xpsi

class CustomData(xpsi.Data):
    """ Custom data container.

    Currently tailored to the NICER light-curve SWG model specification.

    """
    def __init__(self, first, last, counts, phases, exposure_time):
        """
        :param counts: A :class:`numpy.ndarray` of count numbers. The rows of
                       the array must map to a contiguous subset of instrument
                       output channels, with the zeroth row corresponding to
                       the :attr:`first` channel, and the last row
                       corresponding to the channel :attr:`last` minus one.
                       The columns must map to the phases given by
                       :obj:`phases`.

        :param phases: A :class:`numpy.ndarray` of phase *edges* of intervals
                       in which the *synthetic* photons arrive.

        :param exposure_time: The total exposure time in seconds.

        """
        # Execute parent initialisation code
        super(CustomData, self).__init__(first, last)

        try:
            assert isinstance(counts, np.ndarray)
        except AssertionError:
            raise TypeError('Counts object is not a ``numpy.ndarray``.')
        else:
            self._counts = counts

        try:
            assert self._counts.shape[0] == self._last - self._first
        except AssertionError:
            raise AssertionError('The number of rows must be compatible '
                                 'with the first and last output channel '
                                 'numbers.')

        try:
            assert isinstance(phases, np.ndarray)
        except AssertionError:
            raise TypeError('Phases object is not a ``numpy.ndarray``.')
        else:
            self._phases = phases

        self._exposure_time = exposure_time

    @property
    def exposure_time(self):
        """ Get the total exposure time in seconds. """
        return self._exposure_time

    @property
    def counts(self):
        """ Get the photon count data. """
        return self._counts

    @property
    def phases(self):
        """ Get the phases. """
        return self._phases

    @classmethod
    def from_SWG(cls, path, *args):
        """ Constructor which loads photon data from a .txt file.

        :param str path: Path to .txt file which is converted into a
                         two-dimensional :class:`numpy.ndarray`.

        """
        try:
            data = np.loadtxt(path, dtype=np.double)
        except (OSError, IOError, TypeError, ValueError):
            print('Data file could not be loaded.')
            raise

        first = 0; last = 181

        phases = np.linspace(0.0, 1.0, 33)

        return cls(first, last, data, phases, *args)

