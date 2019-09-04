from __future__ import print_function, division

import numpy as np
import math

import xpsi

class CustomInstrument(xpsi.Instrument):
    """ Methods and attributes specific to the NICER instrument.

    Currently tailored to the NICER light-curve SWG model specification.

    """
    def __init__(self, PI_channels, chan_edges, *args):
        super(CustomInstrument, self).__init__(*args)
        self._PI_channels = PI_channels
        self._chan_edges = chan_edges

    @property
    def channels(self):
        return self._PI_channels

    @property
    def channel_edges(self):
        """ Get the channel edges. """
        return self._chan_edges

    def _construct_matrix(self, p):
        """ Implement response matrix parameterisation. """

        return self.matrix

    def __call__(self, p, signal, *args):
        """ Overwrite. """

        matrix = self._construct_matrix(p)

        self._folded_signal = np.dot(matrix, signal)

        return self._folded_signal

    @classmethod
    def from_SWG(cls, num_params, bounds,
                 ARF, RMF, max_input, min_input=0, chan_edges=None):
        """ Constructor which converts files into :class:`numpy.ndarray`s.

        :param str ARF: Path to ARF which is compatible with
                                :func:`numpy.loadtxt`.

        :param str RMF: Path to RMF which is compatible with
                                :func:`numpy.loadtxt`.

        :param str chan_edges: Optional path to edges which is compatible with
                                :func:`numpy.loadtxt`.

        """

        if min_input != 0:
            min_input = int(min_input)

        max_input = int(max_input)

        try:
            ARF = np.loadtxt(ARF, dtype=np.double, skiprows=3)
            RMF = np.loadtxt(RMF, dtype=np.double)
            matrix = np.ascontiguousarray(RMF[min_input:max_input,20:201].T, dtype=np.double)
            if chan_edges:
                chan_edges = np.loadtxt(chan_edges, dtype=np.double, skiprows=3)[:,1:]
        except (OSError, IOError, TypeError, ValueError):
            print('A file could not be loaded.')
            raise

        edges = np.zeros(ARF[min_input:max_input,3].shape[0]+1, dtype=np.double)

        edges[0] = ARF[min_input,1]; edges[1:] = ARF[min_input:max_input,2]

        for i in range(matrix.shape[0]):
            matrix[i,:] *= ARF[min_input:max_input,3]

        PI_channels = np.arange(20, 201)

        return cls(PI_channels, chan_edges[20:202,-2],
                   num_params, bounds, matrix, edges)
