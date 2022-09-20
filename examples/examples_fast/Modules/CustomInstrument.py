from __future__ import print_function, division

import numpy as np
import math

import xpsi

from xpsi import Parameter, make_verbose

# Instrument Class

class CustomInstrument(xpsi.Instrument):

    """ Fake telescope response. """



    def __call__(self, signal, *args):
        """ Overwrite base just to show it is possible.

        We loaded only a submatrix of the total instrument response
        matrix into memory, so here we can simplify the method in the
        base class.

        """
        matrix = self.construct_matrix()

        self._folded_signal = np.dot(matrix, signal)

        return self._folded_signal

    @classmethod
    def from_response_files(cls, ARF, RMF, max_input, min_input=0,channel=[1,1500],
                            channel_edges=None):
        """ Constructor which converts response files into :class:`numpy.ndarray`s.
        :param str ARF: Path to ARF which is compatible with
                                :func:`numpy.loadtxt`.
        :param str RMF: Path to RMF which is compatible with
                                :func:`numpy.loadtxt`.
        :param str channel_edges: Optional path to edges which is compatible with
                                  :func:`numpy.loadtxt`.
        """

        if min_input != 0:
            min_input = int(min_input)

        max_input = int(max_input)

        matrix = np.ascontiguousarray(RMF[min_input:max_input,channel[0]:channel[1]].T, dtype=np.double)

        edges = np.zeros(ARF[min_input:max_input,2].shape[0]+1, dtype=np.double)



        edges[0] = ARF[min_input,0]; edges[1:] = ARF[min_input:max_input,1]

        for i in range(matrix.shape[0]):
            matrix[i,:] *= ARF[min_input:max_input,2]


        channels = np.arange(channel[0],channel[1])


        return cls(matrix, edges, channels, channel_edges[channel[0]:channel[1]+1,1])
