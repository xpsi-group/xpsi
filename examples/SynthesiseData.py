from __future__ import print_function, division

import numpy as np
import math

import xpsi

class CustomData(xpsi.Data):
    """ Custom data container.

    Currently tailored to the NICER light-curve SWG model specification.

    """
    def __init__(self, first, last, phases):
        """
        :param phase_edges: A :class:`numpy.ndarray` of phase interval edges.

        """
        # Execute parent initialisation code
        super(CustomData, self).__init__(first, last)

        self._phases = phases

    @property
    def phases(self):
        """ Get the phase edges. """
        return self._phases

