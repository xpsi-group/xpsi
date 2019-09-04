from __future__ import print_function, division

import numpy as np
import math

import xpsi

class CustomSpacetime(xpsi.Spacetime):
    """ A custom spacetime object.

    For the NICER SWG synthetic data parameter recovery exercise, the coordinate
    rotation frequency of the star is fixed.

    """

    def __init__(self, num_params, bounds, S):
        """
        :param int num_params: The number of spacetime parameters.

        :param float S: The coordinate rotation frequency (Hz).

        """
        super(CustomSpacetime, self).__init__(num_params, bounds)

        try:
            self._S = float(S)
        except TypeError:
            raise TypeError('Coordinate spin frequency must be a ``float``.')
        else:
            self._Omega = 2.0 * math.pi * S
