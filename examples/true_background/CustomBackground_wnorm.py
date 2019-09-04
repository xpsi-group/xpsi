from __future__ import print_function, division

import numpy as np
import math

import xpsi

class CustomBackground(xpsi.Background):
    """ Currently tailored to the NICER light-curve SWG model specification.

    NICER parameter recovery from synthetic photon count data.

    The background must be set using the property method defined in the
    parent class, which will perform basic compatibility checks.

    """

    def __init__(self, num_params, bounds):
        super(CustomBackground, self).__init__(num_params, bounds)

    def __call__(self, p, energy_edges, phases):
        """ Evaluate the incident background field. """
        Gamma = p[0]
        norm = 10.0**(p[1])

        temp = np.zeros((energy_edges.shape[0] - 1, phases.shape[0] - 1))

        temp[:,0] = (energy_edges[1:]**(Gamma + 1.0) - energy_edges[:-1]**(Gamma + 1.0)) / (Gamma + 1.0)

        for i in range(temp.shape[1]):
            temp[:,i] = temp[:,0]

        self.background = norm * temp
