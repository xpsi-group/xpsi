

import numpy as np
import math

import xpsi
from xpsi.global_imports import _2pi

class CustomPhotosphere(xpsi.Photosphere):
    """ Implement method for imaging."""

    @property
    def global_variables(self):

        return np.array([self['h__super_colatitude'],
                          self['h__phase_shift'] * _2pi,
                          self['h__super_radius'],
                          self['h__super_temperature']])
