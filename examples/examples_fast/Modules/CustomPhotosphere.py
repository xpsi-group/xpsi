

import numpy as np

import xpsi
from xpsi.global_imports import _2pi

class CustomPhotosphere(xpsi.Photosphere):
    """ Implement method for imaging."""

    @property
    def global_variables(self):

        return np.array([self['h__super_colatitude'],
                          self['h__phase_shift'] * _2pi,
                          self['h__super_radius'],
                          0.0, #self['p__cede_colatitude'],
                          0.0, #self['p__phase_shift'] * 2.0 * math.pi - self['p__cede_azimuth'],
                          0.0, #self['p__cede_radius'],
                          0.0, #self['s__super_colatitude'],
                          0.0, #(self['s__phase_shift'] + 0.5) * 2.0 * math.pi,
                          0.0, #self['s__super_radius'],
                          0.0, #self['s__cede_colatitude'],
                          0.0, #(self['s__phase_shift'] + 0.5) * 2.0 * math.pi - self['s__cede_azimuth'],
                          0.0, #self['s__cede_radius'],
                          self['h__super_temperature'],
                          0.0, #self['p__cede_temperature'],
                          0.0, #self['s__super_temperature'],
                          0.0]) #self['s__cede_temperature']])
