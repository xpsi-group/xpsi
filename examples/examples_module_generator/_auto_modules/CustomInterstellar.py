""" Interstellar module for X-PSI CST+PDT modelling of NICER PSR J0030+0451 event data. """

import numpy as np
import math

import xpsi
from xpsi import Parameter

from scipy.interpolate import Akima1DInterpolator

class CustomInterstellar(xpsi.Interstellar):
    """ Apply interstellar attenuation model tbnew. """

    def __init__(self, energies, attenuation, bounds, values = None):

        if values is None: values = {}

        assert len(energies) == len(attenuation), 'Array length mismatch.'

        self._lkp_energies = energies # for lookup
        self._lkp_attenuation = attenuation # for lookup

        N_H = Parameter('neutral_hydrogen_column_density',
                        strict_bounds = (0.0, 50.0),
                        bounds = bounds.get('neutral_hydrogen_column_density', None),
                        doc = 'Neutral hydrogen column density in units of the fiducial column density',
                        symbol = r'$N_{\rm H}$',
                        value = values.get('neutral_hydrogen_column_density', None),
                        permit_prepend = False)

        self._interpolator = Akima1DInterpolator(self._lkp_energies,
                                                 self._lkp_attenuation)
        self._interpolator.extrapolate = True

        super(CustomInterstellar, self).__init__(N_H)

    def attenuation(self, energies):
        """ Interpolate the attenuation coefficients.

        Useful for post-processing.

        """
        return self._interpolate(energies)**(self['neutral_hydrogen_column_density'])

    def _interpolate(self, energies):
        """ Helper. """
        _att = self._interpolator(energies)
        _att[_att < 0.0] = 0.0
        return _att

    @classmethod
    def load(cls, path,
             energy_column=0,
             attenuation_column=1,
             **kwargs):
        """ Load attenuation file. """

        # check the loading assumptions and comment out the exception throw if they are true
        #raise NotImplementedError('Implement the class method to load the interstellar attenuation table.')

        # template

        temp = np.loadtxt(path, dtype=np.double)

        energies = temp[:,energy_column]

        attenuation = temp[:,attenuation_column]

        return cls(energies, attenuation, **kwargs)