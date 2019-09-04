from __future__ import print_function, division

import numpy as np
import math

import xpsi

from xpsi.likelihoods.Poisson_loglike import eval_loglike
from xpsi.tools import phase_interpolator
from xpsi.tools.phase_integrator import phase_integrator
from xpsi.tools.synthesise import synthesise

class CustomPulse(xpsi.Pulse):
    """ A custom calculation of the logarithm of the likelihood.

    We extend the :class:`xpsi.Pulse.Pulse` class to make it callable.

    We overwrite the body of the __call__ method. The docstring for the
    abstract method is copied.

    """

    def __init__(self, **kwargs):
        """ Initialise. """

        super(CustomPulse, self).__init__(**kwargs)

        try:
            self._data.counts
        except AttributeError:
            print('No data... can synthesise data but cannot evaluate a '
                  'likelihood function.')

    def __call__(self, p, *args, **kwargs):
        """

        Parameter vector:

        * p[0] = phase shift primary (alias for initial azimuth/phase of photosphere)
        * p[1] = phase shift secondary

        """
        self.shift = np.array(p)

        self.loglikelihood, self.expected_counts = \
                eval_loglike(self._data.exposure_time,
                                      self._data.phases,
                                      self._data.counts,
                                      self._pulse,
                                      self._phases,
                                      self._shift,
                                      self._background.folded_background)

        self._background_signal = self._background.folded_background[:,0] * self._data.exposure_time

    __call__.__doc__ = xpsi.Pulse.__call__.__doc__ + __call__.__doc__

    def synthesise(self, p,
                   require_source_counts,
                   require_background_counts,
                   directory='./', **kwargs):
        """ Synthesise data set.

        * p[0] = phase shift primary (alias for initial azimuth/phase of photosphere)
        * p[1] = phase shift secondary

        """
        self.shift = np.array(p)

        self._expected_counts, synthetic = synthesise(self._data.phases,
                                            require_source_counts,
                                            require_background_counts,
                                            self._pulse,
                                            self._phases,
                                            self._background.folded_background,
                                            self._shift)

        np.savetxt(directory + '/synthetic_realisation.dat',
                   synthetic,
                   fmt = '%u')

        self._write(self.expected_counts,
                    filename = directory + '/synthetic_expected_hreadable.dat',
                    fmt = '%.8e')

        self._write(synthetic,
                    filename = directory + '/synthetic_realisation_hreadable.dat',
                    fmt = '%u')


    def _write(self, counts, filename, fmt):
        """ Write to file. """

        rows = len(self._data.phases) - 1
        rows *= self._data.channel_range[1] - self._data.channel_range[0]

        phases = self._data.phases[:-1]
        array = np.zeros((rows, 3))

        for i in range(counts.shape[0]):
            for j in range(counts.shape[1]):
                array[i*len(phases) + j,:] = i+20, phases[j], counts[i,j]

        np.savetxt(filename, array, fmt=['%u', '%.6f'] + [fmt])


