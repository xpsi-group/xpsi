__all__ = ["Interstellar"]

from xpsi.global_imports import *

from abc import abstractmethod
from xpsi.ParameterSubspace import ParameterSubspace

class Interstellar(ParameterSubspace):
    """ Base class for model interstellar X-ray processes (e.g., attenuation)."""

    @abstractmethod
    def attenuation(self, energies):
        """ Return the attenuation factor at a set of energies.

        The attenuation model can depend on fixed or free variables in the
        subspace, and generally requires loading model data from disk.

        :param ndarray[m] energies:
            An array of energies in keV at which attenuation factors are
            requested.

        :returns: An array of attenuation factors, one at each input energy.

        """
        #raise NotImplementedError('Implement the attenuation method.')

    def __call__(self, energies, signal):
        """ Attenuate a (specific) photon flux signal *in-place*.

        :param ndarray[m] energies:
            An array of energies in keV at which attenuation factors need to
            be applied to the corresponding elements of the signal array.

        :param ndarray[m[,n]] signal:
            A signal array to be attenuated, where the second dimension
            (columns) is optional and generally represents time (phase). The
            number of rows must be equal to the number of energies.

        :returns: ``None``.

        .. note::

            It is expected that the operations performed on a column of
            specific fluxes need to be applied identically to all other
            columns. The referenced object :obj:`signal` needs to be
            directly modified *in-place*, and *not* copied.

        """
        if not isinstance(signal, _np.ndarray):
            raise TypeError('Signal must be a numpy.ndarray.')

        if signal.ndim == 1:
            signal[:] *= self.attenuation(energies)
        elif signal.ndim == 2:
            for i in range(signal.shape[1]):
                signal[:,i] *= self.attenuation(energies)
        else:
            raise ValueError('Invalid number of signal array dimensions.')
