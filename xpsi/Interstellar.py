from __future__ import division, print_function

__all__ = ["Interstellar"]

from .global_imports import *
from . import global_imports

from abc import abstractmethod
from .ParameterSubspace import ParameterSubspace

class Interstellar(ParameterSubspace):
    """ Base class for model interstellar X-ray processes (e.g., absorption)."""

    @abstractmethod
    def __call__(self, energies, pulse):
        """ Subclass to write a specialised specific flux modifier.

        :param energies: A :class:`numpy.ndarray` of energies in keV at which
                         the specific fluxes in :obj:`pulse` are calculated.

        :param pulse: A :class:`numpy.ndarray` of specific fluxes, with energy
                      increasing along rows and phase increasing along columns.

        .. note:: It is expected that the operations performed on a column of
                  specific fluxes need to be applied identically to all other
                  columns. The referenced object :obj:`pulse` needs to be
                  directly modified, and *not* copied.

        :return: ``None``.

        """








