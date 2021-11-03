
__all__ = ["Background", "BackgroundError"]

from .global_imports import *
from . import global_imports

from abc import abstractmethod
from .ParameterSubspace import ParameterSubspace

class BackgroundError(xpsiError):
    """ Raised if there is an issue with the incident background . """

class Background(ParameterSubspace):
    """ Base class for astronomical photon X-ray background. """

    @abstractmethod
    def __call__(self, energy_edges, phases):
        """ Evaluate the background and store it as an attribute.

        :param energy_edges:
            The :class:`numpy.ndarray` of energy edges the photon events span.

        :param phases:
            The :class:`numpy.ndarray` of phases (in cycles) the pulse was
            evaluated at.

        .. note::

            Notice that the background can be phase-dependent, and that the
            cached background must be a two-dimensional matrix equal in shape
            to that of the pulse signal incident on the telescope.

        """
        # do some stuff and set:
        # self.background = ...

    @property
    def incident_background(self):
        """ Return the incident background specifc photon flux.

        Note that the area units need to be consistent with the units
        the effective area is expressed in within the response matrix of
        the instrument object.

        """
        return self._incident_background

    @incident_background.setter
    def incident_background(self, obj):
        """ Check and set the incident background. """
        try:
            assert isinstance(obj, _np.ndarray)
            assert obj.ndim == 2
            assert (obj >= 0.0).all()
        except AssertionError:
            raise BackgroundError('Stored background needs to be a '
                                  'two-dimensional ndarray with '
                                  'elements which are zero or positive.')
        else:
            self._incident_background = obj

    @property
    def registered_background(self):
        """ Return the registered background count rate signal. """
        return self._registered_background

    @registered_background.setter
    def registered_background(self, obj):
        """ Check and set the registered background. """
        try:
            assert isinstance(obj, _np.ndarray)
            assert obj.ndim == 2
            assert (obj >= 0.0).all()
        except AssertionError:
            raise BackgroundError('Stored background needs to be a '
                                  'two-dimensional ndarray with '
                                  'elements which are zero or positive.')
        else:
            self._registered_background = obj
