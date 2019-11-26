from __future__ import division, print_function

__all__ = ["Data"]

from .global_imports import *
from . import global_imports

from abc import ABCMeta, abstractmethod

class Data(object):
    """ A photon data container.

    The instrument associated with this data in an instance of
    :class:`~.Pulse.Pulse` must transform incident pulses into a
    structure congruent to that of the data. The attributes and methods of this
    class and any derived classes must therefore store information required
    for this operation.

    We reserve the initialiser to assign attributes which are required for the
    treating the incident specific flux pulses using the model instrument.
    The body of the initialiser must not be changed to ensure inter-module
    compatibility, but can be extended if appropriate using a call to
    ``super().__init__``. Specialist constructors can be defined in a subclass
    using the ``@classmethod`` decorator.

    The working assumption is that the sampling distribution of the data can be
    written in terms of a set of channel-by-channel *count*\ -rate pulses.

    .. note:: You need to subclass in order to tailor the handling of the
              photon event data.

    :param int first: The first instrument *output* channel which photons
                      of this data are associated with.

    :param int last: The last instrument *output* channel which photons
                     of this data are associated with.

    .. note:: For treatment of the incident pulse, it is assumed that
              that photons space a contiguous subset of output channels,
              between the :obj:`first` and :obj:`last` channel.

    """
    __metaclass__ = ABCMeta

    @abstractmethod
    def __init__(self, first, last):
        try:
            self._first = int(first)
            self._last = int(last)
        except TypeError:
            raise TypeError('The first and last channels must be integers.')

        if self._first >= self._last:
            raise ValueError('The first channel number must be lower than the '
                             'the last channel number.')

    @property
    def channel_range(self):
        """ Get a 2-tuple containing the bounding output channels. """
        return (self._first, self._last)







