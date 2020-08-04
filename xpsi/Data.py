from __future__ import division, print_function

__all__ = ["Data"]

from .global_imports import *
from . import global_imports

class Data(object):
    """ A container for event data.

    The working assumption is that the sampling distribution of this event data
    can be written in terms of a set of channel-by-channel *count*\ -rate
    signals. The instrument associated with this data in an instance of
    :class:`~.Signal.Signal` must transform incident signals into a structure
    congruent to that of the event data. The attributes and methods of this
    class and any derived classes must therefore store information required for
    this operation.

    The initialiser to assign attributes which are required for the treating
    the incident specific flux signals using the model instrument. The body of
    the initialiser may be changed, but to ensure inter-module compatibility,
    the :meth:`channel_range` property must expose the same information. The
    initialiser can also be extended if appropriate using a call to
    ``super().__init__``. Specialist constructors can be defined in a subclass
    using the ``@classmethod`` decorator, for instance to load event data from
    disk into a compatible data structure in memory.

    .. note::

        You can subclass in order to tailor the handling of the event data, for
        instance to implement a likelihood functions for unbinned event data.

    :param int first:
        The first (loaded) instrument channel containing events (see note below).

    :param int last:
        The last (loaded) instrument channel containing events (see note below).

    .. note::

        For treatment of the incident signal, it is assumed that that events
        span a contiguous subset of channels, between and including the
        :obj:`first` and :obj:`last` channels. Moreover, the :obj:`first` and
        :obj:`last+1` channels are used to index the instrument response matrix.
        Therefore, if you load only a submatrix of the full instrument response
        matrix, these indices must be appropriate for the loaded submatrix, and
        must not be the true channel numbers (this information is instead
        loaded in the :class:`xpsi.Instrument`).

    :param float exposure_time:
        The exposure time, in seconds, to acquire this set of event data.

    :param ndarray[n,m] counts:
        A :class:`~numpy.ndarray` of count numbers. The rows of
        the array must map to a contiguous subset of instrument channels,
        with the zeroth row corresponding to the :attr:`first` channel,
        and the last row corresponding to the channel :attr:`last` channel.
        The columns must map to the phase intervals given by :obj:`phases`.

    :param ndarray[m+1] phases:
        A :class:`~numpy.ndarray` of phase interval edges, where events are
        binned into these same intervals in each instrument channel.

    """
    def __init__(self, first, last, exposure_time, counts, phases):
        try:
            self._first = int(first)
            self._last = int(last)
        except TypeError:
            raise TypeError('The first and last channels must be integers.')

        if self._first >= self._last:
            raise ValueError('The first channel number must be lower than the '
                             'the last channel number.')

        self._exposure_time = exposure_time

        if not isinstance(counts, _np.ndarray):
            raise TypeError('Counts object is not a ``numpy.ndarray``.')
        else:
            self._counts = counts

        if self._counts.shape[0] != self._last - self._first + 1:
            raise ValueError('The number of rows must be compatible '
                             'with the first and last channel numbers.')

        if not isinstance(phases, _np.ndarray):
            raise TypeError('Phases object is not a ``numpy.ndarray``.')
        else:
            self._phases = phases

    @property
    def exposure_time(self):
        """ Get the total exposure time in seconds. """
        return self._exposure_time

    @property
    def counts(self):
        """ Get the photon count data. """
        return self._counts

    @property
    def phases(self):
        """ Get the phases. """
        return self._phases

    @property
    def channel_range(self):
        """ Get a 2-tuple containing the bounding channels. """
        return (self._first, self._last + 1) # plus one for array indexing
