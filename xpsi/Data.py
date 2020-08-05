from __future__ import division, print_function

__all__ = ["Data"]

from .global_imports import *
from . import global_imports
from . import make_verbose

from .Instrument import ChannelError

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
    the :meth:`index_range` property must expose the same information. The
    initialiser can also be extended if appropriate using a call to
    ``super().__init__``. Specialist constructors can be defined in a subclass
    using the ``@classmethod`` decorator, for instance to load event data from
    disk into a compatible data structure in memory.

    .. note::

        You can subclass in order to tailor the handling of the event data, for
        instance to implement a likelihood functions for unbinned event data.

    :param ndarray[n,m] counts:
        A :class:`~numpy.ndarray` of count numbers. The columns must map to
        the phase intervals given by :obj:`phases`. The rows of the array map
        to some subset of instrument channels.

    :param ndarray[n] channels:
        Instrument channel numbers which must be equal in number to the first
        dimension of the :attr:`matrix`: the number of channels must be
        :math:`p`. These channels will correspond to the nominal response
        matrix and any deviation from this matrix (see above). In common usage
        patterns, the channel numbers will increase monotonically with row
        number, and usually increment by one (but this is not necessary). It is
        advisable that these numbers are the actual instrument channel numbers
        so that plots using these labels are clear.

    :param ndarray[m+1] phases:
        A :class:`~numpy.ndarray` of phase interval edges, where events are
        binned into these same intervals in each instrument channel.

    :param int first:
        The index of the first row of the loaded response matrix containing
        events (see note below).

    :param int last:
        The index of the last row of the loaded response matrix containing
        events (see note below).

    .. note::

        The :obj:`counts` matrix rows  *must* span a contiguous subset of the
        rows of the loaded response matrix, but in general can span an
        arbitrary subset and order of instrument channels. Note that the
        :obj:`first` and :obj:`last+1` numbers are used to index the loaded
        instrument response matrix.  Therefore, if you load only a submatrix of
        the full instrument response matrix, these indices must be appropriate
        for the loaded submatrix, and must not be the true channel numbers
        (this information is instead loaded in the :class:`xpsi.Instrument`).
        Of course, in all sensible usage patterns the order of the instrument
        channels, when mapped to matrix rows, will be such that channel number
        increases with matrix row number monotonically because, then the
        nominal event energy increases monotonically with row number, which is
        important for visualisation of data and model (because spatial order
        matches energy order and correlations in space can be discerned easily).
        However, there could in principle be channel cuts that mean an increment
        of more than one channel between matrix adjacent rows, and the
        response matrix needs to be manipulated before or during a custom
        loading phase such that its rows match the channel numbers assigned to
        the :obj:`counts` matrix rows.

    :param float exposure_time:
        The exposure time, in seconds, to acquire this set of event data.

    """
    def __init__(self, counts, channels, phases, first, last, exposure_time):

        if not isinstance(counts, _np.ndarray):
            raise TypeError('Counts object is not a ``numpy.ndarray``.')
        else:
            self._counts = counts

        self.channels = channels

        if not isinstance(phases, _np.ndarray):
            raise TypeError('Phases object is not a ``numpy.ndarray``.')
        else:
            self._phases = phases

        try:
            self._first = int(first)
            self._last = int(last)
        except TypeError:
            raise TypeError('The first and last channels must be integers.')

        if self._first >= self._last:
            raise ValueError('The first channel number must be lower than the '
                             'the last channel number.')

        if self._counts.shape[0] != self._last - self._first + 1:
            raise ValueError('The number of rows must be compatible '
                             'with the first and last channel numbers.')

        self._exposure_time = exposure_time

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
    def index_range(self):
        """ Get a 2-tuple of the bounding response-matrix row indices. """
        return (self._first, self._last + 1) # plus one for array indexing

    @property
    def channel_range(self):
        """ Deprecated property name. To be removed for v1.0. """
        return self.index_range

    @property
    def channels(self):
        return self._channels

    @make_verbose('Setting channels for event data',
                  'Channels set')
    @channels.setter
    def channels(self, array):
        if not isinstance(array, _np.ndarray):
            try:
                self._channels = _np.array(array)
            except TypeError:
                raise ChannelError('Channel numbers must be in a '
                                   'one-dimensional array, and must all be '
                                   'positive integers including zero.')
        else:
            self._channels = array

        try:
            assert self._channels.ndim == 1
            assert (self._channels >= 0).all()
            assert self._channels.shape[0] == self._counts.shape[0]
        except AssertionError:
            raise ChannelError('Channel numbers must be in a '
                               'one-dimensional array, and must all be '
                               'positive integers including zero.')

        if not (self._channels[1:] - self._channels[:-1] != 1).any():
            yield ('Warning: Channel numbers do not uniformly increment by one.'
                   '\n         Please check for correctness.')

        yield

    @make_verbose('Loading event list and phase binning',
                  'Events loaded and binned')
    @classmethod
    def phase_bin__event_list(cls, path, channels, phases,
                              phase_shift=0.0,
                              phase_column=1,
                              skiprows=1, *args, **kwargs):
        """ Load a phase-folded event list and bin the events in phase.

        :param str path:
            Path to event list file containing two columns, where the first
            column contains phases on the unit interval, and the second
            column contains the channel number.

        :param list channels:
            An (ordered) subset of instrument channels. It is advisable that
            these channels are a contiguous subset of instrument channels, but
            this not a strict requirement if you are comfortable with the
            handling the instrument response matrix and count number matrix to
            match in row-to-channel definitions.

        :param list phases:
            An ordered sequence of phase-interval edges on the unit interval.
            The first and last elements will almost always be zero and unity
            respectively.

        :param float phase_shift:
            A phase-shift in cycles to be applied when binning the events in
            phase.

        :param int phase_column:
            The column in the loaded file containing event phases. Either zero
            or one.

        :param int skiprows:
            The number of top rows to skip when loading the events from file.
            The top row of couple of rows will typically be reserved for
            column headers.

        """
        events = _np.loadtxt(path, skiprows=skiprows)

        channels = list(channels)
        channel_col = int(not phase_column)

        yield 'Total number of events: %i.' % events.shape[0]

        data = _np.zeros((len(channels), len(phases)-1), dtype=_np.int)

        for i in range(events.shape[0]):
            if events[i,channel_col] in channels:
                _temp = events[i,phase_column] + phase_shift
                _temp -= _np.floor(_temp)

                for j in range(phases.shape[0]-1):
                    if phases[j] <= _temp <= phases[j+1]:
                        data[channels.index(int(events[i,channel_col])),j] += 1
                        break

        yield 'Number of events constituting data set: %i.' % _np.sum(data)

        yield cls(data, channels, _np.array(phases), *args, **kwargs)
