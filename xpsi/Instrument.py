from __future__ import division, print_function

__all__ = ["Instrument"]

from .global_imports import *
from . import global_imports

from . import make_verbose

from .ParameterSubspace import ParameterSubspace

class ResponseError(xpsiError):
    """ Raised if there is a problem with the input response matrix. """

class EdgesError(xpsiError):
    """ Raised if there is a problem with the input energy edges. """

class ChannelError(xpsiError):
    """ Raised if there is a problem with the input channel numbers. """

class Instrument(ParameterSubspace):
    """ Base class for astronomical X-ray instruments on-board space telescopes.

    The body of the initialiser must not be changed to ensure inter-module
    compatibility, but can be extended if appropriate using a call to
    ``super().__init__``. Specialist constructors can be defined in a subclass
    using the ``@classmethod`` decorator.

    :param ndarray[p,q] matrix:
        A :math:`p \\times q` matrix which is the product of a redistribution
        matrix and effective area vector. The input energy intervals must
        increase along the columns of :attr:`matrix`, and the output channels
        must increase along the rows of :attr:`matrix`. The *units* of the
        elements must be that of an *effective* area (:math:`cm^2`). Generally
        there will be some available calibration product, and deviations from
        this nominal response model will be parametrised. So here load some
        nominal response matrix.

    :param ndarray[q+1] energy_edges:
        Energy edges in keV of the instrument energy intervals which must be
        congruent to the first dimension of the :attr:`matrix`: the number of
        edges must be :math:`q + 1`. The edges must be monotonically increasing.
        These edges will correspond to the nominal response matrix and any
        deviation from this matrix (see above).

    :param ndarray[p] channels:
        Instrument channel numbers which must be equal in number to the number
        of rows of the :attr:`matrix`. The number of channels must therefore be
        :math:`p`. These channels will correspond to the nominal response
        matrix and any deviation from this matrix (see above). In common usage
        patterns, the channel numbers will increase monotonically with row
        number, and usually increment by one (but this is not necessary).

    .. note::

        That these channel numbers are not used to index the loaded instrument
        (sub)matrix. The :attr:`xpsi.Data.index_range` property returns
        bounding row numbers that index the loaded instrument response
        (sub)matrix in order to operate on an incident signal flux. The
        channel array contained in :attr:`xpsi.Data.channels` must be a
        contiguous (ordered) subset of the channel array loaded here.

    .. note::

        The dimensions of the response matrix need not be equal, but it is
        required that the number of input intervals be greater than or equal to
        the number of output channels -- i.e., :math:`p \leq q`. If :math:`p <
        q` then it is implied that subsets of adjacent output channels are
        effectively grouped together.

    :param ndarray[p+1] channel_edges:
        The channel (energy) edges of the instrument, in keV. The array must
        be congruent to the zeroth dimension of the :attr:`matrix`: the number
        of edges must be :math:`p + 1`. The edges must be monotonically
        increasing. These edges will correspond to the nominal response matrix
        and any deviation from this matrix (see above).

    :param tuple args:
        Container of parameter instances.

    :param dict kwargs:
        If you want to prefix parameters of an instance of this instrument
        subspace with an identifier, pass it as keyword argument and it will
        find its way to the base class.

    """
    def __init__(self, matrix, energy_edges, channels, channel_edges=None,
                 *args, **kwargs):

        self.matrix = matrix
        self.energy_edges = energy_edges
        self.channels = channels
        if channel_edges is not None:
            self.channel_edges = channel_edges

        super(Instrument, self).__init__(*args, **kwargs)

    @property
    def matrix(self):
        """ Get the reference response matrix.

        In common usage patterns there will be some fiducial or nominal
        response matrix that either defines fixed instrument operation or
        is a basis for parametrised deviations. This matrix is usually a
        calibration product distributed by an instrument calibration team.

        A matrix of dimension :math:`p \\times q`. Here :math:`p` must be the
        number of input energy intervals, and :math:`q \geq p` the number of
        output channels.

        .. note::

            The attribute :attr:`matrix` must be assigned, and it must be
            a :class:`numpy.ndarray` for use with :func:`numpy.dot` (even
            if the matrix is sparse to some degree).

        """
        return self._matrix

    @matrix.setter
    def matrix(self, matrix):
        """ Set the matrix. """
        try:
            assert isinstance(matrix, _np.ndarray)
            assert matrix.ndim == 2
            assert matrix.shape[0] <= matrix.shape[1]
            assert (matrix >= 0.0).all()
        except AssertionError:
            raise ResponseError('Input matrix must be a two-dimensional ndarray, with number of channel intervals'
                                '(dimension 0) smaller than the number of energy intervals (dimension 1),'
                                ' and with all matrix elements that are zero or positive.')
        try:
            for i in range(matrix.shape[0]):
                assert matrix[i,:].any()
            for j in range(matrix.shape[1]):
                assert matrix[:,j].any()
        except AssertionError:
            raise ResponseError('Each row and column of the matrix must contain at least '
                                'one positive number.')
        self._matrix = matrix

    def construct_matrix(self):
        """ Construct the response matrix if it is parameterised.

        If customising, do operations to calculate a matrix, and return it.
        You can access parameters (free, fixed, and derived) via the container
        access ``self[<name>]``.

        If the instrument operation is fixed, you might not need to subclass,
        because the default behaviour is to return the nominal response you
        loaded. If for some reason the matrix you loaded is to be modified in
        some fixed manner, possibly as a function of some custom fixed
        parameters that you defined, you would also have to subclass and
        provide the correct implementation of this method.

        """
        return self.matrix

    def __call__(self, signal, irange, orange):
        """ Register an incident signal.

        :param ndarray[m,n] signal:
            An :math:`m \\times n` matrix, where input energy interval
            increments along rows, and phase increases along columns. The
            number of rows, :math:`m`, must equal the number of columns of
            :attr:`matrix`: :math:`m=q`.

        :param array-like irange:
            Indexable object with two elements respectively denoting the
            indices of the first and last *input* intervals. The response
            matrix :attr:`matrix` must be indexable with these numbers, i.e.,
            they must satisfy :math:`indx < q`.

        :param array-like orange:
            Indexable object with two elements respectively denoting the
            indices of the first and last *output* channels. The response
            matrix :attr:`matrix` must be indexable with these numbers, i.e.,
            they must satisfy :math:`indx < p`.

        :return: *ndarray[p,n]* containing the registered signal.

        .. note::

            The product of the most recent operation is stored as the property
            :attr:`cached_signal`.

        """
        matrix = self.construct_matrix()

        self._cached_signal = _np.dot(matrix[orange[0]:orange[1],
                                             irange[0]:irange[1]], signal)

        return self._cached_signal

    @property
    def cached_signal(self):
        """ Get the cached registered signal. """
        return self._cached_signal

    @property
    def energy_edges(self):
        """ Get the energy edges of the instrument, in keV.

        A :class:`numpy.ndarray` of edges of the input energy intervals which
        map to channels defined in the data space.

        """
        return self._energy_edges

    @energy_edges.setter
    def energy_edges(self, energy_edges):
        """ Set the energy edges in keV. """

        if not isinstance(energy_edges, _np.ndarray):
            try:
                energy_edges = _np.array(energy_edges)
            except TypeError:
                raise EdgesError('Energy edges must be in a one-dimensional array of positive increasing values.')

        try:
            assert energy_edges.ndim == 1
            assert (energy_edges >= 0.0).all()
            assert energy_edges.shape[0] == self._matrix.shape[1] + 1
            assert not (energy_edges[1:] <= energy_edges[:-1]).any()
        except AssertionError:
            raise EdgesError('Energy edges must be in a one-dimensional array of positive increasing values, with a '
                             'length equal to the length of the response matrix + 1.')

        self._energy_edges = energy_edges

    @property
    def channel_edges(self):
        """ Get the channel (energy) edges of the instrument, in keV.

        A :class:`numpy.ndarray` of edges of the registered energy intervals
        labelled as channels defined in the data space. This is relevant when
        there is a detector-by-detector gain scale applied to event data (such
        as for NICER instrument calibration products), meaning that the
        redistribution matrix is effectively shared by detectors and the
        channels across detectors can share an energy scale definition.

        An incident photon of given energy then has a registered-energy
        distribution that generally peaks in the vicinity of the true photon
        energy. The resdistribution matrix will have some energy resolution
        (along with other features such as shelves). With thanks to Paul S.
        Ray for explaining the choice to calibrate in this manner.

        .. note::

            If you made a channel cut that results in a non-contiguous subset
            of channels, you will need to overwrite the setter method because
            the checks will fail.

        """
        return self._channel_edges

    @channel_edges.setter
    def channel_edges(self, channel_edges):
        """ Set the channel (energy) edges in keV. """

        if not isinstance(channel_edges, _np.ndarray):
            try:
                channel_edges = _np.array(channel_edges)
            except TypeError:
                raise EdgesError('Channel edges must be in a one-dimensional array of positive increasing values.')

        try:
            assert channel_edges.ndim == 1
            assert (channel_edges >= 0.0).all()
            assert channel_edges.shape[0] == self._matrix.shape[0] + 1
            assert not (channel_edges[1:] <= channel_edges[:-1]).any()
        except AssertionError:
            raise EdgesError('Channel edges must be in a one-dimensional array of positive increasing values, with a '
                             'length equal to the length of the response matrix + 1.')

        self._channel_edges = channel_edges

    @property
    def channels(self):
        """ Get the array of channels corresponding to rows of the matrix.

        The matrix being the loaded instrument response (sub)matrix.

        """
        return self._channels

    @channels.setter
    @make_verbose('Setting channels for loaded instrument response (sub)matrix',
                  'Channels set')
    def channels(self, channel_array):
        if not isinstance(channel_array, _np.ndarray):
            try:
                channel_array = _np.array(channel_array)
            except TypeError:
                raise ChannelError('Channel numbers must be in a one-dimensional array of positive integers'
                                   '(including zero).')

        try:
            assert channel_array.ndim == 1
            assert (channel_array >= 0).all()
            assert channel_array.shape[0] == self._matrix.shape[0]
        except AssertionError:
            raise ChannelError('Channel numbers must be in a one-dimensional array of positive integers'
                               '(including zero).')

        if (channel_array[1:] - channel_array[:-1] != 1).any():
            yield ('Warning: Channel numbers do not uniformly increment by one.\n'
                   '         Please check for correctness.')

        self._channels = channel_array

        yield
