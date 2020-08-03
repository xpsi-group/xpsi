from __future__ import division, print_function

__all__ = ["Instrument"]

from .global_imports import *
from . import global_imports

from abc import abstractmethod
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
        Instrument channel numbers which must be equal in number to the first
        dimension of the :attr:`matrix`: the number of channels must be
        :math:`p`. The channel numbers must be monotonically increasing.
        These channels will correspond to the nominal response matrix and any
        deviation from this matrix (see above).

    .. note:: The dimensions of the response matrix need not be equal, but
              it is required that the number of input intervals be greater
              than or equal to the number of output channels -- i.e.,
              :math:`p \leq q`. If :math:`p < q` then it is implied that
              subsets of adjacent output channels are effectively grouped
              together.

    :param tuple args:
        Container of parameter instances.

    :param dict kwargs:
        If you want to prefix parameters of an instance of this instrument
        subspace with an identifier, pass it as keyword argument and it will
        find its way to the base class.

    """
    def __init__(self, matrix, energy_edges, channels, *args, **kwargs):

        self.matrix = matrix
        self.energy_edges = energy_edges
        self.channels = channels

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
            raise ResponseError('Input matrix must be a two-dimensional '
                                'ndarray awith elements that are zero '
                                'or positive.')
        else:
            try:
                for i in range(matrix.shape[0]):
                    assert matrix[i,:].any()
                for j in range(matrix.shape[1]):
                    assert matrix[:,j].any()
            except AssertionError:
                raise ResponseError('Each row and column must contain at least '
                                    'one positive number.')
            else:
                self._matrix = matrix

    def construct_matrix(self):
        """ Construct the response matrix if it is parametrised.

        If customising, do operations to calculate a matrix, and return it.
        You can access parameters (free, fixed, and derived) via the container
        access self[<name>].

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
        map to output channels defined in the data space.

        """
        return self._energy_edges

    @energy_edges.setter
    def energy_edges(self, energy_edges):
        """ Set the energy edges in keV. """

        if not isinstance(energy_edges, _np.ndarray):
            try:
                self._energy_edges = _np.array(energy_edges)
            except TypeError:
                raise EdgesError('Energy edges must be in a one-dimensional '
                                 'array, and must all be postive.')
        else:
            self._energy_edges = energy_edges

        try:
            assert self._energy_edges.ndim == 1
            assert (self._energy_edges >= 0.0).all()
            assert self._energy_edges.shape[0] == self._matrix.shape[1] + 1
        except AssertionError:
            raise EdgesError('Energy edges must be in a one-dimensional '
                             'array, and must be postive.')

    @property
    def channels(self):
        return self._channels

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
            assert self._channels.shape[0] == self._matrix.shape[0]
        except AssertionError:
            raise ChannelError('Channel numbers must be in a '
                               'one-dimensional array, and must all be '
                               'positive integers including zero.')
