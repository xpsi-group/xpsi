from __future__ import division, print_function

__all__ = ["Instrument", "ResponseError", "EdgesError"]

from .global_imports import *
from . import global_imports

from abc import abstractmethod
from .ParameterSubspace import ParameterSubspace

class ResponseError(xpsiError):
    """ Raised if there is a problem with the input response matrix. """

class EdgesError(xpsiError):
    """ Raised if there is a problem with the input energy edges. """

class Instrument(ParameterSubspace):
    """ Base class for astronomical X-ray instruments such as NICER.

    The body of the initialiser must not be changed to ensure inter-module
    compatibility, but can be extended if appropriate using a call to
    ``super().__init__``. Specialist constructors can be defined in a subclass
    using the ``@classmethod`` decorator.

    :param ndarray[p,q] matrix:
        A :math:`p \\times q` matrix which is the
        product of a redistribution matrix and effective area
        vector. The input energy intervals must increase along
        the columns of :attr:`matrix`, and the output channels
        must increase along the rows of :attr:`matrix`. The
        *units* of the elements must be that of an *effective*
        area (:math:`cm^2`).

    :param ndarray[q+1] energy_edges:
        Energy edges of the instrument energy intervals which
        must be congruent to the first dimension of the
        :attr:`matrix`: the number of edges must
        be :math:`q + 1`. The edges must be monotonically
        increasing.

    .. note:: The dimensions of the response matrix need not be equal, but
              it is required that the number of input intervals be greater
              than or equal to the number of output channels -- i.e.,
              :math:`p \leq q`. If :math:`p < q` then it is implied that subsets of
              adjacent output channels are actually grouped together.

    """
    def __init__(self, num_params, bounds, matrix, energy_edges):
        super(Instrument, self).__init__(num_params, bounds)

        self.matrix = matrix
        self.energy_edges = energy_edges

    @property
    def matrix(self):
        """ Get the response matrix.

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
            raise ResponseError('Input matrix must be a two-dimensional'
                                '``numpy.ndarray`` with elements '
                                'which are zero or positive.')
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

    @property
    def energy_edges(self):
        """ Get the energy edges of the instrument.

        A :class:`numpy.ndarray` of edges of the input energy
        intervals which map to output channels defined in the
        data space.

        """
        return self._energy_edges

    @energy_edges.setter
    def energy_edges(self):
        """ Set the energy edges. """

        try:
            assert isinstance(energy_edges, _np.ndarray)
        except AssertionError:
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

    def __call__(self, p, signal, irange, orange):
        """ Fold an incident signal.

        :param ndarray[m,n] signal:
            An :math:`m \\times n` matrix, where input energy interval
            increments along rows, and phase increases along columns.
            The number of rows, :math:`m`, must equal the number of columns of
            :attr:`matrix`: :math:`m=q`.

        :param array-like irange:
            Indexable object with two elements respectively denoting
            the indices of the first and last *input* intervals. The
            response matrix :attr:`matrix` must be indexable with
            these numbers, i.e., they must satisfy :math:`i < q`.

        :param array-like orange:
            Indexable object with two elements respectively denoting
            the indices of the first and last *output* channels. The
            response matrix :attr:`matrix` must be indexable with
            these numbers, i.e., they must satisfy :math:`i < p`.

        :return: *ndarray[p,n]* containing the folded signal.

        .. note::

            The profile most recently operated on is stored as the property
            :attr:`folded_signal`.

        """

        self.construct_matrix(p)

        self._folded_signal = _np.dot(self._matrix[orange[0]:orange[1],
                                                   irange[0]:irange[1]], signal)

        return self._folded_signal

    @property
    def folded_signal(self):
        """ Get the cached folded signal. """
        return self._folded_signal

    @abstractmethod
    def _construct_matrix(self, p):
        """ Construct the response matrix. """
