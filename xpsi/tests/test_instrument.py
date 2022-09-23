from __future__ import division, print_function
import numpy as np
from xpsi.Instrument import Instrument
from xpsi.Instrument import ResponseError
from xpsi.Instrument import EdgesError
from xpsi.Instrument import ChannelError
import pytest


class TestInstrument(object):

    def setup_class(cls):
        cls.matrix = np.ones((3,3))
        cls.energy_edges = np.arange(1,5,1)
        cls.channels = [0,1,2]
        cls.channel_edges = None

    def test_instrument_class_initializes(self):
        st = Instrument(self.matrix, self.energy_edges, self.channels, self.channel_edges)


class TestInstrumentMatrix(object):

    def setup_class(cls):
        cls.matrix = np.ones((3,3))
        cls.energy_edges = np.arange(1,5,1)
        cls.channels = [0,1,2]
        cls.channel_edges = None
        cls.st = Instrument(cls.matrix, cls.energy_edges, cls.channels, cls.channel_edges)

    def test_instrument_matrix_dimension(self):

        # Test if matrix is 1-D
        matrix = np.ones((3))
        with pytest.raises(ResponseError):
            st = Instrument(matrix, self.energy_edges, self.channels, self.channel_edges)

        # Test if matrix is 3-D
        matrix = np.ones((3,3,3))
        with pytest.raises(ResponseError):
            st = Instrument(matrix, self.energy_edges, self.channels, self.channel_edges)

    def test_instrument_matrix_shape(self):

        # Tests if matrix has the wrong shape (i.e., dim[0]>dim[1])
        matrix = np.ones((4,3))
        with pytest.raises(ResponseError):
            st = Instrument(matrix, self.energy_edges, self.channels, self.channel_edges)

    def test_instrument_matrix_values(self):

        # Tests if matrix has >= 0 values
        matrix = -1*np.ones((3,3))
        with pytest.raises(ResponseError):
            st = Instrument(matrix, self.energy_edges, self.channels, self.channel_edges)

    def test_instrument_matrix_content(self):

        # Tests if matrix has no positive values per row/column
        matrix = np.array([[0,0,0],[1,1,1],[1,1,1]])
        with pytest.raises(ResponseError):
            st = Instrument(matrix, self.energy_edges, self.channels, self.channel_edges)
        matrix = np.array([[0,1,1],[0,1,1],[0,1,1]])
        with pytest.raises(ResponseError):
            st = Instrument(matrix, self.energy_edges, self.channels, self.channel_edges)

    def test_matrix(self):
        self.st._matrix     ## from matrix()
        self.st.matrix      ## from construct_matrix()


class TestInstrumentEnergyEdges(object):

    def setup_class(cls):
        cls.matrix = np.ones((3,3))
        cls.energy_edges = np.arange(1,5,1)
        cls.channels = [0,1,2]
        cls.channel_edges = None
        cls.st = Instrument(cls.matrix, cls.energy_edges, cls.channels, cls.channel_edges)

    def test_instrument_energyedges_values(self):
        # Tests if energy_edges has >=0 values
        energy_edges = np.arange(-1,3,1)
        with pytest.raises(EdgesError):
            st = Instrument(self.matrix, energy_edges, self.channels, self.channel_edges)

    def test_instrument_energyedges_size(self):
        # Tests if energy_edges has the right length (nb of energy bins of matrix+1)
        energy_edges = np.ones((3))
        with pytest.raises(EdgesError):
            st = Instrument(self.matrix, energy_edges, self.channels, self.channel_edges)

    def test_instrument_energyedges_content(self):
        # Tests if energy_edges has increasing values
        energy_edges = np.flip(np.arange(1,5,1))
        with pytest.raises(EdgesError):
            st = Instrument(self.matrix, energy_edges, self.channels, self.channel_edges)

    def test_energy_edges(self):
        self.st._energy_edges
        self.st.energy_edges



class TestInstrumentChannel(object):

    def setup_class(cls):
        cls.matrix = np.ones((3,3))
        cls.energy_edges = np.arange(1,5,1)
        cls.channels = [0,1,2]
        cls.channel_edges = None
        cls.st = Instrument(cls.matrix, cls.energy_edges, cls.channels, cls.channel_edges)

    def test_instrument_channels_values(self):

        # Tests if channels has >= 0 values
        channels = np.array([-1,0,1])
        with pytest.raises(ChannelError):
            st = Instrument(self.matrix, self.energy_edges, channels, self.channel_edges)

    def test_instrument_channels_size(self):

        # Tests if channels has the right length (nb of channel bins of matrix)
        channels = np.array([0, 1, 2, 3])
        with pytest.raises(ChannelError):
            st = Instrument(self.matrix, self.energy_edges, channels, self.channel_edges)
        channels = np.array([0, 1])
        with pytest.raises(ChannelError):
            st = Instrument(self.matrix, self.energy_edges, channels, self.channel_edges)

    def test_channel(self):
        self.st._channels
        self.st.channels



class TestInstrumentChannelEdges(object):

    def setup_class(cls):
        cls.matrix = np.ones((3,3))
        cls.energy_edges = np.arange(1,5,1)
        cls.channels = [1,2,3]
        cls.channel_edges = [0,1,2,3]
        cls.st = Instrument(cls.matrix, cls.energy_edges, cls.channels, cls.channel_edges)

    def test_instrument_channelsedges_size(self):
        # Tests if channels_edges has the right length (nb of channel bins of matrix+1)
        channel_edges = np.array([0, 1, 2])
        with pytest.raises(EdgesError):
            st = Instrument(self.matrix, self.energy_edges,  self.channels, channel_edges)
        channel_edges = np.array([0, 1, 2, 3, 4])
        with pytest.raises(EdgesError):
            st = Instrument(self.matrix, self.energy_edges,  self.channels, channel_edges)

    def test_instrument_channelsedges_values(self):
        # Tests if channel_edges has >= 0 values
        channel_edges = np.array([-1, 1, 2, 3])
        with pytest.raises(EdgesError):
            st = Instrument(self.matrix, self.energy_edges, self.channels, channel_edges)

    def test_instrument_channelsedges_content(self):
        # Tests if energy_edges has increasing values
        channel_edges = np.flip(np.array([0,1,2,3]))
        with pytest.raises(EdgesError):
            st = Instrument(self.matrix, self.energy_edges, self.channels, channel_edges)

    def test_channeledges(self):
        self.st._channel_edges
        self.st.channel_edges
