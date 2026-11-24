
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


def _write_rmf(path, tlmin=1, write_tlmin=True, extname='MATRIX', detchans=5):
    """ Write a minimal OGIP-like diagonal RMF with ``detchans`` channels and
    energies, channel numbers starting at ``tlmin``. """
    from astropy.io import fits

    n = detchans
    energ_lo = np.arange(1.0, 1.0 + n)
    matrix = np.zeros((n, 2), dtype=np.float32)
    matrix[:, 0] = 1.0
    matrix_cols = [
        fits.Column(name='ENERG_LO', array=energ_lo, format='E'),
        fits.Column(name='ENERG_HI', array=energ_lo + 1.0, format='E'),
        fits.Column(name='N_GRP', array=np.ones(n, dtype=np.int16), format='I'),
        fits.Column(name='F_CHAN', array=np.arange(tlmin, tlmin + n, dtype=np.int32), format='J'),
        fits.Column(name='N_CHAN', array=np.ones(n, dtype=np.int32), format='J'),
        fits.Column(name='MATRIX', array=matrix, format='2E'),
    ]
    matrix_hdu = fits.BinTableHDU.from_columns(matrix_cols, name=extname)
    matrix_hdu.header['INSTRUME'] = 'TEST'
    matrix_hdu.header['DETCHANS'] = n
    if write_tlmin:
        matrix_hdu.header['TLMIN4'] = tlmin
        matrix_hdu.header['TLMAX4'] = tlmin + n - 1

    ebounds_cols = [
        fits.Column(name='CHANNEL', array=np.arange(tlmin, tlmin + n, dtype=np.int32), format='J'),
        fits.Column(name='E_MIN', array=energ_lo, format='E'),
        fits.Column(name='E_MAX', array=energ_lo + 1.0, format='E'),
    ]
    ebounds_hdu = fits.BinTableHDU.from_columns(ebounds_cols, name='EBOUNDS')
    fits.HDUList([fits.PrimaryHDU(), matrix_hdu, ebounds_hdu]).writeto(path, overwrite=True)


@pytest.mark.parametrize('tlmin, write_tlmin, extname', [
    (0, True, 'MATRIX'),
    (1, True, 'MATRIX'),
    (1, False, 'MATRIX'),
    (1, True, 'SPECRESP MATRIX'),
])


@pytest.mark.parametrize('twod_edges', [False, True])
def test_from_ogip_fits_channel_range(tmp_path, tlmin, write_tlmin, extname, twod_edges):
    'All DETCHANS channels are loaded, with or without TLMIN/TLMAX and for RMF or RSP extension names.'

    path = tmp_path / 'response.rmf'
    _write_rmf(str(path), tlmin=tlmin, write_tlmin=write_tlmin, extname=extname)

    instr = Instrument.from_ogip_fits(str(path), twod_edges=twod_edges)

    assert np.array_equal(instr.channels, np.arange(tlmin, tlmin + 5))
    assert instr.TLMIN == tlmin
    assert np.allclose(instr.matrix, np.eye(5))
    if twod_edges:
        assert instr.channel_edges.shape == (2, 5)
        assert np.allclose(instr.channel_edges[0], np.arange(1.0, 6.0))
    else:
        assert np.allclose(instr.channel_edges, np.arange(1.0, 7.0))
        
        
class TestInstrumentParameterArguments(object):
    """ Parameters can be passed positionally after ``channel_edges``.

    Custom instruments (and modules made by the module generator) typically
    end with ``cls(matrix, energy_edges, channels, channel_edges, alpha, **kwargs)``.
    Regression test: in X-PSI v3.4.0, ``alpha`` was taken as ``twod_edges``.
    """

    class ScaledInstrument(Instrument):
        """ Minimal custom instrument with an energy-independent scaling factor. """

        def construct_matrix(self):
            return self['alpha'] * self.matrix

        @classmethod
        def from_arrays(cls, matrix, energy_edges, channels, channel_edges, **kwargs):
            from xpsi.Parameter import Parameter
            alpha = Parameter('alpha',
                              strict_bounds=(0.1, 1.9),
                              bounds=(0.5, 1.5),
                              doc='Energy-independent scaling factor',
                              value=1.2)
            return cls(matrix, energy_edges, channels, channel_edges, alpha, **kwargs)

    def setup_class(cls):
        cls.matrix = np.ones((3,3))
        cls.energy_edges = np.arange(1,5,1)
        cls.channels = [0,1,2]
        cls.channel_edges = np.arange(1,5,1)

    def test_positional_parameter(self):
        st = self.ScaledInstrument.from_arrays(self.matrix, self.energy_edges,
                                               self.channels, self.channel_edges,
                                               prefix='XTI')
        assert st._twod_edges is False
        assert st['alpha'] == 1.2
        assert np.allclose(st.construct_matrix(), 1.2 * self.matrix)

    def test_twod_edges_keyword(self):
        energy_edges = np.array([self.energy_edges[:-1], self.energy_edges[1:]])
        channel_edges = np.array([self.channel_edges[:-1], self.channel_edges[1:]])
        st = self.ScaledInstrument.from_arrays(self.matrix, energy_edges,
                                               self.channels, channel_edges,
                                               twod_edges=True)
        assert st._twod_edges is True
        assert st['alpha'] == 1.2
