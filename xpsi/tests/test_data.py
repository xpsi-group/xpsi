import os
import numpy as np

from xpsi.Data import Data

import pytest

# Make a simple Data to use for tests
COUNTS = np.array( [[ 1 , 2 ],
                    [ 3 , 4 ], 
                    [ 5 , 6 ]] )
CHANNELS = np.array( [ 30 , 31, 32 ] )
PHASES = np.array( [ 0.0 , 0.5, 1.0 ] )
FIRST = 30
LAST = 32
EXPOSURE = 12345.6789

@pytest.fixture
def data():
    return Data(counts=COUNTS, channels=CHANNELS, phases=PHASES, 
                first=FIRST, last=LAST, exposure_time=EXPOSURE)

def compare_data( data1, data2, phase_resolved_vs_averaged=False ):

    assert np.array_equal( data1.channels, data2.channels )
    assert data1._first - data1._last == data2._first - data2._last
    assert data1.exposure_time == data2.exposure_time

    if not phase_resolved_vs_averaged:
        assert np.array_equal( data1.counts, data2.counts )
        assert np.array_equal( data1.phases, data2.phases )
    else:
        phase_averaged_counts = data1.counts.sum(axis=1).reshape(-1,1)
        assert np.array_equal( phase_averaged_counts, data2.counts )

############################################################
### Make some test for loading data in different formats ###
############################################################
class TestDataLoading( object ):
    
    def test_init( self , data ):
        assert isinstance( data, Data )
        assert np.array_equal( data.counts, COUNTS )
        assert np.array_equal( data.channels, CHANNELS )
        assert np.array_equal( data.phases, PHASES )
        assert data._first == FIRST
        assert data._last == LAST
        assert data.exposure_time == EXPOSURE

    # Load a text file corresponding to data()
    def test_txt( self , data ):
        data_txt = Data.bin__event_list(
            path=os.path.join(os.path.dirname(os.path.abspath(__file__)),'data','Data.txt'),
            channels=CHANNELS,
            phases=PHASES,
            channel_column=1,
            phase_column=2,
            skiprows=2,
            first=FIRST,
            last=LAST,
            exposure_time=EXPOSURE
        )
        compare_data( data, data_txt )

    # Load an event file corresponding to data()
    def test_evt( self , data ):
        data_evt = Data.load( 
            os.path.join(os.path.dirname(os.path.abspath(__file__)),'data','Data.evt'),
            n_phases=2,
            channels=CHANNELS
            )
        compare_data( data, data_evt )

    # Load a phase-averaged spectrum corresponding to data()
    def test_pha( self , data ):
        data_pha = Data.load(
            os.path.join(os.path.dirname(os.path.abspath(__file__)),'data','Data.pha'),
            channels=CHANNELS
        )
        compare_data( data, data_pha, phase_resolved_vs_averaged=True )


################################################
### Check data class properties are well set ###
################################################
class TestDataProperties( object ):

    def test_properties( self , data ):
        assert data.counts.shape == (3,2)
        assert data.channels.shape == (3,)
        assert data.phases.shape == (3,)
        assert data._first == FIRST
        assert data._last == LAST
        assert data.exposure_time == EXPOSURE


#################################
### Check wrong input effects ###
#################################
class TestWrongInput( object ):
    
    def test_wrong_counts_shapes( self ):
        with pytest.raises( ValueError ):
            Data(counts=COUNTS.T, channels=CHANNELS, phases=PHASES, 
                 first=FIRST, last=LAST, exposure_time=EXPOSURE)

    def test_3d_counts( self ):
        with pytest.raises( TypeError ):
            Data(counts=COUNTS.reshape(3,2,1), channels=CHANNELS, phases=PHASES, 
                 first=FIRST, last=LAST, exposure_time=EXPOSURE)

    def test_negative_counts( self ):
        with pytest.raises( ValueError ):
            d = Data(counts= -COUNTS, channels=CHANNELS, phases=PHASES, 
                 first=FIRST, last=LAST, exposure_time=EXPOSURE)

    def test_zero_counts( self ):
        with pytest.raises( ValueError ):
            Data(counts= np.zeros_like( COUNTS ), channels=CHANNELS, phases=PHASES, 
                 first=FIRST, last=LAST, exposure_time=EXPOSURE)

    def test_first_larger_than_last( self ):
        with pytest.raises( ValueError ):
            Data(counts=COUNTS, channels=CHANNELS, phases=PHASES,
                 first=LAST, last=FIRST, exposure_time=EXPOSURE)

    def test_negative_exposure_time( self ):
        with pytest.raises( ValueError ):
            Data(counts=COUNTS, channels=CHANNELS, phases=PHASES, 
                 first=FIRST, last=LAST, exposure_time=-EXPOSURE)

#############################################
### Check some objects can not be changed ###
#############################################

class TestImmutableObjects( object ):

    def test_immutable_exposure( self , data ):
        with pytest.raises( AttributeError ):
            data.exposure_time = 1000.0

    def test_immutable_phases( self , data ):
        with pytest.raises( AttributeError ):
            data.phases = PHASES + 1

    def test_immutable_counts( self , data ):
        with pytest.raises( AttributeError ):
            data.counts = COUNTS


#######################################
### Check the trimming is done well ###
#######################################
class TestTrimming( object ):

    @pytest.mark.parametrize('bad', [56.2 , True , 'hello'])
    def test_non_int_channels( self , data, bad ):
        with pytest.raises( AssertionError ):
            data.trim_data( min_channel=bad, max_channel=32 )
            data.trim_data( min_channel=30, max_channel=bad )

    def test_min_greater_than_max( self , data ):
        with pytest.raises( AssertionError ):
            data.trim_data( min_channel=32, max_channel=30 )

    def test_size( self , data ):
        data.trim_data( min_channel=31, max_channel=32 )
        assert data.counts.shape == (2,2)
        assert data.channels.shape == (2,)

    def test_resulting_channels( self , data ):
        data.trim_data( min_channel=31, max_channel=32 )
        np.testing.assert_array_equal( data.channels , np.array([31, 32]) )

    def test_resulting_counts( self , data ):
        data.trim_data( min_channel=31, max_channel=32 )
        np.testing.assert_array_equal( data.counts , COUNTS[1:,:] )

##############################################
### Check the spectrum support computation ###
##############################################
BACKSCAL_SOURCE = 1.0
BACKSCAL_BACKGROUND = 5.0
SUPORT_EXPECTED_VALUES = np.array( [[0.0, 1.327777e-04],
                                    [0.0, 2.419835e-04],
                                    [1.701204e-05, 3.393880e-04]])

@pytest.fixture
def spectra_support( data ):

    # Set up the source and background
    background_source = data
    background_source.backscal = BACKSCAL_BACKGROUND

    # Make the support
    support = background_source.spectra_support( 3.0, BACKSCAL_SOURCE )
    return support

class TestSpectrumSupport( object ):
    
    def test_shape( self , spectra_support ):
        assert spectra_support.shape == (3,2)

    def test_positive_values( self , spectra_support ):
        assert (spectra_support >= 0.0).all()

    def test_computed_values( self , spectra_support ):
        assert np.allclose( spectra_support, SUPORT_EXPECTED_VALUES )