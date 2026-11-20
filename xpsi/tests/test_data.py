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

# Folder containing the test data files
DATA_DIR = os.path.join( os.path.dirname( os.path.abspath(__file__) ), 'data' )

@pytest.fixture
def data():
    return Data(counts=COUNTS, channels=CHANNELS, phases=PHASES, 
                first=FIRST, last=LAST, exposure_time=EXPOSURE)

def compare_data( data1, data2, phase_resolved_vs_averaged=False ):

    assert np.array_equal( data1.channels, data2.channels )
    assert data1.exposure_time == data2.exposure_time

    # The absolute index range may differ (loaders always start at 0), but
    # each must span exactly the rows of its counts matrix
    for d in ( data1, data2 ):
        start, stop = d.index_range
        assert stop - start == d.counts.shape[0]

    if not phase_resolved_vs_averaged:
        assert np.array_equal( data1.counts, data2.counts )
        assert np.array_equal( data1.phases, data2.phases )
    else:
        phase_averaged_counts = data1.counts.sum(axis=1).reshape(-1,1)
        assert np.array_equal( phase_averaged_counts, data2.counts )
        assert np.array_equal( data2.phases, np.array([0.0, 1.0]) )

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
            path=os.path.join( DATA_DIR, 'Data.txt' ),
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
            os.path.join( DATA_DIR, 'Data.evt' ),
            n_phases=2,
            channels=CHANNELS
            )
        compare_data( data, data_evt )
        assert data_evt.index_range == (0, 3)

    # Load a phase-averaged spectrum corresponding to data()
    def test_pha( self , data ):
        data_pha = Data.load(
            os.path.join( DATA_DIR, 'Data.pha' ),
            channels=CHANNELS
        )
        compare_data( data, data_pha, phase_resolved_vs_averaged=True )
        assert data_pha.index_range == (0, 3)


################################################
### Check data class properties are well set ###
################################################
class TestDataProperties( object ):

    def test_index_range( self , data ):
        # Upper bound is last + 1, for array slicing
        assert data.index_range == (FIRST, LAST + 1)

    def test_channel_range_is_deprecated_alias( self , data ):
        assert data.channel_range == data.index_range


#################################
### Check wrong input effects ###
#################################
class TestWrongInput( object ):
    
    def test_wrong_counts_rows( self ):
        # 2 rows, but first/last require 3
        with pytest.raises( ValueError , match='number of rows' ):
            Data(counts=COUNTS.T, channels=CHANNELS, phases=PHASES, 
                 first=FIRST, last=LAST, exposure_time=EXPOSURE)

    def test_wrong_counts_columns( self ):
        # 3 columns, but 3 phase edges give only 2 phase bins
        counts = np.hstack( ( COUNTS, COUNTS[:,:1] ) )
        with pytest.raises( ValueError , match='number of columns' ):
            Data(counts=counts, channels=CHANNELS, phases=PHASES, 
                 first=FIRST, last=LAST, exposure_time=EXPOSURE)

    def test_3d_counts( self ):
        with pytest.raises( TypeError ):
            Data(counts=COUNTS.reshape(3,2,1), channels=CHANNELS, phases=PHASES, 
                 first=FIRST, last=LAST, exposure_time=EXPOSURE)

    @pytest.mark.parametrize('counts', [ -COUNTS,
                                         COUNTS * np.array([[1, 1], [1, -1], [1, 1]]) ],
                             ids=['all_negative', 'one_negative'])
    def test_negative_counts( self , counts ):
        with pytest.raises( ValueError , match='Negative count' ):
            Data(counts=counts, channels=CHANNELS, phases=PHASES, 
                 first=FIRST, last=LAST, exposure_time=EXPOSURE)

    def test_zero_counts( self ):
        with pytest.raises( ValueError ):
            Data(counts= np.zeros_like( COUNTS ), channels=CHANNELS, phases=PHASES, 
                 first=FIRST, last=LAST, exposure_time=EXPOSURE)

    def test_first_larger_than_last( self ):
        with pytest.raises( ValueError ):
            Data(counts=COUNTS, channels=CHANNELS, phases=PHASES,
                 first=LAST, last=FIRST, exposure_time=EXPOSURE)

    @pytest.mark.parametrize('exposure', [ -EXPOSURE, 0.0 ],
                             ids=['negative', 'zero'])
    def test_non_positive_exposure_time( self , exposure ):
        with pytest.raises( ValueError , match='exposure time' ):
            Data(counts=COUNTS, channels=CHANNELS, phases=PHASES, 
                 first=FIRST, last=LAST, exposure_time=exposure)

#############################################
### Check some objects cannot be changed ###
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

    # Values are within the channel range so only the type check can fail.
    # Note: booleans are not tested, since ``isinstance(True, int)`` is True.
    @pytest.mark.parametrize('bad', [31.0 , '31' , None , [31]])
    def test_non_int_min_channel( self , data, bad ):
        with pytest.raises( AssertionError ):
            data.trim_data( min_channel=bad, max_channel=32 )

    @pytest.mark.parametrize('bad', [31.0 , '31' , None , [31]])
    def test_non_int_max_channel( self , data, bad ):
        with pytest.raises( AssertionError ):
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

    def test_resulting_index_range( self , data ):
        # After trimming, the index range restarts at 0 and spans the kept rows
        data.trim_data( min_channel=31, max_channel=32 )
        assert data.index_range == (0, 2)

##############################################
### Check the spectrum support computation ###
##############################################
BACKSCAL_SOURCE = 1.0
BACKSCAL_BACKGROUND = 5.0
N_SIGMA = 3.0

# Expected support: Poisson interval [c - n*sqrt(c), c + n*sqrt(c)] on the
# phase-summed counts c, lower bound clipped at 0, converted to a count rate
# and rescaled by the ratio of source to background BACKSCAL.
# With c = [3, 7, 11], only the last lower bound is non-zero.
_C = COUNTS.sum(axis=1).astype(float)
SUPPORT_EXPECTED_VALUES = ( np.column_stack( ( np.clip( _C - N_SIGMA*np.sqrt(_C), 0.0, None ),
                                               _C + N_SIGMA*np.sqrt(_C) ) )
                            / EXPOSURE * ( BACKSCAL_SOURCE / BACKSCAL_BACKGROUND ) )

@pytest.fixture
def spectra_support():

    # Build a dedicated background Data object, with its BACKSCAL
    background = Data(counts=COUNTS, channels=CHANNELS, phases=PHASES,
                      first=FIRST, last=LAST, exposure_time=EXPOSURE)
    background.backscal = BACKSCAL_BACKGROUND

    # Make the support
    return background.spectra_support( N_SIGMA, BACKSCAL_SOURCE )

class TestSpectrumSupport( object ):
    
    def test_shape( self , spectra_support ):
        assert spectra_support.shape == (3,2)

    def test_positive_values( self , spectra_support ):
        assert (spectra_support >= 0.0).all()

    def test_computed_values( self , spectra_support ):
        assert np.allclose( spectra_support, SUPPORT_EXPECTED_VALUES )