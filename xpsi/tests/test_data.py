import os
import numpy as np
from astropy.io import fits

from xpsi.Data import Data
from xpsi.Instrument import ChannelError

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

def write_pha( path, columns, **header_updates ):
    """ Write a synthetic PHA file, reusing the SPECTRUM header of Data.pha.

    :param columns: dict mapping column name to a (FITS format, array) tuple.
    :param header_updates: header keywords to add or overwrite.
    """
    with fits.open( os.path.join( DATA_DIR, 'Data.pha' ) ) as hdul:
        header = hdul['SPECTRUM'].header.copy()
    header.update( header_updates )

    # Column keywords (TTYPEn, TFORMn, ...) are rebuilt from the new columns
    cols = [ fits.Column( name=name, format=fmt, array=array )
             for name, (fmt, array) in columns.items() ]
    hdu = fits.BinTableHDU.from_columns( cols, header=header )
    fits.HDUList( [ fits.PrimaryHDU(), hdu ] ).writeto( path )
    return str( path )

# CHANNEL and COUNTS columns identical to Data.pha
PHA_COLUMNS = { 'CHANNEL' : ( 'J', CHANNELS ),
                'COUNTS'  : ( 'J', COUNTS.sum(axis=1) ) }

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
        """Constructing Data from arrays stores counts, channels, phases,
        first/last and exposure time unchanged.
        """
        assert isinstance( data, Data )
        assert np.array_equal( data.counts, COUNTS )
        assert np.array_equal( data.channels, CHANNELS )
        assert np.array_equal( data.phases, PHASES )
        assert data._first == FIRST
        assert data._last == LAST
        assert data.exposure_time == EXPOSURE

    def test_init_from_lists( self , data ):
        """Plain Python lists for counts, channels and phases are converted to
        numpy arrays, giving the same Data as the array inputs.
        """
        data_list = Data(counts=COUNTS.tolist(), channels=CHANNELS.tolist(),
                         phases=PHASES.tolist(), first=FIRST, last=LAST,
                         exposure_time=EXPOSURE)
        for attr in ( 'counts', 'channels', 'phases' ):
            assert isinstance( getattr( data_list, attr ), np.ndarray )
        compare_data( data, data_list )

    def test_txt( self , data ):
        """Binning the event list in Data.txt with ``bin__event_list``
        reproduces the reference phase-resolved Data.
        """
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

    def test_evt( self , data ):
        """Loading the event file Data.evt into two phase bins reproduces the
        reference Data, with an index range starting at 0.
        """
        data_evt = Data.load( 
            os.path.join( DATA_DIR, 'Data.evt' ),
            n_phases=2,
            channels=CHANNELS
            )
        compare_data( data, data_evt )
        assert data_evt.index_range == (0, 3)

    def test_pha( self , data ):
        """Loading the spectrum Data.pha gives the phase-summed reference
        counts in a single phase bin, with an index range starting at 0.
        """
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
        """The index range is (first, last + 1), so it can be used directly for
        array slicing.
        """
        # Upper bound is last + 1, for array slicing
        assert data.index_range == (FIRST, LAST + 1)

    def test_channel_range_is_deprecated_alias( self , data ):
        """The deprecated ``channel_range`` property returns the same value as
        ``index_range``.
        """
        assert data.channel_range == data.index_range


######################################
### Check the options of Data.load ###
######################################
class TestLoadOptions( object ):

    def test_unknown_hduclas1( self , tmp_path ):
        """A FITS file whose HDUCLAS1 keyword is neither an event list nor a
        spectrum raises IOError.
        """
        path = write_pha( tmp_path / 'unknown.pha', PHA_COLUMNS, HDUCLAS1='UNKNOWN' )
        with pytest.raises( IOError , match='HDUCLAS1' ):
            Data.load( path )

    def test_evt_all_channels( self , data ):
        """Without a ``channels`` argument, all channels present in the event
        list are used.
        """
        # Without channels, all channels present in the event list are used
        data_evt = Data.load( os.path.join( DATA_DIR, 'Data.evt' ), n_phases=2 )
        compare_data( data, data_evt )

    def test_evt_one_phase_bin( self , data ):
        """With ``n_phases=1``, the event file is binned into a single phase
        bin, giving the phase-summed counts.
        """
        data_evt = Data.load( os.path.join( DATA_DIR, 'Data.evt' ),
                              n_phases=1, channels=CHANNELS )
        compare_data( data, data_evt, phase_resolved_vs_averaged=True )

    @pytest.mark.parametrize('phase_shift, columns', [ (0.5, [1, 0]), (1.0, [0, 1]) ],
                             ids=['half_cycle', 'full_cycle'])
    def test_evt_phase_shift( self , phase_shift, columns ):
        """A half-cycle phase shift swaps the two phase bins of Data.evt, while
        a full-cycle shift leaves them unchanged.
        """
        # No event sits on a bin edge, so a half-cycle shift swaps the two
        # phase bins and a full-cycle shift changes nothing
        data_evt = Data.from_evt( os.path.join( DATA_DIR, 'Data.evt' ),
                                  n_phases=2, channels=CHANNELS,
                                  phase_shift=phase_shift )
        np.testing.assert_array_equal( data_evt.counts, COUNTS[:, columns] )

    def test_txt_phase_averaged( self , data ):
        """With ``phase_averaged=True``, all text-file events go into a single
        phase bin and no phase column is needed.
        """
        # All events go into a single phase bin; no phase column is needed
        data_txt = Data.bin__event_list(
            path=os.path.join( DATA_DIR, 'Data.txt' ),
            channels=CHANNELS,
            phases=np.array([ 0.0, 1.0 ]),
            channel_column=1,
            phase_averaged=True,
            skiprows=2,
            first=FIRST,
            last=LAST,
            exposure_time=EXPOSURE
        )
        compare_data( data, data_txt, phase_resolved_vs_averaged=True )

    def test_txt_phase_shift( self ):
        """A half-cycle phase shift applied to the text event list swaps the
        two phase bins, as for the event file.
        """
        # Same events as Data.evt: a half-cycle shift swaps the two phase bins
        data_txt = Data.bin__event_list(
            path=os.path.join( DATA_DIR, 'Data.txt' ),
            channels=CHANNELS,
            phases=PHASES,
            channel_column=1,
            phase_column=2,
            phase_shift=0.5,
            skiprows=2,
            first=FIRST,
            last=LAST,
            exposure_time=EXPOSURE
        )
        np.testing.assert_array_equal( data_txt.counts, COUNTS[:, [1, 0]] )


######################################
### Check the loading of PHA files ###
######################################
class TestPhaLoading( object ):

    PATH = os.path.join( DATA_DIR, 'Data.pha' )

    def test_channel_subset( self ):
        """Requesting a subset of channels keeps only the matching rows of the
        PHA file, with an index range restarting at 0.
        """
        data_pha = Data.load( self.PATH, channels=CHANNELS[1:] )
        np.testing.assert_array_equal( data_pha.channels, CHANNELS[1:] )
        np.testing.assert_array_equal( data_pha.counts, COUNTS[1:].sum(axis=1).reshape(-1,1) )
        assert data_pha.index_range == (0, 2)

    def test_missing_channel( self ):
        """Requesting a channel absent from the PHA file raises ValueError."""
        with pytest.raises( ValueError , match='Not all channels' ):
            Data.load( self.PATH, channels=np.array([ 31, 33 ]) )

    def test_metadata( self ):
        """The PHA header keywords (BACKSCAL, ANCRFILE, RESPFILE, BACKFILE,
        INSTRUME) are read into the corresponding Data attributes.
        """
        data_pha = Data.load( self.PATH )
        assert data_pha.backscal == 1
        assert data_pha.ancrfile == 'NONE'
        assert data_pha.respfile == 'NONE'
        assert data_pha.backfile == 'NONE'
        assert data_pha.instrument == 'UNKNOWN'

    def test_evt_instrument( self ):
        """The instrument name is also read from the header of an event file.
        """
        data_evt = Data.load( os.path.join( DATA_DIR, 'Data.evt' ),
                              n_phases=2, channels=CHANNELS )
        assert data_evt.instrument == 'UNKNOWN'

    def test_no_backfile_for_background( self , tmp_path ):
        """BACKFILE is only read for a total (source) spectrum, so a background
        spectrum (HDUCLAS2='BKG') has no ``backfile`` attribute.
        """
        # BACKFILE is only read for a total (source) spectrum
        path = write_pha( tmp_path / 'bkg.pha', PHA_COLUMNS, HDUCLAS2='BKG' )
        assert not hasattr( Data.load( path ), 'backfile' )

    @pytest.mark.parametrize('option, column', [ ('withuncs', 'STAT_ERR'),
                                                 ('checkqual', 'QUALITY'),
                                                 ('checkgrp', 'GROUPING') ])
    def test_missing_column( self , option, column ):
        """Requesting uncertainties, quality or grouping from a PHA file
        lacking the STAT_ERR, QUALITY or GROUPING column raises IOError naming
        that column.
        """
        # Data.pha has only CHANNEL and COUNTS columns
        with pytest.raises( IOError , match=column ):
            Data.load( self.PATH, **{ option : True } )

    def test_rate_column( self , tmp_path ):
        """A PHA file with a RATE column instead of COUNTS is loaded as rates,
        with the exposure time set to 1.
        """
        # With RATE instead of COUNTS, the exposure time is set to 1
        rates = np.array([ 0.5, 1.5, 2.5 ])
        path = write_pha( tmp_path / 'rate.pha',
                          { 'CHANNEL' : ( 'J', CHANNELS ), 'RATE' : ( 'E', rates ) } )
        data_pha = Data.load( path )
        assert data_pha.exposure_time == 1.0
        np.testing.assert_array_equal( data_pha.counts, rates.reshape(-1,1) )


###################################################
### Check the quality and grouping of PHA files ###
###################################################
# Six channels; the last one is flagged bad; three groups: 30-31, 32-34, 35
GRP_CHANNELS = np.arange( 30, 36 )
GRP_COUNTS = np.arange( 1, 7 )
GRP_QUALITY = np.array([ 0, 0, 0, 0, 0, 1 ])
GRP_GROUPING = np.array([ 1, -1, 1, -1, -1, 1 ])
GRP_MINBINCNT = 3

@pytest.fixture
def grouped_pha( tmp_path ):
    columns = { 'CHANNEL'  : ( 'J', GRP_CHANNELS ),
                'COUNTS'   : ( 'J', GRP_COUNTS ),
                'STAT_ERR' : ( 'E', np.sqrt( GRP_COUNTS ) ),
                'QUALITY'  : ( 'I', GRP_QUALITY ),
                'GROUPING' : ( 'I', GRP_GROUPING ) }
    return write_pha( tmp_path / 'grouped.pha', columns,
                      TLMAX1=GRP_CHANNELS[-1], DETCHANS=len(GRP_CHANNELS) )

class TestPhaQualityGrouping( object ):

    def test_quality( self , grouped_pha ):
        """With ``checkqual=True``, channels flagged as bad by the QUALITY
        column are dropped and the quality array is stored.
        """
        # The bad channel 35 is dropped
        data_pha = Data.load( grouped_pha, checkqual=True )
        np.testing.assert_array_equal( data_pha.channels, GRP_CHANNELS[:-1] )
        np.testing.assert_array_equal( data_pha.counts, GRP_COUNTS[:-1].reshape(-1,1) )
        np.testing.assert_array_equal( data_pha._quality, GRP_QUALITY )

    def test_grouping( self , grouped_pha ):
        """With ``checkgrp=True``, channels are summed into groups renumbered
        from 0, and uncertainties are summed in quadrature.
        """
        # Groups are renumbered from 0; errors are summed in quadrature
        data_pha = Data.load( grouped_pha, checkgrp=True, withuncs=True,
                              minbincnt=GRP_MINBINCNT )
        np.testing.assert_array_equal( data_pha.channels, [ 0, 1, 2 ] )
        np.testing.assert_array_equal( data_pha.counts, [ [3], [12], [6] ] )
        np.testing.assert_allclose( data_pha.errors, np.sqrt([ [3], [12], [6] ]), rtol=1e-6 )
        np.testing.assert_array_equal( data_pha._grouping, GRP_GROUPING )
        assert data_pha._minbincnt == GRP_MINBINCNT
        assert data_pha._grpdata

    def test_quality_and_grouping( self , grouped_pha ):
        """Bad channels are removed before grouping, so a group made only of
        bad channels disappears.
        """
        # The bad channel 35 is removed before grouping, so its group disappears
        data_pha = Data.load( grouped_pha, checkqual=True, checkgrp=True,
                              minbincnt=GRP_MINBINCNT )
        np.testing.assert_array_equal( data_pha.channels, [ 0, 1 ] )
        np.testing.assert_array_equal( data_pha.counts, [ [3], [12] ] )

    def test_grouping_requires_minbincnt( self , grouped_pha ):
        """Grouping without a ``minbincnt`` value raises AttributeError."""
        with pytest.raises( AttributeError , match='grouping' ):
            Data.load( grouped_pha, checkgrp=True )


#################################
### Check wrong input effects ###
#################################
class TestWrongInput( object ):
    
    def test_wrong_counts_rows( self ):
        """A counts matrix whose number of rows does not match first/last
        raises ValueError.
        """
        # 2 rows, but first/last require 3
        with pytest.raises( ValueError , match='number of rows' ):
            Data(counts=COUNTS.T, channels=CHANNELS, phases=PHASES, 
                 first=FIRST, last=LAST, exposure_time=EXPOSURE)

    def test_wrong_counts_columns( self ):
        """A counts matrix whose number of columns does not match the number of
        phase bins raises ValueError.
        """
        # 3 columns, but 3 phase edges give only 2 phase bins
        counts = np.hstack( ( COUNTS, COUNTS[:,:1] ) )
        with pytest.raises( ValueError , match='number of columns' ):
            Data(counts=counts, channels=CHANNELS, phases=PHASES, 
                 first=FIRST, last=LAST, exposure_time=EXPOSURE)

    def test_3d_counts( self ):
        """A three-dimensional counts array raises TypeError."""
        with pytest.raises( TypeError ):
            Data(counts=COUNTS.reshape(3,2,1), channels=CHANNELS, phases=PHASES, 
                 first=FIRST, last=LAST, exposure_time=EXPOSURE)

    @pytest.mark.parametrize('counts', [ -COUNTS,
                                         COUNTS * np.array([[1, 1], [1, -1], [1, 1]]) ],
                             ids=['all_negative', 'one_negative'])
    def test_negative_counts( self , counts ):
        """Counts with any negative value raise ValueError."""
        with pytest.raises( ValueError , match='Negative count' ):
            Data(counts=counts, channels=CHANNELS, phases=PHASES, 
                 first=FIRST, last=LAST, exposure_time=EXPOSURE)

    def test_zero_counts( self ):
        """A counts matrix containing only zeros raises ValueError."""
        with pytest.raises( ValueError ):
            Data(counts= np.zeros_like( COUNTS ), channels=CHANNELS, phases=PHASES, 
                 first=FIRST, last=LAST, exposure_time=EXPOSURE)

    def test_2d_phases( self ):
        """A two-dimensional phases array raises ValueError, since phase edges
        must be one-dimensional.
        """
        with pytest.raises( ValueError , match='one-dimensional' ):
            Data(counts=COUNTS, channels=CHANNELS, phases=PHASES.reshape(1,-1),
                 first=FIRST, last=LAST, exposure_time=EXPOSURE)

    @pytest.mark.parametrize('first, last', [ (None, LAST), (FIRST, None),
                                              ([FIRST], LAST) ],
                             ids=['first_None', 'last_None', 'first_list'])
    def test_non_integer_first_last( self , first, last ):
        """A ``first`` or ``last`` that cannot be converted to an integer
        (None, list) raises TypeError.
        """
        with pytest.raises( TypeError , match='must be integers' ):
            Data(counts=COUNTS, channels=CHANNELS, phases=PHASES,
                 first=first, last=last, exposure_time=EXPOSURE)

    def test_non_numeric_string_first( self ):
        """A non-numeric string for ``first`` fails in ``int()`` with a
        ValueError, not the TypeError raised for other non-integer types.
        """
        with pytest.raises( ValueError , match='invalid literal' ):
            Data(counts=COUNTS, channels=CHANNELS, phases=PHASES,
                 first='a', last=LAST, exposure_time=EXPOSURE)

    def test_first_larger_than_last( self ):
        """A ``first`` larger than ``last`` raises ValueError."""
        with pytest.raises( ValueError ):
            Data(counts=COUNTS, channels=CHANNELS, phases=PHASES,
                 first=LAST, last=FIRST, exposure_time=EXPOSURE)

    @pytest.mark.parametrize('exposure', [ -EXPOSURE, 0.0 ],
                             ids=['negative', 'zero'])
    def test_non_positive_exposure_time( self , exposure ):
        """A negative or zero exposure time raises ValueError."""
        with pytest.raises( ValueError , match='exposure time' ):
            Data(counts=COUNTS, channels=CHANNELS, phases=PHASES, 
                 first=FIRST, last=LAST, exposure_time=exposure)

#############################################
### Check some objects cannot be changed ###
#############################################

class TestImmutableObjects( object ):

    def test_immutable_exposure( self , data ):
        """The exposure time is read-only after construction."""
        with pytest.raises( AttributeError ):
            data.exposure_time = 1000.0

    def test_immutable_phases( self , data ):
        """The phase edges are read-only after construction."""
        with pytest.raises( AttributeError ):
            data.phases = PHASES + 1

    def test_immutable_counts( self , data ):
        """The counts are read-only after construction."""
        with pytest.raises( AttributeError ):
            data.counts = COUNTS


##################################################
### Check the channels are validated when set ###
##################################################
class TestChannels( object ):

    @pytest.mark.parametrize('channels', [ -CHANNELS, CHANNELS.reshape(1,-1) ],
                             ids=['negative', '2d'])
    def test_invalid_channels_at_init( self , channels ):
        """Negative or two-dimensional channels raise ChannelError at
        construction.
        """
        with pytest.raises( ChannelError ):
            Data(counts=COUNTS, channels=channels, phases=PHASES,
                 first=FIRST, last=LAST, exposure_time=EXPOSURE)

    @pytest.mark.parametrize('channels', [ -CHANNELS, CHANNELS.reshape(1,-1) ],
                             ids=['negative', '2d'])
    def test_setter_revalidates( self , data, channels ):
        """Assigning invalid channels after construction raises ChannelError
        and leaves the existing channels unchanged.
        """
        with pytest.raises( ChannelError ):
            data.channels = channels
        # The failed assignment leaves the channels unchanged
        np.testing.assert_array_equal( data.channels, CHANNELS )

    def test_setter_accepts_list( self , data ):
        """Assigning a plain list of channels stores it as a numpy array."""
        data.channels = [ 40, 41, 42 ]
        assert isinstance( data.channels, np.ndarray )
        np.testing.assert_array_equal( data.channels, [ 40, 41, 42 ] )

    def test_non_contiguous_channels_warning( self , monkeypatch, capsys ):
        """Channels that do not increment by one print a warning when verbosity
        is on.
        """
        # The warning is only printed when verbosity is on
        monkeypatch.setattr( 'xpsi.utils._verbose', True )
        Data(counts=COUNTS, channels=np.array([ 30, 32, 34 ]), phases=PHASES,
             first=FIRST, last=LAST, exposure_time=EXPOSURE)
        assert 'do not uniformly increment by one' in capsys.readouterr().out


#######################################
### Check the trimming is done well ###
#######################################
class TestTrimming( object ):

    # Values are within the channel range so only the type check can fail.
    # Note: booleans are not tested, since ``isinstance(True, int)`` is True.
    @pytest.mark.parametrize('bad', [31.0 , '31' , None , [31]])
    def test_non_int_min_channel( self , data, bad ):
        """A non-integer ``min_channel`` passed to ``trim_data`` raises
        AssertionError (booleans are not tested, since they are ints).
        """
        with pytest.raises( AssertionError ):
            data.trim_data( min_channel=bad, max_channel=32 )

    @pytest.mark.parametrize('bad', [31.0 , '31' , None , [31]])
    def test_non_int_max_channel( self , data, bad ):
        """A non-integer ``max_channel`` passed to ``trim_data`` raises
        AssertionError.
        """
        with pytest.raises( AssertionError ):
            data.trim_data( min_channel=30, max_channel=bad )

    def test_min_greater_than_max( self , data ):
        """A ``min_channel`` larger than ``max_channel`` raises AssertionError.
        """
        with pytest.raises( AssertionError ):
            data.trim_data( min_channel=32, max_channel=30 )

    def test_size( self , data ):
        """Trimming to channels 31-32 reduces counts and channels to two rows.
        """
        data.trim_data( min_channel=31, max_channel=32 )
        assert data.counts.shape == (2,2)
        assert data.channels.shape == (2,)

    def test_resulting_channels( self , data ):
        """Trimming keeps exactly the channels in the requested range."""
        data.trim_data( min_channel=31, max_channel=32 )
        np.testing.assert_array_equal( data.channels , np.array([31, 32]) )

    def test_resulting_counts( self , data ):
        """Trimming keeps the counts rows of the retained channels."""
        data.trim_data( min_channel=31, max_channel=32 )
        np.testing.assert_array_equal( data.counts , COUNTS[1:,:] )

    def test_resulting_index_range( self , data ):
        """After trimming, the index range restarts at 0 and spans the kept
        rows.
        """
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
        """The spectrum support has one (lower, upper) row per channel."""
        assert spectra_support.shape == (3,2)

    def test_positive_values( self , spectra_support ):
        """All bounds of the spectrum support are non-negative, since lower
        bounds are clipped at 0.
        """
        assert (spectra_support >= 0.0).all()

    def test_computed_values( self , spectra_support ):
        """The support matches the expected Poisson interval on the phase-
        summed counts, converted to a rate and rescaled by the BACKSCAL ratio.
        """
        assert np.allclose( spectra_support, SUPPORT_EXPECTED_VALUES )

    def test_channel_selection( self ):
        """With a ``channels`` argument, the support rows are returned in the
        order of the requested channels.
        """
        # Rows are returned in the order of the requested channels
        background = Data(counts=COUNTS, channels=CHANNELS, phases=PHASES,
                          first=FIRST, last=LAST, exposure_time=EXPOSURE)
        background.backscal = BACKSCAL_BACKGROUND
        support = background.spectra_support( N_SIGMA, BACKSCAL_SOURCE,
                                              channels=[ 32, 30 ] )
        assert np.allclose( support, SUPPORT_EXPECTED_VALUES[[2, 0]] )

    @pytest.mark.parametrize('zero_row, donor_row', [ (0, 1), (1, 2), (2, 1) ],
                             ids=['first', 'middle', 'last'])
    def test_zero_counts_upper_bound( self , zero_row, donor_row ):
        """A channel with zero counts keeps a lower bound of 0, and its upper
        bound is copied from the nearest channel (the next one first).
        """
        # A channel with no counts would have a [0, 0] support; its upper
        # bound is copied from the nearest channel, the next one first
        counts = COUNTS.copy()
        counts[zero_row] = 0
        background = Data(counts=counts, channels=CHANNELS, phases=PHASES,
                          first=FIRST, last=LAST, exposure_time=EXPOSURE)
        background.backscal = BACKSCAL_BACKGROUND
        support = background.spectra_support( N_SIGMA, BACKSCAL_SOURCE )
        assert support[zero_row, 0] == 0.0
        assert support[zero_row, 1] == support[donor_row, 1] > 0.0

    def test_missing_backscal( self , data ):
        """Computing the support of a Data without BACKSCAL raises IOError."""
        with pytest.raises( IOError , match='BACKSCAL' ):
            data.spectra_support( N_SIGMA, BACKSCAL_SOURCE )


######################################
### Check the plotting of the data ###
######################################
@pytest.fixture
def agg_backend():
    # Use a non-interactive backend, then restore the previous one
    import matplotlib.pyplot as plt
    backend = plt.get_backend()
    plt.switch_backend( 'Agg' )
    yield
    plt.close( 'all' )
    plt.switch_backend( backend )

class TestPlot( object ):

    def test_plot( self , data, agg_backend ):
        """``Data.plot`` returns a figure at 200 dpi with a 2x2 grid of axes.
        """
        fig, axs = data.plot()
        assert axs.shape == (2, 2)
        assert fig.get_dpi() == 200