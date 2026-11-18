import numpy as np
import pytest
from astropy.io import fits

from xpsi.Data import Data
from xpsi.Signal import Signal
from xpsi.utils import fits_column_keyword


def _write_pha(path, channels, counts, tlmin=None, channel_first=True,
               quality=None, grouping=None):
    """ Write a minimal OGIP-like PHA file.

    :param int | None tlmin:
        Value of the TLMINn keyword of the CHANNEL column, not written if None.

    :param bool channel_first:
        If False, the CHANNEL column is the second column (so its TLMIN
        keyword is TLMIN2), and a decoy TLMIN1 is written for COUNTS.
    """
    channel_col = fits.Column(name='CHANNEL', array=np.asarray(channels), format='1J')
    counts_col = fits.Column(name='COUNTS', array=np.asarray(counts, dtype=float), format='D')
    cols = [channel_col, counts_col] if channel_first else [counts_col, channel_col]
    if quality is not None:
        cols.append(fits.Column(name='QUALITY', array=np.asarray(quality), format='1I'))
    if grouping is not None:
        cols.append(fits.Column(name='GROUPING', array=np.asarray(grouping), format='1I'))

    hdu = fits.BinTableHDU.from_columns(cols, name='SPECTRUM')
    for key, value in {'HDUCLAS1': 'SPECTRUM', 'EXPOSURE': 1000.0, 'BACKSCAL': 1.0,
                       'ANCRFILE': 'none', 'RESPFILE': 'none'}.items():
        hdu.header[key] = value
    if tlmin is not None:
        if channel_first:
            hdu.header['TLMIN1'] = tlmin
        else:
            hdu.header['TLMIN1'] = 999  # decoy, belongs to COUNTS
            hdu.header['TLMIN2'] = tlmin
    fits.HDUList([fits.PrimaryHDU(), hdu]).writeto(path, overwrite=True)


@pytest.mark.parametrize('tlmin, channel_first, expected', [
    (1, True, 1),
    (0, False, 0),
    (None, True, None),
])
def test_fits_column_keyword(tmp_path, tlmin, channel_first, expected):
    'TLMINn is read for the CHANNEL column, whatever its position in the table.'

    path = tmp_path / 'spectrum.pha'
    _write_pha(path, np.arange(5), np.ones(5), tlmin=tlmin, channel_first=channel_first)
    with fits.open(path) as hdul:
        assert fits_column_keyword(hdul['SPECTRUM'], 'CHANNEL', 'TLMIN') == expected
        assert fits_column_keyword(hdul['SPECTRUM'], 'channel', 'tlmin') == expected


def test_fits_column_keyword_missing_column(tmp_path):
    path = tmp_path / 'spectrum.pha'
    _write_pha(path, np.arange(5), np.ones(5), tlmin=0)
    with fits.open(path) as hdul:
        with pytest.raises(KeyError):
            fits_column_keyword(hdul['SPECTRUM'], 'PI', 'TLMIN')


@pytest.mark.parametrize('tlmin, channel_first', [
    (10, True),
    (10, False),
    (None, True),
])
def test_from_pha_quality_channel_numbering(tmp_path, tlmin, channel_first):
    'Quality-filtered channels are renumbered from the first channel of the file.'

    channels = np.arange(10, 20)
    counts = np.arange(100.0, 110.0)
    quality = np.zeros(10, dtype=int)
    quality[2] = 1
    path = tmp_path / 'spectrum.pha'
    _write_pha(path, channels, counts, tlmin=tlmin, channel_first=channel_first, quality=quality)

    data = Data.from_pha(str(path), checkqual=True)

    assert np.array_equal(data.channels, np.arange(10, 19))
    assert np.array_equal(data.counts[:, 0], counts[quality == 0])


@pytest.mark.parametrize('tlmin', [10, None])
def test_from_pha_grouping_with_nonzero_first_channel(tmp_path, tlmin):
    'Grouped counts are correct when channels do not start at 0, with or without TLMIN.'

    channels = np.arange(10, 20)
    counts = np.arange(1.0, 11.0)
    grouping = np.array([1, -1, -1, 1, -1, 1, -1, -1, -1, 1])
    path = tmp_path / 'spectrum.pha'
    _write_pha(path, channels, counts, tlmin=tlmin, grouping=grouping)

    data = Data.from_pha(str(path), checkgrp=True, minbincnt=1)

    starts = np.where(grouping == 1)[0]
    expected = np.add.reduceat(counts, starts)
    assert np.array_equal(data.counts[:, 0], expected)


def test_write_pha_tlmin_roundtrip(tmp_path, capsys):
    'Synthetic PHA files carry TLMIN1/TLMAX1, so Data.from_pha does not need the fallback.'

    channels = np.arange(10, 20)
    counts = np.arange(100.0, 110.0)
    filename = str(tmp_path / 'synthetic')
    # _write_PHA does not use the instance, so it can be called without building a Signal
    Signal._write_PHA(None, channels=channels, counts=counts, filename=filename,
                      instrument='TEST', exposure=1000.0, backscale=1.0)

    with fits.open(filename + '.pha') as hdul:
        assert hdul['SPECTRUM'].header['TLMIN1'] == 10
        assert hdul['SPECTRUM'].header['TLMAX1'] == 19

    data = Data.from_pha(filename + '.pha')
    assert 'TLMIN keyword not found' not in capsys.readouterr().out
    assert np.array_equal(data.channels, channels)
    assert np.array_equal(data.counts[:, 0], counts)


def _write_evt(tmp_path, channels, counts_per_channel, name='events'):
    """ Write an event file with the given counts per channel, split over 2 phase bins. """
    counts_per_channel = np.asarray(counts_per_channel)
    counts = np.stack([counts_per_channel // 2, counts_per_channel - counts_per_channel // 2], axis=1)
    filename = str(tmp_path / name)
    Signal._write_EVT(None, channels=np.asarray(channels), counts=counts, filename=filename,
                      instrument='TEST', exposure=1000.0, backscale=1.0)
    return filename + '.evt'


@pytest.mark.parametrize('evtcnts, expected, last_group_bad', [
    # Last group (channels 18-19) has 2 < 4 counts: flagged as bad and dropped
    ([3, 1, 1, 4, 2, 2, 1, 5, 1, 1], [4, 5, 4, 6], True),
    # Last group has exactly minbincnt counts: kept
    ([3, 1, 1, 4, 2, 2, 1, 5, 3, 1], [4, 5, 4, 6, 4], False),
])
def test_from_evt_minbincnt_grouping(tmp_path, evtcnts, expected, last_group_bad):
    'Grouping of an event list with a minimal number of counts per group.'

    path = _write_evt(tmp_path, np.arange(10, 20), evtcnts)
    data = Data.from_evt(path, n_phases=2, minbincnt=4)

    assert np.array_equal(data.counts.sum(axis=1), expected)
    assert np.array_equal(data.channels, np.arange(len(expected)))

    # Quality and grouping are indexed by channel number, channels outside the event list are bad
    assert len(data._quality) >= 20 and len(data._grouping) == len(data._quality)
    assert np.all(data._quality[:10] == 1) and np.all(data._quality[20:] == 1)
    expected_quality = np.zeros(10)
    if last_group_bad:
        expected_quality[-2:] = 1
    assert np.array_equal(data._quality[10:20], expected_quality)
    assert np.array_equal(data._grouping[10:20], [1, -1, 1, -1, 1, -1, 1, -1, 1, -1])


def test_from_evt_checkdatasrc_evt_source(tmp_path):
    'An event list grouped like a source event list gets the same grouping.'

    source_path = _write_evt(tmp_path, np.arange(10, 20), [3, 1, 1, 4, 2, 2, 1, 5, 1, 1], name='source')
    background_path = _write_evt(tmp_path, np.arange(10, 20), [1, 2, 3, 4, 5, 6, 7, 8, 9, 10], name='background')
    source = Data.from_evt(source_path, n_phases=2, minbincnt=4)
    background = Data.from_evt(background_path, n_phases=2, checkdatasrc=source)

    # Groups of the source: [10,11], [12,13], [14,15], [16,17]; channels 18-19 are bad
    assert np.array_equal(background.counts.sum(axis=1), [1+2, 3+4, 5+6, 7+8])
    assert np.array_equal(background.channels, source.channels)
    assert background._quality is source._quality
    assert background._minbincnt == 4


def test_from_evt_checkdatasrc_pha_source(tmp_path):
    'An event list grouped like a source PHA file with a bad channel inside a group.'

    channels = np.arange(10)
    counts = np.arange(1, 11)
    quality = np.array([0, 0, 0, 0, 1, 0, 0, 0, 0, 0])
    # Channel 4 is bad, so the group starting at channel 2 spans channels 2, 3, 5 and 6
    grouping = np.array([1, -1, 1, -1, 1, -1, -1, 1, -1, -1])
    pha_path = tmp_path / 'source.pha'
    _write_pha(pha_path, channels, counts, tlmin=0, quality=quality, grouping=grouping)
    source = Data.from_pha(str(pha_path), checkqual=True, checkgrp=True, minbincnt=1)

    evt_path = _write_evt(tmp_path, channels, counts)
    data = Data.from_evt(evt_path, n_phases=2, checkdatasrc=source)

    expected = [1+2, 3+4+6+7, 8+9+10]
    assert np.array_equal(source.counts[:, 0], expected)
    assert np.array_equal(data.counts.sum(axis=1), expected)
    assert np.array_equal(data.channels, source.channels)


def test_from_evt_checkdatasrc_channels_beyond_source(tmp_path):
    'A clear error is raised if the event list extends beyond the source quality information.'

    source_path = _write_evt(tmp_path, np.arange(0, 10), np.full(10, 2), name='source')
    source = Data.from_evt(source_path, n_phases=2, minbincnt=4)
    assert len(source._quality) == 10
    path = _write_evt(tmp_path, np.arange(0, 30), np.full(30, 5), name='wide')
    with pytest.raises(ValueError):
        Data.from_evt(path, n_phases=2, checkdatasrc=source)
