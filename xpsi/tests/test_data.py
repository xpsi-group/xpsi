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
