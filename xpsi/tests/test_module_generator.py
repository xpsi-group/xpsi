""" Tests for the modules written by ``xpsi/module_generator.py``. """

import os
import sys
import subprocess
import importlib.util

import numpy as np
import pytest
from astropy.io import fits

import xpsi
from xpsi.Instrument import Instrument

from test_instrument import _write_rmf

ALPHA = 'energy_independent_effective_area_scaling_factor'


def _write_arf(path, detchans=5):
    """ Write a minimal OGIP-like ARF matching the energy grid of ``_write_rmf``. """
    n = detchans
    energ_lo = np.arange(1.0, 1.0 + n)
    cols = [fits.Column(name='ENERG_LO', array=energ_lo, format='E'),
            fits.Column(name='ENERG_HI', array=energ_lo + 1.0, format='E'),
            fits.Column(name='SPECRESP', array=10.0 * np.arange(1.0, 1.0 + n), format='E')]
    hdu = fits.BinTableHDU.from_columns(cols, name='SPECRESP')
    hdu.header['INSTRUME'] = 'TEST'
    fits.HDUList([fits.PrimaryHDU(), hdu]).writeto(path, overwrite=True)


@pytest.fixture(scope='module')
def CustomInstrument(tmp_path_factory):
    """ Run the module generator with ``--use-fits-format`` and import the
    generated ``CustomInstrument`` class. """
    outdir = tmp_path_factory.mktemp('generated')
    generator = os.path.join(os.path.dirname(xpsi.__file__), 'module_generator.py')
    subprocess.run([sys.executable, generator,
                    '--telescope', 'TEST',
                    '--instrument', 'TEST',
                    '--source', 'SRC',
                    '--frequency', '300',
                    '--use-fits-format',
                    '--hot-region-model', 'ST',
                    '--prefix', 'p',
                    '--is-antiphased', 'False',
                    '--hot-atmosphere-model', 'blackbody',
                    '--module-directory-path', str(outdir)],
                   cwd=str(outdir), check=True, capture_output=True)

    spec = importlib.util.spec_from_file_location('generated_CustomInstrument',
                                                  str(outdir / 'CustomInstrument.py'))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module.CustomInstrument


@pytest.mark.parametrize('tlmin, write_tlmin, extname', [
    (0, True, 'MATRIX'),
    (0, False, 'MATRIX'),
    (1, True, 'MATRIX'),
    (1, False, 'MATRIX'),
    (1, True, 'SPECRESP MATRIX'),
])
@pytest.mark.parametrize('twod_edges', [False, True])
def test_generated_from_ogip_fits(CustomInstrument, tmp_path, tlmin, write_tlmin, extname, twod_edges):
    'The generated loader matches Instrument.from_ogip_fits, and adds the scaling factor.'

    rmf = tmp_path / 'response.rmf'
    arf = tmp_path / 'response.arf'
    _write_rmf(str(rmf), tlmin=tlmin, write_tlmin=write_tlmin, extname=extname)
    _write_arf(str(arf))

    ref = Instrument.from_ogip_fits(str(rmf), ARF_path=str(arf), twod_edges=twod_edges)
    # Generated signature is (ARF_path, RMF_path, ...), and is called positionally by users
    instr = CustomInstrument.from_ogip_fits(str(arf), str(rmf), twod_edges=twod_edges)

    assert np.array_equal(instr.channels, np.arange(tlmin, tlmin + 5))
    assert np.array_equal(instr.channels, ref.channels)
    assert np.allclose(instr.matrix, ref.matrix)
    assert np.allclose(instr.matrix, np.diag(10.0 * np.arange(1.0, 6.0)))
    assert np.allclose(instr.channel_edges, ref.channel_edges)
    assert np.allclose(instr.energy_edges, ref.energy_edges)
    assert np.allclose(instr.RMF, ref.RMF)
    assert np.allclose(instr.ARF, ref.ARF)

    assert instr.names == [ALPHA]
    assert instr[ALPHA] == 1.0


def test_generated_from_ogip_fits_parameter(CustomInstrument, tmp_path):
    'The scaling factor can be free, is prefixed, and scales the matrix.'

    rmf = tmp_path / 'response.rmf'
    arf = tmp_path / 'response.arf'
    _write_rmf(str(rmf), tlmin=1)
    _write_arf(str(arf))

    instr = CustomInstrument.from_ogip_fits(str(arf), str(rmf),
                                            bounds={ALPHA: (0.8, 1.2)},
                                            prefix='TEST')

    assert instr.names == ['TEST__' + ALPHA]
    assert not instr.params[0].fixed
    assert instr.params[0].symbol == r'$\alpha_{\rm TEST}$'
    instr['TEST__' + ALPHA] = 1.1
    assert np.allclose(instr.construct_matrix(), 1.1 * instr.matrix)
