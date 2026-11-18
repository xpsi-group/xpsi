import numpy as np
import pytest

import xpsi
from xpsi.Spacetime import Spacetime
from xpsi.HotRegion import HotRegion
from xpsi.Elsewhere import Elsewhere
from xpsi.Everywhere import Everywhere
from xpsi.Photosphere import Photosphere
from xpsi.HotRegion import AtmosError as HotRegionAtmosError
from xpsi.Elsewhere import AtmosError as ElsewhereAtmosError
from xpsi.Everywhere import AtmosError as EverywhereAtmosError


ENERGIES = np.logspace(-0.5, 1.0, 20)


def _spacetime():
    bounds = dict(distance=(0.1, 1.0),
                  mass=(1.0, 3.0),
                  radius=(3.0 * 1.4766, 16.0),
                  cos_inclination=(0.0, 1.0))
    st = Spacetime(bounds=bounds, values=dict(frequency=300.0))
    st['distance'] = 0.3
    st['mass'] = 1.4
    st['radius'] = 12.0
    st['cos_inclination'] = 0.5
    return st


def _hot_region(atm_ext='BB'):
    bounds = dict(super_colatitude=(None, None),
                  super_radius=(None, None),
                  phase_shift=(-0.25, 0.75),
                  super_temperature=(None, None))
    hot = HotRegion(bounds=bounds, values={}, symmetry=True, omit=False,
                    cede=False, concentric=False, sqrt_num_cells=16,
                    min_sqrt_num_cells=10, max_sqrt_num_cells=32,
                    num_leaves=32, num_rays=128, atm_ext=atm_ext,
                    prefix='hot')
    hot['super_colatitude'] = 1.0
    hot['super_radius'] = 0.1
    hot['phase_shift'] = 0.0
    hot['super_temperature'] = 6.2
    return hot


def _everywhere():
    return Everywhere(time_invariant=True, values={'temperature': 6.0},
                      sqrt_num_cells=16, num_rays=64, atm_ext='Num4D')


def _elsewhere():
    return Elsewhere(values={'elsewhere_temperature': 6.0},
                     sqrt_num_cells=16, num_rays=64, atm_ext='Num4D')


def _photosphere(**components):
    ph = Photosphere(values=dict(mode_frequency=300.0), **components)
    xpsi.Star(_spacetime(), [ph]).update()
    return ph


def _fake_atmosphere(logT=(5.0, 5.5), logg=(13.7, 14.0)):
    """ Minimal (logT, logg, mu, logE, buf) tuple. Only the logT and logg
    grids are read by the bounds check in ``Photosphere.integrate``. """
    return (np.array(logT), np.array(logg),
            np.array([0.0, 1.0]), np.array([-1.0, 1.0]), np.zeros(16))


class TestPhotosphereUnloadedAtmosphere(object):
    """ A Num4D component without a loaded atmosphere table must raise the
    component's ``AtmosError``, not an ``IndexError`` from the table-bounds
    check in ``Photosphere.integrate``. """

    def test_everywhere_unloaded_raises_atmos_error(self):
        ph = _photosphere(hot=None, elsewhere=None, everywhere=_everywhere())
        with pytest.raises(EverywhereAtmosError):
            ph.integrate(ENERGIES, 1)

    def test_elsewhere_unloaded_raises_atmos_error(self):
        ph = _photosphere(hot=_hot_region('BB'), elsewhere=_elsewhere())
        with pytest.raises(ElsewhereAtmosError):
            ph.integrate(ENERGIES, 1)

    def test_hotregion_unloaded_raises_atmos_error(self):
        ph = _photosphere(hot=_hot_region('Num4D'), elsewhere=None)
        with pytest.raises(HotRegionAtmosError):
            ph.integrate(ENERGIES, 1)


class TestPhotosphereLoadedAtmosphereBounds(object):
    """ With a table loaded, the bounds check must still run. """

    def test_everywhere_temperature_out_of_bounds(self):
        ph = _photosphere(hot=None, elsewhere=None, everywhere=_everywhere())
        ph._everywhere_atmosphere = _fake_atmosphere(logT=(5.0, 5.5))
        with pytest.raises(Exception,
                           match='outside of the atmosphere table bounds'):
            ph.integrate(ENERGIES, 1)

    def test_elsewhere_temperature_out_of_bounds(self):
        ph = _photosphere(hot=_hot_region('BB'), elsewhere=_elsewhere())
        ph._elsewhere_atmosphere = _fake_atmosphere(logT=(5.0, 5.5))
        with pytest.raises(Exception,
                           match='outside of the atmosphere table bounds'):
            ph.integrate(ENERGIES, 1)

    def test_hotregion_temperature_out_of_bounds(self):
        ph = _photosphere(hot=_hot_region('Num4D'), elsewhere=None)
        ph._hot_atmosphere = _fake_atmosphere(logT=(5.0, 5.5))
        with pytest.raises(Exception,
                           match='outside of the atmosphere table bounds'):
            ph.integrate(ENERGIES, 1)
