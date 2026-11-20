"""Unit tests for xpsi.Interstellar and for the usual CustomInterstellar.

Run with:  pytest -v test_interstellar.py

All parameter values are set at construction time (fixed parameters), which
keeps the tests independent of xpsi's likelihood/parameter-update machinery.
"""

import os
import numpy as np
import pytest

import xpsi

###########################
### Some pre-requisites ###
###########################
class ConcreteInterstellar(xpsi.Interstellar):

    def attenuation(self, energies):
        return None
    
class MultiplicativeInterstellar(xpsi.Interstellar):
    """Minimal concrete subclass of the base class: constant attenuation."""

    def __init__(self, factor=0.5):
        self.factor = factor
        super().__init__()

    def attenuation(self, energies):
        return np.ones_like(energies) * self.factor

class EnergyDependentInterstellar(xpsi.Interstellar):
    " Energy dependent attenuation case "

    def __init__( self , index = 1 ):
        self.index = index
        super().__init__()

    def attenuation(self, energies):
        return energies**(-self.index)
    

#############################################
### Check the ConcreteInterstellar exists ###
#############################################
class TestBaseClassCall:

    def test_is_subclass_of_parameter_subspace(self):
        assert issubclass(xpsi.Interstellar, xpsi.ParameterSubspace)

    def test_base_class_is_abstract(self):
        with pytest.raises(TypeError):
            xpsi.Interstellar()

    def test_subclass_without_attenuation_is_abstract(self):
        class NoAttenuation(xpsi.Interstellar):
            pass
        with pytest.raises(TypeError):
            NoAttenuation()

    def test_concrete_subclass_instantiates(self):
        assert isinstance(ConcreteInterstellar(), xpsi.Interstellar)

############################################
### Check the MultiplicativeInterstellar ###
############################################
class TestMultiplicativeInterstellar:
    def test_1d_signal_attenuated(self):
        sig = np.ones(5)
        MultiplicativeInterstellar(0.25)(np.linspace(1, 2, 5), sig)
        np.testing.assert_allclose(sig, 0.25)

    def test_1d_signal_is_modified_in_place(self):
        sig = np.ones(5)
        ref = sig
        MultiplicativeInterstellar(0.5)(np.ones(5), sig)
        assert sig is ref
        np.testing.assert_allclose(ref, 0.5)

    def test_2d_signal_attenuated(self):
        sig = np.ones((5, 3))
        MultiplicativeInterstellar(0.5)(np.ones(5), sig)
        np.testing.assert_allclose(sig, 0.5)

    def test_2d_signal_is_modified_in_place(self):
        sig = np.ones((5, 3))
        ref = sig
        MultiplicativeInterstellar(0.5)(np.ones(5), sig)
        assert sig is ref

    def test_zero_attenuation_zeroes_signal(self):
        sig = np.ones((4, 2))
        MultiplicativeInterstellar(0.0)(np.ones(4), sig)
        assert np.all(sig == 0.0)
    
    def test_unit_attenuation_leaves_signal_unchanged(self):
        sig = np.arange(8, dtype=float).reshape(4, 2)
        expected = sig.copy()
        MultiplicativeInterstellar(1.0)(np.ones(4), sig)
        np.testing.assert_array_equal(sig, expected)

#############################################
### Check the EnergyDependentInterstellar ###
#############################################
class TestEnergyDependentInterstellar:

    def test_2d_rows_scaled_by_energy_dependent_factor(self):
        energies = np.array([1.0, 2.0, 3.0])
        sig = np.ones((3, 4))
        index = 2
        EnergyDependentInterstellar( index=index )(energies, sig)
        for col in range(4):
            np.testing.assert_allclose(sig[:, col], energies**(-index))

    def test_2d_columns_treated_identically(self):
        energies = np.array([1.0, 2.0, 3.0])
        sig = np.arange(1.0, 10.0).reshape(3, 3)   # distinct columns
        original = sig.copy()
        EnergyDependentInterstellar()(energies, sig)
        for col in range(3):
            np.testing.assert_allclose(sig[:, col], original[:, col] * energies**(-1))

    def test_1d_energy_dependent_values(self):
        energies = np.array([0.5, 1.0, 2.0])
        sig = 2.0 * np.ones( 3 )
        index = 2.0
        EnergyDependentInterstellar( index=index )(energies, sig)
        np.testing.assert_allclose(sig, 2.0 * energies**(-index))


##########################################
### Check the usual interstellar class ###
##########################################

from scipy.interpolate import Akima1DInterpolator
from xpsi.Parameter import Parameter, StrictBoundsError

class CustomInterstellar(xpsi.Interstellar):

    """ Apply interstellar attenuation model tbnew usually used. """

    def __init__(self, energies, attenuation, bounds, values = None):

        if values is None: values = {}
        assert len(energies) == len(attenuation), 'Array length mismatch.'

        self._lkp_energies = energies # for lookup
        self._lkp_attenuation = attenuation # for lookup

        N_H = Parameter('neutral_hydrogen_column_density',
                        strict_bounds = (0.0, 50.0),
                        bounds = bounds.get('neutral_hydrogen_column_density', None),
                        doc = 'Neutral hydrogen column density in units of the fiducial column density',
                        symbol = r'$N_{\rm H}$',
                        value = values.get('neutral_hydrogen_column_density', None),
                        permit_prepend = False)

        self._interpolator = Akima1DInterpolator(self._lkp_energies,
                                                 self._lkp_attenuation)
        self._cache = {}
        self._interpolator.extrapolate = True

        super(CustomInterstellar, self).__init__(N_H)

    def _get_entry(self, energies):
        """
        Per-energy-grid cache keyed by object identity (a few grids).
        Useful when multiple instruments are fitted together.

        If first time with the energy grid, compute the interpolated attenuation.
        """

        # Search the cache for the same energy grid
        entry = self._cache.get(id(energies))

        # If not present, compute it
        if entry is None:
            grid = np.asarray(energies, dtype=np.double)
            att = self._interpolate(grid)

            # Check for division by zero
            with np.errstate(divide='ignore'):
                log_att = np.log(att)                  # att == 0 -> -inf, exp(-inf) = 0

            # Check if not too many entries, ptherwise likely an error
            if len(self._cache) >= 20:                  # safety against unbounded growth
                self._cache.clear()

            # Keep a reference to `energies` so its id cannot be reused by another object
            entry = (energies, log_att, np.empty_like(log_att))
            self._cache[id(energies)] = entry

        return entry

    def _compute(self, energies, out=None):
        """ Compute the attenuation using the cached values if needed. """

        # Get the attenuation reference
        _, log_att, buf = self._get_entry(energies)
        if out is None:
            out = buf

        # Apply the power of the neutral hydrogen column density
        nh = self['neutral_hydrogen_column_density']
        if nh == 0.0:                                  # att**0 == 1, avoids 0 * -inf = nan
            out.fill(1.0)
        else:
            np.multiply(log_att, nh, out=out)
            np.exp(out, out=out)

        return out

    def attenuation(self, energies):
        return self._compute(energies, out=np.empty_like(self._get_entry(energies)[1]))

    def _interpolate(self, energies):
        """ Helper. """
        _att = self._interpolator(energies)
        _att[_att < 0.0] = 0.0
        _att[_att > 1.0] = 1.0
        return _att

    @classmethod
    def load(cls, path,
             energy_column=0,
             attenuation_column=2,
             **kwargs):
        """ Load attenuation file. """
        temp = np.loadtxt(path, dtype=np.double)
        energies = temp[:,energy_column]
        attenuation = temp[:,attenuation_column]
        return cls(energies, attenuation, **kwargs)

def ism_nh(nh):
    inter = CustomInterstellar.load( 
        os.path.join( os.path.dirname(__file__), 'data', 'interstellar_phot_frac.txt' ),
        energy_column=0, 
        attenuation_column=2,
        bounds={'neutral_hydrogen_column_density': (0.0, 50.0)},
        values=None )
    inter['neutral_hydrogen_column_density'] = nh
    return inter

class TestCustomInterstellar:

    @pytest.mark.parametrize('bad', [[1.0, 2.0], (1.0, 2.0), 3.0, None, 'abc'])
    def test_non_ndarray_signal_raises_type_error(self, bad):
        ism = ism_nh(1.0)
        with pytest.raises(TypeError):
            ism(np.ones(2), bad)

    def test_3d_signal_raises_value_error(self):
        ism = ism_nh(1.0)
        with pytest.raises(ValueError):
            ism(np.ones(2), np.ones((2, 2, 2)))

    def test_0d_signal_raises_value_error(self):
        ism = ism_nh(1.0)
        with pytest.raises(ValueError):
            ism(np.ones(1), np.array(1.0))

    def test_row_mismatch_raises(self):
        ism = ism_nh(1.0)
        with pytest.raises(ValueError):
            ism(np.ones(3), np.ones((4, 2)))

    def test_parameter_limits(self):
        with pytest.raises(StrictBoundsError):
            ism_nh(100.0)

    def test_increase_nh_absorbs_more( self ):
        energies = np.linspace(0.3, 2.0, 5)
        sig1 = np.ones( (5 , 2) )
        sig2 = np.ones( (5 , 2) )
        ism_nh(10.0)(energies, sig1)
        ism_nh(20.0)(energies, sig2)
        assert (sig2 < sig1).all()
        
    ### Values against the tabulated model ###
    @staticmethod
    def _table():
        t = np.loadtxt(os.path.join(os.path.dirname(__file__), 'data',
                                    'interstellar_phot_frac.txt'))
        return t[:, 0], t[:, 2]

    def test_nh_one_reproduces_table(self):
        energies, table = self._table()
        np.testing.assert_allclose(ism_nh(1.0).attenuation(energies), table,
                                   rtol=1e-12)

    @pytest.mark.parametrize('nh', [0.5, 2.0, 7.3])
    def test_nh_scales_as_power_of_table(self, nh):
        energies, table = self._table()
        np.testing.assert_allclose(ism_nh(nh).attenuation(energies), table**nh,
                                   rtol=1e-10)

    def test_nh_zero_gives_unit_attenuation(self):
        energies, _ = self._table()
        np.testing.assert_array_equal(ism_nh(0.0).attenuation(energies), 1.0)

    def test_nh_zero_with_zero_in_table_gives_no_nan(self):
        # Table containing an exact zero: log(0) = -inf, and 0 * -inf = nan
        # unless the N_H = 0 case is handled separately.
        energies = np.array([1.0, 2.0, 3.0, 4.0])
        table = np.array([0.0, 0.2, 0.6, 0.9])
        ism = CustomInterstellar(energies, table,
                bounds={'neutral_hydrogen_column_density': (0.0, 50.0)},
                values={'neutral_hydrogen_column_density': 0.0})
        att = ism.attenuation(energies)
        assert not np.isnan(att).any()
        np.testing.assert_array_equal(att, 1.0)

        # And with N_H > 0, the zero entry stays exactly zero.
        ism['neutral_hydrogen_column_density'] = 1.0
        att = ism.attenuation(energies)
        assert not np.isnan(att).any()
        np.testing.assert_allclose(att, table, rtol=1e-12)

    ### Clipping and extrapolation outside the tabulated range (0.1-12 keV) ###
    def test_raw_extrapolation_leaves_unit_interval(self):
        # Sanity check: without clipping, the Akima extrapolation of this table
        # goes below 0 and above 1, so the clipping tests below are meaningful.
        raw = ism_nh(1.0)._interpolator(np.array([0.01, 20.0]))
        assert raw[0] < 0.0
        assert raw[1] > 1.0

    @pytest.mark.parametrize('nh', [0.5, 1.0, 10.0])
    def test_attenuation_within_unit_interval_everywhere(self, nh):
        energies = np.concatenate([np.linspace(0.001, 0.1, 50),
                                   np.linspace(0.1, 12.0, 200),
                                   np.linspace(12.0, 100.0, 50)])
        att = ism_nh(nh).attenuation(energies)
        assert np.all(np.isfinite(att))
        assert np.all((att >= 0.0) & (att <= 1.0))

    def test_below_table_clipped_to_full_absorption(self):
        energies = np.array([0.01, 0.05])
        np.testing.assert_array_equal(ism_nh(1.0).attenuation(energies), 0.0)

    def test_above_table_clipped_to_no_absorption(self):
        energies = np.array([20.0, 50.0, 100.0])
        np.testing.assert_array_equal(ism_nh(1.0).attenuation(energies), 1.0)

    def test_extrapolation_is_enabled_just_outside_table(self):
        # Just below 0.1 keV the extrapolated value is still inside (0, 1):
        # it must be returned as a finite number, not NaN or clipped.
        att = ism_nh(1.0).attenuation(np.array([0.09]))
        assert np.isfinite(att[0])
        assert 0.0 < att[0] < 1.0
