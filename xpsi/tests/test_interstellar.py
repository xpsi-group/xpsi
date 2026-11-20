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
        """Interstellar inherits from ParameterSubspace, so its parameters plug
        into the xpsi parameter machinery.
        """
        assert issubclass(xpsi.Interstellar, xpsi.ParameterSubspace)

    def test_base_class_is_abstract(self):
        """Creating xpsi.Interstellar() directly raises TypeError, because
        attenuation() is abstract.
        """
        with pytest.raises(TypeError):
            xpsi.Interstellar()

    def test_subclass_without_attenuation_is_abstract(self):
        """A subclass that does not define attenuation() cannot be created, so
        the mistake shows up at construction time.
        """
        class NoAttenuation(xpsi.Interstellar):
            pass
        with pytest.raises(TypeError):
            NoAttenuation()

    def test_concrete_subclass_instantiates(self):
        """A subclass that defines attenuation() can be created and is an
        Interstellar instance.
        """
        assert isinstance(ConcreteInterstellar(), xpsi.Interstellar)

############################################
### Check the MultiplicativeInterstellar ###
############################################
class TestMultiplicativeInterstellar:
    def test_1d_signal_attenuated(self):
        """A 1D signal is multiplied by the constant attenuation factor."""
        sig = np.ones(5)
        MultiplicativeInterstellar(0.25)(np.linspace(1, 2, 5), sig)
        np.testing.assert_allclose(sig, 0.25)

    def test_1d_signal_is_modified_in_place(self):
        """The 1D signal array passed in is modified in place, not copied, as
        required by __call__.
        """
        sig = np.ones(5)
        ref = sig
        MultiplicativeInterstellar(0.5)(np.ones(5), sig)
        assert sig is ref
        np.testing.assert_allclose(ref, 0.5)

    def test_2d_signal_attenuated(self):
        """Every element of a 2D (energy x phase) signal is multiplied by the
        factor.
        """
        sig = np.ones((5, 3))
        MultiplicativeInterstellar(0.5)(np.ones(5), sig)
        np.testing.assert_allclose(sig, 0.5)

    def test_2d_signal_is_modified_in_place(self):
        """The 2D signal array is modified in place rather than replaced by a
        copy.
        """
        sig = np.ones((5, 3))
        ref = sig
        MultiplicativeInterstellar(0.5)(np.ones(5), sig)
        assert sig is ref

    def test_zero_attenuation_zeroes_signal(self):
        """An attenuation factor of 0 sets the whole signal to zero."""
        sig = np.ones((4, 2))
        MultiplicativeInterstellar(0.0)(np.ones(4), sig)
        assert np.all(sig == 0.0)
    
    def test_unit_attenuation_leaves_signal_unchanged(self):
        """An attenuation factor of 1 leaves a non-trivial signal exactly as it
        was.
        """
        sig = np.arange(8, dtype=float).reshape(4, 2)
        expected = sig.copy()
        MultiplicativeInterstellar(1.0)(np.ones(4), sig)
        np.testing.assert_array_equal(sig, expected)

#############################################
### Check the EnergyDependentInterstellar ###
#############################################
class TestEnergyDependentInterstellar:

    def test_2d_rows_scaled_by_energy_dependent_factor(self):
        """Each row of a 2D signal is scaled by the factor for its own energy,
        so attenuation is applied along the energy axis.
        """
        energies = np.array([1.0, 2.0, 3.0])
        sig = np.ones((3, 4))
        index = 2
        EnergyDependentInterstellar( index=index )(energies, sig)
        for col in range(4):
            np.testing.assert_allclose(sig[:, col], energies**(-index))

    def test_2d_columns_treated_identically(self):
        """Starting from distinct columns (phases), every column is scaled by
        the same energy-dependent factors.
        """
        energies = np.array([1.0, 2.0, 3.0])
        sig = np.arange(1.0, 10.0).reshape(3, 3)   # distinct columns
        original = sig.copy()
        EnergyDependentInterstellar()(energies, sig)
        for col in range(3):
            np.testing.assert_allclose(sig[:, col], original[:, col] * energies**(-1))

    def test_1d_energy_dependent_values(self):
        """A 1D signal gets the correct energy-dependent factor at each
        energy.
        """
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
        """A signal that is not a numpy array (list, tuple, number, None,
        string) raises TypeError.
        """
        ism = ism_nh(1.0)
        with pytest.raises(TypeError):
            ism(np.ones(2), bad)

    def test_3d_signal_raises_value_error(self):
        """A 3D signal is rejected with ValueError."""
        ism = ism_nh(1.0)
        with pytest.raises(ValueError):
            ism(np.ones(2), np.ones((2, 2, 2)))

    def test_0d_signal_raises_value_error(self):
        """A 0D (scalar) array signal is rejected with ValueError."""
        ism = ism_nh(1.0)
        with pytest.raises(ValueError):
            ism(np.ones(1), np.array(1.0))

    def test_row_mismatch_raises(self):
        """A 2D signal whose number of rows differs from the number of energies
        raises ValueError.
        """
        ism = ism_nh(1.0)
        with pytest.raises(ValueError):
            ism(np.ones(3), np.ones((4, 2)))

    def test_1d_length_mismatch_raises(self):
        """A 1D signal whose length differs from the number of energies raises
        ValueError.
        """
        ism = ism_nh(1.0)
        with pytest.raises(ValueError):
            ism(np.ones(3), np.ones(4))

    @pytest.mark.parametrize('shape', [(4,), (4, 2)])
    def test_single_energy_with_several_rows_raises(self, shape):
        """A single energy with a multi-row signal raises ValueError instead of
        being broadcast silently over every row. The signal is left
        unchanged.
        """
        # A length-1 attenuation would broadcast silently over any number
        # of rows, so this mismatch must be caught by an explicit check.
        ism = ism_nh(1.0)
        sig = np.ones(shape)
        with pytest.raises(ValueError, match='one row per energy'):
            ism(np.ones(1), sig)
        np.testing.assert_array_equal(sig, 1.0)   # left untouched

    def test_constructor_length_mismatch_raises(self):
        """Building the model with energy and attenuation arrays of different
        lengths fails with 'Array length mismatch'.
        """
        with pytest.raises(AssertionError, match='Array length mismatch'):
            CustomInterstellar(np.linspace(0.1, 12.0, 10),
                               np.linspace(0.1, 0.9, 9),
                               bounds={'neutral_hydrogen_column_density': (0.0, 50.0)})

    @pytest.mark.parametrize('nh', [-1.0, -1e-6, 50.1, 100.0])
    def test_parameter_limits(self, nh):
        """N_H values outside the strict bounds [0, 50] raise
        StrictBoundsError, on both sides of the range.
        """
        with pytest.raises(StrictBoundsError):
            ism_nh(nh)

    @pytest.mark.parametrize('nh', [0.0, 50.0])
    def test_parameter_limits_are_inclusive(self, nh):
        """The strict bounds themselves, N_H = 0 and N_H = 50, are accepted."""
        assert ism_nh(nh)['neutral_hydrogen_column_density'] == nh

    def test_increase_nh_absorbs_more( self ):
        """Over 0.3-2 keV, N_H = 20 removes more signal than N_H = 10 at every
        energy.
        """
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
        """With N_H = 1, the attenuation at the tabulated energies equals the
        tabulated values.
        """
        energies, table = self._table()
        np.testing.assert_allclose(ism_nh(1.0).attenuation(energies), table,
                                   rtol=1e-12)

    @pytest.mark.parametrize('nh', [0.5, 2.0, 7.3])
    def test_nh_scales_as_power_of_table(self, nh):
        """For several N_H, the attenuation equals the tabulated values raised
        to the power N_H.
        """
        energies, table = self._table()
        np.testing.assert_allclose(ism_nh(nh).attenuation(energies), table**nh,
                                   rtol=1e-10)

    def test_nh_zero_gives_unit_attenuation(self):
        """N_H = 0 gives an attenuation of exactly 1 at every energy."""
        energies, _ = self._table()
        np.testing.assert_array_equal(ism_nh(0.0).attenuation(energies), 1.0)

    def test_nh_zero_with_zero_in_table_gives_no_nan(self):
        """With a table containing a zero, N_H = 0 still gives 1 rather than
        NaN (0 * log(0)), and N_H = 1 keeps that entry at exactly 0.
        """
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
        """Before clipping, the interpolated table goes below 0 and above 1
        outside its range, so the clipping tests below are meaningful.
        """
        # Sanity check: without clipping, the Akima extrapolation of this table
        # goes below 0 and above 1, so the clipping tests below are meaningful.
        raw = ism_nh(1.0)._interpolator(np.array([0.01, 20.0]))
        assert raw[0] < 0.0
        assert raw[1] > 1.0

    @pytest.mark.parametrize('nh', [0.5, 1.0, 10.0])
    def test_attenuation_within_unit_interval_everywhere(self, nh):
        """From 0.001 to 100 keV and for several N_H, every attenuation value
        is finite and lies in [0, 1].
        """
        energies = np.concatenate([np.linspace(0.001, 0.1, 50),
                                   np.linspace(0.1, 12.0, 200),
                                   np.linspace(12.0, 100.0, 50)])
        att = ism_nh(nh).attenuation(energies)
        assert np.all(np.isfinite(att))
        assert np.all((att >= 0.0) & (att <= 1.0))

    def test_below_table_clipped_to_full_absorption(self):
        """Well below the tabulated range (0.01 and 0.05 keV), the attenuation
        is clipped to exactly 0.
        """
        energies = np.array([0.01, 0.05])
        np.testing.assert_array_equal(ism_nh(1.0).attenuation(energies), 0.0)

    def test_above_table_clipped_to_no_absorption(self):
        """Well above the tabulated range (20-100 keV), the attenuation is
        clipped to exactly 1.
        """
        energies = np.array([20.0, 50.0, 100.0])
        np.testing.assert_array_equal(ism_nh(1.0).attenuation(energies), 1.0)

    def test_extrapolation_is_enabled_just_outside_table(self):
        """Just below the table (0.09 keV), the result is a finite value in (0,
        1), not NaN, showing extrapolation is enabled.
        """
        # Just below 0.1 keV the extrapolated value is still inside (0, 1):
        # it must be returned as a finite number, not NaN or clipped.
        att = ism_nh(1.0).attenuation(np.array([0.09]))
        assert np.isfinite(att[0])
        assert 0.0 < att[0] < 1.0

    ### Per-energy-grid cache ###
    def test_same_grid_reuses_cache_entry(self):
        """Calling twice with the same energy array reuses the same cache entry
        instead of interpolating again.
        """
        ism = ism_nh(1.0)
        energies = np.linspace(0.5, 5.0, 10)
        ism.attenuation(energies)
        entry = ism._cache[id(energies)]
        ism.attenuation(energies)
        assert len(ism._cache) == 1
        assert ism._cache[id(energies)] is entry

    def test_changing_nh_updates_result_with_cached_grid(self):
        """After changing N_H on a cached grid, the result follows the new
        value, so the cache does not hold a result for an old N_H.
        """
        energies, table = self._table()
        ism = ism_nh(1.0)
        np.testing.assert_allclose(ism.attenuation(energies), table, rtol=1e-12)
        ism['neutral_hydrogen_column_density'] = 3.0
        np.testing.assert_allclose(ism.attenuation(energies), table**3,
                                   rtol=1e-10)
        assert len(ism._cache) == 1

    def test_different_grids_get_separate_entries(self):
        """Two different energy grids each get their own cache entry, with the
        right length and values.
        """
        ism = ism_nh(1.0)
        e1 = np.linspace(0.5, 5.0, 10)
        e2 = np.linspace(0.3, 8.0, 7)
        a1 = ism.attenuation(e1)
        a2 = ism.attenuation(e2)
        assert len(ism._cache) == 2
        assert a1.shape == (10,) and a2.shape == (7,)
        np.testing.assert_allclose(a1, ism._interpolate(e1), rtol=1e-12)
        np.testing.assert_allclose(a2, ism._interpolate(e2), rtol=1e-12)

    def test_cache_cleared_after_twenty_grids(self):
        """After 20 grids, a 21st clears the cache and leaves only the new
        entry, so the cache cannot grow without limit.
        """
        ism = ism_nh(1.0)
        grids = [np.linspace(0.5, 5.0, 10) for _ in range(20)]
        for g in grids:
            ism.attenuation(g)
        assert len(ism._cache) == 20
        extra = np.linspace(0.5, 5.0, 10)
        ism.attenuation(extra)
        assert len(ism._cache) == 1
        assert id(extra) in ism._cache

    def test_attenuation_returns_new_array_not_shared_buffer(self):
        """attenuation() returns a new array each time, so a kept result is not
        overwritten by a later call with a different N_H.
        """
        ism = ism_nh(1.0)
        energies = np.linspace(0.5, 5.0, 10)
        first = ism.attenuation(energies)
        kept = first.copy()
        buf = ism._cache[id(energies)][2]
        assert first is not buf

        ism['neutral_hydrogen_column_density'] = 5.0
        second = ism.attenuation(energies)
        assert second is not first
        # The earlier result must not be overwritten by the later call.
        np.testing.assert_array_equal(first, kept)
