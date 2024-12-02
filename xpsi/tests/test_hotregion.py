import numpy as np
from xpsi.HotRegion import HotRegion
import pytest

class TestHotRegion(object):
    
    def setup_class(cls):
        cls.super_colatitude_bounds = (None, None)
        cls.super_radius_bounds = (None, None)
        cls.phase_shift_bounds = (0.0, 0.1)
        cls.super_temperature_bounds = (5.1, 6.8)

        cls.bounds = {'super_colatitude': cls.super_colatitude_bounds,
                      'super_radius': cls.super_radius_bounds,
                      'phase_shift': cls.phase_shift_bounds,
                      'super_temperature': cls.super_temperature_bounds,}
	
        cls.values = {}
        cls.HotRegion = HotRegion(cls.bounds, cls.values)
    
    def test_hotregion_class_initializes(self):
        HotRegion(self.bounds, self.values)
    
    def test_hotregion_breaks_with_bounds_missing(self):
        with pytest.raises(TypeError):
            HotRegion(self.values)
            
    def test_hotregion_breaks_with_values_missing(self):
        with pytest.raises(TypeError):
            HotRegion(self.bounds)
            
    def test_hotregion_breaks_if_phase_shift_bounds_unspecified(self):
        bounds = self.bounds.copy()
        bounds['phase_shift'] = (None, None)
        with pytest.raises(ValueError):
            HotRegion(bounds, self.values)
            
    def test_hotregion_breaks_if_phase_shift_bounds_infinite(self):
        bounds = self.bounds.copy()
        bounds['phase_shift'] = (-np.inf, np.inf)
        with pytest.raises(ValueError):
            HotRegion(bounds, self.values)
            
    def test_hotregion_breaks_if_phase_shift_bounds_more_than_one_cycle(self):
        bounds = self.bounds.copy()
        bounds['phase_shift'] = (0, 1.1)
        with pytest.raises(ValueError):
            HotRegion(bounds, self.values)
            
    def test_hotregion_produces_default_parameters(self):
        self.HotRegion.get_param('phase_shift')
        self.HotRegion.get_param('super_colatitude')
        self.HotRegion.get_param('super_radius')
        self.HotRegion.get_param('super_temperature')
        self.HotRegion.get_param('cede_colatitude')
        self.HotRegion.get_param('cede_radius')
        self.HotRegion.get_param('cede_azimuth')
        self.HotRegion.get_param('omit_colatitude')
        self.HotRegion.get_param('omit_radius')
        self.HotRegion.get_param('omit_azimuth')
        
        
    def test_hotregion_handles_cede(self):
        bounds = self.bounds.copy()
        bounds['cede_temperature'] = (5.1, 6.8)
        HotRegionCede = HotRegion(bounds, self.values, cede = True)
        HotRegionCede.get_param('cede_temperature')
        
    @pytest.mark.parametrize("inplace", [True, False])
    def test_hotregion_sets_symmetry(self, inplace):
        self.HotRegion.symmetry = inplace
        
    @pytest.mark.parametrize("declaration, split, expected_error", [
    (True, False, None),                 # No error, symmetry=True, split=False
    (False, False, None),                # No error, symmetry=False, split=False
    (False, True, "TypeError"),          # Error, split=True, symmetry=False
    ("invalid", False, "TypeError"),     # Error, invalid type for symmetry
])
    def test_symmetry(self, declaration, split, expected_error):
        hot_region = HotRegion(self.bounds, self.values, split=split)
    
        if expected_error:
            with pytest.raises(eval(expected_error)):
                hot_region.symmetry = declaration
        else:
            hot_region.symmetry = declaration
            assert hot_region._integrator is not None
            assert hot_region._integratorIQU is not None
