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
        
    def test_hotregion_class_initializes(self):
        HotRegion(self.bounds, self.values)
    
    def test_hotregion_breaks_with_bounds_missing(self):
        with pytest.raises(TypeError):
            HotRegion(self.values)
            
    def test_hotregion_breaks_with_values_missing(self):
        with pytest.raises(TypeError):
            HotRegion(self.bounds)
            
    def test_hotregion_breaks_if_phase_shift_bounds_unspecified(self):
        self.bounds['phase_shift'] = (None, None)
        with pytest.raises(ValueError):
            HotRegion(self.bounds, self.values)
            
    def test_hotregion_breaks_if_phase_shift_bounds_infinite(self):
        self.bounds['phase_shift'] = (-np.inf, np.inf)
        with pytest.raises(ValueError):
            HotRegion(self.bounds, self.values)
            
    def test_hotregion_breaks_if_phase_shift_bounds_more_than_one_cycle(self):
        self.bounds['phase_shift'] = (0, 1.1)
        with pytest.raises(ValueError):
            HotRegion(self.bounds, self.values)
            
    