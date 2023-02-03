import numpy as np
from xpsi.HotRegion import HotRegion
from xpsi.HotRegions import HotRegions
import pytest

class TestHotRegions(object):
    
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
        
        cls.primary = HotRegion(cls.bounds, cls.values, prefix='p')
        cls.secondary = HotRegion(cls.bounds, cls.values, prefix='s')
        
    def test_hotregions_class_initializes(self):
        hrs = HotRegions((self.primary, self.secondary))
        
    def test_hotregions_breaks_with_one_hotregion(self):
        with pytest.raises(ValueError):
            hrs = HotRegions((self.primary))

