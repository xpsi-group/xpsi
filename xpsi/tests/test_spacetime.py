from __future__ import division, print_function

from xpsi.Spacetime import Spacetime
import pytest

class TestSpacetime(object):

    def setup_class(cls):
        cls.frequency = 300.0
        cls.mass = 1.4
        cls.radius = 10.0
        cls.distance = 3.0
        cls.cos_inclination = 0.6

        cls.frequency_bounds = (0.1, 750.0)
        cls.mass_bounds = (0.01, 2.5)
        cls.radius_bounds = (5.0, 15.0)
        cls.distance_bounds = (0.05, 25.0)
        cls.cos_inclination_bounds = (-0.9, 0.9)

        cls.bounds = {"frequency": cls.frequency_bounds,
                      "mass": cls.mass_bounds,
                      "radius": cls.radius_bounds,
                      "distance": cls.distance_bounds,
                      "cos_inclination": cls.cos_inclination_bounds}

        cls.values = {"frequency": cls.frequency,
                      "mass": cls.mass,
                      "radius": cls.radius,
                      "distance": cls.distance,
                      "cos_inclination": cls.cos_inclination}

    def test_spacetime_class_initializes(self):
        st = Spacetime(self.bounds, self.values)

    def test_spacetime_breaks_with_values_missin(self):
        with pytest.raises(TypeError):
            st = Spacetime(self.values)

    def test_spacetime_breaks_with_bounds_missing(self):
        with pytest.raises(TypeError):
            st = Spacetime(self.bounds)

    def test_spacetime_breaks_with_both_missing(self):
        with pytest.raises(TypeError):
            st = Spacetime()

    def test_spacetime_breaks_if_inputs_switched(self):
       with pytest.raises(TypeError):
           st = Spacetime(self.values, self.bounds)

    def test_failure_when_parameter_out_of_bounds(self):
        mass = 3.5
        values = {"frequency": self.frequency, 
                  "mass": mass, 
                  "radius": self.radius, 
                  "distance": self.distance,
                  "cos_inclination": self.cos_inclination}
        with pytest.raises(TypeError):
            st = Spacetime(values, self.bounds)

    def test_failure_when_parameter_out_of_strict_bounds(self):
        mass_bounds = (0.001, 3.5)
        bounds = {"frequency": self.frequency_bounds,
                  "mass": mass_bounds,
                  "radius": self.radius_bounds,
                  "distance": self.distance_bounds,
                  "cos_inclination": self.cos_inclination_bounds}
        with pytest.raises(ValueError):
            st = Spacetime(bounds, self.values)

    def test_fails_when_input_is_not_dict(self):
        bounds = (self.frequency_bounds, self.mass_bounds,
                  self.radius_bounds, self.distance_bounds,
                  self.cos_inclination_bounds)

        with pytest.raises(TypeError):
            st = Spacetime(bounds, self.values)
