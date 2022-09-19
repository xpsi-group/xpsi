from __future__ import division, print_function

from .Spacetime import Spacetime

def TestSpacetime(object):

    def setup_class(self):
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

        cls.bounds = (self.frequency_bounds,
                      self.mass_bounds,
                      self.radius_bounds,
                      self.distance_bounds,
                      self.cos_inclination_bounds)

        cls.values = (self.frequency,
                      self.mass,
                      self.radius,
                      self.distance,
                      self.cos_inclination)

    def test_spacetime_class_initializes(self):
        st = Spacetime(self.bounds, self.values)

    def test_spacetime_breaks_with_values_missin(self):
        st = Spacetime(self.values)

    def test_spacetime_breaks_with_bounds_missing(self):
        st = Spacetime(self.bounds)

    def test_spacetime_breaks_with_both_missing(self):
        st = Spacetime()

    def test_spacetime_breaks_if_inputs_switched(self):
       st = Spacetime(self.values, self.bounds)

    def test_failure_when_parameter_out_of_bounds(self):
        mass = 3.5
        values = (self.frequency, mass, self.radius, self.distance,
                  self.cos_inclination)

        st = Spacetime(values, self.bounds)

    def test_failure_when_parameter_out_of_strict_bounds(self):
        mass_bounds = (0.001, 3.5)
        bounds = (self.frequency_bounds,
                  mass_bounds,
                  self.radius_bounds,
                  self.distance_bounds,
                  self.cos_inclination_bounds)

        st = Spacetime(bounds, self.values)

    def test_failure_when_input_none(self):
        mass = None
        values = (self.frequency, mass, self.radius, self.distance, 
                  self.cos_inclination)
        st = Spacetime(self.bounds, self.values)


