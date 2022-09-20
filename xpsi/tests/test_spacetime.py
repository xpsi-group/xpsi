from __future__ import division, print_function

import numpy as np
from xpsi.Spacetime import Spacetime
import pytest

class TestSpacetimeInit(object):

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


class TestSpacetimeProperties(object):

    def setup_class(cls):
        cls.frequency = 100.0
        cls.mass = 1.0
        cls.radius = 10.0
        cls.distance = 3.0
        cls.cos_inclination = -1.0

        cls.frequency_bounds = (0.1, 750.0)
        cls.mass_bounds = (0.01, 2.5)
        cls.radius_bounds = (5.0, 15.0)
        cls.distance_bounds = (0.05, 25.0)
        cls.cos_inclination_bounds = (-1.0, 0.9)

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

        cls.st = Spacetime(cls.bounds, cls.values)

        cls.rg_sun = 2.95*10**3
        cls.msun =  1.989*10**30

    def test_M_works(self):
        self.st.M

    def test_M_is_correct(self):
        assert np.isclose(self.msun, self.st.M, rtol=0.001)

    def test_rg_runs(self):
        self.st.r_g

    def test_rg_is_correct(self):
        # XPSI definition of R_g includes a factor of
        # 2 elsewhere, hence the division by 2 below
        rg_sun = self.rg_sun / 2
        assert np.isclose(rg_sun, self.st.r_g, rtol=0.002)

    def test_schwarzschild_radius_computes(self):
        self.st.r_s

    def test_rs_is_correct(self):
        assert np.isclose(self.rg_sun, self.st.r_s , rtol=0.002)

    def test_r_computes(self):
        self.st.R

    def test_r_is_correct(self):
        assert self.st.R == self.radius * 1000

    def test_rrs_computes(self):
        self.st.R_r_s

    def test_rrs_is_correct(self):
        ratio = self.radius * 1000.0 / self.rg_sun
        assert np.isclose(self.st.R_r_s, ratio, rtol=0.002)

    def test_f_computes(self):
        self.st.f

    def test_f_is_correct(self):
        assert self.st.f == self.frequency

    def test_omega_computes(self):
        self.st.Omega

    def test_omega_is_correct(self):
        assert self.st.Omega == 2 * np.pi * self.frequency

    def test_i_computes(self):
        self.st.i

    def test_i_is_correct(self):
        assert self.st.i == np.pi

    def test_d_computes(self):
        self.st.d

    def test_d_is_correct(self):
        dist = self.distance * 3.086 * 10**16 * 1000.0
        assert np.isclose(self.st.d, dist, rtol=0.001)

    def test_dsq_computes(self):
        self.st.d_sq

    def test_dsq_correct(self):
        dist = (self.distance * 3.086 * 10**16 * 1000.0) ** 2.0
        assert np.isclose(self.st.d_sq, dist, rtol=0.001) 

    def test_zeta_gets_calculated(self):
        zeta = (self.rg_sun / 2.) / (self.radius * 1000.0)
        assert np.isclose(zeta, self.st.zeta, rtol=0.002)

    def test_zeta_gets_returned_upon_calculation(self):
        st = Spacetime(self.bounds, self.values)
        zeta_val = 42
        st._zeta = zeta_val

        assert st.zeta == zeta_val

    def test_epsilon_gets_calculated(self):
        omega_sq = (2. * np.pi * self.frequency) ** 2.0
        r_cubed = (self.radius * 1000) ** 3.0
        grav_const = 6.6743e-11
        gm = grav_const  * self.msun
        eps = omega_sq * r_cubed / gm

        assert np.isclose(eps, self.st.epsilon, rtol=0.001) 

    def test_epsilon_gets_returned_upon_calculation(self):
        st = Spacetime(self.bounds, self.values)
        epsilon_val = 42
        st._epsilon = epsilon_val

        assert st.epsilon == epsilon_val

    def test_a_setter(self):
        st = Spacetime(self.bounds, self.values)
        aval = 42
        st.a = aval

        assert st._a == aval

    def test_a_deleter(self):
        st = Spacetime(self.bounds, self.values)
        aval = 42
        st.a = aval

        assert st._a == aval
        del st.a
        with pytest.raises(AttributeError):
            st._a

    def test_q_setter(self):
        st = Spacetime(self.bounds, self.values)
        qval = 42
        st.q = qval

        assert st._q == qval

    def test_q_deleter(self):
        st = Spacetime(self.bounds, self.values)
        qval = 42
        st.q = qval
        
        assert st._q == qval
        del st.q
        with pytest.raises(AttributeError):
            st._q
 
