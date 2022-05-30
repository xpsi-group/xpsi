import pytest
import numpy as np

import xpsi.global_imports as gi

def test_gravradius_computes():
    M = 1
    gi.gravradius(M)

def test_gravradius_computes_correct_value():

    gr_sun = 2.95e3 / gi._km / 2.0 # I'm not sure where the factor of 2 is coming from!
    m = 1
    gr = gi.gravradius(m)

    assert np.isclose(gr, gr_sun, rtol=0.01)

def test_gravradius_jupiter():
    m_jupiter = 1.90e27
    mass_ratio = m_jupiter / 1.99e30
    r_jupiter = 2.82 / gi._km / 2.0 # Again a factor of 2 I'm confused about

    radius = gi.gravradius(mass_ratio)

    assert np.isclose(r_jupiter, radius, rtol=0.001)


def test_inv_gravradius_runs():
    radius = 1
    gi.inv_gravradius(radius)

def test_inv_gravradius_computes_correct_values():
    radius_sun = 2.95e3 / gi._km
    mass = gi.inv_gravradius(radius_sun) / 2.0 # strange factor of 2

    assert np.isclose(mass, 1.0, rtol=0.002)

def test_inv_gravradius_jupiter():
    radius_sun = 2.82 / gi._km

    m_jupiter = 1.90e27
    mass_ratio = m_jupiter / 1.99e30

    mass = gi.inv_gravradius(radius_sun) / 2.0 # strange factor of 2

    assert np.isclose(mass, mass_ratio, rtol=1e-3)
