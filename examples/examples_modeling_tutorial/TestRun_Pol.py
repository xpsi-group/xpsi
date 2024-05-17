'''
Test script with the polarized burst atmosphere extension.
'''

import os
import numpy as np
import math
import time

from matplotlib import pyplot as plt
from matplotlib import rcParams
from matplotlib.ticker import MultipleLocator, AutoLocator, AutoMinorLocator
from matplotlib import gridspec
from matplotlib import cm

import xpsi

from xpsi.global_imports import _c, _G, _dpr, gravradius, _csq, _km, _2pi

np.random.seed(xpsi._rank+10)

print('Rank reporting: %d' % xpsi._rank)

bounds = dict(distance = (0.1, 1.0),                     # (Earth) distance
                mass = (1.0, 3.0),                       # mass
                radius = (3.0 * gravradius(1.0), 16.0),  # equatorial radius
                cos_inclination = (0.0, 1.0))      # (Earth) inclination to rotation axis

spacetime = xpsi.Spacetime(bounds=bounds, values=dict(frequency=400.9752075))

bounds = dict(super_colatitude = (None, None),
              super_radius = (None, None),
              phase_shift = (0.0, 0.1),
              super_temperature = (5.1, 6.8))

primary = xpsi.HotRegion(bounds=bounds,
                        values={},
                        symmetry=True,
                        omit=False,
                        cede=False,
                        concentric=False,
                        sqrt_num_cells=32,
                        min_sqrt_num_cells=10,
                        max_sqrt_num_cells=64,
                        num_leaves=100,
                        num_rays=200,
                        atm_ext="Pol_BB_Burst",
                        image_order_limit=3,
                        prefix='p')

bounds2 = dict(super_colatitude = (None, None),
                        super_radius = (None, None),
                        phase_shift = (0.0, 0.1),
                        super_temperature = (5.1, 6.8),
                        cede_colatitude = (None, None),
                        cede_radius = (None, None),
                        cede_azimuth = (None, None),
                        cede_temperature = (5.1, 6.8))

secondary = xpsi.HotRegion(bounds=bounds2, # can otherwise use same bounds
                            values={},
                            symmetry=True,
                            omit=False,
                            cede=True,
                            concentric=False,
                            sqrt_num_cells=32,
                            min_sqrt_num_cells=10,
                            max_sqrt_num_cells=100,
                            num_leaves=100,
                            num_rays=200,
                            do_fast=False,
                            atm_ext="Pol_BB_Burst",
                            is_antiphased=True,
                            image_order_limit=3,
                            prefix='s')


from xpsi import HotRegions
hot = HotRegions((primary, secondary))
h = hot.objects[0]
hot['p__super_temperature'] = 6.0 # equivalent to ``primary['super_temperature'] = 6.0``

use_elsewhere = False

if use_elsewhere:
    bounds=dict(elsewhere_temperature = (5.2, 6.5))

    elsewhere = xpsi.Elsewhere(bounds=bounds,
                   values={},
                   sqrt_num_cells=512,
                   num_rays=512,
                   atm_ext="BB")
else:
    elsewhere = None

bounds = dict(spin_axis_position_angle = (None, None))
photosphere = xpsi.Photosphere(hot = hot, elsewhere = elsewhere, stokes=True, bounds=bounds,
                                values=dict(mode_frequency = spacetime['frequency']))


photosphere['mode_frequency'] == spacetime['frequency']

star = xpsi.Star(spacetime = spacetime, photospheres = photosphere)

print("Parameters of the star:")
print(star.params, len(star.params))

p = [1.0368513939430604,
     6.087862992320039,
     0.26870812456714116,
     0.39140510783272897,
     0.0,
     0.04346870860640872,
     0.8002010406881243,
     1.1165398710637626,
     5.865655057483478,
     0.07360477761463673,
     2.4602238829718432,
     0.4277092192054918,
     6.0,
     2.3, #cede_colatitude
     0.5, #cede_radius
     0.1, #cede_azimuth
     6.1] #cede_temperature

if use_elsewhere:
    p.append(5.5)

star(p)
star.update()

#start = time.time()

#To get the incident signal before interstellar absorption or operating with the telescope:
#energies = np.logspace(-1.0, np.log10(3.0), 128, base=10.0)
energies = np.logspace(-1.0, np.log10(12.0), 400, base=10.0)
photosphere.integrate(energies, threads=1) # the number of OpenMP threads to use

#end = time.time()
#print("Time spent in integration:",end - start)
#exit()

print("Bolometric profiles for I, Q, and U:")
print("1st spot:")
print(repr(np.sum(photosphere.signal[0][0], axis=0)))
print(repr(np.sum(photosphere.signalQ[0][0], axis=0)))
print(repr(np.sum(photosphere.signalU[0][0], axis=0)))
print("2nd spot:")
print(repr(np.sum(photosphere.signal[1][0], axis=0)))
print(repr(np.sum(photosphere.signalQ[1][0], axis=0)))
print(repr(np.sum(photosphere.signalU[1][0], axis=0)))
print()
print("2nd spot, ceding):")
print(np.sum(photosphere.signal[1][1], axis=0))
print(np.sum(photosphere.signalQ[1][1], axis=0))
print(np.sum(photosphere.signalU[1][1], axis=0))
print()

def get_photosphere_stokes_1spot():
    #Return signal from the 1st spot to ixpeobssim
    return hot.phases_in_cycles[0], energies, photosphere.signal[0][0], photosphere.signalQ[0][0], photosphere.signalU[0][0]
