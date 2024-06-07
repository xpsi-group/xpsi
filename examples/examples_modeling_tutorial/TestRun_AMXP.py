'''
Test script with the polarized 3+2 numerical atmosphere applied to a one-spot model.

Prequisities:
Before running the script, add the atmosphere data to the model_data subdirectory:
Bobrikova_compton_slab_I.npz and Bobrikova_compton_slab_Q.npz. See the
example script in xpsi/examples/produce_atmos_lookuptable for producing these files
from those provided in https://github.com/AnnaBobrikova/ComptonSlabTables.
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

this_directory = os.path.dirname(os.path.abspath(__file__))

from xpsi.global_imports import  _keV, _k_B
k_B_over_keV = _k_B / _keV


from modules.helper_functions import get_T_in_log10_Kelvin



bounds = dict(distance = (0.1, 1.0),                     # (Earth) distance
                mass = (1.0, 3.0),                       # mass
                radius = (3.0 * gravradius(1.0), 16.0),  # equatorial radius
                cos_inclination = (0.0, 1.0))      # (Earth) inclination to rotation axis

spacetime = xpsi.Spacetime(bounds=bounds, values=dict(frequency=400.9752075))

from xpsi.Parameter import Parameter
from modules.CustomHotRegion_Accreting import CustomHotRegion_Accreting


bounds = dict(super_colatitude = (None, None),
              super_radius = (None, None),
              phase_shift = (0.0, 0.1),
              super_tbb = (0.001, 0.003),
              super_tau = (0.5, 3.5),
              super_te = (40.0, 200.0))

primary = CustomHotRegion_Accreting(bounds=bounds,
                                    values={},
                                    symmetry=True,
                                    omit=False,
                                    cede=False,
                                    concentric=False,
                                    sqrt_num_cells=32, #100
                                    min_sqrt_num_cells=10,
                                    max_sqrt_num_cells=64, #100
                                    num_leaves=100,
                                    num_rays=200,
                                    split=True,
                                    image_order_limit=3,
                                    prefix='p')

from xpsi import HotRegions
hot = HotRegions((primary,))

use_elsewhere = False #True

if use_elsewhere:
    bounds=dict(elsewhere_temperature = (3.0, 7.5))

    elsewhere = xpsi.Elsewhere(bounds=bounds,
                   values={},
                   sqrt_num_cells=512,
                   num_rays=512,
                   atm_ext="BB")
else:
    elsewhere = None

from xpsi.HotRegion import HotRegion
from xpsi.Elsewhere import Elsewhere
from xpsi.Everywhere import Everywhere


    
    

# trying to use the new disk
from modules.Disk import Disk, k_disk_derive
k_disk = k_disk_derive()
T_in = get_T_in_log10_Kelvin(0.25) #(0.29)
R_in = 30.0 #55.0
values = {'T_in':T_in,'R_in':R_in,'K_disk': k_disk}
disk = Disk(bounds={}, values=values)

from modules.CustomPhotosphere import CustomPhotosphere_NumA5


bounds = dict(spin_axis_position_angle = (None, None))
photosphere = CustomPhotosphere_NumA5(hot = hot, elsewhere = elsewhere, stokes=True, disk=disk, bounds=bounds,
                                values=dict(mode_frequency = spacetime['frequency']))

photosphere.hot_atmosphere = this_directory+'/model_data/Bobrikova_compton_slab_I.npz'
photosphere.hot_atmosphere_Q = this_directory+'/model_data/Bobrikova_compton_slab_Q.npz'

photosphere['mode_frequency'] == spacetime['frequency']

star = xpsi.Star(spacetime = spacetime, photospheres = photosphere)

k_disk.star = star
k_disk.disk = disk
# disk_flux = disk(energies)

# print("Parameters of the star:")
# print(star.params, len(star.params))

# SAX J1808-like 
mass = 1.4
radius = 11.0
distance = 2.7
inclination = 40.
cos_i = math.cos(inclination*math.pi/180.0)
chi0 = 0.0

# Hotspot
phase_shift =  0.226365126031355196E+00
super_colatitude = 0.175993450466385537E+00
super_radius = 0.156951249537834525E+01

# Compton slab model parameters
tbb=0.103616176435110115E-02
te=0.729440224892133244E+02
tau=0.153014380768402769E+01

#Tbb = 1 keV <=> tbb = 0.002 (roughly)
#Te = 50 keV <=>  te = 100 (roughly)

p = [mass, #grav mass
      radius, #coordinate equatorial radius
      distance, # earth distance kpc
      cos_i, #cosine of earth inclination
      chi0, #spin axis position angle
      phase_shift, #phase of hotregion
      super_colatitude, #colatitude of centre of superseding region
      super_radius,  #angular radius superceding region
      tbb,
      te,
      tau
      ]

# print(len(p))

# elsewhere
elsewhere_T_keV = 0.4 #  keV
elsewhere_T_log10_K = get_T_in_log10_Kelvin(elsewhere_T_keV)

if use_elsewhere:
    p.append(elsewhere_T_log10_K)

star(p)
star.update()
# start = time.time()

#To get the incident signal before interstellar absorption or operating with the telescope:
energies = np.logspace(np.log10(0.15), np.log10(12.0), 100, base=10.0)
photosphere.integrate(energies, threads=1) # the number of OpenMP threads to use

# end = time.time()
# print("Time spent in integration:",end - start)
#exit()

# print("Bolometric profiles for I, Q, and U:")
# print(repr(np.sum(photosphere.signal[0][0], axis=0)))
# print(repr(np.sum(photosphere.signalQ[0][0], axis=0)))
# print(repr(np.sum(photosphere.signalU[0][0], axis=0)))

from modules.CustomInterstellar import CustomInterstellar


StokesI = photosphere.signal[0][0]
StokesQ = photosphere.signalQ[0][0]
StokesU = photosphere.signalU[0][0]
print('stokesI',StokesI)




