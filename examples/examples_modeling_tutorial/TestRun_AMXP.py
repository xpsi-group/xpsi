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

import xpsi

from xpsi.global_imports import _c, gravradius

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

num_leaves=30

primary = CustomHotRegion_Accreting(bounds=bounds,
                                    values={},
                                    symmetry=True,
                                    omit=False,
                                    cede=False,
                                    concentric=False,
                                    sqrt_num_cells=50, #100
                                    min_sqrt_num_cells=10,
                                    max_sqrt_num_cells=64, #100
                                    num_leaves=num_leaves,
                                    num_rays=200,
                                    split=False,
                                    atm_ext='Num5D',
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


from xpsi.ParameterSubspace import ParameterSubspace
from xpsi.Parameter import Parameter, Derive
from xpsi.global_imports import  _keV, _k_B, _h_keV
_c = 2.99792458E8
_c_cgs = _c*1E2
k_B_over_keV = _k_B / _keV
import numpy as np
from scipy.integrate import quad




class Disk(ParameterSubspace):
     
     def __init__(self, bounds=None, values=None, interstellar = None):

         doc = """
         Temperature at inner disk radius in log10 Kelvin.
         """
         inner_temperature = Parameter('T_in',
                                 strict_bounds = (3., 10.),
                                 bounds = bounds.get('T_in', None),
                                 doc = doc,
                                 symbol = r'$T_{in}$',
                                 value = values.get('T_in', None))

         doc = """
         Disk R_in in kilometers.
         """
         inner_radius = Parameter('R_in',
                                 strict_bounds = (0., 1e3),
                                 bounds = bounds.get('R_in', None),
                                 doc = doc,
                                 symbol = r'$R_{in}$',
                                 value = values.get('R_in', None))

         
         doc = """
         Disk normalisation cos_i*R_in^2/D^2 in (km / 10 kpc)^2.
         """
         background_normalisation = Parameter('K_disk',
                                 strict_bounds = (0., 1e8),
                                 bounds = bounds.get('K_disk', None),
                                 doc = doc,
                                 symbol = r'$K_{BB}$',
                                 value = values.get('K_disk', None))
         

         super(Disk, self).__init__(inner_temperature, inner_radius, background_normalisation)

     def __call__(self, energies):
              
         spectral_radiance = self.B_E 
         # distance_m = 3.08567758128e19*distance
         #Converting to photosphere.signal units
         self.disk_flux = self.get_f_disk(energies, spectral_radiance)/energies # *distance_m**2/energies
         return self.disk_flux
         
     def get_f_disk(self, energies, spectral_radiance, attenuate = False):
         """ Evaluate f_disk(energies).
         
         f_disk(E) = 4/3*pi * K_disk * l_disk(b_E/B_E, E)
         Ref: Mitsuda et al. 1984, Makishima et al. 1986
         But note that Mitsuda et al. 1984 has an extra factor here because they
         don't have it in the definition for b_E/B_E.
         
         parameters
         energies[keV]
         spectral_radiance can be: b_E or B_E
         attenuate determines whether to apply interstellar medium attenuation.
         
         returns
         f_disk [photons/s/cm^2/keV] or [keV/s/cm^2/keV]
         
         """
         

         T_in = self['T_in']
         K_disk = self['K_disk']

         # KbT in keV
         T_in_keV = k_B_over_keV * pow(10.0, T_in)
         
         T_out_keV = T_in_keV*1e-1
         
         epsrel = 1e-4

         f_disk_array = np.array([]) #photons/s/cm^2/sr or keV/s/cm^2/sr 
         for energy in energies:
             f_disk_value = self.l_disk(energy, T_in_keV, T_out_keV, spectral_radiance, epsrel) 
             f_disk_array=np.append(f_disk_array,f_disk_value)
         
         # K_disk is cos_i*R_in^2/D^2 in (km / 10 kpc)^2.
         # (1 km / 10 kpc)^2 = 1.0502650e-35 [ cm/cm ]
         
         f_disk_array *=K_disk*4*np.pi/3*1.0502650e-35 # photons/s/cm^2/energy_bin
         
         #print("f_disk_array:",f_disk_array)
         #exit()
         
             # Apply Interstellar if not None
         if attenuate:
             if self.interstellar is not None:
                 self.interstellar(energies, f_disk_array) # bkg is overwritten here
         
         return f_disk_array

     def b_E(self, E, T):
         '''
         photon radiance of a blackbody

         parameters:
             E in keV
             T in keV

         returns:
             b_E in photons/s/keV/cm^2/sr 
         '''

         b = 2*E**2/(_h_keV**3*_c_cgs**2)/(np.exp(E/T)-1)
         return b
         
         
     def B_E(self, E, T):
         '''
         Energy radiance of a blackbody.

         parameters:
             E in keV
             T in keV

         returns:
             B_E in keV/s/keV/cm^2/sr (you will integrate over keV)
         '''
         
         B = 2*E**3/(_h_keV**3*_c_cgs**2)/(np.exp(E/T)-1)
         return B


     def l_disk_integrand(self, T, E, T_in, spectral_radiance):
         '''
         parameters:
             T, T_in in keV
             E in keV

         returns:
             integrand in spectral radiance units/keV. This integrand will 
             be integrated over keV.
         '''

         # print('T: ', T)
         # print('E:', E)
         # print('T_in: ', T_in)
         # print('(T/T_in)**(-11/3)', (T/T_in)**(-11/3))
         # print('spectral_radiance(E, T)/T_in', spectral_radiance(E, T)/T_in)

         integrand = (T/T_in)**(-11/3)*spectral_radiance(E, T)/T_in
         return integrand
     
     def l_disk(self, E, T_in, T_out, spectral_radiance, epsrel):
         '''
         parameters:
             T, T_in in keV
             E in keV

         returns:
             disk luminosity [spectral radiance units]. 
         '''

         disk_luminosity,_= quad(self.l_disk_integrand, T_out, T_in, args=(E, T_in, spectral_radiance), epsrel=epsrel)
         return disk_luminosity
     
from xpsi import Derive
   
def get_k_disk(cos_i, r_in, distance):
   """
   This function calculates the k-disk value for a given set of input parameters.

   Args:
       cos_i: The cosine inclination angle of the disk, can be a scalar or a tuple.
       r_in: The inner radius of the disk in kilometers, can be a scalar or a tuple.
       distance: The distance to the disk in kiloparsecs, can be a scalar or a tuple.

   Returns:
       A tuple containing the k-disk values for each element in the input parameters.

   Raises:
       ValueError: If the input tuples have different lengths.
   """

   if isinstance(cos_i, tuple) and isinstance(r_in, tuple) and isinstance(distance, tuple):
     if len(cos_i) != len(r_in) or len(cos_i) != len(distance):
       raise ValueError("Input tuples must have the same length.")
     # Use a loop instead of recursion
     k_disk_values = []
     for c, r, d in zip(cos_i, r_in, distance):
       k_disk_values.append(c * (r / (d / 10))**2)
     return tuple(k_disk_values)
   else:
     # return cos_i * (r_in / (distance / 10))**2
     
     # scaling k_disk further to match signal units
     distance_m = 3.08567758128e19*distance
     
     return cos_i * (r_in / (distance / 10))**2 * distance_m**2

class k_disk_derive(Derive):
     def __init__(self):
         pass

     def __call__(self, boundto, caller = None):
         # ref is a reference to another hot region object
         return get_k_disk(self.star['cos_inclination'], self.disk['R_in'], self.star['distance'])
   

k_disk = k_disk_derive()
T_in = get_T_in_log10_Kelvin(0.16845756373108872) #(0.29)
R_in = 0.308122224729265000E+02 #55.0
values = {'T_in':T_in,'R_in':R_in,'K_disk': k_disk}
disk = Disk(bounds={}, values=values)

from modules.CustomPhotosphere import CustomPhotosphere_NumA5

stokes = False
bounds = dict(spin_axis_position_angle = (None, None))
photosphere = CustomPhotosphere_NumA5(hot = hot, elsewhere = elsewhere, stokes=stokes, disk=None, bounds=bounds,
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
inclination = 80.
cos_i = math.cos(inclination*math.pi/180.0)
if stokes:
    chi0 = 0.0

# Hotspot
phase_shift = 0.226365126031355196E+00
super_colatitude = 0.175993450466385537E+00
super_radius = 30.*math.pi/180

# Compton slab model parameters
tbb=0.0025
te=100.
tau=2.0

#Tbb = 1 keV <=> tbb = 0.002 (roughly)
#Te = 50 keV <=>  te = 100 (roughly)

if stokes:
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
elif not stokes:
    p = [mass, #grav mass
          radius, #coordinate equatorial radius
          distance, # earth distance kpc
          cos_i, #cosine of earth inclination
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
start = time.time()

#To get the incident signal before interstellar absorption or operating with the telescope:
energies = np.logspace(np.log10(0.15), np.log10(12.0), 40, base=10.0)
multiple_times = 100
for i in range(multiple_times):
    photosphere.integrate(energies, threads=1) # the number of OpenMP threads to use

end = time.time()
print("Time spent in integration:",(end - start)/multiple_times)
#exit()

# print("Bolometric profiles for I, Q, and U:")
# print(repr(np.sum(photosphere.signal[0][0], axis=0)))
# print(repr(np.sum(photosphere.signalQ[0][0], axis=0)))
# print(repr(np.sum(photosphere.signalU[0][0], axis=0)))

StokesI = photosphere.signal[0][0]
# StokesQ = photosphere.signalQ[0][0]
# StokesU = photosphere.signalU[0][0]
print('stokesI',StokesI)

phases = np.linspace(0,1,num_leaves)
from modules.helper_functions import CustomAxes
fig, ax = plt.subplots()
profile = CustomAxes.plot_2D_signal(ax, [StokesI], phases, [0,0], energies, 'energies')
fig.colorbar(profile, ax=ax)


# plt.plot(energies[0:50],np.sum(StokesI,axis=1)[0:50])
# print(energies[0:130])
# print(np.sum(StokesI,axis=1)[0:130])
# plt.ylabel('Flux [?]')
# plt.xlabel('Energy [keV]')
# plt.ylim(0.0,8.0e31)


