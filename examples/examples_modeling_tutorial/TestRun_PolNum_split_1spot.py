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

bounds = dict(distance = (0.1, 1.0),                     # (Earth) distance
                mass = (1.0, 3.0),                       # mass
                radius = (3.0 * gravradius(1.0), 16.0),  # equatorial radius
                cos_inclination = (0.0, 1.0))      # (Earth) inclination to rotation axis

spacetime = xpsi.Spacetime(bounds=bounds, values=dict(frequency=400.9752075))

from xpsi.Parameter import Parameter
class CustomHotRegion_Accreting(xpsi.HotRegion):
    """Custom implementation of HotRegion. Accreting Atmosphere model by 
    Anna Bobrikova. The parameters are ordered I(E < mu < tau < tbb < te).
    
    E is energy.
    mu is cos of zenith angle.
    tau is the optical depth of the comptom slab.
    tbb is the black body temperature.
    te is temperature of the electron gas.
    """

    required_names = ['super_colatitude',
                      'super_radius',
                      'phase_shift',
                      'super_tbb',
                      'super_te',
                      'super_tau']
    optional_names = ['omit_colatitude',
                      'omit_radius',
                      'omit_azimuth',
                      'cede_colatitude',
                      'cede_radius',
                      'cede_azimuth',
                      'cede_tbb',
                      'cede_te',
                      'cede_tau']
    
    def __init__(self,
            bounds,
            values,
            symmetry = 'azimuthal_invariance',
            interpolator = 'split',
            omit = False,
            cede = False,
            concentric = False,
            sqrt_num_cells = 32,
            min_sqrt_num_cells = 10,
            max_sqrt_num_cells = 80,
            num_rays = 200,
            num_leaves = 64,
            num_phases = None,
            phases = None,
            do_fast = False,
            fast_sqrt_num_cells = 16,
            fast_min_sqrt_num_cells = 4,
            fast_max_sqrt_num_cells = 16,
            fast_num_rays = 100,
            fast_num_leaves = 32,
            fast_num_phases = None,
            fast_phases = None,
            is_antiphased = False,
            custom = None,
            image_order_limit = None,
            **kwargs
            ):

        doc = """
        tbb
        """
        super_tbb = Parameter('super_tbb',
                    strict_bounds = (0.001, 0.003), #tbb = Tbb(keV)/511keV
                    bounds = bounds.get('super_tbb', None),
                    doc = doc,
                    symbol = r'tbb',
                    value = values.get('super_tbb', None))
        doc = """
        te
        """
        super_te = Parameter('super_te',
                    strict_bounds = (40., 200.), #te = Te(keV)*1000/511keV
                    bounds = bounds.get('super_te', None),
                    doc = doc,
                    symbol = r'te',
                    value = values.get('super_te', None))
        
        doc = """
        tau
        """
        super_tau = Parameter('super_tau',
                    strict_bounds = (0.5, 3.5),
                    bounds = bounds.get('super_tau', None),
                    doc = doc,
                    symbol = r'tau',
                    value = values.get('super_tau', None))
        

        custom = [super_tbb, super_te, super_tau]

        if cede:
            doc = """
            cede_tbb
            """        
            cede_tbb = Parameter('cede_tbb',
            strict_bounds = (0.001, 0.003),
            bounds = bounds.get('cede_tbb', None),
            doc = doc,
            symbol = r'cede_tbb',
            value = values.get('cede_tbb', None))

            doc = """
            cede_te
            """
            cede_te = Parameter('cede_te',
                        strict_bounds = (40., 200.),
                        bounds = bounds.get('cede_te', None),
                        doc = doc,
                        symbol = r'cede_te',
                        value = values.get('cede_te', None))
            
            doc = """
            cede_tau
            """
            cede_tau = Parameter('cede_tau',
                        strict_bounds = (0.5, 3.5),
                        bounds = bounds.get('cede_tau', None),
                        doc = doc,
                        symbol = r'cede_tau',
                        value = values.get('cede_tau', None))
            
            custom += [cede_tbb,cede_te,cede_tau]

        super(CustomHotRegion_Accreting, self).__init__(
                bounds,
                values,
                symmetry = symmetry,
                interpolator = interpolator,
                omit = omit,
                cede = cede,
                concentric = concentric,
                sqrt_num_cells = sqrt_num_cells,
                min_sqrt_num_cells = min_sqrt_num_cells,
                max_sqrt_num_cells = max_sqrt_num_cells,
                num_rays = num_rays,
                num_leaves = num_leaves,
                num_phases = num_phases,
                phases = phases,
                do_fast = do_fast,
                fast_sqrt_num_cells = fast_sqrt_num_cells,
                fast_min_sqrt_num_cells = fast_min_sqrt_num_cells,
                fast_max_sqrt_num_cells = fast_max_sqrt_num_cells,
                fast_num_rays = fast_num_rays,
                fast_num_leaves = fast_num_leaves,
                fast_num_phases = fast_num_phases,
                fast_phases = fast_phases,
                is_antiphased = is_antiphased,
                custom = custom,
                image_order_limit = image_order_limit,
                **kwargs
                )

    def _HotRegion__compute_cellParamVecs(self):
        self._super_radiates = np.greater(self._super_cellArea, 0.0).astype(np.int32)
        self._super_cellParamVecs = np.ones((self._super_radiates.shape[0],
                                      self._super_radiates.shape[1],
                                      3),
                                     dtype=np.double)

        self._super_cellParamVecs[...,0] *= self['super_te']
        self._super_cellParamVecs[...,1] *= self['super_tbb']
        self._super_cellParamVecs[...,2] *= self['super_tau']

        try:
            self._cede_radiates = np.greater(self._cede_cellArea, 0.0).astype(np.int32)
        except AttributeError:
            pass
        else:
            self._cede_cellParamVecs = np.ones((self._cede_radiates.shape[0],
                                                 self._cede_radiates.shape[1],
                                                 3), dtype=np.double)

            self._cede_cellParamVecs[...,0] *= self['cede_te']
            self._cede_cellParamVecs[...,1] *= self['cede_tbb']
            self._cede_cellParamVecs[...,2] *= self['cede_tau']



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


class CustomPhotosphere_NumA5(xpsi.Photosphere):
    """ A photosphere extension to preload the numerical 5D accretion atmosphere. """

    @xpsi.Photosphere.hot_atmosphere.setter
    def hot_atmosphere(self, path):
        with np.load(path, allow_pickle=True) as data_dictionary:
            NSX = data_dictionary['NSX.npy']
            size_reorderme = data_dictionary['size.npy']

        size = [size_reorderme[3], size_reorderme[4], size_reorderme[2], size_reorderme[1], size_reorderme[0]]

        Energy = np.ascontiguousarray(NSX[0:size[0],0])
        cos_zenith = np.ascontiguousarray([NSX[i*size[0],1] for i in range(size[1])])
        tau = np.ascontiguousarray([NSX[i*size[0]*size[1],2] for i in range(size[2])])
        t_bb = np.ascontiguousarray([NSX[i*size[0]*size[1]*size[2],3] for i in range(size[3])])
        t_e = np.ascontiguousarray([NSX[i*size[0]*size[1]*size[2]*size[3],4] for i in range(size[4])])
        intensities = np.ascontiguousarray(NSX[:,5])

        self._hot_atmosphere = (t_e, t_bb, tau, cos_zenith, Energy, intensities)

    @xpsi.Photosphere.hot_atmosphere_Q.setter
    def hot_atmosphere_Q(self, path):
        with np.load(path, allow_pickle=True) as data_dictionary:
            NSX = data_dictionary['NSX.npy']
            size_reorderme = data_dictionary['size.npy']

        size = [size_reorderme[3], size_reorderme[4], size_reorderme[2], size_reorderme[1], size_reorderme[0]]

        Energy = np.ascontiguousarray(NSX[0:size[0],0])
        cos_zenith = np.ascontiguousarray([NSX[i*size[0],1] for i in range(size[1])])
        tau = np.ascontiguousarray([NSX[i*size[0]*size[1],2] for i in range(size[2])])
        t_bb = np.ascontiguousarray([NSX[i*size[0]*size[1]*size[2],3] for i in range(size[3])])
        t_e = np.ascontiguousarray([NSX[i*size[0]*size[1]*size[2]*size[3],4] for i in range(size[4])])
        intensities = np.ascontiguousarray(NSX[:,5])

        self._hot_atmosphere_Q = (t_e, t_bb, tau, cos_zenith, Energy, intensities)

photosphere = CustomPhotosphere_NumA5(hot = hot, elsewhere = elsewhere, stokes=True,
                                values=dict(mode_frequency = spacetime['frequency']))

photosphere.hot_atmosphere = this_directory+'/model_data/Bobrikova_compton_slab_I.npz'
photosphere.hot_atmosphere_Q = this_directory+'/model_data/Bobrikova_compton_slab_Q.npz'

photosphere['mode_frequency'] == spacetime['frequency']

star = xpsi.Star(spacetime = spacetime, photospheres = photosphere)

print("Parameters of the star:")
print(star.params, len(star.params))

# SAX J1808-like 
mass = 1.4
radius = 12.0
distance = 3.5
inclination = 10.0 #60.0
cos_i = math.cos(inclination*math.pi/180.0)

# Hotspot
phase_shift = 0.0
super_colatitude = 105.0*math.pi/180.0 #45.0*math.pi/180.0
super_radius = 1.0*math.pi/180.0 #15.5*math.pi/180.0

# Compton slab model parameters
tbb=0.002 #0.0012
te=100.0 # 50.0
tau=1.6 #1.0

#Tbb = 1 keV <=> tbb = 0.002 (roughly)
#Te = 50 keV <=>  te = 100 (roughly)

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

print(len(p))

# elsewhere
elsewhere_T_keV = 0.4 #  keV

from xpsi.global_imports import  _keV, _k_B
k_B_over_keV = _k_B / _keV
def get_T_in_log10_Kelvin(T_keV):
    """convert T in 10^x K to T in keV"""
    T_log10_Kelvin = np.log10(T_keV/k_B_over_keV)
    return T_log10_Kelvin

elsewhere_T_log10_K = get_T_in_log10_Kelvin(elsewhere_T_keV)

if use_elsewhere:
    p.append(elsewhere_T_log10_K)

star(p)
star.update()
#start = time.time()

#To get the incident signal before interstellar absorption or operating with the telescope:
energies = np.logspace(-1.0, np.log10(16.0), 400, base=10.0)
photosphere.integrate(energies, threads=1) # the number of OpenMP threads to use

#end = time.time()
#print("Time spent in integration:",end - start)
#exit()

print("Bolometric profiles for I, Q, and U:")
print(repr(np.sum(photosphere.signal[0][0], axis=0)))
print(repr(np.sum(photosphere.signalQ[0][0], axis=0)))
print(repr(np.sum(photosphere.signalU[0][0], axis=0)))

from scipy.interpolate import Akima1DInterpolator

class CustomInterstellar(xpsi.Interstellar):
    """ Apply interstellar attenuation. """

    def __init__(self, energies, attenuation, bounds, value):

        assert len(energies) == len(attenuation), 'Array length mismatch.'

        self._lkp_energies = energies # for lookup
        self._lkp_attenuation = attenuation # for lookup      

        N_H = Parameter('column_density',
                        strict_bounds = (0.0,10.0),
                        bounds = bounds,
                        doc = 'Units of 10^21 cm^-2.',
                        symbol = r'$N_{\rm H}$',
                        value = value)

        self._interpolator = Akima1DInterpolator(self._lkp_energies,
                                                 self._lkp_attenuation)
        self._interpolator.extrapolate = True

        super(CustomInterstellar, self).__init__(N_H)

    def attenuation(self, energies):
        """ Interpolate the attenuation coefficients.

        Useful for post-processing. 

        """
        return self._interpolate(energies)**(self['column_density']/1.4)

    def _interpolate(self, energies):
        """ Helper. """
        _att = self._interpolator(energies)
        _att[_att < 0.0] = 0.0
        return _att

    @classmethod
    def from_SWG(cls, path, **kwargs):
        """ Load attenuation file from the NICER SWG. Should be the 1.4e21 cm^-2 file. """

        temp = np.loadtxt(path, dtype=np.double)

        energies = temp[:,0]

        attenuation = temp[:,2]

        return cls(energies, attenuation, **kwargs)

StokesI = photosphere.signal[0][0]
StokesQ = photosphere.signalQ[0][0]
StokesU = photosphere.signalU[0][0]


#Uncomment the following code (and download the required input table) if want to add interstellar attenuation to the modeled signal:
#column_density = 1.17 #10^21 cm^-2
#interstellar = CustomInterstellar.from_SWG(this_directory+'/model_data/tbnew0.14.txt', bounds=(None, None), value=column_density)
#interstellar(energies, StokesI)
#interstellar(energies, StokesQ)
#interstellar(energies, StokesU)

#plt.plot(energies[0:50],np.sum(StokesI,axis=1)[0:50])
#print(energies[0:130])
#print(np.sum(StokesI,axis=1)[0:130])
#plt.ylabel('Flux [?]')
#plt.xlabel('Energy [keV]')
#plt.ylim(0.0,8.0e31)
#plt.savefig("figs/spectrum_after_ism.png")
#exit()

def get_photosphere_stokes_1spot():
    #Return signal from the spot to ixpeobssim
    return hot.phases_in_cycles[0], energies, StokesI, StokesQ, StokesU

