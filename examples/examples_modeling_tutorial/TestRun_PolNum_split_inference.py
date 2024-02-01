'''
Test script with the polarized 3+2 numerical atmosphere applied to a ST+PDT model.
Example of fitting Stokes signals observed by an instrument will be added here.

Prequisities:
Before running the script, add the atmosphere data to the model_data subdirectory:
Bobrikova_compton_slab_I.npz and Bobrikova_compton_slab_Q.npz. See the
example script in xpsi/examples/produce_atmos_lookuptable for producing these files
from those provided in https://github.com/AnnaBobrikova/ComptonSlabTables.
In addition, the simulated polarization data files need to placed in
xpsi/examples/examples_modeling_tutorial/ixpeobssim/ixpeobssimdata/.
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

import sys
sys.path.append(this_directory+"/ixpeobssim/")

class namespace():
    pass

IXPE_I = namespace()
IXPE_Q = namespace()
IXPE_U = namespace()

from ixpe_read import readData_pcube_ebin

fname_ixpedata = this_directory+"/ixpeobssim/ixpeobssimdata/model_amsp_xpsi"

#In the end, we probably want to use data in PHA format instead of PCUBE, to properly account for the instrument response.
#For now, we just read some PCUBE data (created by ixpeobssim) for testing purposes (binned in 1 energy channel).
phase_IXPE, Idat, qn, un, Iderr_du1, qnerr_du1, unerr_du1, keVdat = readData_pcube_ebin(fname_ixpedata)

IXPE_I.data = xpsi.Data([Idat[:,0]],
                       channels=np.arange(0, 1),
                       phases=phase_IXPE,
                       first=0,
                       last=0,
                       exposure_time=1.0)
IXPE_Q.data = xpsi.Data([qn[:,0]],
                       channels=np.arange(0, 1),
                       phases=phase_IXPE,
                       first=0,
                       last=0,
                       exposure_time=1.0)
IXPE_U.data = xpsi.Data([un[:,0]],
                       channels=np.arange(0, 1),
                       phases=phase_IXPE,
                       first=0,
                       last=0,
                       exposure_time=1.0)

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
                            prefix='p')

bounds2 = dict(super_colatitude = (None, None),
                        super_radius = (None, None),
                        phase_shift = (0.0, 0.1),
                        super_tbb = (0.001, 0.003),
                        super_tau = (0.5, 3.5),
                        super_te = (40.0, 200.0),
                        cede_colatitude = (None, None),
                        cede_radius = (None, None),
                        cede_azimuth = (None, None),
                        cede_tbb = (0.001, 0.003),
                        cede_tau = (0.5, 3.5),
                        cede_te = (40.0, 200.0))

secondary = CustomHotRegion_Accreting(bounds=bounds2, # can otherwise use same bounds
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
                            is_antiphased=True,
                            split=True,
                            prefix='s')


from xpsi import HotRegions
hot = HotRegions((primary, secondary))

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

p = [1.0368513939430604,
     6.087862992320039,
     0.26870812456714116,
     0.39140510783272897,
     0.04346870860640872,
     0.8002010406881243,
     1.1165398710637626,
     0.0015,
     100.0,
     1.0,
     0.07360477761463673,
     2.4602238829718432,
     0.4277092192054918,
     2.3, #cede_colatitude
     0.5, #cede_radius
     0.1, #cede_azimuth
     0.0015,
     100.0,
     1.0,     
     0.0015,
     100.0,
     1.0]
print(len(p))

#Tbb = 1 keV <=> tbb = 0.002 (roughly)
#Te = 50 keV <=>  te = 100 (roughly)

if use_elsewhere:
    p.append(5.5)

star(p)
star.update()
#start = time.time()

#To get the incident signal before interstellar absorption or operating with the telescope:
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

