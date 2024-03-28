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


from xpsi.Parameter import Parameter
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

fname_ixpedata = "/home/xiaotuo/ixpeobssimdata/ixpeobssimdata_24_set1Bas_dt0/model_amsp_xpsi"
#fname_ixpedata = "/home/xiaotuo/ixpeobssimdata/ixpeobssimdata_24_set3_dt0/model_amsp_xpsi"

#In the end, we probably want to use data in PHA format instead of PCUBE, to properly account for the instrument response.
#For now, we just read some PCUBE data (created by ixpeobssim) for testing purposes (binned in 1 energy channel).
phase_IXPE, Idat, qn, un, Iderr, qnerr, unerr, PD, PDerr, keVdat, MDP99 = readData_pcube_ebin(fname_ixpedata)

IXPE_I.data = xpsi.Data([Idat[:,0]],
                       channels=np.arange(0, 1),
                       phases=np.linspace(0,1,len(phase_IXPE)+1),
                       first=0,
                       last=0,
                       exposure_time=1.0)
IXPE_Q.data = xpsi.Data([qn[:,0]],
                       channels=np.arange(0, 1),
                       phases=np.linspace(0,1,len(phase_IXPE)+1),
                       first=0,
                       last=0,
                       exposure_time=1.0)
IXPE_U.data = xpsi.Data([un[:,0]],
                       channels=np.arange(0, 1),
                       phases=np.linspace(0,1,len(phase_IXPE)+1),
                       first=0,
                       last=0,
                       exposure_time=1.0)
IXPE_I.data.errors, IXPE_Q.data.errors, IXPE_U.data.errors = Iderr, qnerr, unerr

from ixpe_read import read_response_IXPE

class CustomInstrument_stokes(xpsi.Instrument):
    """ A model of the NICER telescope response. """

    def __call__(self, signal, *args):
        """ Overwrite base just to show it is possible.

        We loaded only a submatrix of the total instrument response
        matrix into memory, so here we can simplify the method in the
        base class.

        """
        matrix = self.construct_matrix()

        self._folded_signal = np.dot(matrix, signal)

        return self._folded_signal

    @classmethod
    def from_response_files(cls, MRF, RMF, max_input, max_channel, min_input=0, min_channel=0,
                            channel_edges=None):
        """ Constructor which converts response files into :class:`numpy.ndarray`s.
        :param str MRF: Path to MRF which is compatible with
                                :...
        :param str RMF: Path to RMF which is compatible with
                                :...
        :param str channel_edges: Optional path to edges which is compatible with
                                  :func:`numpy.loadtxt`.
        """
        if min_input != 0:
            min_input = int(min_input)
        max_input = int(max_input)
        try:
            matrix, edges, channels, channel_edgesT = read_response_IXPE(MRF,RMF,min_input,max_input,min_channel,max_channel)
            if channel_edges:
                channel_edgesT = np.loadtxt(channel_edges, dtype=np.double, skiprows=3)[:,1:]
        except:
            print('A file could not be loaded.')
            raise
        return cls(matrix, edges, channels, channel_edgesT)

#Let's test modeling using just the instrument files for the first detector unit of IXPE:
IXPE_du1 = CustomInstrument_stokes.from_response_files(MRF = this_directory+'/model_data/ixpe_d1_obssim_v012.mrf',
                                             RMF = this_directory+'/model_data/ixpe_d1_obssim_v012.rmf',
                                             max_input = 275, #175, #275,
                                             max_channel = 200,
                                             min_input = 0, #25, #25 for 2 keV, 75 for 4 keV
                                             min_channel = 50, #50 for 2 keV, 100 for 4 keV
                                             channel_edges = None)

from xpsi.likelihoods._gaussian_likelihood_QnUn import gaussian_likelihood_QnUn
from xpsi.likelihoods._gaussian_likelihood_given_background_IQU import gaussian_likelihood_given_background

from scipy.interpolate import interp1d

class CustomSignal_gaussian(xpsi.Signal):
    """

    A custom calculation of the logarithm of the likelihood.
    We extend the :class:`~xpsi.Signal.Signal` class to make it callable.
    We overwrite the body of the __call__ method. The docstring for the
    abstract method is copied.

    """

    def __init__(self, workspace_intervals = 1000, epsabs = 0, epsrel = 1.0e-8,
                 epsilon = 1.0e-3, sigmas = 10.0, support = None, **kwargs):
        """ Perform precomputation.

        :params ndarray[m,2] support:
            Prior support bounds for background count rate variables in the
            :math:`m` instrument channels, where the lower bounds must be zero
            or positive, and the upper bounds must be positive and greater than
            the lower bound. Alternatively, setting the an upper bounds as
            negative means the prior support is unbounded and the flat prior
            density functions per channel are improper. If ``None``, the lower-
            bound of the support for each channel is zero but the prior is
            unbounded.

        """

        super(CustomSignal_gaussian, self).__init__(**kwargs)

    def __call__(self, *args, **kwargs):
        anegI = (False)
        anegQU = (True)
        background = np.zeros((np.shape(self._data.counts))) #self._background.registered_background
        if self.isI:       
            self.loglikelihood, self.expected_counts = \
                gaussian_likelihood_given_background(self._data.exposure_time,
                                          self._data.phases,
                                          self._data.counts,
                                          self._data.errors,
                                          self._signals,
                                          self._phases,
                                          self._shifts,
                                          background = background,
                                          allow_negative=anegI)
            #print("loglikelihood:",self.loglikelihood)
        else:
            #Note: Signal can be in Qn or Un form if defined so when creating Signal object.
            #In that case: Make sure to have exposure time to 1, background zero, and only 1 hot spot:
            #And also convert signal first to data phase points:
            #(This version is still only for the most simple case)
            sig1 = self._signals[0][0]
            fsig = interp1d(self._phases[0], sig1, kind='linear')
            #signal_dphase = fsig(self._data.phases) #this is wrong
            signal_dphase = fsig(phase_IXPE)

            self.loglikelihood, self.expected_counts = \
                gaussian_likelihood_QnUn(self._data.phases,
                                          self._data.counts,
                                          self._data.errors,
                                          signal_dphase)

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
        #print("att:",_att)
        return _att

    @classmethod
    def from_SWG(cls, path, **kwargs):
        """ Load attenuation file from the NICER SWG. Should be the 1.4e21 cm^-2 file. """

        temp = np.loadtxt(path, dtype=np.double)

        energies = temp[:,0]

        attenuation = temp[:,2]

        return cls(energies, attenuation, **kwargs)

column_density = 1.17 #10^21 cm^-2
interstellar = CustomInterstellar.from_SWG(this_directory+'/model_data/tbnew0.14.txt', bounds=(None, None), value=column_density)

signals = [[],]

signalI = CustomSignal_gaussian(data = IXPE_I.data,
                instrument = IXPE_du1,
                #background = background,
                interstellar = interstellar,
                workspace_intervals = 1000,
                cache = True,
                epsrel = 1.0e-8,
                epsilon = 1.0e-3,
                sigmas = 10.0,
                support = None,
                stokes="I")

signals[0].append(signalI)

signalQ = CustomSignal_gaussian(data = IXPE_Q.data,
                        instrument = IXPE_du1,
                        #background = background,
                        interstellar = interstellar,
                        workspace_intervals = 1000,
                        cache = True,
                        epsrel = 1.0e-8,
                        epsilon = 1.0e-3,
                        sigmas = 10.0,
                        support = None,
                        stokes="Qn")

signals[0].append(signalQ)

signalU = CustomSignal_gaussian(data = IXPE_U.data,
	                instrument = IXPE_du1,
	                #background = background,
	                interstellar = interstellar,
	                workspace_intervals = 1000,
	                cache = True,
	                epsrel = 1.0e-8,
	                epsilon = 1.0e-3,
	                sigmas = 10.0,
	                support = None,
	                stokes="Un")
signals[0].append(signalU)

bounds = dict(distance = (0.1, 1.0),                     # (Earth) distance
                mass = (1.0, 3.0),                       # mass
                radius = (3.0 * gravradius(1.0), 16.0),  # equatorial radius
                cos_inclination = (0.0, 1.0))      # (Earth) inclination to rotation axis

spacetime = xpsi.Spacetime(bounds=bounds, values=dict(frequency=400.9752075))

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

from xpsi import HotRegions
hot = HotRegions((primary,))

use_elsewhere = True

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

#set1
#p = [1.4, 12.0, 3.5, 0.5000000000000001, 0.0, 0.7853981633974483, 0.6981317007977318, 0.002, 100.0, 1.0]

#set1Bas
p = [1.4, 12.0, 3.5, 0.5000000000000001, 0.0, 0.7853981633974483, 0.27052603405912107, 0.0012, 100.0, 1.0]

#set2:
#p = [1.4, 12.0, 3.5, 0.984807753012208, 0.0, 1.8325957145940461, 0.6981317007977318, 0.002, 100.0, 1.6]

#set3
#p = [1.4, 12.0, 3.5, 0.984807753012208, 0.0, 1.8325957145940461, 0.27052603405912107, 0.002, 100.0, 1.6]


print(len(p))
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

#overwriting p with the found maxL solution (+distance= 1.0 or 3.5):
#p = [2.03431421e+00, 1.57761215e+01, 3.5, 4.72596557e-02, 4.77767506e-01,
# 7.49863991e-01, 1.35199302e+00, 1.68339025e-03, 1.16438934e+02,
# 3.30295403e+00, 5.76610075e+00]


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

StokesI = photosphere.signal[0][0]
StokesQ = photosphere.signalQ[0][0]
StokesU = photosphere.signalU[0][0]

interstellar(energies, StokesI)
interstellar(energies, StokesQ)
interstellar(energies, StokesU)

print(StokesI)
print(StokesQ)
print(StokesU)

def find_idx(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx

phase1 = hot.phases_in_cycles[0]
print(hot.phases_in_cycles)

Imod = np.zeros((len(phase1)))
Qmod = np.zeros((len(phase1)))
Umod = np.zeros((len(phase1)))

#print(energies)
for e in range(find_idx(energies,2.0),find_idx(energies,8.0)):
	Imod[:] = Imod[:] + (photosphere.signal[0][0][e,:]+photosphere.signal[0][0][e+1,:])*(energies[e+1]-energies[e])	
	Qmod[:] = Qmod[:] + (photosphere.signalQ[0][0][e,:]+photosphere.signalQ[0][0][e+1,:])*(energies[e+1]-energies[e])
	Umod[:] = Umod[:] + (photosphere.signalU[0][0][e,:]+photosphere.signalU[0][0][e+1,:])*(energies[e+1]-energies[e])

Imod = 1/2*Imod
Qmod = 1/2*Qmod
Umod = 1/2*Umod

########################################
#print(phase1)
#print(phase_IXPE)
ph_bin_len = phase_IXPE[1]-phase_IXPE[0]
Imod_int = np.zeros((len(phase_IXPE)))
Qmod_int = np.zeros((len(phase_IXPE)))
Umod_int = np.zeros((len(phase_IXPE)))
for ipa in range(0,len(phase_IXPE)):
    for ipb in range(find_idx(phase1,phase_IXPE[ipa]-0.5*ph_bin_len),find_idx(phase1,phase_IXPE[ipa]+0.5*ph_bin_len)):
        Imod_int[ipa] = Imod_int[ipa] + (Imod[ipb]+Imod[ipb+1])*(phase1[ipb+1]-phase1[ipb])
        Qmod_int[ipa] = Qmod_int[ipa] + (Qmod[ipb]+Qmod[ipb+1])*(phase1[ipb+1]-phase1[ipb])
        Umod_int[ipa] = Umod_int[ipa] + (Umod[ipb]+Umod[ipb+1])*(phase1[ipb+1]-phase1[ipb])                

Imod_int = 1/2*Imod_int
Qmod_int = 1/2*Qmod_int
Umod_int = 1/2*Umod_int

Inorm = np.zeros((len(phase_IXPE)))
Qnorm = np.zeros((len(phase_IXPE)))
Unorm = np.zeros((len(phase_IXPE)))

for iph in range(0, len(Imod_int[:])):
    if (Imod[iph] < 1e-10):
        Qnorm[iph] = 0.0
        Unorm[iph] = 0.0
    else:
        Qnorm[iph] = Qmod_int[iph] / Imod_int[iph]
        Unorm[iph] = Umod_int[iph] / Imod_int[iph]           
Inorm = Imod_int/np.max(Imod_int)

Imod_int = Inorm
Qmod_int = Qnorm
Umod_int = Unorm
######################################


#######################################
#Inorm = np.zeros((len(phase1)))
#Qnorm = np.zeros((len(phase1)))
#Unorm = np.zeros((len(phase1)))

#for iph in range(0, len(Imod[:])):
#    if (Imod[iph] < 1e-10):
#        Qnorm[iph] = 0.0
#        Unorm[iph] = 0.0
#    else:
#        Qnorm[iph] = Qmod[iph] / Imod[iph]
#        Unorm[iph] = Umod[iph] / Imod[iph]           
#Inorm = Imod/np.max(Imod)

#from scipy.interpolate import interp1d
#Imod_int = interp1d(phase1, Inorm, kind='linear')(phase_IXPE)
#Qmod_int = interp1d(phase1, Qnorm, kind='linear')(phase_IXPE)
#Umod_int = interp1d(phase1, Unorm, kind='linear')(phase_IXPE)

#print(Idat[:,0])
#print(Iderr[:,0])
#exit()
######################################

Idatn = Idat/np.max(Idat[:,0])
Iderrn = Iderr/np.max(Idat[:,0])

fig, ax = plt.subplots(1,3)
ax[0].plot(phase_IXPE, Imod_int,color="darkorange",linewidth=1.0)
ax[0].errorbar(phase_IXPE, Idatn[:,0], yerr=Iderrn[:,0],fmt='.',markersize=0.5,color="blue")
ax[0].set_ylabel('I/Imax')
ax[1].errorbar(phase_IXPE, qn[:,0], yerr=qnerr[:,0],fmt='.')
ax[1].plot(phase_IXPE, Qmod_int)
ax[1].set_ylabel('Q/I')
ax[2].errorbar(phase_IXPE, un[:,0], yerr=unerr[:,0], fmt='.', label='data')
ax[2].plot(phase_IXPE, Umod_int, label='model')
ax[2].set_ylabel('U/I')
for iiii in range(3): ax[iiii].set_xlabel('Phase')
ax[2].legend()
fig.tight_layout()              
plt.savefig("figs/stokes_pulse_comp_xpsi.png",dpi=300.0)

fig, ax = plt.subplots(2,1)
ax[0].errorbar(phase_IXPE, Idatn[:,0]-Imod_int,yerr=Iderrn[:,0],fmt='.',markersize=4.0,color="blue")
ax[1].errorbar(phase_IXPE, Idatn[:,0]/Imod_int,yerr=Iderrn[:,0]/Imod_int,fmt='.',markersize=4.0,color="blue")
ax[0].set_ylabel('I(data)-I(model)')
ax[1].set_ylabel('I(data)/I(model)')
for iiii in range(2): ax[iiii].set_xlabel('Phase')
fig.tight_layout()
plt.savefig("figs/stokes_pulse_comp_xpsi_resid.png",dpi=300.0)



likelihood = xpsi.Likelihood(star = star, signals = signals,
                             num_energies=128,
                             threads=1,
                             externally_updated=False)
                             
p.append(column_density)

#overwriting p with the found maxL solution (+distance= 1.0 or 3.5):

#p = [2.03431421e+00, 1.57761215e+01, 3.5, 4.72596557e-02, 4.77767506e-01,
# 7.49863991e-01, 1.35199302e+00, 1.68339025e-03, 1.16438934e+02,
# 3.30295403e+00, 5.76610075e+00, 8.52169979e+00]

print(likelihood, p)
true_logl = -1.4223095642e+05 #-1.2543696401e+04 #-1.3357765522e+05 #-9.1393983034e+05 #-9.1394271354e+05 #-1.4223099929e+05 #-9.1394280091e+0
likelihood.check(None, [true_logl], 1.0e-6,physical_points=[p],force_update=True)

Isig = np.sum(signals[0][0].signals[0], axis=0)
Isign = Isig/np.max(Isig)
print(np.shape(Isign),np.shape(Imod_int))
print(Isign)

Iexpec = np.sum(signals[0][0].expected_counts, axis=0)
Iexpecn = Iexpec/np.max(Iexpec)
#print(np.shape(signals[0][0].expected_counts))
print(Iexpecn)
print(Idatn)

qexpec = signals[0][1].expected_counts
uexpec = signals[0][2].expected_counts

print(qexpec)
print("uexpec:",uexpec)


fig, ax = plt.subplots(1,3)
ax[0].plot(phase_IXPE, Iexpecn,color="darkorange",linewidth=1.0)
ax[0].errorbar(phase_IXPE, Idatn[:,0], yerr=Iderrn[:,0],fmt='.',markersize=0.5,color="blue")
ax[0].set_ylabel('I/Imax')
ax[1].errorbar(phase_IXPE, qn[:,0], yerr=qnerr[:,0],fmt='.')
ax[1].plot(phase_IXPE, qexpec)
ax[1].set_ylabel('Q/I')
ax[2].errorbar(phase_IXPE, un[:,0], yerr=unerr[:,0], fmt='.', label='data')
ax[2].plot(phase_IXPE, uexpec, label='model')
ax[2].set_ylabel('U/I')
for iiii in range(3): ax[iiii].set_xlabel('Phase')
ax[2].legend()
fig.tight_layout()              
plt.savefig("figs/stokes_signal_comp_xpsi.png",dpi=300.0)


fig, ax = plt.subplots(2,1)
ax[0].errorbar(phase_IXPE, Idatn[:,0]-Iexpecn,yerr=Iderrn[:,0],fmt='.',markersize=4.0,color="blue")
ax[1].errorbar(phase_IXPE, Idatn[:,0]/Iexpecn,yerr=Iderrn[:,0]/Iexpecn,fmt='.',markersize=4.0,color="blue")
ax[0].set_ylabel('I(data)-I(model)')
ax[1].set_ylabel('I(data)/I(model)')
for iiii in range(2): ax[iiii].set_xlabel('Phase')
fig.tight_layout()
plt.savefig("figs/stokes_signal_comp_xpsi_resid.png",dpi=300.0)


#Plot the stokes profiles in the fancier format:
import matplotlib
matplotlib.use('agg')
from pylab import *

rc("text", usetex=True)
fig = figure(figsize=(8,10), dpi=300)
plt.figure(1)
fig.subplots_adjust(hspace=0)
fig.subplots_adjust(bottom=0.15)
fig.subplots_adjust(left=0.18)
lbfontsz = 30
lwidth= 1.0
rc("xtick", labelsize=lbfontsz)
rc("ytick", labelsize=lbfontsz)
rc("axes", linewidth=lwidth)
labelx = -0.16
msize = 2.0

ax2 = plt.subplot(3,1,1)
ax3 = plt.subplot(3,1,2)
ax4 = plt.subplot(3,1,3)                

PD_XPSI = np.sqrt(Qmod**2+Umod**2)/Imod

ax2.tick_params(bottom=True, top=True, left=True, right=True, direction="in", length=10)
ax2.plot(phase1, PD_XPSI,linewidth=lwidth,color="blue")			
ax2.errorbar(phase_IXPE,MDP99[:,0],yerr=0.0,xerr=0.03,fmt='o',capsize=2.0,markersize=msize,color="red",label="MDP99")
ax2.errorbar(phase_IXPE, PD[:,0], yerr=PDerr[:,0], xerr=0.0, fmt='o', color="purple",capsize=2.0,markersize=msize)
ax2.set_xlim(0., 1.)
ax2.set_ylim(0.0, 0.12)
ax2.set_ylabel("$P_{\mathrm{obs}}$",fontsize=lbfontsz)
ax2.set_yticks([0.0,0.05,0.1])
ax2.set_yticklabels(["$0.0$","$0.05$","$0.10$"])
ax2.yaxis.set_label_coords(labelx, 0.5)
                    
plt.setp(ax2.get_xticklabels(), visible=False)

ax3.tick_params(bottom=True, top=True, left=True, right=True, direction="in", length=10)
ax3.plot(phase_IXPE, qexpec,linewidth=lwidth,color="blue")			
ax3.errorbar(phase_IXPE, qn[:,0], yerr=qnerr[:,0], xerr=0.0, fmt='o', color="purple",capsize=2.0,markersize=msize)
ax3.set_xlim(0., 1.)
ax3.set_ylim(-0.08, 0.08)
ax3.set_ylabel("$q$",fontsize=lbfontsz)
ax3.yaxis.set_label_coords(labelx, 0.5)
            
ax4.tick_params(bottom=True, top=True, left=True, right=True, direction="in", length=10)
ax4.plot(phase_IXPE, uexpec,linewidth=lwidth,color="blue")
ax4.errorbar(phase_IXPE, un[:,0], yerr=unerr[0,:], xerr=0.0, fmt='o', color="purple",capsize=2.0,markersize=msize)
ax4.set_xlim(0., 1.)
ax4.set_ylim(-0.08, 0.08)
ax4.set_ylabel("$q$",fontsize=lbfontsz)
ax4.set_xlabel("Phase $\phi/2\pi$",fontsize=lbfontsz)
ax4.yaxis.set_label_coords(labelx, 0.5)

ax3.set_yticks([-0.05,0.0,0.05])
ax3.set_yticklabels(["$-0.05$","$0.0$","$0.05$"])
ax4.set_yticks([-0.05,0.0,0.05])
ax4.set_yticklabels(["$-0.05$","$0.0$","$0.05$"])
ax4.set_xticks([0.0,0.25,0.50,0.75,1.0])
ax4.set_xticklabels(["$0.0$","$0.25$","$0.50$","$0.75$","$1.0$"])

plt.setp(ax3.get_xticklabels(), visible=False)

#    ic = ic+1

fig.savefig("figs/stokes_signalsX.png",dpi=300.0)	
plt.close()


