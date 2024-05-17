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
        anegI = (False, False)
        anegQU = (True, True)
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

        else:
            #Note: Signal can be in Qn or Un form if defined so when creating Signal object.
            #In that case: Make sure to have exposure time to 1, background zero, and only 1 hot spot:
            #And also convert signal first to data phase points:
            #(This version is still only for the most simple case)
            #sig1 = self._signals[0][0]
            #fsig = interp1d(self._phases[0], sig1, kind='linear')
            #signal_dphase = fsig(phase_IXPE)

            #self.loglikelihood, self.expected_counts = \
            #    gaussian_likelihood_QnUn(self._data.phases,
            #                              self._data.counts,
            #                              self._data.errors,
            #                              signal_dphase)

            #This for non-normalized Q and U:
            self.loglikelihood, self.expected_counts = \
                gaussian_likelihood_given_background(self._data.exposure_time,
                                          self._data.phases,
                                          self._data.counts,
                                          self._data.errors,
                                          self._signals,
                                          self._phases,
                                          self._shifts,
                                          background = background,
                                          allow_negative=anegQU)

signals = [[],]

signalI = CustomSignal_gaussian(data = IXPE_I.data,
                instrument = IXPE_du1,
                #background = background,
                interstellar = None,
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
                        interstellar = None,
                        workspace_intervals = 1000,
                        cache = True,
                        epsrel = 1.0e-8,
                        epsilon = 1.0e-3,
                        sigmas = 10.0,
                        support = None,
                        stokes="Q")

signals[0].append(signalQ)

signalU = CustomSignal_gaussian(data = IXPE_U.data,
	                instrument = IXPE_du1,
	                #background = background,
	                interstellar = None,
	                workspace_intervals = 1000,
	                cache = True,
	                epsrel = 1.0e-8,
	                epsilon = 1.0e-3,
	                sigmas = 10.0,
	                support = None,
	                stokes="U")
signals[0].append(signalU)


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
                                    image_order_limit=3,
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

bounds = dict(spin_axis_position_angle = (None, None))
photosphere = CustomPhotosphere_NumA5(hot = hot, elsewhere = elsewhere, stokes=True,
                                values=dict(mode_frequency = spacetime['frequency']),
                                bounds=bounds)

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
     0.0,
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


star = xpsi.Star(spacetime = spacetime, photospheres = photosphere)

likelihood = xpsi.Likelihood(star = star, signals = signals,
                             num_energies=128,
                             threads=1,
                             externally_updated=False)


from scipy.stats import truncnorm
class CustomPrior(xpsi.Prior):
    """ A custom (joint) prior distribution.
    """

    __derived_names__ = ['compactness', 'phase_separation',]

    def __init__(self):
        super(CustomPrior, self).__init__() # not strictly required if no hyperparameters

    def __call__(self, p = None):
        """ Evaluate distribution at ``p``.
        :param list p: Model parameter values.
        :returns: Logarithm of the distribution evaluated at ``p``.

        """
        temp = super(CustomPrior, self).__call__(p)
        if not np.isfinite(temp):
            return temp

        ## based on contemporary EOS theory
        if not self.parameters['radius'] <= 16.0:
            return -np.inf

        ref = self.parameters.star.spacetime # shortcut

        # limit polar radius to be outside the Schwarzschild photon sphere
        R_p = 1.0 + ref.epsilon * (-0.788 + 1.030 * ref.zeta)
        if R_p < 1.505 / ref.R_r_s:
            return -np.inf

        mu = math.sqrt(-1.0 / (3.0 * ref.epsilon * (-0.788 + 1.030 * ref.zeta)))

        # 2-surface cross-section have a single maximum in |z|
        # i.e., an elliptical surface; minor effect on support, if any,
        # for high spin frequenies
        if mu < 1.0:
            return -np.inf

        ref = self.parameters # redefine shortcut

        # enforce order in hot region colatitude
        if ref['p__super_colatitude'] > ref['s__super_colatitude']:
                return -np.inf

        # hot regions cannot overlap
        #To be added...
        return 0.0

    def inverse_sample(self, hypercube=None):
        """ Draw sample uniformly from the distribution via inverse sampling. """

        to_cache = self.parameters.vector

        if hypercube is None:
            hypercube = np.random.rand(len(self))

        # the base method is useful, so to avoid writing that code again:
        _ = super(CustomPrior, self).inverse_sample(hypercube)

        ref = self.parameters # shortcut)

        # flat priors in cosine of hot region centre colatitudes (isotropy)
        # support modified by no-overlap rejection condition
        idx = ref.index('p__super_colatitude')
        a, b = ref.get_param('p__super_colatitude').bounds
        a = math.cos(a); b = math.cos(b)
        ref['p__super_colatitude'] = math.acos(b + (a - b) * hypercube[idx])

        idx = ref.index('s__super_colatitude')
        a, b = ref.get_param('s__super_colatitude').bounds
        a = math.cos(a); b = math.cos(b)
        ref['s__super_colatitude'] = math.acos(b + (a - b) * hypercube[idx])

        # restore proper cache
        for parameter, cache in zip(ref, to_cache):
            parameter.cached = cache

        # it is important that we return the desired vector because it is
        # automatically written to disk by MultiNest and only by MultiNest
        return self.parameters.vector

    def transform(self, p, **kwargs):
        """ A transformation for post-processing. """

        p = list(p) # copy

        # used ordered names and values
        ref = dict(zip(self.parameters.names, p))

        return p

prior = CustomPrior()
likelihood.prior = prior

wrapped_params = [0]*len(likelihood)
wrapped_params[likelihood.index('p__phase_shift')] = 1
wrapped_params[likelihood.index('s__phase_shift')] = 1

runtime_params = {'resume': False,
                  'importance_nested_sampling': False,
                  'multimodal': False,
                  'n_clustering_params': None,
                  'outputfiles_basename': './run/run_PolNum',
                  'n_iter_before_update': 50,
                  'n_live_points': 50,
                  'sampling_efficiency': 0.8,
                  'const_efficiency_mode': False,
                  'wrapped_params': wrapped_params,
                  'evidence_tolerance': 0.5,
                  'seed': 7,
                  'max_iter': 100, # manual termination condition for short test
                  'verbose': True}

likelihood.reinitialise()
likelihood.clear_cache()

true_logl = -1.2738517361e+06

if __name__ == '__main__': # sample from the posterior
    # inform source code that parameter objects updated when inverse sampling
    likelihood.externally_updated = True
    # let's require that checks pass before starting to sample
    check_kwargs = dict(hypercube_points = None,
                    physical_points = p, # externally_updated preserved
                    loglikelihood_call_vals = [true_logl],
                    rtol_loglike = 1.0e-6) # choose a tolerance
    xpsi.Sample.nested(likelihood, prior, check_kwargs, **runtime_params)
