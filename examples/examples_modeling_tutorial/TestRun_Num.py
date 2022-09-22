'''
Test script to check that X-PSI installation is working (with the numerical atmosphere extension).
The run is succesful if it passes the likelihood check (in ~ 1 minute) and sampling finishes without errors (in ~ 5 min).
The model and synthetic data are based on those shown in the Modelling tutorial notebook.

Prequisities:
Before running the script, add the NICER instrument files to the model_data subdirectory:
nicer_v1.01_arf.txt, nicer_v1.01_rmf_energymap.txt, and nicer_v1.01_rmf_matrix.txt, nsx_H_v200804.out
(found from https://doi.org/10.5281/zenodo.7094144).
'''

from __future__ import print_function, division

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

#using an example synthetic data
settings = dict(counts = np.loadtxt('model_data/example_synthetic_realisation.dat', dtype=np.double),
        channels=np.arange(20,201),
        phases=np.linspace(0.0, 1.0, 33),
        first=0, last=180,
        exposure_time=984307.6661)
data = xpsi.Data(**settings)

class CustomInstrument(xpsi.Instrument):
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
    def from_response_files(cls, ARF, RMF, max_input, min_input=0,
                            channel_edges=None):
        """ Constructor which converts response files into :class:`numpy.ndarray`s.
        :param str ARF: Path to ARF which is compatible with
                                :func:`numpy.loadtxt`.
        :param str RMF: Path to RMF which is compatible with
                                :func:`numpy.loadtxt`.
        :param str channel_edges: Optional path to edges which is compatible with
                                  :func:`numpy.loadtxt`.
        """
        link ="https://doi.org/10.5281/zenodo.7094144"
        if min_input != 0:
            min_input = int(min_input)

        max_input = int(max_input)

        try:
            ARF = np.loadtxt(ARF, dtype=np.double, skiprows=3)
        except:
            print("ERROR: You miss the following file: {}\nThe file is found from here: {}".format(ARF, link))
            exit()

        try:
            RMF = np.loadtxt(RMF, dtype=np.double)
        except:
            print("ERROR: You miss the following file: {}\nThe file is found from here: {}".format(RMF, link))
            exit()
        if channel_edges:
            try:
                channel_edges = np.loadtxt(channel_edges, dtype=np.double, skiprows=3)[:,1:]
            except:
                print("ERROR: You miss the following file: {}\nThe file is found from here: {}".format(channel_edges, link))
                exit()

        matrix = np.ascontiguousarray(RMF[min_input:max_input,20:201].T, dtype=np.double)

        edges = np.zeros(ARF[min_input:max_input,3].shape[0]+1, dtype=np.double)

        edges[0] = ARF[min_input,1]; edges[1:] = ARF[min_input:max_input,2]

        for i in range(matrix.shape[0]):
            matrix[i,:] *= ARF[min_input:max_input,3]

        channels = np.arange(20, 201)

        return cls(matrix, edges, channels, channel_edges[20:202,-2])



NICER = CustomInstrument.from_response_files(ARF = 'model_data/nicer_v1.01_arf.txt',
                                     RMF = 'model_data/nicer_v1.01_rmf_matrix.txt',
                                     max_input = 500,
                                     min_input = 0,
                                     channel_edges = 'model_data/nicer_v1.01_rmf_energymap.txt')

bounds = dict(distance = (0.1, 1.0),                     # (Earth) distance
                mass = (1.0, 3.0),                       # mass
                radius = (3.0 * gravradius(1.0), 16.0),  # equatorial radius
                cos_inclination = (0.0, 1.0))      # (Earth) inclination to rotation axis

spacetime = xpsi.Spacetime(bounds=bounds, values=dict(frequency=300.0))

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
	                    prefix='p')

class derive(xpsi.Derive):
    def __init__(self):
        """
        We can pass a reference to the primary here instead
        and store it as an attribute if there is risk of
        the global variable changing.

        This callable can for this simple case also be
        achieved merely with a function instead of a magic
        method associated with a class.
        """
        pass

    def __call__(self, boundto, caller = None):
        # one way to get the required reference
        global primary # unnecessary, but for clarity
        return primary['super_temperature'] - 0.2


bounds['super_temperature'] = None # declare fixed/derived variable

secondary = xpsi.HotRegion(bounds=bounds, # can otherwise use same bounds
	                      values={'super_temperature': derive()},
	                      symmetry=True,
	                      omit=False,
	                      cede=False,
	                      concentric=False,
	                      sqrt_num_cells=32,
	                      min_sqrt_num_cells=10,
	                      max_sqrt_num_cells=100,
	                      num_leaves=100,
	                      num_rays=200,
	                      do_fast=False,
	                      is_antiphased=True,
	                      prefix='s')


from xpsi import HotRegions
hot = HotRegions((primary, secondary))
h = hot.objects[0]
hot['p__super_temperature'] = 6.0 # equivalent to ``primary['super_temperature'] = 6.0``

class CustomPhotosphere_num(xpsi.Photosphere):
    """ A photosphere extension to preload the numerical atmosphere NSX. """

    @xpsi.Photosphere.hot_atmosphere.setter
    def hot_atmosphere_old(self, path):
        try:
            NSX = np.loadtxt(path, dtype=np.double)
        except:
            print("ERROR: You miss the following file:", path)
            print("The file is found from here: https://doi.org/10.5281/zenodo.7094144")
            exit()
        logT = np.zeros(35)
        logg = np.zeros(14)
        mu = np.zeros(67)
        logE = np.zeros(166)

        #reorder_buf = np.zeros((35,11,67,166))
        reorder_buf = np.zeros((35,14,67,166))

        index = 0
        for i in range(reorder_buf.shape[0]):
            for j in range(reorder_buf.shape[1]):
                for k in range(reorder_buf.shape[3]):
                   for l in range(reorder_buf.shape[2]):
                        logT[i] = NSX[index,3]
                        logg[j] = NSX[index,4]
                        logE[k] = NSX[index,0]
                        mu[reorder_buf.shape[2] - l - 1] = NSX[index,1]
                        reorder_buf[i,j,reorder_buf.shape[2] - l - 1,k] = 10.0**(NSX[index,2])
                        index += 1

        buf = np.zeros(np.prod(reorder_buf.shape))

        bufdex = 0
        for i in range(reorder_buf.shape[0]):
            for j in range(reorder_buf.shape[1]):
                for k in range(reorder_buf.shape[2]):
                   for l in range(reorder_buf.shape[3]):
                        buf[bufdex] = reorder_buf[i,j,k,l]; bufdex += 1

        self._hot_atmosphere = (logT, logg, mu, logE, buf)

    @xpsi.Photosphere.hot_atmosphere.setter
    def hot_atmosphere(self,path):
        size=(35, 14, 67, 166)
        path_npy = path[0:len(path)-3]+"npy"
        try:
            NSX = np.load(path_npy)
        except:
            try:
                NSX = np.loadtxt(path, dtype=np.double)
            except:
                print("ERROR: You miss the following file:", path)
                print("The file is found from here: https://doi.org/10.5281/zenodo.7094144")
                exit()
            np.save(path_npy,NSX)

        def reorder_23(array, size):
            new_array=np.zeros(size)
            index=0
            for i in range(size[3]):
                 for j in range(size[2]):
                      new_array[:,:,j,i]=array[:,:,index]
                      index+=1
            return new_array

        _mu_opt = np.ascontiguousarray(NSX[0:size[2],1][::-1])
        logE_opt = np.ascontiguousarray([NSX[i*size[2],0] for i in range(size[3])])
        logT_opt = np.ascontiguousarray([NSX[i*size[1]*size[2]*size[3],3] for i in range(size[0])])
        logg_opt = np.ascontiguousarray([NSX[i*size[2]*size[3],4] for i in range(size[1])])

        reorder_buf_opt=reorder_23(10**NSX[:,2].reshape(size[0],size[1],int(np.prod(size)/(size[0]*size[1]))),size)
        buf_opt=np.ravel(np.flip(reorder_buf_opt,2))

        self._hot_atmosphere = (logT_opt, logg_opt, _mu_opt, logE_opt, buf_opt)


photosphere = CustomPhotosphere_num(hot = hot, elsewhere = None,
                                values=dict(mode_frequency = spacetime['frequency']))

from time import time
start_original = time()

photosphere.hot_atmosphere = 'model_data/nsx_H_v200804.out'
# #nsx_H_v171019.out'

done_original = time()
print('Time taken by photosphere loading =', done_original-start_original)

photosphere['mode_frequency'] == spacetime['frequency']

star = xpsi.Star(spacetime = spacetime, photospheres = photosphere)

print("Parameters of the star:")
print(star.params)

p = [1.0368513939430604,
     6.087862992320039,
     0.26870812456714116,
     0.39140510783272897,
     0.04346870860640872,
     0.8002010406881243,
     1.1165398710637626,
     5.865655057483478,
     0.07360477761463673,
     2.4602238829718432,
     0.4277092192054918]
star(p)
star.update()

#To get the incident signal before interstellar absorption or operating with the telescope:
energies = np.logspace(-1.0, np.log10(3.0), 128, base=10.0)
photosphere.integrate(energies, threads=1) # the number of OpenMP threads to use
print("Bolometric pulse for 1st spot:")
print(repr(np.sum(photosphere.signal[0][0], axis=0)))
print("Bolometric pulse for 2nd spot:")
print(repr(np.sum(photosphere.signal[1][0], axis=0)))

from xpsi.likelihoods.default_background_marginalisation import eval_marginal_likelihood
from xpsi.likelihoods.default_background_marginalisation import precomputation

class CustomSignal(xpsi.Signal):
    """ A custom calculation of the logarithm of the likelihood.

    We extend the :class:`xpsi.Signal.Signal` class to make it callable.

    We overwrite the body of the __call__ method. The docstring for the
    abstract method is copied.

    """

    def __init__(self, workspace_intervals = 1000, epsabs = 0, epsrel = 1.0e-8,
                 epsilon = 1.0e-3, sigmas = 10.0, support = None, *args, **kwargs):
        """ Perform precomputation. """

        super(CustomSignal, self).__init__(*args, **kwargs)

        try:
            self._precomp = precomputation(self._data.counts.astype(np.int32))
        except AttributeError:
            print('No data... can synthesise data but cannot evaluate a '
                  'likelihood function.')
        else:
            self._workspace_intervals = workspace_intervals
            self._epsabs = epsabs
            self._epsrel = epsrel
            self._epsilon = epsilon
            self._sigmas = sigmas

            if support is not None:
                self._support = support
            else:
                self._support = -1.0 * np.ones((self._data.counts.shape[0],2))
                self._support[:,0] = 0.0

    @property
    def support(self):
        return self._support

    @support.setter
    def support(self, obj):
        self._support = obj

    def __call__(self, *args, **kwargs):
        self.loglikelihood, self.expected_counts, self.background_signal, self.background_signal_given_support = \
                eval_marginal_likelihood(self._data.exposure_time,
                                          self._data.phases,
                                          self._data.counts,
                                          self._signals,
                                          self._phases,
                                          self._shifts,
                                          self._precomp,
                                          self._support,
                                          self._workspace_intervals,
                                          self._epsabs,
                                          self._epsrel,
                                          self._epsilon,
                                          self._sigmas,
                                          kwargs.get('llzero'))


signal = CustomSignal(data = data,
                        instrument = NICER,
                        background = None,
                        interstellar = None,
                        workspace_intervals = 1000,
                        cache = True,
                        epsrel = 1.0e-8,
                        epsilon = 1.0e-3,
                        sigmas = 10.0,
                        support = None)

from scipy.stats import truncnorm
class CustomPrior(xpsi.Prior):
    """ A custom (joint) prior distribution.

    Source: Fictitious
    Model variant: ST-U
        Two single-temperature, simply-connected circular hot regions with
        unshared parameters.

    """

    __derived_names__ = ['compactness', 'phase_separation',]
    __draws_from_support__ = 3

    def __init__(self):
        """ Nothing to be done.

        A direct reference to the spacetime object could be put here
        for use in __call__:

        .. code-block::

            self.spacetime = ref

        Instead we get a reference to the spacetime object through the
        a reference to a likelihood object which encapsulates a
        reference to the spacetime object.

        """
        super(CustomPrior, self).__init__() # not strictly required if no hyperparameters

    def __call__(self, p = None):
        """ Evaluate distribution at ``p``.

        :param list p: Model parameter values.

        :returns: Logarithm of the distribution evaluated at ``p``.

        """
        temp = super(CustomPrior, self).__call__(p)
        if not np.isfinite(temp):
            return temp

        # based on contemporary EOS theory
        if not self.parameters['radius'] <= 16.0:
            return -np.inf

        ref = self.parameters.star.spacetime # shortcut

        # limit polar radius to try to exclude deflections >= \pi radians
        # due to oblateness this does not quite eliminate all configurations
        # with deflections >= \pi radians
        R_p = 1.0 + ref.epsilon * (-0.788 + 1.030 * ref.zeta)
        if R_p < 1.76 / ref.R_r_s:
            return -np.inf

        # polar radius at photon sphere for ~static star (static ambient spacetime)
        #if R_p < 1.5 / ref.R_r_s:
        #    return -np.inf

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

        phi = (ref['p__phase_shift'] - 0.5 - ref['s__phase_shift']) * _2pi

        ang_sep = xpsi.HotRegion.psi(ref['s__super_colatitude'],
                                     phi,
                                     ref['p__super_colatitude'])

        # hot regions cannot overlap
        if ang_sep < ref['p__super_radius'] + ref['s__super_radius']:
            return -np.inf

	#print("Calling CustomPrior with these paremeters:",p)

        return 0.0

    def inverse_sample(self, hypercube=None):
        """ Draw sample uniformly from the distribution via inverse sampling. """

        to_cache = self.parameters.vector

        if hypercube is None:
            hypercube = np.random.rand(len(self))

        # the base method is useful, so to avoid writing that code again:
        _ = super(CustomPrior, self).inverse_sample(hypercube)

        ref = self.parameters # shortcut

        idx = ref.index('distance')
        ref['distance'] = truncnorm.ppf(hypercube[idx], -2.0, 7.0, loc=0.3, scale=0.1)

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

        # compactness ratio M/R_eq
        p += [gravradius(ref['mass']) / ref['radius']]

        # phase separation between hot regions
        # first some temporary variables:
        if ref['p__phase_shift'] < 0.0:
            temp_p = ref['p__phase_shift'] + 1.0
        else:
            temp_p = ref['p__phase_shift']

        temp_s = 0.5 + ref['s__phase_shift']

        if temp_s > 1.0:
            temp_s = temp_s - 1.0

        # now append:
        if temp_s >= temp_p:
            p += [temp_s - temp_p]
        else:
            p += [1.0 - temp_p + temp_s]

        return p

prior = CustomPrior()


likelihood = xpsi.Likelihood(star = star, signals = signal,
                             num_energies=128,
                             threads=1,
                             prior=prior,
                             externally_updated=True)


wrapped_params = [0]*len(likelihood)
wrapped_params[likelihood.index('p__phase_shift')] = 1
wrapped_params[likelihood.index('s__phase_shift')] = 1

try:
    os.makedirs("run")
except OSError:
    if not os.path.isdir("run"):
        raise

runtime_params = {'resume': False,
                  'importance_nested_sampling': False,
                  'multimodal': False,
                  'n_clustering_params': None,
                  'outputfiles_basename': './run/run_Num',
                  'n_iter_before_update': 50,
                  'n_live_points': 50,
                  'sampling_efficiency': 0.8,
                  'const_efficiency_mode': False,
                  'wrapped_params': wrapped_params,
                  'evidence_tolerance': 0.5,
                  'seed': 7,
                  'max_iter': 100, # manual termination condition for short test
                  'verbose': True}

# let's require that checks pass before starting to sample
try:
	true_logl = -68147.0113542
	#print(likelihood(p))#Need to print this if not using force_update in the following line.
	likelihood.check(None, [true_logl], 1.0e-6,physical_points=[p],force_update=True)
except:
	print("Likelihood check did not pass. Checking if wrong atmosphere model installed.")
	true_logl = -116504.074
	#print(likelihood(p))#Need to print this if not using force_update in the following line.
	try:
		likelihood.check(None, [true_logl], 1.0e-6,physical_points=[p],force_update=True)
		print("Seems that blacbkody atmosphere extension was used instead of numerical.")
		print("Please re-install X-PSI using numerical atmosphere extension if want to use this test run.")
	except:
		print("Seems that neither of the likelihood checks passed, so something must be wrong.")
		exit()
	exit()

if __name__ == '__main__': # sample from the posterior
    xpsi.Sample.nested(likelihood, prior,**runtime_params)
