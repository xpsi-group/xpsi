
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



#Then we can use that data:

#settings = dict(counts = np.loadtxt('../examples/data_my/new_synthetic_realisation.dat', dtype=np.double),
settings = dict(counts = np.loadtxt('../docs/source/data/new_synthetic_realisation.dat', dtype=np.double),
                channels=np.arange(20,201),
                phases=np.linspace(0.0, 1.0, 33),
                first=0, last=180,
                exposure_time=984307.6661)

data = xpsi.Data(**settings)

rcParams['text.usetex'] = False
rcParams['font.size'] = 14.0

def veneer(x, y, axes, lw=1.0, length=8):
    """ Make the plots a little more aesthetically pleasing. """
    if x is not None:
        if x[1] is not None:
            axes.xaxis.set_major_locator(MultipleLocator(x[1]))
        if x[0] is not None:
            axes.xaxis.set_minor_locator(MultipleLocator(x[0]))
    else:
        axes.xaxis.set_major_locator(AutoLocator())
        axes.xaxis.set_minor_locator(AutoMinorLocator())

    if y is not None:
        if y[1] is not None:
            axes.yaxis.set_major_locator(MultipleLocator(y[1]))
        if y[0] is not None:
            axes.yaxis.set_minor_locator(MultipleLocator(y[0]))
    else:
        axes.yaxis.set_major_locator(AutoLocator())
        axes.yaxis.set_minor_locator(AutoMinorLocator())

    axes.tick_params(which='major', colors='black', length=length, width=lw)
    axes.tick_params(which='minor', colors='black', length=int(length/2), width=lw)
    plt.setp(axes.spines.values(), linewidth=lw, color='black')

def plot_one_pulse(pulse, x, label=r'Counts', cmap=cm.magma, vmin=None, vmax=None):
    """ Plot a pulse resolved over a single rotational cycle. """

    fig = plt.figure(figsize = (7,7))

    gs = gridspec.GridSpec(1, 2, width_ratios=[50,1])
    ax = plt.subplot(gs[0])
    ax_cb = plt.subplot(gs[1])

    profile = ax.pcolormesh(x,
                             data.channels,
                             pulse,
                             vmin = vmin,
                             vmax = vmax,
                             cmap = cmap,
                             linewidth = 0,
                             rasterized = True)

    profile.set_edgecolor('face')

    ax.set_xlim([0.0, 1.0])
    ax.set_yscale('log')
    ax.set_ylabel(r'Channel')
    ax.set_xlabel(r'Phase')

    cb = plt.colorbar(profile,
                      cax = ax_cb)

    cb.set_label(label=label, labelpad=25)
    cb.solids.set_edgecolor('face')

    veneer((0.05, 0.2), (None, None), ax)

    plt.subplots_adjust(wspace = 0.025)
    fig.savefig("figs/dataX.pdf")



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

        if min_input != 0:
            min_input = int(min_input)

        max_input = int(max_input)

        try:
            ARF = np.loadtxt(ARF, dtype=np.double, skiprows=3)
            RMF = np.loadtxt(RMF, dtype=np.double)
            if channel_edges:
                channel_edges = np.loadtxt(channel_edges, dtype=np.double, skiprows=3)[:,1:]
        except:
            print('A file could not be loaded.')
            raise

        matrix = np.ascontiguousarray(RMF[min_input:max_input,20:201].T, dtype=np.double)

        edges = np.zeros(ARF[min_input:max_input,3].shape[0]+1, dtype=np.double)

        edges[0] = ARF[min_input,1]; edges[1:] = ARF[min_input:max_input,2]

        for i in range(matrix.shape[0]):
            matrix[i,:] *= ARF[min_input:max_input,3]

        channels = np.arange(20, 201)

        return cls(matrix, edges, channels, channel_edges[20:202,-2])

NICER = CustomInstrument.from_response_files(ARF = '../examples/model_data/nicer_v1.01_arf.txt',
                                             RMF = '../examples/model_data/nicer_v1.01_rmf_matrix.txt',
                                             max_input = 500,
                                             min_input = 0,
                                             channel_edges = '../examples/model_data/nicer_v1.01_rmf_energymap.txt')
plot_response_and_data=False
if plot_response_and_data:
	fig = plt.figure(figsize = (14,7))

	ax = fig.add_subplot(111)
	veneer((25, 100), (10, 50), ax)

	_ = ax.imshow(NICER.matrix,
		      cmap = cm.viridis,
		      rasterized = True)

	ax.set_ylabel('Channel $-\;20$')
	_ = ax.set_xlabel('Energy interval')
	fig.savefig("figs/responseX.pdf")

	fig = plt.figure(figsize = (7,7))

	ax = fig.add_subplot(111)
	veneer((0.1, 0.5), (50,250), ax)

	ax.plot((NICER.energy_edges[:-1] + NICER.energy_edges[1:])/2.0, np.sum(NICER.matrix, axis=0), 'k-')

	ax.set_ylabel('Effective area [cm$^{-2}$]')
	_ = ax.set_xlabel('Energy [keV]')

	fig.savefig("figs/eff_areaX.pdf")

	plot_one_pulse(data.counts, data.phases)


from xpsi.likelihoods.default_background_marginalisation import eval_marginal_likelihood
from xpsi.likelihoods.default_background_marginalisation import precomputation

class CustomSignal(xpsi.Signal):
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

        super(CustomSignal, self).__init__(**kwargs)

        try:
            self._precomp = precomputation(self._data.counts.astype(np.int32))
        except AttributeError:
            print('Warning: No data... can synthesise data but cannot evaluate a '
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

    def __call__(self, *args, **kwargs):
        self.loglikelihood, self.expected_counts, self.background_signal, extra_element = \
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
                                          kwargs.get('llzero'),
                                          allow_negative=(False, False))


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


spacetime = xpsi.Spacetime.fixed_spin(300.0)
#for p in spacetime:
#    print(p)
bounds = dict(distance = (0.1, 1.0),                     # (Earth) distance
                mass = (1.0, 3.0),                       # mass
                radius = (3.0 * gravradius(1.0), 16.0),  # equatorial radius
                cos_inclination = (0.0, 1.0))      # (Earth) inclination to rotation axis

spacetime = xpsi.Spacetime(bounds=bounds, values=dict(frequency=300.0))

bounds = dict(super_colatitude = (None, None),
              super_radius = (None, None),
              phase_shift = (-0.25, 0.75),
              super_temperature = (None, None))

ceding=False #True

if ceding:
	bounds = dict(super_colatitude=(None,None),
		      super_radius = (None, None),
		      phase_shift = (0.0, 1.0),
		      super_temperature = (None, None),
		      cede_colatitude = (None, None),
		      cede_radius = (None, None),
		      cede_azimuth = (None, None),
		      cede_temperature = (None, None))

# a simple circular, simply-connected spot
primary = xpsi.HotRegion(bounds=bounds,
                            values={}, # no initial values and no derived/fixed
                            symmetry=True,
                            omit=False,
                            cede=ceding, #True, #False,
                            concentric=False,
                            sqrt_num_cells=32,
                            min_sqrt_num_cells=10,
                            max_sqrt_num_cells=64,
                            num_leaves=100,
                            num_rays=200,
                            prefix='p') # unique prefix needed because >1 instance

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
                              values={'super_temperature': derive()}, # create a callable value
                              symmetry=True,
                              omit=False,
                              cede=ceding, #True, #False,
                              concentric=False,
                              sqrt_num_cells=32,
                              min_sqrt_num_cells=10,
                              max_sqrt_num_cells=100,
                              num_leaves=100,
                              num_rays=200,
                              is_antiphased=True,
                              prefix='s') # unique prefix needed because >1 instance

from xpsi import HotRegions
hot = HotRegions((primary, secondary))
h = hot.objects[0]
#h.names
#h.get_param('phase_shift')
hot['p__super_temperature'] = 6.0 # equivalent to ``primary['super_temperature'] = 6.0``
#secondary['super_temperature']

class CustomPhotosphere(xpsi.Photosphere):
    """ Implement method for imaging."""

    @property
    def global_variables(self):

        return np.array([self['p__super_colatitude'],
                          self['p__phase_shift'] * _2pi,
                          self['p__super_radius'],
                          self['p__super_temperature'],
                          self['s__super_colatitude'],
                          (self['s__phase_shift'] + 0.5) * _2pi,
                          self['s__super_radius'],
                          self.hot.objects[1]['s__super_temperature']])

photosphere = CustomPhotosphere(hot = hot, elsewhere = None,
                                values=dict(mode_frequency = spacetime['frequency']))

photosphere['mode_frequency'] == spacetime['frequency']

star = xpsi.Star(spacetime = spacetime, photospheres = photosphere)

likelihood = xpsi.Likelihood(star = star, signals = signal,
                             num_energies=128,
                             threads=1,
                             externally_updated=False)

xpsi.set_phase_interpolant('Akima')

p = [1.4,
     12.5,
     0.2,
     math.cos(1.25),
     0.0,
     1.0,
     0.075,
     6.2,
     0.025,
     math.pi - 1.0,
     0.2]

print(star)

if ceding:
	p = [1.4, #mass
	     12.0, #radius
	     0.2, #distance
	     math.cos(1.25), #cos_inclination
	     0.0, #p__phase_shift
	     1.0, #p__super_colatitude
	     0.075, #p__super_radius
	     6.2, #p__super_temperature
	     0.1, #p__cede_colatitude
	     0.1, #p__cede_radius
	     0.0, #p__cede_azimuth
	     6.2, #p__cede_temperature
	     0.025, #s__phase_shift
	     math.pi - 1.0, #s__super_colatitude
	     0.2, #s__super_radius
	     math.pi-1.0, #s__cede_colatitude ..
	     0.3, #s__cede_radius
	     0.0, #s__cede_azimuth
	     6.2] #s__cede_temperature


likelihood.clear_cache()
t = time.time()
# source code changes since model was applied, so let's be a
# bit lenient when checking the likelihood function
likelihood.check(None, [-26713.6136777], 1.0e-6, #stokes=True,
                 physical_points=[p])
print('time = %.3f s' % (time.time() - t))

print("likelihood.params=",likelihood.params)


#TypeError: __call__() got an unexpected keyword argument 'stokes' !!!

#exit()

# > xpsi.set_phase_interpolant('Akima')
# Checking likelihood and prior evaluation before commencing sampling...
# Cannot import ``allclose`` function from NumPy.
# Using fallback implementation...
# Checking closeness of likelihood arrays:
# -2.67136012e+04 | -2.67136137e+04 .....
# Closeness evaluated.
# Log-likelihood value checks passed on root process.
# Checks passed.
# time = 0.571 s

# > xpsi.set_phase_interpolant('Steffen')
# Checking likelihood and prior evaluation before commencing sampling...
# Cannot import ``allclose`` function from NumPy.
# Using fallback implementation...
# Checking closeness of likelihood arrays:
# -2.67136140e+04 | -2.67136137e+04 .....
# Closeness evaluated.
# Log-likelihood value checks passed on root process.
# Checks passed.
# time = 0.581 s

# > xpsi.set_phase_interpolant('Cubic')
# Checking likelihood and prior evaluation before commencing sampling...
# Cannot import ``allclose`` function from NumPy.
# Using fallback implementation...
# Checking closeness of likelihood arrays:
# -2.67135656e+04 | -2.67136137e+04 .....
# Closeness evaluated.
# Log-likelihood value checks passed on root process.
# Checks passed. (increasing the tolerance to 1e-5)
# time = 0.720 s (consistently slower)



def plot_pulse():
    """ Plot hot region signals before and after telescope operation. """
    fig = plt.figure(figsize=(7,7))
    ax = fig.add_subplot(111)

    ax.set_ylabel('Signal [arbitrary normalisation]')
    ax.set_xlabel('Phase [cycles]')

    temp = np.sum(signal.signals[0], axis=0)
    ax.plot(signal.phases[0], temp/np.max(temp), '-', color='k', lw=0.5)
    temp = np.sum(signal.signals[1], axis=0)
    ax.plot(signal.phases[1], temp/np.max(temp), '-', color='r', lw=0.5)

    temp = np.sum(photosphere.signal[0][0], axis=0)
    ax.plot(signal.phases[0], temp/np.max(temp), 'o-', color='k', lw=0.5, markersize=2)
    temp = np.sum(photosphere.signal[1][0], axis=0)
    ax.plot(signal.phases[1], temp/np.max(temp), 'o-', color='r', lw=0.5, markersize=2)

    veneer((0.05,0.2), (0.05,0.2), ax)
    fig.savefig("figs/signalsX.pdf")


likelihood(p, reinitialise=False)
#_ = plot_pulse()

#The rest of the modelling examples and emcmc are still only in the notebook version.
#Let's still test multinest here:

from scipy.stats import truncnorm
class CustomPrior(xpsi.Prior):
    """ A custom (joint) prior distribution.

    Source: Fictitious
    Model variant: ST-U
        Two single-temperature, simply-connected circular hot regions with
        unshared parameters.

    """

    __derived_names__ = ['compactness', 'phase_separation',]

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

        # polar radius at photon sphere for ~static star (static ambient spacetime)
        #if R_p < 1.5 / ref.R_r_s:
        #    return -np.inf

        # limit polar radius to try to exclude deflections >= \pi radians
        # due to oblateness this does not quite eliminate all configurations
        # with deflections >= \pi radians
        R_p = 1.0 + ref.epsilon * (-0.788 + 1.030 * ref.zeta)
        if R_p < 1.76 / ref.R_r_s:
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

        phi = (ref['p__phase_shift'] - 0.5 - ref['s__phase_shift']) * _2pi

        ang_sep = xpsi.HotRegion.psi(ref['s__super_colatitude'],
                                     phi,
                                     ref['p__super_colatitude'])

        # hot regions cannot overlap
        if ang_sep < ref['p__super_radius'] + ref['s__super_radius']:
            return -np.inf

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
likelihood.prior = prior

wrapped_params = [0]*len(likelihood)
wrapped_params[likelihood.index('p__phase_shift')] = 1
wrapped_params[likelihood.index('s__phase_shift')] = 1

runtime_params = {'resume': False,
                  'importance_nested_sampling': False,
                  'multimodal': False,
                  'n_clustering_params': None,
                  'outputfiles_basename': './run/run', # make ./run directory manually
                  'n_iter_before_update': 50,
                  'n_live_points': 100,
                  'sampling_efficiency': 0.8,
                  'const_efficiency_mode': False,
                  'wrapped_params': wrapped_params,
                  'evidence_tolerance': 0.5,
                  'max_iter': 1000, # manual termination condition for short test
                  'verbose': True}

for h in hot.objects:
    h.set_phases(num_leaves = 100)

likelihood.threads = 3
likelihood.reinitialise()
likelihood.clear_cache()

# inform source code that parameter objects updated when inverse sampling
likelihood.externally_updated = True

p = [1.4,
     12.5,
     0.2,
     math.cos(1.25),
     0.0,
     1.0,
     0.075,
     6.2,
     0.025,
     math.pi - 1.0,
     0.2]

# let's require that checks pass before starting to sample
check_kwargs = dict(hypercube_points = None,
                    physical_points = p, # externally_updated preserved
                    loglikelihood_call_vals = [-26713.613677], # from above
                    rtol_loglike = 1.0e-6) # choose a tolerance

# note that mutual refs are already stored in the likelihood and prior
# objects to facilitate communication externally of the sampling process
xpsi.Sample.nested(likelihood, prior, check_kwargs, **runtime_params)





