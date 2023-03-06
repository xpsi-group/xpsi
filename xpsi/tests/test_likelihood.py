
import numpy as np
import math

import xpsi

np.random.seed(xpsi._rank+10)

class CustomInstrument(xpsi.Instrument):

    """ Fake telescope response. """

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
    def from_response_files(cls, ARF, RMF, channel_edges,max_input,
                            min_input=0,channel=[1,1500],
                            ):
        """ Constructor which converts response files into :class:`numpy.ndarray`s.
        :param str ARF: Path to ARF which is compatible with
                                :func:`numpy.loadtxt`.
        :param str RMF: Path to RMF which is compatible with
                                :func:`numpy.loadtxt`.
        :param str channel_edges: Path to edges which is compatible with
                                  :func:`numpy.loadtxt`.
        """

        if min_input != 0:
            min_input = int(min_input)

        max_input = int(max_input)

        matrix = np.ascontiguousarray(RMF[min_input:max_input,channel[0]:channel[1]].T, dtype=np.double)

        edges = np.zeros(ARF[min_input:max_input,2].shape[0]+1, dtype=np.double)



        edges[0] = ARF[min_input,0]; edges[1:] = ARF[min_input:max_input,1]

        for i in range(matrix.shape[0]):
            matrix[i,:] *= ARF[min_input:max_input,2]


        channels = np.arange(channel[0],channel[1])


        return cls(matrix, edges, channels, channel_edges[channel[0]:channel[1]+1,1])

class CustomPhotosphere(xpsi.Photosphere):
    """ Implement method for imaging."""

    @property
    def global_variables(self):

        return np.array([self['h__super_colatitude'],
                          self['h__phase_shift'] * _2pi,
                          self['h__super_radius'],
                          self['h__super_temperature']])

from xpsi.likelihoods.default_background_marginalisation import eval_marginal_likelihood
from xpsi.likelihoods.default_background_marginalisation import precomputation

class CustomSignal(xpsi.Signal):
    """ A custom calculation of the logarithm of the NICER likelihood.

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
        self.loglikelihood, self.expected_counts, self.background_signal,self.background_given_support = \
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

from scipy.stats import truncnorm
from xpsi.global_imports import _G, _csq, _km, _2pi, gravradius, _dpr
from scipy.interpolate import Akima1DInterpolator
from scipy.stats import truncnorm

class CustomPrior(xpsi.Prior):
    """ A custom (joint) prior distribution.

    Source: Fictitious
    Model variant: ST-

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

        # Compactness limit
        R_p = 1.0 + ref.epsilon * (-0.788 + 1.030 * ref.zeta)
        if R_p < 1.76 / ref.R_r_s:
            return -np.inf

        mu = math.sqrt(-1.0 / (3.0 * ref.epsilon * (-0.788 + 1.030 * ref.zeta)))

        # 2-surface cross-section have a single maximum in |z|
        # i.e., an elliptical surface; minor effect on support, if any,
        # for high spin frequenies
        if mu < 1.0:
            return -np.inf

        ref = self.parameters

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
        ref['distance'] = truncnorm.ppf(hypercube[idx], -5.0, 5.0, loc=1.0, scale=0.1)

        # flat priors in cosine of hot region centre colatitudes (isotropy)
        # support modified by no-overlap rejection condition
        idx = ref.index('hot__super_colatitude')
        a, b = ref.get_param('hot__super_colatitude').bounds
        a = math.cos(a); b = math.cos(b)
        ref['hot__super_colatitude'] = math.acos(b + (a - b) * hypercube[idx])


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

        return p
  
class TestLikelihoodCheck(object):

    def test_likelihood_is_correct(self):
        # Data
        #if __name__ == '__main__':
	    #    data_path = "../Data/xpsi_good_realisation.dat"
        #else:
	    #    data_path = "./Data/xpsi_good_realisation.dat"

        #try:
	    #    data_loaded = np.loadtxt(data_path, dtype=np.double)
        #except:
	    #    print("Loading the data assuming the notebook was run for documentation pages")
	    #    data_loaded = np.loadtxt('../../examples/examples_fast/Data/xpsi_good_realisation.dat', dtype=np.double)

        try:
            data_loaded = np.loadtxt('../../examples/examples_fast/Data/xpsi_good_realisation.dat', dtype=np.double)
        except:
            data_loaded = np.loadtxt('./examples/examples_fast/Data/xpsi_good_realisation.dat', dtype=np.double)

        data = xpsi.Data(data_loaded,
	                         channels=np.arange(10,301),
	                         phases=np.linspace(0.0, 1.0, 33),
	                         first=0,
	                         last=290,
	                         exposure_time=1000.0)

        # # Instrument settings

        channel_number=np.arange(0,1501)    # The channel nnumber
        energy_low=np.arange(0,15.01, 0.01) # Lower bounds of each channel
        energy_high=energy_low+0.01         # Upper bounds of each channel
        channel_edges=np.array([list(channel_number),list(energy_low),list(energy_high)]).T

        # ARF
        arf_energy_low=[0.1]
        arf_energy_high=[0.105]
        arf_val=[1800]

        counter=1
        while arf_energy_low[-1]<=14.995:
	        arf_energy_low.append(arf_energy_low[-1]+0.005)
	        arf_energy_high.append(arf_energy_high[-1]+0.005)
	        arf_val.append(1800)
	        counter +=1


        ARF=np.array([list(arf_energy_low),
	                  list(arf_energy_high),
	                  list(arf_val)]).T

        # RMF
        RMF=np.diag(np.full(counter,1))


        Instrument = CustomInstrument.from_response_files(ARF =ARF,
	                                                 RMF = RMF,
	                                                 channel_edges =channel_edges,
	                                                 max_input = 301,
	                                                 min_input = 10,
	                                                 channel=[10,301])

        # # Signal
        signal = CustomSignal(data = data,
	                          instrument = Instrument,
	                          interstellar = None,
	                          cache = True,
	                          workspace_intervals = 1000,
	                          epsrel = 1.0e-8,
	                          epsilon = 1.0e-3,
	                          sigmas = 10.0)

        # # Space-time
        bounds = dict(distance = (0.5,2),
	                  mass = (1.0,1.6),
	                  radius = (10,13),
	                  cos_inclination = (0,1))


        spacetime = xpsi.Spacetime(bounds,
	                               values=dict(frequency = 314.0))

        # # Hot-spot
        bounds = dict(super_colatitude = (0.001, math.pi/2 - 0.001),
	                  super_radius = (0.001, math.pi/2 - 0.001),
	                  phase_shift = (-0.25, 0.75),
	                  super_temperature = (6., 7.))  # Valery model limit


        hot_spot = xpsi.HotRegion(bounds=bounds,
	                                    values={},
	                                    symmetry=True,
	                                    omit=False,
	                                    cede=False,
	                                    concentric=False,
	                                    sqrt_num_cells=32,
	                                    min_sqrt_num_cells=16,
	                                    max_sqrt_num_cells=64,
	                                    num_leaves=64,
	                                    num_rays=512,
	                                    is_secondary=True,
	                                    image_order_limit=3, # up to tertiary
	                                    prefix='hot')


        # # Photosphere
        photosphere = CustomPhotosphere(hot = hot_spot, elsewhere = None,
	                                    values=dict(mode_frequency = spacetime['frequency']))


        # # Star
        star = xpsi.Star(spacetime = spacetime, photospheres = photosphere)

        # # Prior
        prior = CustomPrior()


        # # Likelihood

        likelihood = xpsi.Likelihood(star = star, signals = signal,
	                                 num_energies = 128,
	                                 threads = 1,
	                                 externally_updated = True,
	                                 prior = prior)

        # Crucial step, if the likelihood check fails, then something went terrible wrong :)
        p=[1.4,10,1.,math.cos(60*np.pi/180),0.0,70*np.pi/180, 0.75,6.8]

        #likelihood.check(None, [-47881.27817666349], 1.0e-5, physical_points=[p])
        assert np.isclose(-47881.27817666349, likelihood(p,force=True,reinitialise=True), rtol=1.0e-5)
