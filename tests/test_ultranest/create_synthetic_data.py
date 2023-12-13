""" Script to automatically create new synthetic data for a simple
ST Model (single temperature). Input parameters for the star and instrument
can be adjusted. 

"""

import argparse
import numpy as np
import math
import os

import xpsi
from xpsi.likelihoods.default_background_marginalisation import eval_marginal_likelihood
from xpsi.likelihoods.default_background_marginalisation import precomputation
from xpsi.global_imports import _2pi, gravradius
from xpsi.tools.synthesise import synthesise_exposure as _synthesise


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
    def from_response_files(cls, ARF, RMF, max_input, min_input=0,channel=[1,1500],
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

        matrix = np.ascontiguousarray(RMF[min_input:max_input,channel[0]:channel[1]].T, dtype=np.double)

        edges = np.zeros(ARF[min_input:max_input,2].shape[0]+1, dtype=np.double)
        
        

        edges[0] = ARF[min_input,0]; edges[1:] = ARF[min_input:max_input,1]

        for i in range(matrix.shape[0]):
            matrix[i,:] *= ARF[min_input:max_input,2]
    

        channels = np.arange(channel[0],channel[1])
    

        return cls(matrix, edges, channels, channel_edges[channel[0]:channel[1]+1,1])


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

class CustomPhotosphere(xpsi.Photosphere):
    """ Implement method for imaging."""

    @property
    def global_variables(self):

        return np.array([self['hot__super_colatitude'],
                        self['hot__phase_shift'] * _2pi,
                        self['hot__super_radius'],
                        self['hot__super_temperature']])
    
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
    
class CustomBackground(xpsi.Background):
    """ The background injected to generate synthetic data. """

    def __init__(self, bounds=None, value=None):

        # first the parameters that are fundemental to this class
        doc = """
        Powerlaw spectral index.
        """
        index = xpsi.Parameter('powerlaw_index',
                                strict_bounds = (1., 3.),
                                bounds = bounds,
                                doc = doc,
                                symbol = r'$\Gamma$',
                                value = value)

        super(CustomBackground, self).__init__(index)

    def __call__(self, energy_edges, phases):
        """ Evaluate the incident background field. """

        G = self['powerlaw_index']

        temp = np.zeros((energy_edges.shape[0] - 1, phases.shape[0]))

        temp[:,0] = (energy_edges[1:]**(G + 1.0) - energy_edges[:-1]**(G + 1.0)) / (G + 1.0)

        for i in range(phases.shape[0]):
            temp[:,i] = temp[:,0]

        self._incident_background= temp

class SynthesiseData(xpsi.Data):
    """ Custom data container to enable synthesis. """

    def __init__(self, channels, phases, first, last):

        self.channels = channels
        self._phases = phases

        try:
            self._first = int(first)
            self._last = int(last)
        except TypeError:
            raise TypeError('The first and last channels must be integers.')
        if self._first >= self._last:
            raise ValueError('The first channel number must be lower than the '
                             'the last channel number.')

def create_instrument():
    """ Create fake custom instrument object. """

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

    Instrument = CustomInstrument.from_response_files(ARF = ARF,
                                                    RMF = RMF,
                                                    max_input = 301,
                                                    min_input = 10,
                                                    channel=[10,301],
                                                    channel_edges =channel_edges) 

    return Instrument     

def create_star():
    """ Create a star object. """

    bounds = dict(distance = (0.5,2.),                  # (Earth) distance
                mass = (1.0, 3.0),                      # mass
                radius = (10., 15.0),                   # equatorial radius
                cos_inclination = (0., 1.))             # (Earth) inclination to rotation axis

    spacetime = xpsi.Spacetime(bounds=bounds, values=dict(frequency=314.0))# Fixing the spin

    bounds = dict(super_colatitude = (0.001, math.pi/2 - 0.001),
              super_radius = (0.001, math.pi/2 - 0.001),
              phase_shift = (-0.25, 0.75),
              super_temperature = (6.5, 7.2))

    # a simple circular, simply-connected spot
    hot_spot = xpsi.HotRegion(bounds=bounds,
                            values={}, # no initial values and no derived/fixed
                            symmetry=True,
                            omit=False,
                            cede=False,
                            concentric=False,
                            sqrt_num_cells=32,
                            min_sqrt_num_cells=10,
                            max_sqrt_num_cells=64,
                            num_leaves=100,
                            num_rays=200,
                            is_antiphased=True, 
                            prefix='hot') # unique prefix needed because >1 instance
      
    photosphere = CustomPhotosphere(hot = hot_spot, elsewhere = None,
                                values=dict(mode_frequency = spacetime['frequency']))
    
    for h in hot_spot.objects:
        h.set_phases(num_leaves = 100)

    star = xpsi.Star(spacetime = spacetime, photospheres = photosphere)

    return star 


def synthesise(self,
               exposure_time,
               expected_background_counts,
               name='synthetic',
               directory='./',
               **kwargs):
    
        """ Synthesise data set.

        """
        self._expected_counts, synthetic, bkg= _synthesise(exposure_time,
                                                             self._data.phases,
                                                             self._signals,
                                                             self._phases,
                                                             self._shifts,
                                                             expected_background_counts,
                                                             self._background.registered_background,
                                                             gsl_seed=42)
        
        try:
            if not os.path.isdir(directory):
                os.mkdir(directory)
        except OSError:
            print('Cannot create write directory.')
            raise

        np.savetxt(os.path.join(directory, name+'_realisation.dat'),
                   synthetic,
                   fmt = '%u')

        self._write(self.expected_counts,
                    filename = os.path.join(directory, name+'_expected_hreadable.dat'),
                    fmt = '%.8e')

        self._write(synthetic,
                    filename = os.path.join(directory, name+'_realisation_hreadable.dat'),
                    fmt = '%u')

def _write(self, counts, filename, fmt):
        """ Write to file in human readable format. """

        rows = len(self._data.phases) - 1
        rows *= len(self._data.channels)

        phases = self._data.phases[:-1]
        array = np.zeros((rows, 3))

        for i in range(counts.shape[0]):
            for j in range(counts.shape[1]):
                array[i*len(phases) + j,:] = self._data.channels[i], phases[j], counts[i,j]

            np.savetxt(filename, array, fmt=['%u', '%.6f'] + [fmt])


def create_synthetic_data(exposure_time, expected_background_counts, name, directory, p_T):
    """ Create synthetic data. """

    star = create_star()
    prior = CustomPrior()

    background = CustomBackground(bounds=(None, None))
    Instrument = create_instrument()
    _data = SynthesiseData(np.arange(10,301), np.linspace(0.0, 1.0, 33), 0, 290 )
    
    CustomSignal.synthesise = synthesise
    CustomSignal._write = _write

    signal = CustomSignal(data = _data,
                        instrument = Instrument,
                        background = background,
                        interstellar = None,
                        cache = True,
                        prefix='Instrument')

    likelihood = xpsi.Likelihood(star = star, signals = signal,
                        num_energies=128, 
                        threads=1,
                        externally_updated=False,
                        prior = prior)       

    print("Processing ...")
    instrument_kwargs = dict(exposure_time=exposure_time, 
                             expected_background_counts=expected_background_counts, 
                             name=name, 
                             directory=directory)
    likelihood.synthesise(p_T, force=True, Instrument=instrument_kwargs) 
    print("Done !")                      


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Creating synthetic data")
    
    parser.add_argument("-d", "--directory", action="store", required=False, dest="directory",
                        default="../synthetic_data/", help="Specify directory for output file\
                        (default: ../synthetic_data/)")

    parser.add_argument("-n", "--name", action="store", required=False, dest="name",
                        default="new_synthetic_data", help="Specify name for output file\
                        (default: new_synthetic_data)")
    
    parser.add_argument("-t", "--exptime", action="store", required=False, dest="exposure_time",
                        default=1000.0, type=float, help="Specify exposure time (default: 1000.0)")
    
    parser.add_argument("-b", "--bkgcounts", action="store", required=False, dest="expected_background_counts",
                        default=10000.0, type=float, help="Specify background counts (default: 10000.0)")

    parser.add_argument("-p", "--parameters", action="store", required=False, dest="p_T",
                        default=[1.4, 12, 1., math.cos(60*np.pi/180), 0.0, 70*np.pi/180, 0.75, 6.7, -2], 
                        type=float, nargs='+', help="Specify star model parameter values as floats like this:\
                        python create_synthetic_data.py -p 1.4 12 1. 0.5 etc. The parameters are: \
                        mass in solar radius, equatorial radius in km, distance in kpc,\
                        cosine of Earth inclination, rotation axis, phase shift,\
                        colatitude of the centre of the superseding region,\
                        angular radius of the (circular) superseding region,\
                        temperature in log 10, background of the spectral index gamma(E^gamma)")
    
    clargs = parser.parse_args()

    directory = clargs.directory 
    name = clargs.name
    exposure_time = clargs.exposure_time
    expected_background_counts = clargs.expected_background_counts
    p_T = clargs.p_T

    create_synthetic_data(exposure_time, expected_background_counts, name, directory, p_T)  