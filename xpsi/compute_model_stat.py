
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

class namespace():
    pass

np.random.seed(10)

#Then we can use that data:

#settings = dict(counts = np.loadtxt('../examples/data_my/new_synthetic_realisation.dat', dtype=np.double),
settings = dict(counts = np.loadtxt('../docs/source/data/new_synthetic_realisation.dat', dtype=np.double),
                channels=np.arange(20,201),
                phases=np.linspace(0.0, 1.0, 33),
                first=0, last=180,
                exposure_time=984307.6661)
                
counts = np.loadtxt('../docs/source/data/new_synthetic_realisation.dat', dtype=np.double)
#print(counts.shape[0],counts.shape[1])

data = xpsi.Data(**settings)



IXPE_du1_I = namespace()
IXPE_du1_Q = namespace()
IXPE_du1_U = namespace()

IXPE_du2_I = namespace()
IXPE_du2_Q = namespace()
IXPE_du2_U = namespace()

IXPE_du3_I = namespace()
IXPE_du3_Q = namespace()
IXPE_du3_U = namespace()

minCH_IXPE = 0
maxCH_IXPE = 1
exposure_time_IXPE = 1.0


from ixpe_read import readData_pcube
fname = "/home/tuomo/polcslab/X-PATAP/x-patap/ad_new_simulations/toy_amsp_hotspot_direct_du1"
phase_IXPE, Idat, qn, un, Iderr_du1, qnerr_du1, unerr_du1, keVdat = readData_pcube(fname)

IXPE_du1_I.data = xpsi.Data(Idat,
                       channels=np.arange(minCH_IXPE, maxCH_IXPE),
                       phases=phase_IXPE,
                       first=0,
                       last=maxCH_IXPE-minCH_IXPE-1,
                       exposure_time=exposure_time_IXPE)

IXPE_du1_Q.data = xpsi.Data(qn,
                       channels=np.arange(minCH_IXPE, maxCH_IXPE),
                       phases=phase_IXPE,
                       first=0,
                       last=maxCH_IXPE-minCH_IXPE-1,
                       exposure_time=exposure_time_IXPE)
                       
IXPE_du1_U.data = xpsi.Data(un,
                       channels=np.arange(minCH_IXPE, maxCH_IXPE),
                       phases=phase_IXPE,
                       first=0,
                       last=maxCH_IXPE-minCH_IXPE-1,
                       exposure_time=exposure_time_IXPE)                                              

#Consider also defining this errors property in a CustomData or directly in Data.py
IXPE_du1_I.data.errors, IXPE_du1_Q.data.errors, IXPE_du1_U.data.errors = Iderr_du1, qnerr_du1, unerr_du1

fname = "/home/tuomo/polcslab/X-PATAP/x-patap/ad_new_simulations/toy_amsp_hotspot_direct_du2"
phase_IXPE, Idat, qn, un, Iderr_du2, qnerr_du2, unerr_du2, keVdat = readData_pcube(fname)

#print(Idat.shape[0],Idat.shape[1])
#exit()

IXPE_du2_I.data = xpsi.Data(Idat,
                       channels=np.arange(minCH_IXPE, maxCH_IXPE),
                       phases=phase_IXPE,
                       first=0,
                       last=maxCH_IXPE-minCH_IXPE-1,
                       exposure_time=exposure_time_IXPE)

IXPE_du2_Q.data = xpsi.Data(qn,
                       channels=np.arange(minCH_IXPE, maxCH_IXPE),
                       phases=phase_IXPE,
                       first=0,
                       last=maxCH_IXPE-minCH_IXPE-1,
                       exposure_time=exposure_time_IXPE)
                       
IXPE_du2_U.data = xpsi.Data(un,
                       channels=np.arange(minCH_IXPE, maxCH_IXPE),
                       phases=phase_IXPE,
                       first=0,
                       last=maxCH_IXPE-minCH_IXPE-1,
                       exposure_time=exposure_time_IXPE)  

IXPE_du2_I.data.errors, IXPE_du2_Q.data.errors, IXPE_du2_U.data.errors = Iderr_du2, qnerr_du2, unerr_du2
                       
fname = "/home/tuomo/polcslab/X-PATAP/x-patap/ad_new_simulations/toy_amsp_hotspot_direct_du3"
phase_IXPE, Idat, qn, un, Iderr_du3, qnerr_du3, unerr_du3, keVdat = readData_pcube(fname)

IXPE_du3_I.data = xpsi.Data(Idat,
                       channels=np.arange(minCH_IXPE, maxCH_IXPE),
                       phases=phase_IXPE,
                       first=0,
                       last=maxCH_IXPE-minCH_IXPE-1,
                       exposure_time=exposure_time_IXPE)

IXPE_du3_Q.data = xpsi.Data(qn,
                       channels=np.arange(minCH_IXPE, maxCH_IXPE),
                       phases=phase_IXPE,
                       first=0,
                       last=maxCH_IXPE-minCH_IXPE-1,
                       exposure_time=exposure_time_IXPE)
                       
IXPE_du3_U.data = xpsi.Data(un,
                       channels=np.arange(minCH_IXPE, maxCH_IXPE),
                       phases=phase_IXPE,
                       first=0,
                       last=maxCH_IXPE-minCH_IXPE-1,
                       exposure_time=exposure_time_IXPE)  

IXPE_du3_I.data.errors, IXPE_du3_Q.data.errors, IXPE_du3_U.data.errors = Iderr_du3, qnerr_du3, unerr_du3

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

def plot_one_pulse_stokes(pulse, x, label=r'Counts', cmap=cm.magma, vmin=None, vmax=None, sname="X",dchan=data.channels,errs=[0]):
    """ Plot a pulse resolved over a single rotational cycle. """

    fig = plt.figure(figsize = (7,7))

    gs = gridspec.GridSpec(1, 2, width_ratios=[50,1])
    ax = plt.subplot(gs[0])
    if errs[0][0]==0:
    	ax.plot(x, pulse[0]/np.max(pulse[0]), '-', color='k', lw=0.5) 
    	#ax.plot(x, pulse[0], '-', color='k', lw=0.5)
    else:
    	ax.errorbar(x, pulse[0], yerr=errs[0], xerr=0.0, fmt='o', color="purple",capsize=2.0,markersize=3.0)

    ax.set_xlim([0.0, 1.0])
    ax.set_ylabel(r'Normalized counts')
    ax.set_ylabel(r'Counts')    
    ax.set_xlabel(r'Phase')

    veneer((0.05, 0.2), (None, None), ax)

    #plt.subplots_adjust(wspace = 0.025)
    fig.savefig("figs/data"+sname+".pdf")


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


#mrf = "essentially the product of the effective area times the modulation factor, and is meant to be used to fit polarimetric models in XSPEC"
                             
                                             
IXPE_du1 = CustomInstrument_stokes.from_response_files(MRF = '/home/tuomo/polcslab/ixpe_sim/ixpeobssim_official/ixpeobssim/ixpeobssim/caldb/bcf/mrf/ixpemcdu1stdcutv006.mrf',
                                             RMF = '/home/tuomo/polcslab/ixpe_sim/ixpeobssim_official/ixpeobssim/ixpeobssim/caldb/cpf/rmf/ixpemcdu1stdcutv006.rmf',
                                             max_input = 175,
                                             max_channel = 200,
                                             min_input = 75,
                                             min_channel = 100,
                                             channel_edges = None)
                                             
IXPE_du2 = CustomInstrument_stokes.from_response_files(MRF = '/home/tuomo/polcslab/ixpe_sim/ixpeobssim_official/ixpeobssim/ixpeobssim/caldb/bcf/mrf/ixpemcdu2stdcutv006.mrf',
                                             RMF = '/home/tuomo/polcslab/ixpe_sim/ixpeobssim_official/ixpeobssim/ixpeobssim/caldb/cpf/rmf/ixpemcdu2stdcutv006.rmf',
                                             max_input = 175,
                                             max_channel = 200,
                                             min_input = 75,
                                             min_channel = 100,
                                             channel_edges = None)
                                             
IXPE_du3 = CustomInstrument_stokes.from_response_files(MRF = '/home/tuomo/polcslab/ixpe_sim/ixpeobssim_official/ixpeobssim/ixpeobssim/caldb/bcf/mrf/ixpemcdu3stdcutv006.mrf',
                                             RMF = '/home/tuomo/polcslab/ixpe_sim/ixpeobssim_official/ixpeobssim/ixpeobssim/caldb/cpf/rmf/ixpemcdu3stdcutv006.rmf',
                                             max_input = 175,
                                             max_channel = 200,
                                             min_input = 75,
                                             min_channel = 100,
                                             channel_edges = None)
                                             
plot_response_and_data=False
if plot_response_and_data:

	all_idus = ["1","2","3"]
	all_dus = [IXPE_du1,IXPE_du2,IXPE_du3]
	
	for idu in range(len(all_dus)):
		fig = plt.figure(figsize = (14,7))

		ax = fig.add_subplot(111)

		_ = ax.imshow(all_dus[idu].matrix,
			      cmap = cm.viridis,
			      rasterized = True)

		ax.set_ylabel('Channel')
		_ = ax.set_xlabel('Energy interval')
		fig.savefig("figs/response_IXPE_du"+all_idus[idu]+"_X.pdf")

		fig = plt.figure(figsize = (7,7))

		ax = fig.add_subplot(111)

		ax.plot((all_dus[idu].energy_edges[:-1] + all_dus[idu].energy_edges[1:])/2.0, np.sum(all_dus[idu].matrix, axis=0), 'k-')

		ax.set_ylabel('Effective area [cm$^{-2}$]')
		_ = ax.set_xlabel('Energy [keV]')

		fig.savefig("figs/eff_area_IXPE_du"+all_idus[idu]+"_X.pdf")					

		stokes_set = ["I","Q","U"]
		for ist in range(len(stokes_set)):
			if stokes_set[ist] == "I":
				if all_idus[idu] == "1":
					counts = IXPE_du1_I.data.counts
					phases = IXPE_du1_I.data.phases
					channels = IXPE_du1_I.data.channels
					errs = Iderr_du1
				elif all_idus[idu] == "2":
					counts = IXPE_du2_I.data.counts
					phases = IXPE_du2_I.data.phases
					channels = IXPE_du2_I.data.channels
					errs = Iderr_du2
				elif all_idus[idu] == "3":
					counts = IXPE_du3_I.data.counts
					phases = IXPE_du3_I.data.phases
					channels = IXPE_du3_I.data.channels
					errs = Iderr_du3
			if stokes_set[ist] == "Q":
				if all_idus[idu] == "1":
					counts = IXPE_du1_Q.data.counts
					phases = IXPE_du1_Q.data.phases
					channels = IXPE_du1_Q.data.channels
					errs = qnerr_du1
				elif all_idus[idu] == "2":
					counts = IXPE_du2_Q.data.counts
					phases = IXPE_du2_Q.data.phases
					channels = IXPE_du2_Q.data.channels
					errs = qnerr_du2
				elif all_idus[idu] == "3":
					counts = IXPE_du3_Q.data.counts
					phases = IXPE_du3_Q.data.phases
					channels = IXPE_du3_Q.data.channels
					errs = qnerr_du3
			if stokes_set[ist] == "U":
				if all_idus[idu] == "1":
					counts = IXPE_du1_U.data.counts
					phases = IXPE_du1_U.data.phases
					channels = IXPE_du1_U.data.channels
					errs = unerr_du1
				elif all_idus[idu] == "2":
					counts = IXPE_du2_U.data.counts
					phases = IXPE_du2_U.data.phases
					channels = IXPE_du2_U.data.channels
					errs = unerr_du2
				elif all_idus[idu] == "3":
					counts = IXPE_du3_U.data.counts
					phases = IXPE_du3_U.data.phases
					channels = IXPE_du3_U.data.channels
					errs = unerr_du3
		
		
			plot_one_pulse_stokes(counts,phases,sname="_IXPE_du"+all_idus[idu]+"_"+stokes_set[ist]+"X",dchan=channels,errs=errs)


from xpsi.likelihoods.default_background_marginalisation import eval_marginal_likelihood
from xpsi.likelihoods.default_background_marginalisation import precomputation

from xpsi.likelihoods._gaussian_likelihood_given_background_IQU import gaussian_likelihood_given_background
#from xpsi.likelihoods._poisson_likelihood_given_background_IQU import poisson_likelihood_given_background
#from xpsi.likelihoods._poisson_likelihood_given_background import poisson_likelihood_given_background

class CustomSignal_poisson(xpsi.Signal):
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

        super(CustomSignal_poisson, self).__init__(**kwargs)

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
        if two_spots:
            anegI = (False, False)
            anegQU = (True, True)
        else:
            anegI = (False)
            anegQU = (True)
    	errors = self._data.counts * 0.1
        if self.isI:
             
            self.loglikelihood, self.expected_counts = \
                gaussian_likelihood_given_background(self._data.exposure_time,
                                          self._data.phases,
                                          self._data.counts,
                                          errors, #self._data.errors,
                                          self._signals,
                                          self._phases,
                                          self._shifts,
                                          background = self._background.registered_background,
                                          allow_negative=anegI)
                                          
        else: #For Q and U use possibly different likelihood evaluation function:
            #TBD: Convert signal to Q/I and U/I before calling likelihood function
            #(since data and errors are currently in that form)
            #signal_ffit = self._signals/ ..... ?
            
            signal_ffit = self._signals
            
            
            self.loglikelihood, self.expected_counts = \
                gaussian_likelihood_given_background(self._data.exposure_time,
                                          self._data.phases,
                                          self._data.counts,
                                          errors, #self._data.errors,
                                          signal_ffit,
                                          self._phases,
                                          self._shifts,
                                          background = self._background.registered_background, 
                                          allow_negative=anegQU)                                                                   
                                          
                                                                                    


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
                                          kwargs.get('llzero'),
                                          allow_negative=(False, False))

class CustomBackground(xpsi.Background):
    """ The background injected to generate synthetic data. """

    def __init__(self, bounds=None, value=None):

        # first the parameters that are fundemental to this class
        doc = """
        Powerlaw spectral index.
        """
        index = xpsi.Parameter('powerlaw_index',
                                strict_bounds = (-3.0, -1.01),
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
            temp[:,i] = temp[:,0]*0.0 #*0.0 if want to have zero bkg

        #self.background = temp
        self.incident_background = temp

#background = CustomBackground(bounds=(None, None)) # use strict bounds, but do not fix/derive

background = CustomBackground(value=-2.0)

signals = [[],]

include_I = True

#TBD: Initialize signals in a loop to avoid repetition.

if include_I:
    signal_du1 = CustomSignal_poisson(data = IXPE_du1_I.data, #data,
                        instrument = IXPE_du1, #NICER,
                        background = background,
                        interstellar = None,
                        workspace_intervals = 1000,
                        cache = True,
                        epsrel = 1.0e-8,
                        epsilon = 1.0e-3,
                        sigmas = 10.0,
                        support = None,
                        stokes="I")
    signals[0].append(signal_du1)

signalQ_du1 = CustomSignal_poisson(data = IXPE_du1_Q.data, 
                        instrument = IXPE_du1, 
                        background = background, 
                        interstellar = None,
                        workspace_intervals = 1000,
                        cache = True,
                        epsrel = 1.0e-8,
                        epsilon = 1.0e-3,
                        sigmas = 10.0,
                        support = None,
                        stokes="Q")

signals[0].append(signalQ_du1)                        
                        
signalU_du1 = CustomSignal_poisson(data = IXPE_du1_U.data, 
                        instrument = IXPE_du1, 
                        background = background, 
                        interstellar = None,
                        workspace_intervals = 1000,
                        cache = True,
                        epsrel = 1.0e-8,
                        epsilon = 1.0e-3,
                        sigmas = 10.0,
                        support = None,
                        stokes="U")                        
signals[0].append(signalU_du1)                        
                        

if include_I:
    signal_du2 = CustomSignal_poisson(data = IXPE_du2_I.data, 
                        instrument = IXPE_du2, 
                        background = background,
                        interstellar = None,
                        workspace_intervals = 1000,
                        cache = True,
                        epsrel = 1.0e-8,
                        epsilon = 1.0e-3,
                        sigmas = 10.0,
                        support = None,
                        stokes="I")
    signals[0].append(signal_du2)

signalQ_du2 = CustomSignal_poisson(data = IXPE_du2_Q.data, 
                        instrument = IXPE_du2, 
                        background = background, 
                        interstellar = None,
                        workspace_intervals = 1000,
                        cache = True,
                        epsrel = 1.0e-8,
                        epsilon = 1.0e-3,
                        sigmas = 10.0,
                        support = None,
                        stokes="Q")

signals[0].append(signalQ_du2)                        
                        
signalU_du2 = CustomSignal_poisson(data = IXPE_du2_U.data, 
                        instrument = IXPE_du2, 
                        background = background, 
                        interstellar = None,
                        workspace_intervals = 1000,
                        cache = True,
                        epsrel = 1.0e-8,
                        epsilon = 1.0e-3,
                        sigmas = 10.0,
                        support = None,
                        stokes="U")                        
signals[0].append(signalU_du2) 


if include_I:
    signal_du3 = CustomSignal_poisson(data = IXPE_du3_I.data, 
                        instrument = IXPE_du3, 
                        background = background,
                        interstellar = None,
                        workspace_intervals = 1000,
                        cache = True,
                        epsrel = 1.0e-8,
                        epsilon = 1.0e-3,
                        sigmas = 10.0,
                        support = None,
                        stokes="I")
    signals[0].append(signal_du3)

signalQ_du3 = CustomSignal_poisson(data = IXPE_du3_Q.data, 
                        instrument = IXPE_du3, 
                        background = background, 
                        interstellar = None,
                        workspace_intervals = 1000,
                        cache = True,
                        epsrel = 1.0e-8,
                        epsilon = 1.0e-3,
                        sigmas = 10.0,
                        support = None,
                        stokes="Q")

signals[0].append(signalQ_du3)                        
                        
signalU_du3 = CustomSignal_poisson(data = IXPE_du3_U.data, 
                        instrument = IXPE_du3, 
                        background = background, 
                        interstellar = None,
                        workspace_intervals = 1000,
                        cache = True,
                        epsrel = 1.0e-8,
                        epsilon = 1.0e-3,
                        sigmas = 10.0,
                        support = None,
                        stokes="U")                        
signals[0].append(signalU_du3) 




#spacetime = xpsi.Spacetime.fixed_spin(300.0)
#bounds = dict(distance = (0.1, 1.0),                     # (Earth) distance
#                mass = (1.0, 3.0),                       # mass
#                radius = (3.0 * gravradius(1.0), 16.0),  # equatorial radius
#                cos_inclination = (0.0, 1.0))      # (Earth) inclination to rotation axis
#values = dict(frequency=300.0)

bounds = dict(cos_inclination = (0.0, 1.0))# (Earth) inclination to rotation axis
values =  dict(frequency = 401.0,mass=1.4,radius=12.0,distance= 1.0)

spacetime = xpsi.Spacetime(bounds=bounds, values=values)

#bounds = dict(super_colatitude = (None, None),
#              super_radius = (None, None),
#              phase_shift = (-0.25, 0.75),
#              super_temperature = (None, None))
#values={}
              
bounds = dict(super_colatitude = (None, None),
              phase_shift = (-0.25, 0.75))

deg2rad = np.pi/180.0
tempkeV = 1.0
tempK = np.log10(tempkeV*11604525.00617)
print("tempK=",tempK)
values = {'super_radius': 1.0*deg2rad,'super_temperature': tempK}              

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
                            values=values, 
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

two_spots = False

if two_spots:
    #bounds['super_temperature'] = None # declare fixed/derived variable
    bounds = dict(super_colatitude = (None, None), phase_shift = (-0.25, 0.75))
    secondary = xpsi.HotRegion(bounds=bounds, # can otherwise use same bounds
	                      #values={'super_temperature': derive()}, # create a callable value
	                      values={'super_radius': 1.0*deg2rad,'super_temperature': derive()},
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
    #hot['p__super_temperature'] = 6.0 # equivalent to ``primary['super_temperature'] = 6.0``
else:
    hot = primary


class CustomPhotosphereBB(xpsi.Photosphere):
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

class CustomPhotosphere(xpsi.Photosphere):
    """ Implement method for imaging."""

    @xpsi.Photosphere.hot_atmosphere.setter    
    def hot_atmosphere(self, path):
        size=(22,281)
        NSX = np.loadtxt(path, dtype=np.double)
        _mu = np.zeros(size[0]) 
        logE = np.zeros(size[1])    

        reorder_buf = np.zeros(size)

        index = 0
        for k in range(reorder_buf.shape[1]):
            for l in range(reorder_buf.shape[0]):            
                logE[k] = NSX[index,0]
                _mu[reorder_buf.shape[0] - l - 1] = NSX[index,1]            
                reorder_buf[reorder_buf.shape[0] - l - 1,k] = 10.0**(NSX[index,2])                                
                index += 1

        buf = np.zeros(np.prod(reorder_buf.shape))

        bufdex = 0
        for k in range(reorder_buf.shape[0]):
            for l in range(reorder_buf.shape[1]):            
                buf[bufdex] = reorder_buf[k,l]; bufdex += 1
        self._hot_atmosphere = (_mu, logE, buf)
        
    @xpsi.Photosphere.hot_atmosphere_Q.setter    
    def hot_atmosphere_Q(self, path):
        size=(22,281)
        NSX = np.loadtxt(path, dtype=np.double)
        _mu = np.zeros(size[0]) 
        logE = np.zeros(size[1])    

        reorder_buf = np.zeros(size)

        index = 0
        for k in range(reorder_buf.shape[1]):
            for l in range(reorder_buf.shape[0]):            
                logE[k] = NSX[index,0]
                _mu[reorder_buf.shape[0] - l - 1] = NSX[index,1] 
                Qsign = NSX[index,3]           
                reorder_buf[reorder_buf.shape[0] - l - 1,k] = Qsign*10.0**(NSX[index,2])                                
                index += 1

        buf = np.zeros(np.prod(reorder_buf.shape))

        bufdex = 0
        for k in range(reorder_buf.shape[0]):
            for l in range(reorder_buf.shape[1]):            
                buf[bufdex] = reorder_buf[k,l]; bufdex += 1
        self._hot_atmosphere_Q  = (_mu, logE, buf)

    @property
    def global_variables(self):
        """ This method is needed if we also want to ivoke the image-plane signal simulator. """

        return np.array([self['p__super_colatitude'],
                          self['p__phase_shift'] * _2pi,
                          self['p__super_radius'],
                          self['p__super_temperature'],
                          self['s__super_colatitude'],
                          (self['s__phase_shift'] + 0.5) * _2pi,
                          self['s__super_radius'],
                          self.hot.objects[1]['s__super_temperature']])


numerical_atmos = True

if numerical_atmos:
    photosphere = CustomPhotosphere(hot = hot, elsewhere = None,
                                values=dict(mode_frequency = spacetime['frequency']))
    photosphere.hot_atmosphere = "/home/tuomo/polcslab/X-PATAP/x-patap/analysis/model/atmos_thomI_corr2.txt"
    photosphere.hot_atmosphere_Q = "/home/tuomo/polcslab/X-PATAP/x-patap/analysis/model/atmos_thomQ_corr2.txt"

else:

    photosphere = CustomPhotosphereBB(hot = hot, elsewhere = None,
                                values=dict(mode_frequency = spacetime['frequency']))

photosphere['mode_frequency'] == spacetime['frequency']

star = xpsi.Star(spacetime = spacetime, photospheres = photosphere)

likelihood = xpsi.Likelihood(star = star, signals = signals,
                             num_energies=128,
                             threads=1,
                             externally_updated=False)

xpsi.set_phase_interpolant('Akima')

#p = [1.4,
#     12.5,
#     0.2,
#     math.cos(1.25),
#     0.0,
#     1.0,
#     0.075,
#     6.2,
#     0.025,
#     math.pi - 1.0,
#     0.2]
     
#For fast geometry simulation using polarimetry
if two_spots:
    p = [math.cos(60.0*deg2rad),
        0.0,
        20.0*deg2rad,
        0.0,
        20.0*deg2rad]
else:
    p = [math.cos(60.0*deg2rad),
        0.0,
        20.0*deg2rad]   

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

true_logl = -4.11344901e+03

likelihood.check(None, [true_logl], 1.0e-6,
                 physical_points=[p])
print('time = %.3f s' % (time.time() - t))

print("likelihood.params=",likelihood.params)


#For a synthetic NICER signal: -2.67136137e+04 #-26713.6136777
#For 3 synthetic NICER signals : -8.01408041e+04
#For 3 Stokes signals convolved with NICER response: -2.78202535e+05
#1) As above but with poisson likelihood given fixed background: -9.75133213e+09 (negative Q and U not yet allowed)
   #a) As above but background=0.0, negative Q and U enabled in likelihood fit: 2.49114467e+07 
   #b) As above but with default likelihood for I (and poisson for Q and U): -6.49958273e+09
#As 1a) but gaussian likelihood for all I, Q, and U: -6.94509688e+05

#For synthetic IXPE_du1 data with re-binned response (and using non-matching parameter values): -1.37016697e+03
#As above but for numerical atmosphere: -1.35961526e+03
#As above but with matching parameter values (1 hot spot): -1.37016381e+03 , or assuming 2 spots: -1.37016303e+03
#As above (1-spot), but with all dus included (and still with artificial errors): -4.11344901e+03

#print("signal I (primary):")
#print(signals[0][0].signals[0])
#print("signal Q (primary):")
#print(signals[0][1].signals[0])
#print("signal U (primary):")
#print(signals[0][2].signals[0])


def plot_pulse_stokes():
    """ Plot hot region signals before and after telescope operation. """
    fig = plt.figure(figsize=(7,7))
    ax = fig.add_subplot(111)

    ax.set_ylabel('Signal [arbitrary normalisation]')
    ax.set_xlabel('Phase [cycles]')

    #print("photosphere.signalI:")
    #print(np.sum(photosphere.signal[0][0], axis=0))

    temp = np.sum(signals[0][0].signals[0], axis=0)
    I1s = temp
    ax.plot(signals[0][0].phases[0], temp/np.max(temp), '-', color='k', lw=0.5)
    if two_spots:
        temp = np.sum(signals[0][0].signals[1], axis=0)    
        I2s = temp
        ax.plot(signals[0][0].phases[1], temp/np.max(temp), '-', color='r', lw=0.5)
    
    temp = np.sum(photosphere.signal[0][0], axis=0)
    I1p = temp
    ax.plot(signals[0][0].phases[0], temp/np.max(temp), 'o-', color='k', lw=0.5, markersize=2)
    if two_spots:
        temp = np.sum(photosphere.signal[1][0], axis=0)    
        I2p = temp
        ax.plot(signals[0][0].phases[1], temp/np.max(temp), 'o-', color='r', lw=0.5, markersize=2)    

    temp = np.sum(signals[0][0].expected_counts, axis=0)
    Iexpect = temp
    data_phases = np.linspace(0.0, 1.0, 10) #(0.0, 1.0, 33)
    #ax.plot(data_phases[0:32], temp/np.max(temp), '--', color='k', lw=0.5)
    ax.plot(data_phases[0:9], temp/np.max(temp), '--', color='k', lw=0.5)    

    ax.errorbar(IXPE_du1_I.data.phases, IXPE_du1_I.data.counts[0]/(np.max(IXPE_du1_I.data.counts[0])), yerr=IXPE_du1_I.data.errors[0]/(np.max(IXPE_du1_I.data.counts[0])), xerr=0.0, fmt='o', color="purple",capsize=2.0,markersize=3.0)

    #veneer((0.05,0.2), (0.05,0.2), ax)
    fig.savefig("figs/signalsIX.pdf")
    
    fig = plt.figure(figsize=(7,7))
    ax = fig.add_subplot(111)

    ax.set_ylabel('Q/I')
    ax.set_xlabel('Phase [cycles]')

    Q1s = np.sum(signals[0][1].signals[0], axis=0)
    Q1sn = np.copy(Q1s)
    for ip in range(len(Q1sn)):
    	if(I1s[ip] > 1e-10):
    		Q1sn[ip] = Q1s[ip]/I1s[ip]
    	else:
    		Q1sn[ip] = 0.0
    ax.plot(signals[0][1].phases[0], Q1sn, '-', color='k', lw=0.5)
    
    if two_spots:
        Q2s = np.sum(signals[0][1].signals[1], axis=0)    
        Q2sn = np.copy(Q2s)
        for ip in range(len(Q2sn)):
            if(I2s[ip] > 1e-10):
                Q2sn[ip] = Q2s[ip]/I2s[ip]
            else:
                Q2sn[ip] = 0.0
        ax.plot(signals[0][1].phases[1], Q2sn, '-', color='r', lw=0.5)
    
    Q1p = np.sum(photosphere.signalQ[0][0], axis=0)
    Q1pn = np.copy(Q1p)
    for ip in range(len(Q1pn)):
    	if(I1p[ip] > 1e-10):
    		Q1pn[ip] = Q1p[ip]/I1p[ip]
    	else:
    		Q1pn[ip] = 0.0
    ax.plot(signalQ_du1.phases[0], Q1pn, 'o-', color='k', lw=0.5, markersize=2)
    
    if two_spots:
        Q2p = np.sum(photosphere.signalQ[1][0], axis=0)
        Q2pn = np.copy(Q2p)
        for ip in range(len(Q1pn)):
            if(I2p[ip] > 1e-10):
                Q2pn[ip] = Q2p[ip]/I2p[ip]
            else:
                Q2pn[ip] = 0.0
        ax.plot(signalQ_du1.phases[1], Q2pn, 'o-', color='r', lw=0.5, markersize=2)   
    
    Qexpect = np.sum(signals[0][1].expected_counts, axis=0)
    data_phases = np.linspace(0.0, 1.0, 10)
    Qexpectn = np.copy(Qexpect)
    for ip in range(len(Qexpectn)):
    	if(Iexpect[ip] > 1e-50):
    		Qexpectn[ip] = Qexpect[ip]/Iexpect[ip]
    	else:
    		Qexpectn[ip] = 0.0 
    ax.plot(data_phases[0:9], Qexpectn, '--', color='k', lw=0.5)
        
    ax.errorbar(IXPE_du1_Q.data.phases, IXPE_du1_Q.data.counts[0], yerr=IXPE_du1_Q.data.errors[0], xerr=0.0, fmt='o', color="purple",capsize=2.0,markersize=3.0)       
        
    veneer((0.05,0.2), (0.05,0.2), ax)
    fig.savefig("figs/signalsQX.pdf")
    
    fig = plt.figure(figsize=(7,7))
    ax = fig.add_subplot(111)

    ax.set_ylabel('U/I')
    ax.set_xlabel('Phase [cycles]')

    U1s = np.sum(signals[0][2].signals[0], axis=0)
    U1sn = np.copy(U1s)
    for ip in range(len(U1sn)):
    	if(I1s[ip] > 1e-10):
    		U1sn[ip] = U1s[ip]/I1s[ip]
    	else:
    		U1sn[ip] = 0.0
    ax.plot(signals[0][2].phases[0], U1sn, '-', color='k', lw=0.5)

    if two_spots:
        U2s = np.sum(signals[0][2].signals[1], axis=0)    
        U2sn = np.copy(U2s)
        for ip in range(len(U2sn)):
            if(I2s[ip] > 1e-10):
                U2sn[ip] = U2s[ip]/I2s[ip]
            else:
                U2sn[ip] = 0.0
        ax.plot(signals[0][2].phases[1], U2sn, '-', color='r', lw=0.5)
    
    U1p = np.sum(photosphere.signalU[0][0], axis=0)
    U1pn = np.copy(U1p)
    for ip in range(len(U1pn)):
    	if(I1p[ip] > 1e-10):
    		U1pn[ip] = U1p[ip]/I1p[ip]
    	else:
    		U1pn[ip] = 0.0  		
    ax.plot(signalU_du1.phases[0], U1pn, 'o-', color='k', lw=0.5, markersize=2)
    
    if two_spots:
        U2p = np.sum(photosphere.signalU[1][0], axis=0)
        U2pn = np.copy(U2p)
        for ip in range(len(U1pn)):
    	    if(I2p[ip] > 1e-10):
                U2pn[ip] = U2p[ip]/I2p[ip]
            else:
                U2pn[ip] = 0.0
        ax.plot(signalU_du1.phases[1], U2pn, 'o-', color='r', lw=0.5, markersize=2)   
    
    Uexpect = np.sum(signals[0][2].expected_counts, axis=0)
    data_phases = np.linspace(0.0, 1.0, 10)
    Uexpectn = np.copy(Uexpect)
    for ip in range(len(Uexpectn)):
    	if(Iexpect[ip] > 1e-10):
    		Uexpectn[ip] = Uexpect[ip]/Iexpect[ip]
    	else:
    		Uexpectn[ip] = 0.0 
    ax.plot(data_phases[0:9], Uexpectn, '--', color='k', lw=0.5)

    ax.errorbar(IXPE_du1_U.data.phases, IXPE_du1_U.data.counts[0], yerr=IXPE_du1_U.data.errors[0], xerr=0.0, fmt='o', color="purple",capsize=2.0,markersize=3.0)  

    veneer((0.05,0.2), (0.05,0.2), ax)
    fig.savefig("figs/signalsUX.pdf")       
    
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
_ = plot_pulse_stokes()

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

        ## based on contemporary EOS theory
        #if not self.parameters['radius'] <= 16.0:
        #    return -np.inf

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
        if two_spots:
            if ref['p__super_colatitude'] > ref['s__super_colatitude']:
                return -np.inf
            phi = (ref['p__phase_shift'] - 0.5 - ref['s__phase_shift']) * _2pi
            ang_sep = xpsi.HotRegion.psi(ref['s__super_colatitude'],
                                     phi,
                                     ref['p__super_colatitude'])
            # hot regions cannot overlap
            #Works only if super_radius is a free parameter... TBD: What do if not?
            #if ang_sep < ref['p__super_radius'] + ref['s__super_radius']:
            #    return -np.inf
        return 0.0

    def inverse_sample(self, hypercube=None):
        """ Draw sample uniformly from the distribution via inverse sampling. """

        to_cache = self.parameters.vector

        if hypercube is None:
            hypercube = np.random.rand(len(self))

        # the base method is useful, so to avoid writing that code again:
        _ = super(CustomPrior, self).inverse_sample(hypercube)

        ref = self.parameters # shortcut

        #idx = ref.index('distance')
        #ref['distance'] = truncnorm.ppf(hypercube[idx], -2.0, 7.0, loc=0.3, scale=0.1)

        # flat priors in cosine of hot region centre colatitudes (isotropy)
        # support modified by no-overlap rejection condition
        idx = ref.index('p__super_colatitude')
        a, b = ref.get_param('p__super_colatitude').bounds
        a = math.cos(a); b = math.cos(b)
        ref['p__super_colatitude'] = math.acos(b + (a - b) * hypercube[idx])

        if two_spots:
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
if two_spots:
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

# let's require that checks pass before starting to sample
check_kwargs = dict(hypercube_points = None,
                    physical_points = p, # externally_updated preserved
                    loglikelihood_call_vals = [true_logl], 
                    rtol_loglike = 1.0e-6) # choose a tolerance

# note that mutual refs are already stored in the likelihood and prior
# objects to facilitate communication externally of the sampling process
xpsi.Sample.nested(likelihood, prior, check_kwargs, **runtime_params)





