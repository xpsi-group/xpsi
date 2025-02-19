
import numpy as np
import math

import xpsi

np.random.seed(xpsi._rank+10)

import sys
import os

sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), "../../examples/examples_fast/Modules/"))

from CustomInstrument import CustomInstrument
from CustomSignal import CustomSignal
from CustomPhotosphere import CustomPhotosphere
from CustomPrior import CustomPrior

class TestMultiNestCheck(object):

    def test_sampling_works(self):
        # Data
        data_loaded = np.loadtxt(os.path.join(os.path.dirname(os.path.abspath(__file__)), "../../examples/examples_fast/Data/xpsi_good_realisation.dat"), dtype=np.double)

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
	                  mass = (1.0,2.0),
	                  radius = (5.0,16.0),
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
                                        is_antiphased=True,
	                                    image_order_limit=3, # up to tertiary
	                                    prefix='p')


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

        wrapped_params = [0] * len(likelihood)
        wrapped_params[likelihood.index('p__phase_shift')] = 1

        #The original (more accurate) run settings shown as commented.
        runtime_params = {'resume': False,
                          'importance_nested_sampling': False,
                          'multimodal': False,
                          'n_clustering_params': None,
                          'outputfiles_basename': '../../examples/examples_fast/Outputs/ST_live_1000_eff_0.3_seed0_v2',
                          'n_iter_before_update': 1, #100,
                          'n_live_points': 10, #1000,
                          'sampling_efficiency': 0.3,
                          'const_efficiency_mode': False,
                          'wrapped_params': wrapped_params,
                          'evidence_tolerance': 0.1,
                          'max_iter': 1, #-1,
                          'seed' : 0, # Fixing the seed
                          'LHS_seed': 42, # Fixing the LHS seed for hypercube fraction estimation
                          'verbose': True}


        prior.__draws_from_support__ = 1
        #Testing that the hypercube volume is estimated correctly:
        assert np.isclose(0.8333333333333334, prior.unit_hypercube_frac(LHS_seed=42), rtol=1.0e-5)

        #Testing sampler:  #To be uncommented once MultiNest installation added to Github workflow: 
        #xpsi.Sample.nested(likelihood, prior, **runtime_params)
        

