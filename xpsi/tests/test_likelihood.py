
import numpy as np
import math

import xpsi

np.random.seed(xpsi._rank+10)

import sys
sys.path.append('../../examples/examples_fast/Modules/')

from CustomInstrument import CustomInstrument
from CustomSignal import CustomSignal
from CustomPhotosphere import CustomPhotosphere
from CustomPrior import CustomPrior

class TestLikelihoodCheck(object):

    def test_likelihood_is_correct(self):
        # Data
        if __name__ == '__main__':
	        data_path = "../Data/xpsi_good_realisation.dat"
        else:
	        data_path = "./Data/xpsi_good_realisation.dat"

        try:
	        data_loaded = np.loadtxt(data_path, dtype=np.double)
        except:
	        print("Loading the data assuming the notebook was run for documentation pages")
	        data_loaded = np.loadtxt('../../examples/examples_fast/Data/xpsi_good_realisation.dat', dtype=np.double)

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
