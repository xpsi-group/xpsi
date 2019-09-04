from __future__ import print_function, division

import numpy as np
import math

import xpsi

print('Rank reporting: %d' % xpsi._rank)

from CustomData import CustomData
from CustomInstrument import CustomInstrument
from CustomPulse import CustomPulse
from CustomSpacetime import CustomSpacetime
from CustomPrior import CustomPrior
from CustomBackground import CustomBackground

data = CustomData.from_SWG('data/synthetic_realisation.dat', 984307.6661)

NICER = CustomInstrument.from_SWG(num_params=0,
                                  bounds=[],
                                  ARF = 'model_data/nicer_v1.01_arf.txt',
                                  RMF = 'model_data/nicer_v1.01_rmf_matrix.txt',
                                  max_input = 500,
                                  min_input = 0,
                                  chan_edges = 'model_data/nicer_v1.01_rmf_energymap.txt')

pulse = CustomPulse(tag = 'all',
                    num_params = 2,
                    bounds = [(-0.25, 0.75), (-0.25, 0.75)],
                    data = data,
                    instrument = NICER,
                    background = None,
                    interstellar = None,
                    energies_per_interval = 0.5,
                    default_energy_spacing = 'logspace',
                    fast_rel_energies_per_interval = 0.5,
                    workspace_intervals = 1000,
                    adaptive_energies = False,
                    store = False,
                    epsrel = 1.0e-8,
                    epsilon = 1.0e-3,
                    sigmas = 10.0)

from xpsi.global_imports import _c, _G, _M_s, _dpr, gravradius

bounds = [(0.1, 1.0),
          (1.0, 3.0),
          (3.0 * gravradius(1.0), 16.0),
          (0.001, math.pi/2.0)]

spacetime = CustomSpacetime(num_params = 4, bounds = bounds, S = 300.0)

bounds = [(0.001, math.pi - 0.001),
          (0.001, math.pi/2.0 - 0.001),
          (5.5, 6.5),
          (0.001, math.pi - 0.001),
          (0.001, math.pi/2.0 - 0.001),
          (5.5, 6.5)]

spot = xpsi.Spots(num_params=(3,3), bounds=bounds,
                    symmetry=True,
                    hole=False,
                    cede=False,
                    concentric=False,
                    antipodal_symmetry=False,
                    sqrt_num_cells=32,
                    min_sqrt_num_cells=10,
                    max_sqrt_num_cells=64,
                    do_fast=False,
                    num_leaves=100,
                    num_rays=200)

photosphere = xpsi.Photosphere(num_params = 0, bounds = [],
                               tag = 'all', spot = spot, elsewhere = None)

star = xpsi.Star(spacetime = spacetime, photospheres = photosphere)

likelihood = xpsi.Likelihood(star = star, pulses = pulse, threads=1)

prior = CustomPrior(bounds=likelihood.bounds, spacetime=spacetime)

likelihood.prior = prior

import time

p = [0.2,
     1.4,
     12.5,
     1.25,
     1.0,
     0.075,
     6.2,
     math.pi - 1.0,
     0.2,
     6.0,
     0.0,
     0.025]

t = time.time()
ll = likelihood(p) # OptiPlex: ll = -26713.613677693844
print('p: ', ll, time.time() - t)

runtime_params = {'resume': False,
                  'importance_nested_sampling': False,
                  'multimodal': False,
                  'n_clustering_params': None,
                  'outputfiles_basename': './run2_nlive1000_eff0.3_noCONST_noMM_noIS_tol-1',
                  'n_iter_before_update': 100,
                  'n_live_points': 1000,
                  'sampling_efficiency': 0.3,
                  'const_efficiency_mode': False,
                  'wrapped_params': [0,0,0,0,0,0,0,0,0,0,1,1],
                  'evidence_tolerance': 0.1,
                  'max_iter': -1,
                  'verbose': True}

xpsi.Sample.MultiNest(likelihood, prior, **runtime_params)
