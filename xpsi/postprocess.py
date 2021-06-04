from __future__ import division

import sys
import os

import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

import numpy as np
import math

from collections import OrderedDict

import xpsi

from xpsi import PostProcessing

# choose a seed for the notebook if you want caching to be useful
# and the notebook exactly reproducible
PostProcessing.set_random_seed(0)

from xpsi.global_imports import gravradius





#Code from Serena to test the convergence:
#data_rough = np.loadtxt('STPST_outputs/run1/run1_nlive1000_eff0.3_noCONST_noMM_noIS_tol-1.txt')
#samples = data_rough
#plt.plot(samples[np.argsort(-0.5*samples[:,1]),0][:])

import matplotlib.pyplot as plt
#data_rough = np.loadtxt('run/runphys_live.points')
data_rough = np.loadtxt('output1/run1_nlive1000_eff0.3_noCONST_noMM_noIS_tol-1.txt')
samples = data_rough
plt.plot(samples[np.argsort(-0.5*samples[:,1]),0][:])
plt.savefig('figs/posteriors/convergence_run_STS.pdf')

#data_rough = np.loadtxt('job1v3/run.txt')
data_rough = np.loadtxt('run_test1f/run.txt')
samples = data_rough
plt.plot(samples[np.argsort(-0.5*samples[:,1]),0][:])
plt.savefig('figs/posteriors/convergence_run_STUmytest1f.pdf')





#from STS_module import main as STS

#from STS_module_yves import main as STS
from STU_modules import main as STU
#from ST_PST_modules import main as ST_PST


print(STU.likelihood)


# names of free parameters ordered as in sample files
# NB: the parameter named cos_inclination in the sample files is
#     actually the inclination, but the name is declared here to
#     link it to the inclination parameter in the current API, which
#     is the cos(inclination). We do this because the inclinations
#     on disk will be transformed to cos(inclination) by the
#     CustomPrior.transform() instance method with the argument
#     old_API=True.
#STU.names = ['distance', 'mass', 'radius', 'cos_inclination',
#             'p__super_colatitude', 'p__super_radius', 'p__super_temperature',
#             's__super_colatitude', 's__super_radius', 's__super_temperature',
##             'column_density',
##             'alpha', 'beta', 'gamma',
#             'p__phase_shift', 's__phase_shift']
             
STU.names = ['mass', 'radius', 'distance', 'cos_inclination', 'p__phase_shift', 
'p__super_colatitude', 'p__super_radius', 'p__super_temperature', 's__phase_shift',
's__super_colatitude', 's__super_radius']


#For STS:
#Free parameters
#---------------
#mass: Gravitational mass [solar masses].
#radius: Coordinate equatorial radius [km].
#cos_inclination: Cosine of Earth inclination to rotation axis.
#p__phase_shift: The phase of the hot region, a periodic parameter [cycles].
#p__super_colatitude: The colatitude of the centre of the superseding region [radians].
#p__super_radius: The angular radius of the (circular) superseding region [radians].
#p__super_temperature: log10(superseding region effective temperature [K]).
#beta: Units of kpc^-2.
#column_density: Units of 10^20 cm^-2.


#STS.names = ['mass', 'radius', 'cos_inclination', 'p__phase_shift', 
#'p__super_colatitude', 'p__super_radius', 'p__super_temperature', 'beta', 'column_density']
#STS.names += ['compactness']#, 's__transformed_phase']


# names of derived variables of interest
#STU.names += ['compactness']#, 's__transformed_phase']

# the hard bounds imposed above for each parameter
# in some cases, depending on how the prior density
# is constructed, automated population of the
# dictionary below might not yield the required bounds,
# whereas manual population of the dictionary allows
# the user to take full control




#Free parameters
#---------------
#mass: Gravitational mass [solar masses].
#radius: Coordinate equatorial radius [km].
#distance: Earth distance [kpc].
#cos_inclination: Cosine of Earth inclination to rotation axis.
#p__phase_shift: The phase of the hot region, a periodic parameter [cycles].
#p__super_colatitude: The colatitude of the centre of the superseding region [radians].
#p__super_radius: The angular radius of the (circular) superseding region [radians].
#p__super_temperature: log10(superseding region effective temperature [K]).
#s__phase_shift: The phase of the hot region, a periodic parameter [cycles].
#s__super_colatitude: The colatitude of the centre of the superseding region [radians].
#s__super_radius: The angular radius of the (circular) superseding region [radians].




#Free parameters
#---------------
#mass: Gravitational mass [solar masses].
#radius: Coordinate equatorial radius [km].
#distance: Earth distance [kpc].
#cos_inclination: Cosine of Earth inclination to rotation axis.
#p__phase_shift: The phase of the hot region, a periodic parameter [cycles].
#p__super_colatitude: The colatitude of the centre of the superseding region [radians].
#p__super_radius: The angular radius of the (circular) superseding region [radians].
#p__super_temperature: log10(superseding region effective temperature [K]).
#s__phase_shift: The phase of the hot region, a periodic parameter [cycles].
#s__super_colatitude: The colatitude of the centre of the superseding region [radians].
#s__super_radius: The angular radius of the (circular) superseding region [radians].



STU.bounds = {'mass': (1.0, 3.0),
              'radius': (3.0 * gravradius(1.0), 16.0),
              #'distance': (0.235, 0.415),
              'distance': (0.14, 0.35),
              'cos_inclination': (0.0, math.cos(0.001)),
              'p__phase_shift': (-0.25,0.75),#(-0.5,0.5),
              'p__super_colatitude': (0.001, math.pi - 0.001),
              'p__super_radius': (0.001, math.pi/2.0 - 0.001),
              'p__super_temperature': (5.1, 6.8),
              's__phase_shift': (-0.25,0.75),#(-0.5,0.5),
              's__super_colatitude': (0.001, math.pi - 0.001),
              's__super_radius': (0.001, math.pi/2.0 - 0.001)}

#STS.bounds = {'mass': (1.0, 3.0),
#              'radius': (3.0 * gravradius(1.0), 16.0),
##              'distance': (0.235, 0.415),
#              'cos_inclination': (0.0, math.cos(0.001)),
#              'p__phase_shift': (-0.5,0.5),
#              'p__super_colatitude': (0.001, math.pi - 0.001),
#              'p__super_radius': (0.001, math.pi/2.0 - 0.001),
#              'p__super_temperature': (5.1, 6.8),
##              's__phase_shift': (-0.5,0.5),
##              's__super_colatitude': (0.001, math.pi - 0.001),
##              's__super_radius': (0.001, math.pi/2.0 - 0.001),
##              'beta': (0.0,1.0),
#              'beta': (0.0,30.0),
#              'column_density': (0.0, 5.0),
#              'compactness': (gravradius(1.0)/16.0, 1.0/3.0)}



# TeX compatible labels
#STU.labels = {'distance': r"D\;\mathrm{[kpc]}",
#              'mass': r"M\;\mathrm{[M}_{\odot}\mathrm{]}",
#              'radius': r"R_{\mathrm{eq}}\;\mathrm{[km]}",
#              'cos_inclination': r"\cos(i)",
#              'p__super_colatitude': r"\Theta_{p}\;\mathrm{[rad]}",
#              'p__super_radius': r"\zeta_{p}\;\mathrm{[rad]}",
#              'p__super_temperature': r"\mathrm{log10}(T_{p}\;[\mathrm{K}])",
#              's__super_colatitude': r"\Theta_{s}\;\mathrm{[rad]}",
#              's__super_radius': r"\zeta_{s}\;\mathrm{[rad]}",
##              's__super_temperature': r"\mathrm{log10}(T_{s}\;[\mathrm{K}])",
##              'column_density': r"N_{\mathrm{H}}\;\mathrm{[10^{20}\;cm^{-2}]}",
##              'alpha': r"\alpha",
##              'beta': r"\beta",
##              'gamma': r"\gamma",
#              'p__phase_shift': r"\phi_{p}\;\mathrm{[cycles]}",}
##              's__phase_shift': r"\phi_{s}\;\mathrm{[cycles]}",}
##              'compactness': r"M/R_{\mathrm{eq}}",
##              's__transformed_phase': r"\phi_{s}\;\mathrm{[cycles]}",}


STU.labels = {'mass': r"M\;\mathrm{[M}_{\odot}\mathrm{]}",
              'radius': r"R_{\mathrm{eq}}\;\mathrm{[km]}",
              'distance': r"D\;\mathrm{[kpc]}",
              'cos_inclination': r"\cos(i)",
              'p__phase_shift': r"\phi_{p}\;\mathrm{[cycles]}",
              'p__super_colatitude': r"\Theta_{p}\;\mathrm{[rad]}",
              'p__super_radius': r"\zeta_{p}\;\mathrm{[rad]}",
              'p__super_temperature': r"\mathrm{log10}(T_{p}\;[\mathrm{K}])",
              's__phase_shift': r"\phi_{p}\;\mathrm{[cycles]}",
              's__super_colatitude': r"\Theta_{s}\;\mathrm{[rad]}",
              's__super_radius': r"\zeta_{s}\;\mathrm{[rad]}"}


#STS.labels = {'mass': r"M\;\mathrm{[M}_{\odot}\mathrm{]}",
#              'radius': r"R_{\mathrm{eq}}\;\mathrm{[km]}",
##              'distance': r"D\;\mathrm{[kpc]}",
#              'cos_inclination': r"\cos(i)",
#              'p__phase_shift': r"\phi_{p}\;\mathrm{[cycles]}",
#              'p__super_colatitude': r"\Theta_{p}\;\mathrm{[rad]}",
#              'p__super_radius': r"\zeta_{p}\;\mathrm{[rad]}",
#              'p__super_temperature': r"\mathrm{log10}(T_{p}\;[\mathrm{K}])",
##              's__phase_shift': r"\phi_{p}\;\mathrm{[cycles]}",
##              's__super_colatitude': r"\Theta_{s}\;\mathrm{[rad]}",
##              's__super_radius': r"\zeta_{s}\;\mathrm{[rad]}",
#              'beta': r"\beta",
#              'column_density': r"N_{\mathrm{H}}\;\mathrm{[10^{20}\;cm^{-2}]}",
#              'compactness': r"M/R_{\mathrm{eq}}"}





getdist_kde_settings = {'ignore_rows': 0,
                         'min_weight_ratio': 1.0e-10,
                         'contours': [0.683, 0.954, 0.997],
                         'credible_interval_threshold': 0.001,
                         'range_ND_contour': 0,
                         'range_confidence': 0.001,
                         'fine_bins': 1024,
                         'smooth_scale_1D': 0.4,
                         'num_bins': 100,
                         'boundary_correction_order': 1,
                         'mult_bias_correction_order': 1,
                         'smooth_scale_2D': 0.4,
                         'max_corr_2D': 0.99,
                         'fine_bins_2D': 512,
                         'num_bins_2D': 40}



#STS.runs = xpsi.Runs.load_runs(ID='ST-S',
#                               #run_IDs=['run 1', 'run 2', 'mode separation'],
#                               #roots=['run1/run1_nlive1000_eff0.3_noCONST_noMM_noIS_tol-3',
#                               #       'run2/run2_nlive1000_eff0.3_noCONST_noMM_noIS_tol-3',
#                               #       'run3/run1_nlive1000_eff0.3_noCONST_MM_tol-1'],
#                               #base_dirs=['../../STU_v2/'] * 3,
#                               #use_nestcheck=[True,True,False],
#                               run_IDs=['run 1'],
#                               #roots=['run/run'],
#                               #roots=['run_test1/run'],
#                               roots=['output1/run1_nlive1000_eff0.3_noCONST_noMM_noIS_tol-1'],
#                               base_dirs=['.'],# * 3
#                               use_nestcheck=[True],#[False],
#                               kde_settings=getdist_kde_settings,
#                               likelihood=STS.likelihood,
#                               names=STS.names,
#                               bounds=STS.bounds,
#                               labels=STS.labels,
#                               implementation='multinest')#, #,
##                               overwrite_transformed=)


STU.runs = xpsi.Runs.load_runs(ID='ST-S',
                               #run_IDs=['run 1', 'run 2', 'mode separation'],
                               #roots=['run1/run1_nlive1000_eff0.3_noCONST_noMM_noIS_tol-3',
                               #       'run2/run2_nlive1000_eff0.3_noCONST_noMM_noIS_tol-3',
                               #       'run3/run1_nlive1000_eff0.3_noCONST_MM_tol-1'],
                               #base_dirs=['../../STU_v2/'] * 3,
                               #use_nestcheck=[True,True,False],
                               run_IDs=['run 1'],
                               #roots=['run/run'],
                               #roots=['job1v2/run'],
                               roots=['run_test1f/run'],
                               base_dirs=['.'],# * 3
                               use_nestcheck=[False],
                               kde_settings=getdist_kde_settings,
                               likelihood=STU.likelihood,
                               names=STU.names,
                               bounds=STU.bounds,
                               labels=STU.labels,
                               implementation='multinest')#, #,
#                               overwrite_transformed=)




pp = xpsi.PostProcessing.CornerPlotter([STU.runs])
_ = pp.plot(
     #params=['radius','compactness','mass'],
     #params=['radius','mass'],
     params=['radius', 'mass', 'distance', 'cos_inclination', 'p__phase_shift', 
'p__super_colatitude', 'p__super_radius', 'p__super_temperature', 's__phase_shift',
's__super_colatitude', 's__super_radius'],
     IDs=OrderedDict([('ST-S', ['run 1',]),]),
     prior_density=True,
     KL_divergence=True,
     ndraws=5e4,
     combine=True, combine_all=True, only_combined=False, overwrite_combined=True,
     param_plot_lims={},#'radius': (7.5,16.0)},#, 'mass': (1.0,3.0)
     bootstrap_estimators=False,
     bootstrap_density=False,
     n_simulate=200,
     crosshairs=False,
     write=True,#False,
     ext='.png',
     directory='figs/posteriors/',
     maxdots=3000,
     root_filename='STUmy_test1f',
     credible_interval_1d=True,
     annotate_credible_interval=True,
     compute_all_intervals=False,
     sixtyeight=True,
     x_label_rotation=45.0,
     num_plot_contours=3,
     subplot_size=4.0,
     legend_corner_coords=(0.675,0.8),
     legend_frameon=False,
     scale_attrs=OrderedDict([('legend_fontsize', 2.0),
                              ('lab_fontsize', 1.35),
                              ('axes_fontsize', 'lab_fontsize'),
                             ]
                            ),
     colormap='Reds',
     shaded=True,
     shade_root_index=-1,
     rasterized_shade=True,
     no_ylabel=True,
     no_ytick=True,
     lw=1.0,
     lw_1d=1.0,
     filled=False,
     normalize=True,
     veneer=True,
     tqdm_kwargs={'disable': False},
     lengthen=2.0,
     embolden=1.0,
     nx=500,
     scale_ymax=1.1)
     

     

