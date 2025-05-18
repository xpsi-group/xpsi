import numpy as np
import math

import xpsi

np.random.seed(xpsi._rank+10)

import time
import os
import sys

this_directory = os.path.dirname(os.path.abspath(__file__))
sys.path.append(this_directory)

# Data
from CustomInstrument import CustomInstrument
from CustomSignal import CustomSignal
from CustomPhotosphere import CustomPhotosphere
from CustomPrior import CustomPrior

data_path = this_directory+"/../Data/xpsi_good_realisation.dat"
data_loaded = np.loadtxt(data_path, dtype=np.double)

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


data_spectrum = np.sum(data.counts,axis=1)


registered_true_background = np.array([12774.19354839, 11975.80645161, 11250.        , 10588.23529412,
        9983.19327731,  9428.57142857,  8918.91891892,  8449.50213371,
        8016.19433198,  7615.38461538,  7243.90243902,  6898.95470383,
        6578.0730897 ,  6279.06976744,  6000.        ,  5739.13043478,
        5494.91211841,  5265.95744681,  5051.02040816,  4848.97959184,
        4658.82352941,  4479.63800905,  4310.59506531,  4150.94339623,
        4000.        ,  3857.14285714,  3721.80451128,  3593.46642468,
        3471.65400351,  3355.93220339,  3245.90163934,  3141.19513485,
        3041.47465438,  2946.42857143,  2855.76923077,  2769.23076923,
        2686.56716418,  2607.55048288,  2531.96930946,  2459.62732919,
        2390.34205231,  2323.94366197,  2260.2739726 ,  2199.18548686,
        2140.54054054,  2084.21052632,  2030.07518797,  1978.02197802,
        1927.94547225,  1879.74683544,  1833.33333333,  1788.61788618,
        1745.51866001,  1703.95869191,  1663.86554622,  1625.17099863,
        1587.81074579,  1551.72413793,  1516.85393258,  1483.14606742,
        1450.54945055,  1419.01576684,  1388.49929874,  1358.95676047,
        1330.34714446,  1302.63157895,  1275.77319588,  1249.73700821,
        1224.48979592,  1200.        ,  1176.23762376,  1153.17414094,
        1130.78241005,  1109.03659447,  1087.91208791,  1067.38544474,
        1047.43431494,  1028.03738318,  1009.17431193,   990.82568807,
         972.97297297,   955.5984556 ,   938.6852086 ,   922.21704704,
         906.1784897 ,   890.55472264,   875.33156499,   860.49543677,
         846.03332859,   831.93277311,   818.18181818,   804.76900149,
         791.68332667,   778.91424076,   766.4516129 ,   754.28571429,
         742.4071991 ,   730.80708661,   719.47674419,   708.4078712 ,
         697.59248385,   687.02290076,   676.69172932,   666.59185277,
         656.71641791,   647.05882353,   637.61270932,   628.37194541,
         619.33062246,   610.48304214,   601.82370821,   593.34731795,
         585.04875406,   576.92307692,   568.96551724,   561.17146906,
         553.53648309,   546.05626034,   538.72664611,   531.54362416,
         524.50331126,   517.6019519 ,   510.83591331,   504.20168067,
         497.69585253,   491.31513648,   485.05634493,   478.9163912 ,
         472.89228565,   466.98113208,   461.18012422,   455.48654244,
         449.89775051,   444.41119258,   439.02439024,   433.73493976,
         428.54050934,   423.43883661,   418.42772612,   413.50504699,
         408.66873065,   403.91676867,   399.24721065,   394.65816225,
         390.14778325,   385.71428571,   381.3559322 ,   377.07103409,
         372.85794991,   368.7150838 ,   364.64088398,   360.6338413 ,
         356.69248784,   352.81539558,   349.00117509,   345.24847428,
         341.55597723,   337.922403  ,   334.34650456,   330.82706767,
         327.36290989,   323.95287958,   320.59585492,   317.29074302,
         314.03647898,   310.83202512,   307.67637004,   304.56852792,
         301.50753769,   298.49246231,   295.52238806,   292.59642382,
         289.71370043,   286.87337004,   284.07460545,   281.31659957,
         278.5985648 ,   275.91973244,   273.27935223,   270.67669173,
         268.11103588,   265.58168649,   263.08796173,   260.62919574,
         258.2047381 ,   255.81395349,   253.4562212 ,   251.13093477,
         248.83750157,   246.57534247,   244.3438914 ,   242.14259508,
         239.97091262,   237.82831518,   235.71428571,   233.62831858,
         231.5699193 ,   229.53860422,   227.53390025,   225.5553446 ,
         223.60248447,   221.67487685,   219.7720882 ,   217.89369429,
         216.03927987,   214.20843851,   212.40077237,   210.61589193,
         208.85341584,   207.11297071,   205.39419087,   203.69671822,
         202.02020202,   200.36429872,   198.7286718 ,   197.11299154,
         195.51693493,   193.94018545,   192.38243296,   190.84337349,
         189.32270916,   187.82014798,   186.33540373,   184.86819582,
         183.41824919,   181.98529412,   180.56906615,   179.16930595,
         177.78575918,   176.41817642,   175.066313  ,   173.72992893,
         172.40878878,   171.1026616 ,   169.81132075,   168.53454391,
         167.27211287,   166.02381352,   164.78943572,   163.56877323,
         162.36162362,   161.16778815,   159.98707175,   158.81928291,
         157.66423358,   156.52173913,   155.39161827,   154.27369296,
         153.16778835,   152.07373272,   150.9913574 ,   149.92049671,
         148.8609879 ,   147.81267108,   146.77538918,   145.74898785,
         144.73331547,   143.728223  ,   142.73356401,   141.74919461,
         140.77497334,   139.81076119,   138.85642153,   137.91182002,
         136.97682463,   136.05130554,   135.13513514,   134.22818792,
         133.33034051,   132.44147157,   131.56146179,   130.69019384,
         129.82755229,   128.97342366,   128.12769629,   127.29026037,
         126.46100786,   125.63983248,   124.82662968,   124.02129659,
         123.22373198,   122.43383626,   121.65151143,   120.87666104,
         120.10919017,   119.34900542,   118.59601485,   117.85012797,
         117.11125569,   116.37931034,   115.65420561])

#print(len(registered_true_background))
#exit()

registered_true_background_rate = registered_true_background/1000.0/32.0 #Divide here with the exposure time!!

support = np.ones((len(registered_true_background_rate), 2), dtype=np.double)
support[:,0] = 0.099*registered_true_background_rate
support[:,1] = 10.91*registered_true_background_rate

#for i in range(support.shape[0]):
#    if support[i,1] == 0.0:
#        for j in range(1, support.shape[0]):
#            if i+j < support.shape[0] and support[i+j,1] > 0.0:
#                support[i,1] = support[i+j,1]
#                break
#            elif i-j >= 0 and support[i-j,1] > 0.0:
#                support[i,1] = support[i-j,1]
#                break
# # Signal
signal = CustomSignal(data = data,
                      instrument = Instrument,
                      support = support,
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
                           values=dict(frequency = 314.0),
                           star_shape="AGM_14")
#Default star shape is an oblate spheroid from AlGendy & Morsink (2014, AGM_14)
#But also a spherical star can be used with:
#spacetime.star_shape = "sphere"

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
                                prefix='p') # For "primary"


# # Photosphere
photosphere = CustomPhotosphere(hot = hot_spot, elsewhere = None,
                                values=dict(mode_frequency = spacetime['frequency']))


# # Star
star = xpsi.Star(spacetime = spacetime, photospheres = photosphere)

# # Prior
prior = CustomPrior()


# # Likelihood

likelihood = xpsi.Likelihood(star = star, signals = signal,
                             num_energies = 64,
                             threads = 1,
                             externally_updated = True,
                             prior = prior)

# Crucial step, if the likelihood check fails, then something went terrible wrong :)
p=[1.4,12,1.,math.cos(60*np.pi/180),0.0,70*np.pi/180, 0.75,6.7]

likelihood.check(None, [-3.1603740790e+04], 1.0e-5, physical_points=[p])

if __name__ == '__main__':

    start = time.time()
    wrapped_params = [0] * len(likelihood)
    wrapped_params[likelihood.index('p__phase_shift')] = 1

    #The original (more accurate) run settings shown as commented.
    runtime_params = {'resume': False,
                      'importance_nested_sampling': False,
                      'multimodal': False,
                      'n_clustering_params': None,
                      'outputfiles_basename': '../Outputs/ST_live_1000_eff_0.3_seed0_v2',
                      'n_iter_before_update': 50, #100,
                      'n_live_points': 50, #1000,
                      'sampling_efficiency': 0.3,
                      'const_efficiency_mode': False,
                      'wrapped_params': wrapped_params,
                      'evidence_tolerance': 0.1,
                      'max_iter': 100, #-1,
                      'seed' : 0, # Fixing the seed
                      'LHS_seed': 42, # Fixing the LHS seed for hypercube fraction estimation
                      'verbose': True}


    xpsi.Sample.run_multinest(likelihood, prior, **runtime_params)

    print('Sampling took', (time.time()-start)/60, 'minutes')
