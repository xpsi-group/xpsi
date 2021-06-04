

#%matplotlib inline

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


freq = 600.0 #800.0 #600.0 #300.0

#spacetime = xpsi.Spacetime.fixed_spin(freq)
#for p in spacetime:
#    print(p)

#Or alternatively :

#xpsi.Spacetime  #? # uncomment to query

#exit()

bounds = dict(distance = (0.1, 1.0),                     # (Earth) distance
                mass = (1.0, 3.0),                       # mass
                radius = (3.0 * gravradius(1.0), 16.0),  # equatorial radius
                cos_inclination = (0.0, 1.0))            # (Earth) inclination to rotation axis

spacetime = xpsi.Spacetime(bounds=bounds, values=dict(frequency=freq))


#print(xpsi.HotRegion) #? # uncomment to query
#print(xpsi.HotRegion.required_names)

#exit()

ceding = False

bounds = dict(super_colatitude = (None, None),
              super_radius = (None, None),
              phase_shift = (0.0, 0.1),
              super_temperature = (None, None))

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
                            #symmetry=False, 
                            omit=False,
                            cede=ceding, 
                            concentric=False,
                            sqrt_num_cells=32,
                            min_sqrt_num_cells=10,
                            max_sqrt_num_cells=64,
                            #num_leaves=100,
                            num_leaves=121,
                            num_rays=200,
                            do_fast=False,
                            prefix='p')

#print(primary.phases_in_cycles)
#print(primary.num_cells)
#print(primary.print_settings)
#exit()

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


print("For a secondary spot:")
bounds['super_temperature'] = None # declare fixed/derived variable

secondary = xpsi.HotRegion(bounds=bounds, # can otherwise use same bounds
                              values={'super_temperature': derive()}, # create a callable value
                              symmetry=True,
                              omit=False,
                              cede=ceding, 
                              concentric=False,
                              sqrt_num_cells=32,
                              min_sqrt_num_cells=10,
                              max_sqrt_num_cells=100,
                              num_leaves=100,
                              num_rays=200,
                              do_fast=False,
                              is_antiphased=True,
                              prefix='s') # unique prefix needed because >1 instance

from xpsi import HotRegions

hot = HotRegions((primary, secondary))


print(hot.objects[0]) # 'p'
print(hot.objects[1]) # 's'

h = hot.objects[0]
print(h.names)
print(h.prefix)

hot['p__super_temperature'] = 6.0 # equivalent to ``primary['super_temperature'] = 6.0``
print(secondary['super_temperature'])


class CustomPhotosphere(xpsi.Photosphere):
    """ Implement method for imaging."""

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

photosphere = CustomPhotosphere(hot = hot, elsewhere = None,
                                values=dict(mode_frequency = spacetime['frequency']))






print(photosphere['mode_frequency'] == spacetime['frequency'])

star = xpsi.Star(spacetime = spacetime, photospheres = photosphere)
print("Printing star:")
print(star)

print("Parameters of the star:")
print(star.params)
#Parameters of the star:
#[Gravitational mass [solar masses], 
#Coordinate equatorial radius [km], 
#Earth distance [kpc], Cosine of Earth inclination to rotation axis, 
#The phase of the hot region, a periodic parameter [cycles], 
#The colatitude of the centre of the superseding region [radians], 
#The angular radius of the (circular) superseding region [radians], 
#log10(superseding region effective temperature [K]) = 6.000e+00, 
#The colatitude of the centre of the ceding region [radians], 
#The angular radius of the (circular) ceding region [radians], 
#The azimuth of the centre of the ceding region relative to the centre of the superseding region [radians], 
#log10(ceding region effective temperature [K]), 
#The phase of the hot region, a periodic parameter [cycles], 
#The colatitude of the centre of the superseding region [radians], 
#The angular radius of the (circular) superseding region [radians], 
#The colatitude of the centre of the ceding region [radians], 
#The angular radius of the (circular) ceding region [radians], 
#The azimuth of the centre of the ceding region relative to the centre of the superseding region [radians], 
#log10(ceding region effective temperature [K])]


#For fixed values:
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

star(p)


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


#Another way to set param values:

#compactness = 0
star['mass'] = 1.4 #2.7088795 #1.4 #0.112*20.0 #1.4
star['radius'] = 12.0 #20.0 #12.0
star['distance'] = 0.2
incl_deg = 40.0 #90.0 #40.0
star['cos_inclination'] = math.cos(math.pi*incl_deg/180.0)#math.cos(2.0)
theta_deg = 60.0
star['p__super_colatitude'] = math.pi*theta_deg/180.0 #0.0 #2.0
rho_deg = 10.0 #10.0 #0.001 #10.0
star['p__super_radius'] = math.pi*rho_deg/180.0
#print("rho[deg]=",math.pi*rho_deg/180.0)
tplanck = 1.0219978 #1.0 #in keV #1 keV -> log10(T[K]) = 7.06 (out of bounds originally)
#print(np.log10(tplanck*11604525.0061657))
star['p__super_temperature'] = np.log10(tplanck*11604525.0061657)

#star['s__phase_shift'] = 0.025
#star['s__super_colatitude'] = 2.142
#star['s__super_radius'] = 0.2


if ceding:
	star['p__cede_colatitude'] = math.pi*theta_deg/180.0 
	star['p__cede_azimuth'] = 0.0 #0.0
	star['p__cede_radius'] = math.pi*rho_deg/180.0+0.02
	star['p__cede_temperature'] = np.log10(tplanck*11604525.0061657)
	star['s__super_radius'] = math.pi*rho_deg/180.0
	star['s__super_colatitude'] = math.pi-math.pi*theta_deg/180.0
	star['s__cede_colatitude'] = math.pi-math.pi*theta_deg/180.0 
	star['s__cede_azimuth'] = 0.1 #math.pi/2.0
	star['s__cede_radius'] = math.pi*rho_deg/180.0+0.02
	star['s__cede_temperature'] = np.log10(tplanck*11604525.0061657)


print("Parameters of the star:")
print(star.params)
#exit()


#photosphere._hot_atmosphere = (logT, logg, mu, logE, buf) ???? -Check the example script and CustomPhotosphere in tutorial


#Now some pulse profiles:

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

def plot_pulse():
    """ Plot hot region signals before telescope operation. """
    fig = plt.figure(figsize=(7,7))
    ax = fig.add_subplot(111)

    ax.set_ylabel('Signal [arbitrary normalisation]')
    ax.set_xlabel('Phase [cycles]')

    temp = np.sum(photosphere.signal[0][0], axis=0)
    #temp = photosphere.signal[0][0][nene-1,:] 
    ax.plot(hot.phases_in_cycles[0], temp/np.max(temp), 'o-', color='k', lw=0.5, markersize=2)
    temp = np.sum(photosphere.signal[1][0], axis=0)
    ax.plot(hot.phases_in_cycles[1], temp/np.max(temp), 'o-', color='r', lw=0.5, markersize=2)

    veneer((0.05,0.2), (0.05,0.2), ax)
    fig.savefig("figs/pulse_profileX.pdf")

nene = 281 #128

def save_pulse(PulsName): #To be continued ...
    """Save the pulse profile in a file. """
    #print("F(E) = ",photosphere.signal[0][0][:,0], len(photosphere.signal[0][0][:,0])) #np.shape
    #print("F(T) = ", photosphere.signal[0][0][0,:], len(photosphere.signal[0][0][0,:]))
    #print("F(E) = ",photosphere.signal[0][0][:,50], len(photosphere.signal[0][0][:,50])) #np.shape
    #print("F(E) = ", photosphere.signal[1][0][:,50], len(photosphere.signal[1][0][:,50]))

    #print("F(E) = ",photosphere.signal[0][1][:,50], len(photosphere.signal[0][1][:,50])) #np.shape
    #print("F(E) = ", photosphere.signal[1][1][:,50], len(photosphere.signal[1][1][:,50]))
    #exit()

    #print(np.shape(photosphere.signal_stokes))
    #print("F(E) = ",photosphere.signal_stokes[0][1][:,0], len(photosphere.signal_stokes[0][0][:,0])) #np.shape
    #print("F(T) = ", photosphere.signal_stokes[0][1][0,:], len(photosphere.signal_stokes[0][0][0,:]))
    #print("F(E) = ",photosphere.signal_stokes[0][2][:,0], len(photosphere.signal_stokes[0][0][:,0])) #np.shape
    #print("F(T) = ", photosphere.signal_stokes[1][0][0,:], len(photosphere.signal_stokes[0][0][0,:]))
    #print("F(E) = ",photosphere.signal_stokes[0][5][:,0], len(photosphere.signal_stokes[0][0][:,0])) #np.shape
    #print("F(T) = ", photosphere.signal_stokes[1][3][0,:], len(photosphere.signal_stokes[0][0][0,:]))

    #print(photosphere.signal)

    #pulse1 = np.sum(photosphere.signal[0][0], axis=0)
    #pulse1 = photosphere.signal[0][0][nene-1,:] #pulse in one energy bin
    #pulse1 = photosphere.signal_stokes[0][0][nene-1,:]
    #pulseQ = photosphere.signal_stokes[0][1][nene-1,:] 
    #pulseU = photosphere.signal_stokes[0][2][nene-1,:] 
    pulse1 = photosphere.signal_stokes[0][0][118,:]
    pulseQ = photosphere.signal_stokes[0][1][118,:] 
    pulseU = photosphere.signal_stokes[0][2][118,:] 




    #pulse2 = np.sum(photosphere.signal[1][0], axis=0)
    phase1 = hot.phases_in_cycles[0]
    #phase2 = hot.phases_in_cycles[1]
    #pulse_tot = pulse1+pulse2
    outF = open(PulsName + '_F.bin','w')
    outf = open(PulsName + '_p.bin','w')
    pulse1.tofile(outF,format="%e") 
    phase1.tofile(outf,format="%e")
    outQ = open(PulsName + '_Q.bin','w')
    pulseQ.tofile(outQ,format="%e") 
    outU = open(PulsName + '_U.bin','w')
    pulseU.tofile(outU,format="%e") 


from numpy import logspace
from numpy import log
evere=0.5109989e6 # electron volts in elecron rest energy 
x_l, x_u = -3.7 , .3 # lower and upper bounds of the log_10 energy span
IntEnergy = logspace(x_l,x_u,nene), log(1e1)*(x_u-x_l)/(nene-1.) # sample poi
x, x_weight = IntEnergy
energies = (x*evere)/1e3

#energies = np.logspace(-1.0, np.log10(3.0), 128, base=10.0)
#energies = np.logspace(-1.0, np.log10(4.94), nene, base=10.0)
#energies = np.logspace(-1.0, np.log10(4.86), nene, base=10.0) #this seems to match better

#print("energies (keV) =",energies)

star.update() #Calculating the space-time integrals etc. 

#NOTE: Atmosphere evaluation is only applied during the following integration:
#print("Now the actual computation!:")


import time

primary.image_order_limit = 1


start = time.time()
photosphere.integrate(energies, threads=1, stokes=True) # the number of OpenMP threads to use
end = time.time()
print("Time spent in integration:",end - start)

#exit()

#print("F(E) = ",photosphere.signal_stokes[0][1][:,0], len(photosphere.signal[0][0][:,0])) #np.shape

#_ = plot_pulse()
#print("Light curve finished,")
#print(photosphere.signal[0][0])


#plot_pulse()
#save_pulse("pulses/pulse_f800r20")
#save_pulse("pulses/pulse_io2")
#save_pulse("pulses/pulse7i_rho10f600m27_Tc_io1")
#save_pulse("pulses/pulse7i_rho10f600_Tc_IQUi2")
save_pulse("pulses/pulse_lisa0")





