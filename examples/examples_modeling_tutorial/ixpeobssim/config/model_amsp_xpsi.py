#!/usr/bin/env python
#
# Copyright (C) 2018, the ixpeobssim team.
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

"""
This is an example, illustrating how to fetch an accreting millisecond pulsar
model from X-PSI to ixpeobssim.
"""

from __future__ import print_function, division

import numpy
import matplotlib
import os

from ixpeobssim.config import file_path_to_model_name
from ixpeobssim.utils.matplotlib_ import plt, setup_gca
from ixpeobssim.srcmodel.roi import xPeriodicPointSource, xROIModel
from ixpeobssim.srcmodel.ephemeris import xEphemeris
from ixpeobssim.srcmodel.spectrum import power_law
from ixpeobssim.srcmodel.polarization import constant
from ixpeobssim.utils.fmtaxis import fmtaxis
import numpy as np
from numpy import logspace, zeros, fromfile
from numpy import pi, exp, log, sqrt, sin, cos, arccos, arctan2
from scipy.interpolate import interp2d
from astropy.io import fits
from pylab import *
from ixpeobssim import IXPEOBSSIM_SRCMODEL
from ixpeobssim.srcmodel.ephemeris import xOrbitalEphemeris

import xpsi

from xpsi.global_imports import _c, _G, _dpr, gravradius, _csq, _km, _2pi

import sys
sys.path.append('../../')
sys.path.append('../')
from TestRun_PolNum_split_1spot import get_photosphere_stokes_1spot
#from TestRun_Pol import get_photosphere_stokes_1spot

phase, energies, photosphere_I, photosphere_Q, photosphere_U = get_photosphere_stokes_1spot()

print(np.shape(photosphere_I))

NEnergy, NPhase = np.shape(photosphere_I)[0], np.shape(photosphere_I)[1]

Imod=numpy.zeros((NPhase,NEnergy))
Qmod=numpy.zeros((NPhase,NEnergy))
Umod=numpy.zeros((NPhase,NEnergy))
PAobs=numpy.zeros((NPhase,NEnergy))
PDobs=numpy.zeros((NPhase,NEnergy))

for j in range(NPhase):
	for i in range(NEnergy):
		Imod[j,i],Qmod[j,i],Umod[j,i]=photosphere_I[i,j], photosphere_Q[i,j], photosphere_U[i,j]

chi = 0.0 #pulsar rotation axis position angle
chi_rad = chi*pi/180.0

distance_kpc = 3.5
distance_m = 3.08567758128e19*distance_kpc

#Note that if want to obtain the specific flux in units of photons/cm/s/keV instead, the output of photosphere.signal needs to be divided by distance squared, where distance is measured in meters

for ii in range(0,NEnergy):

	Q_obs = cos(2*chi_rad)*Qmod[:,ii]-sin(2*chi_rad)*Umod[:,ii]
	U_obs = sin(2*chi_rad)*Qmod[:,ii]+cos(2*chi_rad)*Umod[:,ii]

	PAobs[:,ii] = arctan2(-U_obs,-Q_obs)*90/pi+90
	for jj in range(0,NPhase):
		if(Imod[jj,ii] > 1e-30):
			PDobs[jj,ii] = sqrt(U_obs[jj]**2+Q_obs[jj]**2)/Imod[jj,ii]
		else:
			PDobs[jj,ii] = 0.0
									
	#Unit change to cm^-2 s^-1 keV^-1 from photosphere.signal units 
	Imod[:,ii] = Imod[:,ii]/distance_m**2
	
brightn = 100 # Target's brightness in mCrab
highest_I = np.max(Imod)
Icrab = highest_I/0.0019986928

#Using now a rough magic scaling relation here.
#In a more realistic case, no scaling should be done here.
print("Initial brightness in mCrab:",Icrab)
print("Correction factor:",brightn*0.0019986928/highest_I)
Imod = Imod/(highest_I)*brightn*0.0019986928

highest_I = np.max(Imod)
Icrab = highest_I/0.0019986928
print("Brightness in mCrab:",Icrab)

kev2erg = 1.6021766339999e-9
Isum = 0
for ie in range(len(energies)):
        if (2.0 < energies[ie] < 10.0):
                Iplus = energies[ie]*np.sum(Imod[0:NPhase,ie])/NPhase + energies[ie+1]*np.sum(Imod[0:NPhase,ie+1])/NPhase
                Isum = Isum + Iplus*(energies[ie+1]-energies[ie])/2.0
print("Self-integrated flux:",Isum*kev2erg)

Imodspec = np.sum(Imod,axis=0)/NPhase
from scipy import integrate
Imod_intg = integrate.simpson(energies[236:363]*Imodspec[236:363], energies[236:363])
print("Scipy-integrated flux:",Imod_intg*kev2erg)

print("Brightness in mCrab (more accurate):", 1000.0*(Isum*kev2erg)/2.4e-8)
#Crab: 2- 10 keV: 2.4e-8 erg/cm2/s
#flux ~130 mCrab (2.8×10−9 erg/cm2/s) in 2-8 keV  (https://ui.adsabs.harvard.edu/abs/2023IAUS..363..329D/abstract).
#https://en.wikipedia.org/wiki/Crab_(unit)

phi, energy_keV, St_I, PD_tot, PA_tot = phase, energies, Imod, PDobs, PAobs

__model__ = file_path_to_model_name(__file__)
#ra, dec = 45., 45.
ra, dec = 272.11475, -36.97897
#http://simbad.u-strasbg.fr/simbad/sim-id?Ident=SAX+J1808.4-3658&NbIdent=1&Radius=2&Radius.unit=arcmin&submit=submit+id

#E and phase will be tuple when called from ixpeobssim
#'''
def spec(E, phase):
	"""Definition of the energy spectrum (cm^-2 s^-1 keV^-1).
	"""
	Stokes_I_interp = interp2d(energy_keV, phi, St_I, kind='linear')

	try:
		return Stokes_I_interp(E,phase)
	except: #in case of ixpeobssim
		len_E = len(E[:,0])
		len_phase = len(E[0,:])
		flux = np.zeros((len_E,len_phase))
		for ii in range(len_E):
			for jj in range(len_phase):
				flux[ii,jj] = Stokes_I_interp(E[ii,0],phase[0,jj])
		return flux

def pol_deg(E, phase, ra=None, dec=None):
	"""Polarization degree as a function of the dynamical variables.

	Since we're dealing with a point source the sky direction (ra, dec) is 
	irrelevant and, as they are not used, defaulting the corresponding arguments
	to None allows to call the function passing the energy and phase only.
	"""
	PD_interp = interp2d(energy_keV, phi, PD_tot, kind='linear')
	is_array=True
	try:
		#the following lengths should be same:
		len_E=len(E)
		len_phase=len(phase)
		if(len_E != len_phase):
			print("len_E has to be same as len_phase")
	except:#if calling with a single float
		is_array = False
	if is_array: 
		PD_all = np.zeros((len_E))
		for ii in range(len_E):
			PD_all[ii] = PD_interp(E[ii],phase[ii])
	else:
		PD_all = PD_interp(E, phase)
	return PD_all


def pol_ang(E, phase, ra=None, dec=None):
	"""Definition of the polarization angle (in radians).
	"""
	PA_interp = interp2d(energy_keV, phi, PA_tot, kind='linear')
	is_array=True
	try:
		#the following lengths should be same:
		len_E=len(E)
		len_phase=len(phase)
		if(len_E != len_phase):
			print("len_E has to be same as len_phase")
	except:#if testing with calling with only a single float
		is_array = False
	if is_array: 
		PA_all = np.zeros((len_E))
		for ii in range(len_E):
			PA_all[ii] = PA_interp(E[ii],phase[ii])
	else:
		PA_all = PA_interp(E, phase)
	return np.deg2rad(PA_all)


start_date = '2023-10-06'

nu0 = 400.9752075
ephemeris = xEphemeris(0., nu0)
src = xPeriodicPointSource('AMSP', ra, dec, spec, pol_deg, pol_ang,
                           ephemeris)

ROI_MODEL = xROIModel(ra, dec, src)


def display(emin=1., emax=12.):
    """Display the source model.
    """
    energy = numpy.linspace(emin, emax, 100)
    phase = numpy.linspace(0., 1., 100)    
    # Pulse profile: polarization degree.
    plt.figure('%s polarization degree' % __model__)
    for E in [2., 5., 8.]:
        plt.plot(phase, pol_deg(E, phase), label='Energy = %.2f keV' % E)
    setup_gca(ymin=0., ymax=0.2, legend=True, **fmtaxis.pp_pol_deg)

    # Pulse profile: polarization angle.
    plt.figure('%s polarization angle' % __model__)
    for E in [2., 5., 8.]:
        plt.plot(phase, np.rad2deg(pol_ang(E, phase)), label='Energy = %.2f keV' % E)
    setup_gca(ymin=-90., ymax=180.0, legend=True, **fmtaxis.pp_pol_ang)
    
    # Pulse profile: Flux.
    plt.figure('%s pulse profile' % __model__)
    for E in [2., 5., 8.]:
        #print("E,spec",E,spec(E, phase))
        plt.plot(phase, spec(E, phase), label='Energy = %.2f keV' % E)
    setup_gca(ymin=0., legend=True, **fmtaxis.spec)

    # Energy spectrum at different phase values.
    plt.figure('%s spectrum' % __model__)
    for p in [0.0,0.2,0.4,0.6,0.8,1.0]:
        plt.plot(energy, spec(energy, p), label='Phase = %.2f' % p)
    setup_gca(xmin=emin, xmax=emax, logx=True, logy=True, legend=True,
              grids=True, **fmtaxis.spec)



if __name__ == '__main__':
    print("PA at 5 kev, t= 0.5ph",pol_ang(5.0,0.5))
    from ixpeobssim.config import bootstrap_display
    bootstrap_display()
