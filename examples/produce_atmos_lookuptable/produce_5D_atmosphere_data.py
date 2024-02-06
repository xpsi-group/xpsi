#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 16 10:11:26 2022

@author: bas

Modified by T.Salmi on Dec 7 2023 to save also the Stokes Q parameter.

This script can be used to convert the atmosphere data files from
https://github.com/AnnaBobrikova/ComptonSlabTables
to a format that can be used in X-PSI.
"""

from numpy import linspace, logspace, empty, zeros, ones, array, fromfile
from numpy import pi, exp, log, sqrt, sin, cos, arccos, arctan2
from numpy import absolute, sign, floor, ceil, argmin
from numpy.polynomial.laguerre import laggauss
from numpy.polynomial.legendre import leggauss
from scipy.interpolate import interp1d
from scipy.interpolate import CubicSpline 
from scipy.interpolate import interp2d
from scipy.special import kn
from bisect import bisect
import numpy as np
from astropy.io import fits
from numpy import ones, zeros
from scipy import interpolate

t__e = np.arange(40.0, 202.0, 4.0) #te = Te(keV)*1000/511keV
t__bb = np.arange(0.001, 0.0031, 0.0002) #tbb = Tbb(keV)/511keV
tau__t = np.arange(0.5, 3.55, 0.1) 

NEnergy_i = 150
NZenith_i = 9

I_mighty = ones((len(t__e), len(t__bb), len(tau__t), NEnergy_i, NZenith_i)) #5D Intensity tensor[Te,Tbb,tau_t, Energy, Mu]
Q_mighty = ones((len(t__e), len(t__bb), len(tau__t), NEnergy_i, NZenith_i)) #5D Q-parameter tensor[Te,Tbb,tau_t, Energy, Mu]



def reading_table(): #routine that just reads fits tables to I_mighty and Q_mighty
    path='../examples_modeling_tutorial/model_data/ComptonSlabTables/'
    global I_mighty
    global Q_mighty
    p=0 #as we will not go over length_Te but just over values within t__e array, we'll need a counter for index corr. to Te
    for i in t__e:
        hdu = fits.open(path+'CompSlab_%s.fits' %(i), mode='update')
        print('With CompSlab_%s.fits still works' %(i)) #just to see sth on the screen while the files are being opened
        for ii in range(0,len(t__bb)): #Tbb
            for iii in range(0,len(tau__t)): #that's tau_t
                science_data = hdu[(iii)*len(t__bb)+ii+1].data #this messy index leads to the correct "list" within FITS file
                data = np.transpose(science_data.field(1)) #we need to transpose them, because tables were generated by julia and read in python, and there languages have different "majoring" in rows or columns, afaik
                data2 = np.transpose(science_data.field(0))
                for kk in range(0, NEnergy_i): #energy
                    for kkk in range(0,NZenith_i): #zenith angles
                        I_mighty[p, ii, iii, kk, kkk] = data[kk, kkk]#no +1 anymore
                        Q_mighty[p, ii, iii, kk, kkk] = data2[kk, kkk]
        p +=1
        
        
reading_table()


x_l, x_u = -3.7, .3 # lower and upper bounds of the log_10 energy span
NEnergy = 150 # 50# 101 # number of energy points (x)
IntEnergy = np.logspace(x_l,x_u,NEnergy), np.log(1e1)*(x_u-x_l)/(NEnergy-1.) # sample points and weights for integrations over the spectrum computing sorce function
Energy,x_weight=IntEnergy

from numpy.polynomial.legendre import leggauss

def init_mu(n = 3):
        NMu = n # number of propagation zenith angle cosines (\mu) [0,1]
        NZenith = 2*NMu # number of propagation zenith angles (z) [0,pi]
        mu = np.empty(NZenith)
        m2,mw = leggauss(NMu)
        mu[:NMu] = (m2 - 1.)/2
        mu[NMu:NZenith] = (m2 + 1.)/2
        return mu[NMu:NZenith]

mu = init_mu(9)


# I(E < mu < tau < tbb < te)

INDEX = 0
intensities_size = 1
shapes = list(I_mighty.shape)
for shape in shapes:
    intensities_size *= shape

NSX = np.empty((intensities_size,6))

print(len(t__e))
print(len(t__bb))
print(len(tau__t))
print(len(mu))
print(len(Energy))

for m in range(len(t__e)):
    for l in range(len(t__bb)):
        for k in range(len(tau__t)):
            for j in range(len(mu)):
                for i in range(len(Energy)):
                    NSX[INDEX,0] = Energy[i]
                    NSX[INDEX,1] = mu[j]
                    NSX[INDEX,2] = tau__t[k]
                    NSX[INDEX,3] = t__bb[l]
                    NSX[INDEX,4] = t__e[m]
                    NSX[INDEX,5] = I_mighty[m,l,k,i,j]
                    INDEX += 1
path_save='../examples_modeling_tutorial/model_data/'
np.savez_compressed(path_save+'Bobrikova_compton_slab_I.npz', NSX=NSX, size=shapes)


# Q(E < mu < tau < tbb < te)

INDEX = 0
intensities_size = 1
shapes = list(I_mighty.shape)
for shape in shapes:
    intensities_size *= shape

NSX = np.empty((intensities_size,6))

print(len(t__e))
print(len(t__bb))
print(len(tau__t))
print(len(mu))
print(len(Energy))

for m in range(len(t__e)):
    for l in range(len(t__bb)):
        for k in range(len(tau__t)):
            for j in range(len(mu)):
                for i in range(len(Energy)):
                    NSX[INDEX,0] = Energy[i]
                    NSX[INDEX,1] = mu[j]
                    NSX[INDEX,2] = tau__t[k]
                    NSX[INDEX,3] = t__bb[l]
                    NSX[INDEX,4] = t__e[m]
                    NSX[INDEX,5] = Q_mighty[m,l,k,i,j]
                    INDEX += 1
np.savez_compressed(path_save+'Bobrikova_compton_slab_Q.npz', NSX=NSX, size=shapes)

