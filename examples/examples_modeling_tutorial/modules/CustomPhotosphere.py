#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  6 16:52:43 2024

@author: bas
"""

import xpsi
from xpsi import Everywhere, Elsewhere, HotRegion, Parameter
import numpy as np


class CustomPhotosphere_NumA5(xpsi.Photosphere):
    """ A photosphere extension to preload the numerical 5D accretion atmosphere. """
    
    
    def __init__(self,
                 hot = None, elsewhere = None,
                 everywhere = None,
                 bounds = None, values = None,
                 stokes=False,
                 custom = None,
                 disk = None,
                 disk_blocking = True, #override, set to false iff disk but no blocking (unphysical but ok, just to test)
                 **kwargs):

        if everywhere is not None:
            if hot or elsewhere is not None:
                raise ValueError('Cannot use hot region nor elsewhere '
                                 'functionality if constructing the '
                                 'radiation field everywhere.')
            if not isinstance(everywhere, Everywhere):
                raise TypeError('Invalid type for everywhere object.')
        elif hot is None and elsewhere is None:
            pass # can call image-plane extensions
        else:
            if elsewhere is not None:
                if not isinstance(elsewhere, Elsewhere):
                    raise TypeError('Invalid type for an elsewhere object.')

                if hot is None:
                    raise ValueError('Hot region object(s) must be used in '
                                     'conjuction with an elsewhere object.')

            self._elsewhere_atmosphere = ()
                                              # including derived classes
            if hot is not None and hot is not isinstance(hot, HotRegion):
                if hasattr(hot, 'objects'):
                    for obj in getattr(hot, 'objects'):
                        if not isinstance(obj, HotRegion):
                            raise TypeError('Invalid object for the hot '
                                            'region(s).')
                else:
                    raise TypeError('Invalid object for the hot region(s).')

        self._hot = hot
        self._hot_atmosphere = ()
        self._hot_atmosphere_Q = ()
        self._elsewhere = elsewhere
        self._everywhere = everywhere
        self._stokes = stokes
        self._disk_blocking = disk_blocking # disk occultation

        if disk is not None:
            self._disk = disk
        else:
            self._disk = None

        if hot is not None:
            self._surface = self._hot
        else:
            self._surface = self._everywhere
            self.surface.objects = [self.surface]

        if bounds is None: bounds = {}
        if values is None: values = {}

        doc = """
        Coordinate frequency of the mode of radiative asymmetry in the
        photosphere that is assumed to generate the pulsed signal [Hz].
        """
        mode_frequency = Parameter('mode_frequency',
                                   strict_bounds = (0.0, 2000.0),
                                   bounds = bounds.get('mode_frequency', None),
                                   doc = doc,
                                   symbol = r'$f_{\rm mode}$',
                                   value = values.get('mode_frequency', None))

        custom = []
        if disk is not None:
            custom.append(disk)
        
        if stokes:
            doc = """
            Spin axis position angle measured from the north counterclock-
            wise to the projection of the rotation axis on the plane of the
            sky [in radians].
            """
            spin_axis_position_angle = Parameter('spin_axis_position_angle',
                                       strict_bounds = (-np.pi/2.0, np.pi/2.0),
                                       bounds = bounds.get('spin_axis_position_angle', None),
                                       doc = doc,
                                       symbol = r'$\chi_{0}$',
                                       value = values.get('spin_axis_position_angle', None))

            # print('everywhere:', everywhere)
            # print('hotregion:', hot)
            
            super(CustomPhotosphere_NumA5, self).__init__(mode_frequency=mode_frequency, spin_axis_position_angle=spin_axis_position_angle,
                                              hot=hot, elsewhere=elsewhere, everywhere=everywhere,
                                              bounds=bounds, values=values,
                                              stokes=stokes,
                                              custom=custom,
                                              **kwargs)
        else:
            super(CustomPhotosphere_NumA5, self).__init__(mode_frequency=mode_frequency,
                                              hot=hot, elsewhere=elsewhere, everywhere=everywhere,
                                              bounds=bounds, values=values,
                                              custom=custom,
                                              **kwargs)

    @property
    def disk(self):
        """ Get the instance of :class:`~.Disk.Disk`. """
        return self._disk

    @xpsi.Photosphere.hot_atmosphere.setter
    def hot_atmosphere(self, path):
        with np.load(path, allow_pickle=True) as data_dictionary:
            NSX = data_dictionary['NSX.npy']
            size_reorderme = data_dictionary['size.npy']

        size = [size_reorderme[3], size_reorderme[4], size_reorderme[2], size_reorderme[1], size_reorderme[0]]

        Energy = np.ascontiguousarray(NSX[0:size[0],0])
        cos_zenith = np.ascontiguousarray([NSX[i*size[0],1] for i in range(size[1])])
        tau = np.ascontiguousarray([NSX[i*size[0]*size[1],2] for i in range(size[2])])
        t_bb = np.ascontiguousarray([NSX[i*size[0]*size[1]*size[2],3] for i in range(size[3])])
        t_e = np.ascontiguousarray([NSX[i*size[0]*size[1]*size[2]*size[3],4] for i in range(size[4])])
        intensities = np.ascontiguousarray(NSX[:,5])

        self._hot_atmosphere = (t_e, t_bb, tau, cos_zenith, Energy, intensities)

    @xpsi.Photosphere.hot_atmosphere_Q.setter
    def hot_atmosphere_Q(self, path):
        with np.load(path, allow_pickle=True) as data_dictionary:
            NSX = data_dictionary['NSX.npy']
            size_reorderme = data_dictionary['size.npy']

        size = [size_reorderme[3], size_reorderme[4], size_reorderme[2], size_reorderme[1], size_reorderme[0]]

        Energy = np.ascontiguousarray(NSX[0:size[0],0])
        cos_zenith = np.ascontiguousarray([NSX[i*size[0],1] for i in range(size[1])])
        tau = np.ascontiguousarray([NSX[i*size[0]*size[1],2] for i in range(size[2])])
        t_bb = np.ascontiguousarray([NSX[i*size[0]*size[1]*size[2],3] for i in range(size[3])])
        t_e = np.ascontiguousarray([NSX[i*size[0]*size[1]*size[2]*size[3],4] for i in range(size[4])])
        intensities = np.ascontiguousarray(NSX[:,5])

        self._hot_atmosphere_Q = (t_e, t_bb, tau, cos_zenith, Energy, intensities)
        
    
    def integrate(self, energies, threads):
        """ Integrate over the photospheric radiation field.

        :param energies:
            A one-dimensional :class:`numpy.ndarray` of energies in keV.

        :param int threads:
            Number of ``OpenMP`` threads to spawn for signal integration.

        :param bool stokes:
            If activated, a full Stokes vector is computed and stored in signal, signalQ, and signalU.

        """
        if self._everywhere is not None:
            if self._stokes:
                raise NotImplementedError('Stokes option for everywhere not implmented yet.')      
            spectrum = self._everywhere.integrate(self._spacetime,
                                                   energies,
                                                   threads,
                                                   self._hot_atmosphere)
            if spectrum.ndim == 1:
                self._signal = ((spectrum.reshape(-1,1),),)
            else:
                self._signal = ((spectrum,),)
        else:
            if self._elsewhere is not None:
                spectrum = self._elsewhere.integrate(self._spacetime,
                                                     energies,
                                                     threads,
                                                     *self._elsewhere_atmosphere)

            if self._hot is not None:
                try:
                    else_atm_ext = self._elsewhere.atm_ext
                except:
                    else_atm_ext = None

                R_in = self.disk['R_in'] * 1000
                if self._stokes:
                    if self._disk_blocking:
                        self._signal, self._signalQ, self._signalU  = self._hot.integrate_stokes(self._spacetime,
                                                    energies,
                                                    threads,
                                                    self._hot_atmosphere,
                                                    self._hot_atmosphere_Q,
                                                    self._elsewhere_atmosphere,
                                                    else_atm_ext,
                                                    R_in=R_in)
                    elif not self._disk_blocking:
                        self._signal, self._signalQ, self._signalU  = self._hot.integrate_stokes(self._spacetime,
                                                        energies,
                                                        threads,
                                                        self._hot_atmosphere,
                                                        self._hot_atmosphere_Q,
                                                        self._elsewhere_atmosphere,
                                                        else_atm_ext)
                    if not isinstance(self._signal[0], tuple):
                        self._signal = (self._signal,)
                    if not isinstance(self._signalQ[0], tuple):
                        self._signalQ = (self._signalQ,)
                    if not isinstance(self._signalU[0], tuple):
                        self._signalU = (self._signalU,)
                    #Rotate the Stokes parameters based on position of the spin axis:
                    chi_rad = self["spin_axis_position_angle"]
                    tempQ = [list(x) for x in self._signalQ]
                    tempU = [list(x) for x in self._signalU]
                    for ih in range(0,len(self._signalQ)):
                        for ic in range(0,len(self._signalQ[ih][:])):
                            tempQ[ih][ic] = np.cos(2.0*chi_rad) * tempQ[ih][ic] - np.sin(2.0*chi_rad) * tempU[ih][ic]
                            tempU[ih][ic] = np.sin(2.0*chi_rad) * tempQ[ih][ic] + np.cos(2.0*chi_rad) * tempU[ih][ic]
                    self._signalQ = tuple(map(tuple, tempQ))
                    self._signalU = tuple(map(tuple, tempU))
                else:
                    if self._disk_blocking:
                        self._signal = self._hot.integrate(self._spacetime,
                                                    energies,
                                                    threads,
                                                    self._hot_atmosphere,
                                                    self._elsewhere_atmosphere,
                                                    else_atm_ext,
                                                    R_in=R_in)
                    elif not self._disk_blocking:
                         self._signal = self._hot.integrate(self._spacetime,
                                                    energies,
                                                    threads,
                                                    self._hot_atmosphere,
                                                    self._elsewhere_atmosphere,
                                                    else_atm_ext)
                    if not isinstance(self._signal[0], tuple):
                        self._signal = (self._signal,)

                # add time-invariant component to first time-dependent component
                if self._elsewhere is not None:
                    for i in range(self._signal[0][0].shape[1]):
                        self._signal[0][0][:,i] += spectrum    
    
            if self._disk is not None: 
                self.disk_spectrum = self._disk(energies)
                for i in range(self._signal[0][0].shape[1]):
                    self._signal[0][0][:,i] += self.disk_spectrum    