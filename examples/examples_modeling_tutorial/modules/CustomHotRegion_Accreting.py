#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  6 16:46:36 2024

@author: bas
"""

import xpsi
from xpsi import Parameter
import numpy as np
from xpsi.Photosphere import Photosphere

class CustomHotRegion_Accreting(xpsi.HotRegion):
    """Custom implementation of HotRegion. Accreting Atmosphere model by 
    Anna Bobrikova. The parameters are ordered I(E < mu < tau < tbb < te).
    OR LINEAR MODEL
    
    E is energy.
    mu is cos of zenith angle.
    tau is the optical depth of the comptom slab.
    tbb is the black body temperature.
    te is temperature of the electron gas.
    """

    required_names = ['super_colatitude',
                      'super_radius',
                      'phase_shift']
    optional_names = ['omit_colatitude',
                      'omit_radius',
                      'omit_azimuth',
                      'cede_colatitude',
                      'cede_radius',
                      'cede_azimuth',
                      'cede_tbb',
                      'cede_te',
                      'cede_tau',
                      'super_tbb',
                      'super_te',
                      'super_tau',
                      'super_h']
    
    def __init__(self,
            bounds,
            values,
            symmetry = True,
            interpolator = 'split',
            omit = False,
            cede = False,
            concentric = False,
            sqrt_num_cells = 32,
            min_sqrt_num_cells = 10,
            max_sqrt_num_cells = 80,
            num_rays = 200,
            num_leaves = 64,
            num_phases = None,
            phases = None,
            do_fast = False,
            fast_sqrt_num_cells = 16,
            fast_min_sqrt_num_cells = 4,
            fast_max_sqrt_num_cells = 16,
            fast_num_rays = 100,
            fast_num_leaves = 32,
            fast_num_phases = None,
            fast_phases = None,
            is_antiphased = False,
            custom = None,
            image_order_limit = None,
            atm_ext='Num5D',
            **kwargs
            ):
 

        

        if atm_ext == 'Num5D':
            doc = """
            tbb
            """
            super_tbb = Parameter('super_tbb',
                        strict_bounds = (0.001, 0.003), #tbb = Tbb(keV)/511keV
                        bounds = bounds.get('super_tbb', None),
                        doc = doc,
                        symbol = r'tbb',
                        value = values.get('super_tbb', None))
            doc = """
            te
            """
            super_te = Parameter('super_te',
                        strict_bounds = (40., 200.), #te = Te(keV)*1000/511keV
                        bounds = bounds.get('super_te', None),
                        doc = doc,
                        symbol = r'te',
                        value = values.get('super_te', None))
            
            doc = """
            tau
            """
            super_tau = Parameter('super_tau',
                        strict_bounds = (0.5, 3.5),
                        bounds = bounds.get('super_tau', None),
                        doc = doc,
                        symbol = r'tau',
                        value = values.get('super_tau', None))
            

            custom = [super_tbb, super_te, super_tau]
        elif atm_ext == 'user':
            doc = """
            h
            """
            super_h = Parameter('super_h',
                        strict_bounds = (-1, 1),  
                        bounds = bounds.get('super_h', None),
                        doc = doc,
                        symbol = r'h',
                        value = values.get('super_h', None))

            custom = [super_h]

        if cede:
            doc = """
            cede_tbb
            """        
            cede_tbb = Parameter('cede_tbb',
            strict_bounds = (0.001, 0.003),
            bounds = bounds.get('cede_tbb', None),
            doc = doc,
            symbol = r'cede_tbb',
            value = values.get('cede_tbb', None))

            doc = """
            cede_te
            """
            cede_te = Parameter('cede_te',
                        strict_bounds = (40., 200.),
                        bounds = bounds.get('cede_te', None),
                        doc = doc,
                        symbol = r'cede_te',
                        value = values.get('cede_te', None))
            
            doc = """
            cede_tau
            """
            cede_tau = Parameter('cede_tau',
                        strict_bounds = (0.5, 3.5),
                        bounds = bounds.get('cede_tau', None),
                        doc = doc,
                        symbol = r'cede_tau',
                        value = values.get('cede_tau', None))
            
            custom += [cede_tbb,cede_te,cede_tau]

        super(CustomHotRegion_Accreting, self).__init__(
                bounds,
                values,
                symmetry = symmetry,
                interpolator = interpolator,
                omit = omit,
                cede = cede,
                concentric = concentric,
                sqrt_num_cells = sqrt_num_cells,
                min_sqrt_num_cells = min_sqrt_num_cells,
                max_sqrt_num_cells = max_sqrt_num_cells,
                num_rays = num_rays,
                num_leaves = num_leaves,
                num_phases = num_phases,
                phases = phases,
                do_fast = do_fast,
                fast_sqrt_num_cells = fast_sqrt_num_cells,
                fast_min_sqrt_num_cells = fast_min_sqrt_num_cells,
                fast_max_sqrt_num_cells = fast_max_sqrt_num_cells,
                fast_num_rays = fast_num_rays,
                fast_num_leaves = fast_num_leaves,
                fast_num_phases = fast_num_phases,
                fast_phases = fast_phases,
                is_antiphased = is_antiphased,
                custom = custom,
                image_order_limit = image_order_limit,
                atm_ext=atm_ext,
                **kwargs
                )

    def _HotRegion__compute_cellParamVecs(self):


        self._super_radiates = np.greater(self._super_cellArea, 0.0).astype(np.int32)
        self._super_cellParamVecs = np.ones((self._super_radiates.shape[0],
                                              self._super_radiates.shape[1],
                                              3),
                                             dtype=np.double)
        if self.atm_ext == 5: # user
            self._super_cellParamVecs[...,:-1] *= self['super_h']

 
      
        elif self.atm_ext == 6: # Bobrikova atmosphere
            self._super_cellParamVecs[...,0] *= self['super_te']
            self._super_cellParamVecs[...,1] *= self['super_tbb']
            self._super_cellParamVecs[...,2] *= self['super_tau']

        else:
            print('error in cellparamvecs')

        #try:
        #    self._cede_radiates = np.greater(self._cede_cellArea, 0.0).astype(np.int32)
        #except AttributeError:
        #    pass
        #else:
        #    self._cede_cellParamVecs = np.ones((self._cede_radiates.shape[0],
        #                                         self._cede_radiates.shape[1],
        #                                         3), dtype=np.double)
        #
        #    self._cede_cellParamVecs[...,0] *= self['cede_te']
        #    self._cede_cellParamVecs[...,1] *= self['cede_tbb']
        #    self._cede_cellParamVecs[...,2] *= self['cede_tau']

    # @property
    # def symmetry(self):
    #     """ Get the symmetry declaration (controls integrator invocation). """
    #     return self._symmetry

    # @symmetry.setter
    # def symmetry(self, declaration):
    #     if not isinstance(declaration, bool):
    #         raise TypeError('Declare symmetry existence with a boolean.')

    #     self._symmetry = declaration
    #     # find the required integrator
    #     if declaration: # can we safely assume azimuthal invariance?
    #         if self._split:
    #             if not self._disk_blocking:
    #                 from xpsi.cellmesh.integrator_for_azimuthal_invariance_split import integrate as _integrator
    #                 from xpsi.cellmesh.integratorIQU_for_azimuthal_invariance_split import integrate as _integratorIQU 
    #             elif self._disk_blocking:
    #                 from xpsi.cellmesh.integrator_for_azimuthal_invariance_split_disk import integrate as _integrator
    #                 # print("xpsi.cellmesh.integrator_for_azimuthal_invariance_split_disk import integrate as _integrator")
    #                 from xpsi.cellmesh.integratorIQU_for_azimuthal_invariance_split_disk import integrate as _integratorIQU
    #                 # print("xpsi.cellmesh.integratorIQU_for_azimuthal_invariance_split_disk import integrate as _integrator")

    #         else:
    #             if not self._disk_blocking:
    #                 from xpsi.cellmesh.integrator_for_azimuthal_invariance import integrate as _integrator
    #                 from xpsi.cellmesh.integratorIQU_for_azimuthal_invariance import integrate as _integratorIQU
    #             elif self._disk_blocking:
    #                 from xpsi.cellmesh.integrator_for_azimuthal_invariance_disk import integrate as _integrator
    #                 from xpsi.cellmesh.integratorIQU_for_azimuthal_invariance_disk import integrate as _integratorIQU

    #     else: # more general purpose
    #         if self._split:
    #             raise TypeError("Split version of the integrator has not been implemented for symmetry=False.")
    #         from xpsi.cellmesh.integrator import integrate as _integrator
    #         from xpsi.cellmesh.integratorIQU import integrate as _integratorIQU
    #     self._integrator = _integrator
    #     self._integratorIQU = _integratorIQU


    # def integrate(self, st, energies, threads, # here photosphere should also send the disk
    #               hot_atmosphere, elsewhere_atmosphere, atm_ext_else, R_in = None):
    #     # print("R_in = ", R_in)
    #     """ Integrate over the photospheric radiation field.

    #     Calls the CellMesh integrator, with or without exploitation of
    #     azimuthal invariance of the radiation field of the hot region.

    #     :param st: Instance of :class:`~.Spacetime.Spacetime`.

    #     :param energies: A one-dimensional :class:`numpy.ndarray` of energies
    #                      in keV.

    #     :param int threads: Number of ``OpenMP`` threads for pulse
    #                         integration.

    #     """
    #     if self.fast_mode and not self.do_fast:
    #         try:
    #             if self.cede:
    #                 return (None, None)
    #         except AttributeError:
    #             return (None,)

    #     leaves = self._fast_leaves if self.fast_mode else self._leaves
    #     phases = self._fast_phases if self.fast_mode else self._phases
    #     num_rays = self._fast_num_rays if self.fast_mode else self._num_rays

    #     if isinstance(energies, tuple):
    #         try:
    #             super_energies, cede_energies = energies
    #         except ValueError:
    #             super_energies = energies[0]
    #             try:
    #                 self._cede_cellArea
    #             except AttributeError:
    #                 pass
    #             else:
    #                 cede_energies = super_energies
    #     else:
    #         super_energies = cede_energies = energies

    #     if self.atm_ext==2 or self.atm_ext==4 or self._split:
    #         if hot_atmosphere == ():
    #             raise AtmosError('The numerical atmosphere data were not preloaded, '
    #                              'even though that is required by the current atmosphere extension.')
            
    #     if self._disk_blocking:
    #          super_pulse = self._integrator(threads,    
    #                                         R_in,    # here R_in is added                   
    #                                         st.R, 
    #                                         st.Omega,
    #                                         st.r_s,
    #                                         st.i,
    #                                         self._super_cellArea,
    #                                         self._super_r,
    #                                         self._super_r_s_over_r,
    #                                         self._super_theta,
    #                                         self._super_phi,
    #                                         self._super_cellParamVecs,
    #                                         self._super_radiates,
    #                                         self._super_correctionVecs,
    #                                         num_rays,
    #                                         self._super_deflection,
    #                                         self._super_cos_alpha,
    #                                         self._super_lag,
    #                                         self._super_maxDeflection,
    #                                         self._super_cos_gamma,
    #                                         super_energies,
    #                                         leaves,
    #                                         phases,
    #                                         hot_atmosphere,
    #                                         elsewhere_atmosphere,
    #                                         self.atm_ext,
    #                                         atm_ext_else,
    #                                         self.beam_opt,
    #                                         self._image_order_limit)
    #     elif not self._disk_blocking:
    #          super_pulse = self._integrator(threads,                       
    #                                         st.R,
    #                                         st.Omega,
    #                                         st.r_s,
    #                                         st.i,
    #                                         self._super_cellArea,
    #                                         self._super_r,
    #                                         self._super_r_s_over_r,
    #                                         self._super_theta,
    #                                         self._super_phi,
    #                                         self._super_cellParamVecs,
    #                                         self._super_radiates,
    #                                         self._super_correctionVecs,
    #                                         num_rays,
    #                                         self._super_deflection,
    #                                         self._super_cos_alpha,
    #                                         self._super_lag,
    #                                         self._super_maxDeflection,
    #                                         self._super_cos_gamma,
    #                                         super_energies,
    #                                         leaves,
    #                                         phases,
    #                                         hot_atmosphere,
    #                                         elsewhere_atmosphere,
    #                                         self.atm_ext,
    #                                         atm_ext_else,
    #                                         self.beam_opt,
    #                                         self._image_order_limit)


    #     if super_pulse[0] == 1:
    #         raise PulseError('Fatal numerical error during superseding-'
    #                          'region pulse integration.')

    #     try:
    #         cede_pulse = self._integrator(threads,
    #                                       st.R,
    #                                       st.Omega,
    #                                       st.r_s,
    #                                       st.i,
    #                                       self._cede_cellArea,
    #                                       self._cede_r,
    #                                       self._cede_r_s_over_r,
    #                                       self._cede_theta,
    #                                       self._cede_phi,
    #                                       self._cede_cellParamVecs,
    #                                       self._cede_radiates,
    #                                       self._cede_correctionVecs,
    #                                       num_rays,
    #                                       self._cede_deflection,
    #                                       self._cede_cos_alpha,
    #                                       self._cede_lag,
    #                                       self._cede_maxDeflection,
    #                                       self._cede_cos_gamma,
    #                                       cede_energies,
    #                                       leaves,
    #                                       phases,
    #                                       hot_atmosphere,
    #                                       elsewhere_atmosphere,
    #                                       self.atm_ext,
    #                                       atm_ext_else,
    #                                       self.beam_opt,
    #                                       self._image_order_limit)
    #     except AttributeError:
    #         pass
    #     else:
    #         if cede_pulse[0] == 1:
    #             raise PulseError('Fatal numerical error during ceding-region '
    #                              'pulse integration.')
    #         else:
    #             return (super_pulse[1], cede_pulse[1])
            
    
    #     return (super_pulse[1],)
    


    # def integrate_stokes(self, st, energies, threads,
    #               hot_atmosphere_I, hot_atmosphere_Q, elsewhere_atmosphere, atm_ext_else, R_in = None):
    #     """ Integrate Stokes parameters over the photospheric radiation field.

    #     Calls the CellMesh Stokes integrators, with or without exploitation of
    #     azimuthal invariance of the radiation field of the hot region.

    #     :param st: Instance of :class:`~.Spacetime.Spacetime`.

    #     :param energies: A one-dimensional :class:`numpy.ndarray` of energies
    #                      in keV.

    #     :param int threads: Number of ``OpenMP`` threads for pulse
    #                         integration.

    #     """
    #     if self.fast_mode and not self.do_fast:
    #         try:
    #             if self.cede:
    #                 return (None, None)
    #         except AttributeError:
    #             return (None,)

    #     leaves = self._fast_leaves if self.fast_mode else self._leaves
    #     phases = self._fast_phases if self.fast_mode else self._phases
    #     num_rays = self._fast_num_rays if self.fast_mode else self._num_rays

    #     if isinstance(energies, tuple):
    #         try:
    #             super_energies, cede_energies = energies
    #         except ValueError:
    #             super_energies = energies[0]
    #             try:
    #                 self._cede_cellArea
    #             except AttributeError:
    #                 pass
    #             else:
    #                 cede_energies = super_energies
    #     else:
    #         super_energies = cede_energies = energies

    #     if self.atm_ext==2 or self.atm_ext==4 or self._split:
    #         if hot_atmosphere_I == () or hot_atmosphere_Q == ():
    #             raise AtmosError('The numerical atmosphere data were not preloaded, '
    #                              'even though that is required by the current atmosphere extension.')
            
    #     if self._disk_blocking:
    #         all_pulses = self._integratorIQU(threads,
    #                                 R_in,
    #                                 st.R,
    #                                 st.Omega,
    #                                 st.r_s,
    #                                 st.i,
    #                                 self._super_cellArea,
    #                                 self._super_r,
    #                                 self._super_r_s_over_r,
    #                                 self._super_theta,
    #                                 self._super_phi,
    #                                 self._super_cellParamVecs,
    #                                 self._super_radiates,
    #                                 self._super_correctionVecs,
    #                                 num_rays,
    #                                 self._super_deflection,
    #                                 self._super_cos_alpha,
    #                                 self._super_lag,
    #                                 self._super_maxDeflection,
    #                                 self._super_cos_gamma,
    #                                 super_energies,
    #                                 leaves,
    #                                 phases,
    #                                 hot_atmosphere_I,
    #                                 hot_atmosphere_Q,
    #                                 elsewhere_atmosphere,
    #                                 self.atm_ext,
    #                                 atm_ext_else,
    #                                 self.beam_opt,
    #                                 self._image_order_limit)
                                    
            
    #     if not self._disk_blocking:
    #         all_pulses = self._integratorIQU(threads,
    #                                 st.R,
    #                                 st.Omega,
    #                                 st.r_s,
    #                                 st.i,
    #                                 self._super_cellArea,
    #                                 self._super_r,
    #                                 self._super_r_s_over_r,
    #                                 self._super_theta,
    #                                 self._super_phi,
    #                                 self._super_cellParamVecs,
    #                                 self._super_radiates,
    #                                 self._super_correctionVecs,
    #                                 num_rays,
    #                                 self._super_deflection,
    #                                 self._super_cos_alpha,
    #                                 self._super_lag,
    #                                 self._super_maxDeflection,
    #                                 self._super_cos_gamma,
    #                                 super_energies,
    #                                 leaves,
    #                                 phases,
    #                                 hot_atmosphere_I,
    #                                 hot_atmosphere_Q,
    #                                 elsewhere_atmosphere,
    #                                 self.atm_ext,
    #                                 atm_ext_else,
    #                                 self.beam_opt,
    #                                 self._image_order_limit)            

    #     super_pulse = all_pulses[0], all_pulses[1]
    #     super_pulse_Q = all_pulses[0], all_pulses[2]
    #     super_pulse_U = all_pulses[0], all_pulses[3]


    #     if super_pulse[0] == 1:
    #         raise PulseError('Fatal numerical error during superseding-'
    #                          'region pulse integration.')

    #     try:
  
    #         all_pulses = self._integratorIQU(threads,
    #                                    st.R,
    #                                    st.Omega,
    #                                    st.r_s,
    #                                    st.i,
    #                                    self._cede_cellArea,
    #                                    self._cede_r,
    #                                    self._cede_r_s_over_r,
    #                                    self._cede_theta,
    #                                    self._cede_phi,
    #                                    self._cede_cellParamVecs,
    #                                    self._cede_radiates,
    #                                    self._cede_correctionVecs,
    #                                    num_rays,
    #                                    self._cede_deflection,
    #                                    self._cede_cos_alpha,
    #                                    self._cede_lag,
    #                                    self._cede_maxDeflection,
    #                                    self._cede_cos_gamma,
    #                                    cede_energies,
    #                                    leaves,
    #                                    phases,
    #                                    hot_atmosphere_I,
    #                                    hot_atmosphere_Q, 
    #                                    elsewhere_atmosphere,
    #                                    self.atm_ext,
    #                                    atm_ext_else,
    #                                    self.beam_opt,
    #                                    self._image_order_limit)
    #         cede_pulse = all_pulses[0], all_pulses[1] #success and flux
    #         cede_pulse_Q = all_pulses[0], all_pulses[2]
    #         cede_pulse_U = all_pulses[0], all_pulses[3]

   
    #     except AttributeError:
    #         pass
    #     else:
    #         if cede_pulse[0] == 1:
    #             raise PulseError('Fatal numerical error during ceding-region '
    #                              'pulse integration.')
    #         else:
    #             return (super_pulse[1], cede_pulse[1]), (super_pulse_Q[1], cede_pulse_Q[1]), (super_pulse_U[1], cede_pulse_U[1])
    #     return (super_pulse[1],), (super_pulse_Q[1],), (super_pulse_U[1],)
