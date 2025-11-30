#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  6 16:46:36 2024

@author: bas
"""

import xpsi
from xpsi import Parameter
import numpy as np

class CustomHotRegion_Accreting(xpsi.HotRegion):
    """Custom implementation of HotRegion. Accreting Atmosphere model by 
    Anna Bobrikova if atm_ext = 'Num5D'. The parameters are ordered 
    I(E < mu < tau < tbb < te).
    
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
