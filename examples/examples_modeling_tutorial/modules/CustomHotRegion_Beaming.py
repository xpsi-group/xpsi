
import xpsi
from xpsi.global_imports import *
from xpsi.Parameter import Parameter

# Your derived class can have an initializer that sets up your new parameters and the temperature parameters, and then passes them to the base class initializer using the super(<derived class>, self).__init__(*args, **kwargs) call, where kwargs includes custom parameter

class CustomHotRegion(xpsi.HotRegion):
    """ Implement method for HotRegion with beaming parameters."""
    required_names = ['super_colatitude',
                      'super_radius',
                      'phase_shift',
                      'super_temperature (if no custom specification)']

    optional_names = ['super_abb',
                      'super_bbb',
                      'super_cbb',
                      'super_dbb',    
                      'omit_colatitude',
                      'omit_radius',
                      'omit_azimuth',
                      'cede_colatitude',
                      'cede_radius',
                      'cede_azimuth',
                      'cede_temperature']

    def __init__(self,
                 bounds,
                 values,
                 symmetry = True,
                 use_disk=True,
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
                 fbeam = False,
                 **kwargs):
        
        self.use_disk=use_disk
      

        doc = """
        log10(superseding region effective temperature [K])
        """
        super_temp = Parameter('super_temperature',
	          strict_bounds = (3.0, 7.0), # very cold --> very hot
	          bounds = bounds.get('super_temperature', None),
	          doc = doc,
	          symbol = r'$\log_{10}(T\;[\rm{K}])$',
	          value = values.get('super_temperature', None))
        if cede:
            doc = """
            log10(ceding region effective temperature [K])
            """
            cede_temp = Parameter('cede_temperature',
                      strict_bounds = (3.0, 7.0), # same story
                      bounds = bounds.get('cede_temperature', None),
                      doc = doc,
                      symbol = r'$\log_{10}(\mathcal{T}\;[\rm{K}])$',
                      value = values.get('cede_temperature', None))
        else:
            cede_temp = None

                
                
            if not fbeam: 
                bounds['super_abb'] = None
                values['super_abb'] = 0.0
                no_verb['super_abb'] = True
                bounds['super_bbb'] = None
                values['super_bbb'] = 0.0
                no_verb['super_bbb'] = True
                bounds['super_cbb'] = None
                values['super_cbb'] = 0.0
                no_verb['super_cbb'] = True
                bounds['super_dbb'] = None
                values['super_dbb'] = 0.0
                no_verb['super_dbb'] = True            
                
            doc = """
            abb
            """
            super_abb = Parameter('super_abb',
                     strict_bounds = (-3.0, 3.0),
                     bounds = bounds.get('super_abb', None),
                     doc = doc,
                     symbol = r'$abb$',
                     value = values.get('super_abb', None))
            doc = """
            bbb
            """
            super_bbb = Parameter('super_bbb',
                     strict_bounds = (-3.0, 3.0),
                     bounds = bounds.get('super_bbb', None),
                     doc = doc,
                     symbol = r'$\bbb$',
                     value = values.get('super_bbb', None))   
            doc = """
            cbb
            """
            super_cbb = Parameter('super_cbb',
                     strict_bounds = (-3.0, 3.0),
                     bounds = bounds.get('super_cbb', None),
                     doc = doc,
                     symbol = r'$\cbb$',
                     value = values.get('super_cbb', None))  
            doc = """
            dbb
            """
            super_dbb = Parameter('super_dbb',
                     strict_bounds = (-3.0, 3.0),
                     bounds = bounds.get('super_dbb', None),
                     doc = doc,
                     symbol = r'$\dbb$',
                     value = values.get('super_dbb', None))     



        if cede:
            custom = [super_temp,cede_temp,super_abb, super_bbb, super_cbb, super_dbb]
        else:
            custom = [super_temp, super_abb, super_bbb, super_cbb, super_dbb]	

            super(CustomHotRegion, self).__init__(
                     bounds,
                     values,
                     symmetry = symmetry,
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
                     **kwargs)

    def _HotRegion__compute_cellParamVecs(self):
        """
        Precompute photospheric source radiation field parameter vectors
        cell-by-cell. Free model parameters and derived (fixed) variables can
        be transformed into local comoving radiation field variables.

        Subclass and overwrite with custom functionality if you desire.

        Designed here simply for uniform effective temperature superseding (with free beaming parameters)
        and ceding regions.

        """
        self._super_radiates = _np.greater(self._super_cellArea, 0.0).astype(_np.int32)
        self._super_cellParamVecs = _np.ones((self._super_radiates.shape[0],
                                              self._super_radiates.shape[1],
                                              7),
                                             dtype=_np.double)

        self._super_cellParamVecs[...,0] *= self['super_temperature']

        for i in range(self._super_cellParamVecs.shape[1]):
            self._super_cellParamVecs[:,i,1] *= self._super_effGrav

        self._super_cellParamVecs[...,2] *= self['super_abb']
        self._super_cellParamVecs[...,3] *= self['super_bbb']
        self._super_cellParamVecs[...,4] *= self['super_cbb']
        self._super_cellParamVecs[...,5] *= self['super_dbb']

        #Number of mus in re-normalization integral:
        self._super_cellParamVecs[...,6] *= 10 #100

        try:
            self._cede_radiates = _np.greater(self._cede_cellArea, 0.0).astype(_np.int32)
        except AttributeError:
            pass
        else:
            self._cede_cellParamVecs = _np.ones((self._cede_radiates.shape[0],
                                                 self._cede_radiates.shape[1],
                                                 2), dtype=_np.double)

            self._cede_cellParamVecs[...,:-1] *= self['cede_temperature']

            for i in range(self._cede_cellParamVecs.shape[1]):
                self._cede_cellParamVecs[:,i,-1] *= self._cede_effGrav
        #print(self._super_cellParamVecs)
        #exit()
