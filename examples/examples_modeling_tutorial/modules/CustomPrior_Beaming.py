from __future__ import print_function, division

import numpy as np
import math
from scipy.stats import truncnorm

import xpsi
from xpsi.global_imports import _G, _csq, _km, _2pi, gravradius, _dpr
from xpsi import Parameter

from xpsi.PostProcessing import fix_random_seed
xpsi.PostProcessing.set_random_seed(0) # prevent noise during prior generation

from scipy.interpolate import Akima1DInterpolator

class CustomPrior(xpsi.Prior):
    """ A custom (joint) prior distribution.

    Source: Fictitious
    Model variant: ST-U
        Two single-temperature, simply-connected circular hot regions with mostly
        unshared parameters.
    """

    __derived_names__ = ['compactness','p__phase_shift_shifted','s__phase_shift_shifted']
    __draws_from_support__ = 3

    @fix_random_seed
    def __init__(self):#,photosphere=None):
        
        #Saving the photoshpere for beaming priors if needed:
        #self.photosphere = photosphere
        
        super(CustomPrior, self).__init__()

    def __call__(self, p = None):
        """ Evaluate distribution at ``p``.

        :param list p: Model parameter values.

        :returns: Logarithm of the distribution evaluated at ``p``.

        """            
        temp = super(CustomPrior, self).__call__(p)
        if not np.isfinite(temp):
            return temp

        # based on contemporary EOS theory
        if not self.parameters['radius'] <= 16.0:
            return -np.inf

        ref = self.parameters.star.spacetime # shortcut

        # limit polar radius tobe outside the Schwarzschild photon sphere
        R_p = 1.0 + ref.epsilon * (-0.788 + 1.030 * ref.zeta)
        if R_p < 1.5 / ref.R_r_s:
            return -np.inf

        mu = math.sqrt(-1.0 / (3.0 * ref.epsilon * (-0.788 + 1.030 * ref.zeta)))

        # 2-surface cross-section have a single maximum in |z|
        # i.e., an elliptical surface; minor effect on support, if any,
        # only for high spin frequencies
        if mu < 1.0:
            return -np.inf

        # check effective gravity at pole (where it is maximum) and
        # at equator (where it is minimum) are in NSX limits
        grav = xpsi.surface_radiation_field.effective_gravity(np.array([1.0, 0.0]),
                                                               np.array([ref.R] * 2 ),
                                                               np.array([ref.zeta] * 2),
                                                               np.array([ref.epsilon] * 2))
        for g in grav:
            if not 13.7 <= g <= 15.0:
                return -np.inf

        ref = self.parameters # redefine shortcut

        # enforce order in hot region colatitude
        if ref['p__super_colatitude'] > ref['s__super_colatitude']:
            return -np.inf

        phi = (ref['p__phase_shift'] - 0.5 - ref['s__phase_shift']) * _2pi

        ang_sep = xpsi.HotRegion.psi(ref['s__super_colatitude'],
                                     phi,
                                     ref['p__super_colatitude'])

        # hot regions cannot overlap
        if ang_sep < ref['p__super_radius'] + ref['s__super_radius']:
            return -np.inf
            
        # Prevent negative intensities with free beaming correction:
        abb = self.parameters['p__super_abb']
        bbb = self.parameters['p__super_bbb']
        try:
            cbb = self.parameters['p__super_cbb']
            dbb = self.parameters['p__super_dbb']
        except:
            cbb = 0.0
            dbb = 0.0 
		
        logg = grav[0] #grav[0] = equatorial radius, grav[1] = polar radius         

        K2keV =  8.61732814974493e-08
        if (cbb == 0.0 and dbb == 0.0):
            Elist = np.array([1.0])
            tempps = [6.0]
        else:		
            Elist = np.linspace(-1.5,2.2,num=100) #-1.3,2.0
            tempps = [self.parameters['p__super_temperature'], self.parameters['s__super_temperature']]
	
        for isp in range(0,len(tempps)):
            Ekevp = (10**Elist)*K2keV*(10**tempps[isp])	
            for ie in range(0,len(Ekevp)):
                #Check first if beaming function would be negative, giving negative intensities:
                a = bbb*(Ekevp[ie])**dbb #now named as in f(x)=ax^2+bx+c 
                b = abb*(Ekevp[ie])**cbb
                c = 1.0
                dkr = b**2-4.0*a*c
                if dkr >= 0:
                    if abs(a) < 1e-8:
                        if abs(b) > 1e-8:
                            zp_mu = -c/b
                            if(zp_mu > 0.0 and zp_mu < 1.0):
                                return -np.inf
                    else:
                        zp1_mu = (-b+np.sqrt(dkr))/(2.0*a)
                        zp2_mu = (-b-np.sqrt(dkr))/(2.0*a)
                        if((zp1_mu > 0.0 and zp1_mu < 1.0) or (zp2_mu > 0.0 and zp2_mu < 1.0)):
                            return -np.inf
										
            #Check also if difference compared to exact numerical is higher than we want:
            mus = np.linspace(0.0,1.0,100)
            #enes = np.ones(len(mus))*Ekevp[ie]

            #if using opt=3, err_max should be calculated numerically:
            #atm = self.photosphere.hot_atmosphere
            #opt=0
            #local_vars = np.array([[tempps[isp], abb, bbb, cbb, dbb, opt, 0, logg]]*len(mus))
            #H_nsx = xpsi.surface_radiation_field.intensity(enes, mus, local_vars,
                #                                             atmosphere=atm,
                #                                             extension='hot',
                #                                             numTHREADS=1)
            #opt=3 #This should be the same as in CustomHotRegion.
            #local_vars = np.array([[tempps[isp], abb, bbb, cbb, dbb, opt, 5, logg]]*len(mus))
            #H_nsx_bf = xpsi.surface_radiation_field.intensity(enes, mus, local_vars,
                #                                             atmosphere=atm,
                #                                             extension='hot',
                #                                             numTHREADS=1)  
            #err_max = max(abs(H_nsx-H_nsx_bf)/H_nsx)

            #However, this works for opt=2:
            anorm = 0.5/(0.5+(1.0/3.0)*abb*Ekevp[ie]**cbb+(1.0/4.0)*bbb*Ekevp[ie]**dbb)
            beam_func = anorm*(1.0+abb*(Ekevp[ie]**cbb)*mus+bbb*(Ekevp[ie]**dbb)*mus**2)
            err_max = max(abs(beam_func-1.0))

            if( err_max > 0.1):
                return -np.inf
        return 0.0

    def inverse_sample(self, hypercube=None):
        """ Draw sample uniformly from the distribution via inverse sampling. """

        to_cache = self.parameters.vector

        if hypercube is None:
            hypercube = np.random.rand(len(self))

        # the base method is useful, so to avoid writing that code again:
        _ = super(CustomPrior, self).inverse_sample(hypercube)

        ref = self.parameters # redefine shortcut


        # flat priors in cosine of hot region centre colatitudes (isotropy)
        # support modified by no-overlap rejection condition
        idx = ref.index('p__super_colatitude')
        a, b = ref.get_param('p__super_colatitude').bounds
        a = math.cos(a); b = math.cos(b)
        ref['p__super_colatitude'] = math.acos(b + (a - b) * hypercube[idx])
        
        idx = ref.index('s__super_colatitude')
        a, b = ref.get_param('s__super_colatitude').bounds
        a = math.cos(a); b = math.cos(b)
        ref['s__super_colatitude'] = math.acos(b + (a - b) * hypercube[idx])
        
        
        ###############################
        _scale = 0.3
        idx = ref.index('p__super_abb')        
        ref['p__super_abb'] = truncnorm.ppf(hypercube[idx], -5.0, 5.0,
                                         loc=0.0,
                                         scale=_scale)        

        _scale = 0.1 #0.01 
        idx = ref.index('p__super_bbb')
        abb_val = ref['p__super_abb']
        ref['p__super_bbb'] = truncnorm.ppf(hypercube[idx], -5.0, 5.0,
                                         loc=-abb_val,
                                         scale=_scale) 
        ###############################

        # restore proper cache
        for parameter, cache in zip(self.parameters, to_cache):
            parameter.cached = cache

        # it is important that we return the desired vector because it is
        # automatically written to disk by MultiNest and only by MultiNest
        return self.parameters.vector

    def transform(self, p,**kwargs):
        """ A transformation for post-processing. """

        p = list(p) # copy

        # used ordered names and values
        ref = dict(zip(self.parameters.names, p))

        # compactness ratio M/R_eq
        p += [gravradius(ref['mass']) / ref['radius']]
        
        for phase_shift in ['p__phase_shift', 's__phase_shift']:
            if ref[phase_shift] > 0.5:
                p += [ref[phase_shift] - 1.0]
            else:
                p += [ref[phase_shift]]


        return p
