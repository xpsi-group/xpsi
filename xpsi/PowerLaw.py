from xpsi.global_imports import *

from xpsi.Parameter import Parameter
from xpsi.EmissionModels import EmissionModel

class PowerLaw( EmissionModel ):

    required_names = ['norm',
                      'gamma']

    optional_names = ['oscillation_amplitude',
                      'phase_shift']
    
    def __init__(self, 
                 bounds,
                 values,
                 pulsed=False,
                 **kwargs):

        # Initiate the necessary parameters
        norm = Parameter('norm',
                        strict_bounds = (0.0, _np.infty),
                        bounds = bounds.get('norm', None),
                        doc = "Normalization of the PowerLaw",
                        symbol = 'norm',
                        value = values.get('norm', None))
        
        gamma = Parameter('gamma',
                strict_bounds = (0.0, _np.infty),
                bounds = bounds.get('gamma', None),
                doc = "Power of the PowerLaw",
                symbol = r'$\gamma$',
                value = values.get('gamma', None))
        
        # Add more if pulsed
        self.is_pulsed = pulsed
        if self.is_pulsed:

            oscillation_amplitude = Parameter('oscillation_amplitude',
                    strict_bounds = (0. , _np.infty),
                    bounds = bounds.get('oscillation_amplitude', None),
                    doc = "Amplitude of oscillation of the PowerLaw",
                    symbol = 'A',
                    value = values.get('oscillation_amplitude', None))
            phase_shift_bounds = bounds.get('phase_shift', None)
            phase_shift_value = values.get('phase_shift', None)

        else:
            oscillation_amplitude = None
            phase_shift_bounds = None
            phase_shift_value = 0.

        phase_shift = Parameter('phase_shift',
                strict_bounds = (0. , 1.),
                bounds = phase_shift_bounds,
                doc = "Phase shift of the PowerLaw",
                symbol = r'$\phi_0$',
                value = phase_shift_value)

        # Initiate the parent class
        super(PowerLaw, self).__init__( norm, gamma, oscillation_amplitude, phase_shift, **kwargs ) #prefix in kwargs
        
    # def intensity(self, E, phase ):
    #     if self.is_pulsed:
    #         return ( self['norm'] + self['A'] * _np.sin( phase - self['phi0'] ) ) * E^( -self['gamma'] )
    #     else:
    #         return self['norm'] * E^( -self['gamma'] )
    

    def integrate(self, energies, threads, *args, **kwargs):

        table = self['norm'] * _np.ones((self._num_leaves, len(energies))) * energies**( -self['gamma'] )
        self.signal = (table,)
        print('Integrated PowerLaw')

    # def integrate_with_bin_edges(self,
    #             energy_edges,
    #             phases_edges,
    #             is_pulsed = False):

    #     # Energy term
    #     energy_term = (1. - self['gamma']) * ( energy_edges[1:]**(1. - self['gamma']) - energy_edges[:-1]**(1. - self['gamma']) )
    #     phases_width = phases_edges[1:] - phases_edges[:-1]
    #     normalization_term = self['norm'] * 2. * _np.pi * phases_width
    #     if is_pulsed:
    #         normalization_term += self['A'] * ( _np.sin( 2.*_np.pi*(phases_edges[1:]-self['phi0']) ) - _np.sin( 2.*_np.pi*(phases_edges[:-1]-self['phi0']) ) )

    #     # Get the integrated powerlaw
    #     integrated_powerlaw = _np.dot( normalization_term[:,_np.newaxis], energy_term[_np.newaxis,:] )
    #     return integrated_powerlaw
