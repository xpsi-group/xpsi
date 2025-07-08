from xpsi.global_imports import *

from xpsi.Parameter import Parameter
from xpsi.EmissionModels import EmissionModel

class PowerLaw( EmissionModel ):

    required_names = ['powerlaw_norm',
                      'powerlaw_gamma']

    optional_names = ['powerlaw_oscillation_amplitude',
                      'powerlaw_phase_shift']
    
    def __init__(self, 
                 bounds,
                 values,
                 pulsed=False,
                 *args, **kwargs):

        # Initiate the necessary parameters
        norm = Parameter('powerlaw_norm',
                        strict_bounds = (0.0, _np.infty),
                        bounds = bounds.get('powerlaw_norm', None),
                        doc = "Normalization of the PowerLaw",
                        symbol = 'norm',
                        value = values.get('powerlaw_norm', None))
        
        gamma = Parameter('powerlaw_gamma',
                strict_bounds = (0.0, _np.infty),
                bounds = bounds.get('powerlaw_gamma', None),
                doc = "Power of the PowerLaw",
                symbol = r'$\Gamma$',
                value = values.get('powerlaw_gamma', None))
        
        # Add more if pulsed
        self.is_pulsed = pulsed
        if self.is_pulsed:

            oscillation_amplitude = Parameter('powerlaw_oscillation_amplitude',
                    strict_bounds = (0. , _np.infty),
                    bounds = bounds.get('powerlaw_oscillation_amplitude', None),
                    doc = "Amplitude of oscillation of the PowerLaw",
                    symbol = 'A',
                    value = values.get('powerlaw_oscillation_amplitude', None))
            
            phase_shift = Parameter('powerlaw_phase_shift',
                    strict_bounds = (0. , 1.),
                    bounds = bounds.get('powerlaw_phase_shift', None),
                    doc = "Phase shift of the PowerLaw",
                    symbol = r'$\phi_0$',
                    value = values.get('powerlaw_phase_shift', None))

        else:
            oscillation_amplitude = None
            phase_shift = None

        # Initiate the parent class
        super(PowerLaw, self).__init__( norm, gamma, oscillation_amplitude, phase_shift,  *args, **kwargs )
        
    def intensity(self, E, phase ):
        if self.is_pulsed:
            return ( self['norm'] + self['A'] * _np.sin( phase - self['phi0'] ) ) * E^( -self['Gamma'] )
        else:
            return self['norm'] * E^( -self['Gamma'] )
        
    def integrate(self,
                energy_edges,
                phases_edges,
                is_pulsed = False):

        # Energy term
        energy_term = (1. - self['Gamma']) * ( energy_edges[1:]**(1. - self['Gamma']) - energy_edges[:-1]**(1. - self['Gamma']) )
        phases_width = phases_edges[1:] - phases_edges[:-1]
        normalization_term = self['norm'] * 2. * _np.pi * phases_width
        if is_pulsed:
            normalization_term += self['A'] * ( _np.sin( 2.*_np.pi*(phases_edges[1:]-self['phi0']) ) - _np.sin( 2.*_np.pi*(phases_edges[:-1]-self['phi0']) ) )

        # Get the integrated powerlaw
        integrated_powerlaw = _np.dot( normalization_term[:,_np.newaxis], energy_term[_np.newaxis,:] )
        return integrated_powerlaw
    
    # No need to embed but need to override
    def embed(self):
        pass
