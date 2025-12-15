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

            # Define oscillation amplitude parameter
            oscillation_amplitude = Parameter('oscillation_amplitude',
                    strict_bounds = (0. , _np.infty),
                    bounds = bounds.get('oscillation_amplitude', None),
                    doc = "Amplitude of oscillation of the PowerLaw",
                    symbol = 'A',
                    value = values.get('oscillation_amplitude', None))
            
            # Define phase shift parameter values
            phase_shift_bounds = bounds.get('phase_shift', None)
            phase_shift_value = values.get('phase_shift', None)

            # Define the integration method
            self.integrate = self.integrate_pulsed

        else:

            # Define oscillation amplitude and phase shift parameters values
            oscillation_amplitude = None
            phase_shift_bounds = None
            phase_shift_value = 0.

            # Define the integration method
            self.integrate = self.integrate_phase_independent

        # Define phase shift
        phase_shift = Parameter('phase_shift',
                strict_bounds = (0. , 1.),
                bounds = phase_shift_bounds,
                doc = "Phase shift of the PowerLaw",
                symbol = r'$\phi_0$',
                value = phase_shift_value)
        
        # Set the value of num_leaves
        self._num_leaves = kwargs.pop('num_leaves', 32) if pulsed else 5 # Default to 5 so that energy integrator works
        self._phases =  _np.linspace(0.0, 1.0, int(self._num_leaves))

        # Initiate the parent class
        super(PowerLaw, self).__init__( norm, gamma, oscillation_amplitude, phase_shift, **kwargs ) #prefix in kwargs


    def integrate_phase_independent(self, energies, threads, *args, **kwargs):

        table = self['norm'] * _np.ones((self._num_leaves, len(energies))) * energies**( -self['gamma'] )
        self.signal = (_np.ascontiguousarray(table.T),)

    def integrate_pulsed(self, energies, threads, *args, **kwargs):

        # Make a table of the phases
        phases_array = 1 + self['oscillation_amplitude'] * _np.cos( _2pi * ( self._phases - self['phase_shift'] ) )
        powerlaw_array = energies**( -self['gamma'] )
        table = self['norm'] * _np.outer(phases_array, powerlaw_array)
        self.signal = (_np.ascontiguousarray(table.T),)
