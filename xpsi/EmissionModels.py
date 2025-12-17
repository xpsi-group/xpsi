from xpsi.global_imports import *

from xpsi.Parameter import Parameter
from xpsi.ParameterSubspace import ParameterSubspace
from abc import ABCMeta

class EmissionModel( ParameterSubspace , metaclass=ABCMeta ):
    """Base class for any emission model.

    An emission model needs to have an ``integrate`` method defined to compute the flux in [photon/cm^2/s].
    An ``embed`` method to embed the model into the ambient spacetime can be setup if needed.
    Otherwise it can be left blank.

    """
    def __init__(self, *args, **kwargs):

        # First load the parameters
        super( EmissionModel, self ).__init__( *args, **kwargs )

        # Make the phases if not done before
        if not hasattr(self, '_phases'):
            self._num_leaves = kwargs.pop('num_leaves', 32)
            self._phases =  _np.linspace(0.0, 1.0, int(self._num_leaves))

        # Check if everything needed is well defined
        assert hasattr(self, '_num_leaves'), 'Number of leaves must be defined.'
        assert hasattr(self, '_phases'), 'Array of phases used for computations must be defined.'
        assert isinstance(self._num_leaves, int), 'Number of leaves must be an integer.' 
        assert self._num_leaves > 5, 'Number of leaves must be greater than 5 for energy integrator.'

    def integrate(self, energies, threads, *args, **kwargs):
        pass

    def embed(self, spacetime, threads):
        pass

    @property
    def signal(self):
        return self._signal

    @signal.setter
    def signal(self, value):
        self._signal = value

    @property
    def phases(self):
        return self._phases

    @phases.setter
    def phases(self, value):
        assert isinstance(value, _np.ndarray), 'Invalid type for phases. Must be a numpy array.'
        assert value.ndim == 1, 'Phases must form a one-dimensional sequence.'
        assert value.shape[0] == self._num_leaves, 'Number of phases must match number of leaves.'
        self._phases = value

    @property
    def parameters(self):
        """ Like Prior, this is to access other parameters such as Spacetime parameters if needed. """
        return self._parameters

    @parameters.setter
    def parameters(self, obj):
        if not isinstance(obj, ParameterSubspace):
            raise TypeError('A ParameterSubspace object is required.')
        self._parameters = obj

class EmissionModels( ParameterSubspace ):
    """ A collection of emission models, used for interfacing with the Likelihood and signal registering. 

    :param list models:
        A list of emission models or an emission model.

    """
    def __init__( self, 
                 models ):
        
        # Everything must be an EmissionModel
        try:
            if isinstance( models, EmissionModel ):
                models = [models]
            elif isinstance( models, tuple ) or isinstance( models, list):
                for model in models:
                    assert isinstance(model, EmissionModel), 'Invalid type for emission model. Must be EmissionModel.'
        except:
            raise TypeError('models must be an iterable of EmissionModel.')
        
        # Save it all
        self.models = models

        # Initialize
        super( EmissionModels, self ).__init__( models )

    def integrate(self, energies, threads,
                  *args, **kwargs):
        """ Integrate over the emission models (same as for HotRegions)"""

        for model in self.models:
            model.integrate( energies,
                            threads,
                            *args, **kwargs )
    
    def update(self, threads=1, force_update=False):
        """ Embed the emission models into the ambient spacetime. """

        # TODO link this to the spacetime
        for model in self.models:
            if model.needs_update or force_update:
                model.embed( 1, threads )

    @property
    def signal(self):
        return tuple( model.signal for model in self.models )
    
    def __iter__(self):
        return iter(self.models)

class PowerLaw( EmissionModel ):

    """ Simple powerlaw emission model, which can be pulsed or not.

    $$F_E = K E^{-\Gamma} \left( 1 + A\cos(2\pi (\phi - \phi_0) ) \right)$$ 

    where F_E is the flux in [photon/cm^2/s] at energy E, K is the normalization
    (value of the unabsorbed powerlaw in [photon/cm^2/s] at 1keV), $\Gamma$ is the
    power of the powerlaw, A is the pulsation amplitude, and $\phi_0$ is the pulsation phase.

    :param dict bounds:
        Hard prior parameter bounds for the free parameters. The dictionary
        keys must match the required parameter names, at least. If a required
        name is omitted as a key, the parameter is interpreted as *fixed* or
        *derived*. A key-value pair can take the following forms:

            * ``'name': None``
            * ``'name': (None, None)``, ``(None, x)``, ``(x, None)``
            * ``'name': (x, y)``

        where if a bound is ``None`` that bound is set equal to a strict
        hard-coded bound. We note that the bounds for parameters used in the
        atmosphere model should be restricted (by the user) to be within the
        tabulated values, in case a numerical atmosphere extension is used.

    :param dict values:
        Initial values of *free* parameters, fixed values of *fixed* parameters,
        and callables for *derived* parameters. If a key is omitted for a free
        parameter, the initial value is ``None`` by default. A key cannot be
        omitted for a required name that appears in the ``bounds`` dictionary
        with value ``None``, or a required name that is omitted from the bounds
        dictionary.

    :param bool pulsed:
        Choose whether the emission is pulsed or not. Setting to ``True`` will
        add the oscillation_amplitude parameter and use a cosine profile.

    """
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
        super(PowerLaw, self).__init__( norm, gamma, oscillation_amplitude, phase_shift, **kwargs ) # prefix in kwargs


    def integrate_phase_independent(self, energies, threads, *args, **kwargs):

        table = self['norm'] * _np.ones((self._num_leaves, len(energies))) * energies**( -self['gamma'] ) / self._num_leaves
        self.signal = (_np.ascontiguousarray(table.T),)

    def integrate_pulsed(self, energies, threads, *args, **kwargs):

        # Make a table of the phases
        phases_array = 1 + self['oscillation_amplitude'] * _np.cos( _2pi * ( self._phases - self['phase_shift'] ) )
        powerlaw_array = energies**( -self['gamma'] )
        table = self['norm'] * _np.outer(phases_array, powerlaw_array) / self._num_leaves
        self.signal = (_np.ascontiguousarray(table.T),)
