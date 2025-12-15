from xpsi.global_imports import *

from xpsi.ParameterSubspace import ParameterSubspace
from abc import ABCMeta, abstractmethod

class EmissionModel( ParameterSubspace , metaclass=ABCMeta ):

    def __init__(self, *args, **kwargs):

        # First load the parameters
        super( EmissionModel, self ).__init__( *args, **kwargs )

        # Get some default values if needed
        self._num_leaves = kwargs.pop('num_leaves', None)
        if self._num_leaves is not None:
            self._phases =  _np.linspace(0.0, 1.0, int(self._num_leaves))


    # @abstractmethod
    # def integrate(self, **kwargs):
    #     pass

    # @abstractmethod
    def embed(self, spacetime, threads):
        pass


class EmissionModels( ParameterSubspace ):

    def __init__( self, 
                 models, 
                 *args, **kwargs ):
        
        # Everything must be an EmissionModel
        try:
            if isinstance( models, EmissionModel ):
                models = [models]
            for model in models:
                assert isinstance(model, EmissionModel), 'Invalid type for emission model. Must be EmissionModel.'
        except:
            raise TypeError('emission_models must be an interable of EmissionModel.')
        
        # Save it all
        self.models = models

        # Initialize
        super( EmissionModels, self ).__init__( models, *args, **kwargs )

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