from xpsi.global_imports import *

from xpsi.ParameterSubspace import ParameterSubspace
from abc import ABCMeta, abstractmethod

class EmissionModel( ParameterSubspace , metaclass=ABCMeta ):

    def __init__(self, *args, **kwargs):
        super( EmissionModel, self ).__init__( *args, **kwargs )

    # @abstractmethod
    # def integrate(self, **kwargs):
    #     pass

    # @abstractmethod
    def embed(self, spacetime, threads):
        pass


class EmissionModels( ParameterSubspace ):

    def __init__( self, 
                 components, 
                 *args, **kwargs ):
        
        # Everything must be an EmissionModel
        try:
            if isinstance( components, EmissionModel ):
                components = [components]
            for model in components:
                assert isinstance(model, EmissionModel), 'Invalid type for emission model. Must be EmissionModel.'
        except:
            raise TypeError('emission_models must be an interable of EmissionModel.')
        
        # Save it all
        self.components = components
        super( EmissionModels, self ).__init__( components, *args, **kwargs )

    def integrate(self, energies, threads,
                  *args, **kwargs):
        """ Integrate over the emission models (same as for HotRegions)"""

        signals = []
        for model in self.components:
            signals.append(model.integrate( energies,
                                         threads,
                                         *args, **kwargs ))

        return tuple(signals)
    
    def update(self, threads=1, force_update=False):
        """ Embed the emission models into the ambient spacetime. """

        # TODO link this to the spacetime
        for model in self.components:
            if model.needs_update or force_update:
                model.embed( 1, threads )