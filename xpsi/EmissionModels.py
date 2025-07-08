from xpsi.global_imports import *

from xpsi.ParameterSubspace import ParameterSubspace
from abc import ABCMeta, abstractmethod

class EmissionModel( ParameterSubspace , metaclass=ABCMeta ):

    def __init__(self, *args, **kwargs):
        super( EmissionModel, self ).__init__( *args, **kwargs )

    @abstractmethod
    def integrate(self, st, energies, threads, *args, **kwargs):
        pass

    @abstractmethod
    def embed(self, spacetime, photosphere, fast_total_counts, threads, *args):
        pass


class EmissionModels( ParameterSubspace ):
    pass

    def __init__( self, 
                 emission_models, 
                 *args, **kwargs ):
        
        # Everything must be an EmissionModel
        try:
            for model in emission_models:
                assert isinstance(model, EmissionModel), 'Invalid type for emission model. Must be EmissionModel.'
        except:
            raise TypeError('emission_models must be an interable of EmissionModel.')
        
        # Save it all
        self.emission_models = emission_models
        super( EmissionModels, self ).__init__( *args, **kwargs )

    def integrate(self, st, energies, threads,
                  *args, **kwargs):
        """ Integrate over the emission models (same as for HotRegions)"""
        if not isinstance(energies, tuple):
            energies = [energies] * len(self)

        signals = []
        for model, E in zip(self._emission_models, energies):
            signals.append(model.integrate(st, E,
                                         threads,
                                         *args, **kwargs ))

        return tuple(signals)
    
    def embed(self, spacetime, photosphere, fast_total_counts, threads, *args):
        """ Embed the emission models. """

        if fast_total_counts is None:
            fast_total_counts = [None] * len(self)

        for obj, fast in zip(self._objects, fast_total_counts):
            obj.embed(spacetime,
                      photosphere,
                      fast,
                      threads, *args)
