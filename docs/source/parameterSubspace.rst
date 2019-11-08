.. module:: xpsi.ParameterSubspace

.. _parameterSubspace:

ParameterSubspace
=================

Instances of :class:`~.ParameterSubspace.ParameterSubspace` represent *abstract*
subspaces of the model global parameter space. The class itself is an
Abstract Base Class for derived classes.

.. autoclass:: xpsi.ParameterSubspace.ParameterSubspace
    :members: __init__, num_params, bounds

.. autoclass:: xpsi.ParameterSubspace.BoundsError
    :show-inheritance:
