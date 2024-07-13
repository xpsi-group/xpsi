.. module:: xpsi.Signal

.. _signal:

Signal
======

Instances of :class:`~.Signal.Signal` are objects which define the
generative relationship between the pulse and X-ray photon data.
That is, the parametrised sampling distribution of the photon data is defined.

.. autoclass:: xpsi.Signal.Signal
    :members:
    :show-inheritance:
    :special-members: __call__

.. autoclass:: xpsi.Signal.LikelihoodError
    :show-inheritance:

.. autofunction:: xpsi.Signal.construct_energy_array
