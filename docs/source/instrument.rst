.. module:: xpsi.Instrument

.. _instrument:

Instrument
==========

Instances of :class:`~.Instrument.Instrument` are objects representing a model
instrument.

.. autoclass:: xpsi.Instrument.Instrument
    :members: matrix, energy_edges, cached_signal, construct_matrix
    :special-members: __call__

.. autoclass:: xpsi.Instrument.ResponseError
    :show-inheritance:

.. autoclass:: xpsi.Instrument.EdgesError
    :show-inheritance:

