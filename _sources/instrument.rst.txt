.. module:: xpsi.Instrument

.. _instrument:

Instrument
==========

Instances of :class:`~.Instrument.Instrument` are objects representing a model
instrument.

.. autoclass:: xpsi.Instrument.Instrument
    :members: matrix, energy_edges, channels, cached_signal, construct_matrix,
              channel_edges
    :special-members: __call__
    :show-inheritance:

.. autoclass:: xpsi.Instrument.ResponseError
    :show-inheritance:

.. autoclass:: xpsi.Instrument.EdgesError
    :show-inheritance:

.. autoclass:: xpsi.Instrument.ChannelError
    :show-inheritance:
