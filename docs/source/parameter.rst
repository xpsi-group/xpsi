.. module:: xpsi.Parameter

.. _parameter:

Parameter
=========

Instances of :class:`~.Parameter.Parameter` represent *abstract*
model parameters, both free and fixed/derived variables.

.. autoclass:: xpsi.Parameter.Parameter
    :members:
    :special-members: __call__, __len__, __str__, __repr__
    :exclude-members: __init__

    .. automethod:: __init__(name, strict_bounds, bounds=(None,None), doc=None, symbol=r'', value=None, permit_prepend=True)

.. autoclass:: xpsi.Parameter.Derive
    :members:
    :special-members: __call__
