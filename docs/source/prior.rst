.. module:: xpsi.Prior

.. _prior:

Prior
=====

Instances of :class:`Prior` are called by worker sampling processes for
evaluation of a joint prior distribution.

.. autoclass:: xpsi.Prior.Prior
    :members: __init__, __call__, inverse_sample, draw,
              estimate_hypercube_frac, unit_hypercube_frac
