.. module:: xpsi.Prior

.. _prior:

Prior
=====

Instances of :class:`Prior` are called by worker sampling processes for
evaluation of a joint prior distribution.

.. autoclass:: xpsi.Prior.Prior
    :members: __call__, inverse_sample, draw, inverse_sample_and_transform,
              transform, estimate_hypercube_frac, unit_hypercube_frac
