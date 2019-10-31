.. module:: xpsi

.. _likelihood:

Likelihood
==========

Instances of :class:`~.Likelihood.Likelihood` are called by worker sampling
processes for evaluation of a likelihood function.

.. autoclass:: xpsi.Likelihood.Likelihood
    :members: __init__, __call__, __str__, threads, num_params, theta
    :show-inheritance:

.. autoclass:: xpsi.Likelihood.TagError
    :show-inheritance:
