.. module:: xpsi.Likelihood

.. _likelihood:

Likelihood
==========

Instances of :class:`~.Likelihood.Likelihood` are called by worker sampling
processes for evaluation of a likelihood function.

.. autoclass:: xpsi.Likelihood.Likelihood
    :members: __call__, __str__, threads,
              synthesise, check, llzero, random_near_llzero, less_than_llzero,
              star, prior, signals, prior
    :show-inheritance:
