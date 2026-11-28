.. module:: xpsi.Star

.. _star:

Star
====

Instances of :class:`~.Star.Star` are objects representing a model star.
A model star is constituted by:

* an ambient Schwarzschild spacetime solution;
* a discrete representation (on leaves of a spacetime foliation) of an embedded
  stellar 2-surface exterior to which the ambient spacetime solution is assumed
  valid;
* a discrete representation of the local comoving source radiation field in
  time and photon phase-space over the stellar 2-surface;
* a discrete representation of the time-invariant null mapping from the stellar
  2-surface to effective infinity;
* and low-level integrators of the radiation field incident in the
  neighbourhood of a point at effective infinity to obtain pulses.

.. autoclass:: xpsi.Star.Star
    :members: photospheres, update
    :show-inheritance:
