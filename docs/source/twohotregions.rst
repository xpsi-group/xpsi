.. module:: xpsi

.. _twohotregions:

TwoHotRegions
=============

Instances of :class:`~.TwoHotRegions.TwoHotRegions` are represent pair of
radiatively intense regions of the source photosphere. The class inherits
from :class:`~.HotRegion.HotRegion`, and applies the methods to compute photon
pulses from a pair of spots which may be related under some model; useful for
pulsars under the assumption of a dominantly dipolar magnetic field with two
magnetic polar caps.

.. autoclass:: xpsi.TwoHotRegions.TwoHotRegions
    :members: __init__, embed, integrate, print_settings, num_params, bounds, cellArea
    :show-inheritance:

.. autoclass:: xpsi.TwoHotRegions.PulseError
    :show-inheritance:
