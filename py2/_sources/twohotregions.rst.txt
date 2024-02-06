.. module:: xpsi.TwoHotRegions

.. _twohotregions:

TwoHotRegions
=============

Instances of :class:`~.TwoHotRegions.TwoHotRegions` are represent pair of
radiatively intense regions of the source photosphere. The class inherits
from :class:`~.HotRegion.HotRegion`, and applies the methods to compute photon
pulses from a pair of spots where one is derived from the other.
Specifically, the regions are related via antipodal reflection symmetry.
Useful for pulsars under the assumption of a dominantly dipolar magnetic
field with two magnetic polar caps.

The same configuration can be achieved manually via application of the
:class:`~.HotRegions.HotRegions` class, with very slightly higher likelihood
expense (the difference is not significant enough to concern oneself with).

.. autoclass:: xpsi.TwoHotRegions.TwoHotRegions
    :members: embed, integrate, print_settings, cellArea
    :show-inheritance:

.. autoclass:: xpsi.TwoHotRegions.PulseError
    :show-inheritance:
