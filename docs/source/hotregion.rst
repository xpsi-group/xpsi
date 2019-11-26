.. module:: xpsi.HotRegion

.. _hotregion:

HotRegion
=========

Instances of :class:`~.HotRegion.HotRegion` are objects representing radiatively
intense regions of the source photosphere.

.. autoclass:: xpsi.HotRegion.HotRegion
    :members: _HotRegion__construct_cellMesh, _HotRegion__compute_rays,
              _HotRegion__compute_cellParamVecs, embed,
              _psi, integrate, num_rays, sqrt_num_cells, leaves, phases,
              set_phases, phases_in_cycles, print_settings, num_cells,
              num_params, bounds, cede, concentric
    :show-inheritance:

.. autoclass:: xpsi.HotRegion.RayError
    :show-inheritance:

.. autoclass:: xpsi.HotRegion.PulseError
    :show-inheritance:
