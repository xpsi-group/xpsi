.. module:: xpsi.Everywhere

.. _everywhere:

Everywhere
==========

Instances of :class:`~.Everywhere.Everywhere` represents the global photosphere.

.. autoclass:: xpsi.Everywhere.Everywhere
    :members: _construct_cellMesh, _compute_rays, _compute_cellParamVecs,
              integrate, num_cells, sqrt_num_cells,
              num_rays, print_settings, set_phases, phases_in_cycles, embed
    :show-inheritance:

.. autoclass:: xpsi.Everywhere.RayError
    :show-inheritance:

.. autoclass:: xpsi.Everywhere.IntegrationError
    :show-inheritance:
