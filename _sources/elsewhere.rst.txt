.. module:: xpsi.Elsewhere

.. _elsewhere:

Elsewhere
=========

Instances of :class:`~.Elsewhere.Elsewhere` represent the region of the
photosphere in which smaller radiative features (such as spots) are embedded.

.. autoclass:: xpsi.Elsewhere.Elsewhere
    :members: _construct_cellMesh, _compute_rays, _compute_cellParamVecs,
              integrate, num_cells, sqrt_num_cells,
              num_rays, print_settings, embed

    :show-inheritance:

.. autoclass:: xpsi.Elsewhere.RayError
    :show-inheritance:

.. autoclass:: xpsi.Elsewhere.IntegrationError
    :show-inheritance:
