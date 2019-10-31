.. module:: xpsi

.. _elsewhere:

Elsewhere
=========

Instances of :class:`~.Elsewhere.Elsewhere` represent the region of the
photosphere in which smaller radiative features (such as spots) are embedded.

.. autoclass:: xpsi.Elsewhere.Elsewhere
    :members: __init__, _construct_cellMesh, _compute_rays, _compute_cellParamVecs, eval_srcRadFieldParamVectors, integrate, num_cells, sq_num_cells, num_rays, print_settings, num_params, bounds
    :show-inheritance:

.. autoclass:: xpsi.Elsewhere.RayError
    :show-inheritance:

.. autoclass:: xpsi.Elsewhere.IntegrationError
    :show-inheritance:
