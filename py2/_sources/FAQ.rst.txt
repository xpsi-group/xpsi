.. _faq:

FAQs and common problems
========================


Installation
^^^^^^^^^^^^

.. rubric:: Do I need to edit the package setup script?

You may well have to edit the setup script depending on the target system.
This includes editing compiler flags (see below for example regarding
instruction sets).

.. rubric:: Does it matter what compiler I use?

The Intel compiler collection has been used successfully for X-PSI and
dependencies (namely GSL, MultiNest). We recommend first trying to use Intel
in a context where performance matters.

.. rubric:: What Intel instruction sets should I use?

If you want to test the binaries on a login node, note that you can
compile with multiple instruction sets for auto-dispatch using the ``-x`` and
``-ax`` flags. See the :ref:`hpcsystems` page for examples.

.. rubric:: What atmosphere extension should I use?

If you want to do some quick calculations and run the Modeling tutorial, you should use the default blackbody atmosphere extension. But if you want to use similar atmosphere models as in the published X-PSI applications so far, you should use to the numerical atmosphere extension when installing X-PSI (see instructions in :ref:`install`).

Note that using the model scripts with an unintended atmosphere extension may lead to a segmentation fault error if trying to use numerical atmosphere without providing numerical atmosphere data. Or if using the scripts including the numerical atmosphere data with the blackbody extension, you can get unexpected results, printing though a warning that numerical atmosphere data were preloaded but not used. Examples of numerical atmosphere data, which are required by the numerical atmosphere extension, can be found here: `doi:10.5281/zenodo.7094144`__. Examples of how to use the numerical atmospheres are shown e.g. in Surface radiation field tools -tutorial and in :ref:`example_script`.

.. _Zenodo: https://doi.org/10.5281/zenodo.7094144

__ Zenodo_


Sampling
^^^^^^^^

.. rubric:: Is I/O or disk storage a concern, or are all the files small?

I/O not a concern for likelihood calculation.

Nested sampling writes to disk at user-specified cadence
(so many nested sampling iterations).

Model data such as a four-dimensional atmosphere table can be reasonably
large for I/O.
We recommend loading, at the outset of the run (or a resumed run),
such a table into a contiguous chunk of memory
for each of the Python processes running on one node.
That table is pointed to for access where needed from compiled modules
(C extensions to Python): it is not loaded from disk per likelihood call.
We provide an example custom Python class that handles this loading (as used
in `Riley et al. 2019 <https://ui.adsabs.harvard.edu/abs/2019ApJ...887L..21R/abstract>`_, hereafter R19).

Disk storage required is indeed small: up to :math:`\mathcal{O}(100)` Mbytes for
applications thus far (e.g., R19). There is a variant of MultiNest nested sampling
that is much more memory and disk intensive, but we do not use it.  This is
because importance nested sampling is not compatible with the alternative options
(read: hacks) for prior implementation (see `Riley, PhD thesis <https://hdl.handle.net/11245.1/aa86fcf3-2437-4bc2-810e-cf9f30a98f7a>`_).


Common problems and errors
^^^^^^^^^^^^^^^^^^^^^^^^^^

.. rubric:: How to avoid errors in post-processing?

Do not use X-PSI PostProcessing tools for runs which have not converged yet or have not enough samples. Also, when post-processing, make sure to check the data and output file paths, use ``cache=True`` if plotting the signal, and perform a likelihood check to be sure that the imported model is the same as in the run.

.. rubric:: *AttributeError: ’NestedBackend’ object has no attribute ’\ :math:`\_nc\_bcknd`\ ’*

Can happen in PostProcessing for runs with ``use_nestcheck=[False]`` (e.g. importance sampling). Solution is to turn ``bootstrap_estimators=False``, or alternatively, set ``use_nestcheck=[True]``.

.. rubric:: Why does my skymap show many annular images like this:

.. container:: figure*

   .. image:: _static/ST_PST__NICER__skymap_phase_averaged_run1.png
      :alt: image
      :width: 10cm

The problem is the ``xpsi/xpsi/surface_radiation_field/local_variable.pyx`` file which should be overwritten by ``xpsi/xpsi/surface_radiation_field/archive/local_variables/PST_U.pyx`` or ``xpsi/xpsi/surface_radiation_field/archive/local_variables/two_spots.pyx`` (depending on the model) and then re-install X-PSI.

.. rubric:: *ImportError: No module named tools*

You are running X-PSI from its main directory (the directory where the ``setup.py`` file is). Exit that directory and run it again.

.. rubric:: *<path/to/run/output>dead-birth.txt not found.*

Set ``use_nestcheck=[False]`` or check that nestcheck is installed exactly as instructed in :ref:`install` (by cloning it from ``https://github.com/ThomasEdwardRiley/nestcheck.git``).

.. rubric:: *Invalid caching targets.*

Set ``cache=True`` for the signal.

.. rubric:: *Each row and column must contain at least one positive number.*

There are some rows and/or column in the instrument response that contain only zeros. Solution is to increase the number of channels or decrease the number of energy intervals.

.. rubric:: *Warning: Using native nestcheck KDE instead of GetDist KDE.*

Make sure to to install nestcheck and GetDist packages using the corresponding github repositories as instructed in :ref:`install`.

.. rubric:: *ValueError: There is more than one signal instance.*

Typically occurs when post-processing joint NICER and XMM results, if not setting ``model.likelihood.signals = model.likelihood.signals[0][0]`` (when plotting the inferred NICER signal).
