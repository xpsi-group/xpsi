.. _faq:

(FA)Q
=====

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

If you want to do some quick calculations and run the Modeling tutorial, you should use the default blackbody atmosphere extension. But if you want to use similar atmosphere models as in the published X-PSI applications so far, you should switch to the numerical atmosphere extension before installing X-PSI (see instructions in :ref:`install`).

Note that using the model scripts with an unintended atmosphere extension may lead to a segmentation fault error if trying to use numerical atmosphere without providing numerical atmosphere data. Or if using the scripts including the numerical atmosphere data with the blackbody extension, you can get unexpected results, printing though a warning that numerical atmosphere data were preloaded but not used. Examples of numerical atmosphere data, which are required by the numerical atmosphere extension, can be found e.g. in the Zenodo repository of Riley et al. 2021: `doi:10.5281/zenodo.4697625`__. Examples of how to use the numerical atmospheres are shown e.g. in Surface radiation field tools -tutorial and in :ref:`example_script`.

.. _Zenodo: https://zenodo.org/record/4697625

__ Zenodo_


Model setup
^^^^^^^^^^^

Future questions and answers will be archived here.

Batch usage
^^^^^^^^^^^

Future questions and answers will be archived here.

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
in :ref:`R19`, hereafter R19).

Disk storage required is indeed small: up to :math:`\mathcal{O}(100)` Mbytes for
applications thus far (e.g., R19). There is a variant of MultiNest nested sampling
that is much more memory and disk intensive, but we do not use it.  This is
because importance nested sampling is not compatible with the alternative options
(read: hacks) for prior implementation (see `Riley, PhD thesis <https://hdl.handle.net/11245.1/aa86fcf3-2437-4bc2-810e-cf9f30a98f7a>`_).
