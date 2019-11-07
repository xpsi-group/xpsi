.. _faq:

(FA)Q
-----

.. rubric:: Do I need to edit the package setup script?

You may well have to edit the setup script depending on the target system.

.. rubric:: Does it matter what compiler I use?

The Intel compiler collection has been used successfully for X-PSI and
dependencies (GSL, MultiNest). We recommend first trying to use Intel.

.. rubric:: What Intel instruction sets should I use?

If you want to test the binaries on a login node, note that you may need to
compile with multiple instruction sets.

.. rubric:: Is I/O or disk storage a concern, or are all the files small?

I/O not a concern for likelihood calculation.

Nested sampling writes to disk at user-specified cadence
(so many nested sampling iterations).

Model data such as an four-dimensional atmosphere table can be reasonably
large for I/O.
We recommend loading, at the outset of the run (or a resumed run),
such a table into a contiguous chunk of memory
for each of the Python processes running on one node.
That table is pointed to for access where needed from compiled modules
(C extensions to Python): it is not loaded from disk per likelihood call.
We provide an example custom Python class that handles this loading (as used
in :ref:`R19`, hereafter R19).

Disk storage required is indeed small: up to :math:`\mathcal{O}(100)` Mbytes for
applications thus far (e.g., R19).
There is a variant of MultiNest nested sampling that is much more memory and
disk intensive, but we do not use it.

.. rubric:: Why?

Imporance nested sampling (Feroz et al. 2013) is not compatible with the
alternative options (read: hacks) for prior implementation.
