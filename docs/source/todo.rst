.. _TODO:

Future
------

Below we make (non-exhaustively) note down ideas for future development.
Some features may be added in a backwards compatible manner and thus in
a minor release, whilst others will require inclusion in a major release.
Where we have a concrete expectation, we note which type of release each
listed feature should be included in.

Priority
^^^^^^^^

* Port to a Python 2/3 compatible state.
* Unit testing. At present we are relying on the tutorial
  notebooks and examples to flag problems. *Target: minor release v1.x*.

Tentative
^^^^^^^^^

* Implement a simpler switch between atmosphere model extensions (e.g.,
  blackbody to numerical lookup table), rather than user having to remember to
  modify the relevant ``.pyx`` source file (e.g., by replacing function bodies
  and custom structs with code from the ``surface_radiation_field/archive``)
  and then recompile. Perhaps look into C function pointers passed to Cython for
  runtime specification of shared object.
* Extension to interpolate in arbitrary number of dimensions (currently hard-
  coded four-dimensional cubic polynomial interpolation for, meaning two
  variables in addition to energy and zenith angle).




