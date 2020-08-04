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

The alphabetic version tags below such as v1.a and v1.b give a loose indication
of priorty within a major release cycle, and are not generaly unique. That is,
v1.a might ultimately not contain a feature next to which it is listed, and/or
it might be identical to v1.b.

* Compute secondary images, and more generally images up to some optional order
  using surface discretisation. *Target: minor release v0.a*.
* Port to a Python 2/3 compatible state. *Target: minor release v1.a or major \
  release v2.0*.
* Unit testing. At present we are relying on the tutorial
  notebooks and examples as test beds to flag problems.
  *Target: minor release v1.b*.
* Post-processing option to compute highest-density credible intervals
  (appropriate for multi-modal marginal posterior). *Target: minor release v1.c*.

Prospective
^^^^^^^^^^^

Unordered.

* Make a dedicated directoty for the tutorial notebooks in the documentation
  pages, link to them on GitHub, and state what is necessary to run those
  notebooks completely.
* Implement a simpler switch between atmosphere model extensions (e.g.,
  blackbody to numerical lookup table), rather than user having to remember to
  modify the relevant ``.pyx`` source file (e.g., by replacing function bodies
  and custom structs with code from the ``surface_radiation_field/archive``)
  and then recompile. Perhaps look into C function pointers passed to Cython for
  runtime specification of shared object.
* Extension to interpolate in arbitrary number of dimensions (currently hard-
  coded four-dimensional cubic polynomial interpolation for, meaning two
  variables in addition to energy and zenith angle).
* Add a plot class to render the posterior instrument effective area curves,
  optionally as a set of curves or as conditional posterior bands; relevant
  only for parameterised instrument models.
* Signal plotting tools not associated with post-processing. E.g., simple
  functions to plot single signals (instead of many signals, each associated
  with a posterior sample) cached when the likelihood function is evaluated,
  as demonstrated in the various tutorial notebooks (thus some such functions
  already prototyped). Which module(s) to add these to?
* Support for sensitivity analysis via importance sampling when post-processing
  posterior samples.
* Support for specifying which subset of energies is used for calculating
  signals from which surface components.
* Module containing a class for arbitrary likelihood factor that might be a
  function of parameters defined in X-PSI, such as the mass and distance. It
  would plug into an instance of the existing likelihood class. Currently,
  support for this is provided by the `Prior` class, so it is to be decided
  if a distinct class for arbitrary likelihood factors is of any further use.
* Develop additional extensions for the archive that transform global variables
  and spacetime coordinates into local variables to evaluate the local specific
  intensity emergent from the photosphere along a ray. These archived
  extensions need to implement two overlapping circular regions constituting a
  contiguous surface hot region, leading to morphologies such as
  single-temperature rings and crescents, and two-temperature hot regions that
  are implemented for signal integration via surface discretisation and thus
  for likelihood evaluations. It is very useful to visualise the self-lensed
  images of the star (and the star-receiver configuration) resolved in sky
  direction (i.e., specific intensity, and intensity sky maps, phase resolved
  and phase averaged). The existing extensions can handle two
  single-temperature circular spots, and other complexities, but do not
  precisely implement the aforementioned hot region models. It might be
  possible to develop one complex extension that is all-encompassing in this
  respect, but it would require a user to learn how to deactivate certain
  complexities in order to match the models implemented for surface
  discretisation.
* Develop animated photon incident flux pulse-profile and phase-averaged photon
  incident specific flux spectrum plots, optionally with hot region components
  shown independently, to optionally display alongside sky maps.
* Develop helper functions to wrap the atmosphere checking tools and generate
  standard plot types so a user does not need to handle matplotlib objects as
  much. The relevant tutorial notebook shows how it can be done for some simple
  plots. It is also to be determined what types of plots are of most interest
  and useful, or are considered standard.
* Handle joint modeling of X-ray and Far-UV data with multiple instruments. The
  Far-UV instrument operation and likelihood function form needs to be worked
  out and implemented, but phase-averaged signals, the notion of ``Elsewhere``,
  and different atmospheres (e.g., ionized hot regions + partially-ionized
  elsewhere) currently supported.
* Add customisable method somewhere to transform raw sample file vectors to
  parameter values compatible with current API, in case of backwards
  incompatible changes, to avoid user having to figure out how to do this
  safely.
* Check to see where the extension tools can be called from other extensions
  (e.g., likelihood extensions calling phase integrator tool) where appropriate
  instead of having similar operations coded twice, for maintainability.
* Decide whether to allow user to define their photon energy and effective area
  units (requires a minor tweak and explanation in docstrings/comments).
