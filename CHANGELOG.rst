History
-------

All notable changes to this project will be documented in this file.

The format is based on
`Keep a Changelog <http://keepachangelog.com/en/1.0.0/>`_
and this project adheres to
`Semantic Versioning <http://semver.org/spec/v2.0.0.html>`_.


[Unreleased]
~~~~~~~~~~~~

Summary
^^^^^^^

Fixed
^^^^^

Added
^^^^^

Changed
^^^^^^^

Deprecated
^^^^^^^^^^

Removed
^^^^^^^

Attribution
^^^^^^^^^^^


[v0.5.1] - 2020-08-07
~~~~~~~~~~~~~~~~~~~~~

Fixed
^^^^^

* Bug when plotting intensity sky maps because a line was inadvertently
  removed.
* Some mutable defaults in :class:`xpsi.Elsewhere` and :class:`xpsi.Everywhere`.
* Conditional statement in :meth:`xpsi.Photosphere.embed`.

Added
^^^^^

* Capability to add custom parameters when instantiating
  :class:`xpsi.Photosphere`, which is useful for calling image plane extensions
  whilst passing global variables, without having to instantiate
  surface-discretisation classes and without having to handle global variable
  values at compile time or from disk for runtime access.


[v0.5.0] - 2020-08-06
~~~~~~~~~~~~~~~~~~~~~

Summary
^^^^^^^

* The major change is an update and refactoring of the post-processing module
  to work again with past API changes. (The module was not being kept up to date
  with previous releases listed below because it wasn't to our knowledge
  being used by anyone yet, and thus we focussed on other features.) The module
  has been refactored to be more modular, flexible, and extensible. For
  instance, posterior signal-plot classes can be added by the user and
  complex plotting routines can thus be developed, as demonstrated in the
  concrete classes such as :class:`xpsi.PostProcessing.PulsePlot`. The plot
  classes have been used to reproduce (with improved functionality and
  performance) the relevant signal plots from :ref:`R19`, as demonstrated
  in the post-processing tutorial notebook and embedded in the class docstrings
  for reference.
* Development of online documentation pages, including project organisation
  pages and a Code of Conduct (please read), and development of docstrings.
  Note that some snippets of documentation look forward to v1.0 (e.g., release
  of technical notes in the repo itself).

Fixed
^^^^^

* The :class:`xpsi.Data` docstring explanations have been improved for clarity,
  mainly regarding the instrument channel definitions. The explanation is of how
  the information contained in a :class:`xpsi.Data` instance pertains to the
  *loaded* instrument response (sub)matrix.
* The :class:`xpsi.Instrument` docstrings have also been improved for clarity,
  explaining the relationship to :class:`xpsi.Data` in more detail.
* Update extension module for background marginalisation to take distinct phase
  sets associated with hot regions.
* The constructor :meth:`xpsi.Spacetime.fixed_spin` inclination upper bound
  is :math:`\pi/2` radians to eliminate degeneracy due to equatorial-reflection
  symmetry in the default prior on source-receiver geometric configuration.
* Tweak caching (memoization) so that cache and current vectors are equal at
  the end of likelihood evaluation routine.
* Generally clean up naming and docstrings for extension modules. Add return
  types.
* Bug was fixed for transforming posterior sample sets and prior samples when
  parameter orders different in sample files and a prior object due to API
  updates. Whether this solution is to be long-term is to be decided; more
  generally need to figure out how to elegantly handle derived parameters that
  are not needed for likelihood evaluation (those derived parameters are
  instances of :class:`xpsi.Parameter`) but are of interest for post-processing.
* Handle ``param_plot_lims=None`` correctly in
  :class:`xpsi.PostProcessing.CornerPlotter`.
* Checked for unintended mutable defaults package-wide, and fixed as
  appropriate.
* Fix bugs in ``CustomPrior`` class (:ref:`example_script`; these example
  modules were not run at the time of translation between past API versions, so
  only found bugs when making post-processing tutorial for this release).
* The formatting of annotated credible intervals in
  :class:`xpsi.PostProcessing.CornerPlotter` has been improved by inferring the
  largest number of decimal places needed for two non-zero decimal digits, and
  then formatting the median and quantile differences to this shared decimal
  precision above the on-diagonal panels. If the numbers cannot be well-
  represented by this scheme, the user could try a unit transformation.
* Tried to tweak automated margins for intensity sky map multi-panel plots,
  so as not to sometimes partially cut an axis label.
* Bug that prevented animation of sky map frames written to disk because the
  frames were not cached in memory by reimaging.

Added
^^^^^

* The :class:`xpsi.Data` is now concrete in implementation, such that in common
  usage patterns, it does not need to be subclassed.
* A constructor to :class:`xpsi.Data` to load a phase-folded event list and
  phase-bin the events in a subset of selected channels.
* A :meth:`xpsi.Data.channels` property that holds the instrument channels
  to be checked by a :class:`xpsi.Signal` instance against those declared for
  the loaded instrument response (sub)matrix. This property as also required by
  the post-processing module (namely, :class:`xpsi.PostProcessing.ResidualPlot`
  and the other :class:`xpsi._signalplot.SignalPlot` subclasses).
* A :meth:`xpsi.Instrument.channels` property that holds the instrument
  channels to be checked by a :class:`xpsi.Signal` instance against those
  declared for the event data matrix.
* Support for multiple instruments operating on the same incident signal due to
  assumed effective time-invariance of the signal generated during one
  rotational cycle of the surface radiation field.
* Module :mod:`xpsi.surface_radiation_field` to call atmosphere extensions
  directly (without the calls being embedded in integration algorithms), for
  checking implementation of complicated atmospheres such as those requiring
  interpolation with respect to a numerical lookup table.
* Support for the extension module for calculating the local surface radiation
  field variables to read in numerical model data. An example extension module
  designed to execute nearest-neighbour lookup amonst an general unstructured
  array of points of the openness of magnetic field lines has been developed.
* Add simple energy annotation option to photon specific intensity sky-map
  panels.
* State the energy units (keV) that the :class:`xpsi.Instrument` must comply
  with when energy interval bounds are specified.
* State the units of variables such as energy and specific intensity in the
  surface radiation field extension module. These requirements may be found in
  function body comments.
* Explain in :class:`xpsi.PostProcessing.CornerPlotter` docstring the order in
  which posteriors are plotted given the input order.
* Post-processing switches to overwrite transformed-sample files and
  combined-run files on disk.
* Workaround to handle the case where due to API changes, the relationship
  between sample parameter vectors on disk and the parameter vector in the
  current API are related not just by reordering, but transformations. This
  is demonstrated in the post-processing tutorial instead of transforming the
  original sample files on disk in place, the transformed files written to disk
  contain both the transformed vector (same number of elements) to match the
  parameters defined under the current API (the order of the vector can be
  different between the :class:`xpsi.ParameterSubspace` underlying with a
  :class:`xpsi.Likelihood` instance and the files on disk containing the
  transformed samples), and the additional derived parameters.
* Attempt to free up memory when :meth:`xpsi.Photosphere.images` is no longer
  needed, but memory-intensive operations need to be performed.
* Attempt to free memory properly after animating a sky-map phase sequence.

Changed
^^^^^^^

* Change (Earth) inclination parameter :math:`i` to :math:`\cos(i)` so that the
  default prior density function is isotropic.
* The object formerly named ``xpsi.Pulse`` has had its name changed to
  :class:`xpsi.Signal`, and across the package, names that were ``pulse`` are
  apart from potential corner cases or documentation instances of the word,
  are now ``signal``, because when support joint likelihood functions over
  multiple instruments, some data sets are phase averaged. Moreover, *signal*
  is arguably clearer in meaning than *pulse*, once it has been established
  that the signals the package focuses on are *pulsed* but depending on
  the instrument, the data we confront the model with has some degree of phase
  (timing) resolution that might be insufficient for phase-resolved
  observations.
* The :class:`xpsi.Data` definition of the ``last`` channel has changed to be
  the index of the last row in the loaded instrument response (sub)matrix,
  instead of being the index of the last row plus one; this means that the
  value exposed via a property is ``last+1``.
* For numerical atmospheres of same number of grid dimensions, improved
  extension ``surface_radiation_field/archive/{hot,elsewhere}/numerical.pyx``
  module to infer grid size for memory allocation and interpolation searches
  (implemented automatic inference of grid size, but hard-coded
  four-dimensional cubic polynomial interpolation persistent). Different
  those atmospheres can be loaded simply via a Python subclass without
  the relevant extension module being recompiled.
* The :class:`xpsi.Photosphere` class sometimes does no surface discretisation,
  so allow no hot regions, elsewhere, or everywhere objects; then image-plane
  discretisation can be accessed without dummy object creation.
* Tweak :class:`xpsi.SpectrumPlot` settings to print a warning statement that
  spectrum plot works best with logarithmic spacing, and the user has to shadow
  class attribute with ``logspace_y=False``.
* Do not print :class:`xpsi.HotRegion` instance parameter properties upon
  creation if fixed at boundary value so that the region is fully described by
  fewer parameters.
* Merged energy integration extension modules into one.
* Made phase shift parameters (strictly) unbounded; remember however that for a
  sensible prior, bound the phase shifts on a unit interval, and thus it is
  required that phase bounds are specified and finite.
* In extensions, modified phase shifting so that a shift permitted by unbounded
  phase parameter does not require many iterations to decrement or increment to
  unit interval (achieved simply with floor operation).

Deprecated
^^^^^^^^^^

* The :meth:`xpsi.Data.channel_range` property has been renamed to
  :meth:`xpsi.Data.index_range` so as to avoid confusion between these numbers
  and the true instrument channels. *The old property will be removed for
  release v1.0*.

Removed
^^^^^^^

* The ensemble MCMC sample backend for post-processing because we do not expect
  it to be useful in the immediate future, but requires some non-trivial
  development work to meld properly with the current post-processing module
  which is focussed on nested sampling. This functionality will be reintroduced
  in a future release (refer to :ref:`todo`). The ensemble sampler can still be
  run, however, and the native backend for accessing sample information on disk
  is demonstrated in a tutorial notebook. However, the runs cannot be processed
  for posterior integrals and visualisation using the same tools as available
  for nested sampling runs.

Attribution
^^^^^^^^^^^

* With thanks to Sebastien Guillot (testing and feedback),
  Devarshi Choudhury (testing and feedback),
  Sam Geen & Bob de Witte (Windows installation advice),
  and Anna L. Watts (documentation patches and feedback).


[v0.4.1] - 2020-06-03
~~~~~~~~~~~~~~~~~~~~~

Fixed
^^^^^

* Function signatures to match header declarations in atmosphere extensions:
  ``xpsi/surface_radiation_field/archive/elsewhere/numerical.pyx`` to match
  ``xpsi/surface_radiation_field/elsewhere_radiation_field.pxd``.
  With thanks to Sebastien Guillot.


[v0.4.0] - 2020-02-14
~~~~~~~~~~~~~~~~~~~~~

Summary
^^^^^^^

* Mainly new features.
* Backwards compatible (apart from possible corner cases).

Fixed
^^^^^

* Removed a spurious geometric factor in the integrator that discretises the
  surface with a static mesh. This integrator was called by the ``Elsewhere``
  class. The error when this factor is included is O(1%) at 600 Hz for soft
  emission from the entire stellar disk, and then scales with spin and energy
  beyond this. To reproduce the bug, find the commented out ``/ superlum`` in
  file ``xpsi/cellmesh/integrator_for_time_invariance.pyx`` (line 251) and
  uncomment it. Then reinstall the package. When this factor is included, the
  mesh itself is moving in the context of the images subtended by its
  constituent elements on our sky. We want the mesh to be static so that this
  integrator can be used for faster calculation of time-invariant signals.
* Bug in which the prior density factor is incorporated twice if a ``Likelihood``
  instance held a reference to a ``Prior`` object and these are merged into
  a ``Posterior`` object which is fed to the ensemble sampler. If the prior
  density was *flat*, this bug will have had no effect on posterior
  distributions.

Added
^^^^^

* New features are the simulation of signals from more general surface
  radiation fields that globally span the stellar surface. This can be
  done with several types of integrator.
* The new image-plane discretisation integrator supports imaging of a star,
  and Python functionality has been added to automate plotting and animation
  of intensity sky maps.
* A new tutorial to the documentation to demonstrate these new features and
  an internal cross-check of distinct integration algorithms.
* A visual introduction to the documentation pages with some animated sky maps.


[v0.3.6] - 2020-01-24
~~~~~~~~~~~~~~~~~~~~~

Fixed
^^^^^

* Some code snippets in documentation examples of prior implementation
  with the latest API minor version (v0.3).

Changed
^^^^^^^

* Modify the ``HotRegions`` class to function with two *or more* hot region
  objects.


[v0.3.5] - 2020-01-22
~~~~~~~~~~~~~~~~~~~~~

Summary
^^^^^^^

* Docstring edits and backwards compatible changes to several class
  initialisation arguments.

Attribution
^^^^^^^^^^^

* Based mostly on discussion with and feedback from Devarshi Choudhury.

Fixed
^^^^^

* Some docs formatting problems.
* Some corrections to example scripts/modules updated in v0.3.4 to use
  current API.

Changed
^^^^^^^

* The photospheric mode frequency parameter is not converted to an angular
  frequency until it is used, so the cached value matches the docstring
  description.

Deprecated
^^^^^^^^^^

* The ``is_secondary`` argument of the ``HotRegion`` class. Use ``is_antiphased`` instead
  to ensure future compatibility.
* The ``store`` argument of the ``Pulse`` class. Use ``cache`` instead to ensure future
  compatibility.


[v0.3.4] - 2020-01-20
~~~~~~~~~~~~~~~~~~~~~

Summary
^^^^^^^

* A few patches including backwards compatible improvements.
* Various docstring/comment/doc edits.
* Update docs example model to use v0.3.4 API.

Fixed
^^^^^

* Ensure consistency between input parameter ``bounds`` and ``values`` by
  always requiring dictionaries. Fix applies to ``Elsewhere`` and 
  ``Photosphere``. Courtesy Sebastien Guillot.
* Gravitational mass doc typo fix.

Changed
^^^^^^^

* Add input argument checks to ``Likelihood.check`` method.
* Add default ``hypercube=None`` to ``Prior.inverse_sample_and_transform``
  method.
* If derived parameters found in subspace, assume an update is needed because
  cache mechanism not in place. (WIP.)


[v0.3.3] - 2020-01-20
~~~~~~~~~~~~~~~~~~~~~

Fixed
^^^^^

* At several places in the ``Likelihood`` class, calls were place to ``self``,
  forgetting that ``Likelihood`` overwrites ``ParameterSubspace.__call__``.
  Now calls are ``super(Likelihood, self).__call__()`` to obtain the current
  parameter vector.

[v0.3.2] - 2020-01-16
~~~~~~~~~~~~~~~~~~~~~

Summary
^^^^^^^

* Bug fixes. Backwards compatible.
* When initializing the ensemble-MCMC chains using an nd-ball, the inclusion
  in the prior support was checked by passing a vector to ``Prior.__call__`` but
  that code assumed that the parameter vector had already been assigned and
  can be accessed through the ``ParameterSubspace``. As a result either an
  exception would be thrown (if parameter objects have no value set) or the
  support condition would be evaluated for some preset vector that does not
  change has we iterate through chains.
* The ``Likelihood.check`` method now has a fallback implementation given that
  the NumPy ``allclose`` function in v1.17 does not support Python 2.7.

Attribution
^^^^^^^^^^^

* Based on testing by Sebastien Guillot.

Fixed
^^^^^

* The ``EnsembleSampler`` so that it does not rely on the ``CustomPrior.__call__``
  implementation to handle a vector argument. Chains should now be in
  prior support from the start and never leave.
* The ``Likelihood.check`` method so that a call to a ``Likelihood`` instance
  updates the parameters with a vector if the physical points are passed
  for value checking.
* The ``Likelihood.check`` method error error handling and if/else branching
  has been fixed.
* Some typographic errors in docs.

Changed
^^^^^^^

* The way ``EnsembleSampler`` accesses the prior object.


[v0.3.1] - 2019-12-12
~~~~~~~~~~~~~~~~~~~~~

Fixed
^^^^^

* Some docstring and Sphinx-related formatting.


[v0.3.0] - 2019-12-10
~~~~~~~~~~~~~~~~~~~~~

Summary
^^^^^^^

* Not backwards compatible.
* The main feature is a more sophisticated backend for handling parameters,
  parameter subspaces, and the object hierarchy that forms the modelling
  language. Notably, the parameter objects can be accessed everywhere more
  readily, with dictionary-like functionality that alleviates the problem
  of remembering the imposed order of parameters in a vector. Resultantly,
  there is much more freedom when a user constructs a model and interfaces
  it with sampling software.
* Model parameters can either be *free*, *fixed/frozen* at some scalar value,
  or *derived* deterministically from other model parameters.
* The docs and tutorials have also been updated to reflect these developments.

Attribution
^^^^^^^^^^^

* Feedback and ideas for the above development were discussed at an X-PSI
  workshop in Amsterdam, November 25-29 2019:
  Sebastien Guillot, Emma van der Wateren, Devarshi Choudhury, Pushpita Das,
  Anna Bilous, and Anna Watts.

Added
^^^^^

* A new class ``xpsi.Parameter`` of which every model parameter is an instance.

Changed
^^^^^^^

* The ``xpsi.ParameterSubspace`` class, which has far more sophisticated behaviours
  as a parameter container. The class, upon initialisation with arguments,
  also merges parameters and subspaces into a higher-dimensional (sub)space.
  Most other classes in the modelling language *inherit* from the
  ``xpsi.ParameterSubspace`` class.
* The ``xpsi.TwoHotRegions`` class is now dedicated to representing antipodally
  reflection-symmetric configurations only to simplify the choice of which
  class to use between ``xpsi.HotRegions`` and ``xpsi.TwoHotRegions``. However,
  antipodally reflection-symmetric models can also be constructed using
  just ``xpsi.HotRegions`` because of the new *derived* parameter support. The
  may be a minor speed difference: ``xpsi.TwoHotRegions``
  should be very slightly faster, but it might be imperceptible. Future
  warning: in the future ``xpsi.TwoHotRegions`` might removed altogther for
  simplication.
* The ``xpsi.Photosphere`` class can be instantiated to encapsulate only a
  reference to an ``xpsi.Elsewhere`` instance, and no ``xpsi.HotRegion`` instances.
  An ``xpsi.Elsewhere`` instance can by definition only generate a
  phase-invariant signal. However, further development is needed to handle
  this phase-invariant signal efficiently for likelihood functionality,
  given that operations with respect to phase are not required. Instead
  likelihood functions would be defined only with respect to energy.

Removed
^^^^^^^

* The ``xpsi.ParameterSpace`` module. The global model parameter space is also
  simply an intance of the ``xpsi.ParameterSubspace`` class.
