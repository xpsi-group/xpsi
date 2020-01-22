# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](http://keepachangelog.com/en/1.0.0/)
and this project adheres to
[Semantic Versioning](http://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Fixed

### Added

### Changed

### Deprecated

### Removed


## [v0.3.5] - 2020-01-22

### Summary

* Docstring edits and backwards compatible changes to several class
    initialisation arguments.

### Attribution

* Based mostly on discussion with and feedback from Devarshi Choudhury.

### Fixed

* Some docs formatting problems.

* Some corrections to example scripts/modules updated in v0.3.4 to use
    current API.

### Changed

* The photospheric mode frequency parameter is not converted to an angular
    frequency until it is used, so the cached value matches the docstring
    description.

### Deprecated

* `is_secondary` argument of the `HotRegion` class. Use `is_antiphased` instead
    to ensure future compatibility.

* `store` argument of the `Pulse` class. Use `cache` instead to ensure future
    compatibility.


## [v0.3.4] - 2020-01-20

### Summary

* A few hotfixes and backwards compatible improvements.

* Various docstring/comment/doc edits.

* Update docs example model to use v0.3.4 API.

### Fixed

* Ensure consistency between input parameter `bounds` and `values` by
    always requiring dictionaries. Fix applies to `Elsewhere` and 
    `Photosphere`. Courtesy Sebastien Guillot.

* Gravitational mass doc typo fix.

### Changed

* Add input argument checks to `Likelihood.check` method.

* Add default `hypercube=None` to `Prior.inverse_sample_and_transform method.`

* If derived parameters found in subspace, assume an update is needed because
    cache mechanism not in place. (WIP.)


## [v0.3.3] - 2020-01-20

### Fixed

* At several places in the `Likelihood` class, calls were place to `self`,
    forgetting that `Likelihood` overwrites `ParameterSubspace.__call__`.
    Now calls are `super(Likelihood, self).__call__()` to obtain the current
    parameter vector.

## [v0.3.2] - 2020-01-16

### Summary

* Bug fixes. Backwards compatible.

* When initializing the ensemble-MCMC chains using an nd-ball, the inclusion
    in the prior support was checked by passing a vector to `Prior.__call__` but
    that code assumed that the parameter vector had already been assigned and
    can be accessed through the `ParameterSubspace`. As a result either an
    exception would be thrown (if parameter objects have no value set) or the
    support condition would be evaluated for some preset vector that does not
    change has we iterate through chains.

* The `Likelihood.check` method now has a fallback implementation given that
    the NumPy `allclose` function in v1.17 does not support Python 2.7.

### Attribution

* Based on testing by Sebastien Guillot.

### Fixed

* The `EnsembleSampler` so that it does not rely on the `CustomPrior.__call__`
    implementation to handle a vector argument. Chains should now be in
    prior support from the start and never leave.

* The `Likelihood.check` method so that a call to a `Likelihood` instance
    updates the parameters with a vector if the physical points are passed
    for value checking.

* The `Likelihood.check` method error error handling and if/else branching
    has been fixed.

* Some typographic errors in docs.

### Changed

* The way `EnsembleSampler` accesses the prior object.


## [v0.3.1] - 2019-12-12

### Fixed

* Some docstring and Sphinx-related formatting.


## [v0.3.0] - 2019-12-10

### Summary

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

### Attribution

* Feedback and ideas for the above development were discussed at an X-PSI
    workshop in Amsterdam, November 25-29 2019:
    Sebastien Guillot, Emma van der Wateren, Devarshi Choudhury, Pushpita Das,
    Anna Bilous, and Anna Watts.

### Added

* A new class `xpsi.Parameter` of which every model parameter is an instance.

### Changed

* The `xpsi.ParameterSubspace` class, which has far more sophisticated behaviours
    as a parameter container. The class, upon initialisation with arguments,
    also merges parameters and subspaces into a higher-dimensional (sub)space.
    Most other classes in the modelling language *inherit* from the
    `xpsi.ParameterSubspace` class.

* The `xpsi.TwoHotRegions` class is now dedicated to representing antipodally
    reflection-symmetric configurations only to simplify the choice of which
    class to use between `xpsi.HotRegions` and `xpsi.TwoHotRegions`. However,
    antipodally reflection-symmetric models can also be constructed using
    just `xpsi.HotRegions` because of the new *derived* parameter support. The
    may be a minor speed difference: `xpsi.TwoHotRegions`
    should be very slightly faster, but it might be imperceptible. Future
    warning: in the future `xpsi.TwoHotRegions` might removed altogther for
    simplication.

* The `xpsi.Photosphere` class can be instantiated to encapsulate only a
    reference to an `xpsi.Elsewhere` instance, and no `xpsi.HotRegion` instances.
    An `xpsi.Elsewhere` instance can by definition only generate a
    phase-invariant signal. However, further development is needed to handle
    this phase-invariant signal efficiently for likelihood functionality,
    given that operations with respect to phase are not required. Instead
    likelihood functions would be defined only with respect to energy.

### Removed

* The `xpsi.ParameterSpace` module. The global model parameter space is also
    simply an intance of the `xpsi.ParameterSubspace` class.

