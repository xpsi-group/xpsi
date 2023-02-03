History
-------

All notable changes to this project will be documented in this file.

The format is based on
`Keep a Changelog <http://keepachangelog.com/en/1.0.0/>`_
and this project adheres to
`Semantic Versioning <http://semver.org/spec/v2.0.0.html>`_.

.. REMOVE THE DOTS BELOW TO UNCOMMENT
.. ..[Unreleased]
.. ~~~~~~~~~~~~

.. Summary
.. ^^^^^^^

.. Fixed
.. ^^^^^

.. Added
.. ^^^^^

.. Changed
.. ^^^^^^^

.. Deprecated
.. ^^^^^^^^^^

.. Removed
.. ^^^^^^^

.. Attribution
.. ^^^^^^^^^^^

[v2.0.0] - 2023-02-03
~~~~~~~~~~~~~~~~~~~~~

Summary
^^^^^^^

* This major release candidate migrates X-PSI from Python2 (X-PSI v1.2.1 or lower) to Python3 (X-PSI v2.0 and higher), with corresponding updates and improvements to all documentation and tutorials.  

Fixed
^^^^^

* Debugging projection tool

Added
^^^^^

* Multi-version documentation so that users can view documentation/tutorials for either Python2 or Python3 (with warning on main page)
* Post-processing - adding names of parameters across diagonal in corner plots
* Extra yticks options for plotting functions in the tutorials
* `--noopenmp` install option for Mac Users
* Added option to fix the random seed for the synthetic data generation in Python3 version

Changed
^^^^^^^

* Modified all X-PSI routines to work in Python3
* General Documentation (Applications, Team and Acknowledgements, Citation, Future pages) updated for both Python2 and Python3 documentation branches.
* Installation and tutorial pages modified for Python3
* Module generator updated for Python3 and documentation added
* Projection tool updated for Python3 and documentation added
* Github actions modified to work in Python3
* Github actions modified to use mamba with install commands on one line to improve speed
* Updated references in the documentation and tutorial notebooks
* CustomInstrument channel_edges argument now changed to mandatory in tutorial notebooks and examples
* X-PSI Postprocessing now supports up-to-date versions of Nestcheck, Getdist. 

Deprecated
^^^^^^^^^^

*The Python2 version of X-PSI (v1.2.1) is now considered deprecated, although documentation and tutorials are still available.

Removed
^^^^^^^

* Removed requirement of FFMPEG for Animations in tutorials
* Suppressed printf() statements from c code in tutorial notebooks

Attribution
^^^^^^^^^^^

Devarshi Choudhury,
Bas Dorsman,
Sebastien Guillot,
Daniela Huppenkothen,
Yves Kini,
Tuomo Salmi,
Serena Vinciguerra,
Anna Watts


[v1.2.1] - 2022-12-12
~~~~~~~~~~~~~~~~~~~~~

Summary
^^^^^^^

* Hard upper bound for temperature increased from 7.0 to 7.6, allowed user an option to adjust when the exact likelihood calculation is skipped because of too bright signal, and license information updated.

Changed
^^^^^^^

* Strict bounds for temperature changed in ``xpsi/HotRegion.py``, ``xpsi/Everywhere.py``, and ``xpsi/xpsi/Elsewhere.py`` to allow analysis for hotter neutron stars.

* Added mention in ``xpsi/HotRegion.py``, ``xpsi/Everywhere.py``, and ``xpsi/xpsi/Elsewhere.py`` that the user should set the parameter bounds to be within the values given in the numerical atmosphere table.

* Added a new input parameter ``slim`` to ``xpsi/likelihoods/default_background_marginalisation.pyx``, which can be used to adjust when the exact likelihood calculation is skipped because of the signal being too bright compared to the data. The default value of this parameter is set to the same value as in the code before (20.0).

* Made the warning in synthesise function in ``xpsi/Likelihood.py`` more accurate.

* Fetched the prior to likelihood object in ``examples/examples_fast/Synthetic_data.ipynb`` to make sure prior bounds are checked when synthesising data.

* License of X-PSI was changed from MIT to GPLv3.

Attribution
^^^^^^^^^^^

Tuomo Salmi,
Yves Kini,
Sebastien Guillot,
Anna Watts


[v1.2.0] - 2022-12-05
~~~~~~~~~~~~~~~~~~~~~

Summary
^^^^^^^

* Included a new numerical atmosphere extension in a ``xpsi/surface_radiation_field/archive/hot/`` directory allowing freedom in the predicted atmospheric beaming pattern.

Added
^^^^^

* ``xpsi/surface_radiation_field/archive/hot/numerical_fbeam.pyx``: New numerical atmosphere extension with additional beaming parameters.

* ``examples/examples_modeling_tutorial/TestRun_NumBeam.py``: An example run using the new atmosphere extension.

* ``examples/examples_modeling_tutorial/modules``: Additional modules (e.g. a CustomHotRegion) needed by the new example run.

Changed
^^^^^^^

* ``Setup.py`` file changed to include the option for installing with new atmosphere extension.

* Documentation page for "Example script and modules" updated to include the new example. 

Attribution
^^^^^^^^^^^

Tuomo Salmi


[v1.1.0] - 2022-11-14
~~~~~~~~~~~~~~~~~~~~~

Summary
^^^^^^^

* Additional tools included in a ``xpsi/utilities`` directory for plotting hot regions on a sphere and performing importance sampling in X-PSI. Documentation for these tools is to be appended later. In addition, the internet documentation compilation was automated using GitHub actions for every merged pull request.

Added
^^^^^

* ``xpsi/utilities/ProjectionTool.py``: Tool for projecting hot regions.

* ``xpsi/utilities/ImportanceSample.py``: Tool for calling X-PSI importance sampling.

Changed
^^^^^^^

* ``Setup.py`` file changed to include the new utilities directory.

* Documentation is now compiled automatically using ``.github/workflows/build_docs.yml`` every time merging a pull request into the main branch.

Attribution
^^^^^^^^^^^

Serena Vinciguerra,
Daniela Huppenkothen,
Tuomo Salmi,
Devarshi Choudhury


[v1.0.0] - 2022-09-26
~~~~~~~~~~~~~~~~~~~~~

Summary
^^^^^^^

* This major release contains minor bug fixes, improved error messages, as well as improved documentation and tutorials (jupyter notebooks).  This release coincided with the submission of an X-PSI article to the `Journal of Open Source Software <https://joss.theoj.org/>`_


Fixed
^^^^^

Added
^^^^^

* A modification of the ``setup.py`` with flags (``--NumHot`` and ``--NumElse``) now facilitates switching between surface emission models.

* The post-processing module has now an option to show the credible intervals of each parameter and run (above the 1D distribution of the corner plot) when multiple runs are plotted in the same figure (but not working for multiple models yet). The appropriate tutorial notebook is also provided.

* Some unit tests and continuous integration.

* A tutorial landing page and a link to a dedicated Zenodo repository for large files needed to run the tutorials. 

Changed
^^^^^^^

* The general documentation has been improved, reorganized and clarified.  More details are provided for the installation, locally and on HPC systems.

* The messages of several possible errors have been clarified and detailed to help the user resolve them.

* A small modification now allows production runs without importing matplotlib.

* All tutorials have been updated and improved.

Deprecated
^^^^^^^^^^

Removed
^^^^^^^

* Method ``fixed_spin`` of ``spacetime.py`` module.  A spacetime with fixed spin can be created by specifying a spin frequency ``value`` and omitting the spin frequency ``bounds``

Attribution
^^^^^^^^^^^

Devarshi Choudhury,
Bas Dorsman,
Sebastien Guillot,
Daniela Huppenkothen,
Yves Kini,
Tuomo Salmi,
Serena Vinciguerra,
Anna Watts

[v0.7.12] - 2022-09-15
~~~~~~~~~~~~~~~~~~~~~~~~~~

Summary
^^^^^^^

* Since version 0.7.11. a few changes have been made including updates to the documentation and the handling of numerical problems in ray tracing. The latter fix can potentially have a small effect on the calculated pulse profiles and likelihood values for some parameter vectors, but according to testing that effect is very minor at most.


Fixed
^^^^^

* Numerical problem in  ``xpsi/cellmesh/rays.pyx`` for certain paramaters causing sporadic warnings in later computation. This is prevented by allowing small rounding errors when checking if sin_alpha parameter is unity, and in case NaNs still occur, replacing them with zero (T.S.).

* Comment about returned variables updated to include the best-fitting background limited by the support in ``xpsi/likelihoods/default_background_marginalisation.pyx`` (T.S.).

* The photosphere object validity check in ``xpsi/Star.py`` which incorrectly failed if all photosphere parameters were fixed (D.C., Y.K., T.S.).

Added
^^^^^

* Added more information and warnings about about switching between the blackbody and numerical atmosphere extensions in the documentation for Installation, Surface radiation field tools and (FA)Q pages. Added also a links to the Zenodo publication of Riley+2021 from where the numerical atmosphere data can be obtained (T.S.).

* Added a new kwargs ("prior_samples_fnames") used in ``xpsi/PostProcessing/_corner.py`` to allow user to set the name of file from where the prior samples are read/saved (T.S.).

* Added comments about the new kwargs (introduced already in version 0.7.11) in the function descriptions used in ``xpsi/PostProcessing/_corner.py`` visible also for the documentation (T.S.).

* Added an option to force update ``xpsi/Star.py`` to avoid errors, for example, when all paremeters are fixed and X-PSI thinks otherwise that updating can be skipped (T.S., D.C., Y.K.).

* Added options allowing the user to truly force update the likelihood in ``xpsi/Likelihood.py`` and avoid errors caused by the automatic need-update-checks not working for all the possible cases. Added also an error message suggesting to use those options if the usual "AttributeError: 'CustomSignal' object has no attribute '_loglikelihood'" would be encountered (T.S.).

Changed
^^^^^^^

Deprecated
^^^^^^^^^^

Removed
^^^^^^^

Attribution
^^^^^^^^^^^

* Tuomo Salmi (T.S.), Devarshi Choudhury (D.C.), and Yves Kini (Y.K.)


[v0.7.11] - 2022-08-22
~~~~~~~~~~~~~~~~~~~~~~

Summary
^^^^^^^

* Since version 0.7.10, a few bugs have been fixed in the module generator, error handling and postprocessing. Also, new error/warning messages are given if trying to use wrong atmosphere extension. In addition, some improvements have also been added to the postprocessing (possibility to e.g. save and read the drawn priors to produce corner plots much faster), without mentioning them in the documentation yet.


Fixed
^^^^^

* Bug in ``xpsi/EnsembleSampler.py`` when initializing walkers. Need to use "self._prior" instead of "prior" (Y.K.).

* Bug (typo) in ``xpsi/PostProcessing/_pulse.py`` when plotting the true signal. Need to use "component" instead of "eomponent" (G.L.).

* Several bugs (typos) in ``xpsi/PostProcessing/_spectrum.py`` when plotting the true signal (T.S., G.L.).

* Issues with ``xpsi/PostProcessing/_corner.py`` not being able to plot the cross hairs for true parameter values in the corner plot if only a subset of model parameters chosen for the figure (T.S., Y.K.).

* Error handling in ``xpsi/Signal.py`` when the number of event data channels does not match the number of the instrument data channels (S.G.).

* Fixed reference to incident_background in the modeling tutorial (B.D.).

* Several bug fixes in ``xpsi/module_generator.py`` (D.C.).

Added
^^^^^

* Added a warning message in the blackbody atmosphere extension  ``xpsi/surface_radiation_field/hot.pyx`` if providing numerical atmosphere data (T.S.).

* Added an error message in the numerical atmosphere extension  ``xpsi/surface_radiation_field/archive/hot/numerical.pyx`` before a segmentation fault error caused by not loading the numerical atmosphere data (T.S.).

* Added a warning when trying to synthetize data in ``xpsi/Likelihood.py`` with input parameters outside of the defined prior bounds, finishing without errors but with no data produced (Y.K. & T.S.).

* Added option for the user to set the line colors for different runs in ``xpsi/PostProcessing/_corner.py`` using kwargs (T.S.).

* Added possibility to save and read the previously drawn prior samples in ``xpsi/PostProcessing/_corner.py`` using "force_draw" kwargs (T.S.).

* Added possibility to plot the priors only for the first run in ``xpsi/PostProcessing/_corner.py`` using "priors_identical" kwargs, if known that priors are the same for all runs (T.S.).

* Saved credible intervals in numerical format that can be accessed after plotting the corner plot (see "val_cred" in ``xpsi/PostProcessing/_corner.py`` and ``xpsi/PostProcessing/_postprocessor.py``) (Y.K., T.S.).

Changed
^^^^^^^

Deprecated
^^^^^^^^^^

Removed
^^^^^^^

Attribution
^^^^^^^^^^^

* Tuomo Salmi (T.S.), Yves Kini (Y.K.), Devarshi Choudhury (D.C.), Bas Dorsman (B.D.), Gwénaël Loyer (G.L.), and Sebastien Guillot (S.G.)


[v0.7.10] - 2022-02-10
~~~~~~~~~~~~~~~~~~~~~~

Summary
^^^^^^^

* Since version 0.7.9, several bugs have been fixed. For example, fixing the double counting of the second component of a dual temperature region when caching turned on. Also, documentation and example scripts have been updated.


Fixed
^^^^^

* Bug in ``xpsi/Signal.py`` when looping over dual temperature components while using caching (D.C., T.S, S.V.). 

* Bug in ``xpsi/Signal.py`` merging the new phase-shift parameter to the parameter subspace (T.S. & D.C.).

* Missing global argument added in ``xpsi/module_generator.py`` (D.C.).

* Documentation and example scripts updated and fixed to work with newest X-PSI versions (S.G.).

* Bug in ``xpsi/PostProcessing/_corner.py`` not showing true values correctly in corner plots for simulated data (T.S. & Y.K.).

* Corrected the link to the documentation pages when importing X-PSI (D.C. & T.S.).

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

* Devarshi Choudhury (D.C.), Tuomo Salmi (T.S.), Serena Vinciguerra (S.V.), Sebastien Guillot (S.G.), and Yves Kini (Y.K.)


[v0.7.9] - 2021-11-26
~~~~~~~~~~~~~~~~~~~~~

Summary
^^^^^^^

* New program that automates generation of model modules for common usage
  patterns, in particular the NICER modelling workflow. The program may be
  located at ``xpsi/module_generator.py`` and executed as
  ``python module_generator.py -h`` to see the usage help.


Fixed
^^^^^

* The :class:`~.Background` call method body template and fixed the
  :class:`~.Signal` class to access the correct property of the background
  instance.

* Documentation URLs to reference the organisation repository. (D.H.)


Added
^^^^^

* Functionality to the :class:`~.Data` class method for event handling so that
  it can load events from file when the energy in eV is given.

* Optional maximum energy to use for ray-tracing simulations. Useful if there
  is a background component such as a powerlaw that is jointly modelled with
  higher-energy event data using a subset of instruments.

* A phase-shift parameter for each :class:`~.Signal` instance. If there are
  two or more phase-resolved data-sets, there may be a need to have a phase-
  shifting parameter for each signal. For phase-summed data sets, the phase-
  shift can be arbitrarily fixed. Phase-shifts can be derived from other
  phase-shifts, and one signal's phase-shift can always be fixed as zero and
  thus locked to the phase shifts of the hot regions.


Attribution
^^^^^^^^^^^

* Daniela Huppenkothen (D.H.).


[v0.7.8] - 2021-09-22
~~~~~~~~~~~~~~~~~~~~~

Fixed
^^^^^

* Correction in the importance sampling function. If the number of MPI
  processes is a factor of the number of samples reweighted, a subset of
  samples, with cardinality equal to the size of the MPU world, was not
  reweighted but is included for renormalisation with the same weight as the
  input weight. E.g., if there is one MPI process, then the last sample is not
  reweighted, so the output weight is equal to the input weight. (S.V.)
* Correction of the image appearing on the :mod:`~.HotRegion` page. (S.V.)
* Minor typos corrected. (T.S. & Y.K.)

Changed
^^^^^^^

* Updated the :func:`~.tools.synthesise_exposure` and
  :func:`~.tools.synthesise_given_total_count_number` functions to handle zero
  background and make sure that the input background memory buffer does not get
  modified by the synthesis routines. (T.S. & Y.K.)
* Added a keyword argument to the default background marginalisation function
  to enable passing of a background signal in the form of a channel-phase
  interval buffer. The background should already be averaged over phase
  intervals, having units of counts/s. Useful for phase-dependent backgrounds,
  or a phase-independent background if the channel-by-channel background
  variable prior support is restricted.

Added
^^^^^

* Updates to the project acknowledgements page of the documentation.

Attribution
^^^^^^^^^^^

* Serena Vinciguerra (S.V.), Yves Kini (Y.K.), and Tuomo Salmi (T.S.).


[v0.7.7] - 2021-06-24
~~~~~~~~~~~~~~~~~~~~~

Fixed
^^^^^

* Bugs in mesh cell allocation routine. These bugs occur for some specific
  subset of hot regions with both a superseding member region and a ceding
  member region and both radiate. This bug does not affect any production
  analyses to date, but was encountered by D.C. when preparing a model with
  such a hot region for posterior sampling.
* Importance sampling bug when reweighting the likelihood function.

Added
^^^^^

* Guidelines to the documentation for dependency citation.
* Tips for installing X-PSI on a macOS in the documentation (S.V. & D.C.).
* Some additional lines to install X-PSI on SURFsara's Cartesius (S.V.).
* Instructions to install X-PSI on SURFsara's Lisa (T.S.).

Attribution
^^^^^^^^^^^

* With thanks to Devarshi Choudhury (D.C.) for noticing and investigating
  potentially buggy mesh construction behaviour that was, indeed, buggy.
* With thanks to Serena Vinciguerra for noticing and investigating
  potentially buggy importance sampling behaviour that was, indeed, buggy.
* With thanks to Serena Vinciguerra (S.V.), D.C., and
  Tuomo Salmi (T.S.) for patches to documentation install instructions.

[v0.7.6] - 2021-05-16
~~~~~~~~~~~~~~~~~~~~~

Summary
^^^^^^^

* NB: This patch is unfortunately not backwards compatible. This patch has been
  pushed nevertheless to comply with a NICER collaboration publication which
  uses X-PSI v0.7 with some features from a development version. The analysis is
  open-source, so the development features used have been pushed in this patch.
  The next minor release will officially include these tested features together
  with documentation.

* New skymap plotting functionality and an MPI-capable importance sampling
  method that can handle likelihood function and prior PDF changes. New
  documentation and examples will be made available in the future.

Changed
^^^^^^^

* The extension module for default background marginalisation returns a tuple
  with an extra element. This is probably backwards incompatible with custom
  subclasses of the :class:`~.Signal` class.

Added
^^^^^

* Skymap plotting functionality. Examples will be added to the documentation
  in a future patch. The most useful feature is plotting a skymap time-series
  so that the image of the model surface hot regions rotates across and down
  a static figure. This is useful for papers to summarise an animated figure.
  This feature is functional but still being tested and developed.

* An MPI-capable importance sampling method that can handle likelihood function
  and prior PDF changes. This is useful to save computation time. This feature
  is being tested and developed.

Fixed
^^^^^

* A bug in :meth:`~.Likelihood.Likelihood.check` that prevented checking
  the likelihood function for more than one point.

Attribution
^^^^^^^^^^^

* With thanks to Serena Vinciguerra (S.V.) for testing importance sampling.


[v0.7.5] - 2021-02-10
~~~~~~~~~~~~~~~~~~~~~

Fixed
^^^^^

* Corner-case stability improvements for default background marginalisation.
* If likelihood function is below :attr:`~.Likelihood.Likelihood.llzero` after
  evaluation, the parameter vector is included in the prior support as
  intended.
* Typo in ``_precision`` function in ``xpsi/PostProcessing/__init__.py``. (S.V.)
* Math typo on the :mod:`~.HotRegion` page. (S.V.)
* Explanatory text in the multiple-imaging tutorial. (T.S.)

Changed
^^^^^^^

* A few image components appearing on the :mod:`~.HotRegion` page. (S.V.)
* Bounds exception now prints the name of the offending parameter in
  :class:`~.Parameter.Parameter`. (S.V.)

Added
^^^^^

* An extension module for calculating hot region local variables from global
  variables for hot region configurations under the umbrella of the PST-U model
  introduced in `Riley et al. (2019) <https://ui.adsabs.harvard.edu/abs/2019ApJ...887L..21R/abstract>`_.

Attribution
^^^^^^^^^^^

* With thanks to Serena Vinciguerra (S.V.) and Tuomo Salmi (T.S.).


[v0.7.4] - 2021-01-26
~~~~~~~~~~~~~~~~~~~~~

Fixed
^^^^^

* Missing packages in ``setup.py`` causing errors when importing xpsi.
* A few typos in the documentation.

Added
^^^^^

* A few images in the documentation.

Attribution
^^^^^^^^^^^

* Serena Vinciguerra, Yves Kini, Devarshi Choudhury.


[v0.7.3] - 2020-11-12
~~~~~~~~~~~~~~~~~~~~~

Fixed
^^^^^

* Phase-averaging issue that can sometimes occur due to numerical effects when
  comparing two numbers that should be the same but can differ by tiny degrees
  at machine precision level.
* Some documentation typographic errors.


[v0.7.2] - 2020-11-04
~~~~~~~~~~~~~~~~~~~~~

Fixed
^^^^^

* Error raised while running ``setup.py`` for linking rayXpanda with
  clang compiler.

Attribution
^^^^^^^^^^^

* Serena Vinciguerra.


[v0.7.1] - 2020-10-01
~~~~~~~~~~~~~~~~~~~~~

Fixed
^^^^^

* An ``AttributeError`` raised during runtime linking to the fallback rayXpanda
  implementation.

Attribution
^^^^^^^^^^^

* With thanks to Devarshi Choudhury for bug testing.


[v0.7.0] - 2020-09-30
~~~~~~~~~~~~~~~~~~~~~

Summary
^^^^^^^

* New plotting functionality.
* Should be backwards compatible, but some small internal tweaks or default
  behaviour changes could result in small differences in plots that might not
  even be discernable.

Added
^^^^^

* Option to specify only the number of phases per cycle when calling
  :meth:`~.Photosphere.Photosphere.image`, instead of having to supply the
  phase set.
* New plot type for animated photon specific intensity skymaps with their
  associated photon specific flux pulse-profiles and the photon specific flux
  spectrum that connects the signals at those energies. See the documentation
  of the :meth:`~.Photosphere.Photosphere.image` method for options, details,
  and an example.
* Example plots to the :class:`~.Photosphere.Photosphere` documentation.
* New helper methods :meth:`~.Photosphere.Photosphere.write_image_data`
  and :meth:`~.Photosphere.Photosphere.load_image_data` to write ray map data,
  photon specific intensity image data, and photon specific flux signal data to
  disk, and then read the data back into memory as attributes so that the data
  can be reused to accelerate calls to calculate images and generate static and
  animated plots.
* Option to :meth:`~.Photosphere.Photosphere._plot_sky_maps`,
  ``add_zero_intensity_level``, that applies a colormap such that zero intensity
  corresponds to the lowest colour. In this case a non-radiating part of the
  stellar surface, and the background sky, have well-defined colour. If lowest
  colour in the colormap is instead associated with the lowest finite intensity
  in the skymap panel, then the background sky for instance is assigned the same
  colour so that the least bright part of the image merges with the background
  sky colour. The latter choice resolves the variation in the intensity as a
  function of phase and sky direction better with colour, but the former might
  give more of an indication of the magnitude of the variation in intensity
  as a function of phase and sky direction relative to the background sky.

Changed
^^^^^^^

* A phase set supplied to :meth:`~.Photosphere.Photosphere.image` can have
  units of cycles, not radians as was previously the requirement, by setting
  the ``phase_in_cycles`` keyword argument to ``True`` if the supplied phase
  array as units of cycles.
* The photon specific flux can be calculated with
  :meth:`~.Photosphere.Photosphere.image` at far more energies than photon
  specific intensities are cached at, by using the :obj:`cache_energy_indices`
  keyword to supply and array of integers to index the energy array. This
  saves memory and means that imaging with an extension module can be executed
  once to generate both skymaps (which require cached intensities but only
  typically at a few representative energies) and the photon specific flux
  (which does not require cached intensities, but typically is computed for
  a much finer energy array).

Attribution
^^^^^^^^^^^

* With thanks to Anna Bilous and Serena Vinciguerra for helpful suggestions
  about the new animated plot type.


[v0.6.3] - 2020-10-01
~~~~~~~~~~~~~~~~~~~~~

Fixed
^^^^^

* An ``AttributeError`` raised during runtime linking to the fallback rayXpanda
  implementation.

Attribution
^^^^^^^^^^^

* With thanks to Devarshi Choudhury for bug testing.


[v0.6.2] - 2020-09-28
~~~~~~~~~~~~~~~~~~~~~

Fixed
^^^^^

* Bug in :func:`~.Sample.nested` when initialisation of nested sampler class
  tries to call ``set_default`` dictionary method instead of the correct
  ``setdefault`` method.
* Import errors associated with the :mod:`~.PostProcessing` module.

Changed
^^^^^^^

* The :attr:`~.Parameter.Parameter.cached` property of a
  :class:`~.Parameter.Parameter` instance can be set to ``None``.
* The :class:`~.ParameterSubspace.ParameterSubspace` initialiser is decorated
  to avoid verbose output by every MPI process.
* The :class:`~.Prior.Prior` uses the class attribute
  ``__draws_from_support__`` to set the number of Monte Carlo draws from the
  joint prior support to require to set the MultiNest hypervolume expansion
  factor appropriately. The default value is ``5``, which means :math:`10^5`
  draws from the joint prior support.
* Checks if an instance of  ``six.string_types`` in
  :class:`~.PostProcessing._metadata.Metadata`, e.g., to allow unicode strings
  in posterior ID labels.


[v0.6.1] - 2020-09-14
~~~~~~~~~~~~~~~~~~~~~

Fixed
^^^^^

* Bug wherein multiple :class:`~.Signal.Signal` instances passed to a
  :class:`~.Likelihood.Likelihood` instance do not have references stored.
* The :mod:`~.tools` synthesis functions adhering to the global phase
  interpolant switch, and updated tutorial accordingly.

Changed
^^^^^^^

* The :meth:`~.Data.Data.phase_bin__event_list` constructor signature, so that
  the phase and channel columns can be arbitrary.

Removed
^^^^^^^

* An unused prototype extension module.


[v0.6.0] - 2020-09-05
~~~~~~~~~~~~~~~~~~~~~

Summary
^^^^^^^

* Backwards compatible for most use cases, but possible corner cases.
* Includes a non-critical, but important patch for animating intensity skymaps,
  and updates to the environment file for cloning.
* The new feature is support for higher-order images when invoking an integrator
  that discretises the surface (with a regular mesh). Secondary images can
  be very important, whilst tertiary images less so. Quaternary, quinary, and
  possibly senary images can sometimes be detected and included too, with
  accuracy that decreases with order. Fortunately, the contribution to the
  photon specific flux generally decays rapidly with image order beyond the
  secondary or tertiary images. The computational cost scales almost
  linearly with order *if* an appreciable fraction of every iso-latitudinal ring
  on the surface is multiply-imaged at each order. Note that multiple-imaging
  manifests entirely naturally when an image-plane is discretised in such away
  that the regular mesh resolves the stellar limb sufficiently well, where
  higher-order images get insanely squeezed.

Fixed
^^^^^

* The memory consumption problem of the animator method in
  :class:`~.Photosphere.Photosphere`. Now animation should generally require
  an entirely tracable amount of memory.

Added
^^^^^
.. _rayXpanda: <https://github.com/ThomasEdwardRiley/rayXpanda>

* Multiple-imaging support including an option to specify the maximum image
  order to iterate up to, with automatic truncation when no image at a given
  order is detected. If no limit is specified (the default), then images are
  included as far as they can be detected given the numerical resolution
  settings, which is typically between quaternary and senary images.
* A multiple-imaging tutorial.
* A global switch for changing phase and energy interpolants without
  recompilation of extensions. To change interpolants, you can use top-level
  functions :func:`xpsi.set_phase_interpolant` and
  :func:`xpsi.set_energy_interpolant`. Generally computations are more
  sensitive to the phase interpolants, of which the options from GSL are:
  Steffen spline (pre-v0.6 choice), Akima periodic spline, and cubic periodic
  spline. The default choice is now an Akima periodic spline in an attempt to
  improve interpolation accuracy of the interpolant at function maxima, where
  the accuracy is generally most important in the context of likelihood
  evaluations.  Note that in some corner cases, the signal from a hot region is
  negative in specific flux because there is a correction computed to yield the
  intended signal from :class:`~.Elsewhere.Elsewhere` when it is partially
  masked by hot regions. In this case, when using phase interpolant tools from
  the :mod:`~.tools` and :mod:`~.likelihood` modules it is necessary to use a
  ``allow_negative`` option when calling the tools to specify that a negative
  interpolant is permitted.
* Automatic linking of the package rayXpanda_ for calculation of the inverse of
  the deflection integral, and it's derivative via a high-order symbolic
  expansion, for a subset of primary images. The purpose is to mainly as an
  orthogonal validation of a subset of integrals executed via numerical
  quadrature and inversion via spline interpolation.  The other reason is
  because to support multiple-imaging with the surface-discretisation
  integrators this aforementioned interpolation had to change due to
  non-injectivity of functions when interpolating with respect to the cosine of
  the deflection angle. However, to calculate the convergence derivative
  sufficiently accurately, interpolating with respect to the cosine of the
  deflection seems necessary. Therefore rayXpanda_ can be linked in, if it is
  available, for low deflection angles instead of avoid having to allocate
  additional memory and construct splines specifically for low-deflection
  primary images. Simple testing suggests there are no valuable speed gains,
  however, possibly because the high-order expansion and simultaneous evaluation
  of the polynomial and it's derivate with a nested Horner scheme itself
  requires a substantial number of floating point operations.
* A helper method :meth:`~.ParameterSubspace.ParameterSubspace.merge` that
  merges a set of parameters, or a parameter subspace, or a set of subspaces,
  into a subspace that has already been instantiated.

Changed
^^^^^^^

* Updated the Conda ``environment.yml`` file for replication of the development
  environment. The ``basic_environment.yml`` file was also updated in an
  earlier release in an additional necessary package, ``wrapt``.

Deprecated
^^^^^^^^^^

* The ``repeat``, ``repeat_delay``, and ``ffmpeg_path`` keyword arguments for
  the animator method in :class:`~.Photosphere.Photosphere`. These were
  ultimately not effective. To repeat the animation intrinsically, set the
  number of ``cycles``, and extrinsically, this can be looped when embedded in
  another environment.


[v0.5.4] - 2020-09-01
~~~~~~~~~~~~~~~~~~~~~

Fixed
^^^^^

* Bug due to local variable ``NameError`` when setting instrument channel
  energy edges.
* Bug that prevented a hot region phase parameter from being a fixed or derived
  variable.

Attribution
^^^^^^^^^^^

* With thanks to Devarshi Choudhury.


[v0.5.3] - 2020-08-14
~~~~~~~~~~~~~~~~~~~~~

Summary
^^^^^^^

* Improvement patches. Deliberately backwards incompatible for safety in
  memory allocation.

Fixed
^^^^^

* Add try-except block to :attr:`~.Photosphere.Photosphere.global_to_local_file`
  property so that explicit setting of ``None`` by user is not required if
  file I/O is not needed in the extension module. Actually, ``None`` could
  not be set for the property anyway due to type checking.
* Bug when declaring that sky maps should be animated and memory freed
  beforehand.

Added
^^^^^

* The surface to image-plane ray map is cached in Python process memory so it
  can be efficiently reused for same spacetime configuration and ray map
  resolution settings. Explicit support for writing the ray map to disk and
  loading it is not included, but this should be entirely possible to achieve
  manually. Backwards compatible except for corner cases, such as not using
  keyword arguments when calling :meth:`~.Photosphere.Photosphere.image`, or if
  resolution settings changed between calls to the imager but a ray map
  otherwise exists in Python process memory and the spacetime configuration has
  not been changed.
* A secret keyword argument to :meth:`~.Photosphere.Photosphere.image`,
  :obj:`_OVERRIDE_MEM_LIM`, which can be used to change an internal hard limit
  on the intensity cache size. This setting is for safety and designed so that
  higher memory consumption is deliberate or if something goes awry, it is
  deemed the responsibilty of the user to have read method docstring carefully.
  The tutorials will not use this secret keyword, so if the user tries to run
  them and encounters an exception, they will need to investigate the docstring
  and either adapt the resolution to their system or take the responsibility of
  setting the cache size limit for their system to accomodate the resolution
  settings in the tutorial.
* Optional argument to :meth:`~.Photosphere.Photosphere.image`,
  :obj:`single_precision_intensities`, which flags whether or not to *cache*
  the intensities in single precision do halve intensity cache memory
  requirements. The default is to cache in single precision.
* Verbosity to :meth:`~.Photosphere.Photosphere.image` because execution
  can take many minutes depending on settings chosen. The verbosity
  can be deactivated via a keyword argument (see the method docstring).

Changed
^^^^^^^

* The usage of the :meth:`~.Photosphere.Photosphere.image` argument
  :obj:`cache_intensities`. Instead of simply activating intensity caching
  with boolean, the user must specify a cache size limit that is adhered to.
  If the required cache size given the resolution settings is larger than
  the limit, imaging does not proceed. If the cache size limit is zero or
  equivalent, then imaging safely proceeds without caching the intensities.
* Intensities are by default *cached* in single precision to reduce cache memory
  requirements.


[v0.5.2] - 2020-08-12
~~~~~~~~~~~~~~~~~~~~~

Summary
^^^^^^^

* Python API: small backwards compatible patches to add useful features.
* C API: small backwards incompatible patch to support Python API patch.

Added
^^^^^

* Support for hyperparameters (i.e., parameters of the prior distribution),
  by making :class:`~.Prior.Prior` inherit from
  :class:`~.ParameterSubspace.ParameterSubspace`. Custom hyperparameters can
  then be defined in a subclass initiliser, or otherwise. The hyperparameters
  are merged into the :class:`~.Likelihood.Likelihood` parameter subspace as
  mostly normal parameters (with small caveat in the form of property
  :attr:`~.Parameter.Parameter.is_hyperparameter`) and can have their own
  prior (the hyperprior) implemented in a :class:`~.Prior.Prior` subclass along
  with the other free parameters in the model. A tutorial will be delivered in
  due course. These modifications are backwards compatible.
* Simple support for transforming from global to local variables (for image-
  plane calculations) with the help of a file on disk, whose path can be
  specified dynamically in Python and relayed to the relevant extension where a
  custom model implemention can do I/O with the file. This is useful if one has
  a set of files containing precomputed data, but understandably does not want
  to do filesystem acrobatics or recompile an extension every time the file
  path changes. Setting the file path dynamically in this way is akin to
  changing the value of some discrete variable in the mapping between global
  and local variables. With thanks to Anna Bilous for the suggestion. A tutorial
  will be delivered when possible.
* Added :attr:`~.Instrument.Instrument.channel_edges` property, and updated
  tutorials to reflect this new concrete implementation.

Changed
^^^^^^^

* The ``init_local_variables`` function signature in the header
  ``xpsi/surface_radiation_field/local_variables.pxd``, and in the
  corresponding ``xpsi/surface_radiation_field/archive/local_variables``
  extensions. You would have to modify a custom extension module manually to
  match the function signature declared in the header.

Fixed
^^^^^

* Removed remnant manual Sphinx method signatures; the decorator now preserves
  the method signature so automated Sphinx doc works on those decorated methods.
* Updated package docstring to reflect name change.
* Uses of ``xpsi.Data.channel_range`` property to adhere to future deprecation.


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
  performance) the relevant signal plots from `Riley et al. (2019) <https://ui.adsabs.harvard.edu/abs/2019ApJ...887L..21R/abstract>`_, as demonstrated
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
