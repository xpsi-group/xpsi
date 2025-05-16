.. _landing_page_tutorials:

==========
Start here
==========

All tutorials listed below may be found as Jupyter notebooks under ``docs/source/``. We encourage running these yourself as you go through the tutorials as a means of practice. Any external data files needed such as instrument response files and numerical atmosphere data are available on `Zenodo <https://doi.org/10.5281/zenodo.7094144>`_. These data files should be saved in ``examples/examples_modeling_tutorial/model_data/``.

**Overview of the available tutorials:**

* :doc:`X-PSI 101<XPSI_101>` is a basic tutorial introducing the key aspects of X-PSI.  In it you will build a very simple model of a single hot spot on a neutron star surface, simulate a pulse profile, and perform parameter estimation.  Start here if you have never used X-PSI before!  

* :doc:`Modeling<Modeling>` is an in-depth tutorial that first covers the structure of X-PSI followed by basic usage of X-PSI. Basic usage includes constructing your first model (e.g. instrument, star, and atmosphere) and constructing a likelihood for nested sampling.

* :doc:`Instrument synergy<Instrument_synergy>` shows how to construct a joint likelihood with data from two instruments.

* :doc:`Hot region complexity<Hot_region_complexity>` goes into more detail about multiple and especially multiple overlapping hot regions. We note that :doc:`hotregion<hotregion>` is a useful page that contains overview figures of various overlapping cases.

* :doc:`Global surface emission<Global_surface_emission>` concerns surface emission from fields that span the full star and showcases three different options for signal integration. The first is a globally uniform temperature field, which allows for a time invariant signal integration. In the second case the temperature field is phase dependent and thus requires time dependent integration. The third case is general purpose and discretises a distant image plane instead of the stellar surface.

* :doc:`Surface radiation field tools<Surface_radiation_field_tools>` demonstrates the usage of the default (blackbody) and an alternative (atmosphere interpolated from precomputed data) surface radiation field module to compute photon specific intensities. It also shows how to do beaming pattern and spectrum plots for the radiation fields.

* :doc:`Modeling (without statistics)<Modeling_without_statistics>` is similar to :doc:`Modeling<Modeling>` but omits any statistical inference and adds various plots of signals. This is a useful tutorial if you are only interested to use X-PSI to create synthetic data.

* :doc:`Polarization<Polarization>` is a tutorial that shows how to model polarized X-rays in X-PSI.

* :doc:`Post-processing<Post-processing>` is a simple exercise on how to use the main post-processing tools provided by X-PSI.

* :doc:`Emitting patterns 2D projection<Emitting_patterns_2Dprojection>` is a tutorial showing how to use the 2D projection tool for the hot emitting regions on the star's surface.

* :doc:`Accretion disk<Accretion_disk>` is a tutorial showing how to set up an accretion disk and add its emission to the emission of a star.

* :doc:`Multiple imaging<Multiple_imaging>` is a tutorial studying the effects of multiple imaging.

* :doc:`Importance sampling<Importance_sampling>` is a tutorial for importance sampling.

* :doc:`Module generator tutorial<Module_generator_tutorial>` provides instructions on how to generate Python modules in an automated way to run X-PSI.

* :doc:`Example script and modules<Example_script_and_modules>` shows a couple of simple example scripts for pulse shape computation and nested sampling.

* :doc:`Example job<Example_job>` contains example job scripts for computation on clusters.

* :doc:`Posterior inference using SBI<x_p_sbi>` is a tutorial showing how to use Simulation-Based Inference (SBI) to obtain posteriors.
