.. _landing_page_tutorials:

==========
Start here
==========

All tutorials listed below may be found as Jupyter notebooks under ``docs/source/``. We encourage running these yourself as you go through the tutorials as a means of practice. Any external data files needed such as instrument response files and numerical atmosphere data are available on `Zenodo <https://doi.org/10.5281/zenodo.7094145>`_.

**Overview of the available tutorials:**

* :doc:`Modeling<Modeling>` is an in-depth tutorial that first covers the structure of X-PSI followed by basic usage of X-PSI. Basic usage includes constructing your first model (e.g. instrument, star, and atmosphere) and construct a likelihood for nested sampling.

* :doc:`Instrument synergy<Instrument_synergy>` shows how to construct a joint likelihood with data from two instruments.

* :doc:`Hot region complexity<Hot_region_complexity>` goes into more detail about multiple and especially multiple overlapping hot regions. We note that :doc:`hotregion<hotregion>` is a useful page that contains overview figures of various overlapping cases.

* :doc:`Surface radiation field tools<Surface_radiation_field_tools>` demonstrates the usage of the default (blackbody) and an alternative (atmosphere interpolated from precomputed data) surface radiation field module to compute photon specific intensities. It also shows how to do beaming pattern and spectrum plots for the radiation fields.

* :doc:`Modeling (without statistics)<Modeling_without_statistics>` is similar to :doc:`Modeling<Modeling>` but omits any statistical inference and adds various plots of signals. This is a useful tutorial if you are only interested to use X-PSI to create synthetic data.

* :doc:`Post-processing<Post-processing>` is a simple exercise on how to use the main post-processing tools provided by X-PSI.

* :doc:`Example script and modules<Example_script_and_modules>` shows a couple of simple example scripts for pulse shape computation and nested sampling.

* :doc:`Example job<Example_job>` contains example job scripts for computation on clusters.
