.. _landing_page_tutorials:

==========
Start here
==========

All tutorials listed below may be found as Jupyter notebooks under ``docs/source/``. We encourage running these yourself as you go through the tutorials as a means of practice. Any external data files needed such as instrument response files and numerical atmosphere data are available on `Zenodo <https://doi.org/10.5281/zenodo.7094145>`_.

**Overview of the available tutorials:**

* :doc:`Modeling<Modeling>` is an in-depth tutorial that first covers the structure of X-PSI followed by basic usage of X-PSI. Basic usage includes constructing your first model (e.g. instrument, star, and atmosphere) and construct a likelihood for nested sampling.

* :doc:`Instrument synergy<Instrument synergy>` shows how to construct a joint likelihood with data from two instruments.

* :doc:`Hot region complexity<Hot region complexity>` goes into more detail about multiple and especially multiple overlapping hot regions. We note that :doc:`hotregion<hotregion>` is a useful page that contains overview figures of various overlapping cases.

* :doc:`Global surface emission<Global surface emission>` concerns surface emission from fields that span the full star and showcases three different options for signal integration. The first is a globally uniform temperature field, which allows for a time invariant signal integration. In the second case the temperature field is phase dependent and thus requires time dependent integration. The third case is general purpose and discretises a distant image plane instead of the stellar surface.

* :doc:`Surface radiation field tools<Surface radiation field tools>` demonstrates the usage of the default (blackbody) and an alternative (atmosphere interpolated from precomputed data) surface radiation field module to compute photon specific intensities. It also shows how to do beaming pattern and spectrum plots for the radiation fields.

* :doc:`Modeling (without statistics)<Modeling (without statistics)>` is similar to:doc:`Modeling<Modeling>` but omits any statistical inference and adds various plots of signals. This is a useful tutorial if you are only interested to use X-PSI to create synthetic data.

* :doc:`Multiple imaging<Multiple imaging>` is an exercise in image validation by using different integrators (largely following :doc:`Global surface emission<Global surface emission>`) and additionally using another integrator module: rayXpanda, which is capable of computing higher order lensing. In the process of this comparison, this is also provides examples of generating multiple image figures.

* :doc:`Post-processing<Post-processing>` is under construction.

* :doc:`Example script and modules<Example script and modules>` is an (advanced) example script for nested sampling from :ref:`R19`. It employs custom photosphere, instrument, interstellar, signal, and prior.

* :doc:`Example job<Example job>` contains and example job script and resume job script for computation on a cluster, also from :ref:`R19`.


