.. _readme:

********************************************
X-ray Pulse Simulation and Inference (X-PSI)
********************************************

**An open-source package for neutron star**
**\ X-ray Pulse Simulation and Inference.** 

X-PSI is designed to simulate surface X-ray emission from rotating neutron stars, particularly for cases where the surface emission is not uniform, and the resulting emission can therefore be pulsed. It takes into account the effects of relativity on the emitted radiation. This can then be used to perform Bayesian statistical inference on real or simulated astronomical data sets. Model parameters of
interest may include neutron star mass and radius (useful to constrain the
properties of ultradense nuclear matter) or the system geometry and properties
of the hot emitting regions on the neutron star surface. To achieve this, X-PSI couples code for likelihood functionality (simulation) with open-source software for
posterior sampling (inference).

.. image:: _static/pulse_profile_example.png

X-PSI has been used extensively to model pulse profiles like the one shown above.  Pulse profiles indicate how the emission (flux and spectrum) of neutron stars varies as a function of their rotational phase.  X-PSI has been used extensively for Pulse Profile Modeling (PPM) of rotation-powered millisecond X-ray pulsars, in particular using data from NASA's `Neutron Star Interior Composition Explorer (NICER) <https://www.nasa.gov/nicer>`_.  It can also be used to model pulse profile data for accreting neutron stars, including accretion-powered millisecond X-ray pulsars and thermonuclear burst oscillation sources.  It has the facility to model polarized X-ray emission, and can also be used to simulate and model phase-averaged X-ray spectra.  For more details see :ref:`applications`.      

Getting started
*****************

If you are new to X-PSI, you'll find full :ref:`install` instructions here, and to work through the tutorials we suggest you :ref:`start here<landing_page_tutorials>`.

