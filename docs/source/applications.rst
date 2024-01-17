.. _applications:

Publications
------------

X-PSI has been applied in the following studies. In addition to giving an
idea of the scientific applications, these may also
be useful for performance benchmarking and planning 
of computational resource consumption. 

If you have used X-PSI for a project and would like to link it here, please
contact the X-PSI team and/or submit a pull-request on GitHub.


Early X-PSI development
***********************

X-PSI was initiated by Thomas Riley as part of his PhD work at the University of Amsterdam. 
Chapter 3 of his `PhD thesis (Riley 2019) <https://hdl.handle.net/11245.1/aa86fcf3-2437-4bc2-810e-cf9f30a98f7a>`_ 
provides an extended technical description of ``v0.1`` of X-PSI and contains
results of some parameter recovery tests using synthetic data.  


NICER papers
************

X-PSI has been used in several Pulse Profile Modeling analysis papers using *NICER* data of rotation-powered millisecond pulsars. These are typically published with a Zenodo repository containing data, model files, X-PSI scripts, samples and post-processing notebooks to enable full reproduction and independent analysis. 

**Vinciguerra et al. 2024** `(ApJ accepted) <https://ui.adsabs.harvard.edu/abs/2023arXiv230809469V/abstract>`_ used  ``v0.7.3`` to ``v2.0.0`` of X-PSI to carry out a much more in-depth and updated analysis of the PSR J0030+0451 *NICER* data set analysed in Riley et al. (2019). See also the associated `Zenodo repository`__.   

.. _Zenodo24a: https://doi.org/10.5281/zenodo.8239000
__ Zenodo24a_

**Vinciguerra et al. 2023** `(ApJ, 959, 55) <https://ui.adsabs.harvard.edu/abs/2023ApJ...959...55V/abstract>`_ used  ``v0.7.9`` to ``v2.0.0`` of X-PSI to carry out parameter recovery simulations for surface temperature maps inspired by PSR J0030+0451.  See also the associated `Zenodo repository`__.

.. _Zenodo23b: https://doi.org/10.5281/zenodo.7646352
__ Zenodo23b_

**Salmi et al. 2023** `(ApJ, 956, 138) <https://ui.adsabs.harvard.edu/abs/2023ApJ...956..138S/abstract>`_ used  ``v0.7.3`` to ``v1.2.1`` of X-PSI to investigate astmospheric effects on neutron star parameter contraints derived from *NICER* observations.  See also the associated `Zenodo repository`__.

.. _Zenodo23a: https://doi.org/10.5281/zenodo.7449785
__ Zenodo23a_

**Salmi et al. 2022** `(ApJ, 941, 150) <https://ui.adsabs.harvard.edu/abs/2022ApJ...941..150S/abstract>`_ used  ``v0.7.10`` of X-PSI to model *NICER* observations of the rotation-powered millisecond X-ray pulsar PSR J0740+6620 using *NICER* background estimates.  See also the associated `Zenodo repository`__.

.. _Zenodo22: https://doi.org/10.5281/zenodo.6827536
__ Zenodo22_


**Riley et al. 2021**  `(ApJL, 918, L27) <https://ui.adsabs.harvard.edu/abs/2021ApJ...918L..27R/abstract>`_ used ``v0.7.6`` of X-PSI to model *NICER* observations of the rotation-powered millisecond X-ray pulsar PSR J0740+6620. See also the associated `Zenodo repository`__.

.. _Zenodo21: https://doi.org/10.5281/zenodo.4697624
__ Zenodo21_

**Bogdanov et al. 2021**  `(ApJL, 914, L15) <https://ui.adsabs.harvard.edu/abs/2021ApJ...914L..15B/abstract>`_ provides additional details of the models used for *NICER* Pulse Profile Modeling, reports ray-tracing cross-tests for X-PSI and other codes for very compact stars where multiple imaging is important, and reports some parameter recovery simulations for synthetic data.  

**Bogdanov et al. 2019** `(ApJL, 887, L26) <https://ui.adsabs.harvard.edu/abs/2019ApJ...887L..26B/abstract>`_ reports the results of ray-tracing cross-tests for several codes in use in the *NICER* team including X-PSI.

**Riley et al. 2019** `(ApJL, 887, L21) <https://ui.adsabs.harvard.edu/abs/2019ApJ...887L..21R/abstract>`_ used 
``v0.1`` of X-PSI to model *NICER* observations of the rotation-powered millisecond X-ray pulsar PSR J0030+0451. See also the associated `Zenodo repository`__.

.. _Zenodo: https://doi.org/10.5281/zenodo.3386448

__ Zenodo_


Other papers
************

**Kini et al. 2024** `(MNRAS, 527, 8118) <https://academic.oup.com/mnras/article/527/3/8118/7440002>`_ used ``v0.7.9`` of X-PSI (with small modifications described in the paper) to explore strategies for mitigating variability when doing pulse profile modelling for thermonuclear burst oscillations.  See also the associated `Zenodo repository`__.

.. _Zenodo24kini: http://dx.doi.org/10.5281/zenodo.8033527
__ Zenodo24kini_

**Kini et al. 2023** `(MNRAS, 522, 3389) <https://ui.adsabs.harvard.edu/abs/2023MNRAS.522.3389K/abstract>`_ used  ``v0.7.9`` of X-PSI (with small modifications described in the paper) to model simulated *RXTE* observations of thermonuclear X-ray burst oscillations when ignoring time variability in the properties of the emitting regions.  See also the associated `Zenodo repository`__.

.. _Zenodo23kini: http://dx.doi.org/10.5281/zenodo.7665653
__ Zenodo23kini_


Settings
********

Below we give details of some settings for the **Riley et al. 2019** paper.  For later papers some changes
were made - see papers for details. 

*Resource consumption*:  The calculations reported were performed on the Dutch national supercomputer
Cartesius, mostly on the *Broadwell* nodes (i.e., using CPUs only, usually
*Broadwell* architecture).
A number of models were applied, with increasing complexity.
Typically, ~1000+ cores were used for posterior sampling, when the
parallelisation is weighted by time consumed.
In total, the sampling problems in R19 required ~450,000 core hours.

*Likelihood function settings*:  On compute nodes, the likelihood function evaluation times ranged from ~1 to
~3 seconds (single-threaded), depending on the model and point in parameter
space.\ [#]_ The parallelisation was purely via MPI, with all but one process
receiving likelihood function evaluation requests.

*Resolution settings*: Number of surface cells/elements per hot (closed) region:\ [#]_ 24x24; 
Number of phases (linearly spaced) = 100; Number of energies (logarithmically spaced) = 175; 
Number of rays (per parallel; linear in cosine of ray angle alpha):\ [#]_ 200

*Interpolation settings*:  Steffen splines (GSL) were used everywhere apart from for the atmosphere, and 
for the atmosphere cubic polynomial interpolants were used in four dimensions.


*Compilation*:  Intel compiler collection,\ [#]_ with the CORE-AVX2 instruction set, for X-PSI
and dependencies (apart from :mod:`numpy`, which was centrally installed).



.. rubric:: Footnotes

.. [#] Variation between models on the same processor results from a
       combination of model complexity and resolution settings. Variation
       as a function of parameters also occurs due to issues like mesh
       construction on the surface subject to the resolution settings.

.. [#] This is the approximate number of cells that have centres lying
       within a hot region. For models where a hot region has two temperature
       components these cells are split between the subregions according to
       area only such that the surface density of cells is commensurate---at
       least for subregions that are commensurate in shape. The *fast*
       precomputation mode was deactivated but if it were activated the
       distribution of cells would further be weighted by the approximate
       total number of counts generated by each subregion.

.. [#] The rays were integrated for each likelihood function call instead of
       loading lookup tables from disk, or using an analytic treatment such as
       a high-order expansion (see, e.g.,
       `rayXpanda <https://github.com/ThomasEdwardRiley/rayXpanda>`_).
       Only primary images were included.

.. [#] On Cartesius (and Lisa), one can simply execute
       ``module load intel/2017b`` to access these compilers from the toolchain.


