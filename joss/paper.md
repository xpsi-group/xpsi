---
title: 'X-PSI: A Python package for neutron star X-ray pulse simulation and inference'
tags:
    - Python
    - astrostatistics
    - neutron stars
authors:
    - name: Thomas E. Riley
      orcid: 0000-0001-9313-0493
      affiliation: 1
    - name: Devarshi Choudhury
      orcid: 0000-0002-2651-5286
      affiliation: 1
    - name: Tuomo Salmi
      orcid: 0000-0001-6356-125X
      affiliation: 1
    - name: Serena Vinciguerra
      orcid: 0000-0003-3068-6974
      affiliation: 1
    - name: Yves Kini
      orcid: 0000-0002-0428-8430
      affiliation: 1
    - name: Bas Dorsman
      orcid:  0000-0002-9407-0733
      affiliation: 1
    - name: Anna L. Watts
      orcid: 0000-0002-1009-2354
      affiliation: 1
    - name: Daniela Huppenkothen
      orcid: 0000-0002-1169-7486
      affiliation: 2
    - name:  Sebastien Guillot
      orcid: 0000-0002-6449-106X
      affiliation: 3
affiliations:
   - name: Anton Pannekoek Institute for Astronomy, University of Amsterdam, Science Park 904, 1090GE Amsterdam, The Netherlands
     index: 1
   - name: SRON Netherlands Institute for Space Research, Niels Bohrweg 4, NL-2333 CA Leiden, the Netherlands
     index: 2
   - name: Institut de Recherche en Astrophysique et Planétologie, UPS-OMP, CNRS, CNES, 9 avenue du Colonel Roche, BP 44346, F-31028 Toulouse Cedex 4, France
     index: 3
date: 26 September 2022
bibliography: xpsijoss.bib
---


# Summary

The X-ray Pulse Simulation and Inference (X-PSI) package is a software package designed to simulate rotationally-modulated 
surface X-ray emission from neutron stars and to perform Bayesian 
statistical inference on real or simulated pulse profile data sets. Model parameters 
of interest include neutron star mass and radius and the system geometry and 
properties of the hot emitting surface regions. 

# Statement of need

Pulsed X-ray signals from neutron stars
can be modeled to statistically estimate parameters such as stellar mass and
radius, and properties of the surface radiation field such as a map of
temperature. The mass and radius of a neutron star are a function of the
equation of state of internal matter, especially the dense matter in the core,
and the formation history of the star, which determines the central energy
density and the spin frequency. The state of the surface radiation field is
the product of a potentially long and complex stellar evolutionary history,
especially that of the stellar magnetosphere. Such parameter estimation
requires relativistic tracing of radiation as it propagates from surface to a
distant telescope. Pulse-profile modelling to infer neutron star parameters
is a major science goal for both current X-ray telescopes such as the Neutron
Star Interior Composition ExploreR [NICER, @Gendreau2016] and proposed future telescopes such as eXTP and STROBE-X
[@Watts2019;@Ray2019].

While there are some open-source libraries for simulating the X-ray
signals from rapidly spinning neutron stars and more generally from the
vicinity of general relativistic compact objects including black holes [@Nattila:2016;@Pihajoki:2018] the scope
of these projects does not include statistical modeling, which
necessitates tractable parametrised models and a modular framework for
constructing those models.  X-PSI addresses this need, coupling code for likelihood
functionality (simulation) with existing open-source software for posterior sampling (inference).

# The X-PSI package and science use

X-PSI is an open-source Python package for Bayesian modeling of time- and
energy-resolved X-ray pulsations. X-PSI provides a framework for the
implementation of custom models, including likelihood and prior functions, and for
feeding those models to open-source statistical sampling software for use on
high-performance computing systems. X-PSI supplies modules for post-processing
posterior sample sets, and supplies tools for visualisation. For example, one
can generate time- and energy-resolved images and animations of model X-ray pulsars
by tracing radiation from the stellar surface to a distant observer, such as a
space telescope, as it propagates through spacetime; a snapshot of such an
animation may be found in \autoref{fig:animation_snapshot}. Posterior summaries
can be plotted of the time- and energy-domain signals that a model pulsar is
inferred to generate, conditioned on observational data.

![A snapshot from a time- and energy-resolved animation of a toy neutron star
that generates X-ray pulsations. The top three panels are sky maps of photon
specific intensity at three representative photon energies, and the bottom
panels show instantaneous integrals over the sky maps, together with sky map
integrals at past times normalized to the maximum at each photon energy
(bottom-left panel) and additional photon energies (bottom-right
panel).\label{fig:animation_snapshot}](fig1.png){width=100%}

X-PSI is being used by the NICER collaboration for pulse-profile modeling of X-ray emission from rotation-powered
millisecond pulsars [@Riley:2019;@Riley:2021].  Many more papers have used the accompanying
open-source analysis pipeline and products published on Zenodo [@Riley:2019:Zenodo;
@Riley:2021:Zenodo]. The first were @Raaijmakers:2019 and
@Bilous:2019, respectively on the topics of dense matter inference and
multipolar magnetic fields.

The numerical likelihood routines native to X-PSI are written in Cython
[@cython2011], and are dependent on the GNU Scientific Library [GSL,
@Gough:2009]. High-level object-oriented model construction is performed by a
user in the Python language, as is the interfacing with sampling software.
Low-level customisation is encouraged in the extensions, either directly in
Cython or via calls to external C libraries.  X-PSI is Unix source code
compatible, and release versions are freely available on GitHub under the GNU General Public License.  Extensive documentation, step-by-step tutorials, and reproduction
code for existing data analyses, are available
via the GitHub repository, along with a growing suite of unit tests.  Future plans
include migration to Python 3, further improvements to post-processing software,
 and the implementation of an expanded suite of atmosphere models.



*Software:* Python/C language [@Python2007], GNU Scientific Library [GSL,
@Gough:2009], NumPy [@Numpy2011], Cython [@cython2011], OpenMP [@openmp], MPI
for Python [@mpi4py], Matplotlib [@Hunter:2007; @matplotlibv2], IPython
[@IPython2007], Jupyter [@Kluyver:2016aa], MultiNest [@MultiNest_2009],
PyMultiNest [@PyMultiNest], GetDist [@Lewis19], nestcheck
[@higson2018nestcheck;@higson2018sampling;@higson2019diagnostic], fgivenx
[@fgivenx], emcee [@emcee]. 

# Acknowledgements

All University of Amsterdam co-authors acknowledge
support from ERC Consolidator grant No. 865768 AEONS (PI: ALW).  DH is supported by the 
Women In Science Excel (WISE) programme of the Netherlands Organisation for 
Scientific Research (NWO). SG acknowledges the support of the CNES. More detailed acknowledgements are written in the project
documentation [hosted on GitHub](https://xpsi-group.github.io/xpsi/acknowledgements.html).

# References
