---
title: 'X-PSI: A Python package for neutron star X-ray Pulse Simulation and Inference'
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
      orcid: 0000 0003 3068 6974
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
date: 23 September 2022
bibliography: xpsijoss.bib
---

# Summary

Stars play host to exotic environments that cannot be simulated in terrestrial
laboratories. The focus of this work is neutron stars, thought to be the most
compact extended objects in the Universe. Observable radiation from a neutron
star encodes information about fundamental physics (gravity, electromagnetism,
and nuclear forces) and astrophysical processes (such as the state and evolution
of a stellar magnetosphere). Neutron stars are often detected by astronomers
because their radiative signal - in radio, X-ray and gamma-rays -  is pulsed due to bulk stellar rotation.  By modeling the physical process that generates data registered by telescopes, astronomers and astrostatisticians can make inferential statements
about the nature of the extreme Universe.


# Statement of need

There exist open-source libraries and packages to support a subset of the
modeling treated in astrophysical literature. They provide frameworks,
toolsets, model implementations, and so on. One sub-field for which there does
not exist such an open-source project, is for the *statistical modeling* of
X-ray signals that pulse due to rotational modulation of asymmetric emission
from the surface of a neutron star. Pulsing X-ray signals from neutron stars
are modeled to statistically estimate parameters such as stellar mass and
radius, and properties of the surface radiation field such as a map of
temperature. The mass and radius of a neutron star are a function of the
equation of state of internal matter (especially the dense matter in the core)
and the formation history of the star (that determines the central energy
density and the spin frequency); the state of the surface radiation field is
the product of a potentially long and complex stellar evolutionary history,
especially that of the stellar magnetosphere. Such parameter estimation
requires relativistic tracing of radiation as it propagates from surface to a
distant telescope. Pulse-profile modelling to infer neutron star parameters
is a major science goal for both current X-ray telescopes such as the Neutron
Star Interior Composition ExploreR (NICER, @Gendreau2016) and future telescopes
@Watts2019.

There are a small number of open-source libraries for simulating the X-ray
signals from rapidly spinning neutron stars and more generally from the
vicinity of general relativistic compact objects (including black holes).  
An example is the Arcmancer
library of @Pihajoki:2018, a general purpose toolbox that is wrapped by an
updated version of the bender library of @Nattila:2016 for the purpose of
simulating X-ray signals from hot regions on the surfaces of rapidly rotating
neutron stars.   However the scope of these projects does not include statistical modeling, which
necessitates tractable parametrised models and a modular framework for
constructing those models.  X-PSI provides this functionality. 

#The X-PSI package and science use

X-PSI is an open-source Python package for Bayesian modeling of time- and
energy-resolved X-ray pulsations. X-PSI provides a framework for the
implementation of custom models (likelihood and prior functions) and for
feeding those models to open-source statistical sampling software for use on
high-performance computing systems. X-PSI supplies modules for post-processing
posterior sample sets, and supplies tools for visualisation. For example, one
can generate time- and energy-resolved (animated) images of model X-ray pulsars
by tracing radiation from the stellar surface to a distant observer (such as a
space telescope) as it propagates through spacetime; a snapshot of such an
animation may be found in \autoref{fig:animation_snapshot}. Posterior summaries
can be plotted of the time- and energy-domain signals that a model pulsar is
inferred to generate, conditioned on observational data.

![A snapshot from a time- and energy-resolved animation of a toy neutron star
that generates X-ray pulsations. The top three panels are sky maps of photon
specific intensity at three representative photon energies, and the bottom
panels show instantaneous integrals over the sky maps, together with sky map
integrals at past times normalized to the maximum at each photon energy
(bottom-left panel) and additional photon energies (bottom-right
panel).\label{fig:animation_snapshot}](_skymap_with_pulse_profile_and_spectrum_plot.png){width=100%}

X-PSI is being used by the NICER collaboration for pulse-profile modeling of X-ray emission from rotation-powered
millisecond pulsars. Publications that have directly applied the X-PSI package thus far are
@Riley:2019 and @Riley:2021. Many more papers have used the accompanying
open-source analysis pipeline and products published on Zenodo; these products
may be accessed using the linked DOIs of @Riley:2019:Zenodo and
@Riley:2021:Zenodo. The first were @Raaijmakers:2019 and
@Bilous:2019, respectively on the topics of dense matter inference and
multipolar magnetic fields.

The numerical likelihood routines native to X-PSI are written in Cython
[@cython2011], and are dependent on the GNU Scientific Library (GSL;
@Gough:2009). High-level object-oriented model construction is performed by a
user in the Python language, as is the interfacing with sampling software.
Low-level customisation is encouraged in the extensions (either directly in
Cython or via calls to external C libraries).  X-PSI is Unix source code
compatible, and release versions are freely available on GitHub under the MIT
license.  Extensive documentation, step-by-step tutorials, and reproduction
code for existing data analyses, are available
via the project github, along with a growing suite of unit tests.  Future plans
include migration to Python 3, further improvements to post-processing software,
 and the implementation of an expanded suite of atmosphere models.



*Software:* Python/C language [@Python2007], GNU Scientific Library (GSL;
@Gough:2009), NumPy [@Numpy2011], Cython [@cython2011], OpenMP [@openmp], MPI
for Python [@mpi4py], Matplotlib [@Hunter:2007; @matplotlibv2], IPython
[@IPython2007], Jupyter [@Kluyver:2016aa], MultiNest [@MultiNest_2009],
PyMultiNest [@PyMultiNest], GetDist [@Lewis19], nestcheck
[@higson2018nestcheck;@higson2018sampling;@higson2019diagnostic], fgivenx
[@fgivenx].

# Acknowledgements

All University of Amsterdam co-authors acknowledge
support from ERC Consolidator grant No. 865768 AEONS (PI: ALW).  DH is supported by the 
Women In Science Excel (WISE) programme of the Netherlands Organisation for 
Scientific Research (NWO). SG acknowledges the support of the CNES. More detailed acknowledgements are written in the project
documentation [hosted](https://xpsi-group.github.io/xpsi/acknowledgements.html)
on GitHub.

# References
