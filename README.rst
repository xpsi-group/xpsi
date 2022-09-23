.. _readme:

| Build Status Main | Docs |

X-PSI
=====

**An open-source package for neutron star**
**\ X-ray Pulse Simulation and Inference.**


X-PSI is designed to simulate rotationally-modified (pulsed) surface 
X-ray emission from neutron stars, taking into account relativistic 
effects on the emitted radiation. This can then be used to perform 
Bayesian statistical inference on real or simulated astronomical data 
sets. Model parameters of interest may include neutron star mass and 
radius (useful to constrain the properties of ultradense nuclear matter) 
or the system geometry and properties of the hot emitting surface-regions. 
To achieve this, X-PSI couples code for likelihood functionality (simulation) 
with existing open-source software for posterior sampling (inference).

It provides the following functionality:
* 
* 
* 
* 

For more details on current and planned capabilities, check out the 
`XPSI documentation <https://xpsi-group.github.io/xpsi/index.html>`_.

Installation and Testing
------------------------

X-PSI has a complex set of dependencies, and is therefore currently best 
installed from source. The documentation provides
`step-by-step installation instructions <https://xpsi-group.github.io/xpsi/install.html>`_
for Linux and for limited MacOS systems.

Documentation
-------------

The documentation for XPSI, including a wide range of tutorials and scripts for 
running XPSI on HPC systems, can be found at `https://xpsi-group.github.io/xpsi/ <https://xpsi-group.github.io/xpsi/>`_.

How to get in touch or get involved
-----------------------------------

We always welcome contributions and feedback! We are especially interested in 
hearing from you if
* something breaks,
* you spot bugs, 
* if there is missing functionality, or you have suggestions for future development,

To get in touch, please `open an issue <https://github.com/xpsi-group/xpsi/issues>`_.
Even better, if you have code you'd be interested in contributing, please send a 
`pull request <https://github.com/xpsi-group/xpsi/pulls>`_ (or get in touch 
and we'll help guide you through the process!). 

For more information, you can take a look at the documentation's 
`Contributing page <https://xpsi-group.github.io/xpsi/contributing.html>`_. Please also 
make sure you take a look at the `Code of Conduct <https://xpsi-group.github.io/xpsi/code_of_conduct.html>`_. 


Citing XPSI
-----------
If you find this package useful in your research, please provide the appropriate acknowledgment 
and citation. `Our documentation <https://xpsi-group.github.io/xpsi/index.html#citation>`_ provides 
more detail, including links to appropriate papers and BibTeX entries.

Copyright and Licensing
-----------------------
All content © 2016-2022 the authors. 
The code is distributed under the MIT license; see `LICENSE.rst <LICENSE.rst>`_ for details.

Legacy
------ 
An earlier version (pre-v0.5) of this project was named:
A prototype open-source package for neutron star X-ray Pulsation Simulation
and Inference.

.. |Build Status Main| image:: https://github.com/xpsi-group/xpsi/workflows/CI%20Tests/badge.svg
   :target: https://github.com/xpsi-group/xpsi/actions/
.. |Docs| image:: https://img.shields.io/badge/docs-latest-brightgreen.svg?style=flat
   :target: https://xpsi-group.github.io/xpsi/index.html

