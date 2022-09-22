.. _faq:

FAQs and common problems
========================


Installation
^^^^^^^^^^^^

.. rubric:: Do I need to edit the package setup script?

You may well have to edit the setup script depending on the target system.
This includes editing compiler flags (see below for example regarding
instruction sets).

.. rubric:: Does it matter what compiler I use?

The Intel compiler collection has been used successfully for X-PSI and
dependencies (namely GSL, MultiNest). We recommend first trying to use Intel
in a context where performance matters.

.. rubric:: What Intel instruction sets should I use?

If you want to test the binaries on a login node, note that you can
compile with multiple instruction sets for auto-dispatch using the ``-x`` and
``-ax`` flags. See the :ref:`hpcsystems` page for examples.

.. rubric:: What atmosphere extension should I use?

If you want to do some quick calculations and run the Modeling tutorial, you should use the default blackbody atmosphere extension. But if you want to use similar atmosphere models as in the published X-PSI applications so far, you should switch to the numerical atmosphere extension before installing X-PSI (see instructions in :ref:`install`).

Note that using the model scripts with an unintended atmosphere extension may lead to a segmentation fault error if trying to use numerical atmosphere without providing numerical atmosphere data. Or if using the scripts including the numerical atmosphere data with the blackbody extension, you can get unexpected results, printing though a warning that numerical atmosphere data were preloaded but not used. Examples of numerical atmosphere data, which are required by the numerical atmosphere extension, can be found e.g. in the Zenodo repository of Riley et al. 2021: `doi:10.5281/zenodo.4697625`__. Examples of how to use the numerical atmospheres are shown e.g. in Surface radiation field tools -tutorial and in :ref:`example_script`.

.. _Zenodo: https://zenodo.org/record/4697625

__ Zenodo_


Sampling
^^^^^^^^

.. rubric:: Is I/O or disk storage a concern, or are all the files small?

I/O not a concern for likelihood calculation.

Nested sampling writes to disk at user-specified cadence
(so many nested sampling iterations).

Model data such as a four-dimensional atmosphere table can be reasonably
large for I/O.
We recommend loading, at the outset of the run (or a resumed run),
such a table into a contiguous chunk of memory
for each of the Python processes running on one node.
That table is pointed to for access where needed from compiled modules
(C extensions to Python): it is not loaded from disk per likelihood call.
We provide an example custom Python class that handles this loading (as used
in :ref:`R19`, hereafter R19).

Disk storage required is indeed small: up to :math:`\mathcal{O}(100)` Mbytes for
applications thus far (e.g., R19). There is a variant of MultiNest nested sampling
that is much more memory and disk intensive, but we do not use it.  This is
because importance nested sampling is not compatible with the alternative options
(read: hacks) for prior implementation (see `Riley, PhD thesis <https://hdl.handle.net/11245.1/aa86fcf3-2437-4bc2-810e-cf9f30a98f7a>`_).


Common problems and errors
^^^^^^^^^^^^^^^^^^^^^^^^^^

- To avoid errors in post-processing, do not use standard PostProcessing scripts 
  for runs which have not converged, and when you do the Postprocessing:
   - make sure to change the paths in the config.ini compared to the ones used in the cluster;
   - make sure to change the cache to cache = True in the signal definition in the main;
   - to have the likelihood check, to be use you are using the same x-psi configuration as in the run.

-  | *AttributeError: ’NestedBackend’ object has no attribute
     ’\ :math:`\_nc\_bcknd`\ ’* in PostProcessing for runs with
     use\ :math:`\_`\ nestcheck=[False] (e.g. importance sampling).
   | SOLUTION: turn bootstrap\ :math:`\_`\ estimators=False,
     or alternatively, set use\ :math:`\_`\ nestcheck=[True].

- Skymap animation or average; if you see many annular images, the
  problem is in the local\ :math:`\_`\ variable.pyx file (which should
  have been subsituted to the archive/local/ :math:`\_`\ variables/PST\ :math:`\_`\ U.pyx)

   .. container:: figure*

      .. image:: images/ST_PST__NICER__skymap_phase_averaged_run1.png
         :alt: image
         :width: 10cm


-  | *Traceback (most recent call last):*
   | *File " :math:`<stdin>`", line 1, in :math:`<module>`*
   | *File "xpsi/:math:`_{\_\_}`\ init\ :math:`_{\_\_}`.py", line 135, in <module>*
   | *from tools import set_phase_interpolant set_energy_interpolant*
   | *ImportError: No module named tools*
   | SOLUTION: You are running X-PSI from its main directory ( the
     directory where the* **setup.py** *file is). Exit that directory and
     run it again.

-  | *transform() got an unexpected keyword argument ’\ :math:`old\_API'`*
   | SOLUTION: Double check the **names** and **bounds** if everything
     correct then add **\**kwargs** to the **transform** method in your
     CustomPrior Class

-  | *<path/to/run/output>dead-birth.txt not found.
   | SOLUTION: Set use_nestcheck to* **False** *(
     use_nestcheck=[*\ **False**\ *]*

-  | *ImportError: libgsl.so.23: cannot open shared object file: No such
     file or directory* after while trying to import X-PSI
   | SOLUTION: Somthing has gone wrong but no need to reinstall
     everything ( gsl, multinest and so on). Just clean everything ( rm
     -r build dist \*egg\* xpsi/*/*.c) and reinstall only X-PSI.

-  | *Invalid caching targets.*
   | SOLUTION 1: Turn cache in Signal in main.py to **True**, if it isn’t
     already; call the likelihood or even better the likelihood check.
   | SOLUTION 2: Set **STU.signal.cache = True** or
     **STU.signals[0][i].cache = True** (where i is the index of the
     instrument e.g. i=0 for NICER in a joint NICER+XMM analysis) in your
     postprocessing script (before signal plotting) and then re-run the
     script (then no modifications to the main.py needed). Replace ’STU’
     with the name of the imported model, if using different model name.
     And remember to call the likelihood check also in this case (in the
     postprocessing script).

-  | *Each row and column must contain at least one positive number.*
   | PROBLEMS: there are some rows and/or column in the instrument
     response that contain only zeros. 
   | SOLUTION: increase the number of
     channels or decrease the number of energy intervals.

-  | *kwargs["sampling_efficiency" = self.\ :math:`\
     _`\ prior.unit_hypercube_frac TypeError: unsupported operand
     type(s) for /=: "float" and "NoneType"*
   | PROBLEM: Likely something need to be corrected in your CustomPrior
     or in the way you are initiating or using it in your main. You can
     try to debug and find the actual error by removing the error
     handling in the end of Prior.py (where printing now just ’Cannot
     locate method for estimating fraction.’) and re-installing X-PSI.

-  | **post-processing**: PROBLEM: Invalid index (IndexError or
     JoblibValueError):
   | SOLUTIONS turn KL_divergence to False. Seems to be a problem for
     the circular parameters, i.e. the phases. In this case it is always
     possible to take the phases out.

-  | **post-processing**: PROBLEM ERROR when using KL_divergence to
     True. Errors can arise because the bounds of one or more
     parameter(s) are None.
   | SOLUTION: explicitly write bound values for all parameters
     (geometrical strict bounds can be found in HotRegion.py);
     alternatively put KL_divergence to False.

-  | **post-processing**: PROBLEM **Warning: Using native nestcheck KDE
     instead of GetDist KDE.** And possibly errors later in the
     postprocessing.
   | SOLUTION: Make sure to to install nestcheck and GetDist packages
     using the github repositories (as instructed in the footnotes of
     the installation tutorial) and not using pip.

Common errors and solution NICER and XMM:
=========================================

-  | Error message ``from STU_NICER_XMM import main as STU``:
   | *main.py in <module>()*
   | *–> 420 prior = prior)*
   | *103 energies =*
   | *construct\ :math:`\_`\ energy\ :math:`\_`\ array(num\ :math:`\_`\ energies,*
   | *–> 104 list(signals))*
   | *–> 647 MAX = ordered[0].energy\ :math:`\_`\ edges[-1]*
   | *IndexError: list index out of range*
   | SOLUTION:
     ``args = parser.parse``\ :math:`\_`\ ``args([’@STU_NICER``\ :math:`\_`\ ``XMM/config.ini’,’–NICER’,’–XMM’])``

-  | **PostProcessing:** *—> 13 plots = ’ST-U’: xpsi.ResidualPlot())*
   | *–> 107 output = func(args, kwargs)*
   | *52 if :math:`\_`\ isgeneratorfunction(func):*
   | *—> 53 for msg in func(args, kwargs):*
   | *54 if :math:`\_`\ verbose and not*
   | *deactivate\ :math:`\_`\ verbosity:*
   | *–> 162 likelihood.signal.caching\ :math:`\_`\ targets =*
   | *caching\ :math:`\_`\ targets*
   | *188 if len(self.\ :math:`\_`\ signals) or*
   | *len(self.\ :math:`\_`\ signals[0]) :*
   | *–> 189 raise ValueError(’There is more than one signal instance.’)*
   | *190 return self.\ :math:`\_`\ signals[0][0]*
   | *ValueError: There is more than one signal instance.*
   | SOLUTION: ``STU.likelihood.signals = STU.likelihood.signals[0][0]``

-  | **KL-divergence = NaN**: KDE on prior and posterior samples results 
     inaccurately in zero prior density where there is a posterior sample. 
     To avoid, need more prior samples  Alternatively, you may need to check that the bounds given in
     the postprocessing notebook/script correspond to the actual hard bounds, 
     to make sure that no posterior samples are out of the bounds
     (defined in the postprocessing).

-  | **the role of bootstrap\ :math:`\_`\ estimators** The 1D KL
     divergence has an error, and the 1D credible interval is based on
     distribution of credible intervals from bootstrapped realisations.
     The ends of the credible region will have thin darker vertical
     bands indicating their error.

-  **kernel or terminal dies:** check that the installed hot.pyx is
   consistent with the atmosphere model used in the adopted xpsi model.


