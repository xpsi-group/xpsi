.. _example_script:

Example script and modules
==========================

A few example scripts and modules are found here: `examples`__.

.. _examples: https://github.com/xpsi-group/xpsi/tree/main/examples

__ examples_

The script called `examples/examples_fast/Modules/main.py`__ is a quick test run not requiring any input files. It assumes a single hot region with a uniform temperature using the default blackbody atmosphere extension and a fake instrument matrix.

.. _t1: https://github.com/xpsi-group/xpsi/tree/main/examples/examples_fast/Modules/main.py

__ t1_

Second script called `examples/examples_modeling_tutorial/TestRun_BB.py`__ is a test run, which is similar to the :doc:`Modeling<Modeling>` tutorial, requiring the same input files as in the corresponding tutorial and using the default blackbody atmosphere extension.

.. _t2: https://github.com/xpsi-group/xpsi/tree/main/examples/examples_modeling_tutorial/TestRun_BB.py

__ t2_

Third script called `examples/examples_modeling_tutorial/TestRun_Num.py`__ is also a test run, which is similar to the :doc:`Modeling<Modeling>` tutorial, except using the numerical atmosphere extension.

.. _t3: https://github.com/xpsi-group/xpsi/tree/main/examples/examples_modeling_tutorial/TestRun_Num.py

__ t3_

Fourth script called `examples/examples_modeling_tutorial/TestRun_NumBeam.py`__ is a similar test run as the previous, except using a numerical atmosphere extension with extra parameters allowing deviations to the atmospheric beaming pattern.

.. _t4: https://github.com/xpsi-group/xpsi/tree/main/examples/examples_modeling_tutorial/TestRun_NumBeam.py

__ t4_

Fifth script called `examples/examples_modeling_tutorial/TestRun_Pol.py`__ is an example script for producing polarized pulses in X-PSI using a simple analytical model for the atmosphere emission.

.. _t5: https://github.com/xpsi-group/xpsi/tree/main/examples/examples_modeling_tutorial/TestRun_Pol.py

__ t5_

Sixth script called `examples/examples_modeling_tutorial/TestRun_PolNum_split_1spot.py`__ is an example script for producing polarized pulses in X-PSI using a numerical atmosphere model with 3+2 dimensional interpolation.

.. _t6: https://github.com/xpsi-group/xpsi/tree/main/examples/examples_modeling_tutorial/TestRun_PolNum_split_1spot.py

__ t6_

Seventh script called `examples/examples_modeling_tutorial/TestRun_PolNum_split_inference.py`__ is a preliminary example script for producing and analyzing polarized pulses in X-PSI using a numerical atmosphere model with 3+2 dimensional interpolation and ST+PDT spot model.

.. _t7: https://github.com/xpsi-group/xpsi/tree/main/examples/examples_modeling_tutorial/TestRun_PolNum_split_inference.py

__ t7_

Eigth script called `examples/examples_modeling_tutorial/TestRun_AMXP.py`__ is an example script of an Accreting Millisecond Pulsar with accretion disk implemented through the CustomPhotosphere.py, and 3+2 dimensional interpolation. Note that TestRun_PolNum_split_1spot similarly includes an accretion disk implementation at the end, but there it is a background signal. The results are equivalent but the implementation in this example here is faster and thus recommended for sampling.

.. _t8: https://github.com/xpsi-group/xpsi/tree/main/examples/examples_modeling_tutorial/TestRun_AMXP.py

__ t8_

