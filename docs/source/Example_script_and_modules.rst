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

Third script called `examples/examples_modeling_tutorial/TestRun_Num.py`__ is also a test run, which is similar to the :doc:`Modeling<Modeling>` tutorial, except using the numerical atmosphere extension (X-PSI installed using the flag ``--NumHot``).

.. _t3: https://github.com/xpsi-group/xpsi/tree/main/examples/examples_modeling_tutorial/TestRun_Num.py

__ t3_

Fourth script called `examples/examples_modeling_tutorial/TestRun_NumBeam.py`__ is a similar test run as the previous, except using a numerical atmosphere extension with extra parameters allowing deviations to the atmospheric beaming pattern (X-PSI installed using the flag ``--NumHotBeam``).

.. _t4: https://github.com/xpsi-group/xpsi/tree/main/examples/examples_modeling_tutorial/TestRun_NumBeam.py

__ t4_
