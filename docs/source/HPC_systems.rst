.. _hpcsystems:

HPC systems
================

The information provided in this page is for users who intend to work on 
High-Performance Computing (HPC) systems. These installation instructions are 
system-specific. X-PSI has already been used on different systems, for some of
which we provide the instructions below. This information may also be
translated to other systems by users looking for guidance.


Snellius (SURF)
-------------------

`Snellius <https://servicedesk.surf.nl/wiki/display/WIKI/Snellius>`_ is the 
Dutch National Supercomputer.

Installation
^^^^^^^^^^^^

All of the following must be performed on a login node, in your ``$HOME`` file
system.

Start by cleaning your home file system of existing versions of dependencies
and move anything else to some archive in ``$HOME``. Clean ``.bashrc`` and
``.bash_profile`` of anything related to this software. Clean ``.bashrc`` and
``.bash_profile`` of environment variables such as: ``LD_LIBRARY_PATH``,
``LD_PRELOAD``, ``RUN_PATH``, ``PATH``, and ``PYTHONPATH``. Then logout and
log back in order to get a clean environment.

To be additionally safe, run:

.. code-block:: bash

    module purge

Load environment module and modify clean environment with foss toolchain
information and the needed modules:

.. code-block:: bash

    module load 2022
    module load foss/2022a
    module load SciPy-bundle/2022.05-foss-2022a
    module load wrapt/1.15.0-foss-2022a
    module load matplotlib/3.5.2-foss-2022a
    module load CMake/3.23.1-GCCcore-11.3.0

Prepare a new Python virtual environment for X-PSI (named for example "xpsi_py3") in case the possibility of having several co-existing X-PSI and/or PyMultiNest versions is wished (otherwise proceed to MultiNest installation):

.. code-block:: bash

    mkdir venvs
    python -m venv ./venvs/xpsi_py3

To access all the loaded site packages when activating the virtual environment, one needs to modify the file ``./venvs/xpsi_py3/pyvenv.cfg`` (using e.g. ``vim`` or ``emacs`` text editor) to change "false" into "true":

.. code-block:: bash

    Include system site packages = true

Now the environment can be activated with

.. code-block:: bash

    source ./venvs/xpsi_py3/bin/activate


To prepare `MultiNest <https://github.com/farhanferoz/MultiNest>`_ from
``$HOME``:

.. code-block:: bash

    git clone https://github.com/farhanferoz/MultiNest.git ~/multinest
    cd ~/multinest/MultiNest_v3.12_CMake/multinest
    mkdir build; cd build
    cmake -DCMAKE_{C,CXX}_FLAGS="-O3 -march=znver2 -funroll-loops" -DCMAKE_Fortran_FLAGS="-O3 -march=znver2 -funroll-loops" ..; make
    ls ../lib/

Use the last command to check for the presence of shared objects.

We also need to set the environment variable for library path to point at
MultiNest:

.. code-block:: bash

    export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$HOME/multinest/MultiNest_v3.12_CMake/multinest/lib/

Now you need the Python interface to MultiNest, starting from ``$HOME``:

.. code-block:: bash

    git clone https://github.com/JohannesBuchner/PyMultiNest.git ~/pymultinest
    cd ~/pymultinest
    python setup.py install

.. note::

    If not using a Python virtual environment, you should add ``--user`` flag when installing PyMultiNest.

To test the installation of MultiNest and PyMultiNest on the login node:

.. code-block:: bash

    mpiexec -n 2 python pymultinest_demo.py

Do you obtain parameter values and evidences?

.. note::

    Without knowing exactly the reason, we currently get typically this message ``Open MPI failed an OFI Libfabric library call (fi_domain).  This is highly unusual; your job may behave unpredictably (and/or abort) after this.`` when doing this test. However, the test works otherwise as expected, and this message seem not to appear when submitting jobs in the cluster instead of using the login node.

.. note::

    We assumed above that nested sampling with `MultiNest`_ is desired. If
    ensemble-MCMC with ``emcee`` is desired, you need to install the Python
    packages ``emcee`` and ``schwimmbad``. We assume the user can infer how to
    do this using the information above and on the :ref:`install` page.

For `GSL <https://www.gnu.org/software/gsl/>`_ we can use the default 2.5
version already provided in Snellius. Thus, to prepare X-PSI from ``$HOME``, we
only need:

.. code-block:: bash

    git clone https://github.com/xpsi-group/xpsi.git
    cd ~/xpsi
    LDSHARED="gcc -shared" CC=gcc python setup.py install

.. note::

    If not using a Python virtual environment, you should add ``--user`` flag when installing X-PSI.

If you ever need to reinstall, first clean to recompile C files:

.. code-block:: bash

    rm -r build dist *egg* xpsi/*/*.c xpsi/include/rayXpanda/*.o

.. note::

    We typically do not use the :mod:`~xpsi.PostProcessing` module, but
    instead ``rsync`` output files to a local system to perform plotting. This
    circumvents any potential backend problems and permits straightforward use
    of IPython for interactive plotting. However, if one wishes to use it on an
    HPC, it would require the installation of `GetDist` and `Nestcheck`. See
    :ref:`install` page for relevant details.


Batch usage
^^^^^^^^^^^

For an example job script, refer to :ref:`example_job`.

Lisa (SURF)
-----------

`Lisa <https://servicedesk.surf.nl/wiki/display/WIKI/Lisa>`_ follows mostly the installation instructions as that of Snellius. Small differences in the installation procedure are still to be studied.

Helios (API)
------------

Helios is a cluster of the Anton Pannekoek Institute for Astronomy. 

Installation
^^^^^^^^^^^^

Let's start by loading the necessary modules and creating a conda environment. At the moment, the installation is known to be working only for the specific python 3.10.6 version, and when conda installing the required python packages separately, as followed:

.. code-block:: bash

   module load anaconda3/2021-05
   module load openmpi/3.1.6
   git clone https://github.com/xpsi-group/xpsi.git
   cd xpsi
   conda create -n xpsi_py3 python=3.10.6
   conda activate xpsi_py3
   conda install -c conda-forge mpi4py
   conda install cython~=0.29
   conda install scipy
   conda install matplotlib
   conda install wrapt   
     
Let's then test if mpi4py works:

.. code-block:: bash

   cd; wget https://github.com/mpi4py/mpi4py/releases/download/3.1.5/mpi4py-3.1.5.tar.gz
   tar zxvf mpi4py-3.1.5.tar.gz
   cd mpi4py-3.1.5
   mpiexec -n 4 python demo/helloworld.py
   
Let's then install MultiNest and PyMultiNest:
   
.. code-block:: bash
   
   cd; git clone https://github.com/farhanferoz/MultiNest.git multinest
   cd multinest/MultiNest_v3.12_CMake/multinest
   mkdir build
   cd build
   CC=gcc FC=/zfs/helios/filer0/sw-astro/api/openmpi/3.1.6/bin/mpif90 CXX=g++ cmake -DCMAKE_{C,CXX}_FLAGS="-O3 -march=native -funroll-loops" -DCMAKE_Fortran_FLAGS="-O3 -march=native -funroll-loops" ..
   make
   
.. code-block:: bash

   cd; git clone https://github.com/JohannesBuchner/PyMultiNest.git pymultinest
   cd pymultinest
   python setup.py install   
   
We can then check, if the PyMultiNest installation works:

.. code-block:: bash

   export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$HOME/multinest/MultiNest_v3.12_CMake/multinest/lib/
   mpiexec -n 2 python pymultinest_demo.py

Let's then install GSL:

.. code-block:: bash

   cd; wget -v http://mirror.koddos.net/gnu/gsl/gsl-latest.tar.gz
   tar -xzvf gsl-latest.tar.gz
   cd gsl-{latest} 
   ./configure CC=gcc --prefix=$HOME/gsl
   make
   make check
   make install
   make installcheck
   make clean
   export PATH=$HOME/gsl/bin:$PATH

where ``gsl-{latest}`` should be replaced with the latest version number. Let's
then finally install X-PSI and test that it works:
   
.. code-block:: bash

   cd; cd xpsi;        
   CC=gcc python setup.py install
   cd examples/examples_fast/Modules/
   python main.py

Batch usage
^^^^^^^^^^^

For example job scripts, see the Helios example in :ref:`example_job`.

.. _CALMIPsystem:

CALMIP
------------------------------------

`CALMIP <https://www.calmip.univ-toulouse.fr>`_ is the supercomputer of `Université Fédérale de Toulouse <https://www.univ-toulouse.fr>`_

Installation
^^^^^^^^^^^^

In your ``$HOME`` file system, from the login node, start by loading the necessary modules:

.. code-block:: bash

    module purge
    module load conda
    module load cmake
    module load intel/18.2.199
    module load intelmpi/18.2
    module load gsl/2.5-icc

Then, create the conda environnnement and Install python packages with conda (or pip):

.. code-block:: bash

    conda create -n xpsi --clone base
    conda activate xpsi
    conda install numpy scipy matplotlib wrapt
    conda install cython~=0.29
    conda install h5py
    conda install -c conda-forge fgivenx
    pip install schwimmbad --user

Point to the Intel compilers

.. code-block:: bash

    export FC=ifort
    export CC=icc
    export CXX=icpc

Install mpi4py in your ``$HOME`` (e.g. in ``~/Softwares``):

.. code-block:: bash

    mkdir Softwares
    cd Softwares
    wget https://github.com/mpi4py/mpi4py/releases/download/3.1.5/mpi4py-3.1.5.tar.gz
    tar zxvf mpi4py-3.1.5.tar.gz
    cd mpi4py-3.1.5
    python setup.py build
    python setup.py install
    # Test on login node:
    mpiexec -n 4 python demo/helloworld.py


Download and Install the MultiNest package in your ``$HOME`` (e.g. in ``~/Softwares``:

.. code-block:: bash

    cd ~/Softwares
    git clone https://github.com/farhanferoz/MultiNest.git  ./MultiNest
    cd MultiNest/MultiNest_v3.12_CMake/multinest/
    mkdir build
    cd build
    cmake -DCMAKE_INSTALL_PREFIX=~/Softwares/MultiNest \
                -DCMAKE_{C,CXX}_FLAGS="-O3 -xCORE-AVX512 -mkl" \
                -DCMAKE_Fortran_FLAGS="-O3 -xCORE-AVX512 -mkl" \
                -DCMAKE_C_COMPILER=mpiicc    \
                -DCMAKE_CXX_COMPILER=mpiicpc \
                -DCMAKE_Fortran_COMPILER=mpiifort  ..
    make

    ## Check that libraries have been compiled and are present
    ls ../lib

Install pymultinest in your ``$HOME`` (e.g. in ``~/Softwares``:

.. code-block:: bash

    cd ~/Softwares
    git clone https://github.com/JohannesBuchner/PyMultiNest.git ./pymultinest
    cd pymultinest
    python setup.py install

    # Add MultiNest to Library Path to test PyMultiNest (action to do for every job to run)
    export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$HOME/Softwares/MultiNest/MultiNest_v3.12_CMake/multinest/lib

    # Test pymultinest
    mpiexec -n 2 python pymultinest_demo.py


Clone and Install X-PSI in ~/Softwares/

.. code-block:: bash

    cd ~/Softwares
    git clone https://github.com/xpsi-group/xpsi.git
    cd xpsi/
    LDSHARED="icc -shared" CC=icc python setup.py install

    # Test installation
    cd ~/
    python -c "import xpsi"

    ## Ignore the warnings about GetDist, NestCheck, CornerPlotter
    ##  which are only for PostProcessing (not usually performed on HPC systems).


Set up your library paths:

.. code-block:: bash

    export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$HOME/Softwares/MultiNest/MultiNest_v3.12_CMake/multinest/lib
    export LD_PRELOAD=$MKLROOT/lib/intel64/libmkl_core.so:$MKLROOT/lib/intel64/libmkl_sequential.so

Note that the ``module`` commands, and the library path ``commands`` above will have to be added in your SBATCH script (see :ref:`example_job`) to execute a run.
