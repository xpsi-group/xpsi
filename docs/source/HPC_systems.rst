.. _hpcsystems:

HPC systems
================

The information provided in this page is for users who intend to work on High-Performance Computing (HPC) systems. These installation instructions are system-specific. X-PSI has already been used on different systems, for some of which, we provide the instructions below. This information may also be translated to other systems by users looking for guidance.


Snellius (SURF)
-------------------

`Snellius <https://servicedesk.surf.nl/wiki/display/WIKI/Snellius>`_ is the Dutch National Supercomputer.

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

Load environment module and modify clean environment with Intel toolchain information:

.. code-block:: bash

    module load 2022
    module load intel/2022a

Prepare a conda environment for X-PSI:

.. code-block:: bash

    module load Anaconda3/2022.05
    git clone https://github.com/xpsi-group/xpsi.git
    cd xpsi
    conda env create -f basic_environment.yml
    conda init
    
For changes to take effect, close and re-open the current shell. After that load the modules again, and activate the environment:  

.. code-block:: bash

    module load 2022
    module load intel/2022a
    module load Anaconda3/2022.05
    conda activate xpsi_py3
    
Next, we point to Intel compilers:

.. code-block:: bash

    export FC=ifort
    export CC=icc
    export CXX=icpc

Below we explicitly specify flag arguments to select intel instruction sets that are compatible with the AMD processors present in Snellius.

Load ``cmake`` module:

.. code-block:: bash

    module load CMake/3.20.1-GCCcore-10.3.0

To prepare MPI from ``$HOME``:

.. code-block:: bash

    cd; wget https://bitbucket.org/mpi4py/mpi4py/downloads/mpi4py-3.1.4.tar.gz
    tar zxvf mpi4py-3.1.4.tar.gz
    cd mpi4py-3.1.4
    python setup.py build   --mpicc=/sw/arch/RHEL8/EB_production/2022/software/impi/2021.6.0-intel-compilers-2022.1.0/mpi/2021.6.0/bin/mpicc
    python setup.py install

To test on the login node:

.. code-block:: bash

    mpiexec -n 4 python demo/helloworld.py

Do you see ranks 0 through 3 reporting for duty?

.. note::

    If MPI raises a warning about missing hydra process manager, run the following code-block:

    .. code-block:: bash

        unset I_MPI_PMI_LIBRARY
        export I_MPI_JOB_RESPECT_PROCESS_PLACEMENT=0


To prepare `MultiNest <https://github.com/farhanferoz/MultiNest>`_ from
``$HOME``:

.. code-block:: bash

    git clone https://github.com/farhanferoz/MultiNest.git ~/multinest
    cd ~/multinest/MultiNest_v3.12_CMake/multinest
    mkdir build
    cd build
    cmake -DCMAKE_{C,CXX}_FLAGS="-O3 -xAVX -axCORE-AVX2 -funroll-loops" -DCMAKE_Fortran_FLAGS="-O3 -xAVX -axCORE-AVX2 -funroll-loops" ..; make
    ls ../lib/

Use the last command to check for the presence of shared objects.

.. note::

    In case the Intel compilers on Snellius run into issues with Intel Math Kernal Library (MKL) due to static linkage, you can solve the problem by setting the appropriate paths to the environment variable for the pre-load libs:

    .. code-block:: bash

        export LD_PRELOAD=/sw/arch/Centos8/EB_production/2021/software/imkl/2021.2.0-iimpi-2021a/mkl/2021.2.0/lib/intel64/libmkl_def.so.1:/sw/arch/Centos8/EB_production/2021/software/imkl/2021.2.0-iimpi-2021a/mkl/2021.2.0/lib/intel64/libmkl_avx2.so.1:/sw/arch/Centos8/EB_production/2021/software/imkl/2021.2.0-iimpi-2021a/mkl/2021.2.0/lib/intel64/libmkl_core.so:/sw/arch/Centos8/EB_production/2021/software/imkl/2021.2.0-iimpi-2021a/mkl/2021.2.0/lib/intel64/libmkl_intel_lp64.so:/sw/arch/Centos8/EB_production/2021/software/imkl/2021.2.0-iimpi-2021a/mkl/2021.2.0/lib/intel64/libmkl_intel_thread.so:/sw/arch/Centos8/EB_production/2021/software/imkl/2021.2.0-iimpi-2021a/compiler/2021.2.0/linux/compiler/lib/intel64_lin/libiomp5.so

    Further details on MKL issues can be found in this `thread <https://community.intel.com/t5/Intel-oneAPI-Math-Kernel-Library/mkl-fails-to-load/m-p/1155538>`_

We also need to the set the environment variable for library path to point at MultiNest:

.. code-block:: bash

    export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$HOME/multinest/MultiNest_v3.12_CMake/multinest/lib/

Now you need the Python interface to MultiNest, starting from ``$HOME``:

.. code-block:: bash

    git clone https://github.com/JohannesBuchner/PyMultiNest.git ~/pymultinest
    cd ~/pymultinest
    python setup.py install

To test the installation of MultiNest and PyMultiNest on the login node:

.. code-block:: bash

    mpiexec -n 2 python pymultinest_demo.py

Do you obtain parameter values and evidences?

.. note::

    We assumed above that nested sampling with `MultiNest`_ is desired. If
    ensemble-MCMC with ``emcee`` is desired, you need to install the Python
    packages ``emcee`` and ``schwimmbad``. We assume the user can infer how to
    do this using the information above and on the :ref:`install` page.

Next, we need to load `GSL <https://www.gnu.org/software/gsl/>`_ and set the `PATH` environment variable:

.. code-block:: bash

    module load GSL/2.7-GCC-10.3.0
    export PATH=/sw/arch/Centos8/EB_production/2021/software/GSL/2.7-GCC-10.3.0:$PATH

To prepare X-PSI from ``$HOME``:

.. code-block:: bash

    cd ~/xpsi
    LDSHARED="icc -shared" CC=icc python setup.py install

This ensures that both the compiler and linker are Intel, otherwise gcc linker
would be invoked. Provided the GSL ``<prefix>/bin`` is in your ``PATH``
environment variable, the X-PSI ``setup.py`` script will automatically use the
``gsl-config`` executable to link the shared libraries and give the required
cflags for compilation of the X-PSI extensions. Because the library location
will not change for runtime, we state the runtime linking instructions at
compilation in the ``setup.py`` script.

.. note::

    Since Snellius uses AMD processors and the Intel instruction sets are internally translated, the installation proceeds while repeating `automatic CPU dispatch` and `icc` warnings.
    These warnings are safe to ignore. However, as they get printed, it takes longer for the installation and can exceed the idle time on the login node, resulting in a `broken pipe`. In this case, it would be preferable to direct the output of the installation into an output file, and if required use a `nohup` or similar command.

If you ever need to reinstall, first clean to recompile C files:

.. code-block:: bash

    rm -r build dist *egg* xpsi/*/*.c

.. note::

    We typically do not used the :mod:`~xpsi.PostProcessing` module, but instead
    ``rsync`` output files to a local system to perform plotting.
    This circumvents any potential backend problems and permits straightforward
    use of IPython for interactive plotting. However, if one wishes to use it on a HPC, it would require installation of `GetDist` and `Nestcheck`. See :ref:`install` page for relevant details.


Batch usage
^^^^^^^^^^^

For an example job script, refer to :ref:`example_job`.

Lisa (SURF)
-----------

`Lisa <https://servicedesk.surf.nl/wiki/display/WIKI/Lisa>`_ follows mostly the exact installation instructions as that of Snellius. Small differences are still to be investigated.

Helios (API)
------------

Helios is a cluster of the Anton Pannekoek Institute for Astronomy. 

Installation
^^^^^^^^^^^^

Let's start by loading the necessary modules and creating a conda environment. At the moment, the installation is known to be working only for the specific python 3.10.6 version, and when conda installing the required python packages separately, as followed:

.. code-block:: bash

   git clone https://github.com/xpsi-group/xpsi.git
   cd xpsi
   module load anaconda3/2021-05
   conda create -n xpsi_py3 python=3.10.6
   conda activate xpsi_py3
   conda install -c conda-forge mpi4py
   conda install cython
   conda install scipy
   conda install matplotlib
   conda install wrapt   
   
.. code-block:: bash

   module load openmpi/3.1.6
     
Let's then test if mpi4py works:

.. code-block:: bash

   cd; wget https://bitbucket.org/mpi4py/mpi4py/downloads/mpi4py-3.1.4.tar.gz
   tar zxvf mpi4py-3.1.4.tar.gz
   cd mpi4py-3.1.4
   mpiexec -n 4 python demo/helloworld.py
   
Let's then install MultiNest and PyMultiNest:
   
.. code-block:: bash
   
   cd; git clone https://github.com/farhanferoz/MultiNest.git multinest
   cd multinest/MultiNest_v3.12_CMake/multinest
   mkdir build
   cd build
   CC=gcc FC=mpif90 CXX=g++ cmake -DCMAKE_{C,CXX}_FLAGS="-O3 -march=native -funroll-loops" -DCMAKE_Fortran_FLAGS="-O3 -march=native -funroll-loops" ..
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
   cd gsl-latest 
   mkdir build 
   cd build
   ../configure CC=gcc --prefix=$HOME/gsl
   make
   make check
   make install
   make installcheck
   make clean
   export PATH=$HOME/gsl/bin:$PATH

Let's then finally install X-PSI and test that it works:
   
.. code-block:: bash

   cd; cd xpsi;        
   CC=gcc python setup.py install
   cd examples/examples_fast/Modules/
   python main.py

Batch usage
^^^^^^^^^^^

For example job scripts, see the Helios example in :ref:`example_job`.

.. _CALMIPsystem:

CALMIP (not updated for python3 yet)
------------------------------------

`CALMIP <https://www.calmip.univ-toulouse.fr>`_ is the supercomputer of `Université Fédérale de Toulouse <https://www.univ-toulouse.fr>`_

Installation
^^^^^^^^^^^^

In your ``$HOME`` file system, from the login node, start by loading the necessary modules:

.. code-block:: bash

    module purge
    module load python/2.7.14
    module load cmake
    module load intel/18.2.199
    module load intelmpi/18.2
    module load gsl/2.5-icc

Then, install/update the required python packages:

.. code-block:: bash

    pip install emcee==3.0.2  --user
    pip install --upgrade numpy --user
    pip install --upgrade Cython --user
    pip install schwimmbad --user


Install MPI4PY in your ``$HOME``:

.. code-block:: bash

    wget https://bitbucket.org/mpi4py/mpi4py/downloads/mpi4py-3.0.0.tar.gz
    tar -xvf mpi4py-3.0.0.tar.gz
    cd mpi4py-3.0.0
    python setup.py install --user

Download the MultiNest package in your ``$HOME``:

.. code-block:: bash

    git clone https://github.com/farhanferoz/MultiNest.git ~/multinest
    cd ~/multinest/MultiNest_v3.11_CMake/multinest
    mkdir build
    cd build


Compile MultiNest in your ``$HOME``, following recommendation from CALMIP support:

.. code-block:: bash

    cmake -DCMAKE_INSTALL_PREFIX=~/multiNest \
            -DCMAKE_{C,CXX}_FLAGS="-O3 -xCORE-AVX512 -mkl" \
            -DCMAKE_Fortran_FLAGS="-O3 -xCORE-AVX512 -mkl" \
            -DCMAKE_C_COMPILER=mpiicc    \
            -DCMAKE_CXX_COMPILER=mpiicpc \
            -DCMAKE_Fortran_COMPILER=mpiifort  ..
    make

Set up your library paths:

.. code-block:: bash

    export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$HOME/multiNest/MultiNest_v3.12_CMake/multinest/lib
    export PYTHONPATH=$HOME/.local/lib/python2.7/site-packages/:$PYTHONPATH
    export LD_PRELOAD=$MKLROOT/lib/intel64/libmkl_core.so:$MKLROOT/lib/intel64/libmkl_sequential.so


Note that the ``module`` commands, and the library path ``commands`` above will have to be added in your SBATCH script (see :ref:`example_job`) to execute a run.
