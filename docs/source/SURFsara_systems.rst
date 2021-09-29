.. _surfsystems:

SURFsara systems
================

The information provided below is intended for users who intend to work on the
SURFsara systems Cartesius and/or Lisa. These specific systems are typically
used by members of the Anton Pannekoek Institute for Astronomy, where X-PSI
was first developed. However, this information may also be translated to
other systems by users looking for guidance.

The X-PSI project is yet to develop a catalogue of tips/guidelines/instructions
for different types of systems or specific systems. Such information should be
gained to some degree as the package is applied by users.

Cartesius
---------

`Cartesius <https://userinfo.surfsara.nl/systems/cartesius>`_ is the Dutch National
Supercomputer.

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

Modify clean environment with Intel toolchain information:

.. code-block:: bash

    module load intel/2017b

Load the module environment:

.. code-block:: bash

    module load pre2019

Point to Intel compilers:

.. code-block:: bash

    export FC=ifort
    export CC=icc
    export CXX=icpc

When compiling, use the ``-ax`` flag to select instruction sets that target
Intel Broadwell/Haswell processors on batch nodes via auto-dispatch.
The ``-x`` instruction set is for sneaky testing purposes on the login node.
Below we write these flags arguments explicitly where required.

The default ``cmake`` module should be sufficient:

.. code-block:: bash

    module load cmake

Prepare Python environment:

.. code-block:: bash

    module load python/2.7.9

We will now install the various python packages we require. We use the module
``/sara/sw/python-2.7.9/`` and its ``pip`` package manager to install packages
locally in ``$HOME/.local/lib/python2.7/site-packages/`` if they are not
installed globally or are outdated.

.. code-block:: bash

    pip install --user Cython==0.27.3
    pip install --user wrap
    pip install --user tools

We set the library and python paths: 

.. code-block:: bash
    
    export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$HOME/multinest/MultiNest_v3.12_CMake/multinest/lib/
    export PYTHONPATH=$HOME/.local/lib/python2.7/site-packages/:$PYTHONPATH

To prepare MPI from ``$HOME``:

.. code-block:: bash

    wget https://bitbucket.org/mpi4py/mpi4py/downloads/mpi4py-3.0.0.tar.gz
    tar -xvf mpi4py-3.0.0.tar.gz
    cd mpi4py-3.0.0
    python setup.py install --user

To test on the login node:

.. code-block:: bash

    mpiexec -n 8 python demo/helloworld.py

Do you see ranks 0 through 7 reporting for duty?

To prepare `MultiNest <https://github.com/farhanferoz/MultiNest>`_ from
``$HOME``:

.. code-block:: bash

    git clone https://github.com/farhanferoz/MultiNest.git ~/multinest
    cd ~/multinest/MultiNest_v3.11_CMake/multinest
    mkdir build
    cd build
    cmake -DCMAKE_{C,CXX}_FLAGS="-O3 -xAVX -axCORE-AVX2 -funroll-loops" -DCMAKE_Fortran_FLAGS="-O3 -xAVX -axCORE-AVX2 -funroll-loops" ..
    make
    ls ../lib/

Use the last command to check for the presence of shared objects.

Now you need the Python interface to MultiNest, starting from ``$HOME``:

.. code-block:: bash

    git clone https://github.com/JohannesBuchner/PyMultiNest.git pymultinest
    cd pymultinest
    python setup.py install --user

.. note::

    We assumed above that nested sampling with `MultiNest`_ is desired. If
    ensemble-MCMC with ``emcee`` is desired, you need to install the Python
    packages ``emcee`` and ``schwimmbad``. We assume the user can infer how to
    do this using the information above and on the :ref:`install` page.

To build and install `GSL <https://www.gnu.org/software/gsl/>`_ from ``$HOME``:

.. code-block:: bash

    wget -v http://mirror.koddos.net/gnu/gsl/http://mirror.koddos.net/gnu/gsl/gsl-latest.tar.gz
    tar -xvf gsl-latest.tar.gz
    mkdir gsl-latest/build
    cd gsl-latest/build
    ./configure FC=ifort CC=icc CFLAGS='-O3 -xAVX -axCORE-AVX2 -mieee-fp -funroll-loops' --prefix=$HOME/gsl
    make

Optionally ``make check`` can be executed next, but should fail on linear
algebra (linalg) checks because precision checks designed for GNU compiler
collection, not Intel. Now:

.. code-block:: bash

    make install

You can check the prefix (which should be ``$HOME/gsl``) and version of GSL
on your path:

.. code-block:: bash

    gsl-config --version
    gsl-config --prefix

Note that if you need to restart installation for some reason, first execute:

.. code-block:: bash

    make clean; make distclean

To prepare X-PSI from ``$HOME``:

.. code-block:: bash

    git clone https://github.com/ThomasEdwardRiley/xpsi.git
    cd xpsi
    LDSHARED="icc -shared" CC=icc python setup.py install --user

This ensures that both the compiler and linker are Intel, otherwise gcc linker
would be invoked. Provided the GSL ``<prefix>/bin`` is in your ``PATH``
environment variable, the X-PSI ``setup.py`` script will automatically use the
``gsl-config`` executable to link the shared libraries and give the required
cflags for compilation of the X-PSI extensions. Because the library location
will not change for runtime, we state the runtime linking instructions at
compilation in the ``setup.py`` script.

If you ever need to reinstall, first clean to recompile C files:

.. code-block:: bash

    rm -r build dist *egg* xpsi/*/*.c

.. note::

    We will not use the :mod:`~xpsi.PostProcessing` module, but instead
    ``scp`` output files to a local system to perform plotting.
    This circumvents any potential backend problems and permits straightforward
    use of IPython for interactive plotting. See also the :ref:`install` page.

Environment variables
^^^^^^^^^^^^^^^^^^^^^

The following environment variables need to be exported in your job script
script so that all relevant libraries can be located at *runtime* by the
dynamic loader (ensure that the environment variables are only extended, and
not overwritten because module loading modifies these variables).

Set runtime linking path for MultiNest:

.. code-block:: bash

    export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$HOME/multinest/Multinest_v3.11_CMake/multinest/lib

We want to ensure that your locally installed Python packages take
precedence over globally installed packages:

.. code-block:: bash

    export PYTHONPATH=$HOME/.local.lib/python2.7/site-packages/:$PYTHONPATH

If you are to perform small tests on login nodes in your login shell, these
environment variables need to be exported in your ``.bash_profile`` script, or
in your ``.bash.rc`` script which can be sourced by your ``.bash_profile``
script (the default default behaviour).

The ``/sara/sw/python-2.7.9/`` Python distribution does not
seem to have :mod:`numpy` linked against the Intel MKL library. Instead it
uses the open-source, multithreaded OpenBLAS library which still offers an
optimised interface to BLAS and LAPACK. However for our purposes on distributed
memory architectures, we  wish to export the following environment variables
in our batch job script if we do not want multithreaded libraries to spawn
worker (OpenMP or POSIX) threads:

.. code-block:: bash

    export OMP_NUM_THREADS=1
    export GOTO_NUM_THREADS=1
    export OPENBLAS_NUM_THREADS=1
    export MKL_NUM_THREADS=1

If we instruct our likelihood evaluation object to OpenMP multithread, local
multithreading regions are used which do not take instructions from the
``OMP_NUM_THREADS`` environment variable, so we can invariantly ``export`` it as
above.
However, the ``MKL_NUM_THREADS`` environment variable should either not be
exported (in which case the ``OMP_NUM_THREADS`` variable is used), or increased
so that :mod:`numpy` can multithread outside of the local multithreading
regions in the X-PSI extension modules.

Note that OpenBLAS may not be compiled against the OpenMP library but instead
use Pthreads. If :mod:`numpy` *is* linked against MKL, we have covered all
possibilities because MKL whilst uses OpenMP threading but the
``MKL_NUM_THREADS`` environment variable takes precedence if set and thus we
ensure it is set to one.

The GSL library we installed (see above) is not a parallel library itself,
and actually supplies a low-level layer of its own as a CBLAS implementation.
This may be replaced with an optimised implementation, in which case the
question of nested multithreading arises. The OpenBLAS and MKL implementations
can detect whether library calls are made within OpenMP-parallel regions of
the X-PSI source code provided the same threading library is used: e.g.,
OpenBLAS compiled with ``USE_OPENMP=1``, or X-PSI compiled with an Intel
compiler and linked against MKL.

Batch usage
^^^^^^^^^^^

For an example job script, refer to :ref:`example_script`.


Lisa
----

The following are the instructions for the
`Lisa <https://userinfo.surfsara.nl/systems/lisa>`_ Cluster, that mostly correspond to those 
of Cartesius instructions given above but with few exceptions.
 
After cleaning your home file system of existing versions of dependencies (as explained for Cartesius), we need to make the environment with Intel toolchain information:

.. code-block:: bash

    module load 2019
    module load intel/2019b
    
The default ``cmake`` module can be obtained from:

.. code-block:: bash

    module load Cmake

To prepare the correct Python environment:

.. code-block:: bash

    module load Python/2.7.15-intel-2019b

Next, we need to use ``pip`` to install packages
locally in ``$HOME/.local/lib/python2.7/site-packages/``:

.. code-block:: bash

    pip install --user Cython==0.28; pip install --user wrapt
    
To prepare MPI and MultiNest from ``$HOME``:

.. code-block:: bash

    wget https://bitbucket.org/mpi4py/mpi4py/downloads/mpi4py-3.0.0.tar.gz
    tar -xvf mpi4py-3.0.0.tar.gz
    cd mpi4py-3.0.0
    python setup.py install --user
    
.. code-block:: bash
    
    git clone https://github.com/farhanferoz/MultiNest.git ~/multinest
    cd ~/multinest/MultiNest_v3.12_CMake/multinest
    mkdir build
    cd build
    cmake -DCMAKE_{C,CXX}_FLAGS="-O3 -xAVX -axCORE-AVX2 -funroll-loops" -DCMAKE_Fortran_FLAGS="-O3 -xAVX -axCORE-AVX2 -funroll-loops" ..
    make   

In case of having problems, the ``-xAVX`` and ``-axCORE-AVX2`` compiler flags can be omitted. 

Finally, we need the Python interface to MultiNest, GSL, and X-PSI, all starting from ``$HOME``:

.. code-block:: bash

    cd; git clone https://github.com/JohannesBuchner/PyMultiNest.git pymultinest
    cd pymultinest
    python setup.py install --user    

.. code-block:: bash

    cd; wget -v http://mirror.koddos.net/gnu/gsl/http://mirror.koddos.net/gnu/gsl/gsl-latest.tar.gz
    tar -xvf gsl-latest.tar.gz
    cd gsl-2.6
    mkdir build; cd build
    ../configure CC=icc CFLAGS='-O3 -xAVX -axCORE-AVX2 -mieee-fp -funroll-loops' --prefix=$HOME/gsl
    make
    make install
    export PATH=$HOME/gsl/bin:$PATH

.. code-block:: bash

    cd; git clone https://github.com/ThomasEdwardRiley/xpsi.git
    cd xpsi
    LDSHARED="icc -shared" CC=icc python setup.py install --user    

The success of different installation steps can be checked as described in the Cartesius instructions, as well as the instructions for possible re-installations. 

The structure of the Lisa filesystem for batch jobs is different to Cartesius. Here is given an example job: 

.. code-block:: bash

   #!/bin/bash
   #SBATCH -N 4
   #SBATCH --tasks-per-node=8
   #SBATCH -t 3-00:00:00
   #SBATCH -p normal
   #SBATCH --job-name=run1

   echo start of job in directory $SLURM_SUBMIT_DIR
   echo number of nodes is $SLURM_JOB_NUM_NODES
   echo the allocated nodes are:
   echo $SLURM_JOB_NODELIST

   module load 2019
   module load intel/2019b
   module load Python/2.7.15-intel-2019b

   cp -r $HOME/NICER_analyses/J0030_ST_PST $TMPDIR

   cd $TMPDIR/J0030_ST_PST

   export PYTHONPATH=$HOME/.local/lib/python2.7/site-packages/:$PYTHONPATH

   export OMP_NUM_THREADS=1
   export OPENBLAS_NUM_THREADS=1
   export GOTO_NUM_THREADS=1
   export MKL_NUM_THREADS=1
   export LD_LIBRARY_PATH=$HOME/multinest/MultiNest_v3.12_CMake/multinest/lib/:$LD_LIBRARY_PATH
   export PATH=$HOME/gsl/bin:$PATH

   srun python main_run1.py > out_run1 2> err_run1

   cp run1* out_run1 err_run1 $HOME/NICER_analyses/J0030_ST_PST/.
   #end of job file

.. todo::

    Things here.
