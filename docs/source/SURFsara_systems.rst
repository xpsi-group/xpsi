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

Start by cleaning your home file system of existing versions of dependencies
and move anything else to some archive in ``$HOME``. Clean ``.bashrc`` and
``.bash_profile`` of anything related to this software. Clean ``.bashrc`` and
``.bash_profile`` of environment variables such as: ``LD_LIBRARY_PATH``,
``LD_PRELOAD``, ``RUN_PATH``, ``PATH``, and ``PYTHONPATH``. Then logout and
log back in order to get a clean environment.

Modify clean environment with intel toolchain information:

.. code-block:: bash

    module load intel/2017b

Point to intel compilers:

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

To prepare MultiNest from ``$HOME``:

.. code-block:: bash

    git clone https://github.com/farhanferoz/MultiNest.git ~/multinest
    cd ~/multinest/MultiNest_v3.11_CMake/multinest
    mkdir build
    cd build
    cmake -DCMAKE_{C,CXX}_FLAGS="-O3 -xAVX -axCORE-AVX2 -funroll-loops" -DCMAKE_Fortran_FLAGS="-O3 -xAVX -axCORE-AVX2 -funroll-loops" ..
    make
    ls ../lib/

Use the last command to check for the presence of shared objects.

Set runtime linking path for MultiNest:

.. code-block:: bash

    export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$HOME/multinest/Multinest_v3.11_CMake/multinest/lib

Now you need the Python interface to MultiNest, starting from ``$HOME``:

.. code-block:: bash

    git clone https://github.com/JohannesBuchner/PyMultiNest.git pymultinest
    cd pymultinest
    python setup.py install --user

To build and install GSL_ from ``$HOME``:

.. code-block:: bash

    wget -v http://mirror.koddos.net/gnu/gsl/http://mirror.koddos.net/gnu/gsl/gsl-latest.tar.gz
    tar -xvf gsl-latest.tar.gz
    cd gsl-latest/build
    ./configure FC=ifort CC=icc CFLAGS='-O3 -xAVX -axCORE-AVX2 -mieee-fp -funroll-loops' --prefix=$HOME/gsl
    make

Optionally ``make check`` can be executed next, but should fail on linear
algebra (linalg) checks because precision checks designed for GNU compiler
collection, not intel. Now:

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
would be invoked

If you ever need to reinstall, first clean to recompile C files:

.. code-block:: bash

    rm -r build dist *egg* xpsi/*/*.c

Lisa
----

The following are *system-specific* instructions for the
`Lisa <https://userinfo.surfsara.nl/systems/lisa>`_ Cluster, *that differ from
the Cartesius instructions given above*.

The instruction sets for targeting intel ivybridge processors (Lisa normal nodes)

To get started, X-PSI and all package and library dependencies need to be
installed. The necessary compilers, wrappers, and low-level parallelisation
libraries are already globally installed on Lisa.

Note that all of the following must be performed on a login node in your
home directory ``$HOME``.

.. _GSL: https://www.gnu.org/software/gsl/

Let's start with GSL_. Assuming you are on your home file system on a login
node, `cd` to the package source code directory (e.g., ``$HOME/src``).
We need to install the library in our home file system, so we give a prefix to
the configure script, 

.. code-block:: bash

    module load gcc
    ./configure CC=gcc --prefix=$HOME/gsl
    make
    make check
    make install
    make installcheck
    make clean

We will now install the various python packages we require. We use the module
``/sara/sw/python-2.7.9/`` and its ``pip`` package manager to install packages
locally in ``$HOME/.local/lib/python2.7/site-packages/`` if they are not
installed globally or are outdated. For emcee_ we want the bleeding-edge
version, so we install from source.

.. _emcee: http://emcee.readthedocs.io/en/latest/

.. code-block:: bash

    module load python/2.7.9
    module load gcc

    export CC=gcc

    pip install --user Cython==0.27.3
    pip install --user mpi4py==2.0.0
    #pip install --user schwimmbad

    git clone https://github.com/dfm/emcee.git
    git cd emcee
    python setup.py install --user
    py.test -v tests
    cd ..
    rm -r emcee

    cd XPSI/src
    python build.py install --user --use-cython
    cd $HOME

Provided the GSL prefix is in your ``PATH`` environment variable (see below for
environment variables), the ``XPSI`` setup script will automatically use
the ``gsl-config`` executable script to link the shared libraries and give the
required cflags for compilation of the ``XPSI`` source code.

We will not use the :mod:`~xpsi.PostProcessing` module, but instead
`scp` output files to a local system to perform plotting.
This circumvents any potential backend problems and permits straightforward
use of IPython for interactive plotting.

.. We will now install `PolyChord`_. Untar the source code archive and `cd` into
    it. Edit the ``PyPolyChord`` target in the ``Makefile``:
    .. code-block:: bash
        PyPolyChord: environment $(LIB_DIR)/libchord.so
            python setup.py install --user
    .. code-block:: bash
        module load python/2.7.9
        module load openmpi/gnu
        #optionally DEBUG=1
        make PyPolyChord MPI=1 COMPILER_TYPE=gnu
        make clean

The following environment variables need to be exported in your job script
script so that all relevant libraries can be located at *runtime* by the
dynamic loader (ensure that the environment variables are only extended, and
not overwritten because module loading modifies these variables):

.. code-block:: bash

    # if you want to ensure that your locally installed packages take
    # precedence over globally installed packages:
    #export PYTHONPATH=$HOME/.local.lib/python2.7/site-packages/:$PYTHONPATH

    # we point the dynamic loader to the runtime path for the GSL library
    # when we link the XPSI binaries into an executable, so we do not require
    # it here:
    export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$HOME/PolyChord/lib

    # if you intend to use PolyChord, the authors require that the dynamic
    # loader imports the MPI library before all others:
    #export LD_PRELOAD=/sara/sw/openmpi-gnu-1.6.5-x/lib/libmpi.so.1:$LD_PRELOAD

If you are to perform small tests on login nodes in your login shell, these
environment variables need to be exported in your ``.bash_profile`` script, or
in your ``.bash.rc`` script which can be sourced by your ``.bash_profile``
script. NB: this is a default behaviour on Lisa.

Unfortunately, the ``/sara/sw/python-2.7.9/`` Python distribution does not
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
multithreading regions are used which do not use the ``OMP_NUM_THREADS``
environment variable, so we can invariantly export it as above. However, the
``MKL_NUM_THREADS`` environment variable should either not be exported (in
which case the ``OMP_NUM_THREADS`` variable is used) or increased so that 
:mod:`numpy` can multithread outside of the our local multithreading regions
in the low-level ``XPSI`` source code.

Note that OpenBLAS may not be compiled against the OpenMP library but use
Pthreads. If :mod:`numpy` *is* linked against MKL, we have covered all
possibilities because MKL whilst uses OpenMP threading but the
``MKL_NUM_THREADS`` environment variable takes precedence if set and thus we
ensure it is set to one.

The GSL library we installed (see above) is not a parallel library itself,
and actually supplies a low-level layer of its own as a CBLAS implementation.
This may be replaced with an optimised implementation, in which case the
question of nested multithreading arises. The OpenBLAS and MKL implementations
can detect whether library calls are made within OpenMP-parallel regions of
the ``XPSI`` source code provided the same threading library is used: e.g.,
OpenBLAS compiled with ``USE_OPENMP=1``, or ``XPSI`` compiled with an Intel
compiler and linked against MKL.


