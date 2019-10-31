.. _install:

Installation
============

X-PSI has a number of dependencies which need to be carefully installed.

X-PSI has been tested in Python 2.7, and the following Python packages are
required:

* :mod:`numpy`
* :mod:`cython`
* :mod:`emcee` (v3.0.1)
* :mod:`h5py`
* :mod:`mpi4py`
* :mod:`schwimmbad`
* `tqdm <https://pypi.python.org/pypi/tqdm>`_
* :mod:`getdist`
* `Matplotlib <https://matplotlib.org/>`_

.. _emcee: http://emcee.readthedocs.io/en/latest/

For likelihood evaluation, you require the GNU Scientific (no reason not to
use v2.5) Library (`GSL <https://www.gnu.org/software/gsl/>`_). You also
require an `OpenMP`_-enabled C compiler (tested with icc, gcc, clang).

.. _OpenMP: http://www.openmp.org

To use `MultiNest`_, you require Version 1.14, which is free for
academic use. The source code can be downloaded once you have registered with
CCPForge as an academic. To build the MultiNest library, you require an
MPI-wrapped Fortran compiler (e.g., mpifort in Open MPI v1.7+).

.. _MultiNest: https://github.com/farhanferoz/MultiNest

.. _source:

From source
-----------

You first need to build GSL_:

.. code-block:: bash

    ./configure CC=<path/to/compiler/executable>
    make
    make check
    make install
    make installcheck
    make clean


To build and install ``xpsi``:

.. code-block:: bash

    python build.py install --user

Alternatively, to build in-place:

.. code-block:: bash

    python build.py build_ext -i

This will build extension modules in the source code directory. You must in
this case ensure that the source code directory is on your ``PYTHONPATH``
environment variable, or inserted into ``sys.path`` within a calling module.

Lisa
----

The following are *system-specific* instructions for the SURFsara
`Lisa <https://userinfo.surfsara.nl/systems/lisa>`_ Cluster.

To get started, ``XPSI`` and all package and library dependencies need to be
installed. The necessary compilers, wrappers, and low-level parallelisation
libraries are already globally installed on Lisa.

Note that all of the following must be performed on a login node in your
home directory ``$HOME``.

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

We will not use :mod:`getdist` or Matplotlib_ on Lisa, but instead `scp` output
files to a local system to perform plotting. This circumvents any potential
backend problems and permits straightforward use of IPython for interactive
plotting.

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

Documentation
-------------

If you wish to compile the documentation, and you are in the ``src`` directory:

.. code-block:: bash

    cd docs

    #optionally:
    #make clean

    make html

The ``.html`` files can then found in ``src/docs/build/html``, along with the
notebooks for the tutorials in this documentation. The ``.html`` files can 
naturally be opened in a browser. To do this you need :mod:`sphinx` and the
relevant extensions and the ``sphinx_rtd_theme``. Customisation can be made
in the ``src/docs/source/conf.py`` script.

Note that if you require links to the source code in the HTML files, you need
to ensure Sphinx imports the ``XPSI`` package from the *source* directory
instead of from the ``~/.local/lib`` directory of the user. To enforce this,
insert the path to the source directory into ``sys.path`` in the ``conf.py``
script. Then make sure the extension modules are inside the source directory
-- i.e., the package is built in-place (see above).





