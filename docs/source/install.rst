.. _install:

Installation
============

X-PSI is an open-source software package that is available on `GitHub
<https://github.com/xpsi-group/xpsi.git>`_ and can be cloned as:

.. code-block:: bash

    git clone https://github.com/xpsi-group/xpsi.git </path/to/xpsi>

In this page, we lay down the instructions for installing X-PSI and all the
necessary prerequisites on your local self-administered system.

.. note::

    For installation on a high-performance computing system, we direct the 
    reader to the :ref:`hpcsystems` page for guidance since the instructions 
    on this page are either not applicable or do not target performance.

.. _dev_env:

Prerequisite Python Packages
----------------------------

X-PSI was originally developed in Python 2.7 and was ported to Python 3 as of 
X-PSI v2.0. We recommend creating a conda virtual environment with anaconda3 as
per instructions below so as to not disrupt your Python ecosystem.

.. _basic_env:

Conda Environment
^^^^^^^^^^^^^^^^^

In the source directory we provide a dependency file ``basic_environment.yml`` that
installs the Python packages required for basic functionality of X-PSI. Its
contents are:

.. code-block:: bash

    name: xpsi_py3
    channels:
        - defaults
    dependencies:
        - numpy
        - cython
        - matplotlib
        - scipy
        - wrapt


The core packages required for likelihood functionality are
`numpy <https://docs.scipy.org/doc/numpy/index.html>`_,
`cython <http://cython.readthedocs.io/en/latest>`_,
`matplotlib <https://matplotlib.org/stable/index.html>`_,
`scipy <https://docs.scipy.org/doc//scipy/index.html>`_, and
`wrapt <https://wrapt.readthedocs.io/en/latest/>`_. 


To create a virtual environment from this file:

.. code-block:: bash

     conda env create -f <path/to/xpsi>/basic_environment.yml

If conda does not solve the environment dependencies, you may need to create
an environment manually via

.. code-block:: bash

     conda create -n xpsi_py3

and then install the core dependencies listed in ``basic_environment.yml`` via
conda.

Activate the environment as:

.. code-block:: bash

    conda activate xpsi_py3

.. note::

    **ALL THE FOLLOWING STEPS SHOULD BE PERFORMED IN THIS NEWLY CREATED
    ENVIRONMENT.** Pay special attention to reactivate the environment if you
    ever have to restart the kernel.
    
Next, install
`mpi4py <https://bitbucket.org/mpi4py/mpi4py/downloads/>`_ which is required for 
nested sampling:

.. code-block:: bash

    conda install -c conda-forge mpi4py


We also need `PyMultiNest <https://github.com/JohannesBuchner/PyMultiNest>`_
(the interface to the MultiNest library) for nested sampling.
However, `conda install -c conda-forge pymultinest` might try
to install dependencies in the environment,
including binaries for MPI, BLAS/LAPACK, and a Fortran compiler,
all in order to install MultiNest. Moreover, the MultiNest version
listed is a minor release too low to satisfy all our needs.
Thus, see the PyMultiNest instructions below.

Then, install optional packages
`getdist <https://getdist.readthedocs.io/en/latest/>`_,
`h5py <https://docs.h5py.org/en/stable/index.html>`_,
`nestcheck <https://nestcheck.readthedocs.io/en/latest/>`_, and
`fgivenx <https://fgivenx.readthedocs.io/en/latest/>`_ which are required for
post-processing:

.. code-block:: bash

    conda install -c conda-forge getdist h5py nestcheck fgivenx

.. note::

    However, to get the most updated versions of getdist and nestcheck (which may be needed by
    some of the X-PSI post-processing features), they should be installed from the source
    (https://github.com/cmbant/getdist and https://github.com/ejhigson/nestcheck)
    by cloning the repositories and running ``python setup.py install`` in them.

In addition, some optional miscellaneous packages are:

#. `jupyter <https://jupyter-notebook.readthedocs.io/en/stable/>`_ if you want to run X-PSI in a notebook.
#. `pytest <https://docs.pytest.org/en/7.2.x/>`_ if you want to run functionality tests for X-PSI.
#. `emcee <https://emcee.readthedocs.io/en/latest/>`_ for optional ensemble-MCMC functionality.


.. _nonpython:

Prerequisite Non-Python Packages and PyMultiNest
------------------------------------------------

X-PSI has several dependencies that are not Python packages,
or which are Python packages but need to be installed from source (PyMultiNest).
Build and install guidelines are given below.

GSL
^^^

GSL is the GNU Scientific Library. To obtain the latest 
`GSL <https://www.gnu.org/software/gsl/>`_ source code (otherwise ``v2.5`` 
works):

.. code-block:: bash

   wget -v http://mirror.koddos.net/gnu/gsl/gsl-latest.tar.gz

.. note::

    The next steps require an `OpenMP <http://www.openmp.org>`_-enabled C 
    compiler (known compatibility with ``icc``, ``gcc``, and ``clang``). Most 
    linux systems come with `GCC <https://gcc.gnu.org>`_ built-in. To find out
    the GCC path-executable on your system, run ``which gcc``.

Untar, navigate to the directory (e.g., ``cd gsl-latest``), and
then build and install:

.. code-block:: bash

    ./configure CC=<path/to/compiler/executable> --prefix=$HOME/gsl
    make
    make check
    make install
    make installcheck
    make clean
    
This will install the library in your ``$HOME``, as an example. Next, add GSL
to your path by adding the following line to ``~/.bashrc``:

.. code-block:: bash

    export PATH=$HOME/gsl/bin:$PATH

You can check the prefix and version of GSL on your path:

.. code-block:: bash

    gsl-config --version
    gsl-config --prefix


.. _multinest:

MultiNest
^^^^^^^^^

Although production sampling runs need to be performed on a high-performance 
system and X-PSI can be installed locally without sampling functionality, it is
advisable to install MultiNest on your personal machine to gain experience in
application to inexpensive test problems. In addition, to leverage some
capabilities of sample post-processing software you 
`require MultiNest <https://github.com/farhanferoz/MultiNest>`_ ``v3.12``.
To build the MultiNest library, you require an MPI-wrapped Fortran compiler
(e.g.,  `openmpi-mpifort <https://anaconda.org/conda-forge/openmpi-mpifort>`_
from Open MPI).

.. note::

    The following assumes you have installed mpi4py. If you
    have not already installed it through the ``environment.yml`` file, you may
    do so e.g. via ``conda install -c conda-forge mpi4py``.

Prerequisites for MultiNest are c and fortran
compilers (e.g. ``gcc`` and ``gfortran``), ``cmake``, ``blas``, ``lapack``, and
``atlas``. In case missing them, they can be installed by:

.. code-block:: bash

    sudo apt-get install cmake libblas-dev liblapack-dev libatlas-base-dev

Assuming these libraries are available, first clone the repository,
then navigate to it and build:

.. code-block:: bash

    git clone https://github.com/farhanferoz/MultiNest.git <path/to/clone>/multinest
    cd <path/to/clone>/multinest/MultiNest_v3.12_CMake/multinest/
    mkdir build
    cd build
    CC=gcc FC=mpif90 CXX=g++ cmake -DCMAKE_{C,CXX}_FLAGS="-O3 -march=native -funroll-loops" -DCMAKE_Fortran_FLAGS="-O3 -march=native -funroll-loops" ..
    make
    ls ../lib/

Now you need the Python interface to MultiNest:

.. code-block:: bash

    git clone https://github.com/JohannesBuchner/PyMultiNest.git <path/to/clone>/pymultinest
    cd <path/to/clone>/pymultinest
    python setup.py install [--user]

The package will be installed in your conda environment, if the environment is activated.
In that case, the optional ``--user`` flag should be omitted.
We also need PyMultinest to interface with MultiNest. To do so, add the
following line to ``~/.bashrc``:

.. code-block:: bash

    export LD_LIBRARY_PATH=/my/directory/MultiNest/lib/:$LD_LIBRARY_PATH
    
It's also good to check whether this has worked. In a new kernel, try 

.. code-block:: bash

    python -c 'import pymultinest'
    
which should import without any errors. If you get ``ERROR:   Could not load
MultiNest library "libmultinest.so"``, that means either MultiNest was not
succesfully installed or could not be found.  While X-PSI will run properly,
the nested-sampling capabilities (requiring MultiNest) will crash. The user can
use emcee as the back-up sampler (see example in :doc:`Modeling<Modeling>`).
Note however that the post-processing turorials have only been implemented
for the outputs of MultiNest.


X-PSI
-----

Finally, to build and install from the X-PSI clone root, execute:

.. code-block:: bash

    CC=<path/to/compiler/executable> python setup.py install [--user]

The ``--user`` flag is optional and specifies where the package is installed;
if you want to install the package in a virtual environment (as recommended), omit this flag.

For ``icc``, you may need to prepend this command with
``LDSHARED="icc -shared"``. This ensures that both the compiler and linker
are Intel, otherwise the ``gcc`` linker might be invoked.

Provided the GSL ``<prefix>/bin`` is in your ``PATH``
environment variable, the X-PSI ``setup.py`` script will automatically use the
``gsl-config`` executable to link the shared libraries and give the required
C flags for compilation of the X-PSI extensions. Because the library location
will not change for runtime, we state the runtime linking instructions at
compilation in the ``setup.py`` script.

A quick check of the X-PSI installation can be done with ``import xpsi``, which
should print to screen something like the following:

.. code-block:: bash

    /=============================================\
    | X-PSI: X-ray Pulse Simulation and Inference |
    |---------------------------------------------|
    |                Version: 2.0.0               |
    |---------------------------------------------|
    |      https://xpsi-group.github.io/xpsi      |
    \=============================================/

    Imported GetDist version: 1.4
    Imported nestcheck version: 0.2.1


.. note::

   Importing X-PSI should not be done in the X-PSI root directory (where the ``setup.py`` file locates).
   Otherwise, a following type of error is expected:
   ``ImportError: cannot import name 'set_phase_interpolant' from 'xpsi.tools' (unknown location)``

For a more complete verification of the X-PSI installation, you can execute
the following:

.. code-block:: bash

    cd examples/examples_fast/Modules/
    python main.py

This module performs a ``likelihood check``. If the likelihood value calculated
matches the given value, X-PSI is functioning as expected, else it will raise
an error message.  The following part of this module requires a functioning
MultiNest installation. It initiate sampling using MultiNest, and given the
settings, it should take ~5 minutes. To cancel mid-way press ``ctrl + C``.

If you ever need to reinstall, first clean to recompile the C files:

.. code-block:: bash

    rm -r build dist *egg* xpsi/*/*.c

Alternatively, to build X-PSI in-place:

.. code-block:: bash

    CC=<path/to/compiler/executable> python setup.py build_ext -i

This will build extension modules in the source code directory. You must in
this case ensure that the source code directory is on your ``PYTHONPATH``
environment variable, or inserted into ``sys.path`` within a calling module.

Documentation
-------------


If you wish to compile the documentation you require 
`Sphinx <http://www.sphinx-doc.org/en/master>`_ and extensions. To install
these, run the following command:

.. code-block:: bash

    conda install sphinx
    conda install -c conda-forge nbsphinx
    conda install decorator
    conda install sphinxcontrib-websupport
    conda install sphinx_rtd_theme

Now the documentation can be compiled using:

.. code-block:: bash

    cd xpsi/docs; [make clean;] make html

To rebuild the documentation after a change to source code docstrings:

.. code-block:: bash

    [CC=<path/to/compiler/executable>] python setup.py install [--user]; cd
    docs; make clean; make html; cd ..

The ``.html`` files can then found in ``xpsi/docs/build/html``, along with the
notebooks for the tutorials in this documentation. The ``.html`` files can
naturally be opened in a browser, handily via a Jupyter session (this is
particularly useful if the edits are to tutorial notebooks).

Note that if you require links to the source code in the HTML files, you need
to ensure Sphinx imports the ``xpsi`` package from the source directory
instead of from the ``~/.local/lib`` directory of the user. To enforce this,
insert the path to the source directory into ``sys.path`` in the ``conf.py``
script. Then make sure the extension modules are inside the source directory
-- i.e., the package is built in-place (see above).

.. note::

   To build the documentation, all modules need to be imported, and the
   dependencies that are not resolved will print warning messages.

Tips for installing on Mac OS
-----------------------------

Most of the aforementioned instructions for linux are also applicable for Mac
OS. Here we note some of the changes required.

After creating the environment using the ``environment.yml`` file, 
install ``xcode`` or ``xcode tools``. Be mindful of the sequence of programs to
be installed hereafter. Use ``pip install`` to download and install ``h5py``
and ``emcee`` (and ``maplotlib``, ``numpy``, ``scipy`` and ``cython`` if not
using the ``environment.yml``. You may use the file as a reference of the
packages required).

On Mac OS, it's preferable to use ``llvm clang`` rather than ``gcc``.  The
``homebrew`` version of ``clang`` works, but some users may face potential
issues (see below for the MacOS native ``clang``).  To use ``homebrew`` version
of ``clang``, first install  ``homebrew``:

.. code-block:: bash

    /usr/bin/ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"

Install ``llvm`` with homebrew, even if weird messages appear, saying llvm is
already present in the Mac OS:

.. code-block:: bash

    brew install llvm 
   
Install ``GSL`` (see above).

Install ``fortran`` before ``MPI``. If faced with issues when specifying or
using gfortran (and it "does not pass simple tests") specify the compiler as
being gfortran in the ``mpif90`` wrapper files and delete the files that were
already in the build directory. Once ``MPI`` is installed, export the following
environment variables:

.. code-block:: bash

    export LD_LIBRARY_PATH="/Users/<your_path>/openmpi/lib:$LD_LIBRARY_PATH"
    export PATH=$PATH:/Users/<your_path>/mpi/bin/
    export LDFLAGS="-L/usr/local/opt/llvm/lib"
    export CPPFLAGS="-I/usr/local/opt/llvm/include"
    export KMP_DUPLICATE_LIB_OK=TRUE


Consider adding these lines directly in your bashrc (or equivalent file for a
different shell e.g. zshrc).

Install ``X-PSI`` using:

.. code-block:: bash

    CC=/usr/local/opt/llvm/bin/clang python setup.py install [--user] 


If you are facing problem with this installation (e.g., linker problems, or
--fopenmp libraries missing), you may try the following:

.. code-block:: bash

    CC=/usr/local/opt/llvm/bin/clang python setup.py install --noopenmp [--user] 


You may also try to use the MacOS native version of ``clang``:

.. code-block:: bash

    CC=/usr/bin/clang python setup.py install --noopenmp [--user] 



If you encounter any problems with permissions when installing X-PSI, use the
``--user`` option (This will install X-PSI globally, and not just within your
virtual environment).

Tips for installing on Windows
------------------------------

.. note::

    We do not recommend installing and running X-PSI on windows. However, if
    you must, this section details some of the relevant procedures.


X-PSI was successfully installed and run on Windows in the year 2020, at least
for the purpose of likelihood functionality, using the following 
user-contributed procedure.

* Clone the X-PSI repository to a directory on your Windows computer (see above).
* Download `Ubuntu <https://www.windowscentral.com/install-windows-subsystem-linux-windows-10>`_ for Windows.
* Install a Anaconda or Miniconda  virtual Python environment in an Ubuntu shell.
* Follow the instructions of this page to install all the python and non-python packages.
