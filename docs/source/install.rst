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

Note that ``python >= 3.9.0`` is now required because of the version of
``matplotlib`` needed.

.. _basic_env:

Conda Environment
^^^^^^^^^^^^^^^^^

In the source directory we provide a dependency file ``environment.yml`` that
installs the Python packages required for basic functionality of X-PSI. Note
the script requirement for ``matplotlib``.  For now, we leave out the
multinest/ultranest samplers, because of issues on MacOS.  Details for installation
of those are below.

The content of the ``environment.yml`` are:

.. code-block:: bash

    name: xpsi_py3
    channels:
        - defaults
        - conda-forge
    dependencies:
        - python >= 3.9.0
        - numpy < 2.0.0
        - matplotlib == 3.9.2       # STRICT REQUIREMENT FROM FGIVENX
        - scipy
        - wrapt
        - gsl                       # GNU Science library
        - pytest                    # running functionality self-tests
        - getdist                   # posterior KDE corner plotting
        - tqdm                      # progress bar package
        - h5py                      # storage of X-ray signals computed from posterior samples
        - nestcheck                 # posterior error analysis, plotting, run combination, etc.
        - fgivenx                   # conditional posterior plotting; also required by nestcheck
        - astropy >= 5.2, < 7.0.0   # reading FITS files
        - emcee                     # MCMC sammpler


The core packages required for likelihood functionality are
`numpy <https://docs.scipy.org/doc/numpy/index.html>`_,
`cython <http://cython.readthedocs.io/en/latest>`_,
`matplotlib <https://matplotlib.org/stable/index.html>`_,
`scipy <https://docs.scipy.org/doc//scipy/index.html>`_, and
`wrapt <https://wrapt.readthedocs.io/en/latest/>`_.


To create a virtual environment from this file:

.. code-block:: bash

     conda env create -f <path/to/xpsi>/environment.yml

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



We now install
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

#. `jupyter <https://jupyter-notebook.readthedocs.io/en/stable/>`_ if you want to run X-PSI in a notebook. You may also need the ``ipywidgets`` that can be installed with ``conda install -c conda-forge ipywidgets``.
#. `pytest <https://docs.pytest.org/en/7.2.x/>`_ if you want to run functionality tests for X-PSI.
#. `emcee <https://emcee.readthedocs.io/en/latest/>`_ for optional ensemble-MCMC functionality.
#. `UltraNest <https://johannesbuchner.github.io/UltraNest/readme.html>`_ as alternative sampler.
#. `sbi <https://sbi-dev.github.io/sbi/latest/>`_ for performing Simulation-Based Inference (SBI). This additionally requires installation of `torch and torchvision <https://pytorch.org/get-started/locally/>`_.


.. _nonpython:

Prerequisite Non-Python Packages and PyMultiNest
------------------------------------------------

X-PSI has dependencies that are not Python packages,
or which are Python packages but need to be installed from source (PyMultiNest).
Build and install guidelines are given below.

.. note::

    The next steps require an `OpenMP <http://www.openmp.org>`_-enabled C 
    compiler (known compatibility with ``icc``, ``gcc``, and ``clang``). Most 
    linux systems come with `GCC <https://gcc.gnu.org>`_ built-in. To find out
    the GCC path-executable on your system, run ``which gcc``.


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

Prerequisites for MultiNest are c and fortran
compilers (e.g. ``gcc`` and ``gfortran``), ``cmake``, ``blas``, ``lapack``, and
``atlas``. In case missing them, they can be installed by:

.. code-block:: bash

    sudo apt-get install cmake libblas-dev liblapack-dev libatlas-base-dev

To have MPI-wrapped compilers, one should also install ``mpich`` if not installed already:

.. code-block:: bash

    sudo apt install mpich

Assuming these libraries are available, first clone the repository,
then navigate to it and build:

.. code-block:: bash

    git clone https://github.com/farhanferoz/MultiNest.git <path/to/clone>/multinest
    cd <path/to/clone>/multinest/MultiNest_v3.12_CMake/multinest/
    mkdir build
    cd build
    CC=gcc FC=<path/to/working/mpifortran/compiler/>mpif90 CXX=g++ cmake -DCMAKE_{C,CXX}_FLAGS="-O3 -march=native -funroll-loops" -DCMAKE_Fortran_FLAGS="-O3 -march=native -funroll-loops" ..
    make
    ls ../lib/

.. note::

   We note that new default mpif90 created by mpi4py conda installation may not work here. Thus, one needs to point the path to the native mpif90 compiler of the system (e.g. ``CC=gcc FC=/usr/bin/mpif90 CXX=g++ ...``) or install mpi4py only after MultiNest has been installed and use then ``FC=mpif90``.

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
successfully installed or could not be found.  While X-PSI will run properly,
the nested-sampling capabilities (requiring MultiNest) will crash. The user can
use emcee as the back-up sampler (see example in :doc:`Modeling<Modeling>`).
Note however that the post-processing tutorials have only been implemented
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

A quick check of the X-PSI installation can be done with ``import xpsi``, which
should print to screen something like the following:

.. code-block:: bash

    /=============================================\
    | X-PSI: X-ray Pulse Simulation and Inference |
    |---------------------------------------------|
    |                Version: 3.0.6               |
    |---------------------------------------------|
    |      https://xpsi-group.github.io/xpsi      |
    \=============================================/

    Imported emcee version: 3.1.6
    Warning: Cannot import torch and test SBI_wrapper.
    Imported GetDist version: 1.6.4
    Imported nestcheck version: 0.2.1


Some warnings may appear if you are missing Multinest/Ultranest packages.


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
MultiNest installation. It initiates sampling using MultiNest, and given the
settings, it should take ~5 minutes. To cancel mid-way press ``ctrl + C``.

.. note::

   Note that in X-PSI versions before 2.1.0 the selection of the atmosphere
   extension needed to be done when installing X-PSI using appropriate flags:

   .. code-block:: bash

      CC=<path/to/compiler/executable> python setup.py --help
      CC=<path/to/compiler/executable> python setup.py install [--NumHot] [--NumElse] [--user]

   This installed the numerical atmosphere for the hot regions and/or for
   the rest of the surface (``elsewhere``). To (re-) install the default
   blackbody surface emission model, the following command without the flags
   was used:

   .. code-block:: bash

      CC=<path/to/compiler/executable> python setup.py install [--user]

   For X-PSI versions newer than 2.1.0 atmosphere selection is done without
   reinstalling X-PSI.

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

    conda install "sphinx<7.0"
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

The ``.html`` files can then be found in ``xpsi/docs/build/html``, along with the
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
and ``emcee`` (and ``maplotlib``, ``numpy``, ``scipy`` and ``cython ~= 3.0.11``
if not using the ``environment.yml``. You may use the file as a reference of the
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


If you are facing problems with this installation (e.g., linker problems, or
--fopenmp libraries missing), you may try the following:

.. code-block:: bash

    CC=/usr/local/opt/llvm/bin/clang python setup.py install --noopenmp [--user] 


You may also try to use the MacOS native version of ``clang``:

.. code-block:: bash

    CC=/usr/bin/clang python setup.py install --noopenmp [--user] 



If you encounter any problems with permissions when installing X-PSI, use the
``--user`` option (This will install X-PSI globally, and not just within your
virtual environment).

.. note::

    We are encountering issues with installing MultiNest on Mac and we are working on proposing a solution.


.. note::

   See the :ref:`faq` page for issues that might arise if you are trying to install on Mac using a non-native ``gcc`` compiler.  

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
