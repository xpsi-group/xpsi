.. _install:

Installation
============

.. _dev_env:

Python environment
------------------

X-PSI was developed in Python 2.7, and has not yet been ported to Python 3.
Fortunately, there are several ways to create a virtual environment with a
different version of Python, without disrupting your Python ecosystem.

Clone X-PSI:

.. code-block:: bash

    git clone https://github.com/ThomasEdwardRiley/xpsi.git <path/to/xpsi>

Basic Conda environment
-----------------------

In the source directory we provide a basic dependency file that installs
the core Python packages required for *likelihood* functionality. These
packages are:

* `NumPy <https://docs.scipy.org/doc/numpy/index.html>`_
* `Cython <http://cython.readthedocs.io/en/latest>`_

For likelihood evaluation, you also require the GNU Scientific Library
(`GSL <https://www.gnu.org/software/gsl/>`_). We have included this in the
Conda environment file, but we give installation
instructions from `source`_ below; in the latter case, you can remove the
GSL entry from the environment file prior to creation.

If you want to run X-PSI in a
`Jupyter <https://jupyter-notebook.readthedocs.io/en/stable/>`_
notebook, you can add this as an entry (e.g., ``- jupyter=1.0``) to the
environment file or you can install it via Conda (or pip) after environment
creation.

To create a virtual environment from file:

.. code-block:: bash

     conda env create -f <path/to/xpsi>/basic_environment.yml

If Conda does not solve the environment dependencies, you may need to create
an environment manually via

.. code-block:: bash

     conda create -n xpsi python=2.7

and then install the core dependencies `NumPy`_ and `Cython`_, and also `GSL`_

Conda environment duplication
-----------------------------

In the source repository we provide a dependency file that can facilitate
exact duplication of the environment from which X-PSI ``v0.2.0-alpha`` was
released. This information may be useful if trying to diagnose installation
problems, but can only be expected to be compatible with the same platform.

The development environment:

    * Ubuntu 14.04
    * Installed globally via ``apt``:
        * GCC 4.8.4
        * Open MPI 1.6.5 ("ancient")
        * BLAS, LAPACK, ATLAS
    * `Miniconda2 <https://docs.conda.io/en/latest/miniconda.html>`_
      (Python 2.7; 64-bit)
    * Conda environment exported to ``xpsi/environment.yml``

When inspecting the ``xpsi/environment.yml`` file, note that most of the
entries were installed via automatic resolution of a strict dependency chain
when core packages were specified. Also note the packages that
were installed into a Conda environment via pip. There are a few reasons
for these choices, but the main one is that pip is purely for Python
packages and will not install unwanted non-Python libraries. To be clear, such
libraries would be dependencies that could have been installed via Conda,
if we had not already satisfied them as listed above in this instance.

The Python packages below can be installed straightforwardly from source
or via a package manager (Conda, pip, or a combination), via the instructions
native to the packages. When searching for an open-source package you may need
to add *conda-forge* package channel.

.. note::

    The specifications on this page regard the development environment:
    you are free to set up alternative environment. For installation on a
    high-performance system, instructions on this page, which tailor to a
    self-administered machine, are either not applicable or do not target
    performance. We direct the reader to the :ref:`surfsystems` page for
    guidance.

To duplicate from file:

.. code-block:: bash

     conda env create -f <path/to/xpsi>/environment.yml

Dependencies
------------

The following Python packages are required for nested sampling:

* `PyMultiNest <https://github.com/JohannesBuchner/PyMultiNest>`_
  (the interface to the MultiNest library)
* `mpi4py <https://bitbucket.org/mpi4py/mpi4py/downloads/>`_
  (for parallelisation)
* `SciPy <https://docs.scipy.org/doc/scipy/reference/>`_
  (optional core package useful for, e.g., inverse prior sampling)

.. note::

    That ``conda install -c conda-forge pymultinest`` might try to install
    dependencies in the environment, including binaries for MPI, BLAS/LAPACK,
    and a Fortran compiler, all in order to install MultiNest. Moreover, the
    MultiNest version listed is a minor release too low to satisfy all our
    needs. Although production sampling runs need to be performed on a
    high-performance system, it is advisable to install MultiNest on your
    personal machine to gain experience on application to inexpensive test
    problems. Below we offer `from source`__ instructions.

The following Python packages are required for full functionality of the
post-processing module:

* `Matplotlib <https://matplotlib.org/>`_
* `GetDist <https://getdist.readthedocs.io/en/latest/>`_
  (posterior KDE corner plotting)\ [#]_
* `h5py <http://docs.h5py.org/en/stable/>`_
  (storage of X-ray signals computed from posterior samples; also used by
  emcee_)
* `nestcheck <https://nestcheck.readthedocs.io/en/latest/>`_
  (posterior error analysis, plotting, run combination, etc.)\ [#]_
* `fgivenx <https://fgivenx.readthedocs.io/en/latest/>`_
  (conditional posterior plotting; also required by nestcheck)

Note that post-processing can generally be done on a desktop computer and thus
these packages are not necessary for running sampling processes on a
high-performance system. If they are not installed, a warning message is
printed or an exception is raised (by the root process if MPI world size >1).

The `emcee <https://emcee.readthedocs.io/en/latest/>`_ Python package for
ensemble-MCMC is optional.

.. note::

    That ``conda install -c conda-forge emcee`` will handle dependencies
    recursively to the extent that MPI would be installed if you accept.

.. rubric:: Footnotes

.. [#] The version of GetDist_ currently compatible with X-PSI, and used in
       :ref:`R19`, is v0.3.1. It may be cloned as follows:

       .. code-block:: bash

          git clone [--single-branch] -b customisation \
          https://github.com/ThomasEdwardRiley/getdist.git

.. [#] The version of nestcheck_ currently compatible with X-PSI, and used in
       :ref:`R19`, is v0.2.0. It may be cloned as follows:

       .. code-block:: bash

          git clone [--single-branch] -b feature/getdist_kde \
          https://github.com/ThomasEdwardRiley/nestcheck.git

__ source_

.. _source:

From source
-----------

X-PSI has several dependencies that are not Python packages. Build and
install guidelines are given below.

GSL
^^^

To obtain the latest GSL_ source code (otherwise ``v2.5`` works):

.. code-block:: bash

   wget -v http://mirror.koddos.net/gnu/gsl/gsl-latest.tar.gz

Untar, navigate to the build directory (e.g., ``cd gsl-latest/build``), and
then build and install:

.. code-block:: bash

    ../configure CC=<path/to/compiler/executable> --prefix=$HOME/gsl
    make
    make check
    make install
    make installcheck
    make clean

This will install the library in your ``$HOME``, as an example. You can check
the prefix and version of GSL on your path:

.. code-block:: bash

    gsl-config --version
    gsl-config --prefix


MultiNest
^^^^^^^^^

To leverage some capabilities of sample post-processing software you require
`MultiNest`_ ``v3.11``. To build the MultiNest library,
you require an MPI-wrapped Fortran compiler (e.g., ``mpifort`` from Open MPI).

.. _MultiNest: https://github.com/farhanferoz/MultiNest

.. note::

    The following assumes an environment similar to that summarised in
    the in the :ref:`dev_env` section above, specifically to emphasise where an
    MPI compiler wrapper is required.

First clone the repository, then navigate to it and build:

.. code-block:: bash

    git clone https://github.com/farhanferoz/MultiNest.git <path/to/clone>/multinest
    cd <path/to/clone>/multinest/MultiNest_v3.11_CMake/multinest
    mkdir build
    cd build
    CC=gcc FC=mpif90 CXX=g++ cmake -DCMAKE_{C,CXX}_FLAGS="-O3 -march=native -funroll-loops" -DCMAKE_Fortran_FLAGS="-O3 -march=native -funroll-loops" ..
    make
    ls ../lib/

Use the last command to check for the presence of shared objects. There is
*no* need to ``make install`` as suggested in the source code documentation.

If you have not already installed mpi4py using pip (or Conda assuming a
different environment setup to that summarised in :ref:`dev_env`), then here
is how to do it from source (e.g., on some path such as ``$HOME``):

.. code-block:: bash

    wget https://bitbucket.org/mpi4py/mpi4py/downloads/mpi4py-3.0.0.tar.gz

    tar -xf mpi4py-3.0.0.tar.gz

    python setup.py build --mpicc=mpicc

    python setup.py install


The package will be installed in your Conda environment (if activated).

To test:

.. code-block:: bash

    mpiexec -n 4 python demo/helloworld.py

Do you see ranks 0 through 3 reporting for duty? The number of MPI processes
might be best set to somewhere between the number of physical cores and
logical cores in your machine for test sampling applications. For a typical
laptop that might be up to ``-n 4``.

Now you need the Python interface to MultiNest:

.. code-block:: bash

    git clone https://github.com/JohannesBuchner/PyMultiNest.git <path/to/clone>/pymultinest
    cd <path/to/clone>/pymultinest
    python setup.py install --user

The package will be installed in your Conda environment (if activated).

.. note::

    Here we clone the PyMultiNest repository. However, for :ref:`R19`,
    working with X-PSI ``v0.1``, we used the repository as frozen in a *fork*.
    To clone this version instead:

    .. code-block:: bash

        git clone https://github.com/ThomasEdwardRiley/PyMultiNest.git <path/to/clone>

    and then simply follow the same installation procedure.

X-PSI
^^^^^

.. _OpenMP: http://www.openmp.org

To build and install from the X-PSI clone root, you require an
`OpenMP`_-enabled C compiler (known compatibility with ``icc``, ``gcc``, and
``clang``):

.. code-block:: bash

    CC=<path/to/compiler/executable> python setup.py install [--user]

For ``icc``, you may need to prepend this command with
``LDSHARED="icc -shared"``. This ensures that both the compiler and linker
are Intel, otherwise the ``gcc`` linker might be invoked.

Provided the GSL ``<prefix>/bin`` is in your ``PATH``
environment variable, the X-PSI ``setup.py`` script will automatically use the
``gsl-config`` executable to link the shared libraries and give the required
C flags for compilation of the X-PSI extensions. Because the library location
will not change for runtime, we state the runtime linking instructions at
compilation in the ``setup.py`` script.

.. note::

   To install X-PSI on Mac OS, you can use ``llvm clang`` rather than ``gcc``.
   First install ``homebrew`` and use that to install ``llvm``:

   .. code-block:: bash

      /usr/bin/ruby -e 
      "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"

      brew install llvm

   Modify your ``.profile`` file as follows:

   .. code-block:: bash

      export PATH=/usr/local/opt/llvm/bin:$PATH
      export LDFLAGS="-L/usr/local/opt/llvm/lib"
      export CPPFLAGS="-I/usr/local/opt/llvm/include"
      export KMP_DUPLICATE_LIB_OK=TRUE

   Install X-PSI using

   .. code-block:: bash

      CC=/usr/local/opt/llvm/bin/clang python setup.py install [--user]


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

.. _Sphinx: http://www.sphinx-doc.org/en/master

If you wish to compile the documentation you require `Sphinx`_:

.. code-block:: bash

    cd xpsi/docs; make html

The ``.html`` files can then found in ``xpsi/docs/build/html``, along with the
notebooks for the tutorials in this documentation. The ``.html`` files can
naturally be opened in a browser. You need the relevant extensions (such as 
``nbsphinx``, which you will be prompted to install) and atheme such as the 
Sphinx `Read the Docs theme`__. Customisation can be made in the 
``xpsi/docs/source/conf.py`` script.

__ https://sphinx-rtd-theme.readthedocs.io/en/latest/

Note that if you require links to the source code in the HTML files, you need
to ensure Sphinx imports the ``xpsi`` package from the source directory
instead of from the ``~/.local/lib`` directory of the user. To enforce this,
insert the path to the source directory into ``sys.path`` in the ``conf.py``
script. Then make sure the extension modules are inside the source directory
-- i.e., the package is built in-place (see above).

.. note::

   To build the documentation, all modules need to be imported, and the
   dependencies that are not resolved will print warning messages.

Installing on Windows
-------------------------------

X-PSI has been successfully installed and run on Windows using the following
procedure.  

.. _Ubuntu: https://www.windowscentral.com/install-windows-subsystem-linux-windows-10

.. _Python 2.7: https://help.dreamhost.com/hc/en-us/articles/115000218612-Installing-a-custom-version-of-Python

.. _virtual Python environment: https://help.dreamhost.com/hc/en-us/articles/215489338-Installing-and-using-virtualenv-with-Python-2

* Clone the X-PSI repository to a directory on your Windows computer (see above).
* Download `Ubuntu`_ for Windows.
* Install `Python 2.7`_.
* Create a `virtual Python environment`_ in an Ubuntu shell.
* Install supporting packages ``pip install matplotlib numpy scipy 
  pymultinest cython`` followed by ``sudo apt-get install libgsl-dev``.
* Ensure you are in the X-PSI directory and install X-PSI 
  ``CC=gcc python setup.py install``.
* Install any missing packages e.g. ``pip install h5py``.
* Install Jupyter notebook using ``pip install notebook``.
* Start the kernel with the command ``Jupyter notebook` and import X-PSI 
  (see tutorials)``.
