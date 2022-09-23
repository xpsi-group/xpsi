.. _install:

Installation
============

X-PSI is an open-source software package that is available on `GitHub
<http://github.com/>`_ and can be cloned as:

.. code-block:: bash

    git clone https://github.com/xpsi-group/xpsi.git </path/to/xpsi>

In this page, we lay down the instructions for installing X-PSI and all the  necessary pre-requisites on your local self-administered system.

.. note::

    For installation on a high-performance computing system, we direct the reader to the :ref:`hpcsystems` page for guidance since the instructions on this page are either not applicable or do not target performance.

.. _dev_env:

Python environment
------------------

X-PSI was developed in Python 2.7, and is in the process of being ported to Python 3.
Fortunately, there are several ways to create a virtual environment with a
different version of Python, without disrupting your Python ecosystem.

This section is divided into two subsections. We recommend that the user follow the instructions in the :ref:`basic_env` subsection to begin with. If faced with installation issues, the user may refer to the :ref:`diagnosis_env` subsection.



.. _basic_env:

Basic Conda environment
^^^^^^^^^^^^^^^^^^^^^^^

In the source directory we provide a basic dependency file that installs
the core Python packages required for *likelihood* functionality. These
packages are:

* `NumPy <https://docs.scipy.org/doc/numpy/index.html>`_
* `Cython <http://cython.readthedocs.io/en/latest>`_

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

and then install the core dependencies listed in `basic_environment.yml`,
such as `NumPy`_ and `Cython`_.

All the following steps need to be performed in this newly created environment which can be activated as:

.. code-block:: bash

    conda activate xpsi

.. _diagnosis_env:

Environment duplication for diagnosis
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In the source repository we provide another dependency file that can facilitate
exact duplication of the environment from which X-PSI ``v0.6`` was
released. This information may be useful if trying to diagnose installation
problems, but can only be expected to be compatible with the same platform. We therefore recommend users to instead try and follow the steps mentioned in the previous subsection, and use this section as more of a guidance if required.

The development environment:

    * Ubuntu 14.04
    * Installed globally via ``apt``: GCC 4.8.4; Open MPI 1.6.5 ("ancient");
      BLAS; LAPACK; ATLAS.
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

To duplicate from file:

.. code-block:: bash

     conda env create -f <path/to/xpsi>/environment.yml

Dependencies
------------

.. note::

    For installing X-PSI on a Mac OS or Windows, please look at the tips below before proceeding with the installation of the various depnedencies.

Python dependencies
^^^^^^^^^^^^^^^^^^^

The following Python packages are required for nested sampling:

* `PyMultiNest <https://github.com/JohannesBuchner/PyMultiNest>`_
  (the interface to the MultiNest library)
* `mpi4py <https://bitbucket.org/mpi4py/mpi4py/downloads/>`_
  (for parallelisation)
* `mpifort <https://anaconda.org/conda-forge/openmpi-mpifort>`_
  (MPI-wrapped Fortran compiler for building library)

.. note::

    That ``conda install -c conda-forge pymultinest`` might try to install
    dependencies in the environment, including binaries for MPI, BLAS/LAPACK,
    and a Fortran compiler, all in order to install MultiNest. Moreover, the
    MultiNest version listed is a minor release too low to satisfy all our
    needs. Although production sampling runs need to be performed on a
    high-performance system and X-PSI can locally be installed without sampling
    functionality, it is advisable to install MultiNest on your
    personal machine to gain experience on application to inexpensive test
    problems. Below we offer `from source`__ instructions.

Running the tests requires:
* `Pytest <http://pytest.org>`_

The following Python packages are required for full functionality of the
post-processing module:

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

    That ``pip install emcee==3.0.2  [--user]`` installs a version working with Python 2.

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
^^^^^^^^^^^

X-PSI has several dependencies that are not Python packages. Build and
install guidelines are given below.

GSL
```

To obtain the latest GSL_ source code (otherwise ``v2.5`` works):

.. code-block:: bash

   wget -v http://mirror.koddos.net/gnu/gsl/gsl-latest.tar.gz

.. note::

    The next steps require an `OpenMP`_-enabled C compiler (known compatibility with ``icc``, ``gcc``, and
    ``clang``). Most linux systems come with `GCC <https://gcc.gnu.org>`_ built-in. To find out the GCC path-executable on your system, run ``which gcc``.

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
`````````

To leverage some capabilities of sample post-processing software you require
`MultiNest`_ ``v3.12``. To build the MultiNest library,
you require an MPI-wrapped Fortran compiler (e.g., ``mpifort`` from Open MPI).

.. _MultiNest: https://github.com/farhanferoz/MultiNest

.. note::

    The following assumes an environment similar to that summarised in
    the in the :ref:`dev_env` section above, specifically to emphasise where an
    MPI compiler wrapper is required.

First clone the repository, then navigate to it and build:

.. code-block:: bash

    git clone https://github.com/farhanferoz/MultiNest.git <path/to/clone>/multinest
    cd <path/to/clone>/multinest/MultiNest_v3.12_CMake/multinest/
    mkdir build
    cd build
    CC=gcc FC=mpif90 CXX=g++ cmake -DCMAKE_{C,CXX}_FLAGS="-O3 -march=native -funroll-loops" -DCMAKE_Fortran_FLAGS="-O3 -march=native -funroll-loops" ..
    make
    ls ../lib/

Use the last command to check for the presence of shared objects. There is
*no* need to ``make install`` as suggested in the source code documentation.

.. note::

    If prompted about missing ``cmake`` and ``gfortran``, they can simply be installed as ``sudo apt-get install cmake gfortran``

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
    python setup.py install [--user]

The package will be installed in your Conda environment (if activated).

.. note::

    Here we clone the PyMultiNest repository. However, for :ref:`R19`,
    working with X-PSI ``v0.1``, we used the repository as frozen in a *fork*.
    To clone this version instead:

    .. code-block:: bash

        git clone https://github.com/ThomasEdwardRiley/PyMultiNest.git <path/to/clone>

    and then simply follow the same installation procedure.

X-PSI
-----

.. _OpenMP: http://www.openmp.org

Finally, to build and install from the X-PSI clone root, execute:

.. code-block:: bash

    CC=<path/to/compiler/executable> python setup.py install [--user]

The ``--user`` flag is optional and specifies where the package is installed;
if you want to install the package in a virtual environment, omit this flag.

For ``icc``, you may need to prepend this command with
``LDSHARED="icc -shared"``. This ensures that both the compiler and linker
are Intel, otherwise the ``gcc`` linker might be invoked.

Provided the GSL ``<prefix>/bin`` is in your ``PATH``
environment variable, the X-PSI ``setup.py`` script will automatically use the
``gsl-config`` executable to link the shared libraries and give the required
C flags for compilation of the X-PSI extensions. Because the library location
will not change for runtime, we state the runtime linking instructions at
compilation in the ``setup.py`` script.

To check whether installation proceeded correctly and the software is functioning as expected,
execute the following:

.. code-block:: bash

    cd examples/examples_fast/Modules/
    python main.py

This module performs a ``likelihood check``. If the likelihood value calculated matches
the given value, X-PSI is functioning as expected, else it will raise an error message.
The module will then initiate sampling using MultiNest (assuming that it's installed),
and given the settings, it should take ~5 minutes. To cancel mid-way press ``ctrl + C``.

.. note::

   The default X-PSI is installed with an analytical blackbody surface emission model extension. If you want to use alternative models for the surface radiation field, you will need to (re-)install / (re-)compile XPSI with the appropriate flags:

   .. code-block:: bash

      CC=<path/to/compiler/executable> python setup.py --help
      CC=<path/to/compiler/executable> python setup.py install [--NumHot] [--NumElse] [--user]

   This will install the numerical atmosphere for the hot regions and/or for the rest of the surface (``elsewhere``). To (re-) install the default blackbody surface emission model, run the command again without the flags:

   .. code-block:: bash

      CC=<path/to/compiler/executable> python setup.py install [--user]

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

To install sphinx, run the following command in the X-PSI environment:

.. code-block:: bash

    conda install sphinx=1.8.5

You then need the relevant extensions and need to ensure versions compatible with python2.
Make sure to run each line individually and not copy-paste the whole block into your terminal for proper installation.

.. code-block:: bash

    conda install -c conda-forge nbsphinx=0.5.1
    conda install decorator=4.4.1
    pip install sphinxcontrib-websupport==1.1.2
    pip install sphinx_rtd_theme==0.4.3

Now the documentation can be compiled using:

.. code-block:: bash

    cd xpsi/docs; [make clean;] make html

To rebuild the documentation after a change to source code docstrings:

.. code-block:: bash

    [CC=<path/to/compiler/executable>] python setup.py install [--user]; cd docs; make clean; make html; cd ..

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

Most of the aforementioned instructions for linux are also applicable for Mac OS.
Here we note some of the changes required.

After creating the environment using the ``basic_environment.yml`` file, install ``xcode`` or ``xcode tools``. Be mindful of the sequence of programs to be installed hereafter.
Use ``pip install`` to download and install ``h5py`` and ``emcee`` (and ``maplotlib``, ``numpy``, ``scipy`` and ``cython`` if not using the ``basic_environment.yml``. You may use the file as a reference of the packages required).

On Mac OS, it's preferable to use ``llvm clang`` rather than ``gcc``. In order to do so, first install  ``homebrew``:

.. code-block:: bash

   /usr/bin/ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"

Install ``llvm`` with homebrew, even if weird messages appear, saying llvm is already present in the Mac OS:

.. code-block:: bash

   brew install llvm

Install ``GSL`` (see above).

Install ``fortran`` before ``MPI``.
If faced with issues when specifying or using gfortran (and it "does not pass simple tests") specify the compiler as being gfortran in the ``mpif90`` wrapper files and delete the files that were already in the build directory.
Once ``MPI`` is installed,
export the following environment variables:

.. code-block:: bash

   export LD_LIBRARY_PATH="/Users/<your_path>/openmpi/lib:$LD_LIBRARY_PATH"
   export PATH=$PATH:/Users/<your_path>/mpi/bin/
   export LDFLAGS="-L/usr/local/opt/llvm/lib"
   export CPPFLAGS="-I/usr/local/opt/llvm/include"
   export KMP_DUPLICATE_LIB_OK=TRUE

Consider adding these lines directly in your bashrc (or equivalent file for a different shell e.g. zshrc).

Install ``X-PSI`` using:

.. code-block:: bash

   CC=/usr/local/opt/llvm/bin/clang python setup.py install [--user]

If it gives problem, remove the ``tools`` and ``surface_radiation_field`` entires from ``setup.py`` of ``X-PSI``.
The line in the setup.py file would then look like:

.. code-block:: bash

   packages = ['xpsi', 'xpsi/PostProcessing']

If you encounter any problems with permissions when installing X-PSI, use the ``--user`` option (This will install X-PSI globally, and not just within your virtual environment).

For compatibility, install the specified ``fgivenx``, ``GetDist`` and ``nestcheck`` (see above).


Tips for installing on Windows
------------------------------

.. note::

    We do not recommend installing and running X-PSI on windows. However, if you must, this section details some of the relevant procedures.

X-PSI was successfully installed and run on Windows in the year 2020, at least for the purpose of likelihood functionality, using the following user-contributed procedure.

.. _Ubuntu: https://www.windowscentral.com/install-windows-subsystem-linux-windows-10

.. _Python 2.7: https://help.dreamhost.com/hc/en-us/articles/115000218612-Installing-a-custom-version-of-Python

.. _virtual Python environment: https://help.dreamhost.com/hc/en-us/articles/215489338-Installing-and-using-virtualenv-with-Python-2

* Clone the X-PSI repository to a directory on your Windows computer (see above).
* Download `Ubuntu`_ for Windows.
* Install `Python 2.7`_.
* Create a `virtual Python environment`_ in an Ubuntu shell.
* Install supporting packages ``pip install matplotlib numpy cython scipy``
  followed by ``sudo apt-get install libgsl-dev``.
* Ensure you are in the X-PSI directory and install X-PSI
  ``CC=gcc python setup.py install``.
* Install any missing packages that you need, e.g., ``pip install h5py`` for
  post-processing functionality if you have posterior sample sets available.
* Install Jupyter notebook using ``pip install notebook``.
* Start the kernel with the command ``Jupyter notebook``.
