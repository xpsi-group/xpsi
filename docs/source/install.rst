.. _install:

Installation
============

X-PSI is an open-source software package that is available on `GitHub
<https://github.com/xpsi-group/xpsi.git>`_ and can be cloned as:

.. code-block:: bash

    git clone https://github.com/xpsi-group/xpsi.git </path/to/xpsi>

In this page, we lay down the instructions for installing X-PSI and all the
necessary pre-requisites on your local self-administered system.

.. note::

    For installation on a high-performance computing system, we direct the 
    reader to the :ref:`hpcsystems` page for guidance since the instructions 
    on this page are either not applicable or do not target performance.

.. _dev_env:

Python environment
------------------

X-PSI was originally developed in Python 2.7 and was ported to Python 3 as of 
X-PSI 2. We recommend creating a conda virtual environment with anaconda 3 as
per instructions below so as to not disrupt your Python ecosystem.

.. _basic_env:

Conda Environment
^^^^^^^^^^^^^^^^^

In the source directory we provide a dependency file that installs
the Python packages required for full functionality of X-PSI:

.. code-block:: 
    name: xpsi2
    channels:
      - defaults
      - conda-forge
    dependencies:
      - numpy
      - cython
      - matplotlib
      - scipy
      - wrapt
      - pymultinest  # nested sampling
      - mpi4py  # nested sampling
      - getdist  # posterior KDE corner plotting
      - h5py  # storage of X-ray signals computed from posterior samples
      - nestcheck  # posterior error analysis, plotting, run combination, etc.
      - fgivenx  # conditional posterior plotting; also required by nestcheck

The core packages required for likelihood functionality are:

* `NumPy <https://docs.scipy.org/doc/numpy/index.html>`_
* `Cython <http://cython.readthedocs.io/en/latest>`_
* `Matplotlib <https://matplotlib.org/stable/index.html>`_
* `Scipy <https://docs.scipy.org/doc//scipy/index.html>`_
* `Wrapt <https://wrapt.readthedocs.io/en/latest/>`_

Then, optional packages required for nested sampling are: 

* `pyMultiNest <https://johannesbuchner.github.io/PyMultiNest/>`_
* `mpi4py <http://cython.readthedocs.io/en/latest>`_

Note that pyMultiNest expects a MultiNest, which can be installed 
afterwards. Installation instructions for `multinest`__ are given below. 
Finally, optional packages required for post-processing are:

* `getdist <https://getdist.readthedocs.io/en/latest/>`_
* `h5py <https://docs.h5py.org/en/stable/index.html>`_
* `nestcheck <https://nestcheck.readthedocs.io/en/latest/>`_
* `fgivenx <https://fgivenx.readthedocs.io/en/latest/>`_

If you want to run X-PSI in a
`Jupyter <https://jupyter-notebook.readthedocs.io/en/stable/>`_ notebook, you 
can add this as an entry (``- jupyter``) to the environment file. If
you want to run the functionality tests you require 
`PyTest <https://docs.pytest.org/en/7.2.x/>`_ (``- pytest``). The
`emcee <https://emcee.readthedocs.io/en/latest/>`_ (``- emcee``) package for
ensemble-MCMC is also optional. Alternatively, you can install these 
via Conda (or pip) after environment creation.

To create a virtual environment from this file:

.. code-block:: bash

     conda env create -f <path/to/xpsi>/environment.yml

If Conda does not solve the environment dependencies, you may need to create
an environment manually via

.. code-block:: bash

     conda create -n xpsi

and then install the core dependencies listed in `environment.yml`.

*ALL the following steps should to be performed in this newly created 
environment* which can be activated as:

.. code-block:: bash

    conda activate xpsi



__ source_

.. _source:

From source
^^^^^^^^^^^

X-PSI has several dependencies that are not Python packages. Build and
install guidelines are given below.

GSL
```

To obtain the latest `GSL <https://www.gnu.org/software/gsl/>`_ source code (otherwise ``v2.5`` works):

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

__ multinest_

.. _multinest:

MultiNest
`````````

Although production sampling runs need to be performed on a high-performance 
system and X-PSI can locally be installed without sampling functionality, it is
advisable to install MultiNest on your personal machine to gain experience on
application to inexpensive test problems.

To leverage some capabilities of sample post-processing software you require
`MultiNest`_ ``v3.12``. To build the MultiNest library,
you require an MPI-wrapped Fortran compiler (e.g., 
`openmpi-mpifort <https://anaconda.org/conda-forge/openmpi-mpifort>`_ from 
Open MPI).

.. _MultiNest: https://github.com/farhanferoz/MultiNest

.. note::

    The following assumes an environment similar to that summarised in
    the in the :ref:`dev_env` section above.

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

If you have not already installed mpi4py (e.g. through the
environment file as listed in :ref:`dev_env`), then here is how to do it from
source (e.g., on some path such as ``$HOME``):

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

Now you need the Python to interface with MultiNest. For now, we recommend
looking at: `installing MultiNest 
<https://johannesbuchner.github.io/PyMultiNest/install.html>`.



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
