.. _install:

Installation
============

Conda environment duplication
-----------------------------

In the source repository we provide dependency files that can facilitate
the duplication of the environment from which X-PSI ``v0.1`` was released.
This information may be useful if trying to diagnose installation problems.

The development environment is summarised as follows:

    * Ubuntu 14.04
    * GCC 4.8.4
    * Open MPI 1.6.5
    * BLAS, LAPACK, ATLAS installed globally via ``apt``
    * Miniconda (Python 2.7; 64-bit)
    * conda environment exported to ``xpsi/environment.yml``

When inspecting the ``xpsi/environment.yml`` file, note the packages that
where installed into a conda environment via pip. There are a few reasons
for these choices, but one main one is that pip is purely for Python
packages and will not install unwanted non-Python libraries. To be clear, such
libraries would be dependencies, if we had not already satisfied them as listed
above.

The Python packages below can be installed straightforwardly from source
or via a package manager (conda, pip, or a combination), via the instructions
native to the packages. When searching for an open-source package you may need
to add *conda-forge* package channel.

.. note::

    The specifications on this page regard the development environment:
    you are free to set up alternative environment. For installation on a
    high-performance system, instructions on this page, which tailor to a
    self-administered machine, may not be applicable. We direct the reader to
    the :ref:`surfsystems` page for guidance.

Dependencies
------------

X-PSI was developed in Python 2.7 environments. The following
Python packages are required for strictly for likelihood functionality:

* `NumPy <https://docs.scipy.org/doc/numpy/index.html>`_
* `Cython <http://cython.readthedocs.io/en/latest>`_

The following Python packages are required for nested sampling:

* `PyMultiNest <https://github.com/JohannesBuchner/PyMultiNest>`_
* `mpi4py <https://bitbucket.org/mpi4py/mpi4py/downloads/>`_

.. note::

    That ``conda install -c conda-forge pymultinest`` might try to install
    many dependencies, including MPI, BLAS/LAPACK, and a Fortran compiler
    binaries in order to install MultiNest. Moreover, the MultiNest version
    listed is a minor release too low to satisfy all our needs. It might not
    be worth installing MultiNest on your personal machine unless you want
    to apply it more generally or gain experience on test problems. In the
    latter case we provide offer `from source`__ instructions.

The following Python packages are required for full functionality of the
post-processing module:

* `Matplotlib <https://matplotlib.org/>`_
* `getdist <https://getdist.readthedocs.io/en/latest/>`_
  (posterior KDE corner plotting)\ [#]_
* `h5py <http://docs.h5py.org/en/stable/>`_
  (storage of X-ray signals computed from posterior samples)
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

.. [#] The getdist_ software used in :ref:`R19` and which X-PSI ``v0.1``
       interfaces may be cloned as follows:

       .. code-block:: bash

            git clone [--single-branch] -b customisation https://github.com/ThomasEdwardRiley/getdist

.. [#] The nestcheck_ software used in :ref:`R19` and which X-PSI ``v0.1``
       interfaces may be cloned as follows:

       .. code-block:: bash

            git clone [--single-branch] -b feature/getdist_kde https://github.com/ThomasEdwardRiley/nestcheck

.. note::

    X-PSI has several dependencies that are not Python packages. Build and
    install guidelines are given below.

For likelihood evaluation, you require the GNU Scientific Library
(`GSL <https://www.gnu.org/software/gsl/>`_). You also
require an `OpenMP`_-enabled C compiler (known compatibility with icc, gcc,
clang).

.. _OpenMP: http://www.openmp.org

To use `MultiNest`_, you require Version 3.11. To build the MultiNest library,
you require an MPI-wrapped Fortran compiler (e.g., mpifort in Open MPI v1.7+).

.. _MultiNest: https://github.com/farhanferoz/MultiNest

__ source_

.. _source:

From source
-----------

GSL
^^^

To obtain the latest GSL_ source code (v2.5 as of writing):

.. code-block:: bash

   wget -v http://mirror.koddos.net/gnu/gsl/gsl-latest.tar.gz

Untar and navigate to the build directory (e.g. ``cd gsl-latest/build``) and
then build and install:

.. code-block:: bash

    ./configure CC=<path/to/compiler/executable>
    make
    make check
    make install
    make installcheck
    make clean

MultiNest
^^^^^^^^^


X-PSI
^^^^^

Clone X-PSI:

.. code-block:: bash

    git clone https://github.com/ThomasEdwardRiley/xpsi.git <path/to/xpsi>

To build and install ``xpsi`` from the clone root, you require a C compiler:

.. code-block:: bash

    CC=<path/to/compiler/executable> python setup.py install --user

For ``icc``, You may need to prepend this command with
``LDSHARED="icc -shared"``.

Alternatively, to build in-place:

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
naturally be opened in a browser. You need the relevant extensions and a
theme such as `sphinx_rtd_theme`_. Customisation can be made
in the ``xpsi/docs/source/conf.py`` script.

.. _sphinx_rtd_theme: https://sphinx-rtd-theme.readthedocs.io/en/latest/

Note that if you require links to the source code in the HTML files, you need
to ensure Sphinx imports the ``xpsi`` package from the *source* directory
instead of from the ``~/.local/lib`` directory of the user. To enforce this,
insert the path to the source directory into ``sys.path`` in the ``conf.py``
script. Then make sure the extension modules are inside the source directory
-- i.e., the package is built in-place (see above).

.. note::

   To build the documentation, all modules need to be imported, and the
   dependencies that are not resolved will print warning messages.
