"""
Build hook for X-PSI Cython/C extensions.

This file is intentionally kept as a thin shim so that the imperative logic
required to detect GSL, the compiler type, and OpenMP can still live in
Python rather than in a static pyproject.toml.  The actual package metadata
and build-system declaration have moved to pyproject.toml.

Usage
-----
Standard install (replaces the old `python setup.py install`):

    CC=gcc pip install .  #If compiling & linking with gcc, which is on $PATH
    CC=icc LDSHARED="icc -shared" pip install .  #If compiling and linking with Intel icc, with is on $PATH
    CC=clang pip install .  #If compiling and linking with llvm clang, which is on $PATH

No-OpenMP install (replaces the old `--noopenmp` flag):

    XPSI_NO_OPENMP=1 CC=gcc pip install .

In-place build for development (replaces `python setup.py build_ext -i`):

    CC=gcc pip install --no-build-isolation -e .

All old `python setup.py install` invocations now raise an error that
explains the migration.
"""

import os
import sys

# -----------------------------------------------------------------------
# Hard-stop: prevent the old egg-based install path entirely.
# Setuptools >= 80 removed it anyway; this gives a clear message on older
# setuptools too.
# -----------------------------------------------------------------------
if "install" in sys.argv:
    print(
        "\n"
        "ERROR: `python setup.py install` is no longer supported.\n"
        "\n"
        "Please use pip instead:\n"
        "    CC=gcc pip install .\n"
        "\n"
        "To disable OpenMP, set the environment variable before calling pip:\n"
        "    XPSI_NO_OPENMP=1 CC=gcc pip install .\n"
        "\n"
        "For an editable / in-place development install:\n"
        "    CC=gcc pip install --no-build-isolation -e .\n",
        file=sys.stderr,
    )
    sys.exit(1)

# ---------------------------------------------------------------------------
# Everything below is only reached by setuptools/pip during the *build* phase
# (i.e. when pip calls `python setup.py build_ext` internally via the
# setuptools.build_meta backend declared in pyproject.toml).
# ---------------------------------------------------------------------------

from setuptools import Extension, setup  # noqa: E402 — must stay after guard

# ----------------------------------------------------------------
# OpenMP opt-out: honour the XPSI_NO_OPENMP environment variable.
# The old CLI flag `--noopenmp` cannot be used with `pip install`.
# ----------------------------------------------------------------
noopenmp = bool(os.environ.get("XPSI_NO_OPENMP", ""))
print("NOOPENMP =", noopenmp)

# --------------
# Platform check
# --------------
OS = sys.platform
if "darwin" not in OS and "linux" not in OS:
    raise RuntimeError(
        "Unsupported operating system '%s'. "
        "Manually inspect and modify setup.py." % OS
    )
print("Operating system:", OS)

# ----------
# Locate GSL
# ----------
try:
    import subprocess as sub
except ImportError:
    print("The subprocess module is required to locate the GSL library.")
    raise

try:
    gsl_version = sub.check_output(["gsl-config", "--version"])[:-1].decode("UTF-8")
    gsl_prefix  = sub.check_output(["gsl-config", "--prefix"])[:-1].decode("UTF-8")
except Exception:
    print("GNU Scientific Library cannot be located.")
    raise
else:
    print("GSL version:", gsl_version)
    import numpy  # noqa: E402
    _src_dir = os.path.dirname(os.path.abspath(__file__))
    libraries    = ["gsl", "gslcblas", "m"]
    library_dirs = [gsl_prefix + "/lib"]
    include_dirs = [
        gsl_prefix + "/include",
        numpy.get_include(),
        os.path.join(_src_dir, "xpsi/include"),
    ]
    # point to shared library at compile time so runtime resolution
    # is not affected by environment variables, but is determined
    # by the binary itself
    extra_link_args = ["-Wl,-rpath,%s" % (gsl_prefix + "/lib")]

# -----------------------
# Compiler-specific flags
# -----------------------
try:
    print("NOOPENMP =", noopenmp)
    if not noopenmp:
        if "gcc" in os.environ["CC"]:
            extra_compile_args = ["-fopenmp", 
                                  "-march=native", 
                                  "-O3", 
                                  "-funroll-loops",
                                  "-Wno-unused-function", 
                                  "-Wno-uninitialized", 
                                  "-Wno-cpp"]
            extra_link_args.append("-fopenmp")
        elif "icc" in os.environ["CC"]:
            # on high-performance systems using Intel processors
            # on compute nodes, it is usually recommended to select the
            # instruction set (extensions) optimised for a given processor
            extra_compile_args = ["-qopenmp", 
                                  "-O3", 
                                  "-xHOST",
                                  # alternative instruction set
                                  "-axCORE-AVX2,AVX",
                                  "-funroll-loops", 
                                  "-Wno-unused-function"]
            extra_link_args.append("-qopenmp")
        elif "clang" in os.environ["CC"]:
            extra_compile_args = ["-fopenmp",
                                  "-Wno-unused-function", 
                                  "-Wno-uninitialized",
                                  "-Wno-#warnings", 
                                  "-Wno-error=format-security"]
            extra_link_args.append("-fopenmp")
            # you might need these lookup paths for llvm clang on macOS
            # or you might need to edit these paths for your compiler
            #library_dirs.append('/usr/local/opt/llvm/lib')
            #include_dirs.append('/usr/local/opt/llvm/include')
    else:
        if "gcc" in os.environ["CC"]:
            extra_compile_args = ["-march=native", 
                                  "-O3", 
                                  "-funroll-loops",
                                  "-Wno-unused-function", 
                                  "-Wno-uninitialized", 
                                  "-Wno-cpp"]
        elif "icc" in os.environ["CC"]:
            # on high-performance systems using Intel processors
            # on compute nodes, it is usually recommended to select the
            # instruction set (extensions) optimised for a given processor
            extra_compile_args = ["-O3", 
                                  "-xHOST",
                                  # alternative instruction set
                                  "-axCORE-AVX2,AVX",
                                  "-funroll-loops", 
                                  "-Wno-unused-function"]
        elif "clang" in os.environ["CC"]:
            extra_compile_args = ["-Wno-unused-function", 
                                  "-Wno-uninitialized",
                                  "-Wno-#warnings", 
                                  "-Wno-error=format-security"]
            # you might need these lookup paths for llvm clang on macOS
            # or you might need to edit these paths for your compiler
            #library_dirs.append('/usr/local/opt/llvm/lib')
            #include_dirs.append('/usr/local/opt/llvm/include')
except KeyError:
    print('Export CC environment variable to "icc" or "gcc" or '
          '"clang", or modify the setup script for your compiler.')
    raise

# -------------------------------
# Cython vs pre-generated C files
# -------------------------------
cmdclass = {}
try:
    import Cython
    print("Cython.__version__ == %s" % Cython.__version__)
    from Cython.Distutils import build_ext
except ImportError:
    print("Cannot use Cython. Trying to build extension from C files...")
    try:
        from distutils.command import build_ext
    except ImportError:
        print("Cannot import build_ext from distutils...")
        raise
    else:
        cmdclass["build_ext"] = build_ext
        file_extension = ".c"
else:
    print("Using Cython to build extension from .pyx files...")
    file_extension = ".pyx"
    cmdclass["build_ext"] = build_ext

# --------------
# Extension list
# --------------
modnames = [
    "xpsi.surface_radiation_field.effective_gravity_universal",
    "xpsi.cellmesh.mesh_tools",
    "xpsi.cellmesh.mesh",
    "xpsi.cellmesh.polar_mesh",
    "xpsi.cellmesh.global_mesh",
    "xpsi.cellmesh.rays",
    "xpsi.tools.energy_interpolator",
    "xpsi.tools.energy_integrator",
    "xpsi.tools.phase_integrator",
    "xpsi.tools.phase_interpolator",
    "xpsi.tools.synthesise",
    "xpsi.tools.core",
    "xpsi.likelihoods.default_background_marginalisation",
    "xpsi.likelihoods._poisson_likelihood_given_background",
    "xpsi.likelihoods._gaussian_likelihood_given_background_IQU",
    "xpsi.likelihoods._gaussian_likelihood_QnUn",
    "xpsi.likelihoods.compute_expected_counts",
    "xpsi.surface_radiation_field.core",
    "xpsi.surface_radiation_field.preload",
    "xpsi.surface_radiation_field.hot_user",
    "xpsi.surface_radiation_field.hot_BB",
    "xpsi.surface_radiation_field.hot_Num4D",
    "xpsi.surface_radiation_field.hot_BB_burst",
    "xpsi.surface_radiation_field.hot_Num2D",
    "xpsi.surface_radiation_field.hot_Num2D_split",
    "xpsi.surface_radiation_field.hot_Num5D_split",
    "xpsi.surface_radiation_field.hot_wrapper",
    "xpsi.surface_radiation_field.elsewhere_user",
    "xpsi.surface_radiation_field.elsewhere_wrapper",
    "xpsi.cellmesh.integrator",
    "xpsi.cellmesh.integratorIQU",
    "xpsi.cellmesh.integrator_for_azimuthal_invariance",
    "xpsi.cellmesh.integrator_for_azimuthal_invariance_split",
    "xpsi.cellmesh.integratorIQU_for_azimuthal_invariance",
    "xpsi.cellmesh.integratorIQU_for_azimuthal_invariance_split",
    "xpsi.cellmesh.integrator_for_time_invariance",
    "xpsi.pixelmesh.METRIC_qK",
    "xpsi.pixelmesh.RODES_qK",
    "xpsi.pixelmesh.BOUNDARY_CONDITIONS",
    "xpsi.pixelmesh.surfaceBisection",
    "xpsi.pixelmesh.coordinateTransformation",
    "xpsi.pixelmesh.RK_IP2S_tracer",
    "xpsi.pixelmesh.get_IP_radius",
    "xpsi.pixelmesh.globalRayMap",
    "xpsi.surface_radiation_field.local_variables",
    "xpsi.pixelmesh.integrator",
]


def make_extension(modname):
    pathname = modname.replace(".", os.path.sep)
    ext = Extension(modname,
                    [pathname + file_extension],
                    language="c",
                    libraries=libraries,
                    library_dirs=library_dirs,
                    include_dirs=include_dirs,
                    extra_compile_args=extra_compile_args,
                    extra_link_args=extra_link_args,
                )
    ext.cython_directives = {"language_level": "3"}
    return ext

extensions = [make_extension(m) for m in modnames]

# -----------------------------------------------------------
# setup() — metadata is minimal here; the canonical source is
# pyproject.toml.  Only ext_modules and cmdclass are needed.
# -----------------------------------------------------------
setup(
    ext_modules=extensions,
    cmdclass=cmdclass,
)
