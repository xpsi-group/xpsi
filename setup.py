"""
To compile to C code, and then compile the C code and link libraries:
    --> CC=<path/to/compiler/executable> python setup.py install [--user]

If compiling and linking with gcc, which is on $PATH:
    --> CC=gcc python setup.py install [--user]

If compiling and linking with Intel icc, with is on $PATH:
    --> LDSHARED="icc -shared" CC=icc python setup.py install [--user]
"""

import os

from setuptools import setup, Extension

if __name__ == '__main__':
    import numpy
    import sys
    OS = sys.platform
    import os
    from os.path import join

    if 'darwin' in OS or 'linux' in OS:
        print('Operating system: ' + OS)

        try:
            import subprocess as sub
        except ImportError:
            print('The subprocess module is required to locate the GSL library.')
            raise

        try:
            gsl_version = sub.check_output(['gsl-config','--version'])[:-1]
            gsl_prefix = sub.check_output(['gsl-config','--prefix'])[:-1]
        except Exception:
            print('GNU Scientific Library cannot be located.')
            raise
        else:
            print('GSL version: ' + gsl_version)
            libraries = ['gsl','gslcblas','m'] # default BLAS interface for gsl
            library_dirs = [gsl_prefix + '/lib']
            _src_dir = os.path.dirname(os.path.abspath(__file__))
            include_dirs = [gsl_prefix + '/include',
                            numpy.get_include(),
                            join(_src_dir, 'xpsi/include')]

            # point to shared library at compile time so runtime resolution
            # is not affected by environment variables, but is determined
            # by the binary itself
            extra_link_args = ['-Wl,-rpath,%s'%(gsl_prefix+'/lib')]

        # try to get the rayXpanda library:
        # please modify these compilation steps it does not work for your
        # environment; this specification of the (shared) object files
        # seems to work fine for gcc and icc compilers at least
        try:
            import rayXpanda
        except ImportError:
            print('Warning: the rayXpanda package cannot be imported. '
                  'Using fallback implementation.')
            CC = os.environ['CC']
            sub.call(['%s'%CC,
                      '-c',
                      join(_src_dir, 'xpsi/include/rayXpanda/inversion.c'),
                      '-o',
                      join(_src_dir, 'xpsi/include/rayXpanda/inversion.o')])
            sub.call(['%s'%CC,
                      '-c',
                      join(_src_dir, 'xpsi/include/rayXpanda/deflection.c'),
                      '-o',
                      join(_src_dir, 'xpsi/include/rayXpanda/deflection.o')])
            use_rayXpanda = False
        else:
            use_rayXpanda = True

        if use_rayXpanda:
            if 'clang' in os.environ['CC']:
                libraries += ['inversion.so', 'deflection.so']
            else:
                libraries += [':inversion.so', ':deflection.so']
            library_dirs += [rayXpanda.__path__[0]]
            extra_link_args += ['-Wl,-rpath,%s'%rayXpanda.__path__[0]]
        else: # get the native dummy interface
            if 'clang' in os.environ['CC']:
                libraries += ['inversion.o', 'deflection.o']
            else:
                libraries += [':inversion.o', ':deflection.o']
            library_dirs += [join(_src_dir, 'xpsi/include/rayXpanda')]
            extra_link_args += ['-Wl,-rpath,%s'%join(_src_dir,
                                                 'xpsi/include/rayXpanda')]

        try:
            if 'gcc' in os.environ['CC']:
                extra_compile_args=['-fopenmp',
                                    '-march=native',
                                    '-O3',
                                    '-funroll-loops',
                                    '-Wno-unused-function',
                                    '-Wno-uninitialized',
                                    '-Wno-cpp']
                extra_link_args.append('-fopenmp')
            elif 'icc' in os.environ['CC']:
                # on high-performance systems using Intel processors
                # on compute nodes, it is usually recommended to select the
                # instruction set (extensions) optimised for a given processor
                extra_compile_args=['-qopenmp',
                                    '-O3',
                                    '-xHOST',
                                    # alternative instruction set
                                    '-axCORE-AVX2,AVX',
                                    '-funroll-loops',
                                    '-Wno-unused-function']
                extra_link_args.append('-qopenmp')
            elif 'clang' in os.environ['CC']:
                extra_compile_args=['-fopenmp',
                                    '-Wno-unused-function',
                                    '-Wno-uninitialized',
                                    '-Wno-#warnings',
                                    '-Wno-error=format-security']
                extra_link_args.append('-fopenmp')
                # you might need these lookup paths for llvm clang on macOS
                # or you might need to edit these paths for your compiler
                #library_dirs.append('/usr/local/opt/llvm/lib')
                #include_dirs.append('/usr/local/opt/llvm/include')
        except KeyError:
            print('Export CC environment variable to "icc" or "gcc" or '
                  '"clang", or modify the setup script for your compiler.')
            raise
    else:
        print('Unsupported operating system. Manually inspect and modify '
              'setup.py script.')
        raise Exception

    cmdclass = {}

    try:
        import Cython
        print('Cython.__version__ == %s' % Cython.__version__)
        from Cython.Distutils import build_ext
    except ImportError:
        print('Cannot use Cython. Trying to build extension from C files...')
        try:
            from distutils.command import build_ext
        except ImportError:
            print('Cannot import build_ext from distutils...')
            raise
        else:
            cmdclass['build_ext'] = build_ext
            file_extension = '.c'
    else:
        print('Using Cython to build extension from .pyx files...')
        file_extension = '.pyx'
        cmdclass['build_ext'] = build_ext

    def EXTENSION(modname):

        pathname = modname.replace('.', os.path.sep)

        return Extension(modname,
                         [pathname + file_extension],
                         language = 'c',
                         libraries = libraries,
                         library_dirs = library_dirs,
                         include_dirs = include_dirs,
                         extra_compile_args = extra_compile_args,
                         extra_link_args = extra_link_args)

    modnames = ['xpsi.surface_radiation_field.effective_gravity_universal',
                'xpsi.cellmesh.mesh_tools',
                'xpsi.cellmesh.mesh',
                'xpsi.cellmesh.polar_mesh',
                'xpsi.cellmesh.global_mesh',
                'xpsi.cellmesh.rays',
                'xpsi.tools.energy_interpolator',
                'xpsi.tools.energy_integrator',
                'xpsi.tools.phase_integrator',
                'xpsi.tools.phase_interpolator',
                'xpsi.tools.synthesise',
                'xpsi.tools.__init__',
                'xpsi.likelihoods.default_background_marginalisation',
                'xpsi.likelihoods._poisson_likelihood_given_background',
                'xpsi.surface_radiation_field.__init__',
                'xpsi.surface_radiation_field.preload',
                'xpsi.surface_radiation_field.hot',
                'xpsi.surface_radiation_field.elsewhere',
                'xpsi.cellmesh.integrator',
                'xpsi.cellmesh.integrator_for_azimuthal_invariance',
                'xpsi.cellmesh.integrator_for_time_invariance',
                'xpsi.pixelmesh.METRIC_qK',
                'xpsi.pixelmesh.RODES_qK',
                'xpsi.pixelmesh.BOUNDARY_CONDITIONS',
                'xpsi.pixelmesh.surfaceBisection',
                'xpsi.pixelmesh.coordinateTransformation',
                'xpsi.pixelmesh.RK_IP2S_tracer',
                'xpsi.pixelmesh.get_IP_radius',
                'xpsi.pixelmesh.globalRayMap',
                'xpsi.surface_radiation_field.local_variables',
                'xpsi.pixelmesh.integrator']

    extensions = []

    for mod in modnames:
        extensions.append(EXTENSION(mod))

    setup(
        name = 'xpsi',
        version = '0.7.8',
        author = 'Thomas Edward Riley',
        author_email = 't.e.riley@uva.nl; t.riley.phd@gmail.com',
        url = 'https://github.com/ThomasEdwardRiley/xpsi',
        license = 'MIT',
        description = """X-PSI: An open-source package for
                         neutron star X-ray Pulse Simulation and Inference.""",
        long_description = open('README.rst').read(),
        packages = ['xpsi',
                    'xpsi/PostProcessing',
                    'xpsi/cellmesh',
                    'xpsi/tools',
                    'xpsi/surface_radiation_field',
                    'xpsi/likelihoods',
                    'xpsi/pixelmesh'],
        install_requires = ['numpy'],
        setup_requires = ['cython'],
        package_data = {'': ['README.rst', 'LICENSE']},
        include_package_data = True,
        ext_modules = extensions,
        cmdclass = cmdclass,
        classifiers = ['Development Status :: 3 - Alpha',
                       'Intended Audience :: Science/Research',
                       'Operating System :: Linux, macOS',
                       'License :: OSI Approved :: MIT License',
                       'Programming Language :: Python'],
        zip_safe = False,
    )

else:
    pass
