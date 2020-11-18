""" X-PSI: X-ray Pulse Simulation and Inference

A open-source package for neutron star astrostatistics.

"""
from __future__ import print_function
__version__ = "0.7.4"
__author__ = "Thomas E. Riley"

try:
    __XPSI_SETUP__
except NameError:
    __XPSI_SETUP__ = False

if not __XPSI_SETUP__:

    try:
        from mpi4py import MPI
    except ImportError:
        _comm = None
        _rank = 0
        _size = 1
        _verbose = True
    else:
        _comm = MPI.COMM_WORLD
        _rank = MPI.COMM_WORLD.rank
        _size = MPI.COMM_WORLD.size
        if _rank == 0:
            _verbose = True
        else:
            _verbose = False

    import six as _six
    from inspect import isgeneratorfunction as _isgeneratorfunction
    import wrapt

    def make_verbose(enter_msg='', exit_msg=''):
        """ Decorator factory for a decorator that controls verbosity. """
        @wrapt.decorator
        def decorator(func, instance, args, kwargs):
            deactivate_all_verbosity = kwargs.pop('deactivate_all_verbosity',
                                                  False)
            if deactivate_all_verbosity:
                deactivate_verbosity = True
                _ =  kwargs.setdefault('deactivate_verbosity', True)
            else:
                deactivate_verbosity = kwargs.pop('deactivate_verbosity', False)
            if _verbose and not deactivate_verbosity:
                if enter_msg and isinstance(enter_msg, _six.string_types):
                    msg = enter_msg
                    print(msg + ('...' if enter_msg[-1] != ':' else ''))
            if _isgeneratorfunction(func):
                for msg in func(*args, **kwargs):
                    if _verbose and not deactivate_verbosity:
                        if msg and isinstance(msg, _six.string_types):
                            print(msg + ('...' if msg[-1] != '.' else ''))
                    final = msg # catch last yield if generator
            else:
                final = func(*args, **kwargs)
            if _verbose and not deactivate_verbosity:
                if exit_msg and isinstance(exit_msg, _six.string_types):
                    if exit_msg == '\n':
                        print(exit_msg)
                    else:
                        print(exit_msg + ('.' if exit_msg[-1] != '.' else ''))
            return final
        return decorator

    class verbose(object):
        """ Context manager for verbosity. """
        def __init__(self, condition, enter_msg, exit_msg):
            self.condition = condition
            self.enter_msg = enter_msg
            self.exit_msg = exit_msg

        def __enter__(self):
            if self.condition and _verbose:
                print(self.enter_msg + '...')
            return self.condition

        def __exit__(self, *args, **kwargs):
            if self.condition and _verbose:
                print(self.exit_msg + '.')

    class fragile(object):
        """ A solution straight from Stack Overflow.

        Reference: questions/11195140/break-or-exit-out-of-with-statement

        """
        class Break(Exception):
            """ Break out of the with statement. """

        def __init__(self, value):
            self.value = value

        def __enter__(self):
            return self.value.__enter__()

        def __exit__(self, etype, value, traceback):
            error = self.value.__exit__(etype, value, traceback)
            if etype == self.Break:
                return True
            return error

    if _verbose:
        vstring = "Version: %s" % __version__
        name = "X-PSI: X-ray Pulse Simulation and Inference"
        rtds = "https://thomasedwardriley.github.io/xpsi/"
        print("/=============================================\\")
        print("| " + name + " |")
        print("|---------------------------------------------|")
        print("|" + vstring.center(len(name)+2) + "|")
        print("|---------------------------------------------|")
        print("|" + rtds.center(len(name)+2) + "|")
        print("\\=============================================/\n")

    def _warning(msg):
        if _verbose:
            _ends = ['.',',',':',';','!','?']
            _append = ('.' if msg[-1] not in _ends else '')
            print('Warning: ' + msg + _append)

    import global_imports

    try:
        import rayXpanda
    except ImportError:
        __rayXpanda_installed__ = False
    else:
        __rayXpanda_installed__ = True
        from cellmesh import set_rayXpanda_deflection_limit
        set_rayXpanda_deflection_limit(global_imports._pi/2.0)

    from tools import set_phase_interpolant, set_energy_interpolant
    set_phase_interpolant('Akima')
    set_energy_interpolant('Steffen')

    from .Parameter import Parameter, Derive
    from .ParameterSubspace import ParameterSubspace

    from .Likelihood import Likelihood
    from .Data import Data
    from .Instrument import Instrument
    from .Signal import Signal

    from .Star import Star
    from .Spacetime import Spacetime
    from .Photosphere import Photosphere
    from .HotRegion import HotRegion
    from .TwoHotRegions import TwoHotRegions
    from .HotRegions import HotRegions
    from .Elsewhere import Elsewhere
    from .Everywhere import Everywhere

    from .Background import Background
    from .Interstellar import Interstellar

    from .Sample import nested
    from .Prior import Prior

    if global_imports._size == 1:
        try:
            from .PostProcessing import *
        except ImportError:
            pass

    __XPSI_SETUP__ = True
