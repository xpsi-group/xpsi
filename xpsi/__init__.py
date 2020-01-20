""" X-PSI: A prototype open-source package for neutron star
    X-ray Pulse Simulation and Inference. """
from __future__ import print_function
__version__ = "0.3.3"
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
    from functools import wraps

    def make_verbose(enter_msg='', exit_msg=''):
        def decorator(func):
            @wraps(func)
            def wrapper(*args, **kwargs):
                if _verbose:
                    if enter_msg and isinstance(enter_msg, _six.string_types):
                        msg = enter_msg
                        print(msg + ('...' if enter_msg[-1] != ':' else ''))
                if _isgeneratorfunction(func):
                    for msg in func(*args, **kwargs):
                        if _verbose:
                            if msg and isinstance(msg, _six.string_types):
                                print(msg + ('...' if msg[-1] != '.' else ''))
                        _ = msg # catch last yield if generator
                else:
                    _ = func(*args, **kwargs)
                if _verbose:
                    if exit_msg and isinstance(exit_msg, _six.string_types):
                        if exit_msg == '\n':
                            print(exit_msg)
                        else:
                            print(exit_msg + ('.' if exit_msg[-1] != '.' else ''))
                return _
            return wrapper
        return decorator

    class verbose(object):
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
        class Break(Exception):
            """Break out of the with statement"""

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

    import global_imports

    from .Parameter import Parameter, Derive
    from .ParameterSubspace import ParameterSubspace

    from .Likelihood import Likelihood
    from .Data import Data
    from .Instrument import Instrument
    from .Pulse import Pulse

    from .Star import Star
    from .Spacetime import Spacetime
    from .Photosphere import Photosphere
    from .HotRegion import HotRegion
    from .TwoHotRegions import TwoHotRegions
    from .HotRegions import HotRegions
    from .Elsewhere import Elsewhere

    from .Background import Background
    from .Interstellar import Interstellar

    from .Sample import nested
    from .Prior import Prior

    if global_imports._size == 1:
        from .PostProcessing import *

