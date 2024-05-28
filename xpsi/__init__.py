""" X-PSI: X-ray Pulse Simulation and Inference

An open-source package for neutron star astrostatistics.

"""

__version__ = "2.2.4"
__author__ = "The X-PSI Core Team"

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


    if _verbose:
        vstring = "Version: %s" % __version__
        name = "X-PSI: X-ray Pulse Simulation and Inference"
        rtds = "https://xpsi-group.github.io/xpsi"
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

    from xpsi import global_imports

    try:
        import rayXpanda
    except ImportError:
        __rayXpanda_installed__ = False
    else:
        __rayXpanda_installed__ = True
        from .cellmesh import set_rayXpanda_deflection_limit
        set_rayXpanda_deflection_limit(global_imports._pi/2.0)

    from xpsi.tools import set_phase_interpolant, set_energy_interpolant
    set_phase_interpolant('Akima')
    set_energy_interpolant('Steffen')

    from xpsi.Parameter import Parameter, Derive
    from xpsi.ParameterSubspace import ParameterSubspace

    from xpsi.Likelihood import Likelihood
    from xpsi.Data import Data
    from xpsi.Instrument import Instrument
    from xpsi.Signal import Signal

    from xpsi.Star import Star
    from xpsi.Spacetime import Spacetime
    from xpsi.Photosphere import Photosphere
    from xpsi.HotRegion import HotRegion
    from xpsi.TwoHotRegions import TwoHotRegions
    from xpsi.HotRegions import HotRegions
    from xpsi.Elsewhere import Elsewhere
    from xpsi.Everywhere import Everywhere

    from xpsi.Background import Background
    from xpsi.Interstellar import Interstellar

    from xpsi.Sample import nested
    from xpsi.Prior import Prior

    if global_imports._size == 1:
        try:
            from xpsi.PostProcessing import *
        except ImportError:
            pass

    __XPSI_SETUP__ = True
