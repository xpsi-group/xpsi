""" X-PSI: X-ray Pulse Simulation and Inference

An open-source package for neutron star astrostatistics.

"""

__version__ = "2.2.7"
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
        name = "X-PSI Silva: X-ray Pulse Simulation and Inference"
        rtds = "https://xpsi-group.github.io/xpsisilva"
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

    from xpsisilva import global_imports

    try:
        import rayXpanda
    except ImportError:
        __rayXpanda_installed__ = False
    else:
        __rayXpanda_installed__ = True
        from .cellmesh import set_rayXpanda_deflection_limit
        set_rayXpanda_deflection_limit(global_imports._pi/2.0)

    from xpsisilva.tools import set_phase_interpolant, set_energy_interpolant
    set_phase_interpolant('Akima')
    set_energy_interpolant('Steffen')

    from xpsisilva.Parameter import Parameter, Derive
    from xpsisilva.ParameterSubspace import ParameterSubspace

    from xpsisilva.Likelihood import Likelihood
    from xpsisilva.Data import Data
    from xpsisilva.Instrument import Instrument
    from xpsisilva.Signal import Signal

    from xpsisilva.Star import Star
    from xpsisilva.Spacetime import Spacetime
    from xpsisilva.Photosphere import Photosphere
    from xpsisilva.HotRegion import HotRegion
    from xpsisilva.TwoHotRegions import TwoHotRegions
    from xpsisilva.HotRegions import HotRegions
    from xpsisilva.Elsewhere import Elsewhere
    from xpsisilva.Everywhere import Everywhere

    from xpsisilva.Background import Background
    from xpsisilva.Interstellar import Interstellar

    from xpsisilva.Sample import nested
    from xpsisilva.Prior import Prior

    if global_imports._size == 1:
        try:
            from xpsisilva.PostProcessing import *
        except ImportError:
            pass

    __XPSI_SETUP__ = True
