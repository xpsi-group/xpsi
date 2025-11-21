""" X-PSI: X-ray Pulse Simulation and Inference

An open-source package for neutron star astrostatistics.

"""

__version__ = "3.1.2"
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

    from xpsi.tools import set_phase_interpolant, set_energy_interpolant
    set_phase_interpolant('Akima')
    set_energy_interpolant('Steffen')

    from xpsi.Parameter import Parameter, Derive
    from xpsi.ParameterSubspace import ParameterSubspace

    from xpsi.Likelihood import Likelihood
    from xpsi.Data import Data
    from xpsi.Instrument import Instrument, InstrumentPileup
    from xpsi.Signal import Signal

    from xpsi.Star import Star
    from xpsi.Spacetime import Spacetime
    from xpsi.Photosphere import Photosphere
    from xpsi.PhotospherePlotter import PhotospherePlotter   
    from xpsi.HotRegion import HotRegion
    from xpsi.HotRegions import HotRegions
    from xpsi.Elsewhere import Elsewhere
    from xpsi.Everywhere import Everywhere

    from xpsi.Background import Background
    from xpsi.Interstellar import Interstellar

    from xpsi.Sample import nested
    from xpsi.Prior import Prior
    try:
        from xpsi.SBI_wrapper import Custom_SBI_Likelihood
    except ImportError:
        _warning('Cannot import torch and test SBI_wrapper.')
        

    if global_imports._size == 1:
        try:
            from xpsi.PostProcessing import *
        except ImportError:
            pass

    __XPSI_SETUP__ = True
