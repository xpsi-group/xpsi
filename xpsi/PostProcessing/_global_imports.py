__all__ = ["make_verbose",
           "verbose",
           "fragile",
           "_verbose",
           "_warning",
           "xpsiError",
           "random_seed",
           "fix_random_seed",
           "getdist",
           "nestcheck",
           "OrderedDict",
           "ABCMeta",
           "abstractmethod",
           "_six",
           "_os",
           "_sys",
           "_np",
           "AmbiguityError",
           "rc", "rcParams",
           "_mpl",
           "plt",
           "gridspec",
           "cm",
           "_get_default_locator",
           "_get_default_formatter",
           "MultipleLocator",
           "MaxNLocator",
           "AutoMinorLocator",
           "AutoLocator",
           "ScalarFormatter",
           "LogLocator",
           "NullFormatter"]

from abc import ABCMeta, abstractmethod

from xpsi.global_imports import *
from xpsi import _warning, _verbose
from xpsi.utils import make_verbose, verbose, fragile

import wrapt

from collections import OrderedDict

from matplotlib import pyplot as plt
from matplotlib import rc, rcParams
from matplotlib.ticker import MultipleLocator, MaxNLocator, AutoMinorLocator,\
                                AutoLocator, ScalarFormatter, LogLocator,\
                                NullFormatter

def _get_default_locator(prune):
    return MaxNLocator(nbins=5, min_n_ticks=3, prune=prune)

def _get_default_formatter():
    default = ScalarFormatter(useOffset=False)
    default.set_powerlimits(lims = (-2.0,3.0))
    return default

from matplotlib import gridspec
from matplotlib import cm

try:
    import getdist
except ImportError:
   _warning('Cannot import GetDist.')
   getdist = None
else:
    if _verbose:
        print('Imported GetDist version: %s' % getdist.__version__)

    _expected_version = '0.3.1'
    if getdist.__version__ != _expected_version:
        if _verbose:
            _warning('The PostProcessing module is compatible with a specific '
                     'GetDist commit, with version %s, so this module will '
                     'likely not work as intended.' % _expected_version)

    # the following disables getdist.chains.print_load_details
    getdist.chains.print_load_details = False

try:
    import nestcheck
except ImportError:
    _warning('Cannot import nestcheck.')
    nestcheck = None
else:
    if _verbose:
        print('Imported nestcheck version: %s' % nestcheck.__version__)

    _expected_version = '0.2.0'
    if nestcheck.__version__ != _expected_version:
        if _verbose:
            _warning('The PostProcessing module is compatible with a specific '
                     'nestcheck commit, with version %s, so this module will '
                     'likely not work as intended.' % _expected_version)

class AmbiguityError(xpsiError):
    """ Thrown if ambiguous IDs are declared for objects. """

random_seed = None

@wrapt.decorator
def fix_random_seed(func, instance, args, kwargs):
    global random_seed
    state = _np.random.get_state()
    _np.random.seed(random_seed)
    output = func(*args, **kwargs)
    _np.random.set_state(state)
    return output
