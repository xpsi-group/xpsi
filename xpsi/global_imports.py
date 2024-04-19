from xpsi import _size

__all__ = ["_np",
           "_mpl",
           "_m",
           "_sys",
           "_os",
           "_six",
           "_contig",
           "_pi",
           "_2pi",
           "_4pi",
           "_c",
           "_c_cgs",
           "_csq",
           "_km",
           "_kpc",
           "_keV",
           "_G",
           "_GM",
           "_h_keV",
           "_h",
           "_k_B",
           "_dpr",
           "_cos",
           "_sin",
           "_arccos",
           "_arcsin",
           "_arctan",
           "_sqrt",
           "_exp",
           "_log",
           "_time",
           "gravradius",
           "xpsiError"]

import os as _os
import sys as _sys

import numpy as _np
import math as _m

import six as _six

from numpy import ascontiguousarray as _contig
from numpy import cos as _cos
from numpy import sin as _sin
from numpy import arccos as _arccos
from numpy import arcsin as _arcsin
from numpy import arctan as _arctan
from numpy import sqrt as _sqrt
from numpy import log as _log
from numpy import exp as _exp
from time import time as _time

if _size == 1: # else assume production run usage on HPC system
    try:
        import matplotlib as _mpl
    except ImportError:
        print('Warning: could not import a matplotlib library.')
        print('         X-PSI Inference runs should still work.')
        print('         Matpotlib is necessary for the PostProcessing module and producing sky maps.')
        _mpl = None
else:
    _mpl = None

# math + phys constants
_pi = _m.pi
_2pi = 2.0 * _pi
_4pi = 4.0 * _pi
_c = 2.99792458E8
_c_cgs = _c*1E2
_csq = _c * _c
_kpc = 3.08567758e19
_keV = 1.60217662e-16
_GM = 0.5 * 2.95325024e3 # Solar Schwarzschild gravitational radius in SI
_G = 6.6740831e-11
_h_keV = 4.135667662e-18 #Conversion from Hz^{-1} and keV^{-1} units
_h = 6.62607004e-34
_dpr = 180.0 / _pi
_km = 1.0e3
_k_B = 1.38064852e-23

def gravradius(M):
    """ Get the gravitational radius in km given gravitational mass in Solar masses. """
    return M * _GM / _km

def inv_gravradius(R):
    """ Get the gravitational mass in solar masses given gravitational radius in km. """
    return R * _km / _GM

class xpsiError(Exception):
    """ Base exception for xpsi-specific runtime errors. """

