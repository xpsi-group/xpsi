from __future__ import division

from . import _size

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
           "_csq",
           "_km",
           "_kpc",
           "_keV",
           "_G",
           "_M_s",
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
        print('Matpotlib is necessary for use of the PostProcessing module.')
        _mpl = None
else:
    _mpl = None

# math + phys constants
_pi = _m.pi
_2pi = 2.0 * _pi
_4pi = 4.0 * _pi
_c = 2.99792458E8
_csq = _c * _c
_kpc = 3.08567758e19
_keV = 1.60217662e-16
_G = 6.6730831e-11
_M_s = 1.9887724767047002e30
_h_keV = 4.135667662e-18
_h = 6.62607004e-34
_dpr = 180.0 / _pi
_km = 1.0e3
_k_B = 1.38064852e-23

def gravradius(M):
    """ Get the gravitational radius in km given gravitational mass. """
    return _G * M * _M_s/(1000.0 * _c * _c)

def inv_gravradius(R):
    """ Get the gravitational mass in solar masses given gravitational radius. """
    return R * 1000.0 * _c * _c / (_G * _M_s)

class xpsiError(Exception):
    """ Base exception for xpsi-specific runtime errors. """


__interpolants__ = {'Akima': 0, 'Steffen': 1, 'Cubic': 2}

def get_phase_interpolant():
    """ Get the globally-set phase interpolant. """
    try:
        __phase_interpolant__
    except NameError:
        raise xpsiError('Phase interpolant not set.')
    else:
        for key, value in __interpolants__.items():
            if __phase_interpolant__ == value:
                return key

def set_phase_interpolant(interpolant):
    """ Globally set the phase interpolant to invoke from the GSL library.

    :param str interpolant:
        Name of the spline interpolant. Options are "Akima" (default), "Steffen"
        (pre-v0.6 choice), and "Cubic". The first and last have periodic
        boundary conditions.

    """
    global __phase_interpolant__
    global __interpolants__

    if not isinstance(interpolant, _six.string_types):
        raise TypeError('Interpolant declaration must be a string.')

    if interpolant not in __interpolants__.keys():
        raise ValueError('Invalid interpolant name. See the docstring.')

    __phase_interpolant__ = __interpolants__[interpolant]

def get_energy_interpolant():
    """ Get the globally-set energy interpolant. """
    try:
        __energy_interpolant__
    except NameError:
        raise xpsiError('Energy interpolant not set.')
    else:
        for key, value in __interpolants__.items():
            if __energy_interpolant__ == value:
                return key

def set_energy_interpolant(interpolant):
    """ Globally set the energy interpolant to invoke from the GSL library.

    :param str interpolant:
        Name of the spline interpolant. Options are "Akima" (default), "Steffen"
        (pre-v0.6 choice), and "Cubic". All have natural boundary conditions.

    """
    global __energy_interpolant__
    global __interpolants__

    if not isinstance(interpolant, _six.string_types):
        raise TypeError('Interpolant declaration must be a string.')

    if interpolant not in __interpolants__.keys():
        raise ValueError('Invalid interpolant name. See the docstring.')

    __energy_interpolant__ = __interpolants__[interpolant]
