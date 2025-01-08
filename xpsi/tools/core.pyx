""" Tools for signal handling.

The extensions in this module are available for users to implement custom
likelihood functions, and for developers who wish to contribute to the source
code.

"""
from ..global_imports import _six, xpsiError

from .phase_integrator import phase_integrator
from .phase_interpolator import phase_interpolator

from .energy_integrator import energy_integrator
from .energy_interpolator import energy_interpolator

from .synthesise import synthesise_exposure
from .synthesise import synthesise_given_total_count_number

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

cdef const gsl_interp_type* _get_phase_interpolant() except *:
    global __phase_interpolant__

    try:
        __phase_interpolant__
    except NameError:
        __phase_interpolant__ = 0 # default Akima (periodic)
    finally:
        if __phase_interpolant__ == 0:
            return gsl_interp_akima_periodic
        elif __phase_interpolant__ == 1:
            return gsl_interp_steffen
        elif __phase_interpolant__ == 2:
            return gsl_interp_cspline_periodic
        else:
            raise ValueError('Invalid phase interpolant setting.')

cdef const gsl_interp_type* _get_energy_interpolant() except *:
    global __energy_interpolant__

    try:
        __energy_interpolant__
    except NameError:
        __energy_interpolant__ = 1 # default Steffen
    finally:
        if __energy_interpolant__ == 0:
            return gsl_interp_akima
        elif __energy_interpolant__ == 1:
            return gsl_interp_steffen
        elif __energy_interpolant__ == 2:
            return gsl_interp_cspline
        else:
            raise ValueError('Invalid eneergy interpolant setting.')
