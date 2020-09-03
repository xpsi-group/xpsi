from ..global_imports import _pi

def set_rayXpanda_deflection_limit(limit):
    """ Control the domain of applicability of rayXpanda library.

    :param float rayXpanda_defl_lim:
        The limiting deflection angle on the interval (0,pi) radians
        below which the rayXpanda library is called for convenient
        calculation of the lensing factor derivative without need for
        interpolation. Silently ignored if compile-time linking fails
        (due to the rayXpanda package not being installed).

    """
    global __rayXpanda_defl_lim__
    _check_rayXpanda_defl_lim(limit)
    __rayXpanda_defl_lim__ = limit

def _check_rayXpanda_defl_lim(limit):
    if not isinstance(limit, float):
        raise TypeError('rayXpanda deflection limit must be a float.')
    if not 0.0 < limit < _pi:
        raise ValueError('The rayXpanda limit declared is outside the '
                         'domain, [0,pi) radians, of the rayXpanda expansion.')
