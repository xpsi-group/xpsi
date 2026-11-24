
__all__ = ["make_verbose", "verbose", "fragile", "fits_column_keyword"]

from inspect import isgeneratorfunction as _isgeneratorfunction
import wrapt
from xpsi import _verbose

def make_verbose(enter_msg='', exit_msg=''):
    """ Decorator factory for a decorator that controls verbosity. """

    @wrapt.decorator
    def decorator(func, instance, args, kwargs):
        deactivate_all_verbosity = kwargs.pop('deactivate_all_verbosity',
                                              False)
        if deactivate_all_verbosity:
            deactivate_verbosity = True
            _ = kwargs.setdefault('deactivate_verbosity', True)
        else:
            deactivate_verbosity = kwargs.pop('deactivate_verbosity', False)
        if _verbose and not deactivate_verbosity:
            if enter_msg and isinstance(enter_msg, str):
                msg = enter_msg
                print(msg + ('...' if enter_msg[-1] != ':' else ''))
        if _isgeneratorfunction(func):
            for msg in func(*args, **kwargs):
                if _verbose and not deactivate_verbosity:
                    if msg and isinstance(msg, str):
                        print(msg + ('...' if msg[-1] != '.' else ''))
                final = msg  # catch last yield if generator
        else:
            final = func(*args, **kwargs)
        if _verbose and not deactivate_verbosity:
            if exit_msg and isinstance(exit_msg, str):
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


def fits_column_keyword(hdu, column, keyword):
    """ Get a column-indexed keyword (e.g. ``TLMIN``, ``TLMAX``) of a FITS table column.

    In FITS binary tables, keywords such as ``TLMINn`` and ``TLMAXn`` refer to
    the ``n``-th column (1-based), so the keyword name depends on the position
    of the column in the table (e.g. ``TLMIN1`` for ``CHANNEL`` in OGIP PHA
    files, ``TLMIN4`` for ``F_CHAN`` in OGIP RMF files).

    :param hdu:
        A FITS binary table HDU (e.g. ``hdul['SPECTRUM']``).

    :param str column:
        Name of the column (case-insensitive).

    :param str keyword:
        Root of the keyword, without the column number (e.g. ``'TLMIN'``).

    :returns:
        The keyword value, or ``None`` if the keyword is absent from the header.

    :raises KeyError:
        If the column does not exist in the table.
    """
    names = [ name.upper() for name in hdu.columns.names ]
    try:
        n = names.index( column.upper() ) + 1
    except ValueError:
        raise KeyError(f'Column {column} not found in the {hdu.name} extension.')

    return hdu.header.get( f'{keyword.upper()}{n}' )
