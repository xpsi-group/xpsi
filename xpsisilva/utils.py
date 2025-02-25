
__all__ = ["make_verbose", "verbose", "fragile"]

import six as _six
from inspect import isgeneratorfunction as _isgeneratorfunction
import wrapt
from xpsisilva import _verbose

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
            if enter_msg and isinstance(enter_msg, _six.string_types):
                msg = enter_msg
                print(msg + ('...' if enter_msg[-1] != ':' else ''))
        if _isgeneratorfunction(func):
            for msg in func(*args, **kwargs):
                if _verbose and not deactivate_verbosity:
                    if msg and isinstance(msg, _six.string_types):
                        print(msg + ('...' if msg[-1] != '.' else ''))
                final = msg  # catch last yield if generator
        else:
            final = func(*args, **kwargs)
        if _verbose and not deactivate_verbosity:
            if exit_msg and isinstance(exit_msg, _six.string_types):
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
