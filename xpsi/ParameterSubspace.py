""" Abstract parameter subspace :math:`\mathbb{R}^{d}`. """

from __future__ import division, print_function

from .global_imports import *
from . import global_imports

from abc import ABCMeta

from .Parameter import Parameter

class ParameterSubspace(object):
    """ Ordered parameter subspace.

    :param str prefix:
        In common modelling patterns there are multiple instances
        of a particular subspace, so provide a unique naming prefix. The prefix
        will be suffixed with a double-underscore. Note that this prefix will
        prepend prefixes for encapsulated parameters lower in the hierarchy
        which themselves have identification prefixes.

    :param args:
        Instances of :class:`~.Parameter.Parameter`, or iterables over
        instances of :class:`~.Parameter.Parameter` or
        :class:`~.ParameterSubspace.ParameterSubspace`. All parameters
        and subspaces will be merged into this new subspace.

    """

    __metaclass__ = ABCMeta

    def __init__(self, *args, **kwargs):

        prefix = kwargs.get('prefix', None)

        if prefix is not None:
            self.prefix = prefix

        if not args:
            print('No parameters supplied... empty subspace created.')

        self._params = []
        self.merge(*args)

    def merge(self, *args):
        """ Merge parameters into the subspace, thereby expanding it. """
        for obj in args: # siphon off the parameters for merge
            if isinstance(obj, Parameter):
                self._handle(obj)
            elif obj is not None: # some arguments will default to ``None``
                try:
                    iter(obj)
                except TypeError:
                    print('Argument must be parameter or an iterable of parameters or subspaces.')
                    raise
                else:
                    if isinstance(obj, ParameterSubspace):
                        self._handle_prefix(obj)
                    for o in obj:
                        if isinstance(o, Parameter):
                            self._handle(o)
                        elif isinstance(o, ParameterSubspace):
                            self._handle_prefix(o)
                            for param in o:
                                self._handle(param)

    def _handle(self, param):
        """ Check for duplicates and grow in dimensionality if not. """
        try:
            if len(self._params) >= 1:
                for i in range(1, len(self._params) + 1):
                    temp = self._params[len(self._params) - i]
                    if param is temp:
                        return None # quietly do nothing
                    elif param.name == temp.name:
                        raise ValueError('Duplicated parameter name '
                                         '%s.' % param.name)

        except AttributeError:
            pass

        try:
            self._params.append(param)
        except AttributeError:
            self._params = [param]

        try:
            if param.permit_prepend: # check permissions
                self._prepend(param)
        except AttributeError:
            pass # quietly do nothing

    def _handle_prefix(self, subspace):
        try:
            self.prefix
        except AttributeError:
            pass # quietly do nothing
        else:
            try:
                subspace.prefix
            except AttributeError:
                subspace.prefix = self.prefix
            else:
                if self.prefix not in subspace.prefix:
                    subspace.prefix = self.prefix + '__' + subspace.prefix

    @property
    def params(self):
        """ Get the list of parameter objects. """
        return self._params

    def get_param(self, name):
        """ Get a reference to a parameter object by name. """

        if isinstance(name, _six.string_types):
            for param in self._params:
                if name == param.name:
                    return param

            try:
                if self.prefix + '__' not in name:
                    name = self.prefix + '__' + name
            except AttributeError:
                pass
            else:
                for param in self._params:
                    if name == param.name:
                        return param

            raise KeyError('No parameter in subspace with matching name.')

    def index(self, name):
        """ Get the index of a free parameter. """
        for i, param in enumerate(self):
            if name == param.name:
                return i

            try:
                self.prefix
            except AttributeError:
                continue # with iteration
            else:
                if self.prefix + '__' not in name: # safety guard for prefix
                    if self.prefix + '__' + name == param.name:
                        return i

    def __len__(self):
        """ Get the number of free parameters in the subspace. """
        return sum(1 for _ in self)

    def __str__(self):
        """ Get a summary of the parameters constituting the subspace. """

        summary =  'Free parameters\n'
        summary += '---------------\n'
        for param in self:
            summary += param.name + ': ' + str(param) + '\n'

        derived = ''
        for param in self._params:
            if param.fixed:
                derived += param.name + ': ' + str(param) + '\n'

        if derived:
            summary += '\nDerived/fixed parameters\n'
            summary +=   '------------------------\n'
            summary += derived

        return summary

    def __repr__(self):
        """ Redirect to the string representation. """
        return str(self)

    @property
    def names(self):
        """ Get a list of parameter names. """
        return [param.name for param in self._params]

    @property
    def vector(self):
        """ Get all variable values (free parameters and fixed/derived). """
        return [p.evaluate(self) for p in self._params]

    def __getitem__(self, key):
        """ Pass a string or an integer (latter only for free parameters). """

        if isinstance(key, _six.string_types):
            return self.get_param(key).evaluate(self) # "caller" is subspace

        elif isinstance(key, int):
            try: # note only free parameters considered for this variant
                return [param for param in self][key].evaluate(self)
            except IndexError:
                raise

    def __setitem__(self, key, value):
        """ Pass a string or an integer (latter only for free parameters). """

        if isinstance(key, _six.string_types):
            try:
                if self.prefix + '__' not in key:
                    key = self.prefix + '__' + key
            except AttributeError:
                pass

            for param in self._params:
                if key == param.name:
                    param.value = value
                    del value
                    break
            try:
                value
            except NameError:
                pass
            else:
                raise KeyError('No parameter in subspace with matching name.')
        elif isinstance(key, int):
            try: # note only free parameters considered for this variant
                [param for param in self][key].value = value
            except IndexError:
                print('Invalid index to parameter subspace.')
                raise

    def __iter__(self):
        """ Get an iterator over *free* parameters. """
        self._index = -1
        return self

    def __next__(self):
        """ Redirect for Python 3 compatibility. """
        return self.next()

    def next(self):
        """ Get next free parameter. """
        self._index += 1
        if self._index == len(self._params):
            raise StopIteration
        elif self._params[self._index].fixed:
            return self.next()
        else:
            return self._params[self._index]

    def __call__(self, p = None):
        """ Useful callable way to update parameters with a simple vector.

        The order must match the order in which parameters were merged into
        the subspace. Note that this method is often overwritten in subclasses,
        but this parent method can always be accessed if required.

        :param array-like p:
            New point to update to. If ``None``, the current point is returned.
            The object must be ordered to match that of container ``self``.

        :returns: Current point (if the call was not an update).
        :rtype: array-like

        """
        if p is not None:
            if len(p) != len(self):
                raise ValueError('Parameter vector length mismatch.')
            for param, value in zip(self, p):
                param(value)
        else:
            return [param() for param in self]

    @property
    def needs_update(self):
        """ Do cached dependencies need to be updated? """
        for param in self.params:
            if param.needs_update:
                return True

        return False

    @property
    def cached(self):
        """ Get the cached values of all free parameters in the subspace. """
        return [param.cached for param in self]

    def clear_cache(self):
        """ Clear the cache for all free parameters in the subspace.


        .. note::

            That this also deletes the current parameter values, in
            addition to the cached previous values.

        """
        for param in self:
            del param.cached
            param.value = param.cached

    def _prepend(self, param):
        if self.prefix + '__' not in param.name:
            param.name = self.prefix + '__' + param.name

    @property
    def prefix(self):
        """ Get the unique prefix for subspace identification. """
        return self._prefix

    @prefix.setter
    def prefix(self, prefix):
        if isinstance(prefix, _six.string_types):
            self._prefix = prefix
        else:
            raise ValueError('Invalid type for prefix string.')

    @classmethod
    def _update_doc(cls):
        """ Append the docstring with parameter names in class attributes. """
        try:
            cls.required_names
        except AttributeError:
            pass
        else:
            cls.__doc__ += '\n    Required parameter names:'
            for name in cls.required_names:
                cls.__doc__ += '\n        * ' + name

        cls.__doc__ += '\n'

        try:
            cls.optional_names
        except AttributeError:
            pass
        else:
            cls.__doc__ += '\n    Optional parameter names:'
            for name in cls.optional_names:
                cls.__doc__ += '\n        * ' + name

        cls.__doc__ += '\n'
