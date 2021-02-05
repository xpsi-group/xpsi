from __future__ import division, print_function

from .global_imports import *
from . import global_imports

from . import make_verbose

import string

from types import MethodType
from abc import ABCMeta, abstractmethod

class StrictBoundsError(xpsiError):
    """ Raised if the set parameter value lies beyond strict bounds. """

class Derive(object):
    """ Helper class to bind to parameter instances as a method.

    This is a powerful abstract base class for customisting how derived
    parameters are evaluated from references to existing parameter objects.

    :param obj refs:
        Some references to parameter objects or subspace objects that
        need to be stored for deriving variable values. For example,
        a dictionary with key-ref pairs.

    .. note::

        In principle, it might be the case that the calling parameter
        subspace does not have references to other parameter objects
        required, *and* that distinct subspaces require mutual references
        to each other. An example would be two hot regions, each of which
        has one or more parameters that are derived in part from parameters
        of the other hot region. In this case you need to instantiate
        the subspaces first by binding instances of this present
        class to parameters. However, you then need to complete these instances
        of this present class (or more precisely instances of subclasses) with
        the required references. As an example, consider the following:

    .. highlight:: python
    .. code-block:: python

        bounds = dict(super_colatitude = (None, None),
                      super_radius = (None, None),
                      super_temperature = (None, None))

        class derive(xpsi.Derive):
            def __init__(self):
                pass

            def __call__(self, boundto, caller = None):
                # ref is a reference to another hot region object
                return self.ref['phase_shift'] - 0.5

        ref_1 = derive()

        # a simple circular, simply-connected spot
        primary = xpsi.HotRegion(bounds=bounds,
                                    values={'phase_shift': ref_1},
                                    symmetry=True,
                                    omit=False,
                                    cede=False,
                                    concentric=False,
                                    sqrt_num_cells=32,
                                    min_sqrt_num_cells=10,
                                    max_sqrt_num_cells=64,
                                    num_leaves=100,
                                    num_rays=200,
                                    do_fast=False,
                                    prefix='p')

        bounds = dict(omit_colatitude = (None, None),
                      super_radius = (None, None),
                      phase_shift = (None, None),
                      super_temperature = (None, None),
                      omit_radius = (None, None),
                      omit_azimuth = (None, None))

        class derive(xpsi.Derive):
            def __init__(self):
                pass

            def __call__(self, boundto, caller = None):
                return math.pi - self.ref['super_colatitude']

        ref_2 = derive()

        # overlap of an omission region and
        # and a radiating super region
        secondary = xpsi.HotRegion(bounds=bounds,
                                      values={'super_colatitude': ref_2},
                                      symmetry=True,
                                      omit=True,
                                      cede=False,
                                      concentric=False,
                                      sqrt_num_cells=32,
                                      min_sqrt_num_cells=10,
                                      max_sqrt_num_cells=100,
                                      num_leaves=100,
                                      num_rays=200,
                                      do_fast=False,
                                      is_secondary=True,
                                      prefix='s')

        from xpsi import HotRegions

        hot = HotRegions((primary, secondary))

        # the crux: resolve the mutual refs
        ref_1.ref = secondary
        ref_2.ref = primary

    """

    __metaclass__ = ABCMeta

    @abstractmethod
    def __init__(self, refs):
        self.refs = refs

    @abstractmethod
    def __call__(self, boundto, caller = None):
        """ Derive value from some parameters.

        The second argument is the parameter object to which this callable
        is bound. The third argument is a subspace from which the call
        comes, which might be useful or even sufficient for retrieving the
        required information, in which case write an initialiser with ``pass``
        as the body.

        """
        return 0 # calculate something and return

class Parameter(object):
    """ A parameter.

    :param str name:
        A unique parameter name for identification in attribute lookup.

    :param tuple strict_bounds:
        One 2-tuple of hard bounds per parameter. Can be unbounded
        *in principle*, but read the documentation for the
        :class:`~.Prior.Prior` class first.

    :param tuple bounds:
        One 2-tuple of hard bounds per parameter. Can be unbounded
        *in principle*, but read the documentation for the
        :class:`~.Prior.Prior` class first.

    :param bool permit_prepend:
        Allow encapsulating subspaces to prepend the parameter name with
        a prefix? Note that this gives permission recursively to all
        encapsulating subspaces higher in the hierarchy.

    :param bool is_hyperparameter:
        A boolean declaring whether the parameter is a hyperparameter.

    """

    @make_verbose('Creating parameter:')
    def __init__(self, name, strict_bounds, bounds=(None,None),
                 doc=None, symbol=r'', value=None, permit_prepend=True,
                 deactivate_verbosity=False, is_hyperparameter=False):
        """ See the class docstring. """

        self.name = name
        self.strict_bounds = strict_bounds
        self.fixed = True if bounds is None else False
        self.is_hyperparameter = is_hyperparameter
        self.bounds = bounds
        self.doc = doc
        self.symbol = symbol
        if callable(value):
            if not self.fixed:
                raise TypeError('Initial value should be a scalar not callable.')
            if not isinstance(value, Derive):
                raise TypeError('It is recommended to subclass the prototype '
                                'abstract base class ``Derive``.')
            self.evaluate = MethodType(value, self, Parameter)
            self.derived = True
        else:
            self.value = value
            self.derived = False
        self.permit_prepend = permit_prepend

        if self.fixed: # fixed can also encapsulate derived variables
            if callable(value):
                end = 'that is derived from ulterior variables'
            else:
                end = 'with fixed value %.3e' % value
        else:
            bounds = self.bounds # bounds might have been automatically set
            if None in bounds:
                if bounds[0] is not None and bounds[1] is None:
                    bounds = 'lower-bound %.3e' % bounds[0]
                elif bounds[0] is None and bounds[1] is not None:
                    bounds = 'upper-bound %.3e' % bounds[1]
                else:
                    bounds=''
            else:
                bounds  = 'bounds [%.3e, %.3e]' % tuple(bounds)

            if value is None:
                value = ''
            else:
                value = 'initial value %.3e' % value

            if bounds and value:
                end = 'with ' + bounds + ' and ' + value
            elif bounds:
                end = 'with ' + bounds
            elif value:
                end = 'with ' + value
            else:
                end = ''

        yield ('    > Named "%s" %s.' % (name, end) if end
               else  '    > Named "%s".' % name)

        if doc is not None:
            yield '    > %s' % self.doc # get set version

        yield None # initialiser must return NoneType

    @property
    def name(self):
        """ Get the name of the parameter. """
        return self._name

    @name.setter
    def name(self, name):
        if isinstance(name, _six.string_types):
            self._name = name
        else:
            raise TypeError('Name must be a string.')

    @property
    def is_hyperparameter(self):
        """ Is the variable a hyperparameter? """
        return self._is_hyperparameter

    @is_hyperparameter.setter
    def is_hyperparameter(self, is_hyper):
        if not isinstance(is_hyper, bool):
            raise TypeError('A boolean is required to define variable type.')
        self._is_hyperparameter = is_hyper

    @property
    def permit_prepend(self):
        """ Allow subspaces to prepend parameter with prefixes? """
        return self._permit_prepend

    @permit_prepend.setter
    def permit_prepend(self, permit):
        if not isinstance(permit, bool):
            raise TypeError('Provide a boolean to define prepend permissions.')
        self._permit_prepend = permit

    @property
    def doc(self):
        """ Redirect to the magic docstring. """
        return self.__doc__

    @doc.setter
    def doc(self, doc):
        if isinstance(doc, _six.string_types):
            lines = [string.strip(line) for line in doc.splitlines()]
            doc = string.join([line for line in lines if line], '\n')
            if doc[-1] != '.': doc += '.'
            self.__doc__ =  doc
        elif doc is not None:
            raise TypeError('Parameter description must be a string and you '
                            'a description must be provided.')

    @doc.deleter
    def doc(self):
        del self.__doc__

    def __repr__(self):
        """ Get a string summary of the parameter and current value. """
        try:
            val = self.evaluate()
        except (TypeError, AttributeError, NameError):
            msg = ''
        else:
            msg = (' = %.3e' % val if val is not None else '')
        return str(self)[:-1] + msg

    def __str__(self):
        """ Redirect to the magic doctring. """
        return self.__doc__

    @property
    def symbol(self):
        """ Get TeX-compatible symbol."""
        return self._tex

    @symbol.setter
    def symbol(self, symbol):
        if isinstance(symbol, _six.string_types):
            self._tex = symbol
        elif symbol is not None:
            raise TypeError('Invalid type for tex-compatible symbol string.')

    @property
    def strict_bounds(self):
        """ Get the strict bounds of the parameter. """
        return self._strict_bounds

    @strict_bounds.setter
    def strict_bounds(self, bounds):
        try:
            iter(bounds)
        except TypeError:
            raise TypeError('Bounds must be an ordered container with '
                            'two elements.')
        else:
            if len(bounds) != 2:
                raise TypeError('Bounds must be an ordered container with two elements.')
            else:
                if None not in bounds:
                    if bounds[0] >= bounds[1]:
                        raise ValueError('Lower-bound is greater than or equal to upper-bound.')

                if bounds[0] is None: bounds[0] = -_np.inf
                if bounds[1] is None: bounds[1] = _np.inf

                self._strict_bounds = bounds

    @property
    def bounds(self):
        """ Get the hard bounds of the parameter. """
        return self._bounds

    @bounds.setter
    def bounds(self, bounds):
        try:
            iter(bounds)
        except TypeError:
            if bounds is None and self.fixed:
                self._bounds = None
            else:
                raise TypeError('Bounds must be an ordered container with '
                                'two elements if the parameter is free, '
                                'or ``None`` if fixed.')
        else:
            if self.fixed:
                raise TypeError('Check if parameter %s should actually be '
                                'free.' % self._name)
            elif len(bounds) != 2:
                raise TypeError('Bounds must be an ordered container with two elements.')
            else:
                if None not in bounds:
                    if bounds[0] >= bounds[1]:
                        raise ValueError('Lower-bound is greater than or equal to upper-bound.')
                bounds = list(bounds) # make mutable
                for i, bound in enumerate(bounds):
                    if bound is not None:
                        if not self.strict_bounds[0] <= bound <= self.strict_bounds[1]:
                            raise ValueError('Invalid bound for parameter '
                                             'named "%s".' % self.name)
                    else:
                        bounds[i] = self.strict_bounds[i]

                self._bounds = tuple(bounds) # back to immutable

    @property
    def fixed(self):
        """ Is the variable fixed (or derived) or a free parameter? """
        return self._fixed

    @fixed.setter
    def fixed(self, fix):
        if not isinstance(fix, bool):
            raise TypeError('A boolean is required to define variable type.')
        self._fixed = fix

    @property
    def value(self):
        """ Get the current parameter value. """
        return self._value

    @value.setter
    def value(self, value):
        try:
            if not self.strict_bounds[0] <= float(value) <= self.strict_bounds[1]:
                # handle this exception externally if sampling software can
                # make proposals outside of strict bounds
                raise StrictBoundsError('Value of parameter %s is not within the strict bounds.'%self.name)
        except TypeError:
            if self._fixed:
                raise ValueError('Value must be a float.')
            if value is not None:
                print('Value must be a float.')
                raise
            else:
                self._value = None
        else:
            try:
                self._cache() # cache it!
            except AttributeError:
                # first time being set so nothing to cache
                pass
            self._value = float(value)

    def evaluate(self, caller = None):
        """ Symlink to property pending dynamic overwrite.

        :param obj caller:
            An object, such as the calling class itself, used to transfer
            information from higher in the hierarchy.

        Overwrite if value must be explicitly computed from other variables
        and parameters. That is, subclass :class:`~.Derive`, instantiate,
        and pass as a callable to the initialiser of the present class as
        a value of a *derived* parameter. The callable will automatically
        replace the bound variant of this present method and will have access
        to some caller object, plus other references (class and/or instance)
        attributes you define in the subclass.

        """
        return self.value

    @property
    def cached(self):
        """ Get the cached value. """
        try:
            return self._cached
        except AttributeError:
            return None

    @cached.setter
    def cached(self, value):
        """ To clear the cache, use the deleter. """
        try:
            self._cached = float(value)
        except TypeError:
            if value is None:
                self._cached = None
            else:
                raise TypeError('A float is required.')

    @cached.deleter
    def cached(self):
        """ Clear the cache. """
        try:
            del self._cached
        except AttributeError:
            pass # quietly do nothing

    def _cache(self):
        self._cached = self._value

    @property
    def needs_update(self):
        """ Do cached dependencies need to be updated? """
        if self.is_hyperparameter:
            return False # likelihood implicitly dependent on hyperparameters

        if self.derived:
            return True # assume ulterior variables have changed
        elif self.fixed:
            return False

        try:
            return self.cached != self._value
        except AttributeError:
            return True

    def __call__(self, value = None):
        """ Update or get current point if the parameter is *free*.

        :param array-like p:
            New point to update to. If ``None``, the current point is returned.

        :returns: Current point (if the call was not an update).
        :rtype: array-like

        :raises AttributeError:
            If parameter is derived or has no value yet but the argument is
            ``None``.

        """
        if value is not None:
            self.value = value
        else:
            return self.value
