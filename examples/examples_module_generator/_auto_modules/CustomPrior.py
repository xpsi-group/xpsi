""" Prior module for X-PSI CST+PDT modelling of NICER PSR J0030+0451 event data. """

import argparse
import re

class ArgumentParserCustom(argparse.ArgumentParser):
    """A custom implementation of argparse.ArgumentParser for handling arguments specified in a configuration file."""

    def convert_arg_line_to_args(self, arg_line):
        """ Convert a line from a configuration file to a list of arguments.

        :param arg_line (str): Line from the configuration file.
        :return: A list of arguments.
        """
        if (re.match(r'^[\s]*#', arg_line) or   # look for any number of whitespace characters up to a `#` character
            re.match(r'^[\s]*$', arg_line)):    # look for lines containing nothing or just whitespace
            return []
        else:
            try:
                _idx = arg_line.index('#')
            except ValueError:
                pass
            else:
                arg_line = arg_line[:_idx].rstrip()

            return [arg_line]

parser = ArgumentParserCustom(
    description="""
    Prior module for X-PSI CST+PDT modelling of NICER PSR J0030+0451 event data.

    You should import this module.

    For help: python %(prog)s -h

    """,
    fromfile_prefix_chars='@')

class CompileAction(argparse._StoreAction):
    """ Compile arguments for dynamic evaluation. """
    def __call__(self, parser, namespace, values, option_string=None):
        if isinstance(values, list):
            if 'DEFAULT UNIFORM' in values:
                setattr(namespace, self.dest, None)
                return None
        elif values == 'DEFAULT UNIFORM':
            setattr(namespace, self.dest, None)
            return None

        if isinstance(values, list):
            for i, value in enumerate(values[:-1]):
                values[i] = compile(value, '<string>', 'exec')
            values[-1] = compile(values[-1], '<string>', 'eval')
            setattr(namespace, self.dest, values)
        else:
            setattr(namespace, self.dest, compile(values, '<string>', 'eval'))

parser.add_argument('--prior-import-statements',
                    type=str,
                    nargs='*',
                    default=['from scipy.stats import truncnorm', 'import math'],
                    help='Custom import statements needed for evaluation of prior CDFs. Each statement is executed with the ``exec(...)`` builtin function.')

parser.add_argument('--prior-global-statements',
                    type=str,
                    nargs='*',
                    help='Custom assignment statements to be evaluated on the global level that are useful, e.g., for evaluation of prior CDFs. Each statement is executed with the ``exec(...)`` builtin function.')

parser.add_argument('--mass-prior',
                    type=str,
                    nargs='*',
                    action=CompileAction,
                    help='Prior inverse CDF of gravitation mass (solar masses). Supply a function of one variable (the probability mass ``x``), in the form of an expression that can be evaluated with the ``eval(...)`` builtin function, i.e., scipy.stats.truncnorm(x, ...). Note that the prior default PDF is uniform (with compact support), so do not supply a CDF if a uniform prior is desired, or to be explicit, use: DEFAULT UNIFORM. You must use DEFAULT UNIFORM to overwrite a default CDF shown in the auto-generated configuration file, unless the parameter is fixed/derived in which case the prior flag is silently ignored. You can also use the flag more than once: the last usage must be an expression that will be dynamically evaluated using the ``eval(...)`` builtin and must return a float to set as the parameter value; the other usages can be helper statements executed with the ``exec(...)`` builtin, e.g., to set temporary local variables to make the code (and configuration file more readable).')

parser.add_argument('--distance-prior',
                    type=str,
                    nargs='*',
                    action=CompileAction,
                    help='Prior inverse CDF of Earth distance (kpc). Supply a function of one variable (the probability mass ``x``), in the form of an expression that can be evaluated with the ``eval(...)`` builtin function, i.e., scipy.stats.truncnorm(x, ...). Note that the prior default PDF is uniform (with compact support), so do not supply a CDF if a uniform prior is desired, or to be explicit, use: DEFAULT UNIFORM. You must use DEFAULT UNIFORM to overwrite a default CDF shown in the auto-generated configuration file, unless the parameter is fixed/derived in which case the prior flag is silently ignored. You can also use the flag more than once: the last usage must be an expression that will be dynamically evaluated using the ``eval(...)`` builtin and must return a float to set as the parameter value; the other usages can be helper statements executed with the ``exec(...)`` builtin, e.g., to set temporary local variables to make the code (and configuration file more readable).')

parser.add_argument('--cos-inclination-prior',
                    type=str,
                    nargs='*',
                    action=CompileAction,
                    help='Prior inverse CDF of cosine of Earth inclination to stellar spin axis. Supply a function of one variable (the probability mass ``x``), in the form of an expression that can be evaluated with the ``eval(...)`` builtin function, i.e., scipy.stats.truncnorm(x, ...). Note that the prior default PDF is uniform (with compact support), so do not supply a CDF if a uniform prior is desired, or to be explicit, use: DEFAULT UNIFORM. You must use DEFAULT UNIFORM to overwrite a default CDF shown in the auto-generated configuration file, unless the parameter is fixed/derived in which case the prior flag is silently ignored. You can also use the flag more than once: the last usage must be an expression that will be dynamically evaluated using the ``eval(...)`` builtin and must return a float to set as the parameter value; the other usages can be helper statements executed with the ``exec(...)`` builtin, e.g., to set temporary local variables to make the code (and configuration file more readable).')

parser.add_argument('--neutral-hydrogen-column-density-prior',
                    type=str,
                    nargs='*',
                    action=CompileAction,
                    help='Prior inverse CDF of ratio of interstellar neutral hydrogen column density to the fiducial density. Supply a function of one variable (the probability mass ``x``), in the form of an expression that can be evaluated with the ``eval(...)`` builtin function, i.e., scipy.stats.truncnorm(x, ...). Note that the prior default PDF is uniform (with compact support), so do not supply a CDF if a uniform prior is desired, or to be explicit, use: DEFAULT UNIFORM. You must use DEFAULT UNIFORM to overwrite a default CDF shown in the auto-generated configuration file, unless the parameter is fixed/derived in which case the prior flag is silently ignored. You can also use the flag more than once: the last usage must be an expression that will be dynamically evaluated using the ``eval(...)`` builtin and must return a float to set as the parameter value; the other usages can be helper statements executed with the ``exec(...)`` builtin, e.g., to set temporary local variables to make the code (and configuration file more readable).')

parser.add_argument('--p-super-temperature-prior',
                    type=str,
                    nargs='*',
                    action=CompileAction,
                    help='Prior inverse CDF of hot-region p superseding region log10(temperature [K]). Supply a function of one variable (the probability mass ``x``), in the form of an expression that can be evaluated with the ``eval(...)`` builtin function, i.e., scipy.stats.truncnorm(x, ...). Note that the prior default PDF is uniform (with compact support), so do not supply a CDF if a uniform prior is desired, or to be explicit, use: DEFAULT UNIFORM. You must use DEFAULT UNIFORM to overwrite a default CDF shown in the auto-generated configuration file, unless the parameter is fixed/derived in which case the prior flag is silently ignored. You can also use the flag more than once: the last usage must be an expression that will be dynamically evaluated using the ``eval(...)`` builtin and must return a float to set as the parameter value; the other usages can be helper statements executed with the ``exec(...)`` builtin, e.g., to set temporary local variables to make the code (and configuration file more readable).')

parser.add_argument('--XTI-energy-independent-effective-area-scaling-factor-prior',
                    type=str,
                    nargs='*',
                    action=CompileAction,
                    default=[compile('truncnorm.ppf(x, -5.0, 5.0, loc=1.0, scale=0.104)', '<string>', 'eval')],
                    help='Prior inverse CDF of the energy-independent effective area scaling factor. Supply a function of one variable (the probability mass ``x``), in the form of an expression that can be evaluated with the ``eval(...)`` builtin function, i.e., scipy.stats.truncnorm(x, ...). Note that the prior default PDF is uniform (with compact support), so do not supply a CDF if a uniform prior is desired, or to be explicit, use: DEFAULT UNIFORM. You must use DEFAULT UNIFORM to overwrite a default CDF shown in the auto-generated configuration file, unless the parameter is fixed/derived in which case the prior flag is silently ignored. You can also use the flag more than once: the last usage must be an expression that will be dynamically evaluated using the ``eval(...)`` builtin and must return a float to set as the parameter value; the other usages can be helper statements executed with the ``exec(...)`` builtin, e.g., to set temporary local variables to make the code (and configuration file more readable).')

parser.add_argument('--s-super-temperature-prior',
                    type=str,
                    nargs='*',
                    action=CompileAction,
                    help='Prior inverse CDF of hot-region s superseding region log10(temperature [K]). Supply a function of one variable (the probability mass ``x``), in the form of an expression that can be evaluated with the ``eval(...)`` builtin function, i.e., scipy.stats.truncnorm(x, ...). Note that the prior default PDF is uniform (with compact support), so do not supply a CDF if a uniform prior is desired, or to be explicit, use: DEFAULT UNIFORM. You must use DEFAULT UNIFORM to overwrite a default CDF shown in the auto-generated configuration file, unless the parameter is fixed/derived in which case the prior flag is silently ignored. You can also use the flag more than once: the last usage must be an expression that will be dynamically evaluated using the ``eval(...)`` builtin and must return a float to set as the parameter value; the other usages can be helper statements executed with the ``exec(...)`` builtin, e.g., to set temporary local variables to make the code (and configuration file more readable).')

parser.add_argument('--s-cede-temperature-prior',
                    type=str,
                    nargs='*',
                    action=CompileAction,
                    help='Prior inverse CDF of hot-region s ceding region log10(temperature [K]). Supply a function of one variable (the probability mass ``x``), in the form of an expression that can be evaluated with the ``eval(...)`` builtin function, i.e., scipy.stats.truncnorm(x, ...). Note that the prior default PDF is uniform (with compact support), so do not supply a CDF if a uniform prior is desired, or to be explicit, use: DEFAULT UNIFORM. You must use DEFAULT UNIFORM to overwrite a default CDF shown in the auto-generated configuration file, unless the parameter is fixed/derived in which case the prior flag is silently ignored. You can also use the flag more than once: the last usage must be an expression that will be dynamically evaluated using the ``eval(...)`` builtin and must return a float to set as the parameter value; the other usages can be helper statements executed with the ``exec(...)`` builtin, e.g., to set temporary local variables to make the code (and configuration file more readable).')

if __name__ == '__main__':
    args, _ = parser.parse_known_args()
else:
    args, _ = parser.parse_known_args(['@./config.ini'])

import numpy as np
import math

import xpsi
from xpsi.global_imports import _2pi, gravradius, _dpr
from xpsi import Parameter

from xpsi.cellmesh.mesh_tools import eval_cedeCentreCoords as eval_coords_under_rotation
from scipy.interpolate import Akima1DInterpolator

if args.prior_import_statements is not None:
    for import_statement in args.prior_import_statements:
        exec(import_statement)

if args.prior_global_statements is not None:
    for global_statement in args.prior_global_statements:
        exec(global_statement)

class CustomPrior(xpsi.Prior):
    """ A joint prior PDF. """

    __derived_names__ = ['compactness']

    __draws_from_support__ = 4

    def __init__(self):

        super(CustomPrior, self).__init__()

        self.a_f = 0.001
        self.b_f = 1.0
        self.a_zeta = 0.001
        self.b_zeta = math.pi/2.0 - self.a_zeta

        vals = np.linspace(0.0, self.b_zeta, 1000)
        self._interpolator_super_smaller = Akima1DInterpolator(self._vector_super_smaller_radius_mass(vals), vals)
        self._interpolator_super_smaller.extrapolate = True

        self.c_f = 0.001
        self.d_f = 2.0
        self.a_xi = 0.001
        self.b_xi = math.pi/2.0 - self.a_xi

        vals = np.linspace(0.0, self.b_xi, 1000)
        self._interpolator = Akima1DInterpolator(self._vector_super_radius_mass(vals), vals)
        self._interpolator.extrapolate = True

    def __call__(self, p = None):
        """ Evaluate distribution at point ``p``.

        :param list p: Model parameter values.

        :returns: Logarithm of the distribution evaluated at point ``p``.

        """
        temp = super(CustomPrior, self).__call__(p)
        if not np.isfinite(temp):
            return temp

        ref = self.parameters.star.spacetime # shortcut

        # check the prior PDF support conditions below and comment out the exception throw
        # if you want to condition on those model assumptions
        #raise NotImplementedError('Implement the prior __call__ method.')

        # based on contemporary EOS theory
        if not ref['radius'] <= 16.0:
            return -np.inf

        # limit polar radius to be just outside the Schwarzschild photon sphere
        R_p = 1.0 + ref.epsilon * (-0.788 + 1.030 * ref.zeta)
        if R_p < 1.505 / ref.R_r_s:
            return -np.inf

        mu = math.sqrt(-1.0 / (3.0 * ref.epsilon * (-0.788 + 1.030 * ref.zeta)))
        # 2-surface cross-section have a single maximum in |z|
        # i.e., an elliptical surface; minor effect on support, if any,
        # only for high spin frequencies
        if mu < 1.0:
            return -np.inf

        # check effective gravity at pole (where it is maximum) and
        # at equator (where it is minimum) are in NSX limits
        grav = xpsi.surface_radiation_field.effective_gravity(np.array([1.0, 0.0]),
                                                              np.array([ref.R] * 2 ),
                                                              np.array([ref.zeta] * 2),
                                                              np.array([ref.epsilon] * 2))

        for g in grav:
            if not 13.7 <= g <= 15.0: # check that these NSX effective gravity table limits are correct
                return -np.inf

        ref = self.parameters # redefine shortcut

        # require that hot regions do not overlap within the prior support
        if self._overlap(ref, 's', 'p', 'super', 'super', 0.0, 0.0):
            if not self._overlap(ref, 's', 'p', 'super', 'omit', 0.0, 0.0, superset='p'):
                return -np.inf

        if self._overlap(ref, 's', 'p', 'cede', 'super', 0.0, 0.0):
            if not self._overlap(ref, 's', 'p', 'cede', 'omit', 0.0, 0.0, superset='p'):
                return -np.inf

        return 0.0

    @staticmethod
    def _colatitude(ref, z, z_member):
        """ Helper to bypass exception for concentric regions. """
        try:
            return ref[z + '__' + z_member + '_colatitude']
        except KeyError:
            return ref[z + '__super_colatitude']

    @staticmethod
    def _azimuth(ref, z, z_member):
        """ Helper to bypass exception for concentric regions. """
        try:
            return ref[z + '__' + z_member + '_azimuth']
        except KeyError:
            return 0.0

    def _overlap(self, ref, x, y, x_member, y_member, x_antiphase, y_antiphase, superset=None, use_cached=False):
        """ Determine overlap between two spherical circles. """

        if not use_cached:
            _tmp_phase_x = (ref[x + '__phase_shift'] + x_antiphase) * _2pi
            if x_member == 'super':
                _tmp_phase_x -= self._azimuth(ref, x, 'omit')
            elif x_member == 'cede':
                _tmp_phase_x += self._azimuth(ref, x, 'cede')

            _tmp_phase_y = (ref[y + '__phase_shift'] + y_antiphase) * _2pi
            if y_member == 'super':
                _tmp_phase_y -= self._azimuth(ref, y, 'omit')
            elif y_member == 'cede':
                _tmp_phase_y += self._azimuth(ref, y, 'cede')

            _phi = _tmp_phase_x - _tmp_phase_y

            self._ang_sep = xpsi.HotRegion.psi(self._colatitude(ref, x, x_member),
                                               _phi,
                                               self._colatitude(ref, y, y_member))

        if superset is None:
            if self._ang_sep < ref[x + '__' + x_member + '_radius'] + ref[y + '__' + y_member + '_radius']:
                return True

        elif superset == y:
            if self._ang_sep + ref[x + '__' + x_member + '_radius'] < ref[y + '__' + y_member + '_radius']:
                return True

        elif superset == x:
            if ref[x + '__' + x_member + '_radius'] > self._ang_sep + ref[y + '__' + y_member + '_radius']:
                return True

        return False

    def _I_super_smaller(self, x):
        return x * np.log(self.b_zeta/self.a_zeta)

    def _II_super_smaller(self, x):
        return x - self.a_zeta - x*np.log(x/self.b_zeta)

    def _scalar_super_smaller_radius_mass(self, x):
        if x >= self.a_zeta:
            mass = self._II_super_smaller(x)
        else:
            mass = self._I_super_smaller(x)

        return mass

    def _vector_super_smaller_radius_mass(self, x):
        masses = np.zeros(len(x))

        for i, _ in enumerate(x):
            masses[i] = self._scalar_super_smaller_radius_mass(_)

        masses /= (self.b_f - self.a_f)
        masses /= (self.b_zeta - self.a_zeta)

        return masses

    def _inverse_sample_cede_larger_radius(self, x, psi):
        if psi < self.a_zeta:
            return self.a_zeta*np.exp(x * np.log(self.b_zeta/self.a_zeta))
        else:
            return psi*np.exp(x*np.log(self.b_zeta/psi))

    def _I(self, x):
        return x * np.log(self.b_xi/self.a_xi)

    def _II(self, x):
        return 2.0*(x - self.a_xi) - x*np.log(x/self.b_xi)

    def _scalar_super_radius_mass(self, x):
        if x >= self.a_xi:
            mass = self._II(x)
        else:
            mass = self._I(x)

        return mass

    def _vector_super_radius_mass(self, x):
        masses = np.zeros(len(x))

        for i, _ in enumerate(x):
            masses[i] = self._scalar_super_radius_mass(_)

        masses /= (self.d_f - self.c_f)
        masses /= (self.b_xi - self.a_xi)

        return masses

    def _inverse_sample_cede_radius(self, x, psi):
        if psi < self.a_xi:
            return self.a_xi*np.exp(x * np.log(self.b_xi/self.a_xi))
        elif psi >= self.a_xi and x <= 1.0/(1.0 + np.log(self.b_xi/psi)):
            return x*psi*(1.0 + np.log(self.b_xi/psi))
        else:
            return psi*np.exp(x*(1.0 + np.log(self.b_xi/psi)) - 1.0)

    def inverse_sample(self, hypercube=None):
        """ Draw sample uniformly from the distribution via inverse sampling. """

        global args

        to_cache = self.parameters.vector

        if hypercube is None:
            hypercube = np.random.rand(len(self))

        # the base method is useful, so to avoid writing that code again:
        _ = super(CustomPrior, self).inverse_sample(hypercube)

        ref = parameters = self.parameters # redefine shortcut

        try:
            self._modded_names
        except AttributeError:
            self._modded_names = [name.replace('__', '_') for name in ref.names]

        for modded_name, name in list(zip(self._modded_names, ref.names)):
            if getattr(args, modded_name + '_prior', None) is not None:
                idx = ref.index(name)
                x = hypercube[idx]
                for _statement in getattr(args, modded_name + '_prior')[:-1]:
                    exec(_statement)
                ref[name] = eval(getattr(args, modded_name + '_prior')[-1])

        # inverse sample parameters of hot-region p
        idx = ref.index('p__super_colatitude')
        a, b = ref.get_param('p__super_colatitude').bounds
        a = math.cos(a); b = math.cos(b)
        ref['p__super_colatitude'] = math.acos(b + (a - b) * hypercube[idx])

        # radius of superseding region (omit or super code object)
        idx = ref.index('p__omit_radius')
        ref['p__omit_radius'] = float(self._interpolator_super_smaller(hypercube[idx]))

        # radius of ceding region (super or cede code object)
        idx = ref.index('p__super_radius')
        ref['p__super_radius'] = self._inverse_sample_cede_larger_radius(hypercube[idx], ref['p__omit_radius'])

        # inverse sample parameters of hot-region s
        idx = ref.index('s__super_colatitude')
        a, b = ref.get_param('s__super_colatitude').bounds
        a = math.cos(a); b = math.cos(b)
        ref['s__super_colatitude'] = math.acos(b + (a - b) * hypercube[idx])

        # radius of superseding region (omit or super code object)
        idx = ref.index('s__super_radius')
        ref['s__super_radius'] = float(self._interpolator(hypercube[idx]))

        # radius of ceding region (super or cede code object)
        idx = ref.index('s__cede_radius')
        ref['s__cede_radius'] = self._inverse_sample_cede_radius(hypercube[idx], ref['s__super_radius'])

        # coordinates of mask or ceding region (omit or cede code object)
        idx = ref.index('s__cede_colatitude')
        if ref['s__super_radius'] <= ref['s__cede_radius']:
            ref[idx] = hypercube[idx] * (ref['s__cede_radius'] + ref['s__super_radius'])
        else:
            ref[idx] = (  ref['s__super_radius']
                        - ref['s__cede_radius']
                        + 2.0*hypercube[idx]*ref['s__cede_radius'] )

        ref[idx], ref['s__cede_azimuth'] = eval_coords_under_rotation(ref['s__super_colatitude'],
                                                                       ref['s__cede_colatitude'],
                                                                       ref['s__cede_azimuth'])

        # restore proper cache
        for parameter, cache in list(zip(self.parameters, to_cache)):
            parameter.cached = cache

        # it is important that we return the desired vector because it is
        # automatically written to disk by MultiNest and only by MultiNest
        return self.parameters.vector

    def transform(self, p, **kwargs):
        """ A transformation for post-processing. """

        p = list(p) # copy

        # used ordered names and values
        ref = dict(list(zip(self.parameters.names, p)))

        # compactness ratio M/R_eq
        p += [gravradius(ref['mass']) / ref['radius']]

        return p