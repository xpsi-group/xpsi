""" Main module for NICER PSR J0030+0451 <- X-PSI 2.1.1 CST+PDT"""
import os
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

            if xpsi._verbose:
                print(arg_line)
            return [arg_line]

    def add_argument(self, *args, **kwargs):
        """
        Add an argument to the argument parser.
        """
        if kwargs.pop('destined_for_config_file', True) and args[0] != '-h':
            _ = (args[0],
                 kwargs.get('default', None),
                 kwargs.get('nargs', 1) if kwargs.get('action') != 'store_true' else 0,
                 kwargs.pop('comment_line_above', None),
                 kwargs.pop('empty_lines_below', 0),
                 kwargs.pop('comment', False),
                 kwargs.pop('inline_comment', None),
                 kwargs.get('action', None))
            try:
                self._config_file_args
            except AttributeError:
                self._config_file_args = [_]
            else:
                self._config_file_args.append(_)
        else:
            _ = kwargs.pop('comment_line_above', None)
            _ = kwargs.pop('empty_lines_below', 0)
            _ = kwargs.pop('comment', False)
            _ = kwargs.pop('inline_comment', None)

        super(ArgumentParserCustom, self).add_argument(*args, **kwargs)

class CompileAction(argparse._StoreAction):
    """ Compile arguments for dynamic evaluation. """
    def __call__(self, parser, namespace, values, option_string=None):
        if isinstance(values, list):
            for i, value in enumerate(values):
                values[i] = compile(value, '<string>', 'eval')
            setattr(namespace, self.dest, values)
        else:
            setattr(namespace, self.dest, compile(values, '<string>', 'eval'))

class GenerateConfigAction(argparse.Action):
    """ Class that generates a configuration file based on the arguments provided to argparse.

    The class inherits from argparse.Action and overrides the __init__ and __call__ methods to
    implement the configuration file generation.
    """
    def __init__(self, option_strings, dest, **kwargs):
        """ Initialize the class instance.

        :param option_strings (list): A list of command-line option strings.
        :param dest (str): The name of the attribute to be added to the namespace.
        """
        super(GenerateConfigAction, self).__init__(option_strings, dest, nargs=0, **kwargs)

    @staticmethod
    def _typeset(arg, default, nargs, comment_line_above, empty_lines_below, comment, inline_comment, action, newline=True):
        """ Helper method to generate the text for a single argument in the configuration file.

        :param arg (str): The name of the argument.
        :param default (str, list): The default value for the argument.
        :param nargs (int, str): The number of values to take as input for the argument.
        :param comment_line_above (str): Text to place above the argument. Enter 'rule' if you just want to place a
                                         separating line above, or enter a header text describing a group of arguments.
        :param empty_lines_below (int): The number of empty lines to include below the argument.
        :param comment (bool): Whether the argument should be commented out or not.
        :param inline_comment (str): A comment to include next to the argument.
        :param action (str): The action to be performed with the argument.
        :param newline (bool): Whether to include a newline before the argument.

        :return str: The text for the argument in the configuration file.
        """
        entry = '\n' if newline else ''

        if comment_line_above is not None:
            if comment_line_above == 'rule':
               entry += '#' + '-'*78 + '#\n'
            else:
                _ = '## {} ##'.format(comment_line_above)
                entry += '##' + '-' * (len(_) - 4) + '##\n' + _ + '\n##' + '-' * (len(_) - 4) + '##\n'

        _ = ' ## {}'.format(inline_comment) if inline_comment is not None else ''
        if not _ and isinstance(nargs, int) and nargs > 1:
            _ = ' ## enter {} values below, one per empty line'.format(nargs)
        elif not _ and not isinstance(nargs, int):
            _ = ' ## enter code below, one statement per line'

        if isinstance(default, list):
            for i, _default in enumerate(default):
                entry += '{5}{0}{1}{2}{3}{4}'.format('#' if comment else '',
                                                                 '' if (i > 0 and nargs != 1) else arg,
                                                                 '' if nargs != 1 else '=',
                                                                 '' if (i == 0 and nargs != 1) else str(_default),
                                                                 (_ if i == 0 else '') + ('\n{0}{1}'.format('#' if comment else '', str(_default)) if i == 0 and nargs != 1 else ''),
                                                                 '\n' if i > 0 else '')
        else:
            entry += '{0}{1}{2}{3}{4}'.format('#' if comment else '',
                                                        arg,
                                                        '=' if nargs == 1 else '',
                                                        _ if nargs != 1 else (str(default) if default is not None else ''),
                                                        ('\n' + str(default) if default is not None else '') if nargs != 1 else _)

        if action == 'append':
            entry += '\n#{0}='.format(arg)

        if isinstance(nargs, int) and nargs > 1:
            entry += '\n' * nargs
        elif isinstance(nargs, str):
            entry += '\n' * 3

        entry += '\n' * empty_lines_below

        return entry

    def __call__(self, parser, namespace, values, option_string=None):
        """Method that generates the configuration file.

            :param parser (argparse.ArgumentParser): The ArgumentParser object.
            :param namespace (argparse.Namespace): The Namespace object.
            :param values (list): The values for the arguments.
            :param option_string (str, optional): The option string for the argument.

            :returns None:
        """
        for _ in parser._config_file_args:
            try:
                config_file
            except NameError:
                config_file = self._typeset(*_, newline=False)
            else:
                config_file += self._typeset(*_)

        with open('./config.ini', 'w') as file:
            file.write(config_file)

        print('Configuration file generated.')
        parser.exit()

class NullAction(argparse.Action):
    """ Do not store value in namespace. """
    def __call__(self, parser, namespace, values, option_string=None):
        pass

parser = ArgumentParserCustom(
    description="""
    Main module for X-PSI CST+PDT modelling of NICER PSR J0030+0451 event data.

    You can run this module as a script and launch a sampler, optionally
    with a world of MPI processes.

    Alternate usage: mpiexec -n 4 python -m mpi4py %(prog)s [-h] @<config.ini>

    """,
    fromfile_prefix_chars='@')

def str_to_bool(x):
    if x == 'False':
        return False
    elif x == 'True':
        return True

    raise ValueError('Invalid argument where boolean ``True`` or ``False`` is required.')


parser.add_argument('--generate-config-file', default=argparse.SUPPRESS, action=GenerateConfigAction, help='Generate the configuration file template.',
                         destined_for_config_file=False)

parser.add_argument('--main-import-statements',
                    type=str,
                    nargs='*',
                    default=['from xpsi.global_imports import gravradius', 'import math'],
                    help='Custom import statements needed for main module. Each statement is executed with the ``exec(...)`` builtin function. Note that if you pass statements, the default statements are deleted unless you uncomment the defaults in the configuration file.',
                    comment=True,
                    comment_line_above='import statements needed for main module',
                    empty_lines_below=2,
                    inline_comment='e.g., from ... import ... as ...')

parser.add_argument('--main-global-statements',
                    type=str,
                    nargs='*',
                    help='Custom assignment statements to be evaluated on the global level in the main module. Each statement is executed with the ``exec(...)`` builtin function. Note that if you pass statements, the default statements are deleted unless you uncomment the defaults in the configuration file.',
                    comment=True,
                    comment_line_above='global statements needed for main module',
                    empty_lines_below=2,
                    inline_comment='e.g., global_variable = math.pi')

parser.add_argument('--prior-import-statements',
                    type=str,
                    nargs='*',
                    default=['from scipy.stats import truncnorm', 'import math'],
                    action=NullAction,
                    help='Custom import statements needed for evaluation of prior CDFs. Each statement is executed with the ``exec(...)`` builtin function. Note that if you pass statements, the default statements are deleted unless you uncomment the defaults in the configuration file.',
                    comment=True,
                    comment_line_above='import statements needed for prior',
                    empty_lines_below=2,
                    inline_comment='e.g., from ... import ... as ...')

parser.add_argument('--prior-global-statements',
                    type=str,
                    nargs='*',
                    action=NullAction,
                    help='Custom assignment statements to be evaluated on the global level that are useful, e.g., for evaluation of prior CDFs. Each statement is executed with the ``exec(...)`` builtin function. Note that if you pass statements, the default statements are deleted unless you uncomment the defaults in the configuration file.',
                    comment=True,
                    comment_line_above='global statements needed for prior',
                    empty_lines_below=2,
                    inline_comment='e.g., global_variable = math.pi')

parser.add_argument('--XTI-exposure-time',
                    type=float,
                    help='XTI exposure time in seconds.',
                    comment_line_above='XTI configuration flags')

parser.add_argument('--XTI-count-matrix-path', type=str, help='Absolute or relative path to XTI channel-phase count matrix. If the data is a spectrum (phase-averaged), then the file must contain a vector. This path is written to if the file does not exist by processing the event files.')
parser.add_argument('--XTI-count-matrix-type', type=str, default='double', help='XTI count matrix NumPy data type.',
                    comment=True)
parser.add_argument('--XTI-event-path', type=str, help='Absolute or relative path to XTI event list file.')
parser.add_argument('--XTI-number-phase-bins', type=int, help='Number of phases bins for binning XTI event list file.')
parser.add_argument('--XTI-event-file-channel-column', type=int, default=1, help='Channel column in XTI event list file.',
                    comment=True)
parser.add_argument('--XTI-event-file-phase-column', type=int, default=2, help='Phase column in XTI event list file.',
                    comment=True)
parser.add_argument('--XTI-event-file-skiprows', type=int, default=3, help='Number of top rows to skip when loading XTI event list file.',
                    comment=True)
parser.add_argument('--XTI-events-in-eV', action='store_true', help='XTI event list file lists events by energy in eV?',
                    comment=True)
parser.add_argument('--XTI-arf-path', type=str, help='Absolute or relative path to XTI ARF file.',
                    comment_line_above='rule')
parser.add_argument('--XTI-effective-area-scaling-factor', type=str, default='1.0',
                    help='Factor by which to scale the nominal effective area model, as a mathematical expression, e.g., a ratio of integers such as 51.0/52.0.',
                    comment=True)
parser.add_argument('--XTI-arf-skiprows', type=int, default=3,
                    help='Number of header rows to skip when loading ARF file.',
                    comment=True)
parser.add_argument('--XTI-arf-low-column', type=int, default=1,
                    help='Column (zero-indexed) containing the low energy edges in the ARF file.',
                    comment=True)
parser.add_argument('--XTI-arf-high-column', type=int, default=2,
                    help='Column (zero-indexed) containing the high energy edges in the ARF file.',
                    comment=True)
parser.add_argument('--XTI-arf-area-column', type=int, default=3,
                    help='Column (zero-indexed) containing the effective area in the ARF file.',
                    comment=True)

parser.add_argument('--XTI-rmf-path', type=str, help='Absolute or relative path to XTI RMF file.',
                    comment_line_above='rule')
parser.add_argument('--XTI-rmf-skiprows', type=int, default=3,
                    help='Number of header rows to skip when loading RMF file.',
                    comment=True)
parser.add_argument('--XTI-rmf-usecol', type=int, default=-1,
                    help='Column (zero-indexed) containing the flattened redistribution matrix elements in the RMF file.',
                    comment=True)

parser.add_argument('--XTI-channels-path', type=str, help='Absolute or relative path to XTI channel-energy mapping file.',
                    comment_line_above='rule')
parser.add_argument('--XTI-channel-energies-skiprows', type=int, default=0,
                    help='Number of header rows to skip when loading channel-energy mapping file.',
                    comment=True)
parser.add_argument('--XTI-channel-energies-low-column', type=int, default=0,
                    help='Column (zero-indexed) containing the low energy edges in the channel-energy mapping file.',
                    comment=True)

parser.add_argument('--XTI-input-bounds',
                    type=int,
                    nargs=2,
                    help='XTI bounding input energy intervals of instrument response submatrix for use with NumPy slice notation.',
                    comment_line_above='rule')

parser.add_argument('--XTI-channel-bounds',
                    type=int,
                    nargs=2,
                    help='XTI bounding channels of instrument response submatrix for use with NumPy slice notation.')

parser.add_argument('--XTI-energy-independent-effective-area-scaling-factor-bounds',
                    type=str,
                    nargs=2,
                    action=CompileAction,
                    help='Bounds for XTI energy-independent effective area scaling factor parameter. If no bounds are given (``None``), and no value is given (``None``), the parameter value is fixed at unity, and the instrument response model is locked to the nominal response model (unless a custom model is implemented).',
                    comment=True,
                    comment_line_above='rule')

parser.add_argument('--XTI-energy-independent-effective-area-scaling-factor-prior',
                    type=str,
                    nargs='*',
                    default=['truncnorm.ppf(x, -3.0, 3.0, loc=1.0, scale=0.1)'],
                    action=NullAction,
                    help='Prior inverse CDF of the energy-independent effective area scaling factor. Supply a function of one variable (the probability mass ``x``), in the form of an expression that can be evaluated with the ``eval(...)`` builtin function, i.e., scipy.stats.truncnorm(x, ...). Note that the prior default PDF is uniform (with compact support), so do not supply a CDF if a uniform prior is desired, or to be explicit, use: DEFAULT UNIFORM. You must use DEFAULT UNIFORM to overwrite a default CDF shown in the auto-generated configuration file, unless the parameter is fixed/derived in which case the prior flag is silently ignored. You can also use the flag more than once: the last usage must be an expression that will be dynamically evaluated using the ``eval(...)`` builtin and must return a float to set as the parameter value; the other usages can be helper statements executed with the ``exec(...)`` builtin, e.g., to set temporary local variables to make the code (and configuration file more readable).',
                    comment=True,
                    inline_comment='Normal distribution with std. dev. 10%, truncated at +/- 3 std. dev.')

parser.add_argument('--XTI-energy-independent-effective-area-scaling-factor-value',
                    type=str,
                    action=CompileAction,
                    help='Value for XTI energy-independent effective area scaling parameter. Either the name of an instrument to share the parameter with, as a string, or a float. No value means the parameter is free, whilst a value means the parameter is fixed (and the prior may need to be modified manually). If you want the parameter to be derived from other parameters in a complex way, manual modification of the main module is necessary, but we support functions of one parameter in the configuration file. If the parameter is derived from one other parameter, e.g., the temperature of a hot region is derived from the temperature of the other hot region, then the value needs to be written using the following template: lambda x: f(x), "parameter", "space". In this template: f(x) is a function of the parameter x from which the value is derived; "parameter" is the name of the parameter x as a string; and "space" is an object in the global namespace, with name written as a string, which is a (sub)space of parameters from which the current value of parameter x can be accessed via getitem magic using "parameter".',
                    comment=True,
                    empty_lines_below=2)

parser.add_argument('--XTI-phase-shift-bounds',
                    type=str,
                    nargs=2,
                    action=CompileAction,
                    help='Bounds for XTI phase-shift parameter. If no bounds are given (``None``), and no value is given (``None``), the parameter value is fixed at zero, and is therefore locked to the phase of the signal specified by the hot region phases. For one phase-resolving instrument, this default behaviour is advised, and additional phase-resolving instruments can in principle have a different fixed, derived, or free phase-shift parameter. For instruments that phase-average, the phase-shift can be arbitrarily fixed or derived, but not free because the likelihood is not a function of it.',
                    comment=True,
                    comment_line_above='rule')

parser.add_argument('--XTI-phase-shift-value',
                    type=str,
                    action=CompileAction,
                    help='Value for XTI phase-shift parameter. Either the name of an instrument to share the parameter with, as a string, or a float. No value means the parameter is free, whilst a value means the parameter is fixed (and the prior may need to be modified manually). If you want the parameter to be derived from other parameters in a complex way, manual modification of the main module is necessary, but we support functions of one parameter in the configuration file. If the parameter is derived from one other parameter, e.g., the temperature of a hot region is derived from the temperature of the other hot region, then the value needs to be written using the following template: lambda x: f(x), "parameter", "space". In this template: f(x) is a function of the parameter x from which the value is derived; "parameter" is the name of the parameter x as a string; and "space" is an object in the global namespace, with name written as a string, which is a (sub)space of parameters from which the current value of parameter x can be accessed via getitem magic using "parameter".',
                    comment=True,
                    empty_lines_below=2)

parser.add_argument('--XTI-background-prior-support-path', type=str, help='Absolute or relative path to XTI background prior support file. The channel-by-channel lower count-rate limits in the zeroth column, and the upper count-rate limits in the first column. The channels must already match the data.',
                    comment=True,
                    comment_line_above='rule')
parser.add_argument('--XTI-background-skiprows', type=int, default=0, help='Number of top rows to skip when loading XTI background file (prior support or spectrum file).',
                    comment=True)

parser.add_argument('--XTI-background-path', type=str, help='Absolute or relative path to XTI background spectrum file (for imaging telescope).',
                    comment=True)
parser.add_argument('--XTI-background-usecol', type=int, help='Column to use when loading XTI background spectrum file (for imaging telescope).',
                    comment=True)
parser.add_argument('--XTI-background-prior-support-half-width', type=float, help='XTI background prior support half-width (for imaging telescope). The half-width is in units of standard deviation of background count number per instrument channel.',
                    comment=True)
parser.add_argument('--XTI-background-exposure-time',
                    type=float,
                    help='XTI background exposure time in seconds (for imaging telescope).',
                    comment=True)
parser.add_argument('--XTI-background-scaling-factor',
                    type=str,
                    help='XTI background scaling factor, nominally the ratio of on-source CCD extraction area to background CCD extraction area (ideally on same CCD) for imaging telescope. Supply an expression for evaluation by the ``eval(...)`` builtin function.',
                    comment=True,
                    empty_lines_below=2)

parser.add_argument('--attenuation-path', type=str, help='Absolute or relative path to attenuation file.',
                         comment_line_above='attenuation flags')
parser.add_argument('--attenuation-energy-column', type=int, default=0,
                    help='Column (zero-indexed) containing the energies in the attenuation file.',
                    comment=True)
parser.add_argument('--attenuation-column', type=int, default=1,
                    help='Column (zero-indexed) containing the attenuation factors in the attenuation file.',
                    comment=True)

parser.add_argument('--neutral-hydrogen-column-density-bounds',
                    type=str,
                    nargs=2,
                    action=CompileAction,
                    help='Bounds of the neutral hydrogen column density parameter. If no bounds are given (``None``), and no value is given (``None``), bounds default to the source code strict bounds.',
                    comment=True,
                    comment_line_above='rule')

parser.add_argument('--neural-hydrogen-column-density-prior',
                    type=str,
                    nargs='*',
                    action=NullAction,
                    help='Prior inverse CDF of ratio of interstellar neutral hydrogen column density to the fiducial density. Supply a function of one variable (the probability mass ``x``), in the form of an expression that can be evaluated with the ``eval(...)`` builtin function, i.e., scipy.stats.truncnorm(x, ...). Note that the prior default PDF is uniform (with compact support), so do not supply a CDF if a uniform prior is desired, or to be explicit, use: DEFAULT UNIFORM. You must use DEFAULT UNIFORM to overwrite a default CDF shown in the auto-generated configuration file, unless the parameter is fixed/derived in which case the prior flag is silently ignored. You can also use the flag more than once: the last usage must be an expression that will be dynamically evaluated using the ``eval(...)`` builtin and must return a float to set as the parameter value; the other usages can be helper statements executed with the ``exec(...)`` builtin, e.g., to set temporary local variables to make the code (and configuration file more readable).',
                    comment=True)

parser.add_argument('--neutral-hydrogen-column-density-value',
                    type=str,
                    action=CompileAction,
                    help='Value of the neutral hydrogen column density parameter. No value means the parameter is free, whilst a value means the parameter is fixed (and the prior may need to be modified manually). If you want the parameter to be derived from other parameters in a complex way, manual modification of the main module is necessary, but we support functions of one parameter in the configuration file.',
                    comment=True,
                    empty_lines_below=2)


parser.add_argument('--mass-bounds',
                    type=str,
                    nargs=2,
                    action=CompileAction,
                    help='Bounds of gravitational mass (solar masses). If no bounds are given (``None``), and no value is given (``None``), bounds default to the source code strict bounds.',
                    comment=True,
                    comment_line_above='spacetime flags')

parser.add_argument('--mass-prior',
                    type=str,
                    nargs='*',
                    action=NullAction,
                    help='Prior inverse CDF of gravitation mass (solar masses). Supply a function of one variable (the probability mass ``x``), in the form of an expression that can be evaluated with the ``eval(...)`` builtin function, i.e., scipy.stats.truncnorm(x, ...). Note that the prior default PDF is uniform (with compact support), so do not supply a CDF if a uniform prior is desired, or to be explicit, use: DEFAULT UNIFORM. You must use DEFAULT UNIFORM to overwrite a default CDF shown in the auto-generated configuration file, unless the parameter is fixed/derived in which case the prior flag is silently ignored. You can also use the flag more than once: the last usage must be an expression that will be dynamically evaluated using the ``eval(...)`` builtin and must return a float to set as the parameter value; the other usages can be helper statements executed with the ``exec(...)`` builtin, e.g., to set temporary local variables to make the code (and configuration file more readable).',
                    comment=True)

parser.add_argument('--mass-value',
                    type=str,
                    action=CompileAction,
                    help='Value of gravitational mass (solar masses). No value means the parameter is free, whilst a value means the parameter is fixed (and the prior may need to be modified manually). If you want the parameter to be derived from other parameters in a complex way, manual modification of the main module is necessary, but we support functions of one parameter in the configuration file. If the parameter is derived from one other parameter, e.g., the temperature of a hot region is derived from the temperature of the other hot region, then the value needs to be written using the following template: lambda x: f(x), "parameter", "space". In this template: f(x) is a function of the parameter x from which the value is derived; "parameter" is the name of the parameter x as a string; and "space" is an object in the global namespace, with name written as a string, which is a (sub)space of parameters from which the current value of parameter x can be accessed via getitem magic using "parameter".',
                    comment=True,
                    empty_lines_below=2)

parser.add_argument('--radius-bounds',
                    type=str,
                    nargs=2,
                    action=CompileAction,
                    help='Bounds of coordinate equatorial radius (km). If no bounds are given (``None``), and no value is given (``None``), bounds default to the source code strict bounds.',
                    comment=True,
                    comment_line_above='rule')

parser.add_argument('--radius-value',
                    type=str,
                    action=CompileAction,
                    help='Value of coordinate equatorial radius (km). No value means the parameter is free, whilst a value means the parameter is fixed (and the prior may need to be modified manually). If you want the parameter to be derived from other parameters in a complex way, manual modification of the main module is necessary, but we support functions of one parameter in the configuration file. If the parameter is derived from one other parameter, e.g., the temperature of a hot region is derived from the temperature of the other hot region, then the value needs to be written using the following template: lambda x: f(x), "parameter", "space". In this template: f(x) is a function of the parameter x from which the value is derived; "parameter" is the name of the parameter x as a string; and "space" is an object in the global namespace, with name written as a string, which is a (sub)space of parameters from which the current value of parameter x can be accessed via getitem magic using "parameter".',
                    comment=True,
                    empty_lines_below=2)

parser.add_argument('--cos-inclination-bounds',
                    type=str,
                    nargs=2,
                    action=CompileAction,
                    help='Bounds of cosine of Earth colatitude (inclination) w.r.t to stellar rotation axis. If no bounds are given (``None``), and no value is given (``None``), bounds default to the source code strict bounds.',
                    comment=True,
                    comment_line_above='rule')

parser.add_argument('--cos-inclination-prior',
                    type=str,
                    nargs='*',
                    action=NullAction,
                    help='Prior inverse CDF of cosine of Earth inclination to stellar spin axis. Supply a function of one variable (the probability mass ``x``), in the form of an expression that can be evaluated with the ``eval(...)`` builtin function, i.e., scipy.stats.truncnorm(x, ...). Note that the prior default PDF is uniform (with compact support), so do not supply a CDF if a uniform prior is desired, or to be explicit, use: DEFAULT UNIFORM. You must use DEFAULT UNIFORM to overwrite a default CDF shown in the auto-generated configuration file, unless the parameter is fixed/derived in which case the prior flag is silently ignored. You can also use the flag more than once: the last usage must be an expression that will be dynamically evaluated using the ``eval(...)`` builtin and must return a float to set as the parameter value; the other usages can be helper statements executed with the ``exec(...)`` builtin, e.g., to set temporary local variables to make the code (and configuration file more readable).',
                    comment=True)

parser.add_argument('--cos-inclination-value',
                    type=str,
                    action=CompileAction,
                    help='Value of cosine of Earth colatitude (inclination) w.r.t to stellar rotation axis. No value means the parameter is free, whilst a value means the parameter is fixed (and the prior may need to be modified manually). If you want the parameter to be derived from other parameters in a complex way, manual modification of the main module is necessary, but we support functions of one parameter in the configuration file.',
                    comment=True,
                    empty_lines_below=2)

parser.add_argument('--distance-bounds',
                    type=str,
                    nargs=2,
                    action=CompileAction,
                    help='Bounds of distance to source (kpc). If no bounds are given (``None``), and no value is given (``None``), bounds default to the source code strict bounds.',
                    comment=True,
                    comment_line_above='rule')

parser.add_argument('--distance-prior',
                    type=str,
                    nargs='*',
                    action=NullAction,
                    help='Prior inverse CDF of distance to source (kpc). Supply a function of one variable (the probability mass ``x``), in the form of an expression that can be evaluated with the ``eval(...)`` builtin function, i.e., scipy.stats.truncnorm(x, ...). Note that the prior default PDF is uniform (with compact support), so do not supply a CDF if a uniform prior is desired, or to be explicit, use: DEFAULT UNIFORM. You must use DEFAULT UNIFORM to overwrite a default CDF shown in the auto-generated configuration file, unless the parameter is fixed/derived in which case the prior flag is silently ignored. You can also use the flag more than once: the last usage must be an expression that will be dynamically evaluated using the ``eval(...)`` builtin and must return a float to set as the parameter value; the other usages can be helper statements executed with the ``exec(...)`` builtin, e.g., to set temporary local variables to make the code (and configuration file more readable).',
                    comment=True)

parser.add_argument('--distance-value',
                    type=str,
                    action=CompileAction,
                    help='Value of distance to source (kpc). No value means the parameter is free, whilst a value means the parameter is fixed (and the prior may need to be modified manually). If you want the parameter to be derived from other parameters in a complex way, manual modification of the main module is necessary, but we support functions of one parameter in the configuration file. If the parameter is derived from one other parameter, e.g., the temperature of a hot region is derived from the temperature of the other hot region, then the value needs to be written using the following template: lambda x: f(x), "parameter", "space". In this template: f(x) is a function of the parameter x from which the value is derived; "parameter" is the name of the parameter x as a string; and "space" is an object in the global namespace, with name written as a string, which is a (sub)space of parameters from which the current value of parameter x can be accessed via getitem magic using "parameter".',
                    comment=True,
                    empty_lines_below=2)


parser.add_argument('--p-super-colatitude-bounds',
                    type=str,
                    nargs=2,
                    action=CompileAction,
                    help='Bounds of hot region "p" super-member colatitude w.r.t stellar spin axis (radians). If no bounds are given (``None``), and no value is given (``None``), bounds default to the source code strict bounds.',
                    comment_line_above='"p" hot region parameter flags',
                    comment=True)

parser.add_argument('--p-super-colatitude-value',
                    type=str,
                    action=CompileAction,
                    help='Value of hot region "p" super-member colatitude w.r.t stellar spin axis (radians). No value means the parameter is free, whilst a value means the parameter is fixed (and the prior may need to be modified manually). If you want the parameter to be derived from other parameters in a complex way, manual modification of the main module is necessary, but we support functions of one parameter in the configuration file. If the parameter is derived from one other parameter, e.g., the temperature of a hot region is derived from the temperature of the other hot region, then the value needs to be written using the following template: lambda x: f(x), "parameter", "space". In this template: f(x) is a function of the parameter x from which the value is derived; "parameter" is the name of the parameter x as a string; and "space" is an object in the global namespace, with name written as a string, which is a (sub)space of parameters from which the current value of parameter x can be accessed via getitem magic using "parameter".',
                    comment=True,
                    empty_lines_below=2)

parser.add_argument('--p-super-radius-bounds',
                    type=str,
                    nargs=2,
                    action=CompileAction,
                    help='Bounds of hot region "p" super-member angular radius (radians). If no bounds are given (``None``), and no value is given (``None``), bounds default to the source code strict bounds.',
                    comment=True,
                    comment_line_above='rule')

parser.add_argument('--p-super-radius-value',
                    type=str,
                    action=CompileAction,
                    help='Value of hot region "p" super-member angular radius (radians). No value means the parameter is free, whilst a value means the parameter is fixed (and the prior may need to be modified manually). If you want the parameter to be derived from other parameters in a complex way, manual modification of the main module is necessary, but we support functions of one parameter in the configuration file. If the parameter is derived from one other parameter, e.g., the temperature of a hot region is derived from the temperature of the other hot region, then the value needs to be written using the following template: lambda x: f(x), "parameter", "space". In this template: f(x) is a function of the parameter x from which the value is derived; "parameter" is the name of the parameter x as a string; and "space" is an object in the global namespace, with name written as a string, which is a (sub)space of parameters from which the current value of parameter x can be accessed via getitem magic using "parameter".',
                    comment=True,
                    empty_lines_below=2)

parser.add_argument('--p-super-temperature-bounds',
                    type=str,
                    nargs=2,
                    action=CompileAction,
                    help='Bounds of hot region "p" super-member log10(temperature [K]). If no bounds are given (``None``), and no value is given (``None``), bounds default to the source code strict bounds.',
                    comment=True,
                    comment_line_above='rule')

parser.add_argument('--p-super-temperature-prior',
                    type=str,
                    nargs='*',
                    action=NullAction,
                    help='Prior inverse CDF of hot-region p superseding region log10(temperature [K]). Supply a function of one variable (the probability mass ``x``), in the form of an expression that can be evaluated with the ``eval(...)`` builtin function, i.e., scipy.stats.truncnorm(x, ...). Note that the prior default PDF is uniform (with compact support), so do not supply a CDF if a uniform prior is desired, or to be explicit, use: DEFAULT UNIFORM. You must use DEFAULT UNIFORM to overwrite a default CDF shown in the auto-generated configuration file, unless the parameter is fixed/derived in which case the prior flag is silently ignored. You can also use the flag more than once: the last usage must be an expression that will be dynamically evaluated using the ``eval(...)`` builtin and must return a float to set as the parameter value; the other usages can be helper statements executed with the ``exec(...)`` builtin, e.g., to set temporary local variables to make the code (and configuration file more readable).',
                    comment=True)

parser.add_argument('--p-super-temperature-value',
                    type=str,
                    action=CompileAction,
                    help='Value of hot region "p" super-member log10(temperature [K]). No value means the parameter is free, whilst a value means the parameter is fixed (and the prior may need to be modified manually). If you want the parameter to be derived from other parameters in a complex way, manual modification of the main module is necessary, but we support functions of one parameter in the configuration file. If the parameter is derived from one other parameter, e.g., the temperature of a hot region is derived from the temperature of the other hot region, then the value needs to be written using the following template: lambda x: f(x), "parameter", "space". In this template: f(x) is a function of the parameter x from which the value is derived; "parameter" is the name of the parameter x as a string; and "space" is an object in the global namespace, with name written as a string, which is a (sub)space of parameters from which the current value of parameter x can be accessed via getitem magic using "parameter".',
                    comment=True,
                    empty_lines_below=2)

parser.add_argument('--p-omit-radius-bounds',
                    type=str,
                    nargs=2,
                    action=CompileAction,
                    help='Bounds of hot region "p" omit-member angular radius (radians). If no bounds are given (``None``), and no value is given (``None``), bounds default to the source code strict bounds.',
                    comment=True,
                    comment_line_above='rule')

parser.add_argument('--p-omit-radius-value',
                    type=str,
                    action=CompileAction,
                    help='Value of hot region "p" omit-member angular radius (radians). No value means the parameter is free, whilst a value means the parameter is fixed (and the prior may need to be modified manually). If you want the parameter to be derived from other parameters in a complex way, manual modification of the main module is necessary, but we support functions of one parameter in the configuration file. If the parameter is derived from one other parameter, e.g., the temperature of a hot region is derived from the temperature of the other hot region, then the value needs to be written using the following template: lambda x: f(x), "parameter", "space". In this template: f(x) is a function of the parameter x from which the value is derived; "parameter" is the name of the parameter x as a string; and "space" is an object in the global namespace, with name written as a string, which is a (sub)space of parameters from which the current value of parameter x can be accessed via getitem magic using "parameter".',
                    comment=True,
                    empty_lines_below=2)

parser.add_argument('--p-phase-shift-bounds',
                    type=str,
                    nargs=2,
                    action=CompileAction,
                    help='Bounds of hot region "p" phase shift (cycles). If no bounds are given (``None``), and no value is given (``None``), bounds default to the source code strict bounds.',
                    comment=True,
                    comment_line_above='rule')

parser.add_argument('--p-phase-shift-value',
                type=str,
                action=CompileAction,
                help='Value of hot region "p" phase shift (cycles). No value means the parameter is free, whilst a value means the parameter is fixed (and the prior may need to be modified manually). If you want the parameter to be derived from other parameters in a complex way, manual modification of the main module is necessary, but we support functions of one parameter in the configuration file. If the parameter is derived from one other parameter, e.g., the temperature of a hot region is derived from the temperature of the other hot region, then the value needs to be written using the following template: lambda x: f(x), "parameter", "space". In this template: f(x) is a function of the parameter x from which the value is derived; "parameter" is the name of the parameter x as a string; and "space" is an object in the global namespace, with name written as a string, which is a (sub)space of parameters from which the current value of parameter x can be accessed via getitem magic using "parameter".',
                comment=True,
                empty_lines_below=2)

parser.add_argument('--s-super-colatitude-bounds',
                    type=str,
                    nargs=2,
                    action=CompileAction,
                    help='Bounds of hot region "s" super-member colatitude w.r.t stellar spin axis (radians). If no bounds are given (``None``), and no value is given (``None``), bounds default to the source code strict bounds.',
                    comment_line_above='"s" hot region parameter flags',
                    comment=True)

parser.add_argument('--s-super-colatitude-value',
                    type=str,
                    action=CompileAction,
                    help='Value of hot region "s" super-member colatitude w.r.t stellar spin axis (radians). No value means the parameter is free, whilst a value means the parameter is fixed (and the prior may need to be modified manually). If you want the parameter to be derived from other parameters in a complex way, manual modification of the main module is necessary, but we support functions of one parameter in the configuration file. If the parameter is derived from one other parameter, e.g., the temperature of a hot region is derived from the temperature of the other hot region, then the value needs to be written using the following template: lambda x: f(x), "parameter", "space". In this template: f(x) is a function of the parameter x from which the value is derived; "parameter" is the name of the parameter x as a string; and "space" is an object in the global namespace, with name written as a string, which is a (sub)space of parameters from which the current value of parameter x can be accessed via getitem magic using "parameter".',
                    comment=True,
                    empty_lines_below=2)

parser.add_argument('--s-super-radius-bounds',
                    type=str,
                    nargs=2,
                    action=CompileAction,
                    help='Bounds of hot region "s" super-member angular radius (radians). If no bounds are given (``None``), and no value is given (``None``), bounds default to the source code strict bounds.',
                    comment=True,
                    comment_line_above='rule')

parser.add_argument('--s-super-radius-value',
                    type=str,
                    action=CompileAction,
                    help='Value of hot region "s" super-member angular radius (radians). No value means the parameter is free, whilst a value means the parameter is fixed (and the prior may need to be modified manually). If you want the parameter to be derived from other parameters in a complex way, manual modification of the main module is necessary, but we support functions of one parameter in the configuration file. If the parameter is derived from one other parameter, e.g., the temperature of a hot region is derived from the temperature of the other hot region, then the value needs to be written using the following template: lambda x: f(x), "parameter", "space". In this template: f(x) is a function of the parameter x from which the value is derived; "parameter" is the name of the parameter x as a string; and "space" is an object in the global namespace, with name written as a string, which is a (sub)space of parameters from which the current value of parameter x can be accessed via getitem magic using "parameter".',
                    comment=True,
                    empty_lines_below=2)

parser.add_argument('--s-super-temperature-bounds',
                    type=str,
                    nargs=2,
                    action=CompileAction,
                    help='Bounds of hot region "s" super-member log10(temperature [K]). If no bounds are given (``None``), and no value is given (``None``), bounds default to the source code strict bounds.',
                    comment=True,
                    comment_line_above='rule')

parser.add_argument('--s-super-temperature-prior',
                    type=str,
                    nargs='*',
                    action=NullAction,
                    help='Prior inverse CDF of hot-region s superseding region log10(temperature [K]). Supply a function of one variable (the probability mass ``x``), in the form of an expression that can be evaluated with the ``eval(...)`` builtin function, i.e., scipy.stats.truncnorm(x, ...). Note that the prior default PDF is uniform (with compact support), so do not supply a CDF if a uniform prior is desired, or to be explicit, use: DEFAULT UNIFORM. You must use DEFAULT UNIFORM to overwrite a default CDF shown in the auto-generated configuration file, unless the parameter is fixed/derived in which case the prior flag is silently ignored. You can also use the flag more than once: the last usage must be an expression that will be dynamically evaluated using the ``eval(...)`` builtin and must return a float to set as the parameter value; the other usages can be helper statements executed with the ``exec(...)`` builtin, e.g., to set temporary local variables to make the code (and configuration file more readable).',
                    comment=True)

parser.add_argument('--s-super-temperature-value',
                    type=str,
                    action=CompileAction,
                    help='Value of hot region "s" super-member log10(temperature [K]). No value means the parameter is free, whilst a value means the parameter is fixed (and the prior may need to be modified manually). If you want the parameter to be derived from other parameters in a complex way, manual modification of the main module is necessary, but we support functions of one parameter in the configuration file. If the parameter is derived from one other parameter, e.g., the temperature of a hot region is derived from the temperature of the other hot region, then the value needs to be written using the following template: lambda x: f(x), "parameter", "space". In this template: f(x) is a function of the parameter x from which the value is derived; "parameter" is the name of the parameter x as a string; and "space" is an object in the global namespace, with name written as a string, which is a (sub)space of parameters from which the current value of parameter x can be accessed via getitem magic using "parameter".',
                    comment=True,
                    empty_lines_below=2)

parser.add_argument('--s-cede-radius-bounds',
                    type=str,
                    nargs=2,
                    action=CompileAction,
                    help='Bounds of hot region "s" cede-member angular radius (radians). If no bounds are given (``None``), and no value is given (``None``), bounds default to the source code strict bounds.',
                    comment=True,
                    comment_line_above='rule')

parser.add_argument('--s-cede-radius-value',
                    type=str,
                    action=CompileAction,
                    help='Value of hot region "s" cede-member angular radius (radians). No value means the parameter is free, whilst a value means the parameter is fixed (and the prior may need to be modified manually). If you want the parameter to be derived from other parameters in a complex way, manual modification of the main module is necessary, but we support functions of one parameter in the configuration file. If the parameter is derived from one other parameter, e.g., the temperature of a hot region is derived from the temperature of the other hot region, then the value needs to be written using the following template: lambda x: f(x), "parameter", "space". In this template: f(x) is a function of the parameter x from which the value is derived; "parameter" is the name of the parameter x as a string; and "space" is an object in the global namespace, with name written as a string, which is a (sub)space of parameters from which the current value of parameter x can be accessed via getitem magic using "parameter".',
                    comment=True,
                    empty_lines_below=2)

parser.add_argument('--s-cede-temperature-bounds',
                    type=str,
                    nargs=2,
                    action=CompileAction,
                    help='Bounds of hot region "s" cede-member log10(temperature [K]). If no bounds are given (``None``), and no value is given (``None``), bounds default to the source code strict bounds.',
                    comment=True,
                    comment_line_above='rule')

parser.add_argument('--s-cede-temperature-prior',
                    type=str,
                    nargs='*',
                    action=NullAction,
                    help='Prior inverse CDF of hot-region s ceding region log10(temperature [K]). Supply a function of one variable (the probability mass ``x``), in the form of an expression that can be evaluated with the ``eval(...)`` builtin function, i.e., scipy.stats.truncnorm(x, ...). Note that the prior default PDF is uniform (with compact support), so do not supply a CDF if a uniform prior is desired, or to be explicit, use: DEFAULT UNIFORM. You must use DEFAULT UNIFORM to overwrite a default CDF shown in the auto-generated configuration file, unless the parameter is fixed/derived in which case the prior flag is silently ignored. You can also use the flag more than once: the last usage must be an expression that will be dynamically evaluated using the ``eval(...)`` builtin and must return a float to set as the parameter value; the other usages can be helper statements executed with the ``exec(...)`` builtin, e.g., to set temporary local variables to make the code (and configuration file more readable).',
                    comment=True)

parser.add_argument('--s-cede-temperature-value',
                    type=str,
                    action=CompileAction,
                    help='Value of hot region "s" cede-member log10(temperature [K]). No value means the parameter is free, whilst a value means the parameter is fixed (and the prior may need to be modified manually). If you want the parameter to be derived from other parameters in a complex way, manual modification of the main module is necessary, but we support functions of one parameter in the configuration file. If the parameter is derived from one other parameter, e.g., the temperature of a hot region is derived from the temperature of the other hot region, then the value needs to be written using the following template: lambda x: f(x), "parameter", "space". In this template: f(x) is a function of the parameter x from which the value is derived; "parameter" is the name of the parameter x as a string; and "space" is an object in the global namespace, with name written as a string, which is a (sub)space of parameters from which the current value of parameter x can be accessed via getitem magic using "parameter".',
                    comment=True,
                    empty_lines_below=2)

parser.add_argument('--s-cede-colatitude-bounds',
                    type=str,
                    nargs=2,
                    action=CompileAction,
                    help='Bounds of hot region "s" cede-member colatitude w.r.t stellar spin axis (radians). If no bounds are given (``None``), and no value is given (``None``), bounds default to the source code strict bounds.',
                    comment=True,
                    comment_line_above='rule')

parser.add_argument('--s-cede-colatitude-value',
                    type=str,
                    action=CompileAction,
                    help='Value of hot region "s" cede-member colatitude w.r.t stellar spin axis (radians). No value means the parameter is free, whilst a value means the parameter is fixed (and the prior may need to be modified manually). If you want the parameter to be derived from other parameters in a complex way, manual modification of the main module is necessary, but we support functions of one parameter in the configuration file. If the parameter is derived from one other parameter, e.g., the temperature of a hot region is derived from the temperature of the other hot region, then the value needs to be written using the following template: lambda x: f(x), "parameter", "space". In this template: f(x) is a function of the parameter x from which the value is derived; "parameter" is the name of the parameter x as a string; and "space" is an object in the global namespace, with name written as a string, which is a (sub)space of parameters from which the current value of parameter x can be accessed via getitem magic using "parameter".',
                    comment=True,
                    empty_lines_below=2)

parser.add_argument('--s-cede-azimuth-bounds',
                    type=str,
                    nargs=2,
                    action=CompileAction,
                    help='Bounds of hot region "s" cede-member azimuth relative to super-member (radians). If no bounds are given (``None``), and no value is given (``None``), bounds default to the source code strict bounds.',
                    comment=True,
                    comment_line_above='rule')

parser.add_argument('--s-cede-azimuth-value',
                    type=str,
                    action=CompileAction,
                    help='Value of hot region "s" cede-member azimuth relative to super-member (radians). No value means the parameter is free, whilst a value means the parameter is fixed (and the prior may need to be modified manually). If you want the parameter to be derived from other parameters in a complex way, manual modification of the main module is necessary, but we support functions of one parameter in the configuration file. If the parameter is derived from one other parameter, e.g., the temperature of a hot region is derived from the temperature of the other hot region, then the value needs to be written using the following template: lambda x: f(x), "parameter", "space". In this template: f(x) is a function of the parameter x from which the value is derived; "parameter" is the name of the parameter x as a string; and "space" is an object in the global namespace, with name written as a string, which is a (sub)space of parameters from which the current value of parameter x can be accessed via getitem magic using "parameter".',
                    comment=True,
                    empty_lines_below=2)

parser.add_argument('--s-phase-shift-bounds',
                    type=str,
                    nargs=2,
                    action=CompileAction,
                    help='Bounds of hot region "s" phase shift (cycles). If no bounds are given (``None``), and no value is given (``None``), bounds default to the source code strict bounds.',
                    comment=True,
                    comment_line_above='rule')

parser.add_argument('--s-phase-shift-value',
                type=str,
                action=CompileAction,
                help='Value of hot region "s" phase shift (cycles). No value means the parameter is free, whilst a value means the parameter is fixed (and the prior may need to be modified manually). If you want the parameter to be derived from other parameters in a complex way, manual modification of the main module is necessary, but we support functions of one parameter in the configuration file. If the parameter is derived from one other parameter, e.g., the temperature of a hot region is derived from the temperature of the other hot region, then the value needs to be written using the following template: lambda x: f(x), "parameter", "space". In this template: f(x) is a function of the parameter x from which the value is derived; "parameter" is the name of the parameter x as a string; and "space" is an object in the global namespace, with name written as a string, which is a (sub)space of parameters from which the current value of parameter x can be accessed via getitem magic using "parameter".',
                comment=True,
                empty_lines_below=2)

parser.add_argument('--p-sqrt-num-cells',
                    type=int,
                    default=32,
                    help='Target square-root of the number of cells spanning (raditing subset of) hot region "p".',
                    comment_line_above='"p" hot region resolution flags')

parser.add_argument('--p-min-sqrt-num-cells',
                    type=int,
                    default=10,
                    help='Minimum square-root of the number of cells constituting hot region "p" mesh.')

parser.add_argument('--p-max-sqrt-num-cells',
                    type=int,
                    default=80,
                    help='Maximum square-root of the number of cells constituting hot region "p" mesh.')

parser.add_argument('--p-num-leaves',
                    type=int,
                    default=64,
                    help='Number of phases on unit interval at which to compute hot region "p" signal.')

parser.add_argument('--p-num-rays',
                    type=int,
                    default=512,
                    help='Number of rays per iso-latitude mesh subset to trace when computing hot region "p" signal.',
                    empty_lines_below=2)

parser.add_argument('--s-sqrt-num-cells',
                    type=int,
                    default=32,
                    help='Target square-root of the number of cells spanning (raditing subset of) hot region "s".',
                    comment_line_above='"s" hot region resolution flags')

parser.add_argument('--s-min-sqrt-num-cells',
                    type=int,
                    default=10,
                    help='Minimum square-root of the number of cells constituting hot region "s" mesh.')

parser.add_argument('--s-max-sqrt-num-cells',
                    type=int,
                    default=80,
                    help='Maximum square-root of the number of cells constituting hot region "s" mesh.')

parser.add_argument('--s-num-leaves',
                    type=int,
                    default=64,
                    help='Number of phases on unit interval at which to compute hot region "s" signal.')

parser.add_argument('--s-num-rays',
                    type=int,
                    default=512,
                    help='Number of rays per iso-latitude mesh subset to trace when computing hot region "s" signal.',
                    empty_lines_below=2)

parser.add_argument('--hot-atmosphere-path', type=str, help='Absolute or relative path to hot atmosphere file.',
                             comment_line_above='hot atmosphere flags')

parser.add_argument('--hot-atmosphere-size',
                    type=int,
                    nargs=4,
                    action=NullAction,
                    help='Size of each of the four dimensions of the numeric atmosphere table for the hot regions.',
                    empty_lines_below=2)

parser.add_argument('--image-order-limit',
                         type=int,
                         default=3,
                         help='The highest-order image to sum over. Either a positive integer, or do not pass an argument if a hard limit is not desired.',
                         comment_line_above='global resolution flags')

parser.add_argument('--number-energies',
                         type=int,
                         default=128,
                         help='Number of energies, distributed over instrument wavebands, to compute incident photon specific flux at.')

parser.add_argument('--maximum-energy-ray-tracing',
                         type=int,
                         help='Maximum energy for ray tracing. Useful if there is a background component such as a powerlaw that is jointly modelled with higher-energy event data using a subset of instruments.',
                         comment=True,
                         empty_lines_below=2)

parser.add_argument('--openmp-threads',
                    type=int,
                    default=1,
                    help='Number of OpenMP threads to spawn during likelihood function calls.',
                    comment_line_above='miscellaneous flags')

parser.add_argument('--parameters-externally-updated',
                    type=str_to_bool,
                    default=True,
                    help='Are the parameters updated before calling the likelihood object, e.g., in the prior object?',
                    empty_lines_below=2)

parser.add_argument('--multinest', action='store_true', help='Launch MultiNest sampler if module is executed.',
                         comment_line_above='runtime flags')

parser.add_argument('--resume', action='store_true', help='Resume sampling if module is executed.')


parser.add_argument('--sample-files-directory-path',
                    type=str,
                    default='samples/',
                    help='Absolute or relative path to sample file directory. If no path is provided, the default (relative) path is "samples/".')

parser.add_argument('--sample-files-root',
                    type=str,
                    help='The root name of the sample files (i.e., without a file extension) to be generated or already generated by MultiNest. If no path is provided, the sample file root name will be constructed automatically from other sampling process settings.',
                    comment=True)

parser.add_argument('--number-iterations-per-write',
                    type=int,
                    default=100,
                    help='Number of nested replacements per write to disk of the sampling process to enable resumption. Posterior files are generated at a cadence of 10x this number.')

parser.add_argument('--number-live-points',
                    type=int,
                    default=1000,
                    help='Number of live points in nested sampling process.')

parser.add_argument('--hypervolume-expansion-factor',
                    type=float,
                    default=10.0,
                    help='Factor by which to expand the hyperellisoid union that approximately minimally bounds the set of live points.')

parser.add_argument('--constant-efficiency-variant',
                    action='store_true',
                    help='Activate MultiNest constant efficiency sampling variant? Warning: only use this option if computational resources are limited.',
                    comment=True)

parser.add_argument('--mode-separation-variant',
                    action='store_true',
                    help='Activate mode-separation sampling variant? Live point threads (an initial live point and the chain of subsequent replacements) do not migrate between threads by default.',
                    comment=True)

parser.add_argument('--estimated-remaining-log-evidence',
                    type=float,
                    default=0.1,
                    help='Estimated remaining log-evidence for sampling process termination.')

parser.add_argument('--maximum-number-nested-replacement-iterations',
                    type=int,
                    default=-1,
                    help='Maximum number of nested replacements for termination of the nested sampling process. Use negative one (the default) to terminate based on estimated remaining log-evidence instead.')


import xpsi

if __name__ == '__main__':
    if xpsi._verbose:
        print('Parsing configuration file...')
    args, _ = parser.parse_known_args()
    if xpsi._verbose:
        print('Configuration file parsed.')
else:
    if xpsi._verbose:
        print('Parsing configuration file...')
    args, _ = parser.parse_known_args(['@./config.ini'])
    if xpsi._verbose:
        print('Configuration file parsed.')


import os

import numpy as np
import math

from xpsi.Parameter import Derive
from xpsi import HotRegions

print('Rank reporting: %d' % xpsi._rank)
if __name__ == '__main__':
    from CustomInstrument import CustomInstrument
    from CustomSignal import CustomSignal
    from CustomInterstellar import CustomInterstellar

    try:
        from CustomPhotosphere import CustomPhotosphere
    except ImportError:
        from xpsi import Photosphere as CustomPhotosphere

    from CustomPrior import CustomPrior

else:
    from .CustomInstrument import CustomInstrument
    from .CustomSignal import CustomSignal
    from .CustomInterstellar import CustomInterstellar

    try:
        from .CustomPhotosphere import CustomPhotosphere
    except ImportError:
        from xpsi import Photosphere as CustomPhotosphere

    from .CustomPrior import CustomPrior

if args.main_import_statements is not None:
    for import_statement in args.main_import_statements:
        exec(import_statement)

if args.main_global_statements is not None:
    for global_statement in args.main_global_statements:
        exec(global_statement)

class namespace():
    pass

def parse_bounds(bounds, value, default_to_free=True):
    if bounds is not None:
        bounds[0] = eval(bounds[0])
        bounds[1] = eval(bounds[1])
        return tuple(bounds)
    elif default_to_free:
        return None if value is not None else (None, None)

    return None

def derived_parameter(func, parameter, space='caller'):

    class derive(Derive):

        def __init__(self):
            self.space = compile(space, '<string>', 'eval')

        def __call__(self, boundto, caller=None):
            return func(eval(self.space)[parameter])

    return derive()

def parse_value(value):
    if value is not None:
        try:
            return float(eval(value))
        except ValueError:
            return derived_parameter(*eval(value))
    else:
        return None

bounds = dict(neutral_hydrogen_column_density = parse_bounds(args.neutral_hydrogen_column_density_bounds,
                                                              args.neutral_hydrogen_column_density_value))
values = dict(neutral_hydrogen_column_density = parse_value(args.neutral_hydrogen_column_density_value))

interstellar = CustomInterstellar.load(args.attenuation_path,
                                       args.attenuation_energy_column,
                                       args.attenuation_column,
                                       bounds = bounds,
                                       values = values)

signals = [[],]

XTI = namespace()

if args.XTI_energy_independent_effective_area_scaling_factor_value is not None:
    if eval(args.XTI_energy_independent_effective_area_scaling_factor_value) in ['XTI']:
        values = dict(energy_independent_effective_area_scaling_factor = derived_parameter(lambda x: x,
                                              'energy_independent_effective_area_scaling_factor',
                                              eval(args.XTI_energy_independent_effective_area_scaling_factor_value) + '.instrument'))
    else:
        values = dict(energy_independent_effective_area_scaling_factor = parse_value(args.XTI_energy_independent_effective_area_scaling_factor_value))
else:
    values = {}

bounds = dict(energy_independent_effective_area_scaling_factor = parse_bounds(args.XTI_energy_independent_effective_area_scaling_factor_bounds,
                                   args.XTI_energy_independent_effective_area_scaling_factor_value,
                                   default_to_free = False))

XTI.instrument = CustomInstrument.XTI(bounds=bounds,
                                  values=values,
                                  ARF=args.XTI_arf_path,
                                  RMF=args.XTI_rmf_path,
                                  channel_energies=args.XTI_channels_path,
                                  max_input=args.XTI_input_bounds[1],
                                  max_channel=args.XTI_channel_bounds[1],
                                  min_input=args.XTI_input_bounds[0],
                                  min_channel=args.XTI_channel_bounds[0],
                                  effective_area_scaling_factor=eval(args.XTI_effective_area_scaling_factor),
                                  ARF_skiprows=args.XTI_arf_skiprows,
                                  ARF_low_column=args.XTI_arf_low_column,
                                  ARF_high_column=args.XTI_arf_high_column,
                                  ARF_area_column=args.XTI_arf_area_column,
                                  RMF_skiprows=args.XTI_rmf_skiprows,
                                  RMF_usecol=args.XTI_rmf_usecol,
                                  channel_energies_skiprows=args.XTI_channel_energies_skiprows,
                                  channel_energies_low_column=args.XTI_channel_energies_low_column)

try:
    counts = np.loadtxt(args.XTI_count_matrix_path, dtype=np.double)
except ValueError:
    XTI.data = xpsi.Data.bin__event_list(args.XTI_event_path,
                                         channels=XTI.instrument.channels,
                                         phases=np.linspace(0.0, 1.0, args.XTI_number_phase_bins + 1),
                                         channel_column=args.XTI_event_file_channel_column,
                                         phase_column=args.XTI_event_file_phase_column if args.XTI_number_phase_bins > 1 else None,
                                         phase_averaged=True if args.XTI_number_phase_bins == 1 else False,
                                         channel_edges=XTI.instrument.channel_edges,
                                         skiprows=args.XTI_event_file_skiprows,
                                         eV=True if args.XTI_events_in_eV else False,
                                         dtype=getattr(np, args.XTI_count_matrix_type),
                                         first=0,
                                         last=len(XTI.instrument.channels) - 1,
                                         exposure_time=args.XTI_exposure_time)

    np.savetxt(args.XTI_event_path.replace('.txt','_converted_to_counts.txt'), XTI.data.counts)
    print('Counts file saved as: '+args.XTI_event_path.replace('.txt','_converted_to_counts.txt'))
    print('Update configuration file to take in counts file to save computation time.')
else:
    if counts.ndim == 1:
        counts = counts.reshape(-1,1)

    XTI.data = xpsi.Data(counts,
                           channels=XTI.instrument.channels,
                           phases=np.linspace(0.0, 1.0, args.XTI_number_phase_bins + 1),
                           first=0,
                           last=len(XTI.instrument.channels) - 1,
                           exposure_time=args.XTI_exposure_time)

if args.XTI_background_prior_support_path:
    support = np.loadtxt(args.XTI_background_path,
                         skiprows=args.XTI_background_skiprows,
                         dtype=np.double)

elif args.XTI_background_path:
    spectrum = np.loadtxt(args.XTI_background_path,
                          skiprows=args.XTI_background_skiprows,
                          usecols=args.XTI_background_usecol,
                          dtype=np.double)[XTI.instrument.channels]

    support = np.zeros((len(spectrum), 2), dtype=np.double)
    support[:,0] = spectrum - args.XTI_background_prior_support_half_width * np.sqrt(spectrum)
    support[support[:,0] < 0.0, 0] = 0.0
    support[:,1] = spectrum + args.XTI_background_prior_support_half_width * np.sqrt(spectrum)

    for i in range(support.shape[0]):
        if support[i,1] == 0.0:
            for j in range(1, support.shape[0]):
                if i+j < support.shape[0] and support[i+j,1] > 0.0:
                    support[i,1] = support[i+j,1]
                    break
                elif i-j >= 0 and support[i-j,1] > 0.0:
                    support[i,1] = support[i-j,1]
                    break

    support *= (XTI.data.exposure_time / args.XTI_background_exposure_time) * float(eval(args.XTI_background_scaling_factor)) # exposure ratio * scaling

    support /= XTI.data.exposure_time # need count rate, so divide by exposure time
else:
    support = None

XTI.background = None

if args.XTI_phase_shift_value is not None:
    if eval(args.XTI_phase_shift_value) in ['XTI']:
        values = dict(phase_shift = derived_parameter(lambda x: x,
                                                'phase_shift',
                                                eval(args.XTI_phase_shift_value) + '.signal'))
    else:
        values = dict(phase_shift = parse_value(args.XTI_phase_shift_value))
else:
    values = {}

bounds = dict(phase_shift = parse_bounds(args.XTI_phase_shift_bounds,
                                         args.XTI_phase_shift_value,
                                         default_to_free=False))

XTI.signal = CustomSignal(data = XTI.data,
                          instrument = XTI.instrument,
                          interstellar = interstellar,
                          background = XTI.background,
                          cache = False if __name__ == '__main__' else True,
                          bounds = bounds,
                          values = values,
                          workspace_intervals = 1000,
                          epsrel = 1.0e-8,
                          epsilon = 1.0e-3,
                          sigmas = 10.0,
                          support = support,
                          prefix = 'XTI')

signals[0].append(XTI.signal)

bounds = dict(mass = parse_bounds(args.mass_bounds,
                                       args.mass_value),
              radius = parse_bounds(args.radius_bounds,
                                    args.radius_value),
              distance = parse_bounds(args.distance_bounds,
                                      args.distance_value),
              cos_inclination = parse_bounds(args.cos_inclination_bounds,
                                             args.cos_inclination_value))

values = dict(mass = parse_value(args.mass_value),
              radius = parse_value(args.radius_value),
              distance = parse_value(args.distance_value),
              cos_inclination = parse_value(args.cos_inclination_value),
              frequency = 205.0)

spacetime = xpsi.Spacetime(bounds, values)

symmetry = True
omit = True
cede = False
concentric = True

bounds = dict(super_colatitude = parse_bounds(args.p_super_colatitude_bounds,
                                              args.p_super_colatitude_value),
              super_radius = parse_bounds(args.p_super_radius_bounds,
                                          args.p_super_radius_value),
              phase_shift = parse_bounds(args.p_phase_shift_bounds,
                                         args.p_phase_shift_value),
              super_temperature = parse_bounds(args.p_super_temperature_bounds,
                                               args.p_super_temperature_value))

values = dict(super_colatitude = parse_value(args.p_super_colatitude_value),
              super_radius = parse_value(args.p_super_radius_value),
              phase_shift = parse_value(args.p_phase_shift_value),
              super_temperature = parse_value(args.p_super_temperature_value))

bounds['omit_radius'] = parse_bounds(args.p_omit_radius_bounds,
                                     args.p_omit_radius_value)

values['omit_radius'] = parse_value(args.p_omit_radius_value)

primary = xpsi.HotRegion(bounds=bounds,
                            values=values,
                            symmetry=symmetry,
                            omit=omit,
                            cede=cede,
                            concentric=concentric,
                            sqrt_num_cells=args.p_sqrt_num_cells,
                            min_sqrt_num_cells=args.p_min_sqrt_num_cells,
                            max_sqrt_num_cells=args.p_max_sqrt_num_cells,
                            num_leaves=args.p_num_leaves,
                            num_rays=args.p_num_rays,
                            is_antiphased=False,
                            image_order_limit=args.image_order_limit,
                            atm_ext="Num4D" if 'NSX' or 'nsx' in args.hot_atmosphere_model else "BB",
                            prefix='p')

symmetry = True
omit = False
cede = True
concentric = False

bounds = dict(super_colatitude = parse_bounds(args.s_super_colatitude_bounds,
                                              args.s_super_colatitude_value),
              super_radius = parse_bounds(args.s_super_radius_bounds,
                                          args.s_super_radius_value),
              phase_shift = parse_bounds(args.s_phase_shift_bounds,
                                         args.s_phase_shift_value),
              super_temperature = parse_bounds(args.s_super_temperature_bounds,
                                               args.s_super_temperature_value))

values = dict(super_colatitude = parse_value(args.s_super_colatitude_value),
              super_radius = parse_value(args.s_super_radius_value),
              phase_shift = parse_value(args.s_phase_shift_value),
              super_temperature = parse_value(args.s_super_temperature_value))

bounds['cede_radius'] = parse_bounds(args.s_cede_radius_bounds,
                                     args.s_cede_radius_value)

values['cede_radius'] = parse_value(args.s_cede_radius_value)

bounds['cede_temperature'] = parse_bounds(args.s_cede_temperature_bounds,
                                          args.s_cede_temperature_value)

values['cede_temperature'] = parse_value(args.s_cede_temperature_value)

bounds['cede_colatitude'] = parse_bounds(args.s_cede_colatitude_bounds,
                                         args.s_cede_colatitude_value)

values['cede_colatitude'] = parse_value(args.s_cede_colatitude_value)

bounds['cede_azimuth'] = parse_bounds(args.s_cede_azimuth_bounds,
                                      args.s_cede_azimuth_value)

values['cede_azimuth'] = parse_value(args.s_cede_azimuth_value)

secondary = xpsi.HotRegion(bounds=bounds,
                                values=values,
                                symmetry=symmetry,
                                omit=omit,
                                cede=cede,
                                concentric=concentric,
                                sqrt_num_cells=args.s_sqrt_num_cells,
                                min_sqrt_num_cells=args.s_min_sqrt_num_cells,
                                max_sqrt_num_cells=args.s_max_sqrt_num_cells,
                                num_leaves=args.s_num_leaves,
                                num_rays=args.s_num_rays,
                                is_antiphased=False,
                                image_order_limit=args.image_order_limit,
                                atm_ext="Num4D" if 'NSX' or 'nsx' in args.hot_atmosphere_model else "BB",
                                prefix='s')

hot = HotRegions((primary, secondary))

elsewhere = None

photosphere = CustomPhotosphere(hot = hot,
                                     elsewhere = elsewhere,
                                     values = dict(mode_frequency = spacetime['frequency']))

photosphere.hot_atmosphere = args.hot_atmosphere_path

star = xpsi.Star(spacetime = spacetime, photospheres = photosphere)

prior = CustomPrior()

likelihood = xpsi.Likelihood(star = star,
                             signals = signals,
                             num_energies = args.number_energies,
                             threads = args.openmp_threads,
                             externally_updated = args.parameters_externally_updated if args.parameters_externally_updated is not None else True,
                             prior = prior,
                             max_energy = args.maximum_energy_ray_tracing)


if __name__ == '__main__':

    if args.multinest:

        wrapped_params = [0] * len(likelihood)
        for name in likelihood.names:
            if 'phase_shift' or 'azimuth' in name:
                wrapped_params[likelihood.index(name)] = 1

        if args.sample_files_root is None:
            args.sample_files_root = 'nlive{:d}_expf{:.1f}_{}_{}_tol{:.1g}'.format(args.number_live_points,
                                                                                   args.hypervolume_expansion_factor,
                                                                                   'noCONST' if not args.constant_efficiency_variant else 'CONST',
                                                                                   'noMM' if not args.mode_separation_variant else 'MM',
                                                                                   args.estimated_remaining_log_evidence)
        runtime_params = {'resume': args.resume,
                          'importance_nested_sampling': False, # incompatible with xpsi likelihood function
                          'multimodal': args.mode_separation_variant, # this variant, if activated is incompatible with nestcheck
                          'n_clustering_params': None,
                          'outputfiles_basename': os.path.join(args.sample_files_directory_path, args.sample_files_root),
                          'n_iter_before_update': args.number_iterations_per_write,
                          'n_live_points': args.number_live_points,
                          'sampling_efficiency': 1.0 / args.hypervolume_expansion_factor,
                          'const_efficiency_mode': args.constant_efficiency_variant,
                          'wrapped_params': wrapped_params,
                          'evidence_tolerance': args.estimated_remaining_log_evidence,
                          'max_iter': args.maximum_number_nested_replacement_iterations,
                          'verbose': True}

        xpsi.Sample.nested(likelihood, prior, **runtime_params)

else:

    pass
