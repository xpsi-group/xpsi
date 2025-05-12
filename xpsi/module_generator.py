
import os
import sys

import xpsi

def write(filename, module):
    """ Write a module to a file.

    :param filename (str): Name of the file to write to.
    :param module (str): The module to write to the file.
    """
    with open(filename, 'w') as mod:
        _module = ''''''
        for _line in module.splitlines():
            if _module:
                _module += '\n'
            _module += _line.rstrip()
        mod.write(_module)

nindent = '\n    '
indent = '    '

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

        with open('./generate.ini', 'w') as file:
            file.write(config_file)

        print('Configuration file generated.')
        parser.exit()

parser = ArgumentParserCustom(
    description='''
    Script for automated generation of X-PSI model module set.

    Usage: python %(prog)s [-h] @<generate.ini>

    ''',
    fromfile_prefix_chars='@')

parser.add_argument('--generate-config-file', default=argparse.SUPPRESS, action=GenerateConfigAction, help='Generate the meta configuration file template.',
                    destined_for_config_file=False)

parser.add_argument('--telescope',
                    type=str,
                    action='append',
                    help='Telescope name, e.g., NICER. Use argument once per telescope name, and no whitespaces.',
                    comment_line_above='telescope instrument flags')

parser.add_argument('--instrument',
                    type=lambda x: str(x).replace(' ', '_').replace('-', '_'),
                    action='append',
                    help='Name of an instrument on-board a telescope, e.g., XTI. Can use one or more instrument names per telescope name, and no whitespaces.',
                    empty_lines_below=2)

parser.add_argument('--source',
                    type=str,
                    help='The name of the star, e.g., PSR J0740+6620.',
                    comment_line_above='target source flags')

parser.add_argument('--frequency',
                    type=float,
                    required=True,
                    help='The coordinate spin frequency of the star (Hz).',
                    empty_lines_below=2)
                    
parser.add_argument('--use-fits-format',
                    action='store_true',
                    help='Are the source and instrumental files stored in FITS format?',
                    comment_line_above='Input data in FITS format?',
                    empty_lines_below=2,
                    comment=True)

parser.add_argument('--model',
                    type=str,
                    help='A custom model name, e.g., ST-U + NSX-H, otherwise the model name is constructed from other arguments.',
                    comment_line_above='model flags',
                    comment=True)

parser.add_argument('--hot-region-model',
                    type=str,
                    action='append',
                    choices=['ST', 'CST', 'EST', 'PST', 'CDT', 'EDT', 'PDT'],
                    required=True,
                    help='The name of the hot-region model, e.g., ST. Maximum of two argument uses.')

parser.add_argument('--antipodal-reflection-symmetry',
                    action='store_true',
                    help='Are the two hot regions related via antipodal reflection symmetry? E.g., ST-S.',
                    comment=True)

parser.add_argument('--break-hot-region-exchange-degeneracy-with',
                    type=str,
                    default='super_colatitude',
                    help='Hot region parameter name to break hot-region exchange degeneracy with when there are two hot-regions of the same type that are not antipodally reflection-symmetric, e.g., ST+ST (ST-U). An example is e.g., "super_temperature".',
                    comment=True)

def str_to_bool(x):
    if x == 'False':
        return False
    elif x == 'True':
        return True

    raise ValueError('Invalid argument where boolean ``True`` or ``False`` is required.')

parser.add_argument('--is-antiphased',
                    type=str_to_bool,
                    action='append',
                    help='Specify whether the hot regions are anti-phased w.r.t to Earth. If True, the cell mesh shifts by pi radians about the stellar rotation axis for pulse integration and therefore the hot region at phase zero is aligned with the meridian on which the observer’s antipode lies.')

parser.add_argument('--prefix',
                    type=str,
                    action='append',
                    help='Specify the prefixes for hot region parameter naming.')

parser.add_argument('--hot-atmosphere-model',
                    type=str,
                    help='Name of atmosphere model within hot regions, e.g., blackbody or NSX-H.')

parser.add_argument('--hot-atmosphere-load',
                    action='store_true',
                    help='Does a numeric atmosphere table need to be loaded from disk for the hot regions?',
                    comment=True)

parser.add_argument('--elsewhere-atmosphere-model',
                    type=str,
                    help='Name of atmosphere model elsewhere, e.g., blackbody or NSX-H.')

parser.add_argument('--elsewhere-atmosphere-load',
                    action='store_true',
                    help='Does a numeric atmosphere table need to be loaded from disk for elsewhere?',
                    comment=True)

parser.add_argument('--attenuation-model',
                    type=str,
                    help='Name of interstellar attenuation model, e.g., tbnew.',
                    empty_lines_below=2)

parser.add_argument('--background-model',
                    action='store_true',
                    help='Include an incident background component?',
                    comment=True)

parser.add_argument('--background-shared-instance',
                    action='store_true',
                    help='Do all instruments share the same background model instance?')

parser.add_argument('--background-shared-class',
                    action='store_true',
                    help='Do all instrument models share a background class?')

parser.add_argument('--background-parameters',
                    type=lambda x: ( str(x).replace(' ', '_') ).replace('-', '_'),
                    nargs='*',
                    default=['powerlaw_index', 'powerlaw_normalization'],
                    help='Background model parameter names.',
                    comment=True,
                    inline_comment='enter one name per line below',
                    empty_lines_below=2)

parser.add_argument('--print-MPI-rank',
                    action='store_true',
                    help='Print MPI rank from main module?',
                    comment_line_above='miscellaneous flags',
                    empty_lines_below=2)

parser.add_argument('--config-path',
                    type=str,
                    help='If main module is imported, use this argument to specify the relative or absolute path to the configuration file.',
                    comment_line_above='write flags')

parser.add_argument('--module-directory-path',
                    type=str,
                    help='Absolute path to directory to write module files to.')

parser.add_argument('--main-module',
                    type=str,
                    default='main',
                    help='Name of the main module.')

parser.add_argument('--custom-signal-module',
                    type=str,
                    default='CustomSignal',
                    help='Name of the module containing the CustomSignal subclass.')

parser.add_argument('--custom-instrument-module',
                    type=str,
                    default='CustomInstrument',
                    help='Name of the module containing the CustomInstrument subclass.')

parser.add_argument('--custom-photosphere-module',
                    type=str,
                    default='CustomPhotosphere',
                    help='Name of the module containing the CustomPhotosphere subclass.')

parser.add_argument('--custom-interstellar-module',
                    type=str,
                    default='CustomInterstellar',
                    help='Name of the module containing the CustomInterstellar subclass.')

parser.add_argument('--custom-prior-module',
                    type=str,
                    default='CustomPrior',
                    help='Name of the module containing the CustomPrior subclass.')

parser.add_argument('--custom-background-module',
                    type=str,
                    default='CustomBackground',
                    help='Name of the module containing the CustomBackground subclass(es).')

if __name__ == '__main__':
    if xpsi._verbose:
        print('Parsing configuration file...')
    args, _ = parser.parse_known_args()
    if xpsi._verbose:
        print('Configuration file parsed.')
else:
    if xpsi._verbose:
        print('Parsing configuration file...')
    args, _ = parser.parse_known_args(['@generate.ini'])
    if xpsi._verbose:
        print('Configuration file parsed.')

if len(args.hot_region_model) > 2:
    raise ValueError('A maximum of two hot regions are permitted for module autogeneration.')

_telescopes = args.telescope[0]
for _x in args.telescope[1:]:
    _telescopes += ' x {}'.format(_x)

if args.model is None:
    if len(args.hot_region_model) == 2:
        if args.hot_region_model[0] == args.hot_region_model[1]:
            if args.antipodal_reflection_symmetry:
                _tmp = '{}-S'.format(args.hot_region_model[0])
            else:
                _tmp = '{}-U'.format(args.hot_region_model[0])
        else:
            if args.antipodal_reflection_symmetry:
                raise ValueError('Hot region models are not identical, so antipodal reflection symmetry cannot be imposed.')
            _tmp = '{}+{}'.format(args.hot_region_model[0],
                                  args.hot_region_model[1])
    else:
        _tmp = args.hot_region_model[0]

    args.model = '{} + {}'.format(_tmp,
                                  args.hot_atmosphere_model)
    if args.elsewhere_atmosphere_model is not None:
        args.model += ' + {}'.format(args.elsewhere_atmosphere_model)

# Creating Main module
module = (
'''""" Main module for {} {} <- X-PSI {} {}"""'''.format(_telescopes,
                                                         args.source,
                                                         xpsi.__version__,
                                                         args.model)
)

_telescopes = args.telescope[0]
for _x in args.telescope[1:]:
    _telescopes += ' & {}'.format(_x)

module += (
r'''
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
                _ = '## {{}} ##'.format(comment_line_above)
                entry += '##' + '-' * (len(_) - 4) + '##\n' + _ + '\n##' + '-' * (len(_) - 4) + '##\n'

        _ = ' ## {{}}'.format(inline_comment) if inline_comment is not None else ''
        if not _ and isinstance(nargs, int) and nargs > 1:
            _ = ' ## enter {{}} values below, one per empty line'.format(nargs)
        elif not _ and not isinstance(nargs, int):
            _ = ' ## enter code below, one statement per line'

        if isinstance(default, list):
            for i, _default in enumerate(default):
                entry += '{{5}}{{0}}{{1}}{{2}}{{3}}{{4}}'.format('#' if comment else '',
                                                                 '' if (i > 0 and nargs != 1) else arg,
                                                                 '' if nargs != 1 else '=',
                                                                 '' if (i == 0 and nargs != 1) else str(_default),
                                                                 (_ if i == 0 else '') + ('\n{{0}}{{1}}'.format('#' if comment else '', str(_default)) if i == 0 and nargs != 1 else ''),
                                                                 '\n' if i > 0 else '')
        else:
            entry += '{{0}}{{1}}{{2}}{{3}}{{4}}'.format('#' if comment else '',
                                                        arg,
                                                        '=' if nargs == 1 else '',
                                                        _ if nargs != 1 else (str(default) if default is not None else ''),
                                                        ('\n' + str(default) if default is not None else '') if nargs != 1 else _)

        if action == 'append':
            entry += '\n#{{0}}='.format(arg)

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

        with open('{3}', 'w') as file:
            file.write(config_file)

        print('Configuration file generated.')
        parser.exit()

class NullAction(argparse.Action):
    """ Do not store value in namespace. """
    def __call__(self, parser, namespace, values, option_string=None):
        pass

parser = ArgumentParserCustom(
    description="""
    Main module for X-PSI {0} modelling of {1} {2} event data.

    You can run this module as a script and launch a sampler, optionally
    with a world of MPI processes.

    Alternate usage: mpiexec -n 4 python -m mpi4py %(prog)s [-h] @<config.ini>

    """,
    fromfile_prefix_chars='@')
'''.format(args.model,
           _telescopes,
           args.source,
           args.config_path)
)

module += (
'''\ndef str_to_bool(x):
    if x == 'False':
        return False
    elif x == 'True':
        return True

    raise ValueError('Invalid argument where boolean ``True`` or ``False`` is required.')

'''
)

module += (
'''
parser.add_argument('--generate-config-file', default=argparse.SUPPRESS, action=GenerateConfigAction, help='Generate the configuration file template.',
                         destined_for_config_file=False)
'''
)


module += (
'''
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
'''
)

_path = 'Absolute or relative path to'
_bounds_default_notice = 'If no bounds are given (``None``), and no value is given (``None``), bounds default to the source code strict bounds.'
_value_notice = 'No value means the parameter is free, whilst a value means the parameter is fixed (and the prior may need to be modified manually). If you want the parameter to be derived from other parameters in a complex way, manual modification of the main module is necessary, but we support functions of one parameter in the configuration file.'
_derived_notice = 'If the parameter is derived from one other parameter, e.g., the temperature of a hot region is derived from the temperature of the other hot region, then the value needs to be written using the following template: lambda x: f(x), "parameter", "space". In this template: f(x) is a function of the parameter x from which the value is derived; "parameter" is the name of the parameter x as a string; and "space" is an object in the global namespace, with name written as a string, which is a (sub)space of parameters from which the current value of parameter x can be accessed via getitem magic using "parameter".'
_CDF_notice = 'Supply a function of one variable (the probability mass ``x``), in the form of an expression that can be evaluated with the ``eval(...)`` builtin function, i.e., scipy.stats.truncnorm(x, ...). Note that the prior default PDF is uniform (with compact support), so do not supply a CDF if a uniform prior is desired, or to be explicit, use: DEFAULT UNIFORM. You must use DEFAULT UNIFORM to overwrite a default CDF shown in the auto-generated configuration file, unless the parameter is fixed/derived in which case the prior flag is silently ignored. You can also use the flag more than once: the last usage must be an expression that will be dynamically evaluated using the ``eval(...)`` builtin and must return a float to set as the parameter value; the other usages can be helper statements executed with the ``exec(...)`` builtin, e.g., to set temporary local variables to make the code (and configuration file more readable).'

for instrument in args.instrument:
    module += (
    '''
parser.add_argument('--{0}-exposure-time',
                    type=float,
                    help='{0} exposure time in seconds.',
                    comment_line_above='{0} configuration flags')

parser.add_argument('--{0}-count-matrix-path', type=str, help='{1} {0} channel-phase count matrix. If the data is a spectrum (phase-averaged), then the file must contain a vector. This path is written to if the file does not exist by processing the event files.')
parser.add_argument('--{0}-count-matrix-type', type=str, default='double', help='{0} count matrix NumPy data type.',
                    comment=True)
parser.add_argument('--{0}-event-path', type=str, help='{1} {0} event list file.')
parser.add_argument('--{0}-number-phase-bins', type=int, help='Number of phases bins for binning {0} event list file.')
parser.add_argument('--{0}-event-file-channel-column', type=int, default=1, help='Channel column in {0} event list file.',
                    comment=True)
parser.add_argument('--{0}-event-file-phase-column', type=int, default=2, help='Phase column in {0} event list file.',
                    comment=True)
parser.add_argument('--{0}-event-file-skiprows', type=int, default=3, help='Number of top rows to skip when loading {0} event list file.',
                    comment=True)
parser.add_argument('--{0}-events-in-eV', action='store_true', help='{0} event list file lists events by energy in eV?',
                    comment=True)
parser.add_argument('--{0}-arf-path', type=str, help='{1} {0} ARF file.',
                    comment_line_above='rule')
parser.add_argument('--{0}-effective-area-scaling-factor', type=str, default='1.0',
                    help='Factor by which to scale the nominal effective area model, as a mathematical expression, e.g., a ratio of integers such as 51.0/52.0.',
                    comment=True)
parser.add_argument('--{0}-arf-skiprows', type=int, default=3,
                    help='Number of header rows to skip when loading ARF file.',
                    comment=True)
parser.add_argument('--{0}-arf-low-column', type=int, default=1,
                    help='Column (zero-indexed) containing the low energy edges in the ARF file.',
                    comment=True)
parser.add_argument('--{0}-arf-high-column', type=int, default=2,
                    help='Column (zero-indexed) containing the high energy edges in the ARF file.',
                    comment=True)
parser.add_argument('--{0}-arf-area-column', type=int, default=3,
                    help='Column (zero-indexed) containing the effective area in the ARF file.',
                    comment=True)

parser.add_argument('--{0}-rmf-path', type=str, help='{1} {0} RMF file.',
                    comment_line_above='rule')
parser.add_argument('--{0}-rmf-skiprows', type=int, default=3,
                    help='Number of header rows to skip when loading RMF file.',
                    comment=True)
parser.add_argument('--{0}-rmf-usecol', type=int, default=-1,
                    help='Column (zero-indexed) containing the flattened redistribution matrix elements in the RMF file.',
                    comment=True)

parser.add_argument('--{0}-channels-path', type=str, help='{1} {0} channel-energy mapping file.',
                    comment_line_above='rule')
parser.add_argument('--{0}-channel-energies-skiprows', type=int, default=0,
                    help='Number of header rows to skip when loading channel-energy mapping file.',
                    comment=True)
parser.add_argument('--{0}-channel-energies-low-column', type=int, default=0,
                    help='Column (zero-indexed) containing the low energy edges in the channel-energy mapping file.',
                    comment=True)

parser.add_argument('--{0}-input-bounds',
                    type=int,
                    nargs=2,
                    help='{0} bounding input energy intervals of instrument response submatrix for use with NumPy slice notation.',
                    comment_line_above='rule')

parser.add_argument('--{0}-channel-bounds',
                    type=int,
                    nargs=2,
                    help='{0} bounding channels of instrument response submatrix for use with NumPy slice notation.')

parser.add_argument('--{0}-energy-independent-effective-area-scaling-factor-bounds',
                    type=str,
                    nargs=2,
                    action=CompileAction,
                    help='Bounds for {0} energy-independent effective area scaling factor parameter. If no bounds are given (``None``), and no value is given (``None``), the parameter value is fixed at unity, and the instrument response model is locked to the nominal response model (unless a custom model is implemented).',
                    comment=True,
                    comment_line_above='rule')

parser.add_argument('--{0}-energy-independent-effective-area-scaling-factor-prior',
                    type=str,
                    nargs='*',
                    default=['truncnorm.ppf(x, -3.0, 3.0, loc=1.0, scale=0.1)'],
                    action=NullAction,
                    help='Prior inverse CDF of the energy-independent effective area scaling factor. {5}',
                    comment=True,
                    inline_comment='Normal distribution with std. dev. 10%, truncated at +/- 3 std. dev.')

parser.add_argument('--{0}-energy-independent-effective-area-scaling-factor-value',
                    type=str,
                    action=CompileAction,
                    help='Value for {0} energy-independent effective area scaling parameter. Either the name of an instrument to share the parameter with, as a string, or a float. {3} {4}',
                    comment=True,
                    empty_lines_below=2)

parser.add_argument('--{0}-phase-shift-bounds',
                    type=str,
                    nargs=2,
                    action=CompileAction,
                    help='Bounds for {0} phase-shift parameter. If no bounds are given (``None``), and no value is given (``None``), the parameter value is fixed at zero, and is therefore locked to the phase of the signal specified by the hot region phases. For one phase-resolving instrument, this default behaviour is advised, and additional phase-resolving instruments can in principle have a different fixed, derived, or free phase-shift parameter. For instruments that phase-average, the phase-shift can be arbitrarily fixed or derived, but not free because the likelihood is not a function of it.',
                    comment=True,
                    comment_line_above='rule')

parser.add_argument('--{0}-phase-shift-value',
                    type=str,
                    action=CompileAction,
                    help='Value for {0} phase-shift parameter. Either the name of an instrument to share the parameter with, as a string, or a float. {3} {4}',
                    comment=True,
                    empty_lines_below=2)

parser.add_argument('--{0}-background-prior-support-path', type=str, help='{1} {0} background prior support file. The channel-by-channel lower count-rate limits in the zeroth column, and the upper count-rate limits in the first column. The channels must already match the data.',
                    comment=True,
                    comment_line_above='rule')
parser.add_argument('--{0}-background-skiprows', type=int, default=0, help='Number of top rows to skip when loading {0} background file (prior support or spectrum file).',
                    comment=True)

parser.add_argument('--{0}-background-path', type=str, help='{1} {0} background spectrum file (for imaging telescope).',
                    comment=True)
parser.add_argument('--{0}-background-usecol', type=int, help='Column to use when loading {0} background spectrum file (for imaging telescope).',
                    comment=True)
parser.add_argument('--{0}-background-prior-support-half-width', type=float, help='{0} background prior support half-width (for imaging telescope). The half-width is in units of standard deviation of background count number per instrument channel.',
                    comment=True)
parser.add_argument('--{0}-background-exposure-time',
                    type=float,
                    help='{0} background exposure time in seconds (for imaging telescope).',
                    comment=True)
parser.add_argument('--{0}-background-scaling-factor',
                    type=str,
                    help='{0} background scaling factor, nominally the ratio of on-source CCD extraction area to background CCD extraction area (ideally on same CCD) for imaging telescope. Supply an expression for evaluation by the ``eval(...)`` builtin function.',
                    comment=True,
                    empty_lines_below=2)
    '''.format(instrument.replace('_','-'),
               _path,
               _bounds_default_notice,
               _value_notice,
               _derived_notice,
               _CDF_notice)
    )

if args.background_model:
    for i, _parameter in enumerate(args.background_parameters):
        module += (
        '''
parser.add_argument('--{0}-bounds',
                    type=str,
                    nargs=2,
                    action=CompileAction,
                    help='Bounds of the ``{0}`` parameter. {1}',
                    comment=True,
                    comment_line_above='{4}')

parser.add_argument('--{0}-prior',
                    type=str,
                    nargs='*',
                    action=NullAction,
                    help='Prior inverse CDF of the ``{0}`` parameter. {3}',
                    comment=True)

parser.add_argument('--{0}-value',
                    type=str,
                    action=CompileAction,
                    default='{5}',
                    help='Value of the ``{0}`` parameter. {2}',
                    comment=True,
                    inline_comment='{6}',
                    empty_lines_below=2)
        '''.format(_parameter.replace('_', '-'),
                   _bounds_default_notice,
                   _value_notice,
                   _CDF_notice,
                   'background flags' if i == 0 else 'rule',
                   str(1.0) if 'powerlaw_norm' in _parameter else None,
                   'to allow parameter to be free, uncomment and use: None' if 'powerlaw_norm' in _parameter else None)
)

module += (
'''
parser.add_argument('--attenuation-path', type=str, help='{0} attenuation file.',
                         comment_line_above='attenuation flags')
parser.add_argument('--attenuation-energy-column', type=int, default=0,
                    help='Column (zero-indexed) containing the energies in the attenuation file.',
                    comment=True)
parser.add_argument('--attenuation-column', type=int, default=2,
                    help='Column (zero-indexed) containing the attenuation factors in the attenuation file.',
                    comment=True)

parser.add_argument('--neutral-hydrogen-column-density-bounds',
                    type=str,
                    nargs=2,
                    action=CompileAction,
                    help='Bounds of the neutral hydrogen column density parameter. {1}',
                    comment=True,
                    comment_line_above='rule')

parser.add_argument('--neural-hydrogen-column-density-prior',
                    type=str,
                    nargs='*',
                    action=NullAction,
                    help='Prior inverse CDF of ratio of interstellar neutral hydrogen column density to the fiducial density. {3}',
                    comment=True)

parser.add_argument('--neutral-hydrogen-column-density-value',
                    type=str,
                    action=CompileAction,
                    help='Value of the neutral hydrogen column density parameter. {2}',
                    comment=True,
                    empty_lines_below=2)

'''.format(_path,
           _bounds_default_notice,
           _value_notice,
           _CDF_notice)
)

module += (
'''
parser.add_argument('--mass-bounds',
                    type=str,
                    nargs=2,
                    action=CompileAction,
                    help='Bounds of gravitational mass (solar masses). {1}',
                    comment=True,
                    comment_line_above='spacetime flags')

parser.add_argument('--mass-prior',
                    type=str,
                    nargs='*',
                    action=NullAction,
                    help='Prior inverse CDF of gravitation mass (solar masses). {4}',
                    comment=True)

parser.add_argument('--mass-value',
                    type=str,
                    action=CompileAction,
                    help='Value of gravitational mass (solar masses). {2} {3}',
                    comment=True,
                    empty_lines_below=2)

parser.add_argument('--radius-bounds',
                    type=str,
                    nargs=2,
                    action=CompileAction,
                    help='Bounds of coordinate equatorial radius (km). {1}',
                    comment=True,
                    comment_line_above='rule')

parser.add_argument('--radius-value',
                    type=str,
                    action=CompileAction,
                    help='Value of coordinate equatorial radius (km). {2} {3}',
                    comment=True,
                    empty_lines_below=2)

parser.add_argument('--cos-inclination-bounds',
                    type=str,
                    nargs=2,
                    action=CompileAction,
                    help='Bounds of cosine of Earth colatitude (inclination) w.r.t to stellar rotation axis. {1}',
                    comment=True,
                    comment_line_above='rule')

parser.add_argument('--cos-inclination-prior',
                    type=str,
                    nargs='*',
                    action=NullAction,
                    help='Prior inverse CDF of cosine of Earth inclination to stellar spin axis. {4}',
                    comment=True)

parser.add_argument('--cos-inclination-value',
                    type=str,
                    action=CompileAction,
                    help='Value of cosine of Earth colatitude (inclination) w.r.t to stellar rotation axis. {2}',
                    comment=True,
                    empty_lines_below=2)

parser.add_argument('--distance-bounds',
                    type=str,
                    nargs=2,
                    action=CompileAction,
                    help='Bounds of distance to source (kpc). {1}',
                    comment=True,
                    comment_line_above='rule')

parser.add_argument('--distance-prior',
                    type=str,
                    nargs='*',
                    action=NullAction,
                    help='Prior inverse CDF of distance to source (kpc). {4}',
                    comment=True)

parser.add_argument('--distance-value',
                    type=str,
                    action=CompileAction,
                    help='Value of distance to source (kpc). {2} {3}',
                    comment=True,
                    empty_lines_below=2)

'''.format(_path,
           _bounds_default_notice,
           _value_notice,
           _derived_notice,
           _CDF_notice)
)

for _h, _m in list(zip(args.prefix, args.hot_region_model))[:1 if (len(args.hot_region_model) == 1 or args.antipodal_reflection_symmetry) else 2]:
    module += (
    '''
parser.add_argument('--{0}-super-colatitude-bounds',
                    type=str,
                    nargs=2,
                    action=CompileAction,
                    help='Bounds of hot region "{0}" super-member colatitude w.r.t stellar spin axis (radians). {1}',
                    comment_line_above='"{0}" hot region parameter flags',
                    comment=True)

parser.add_argument('--{0}-super-colatitude-value',
                    type=str,
                    action=CompileAction,
                    help='Value of hot region "{0}" super-member colatitude w.r.t stellar spin axis (radians). {2} {3}',
                    comment=True,
                    empty_lines_below=2)

parser.add_argument('--{0}-super-radius-bounds',
                    type=str,
                    nargs=2,
                    action=CompileAction,
                    help='Bounds of hot region "{0}" super-member angular radius (radians). {1}',
                    comment=True,
                    comment_line_above='rule')

parser.add_argument('--{0}-super-radius-value',
                    type=str,
                    action=CompileAction,
                    help='Value of hot region "{0}" super-member angular radius (radians). {2} {3}',
                    comment=True,
                    empty_lines_below=2)

parser.add_argument('--{0}-super-temperature-bounds',
                    type=str,
                    nargs=2,
                    action=CompileAction,
                    help='Bounds of hot region "{0}" super-member log10(temperature [K]). {1}',
                    comment=True,
                    comment_line_above='rule')

parser.add_argument('--{0}-super-temperature-prior',
                    type=str,
                    nargs='*',
                    action=NullAction,
                    help='Prior inverse CDF of hot-region {0} superseding region log10(temperature [K]). {4}',
                    comment=True)

parser.add_argument('--{0}-super-temperature-value',
                    type=str,
                    action=CompileAction,
                    help='Value of hot region "{0}" super-member log10(temperature [K]). {2} {3}',
                    comment=True,
                    empty_lines_below=2)
    '''.format(_h,
               _bounds_default_notice,
               _value_notice,
               _derived_notice,
               _CDF_notice)
    )

    if _m in ['CST', 'EST', 'PST']:
        module += (
        '''
parser.add_argument('--{0}-omit-radius-bounds',
                    type=str,
                    nargs=2,
                    action=CompileAction,
                    help='Bounds of hot region "{0}" omit-member angular radius (radians). {1}',
                    comment=True,
                    comment_line_above='rule')

parser.add_argument('--{0}-omit-radius-value',
                    type=str,
                    action=CompileAction,
                    help='Value of hot region "{0}" omit-member angular radius (radians). {2} {3}',
                    comment=True,
                    empty_lines_below=2)
        '''.format(_h,
                   _bounds_default_notice,
                   _value_notice,
                   _derived_notice)
        )

    if _m in ['EST', 'PST']:
        module += (
        '''
parser.add_argument('--{0}-omit-colatitude-bounds',
                    type=str,
                    nargs=2,
                    action=CompileAction,
                    help='Bounds of hot region "{0}" omit-member colatitude w.r.t stellar spin axis (radians). {1}',
                    comment=True,
                    comment_line_above='rule')

parser.add_argument('--{0}-omit-colatitude-value',
                    type=str,
                    action=CompileAction,
                    help='Value of hot region "{0}" omit-member colatitude w.r.t stellar spin axis (radians). {2} {3}',
                    comment=True,
                    empty_lines_below=2)

parser.add_argument('--{0}-omit-azimuth-bounds',
                    type=str,
                    nargs=2,
                    action=CompileAction,
                    help='Bounds of hot region "{0}" omit-member azimuth relative to super-member (radians). {1}',
                    comment=True,
                    comment_line_above='rule')

parser.add_argument('--{0}-omit-azimuth-value',
                    type=str,
                    action=CompileAction,
                    help='Value of hot region "{0}" omit-member azimuth relative to super-member (radians). {2} {3}',
                    comment=True,
                    empty_lines_below=2)
        '''.format(_h,
                   _bounds_default_notice,
                   _value_notice,
                   _derived_notice)
        )

    if _m in ['CDT', 'EDT', 'PDT']:
        module += (
        '''
parser.add_argument('--{0}-cede-radius-bounds',
                    type=str,
                    nargs=2,
                    action=CompileAction,
                    help='Bounds of hot region "{0}" cede-member angular radius (radians). {1}',
                    comment=True,
                    comment_line_above='rule')

parser.add_argument('--{0}-cede-radius-value',
                    type=str,
                    action=CompileAction,
                    help='Value of hot region "{0}" cede-member angular radius (radians). {2} {3}',
                    comment=True,
                    empty_lines_below=2)

parser.add_argument('--{0}-cede-temperature-bounds',
                    type=str,
                    nargs=2,
                    action=CompileAction,
                    help='Bounds of hot region "{0}" cede-member log10(temperature [K]). {1}',
                    comment=True,
                    comment_line_above='rule')

parser.add_argument('--{0}-cede-temperature-prior',
                    type=str,
                    nargs='*',
                    action=NullAction,
                    help='Prior inverse CDF of hot-region {0} ceding region log10(temperature [K]). {4}',
                    comment=True)

parser.add_argument('--{0}-cede-temperature-value',
                    type=str,
                    action=CompileAction,
                    help='Value of hot region "{0}" cede-member log10(temperature [K]). {2} {3}',
                    comment=True,
                    empty_lines_below=2)
        '''.format(_h,
                   _bounds_default_notice,
                   _value_notice,
                   _derived_notice,
                   _CDF_notice)
        )

    if _m in ['EDT', 'PDT']:
        module += (
        '''
parser.add_argument('--{0}-cede-colatitude-bounds',
                    type=str,
                    nargs=2,
                    action=CompileAction,
                    help='Bounds of hot region "{0}" cede-member colatitude w.r.t stellar spin axis (radians). {1}',
                    comment=True,
                    comment_line_above='rule')

parser.add_argument('--{0}-cede-colatitude-value',
                    type=str,
                    action=CompileAction,
                    help='Value of hot region "{0}" cede-member colatitude w.r.t stellar spin axis (radians). {2} {3}',
                    comment=True,
                    empty_lines_below=2)

parser.add_argument('--{0}-cede-azimuth-bounds',
                    type=str,
                    nargs=2,
                    action=CompileAction,
                    help='Bounds of hot region "{0}" cede-member azimuth relative to super-member (radians). {1}',
                    comment=True,
                    comment_line_above='rule')

parser.add_argument('--{0}-cede-azimuth-value',
                    type=str,
                    action=CompileAction,
                    help='Value of hot region "{0}" cede-member azimuth relative to super-member (radians). {2} {3}',
                    comment=True,
                    empty_lines_below=2)
        '''.format(_h,
                   _bounds_default_notice,
                   _value_notice,
                   _derived_notice)
        )

    module += (
    '''
parser.add_argument('--{0}-phase-shift-bounds',
                    type=str,
                    nargs=2,
                    action=CompileAction,
                    help='Bounds of hot region "{0}" phase shift (cycles). {1}',
                    comment=True,
                    comment_line_above='rule')

parser.add_argument('--{0}-phase-shift-value',
                type=str,
                action=CompileAction,
                help='Value of hot region "{0}" phase shift (cycles). {2} {3}',
                comment=True,
                empty_lines_below=2)
    '''.format(_h,
               _bounds_default_notice,
               _value_notice,
               _derived_notice)
    )

for _h in args.prefix:
    module += (
    '''
parser.add_argument('--{0}-sqrt-num-cells',
                    type=int,
                    default=32,
                    help='Target square-root of the number of cells spanning (raditing subset of) hot region "{0}".',
                    comment_line_above='"{0}" hot region resolution flags')

parser.add_argument('--{0}-min-sqrt-num-cells',
                    type=int,
                    default=10,
                    help='Minimum square-root of the number of cells constituting hot region "{0}" mesh.')

parser.add_argument('--{0}-max-sqrt-num-cells',
                    type=int,
                    default=80,
                    help='Maximum square-root of the number of cells constituting hot region "{0}" mesh.')

parser.add_argument('--{0}-num-leaves',
                    type=int,
                    default=64,
                    help='Number of phases on unit interval at which to compute hot region "{0}" signal.')

parser.add_argument('--{0}-num-rays',
                    type=int,
                    default=512,
                    help='Number of rays per iso-latitude mesh subset to trace when computing hot region "{0}" signal.',
                    empty_lines_below=2)
    '''.format(_h)
    )

if args.hot_atmosphere_load:
    module += (
    '''
parser.add_argument('--hot-atmosphere-path', type=str, help='{0} hot atmosphere file.',
                             comment_line_above='hot atmosphere flags')

parser.add_argument('--hot-atmosphere-size',
                    type=int,
                    nargs=4,
                    action=NullAction,
                    help='Size of each of the four dimensions of the numeric atmosphere table for the hot regions.',
                    empty_lines_below=2)
    '''.format(_path)
    )

if args.elsewhere_atmosphere_model is not None:
    module += (
    '''
parser.add_argument('--elsewhere-temperature-bounds',
                    type=str,
                    nargs=2,
                    action=CompileAction,
                    help='Bounds of log10(temperature [K]) elsewhere. {0}',
                    comment_line_above='elsewhere flags',
                    comment=True)

parser.add_argument('--elsewhere-temperature-value',
                    type=str,
                    action=CompileAction,
                    help='Value of log10(temperature [K]) elsewhere. {1} {2}',
                    comment=True,
                    empty_lines_below=2)

parser.add_argument('--elsewhere-sqrt-num-cells',
                    type=int,
                    default=64,
                    help='Target square-root of the number of cells spanning (raditing subset of) the elsewhere region.')

parser.add_argument('--elsewhere-num-rays',
                    type=int,
                    default=1024,
                    help='Number of rays per iso-latitude mesh subset to trace when computing the signal from elsewhere.',
                    empty_lines_below={3:d})
    '''.format(_bounds_default_notice,
               _value_notice,
               _derived_notice,
               0 if args.elsewhere_atmosphere_load else 2)
    )

    if args.elsewhere_atmosphere_load:
        module += (
        '''
parser.add_argument('--elsewhere-atmosphere-path', type=str, help='{0} Elsewhere atmosphere file.')

parser.add_argument('--elsewhere-atmosphere-size',
                    type=int,
                    nargs=4,
                    action=NullAction,
                    help='Size of each of the four dimensions of the numeric atmosphere table for elsewhere.',
                    empty_lines_below=2)
        '''.format(_bounds_default_notice)
        )

module += (
'''
parser.add_argument('--image-order-limit',
                         type=int,
                         default=None,
                         help='The highest-order image to sum over. Either a positive integer, or do not pass an argument if a hard limit is not desired.',
                         comment_line_above='global resolution flags')
'''
)

module += (
'''
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
'''
)

module += (
'''
parser.add_argument('--multinest', action='store_true', help='Launch MultiNest sampler if module is executed.',
                         comment_line_above='runtime flags')

parser.add_argument('--resume', action='store_true', help='Resume sampling if module is executed.')


parser.add_argument('--sample-files-directory-path',
                    type=str,
                    default='samples/',
                    help='{} sample file directory. If no path is provided, the default (relative) path is "samples/".')

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

'''.format(_path)
)


module += (
'''
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
    args, _ = parser.parse_known_args(['@{}'])
    if xpsi._verbose:
        print('Configuration file parsed.')

'''.format(args.config_path)
)

module += (
'''
import os

import numpy as np
import math

from xpsi.Parameter import Derive
from xpsi import HotRegions
'''
)

if args.print_MPI_rank:
    module += '''\nprint('Rank reporting: %d' % xpsi._rank) '''

module += (
'''
if __name__ == '__main__':
    from {0} import CustomInstrument
    from {1} import CustomSignal
    from {2} import CustomInterstellar

    try:
        from {3} import CustomPhotosphere
    except ImportError:
        from xpsi import Photosphere as CustomPhotosphere

    from {4} import CustomPrior

else:
    from .{0} import CustomInstrument
    from .{1} import CustomSignal
    from .{2} import CustomInterstellar

    try:
        from .{3} import CustomPhotosphere
    except ImportError:
        from xpsi import Photosphere as CustomPhotosphere

    from .{4} import CustomPrior
'''.format(args.custom_instrument_module,
           args.custom_signal_module,
           args.custom_interstellar_module,
           args.custom_photosphere_module,
           args.custom_prior_module)
)

if args.background_shared_instance:
    args.background_shared_class = True
    args.background_shared_parameters = True
elif args.background_shared_class:
    args.background_shared_parameters = True
elif not args.background_shared_class:
    args.background_shared_parameters = False

if args.background_shared_class:
    module += (
    '''
if __name__ == '__main__':
    from {0} import CustomBackground
else:
    from .{0} import CustomBackground
    '''.format(args.custom_background_module)
    )
elif args.background_shared_class:
    for _instrument in args.instruments:
        module += (
        '''
if __name__ == '__main__':
    from {0} import {1}_CustomBackground
        '''.format(args.custom_background_module,
                   _instrument)
        )
        module += (
        '''
else:
    from .{0} import {1}_CustomBackground
        '''.format(args.custom_background_module,
                   _instrument)
        )

module += (
'''
if args.main_import_statements is not None:
    for import_statement in args.main_import_statements:
        exec(import_statement)

if args.main_global_statements is not None:
    for global_statement in args.main_global_statements:
        exec(global_statement)
'''
)

module += (
'''
class namespace():
    pass
'''
)

module += (
'''
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
'''
)

module += (
'''
bounds = dict(neutral_hydrogen_column_density = parse_bounds(args.neutral_hydrogen_column_density_bounds,
                                                              args.neutral_hydrogen_column_density_value))
values = dict(neutral_hydrogen_column_density = parse_value(args.neutral_hydrogen_column_density_value))

interstellar = CustomInterstellar.load(args.attenuation_path,
                                       args.attenuation_energy_column,
                                       args.attenuation_column,
                                       bounds = bounds,
                                       values = values)
'''
)

module += (
'''
signals = [[],]
'''
)

if args.background_model:
    if args.background_shared_instance:
        module += (
        '''
bounds = {}
values = {}
        '''
        )

        for _parameter in args.background_parameters:
            module += (
            '''
bounds['{0}'] = parse_bounds(args.{0}_bounds, args.{0}_value)
values['{0}'] = parse_value(args.{0}_value)
            '''.format(_parameter)
            )

        module += (
        '''
background = CustomBackground(bounds=bounds, values=values)
        '''
    )

for instrument in args.instrument:
    module += (
    '''
{0} = namespace()

if args.{0}_{2}_value is not None:
    if eval(args.{0}_{2}_value) in {1}:
        values = dict({2} = derived_parameter(lambda x: x,
                                              '{2}',
                                              eval(args.{0}_{2}_value) + '.instrument'))
    else:
        values = dict({2} = parse_value(args.{0}_{2}_value))
else:
    values = {{}}

bounds = dict({2} = parse_bounds(args.{0}_{2}_bounds,
                                   args.{0}_{2}_value,
                                   default_to_free = False))

{0}.instrument = CustomInstrument.{0}(bounds=bounds,
                                  values=values,
                                  ARF=args.{0}_arf_path,
                                  RMF=args.{0}_rmf_path,
                                  channel_energies=args.{0}_channels_path,
                                  max_input=args.{0}_input_bounds[1],
                                  max_channel=args.{0}_channel_bounds[1],
                                  min_input=args.{0}_input_bounds[0],
                                  min_channel=args.{0}_channel_bounds[0],
                                  effective_area_scaling_factor=eval(args.{0}_effective_area_scaling_factor),
                                  ARF_skiprows=args.{0}_arf_skiprows,
                                  ARF_low_column=args.{0}_arf_low_column,
                                  ARF_high_column=args.{0}_arf_high_column,
                                  ARF_area_column=args.{0}_arf_area_column,
                                  RMF_skiprows=args.{0}_rmf_skiprows,
                                  RMF_usecol=args.{0}_rmf_usecol,
                                  channel_energies_skiprows=args.{0}_channel_energies_skiprows,
                                  channel_energies_low_column=args.{0}_channel_energies_low_column)

try:
    counts = np.loadtxt(args.{0}_count_matrix_path, dtype=np.double)
except ValueError:
    {0}.data = xpsi.Data.bin__event_list(args.{0}_event_path,
                                         channels={0}.instrument.channels,
                                         phases=np.linspace(0.0, 1.0, args.{0}_number_phase_bins + 1),
                                         channel_column=args.{0}_event_file_channel_column,
                                         phase_column=args.{0}_event_file_phase_column if args.{0}_number_phase_bins > 1 else None,
                                         phase_averaged=True if args.{0}_number_phase_bins == 1 else False,
                                         channel_edges={0}.instrument.channel_edges,
                                         skiprows=args.{0}_event_file_skiprows,
                                         eV=True if args.{0}_events_in_eV else False,
                                         dtype=getattr(np, args.{0}_count_matrix_type),
                                         first=0,
                                         last=len({0}.instrument.channels) - 1,
                                         exposure_time=args.{0}_exposure_time)

    np.savetxt(args.{0}_event_path.replace('.txt','_converted_to_counts.txt'), {0}.data.counts)
    print('Counts file saved as: '+args.{0}_event_path.replace('.txt','_converted_to_counts.txt'))
    print('Update configuration file to take in counts file to save computation time.')
else:
    if counts.ndim == 1:
        counts = counts.reshape(-1,1)

    {0}.data = xpsi.Data(counts,
                           channels={0}.instrument.channels,
                           phases=np.linspace(0.0, 1.0, args.{0}_number_phase_bins + 1),
                           first=0,
                           last=len({0}.instrument.channels) - 1,
                           exposure_time=args.{0}_exposure_time)

if args.{0}_background_prior_support_path:
    support = np.loadtxt(args.{0}_background_path,
                         skiprows=args.{0}_background_skiprows,
                         dtype=np.double)

elif args.{0}_background_path:
    spectrum = np.loadtxt(args.{0}_background_path,
                          skiprows=args.{0}_background_skiprows,
                          usecols=args.{0}_background_usecol,
                          dtype=np.double)[{0}.instrument.channels]

    support = np.zeros((len(spectrum), 2), dtype=np.double)
    support[:,0] = spectrum - args.{0}_background_prior_support_half_width * np.sqrt(spectrum)
    support[support[:,0] < 0.0, 0] = 0.0
    support[:,1] = spectrum + args.{0}_background_prior_support_half_width * np.sqrt(spectrum)

    for i in range(support.shape[0]):
        if support[i,1] == 0.0:
            for j in range(1, support.shape[0]):
                if i+j < support.shape[0] and support[i+j,1] > 0.0:
                    support[i,1] = support[i+j,1]
                    break
                elif i-j >= 0 and support[i-j,1] > 0.0:
                    support[i,1] = support[i-j,1]
                    break

    support *= ({0}.data.exposure_time / args.{0}_background_exposure_time) * float(eval(args.{0}_background_scaling_factor)) # exposure ratio * scaling

    support /= {0}.data.exposure_time # need count rate, so divide by exposure time
else:
    support = None
    '''.format(instrument,
               args.instrument,
               'energy_independent_effective_area_scaling_factor')
    )

    if not args.background_model:
        module += (
        '''
{0}.background = None
        '''.format(instrument)
        )
    else:
        if args.background_shared_instance:
            module += (
            '''
{0}.background = background
            '''.format(instrument)
            )
        elif args.background_shared_class:
            if instrument == args.instrument[0]:
                module += (
                '''
bounds = {}
values = {}
                '''
                )
                for _parameter in args.background_parameters:
                    module += (
                    '''
bounds['{0}'] = parse_bounds(args.{0}_bounds, args.{0}_value)
values['{0}'] = parse_value(args.{0}_value)
                    '''.format(_parameter)
                    )
            else:
                module += (
                '''
bounds = None
values = {}
                '''
                )
                for _parameter in args.background_parameters:
                    module += (
                    '''
values['{0}'] = derived_parameter(lambda x: x,
                                  '{0}',
                                  {1}.background)
                    '''.format(_parameter, args.instrument[0])
                    )

            module += (
            '''
{0}.background = CustomBackground(bounds=bounds, values=values)
            '''.format(instrument)
            )
        elif not args.background_shared_class:
            module += (
            '''
bounds = {}
values = {}
            '''
            )

            for _parameter in args.background_parameters:
                if instrument in _parameter:
                    _parameter = _parameter[len(instrument):].lstrip('_')

                module += (
                '''
bounds['{0}'] = parse_bounds(args.{0}_bounds, args.{0}_value)
values['{0}'] = parse_value(args.{0}_value)
                '''.format(_parameter)
                )

            module += (
            '''
{0}.background = {0}_CustomBackground(bound=bounds, values=values)
            '''
            )

    module += (
    '''
if args.{0}_phase_shift_value is not None:
    if eval(args.{0}_phase_shift_value) in {1}:
        values = dict(phase_shift = derived_parameter(lambda x: x,
                                                'phase_shift',
                                                eval(args.{0}_phase_shift_value) + '.signal'))
    else:
        values = dict(phase_shift = parse_value(args.{0}_phase_shift_value))
else:
    values = {{}}

bounds = dict(phase_shift = parse_bounds(args.{0}_phase_shift_bounds,
                                         args.{0}_phase_shift_value,
                                         default_to_free=False))

{0}.signal = CustomSignal(data = {0}.data,
                          instrument = {0}.instrument,
                          interstellar = interstellar,
                          background = {0}.background,
                          cache = False if __name__ == '__main__' else True,
                          bounds = bounds,
                          values = values,
                          workspace_intervals = 1000,
                          epsrel = 1.0e-8,
                          epsilon = 1.0e-3,
                          sigmas = 10.0,
                          support = support,
                          prefix = '{0}')

signals[0].append({0}.signal)
    '''.format(instrument,
               args.instrument,
               'energy_independent_effective_area_scaling_factor')
    )

module += (
'''
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
              frequency = {0})

spacetime = xpsi.Spacetime(bounds, values)
'''.format(str(args.frequency))
)

def get_bounds_and_values(prefix, model):
    global module

    module += (
    '''
bounds = dict(super_colatitude = parse_bounds(args.{0}_super_colatitude_bounds,
                                              args.{0}_super_colatitude_value),
              super_radius = parse_bounds(args.{0}_super_radius_bounds,
                                          args.{0}_super_radius_value),
              phase_shift = parse_bounds(args.{0}_phase_shift_bounds,
                                         args.{0}_phase_shift_value),
              super_temperature = parse_bounds(args.{0}_super_temperature_bounds,
                                               args.{0}_super_temperature_value))

values = dict(super_colatitude = parse_value(args.{0}_super_colatitude_value),
              super_radius = parse_value(args.{0}_super_radius_value),
              phase_shift = parse_value(args.{0}_phase_shift_value),
              super_temperature = parse_value(args.{0}_super_temperature_value))
    '''.format(prefix)
    )

    if model in ['CST','EST','PST']:
        module += (
        '''
bounds['omit_radius'] = parse_bounds(args.{0}_omit_radius_bounds,
                                     args.{0}_omit_radius_value)

values['omit_radius'] = parse_value(args.{0}_omit_radius_value)
        '''.format(prefix)
        )

    if model in ['EST', 'PST']:
        module += (
        '''
bounds['omit_colatitude'] = parse_bounds(args.{0}_omit_colatitude_bounds,
                                         args.{0}_omit_colatitude_value)

values['omit_colatitude'] = parse_value(args.{0}_omit_colatitude_value)

bounds['omit_azimuth'] = parse_bounds(args.{0}_omit_azimuth_bounds,
                                      args.{0}_omit_azimuth_value)

values['omit_azimuth'] = parse_value(args.{0}_omit_azimuth_value)
        '''.format(prefix)
        )

    if 'DT' in model:
        module += (
        '''
bounds['cede_radius'] = parse_bounds(args.{0}_cede_radius_bounds,
                                     args.{0}_cede_radius_value)

values['cede_radius'] = parse_value(args.{0}_cede_radius_value)

bounds['cede_temperature'] = parse_bounds(args.{0}_cede_temperature_bounds,
                                          args.{0}_cede_temperature_value)

values['cede_temperature'] = parse_value(args.{0}_cede_temperature_value)
        '''.format(prefix)
        )

    if model in ['EDT', 'PDT']:
        module += (
        '''
bounds['cede_colatitude'] = parse_bounds(args.{0}_cede_colatitude_bounds,
                                         args.{0}_cede_colatitude_value)

values['cede_colatitude'] = parse_value(args.{0}_cede_colatitude_value)

bounds['cede_azimuth'] = parse_bounds(args.{0}_cede_azimuth_bounds,
                                      args.{0}_cede_azimuth_value)

values['cede_azimuth'] = parse_value(args.{0}_cede_azimuth_value)
        '''.format(prefix)
        )

def get_member_settings(model):
    global module

    if model == 'ST':
        module += (
        '''
symmetry = True
omit = False
cede = False
concentric = False
        '''
        )
    elif model == 'CST':
        module += (
        '''
symmetry = True
omit = True
cede = False
concentric = True
        '''
        )
    elif model in ['EST','PST']:
        module += (
        '''
symmetry = True
omit = True
cede = False
concentric = False
        '''
        )
    elif model == 'CDT':
        module += (
        '''
symmetry = True
omit = False
cede = True
concentric = True
        '''
        )
    elif model in ['EDT', 'PDT']:
        module += (
        '''
symmetry = True
omit = False
cede = True
concentric = False
        '''
        )

get_member_settings(args.hot_region_model[0])
get_bounds_and_values(args.prefix[0], args.hot_region_model[0])

if 'NSX' in args.hot_atmosphere_model or 'nsx' in args.hot_atmosphere_model:
    module += (
    '''
primary = xpsi.HotRegion(bounds=bounds,
                            values=values,
                            symmetry=symmetry,
                            omit=omit,
                            cede=cede,
                            concentric=concentric,
                            sqrt_num_cells=args.{1}_sqrt_num_cells,
                            min_sqrt_num_cells=args.{1}_min_sqrt_num_cells,
                            max_sqrt_num_cells=args.{1}_max_sqrt_num_cells,
                            num_leaves=args.{1}_num_leaves,
                            num_rays=args.{1}_num_rays,
                            is_antiphased={0},
                            image_order_limit=args.image_order_limit,
                            atm_ext="Num4D",
                            prefix='{1}')
    '''.format(str(args.is_antiphased[0]), args.prefix[0])
    )

else:
    module += (
    '''
primary = xpsi.HotRegion(bounds=bounds,
                            values=values,
                            symmetry=symmetry,
                            omit=omit,
                            cede=cede,
                            concentric=concentric,
                            sqrt_num_cells=args.{1}_sqrt_num_cells,
                            min_sqrt_num_cells=args.{1}_min_sqrt_num_cells,
                            max_sqrt_num_cells=args.{1}_max_sqrt_num_cells,
                            num_leaves=args.{1}_num_leaves,
                            num_rays=args.{1}_num_rays,
                            is_antiphased={0},
                            image_order_limit=args.image_order_limit,
                            atm_ext="BB",
                            prefix='{1}')
    '''.format(str(args.is_antiphased[0]), args.prefix[0])
    )

if len(args.hot_region_model) == 2:

    if args.antipodal_reflection_symmetry:
        bounds = {}

        module += (
        '''
values = dict(super_colatitude = derived_parameter(lambda x: math.pi - x, '{0}__super_colatitude', 'primary'),
              super_radius = derived_parameter(lambda x: x, '{0}__super_radius', 'primary'),
              phase_shift = derived_parameter(lambda x: x, '{0}__phase_shift', 'primary'),
              super_temperature = derived_parameter(lambda x: x, '{0}__super_temperature', 'primary'),
              omit_colatitude = derived_parameter(lambda x: math.pi - x, '{0}__omit_colatitude', 'primary'),
              omit_radius = derived_parameter(lambda x: x, '{0}__omit_radius', 'primary'),
              omit_azimuth = derived_parameter(lambda x: x, '{0}__omit_azimuth', 'primary'),
              cede_colatitude = derived_parameter(lambda x: math.pi - x, '{0}__cede_colatitude', 'primary'),
              cede_radius = derived_parameter(lambda x: x, '{0}__cede_radius', 'primary'),
              cede_azimuth = derived_parameter(lambda x: x, '{0}__cede_azimuth', 'primary'),
              cede_temperature = derived_parameter(lambda x: x, '{0}__cede_temperature', 'primary'))
        '''.format(args.prefix[0])
        )

    else:
        get_member_settings(args.hot_region_model[1])
        get_bounds_and_values(args.prefix[1], args.hot_region_model[1])

    if 'NSX' in args.hot_atmosphere_model or 'nsx' in args.hot_atmosphere_model:
        module += (
        '''
secondary = xpsi.HotRegion(bounds=bounds,
                                values=values,
                                symmetry=symmetry,
                                omit=omit,
                                cede=cede,
                                concentric=concentric,
                                sqrt_num_cells=args.{1}_sqrt_num_cells,
                                min_sqrt_num_cells=args.{1}_min_sqrt_num_cells,
                                max_sqrt_num_cells=args.{1}_max_sqrt_num_cells,
                                num_leaves=args.{1}_num_leaves,
                                num_rays=args.{1}_num_rays,
                                is_antiphased={0},
                                image_order_limit=args.image_order_limit,
                                atm_ext="Num4D",
                                prefix='{1}')
        '''.format(str(not args.is_antiphased[0] if args.antipodal_reflection_symmetry else args.is_antiphased[1]), args.prefix[1])
        )
    
    else:
        module += (
        '''
secondary = xpsi.HotRegion(bounds=bounds,
                                values=values,
                                symmetry=symmetry,
                                omit=omit,
                                cede=cede,
                                concentric=concentric,
                                sqrt_num_cells=args.{1}_sqrt_num_cells,
                                min_sqrt_num_cells=args.{1}_min_sqrt_num_cells,
                                max_sqrt_num_cells=args.{1}_max_sqrt_num_cells,
                                num_leaves=args.{1}_num_leaves,
                                num_rays=args.{1}_num_rays,
                                is_antiphased={0},
                                image_order_limit=args.image_order_limit,
                                atm_ext="BB",
                                prefix='{1}')
        '''.format(str(not args.is_antiphased[0] if args.antipodal_reflection_symmetry else args.is_antiphased[1]), args.prefix[1])
        )        

    module += (
    '''
hot = HotRegions((primary, secondary))
    '''
    )
else:
    module += (
    '''
hot = primary
    '''
    )

if args.elsewhere_atmosphere_model is not None:
    if 'NSX' in args.elsewhere_atmosphere_model or 'nsx' in args.elsewhere_atmosphere_model:
        module += (
        '''
elsewhere = xpsi.Elsewhere(bounds=dict(elsewhere_temperature = parse_bounds(args.elsewhere_temperature_bounds,
                                                                                    args.elsewhere_temperature_value)),
                                    values=dict(elsewhere_temperature = parse_value(args.elsewhere_temperature_value)),
                                    sqrt_num_cells=args.elsewhere_sqrt_num_cells,
                                    num_rays=args.elsewhere_num_rays,
                                    image_order_limit=args.image_order_limit,
                                    atm_ext="Num4D")
        '''
        )
    else:
        module += (
        '''
elsewhere = xpsi.Elsewhere(bounds=dict(elsewhere_temperature = parse_bounds(args.elsewhere_temperature_bounds,
                                                                                    args.elsewhere_temperature_value)),
                                    values=dict(elsewhere_temperature = parse_value(args.elsewhere_temperature_value)),
                                    sqrt_num_cells=args.elsewhere_sqrt_num_cells,
                                    num_rays=args.elsewhere_num_rays,
                                    image_order_limit=args.image_order_limit,
                                    atm_ext="BB")
        '''
        )        
else:
    module += (
    '''
elsewhere = None
    '''
    )


module += (
'''
photosphere = CustomPhotosphere(hot = hot,
                                     elsewhere = elsewhere,
                                     values = dict(mode_frequency = spacetime['frequency']))
'''
)

if args.hot_atmosphere_load:
    module += (
    '''
photosphere.hot_atmosphere = args.hot_atmosphere_path
    '''
    )

if args.elsewhere_atmosphere_model is not None and args.elsewhere_atmosphere_load:
    module += (
    '''
photosphere.elsewhere_atmosphere = args.elsewhere_atmosphere_path
    '''
    )

module += (
'''
star = xpsi.Star(spacetime = spacetime, photospheres = photosphere)

prior = CustomPrior()

likelihood = xpsi.Likelihood(star = star,
                             signals = signals,
                             num_energies = args.number_energies,
                             threads = args.openmp_threads,
                             externally_updated = args.parameters_externally_updated if args.parameters_externally_updated is not None else True,
                             prior = prior,
                             max_energy = args.maximum_energy_ray_tracing)

'''
)

module += (
'''
if __name__ == '__main__':

    if args.multinest:

        wrapped_params = [0] * len(likelihood)
        for name in likelihood.names:
            if ('phase_shift' in name) or ('azimuth' in name):
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

'''
)

if not os.path.isdir(args.module_directory_path):
    os.mkdir(args.module_directory_path)

write(r'{}.py'.format(os.path.join(args.module_directory_path,
                                   args.main_module)), module)

write(r'{}.py'.format(os.path.join(args.module_directory_path, '__init__')), '')

# Creating Signal module
module = (
'''""" Signal module for X-PSI {0} modelling of {1} {2} event data. """

import numpy as np
import math

import xpsi

from xpsi.likelihoods.default_background_marginalisation import eval_marginal_likelihood
from xpsi.likelihoods.default_background_marginalisation import precomputation

class CustomSignal(xpsi.Signal):
    """ A subclass of the :class:`xpsi.Signal.Signal` class to make it callable.

    """

    def __init__(self, workspace_intervals = 1000, epsabs = 0, epsrel = 1.0e-8,
                 epsilon = 1.0e-3, sigmas = 10.0, support = None, *args, **kwargs):
        """ Perform precomputation. """

        super(CustomSignal, self).__init__(*args, **kwargs)

        try:
            self._precomp = precomputation(self._data.counts.astype(np.int32))
        except AttributeError:
            print('No data... can synthesise data but cannot evaluate a '
                  'likelihood function.')
        else:
            self._workspace_intervals = workspace_intervals
            self._epsabs = epsabs
            self._epsrel = epsrel
            self._epsilon = epsilon
            self._sigmas = sigmas

            if support is not None:
                self._support = support
            else:
                self._support = -1.0 * np.ones((self._data.counts.shape[0],2))
                self._support[:,0] = 0.0

    @property
    def support(self):
        return self._support

    @support.setter
    def support(self, obj):
        self._support = obj

    def __call__(self, *args, **kwargs):
        self.loglikelihood, self.expected_counts, self.background_signal, self.background_signal_given_support = \\
                eval_marginal_likelihood(self._data.exposure_time,
                                          self._data.phases,
                                          self._data.counts,
                                          self._signals,
                                          self._phases,
                                          self._shifts,
                                          self._precomp,
                                          self._support,
                                          self._workspace_intervals,
                                          self._epsabs,
                                          self._epsrel,
                                          self._epsilon,
                                          self._sigmas,
                                          kwargs.get('llzero'),
                                          background={3})
'''.format(args.model,
           _telescopes,
           args.source,
           'None' if not args.background_model else 'self.background.registered_background')
)

write(r'{}.py'.format(os.path.join(args.module_directory_path, args.custom_signal_module)), module)

# Creating Photosphere module
module = (
r'''""" Photosphere module for X-PSI {0} modelling of {1} {2} event data. """

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
    Photosphere module for X-PSI {0} modelling of {1} {2} event data.

    You should import this module.

    For help: python %(prog)s -h

    """,
    fromfile_prefix_chars='@')
'''.format(args.model,
           _telescopes,
           args.source)
)

if args.hot_atmosphere_load:
    module += (
    '''
parser.add_argument('--hot-atmosphere-size',
                             type=int,
                             nargs=4,
                             help='Size of each of the four dimensions of the numeric atmosphere table for the hot regions.')
    '''
    )

if args.elsewhere_atmosphere_model is not None:
    if args.elsewhere_atmosphere_load:
        module += (
        '''
parser.add_argument('--elsewhere-atmosphere-size',
                                 type=int,
                                 nargs=4,
                                 help='Size of each of the four dimensions of the numeric atmosphere table for elsewhere.')
        '''
        )

module += (
'''
if __name__ == '__main__':
    args, _ = parser.parse_known_args()
else:
    args, _ = parser.parse_known_args(['@{}'])
'''.format(args.config_path)
)

module += (
'''
import numpy as np
import math

import xpsi

class CustomPhotosphere(xpsi.Photosphere):
'''
)

if args.hot_atmosphere_load:
    module += (
    '''
    @xpsi.Photosphere.hot_atmosphere.setter
    def hot_atmosphere(self, path):
    '''
    )

    if ('NSX' in args.hot_atmosphere_model or 'nsx' in args.hot_atmosphere_model):
        module += (
        '''
        table = np.loadtxt(path, dtype=np.double)
        logT = np.zeros(args.hot_atmosphere_size[0])
        logg = np.zeros(args.hot_atmosphere_size[1])
        mu = np.zeros(args.hot_atmosphere_size[2])
        logE = np.zeros(args.hot_atmosphere_size[3])

        reorder_buf = np.zeros((args.hot_atmosphere_size[0],
                                args.hot_atmosphere_size[1],
                                args.hot_atmosphere_size[2],
                                args.hot_atmosphere_size[3],))

        index = 0
        for i in range(reorder_buf.shape[0]):
            for j in range(reorder_buf.shape[1]):
                for k in range(reorder_buf.shape[3]):
                    for l in range(reorder_buf.shape[2]):
                        logT[i] = table[index,3]
                        logg[j] = table[index,4]
                        logE[k] = table[index,0]
                        mu[reorder_buf.shape[2] - l - 1] = table[index,1]
                        reorder_buf[i,j,reorder_buf.shape[2] - l - 1,k] = 10.0**(table[index,2])
                        index += 1

        buf = np.zeros(np.prod(reorder_buf.shape))

        bufdex = 0
        for i in range(reorder_buf.shape[0]):
            for j in range(reorder_buf.shape[1]):
                for k in range(reorder_buf.shape[2]):
                    for l in range(reorder_buf.shape[3]):
                        buf[bufdex] = reorder_buf[i,j,k,l]; bufdex += 1

        self._hot_atmosphere = (logT, logg, mu, logE, buf)
        '''
        )

    else:
        module += (
        '''
        raise NotImplementedError('You need to implement the setter that loads the hot atmosphere table.')
        '''
        )

if args.elsewhere_atmosphere_model is not None:
    if args.elsewhere_atmosphere_load:
        module += (
        '''
    @xpsi.Photosphere.elsewhere_atmosphere.setter
    def elsewhere_atmosphere(self, path):
        '''
        )

        if ('NSX' in args.elsewhere_atmosphere_model or 'nsx' in args.elsewhere_atmosphere_model):
            module += (
            '''
        table = np.loadtxt(path, dtype=np.double)
        logT = np.zeros(args.elsewhere_atmosphere_size[0])
        logg = np.zeros(args.elsewhere_atmosphere_size[1])
        mu = np.zeros(args.elsewhere_atmosphere_size[2])
        logE = np.zeros(args.elsewhere_atmosphere_size[3])

        reorder_buf = np.zeros((args.elsewhere_atmosphere_size[0],
                                args.elsewhere_atmosphere_size[1],
                                args.elsewhere_atmosphere_size[2],
                                args.elsewhere_atmosphere_size[3]))

        index = 0
        for i in range(reorder_buf.shape[0]):
            for j in range(reorder_buf.shape[1]):
                for k in range(reorder_buf.shape[3]):
                    for l in range(reorder_buf.shape[2]):
                        logT[i] = table[index,3]
                        logg[j] = table[index,4]
                        logE[k] = table[index,0]
                        mu[reorder_buf.shape[2] - l - 1] = table[index,1]
                        reorder_buf[i,j,reorder_buf.shape[2] - l - 1,k] = 10.0**(table[index,2])
                        index += 1

        buf = np.zeros(np.prod(reorder_buf.shape))

        bufdex = 0
        for i in range(reorder_buf.shape[0]):
            for j in range(reorder_buf.shape[1]):
                for k in range(reorder_buf.shape[2]):
                    for l in range(reorder_buf.shape[3]):
                        buf[bufdex] = reorder_buf[i,j,k,l]; bufdex += 1

        self._elsewhere_atmosphere = (logT, logg, mu, logE, buf)
            '''
            )
        else:
            module += (
            '''
        raise NotImplementedError('You need to implement the setter that loads the elsewhere atmosphere table.')
            '''
            )

module += (
'''
    @property
    def global_variables(self):
        """ For interfacing with the image-plane signal simulator.

        The extension module compiled is surface_radiation_field/archive/local_variables/PDT_U.pyx,
        which replaces the contents of surface_radiation_field/local_variables.pyx.

        """
'''
)

if len(args.hot_region_model) == 2:
    module += (
    '''
        ref_p = self.hot.objects[0]
        ref_s = self.hot.objects[1]

        return np.array([ref_p['{0}_colatitude'],
                          (ref_p['phase_shift'] + {4}) * 2.0 * math.pi,
                          ref_p['{0}_radius'],
                          ref_p['{1}_colatitude'],
                          (ref_p['phase_shift'] + {4}) * 2.0 * math.pi {3} ref_p['{2}_azimuth'],
                          ref_p['{1}_radius'],
                          ref_s['{5}_colatitude'],
                          (ref_s['phase_shift'] + {9}) * 2.0 * math.pi,
                          ref_s['{5}_radius'],
                          ref_s['{6}_colatitude'],
                          (ref_s['phase_shift'] + {9}) * 2.0 * math.pi {8} ref_p['{7}_azimuth'],
                          ref_s['{6}_radius'],
                          {10},
                          {11},
                          {12},
                          {13}])
    '''.format('super' if 'DT' in args.hot_region_model[0] else 'omit',
               'cede' if 'DT' in args.hot_region_model[0] else 'super',
               'cede' if 'DT' in args.hot_region_model[0] else 'omit',
               '+' if 'DT' in args.hot_region_model[0] else '-',
               str(0.0) if not args.is_antiphased[0] else str(0.5),
               'super' if 'DT' in args.hot_region_model[1] else 'omit',
               'cede' if 'DT' in args.hot_region_model[1] else 'super',
               'cede' if 'DT' in args.hot_region_model[1] else 'omit',
               '+' if 'DT' in args.hot_region_model[1] else '-',
               str(0.0) if not args.is_antiphased[1] else str(0.5),
               "ref_p['super_temperature']" if 'DT' in args.hot_region_model[0] else 0.0,
               "ref_p['cede_temperature']" if 'DT' in args.hot_region_model[0] else "ref_p['super_temperature']",
               "ref_s['super_temperature']" if 'DT' in args.hot_region_model[1] else 0.0,
               "ref_s['cede_temperature']" if 'DT' in args.hot_region_model[1] else "ref_s['super_temperature']")
    )
else:
    module += (
    '''
        ref = self.hot

        return np.array([ref['{0}_colatitude'],
                          (ref['phase_shift'] + {4}) * 2.0 * math.pi,
                          ref['{0}_radius'],
                          ref['{1}_colatitude'],
                          (ref['phase_shift'] + {4}) * 2.0 * math.pi {3} ref['{2}_azimuth'],
                          ref['{1}_radius'],
                          0.0,
                          0.0,
                          0.0,
                          0.0,
                          0.0,
                          0.0,
                          {5},
                          {6},
                          0.0,
                          0.0])
    '''.format('super' if 'DT' in args.hot_region_model[0] else 'omit',
               'cede' if 'DT' in args.hot_region_model[0] else 'super',
               'cede' if 'DT' in args.hot_region_model[0] else 'omit',
               '+' if 'DT' in args.hot_region_model[0] else '-',
               str(0.0) if not args.is_antiphased[0] else str(0.5),
               "ref_p['super_temperature']" if 'DT' in args.hot_region_model[0] else 0.0,
               "ref_p['cede_temperature']" if 'DT' in args.hot_region_model[0] else "ref_p['super_temperature']")
    )

write(r'{}.py'.format(os.path.join(args.module_directory_path, args.custom_photosphere_module)), module)

# Creating Prior module
module = (
r'''""" Prior module for X-PSI {0} modelling of {1} {2} event data. """

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
    Prior module for X-PSI {0} modelling of {1} {2} event data.

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
'''.format(args.model,
           _telescopes,
           args.source)
)

module += (
'''
parser.add_argument('--mass-prior',
                    type=str,
                    nargs='*',
                    action=CompileAction,
                    help='Prior inverse CDF of gravitation mass (solar masses). {1}')

parser.add_argument('--distance-prior',
                    type=str,
                    nargs='*',
                    action=CompileAction,
                    help='Prior inverse CDF of Earth distance (kpc). {1}')

parser.add_argument('--cos-inclination-prior',
                    type=str,
                    nargs='*',
                    action=CompileAction,
                    help='Prior inverse CDF of cosine of Earth inclination to stellar spin axis. {1}')

parser.add_argument('--neutral-hydrogen-column-density-prior',
                    type=str,
                    nargs='*',
                    action=CompileAction,
                    help='Prior inverse CDF of ratio of interstellar neutral hydrogen column density to the fiducial density. {1}')

parser.add_argument('--{0}-super-temperature-prior',
                    type=str,
                    nargs='*',
                    action=CompileAction,
                    help='Prior inverse CDF of hot-region {0} superseding region log10(temperature [K]). {1}')
'''.format(args.prefix[0],
           _CDF_notice)
)

if args.background_model:
    for _parameter in args.background_parameters:
        module += (
        '''
parser.add_argument('--{0}-prior',
                    type=str,
                    nargs='*',
                    action=CompileAction,
                    help='Prior inverse CDF of the ``{0}`` parameter. {1}')
        '''.format(_parameter, _CDF_notice)
)

for instrument in args.instrument:
    module += (
    '''
parser.add_argument('--{0}-energy-independent-effective-area-scaling-factor-prior',
                    type=str,
                    nargs='*',
                    action=CompileAction,
                    default=[compile('truncnorm.ppf(x, -5.0, 5.0, loc=1.0, scale=0.104)', '<string>', 'eval')],
                    help='Prior inverse CDF of the energy-independent effective area scaling factor. {1}')
    '''.format(instrument,
               _CDF_notice)
    )

def add_cede_temperature_prior(i):
    global args
    global module

    if 'DT' in args.hot_region_model[i]:
        module += (
        '''
parser.add_argument('--{0}-cede-temperature-prior',
                    type=str,
                    nargs='*',
                    action=CompileAction,
                    help='Prior inverse CDF of hot-region {0} ceding region log10(temperature [K]). {1}')
        '''.format(args.prefix[i],
                   _CDF_notice)
        )

add_cede_temperature_prior(0)

if len(args.hot_region_model) == 2 and not args.antipodal_reflection_symmetry:
    module += (
    '''
parser.add_argument('--{0}-super-temperature-prior',
                    type=str,
                    nargs='*',
                    action=CompileAction,
                    help='Prior inverse CDF of hot-region {0} superseding region log10(temperature [K]). {1}')
    '''.format(args.prefix[1],
               _CDF_notice)
    )

    add_cede_temperature_prior(1)


module += (
'''
if __name__ == '__main__':
    args, _ = parser.parse_known_args()
else:
    args, _ = parser.parse_known_args(['@{}'])
'''.format(args.config_path)
)

module += (
'''
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
'''.format(args.model,
           _telescopes,
           args.source)
)

module += (
'''
class CustomPrior(xpsi.Prior):
    """ A joint prior PDF. """

    __derived_names__ = ['compactness']

    __draws_from_support__ = 4

    def __init__(self):

        super(CustomPrior, self).__init__()
'''
)

if (   'CST' in args.hot_region_model
    or 'CDT' in args.hot_region_model
    or 'EST' in args.hot_region_model
    or 'EDT' in args.hot_region_model ):
    module += (
    '''
        self.a_f = 0.001
        self.b_f = 1.0
        self.a_zeta = 0.001
        self.b_zeta = math.pi/2.0 - self.a_zeta

        vals = np.linspace(0.0, self.b_zeta, 1000)
        self._interpolator_super_smaller = Akima1DInterpolator(self._vector_super_smaller_radius_mass(vals), vals)
        self._interpolator_super_smaller.extrapolate = True
    '''
    )

if 'PST' in args.hot_region_model or 'PDT' in args.hot_region_model:
    module += (
    '''
        self.c_f = 0.001
        self.d_f = 2.0
        self.a_xi = 0.001
        self.b_xi = math.pi/2.0 - self.a_xi

        vals = np.linspace(0.0, self.b_xi, 1000)
        self._interpolator = Akima1DInterpolator(self._vector_super_radius_mass(vals), vals)
        self._interpolator.extrapolate = True
    '''
    )

module += (
'''
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
'''
)

if 'NSX' in args.hot_atmosphere_model or 'nsx' in args.hot_atmosphere_model:
    module += (
    '''
        for g in grav:
            if not 13.7 <= g <= 15.0: # check that these NSX effective gravity table limits are correct
                return -np.inf
    '''
    )

module += (
'''
        ref = self.parameters # redefine shortcut
'''
)

if len(args.hot_region_model) == 2:
    if args.hot_region_model[0] == args.hot_region_model[1] and not args.antipodal_reflection_symmetry:
        module += (
        '''
        # hot regions can exchange to yield the exact same configuration
        # break degeneracy:
        if ref['{0}__{2}'] > ref['{1}__{2}']:
            return -np.inf
        '''.format(args.prefix[0],
                   args.prefix[1],
                   args.break_hot_region_exchange_degeneracy_with)
        )

if len(args.hot_region_model) == 2 and not args.antipodal_reflection_symmetry:

    _ST_idx = None
    if args.hot_region_model[0] == 'ST':
        _ST_idx = 0
    elif args.hot_region_model[1] == 'ST':
        _ST_idx = 1

    module += (
    '''
        # require that hot regions do not overlap within the prior support'''
    )

    if _ST_idx is not None:
        if args.hot_region_model[1 - _ST_idx] == 'ST':
            module += (
            '''
        if self._overlap(ref, '{0}', '{1}', 'super', 'super', {2:.1f}, {3:.1f}):
            return -np.inf
            '''.format(args.prefix[_ST_idx],
                       args.prefix[1 - _ST_idx],
                       0.5 if args.is_antiphased[_ST_idx] else 0.0,
                       0.5 if args.is_antiphased[1 - _ST_idx] else 0.0)
            )
        elif args.hot_region_model[1 - _ST_idx] in ['CDT', 'EDT']:
            module += (
            '''
        if self._overlap(ref, '{0}', '{1}', 'super', 'cede', {2:.1f}, {3:.1f}):
            return -np.inf
            '''.format(args.prefix[_ST_idx],
                       args.prefix[1 - _ST_idx],
                       0.5 if args.is_antiphased[_ST_idx] else 0.0,
                       0.5 if args.is_antiphased[1 - _ST_idx] else 0.0)
            )
        elif args.hot_region_model[1 - _ST_idx] == 'PDT':
            module += (
            '''
        if self._overlap(ref, '{0}', '{1}', 'super', 'super', {2:.1f}, {3:.1f}):
            return -np.inf

        if self._overlap(ref, '{0}', '{1}', 'super', 'cede', {2:.1f}, {3:.1f}):
            return -np.inf
            '''.format(args.prefix[_ST_idx],
                       args.prefix[1 - _ST_idx],
                       0.5 if args.is_antiphased[_ST_idx] else 0.0,
                       0.5 if args.is_antiphased[1 - _ST_idx] else 0.0)
            )
        elif args.hot_region_model[1 - _ST_idx] in ['CST', 'EST']:
            module += (
            '''
        if self._overlap(ref, '{0}', '{1}', 'super', 'cede', {2:.1f}, {3:.1f}):
            if not self._overlap(ref, '{0}', '{1}', 'super', 'omit', {2:.1f}, {3:.1f}, superset='{1}'):
                return -np.inf
            '''.format(args.prefix[_ST_idx],
                       args.prefix[1 - _ST_idx],
                       0.5 if args.is_antiphased[_ST_idx] else 0.0,
                       0.5 if args.is_antiphased[1 - _ST_idx] else 0.0)
            )
        elif args.hot_region_model[1 - _ST_idx] == 'PST':
            module += (
            '''
        if self._overlap_outside_mask(ref, '{0}', '{1}', 'super', {2:.1f}, {3:.1f}):
            return -np.inf
            '''.format(args.prefix[_ST_idx],
                       args.prefix[1 - _ST_idx],
                       0.5 if args.is_antiphased[_ST_idx] else 0.0,
                       0.5 if args.is_antiphased[1 - _ST_idx] else 0.0)
            )

    _DT_idx = None
    if _ST_idx is None:
        if args.hot_region_model[0] in ['CDT','EDT']:
            _DT_idx = 0
        elif args.hot_region_model[1] in ['CDT','EDT']:
            _DT_idx = 1

    if _ST_idx is None and _DT_idx is not None:
        if args.hot_region_model[1 - _DT_idx] in ['CDT', 'EDT']:
            module += (
            '''
        if self._overlap(ref, '{0}', '{1}', 'cede', 'cede', {2:.1f}, {3:.1f}):
            return -np.inf
            '''.format(args.prefix[_DT_idx],
                       args.prefix[1 - _DT_idx],
                       0.5 if args.is_antiphased[_DT_idx] else 0.0,
                       0.5 if args.is_antiphased[1 - _DT_idx] else 0.0)
            )
        elif args.hot_region_model[1 - _DT_idx] == 'PDT':
            module += (
            '''
        if self._overlap(ref, '{0}', '{1}', 'cede', 'super', {2:.1f}, {3:.1f}):
            return -np.inf

        if self._overlap(ref, '{0}', '{1}', 'cede', 'cede', {2:.1f}, {3:.1f}):
            return -np.inf
            '''.format(args.prefix[_DT_idx],
                       args.prefix[1 - _DT_idx],
                       0.5 if args.is_antiphased[_DT_idx] else 0.0,
                       0.5 if args.is_antiphased[1 - _DT_idx] else 0.0)
            )
        elif args.hot_region_model[1 - _DT_idx] in ['CST', 'EST']:
            module += (
            '''
        if self._overlap(ref, '{0}', '{1}', 'cede', 'super', {2:.1f}, {3:.1f}):
            if not self._overlap(ref, '{0}', '{1}', 'cede', 'omit', {2:.1f}, {3:.1f}, superset='{1}'):
                return -np.inf
            '''.format(args.prefix[_DT_idx],
                       args.prefix[1 - _DT_idx],
                       0.5 if args.is_antiphased[_DT_idx] else 0.0,
                       0.5 if args.is_antiphased[1 - _DT_idx] else 0.0)
            )
        elif args.hot_region_model[1 - _DT_idx] == 'PST':
            module += (
            '''
        if self._overlap_outside_mask(ref, '{0}', '{1}', 'cede', {2:.1f}, {3:.1f}):
            return -np.inf
            '''.format(args.prefix[_DT_idx],
                       args.prefix[1 - _DT_idx],
                       0.5 if args.is_antiphased[_DT_idx] else 0.0,
                       0.5 if args.is_antiphased[1 - _DT_idx] else 0.0)
            )

    _PDT_idx = None
    if _ST_idx is None and _DT_idx is None:
        if args.hot_region_model[0] == 'PDT':
            _PDT_idx = 0
        elif args.hot_region_model[1] == 'PDT':
            _PDT_idx = 1

    if _ST_idx is None and _DT_idx is None and _PDT_idx is not None:
        if args.hot_region_model[1 - _PDT_idx] == 'PDT':
            module += (
            '''
        if self._overlap(ref, '{0}', '{1}', 'super', 'super', {2:.1f}, {3:.1f}):
            return -np.inf

        if self._overlap(ref, '{0}', '{1}', 'super', 'cede', {2:.1f}, {3:.1f}):
            return -np.inf

        if self._overlap(ref, '{0}', '{1}', 'cede', 'super', {2:.1f}, {3:.1f}):
            return -np.inf

        if self._overlap(ref, '{0}', '{1}', 'cede', 'cede', {2:.1f}, {3:.1f}):
            return -np.inf
            '''.format(args.prefix[_PDT_idx],
                       args.prefix[1 - _PDT_idx],
                       0.5 if args.is_antiphased[_PDT_idx] else 0.0,
                       0.5 if args.is_antiphased[1 - _PDT_idx] else 0.0)
            )
        elif args.hot_region_model[1 - _PDT_idx] in ['CST', 'EST']:
            module += (
            '''
        if self._overlap(ref, '{0}', '{1}', 'super', 'super', {2:.1f}, {3:.1f}):
            if not self._overlap(ref, '{0}', '{1}', 'super', 'omit', {2:.1f}, {3:.1f}, superset='{1}'):
                return -np.inf

        if self._overlap(ref, '{0}', '{1}', 'cede', 'super', {2:.1f}, {3:.1f}):
            if not self._overlap(ref, '{0}', '{1}', 'cede', 'omit', {2:.1f}, {3:.1f}, superset='{1}'):
                return -np.inf
            '''.format(args.prefix[_PDT_idx],
                       args.prefix[1 - _PDT_idx],
                       0.5 if args.is_antiphased[_PDT_idx] else 0.0,
                       0.5 if args.is_antiphased[1 - _PDT_idx] else 0.0)
            )
        elif args.hot_region_model[1 - _PDT_idx] == 'PST':
            module += (
            '''
        if self._overlap_outside_mask(ref, '{0}', '{1}', 'super', {2:.1f}, {3:.1f}):
            return -np.inf
        elif self._overlap_outside_mask(ref, '{0}', '{1}', 'cede', {2:.1f}, {3:.1f}):
            return -np.inf
            '''.format(args.prefix[_PDT_idx],
                       args.prefix[1 - _PDT_idx],
                       0.5 if args.is_antiphased[_PDT_idx] else 0.0,
                       0.5 if args.is_antiphased[1 - _PDT_idx] else 0.0)
            )

    _OST_idx = None
    if _ST_idx is None and _DT_idx is None and _PDT_idx is None:
        if args.hot_region_model[0] in ['CST','EST']:
            _OST_idx = 0
        elif args.hot_region_model[1] in ['CST','EST']:
            _OST_idx = 1

    if _ST_idx is None and _DT_idx is None and _PDT_idx is None and _OST_idx is not None:
        if args.hot_region_model[1 - _OST_idx] in ['CST', 'EST']:
            module += (
            '''
        if self._overlap(ref, '{0}', '{1}', 'cede', 'cede', {2:.1f}, {3:.1f}):
            if not self._overlap(ref, '{0}', '{1}', 'cede', 'omit', {2:.1f}, {3:.1f}, superset='{1}'):
                if not self._overlap(ref, '{0}', '{1}', 'omit', 'cede', {2:.1f}, {3:.1f}, superset='{0}'):
                    return -np.inf
            '''.format(args.prefix[_OST_idx],
                       args.prefix[1 - _OST_idx],
                       0.5 if args.is_antiphased[_OST_idx] else 0.0,
                       0.5 if args.is_antiphased[1 - _OST_idx] else 0.0)
            )
        elif args.hot_region_model[1 - _OST_idx] == 'PST':
            module += (
            '''
        # incomplete prior support:
        # some configurations should be included but are not with the conditions below
        if self._overlap_outside_mask(ref, '{0}', '{1}', 'super', {2:.1f}, {3:.1f}):
            if self._overlap_outside_mask(ref, '{1}', '{0}', 'super', {3:.1f}, {2:.1f}):
                return -np.inf
            '''.format(args.prefix[_OST_idx],
                       args.prefix[1 - _OST_idx],
                       0.5 if args.is_antiphased[_OST_idx] else 0.0,
                       0.5 if args.is_antiphased[1 - _OST_idx] else 0.0)
            )

    if _ST_idx is None and _DT_idx is None and _PDT_idx is None and _OST_idx is None:
        module += (
        '''
        # incomplete prior support:
        # some configurations should be included but are not with the conditions below
        if self._overlap_outside_mask(ref, '{0}', '{1}', 'super', {2:.1f}, {3:.1f}):
            if self._overlap_outside_mask(ref, '{1}', '{0}', 'super', {3:.1f}, {2:.1f}):
                return -np.inf
        '''.format(args.prefix[0],
                   args.prefix[1],
                   0.5 if args.is_antiphased[0] else 0.0,
                   0.5 if args.is_antiphased[1] else 0.0)
        )

module += (
'''
        return 0.0
'''
)

if len(args.hot_region_model) > 1:
    module += (
    '''
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
    '''
    )

if 'PST' in args.hot_region_model and len(args.hot_region_model) > 1:

    module += (
    '''
    def _overlap_outside_mask(self, ref, x, y, x_member, x_antiphase, y_antiphase):
        """ Determine if two spherical circles overlap outisde of a masking spherical circle. """
        # check if x and PST-super overlap
        if not self._overlap(ref, x, y, x_member, 'super', x_antiphase, y_antiphase): # then no overlap
            return False

        # check if PST-super is entirely engulfed by x
        if self._overlap(ref, x, y, x_member, 'super', x_antiphase, y_antiphase,
                             superset=x, use_cached=True): # then overlap
            return True

        phi = ( ref[x + '__phase_shift'] + x_antiphase - ( ref[y + '__phase_shift'] + y_antiphase ) ) * _2pi
        # find the x and PST-super region centre coordinates in rotated coordinate frame
        x__ang_sep, x__azi = eval_coords_under_rotation(-1.0*self._colatitude(ref, y, 'omit'),
                                                        self._colatitude(ref, x, x_member),
                                                        phi)
        y__ang_sep, y__super_azi = eval_coords_under_rotation(-1.0*self._colatitude(ref, y, 'omit'),
                                                              ref[y + '__super_colatitude'],
                                                              -1.0*self._azimuth(ref, y, 'omit'))

        # check if x overlaps with mask
        if x__ang_sep > ref[x + '__super_radius'] + ref[y + '__omit_radius']: # then overlap
            return True

        # check if x is entirely engulfed by mask
        if x__ang_sep + ref[x + '__super_radius'] < ref[y + '__omit_radius']: # then no overlap
            return False

        # check if mask is entirely within PST-super
        # meaning ring-like topology and thus x must overlap with PST-super
        if y__ang_sep + ref[y + '__omit_radius'] < ref[y + '__super_radius']: # then overlap
            return True
        elif x__ang_sep + ref[y + '__omit_radius'] < ref[x + '__super_radius']: # then overlap
            return True

        # finding the terminal half angles in rotated coordinate frame
        x__phi_term = self._phi_calculator(ref[y + '__omit_radius'],
                                           ref[x + '__super_radius'],
                                           x__ang_sep)

        y__phi_term = self._phi_calculator(ref[y + '__omit_radius'],
                                           ref[y + '__super_radius'],
                                           y__ang_sep)

        # find the terminal intervals
        x__term_int = np.array([self._make_periodic(x__azi - x__phi_term), self._make_periodic(x__azi + x__phi_term)])
        y__term_int = np.array([self._make_periodic(y__super_azi - y__phi_term), self._make_periodic(y__super_azi + y__phi_term)])

        # find the widest interval first, and use it as the reference interval
        # the widest interval cannot by definition be within the narrowest interval
        # so if there is an overlap, at least one of the narrowest interval bounds
        # are within the widest interval
        if x__term_int[0] > x__term_int[1]:
            _x__interval_width = x__term_int[1] + (2.0*np.pi - x__term_int[0])
        else:
            _x__interval_width = x__term_int[1] - x__term_int[0]

        if y__term_int[0] > y__term_int[1]:
            _y__interval_width = y__term_int[1] + (2.0*np.pi - y__term_int[0])
        else:
            _y__interval_width = y__term_int[1] - y__term_int[0]

        _widest_interval = x__term_int if _x__interval_width > _y__interval_width else y__term_int
        _narrowest_interval = x__term_int if _y__interval_width < _y__interval_width else y__term_int

        # check if terminal azimuth intervals overlap
        if _widest_interval[0] > _widest_interval[1]:
            if _narrowest_interval[0] < _widest_interval[1] or _narrowest_interval[0] > _widest_interval[0]\\
            or _narrowest_interval[1] < _widest_interval[1] or _narrowest_interval[1] > _widest_interval[0]: # then overlap
                return True
        else:
            if _widest_interval[0] < _narrowest_interval[0] < _widest_interval[1]\\
            or _widest_interval[0] < _narrowest_interval[1] < _widest_interval[1]: # then overlap
                return True

        # check if both regions have a solution to maximum azimuthal extent w.r.t mask centre
        # also need to check if the x region or PST-super region contains the antipode of the mask
        # centre, because in this case there is no solution either
        if x__ang_sep + ref[x + '__super_radius'] > np.pi and y__ang_sep + ref[y + '__super_radius'] > np.pi:
            return True

        if x__ang_sep < ref[x + '__super_radius']:
            x__sep_max_azi = 0.0
        elif x__ang_sep + ref[x + '__super_radius'] < np.pi:
            # find colatitude of the point subtended by the maximum azimuthal points w.r.t mask centre
            x__sep_max_azi = np.arccos(np.cos(x__ang_sep) / np.cos(ref[x + '__super_radius']))
        else:
            x__sep_max_azi = np.pi

        if y__ang_sep < ref[y + '__super_radius']:
            y__sep_max_azi = 0.0
        elif y__ang_sep + ref[y + '__super_radius'] < np.pi:
            y__sep_max_azi = np.arccos(np.cos(y__ang_sep) / np.cos(ref[y + '__super_radius']))
        else:
            y__sep_max_azi = np.pi

        # check if maximum azimuthal points for both x and PST-super regions are greater than mask radius
        if (x__sep_max_azi > ref[y + '__omit_radius']) and (y__sep_max_azi > ref[y + '__omit_radius']): # then overlap
            return True

        # find the angle between the lines joining the terminal points to the region centres
        x__term_ang = self._phi_calculator(ref[y + '__omit_radius'],
                                           x__ang_sep,
                                           ref[x + '__super_radius'])

        y__term_ang = self._phi_calculator(ref[y + '__omit_radius'],
                                           y__ang_sep,
                                           ref[y + '__super_radius'])

        _tmp = np.abs(x__term_ang - math.pi/2) > np.abs(y__term_ang - math.pi/2)

        # check if only the x maximum azimuthal points are outside the mask region
        if x__sep_max_azi > ref[y + '__omit_radius'] and _tmp: # then overlap
            return True

        # check if only the PST-super maximum azimuthal points are outside the mask region
        if y__sep_max_azi > ref[y + '__omit_radius'] and not _tmp: # then overlap
            return True

        # else, the overlap between the regions is a subset of the mask region, and permitted
        return False

    @staticmethod
    def _phi_calculator(psi, zeta, v):
        cos_phi = (np.cos(zeta) -  np.cos(v) * np.cos(psi)) / (np.sin(v) * np.sin(psi))
        return np.arccos(cos_phi)

    @staticmethod
    def _make_periodic(angle):
        if angle < 0:
            angle += _2pi
        elif angle > _2pi:
            angle -= _2pi
        return angle
    '''
    )

if (   'CST' in args.hot_region_model
    or 'CDT' in args.hot_region_model
    or 'EST' in args.hot_region_model
    or 'EDT' in args.hot_region_model ):
    module += (
    '''
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
    '''
    )

if 'PST' in args.hot_region_model or 'PDT' in args.hot_region_model:
    module += (
    '''
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
    '''
    )

module += (
'''
    def inverse_sample(self, hypercube=None):
        """ Draw sample uniformly from the distribution via inverse sampling. """

        global args

        to_cache = self.parameters.vector

        if hypercube is None:
            hypercube = np.random.rand(len(self))

        # the base method is useful, so to avoid writing that code again:
        _ = super(CustomPrior, self).inverse_sample(hypercube)

        ref = parameters = self.parameters # redefine shortcut
'''
)

module += (
'''
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
'''
)

def construct_hot_region_prior(i):

    global args
    global module

    module += (
    '''
        # inverse sample parameters of hot-region {0}
        idx = ref.index('{0}__super_colatitude')
        a, b = ref.get_param('{0}__super_colatitude').bounds
        a = math.cos(a); b = math.cos(b)
        ref['{0}__super_colatitude'] = math.acos(b + (a - b) * hypercube[idx])
    '''.format(args.prefix[i])
    )

    if args.hot_region_model[i] in ['CST','CDT','EST','EDT']:
        module += (
        '''
        # radius of superseding region (omit or super code object)
        idx = ref.index('{0}__{1}_radius')
        ref['{0}__{1}_radius'] = float(self._interpolator_super_smaller(hypercube[idx]))

        # radius of ceding region (super or cede code object)
        idx = ref.index('{0}__{2}_radius')
        ref['{0}__{2}_radius'] = self._inverse_sample_cede_larger_radius(hypercube[idx], ref['{0}__{1}_radius'])
        '''.format(args.prefix[i],
                   'omit'  if 'ST' in args.hot_region_model[i] else 'super',
                   'super' if 'ST' in args.hot_region_model[i] else 'cede')
        )
    elif args.hot_region_model[i] in ['EST','EDT']:
        module += (
        '''
        # coordinates of mask or ceding region (omit or cede code object)
        idx = ref.index('{0}__{3}_colatitude')
        ref[idx] = hypercube[idx] * (ref['{0}__{2}_radius'] - ref['{0}__{1}_radius'])

        ref[idx], ref['{0}__{3}_azimuth'] = eval_coords_under_rotation(ref['{0}__super_colatitude'],
                                                                       ref['{0}__{3}_colatitude'],
                                                                       ref['{0}__{3}_azimuth'])
        '''.format(args.prefix[i],
                   'omit'  if 'ST' in args.hot_region_model[i] else 'super',
                   'super' if 'ST' in args.hot_region_model[i] else 'cede',
                   'omit'  if 'ST' in args.hot_region_model[i] else 'cede')
        )
    elif args.hot_region_model[i] in ['PST','PDT']:
        module += (
        '''
        # radius of superseding region (omit or super code object)
        idx = ref.index('{0}__{1}_radius')
        ref['{0}__{1}_radius'] = float(self._interpolator(hypercube[idx]))

        # radius of ceding region (super or cede code object)
        idx = ref.index('{0}__{2}_radius')
        ref['{0}__{2}_radius'] = self._inverse_sample_cede_radius(hypercube[idx], ref['{0}__{1}_radius'])

        # coordinates of mask or ceding region (omit or cede code object)
        idx = ref.index('{0}__{3}_colatitude')
        if ref['{0}__{1}_radius'] <= ref['{0}__{2}_radius']:
            ref[idx] = hypercube[idx] * (ref['{0}__{2}_radius'] + ref['{0}__{1}_radius'])
        else:
            ref[idx] = (  ref['{0}__{1}_radius']
                        - ref['{0}__{2}_radius']
                        + 2.0*hypercube[idx]*ref['{0}__{2}_radius'] )

        ref[idx], ref['{0}__{3}_azimuth'] = eval_coords_under_rotation(ref['{0}__super_colatitude'],
                                                                       ref['{0}__{3}_colatitude'],
                                                                       ref['{0}__{3}_azimuth'])
        '''.format(args.prefix[i],
                   'omit'  if 'ST' in args.hot_region_model[i] else 'super',
                   'super' if 'ST' in args.hot_region_model[i] else 'cede',
                   'omit'  if 'ST' in args.hot_region_model[i] else 'cede')
        )

construct_hot_region_prior(0)

if len(args.hot_region_model) == 2 and not args.antipodal_reflection_symmetry:
    construct_hot_region_prior(1)

module += (
'''
        # restore proper cache
        for parameter, cache in list(zip(self.parameters, to_cache)):
            parameter.cached = cache

        # it is important that we return the desired vector because it is
        # automatically written to disk by MultiNest and only by MultiNest
        return self.parameters.vector
'''
)

module += (
'''
    def transform(self, p, **kwargs):
        """ A transformation for post-processing. """

        p = list(p) # copy

        # used ordered names and values
        ref = dict(list(zip(self.parameters.names, p)))

        # compactness ratio M/R_eq
        p += [gravradius(ref['mass']) / ref['radius']]

        return p
'''
)

write(r'{}.py'.format(os.path.join(args.module_directory_path, args.custom_prior_module)), module)

# Creating Interstellar module
module = (
'''""" Interstellar module for X-PSI {0} modelling of {1} {2} event data. """

import numpy as np
import math

import xpsi
from xpsi import Parameter

from scipy.interpolate import Akima1DInterpolator

class CustomInterstellar(xpsi.Interstellar):
    """ Apply interstellar attenuation model {3}. """

    def __init__(self, energies, attenuation, bounds, values = None):

        if values is None: values = {{}}

        assert len(energies) == len(attenuation), 'Array length mismatch.'

        self._lkp_energies = energies # for lookup
        self._lkp_attenuation = attenuation # for lookup

        N_H = Parameter('neutral_hydrogen_column_density',
                        strict_bounds = (0.0, 50.0),
                        bounds = bounds.get('neutral_hydrogen_column_density', None),
                        doc = 'Neutral hydrogen column density in units of the fiducial column density',
                        symbol = r'$N_{{\\rm H}}$',
                        value = values.get('neutral_hydrogen_column_density', None),
                        permit_prepend = False)

        self._interpolator = Akima1DInterpolator(self._lkp_energies,
                                                 self._lkp_attenuation)
        self._interpolator.extrapolate = True

        super(CustomInterstellar, self).__init__(N_H)

    def attenuation(self, energies):
        """ Interpolate the attenuation coefficients.

        Useful for post-processing.

        """
        return self._interpolate(energies)**(self['neutral_hydrogen_column_density'])

    def _interpolate(self, energies):
        """ Helper. """
        _att = self._interpolator(energies)
        _att[_att < 0.0] = 0.0
        return _att

    @classmethod
    def load(cls, path,
             energy_column=0,
             attenuation_column=2,
             **kwargs):
        """ Load attenuation file. """

        # check the loading assumptions and comment out the exception throw if they are true
        #raise NotImplementedError('Implement the class method to load the interstellar attenuation table.')

        # template

        temp = np.loadtxt(path, dtype=np.double)

        energies = temp[:,energy_column]

        attenuation = temp[:,attenuation_column]

        return cls(energies, attenuation, **kwargs)
'''.format(args.model,
           _telescopes,
           args.source,
           args.attenuation_model)
)

write(r'{}.py'.format(os.path.join(args.module_directory_path, args.custom_interstellar_module)), module)

_instruments = args.instrument[0]
for _x in args.instrument[1:-1]:
    _instruments += ', {}'.format(_x)
_instruments += ', and {}'.format(args.instrument[-1])

# Creating Instrument module
module = (
'''""" Instrument module for X-PSI {0} modelling of {1} {2} event data. """

import numpy as np
import math

import xpsi

from xpsi import Parameter
from xpsi.utils import make_verbose

class CustomInstrument(xpsi.Instrument):
    """ {3}. """
    def construct_matrix(self):
        """ Implement response matrix parameterisation. """
        matrix = self['energy_independent_effective_area_scaling_factor'] * self.matrix
        matrix[matrix < 0.0] = 0.0

        return matrix

    def __call__(self, signal, *args):
        """ Overwrite. """

        matrix = self.construct_matrix()

        self._cached_signal = np.dot(matrix, signal)

        return self._cached_signal
'''.format(args.model,
           _telescopes,
           args.source,
           _instruments)
)

for instrument in args.instrument:
    module += (
    '''
    @classmethod
    @make_verbose('Loading {0} response matrix',
                  'Response matrix loaded')
    def {0}(cls,
              bounds,
              values,
              ARF,
              RMF,
              channel_energies,
              max_input,
              max_channel,
              min_input=0,
              min_channel=0,
              effective_area_scaling_factor=1.0,
              ARF_skiprows=0,
              ARF_low_column=1,
              ARF_high_column=2,
              ARF_area_column=3,
              RMF_skiprows=0,
              RMF_usecol=-1,
              channel_energies_skiprows=0,
              channel_energies_low_column=0,
              **kwargs):
        """ Load {0} instrument response matrix. """

        alpha = Parameter('{1}',
                          strict_bounds = (0.1,1.9),
                          bounds = bounds.get('{1}', None),
                          doc='{0} energy-independent effective area scaling factor',
                          symbol = r'$\\alpha_{{\\rm {0}}}$',
                          value = values.get('{1}',
                                             1.0 if bounds.get('{1}', None) is None else None))

        # check the loading assumptions and comment out the exception throw if they are true
        #raise NotImplementedError('Implement the class method for loading the {0} instrument.')

        #working template for the github example:
        #if min_input != 0:
        #    min_input = int(min_input)
        #max_input = int(max_input)
        #ARF = np.loadtxt(ARF, dtype=np.double, skiprows=3)
        #RMF = np.loadtxt(RMF, dtype=np.double)
        #channel_energies = np.loadtxt(channel_energies, dtype=np.double, skiprows=3)[:,1:]
        #matrix = np.ascontiguousarray(RMF[min_input:max_input,20:201].T, dtype=np.double)
        #edges = np.zeros(ARF[min_input:max_input,3].shape[0]+1, dtype=np.double)
        #edges[0] = ARF[min_input,1]; edges[1:] = ARF[min_input:max_input,2]
        #for i in range(matrix.shape[0]):
        #    matrix[i,:] *= ARF[min_input:max_input,3]
        #channels = np.arange(20, 201)
        #return cls(matrix, edges, channels, channel_energies[20:202,-2], alpha, **kwargs)


        # another template
        #ARF = np.loadtxt(ARF, dtype=np.double, skiprows=ARF_skiprows)
        #RMF = np.loadtxt(RMF, dtype=np.double, skiprows=RMF_skiprows, usecols=RMF_usecol)
        #channel_energies = np.loadtxt(channel_energies, dtype=np.double, skiprows=channel_energies_skiprows)

        #matrix = np.zeros((channel_energies.shape[0], ARF.shape[0]))

        #for i in range(ARF.shape[0]):
        #    matrix[:,i] = RMF[i*channel_energies.shape[0]:(i+1)*channel_energies.shape[0]]

        #max_input = int(max_input)
        #if min_input != 0:
        #    min_input = int(min_input)

        #edges = np.zeros(max_input - min_input + 1, dtype=np.double)

        #edges[0] = ARF[min_input, ARF_low_column]; edges[1:] = ARF[min_input:max_input, ARF_high_column]

        #RSP = np.zeros((max_channel - min_channel,
        #                max_input - min_input), dtype=np.double)

        #for i in range(RSP.shape[0]):
        #    RSP[i,:] = matrix[i+min_channel, min_input:max_input] * ARF[min_input:max_input, ARF_area_column] * effective_area_scaling_factor

        #channels = np.arange(min_channel, max_channel)

        _RMF_path = RMF
        ARF = np.loadtxt(ARF, dtype=np.double, skiprows=ARF_skiprows)
        RMF = np.loadtxt(_RMF_path, dtype=np.double, skiprows=RMF_skiprows, usecols=RMF_usecol)
        RMF_zerocol = np.loadtxt(_RMF_path, dtype=np.double, skiprows=RMF_skiprows, usecols=0)
        channel_energies = np.loadtxt(channel_energies, dtype=np.double, skiprows=channel_energies_skiprows)

        matrix = np.zeros((channel_energies.shape[0], ARF.shape[0]))

        last = 0
        k = 0
        counter = 0
        for i in range(RMF_zerocol.shape[0]):
            if math.floor(RMF_zerocol[i]) == RMF_zerocol[i] and RMF_zerocol[i] != 0.0:
                counter += 1
                if i == 0: continue
                else:
                    for j in range(i - last):
                        matrix[channel_energies.shape[0] - i + last + j, k] = RMF[last + j] #* ARF[k, ARF_area_column]
                    #if i - last != channel_energies.shape[0]:
                    #    print('arf i=%i'%RMF_zerocol[i], 'i=%i'%i, 'last=%i'%last, 'nchans=%i'%(i - last))
                    last = i
                    k += 1

        max_input = int(max_input)
        if min_input != 0:
            min_input = int(min_input)

        edges = np.zeros(max_input - min_input + 1, dtype=np.double)

        edges[0] = ARF[min_input, ARF_low_column]; edges[1:] = ARF[min_input:max_input, ARF_high_column]

        RSP = np.zeros((max_channel - min_channel,
                        max_input - min_input), dtype=np.double)

        for i in range(RSP.shape[0]):
            RSP[i,:] = matrix[i+min_channel, min_input:max_input] * ARF[min_input:max_input, ARF_area_column] * effective_area_scaling_factor

        channels = np.arange(min_channel, max_channel)

        return cls(RSP,
                   edges,
                   channels,
                   channel_energies[min_channel:max_channel+1,channel_energies_low_column],
                   alpha, **kwargs)
    '''.format(instrument,
               'energy_independent_effective_area_scaling_factor')
    )

write(r'{}.py'.format(os.path.join(args.module_directory_path, args.custom_instrument_module)), module)

if not args.background_model:
    sys.exit(0)

# Creating Background module
module = (
'''""" Background module for X-PSI {0} modelling of {1} {2} event data. """

import numpy as np
import math

import xpsi

from xpsi import Parameter
'''.format(args.model,
           _telescopes,
           args.source)
)

def _write_background_class(instrument=None):

    global module

    module += (
    '''
class {0}CustomBackground(xpsi.Background):
    """ Model for incident background photon specific flux. """

    def __init__(self, bounds=None, value=None):
        # check the model implementation and comment out the exception
        # throw if appropriate
        raise NotImplementedError('Implement the background model.')

        # template is a powerlaw model

        doc = """
        Powerlaw spectral index.
        """
        index = xpsi.Parameter('powerlaw_index',
                               strict_bounds = (-5.0, -1.01),
                               bounds = bounds.get('powerlaw_index', None),
                               doc = doc,
                               symbol = r'$\\Gamma{2}$',
                               value = values.get('powerlaw_index', None),
                               permit_prepend = {1})

        doc = """
        Powerlaw spectral normalization (photons/keV/s/cm^2 @ 1 keV).
        """
        index = xpsi.Parameter('powerlaw_normalization',
                               strict_bounds = (0.0, None),
                               bounds = bounds.get('powerlaw_normalization', None),
                               doc = doc,
                               symbol = r'$N{2}$',
                               value = values.get('powerlaw_normalization', None),
                               permit_prepend = {1})

        super(CustomBackground, self).__init__(index)

    def __call__(self, energy_edges, phases):
        """ Evaluate the incident background field. """

        # check the model implementation and comment out the exception
        # throw if appropriate
        raise NotImplementedError('Implement the background model.')

        # template is a powerlaw model
        G = self['powerlaw_index']

        temp = np.zeros((energy_edges.shape[0] - 1, phases.shape[0]))

        temp[:,0] = (energy_edges[1:]**(G + 1.0) - energy_edges[:-1]**(G + 1.0)) / (G + 1.0)

        for i in range(phases.shape[0]):
            temp[:,i] = temp[:,0]

        self.background = temp * self['powerlaw_normalization']
    '''.format(instrument + '_' if instrument is not None else '',
               'False' if args.background_shared_class else 'True',
               '_{{\\rm {0}}}'.format(instrument) if instrument is not None else '')
)

if args.background_shared_class:
    _write_background_class(None)
else:
    for _instrument in args.instrument:
        _write_background_class(_instrument)

write(r'{}.py'.format(os.path.join(args.module_directory_path, args.custom_background_module)), module)
