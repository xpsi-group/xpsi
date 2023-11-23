""" Photosphere module for X-PSI CST+PDT modelling of NICER PSR J0030+0451 event data. """

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
    Photosphere module for X-PSI CST+PDT modelling of NICER PSR J0030+0451 event data.

    You should import this module.

    For help: python %(prog)s -h

    """,
    fromfile_prefix_chars='@')

parser.add_argument('--hot-atmosphere-size',
                             type=int,
                             nargs=4,
                             help='Size of each of the four dimensions of the numeric atmosphere table for the hot regions.')

if __name__ == '__main__':
    args, _ = parser.parse_known_args()
else:
    args, _ = parser.parse_known_args(['@./config.ini'])

import numpy as np
import math

import xpsi

class CustomPhotosphere(xpsi.Photosphere):

    @xpsi.Photosphere.hot_atmosphere.setter
    def hot_atmosphere(self, path):

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

    @property
    def global_variables(self):
        """ For interfacing with the image-plane signal simulator.

        The extension module compiled is surface_radiation_field/archive/local_variables/PDT_U.pyx,
        which replaces the contents of surface_radiation_field/local_variables.pyx.

        """

        ref_p = self.hot.objects[0]
        ref_s = self.hot.objects[1]

        return np.array([ref_p['omit_colatitude'],
                          (ref_p['phase_shift'] + 0.0) * 2.0 * math.pi,
                          ref_p['omit_radius'],
                          ref_p['super_colatitude'],
                          (ref_p['phase_shift'] + 0.0) * 2.0 * math.pi - ref_p['omit_azimuth'],
                          ref_p['super_radius'],
                          ref_s['super_colatitude'],
                          (ref_s['phase_shift'] + 0.0) * 2.0 * math.pi,
                          ref_s['super_radius'],
                          ref_s['cede_colatitude'],
                          (ref_s['phase_shift'] + 0.0) * 2.0 * math.pi + ref_p['cede_azimuth'],
                          ref_s['cede_radius'],
                          0.0,
                          ref_p['super_temperature'],
                          ref_s['super_temperature'],
                          ref_s['cede_temperature']])
