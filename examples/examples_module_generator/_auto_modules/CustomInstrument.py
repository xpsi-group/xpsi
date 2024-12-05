""" Instrument module for X-PSI CST+PDT modelling of NICER PSR J0030+0451 event data. """

import numpy as np

import xpsi

from xpsi import Parameter
from xpsi.utils import make_verbose

class CustomInstrument(xpsi.Instrument):
    """ XTI, and XTI. """
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

    @classmethod
    @make_verbose('Loading XTI response matrix',
                  'Response matrix loaded')
    def XTI(cls,
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
        """ Load XTI instrument response matrix. """

        alpha = Parameter('energy_independent_effective_area_scaling_factor',
                          strict_bounds = (0.1,1.9),
                          bounds = bounds.get('energy_independent_effective_area_scaling_factor', None),
                          doc='XTI energy-independent effective area scaling factor',
                          symbol = r'$\alpha_{\rm XTI}$',
                          value = values.get('energy_independent_effective_area_scaling_factor',
                                             1.0 if bounds.get('energy_independent_effective_area_scaling_factor', None) is None else None))

        # check the loading assumptions and comment out the exception throw if they are true
        #raise NotImplementedError('Implement the class method for loading the XTI instrument.')

        #working template for the github example:
        if min_input != 0:
            min_input = int(min_input)
        max_input = int(max_input)
        ARF = np.loadtxt(ARF, dtype=np.double, skiprows=3)
        RMF = np.loadtxt(RMF, dtype=np.double)
        channel_energies = np.loadtxt(channel_energies, dtype=np.double, skiprows=3)[:,1:]
        matrix = np.ascontiguousarray(RMF[min_input:max_input,20:201].T, dtype=np.double)
        edges = np.zeros(ARF[min_input:max_input,3].shape[0]+1, dtype=np.double)
        edges[0] = ARF[min_input,1]; edges[1:] = ARF[min_input:max_input,2]
        for i in range(matrix.shape[0]):
            matrix[i,:] *= ARF[min_input:max_input,3]
        channels = np.arange(20, 201)
        return cls(matrix, edges, channels, channel_energies[20:202,-2], alpha, **kwargs)


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

        #_RMF_path = RMF
        #ARF = np.loadtxt(ARF, dtype=np.double, skiprows=ARF_skiprows)
        #RMF = np.loadtxt(_RMF_path, dtype=np.double, skiprows=RMF_skiprows, usecols=RMF_usecol)
        #RMF_zerocol = np.loadtxt(_RMF_path, dtype=np.double, skiprows=RMF_skiprows, usecols=0)
        #channel_energies = np.loadtxt(channel_energies, dtype=np.double, skiprows=channel_energies_skiprows)

        #matrix = np.zeros((channel_energies.shape[0], ARF.shape[0]))

        #last = 0
        #k = 0
        #counter = 0
        #for i in range(RMF_zerocol.shape[0]):
        #    if math.floor(RMF_zerocol[i]) == RMF_zerocol[i] and RMF_zerocol[i] != 0.0:
        #        counter += 1
        #        if i == 0: continue
        #        else:
        #            for j in range(i - last):
        #                matrix[channel_energies.shape[0] - i + last + j, k] = RMF[last + j] #* ARF[k, ARF_area_column]
        #            #if i - last != channel_energies.shape[0]:
        #            #    print('arf i=%i'%RMF_zerocol[i], 'i=%i'%i, 'last=%i'%last, 'nchans=%i'%(i - last))
        #            last = i
        #            k += 1

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

       # return cls(RSP,
       #            edges,
       #            channels,
       #            channel_energies[min_channel:max_channel+1,channel_energies_low_column],
       #            alpha, **kwargs)
