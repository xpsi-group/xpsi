from xpsi.global_imports import *
from xpsi import _verbose
from xpsi.utils import make_verbose
from xpsi.Likelihood import Likelihood
from xpsi.Prior import Prior

import h5py

try:
    import pymultinest
except ImportError:
    print('Check your PyMultiNest installation.')
    raise
else:
    if _verbose:
        print('Imported PyMultiNest.')

class NestedSampler(object):
    """ Extended MultiNest wrapper.

    :param int ndims: Number of parameters.

    :param likelihood: An instance of :class:`~.Likelihood.Likelihood`.

    :param prior: An instance of :class:`~.Prior.Prior`.

    :param kwargs: A dictionary of MultiNest runtime settings, to be
                   handled via the PyMultiNest wrapper.

    """
    def __init__(self,
                 ndims,
                 likelihood,
                 prior):

        try:
            self._ndims = int(ndims)
        except TypeError:
            raise TypeError('Dimensionality must be an integer.')

        if not isinstance(likelihood, Likelihood):
            raise TypeError('Invalid type for likelihood object.')
        else:
            self._likelihood = likelihood

        if not isinstance(prior, Prior):
            raise TypeError('Invalid type for prior object.')
        else:
            self._prior = prior

    @make_verbose('Commencing integration', 'Integration completed')
    def __call__(self, LHS_seed=None, **kwargs):
        """ Integrate.

        :param kwargs: Keyword arguments passed to :func:`pymultinest.solve`.

        """

        kwargs.setdefault('sampling_efficiency', 0.8)
        kwargs['sampling_efficiency'] /= self._prior.unit_hypercube_frac(LHS_seed=LHS_seed)

        yield 'Sampling efficiency set to: %.4f.'%kwargs['sampling_efficiency']

        # ignore the pymultinest output object
        # it has never been used in the author's workflows but change this
        # it is useful to you
        _ = pymultinest.solve(self._likelihood,
                              self._prior.inverse_sample,
                              self._ndims,
                              log_zero = self._likelihood.llzero,
                              **kwargs)
        
    def write_results_as_hdf5(self, output_basename: str):
        """Save the output files of MultiNest as .hdf5 file, which is faster for 
        post-processing.

        :param output_basename: Keyword argument outputfiles_basename originally 
            passed to :func:`pymultinest.solve`.
        """
        output_ext = ['', 'dead-birth', 'phys_live-birth']

        paramnames = self._likelihood.names

        with h5py.File(f'{output_basename}.hdf5', 'w') as f:
            
            for ext in output_ext:
                # load original multinest output file 
                original_output = _np.loadtxt(f"{output_basename}{ext}.txt")

                # set meta depending on output file 
                if ext == '':
                    metadata = ['weights', '-2logl'] + paramnames
                    ext = 'post_unweighted' # rename for hdf5 file
                if ext == 'dead-birth':
                    metadata = paramnames + ['logl', 'logl_birth', 
                                            'log_prior_mass', 'node_no']
                if ext == 'phys_live-birth':
                    metadata = paramnames + ['logl', 'logl_birth', 
                                            'node_no']
                
                # create a group for data, keep original order of datafile 
                data_group = f.create_group(ext, track_order=True)

                # add datasets, one for each column
                for idx, name in enumerate(metadata):
                    if original_output.ndim == 1:
                        data_group.create_dataset(name, data=original_output[idx])
                    else:
                        data_group.create_dataset(name, data=original_output[:,idx])