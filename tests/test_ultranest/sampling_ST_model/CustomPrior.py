

import numpy as np
import math
from scipy.stats import truncnorm

import xpsi
from xpsi.global_imports import _G, _csq, _km, _2pi, gravradius, _dpr
from xpsi import Parameter

from scipy.interpolate import Akima1DInterpolator

from scipy.stats import truncnorm

class CustomPrior(xpsi.Prior):
    """ A custom (joint) prior distribution.

    Source: Fictitious
    Model variant: ST-

    """

    # __derived_names__ = ['compactness', 'phase_separation',]

    def __init__(self):
        """ Nothing to be done.

        A direct reference to the spacetime object could be put here
        for use in __call__:

        .. code-block::

            self.spacetime = ref

        Instead we get a reference to the spacetime object through the
        a reference to a likelihood object which encapsulates a
        reference to the spacetime object.

        """
        super(CustomPrior, self).__init__() # not strictly required if no hyperparameters

    def __call__(self, p = None):
        """ Evaluate distribution at ``p``.

        :param list p: Model parameter values.

        :returns: Logarithm of the distribution evaluated at ``p``.

        """
        temp = super(CustomPrior, self).__call__(p)
        if not np.isfinite(temp):
            return temp

        # based on contemporary EOS theory
        # if not self.parameters['radius'] <= 16.0:
        #     return -np.inf

        # ref = self.parameters.star.spacetime # shortcut

        # # Compactness limit
        # R_p = 1.0 + ref.epsilon * (-0.788 + 1.030 * ref.zeta)
        # if R_p < 1.76 / ref.R_r_s:
        #     return -np.inf

        # mu = math.sqrt(-1.0 / (3.0 * ref.epsilon * (-0.788 + 1.030 * ref.zeta)))

        # 2-surface cross-section have a single maximum in |z|
        # i.e., an elliptical surface; minor effect on support, if any,
        # for high spin frequenies
        # if mu < 1.0:
        #     return -np.inf

        # ref = self.parameters

        return 0.0

    def inverse_sample(self, hypercube=None):
        """ Draw sample uniformly from the distribution via inverse sampling. """

        to_cache = self.parameters.vector

        if hypercube is None:
            hypercube = np.random.rand(len(self))

        # the base method is useful, so to avoid writing that code again:
        _ = super(CustomPrior, self).inverse_sample(hypercube)

        ref = self.parameters # shortcut

        # idx = ref.index('distance')
        # ref['distance'] = truncnorm.ppf(hypercube[idx], -5.0, 5.0, loc=1.0, scale=0.1)

        # flat priors in cosine of hot region centre colatitudes (isotropy)
        # support modified by no-overlap rejection condition
        # idx = ref.index('hot__super_colatitude')
        # a, b = ref.get_param('hot__super_colatitude').bounds
        # a = math.cos(a); b = math.cos(b)
        # ref['hot__super_colatitude'] = math.acos(b + (a - b) * hypercube[idx])


        # restore proper cache
        for parameter, cache in zip(ref, to_cache):
            parameter.cached = cache

        # it is important that we return the desired vector because it is
        # automatically written to disk by MultiNest and only by MultiNest
        return self.parameters.vector

    def transform(self, p, **kwargs):
        """ A transformation for post-processing. """

        p = list(p) # copy

        # used ordered names and values
        ref = dict(zip(self.parameters.names, p))

        # compactness ratio M/R_eq
        # p += [gravradius(ref['mass']) / ref['radius']]

        return p
