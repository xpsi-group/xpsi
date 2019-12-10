from __future__ import division, print_function

from .global_imports import *
from . import global_imports

from .Parameter import Parameter

from .HotRegion import HotRegion, PulseError

class TwoHotRegions(HotRegion):
    """ Two photospheric hot regions, related via antipodal reflection symmetry.

    The *primary* hot region is represented by the class from which the
    :class:`.TwoHotRegions.TwoHotRegions` derives -- i.e., the
    :class:`.HotRegion.HotRegion` class.

    This class differs from the :class:`.HotRegions.HotRegions` class because
    it works *specifically* for two hot regions, via inheritance, where one
    is derived from the other purely via antipodal reflection symmetry.
    The *secondary* hot region is handled by adding behaviour to the parent
    class.

    This means that all parameters describing the secondary are
    derived from the primary: the secondary is antipodal and given by
    an equatorial reflection of the primary, followed by a rotation
    by :math:`\pi` radians about the rotation axis. The mirroring is
    thus with respect to a 2-plane through the coordinate origin
    which is perpendicular to the line through the origin and the
    primary centre.

    :param dict kwargs:
        Keyword arguments passed to :class:`.HotRegion.HotRegion` class.

    """
    def __init__(self, *args, **kwargs):
        # force primary because secondary phase  handled in subclass methods
        _ = kwargs.pop('is_secondary', None)
        super(TwoHotRegions, self).__init__(is_secondary = False,
                                            *args, **kwargs)

    @property
    def cellArea(self):
        """ Get the areas of cells in the secondary mesh. """
        try:
            return (self.__super_cellArea, self.__cede_cellArea)
        except AttributeError:
            return (self.__super_cellArea, None)

    def embed(self, spacetime, photosphere, fast_total_counts, threads, *args):
        """ Embed the hot regions. """

        if fast_total_counts is not None:
            fast_primary_total_counts = fast_total_counts[0]
            fast_secondary_total_counts = fast_total_counts[1]
        else:
            fast_primary_total_counts = None
            fast_secondary_total_counts = None

        # embed secondary first using methods of parent class, using
        # derived parameters as required
        super(TwoHotRegions, self).embed(spacetime,
                                         photosphere,
                                         fast_secondary_total_counts,
                                         threads, *args)

        # reflect in equatorial plane
        self.__super_theta = _pi - self._super_theta
        try:
            self.__cede_theta = _pi - self._cede_theta
        except AttributeError:
            pass

        # bind NumPy arrays to new name with mangle; no shallow/deep copy
        # required. The data in these arrays is required for integration of
        # secondary pulse.
        self.__super_phi = self._super_phi + _pi
        self.__super_r = self._super_r
        self.__super_cellArea = self._super_cellArea
        self.__super_cos_gamma = self._super_cos_gamma
        self.__super_r_s_over_r = self._super_r_s_over_r
        self.__super_lag = self._super_lag
        self.__super_deflection = self._super_deflection
        self.__super_cos_alpha = self._super_cos_alpha
        self.__super_maxDeflection = self._super_maxDeflection
        self.__super_radiates = self._super_radiates
        self.__super_cellParamVecs = self._super_cellParamVecs
        self.__super_correctionVecs = self._super_correctionVecs

        try:
            self.__cede_phi = self._cede_phi + _pi
        except AttributeError:
            pass
        else:
            self.__cede_r = self._cede_r
            self.__cede_cellArea = self._cede_cellArea
            self.__cede_cos_gamma = self._cede_cos_gamma
            self.__cede_r_s_over_r = self._cede_r_s_over_r
            self.__cede_lag = self._cede_lag
            self.__cede_deflection = self._cede_deflection
            self.__cede_cos_alpha = self._cede_cos_alpha
            self.__cede_maxDeflection = self._cede_maxDeflection
            self.__cede_radiates = self._cede_radiates
            self.__cede_cellParamVecs = self._cede_cellParamVecs
            self.__cede_correctionVecs = self._cede_correctionVecs

        if fast_total_counts is not None:
            # embed primary; names (attributes) in parent class rebound
            super(TwoHotRegions, self).embed(spacetime, photosphere,
                                             fast_primary_total_counts,
                                             threads, *args)

    def integrate(self, st, energies, threads,
                  hot_atmosphere, elsewhere_atmosphere):
        """ Integrate over the photospheric radiation field.

        Calls the CellMesh integrator, with or without exploitation of
        azimuthal invariance of the radiation field of the hot region.

        :param st: Instance of :class:`~.Spacetime.Spacetime`.

        :param energies:
            A one-dimensional :class:`numpy.ndarray` of energies in keV.

        :param int threads:
            Number of ``OpenMP`` threads for pulse integration.

        """
        if isinstance(energies, tuple):
            primary_energies = energies[0]
            secondary_energies = energies[1]
        else:
            primary_energies = secondary_energies = energies

        if isinstance(secondary_energies, tuple):
            try:
                super_energies, cede_energies = secondary_energies
            except ValueError:
                super_energies = secondary_energies[0]
                try:
                    self.__cede_cellArea
                except AttributeError:
                    pass
                else:
                    cede_energies = super_energies
        else:
            super_energies = cede_energies = secondary_energies

        primary = super(TwoHotRegions, self).integrate(st, primary_energies,
                                                       threads,
                                                       hot_atmosphere,
                                                       elsewhere_atmosphere)

        leaves = self._fast_leaves if self.fast_mode else self._leaves
        phases = self._fast_phases if self.fast_mode else self._phases
        num_rays = self._fast_num_rays if self.fast_mode else self._num_rays

        super_pulse = self._integrator(threads,
                                          st.M,
                                          st.R,
                                          st.Omega,
                                          st.r_s,
                                          st.zeta,
                                          st.epsilon,
                                          st.i,
                                          self.__super_cellArea,
                                          self.__super_r,
                                          self.__super_r_s_over_r,
                                          self.__super_theta,
                                          self.__super_phi,
                                          self.__super_cellParamVecs,
                                          self.__super_radiates,
                                          self.__super_correctionVecs,
                                          num_rays,
                                          self.__super_deflection,
                                          self.__super_cos_alpha,
                                          self.__super_lag,
                                          self.__super_maxDeflection,
                                          self.__super_cos_gamma,
                                          super_energies,
                                          leaves,
                                          phases,
                                          hot_atmosphere,
                                          elsewhere_atmosphere)

        if super_pulse[0] == 1:
            raise PulseError('Fatal numerical error during secondary '
                             'superseding-region pulse integration.')

        try:
            cede_pulse = self._integrator(threads,
                                             st.M,
                                             st.R,
                                             st.Omega,
                                             st.r_s,
                                             st.zeta,
                                             st.epsilon,
                                             st.i,
                                             self.__cede_cellArea,
                                             self.__cede_r,
                                             self.__cede_r_s_over_r,
                                             self.__cede_theta,
                                             self.__cede_phi,
                                             self.__cede_cellParamVecs,
                                             self.__cede_radiates,
                                             self.__cede_correctionVecs,
                                             num_rays,
                                             self.__cede_deflection,
                                             self.__cede_cos_alpha,
                                             self.__cede_lag,
                                             self.__cede_maxDeflection,
                                             self.__cede_cos_gamma,
                                             cede_energies,
                                             leaves,
                                             phases,
                                             hot_atmosphere,
                                             elsewhere_atmosphere)
        except AttributeError:
            pass
        else:
            if cede_pulse[0] == 1:
                raise PulseError('Fatal numerical error during pulse integration.')
            else:
                return (primary, (super_pulse[1], cede_pulse[1]))

        return (primary, (super_pulse[1],))
