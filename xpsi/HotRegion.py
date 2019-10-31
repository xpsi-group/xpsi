from __future__ import division, print_function

from .global_imports import *
from . import global_imports

from .cellmesh.mesh_tools import allocate_cells as _allocate_cells
from .cellmesh.mesh import construct_spot_cellMesh as _construct_spot_cellMesh
from .cellmesh.polar_mesh import construct_polar_cellMesh as _construct_polar_cellMesh
from .cellmesh.rays import compute_rays as _compute_rays

from .ParameterSubspace import ParameterSubspace, BoundsError

class RayError(xpsiError):
    """ Raised if a numerical problems encountered during ray integration. """

class PulseError(xpsiError):
    """ Raised if a numerical problems encountered during integration. """

class HotRegion(ParameterSubspace):
    """ A photospheric hot region that is contiguous.

    Instances of this class can take the following forms:

        * a circular, simply-connected region;
        * a ceding region with a non-radiating superseding region;
        * a ceding region with a radiating superseding region.

    For the first option, set ``hole=False`` and ``cede=False``.
    For the third option, set ``hole=False`` and ``cede=True``.

    For the second option, set ``hole=True`` and ``cede=False``. This is
    counter-intuitive. These settings produce a superseding region with a
    non-radiating hole, which is equivalent to a ceding region with a
    non-radiating superseding region. The *hole* may not strictly always be
    a hole, because it behaves as a superseding region that does not radiate.
    A superseding region can exist partially outside its relative ceding
    region, so the non-superseded subset of the ceding region is
    simply-connected (i.e., does not have a hole).

    .. note:: Setting ``hole=True`` and ``cede=True`` ignores the hole.

    The ``concentric`` setting defines whether the superseding region is
    concentric with the ceding region or a hole, supposing either
    ``hole=True`` or ``cede=True``.

    """

    def __init__(self,
                 num_params,
                 bounds,
                 symmetry = True,
                 hole = False,
                 cede = False,
                 concentric = True,
                 sqrt_num_cells = 32,
                 min_sqrt_num_cells = 10,
                 max_sqrt_num_cells = 80,
                 num_rays = 200,
                 num_leaves = 64,
                 num_phases = None,
                 phases = None,
                 do_fast = False,
                 fast_sqrt_num_cells = 16,
                 fast_min_sqrt_num_cells = 4,
                 fast_max_sqrt_num_cells = 16,
                 fast_num_rays = 100,
                 fast_num_leaves = 32,
                 is_secondary = False):
        """
        :param num_params: Number of parameters for hot region model. For a
                           circular hot-spot, this includes the
                           the spot centre colatitude and spot angular radius.

        :param list bounds: Hard parameter bounds for the instance of
                            :class:`.ParameterSubspace.ParameterSubspace`.

        .. note::

            The order of the list of the parameter bounds must follow an order
            of precedence:

            * super colatitude
            * super angular radius
            * cede or hole colatitude
            * cede or hole angular radius
            * cede or hole relative azimuth
            * parameters controlling the local comoving radiation field over
              the photospheric 2-surface (entirely the user's responsibility)

            These bounds might not be actually used, depending the user's
            implementation of the joint prior, and the user can in that case
            specify ``[None,None]`` for bounds pertaining to the ceding region
            or the hole.

        :param bool symmetry:
            Is the radiation field axisymmetric (w.r.t the stellar rotation
            axis) within a hot region with only a superseding region?
            This determines which ``CellMesh`` integrator to
            deploy. Only set to ``False`` if you overwrite
            :meth:`_eval_srcRadFieldParamVectors` in a custom subclass.

        :param sqrt_num_cells:
            Number of cells in both colatitude and azimuth which form a
            regular mesh on a curved 2-surface (a spacelike leaf of the
            spacetime foliation). This is square-root of the approximate
            number of cells whose centres should lie within a hot region.

        :param min_sqrt_num_cells:
            Sets the minimum number of cells per *subregion*, discretised to
            the same degree in colatitude and azimuth. This setting has an
            effect only when the hot region has two temperature components, in
            the form of a superseding region and a ceding region.

        :param min_sqrt_num_cells:
            Sets the maximum number of cells per *subregion*, discretised to
            the same degree in colatitude and azimuth. This setting has an
            effect even when there is only one temperature component.

        :param num_rays: Number of rays to trace (integrate) at each colatitude,
                         distributed in angle subtended between ray tangent
                         4-vector and radial outward unit vector w.r.t a local
                         orthonormal tetrad.

        :param int num_leaves: Number of leaves mesh motion is discretised into.

        :param int num_phases: Number of phases in a discrete representation of
                               the specific flux pulses on the interval
                               :math:`[0,1]`.

        :param phases: If not ``None``, a :class:`numpy.ndarray` of phases
                       for a discrete representation of the specific flux
                       pulses on the interval :math:`[0,1]`. If ``None``
                       (default), the :obj:`num_phases` argument is utilised.

        :param bool do_fast:
            Activate fast precomputation to guide cell distribution between
            two radiating regions at distinct temperatures.

        .. note:: For each of the resolution parameters listed above, there is
                  a corresponding parameter whose value is respected if the
                  fast precomputation mode is activated.

        :param bool is_secondary:
            If ``True``, shifts the cell mesh by :math:`\pi` radians about
            the stellar rotation axis for pulse integration.

        """
        super(HotRegion, self).__init__(num_params, bounds)

        if self._num_params < 3:
            raise BoundsError('A hot region requires at least three '
                              'sets of parameter bounds.')

        for colatitude in self._bounds[0]:
            if not 0.0 < colatitude < _pi:
                raise BoundsError('Invalid superseding region colatitude bound.')

        for radius in self._bounds[1]:
            if not 0.0 < radius <= _pi/2.0:
                raise BoundsError('Invalid superseding region angular radius bound.')

        self._is_secondary = is_secondary

        self.do_fast = do_fast

        self.set_num_rays(num_rays, fast_num_rays)
        self.set_num_cells(sqrt_num_cells,
                           min_sqrt_num_cells, max_sqrt_num_cells,
                           fast_sqrt_num_cells,
                           fast_min_sqrt_num_cells, fast_max_sqrt_num_cells)
        self.set_phases(num_leaves, fast_num_leaves, num_phases, phases)

        if symmetry:
            from .cellmesh.integrator_for_azimuthal_invariance import integrate_radField as _integrator
            self._hole = hole
            self._concentric = concentric
        elif cede:
            self._cede = cede
            self._concentric = concentric
            from .cellmesh.integrator_for_azimuthal_invariance import integrate_radField as _integrator
        else:
            from .cellmesh.integrator import integrate_radField as _integrator
        self._integrator = _integrator

    def set_num_rays(self, num_rays, fast_num_rays):
        self.num_rays = num_rays
        self._fast_num_rays = fast_num_rays

    @property
    def num_rays(self):
        """ Get the number of rays integrated per parallel. """
        return self._num_rays

    @num_rays.setter
    def num_rays(self, n):
        """ Set the number of rays integrated per parallel. """
        try:
            self._num_rays = int(n)
        except TypeError:
            raise TypeError('Number of rays must be an integer.')

    def set_num_cells(self, sqrt_num_cells,
                      min_sqrt_num_cells, max_sqrt_num_cells,
                      fast_sqrt_num_cells,
                      fast_min_sqrt_num_cells, fast_max_sqrt_num_cells):

        self.sqrt_num_cells = sqrt_num_cells
        self._num_cells = int(self.sqrt_num_cells**2)
        self._min_sqrt_num_cells = int(min_sqrt_num_cells)
        self._max_sqrt_num_cells = int(max_sqrt_num_cells)

        assert self._min_sqrt_num_cells <= self._max_sqrt_num_cells,\
               ('Minimum number of cells must be less than or equal to '
                'maximum number of cells.')

        self._fast_num_cells = int(fast_sqrt_num_cells**2)
        self._fast_min_sqrt_num_cells = int(fast_min_sqrt_num_cells)
        self._fast_max_sqrt_num_cells = int(fast_max_sqrt_num_cells)

        assert self._fast_min_sqrt_num_cells <= self._fast_max_sqrt_num_cells,\
               ('Minimum number of cells must be less than or equal to '
                'maximum number of cells.')

    @property
    def sqrt_num_cells(self):
        """ Get the number of cell parallels. """
        return self._sqrt_num_cells

    @sqrt_num_cells.setter
    def sqrt_num_cells(self, n):
        """ Set the number of cell parallels. """
        try:
            self._sqrt_num_cells = int(n)
        except TypeError:
            raise TypeError('Number of cells must be an integer.')

    @property
    def leaves(self):
        """ Get the leaves of the photospheric foliation. """
        return self._leaves

    @property
    def phases(self):
        """ Get the leaves of the photospheric foliation. """
        return self._phases

    def set_phases(self, num_leaves, fast_num_leaves=None,
                   num_phases=None, phases=None):
        """ Construct the set of interpolation phases and foliation leaves.

        :param int num_leaves: Number of leaves mesh motion is discretised into.

        :param int num_phases: Number of phases in a discrete representation of
                               the specific flux pulses on the interval
                               :math:`[0,1]`. If ``None``, the number of phases
                               interpolated at is set equal to the number of
                               leaves.

        :param phases: If not ``None``, a :class:`numpy.ndarray` of phases
                       for a discrete representation of the specific flux
                       pulses on the interval :math:`[0,1]`. If ``None``
                       (default), the :obj:`num_phases` argument is utilised.
        """
        if num_phases is None:
            num_phases = num_leaves
            fast_num_phases = fast_num_leaves

        if phases is None:
            try:
                self._phases_cycles = _np.linspace(0.0, 1.0, int(num_phases))
                self._fast_phases_cycles = _np.linspace(0.0, 1.0, int(fast_num_phases))
            except TypeError:
                raise TypeError('Number of phases must be an integer.')
        else:
            try:
                assert isinstance(phases, _np.ndarray)
                assert phases.ndim == 1
                assert (phases >= 0.0).all() & (phases <= 1.0).all()
                for i in range(phases.shape[0] - 1):
                    assert phases[i] < phases[i+1]
            except AssertionError:
                raise TypeError('Phases must be a one-dimensional '
                                '``numpy.ndarray`` with monotonically '
                                'increasing elements on the interval [0,1].')
            else:
                self._phases_cycles = phases

        self._phases = _2pi * self._phases_cycles
        self._fast_phases = _2pi * self._fast_phases_cycles

        try:
            self._leaves =  _np.linspace(0.0, _2pi, int(num_leaves))
            self._fast_leaves =  _np.linspace(0.0, _2pi, int(fast_num_leaves))
        except TypeError:
            raise TypeError('Number of leaves must be an integer.')

    @property
    def phases_in_cycles(self):
        """ Get the phases (in cycles) the pulse is interpolated at. """
        return self._phases_cycles

    @property
    def fast_phases_in_cycles(self):
        return self._fast_phases_cycles

    @property
    def num_cells(self):
        """ Get the total number of cells in the mesh. """
        return self._num_cells

    def print_settings(self):
        """ Print numerical settings. """
        print('Base number of cell parallels: ', self.sqrt_num_cells)
        print('Number of rays per parallel: ', self.num_rays)
        print('Number of photospheric leaves: ', self.leaves.shape[0])
        print('Number of interpolation phases: ', self.phases.shape[0])

    def __construct_cellMesh(self, st, superColatitude, superRadius,
                             cedeColatitude, cedeRadius, cedeAzimuth,
                             holeColatitude, holeRadius, holeAzimuth,
                             fast_total_counts, threads):
        """ Call a low-level routine to construct a mesh representation.

        :param st: Instance of :class:`~.Spacetime.Spacetime`.

        :param float superColatitude: Superseding region colatitude (radians).

        :param float cedeRadius: Ceding region angular radius (radians).

        :param float superRadius: Superseding region angular radius (radians).

        :param float cedeOffset: Ceding region angular offset (radians).

        :param int threads: Number of ``OpenMP`` threads for mesh construction.

        """
        num_cells = self._fast_num_cells if self.fast_mode \
                                                    else self._num_cells
        min_sqrt_num_cells = self._fast_min_sqrt_num_cells if self.fast_mode \
                                                else self._min_sqrt_num_cells
        max_sqrt_num_cells = self._fast_max_sqrt_num_cells if self.fast_mode \
                                                else self._max_sqrt_num_cells

        (self._super_sqrt_num_cells,
         self._super_num_cells,
         self._cede_sqrt_num_cells,
         self._cede_num_cells) = _allocate_cells(num_cells,
                                                      min_sqrt_num_cells,
                                                      max_sqrt_num_cells,
                                                      st.M,
                                                      st.r_s,
                                                      st.R,
                                                      st.Omega,
                                                      st.zeta,
                                                      st.epsilon,
                                                      superRadius,
                                                      cedeRadius,
                                                      superColatitude,
                                                      cedeColatitude,
                                                      -cedeAzimuth,
                                                      holeColatitude,
                                                      holeRadius,
                                                      holeAzimuth,
                                                      fast_total_counts)

        #print('----------------')
        #print('Cell allocation:')
        #print('super region: %i' % self._super_sqrt_num_cells)
        #print('cede region: %i' % self._cede_sqrt_num_cells)
        #print('----------------')

        if superColatitude - superRadius < 0.0 or superColatitude + superRadius > _pi:
            mesh_func = _construct_polar_cellMesh
        else:
            mesh_func = _construct_spot_cellMesh

        (self._super_theta,
         self._super_phi,
         self._super_r,
         self._super_cellArea,
         self._super_maxAlpha,
         self._super_cos_gamma,
         self._super_effGrav) = mesh_func(threads,
                                          self._super_num_cells,
                                          self._super_sqrt_num_cells,
                                          st.M,
                                          st.r_s,
                                          st.R,
                                          st.zeta,
                                          st.epsilon,
                                          superRadius,
                                          superColatitude,
                                          holeRadius,
                                          holeColatitude,
                                          holeAzimuth)

        self._super_phi -= holeAzimuth

        if self._is_secondary:
            self._super_phi += _pi

        if cedeRadius > 0.0:
            if cedeColatitude - cedeRadius < 0.0 or cedeColatitude + cedeRadius > _pi:
                mesh_func = _construct_polar_cellMesh
            else:
                mesh_func = _construct_spot_cellMesh

            (self._cede_theta,
             self._cede_phi,
             self._cede_r,
             self._cede_cellArea,
             self._cede_maxAlpha,
             self._cede_cos_gamma,
             self._cede_effGrav) = mesh_func(threads,
                                             self._cede_num_cells,
                                             self._cede_sqrt_num_cells,
                                             st.M,
                                             st.r_s,
                                             st.R,
                                             st.zeta,
                                             st.epsilon,
                                             cedeRadius,
                                             cedeColatitude,
                                             superRadius,
                                             superColatitude,
                                             -cedeAzimuth)

            self._cede_phi += cedeAzimuth

            if self._is_secondary:
                self._cede_phi += _pi

    def __calibrate_lag(self, st):
        """ Calibrate lag for cell mesh and normalise by oscillation period. """

        R_i = st.R_r_s

        C = (1.0 / self._super_r_s_over_r - R_i)
        C += _log((1.0 / self._super_r_s_over_r - 1.0) / (R_i - 1.0))
        C *= st.r_s / _c
        for j in range(self._super_lag.shape[1]):
            self._super_lag[:,j] -= C
        self._super_lag *= st.Omega

        try:
            C = (1.0 / self._cede_r_s_over_r - R_i)
        except AttributeError:
            pass
        else:
            C += _log((1.0 / self._cede_r_s_over_r - 1.0) / (R_i - 1.0))
            C *= st.r_s / _c
            for j in range(self._cede_lag.shape[1]):
                self._cede_lag[:,j] -= C
            self._cede_lag *= st.Omega

    def __compute_rays(self, st, threads):
        """ Integrate rays for null mapping.

        These rays represent a null mapping from photosphere to a point at
        some effective infinity.

        :param st: Instance of :class:`~.Spacetime.Spacetime`.

        :param int threads: Number of ``OpenMP`` threads for ray integration.

        """
        self._super_r_s_over_r = _contig(st.r_s / self._super_r, dtype = _np.double)

        num_rays = self._fast_num_rays if self.fast_mode else self._num_rays

        (terminate_calculation,
         self._super_deflection,
         self._super_cos_alpha,
         self._super_lag,
         self._super_maxDeflection) = _compute_rays(threads,
                                                    self._super_sqrt_num_cells,
                                                    st.r_s,
                                                    self._super_r_s_over_r,
                                                    self._super_maxAlpha,
                                                    num_rays)

        if terminate_calculation == 1:
            raise RayError('Fatal numerical problem during ray integration.')

        try:
            self._cede_r_s_over_r = _contig(st.r_s / self._cede_r, dtype = _np.double)
        except AttributeError:
            pass
        else:
            (terminate_calculation,
             self._cede_deflection,
             self._cede_cos_alpha,
             self._cede_lag,
             self._cede_maxDeflection) = _compute_rays(threads,
                                                       self._cede_sqrt_num_cells,
                                                       st.r_s,
                                                       self._cede_r_s_over_r,
                                                       self._cede_maxAlpha,
                                                       num_rays)

            if terminate_calculation == 1:
                raise RayError('Fatal numerical problem during ray integration.')

        self.__calibrate_lag(st)

    def __compute_cellParamVecs(self, p):
        """
        Precomputes photospheric source radiation field parameter vectors
        cell-by-cell.

        Designed here simply for uniform effective temperature superseding
        and ceding regions.

        :params p: Arguments for :meth:`_eval_srcRadFieldParamVectors`.

        """
        self._super_radiates = _np.greater(self._super_cellArea, 0.0).astype(_np.int32)
        self._super_cellParamVecs = _np.ones((self._super_radiates.shape[0],
                                              self._super_radiates.shape[1],
                                              2), dtype=_np.double)
        self._super_cellParamVecs[...,:-1] *= _np.array(p[0])
        for i in range(self._super_cellParamVecs.shape[1]):
            self._super_cellParamVecs[:,i,-1] *= self._super_effGrav

        try:
            self._cede_radiates = _np.greater(self._cede_cellArea, 0.0).astype(_np.int32)
        except AttributeError:
            pass
        else:
            self._cede_cellParamVecs = _np.ones((self._cede_radiates.shape[0],
                                                 self._cede_radiates.shape[1],
                                                 2), dtype=_np.double)
            self._cede_cellParamVecs[...,:-1] *= _np.array(p[1])
            for i in range(self._cede_cellParamVecs.shape[1]):
                self._cede_cellParamVecs[:,i,-1] *= self._cede_effGrav

    @staticmethod
    def _psi(theta, phi, colatitude):
        """ Coordinates of cell centres in rotated spherical coordinate system.

        Transformation is anticlockwise rotation about y-axis of Cartesian
        basis.

        """
        return _arccos(_cos(colatitude)
                       * _cos(theta)
                       + _sin(colatitude)
                       * _sin(theta)
                       * _cos(phi))

    def embed(self, spacetime, p, fast_total_counts, threads, *args):
        """ Embed the hot region.

        :param bool correction: Correct the integral over the radiation field
                                *elsewhere* by accounting for the time-dependent
                                component arising from the presence of the
                                hot region.

        :param cellParamVecs: A :class:`numpy.ndarray` of ``float``\ s. If a
                              :obj:`correction` is to be made, this array
                              contains a parameter vector for each cell of the
                              hot region mesh. The parameter vectors should be
                              identical and equivalent to the parameter vector
                              used for the radiation field *elsewhere*\ . If no
                              correction is to be made, the array needs to be
                              passed but can be a zero array.

        """
        if self.fast_mode and not self.do_fast:
            return None

        try:
            self._cede
        except AttributeError:
            if not self._hole:
                self.__construct_cellMesh(spacetime,
                                          p[0],
                                          p[1],
                                          p[0],
                                          0.0,
                                          0.0,
                                          0.0,
                                          0.0,
                                          0.0,
                                          fast_total_counts,
                                          threads)
            else:
                if self._concentric:
                    self.__construct_cellMesh(spacetime,
                                              p[0],
                                              p[1],
                                              p[0],
                                              0.0,
                                              0.0,
                                              p[0],
                                              p[2],
                                              0.0,
                                              fast_total_counts,
                                              threads)
                else:
                    self.__construct_cellMesh(spacetime,
                                              p[0],
                                              p[1],
                                              p[0],
                                              0.0,
                                              0.0,
                                              p[2],
                                              p[3],
                                              p[4],
                                              fast_total_counts,
                                              threads)
        else:
            if self._concentric:
                self.__construct_cellMesh(spacetime,
                                          p[0], # super colatitude
                                          p[1], # super radius
                                          p[0], # cede colatitude
                                          p[2], # cede radius
                                          0.0, # cede relative azimuth
                                          p[0], # super hole colatitude
                                          0.0, # super hole radius
                                          0.0, # super hole azimuth
                                          fast_total_counts,
                                          threads)
            else:
                self.__construct_cellMesh(spacetime,
                                          p[0],
                                          p[1],
                                          p[2],
                                          p[3],
                                          p[4],
                                          p[0],
                                          0.0,
                                          0.0,
                                          fast_total_counts,
                                          threads)

        self.__compute_rays(spacetime, threads)

        try:
            if self._cede:
                self.__compute_cellParamVecs(p[3:5])
        except AttributeError:
            if self._hole:
                if self._concentric:
                    self.__compute_cellParamVecs(p[3:4])
                else:
                    self.__compute_cellParamVecs(p[5:6])
            else:
                self.__compute_cellParamVecs(p[2:3])
        else:
            if not self._cede:
                self.__compute_cellParamVecs(p[5:7])

        if args:
            self._super_correctionVecs = args[0](self._super_theta.shape, args[1])
            for i in range(self._super_theta.shape[1]):
                self._super_correctionVecs[:,i,-1] *= self._super_effGrav

            try:
                self._cede_correctionVecs = args[0](self._cede_theta.shape, args[1])
            except AttributeError:
                pass
            else:
                for i in range(self._cede_theta.shape[1]):
                    self._cede_correctionVecs[:,i,-1] *= self._cede_effGrav
        else:
            self._super_correctionVecs = None
            self._cede_correctionVecs = None

    @property
    def do_fast(self):
        return self._do_fast

    @do_fast.setter
    def do_fast(self, activate):
        self._do_fast = activate

    @property
    def fast_mode(self):
        return self._fast_mode

    @fast_mode.setter
    def fast_mode(self, activate):
        self._fast_mode = activate

    def integrate(self, st, energies, threads,
                  hot_atmosphere, elsewhere_atmosphere):
        """ Integrate over the photospheric radiation field.

        Calls the CellMesh integrator, with or without exploitation of
        azimuthal invariance of the radiation field of the hot region.

        :param st: Instance of :class:`~.Spacetime.Spacetime`.

        :param energies: A one-dimensional :class:`numpy.ndarray` of energies
                         in keV.

        :param int threads: Number of ``OpenMP`` threads for pulse
                            integration.

        """
        if self.fast_mode and not self.do_fast:
            try:
                if self.cede:
                    return (None, None)
            except AttributeError:
                return (None,)

        leaves = self._fast_leaves if self.fast_mode else self._leaves
        phases = self._fast_phases if self.fast_mode else self._phases
        num_rays = self._fast_num_rays if self.fast_mode else self._num_rays

        if isinstance(energies, tuple):
            try:
                super_energies, cede_energies = energies
            except ValueError:
                super_energies = energies[0]
                try:
                    self._cede_cellArea
                except AttributeError:
                    pass
                else:
                    cede_energies = super_energies
        else:
            super_energies = cede_energies = energies

        #print('Primary super energies', super_energies)
        #print('Primary cede energies', cede_energies)

        # change low-level code so as to require only a parameter vector
        # for the correction, and not a larger array of the same parameter
        # vector? Effective gravity varies with colatitude however.
        super_pulse = self._integrator(threads,
                                       st.M,
                                       st.R,
                                       st.Omega,
                                       st.r_s,
                                       st.zeta,
                                       st.epsilon,
                                       st.i,
                                       self._super_cellArea,
                                       self._super_r,
                                       self._super_r_s_over_r,
                                       self._super_theta,
                                       self._super_phi,
                                       self._super_cellParamVecs,
                                       self._super_radiates,
                                       self._super_correctionVecs,
                                       num_rays,
                                       self._super_deflection,
                                       self._super_cos_alpha,
                                       self._super_lag,
                                       self._super_maxDeflection,
                                       self._super_cos_gamma,
                                       super_energies,
                                       leaves,
                                       phases,
                                       hot_atmosphere,
                                       elsewhere_atmosphere)

        if super_pulse[0] == 1:
            raise PulseError('Fatal numerical error during superseding-'
                             'region pulse integration.')

        try:
            cede_pulse = self._integrator(threads,
                                          st.M,
                                          st.R,
                                          st.Omega,
                                          st.r_s,
                                          st.zeta,
                                          st.epsilon,
                                          st.i,
                                          self._cede_cellArea,
                                          self._cede_r,
                                          self._cede_r_s_over_r,
                                          self._cede_theta,
                                          self._cede_phi,
                                          self._cede_cellParamVecs,
                                          self._cede_radiates,
                                          self._cede_correctionVecs,
                                          num_rays,
                                          self._cede_deflection,
                                          self._cede_cos_alpha,
                                          self._cede_lag,
                                          self._cede_maxDeflection,
                                          self._cede_cos_gamma,
                                          cede_energies,
                                          leaves,
                                          phases,
                                          hot_atmosphere,
                                          elsewhere_atmosphere)
        except AttributeError:
            pass
        else:
            if cede_pulse[0] == 1:
                raise PulseError('Fatal numerical error during ceding-region '
                                 'pulse integration.')
            else:
                return (super_pulse[1], cede_pulse[1])

        return (super_pulse[1],)
