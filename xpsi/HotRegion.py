from __future__ import division, print_function

from .global_imports import *
from . import global_imports

from .cellmesh.mesh_tools import allocate_cells as _allocate_cells
from .cellmesh.mesh import construct_spot_cellMesh as _construct_spot_cellMesh
from .cellmesh.polar_mesh import construct_polar_cellMesh as _construct_polar_cellMesh
from .cellmesh.rays import compute_rays as _compute_rays

from .Parameter import Parameter, Derive
from .ParameterSubspace import ParameterSubspace

class RayError(xpsiError):
    """ Raised if a problem was encountered during ray integration. """

class PulseError(xpsiError):
    """ Raised if a problem was encountered during signal integration. """

class HotRegion(ParameterSubspace):
    """ A photospheric hot region that is contiguously radiating.

    Instances of this class can take the following forms:

        * a circular, simply-connected region;
        * a ceding region with a non-radiating superseding region;
        * a ceding region with a radiating superseding region.

    For the first option, set ``omit=False`` and ``cede=False``.
    For the third option, set ``omit=False`` and ``cede=True``.

    For the second option, set ``omit=True`` and ``cede=False``. Yes, this is
    counter-intuitive. These settings produce a superseding region with a
    non-radiating omission region, which is *equivalent* to a ceding region with
    a non-radiating superseding region. The *omission* region may not strictly
    always be a hole, because it behaves as a superseding region that does
    not radiate. A superseding region can exist partially outside its relative
    ceding region, so the non-superseded subset of the ceding region is
    then simply-connected (i.e., does not have a hole).

    .. note::

        Setting ``omit=True`` and ``cede=True`` ignores the omission setting.

    The ``concentric`` setting defines whether the superseding region is
    concentric with the ceding region or a omission region, supposing either
    ``omit=True`` or ``cede=True``.

    These helper settings simply set up the optional parameter definitions
    so you do not have to do so more manually by providing keys-value pairs
    in the ``bounds`` and ``values`` dictionaries (see below).

    :param dict bounds:
        Hard prior parameter bounds for the free parameters.

    :param bool symmetry:
        Is the radiation field axisymmetric (w.r.t the stellar rotation
        axis) within superseding (and ceding) member regions? This determines
        which ``CellMesh`` integrator to deploy based on this safety check.

    :param sqrt_num_cells:
        Number of cells in both colatitude and azimuth which form a
        regular mesh on a curved 2-surface (a spacelike leaf of the
        spacetime foliation). This is the square-root of the approximate
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
        the stellar rotation axis for pulse integration. This is merely
        a choice that can be made, and is not crucial. Note that the
        (fast) phase-shifting applied near the end of the likelihood
        evaluation is related to this choice and thus phase-shift parameter
        prior support can be chosen accordingly.

    .. note::

        The parameters are as follows:

            + super colatitude
            + super angular radius
            + cede *or* omit colatitude
            + cede *or* omit angular radius
            + cede *or* omit relative azimuth
            + parameters controlling the local comoving radiation field over
              the photospheric 2-surface (entirely the user's responsibility
              unless defaulting with ``custom is None``.

        If the ceding region or the omission region are concentric with the
        superseding region, the colatitude and relative azimuth of the
        ceding region or omission region are not parameters.

        These bounds might not be actually used, depending the user's
        implementation of the joint prior, and the user can in that case
        specify ``(None,None)`` for bounds pertaining to the ceding region
        or the omission region. These bounds will be set to strict bounds,
        and the user can implement more complicated prior support boundaries
        in a :class:`~.Prior.Prior` subclass instance.

    :param iterable custom:
        Iterable over :class:`~.Parameter.Parameter` instances. If you
        supply custom parameter definitions, you need to overwrite the
        :func:`~.HotRegion.HotRegion.__compute_cellParamVecs` method to
        handle your custom behaviour. This custom behaviour can only
        target the radiation field within the hot regions as defined
        above.

    """
    required_names = ['super_colatitude',
                      'super_radius',
                      'phase_shift',
                      'super_temperature (if no custom specification)']

    optional_names = ['omit_colatitude',
                      'omit_radius',
                      'omit_azimuth',
                      'cede_colatitude',
                      'cede_radius',
                      'cede_azimuth',
                      'cede_temperature']

    def __init__(self,
                 bounds,
                 values,
                 symmetry = True,
                 omit = False,
                 cede = False,
                 concentric = False,
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
                 fast_num_phases = None,
                 fast_phases = None,
                 is_secondary = False,
                 custom = None,
                 **kwargs):

        self.is_secondary = is_secondary

        self.do_fast = do_fast

        self.set_num_rays(num_rays, fast_num_rays)

        self.set_num_cells(sqrt_num_cells,
                           min_sqrt_num_cells, max_sqrt_num_cells,
                           fast_sqrt_num_cells,
                           fast_min_sqrt_num_cells, fast_max_sqrt_num_cells)

        self.set_phases(num_leaves, num_phases, phases,
                        fast_num_leaves, fast_num_phases, fast_phases)

        # find the required integrator
        if symmetry: # can we safely assume azimuthal invariance?
            from .cellmesh.integrator_for_azimuthal_invariance import integrate_radField as _integrator
        else: # more general purpose
            from .cellmesh.integrator import integrate_radField as _integrator
        self._integrator = _integrator

        # first the parameters that are fundemental to this class
        doc = """
        The colatitude of the centre of the superseding region [radians].
        """
        super_colat = Parameter('super_colatitude',
                                strict_bounds = (0.0, _pi),
                                bounds = bounds.get('super_colatitude', None),
                                doc = doc,
                                symbol = r'$\theta$',
                                value = values.get('super_colatitude', None))

        doc = """
        The angular radius of the (circular) superseding region [radians].
        """
        super_radius = Parameter('super_radius',
                                 strict_bounds = (0.0, _pi/2.0),
                                 bounds = bounds.get('super_radius', None),
                                 doc = doc,
                                 symbol = r'$\psi$',
                                 value = values.get('super_radius', None))

        doc = """
        The phase of the hot region, a periodic parameter [cycles].
        """
        phase_shift = Parameter('phase_shift',
                                strict_bounds = (-0.5, 0.5),
                                bounds = bounds.get('phase_shift', None),
                                doc = doc,
                                symbol = r'$\phi$',
                                value = values.get('phase_shift', None))

        self.cede = cede # takes precedence below over omission setting
        self.omit = omit
        self.concentric = concentric

        # helper callable
        class BindMe(Derive):
            def __init__(self):
                pass

            def __call__(self, boundto, caller):
                return caller['super_colatitude']

        bindme = BindMe() # to parameter instances

        if not self.cede: # means just super region, possibly with omission
            bounds['cede_radius'] = None
            values['cede_radius'] = 0.0

            bounds['cede_colatitude'] = None
            values['cede_colatitude'] = bindme

            bounds['cede_azimuth'] = None
            values['cede_azimuth'] = 0.0

            if not self.omit: # else take free settings from input
                bounds['omit_radius'] = None
                values['omit_radius'] = 0.0

            if self.concentric or not self.omit:
                bounds['omit_colatitude'] = None
                values['omit_colatitude'] = bindme

                bounds['omit_azimuth'] = None
                values['omit_azimuth'] = 0.0
        else: # means super + cede regions, no omission
            if self.concentric:
                bounds['cede_colatitude'] = None
                values['cede_colatitude'] = bindme

                bounds['cede_azimuth'] = None
                values['cede_azimuth'] = 0.0

            bounds['omit_radius'] = None
            values['omit_radius'] = 0.0

            bounds['omit_colatitude'] = None
            values['omit_colatitude'] = bindme

            bounds['omit_azimuth'] = None
            values['omit_azimuth'] = 0.0

        doc = """
        The colatitude of the centre of the omission region [radians].
        """
        omit_colat = Parameter('omit_colatitude',
                               strict_bounds = (0.0, _pi),
                               bounds = bounds.get('omit_colatitude', None),
                               doc = doc,
                               symbol = r'$\upsilon$',
                               value = values.get('omit_colatitude', None))

        doc = """
        The angular radius of the (circular) omission region [radians].
        """
        omit_radius = Parameter('omit_radius',
                                strict_bounds = (0.0, _pi/2.0),
                                bounds = bounds.get('omit_radius', None),
                                doc = doc,
                                symbol = r'$\Delta$',
                                value = values.get('omit_radius', None))

        doc = """
        The azimuth of the centre of the omission region relative to the
        centre of the superseding region [radians].
        """
        omit_azi = Parameter('omit_azimuth',
                             strict_bounds = (-_pi, _pi),
                             bounds = bounds.get('omit_azimuth', None),
                             doc = doc,
                             symbol = r'$\Phi$',
                             value = values.get('omit_azimuth', None))

        doc = """
        The colatitude of the centre of the ceding region [radians].
        """
        cede_colat = Parameter('cede_colatitude',
                               strict_bounds = (0.0, _pi),
                               bounds = bounds.get('cede_colatitude', None),
                               doc = doc,
                               symbol = r'$\Theta$',
                               value = values.get('cede_colatitude', None))

        doc = """
        The angular radius of the (circular) ceding region [radians].
        """
        cede_radius = Parameter('cede_radius',
                                strict_bounds = (0.0, _pi/2.0),
                                bounds = bounds.get('cede_radius', None),
                                doc = doc,
                                symbol = r'$\zeta$',
                                value = values.get('cede_radius', None))

        doc = """
        The azimuth of the centre of the ceding region relative to the
        centre of the superseding region [radians].
        """
        cede_azi = Parameter('cede_azimuth',
                             strict_bounds = (-_pi, _pi),
                             bounds = bounds.get('cede_azimuth', None),
                             doc = doc,
                             symbol = r'$\Phi$',
                             value = values.get('cede_azimuth', None))

        if not custom: # setup default temperature parameters
            doc = """
            log10(superseding region effective temperature [K])
            """
            super_temp = Parameter('super_temperature',
                          strict_bounds = (3.0, 7.0), # very cold --> very hot
                          bounds = bounds.get('super_temperature', None),
                          doc = doc,
                          symbol = r'$\log_{10}(T\;[\rm{K}])$',
                          value = values.get('super_temperature', None))
            if cede:
                doc = """
                log10(ceding region effective temperature [K])
                """
                cede_temp = Parameter('cede_temperature',
                              strict_bounds = (3.0, 7.0), # same story
                              bounds = bounds.get('cede_temperature', None),
                              doc = doc,
                              symbol = r'$\log_{10}(\mathcal{T}\;[\rm{K}])$',
                              value = values.get('cede_temperature', None))
            else:
                cede_temp = None

        else: # let the custom subclass handle definitions; ignore bounds
            super_temp = cede_temp = None

        super(HotRegion, self).__init__(phase_shift, super_colat, super_radius,
                                        super_temp, cede_colat, cede_radius,
                                        cede_azi, cede_temp,
                                        omit_colat, omit_radius, omit_azi,
                                        custom, **kwargs) # prefix in kwargs

    @property
    def objects(self):
        """ Return self for uniform interface with other classes. """
        return [self]

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

    def set_num_cells(self,
                      sqrt_num_cells = 32,
                      min_sqrt_num_cells = 10,
                      max_sqrt_num_cells = 80,
                      fast_sqrt_num_cells = 16,
                      fast_min_sqrt_num_cells = 4,
                      fast_max_sqrt_num_cells = 16):

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

    def set_phases(self, num_leaves,
                   num_phases = None,
                   phases = None,
                   fast_num_leaves = None,
                   fast_num_phases = None,
                   fast_phases = None):
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

        if phases is None:
            try:
                self._phases_cycles = _np.linspace(0.0, 1.0, int(num_phases))
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

        try:
            self._leaves =  _np.linspace(0.0, _2pi, int(num_leaves))
        except TypeError:
            raise TypeError('Number of leaves must be an integer.')

        if self.do_fast:
            if fast_num_phases is None:
                fast_num_phases = fast_num_leaves

            if fast_phases is None:
                try:
                    self._fast_phases_cycles = _np.linspace(0.0, 1.0, int(fast_num_phases))
                except TypeError:
                    raise TypeError('Number of phases must be an integer.')
            else:
                try:
                    assert isinstance(fast_phases, _np.ndarray)
                    assert fast_phases.ndim == 1
                    assert (fast_phases >= 0.0).all() & (fast_phases <= 1.0).all()
                    for i in range(fast_phases.shape[0] - 1):
                        assert fast_phases[i] < fast_phases[i+1]
                except AssertionError:
                    raise TypeError('Phases must be a one-dimensional '
                                    '``numpy.ndarray`` with monotonically '
                                    'increasing elements on the interval [0,1].')
                else:
                    self._fast_phases_cycles = fast_phases

            self._fast_phases = _2pi * self._fast_phases_cycles

            try:
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

    @property
    def is_secondary(self):
        """ Shift the hot region by half a rotational cycle? """
        return self._is_secondary

    @is_secondary.setter
    def is_secondary(self, is_secondary):
        if not isinstance(is_secondary, bool):
            raise TypeError('Use a boolean to specify whether or not the '
                            'hot region should be shifted by half a cycle.')
        else:
            self._is_secondary = is_secondary

    def print_settings(self):
        """ Print numerical settings. """
        print('Base number of cell parallels: ', self.sqrt_num_cells)
        print('Number of rays per parallel: ', self.num_rays)
        print('Number of photospheric leaves: ', self.leaves.shape[0])
        print('Number of interpolation phases: ', self.phases.shape[0])

    def __construct_cellMesh(self, st, fast_total_counts, threads):
        """ Call a low-level routine to construct a mesh representation.

        :param st: Instance of :class:`~.Spacetime.Spacetime`.
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
                                                 self['super_radius'],
                                                 self['cede_radius'],
                                                 self['super_colatitude'],
                                                 self['cede_colatitude'],
                                                 -self['cede_azimuth'],
                                                 self['omit_colatitude'],
                                                 self['omit_radius'],
                                                 self['omit_azimuth'],
                                                 fast_total_counts)

        if (self['super_colatitude'] - self['super_radius'] < 0.0
            or self['super_colatitude'] + self['super_radius'] > _pi):
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
                                          self['super_radius'],
                                          self['super_colatitude'],
                                          self['omit_radius'],
                                          self['omit_colatitude'],
                                          self['omit_azimuth'])

        self._super_phi -= self['omit_azimuth']

        if self._is_secondary:
            self._super_phi += _pi

        if self['cede_radius'] > 0.0:
            if (self['cede_colatitude'] - self['cede_radius'] < 0.0
                or self['cede_colatitude'] + self['cede_radius'] > _pi):
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
                                             self['cede_radius'],
                                             self['cede_colatitude'],
                                             self['super_radius'],
                                             self['super_colatitude'],
                                             -self['cede_azimuth'])

            self._cede_phi += self['cede_azimuth']

            if self._is_secondary:
                self._cede_phi += _pi

    @property
    def __cellArea(self):
        """ Get the areas of cells in the secondary-spot mesh. """
        try:
            return (self._super_cellArea, self._cede_cellArea)
        except AttributeError:
            return (self._super_cellArea, None)

    def __calibrate_lag(self, st, photosphere):
        """ Calibrate lag for cell mesh and normalise by pulse period. """

        R_i = st.R_r_s

        C = (1.0 / self._super_r_s_over_r - R_i)
        C += _log((1.0 / self._super_r_s_over_r - 1.0) / (R_i - 1.0))
        C *= st.r_s / _c
        for j in range(self._super_lag.shape[1]):
            self._super_lag[:,j] -= C
        self._super_lag *= photosphere['mode_frequency']

        try:
            C = (1.0 / self._cede_r_s_over_r - R_i)
        except AttributeError:
            pass
        else:
            C += _log((1.0 / self._cede_r_s_over_r - 1.0) / (R_i - 1.0))
            C *= st.r_s / _c
            for j in range(self._cede_lag.shape[1]):
                self._cede_lag[:,j] -= C
            self._cede_lag *= photosphere['mode_frequency']

    def __compute_rays(self, st, photosphere, threads):
        """ Integrate rays.

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

        self.__calibrate_lag(st, photosphere)

    def __compute_cellParamVecs(self):
        """
        Precompute photospheric source radiation field parameter vectors
        cell-by-cell. Free model parameters and derived (fixed) variables can
        be transformed into local comoving radiation field variables.

        Subclass and overwrite with custom functionality if you desire.

        Designed here simply for uniform effective temperature superseding
        and ceding regions.

        """
        self._super_radiates = _np.greater(self._super_cellArea, 0.0).astype(_np.int32)
        self._super_cellParamVecs = _np.ones((self._super_radiates.shape[0],
                                              self._super_radiates.shape[1],
                                              2),
                                             dtype=_np.double)

        self._super_cellParamVecs[...,:-1] *= self['super_temperature']

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

            self._cede_cellParamVecs[...,:-1] *= self['cede_temperature']

            for i in range(self._cede_cellParamVecs.shape[1]):
                self._cede_cellParamVecs[:,i,-1] *= self._cede_effGrav

    @staticmethod
    def psi(theta, phi, colatitude):
        """ Coordinates of cell centres in rotated spherical coordinate system.

        Transformation is anticlockwise rotation about y-axis of Cartesian
        basis.

        """
        return _arccos(_cos(colatitude)
                       * _cos(theta)
                       + _sin(colatitude)
                       * _sin(theta)
                       * _cos(phi))

    def embed(self, spacetime, photosphere, fast_total_counts, threads, *args):
        """ Embed the hot region of the photosphere into the ambient spacetime.

        :param tuple args:
            Correct the integral over the radiation field *elsewhere* by
            accounting for the time-dependent component arising from the
            presence of the hot region.

        """
        if self.fast_mode and not self.do_fast:
            return None
        elif not self.needs_update: # dynamically evaluate if stuff to do
            return None # what about change in mesh and ray resolution?

        self.__construct_cellMesh(spacetime,
                                  fast_total_counts,
                                  threads)

        self.__compute_rays(spacetime, photosphere, threads)
        self.__compute_cellParamVecs()

        if args:
            self._super_correctionVecs = args[0](self._super_theta)
            for i in range(self._super_theta.shape[1]):
                self._super_correctionVecs[:,i,-1] *= self._super_effGrav

            try:
                self._cede_correctionVecs = args[0](self._cede_theta)
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

    @property
    def cede(self):
        """ Does the hot region have a ceding member? """
        return self._cede

    @cede.setter
    def cede(self, cede):
        if not isinstance(cede, bool):
            raise TypeError('A boolean is required to activate or deactivate '
                            'the ceding region.')
        else:
            self._cede = cede

    @property
    def omit(self):
        """ Does the hot region defined through an omission region? """
        return self._omit

    @omit.setter
    def omit(self, omit):
        if not isinstance(omit, bool):
            raise TypeError('A boolean is required to activate or deactivate '
                            'the omission region.')
        else:
            self._omit = omit

    @property
    def concentric(self):
        """ Is the superseding region concentric with ceding member?

        If not, is the superseding region concentric with an omission member?

        """
        return self._concentric

    @concentric.setter
    def concentric(self, concentric):
        if not isinstance(concentric, bool):
            raise TypeError('A boolean is required to activate or deactivate '
                            'concentricity.')
        else:
            self._concentric = concentric

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

HotRegion._update_doc()
