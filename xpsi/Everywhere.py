from __future__ import division, print_function

from .global_imports import *
from . import global_imports

from .cellmesh.global_mesh import construct_closed_cellMesh as _construct_closed_cellMesh
from .cellmesh.rays import compute_rays as _compute_rays

from .Parameter import Parameter
from .ParameterSubspace import ParameterSubspace

class RayError(xpsiError):
    """ Raised if a problem was encountered during ray integration. """

class IntegrationError(xpsiError):
    """ Raised if a problem was encountered during signal integration. """

class Everywhere(ParameterSubspace):
    """ The photospheric radiation field represented by a global mesh.

    An instance of this class *cannot* be used in conjuction with hot region
    objects nor with an elsewhere instance.

    The local comoving radiation field properties are *not* assumed to be
    azimuthally invariant, so a time-dependent signal can be generated upon
    integration over the image.

    However, if you want to generate a time-dependent signal, you need
    to customise the parameters that control the radiation field w.r.t. surface
    local comoving frames. By default a single temperature parameter is
    defined which is globally invariant and thus the radiation field is
    azimuthally invariant and the image subtended on a distant observer's sky
    is phase-invariant.

    :param bool time_invariant:
        See above. Is the radiation field azimuthally invariant, thus
        generating a phase-invariant signal? Although the default configuration
        is simple and azimuthally invariant, we do not provide a default
        value for this argument to force the user to verify their intentions.

    :param int num_rays:
        Number of rays to trace (integrate) at each colatitude, distributed
        in angle subtended between ray tangent 4-vector and radial outward unit
        vector w.r.t a local orthonormal tetrad.

    :param int sqrt_num_cells:
        Number of cells in both colatitude and azimuth which form a regular
        mesh on the surface. The total number of cells is the square of this
        argument value. The mesh suffers from squeezing in the polar regions,
        leading to a high degree of non-congruence in cell shape over the
        surface.

    :param iterable custom:
        Iterable over :class:`~.Parameter.Parameter` instances. If you
        supply custom parameter definitions, you need to overwrite the
        :func:`~.Everywhere.Everwhere._compute_cellParamVecs` method to
        implement your custom behaviour.

    :param dict bounds:
        If ``custom is None``, these bounds are supplied for instantiation
        of a temperature parameter. The parameter name
        ``'temperature'`` must be a key in the dictionary unless the
        parameter is *fixed* or *derived*. If a bound is ``None`` that bound
        is set equal to a strict hard-coded bound.

    :param dict values:
        Either the fixed value of the temperature, a callable if the
        temperature is *derived*, or a value upon initialisation if the
        temperature is free. The dictionary must have a key with name
        ``'temperature'`` if it is *fixed* or *derived*.

    """
    required_names = ['temperature (if no custom specification)']

    def __init__(self,
                 time_invariant,
                 bounds = {},
                 values = {},
                 sqrt_num_cells = 64,
                 num_rays = 200,
                 num_leaves = 100,
                 num_phases = None,
                 phases = None,
                 custom = None):

        self.num_rays = num_rays

        self.sqrt_num_cells = sqrt_num_cells

        self.set_phases(num_leaves, num_phases, phases)

        if not custom: # setup default temperature parameter
            T = Parameter('temperature',
                          strict_bounds = (3.0, 7.0), # very cold --> very hot
                          bounds = bounds.get('temperature', None),
                          doc = 'log10(effective temperature [K] everywhere)',
                          symbol = r'$\log_{10}(T\;[\rm{K}])$',
                          value = values.get('temperature', None))
        else: # let the custom subclass handle definitions; ignore bounds
            T = None

        super(Everywhere, self).__init__(T, custom)

        self.time_invariant = time_invariant

    @property
    def time_invariant(self):
        """ Get the declared state of time-invariance. """
        return self._time_invariant

    @time_invariant.setter
    def time_invariant(self, invariant):
        """ Declare whether the signal is time-invariant. """
        # find the required integrator
        if invariant: # can we safely assume azimuthal invariance?
            self._time_invariant = True
            from .cellmesh.integrator_for_time_invariance import integrate as _integrator
        else: # more general purpose
            self._time_invariant = False
            from .cellmesh.integrator import integrate as _integrator
        self._integrator = _integrator

    @property
    def num_rays(self):
        """ Get the number of rays integrated per colatitude. """
        return self._num_rays

    @num_rays.setter
    def num_rays(self, n):
        """ Set the number of rays integrated per colatitude. """
        try:
            self._num_rays = int(n)
        except TypeError:
            raise TypeError('Number of rays must be an integer.')

    @property
    def sqrt_num_cells(self):
        """ Get the number of cell colatitudes. """
        return self._sqrt_num_cells

    @sqrt_num_cells.setter
    def sqrt_num_cells(self, n):
        """ Set the number of cell colatitudes. """
        try:
            self._sqrt_num_cells = int(n)
        except TypeError:
            raise TypeError('Number of cells must be an integer.')
        else:
            self._num_cells = n**2

    @property
    def num_cells(self):
        """ Get the total number of cells in the mesh. """
        return self._num_cells

    def print_settings(self):
        """ Print numerical settings. """
        print('Base number of cell parallels: ', self.sqrt_num_cells)
        print('Number of rays per parallel: ', self.num_rays)
        print('Number of photospheric leaves: ', len(self.leaves))
        print('Number of interpolation phases: ', len(self.phases))

    def set_phases(self, num_leaves,
                   num_phases = None,
                   phases = None):
        """ Construct the set of interpolation phases and foliation leaves.

        :param int num_leaves:
            Number of leaves mesh motion is discretised into.

        :param int num_phases:
            Number of phases in a discrete representation of the specific flux
            pulses on the interval :math:`[0,1]`. If ``None``, the number of
            phases interpolated at is set equal to the number of leaves.

        :param phases:
            If not ``None``, a :class:`numpy.ndarray` of phases for a discrete
            representation of the specific flux pulses on the interval
            :math:`[0,1]`. If ``None`` (default), the :obj:`num_phases`
            argument is utilised.

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

    @property
    def phases_in_cycles(self):
        """ Get the phases (in cycles) the pulse is interpolated at. """
        return self._phases_cycles

    def _construct_cellMesh(self, st, threads):
        """ Call a low-level routine to construct a mesh representation.

        :param st: Instance of :class:`~.Spacetime.Spacetime`.

        :param int threads: Number of ``OpenMP`` threads for mesh construction.

        """
        (self._theta,
         self._phi,
         self._r,
         self._cellArea,
         self._maxAlpha,
         self._cos_gamma,
         self._effGrav) = _construct_closed_cellMesh(threads,
                                                     self._sqrt_num_cells,
                                                     self._num_cells,
                                                     st.M,
                                                     st.r_s,
                                                     st.R,
                                                     st.zeta,
                                                     st.epsilon)

    def _calibrate_lag(self, st, photosphere):
        """ Calibrate lag for cell mesh and normalise by pulse period. """

        R_i = st.R_r_s

        C = (1.0 / self._r_s_over_r - R_i)
        C += _log((1.0 / self._r_s_over_r - 1.0) / (R_i - 1.0))
        C *= st.r_s / _c
        for j in range(self._lag.shape[1]):
            self._lag[:,j] -= C
        self._lag *= photosphere['mode_frequency'] * _2pi

    def _compute_rays(self, st, photosphere, threads):
        """ Trace (integrate) a set of rays.

        These rays represent a null mapping from photosphere to a point at
        some effective infinity.

        :param st: Instance of :class:`~.Spacetime.Spacetime`.

        :param int threads: Number of ``OpenMP`` threads for ray integration.

        """
        self._r_s_over_r = _contig(st.r_s / self._r, dtype = _np.double)

        (terminate_calculation,
         self._deflection,
         self._cos_alpha,
         self._lag,
         self._maxDeflection) = _compute_rays(threads,
                                              self._sqrt_num_cells,
                                              st.r_s,
                                              self._r_s_over_r,
                                              self._maxAlpha,
                                              self._num_rays)

        self._calibrate_lag(st, photosphere)

        if terminate_calculation == 1:
            raise RayError('Fatal numerical problem during ray integration.')

    def _compute_cellParamVecs(self):
        """
        Precompute photospheric source radiation field parameter vectors
        cell-by-cell. Free model parameters and derived (fixed) variables can
        be transformed into local comoving radiation field variables.

        Subclass and overwrite with custom functionality if you desire.

        """
        # all radiate, but can be changed with overwrite
        self._cellRadiates = _np.ones(self._theta.shape, dtype=_np.int32)

        self._cellParamVecs = _np.ones((self._theta.shape[0],
                                        self._theta.shape[1],
                                        len(self.vector)+1),
                                       dtype=_np.double)

        self._cellParamVecs[...,:-1] *= _np.array(self.vector)

        for i in range(self._cellParamVecs.shape[1]):
            self._cellParamVecs[:,i,-1] *= self._effGrav

    def embed(self, spacetime, photosphere, threads):
        """ Embed the photosphere everywhere. """

        self._construct_cellMesh(spacetime, threads)
        self._compute_rays(spacetime, photosphere, threads)
        self._compute_cellParamVecs()

    def integrate(self, st, energies, threads, atmosphere):
        """ Integrate over the photospheric radiation field.

        :param st: Instance of :class:`~.Spacetime.Spacetime`.

        :param energies: A one-dimensional :class:`numpy.ndarray` of energies
                         in keV.

        :param int threads: Number of ``OpenMP`` threads the integrator is
                            permitted to spawn.

        """
        if self._time_invariant:
            out = self._integrator(threads,
                                   st.R,
                                   st.Omega,
                                   st.r_s,
                                   st.i,
                                   self._sqrt_num_cells,
                                   self._cellArea,
                                   self._r,
                                   self._r_s_over_r,
                                   self._theta,
                                   self._phi,
                                   self._cellParamVecs,
                                   self._num_rays,
                                   self._deflection,
                                   self._cos_alpha,
                                   self._maxDeflection,
                                   self._cos_gamma,
                                   energies,
                                   atmosphere)
        else:
            out = self._integrator(threads,
                                   st.R,
                                   st.Omega,
                                   st.r_s,
                                   st.i,
                                   _np.ones(self._theta.shape) * self._cellArea,
                                   self._r,
                                   self._r_s_over_r,
                                   self._theta,
                                   self._phi,
                                   self._cellParamVecs,
                                   self._cellRadiates,
                                   None, # no intensity correction required
                                   self._num_rays,
                                   self._deflection,
                                   self._cos_alpha,
                                   self._lag,
                                   self._maxDeflection,
                                   self._cos_gamma,
                                   energies,
                                   self._leaves,
                                   self._phases,
                                   atmosphere,
                                   None) # no other atmosphere needed
        if out[0] == 1:
            raise IntegrationError('Fatal numerical error during integration.')

        return out[1]

Everywhere._update_doc()
