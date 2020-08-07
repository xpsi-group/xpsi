from __future__ import division, print_function

from .global_imports import *
from . import global_imports

from .cellmesh.global_mesh import construct_closed_cellMesh as _construct_closed_cellMesh
from .cellmesh.rays import compute_rays as _compute_rays
from .cellmesh.integrator_for_time_invariance import integrate as _integrator

from .Parameter import Parameter
from .ParameterSubspace import ParameterSubspace

class RayError(xpsiError):
    """ Raised if a problem was encountered during ray integration. """

class IntegrationError(xpsiError):
    """ Raised if a problem was encountered during signal integration. """

class Elsewhere(ParameterSubspace):
    """ The photospheric radiation field *elsewhere*.

    This means the radiation field exterior to the hot regions. The local
    comoving radiation field properties are *assumed* (for now) to be
    azimuthally invariant but can in principle vary colatitudinally.

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
        :func:`~.Elsewhere.Elsewhere._compute_cellParamVecs` method to
        handle your custom behaviour.

    :param dict bounds:
        If ``custom is None``, these bounds are supplied for instantiation
        of a temperature parameter. The parameter name
        ``'elsewhere_temperature'`` must be a key in the dictionary unless the
        parameter is *fixed* or *derived*. If a bound is ``None`` that bound
        is set equal to a strict hard-coded bound.

    :param dict values:
        Either the fixed value of the temperature elsewhere, a callable if the
        temperature is *derived*, or a value upon initialisation if the
        temperature is free. The dictionary must have a key with name
        ``'elsewhere_temperature'`` if it is *fixed* or *derived*.

    """
    required_names = ['elsewhere_temperature (if no custom specification)']

    def __init__(self,
                 sqrt_num_cells = 64,
                 num_rays = 1000,
                 bounds = None,
                 values = None,
                 custom = None):

        self.sqrt_num_cells = sqrt_num_cells
        self.num_rays = num_rays

        if bounds is None: bounds = {}
        if values is None: values = {}

        if not custom: # setup default temperature parameter
            T = Parameter('elsewhere_temperature',
                          strict_bounds = (3.0, 7.0), # very cold --> very hot
                          bounds = bounds.get('elsewhere_temperature', None),
                          doc = 'log10 of the effective temperature elsewhere',
                          symbol = r'$\log_{10}(T_{\rm EW}\;[\rm{K}])$',
                          value = values.get('elsewhere_temperature', None))
        else: # let the custom subclass handle definitions; ignore bounds
            T = None

        super(Elsewhere, self).__init__(T, custom)

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
        print('Number of cell colatitudes: ', self.sqrt_num_cells)
        print('Number of rays per colatitude: ', self.num_rays)

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

    def _compute_rays(self, st, threads):
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

        if terminate_calculation == 1:
            raise RayError('Fatal numerical problem during ray integration.')

    def _compute_cellParamVecs(self, *args):
        """
        Precompute photospheric source radiation field parameter vectors
        cell-by-cell. Free model parameters and derived (fixed) variables can
        be transformed into local comoving radiation field variables.

        Subclass and overwrite with custom functionality if you desire.

        :param tuple args:
            An *ndarray[n,n]* of mesh-point colatitudes.

        """
        if args: # hot region mesh shape information
            cellParamVecs = _np.ones((args[0].shape[0],
                                      args[0].shape[1],
                                      len(self.vector)+1),
                                     dtype=_np.double)

            # get self.vector because there may be fixed variables
            # that also need to be directed to the integrators
            # for intensity evaluation
            cellParamVecs[...,:-1] *= _np.array(self.vector)

            return cellParamVecs

        else:
            self._cellParamVecs = _np.ones((self._theta.shape[0],
                                            self._theta.shape[1],
                                            len(self.vector)+1),
                                           dtype=_np.double)

            self._cellParamVecs[...,:-1] *= _np.array(self.vector)

            for i in range(self._cellParamVecs.shape[1]):
                self._cellParamVecs[:,i,-1] *= self._effGrav

    def embed(self, spacetime, threads):
        """ Embed the photosphere elsewhere. """

        self._construct_cellMesh(spacetime, threads)
        self._compute_rays(spacetime, threads)
        self._compute_cellParamVecs()

    def integrate(self, st, energies, threads, *atmosphere):
        """ Integrate over the photospheric radiation field.

        :param st: Instance of :class:`~.Spacetime.Spacetime`.

        :param energies: A one-dimensional :class:`numpy.ndarray` of energies
                         in keV.

        :param int threads: Number of ``OpenMP`` threads the integrator is
                            permitted to spawn.

        """
        if isinstance(energies, tuple): # resolve energy container type
            if not isinstance(energies[0], tuple):
                _energies = energies[0]
            else:
                _energies = energies[0][0]
        else:
            _energies = energies

        out = _integrator(threads,
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
                           _energies,
                           atmosphere)
        if out[0] == 1:
            raise IntegrationError('Fatal numerical error during elsewhere integration.')

        return out[1]

Elsewhere._update_doc()
