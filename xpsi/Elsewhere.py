from xpsi.global_imports import *

from xpsi.cellmesh.global_mesh import construct_closed_cellMesh as _construct_closed_cellMesh
from xpsi.cellmesh.rays import compute_rays as _compute_rays
from xpsi.cellmesh.integrator_for_time_invariance import integrate as _integrator

from xpsi.Parameter import Parameter
from xpsi.ParameterSubspace import ParameterSubspace

class AtmosError(xpsiError):
    """ Raised if the numerical atmosphere data were not preloaded. """

class RayError(xpsiError):
    """ Raised if a problem was encountered during ray integration. """

class IntegrationError(xpsiError):
    """ Raised if a problem was encountered during signal integration. """

class Elsewhere(ParameterSubspace):
    """ The photospheric radiation field *elsewhere*.

    This means the radiation field exterior to the hot regions. The local
    comoving radiation field properties are *assumed* (for now) to be
    azimuthally invariant but can in principle vary colatitudinally.

    :param int sqrt_num_cells:
        Number of cells in both colatitude and azimuth which form a regular
        mesh on the surface. Must be an even number such that half of the cells
        are exactly in one hemisphere. The total number of cells is the square
        argument value. The mesh suffers from squeezing in the polar regions,
        leading to a high degree of non-congruence in cell shape over the
        surface.

    :param int num_rays:
        Number of rays to trace (integrate) at each colatitude, distributed
        in angle subtended between ray tangent 4-vector and radial outward unit
        vector w.r.t a local orthonormal tetrad.

    :param dict bounds:
        If ``custom is None``, these bounds are supplied for instantiation
        of a temperature parameter. The parameter name
        ``'elsewhere_temperature'`` must be a key in the dictionary unless the
        parameter is *fixed* or *derived*. If a bound is ``None`` that bound
        is set equal to a strict hard-coded bound. We note that the bounds for
        parameters used in the atmosphere model should be restricted (by the user)
        to be within the tabulated values, in case a numerical atmosphere extension is used.

    :param dict values:
        Either the fixed value of the temperature elsewhere, a callable if the
        temperature is *derived*, or a value upon initialisation if the
        temperature is free. The dictionary must have a key with name
        ``'elsewhere_temperature'`` if it is *fixed* or *derived*.

    :param str atm_ext:
        Used to determine which atmospheric extension to use.
        Options at the moment:
        "BB": Analytical blackbody (default)
        "Num4D": Numerical atmosphere using 4D-interpolation from the provided
        atmosphere data
        "user": A user-provided extension which can be set up by replacing the contents of 
        the file hot_user.pyx (and elsewhere_user.pyx if needed) and re-installing X-PSI
        (if not changed, "user" is the same as "BB").

    :param iterable custom:
        Iterable over :class:`~.Parameter.Parameter` instances. If you
        supply custom parameter definitions, you need to overwrite the
        :func:`~.Elsewhere.Elsewhere._compute_cellParamVecs` method to
        handle your custom behaviour.

    :param int image_order_limit:
        The highest-order image to sum over. A value of *one* means primary
        images only (deflections :math:`<\pi`) whilst a value of *two* means
        primary and secondary images (deflections :math:`<2pi`) where visible,
        and so on. If ``None`` (the default), there is no hard limit. In this
        case the limit is determined quasi-naturally for each mesh element,
        meaning that images will be summed over until higher order images are
        not visible or the visibility limit is truncated due to lack of
        numerical precision (e.g. for rays that orbit very close to the
        Schwarzschild photon sphere three times or more).  Higher-order images
        generally contribute less and less due to geometric projection effects
        (higher-order images become more tangential), and the images of
        elements get squeezed in solid angle at the stellar limb. In principle,
        effects such as relativistic beaming can counter this effect to a
        degree for certain source-receiver configurations, by increasing
        brightness whilst solid angle decreases, and thus the flux contribution
        relative to that from a primary image can be greater than suggested
        simply by geometric project effects. Nevertheless, inclusion of these
        images is more computationally expensive. If, when iterating through
        image orders, an image is not visible because the deflection required
        is greater than the highest deflection permitted at a given colatitude
        on a surface (accounting for the surface tilt due to rotation), then
        the iteration over image orders terminates.

    """
    required_names = ['elsewhere_temperature (if no custom specification)']

    def __init__(self,
                 sqrt_num_cells = 64,
                 num_rays = 1000,
                 bounds = None,
                 values = None,
                 atm_ext="BB",
                 custom = None,
                 image_order_limit = None):

        self.sqrt_num_cells = sqrt_num_cells
        self.num_rays = num_rays
        self.image_order_limit = image_order_limit
        self.atm_ext = atm_ext

        if bounds is None: bounds = {}
        if values is None: values = {}

        if not custom: # setup default temperature parameter
            T = Parameter('elsewhere_temperature',
                          strict_bounds = (3.0, 7.6), # very cold --> very hot
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
             _n = int(n)
        except TypeError:
            raise TypeError('Number of cells must be an integer.')
        else:
            if not _n >= 10 or _n%2 != 0:
                raise ValueError('Number of cells must be a positive even '
                                 'integer greater than or equal to ten.')
        self._sqrt_num_cells = _n
        self._num_cells = _n**2

    @property
    def num_cells(self):
        """ Get the total number of cells in the mesh. """
        return self._num_cells

    @property
    def atm_ext(self):
        """ ... """
        return self._atm_ext

    @atm_ext.setter
    def atm_ext(self,extension):
        if extension=="BB":
            self._atm_ext = 1
        elif extension=="Num4D":
            self._atm_ext = 2
        elif extension=="user":
            self._atm_ext  = 3
        else:
            raise TypeError('Got an unrecognised atm_ext argument. Note that the only allowed '
                            'atmosphere options are at the moment "BB", "Num4D", and "user".')

    @property
    def image_order_limit(self):
        """ Get the image order limit. """
        return self._image_order_limit

    @image_order_limit.setter
    def image_order_limit(self, limit):
        """ Set the image order limit. """
        if limit is not None:
            if not isinstance(limit, int):
                raise TypeError('Image order limit must be an positive integer '
                                'if not None.')
        self._image_order_limit = limit

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

        if self._atm_ext==2:
            if atmosphere == ():
                raise AtmosError('The numerical atmosphere data were not preloaded, '
                                 'even though that is required by the current atmosphere extension.')

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
                           atmosphere,
                           self._atm_ext,
                           self._image_order_limit)
        if out[0] == 1:
            raise IntegrationError('Fatal numerical error during elsewhere integration.')

        return out[1]

Elsewhere._update_doc()
