from __future__ import division, print_function

from .global_imports import *
from . import global_imports

from .cellmesh.global_mesh import construct_closed_cellMesh as _construct_closed_cellMesh
from .cellmesh.rays import compute_rays as _compute_rays
from .cellmesh.integrator_for_time_invariance import integrate_radField as _integrator

from .ParameterSubspace import ParameterSubspace

class RayError(xpsiError):
    """ Raised if a numerical problems encountered during ray integration. """

class IntegrationError(xpsiError):
    """ Raised if a numerical problems encountered during integration. """

class Elsewhere(ParameterSubspace):
    """ The photospheric radiation field elsewhere (other than the spot).

    :param int num_params: Number of parameters for model of photospheric
                           radiation field elsewhere.

    :param list bounds: Hard parameter bounds for the instance of
                        :class:`.ParameterSubspace.ParameterSubspace`.

    :param int num_rays: Number of rays to trace (integrate) at each
                         colatitude, distributed in angle subtended between
                         ray tangent 4-vector and radial outward unit
                         vector w.r.t a local orthonormal tetrad.

    :param int sq_num_cells: Number of cells in both colatitude and azimuth
                             which form a regular mesh on a closed curved
                             2-surface (which is a compact subset of a
                             spacelike leaf of the spacetime foliation).

    """

    def __init__(self,
                 num_params,
                 bounds,
                 num_rays = 1000,
                 sq_num_cells = 64):
        super(Elsewhere, self).__init__(num_params, bounds)

        self.num_rays = num_rays
        self.sq_num_cells = sq_num_cells

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

    @property
    def sq_num_cells(self):
        """ Get the number of cell parallels. """
        return self._num_rays

    @sq_num_cells.setter
    def sq_num_cells(self, n):
        """ Set the number of cell parallels. """
        try:
            self._sq_num_cells = int(n)
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
        print('Number of cell parallels: ', self.sq_num_cells)
        print('Number of rays per parallel: ', self.num_rays)

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
                                                     self._sq_num_cells,
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
                                              self._sq_num_cells,
                                              st.r_s,
                                              self._r_s_over_r,
                                              self._maxAlpha,
                                              self._num_rays)

        if terminate_calculation == 1:
            raise RayError('Fatal numerical problem during ray integration.')

    def _compute_cellParamVecs(self, p):
        """
        Precomputes photospheric source radiation field parameter vectors
        cell-by-cell.

        :param list p: Arguments for the
                       :meth:`_eval_srcRadFieldParamVectors` method.

        :param optional spot: An instance of :class:`.Spot`. If ``None``
                              the ambient source radiation field is evaluated
                              in terms of parameter vectors at the centres of
                              cells defined in the ambient cell mesh. If
                              :obj:`spot` is an instance, the parameter vectors
                              are evaluated at the centres of cells defined
                              in the spot cell mesh.
        """
        self._cellParamVecs = self.eval_srcRadFieldParamVectors(self._theta.shape, p)

        for i in range(self._cellParamVecs.shape[1]):
            self._cellParamVecs[:,i,-1] *= self._effGrav

    @staticmethod
    def eval_srcRadFieldParamVectors(shape, p):
        """
        .. note:: Default. Can be overwritten in a custom subclass, *if you know
                  what you are doing*\ .

        The ambient field is uniform with respect to local Lagrangian tetrads.
        Parameters are matrices, and a matrix must be returned.

        Function to evaluate the source radiation field parameter vector at
        points on a closed 2-surface (a spacelike leaf of the spacetime
        foliation).

        :param tuple shape: Dimensions of a 2D array of points.

        :param list p: The local radiation field parameter vector, which
                       is to be passed to the low-level routines.

        """
        cellParamVecs = _np.ones((shape[0], shape[1], len(p)+1),
                                 dtype=_np.double)

        cellParamVecs[...,:-1] *= _np.array(p)

        return cellParamVecs

    def embed(self, spacetime, p, threads):
        """ Embed the photosphere elsewhere. """

        self._construct_cellMesh(spacetime, threads)
        self._compute_rays(spacetime, threads)
        self._compute_cellParamVecs(p)

    def integrate(self, st, energies, threads, *atmosphere):
        """ Integrate over the photospheric radiation field.

        :param st: Instance of :class:`~.Spacetime.Spacetime`.

        :param energies: A one-dimensional :class:`numpy.ndarray` of energies
                         in keV.

        :param int threads: Number of ``OpenMP`` threads the integrator is
                            permitted to spawn. 

        """
        out = _integrator(threads,
                           st.R,
                           st.Omega,
                           st.r_s,
                           st.i,
                           self._sq_num_cells,
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
                           *atmosphere)
        if out[0] == 1:
            raise IntegrationError('Fatal numerical error during elsewhere integration.')

        return out[1]
