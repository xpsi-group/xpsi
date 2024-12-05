#cython: cdivision=True
#cython: boundscheck=False
#cython: nonecheck=False
#cython: wraparound=False
#cython: embedsignature=True

""" Extensions to compute specific intensities in a local comoving frame.

Overview
--------

We require the surface radiation field model to be implemented in C extensions
to Python if signal computation is to be sufficiently fast. Signal computation
includes:

    * integrating over the image of a star on a distant observer's sky,
      yielding a phase-energy resolved signal for likelihood function
      evaluation;
    * resolving the image of star on a distant observer's sky as a function of
      phase and energy for visualisation purposes, and for more general
      simulation purposes.

The first type of signal computation is enabled by surface discretization,
whilst the latter is enabled via image-plane discretization. For surface
discretization there are several extension modules for users and developers to
work with; ``hot_BB.pyx``, and ``hot_Num4D.pyx`` are the built-in options for
isotropic blackbody radiation field and numerical 4D-interpolation from a preloaded
atmosphere table (they can be selected using appropriate flags when creating
instances of the radiating regions); ``hot_user.pyx`` and ``elsewhere_user.pyx``
are additional extensions which can be replaced by user-developed versions, and
used after re-installing X-PSI (by default they are the same as ``hot_BB.pyx``).
These extension modules contain functions whose bodies enable evaluation of
specific intensities with respect to local comoving frames at the stellar surface.
The ``hot_*.pyx`` modules define the radiation field inside surface hot regions
represented by :class:`~xpsi.HotRegion.HotRegion` instances; however, the same
modules are used by an instance of the :class:`~xpsi.Everywhere.Everywhere` class.
The ``elsewhere_user.pyx`` extension module defines the customized radiation
field *exterior* to the :class:`~xpsi.HotRegion.HotRegion` instances on the surface.
In case no customized module for the exterior surface is needed, the surface
radiation field for exterior surface can be selected from one of the built-in
options for the hot regions.

For image-plane discretization, the ``hot_*.pyx`` extensions define the local
surface radiation field properties with respect to a comoving frame. However,
to compute the variables required by this module to evaluate specific intensity,
another extension module is required with transforms a set of *global*
variables into local variables given spacetime coordinates of an event
belonging to a null geodesic connecting the surface to a radiation-detecting
distant instrument. Users and developers must customize the function bodies
in ``local_variables.pyx`` in order to implement more advanced physics models.

For all of these modules, numerical data can be read from disk to define the
functions for specific intensity evaluation given local variables, and for
transformation of global variables to local variables. A common usage pattern
is to statistically model data assuming some numerical physics model for the
radiation emergent from an atmosphere. Statistical modeling requires many
likelihood function evaluations, and therefore to avoid unnecessarily reading
numerical data from disk every time the likelihood function is called, users
can preload the numerical atmosphere data in an encapsulation Python process
and point to the data structures in memory when computing signals; developers
can develop extension modules to work with preloaded numerical atmosphere data.
When preloading an atmosphere, the number of grid points in each dimension
is dynamically inferred, but the existing source code for multi-dimensional
interpolation (``hot_Num4D.pyx``) is hard-coded for four dimensions. With
development, this can be made more sophisticated.

Developers can contribute extension modules to the
``surface_radiation_field/archive``. Extensions from this archive must replace
the contents of the ``hot_user.pyx``, ``elsewhere_user.pyx``, and ``local_variables.pyx``
modules, adhering to the function prototypes defined in the associated
header files and using existing examples in the archive for guidance.

Use the Python-exposed functions below to test the C implementations of
surface radiation field models.

"""

from cython.parallel cimport *

import numpy as np
cimport numpy as np

from ..global_imports import _keV

cdef double keV = _keV

from .preload cimport (_preloaded,
                       init_preload,
                       free_preload)

from .hot_wrapper cimport (init_hot,
                   eval_hot_I,
                   eval_hot_Q,
                   eval_hot_norm,
                   free_hot)

from .elsewhere_wrapper cimport (init_elsewhere,
                         free_elsewhere,
                         eval_elsewhere,
                         eval_elsewhere_norm)

from .local_variables cimport (storage,
                               HIT_or_MISS,
                               init_local_variables,
                               free_local_variables,
                               eval_local_variables)

from .effective_gravity_universal cimport effectiveGravity

from xpsi.pixelmesh.geometricConfiguration cimport _GEOM

ctypedef void* (*fptr_init)(size_t, const _preloaded *const, size_t) nogil

ctypedef int (*fptr_free)(size_t, void *const) nogil

ctypedef double (*fptr_eval)(size_t,
                             double,
                             double,
                             const double *const,
                             void *const,
                             size_t) nogil

ctypedef double (*fptr_norm)() nogil

def intensity(double[::1] energies,
              double[::1] mu,
              double[:,::1] local_variables,
              atmosphere = None,
              int stokesQ = 0,
              region_extension = 'hot',
              atmos_extension = "BB",
              beam_opt = 0,
              size_t numTHREADS = 1):
    """ Evaluate the photon specific intensity using an extension module.

    The intensities are computed directly from local variable values.

    :param double[::1] energies:
        A 1D :class:`numpy.ndarray` of energies in keV to evaluate photon
        specific intensity at, in the local comoving frame.

    :param double[::1] mu:
        A 1D :class:`numpy.ndarray` of angles at which to evaluate the photon
        specific intensity at. Specifically, the angle is the cosine of the
        angle of the ray direction to the local surface normal, in the local
        comoving frame.

    :param double[:,::1] local_variables:
        A 2D :class:`numpy.ndarray` (with beam_opt=0) or 7D :class:`numpy.ndarray`
        (with beam_opt > 0) of local variables such as temperature, effective gravity,
        and beaming parameters. Rows correspond to the sequence of points in the
        space of energy and angle specified above, and columns contain the
        required variable values, one variable per column, in the expected
        order.

    :param tuple atmosphere:
        A numerical atmosphere preloaded into buffers. The buffers must be
        supplied in an :math:`n`-tuple whose :math:`n^{th}` element is an
        :math:`(n-1)`-dimensional array flattened into a one-dimensional
        :class:`numpy.ndarray`. The first :math:`n-1`
        elements of the :math:`n`-tuple must each be an ordered one-dimensional
        :class:`numpy.ndarray` of parameter values for the purpose of
        multi-dimensional interpolation in the :math:`n^{th}` buffer. The
        first :math:`n-1` elements must be ordered to match the index
        arithmetic applied to the :math:`n^{th}` buffer. An example would be
        ``(logT, logg, mu, logE, buf)``, where:
        ``logT`` is a logarithm of local comoving effective temperature;
        ``logg`` is a logarithm of effective surface gravity;
        ``mu`` is the cosine of the angle from the local surface normal;
        ``logE`` is a logarithm of the photon energy; and
        ``buf`` is a one-dimensional buffer of intensities of size given by
        the product of sizes of the first :math:`n-1` tuple elements.

        It is highly recommended that buffer preloading is used, instead
        of loading from disk in the customisable radiation field extension
        module, to avoid reading from disk for every signal
        (likelihood) evaluation. This can be a non-negligible waste of compute
        resources. By preloading in Python, the memory is allocated and
        references to that memory are not in general deleted until a sampling
        script exits and the kernel stops. The likelihood callback accesses
        the same memory upon each call without I/O.

    :param int StokesQ:
        If StokesQ=1, Stokes Q intensities will be returned. Otherwise Stokes I
        intensities. Default is StokesQ=0.

    :param str region_extension:
        Specify the radiating region extension module to invoke. Options are ``'hot'`` and
        ``'elsewhere'``.

    :param str atmos_extension:
        Used to determine which atmospheric extension to use.
        Options at the moment:
        "BB": Analytical blackbody,
        "Num4D": Numerical atmosphere using 4D-interpolation from the provided
        atmosphere data,
        "Pol_BB_burst": Polarized analytical blackbody+burst approximation,
        "Pol_Num2D": Polarized numerical atmosphere using 2D-interpolation from the provided
        atmosphere data,
        "user": A user-provided extension which can be set up by replacing the contents of
        the file hot_user.pyx (and elsewhere_user.pyx if needed) and re-installing X-PSI
        (if not changed, "user" is the same as "BB").

    :param int beam_opt:
        Used to determine which atmospheric beaming modification model to use.
        Options at the moment:
        0: No modification (default)
        1: Original*beaming_correction without re-normalization
        2: Original*beaming_correction with analytical re-normalization estimate
        3: Original*beaming_correction with numerical re-normalization

    :param int numTHREADS:
        Number of OpenMP threads to launch.

    :returns:
        A 1D :class:`numpy.ndarray` of the photon specific intensities in
        units of photons/s/keV/cm^2/sr.

    """
    cdef size_t _beam_opt = beam_opt
    cdef size_t _atmos_extension

    cdef fptr_init init_ptr = NULL
    cdef fptr_free free_ptr = NULL
    cdef fptr_eval eval_ptr = NULL
    cdef fptr_norm norm_ptr = NULL

    if region_extension == 'hot':
        init_ptr = init_hot
        free_ptr = free_hot
        eval_ptr_I = eval_hot_I
        eval_ptr_Q = eval_hot_Q
        norm_ptr = eval_hot_norm
    elif region_extension == 'elsewhere':
        init_ptr = init_elsewhere
        free_ptr = free_elsewhere
        eval_ptr_I = eval_elsewhere
        norm_ptr = eval_elsewhere_norm
        if stokesQ == 1:
            raise ValueError("StokesQ option is not allowed for the elsewhere extension.")
    else:
        raise ValueError("Region extension module must be 'hot' or 'elsewhere'.")

    # initialise the source radiation field
    cdef _preloaded *preloaded = NULL
    cdef void *data = NULL

    if atmos_extension == "BB":
        _atmos_extension = 1
    elif atmos_extension == "Num4D":
        _atmos_extension = 2
        if atmosphere == None:
            raise ValueError("Atmosphere data must be loaded if using numerical atmosphere extension.")
    elif atmos_extension == "Pol_BB_Burst":
        _atmos_extension = 3
        if region_extension == 'elsewhere':
            raise ValueError("'Pol_BB_Burst' is not supported by the elsewhere extension")
    elif atmos_extension == "Pol_Num2D":
        _atmos_extension = 4
        if region_extension == 'elsewhere':
            raise ValueError("'Pol_Num2D' is not supported by the elsewhere extension")
        if atmosphere == None:
            raise ValueError("Atmosphere data must be loaded if using numerical atmosphere extension.")
    elif atmos_extension == "user":
        _atmos_extension = 5
        if region_extension == 'elsewhere':
            _atmos_extension = 3
    else:
        raise ValueError("Atmosphere extension module must be 'BB', 'Num4D', 'Pol_BB_Burst', 'Pol_Num2D', or 'user'.")

    if atmosphere:
        preloaded = init_preload(atmosphere)
        data = init_ptr(numTHREADS, preloaded, _atmos_extension)
    else:
        data = init_ptr(numTHREADS, NULL, _atmos_extension)

    cdef double[::1] intensities = np.zeros(energies.shape[0],
                                            dtype = np.double)

    cdef size_t i, T
    cdef signed int ii
    for ii in prange(<signed int>energies.shape[0],
                     nogil = True,
                     schedule = 'static',
                     num_threads = <size_t> numTHREADS,
                     chunksize = 1):

        T = threadid()
        i = <size_t> ii

        if stokesQ == 1:
            intensities[i] = eval_ptr_Q(T,
                                  energies[i],
                                  mu[i],
                                  &(local_variables[i,0]),
                                  data,
                                  _beam_opt)
        else:
            intensities[i] = eval_ptr_I(T,
                                  energies[i],
                                  mu[i],
                                  &(local_variables[i,0]),
                                  data,
                                  _beam_opt)
        # get photon specific intensity
        intensities[i] *= norm_ptr() / (energies[i] * keV)

    if atmosphere:
        free_preload(preloaded)

    free_ptr(numTHREADS, data)

    return np.asarray(intensities, dtype = np.double, order = 'C')


def intensity_from_globals(double[::1] energies,
                           double[::1] mu,
                           double[::1] colatitude,
                           double[::1] azimuth,
                           double[::1] phase,
                           double[::1] global_variables,
                           double R_eq,
                           double zeta,
                           double epsilon,
                           atmosphere = None,
                           atmos_extension = "BB",
                           size_t numTHREADS = 1):
    """ Evaluate the photon specific intensity using extension modules.

    The local variables are computed from global variables and spacetime
    coordinates of events at the surface.

    :param double[::1] energies:
        A 1D :class:`numpy.ndarray` of energies in keV to evaluate photon
        specific intensity at, in the local comoving frame.

    :param double[::1] mu:
        A 1D :class:`numpy.ndarray` of angles at which to evaluate the photon
        specific intensity at. Specifically, the angle is the cosine of the
        angle of the ray direction to the local surface normal, in the local
        comoving frame.

    :param double[::1] colatitude:
        A 1D :class:`numpy.ndarray` of colatitudes in radians at which to
        evaluate photon specific intensity, in the local comoving frame. This
        is the coordinate of a fixed point relative to a distant static
        observer.

    :param double[::1] azimuth:
        A 1D :class:`numpy.ndarray` of azimuths in radians at which to
        evaluate photon specific intensity, in the local comoving frame. This
        is the coordinate of a fixed point relative to a distant static
        observer.

    :param double[::1] phase:
        A 1D :class:`numpy.ndarray` of rotational phases in radians at which to
        evaluate photon specific intensity, in the local comoving frame. The
        phase describes the rotation of the surface radiation field through
        the fixed points alluded to above. A localised hot region for instance
        rotates so that it's azimuth changes relative to a point with fixed
        azimuth described above. Specifying zero phase means that the intensity
        is evaluated at the fixed points for the initial configuration of the
        hot region (or more general radiation field).

    :param double[::1] global_variables:
        A 1D :class:`numpy.ndarray` of global variables controlling the surface
        radiation field.

    :param double R_eq:
        The equatorial radius in metres. Needed for effective gravity universal
        relation.

    :param double zeta:
        See the :class:`~xpsi.Spacetime` class property. Needed for
        effective gravity universal relation. It is recommended to supply this
        dimensionless variable by using an instance of :class:`~xpsi.Spacetime`.

    :param double epsilon:
        See the :class:`~xpsi.Spacetime` class property. Needed for
        effective gravity universal relation. It is recommended to supply this
        dimensionless variable by using an instance of :class:`~xpsi.Spacetime`.

    :param tuple atmosphere:
        A numerical atmosphere preloaded into buffers. The buffers must be
        supplied in an :math:`n`-tuple whose :math:`n^{th}` element is an
        :math:`(n-1)`-dimensional array flattened into a one-dimensional
        :class:`numpy.ndarray`. The first :math:`n-1`
        elements of the :math:`n`-tuple must each be an ordered one-dimensional
        :class:`numpy.ndarray` of parameter values for the purpose of
        multi-dimensional interpolation in the :math:`n^{th}` buffer. The
        first :math:`n-1` elements must be ordered to match the index
        arithmetic applied to the :math:`n^{th}` buffer. An example would be
        ``(logT, logg, mu, logE, buf)``, where:
        ``logT`` is a logarithm of local comoving effective temperature;
        ``logg`` is a logarithm of effective surface gravity;
        ``mu`` is the cosine of the angle from the local surface normal;
        ``logE`` is a logarithm of the photon energy; and
        ``buf`` is a one-dimensional buffer of intensities of size given by
        the product of sizes of the first :math:`n-1` tuple elements.

        It is highly recommended that buffer preloading is used, instead
        of loading from disk in the customisable radiation field extension
        module, to avoid reading from disk for every signal
        (likelihood) evaluation. This can be a non-negligible waste of compute
        resources. By preloading in Python, the memory is allocated and
        references to that memory are not in general deleted until a sampling
        script exits and the kernel stops. The likelihood callback accesses
        the same memory upon each call without I/O.

    :param str atmos_extension:
        Used to determine which atmospheric extension to use.
        Options at the moment:
        "BB": Analytical blackbody,
        "Num4D": Numerical atmosphere using 4D-interpolation from the provided
        atmosphere data,
        "user": A user-provided extension which can be set up by replacing the contents of
        the file hot_user.pyx (and elsewhere_user.pyx if needed) and re-installing X-PSI
        (if not changed, "user" is the same as "BB").

    :param int numTHREADS:
        Number of OpenMP threads to launch.

    :returns:
        A 1D :class:`numpy.ndarray` of the photon specific intensities in
        units of photons/s/keV/cm^2/sr.

    """
    cdef size_t _atmos_extension
    cdef _GEOM GEOM

    GEOM.R_eq = R_eq
    GEOM.epsilon = epsilon
    GEOM.zeta = zeta
    GEOM.star_shape_ind = 0 #assuming only oblate stars (Algendy & Morsink 2014) for global variables.

    cdef fptr_init init_ptr = init_hot
    cdef fptr_free free_ptr = free_hot
    cdef fptr_eval eval_ptr = eval_hot_I
    cdef fptr_norm norm_ptr = eval_hot_norm

    # initialise the source radiation field
    cdef _preloaded *preloaded = NULL
    cdef void *data = NULL

    if atmos_extension == "BB":
        _atmos_extension = 1
    elif atmos_extension == "Num4D":
        _atmos_extension = 2
        if atmosphere == None:
            raise ValueError("Atmosphere data must be loaded if using numerical atmosphere extension.")
    elif atmos_extension == "user":
        _atmos_extension = 3
    else:
        raise ValueError("Atmosphere extension module must be 'BB', 'Num4D', or 'user'.")

    if atmosphere:
        preloaded = init_preload(atmosphere)
        data = init_ptr(numTHREADS, preloaded, _atmos_extension)
    else:
        data = init_ptr(numTHREADS, NULL, _atmos_extension)

    cdef storage *local_vars_buf = init_local_variables(numTHREADS, NULL)

    cdef double[::1] intensities = np.zeros(energies.shape[0],
                                            dtype = np.double)

    cdef size_t i, T
    cdef int HIT
    cdef signed int ii
    for ii in prange(<signed int>energies.shape[0],
                     nogil = True,
                     schedule = 'static',
                     num_threads = <size_t> numTHREADS,
                     chunksize = 1):

        T = threadid()
        i = <size_t> ii

        # check whether finite quantity of radiation emergent from the
        # neighbourhood of the point
        HIT = HIT_or_MISS(colatitude[i],
                          azimuth[i],
                          phase[i],
                          &(global_variables[0]),
                          local_vars_buf)

        if HIT == 1:
            eval_local_variables(colatitude[i],
                                 azimuth[i],
                                 phase[i],
                                 &GEOM,
                                 &(global_variables[0]),
                                 local_vars_buf,
                                 T)

            intensities[i] = eval_ptr(T,
                                      energies[i],
                                      mu[i],
                                      local_vars_buf.local_variables[T],
                                      data,
                                      0)
        else:
            intensities[i] = 0.0

        # get photon specific intensity
        intensities[i] *= norm_ptr() / (energies[i] * keV)

    if atmosphere:
        free_preload(preloaded)

    free_local_variables(numTHREADS, local_vars_buf)
    free_ptr(numTHREADS, data)

    return np.asarray(intensities, dtype = np.double, order = 'C')


def effective_gravity(double[::1] cos_colatitude,
                      double[::1] R_eq,
                      double[::1] zeta,
                      double[::1] epsilon,
                      str star_shape = "AGM14"):
    """ Approximate local effective gravity using a universal relation.
    (or a spherical star if setting star_shape="Sphere")

    See Morsink et al. (2007), AlGendy & Morsink (2014), and Bogdanov et al.
    (2019, ApJL, 887, L26).

    :param double[::1] cos_colatitude:
        A 1D :class:`numpy.ndarray` of colatitudes, with respect to *a*
        rotational pole, at which to evaluate the local effective gravity.
        on the surface. Specifically, the angle is the cosine of the
        colatitude of the ray direction to the local surface normal, in the
        local comoving frame.

    :param double[::1] R_eq:
        A 1D :class:`numpy.ndarray` of surface coordinate equatorial radii in
        metres.

    :param double[::1] zeta:
        A 1D :class:`numpy.ndarray` of a dimensionless function of stellar
        properties. See :attr:`~xpsi.Spacetime.Spacetime.zeta`.

    :param double[::1] epsilon:
        A 1D :class:`numpy.ndarray` of a dimensionless function of stellar
        properties. See :attr:`~xpsi.Spacetime.Spacetime.epsilon`.

    :param double[::1] star_shape:
        A string defining the shape model of the star ("AGM14" or "Sphere").
        See :attr:`~xpsi.Spacetime.Spacetime.star_shape`.

    :returns:
        A 1D :class:`numpy.ndarray` of the logarithms (base 10) of the
        effective gravities in units of cm/s^2.

    """

    cdef double[::1] gravity = np.zeros(cos_colatitude.shape[0],
                                        dtype=np.double)

    cdef unsigned int i

    allowed_models = ["AGM14","Sphere"]
    cdef unsigned int star_shape_ind = allowed_models.index(star_shape)

    for i in range(<size_t>gravity.shape[0]):
        gravity[i] = effectiveGravity(cos_colatitude[i],
                                      R_eq[i],
                                      zeta[i],
                                      epsilon[i],
                                      star_shape_ind)

    return np.asarray(gravity, dtype = np.double, order = 'C')
