from xpsi.global_imports import *
from xpsi import _warning, _verbose
from xpsi.utils import  make_verbose, verbose

from os.path import join as _join

import xpsi
from xpsi.global_imports import _2pi
from xpsi.pixelmesh.integrator import integrate as _integrate
from xpsi.tools.energy_integrator import energy_integrator
from xpsi.tools.phase_integrator import phase_integrator

if _mpl is not None:
    import matplotlib
    from matplotlib import pyplot as plt
    from matplotlib.figure import Figure
    from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
    from matplotlib import rcParams
    from matplotlib.ticker import MultipleLocator, AutoLocator, AutoMinorLocator
    from matplotlib import gridspec
    from matplotlib import cm
    from matplotlib import animation
    import matplotlib.image as mgimg
else:
    pass

from xpsi.Photosphere import Photosphere

class PhotospherePlotter( Photosphere ):
    """ A photosphere plotter to do all sorts of cool images.

    """
    def __init__ (self, Photosphere) :
        self.photosphere = Photosphere
        print(self.photosphere._params)


    @property
    def global_variables(self):
        """ Get a vector of global surface radiation field variables.

        :returns: An *ndarray[n]* of scalars required to evaluate variables
                  that control the radiation field w.r.t local comoving frames
                  across the stellar surface.

        The following code block is how one would pass the properties of a
        single-temperature circular ``HotRegion`` to the extension modules. If
        you have more than one ``HotRegion`` object merged into the subspace
        associated with the ``Photosphere`` object, they may each be prefixed,
        meaning that the set of parameter names below would need to be prefixed
        at the least, and unless you only want to image one ``HotRegion``, the
        parameters of the ``HotRegions`` object are required.

        .. highlight:: python
        .. code-block:: python

            return _np.array([self['super_colatitude'],
                              self['phase_shift'] * _2pi,
                              self['super_radius'],
                              self['super_temperature']])

        The phase shift controls the initial rotational phase of the
        ``HotRegion`` when imaging commences.

        """
        try:
            return _np.array([self['temperature']])
        except KeyError:
            raise NotImplementedError('Subclass and provide an implementation.')

    @property
    def global_to_local_file(self):
        try:
            return self._global_to_local_file
        except AttributeError:
            return None

    @global_to_local_file.setter
    def global_to_local_file(self, filepath):
        if not isinstance(filepath, _six.string_types):
            raise TypeError('File path must be a string.')
        elif not _os.path.isfile(filepath):
            raise IOError('File does not exist.')
        self._global_to_local_file = filepath

    @property
    def images(self):
        """ Get the precomputed image information. """
        return self._images

    @images.setter
    def images(self, images):
        """ Store an *ndarray[i,j,k]* of images and associated information. """
        try:
            for i, obj in enumerate(images):
                if not isinstance(obj, _np.ndarray):
                    if i < len(images) - 3:
                        raise TypeError('Image information was expected to be '
                                        'contained in an ndarray.')
                    elif obj is not None and not isinstance(obj, float):
                        raise TypeError('Unexpected type for image information.')
        except TypeError:
            raise TypeError('An iterable of objects containing image '
                            'information must be supplied.')

        if len(images) != 13:
            raise ValueError('There must be six ndarray objects specifing '
                             'image information.')

        msg = 'Image information element %i must have %i dimensions.'

        # tuple elements:
        #   energy-phase resolved signal (2D array)
        #   x coordinate on image plane (1D array)
        #   y coordinate on image plane (1D array)
        #   colatitude mapped to point (x,y) on image plane (1D array)
        #   azimuth mapped to point (x,y) on image plane (1D array)
        #   radial coord mapped to point (x,y) on image plane (1D array)
        #   phase lag
        #   redshift
        #   aberrated ray angle to local surface normal
        #   elliptical image-plane radial array
        #   elliptical image-plane semi-major axis
        #   elliptical image-plane semi-minor axis
        #   energy-phase resolved specific intensity sky maps (3D array)
        # the last element is None if intensities not cached
        assert images[0].ndim == 2, msg % (0, 2)
        assert images[1].ndim == 1, msg % (1, 1)
        assert images[2].ndim == 1, msg % (2, 1)
        assert images[3].ndim == 1, msg % (3, 1)
        assert images[4].ndim == 1, msg % (4, 1)
        assert images[5].ndim == 1, msg % (5, 1)
        assert images[6].ndim == 1, msg % (6, 1)
        assert images[7].ndim == 1, msg % (7, 1)
        assert images[8].ndim == 1, msg % (8, 1)
        assert images[9].ndim == 1, msg % (9, 1)
        if images[12] is not None:
            assert images[12].ndim == 3, msg % (12, 3)

        _num_rays = len(images[1])
        for i in range(2,9):
            assert len(images[i] == _num_rays),\
                ('Ray map: array length mismatch (array at tuple index %i is '
                 'not equal in length to array at tuple index 1).' % i)

        assert int( _m.sqrt( _num_rays - 1 ) ) == len(images[9]),\
                ('Ray map: array length mismatch for image-plane radial '
                 'coordinate array (array at tuple index 9).')

        if images[12] is not None:
            assert images[12].shape[0] == images[0].shape[1],\
                    ('Intensity cache dimension 0 does not match the length of '
                     'dimension 1 of the specific flux array '
                     '(at tuple index 1), meaning the number of phases '
                     'is mismatched.')
            assert images[12].shape[2] == _num_rays,\
                    ('Intensity cache dimension 2 does not match the length of '
                     'ray map arrays (e.g., array at tuple index 1), meaning '
                     'the number of rays is mismatched.')

        self._images = images

    @images.deleter
    def images(self):
        del self._images

    def load_image_data(self, directory):
        """ Load imaging data from disk.

        :param str directory:
            Path to directory to load files from. Should contain files written
            to disk by :meth:`write_image_data`.

        """
        _d = directory

        photon_specific_flux = _np.load(_join(_d, 'photon_specific_flux.npy'))

        x_coordinate = _np.load(_join(_d, 'x_coordinate.npy'))

        y_coordinate = _np.load(_join(_d, 'y_coordinate.npy'))

        colatitude = _np.load(_join(_d, 'colatitude.npy'))

        azimuth = _np.load(_join(_d, 'azimuth.npy'))

        radial_coord = _np.load(_join(_d, 'radial_coord.npy'))

        phase_lag = _np.load(_join(_d, 'phase_lag.npy'))

        redshift = _np.load(_join(_d, 'redshift.npy'))

        abberated_angle = _np.load(_join(_d, 'abberated_angle.npy'))

        IP_radial_array = _np.load(_join(_d, 'IP_radial_array.npy'))

        IP_ellipse_axes = _np.load(_join(_d, 'IP_ellipse_axes.npy'))

        intensity = _np.load(_join(_d, 'intensity.npy'))

        self.images = [photon_specific_flux,
                       x_coordinate,
                       y_coordinate,
                       colatitude,
                       azimuth,
                       radial_coord,
                       phase_lag,
                       redshift,
                       abberated_angle,
                       IP_radial_array,
                       IP_ellipse_axes[0],
                       IP_ellipse_axes[1],
                       intensity]

    def write_image_data(self, directory):
        """ Write imaging data to disk.

        :param str directory:
            Path to directory to write to. Must exist.

        """
        _d = directory

        _np.save(_join(_d, 'photon_specific_flux.npy'), self.images[0])

        _np.save(_join(_d, 'x_coordinate.npy'), self.images[1])

        _np.save(_join(_d, 'y_coordinate.npy'), self.images[2])

        _np.save(_join(_d, 'colatitude.npy'), self.images[3])

        _np.save(_join(_d, 'azimuth.npy'), self.images[4])

        _np.save(_join(_d, 'radial_coord.npy'), self.images[5])

        _np.save(_join(_d, 'phase_lag.npy'), self.images[6])

        _np.save(_join(_d, 'redshift.npy'), self.images[7])

        _np.save(_join(_d, 'abberated_angle.npy'), self.images[8])

        _np.save(_join(_d, 'IP_radial_array.npy'), self.images[9])

        _np.save(_join(_d, 'IP_ellipse_axes.npy'),
                 _np.array(self.images[10:12], dtype=_np.double))

        _np.save(_join(_d, 'intensity.npy'), self.images[12])

    @property
    def photon_specific_flux(self):
        """ Get the photon specific flux as a function of phase and energy.

        :return: A two-dimensional :class:`numpy.ndarray`, where photon energy
                 varies with row number, and phase varies with column number.

        """
        return self._images[0]

    @property
    def photon_specific_intensity(self):
        """ Get the photon specific intensity.

        Function of phase, energy and sky direction.

        :return: A three-dimensional :class:`numpy.ndarray`, where the first
                 dimension is phase, the second dimension is photon energy,
                 and the third dimension is sky direction (flattened from
                 two-dimensional sky coordinates to one dimension).

        """
        return self._images[12]

    @make_verbose('Imaging the star', 'Star imaged')
    def image(self,
              reimage = False,
              reuse_ray_map = True,
              energies = None,
              num_phases = None,
              phases = None,
              phases_in_cycles = False,
              sqrt_num_rays = 100,
              epsabs_ray = 1.0e-12,
              epsrel_ray = 1.0e-12,
              max_steps = 100000,
              init_step = 0.1,
              image_plane_radial_increment_power = 1.0 / 2.0,
              threads = 1,
              cache_intensities = False,
              cache_energy_indices = None,
              cache_phase_indices = None,
              single_precision_intensities = True,
              plot_sky_maps = False,
              sky_map_kwargs = None,
              animate_sky_maps = False,
              free_memory = True,
              atm_ext="BB",
              animate_kwargs = None,
              **kwargs):
        r""" Image the star as a function of phase and energy.

        :param bool reimage:
            (Re)image the star. If ``False``, but the spacetime configuration
            has been updated or the photosphere parameters have been updated,
            a warning will be generated. In principle, one might want to plot
            sky maps using cached imaging information, or animate sky maps
            using images on disk, so reimaging is not forced if (non-fixed)
            parameters have been changed.

        :param bool reuse_ray_map:
            Reuse a precomputed ray map from the stellar surface to the image
            plane. If the spacetime configuration has changed (non-fixed
            parameters have changed), a cached ray map will *not* be reused. If
            the spacetime configuration is unchanged, but resolution settings
            have changed for ray tracing, pass ``False`` to adhere to the new
            resolution settings.

        :param ndarray[n] energies:
            Energies in keV to evaluate incident specific intensities at.

        :param int num_phases:
            The number of phases spanning the unit interval (zero and unity
            inclusive) to image at.

        :param ndarray[m] phases:
            Phases in *radians* or *cycles* at which to evaluate incident
            specific intensities at. If not ``None``, takes precedence over
            :obj:`num_phases`. The units need to be specified with the
            :obj:`phases_in_cycles` keyword argument: if ``False``, give
            the phase array in *radians*.

        :param bool phases_in_cycles:
            Is the phase array, if not ``None``, in units of rotational cycles?

        :param int sqrt_num_rays:
            Square-root of the number of rays. This is the level of
            discretisation in both a radial coordinate and a polar coordinate
            on an elliptical image plane.

        .. note::

            When the spacetime is static or extremely close to being static in
            a numerical context, at the resolutions we are interested in, we
            need to mitigate problems with rays that graze the pole
            infinitesimally close to the polar axis. In the vicinity of the
            polar coordinate singularity the ODE system is stiff and the
            solution is unstable. The most straightforward way to mitigate this
            is to perform a fallback forward Euler step for a ray that passes
            exactly through the pole, and use that ray as an approximation for
            the grazing ray that started very nearby on the image plane.
            Internally, if a ray intersects the image plane at
            :math:`x`-coordinate that is numerically very close to, but not
            exactly, zero (which would mean alignment to the rotational axis),
            it is approximated by a ray that intersects :math:`x=0`.
            Image-plane interpolation of quantities (such as intensity) for
            the purpose of visualisation will then smooth out any such
            artefacts.

            Moreover, as an additional measure against artefacts in the sky
            maps in the vicinity of the rotational pole, rays are distributed
            accordingingly. For example, if we request :math:`n=400` rays per
            dimension, a maximal spacing of the rays from the rotational axis
            is achieved by rotating the *spokes* of rays (by up to
            :math:`\pm\pi/n`) so that no spoke is aligned (or anti-aligned)
            with the :math:`y`-direction.

        :param float epsabs_ray:
            Absolute error tolerance per ray to adhere to during numerical
            integration.

        :param float epsrel_ray:
            Relative error tolerance per ray to adhere to during numerical
            integration.

        :param int max_steps:
            Maximum number of steps to permit per ray before forced termination
            of integration.

        :param float init_step:
            The initial *suggested* step size at the image plane for the
            affine parameter for each ray.

        :param float image_plane_radial_increment_power:
            Controls the behaviour of the radial discretisation.
            Higher values towards unity result in linear spacing of rays with
            the radial coordinate. Lower values towards zero squeeze the rays
            towards the visible stellar limb, which is necessary for resolving
            images of extended radiating elements. Values above unity are not
            recommended, and would squeeze rays towards the image-plane origin,
            compromising resolution at the stellar limb.

        :param int threads:
            Number of OpenMP threads to spawn for parallel blocks of code.
            Parallel blocks include ray integration to generate a global ray
            map from image plane to surface; and image calculation at a
            sequence of rotational phases.

        :param float cache_intensities:
            Cache the photon specific intensity sky maps in memory, as a
            function of phase and energy? The type must be a float (greater than
            or equal to zero) or ``False``. The value represents the limiting
            size in GB that can be allocated for the intensity cache. Defaults
            to zero because this dominates memory consumption. You need to
            activate this option if you want to plot the sky maps (see below).
            To activate, supply a limit. A hard limit of 2 GB is imposed for
            safety. To override, use the secret :obj:`_OVERRIDE_MEM_LIM`
            keyword argument to supply a positive limit in GB.

        :param ndarray[m] cache_phase_indices:
            A one-dimensional :class:`numpy.ndarray` of ``dtype=numpy.int32``,
            specifying the phase-array indices to cache intensities at. This
            is useful to save memory when you want to plot specific intensity
            skymaps but also compute the specific flux at many more phases. If
            ``None``, intensities will be cached at all phases subject to
            memory constraints. Note that the order of the list matters for
            plotting order, so the indices should generally increase, as should
            the phases themselves. If plotting the pulse-profile and spectrum,
            then this is a case where many more phases are useful for the
            resolving specific flux pulse-profile than are needed to plot
            specific intensity skymaps and specific flux spectra at three
            representative phases.

        :param ndarray[m] cache_energy_indices:
            A one-dimensional :class:`numpy.ndarray` of ``dtype=numpy.int32``,
            specifying the energy-array indices to cache intensities at. This
            is useful to save memory when you want to plot specific intensity
            skymaps but also compute the specific flux at many more energies. If
            ``None``, intensities will be cached at all energies subject to
            memory constraints. Note that the order of the list matters for
            plotting order, so the indices should generally increase, as should
            the energies themselves. If plotting the pulse-profile and spectrum,
            then this is a case where many more energies are useful for the
            resolving specific flux spectrum than are needed to plot specific
            intensity skymaps and specific flux pulse-profiles at three
            representative energies.

        :param bool single_precision_intensities:
            Cache the intensities in single precision? In most use cases,
            double precision is simply unnecessary, and because memory
            consumption can be high, choosing single precision can reduce
            memory requirements by a factor of two. Note that this only applies
            to the caching of intensities, not the calculation of intensities,
            which is done safely in double precision; only the final caching
            operation is a demotion cast to single precision. The default
            is single precision caching. Option ignored if intensities are not
            cached.

        :param bool plot_sky_maps:
            Plot (specific) intensity sky maps at a sequence of phases, or
            by averaging over phase. Maps can be made at one more energies
            or energy intervals. The images will be written to disk and
            can be used as frames in an animated sequence.

        :param dict sky_map_kwargs:
            Dictionary of keyword arguments passed to
            :meth:`~Photosphere._plot_sky_maps`. Refer to the associated
            method docstring for available options.

        :param bool animate_sky_maps:
            Compile images from disk into an animated sequence.

        :param bool free_memory:
            Try to free the imaging information before animating a sequence of
            sky maps written to disk, to try to avoid high memory usage. For
            safety the default is to free the memory, so deactivate this at your
            own risk. If there are other non-weak references created to the
            underlying objects, the memory may fail to be freed. In the
            methods below, the aim is that the native garbage collection cleans
            up the references because they only exist in the method local scope
            (no closures or globals).

        .. note::

            Memory used for plotting the sky maps and loading the images from
            disk to animate a phase sequence might not be straightforwardly
            freed despite efforts to do so, because of non-weak references
            covertly held by the matplotlib module.

        :param str atm_ext:
            Used to determine which atmospheric extension to use.
            Options at the moment:
            "BB": Analytical blackbody
            "Num4D": Numerical atmosphere using 4D-interpolation from the provided
            atmosphere data
            "user": A user-provided extension which can be set up by replacing the contents of 
            the file hot_user.pyx (and elsewhere_user.pyx if needed) and re-installing X-PSI
            (if not changed, "user" is the same as "BB").

        :param dict animate_kwargs:
            Dictionary of keyword arguments passed to
            :meth:`~Photosphere._animate`. Refer to the associated method
            docstring for available options.

        :param bool deactivate_all_verbosity:
            Deactivate the verbose output? Note that despite this keyword
            argument not appearing in the method signature, it is a valid
            switch.

        """
        if atm_ext=="BB":
            atmosphere_extension = 1
        elif atm_ext=="Num4D":
            atmosphere_extension = 2
        elif atm_ext=="user":
            atmosphere_extension = 3
        else:
            raise TypeError('Got an unrecognised atm_ext argument. Note that the only allowed '
                            'atmosphere options are at the moment "BB", "Num4D", and "user".')
        ref = self._spacetime # geometry shortcut saves characters
        try:
            _DV = deactivate_verbosity
        except NameError:
            _DV = False

        _exc = ValueError('You need to cache intensity sky maps if you '
                          'want to plot them.')
        try:
            self.images
        except AttributeError:
            if not reimage:
                if plot_sky_maps:
                    raise _exc
                else:
                    yield ('Warning: star will not be reimaged... assuming '
                           'images exist on disk.')
        else:
            if not reimage and plot_sky_maps and self.images[-1] is None:
                raise _exc

        if phases is not None and not isinstance(phases, _np.ndarray):
            raise TypeError('Imaging phases must be in a 1D ndarray.')
        elif isinstance(phases, _np.ndarray):
            if phases_in_cycles:
                if phases[0] != 0.0 or phases[-1] != 1.0:
                    _warning('Phase array does not span the unit interval.')
                phases *= _2pi
        elif phases is None:
            if num_phases is None or not isinstance(num_phases, int):
                raise TypeError('Integer number of phases required.')
            phases = _np.linspace(0.0, 1.0, num_phases) * _2pi

        if not isinstance(energies, _np.ndarray):
            raise TypeError('Imaging energies must be in a 1D ndarray.')

        if sky_map_kwargs is not None:
            time_is_space = sky_map_kwargs.get('time_is_space', False)

        if reimage:
            if plot_sky_maps and not cache_intensities:
                raise _exc

            if cache_intensities:
                _override_mem_lim = kwargs.get('_OVERRIDE_MEM_LIM', 1.0)
                if not isinstance(_override_mem_lim, float):
                    raise TypeError('Intensity cache limit override must be a '
                                    'float.')
                elif _override_mem_lim < 0.0:
                    raise ValueError('Intensity cache limit override must be '
                                     'positive or zero.')
                if not isinstance(cache_intensities, float):
                    raise TypeError('Intensity cache limit must be a float.')
                elif not 0.0 <= cache_intensities <= _override_mem_lim:
                    raise ValueError('Intensity cache limit must be positive '
                                     'and less than the safety limit, which '
                                     'in turn can be overridden as described '
                                     'in the method docstring.')

                if cache_energy_indices is None:
                    cache_energy_indices = _np.arange(len(energies),
                                                      dtype=_np.int32)
                elif not isinstance(cache_energy_indices, _np.ndarray):
                    raise TypeError('Energy indices for intensity caching '
                                    'must be supplied in a 1D numpy.ndarray.')
                elif cache_energy_indices.dtype != _np.int32:
                    raise TypeError('Energy indices for intensity caching '
                                    'must be integers.')
                elif time_is_space and len(cache_energy_indices) != len(energies):
                    raise TypeError('Sky maps must be cached at all energies.')

                if cache_phase_indices is None:
                    cache_phase_indices = _np.arange(len(phases),
                                                      dtype=_np.int32)
                elif not isinstance(cache_phases_indices, _np.ndarray):
                    raise TypeError('Phase indices for intensity caching '
                                    'must be supplied in a 1D numpy.ndarray.')
                elif cache_phase_indices.dtype != _np.int32:
                    raise TypeError('Phase indices for intensity caching '
                                    'must be integers.')
                elif not time_is_space and len(cache_phases_indices) != len(phases):
                    raise TypeError('Sky maps must be cached at all phases.')

                _req_size = 4.0 if single_precision_intensities else 8.0
                _req_size *= len(cache_phase_indices) * len(cache_energy_indices) # bytes
                _req_size *= sqrt_num_rays**2.0 # + 1.0 # origin ray negligible

                if _req_size/1.0e9 >= cache_intensities:
                    raise MemoryError('Too much memory would be required to '
                                      'cache the intensities at this '
                                      'resolution. Try decreasing the number '
                                      'of rays, energies, and/or phases, or '
                                      'override the cache size limit if '
                                      'safe.')
                cache_intensities = True
            else:
                cache_intensities = False

            try:
                self.images
            except AttributeError:
                if reuse_ray_map:
                    yield ('Warning: a ray map has not been cached... '
                           'tracing new ray set')
            else:
                # if spacetime configuration was updated
                if ref.needs_update or not reuse_ray_map:
                    # try to free up memory; CPython reference counting means
                    # this should have immediate effect
                    del self.images
                else:
                    # del self.images[0] # doesn't require much memory
                    del self.images[-1]  # requires far more memory

            try:
                _ray_map = tuple(self.images[1:])
                yield 'Cached ray set to be reused... commencing imaging'
            except AttributeError:
                _ray_map = None
                yield 'Commencing ray tracing and imaging'

            images = _integrate(threads,
                                ref.r_s,
                                ref.R,
                                ref.Omega,
                                self['mode_frequency'],
                                ref.zeta,
                                ref.epsilon,
                                ref.a, # dimensionless spin
                                ref.q, # mass quadrupole
                                ref.d,
                                ref.i,
                                sqrt_num_rays,
                                epsabs_ray,
                                epsrel_ray,
                                max_steps,
                                init_step,
                                image_plane_radial_increment_power,
                                self.global_variables,
                                energies,
                                phases,
                                cache_intensities,
                                cache_energy_indices,
                                cache_phase_indices,
                                single_precision_intensities,
                                _ray_map,
                                self.global_to_local_file,
                                self._hot_atmosphere,
                                atmosphere_extension)

            if images[0] == 1:
                raise Exception('A numerical error arose during imaging '
                                'computation... terminating simulation.')
            elif _ray_map is not None: # only recalculated info is returned
                # tuple elements:
                #   energy-phase resolved signal (2D array)
                #   energy-phase resolved specific intensity sky maps (3D array)
                # the last element is None if intensities not cached
                # transpose so signal phase increments along columns
                self.images[0] = images[1].T
                self.images.append(images[2])
            else: # the ray map is also returned
                # tuple elements:
                #   energy-phase resolved signal (2D array)
                #   x coordinate on image plane (1D array)
                #   y coordinate on image plane (1D array)
                #   colatitude mapped to point (x,y) on image plane (1D array)
                #   azimuth mapped to point (x,y) on image plane (1D array)
                #   radial coord mapped to point (x,y) on image plane (1D array)
                #   phase lag
                #   redshift
                #   aberrated ray angle to local surface normal
                #   elliptical image-plane radial array
                #   elliptical image-plane semi-major axis
                #   elliptical image-plane semi-minor axis
                #   energy-phase resolved specific intensity sky maps (3D array)
                # the last element is None if intensities not cached
                # transpose so signal phase increments along columns
                self.images = [images[1].T] + list(images[2:])
                yield 'Ray tracing complete.'
                yield 'Ray set cached.'

            if cache_intensities:
                yield 'Intensity caching complete.'
            else:
                if len(phases) > 1:
                    yield 'Phase-resolved specific flux integration complete.'
                else:
                    yield 'Specific flux integration complete.'

            # memoization
            self._spacetime([param.value for param in self._spacetime])

        if sky_map_kwargs is None: sky_map_kwargs = {}
        if animate_kwargs is None: animate_kwargs = {}

        if plot_sky_maps or animate_sky_maps:
            if cache_phase_indices is None:
                cache_phase_indices = _np.arange(len(phases),
                                                  dtype=_np.int32)
            if cache_energy_indices is None:
                cache_energy_indices = _np.arange(len(energies),
                                                  dtype=_np.int32)

            root_dir = sky_map_kwargs.pop('root_dir', './images')
            file_root = sky_map_kwargs.pop('file_root', 'skymap')
            file_root = _os.path.join(root_dir, file_root)

            phase_average = sky_map_kwargs.get('phase_average', False)
            if phase_average and time_is_space:
                raise ValueError('Cannot phase average sky maps when spatial '
                                 'dimensions are used to render time.')
            if phase_average and animate_sky_maps:
                raise ValueError('Phase averaged sky maps cannot be animated.')

            bolometric = sky_map_kwargs.get('bolometric', False)
            if bolometric and not time_is_space:
                raise ValueError('Cannot energy-integrate sky maps when spatial '
                                 'dimensions are used to render energy.')
            if bolometric and animate_sky_maps:
                raise ValueError('Bolometric sky maps cannot be animated.')

        if plot_sky_maps:
            if not _os.path.isdir(root_dir):
                _os.mkdir(root_dir)
            elif _os.path.isfile(file_root + '_0.png'):
                yield ('\nWarning: at least one image file exists '
                      'in ``%s``.' % root_dir)
                yield ('Attempting to move image files to a subdirectory '
                      'of ``%s``.' % root_dir)
                try: # to archive the existing image files
                    from datetime import datetime
                    obj = datetime.now()
                    temp = '__datetime__%i.%i.%i__%i.%i.%i' % (obj.day,
                                                               obj.month,
                                                               obj.year,
                                                               obj.hour,
                                                               obj.minute,
                                                               obj.second)
                    temp = _os.path.join(root_dir, 'archived_%s' % temp)
                    _os.mkdir(temp)

                    image_files = _os.listdir(root_dir)

                    for image in image_files:
                        if '.png' in image:
                            _os.rename(_os.path.join(root_dir, image),
                                       _os.path.join(temp, image))
                except Exception as e:
                    raise Exception('Aborting: image files would be '
                                    'overwritten. %s' %  str(e))
                else:
                    yield 'Image files archived in subdirectory ``%s``.' % temp

            figsize, dpi, num_frames = self._plot_sky_maps(file_root,
                                               _phases = phases,
                                               _energies = energies,
                                               _c_idxs = cache_energy_indices,
                                               _c_pidxs = cache_phase_indices,
                                               _redraw = True,
                                               deactivate_verbosity = _DV,
                                               **sky_map_kwargs)
        elif animate_sky_maps:
            if reimage:
                raise ValueError('Star was reimaged but sky maps were not '
                                 'plotted... aborting animation.')

            figsize, dpi, num_frames = self._plot_sky_maps(file_root,
                                               _phases = phases,
                                               _energies = energies,
                                               _c_idxs = cache_energy_indices,
                                               _c_pidxs = cache_phase_indices,
                                               _redraw = False,
                                               deactivate_verbosity = _DV,
                                               **sky_map_kwargs)

        if animate_sky_maps:
            if not _os.path.isfile(file_root + '_0.png'):
                raise IOError('No images located for animation.')

            if num_frames is None and reimage:
                if not time_is_space:
                    num_frames = self.images[-1].shape[0]
                else:
                    num_frames = self.images[-1].shape[1]
            elif num_frames is None:
                if not time_is_space:
                    try:
                        num_frames = len(phases)
                    except TypeError:
                        raise TypeError('You need to declare the image phases '
                                        'in order to include all images from disk.')
                else:
                    try:
                        num_frames = len(energies)
                    except TypeError:
                        raise TypeError('You need to declare the image energies '
                                        'in order to include all images from disk.')
            if free_memory:
                try:
                    del self.images # try to free up memory
                except AttributeError:
                    pass

            self._animate(file_root, num_frames,
                          figsize, dpi,
                          deactivate_verbosity = _DV,
                          **animate_kwargs)

        yield None

    @make_verbose('Plotting intensity sky maps', 'Intensity sky maps plotted')
    def _plot_sky_maps(self,
                       _file_root,
                       _phases,
                       _energies,
                       _c_idxs,
                       _c_pidxs,
                       _redraw,
                       threads = 1,
                       with_pulse_profile_and_spectrum = False,
                       time_is_space = False,
                       panel_layout = None,
                       panel_indices = None,
                       cycles = 1,
                       phase_average = False,
                       bolometric = False,
                       energy_bounds = None,
                       phase_bounds = None,
                       num_levels = 100,
                       add_zero_intensity_level = True,
                       normalise_each_panel = True,
                       invert = False,
                       annotate_energies=False,
                       annotate_phases=False,
                       energy_annotation_format='[%.1f keV]',
                       phase_annotation_format='[%.1f cycles]',
                       annotate_location=(0.05,0.05),
                       colormap = None,
                       figsize = (10,10),
                       usetex = False,
                       fontsize_scale = 1.0,
                       tick_spacing = (0.2,1.0),
                       tick_length_scaling = 1.0,
                       dpi_scale = 1.0,
                       **kwargs):
        r""" Helper method for specific intensity sky map visualization.

        Uses Delaunay triangulation to create an irregular sky mesh and
        calculate photon (specific) intensity contours at a sequence of phases.

        Each figure generated contains a sequence of panels arranged in
        one or two spatial dimensions. Each panel is an intensity sky map,
        either at a particular energy or integrated over a finite energy
        interval. Panels cannot mix specific and integrated intensities. Only
        a some sequence of energy (intervals), in any order, can be identified
        as labelling panels in one or two spatial dimensions. Time (rotational
        phase), whilst it could be defined to label a sequence of panels, is
        only identified as a labelling a sequence of *figures*.

        Similarly, sequence of energies could be identified as a variable
        labelling a sequence of figures, but is not. Moreover, energy and time
        could label panels in to spatial dimensions, but such mixing is not
        permitted. Finally, variables controlling the source-receiver system
        could be identified as labels of panels and/or figures, but this
        is also not supported. More complicated rendering patterns may be
        supported in future versions, but for now can be achieved via
        custom extensions building off of the current functionality.

        :param str _file_root:
            Relative or absolute path to parent directory for images,
            extended with the root name for the image files. E.g., the
            default is ``./images/skymap``. You do not need to change this
            unless you wish to, and is otherwise reserved for internal use.
            You may supply a custom file path via keywords ``root_dir`` and
            ``file_root`` upon calling :meth:`~Photosphere.image`, which are
            concatenated appropriately.

        :param ndarray[n] _phases:
            The phases at which the star was imaged. This is handled
            internally, so do *not* pass a keyword argument. If phase averaging,
            the minimum and maximum phases must be zero and :math:`2\pi`
            radians (i.e., zero and one cycles).

        :param ndarray[n] _energies:
            The energies at which the star was imaged. This is handled
            internally, so do *not* pass a keyword argument.

        :param ndarray[n] _c_idxs:
            The energy indices for which the intensity maps were cached for
            memory-efficieny plotting. This is handled internally, so do *not*
            pass a keyword argument.

        :param ndarray[n] _c_pidxs:
            The energy indices for which the intensity maps were cached for
            memory-efficieny plotting. This is handled internally, so do *not*
            pass a keyword argument.

        :param bool _redraw:
            Redraw the sky maps? This is handled internally, so do *not* pass
            a keyword argument.

        :param int threads:
            Number of OpenMP threads to spawn.

        :param bool with_pulse_profile_and_spectrum:
            A setting that fundamentally changes some behaviours. If
            deactivated (the default), only photon (specific) intensity skymaps
            are plotted. The following frame is an example:

        .. image:: _static/_skymap_plot.png

        :param bool with_pulse_profile_and_spectrum:
            If deactivated, the following keyword arguments do not have a use:
            :obj:`cycles` and :obj:`colormap`. If *activated*, photon specific
            intensity skymaps at three energies are plotted in each frame,
            together with their associated photon specific flux pulse-profiles,
            and also the photon specific flux spectrum at a finer array of
            energies. Use :obj:`panel_indices` to select the energies. The
            pulse-profiles are each normalised to their respective maxima,
            and the spectrum shows the relative orders of magnitude of the
            specific flux signals. The following frame is an example:

        .. image:: _static/_skymap_with_pulse_profile_and_spectrum_plot.png

        :param bool with_pulse_profile_and_spectrum:
            If activated, a subset of other keyword arguments are ignored:
            :obj:`energy_bounds`, :obj:`phase_average`, :obj:`panel_layout`,
            and :obj:`invert`. The panel layout is rigid (not customisable) in
            order to focus on the plot quality. If more energies were added, the
            information density in the plot-space might become too high without
            adding much more new information.

        :param bool time_is_space:
            Each image is at constant energy (or is a spectral trace up to
            that energy) instead of being at constant phase (or a pulse-profile
            trace up to that phase).

        :param tuple[int,int] panel_layout:
            Two elements: the number of rows and columns of panels. If ``None``,
            a layout is automatically determined based on the number of
            images to be plotted.

        :param iterable panel_indices:
            These ordered integers will be used to select intensity information
            by indexing the energy dimension of a 3D intensity array. If
            specific intensites are plotted, these integers should index
            a subset of energies at which the star was imaged. If intensities
            are plotted, these integers should index a subset of the energy
            intervals over which specific intensities are integrated. See
            the :obj:`energy_bounds` keyword argument. If the flux is
            calculated at more energies than specific intensities are cached
            at, then these integers need to index the
            :obj:`cache_energy_indices` array appropriately.

        :param int cycles:
            Nuber of cycles to generate images for. Only relevant if one cycle
            is different to the next in terms of the frame, e.g., most commonly
            if plotting the pulse-profile traces over more than one cycle. If
            frames separated by one cycle are identical, declare the number
            of cycles to the animator instead.

        :param bool phase_average:
            Average each sky map over one revolution (cycle) in phase?
            Note that the resulting image is incompatible with the currently
            supported animation mode. The following image is an example:

        .. image:: _static/_skymap_phaseaveraged.png

        :param bool bolometric:
            Integrate each sky map over energy? Note that the resulting image
            is incompatible with the currently supported animation mode.

        :param iterable energy_bounds:
            A set of two-element containers. Each container has an ordered pair
            of energies which delimit an integration domain. Specific intensity
            is integrated along each sky direction, at each phase, between
            these energy bounds. The bounds must be between the minimum and
            maximum energies at which the star was imaged. If ``None``,
            specific intensity sky maps will be plotted (the default). This
            option is ignored if energy is defined as the time dimension,
            meaning the sky maps are animated with respect to energy.

        :param iterable phase_bounds:
            A set of two-element containers. Each container has an ordered pair
            of energies which delimit an integration domain. Specific intensity
            is integrated along each sky direction, at each phase, between
            these energy bounds. The bounds must be between the minimum and
            maximum energies at which the star was imaged. If ``None``,
            specific intensity sky maps will be plotted (the default). This
            option is ignored if phase is defined as the time dimension,
            meaning the sky maps are animated with respect to phase.

        .. note::

            To use this functionality, the specific intensities must have
            been cached at all energies the specific flux is calculated at.

        :param int num_levels:
            Number of contour levels in (specific) intensity, distributed
            between minimum finite, and maximum values per panel, or over all
            panels. See :obj:`normalise_each_panel` keyword argument.

        :param bool add_zero_intensity_level:
            Add a contour level at zero intensity such that the colormap
            minimum corresponds to zero intensity? If ``True`` (the default),
            then the background sky, where there is by definition zero model
            intensity, has the same colour only as the subset of the image of
            the surface that is not radiating in the model. The disadvantage of
            this choice is that the intensity structure of the image as a
            function of phase and sky direction is generally not as well-
            resolved by the colour and greyscale variation. In the limit
            that the minimum finite intensity of the image is far smaller than
            the maximum, then the intensity resolution by colour and greyscale
            values is highest. If ``False``, then the minimum colour is
            assigned to the minimum finite intensity as a function of phase
            and sky direction. This also maximally resolves the intensity by
            colour and greyscale values, which is useful for models wherein the
            surface radiation field is constructed, for instance, from
            uniform-temperature localised hot regions. However, in this case
            the background sky colour is undefined; the background sky colour
            is thus set to the minimum colour in the colormap, meaning that the
            fainest subset of the image over phase and sky direction merges
            with the background sky in terms of colour and greyscale values.

        :param bool normalise_each_panel:
            Normalise the contour colormap to each skymap panel uniquely, or
            globally over all panels? The former yields relative intensity
            as function of phase and sky direction for an energy or energy
            interval, whilst the latter offers more spectral information but
            emission in some panels may not be discernable.

        :param bool invert:
            Invert the greyscale to show bright pixels as dark on a white
            plot background. If a colormap is manually supplied, this just
            controls the plot background colour. Inversion is recommended
            for printed format, whilst a black background is more intuitive
            when in digital format.

        :param obj colormap:
            Usage dependent on other settings. If not plotting the pulse-profile
            and spectrum, then this is simply a (matplotlib) colormap object.
            Choose something appropriate and *accessible* for a non-negative
            scalar field (sky intensities). If plotting the pulse-profile and
            spectrum too, then :obj:`colormap` can be the string
            ``'RedGreyBlue'`` to invoke the default colour scheme which is reds
            for the lowest energy intensity skymap and pulse-profile; pure
            greyscale for the intermediate energy; and blues for the highest
            energy.  Alternatively, if :obj:`colormap` is simply ``None``, the
            default greyscale will be used for all energies, with all
            pulse-profiles in black. Lastly, you can supply a three-element
            list or tuple of colormap objects, ordered from lowest to highest
            energy; the pulse-profile line colours will be retrieved as the
            midpoint of the colormap. Note that the background sky colour will
            be set to the lowest colour in each colourmap.

        :param tuple(int,int) figsize:
            The figure size (width, height) in *inches*. If the dimensions are
            inconsistent with the aspect ratio suggested by the
            :obj`panel_layout` settings, the height of the figure will be
            automatically rescaled to achieve congruence, meaning each panel is
            approximately square.

        :param bool usetex:
            Use TeX backend for figure text.

        :param float fontsize_scale:
            Use this argument to scale the font size of figure text relative
            to the default font size that is automatically determined based
            on the approximate panel size, which is in turn based on the
            figure size and the panel layout.

        :param tuple[float,float] tick_spacing:
            A two-element container. The first element is the minor tick
            spacing, and the second is the major tick spacing, for both
            the :math:`x` and :math:`y` directions on the image plane. The
            units are both the maximum possible angular size of the image of the
            surface in an ambient Schwarzschild spacetime,
            :math:`R_{\\rm eq}/\\sqrt{1-r_{\\rm s}/R_{\\rm eq}}`.

        :param float tick_length_scaling:
            Use this argument to scale the axis tick lengths relative to the
            default lengths that are automatically determined based on the
            panel size.

        :param float dpi_scale:
            Use this argument to scale the dots per inch of the figure,
            relative to the default that is automatically determined based
            on the panel size.

        """
        file_root = _file_root
        phases = _phases
        energies = _energies
        redraw = _redraw

        if with_pulse_profile_and_spectrum:
            if len(panel_indices) != 3:
                raise ValueError('Selected plot type designed for showcasing '
                                 'the specific photon intensity skymaps and '
                                 'their associated specific photon flux '
                                 'pulse-profiles or spectra specifically at '
                                 'at three energies or phases to avoid '
                                 'excessive information density.')
            panel_layout = (2,3)

            if not isinstance(cycles, int):
                raise TypeError('Declare the number of cycles with an integer.')
            elif cycles < 1:
                cycles = 1 # quietly ignore input
                if not time_is_space:
                    num_frames = len(phases)
                else:
                    num_frames = len(energies)
            elif cycles > 1:
                if not time_is_space:
                    num_frames = len(phases) + (cycles - 1) * (len(phases) - 1)
                else:
                    num_frames = len(energies)
        elif panel_layout is None:
            x = int(_m.ceil(_m.sqrt(len(panel_indices))))
            if x * (x - 1) >= len(panel_indices):
                panel_layout = (x, x - 1)
            else:
                panel_layout = (x, x)
            if not time_is_space:
                num_frames = len(phases)
            else:
                num_frames = len(energies)

        # try to improve the aspect ratio so that each panel is
        # approximately square
        width = panel_layout[1] + (panel_layout[1] - 1) * 0.2
        _hspace = 0.25 if with_pulse_profile_and_spectrum else 0.2
        height = panel_layout[0] + (panel_layout[0] - 1) * _hspace
        aspect_ratio = height/float(width)
        if aspect_ratio != 1.0:
            if _np.abs(figsize[1]/float(figsize[0]) - aspect_ratio)/aspect_ratio > 0.1:
               figsize = (figsize[0], figsize[0] * aspect_ratio)

        # calculate an appropriate dpi to resolve each panel adequately
        dpi = (max(panel_layout) / 2.0) * 150.0 * dpi_scale

        if redraw:
            rcParams['text.usetex'] = usetex

            try:
                iter(panel_indices)
            except TypeError:
                raise TypeError('Panel indices object must be iterable.')

            if _np.product(panel_layout) < len(panel_indices):
                raise ValueError('There are too few panels for the requested '
                                 'number of intensity sky maps.')

            # some scaling for appropriate fontsize
            panel_size = max(figsize[1]/float(panel_layout[0]),
                             figsize[0]/float(panel_layout[1]))
            fontsize = (panel_size/5.0) * 14.0 * fontsize_scale
            rcParams['font.size'] = int(fontsize)

            tick_length = int((panel_size/5.0) * 8 * tick_length_scaling)

            # get coordinates of irregular set of points for triangulation
            X = self.images[1]
            Y = self.images[2]

            #if not isinstance(energies, _np.ndarray):
            #    raise TypeError('Imaging energies must be in an ndarray.')
            #if not isinstance(phases, _np.ndarray):
            #    raise TypeError('Imaging phases must be in an ndarray.')

            images = self.images[-1]
            if time_is_space: # transpose dimensions
                _images_ = _np.zeros((images.shape[1],
                                     images.shape[0],
                                     images.shape[-1]), dtype=_np.double)
                for i in range(images.shape[-1]):
                    _images_[:,:,i] = images[:,:,i].T

                images = _images_

            if with_pulse_profile_and_spectrum:
                if not time_is_space:
                    flux = self.images[0]
                else:
                    flux = self.images[0].T

            if not with_pulse_profile_and_spectrum and not time_is_space and energy_bounds:
                with verbose(True,
                         'Integrating specific intensity over energy intervals',
                         'Integrated specific intensity over energy intervals'):
                    for bounds in energy_bounds:
                        if bounds[0] > bounds[1]:
                            raise ValueError('Energy bounds in a tuple must be '
                                             'ordered.')
                        for bound in bounds:
                            if not energies[0] <= bound <= energies[-1]:
                                raise ValueError('Extrapolation would be required.')

                    if len(panel_indices) < len(energy_bounds):
                        yield 'Warning: fewer panels than energy intervals.'

                    integrated = _np.zeros((images.shape[0],
                                            len(energy_bounds),
                                            images.shape[2]), dtype=_np.double)

                    intensities = _np.zeros((images.shape[1],
                                             images.shape[2]), dtype=_np.double)

                    for i in range(images.shape[0]): # phases
                        intensities[...] = images[i,...] # sky directions

                        for k in range(len(energy_bounds)):
                            bounds = _np.log10( _np.array(energy_bounds[k]) )
                            _integrated = energy_integrator(threads,
                                                            intensities,
                                                            _np.log10(energies),
                                                            bounds)

                            integrated[i,k,:] = _integrated[0,:]

                    images = integrated
            elif not with_pulse_profile_and_spectrum and not time_is_space:
                if len(panel_indices) != len(_c_idxs):
                    yield ('Warning: number of panels not equal to number of '
                            'phases.')
            elif not with_pulse_profile_and_spectrum and time_is_space:
                if len(panel_indices) != len(_c_pidxs):
                    yield ('Warning: number of panels not equal to number of '
                           'phases.')

            if not with_pulse_profile_and_spectrum and not time_is_space and phase_average:
                with verbose(True,
                         'Averaging (specific) intensity over rotational phase',
                         'Averaged (specific) intensity over rotational phase'):
                    if phases[0] != 0.0 or phases[-1] != _2pi:
                        raise ValueError('Minimum and maximum phases at which '
                                         'star is imaged must be zero and unity '
                                         'if you are phase averaging.')

                    averaged = _np.zeros((1,
                                          images.shape[1],
                                          images.shape[2]), dtype = _np.double)

                    intensities = _np.zeros((images.shape[2],
                                             images.shape[0]), dtype = _np.double)

                    for i in range(images.shape[1]): # energies
                        for j in range(images.shape[2]): # sky directions
                            intensities[j,:] = images[:,i,j]

                        _averaged = phase_integrator(1.0, # exposure time
                                                     _np.array([0.0, 1.0]),
                                                     intensities,
                                                     phases / _2pi,
                                                     0.0) # phase shift

                        for j in range(images.shape[2]):
                            averaged[:,i,j] = _averaged[j,:]

                    images = averaged

            if not with_pulse_profile_and_spectrum and time_is_space and phase_bounds:
                with verbose(True,
                         'Integrating specific intensity over phase intervals',
                         'Integrated specific intensity over phase intervals'):
                    for bounds in phase_bounds:
                        if bounds[0] > bounds[1]:
                            raise ValueError('Phase bounds in a tuple must be '
                                             'ordered.')
                        for bound in bounds:
                            if not phases[0] <= bound <= phases[-1]:
                                raise ValueError('Extrapolation would be required.')

                    if len(panel_indices) < len(phase_bounds):
                        yield 'Warning: fewer panels than phase intervals.'

                    integrated = _np.zeros((images.shape[0],
                                            len(phase_bounds),
                                            images.shape[2]), dtype=_np.double)

                    intensities = _np.zeros((images.shape[2],
                                             images.shape[1]), dtype=_np.double)

                    for i in range(images.shape[0]): # energies
                        intensities[...] = images[i,...].T # sky directions

                        for k in range(len(phase_bounds)):
                            bounds = _np.array(phase_bounds[k])
                            _integrated = phase_integrator(1.0,
                                                           bounds,
                                                           intensities,
                                                           phases / _2pi,
                                                           0.0)

                            integrated[i,k,:] = _integrated[:,0]

                    images = integrated

            if not with_pulse_profile_and_spectrum and time_is_space and bolometric:
                with verbose(True,
                         'Integrating bolometric intensity',
                         'Averaged bolometric intensity'):
                    if phases[0] != 0.0 or phases[-1] != _2pi:
                        raise ValueError('Minimum and maximum phases at which '
                                         'star is imaged must be zero and unity '
                                         'if you are phase averaging.')

                    integrated = _np.zeros((images.shape[0],
                                            1,
                                            images.shape[2]), dtype = _np.double)

                    intensities = _np.zeros((images.shape[1],
                                             images.shape[2]), dtype = _np.double)

                    for i in range(images.shape[0]): # phases
                        for j in range(images.shape[2]): # sky directions
                            intensities[:,j] = images[i,:,j]

                        bounds = _np.log10( _np.array([energies[0], energies[-1]]) )
                        _integrated = energy_integrator(threads,
                                                        intensities,
                                                        _np.log10(energies),
                                                        bounds)

                        for j in range(images.shape[2]):
                            integrated[i,:,j] = _integrated[:,j]

                    images = integrated

            if normalise_each_panel:
                with verbose(True,
                             'Normalising each sky map panel separately',
                             'Normalised sky map panels separately'):
                    # normalise intensity for each individual panel
                    levels = []

                    if not time_is_space:
                        for j in range(images.shape[1]): # at each energy
                            # find extreme intensities over discrete set of image
                            # phases and sky directions
                            MIN = _np.min(images[:,j,:][images[:,j,:] > 0.0])
                            MAX = _np.max(images[:,j,:])
                            levels.append(_np.linspace(MIN, MAX, num_levels))

                            if add_zero_intensity_level:
                                levels[-1] = _np.array([0.0, 0.001*MIN] + list(levels[-1]))
                    else:
                        for j in range(images.shape[0]): # at each energy
                            # find extreme intensities over discrete set of image
                            # phases and sky directions
                            MIN = _np.min(images[j,:,:][images[j,:,:] > 0.0])
                            MAX = _np.max(images[j,:,:])
                            levels.append(_np.linspace(MIN, MAX, num_levels))

                            if add_zero_intensity_level:
                                levels[-1] = _np.array([0.0, 0.001*MIN] + list(levels[-1]))
            else:
                with verbose(True,
                             'Normalising sky map panels globally',
                             'Normalised sky map panels globally'):
                    MIN = _np.min(images[:,:,:][images[:,:,:] > 0.0])
                    MAX = _np.max(images[:,:,:])
                    levels = _np.linspace(MIN, MAX, num_levels)

                    if add_zero_intensity_level:
                        levels = _np.array([0.0, 0.001*MIN] + list(levels))

            # because of default tick formatting and a minus sign,
            # the left and bottom margins need to be different
            left = 0.09 * (fontsize/14.0)
            bottom = 0.11 * (fontsize/14.0)
            right = 0.975
            top = bottom + (right - left)

            ref = self._spacetime

            fig = Figure(figsize = figsize)
            canvas = FigureCanvas(fig)

            if with_pulse_profile_and_spectrum:
                if not time_is_space and colormap == 'RedGreyBlue':
                    cmap = [cm.Reds_r, cm.Greys_r, cm.Blues_r]
                    _line_colors = [cm.Reds_r(0.25),
                                    cm.Greys_r(0.0),
                                    cm.Blues_r(0.25)]
                elif not isinstance(colormap, (list, tuple)):
                    cmap = [cm.Greys_r] * 3
                    _line_colors = [cm.Greys_r(0.0)] * 3
                else:
                    cmap = colormap
                    _line_colors = [cmap[j](0.5) for j in range(len(cmap))]

                gs = gridspec.GridSpec(panel_layout[0],
                                       panel_layout[1],
                                       left=left, right=right,
                                       bottom=bottom, top=top,
                                       wspace=0.2, hspace=_hspace)

                axes = [fig.add_subplot(gs[j]) for j in range(len(panel_indices))]
                pp_ax = fig.add_subplot(gs[len(panel_indices):-1])
                spec_ax = fig.add_subplot(gs[-1])

                line_styles = ['-', '--', '-.']
            else:
                cmap = colormap or (cm.Greys if invert else cm.Greys_r)
                gs = gridspec.GridSpec(panel_layout[0],
                                       panel_layout[1],
                                       left=left, right=right,
                                       bottom=bottom, top=top,
                                       wspace=0.2, hspace=_hspace)

                axes = [fig.add_subplot(gs[j]) for j in range(len(panel_indices))]

            if with_pulse_profile_and_spectrum:
                if not time_is_space: annotate_energies = True
                else: annotate_phases = True

            _I = 10
            for i in range(images.shape[0]):
                if phase_average:
                    yield 'Rendering phase-averaged images'
                elif bolometric:
                    yield 'Rendering bolometric images'
                elif i == 0 and images.shape[0] < 10:
                    yield 'Rendering images'
                elif i == 0 and images.shape[0] >= 10:
                    yield 'Rendering image numbers [%i, %i]'%(i+1, i+_I)
                elif i%_I == 0:
                    yield 'Rendering image numbers (%i, %i]'%(i, i+_I)

                for j, idx in enumerate(panel_indices):
                    ax = axes[j]
                    if with_pulse_profile_and_spectrum:
                        _cmap = cmap[j]
                    else:
                        _cmap = cmap
                    if (with_pulse_profile_and_spectrum or
                        _np.product(panel_layout) - j - 1 < panel_layout[1]):
                        if ref.R < 1.5 * ref.r_s:
                            ax.set_xlabel(r'$(2x/(3\sqrt{3}r_{\rm s}))$')
                        else:
                            ax.set_xlabel(r'$(x/R_{\rm eq})\sqrt{1-r_{\rm s}/R_{\rm eq}}$')

                    if j % panel_layout[1] == 0:
                        if ref.R < 1.5 * ref.r_s:
                            ax.set_ylabel(r'$(2y/(3\sqrt{3}r_{\rm s}))$')
                        else:
                            ax.set_ylabel(r'$(y/R_{\rm eq})\sqrt{1-r_{\rm s}/R_{\rm eq}}$')

                    _veneer(tick_spacing, tick_spacing, ax,
                            length = tick_length)
                    if with_pulse_profile_and_spectrum:
                        ax.set_facecolor(_cmap(0.0))
                    else:
                        ax.set_facecolor('white' if invert else 'black')

                    if not time_is_space:
                        lvls = levels if isinstance(levels, _np.ndarray) else levels[idx]
                    else:
                        lvls = levels if isinstance(levels, _np.ndarray) else levels[i]

                    ax.tricontourf(X,
                                   Y,
                                   images[i,idx,:],
                                   cmap = _cmap,
                                   levels = lvls)

                    # correct the aspect ratio
                    x_view = ax.xaxis.get_view_interval()
                    diff = x_view[1] - x_view[0]
                    ax.xaxis.set_view_interval(x_view[0] - diff * 0.025,
                                               x_view[1] + diff * 0.025)
                    y_view = ax.yaxis.get_view_interval()
                    ax.yaxis.set_view_interval(y_view[1] - diff * 1.025,
                                               y_view[1] + diff * 0.025)

                    if not time_is_space and annotate_energies:
                        ax.text(annotate_location[0], annotate_location[1],
                           s=energy_annotation_format % energies[_c_idxs[idx]],
                           fontdict={'color': 'black' if invert else 'white'},
                           transform=ax.transAxes)

                    if time_is_space and annotate_phases:
                        ax.text(annotate_location[0], annotate_location[1],
                           s=phase_annotation_format % (phases[_c_pidxs[idx]]/_2pi),
                           fontdict={'color': 'black' if invert else 'white'},
                           transform=ax.transAxes)

                if with_pulse_profile_and_spectrum: # plot summaries
                    _upper_view_lim = _np.max(flux) * 100.0
                    _lower_view_lim = _np.min(flux)
                    _view_lim_diff = _np.log10(_upper_view_lim)
                    _view_lim_diff -= _np.log10(_lower_view_lim)
                    for j, idx in enumerate(panel_indices):
                        _idx = _c_idxs[idx]
                        _diff = _np.log10(flux[_idx,i])
                        _diff -= _np.log10(_lower_view_lim)
                        spec_ax.axvline(energies[_idx],
                                 0.0,
                                 _diff/_view_lim_diff,
                                 color=_line_colors[j],
                                 ls=line_styles[j], lw=1.0)

                    spec_ax.set_xscale('log')
                    spec_ax.set_yscale('log')
                    spec_ax.set_xlim(energies[0], energies[-1])
                    spec_ax.set_ylim(_lower_view_lim, _upper_view_lim)
                    _veneer(None, None, spec_ax, length = tick_length,
                            log=(True, True))
                    spec_ax.set_xlabel('Photon energy [keV]')

                    if not time_is_space:
                        spec_ax.plot(energies, flux[:,i], 'k-', lw=1.0)
                    else:
                        for j, idx in enumerate(panel_indices):
                            _idx = _c_pidxs[idx]
                            spec_ax.plot(energies, flux[_idx,:i],
                                         color=_line_colors[j],
                                         ls=line_styles[j], lw=1.0,
                                         label=phase_annotation_format % phases[_idx])

                    if not time_is_space:
                        for _cycle in range(cycles): # plot pulse-profiles
                            _ext_phases = []
                            for _i in range(_cycle):
                                _ext_phases += list(phases[1 if _i > 0  else 0:] + _i * _2pi)
                            _ext_phases += list(phases[1 if _cycle > 0 else 0:i+1] + _cycle * _2pi)
                            _ext_phases = _np.array(_ext_phases)
                            for j, idx in enumerate(panel_indices):
                                _idx = _c_idxs[idx]
                                _max = _np.max(flux[_idx,:])
                                _ext_flux = []
                                for _i in range(_cycle):
                                    _ext_flux += list(flux[_idx, 1 if _i > 0 else 0:])
                                _ext_flux += list(flux[_idx, 1 if _cycle > 0 else 0:i+1])
                                _ext_flux = _np.array(_ext_flux)
                                pp_ax.plot(_ext_phases/_2pi,
                                       _ext_flux/_max,
                                       color=_line_colors[j],
                                       ls=line_styles[j], lw=1.0,
                                       label=energy_annotation_format % energies[_idx])

                            pp_ax.set_xlim(0.0, float(cycles))
                            pp_ax.set_ylim(0.0,1.2)
                            _veneer((0.1,0.5), (0.05,0.2), pp_ax, length = tick_length)
                            pp_ax.legend(loc='upper center', ncol=3, mode='expand',
                                         handlelength=4.0,
                                         frameon=False, fancybox=False)
                            pp_ax.set_xlabel(r'Phase [$2\pi$ radians]')
                            pp_ax.set_ylabel(r'photons/cm$^2$/s/keV')

                            fig.savefig(file_root + '_%i.png' % (len(_ext_phases) - 1),
                                        dpi=dpi)
                            pp_ax.clear()
                        spec_ax.clear()
                    else:
                        pp_ax.plot(phases/_2pi,
                                 flux[i,:]/_np.max(flux[i,:]), 'k-', lw=1.0)

                        pp_ax.set_xlim(0.0, 1.0)
                        pp_ax.set_ylim(0.0,1.2)
                        _veneer((0.1,0.5), (0.05,0.2), pp_ax, length = tick_length)
                        pp_ax.legend(loc='upper center', ncol=3, mode='expand',
                                     handlelength=4.0,
                                     frameon=False, fancybox=False)
                        pp_ax.set_xlabel(r'Phase [$2\pi$ radians]')
                        pp_ax.set_ylabel(r'photons/cm$^2$/s/keV')

                        fig.savefig(file_root + '_%i.png' % i, dpi=dpi)
                else:
                    fig.savefig(file_root + '_%i.png' % i, dpi=dpi)

                for ax in axes:
                    ax.clear()

            for ax in axes:
                ax.cla()
            plt.close(fig)

        yield figsize, dpi, num_frames

    @staticmethod
    @make_verbose('Animating intensity sky maps', 'Intensity sky maps animated')
    def _animate(_file_root, _num_frames, _figsize, _dpi,
                 cycles = 1, fps = None, **kwargs):
        """ Helper method to animate the intensity sky maps.

        :param str _file_root:
            Reserved for internal use.

        :param int _num_frames:
            Reserved for internal use.

        :param int _figsize:
            Reserved for internal use.

        :param int _dpi:
            Reserved for internal use.

        :param int cycles:
            Number of explicit cycles to generate frames for. The frames from
            the principal cycle are reused, so the images are periodic in
            their phase evolution. There is no delay between cycles in this
            type of repitition.

        :param int fps:
            Frames per second. If ``None``, then one cycle (assuming images
            have been precomputed for a complete cycle), consisting of
            so many frames, will exhibit a period of one second.

        :param bool repeat:
            Inform *ffmpeg* to enter a loop when video play back commences.
            Deprecated.

        :param repeat_delay:
            Delay between repeats in milliseconds. Deprecated.

        :param str ffmpeg_path:
            Absolute path to *ffmpeg* executable. If ``None``, defaults
            to matplotlib rcParams settings, but no guarantee that the package
            will be found even if installed on system. Deprecated.

        """
        file_root = _file_root
        num_frames = _num_frames
        figsize = _figsize
        dpi = _dpi

        fig = plt.figure(figsize = figsize)
        ax = plt.subplot(111)
        plt.axis('off')
        fig.subplots_adjust(left=0, bottom=0, right=1, top=1,
                            wspace=None, hspace=None)

        # animation code based on:
        # http://jakevdp.github.io/blog/2012/08/18/matplotlib-animation-tutorial/
        filename = file_root + '_0.png'
        img = mgimg.imread(filename)
        imgplot = ax.imshow(img, aspect='auto')

        cycles = int(cycles)

        if cycles < 1:
            cycles = 1
        elif cycles > 1:
            num_frames += (cycles - 1) * (_num_frames - 1)

        class _context: # mutable nonlocal namespace in closure
            _cycle_idx = 0 # track cycle
            _j = -1        # track image index
            _last = -1

        def _update(i): # load phase-ordered set of images
            if _context._last == i:
                return imgplot,

            if _context._cycle_idx == 0 and i == _num_frames:
                _context._j = 1
                _context._cycle_idx = 1
            elif i == _num_frames + _context._cycle_idx * (_num_frames - 1):
                _context._j = 1
                _context._cycle_idx += 1
            else:
                _context._j += 1

            _context._last = i

            filename = file_root + '_%i.png' % _context._j
            img = mgimg.imread(filename)
            imgplot.set_data(img)
            return imgplot,

        ani = animation.FuncAnimation(fig, _update, frames=num_frames, blit=True)

        if fps is None:
            fps = num_frames # all frames span one second

        # secret keyword argument; not clear whether should be exposed to user
        bitrate = kwargs.get('bitrate', -1) # default is let _mpl choose

        filename = file_root + '_animated.mp4'
        yield 'Writing to disk: %s' % filename
        ani.save(filename,
                 dpi = dpi, fps = fps, bitrate = bitrate,
                 )

        fig.clf() # this or ax.cla() needed to free memory
        plt.close(fig)

        yield None

Photosphere._update_doc()

def _veneer(x, y, axes, lw=1.0, length=8, log=(False, False)):
    """ Make the plots a little more aesthetically pleasing. """
    if x is not None:
        if x[1] is not None:
            axes.xaxis.set_major_locator(MultipleLocator(x[1]))
        if x[0] is not None:
            axes.xaxis.set_minor_locator(MultipleLocator(x[0]))
    elif not log[0]:
        axes.xaxis.set_major_locator(AutoLocator())
        axes.xaxis.set_minor_locator(AutoMinorLocator())

    if y is not None:
        if y[1] is not None:
            axes.yaxis.set_major_locator(MultipleLocator(y[1]))
        if y[0] is not None:
            axes.yaxis.set_minor_locator(MultipleLocator(y[0]))
    elif not log[1]:
        axes.yaxis.set_major_locator(AutoLocator())
        axes.yaxis.set_minor_locator(AutoMinorLocator())

    axes.tick_params(which='major', colors='black', length=length, width=lw)
    axes.tick_params(which='minor', colors='black', length=int(length/2), width=lw)
    plt.setp(list(axes.spines.values()), linewidth=lw, color='black')

