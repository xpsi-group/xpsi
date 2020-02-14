from __future__ import division, print_function

from .global_imports import *
from . import global_imports

from .Spacetime import Spacetime
from .HotRegion import HotRegion
from .Elsewhere import Elsewhere
from .Everywhere import Everywhere

from .Parameter import Parameter
from .ParameterSubspace import ParameterSubspace

from .pixelmesh.integrator import integrate as _integrate
from .tools.energy_integrator import energy_integrator
from .tools.phase_integrator import phase_integrator

try:
    _mpl
except NameError:
    pass
else:
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

class Photosphere(ParameterSubspace):
    """ A photosphere embedded in an ambient Schwarzschild spacetime.

    :param obj hot:
        An instance of :class:`~.HotRegion.HotRegion` (or a
        derived class). This objects represents the hot
        regions of the surface that in most use-cases will be
        assumed to contain radiating material that is hotter
        than that *elsewhere*.

    :param obj elsewhere:
        An instance of :class:`~.Elsewhere.Elsewhere` (or a derived class).

    :param obj everywhere:
        An instance of :class:`~.Everywhere.Everywhere` (or a derived class).
        Note that if you use construct the surface radiation field in this
        way, you should use the :attr:`~.Photosphere.Photosphere.hot_atmosphere`
        property to pass a buffer of numerical data to the integrator
        routines. You then need to ensure that the extension modules
        ``xpsi/surface_radiation_field/hot_radiation_field.pyx`` and
        ``xpsi/surface_radiation_field/elsewhere_radiation_field.pyx`` match.

    .. note::

        You cannot specify the surface radiation field *everywhere* if you
        use hot regions (the latter usage may also include specification of
        the radiation field *elsewhere*).

    :param dict bounds:
        Bounds are supplied for instantiation of a frequency parameter.
        The parameter name ``'mode_frequency'`` must be a key in the
        dictionary unless the parameter is *fixed* or *derived*. If a bound
        is ``None`` that bound is set equal to a strict hard-coded bound.
        If ``None``, lock the coordinate rotation frequency of a mode of
        asymmetry in the photosphere to a fixed frequency, e.g., the stellar
        rotation frequency. If bounds are passed, the frequency is interpreted
        as a free parameter.

    :param dict values:
        Either the fixed value of the mode frequency, a callable if the
        frequency is *derived*, or a value upon initialisation if the
        frequency is free. The dictionary must have a key with name
        ``'mode_frequency'`` if it is *fixed* or *derived*.
        If the asymmetry is locked to the stellar spin, then you need to pass
        the spin frequency. If fixed but different to the spin frequency, this
        value needs to be passed instead. In the hot region base class this
        mode frequency is applied to normalise the ray lags instead of the
        stellar rotation frequency.

    .. note::

        In basic modelling patterns the frequency is the spin frequency,
        and thus you only need to explicitly pass the spin as ``value`` whilst
        leaving ``bounds`` to default. If the spin frequency happens to be a
        free parameter (perhaps with informative prior information), then
        pass a callable instead that can be used to get the spin frequency
        dynamically when the derived mode frequency variable is called for.

    """
    required_names = ['mode_frequency']

    def __init__(self,
                 hot = None, elsewhere = None,
                 everywhere = None,
                 bounds = {}, values = {},
                 **kwargs):

        if everywhere is not None:
            if hot or elsewhere is not None:
                raise ValueError('Cannot use hot region nor elsewhere '
                                 'functionality if constructing the '
                                 'radiation field everywhere.')
            if not isinstance(everywhere, Everywhere):
                raise TypeError('Invalid type for everywhere object.')
            self._everywhere = everywhere
        else:
            self._everywhere = None

            if elsewhere is not None:
                if not isinstance(elsewhere, Elsewhere):
                    raise TypeError('Invalid type for an elsewhere object.')
                else:
                    self._elsewhere = elsewhere
            else:
                self._elsewhere = None
                if hot is None:
                    raise ValueError('The photosphere must radiate.')

                                              # including derived classes
            if hot is not None and hot is not isinstance(hot, HotRegion):
                if hasattr(hot, 'objects'):
                    for obj in getattr(hot, 'objects'):
                        if not isinstance(obj, HotRegion):
                            raise TypeError('Invalid object for the hot '
                                            'region(s).')
                else:
                    raise TypeError('Invalid object for the hot region(s).')

            self._hot = hot

            self._elsewhere_atmosphere = ()

        self._hot_atmosphere = ()

        doc = """
        Coordinate frequency of the mode of radiative asymmetry in the
        photosphere that is assumed to generate the pulsed signal [Hz].
        """
        mode_frequency = Parameter('mode_frequency',
                                   strict_bounds = (0.0, 2000.0),
                                   bounds = bounds.get('mode_frequency', None),
                                   doc = doc,
                                   symbol = r'$f_{\rm mode}$',
                                   value = values.get('mode_frequency', None))

        super(Photosphere, self).__init__(mode_frequency,
                                          hot, elsewhere, everywhere,
                                          **kwargs)

    @property
    def hot_atmosphere(self):
        """ Get the numerical atmosphere buffers for hot regions if used.

        To preload a numerical atmosphere into a buffer, subclass and
        overwrite the setter. The underscore attribute set by the setter
        must be an :math:`n`-tuple whose :math:`n^{th}` element is an
        :math:`(n-1)`-dimensional array flattened into a one-dimensional
        :class:`numpy.ndarray`. The first :math:`n-1`
        elements of the :math:`n`-tuple must each be an ordered one-dimensional
        :class:`numpy.ndarray` of parameter values for the purpose of
        multi-dimensional interpolation in the :math:`n^{th}` buffer. The
        first :math:`n-1` elements must be ordered to match the index
        arithmetic applied to the :math:`n^{th}` buffer. An example would be
        ``self._hot_atmosphere = (logT, logg, mu, logE, buf)``, where:
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

        """
        return self._hot_atmosphere

    @hot_atmosphere.setter
    def hot_atmosphere(self, path):
        """ Implement if required. """
        raise NotImplementedError('Implement setter if required.')

    @property
    def elsewhere_atmosphere(self):
        """ Get the numerical atmosphere buffers for elsewhere if used.

        To preload a numerical atmosphere into a buffer, subclass and
        overwrite the setter. The underscore attribute set by the setter
        must be an :math:`n`-tuple whose :math:`n^{th}` element is an
        :math:`(n-1)`-dimensional array flattened into a one-dimensional
        :class:`numpy.ndarray`. The first :math:`n-1`
        elements of the :math:`n`-tuple must each be an ordered one-dimensional
        :class:`numpy.ndarray` of parameter values for the purpose of
        multi-dimensional interpolation in the :math:`n^{th}` buffer. The
        first :math:`n-1` elements must be ordered to match the index
        arithmetic applied to the :math:`n^{th}` buffer. An example would be
        ``self._hot_atmosphere = (logT, logg, mu, logE, buf)``, where:
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

        """
        return self._elsewhere_atmosphere

    @elsewhere_atmosphere.setter
    def elsewhere_atmosphere(self, path):
        """ Implement if required. """
        raise NotImplementedError('Implement setter if required.')

    @property
    def hot(self):
        """ Get the instance of :class:`~.HotRegion.HotRegion`. """
        return self._hot

    @property
    def elsewhere(self):
        """ Get the instance of :class:`~.Elsewhere.Elsewhere`. """
        return self._elsewhere

    @property
    def everywhere(self):
        """ Get the instance of :class:`~.Everywhere.Everywhere`. """
        return self._everywhere

    @property
    def spacetime(self):
        """ Return instance of :class:`~.Spacetime.Spacetime`. """
        return self._spacetime

    @spacetime.setter
    def spacetime(self, obj):
        if not isinstance(obj, Spacetime):
            raise TypeError('Invalid type for spacetime object.')
        # otherwise store a reference to the spacetime object
        self._spacetime = obj

    def embed(self, fast_total_counts, threads):
        """ Embed the photosphere in an ambient Schwarzschild spacetime.

        In other words, generate a discrete representation of the photospheric
        radiation field and the null mapping from the photosphere to infinity,
        for use in flux integrators called by distant observers.

        """
        if self._everywhere is not None:
            self._everywhere.embed(self._spacetime,
                                   self,
                                   threads)
        else:
            if self._elsewhere is not None:
                self._elsewhere.embed(self._spacetime, threads)

                if self._hot is not None:
                    self._hot.embed(self._spacetime,
                                    self,
                                    fast_total_counts,
                                    threads,
                                    self._elsewhere._compute_cellParamVecs)
            else:
                self._hot.embed(self._spacetime,
                                self,
                                fast_total_counts,
                                threads)

    def integrate(self, energies, threads):
        """ Integrate over the photospheric radiation field.

        :param energies:
            A one-dimensional :class:`numpy.ndarray` of energies in keV.

        :param int threads:
            Number of ``OpenMP`` threads to spawn for pulse integration.

        """
        if self._everywhere is not None:
            # temp var
            t = self._everywhere.integrate(self._spacetime,
                                           energies,
                                           threads,
                                           self._hot_atmosphere)
            if t.ndim == 1:
                self._time_invariant = t
                self._pulse = ((_np.zeros((len(t),
                                    len(self._everywhere.phases_in_cycles))),),)
                for i in range(self._pulse[0][0].shape[1]):
                    self._pulse[0][0][:,i] += self._time_invariant
            else:
                self._pulse = ((t,),)
        else:
            if self._elsewhere is not None:
                if isinstance(energies, tuple): # resolve energy container type
                    if not isinstance(energies[0], tuple):
                        _energies = energies[0]
                    else:
                        _energies = energies[0][0]
                else:
                    _energies = energies
                self._time_invariant = self._elsewhere.integrate(self._spacetime,
                                                       _energies,
                                                       threads,
                                                       *self._elsewhere_atmosphere)

            if self._hot is not None:
                self._pulse = self._hot.integrate(self._spacetime,
                                                  energies,
                                                  threads,
                                                  self._hot_atmosphere,
                                                  self._elsewhere_atmosphere)

                if not isinstance(self._pulse[0], tuple):
                    self._pulse = (self._pulse,)

                # add time-invariant component to first time-dependent component
                if self._elsewhere is not None:
                    for i in range(self._pulse[0][0].shape[1]):
                        self._pulse[0][0][:,i] += self._time_invariant

    @property
    def pulse(self):
        """ Get the stored pulse.

        :returns:
            *ndarray[m,n]*, where :math:`m` is the number of energies, and
            :math:`n` is the number of phases. Units are photon/s/keV; the
            distance is a fast parameter so the areal units are not yet
            factored in.

        """
        return self._pulse

    @property
    def time_invariant(self):
        """ Get the cached time-invariant signal.

        :returns:
            *ndarray[n]* containing the time-invariant signal at
            :math:`n` energies. Units are photon/s/keV. The distance is a
            fast parameter so the areal units are not yet factored in.

        """
        return self._time_invariant

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
    def images(self):
        """ Get the precomputed image information. """
        return self._images

    @images.setter
    def images(self, images):
        """ Store an *ndarray[i,j,k]* of images. """
        try:
            for i, obj in enumerate(images):
                if not isinstance(obj, _np.ndarray):
                    if i < len(images) - 1:
                        raise TypeError('All image information must be '
                                        'contained in ndarrays.')
                    elif obj is not None:
                        raise TypeError('All image information must be '
                                        'contained in ndarrays.')
        except TypeError:
            raise TypeError('An iterable of objects containing image '
                            'information must be supplied.')

        if len(images) != 9:
            raise ValueError('There must be six ndarray objects specifing '
                             'image information.')

        msg = 'Image information element %i must have %i dimensions.'

        assert images[0].ndim == 2, msg % (0, 2)
        assert images[1].ndim == 1, msg % (1, 1)
        assert images[2].ndim == 1, msg % (2, 1)
        assert images[3].ndim == 1, msg % (3, 1)
        assert images[4].ndim == 1, msg % (4, 1)
        assert images[5].ndim == 1, msg % (5, 1)
        assert images[6].ndim == 1, msg % (6, 1)
        assert images[7].ndim == 1, msg % (7, 1)
        if images[8] is not None:
            assert images[8].ndim == 3, msg % (8, 3)

        self._images = images

    def image(self,
              reimage = False,
              energies = None,
              phases = None,
              sqrt_num_rays = 100,
              epsabs_ray = 1.0e-12,
              epsrel_ray = 1.0e-12,
              max_steps = 100000,
              init_step = 0.1,
              image_plane_radial_increment_power = 1.0 / 2.0,
              threads = 1,
              cache_intensities = False,
              plot_sky_maps = False,
              sky_map_kwargs = {},
              animate_sky_maps = False,
              animate_kwargs = {}):
        """ Image the star as a function of phase and energy.

        :param bool reimage:
            Image the star. Ignored if there is no precomputed image information
            to plot.

        :param ndarray[n] energies:
            Energies in keV to evaluate incident specific intensities at.

        :param ndarray[m] phases:
            Phases in radians at which to evaluate incident specific
            intensities at.

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

            Moreover, as an additional measure against
            artefacts in the sky maps in the vicinity of the rotational pole,
            rays are distributed accordingingly. For example, if we request
            :math:`n=400` rays per dimension, a maximal spacing of the rays
            from the rotational axis is achieved by rotating the *spokes*
            of rays (by up to :math:`\pm\pi/n`) so that no spoke is
            aligned (or anti-aligned) with the :math:`y`-direction.

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
            sequence of rotational phases..

        :param bool cache_intensities:
            Cache the photon specific intensity sky maps in memory, as a
            function of phase and energy? Defaults to ``False`` because this
            dominates memory consumption. You need to activate this option
            if you want to plot the sky maps (see below).

        :param bool plot_sky_maps:
            Plot (specific) intenity sky maps at a sequence of phases, or
            by averaging over phase. Maps can be made at one more energies
            or energy intervals. The images will be written to disk and
            can be used as frames in an animated sequence.

        :param dict sky_map_kwargs:
            Dictionary of keyword arguments passed to
            :meth:`~Photosphere._plot_sky_maps`. Refer to the associated
            method docstring for available options.

        :param bool animate_sky_maps:
            Compile images from disk into an animated sequence.

        :param dict animate_kwargs:
            Dictionary of keyword arguments passed to
            :meth:`~Photosphere._animate`. Refer to the associated
            method docstring for available options.

        """
        ref = self._spacetime # geometry shortcut saves characters

        try:
            self.images
        except AttributeError:
            if not reimage:
                print('Warning: star will not be reimaged... assuming images '
                      'exist on disk.')

        if reimage:
            if not isinstance(phases, _np.ndarray):
                raise TypeError('Imaging phases must be form an ndarray.')

            if not isinstance(energies, _np.ndarray):
                raise TypeError('Imaging energies must be form an ndarray.')

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
                                self._hot_atmosphere)

            if images[0] == 1:
                raise Exception('A numerical error arose during imaging '
                                'computation... terminating simulation.')
            else:
                # tuple elements:
                #   energy-phase resolved pulse (2D array)
                #   x coordinate on image plane (1D array)
                #   y coordinate on image plane (1D array)
                #   colatitude mapped to point (x,y) on image plane (1D array)
                #   azimuth mapped to point (x,y) on image plane (1D array)
                #   phase lag
                #   redshift
                #   aberrated ray angle to local surface normal
                #   energy-phase resolved specific intensity sky maps (3D array)
                # the last element is None if you do not cache intensities
                self.images = list(images[1:])

                # transpose so pulse phase increments along columns
                self.images[0] = self.images[0].T

        if plot_sky_maps or animate_sky_maps:
            if self.images[-1] is None:
                raise ValueError('You need to cache intensity sky maps if you '
                                 'want to plot them.')
            root_dir = sky_map_kwargs.pop('root_dir', './images')
            file_root = sky_map_kwargs.pop('file_root', 'skymap')
            file_root = _os.path.join(root_dir, file_root)

            phase_average = sky_map_kwargs.get('phase_average', False)
            if phase_average and animate_sky_maps:
                raise ValueError('Phase averaged sky maps cannot be animated.')

        if plot_sky_maps:
            if not _os.path.isdir(root_dir):
                _os.mkdir(root_dir)
            elif _os.path.isfile(file_root + '_0.png'):
                print('\nWarning: at least one image file exists '
                      'in ``%s``.' % root_dir)
                print('Attempting to move image files to a subdirectory '
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
                    print('Image files archived in subdirectory ``%s``.' % temp)

            figsize, dpi = self._plot_sky_maps(file_root,
                                               _phases = phases,
                                               _energies = energies,
                                               _redraw = True,
                                               **sky_map_kwargs)

        elif animate_sky_maps:
            if reimage:
                raise ValueError('Star was reimaged but sky maps were not '
                                 'plotted... aborting animation.')

            figsize, dpi = self._plot_sky_maps(file_root,
                                               _energies = energies,
                                               _redraw = False,
                                               threads = threads,
                                               **sky_map_kwargs)

        if animate_sky_maps:
            if not _os.path.isfile(file_root + '_0.png'):
                raise IOError('No images located for animation.')

            if reimage:
                num_frames = self.images[-1].shape[0]
            else:
                try:
                    num_frames = len(phases)
                except TypeError:
                    raise TypeError('You need to declare the image phases '
                                    'in order to include all images from disk.')

            self._animate(file_root, num_frames,
                          figsize, dpi, **animate_kwargs)

    def _plot_sky_maps(self,
                       _file_root,
                       _phases,
                       _energies,
                       _redraw,
                       threads = 1,
                       panel_layout = None,
                       panel_indices = None,
                       phase_average = False,
                       energy_bounds = [],
                       num_levels = 100,
                       normalise_each_panel = True,
                       invert = False,
                       colormap = None,
                       figsize = (10,10),
                       usetex = False,
                       fontsize_scale = 1.0,
                       tick_spacing = (0.2,1.0),
                       tick_length_scaling = 1.0,
                       dpi_scale = 1.0,
                       **kwargs):
        """ Helper method for specific intensity sky map visualization.

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

        :param bool _redraw:
            Redraw the sky maps? This is handled internally, so do *not* pass
            a keyword argument.

        :param tuple[int,int] panel_layout:
            Two elements: the number of rows and columns of panels. If ``None``
            a layout is automatically determined based on the number of
            images to be plotted.

        :param iterable panel_indices:
            These ordered integers will be used to select intensity information
            by indexing the energy dimension of a 3D intensity array. If
            specific intensites are plotted, these integers should index
            a subset of energies at which the star was imaged. If intensities
            are plotted, these integers should index a subset of the energy
            intervals over which specific intensities are integrated. See
            the ``energy_bounds`` keyword argument.

        :param bool phase_average:
            Average each sky map over one revolution (cycle) in phase?
            Note that the resulting image is incompatible with the currently
            supported animation mode.

        :param iterable energy_bounds:
            A set of two-element containers. Each container has an ordered pair
            of energies which delimit an integration domain. Specific intensity
            is integrated along each sky direction, at each phase, between
            these energy bounds. The bounds must be between the minimum and
            maximum energies at which the star was imaged. If ``None``,
            specific intensity sky maps will be plotted (the default).

        :param int num_levels:
            Number of contour levels in (specific) intensity, distributed
            between minimum finite, and maximum values per panel, or over all
            panels. See ``normalise_each_panel`` keyword argument.

        :param bool normalise_each_panel:
            Normalise the contour grey scale to each panel uniquely, or
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
            A matplotlib colormap object. Choose something appropriate and
            *accessible* for a non-negative scalar field (sky intensities).

        :param tuple(int,int) figsize:
            The figure size (width, height) in *inches*. If the dimensions are
            inconsistent with the aspect ratio suggested by the ``panel_layout``
            settings, the height of the figure will be automatically rescaled
            to achieve congruence, meaning each panel is approximately square.

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
            units are gravitational radii (:math:`GM/c^2`).

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

        if panel_layout is None:
            x = int(_m.ceil(_m.sqrt(len(panel_indices))))
            if x * (x - 1) >= len(panel_indices):
                panel_layout = (x, x - 1)
            else:
                panel_layout = (x, x)

        # try to improve the aspect ratio so that each panel is
        # approximately square
        width = panel_layout[1] + (panel_layout[1] - 1) * 0.2
        height = panel_layout[0] + (panel_layout[0] - 1) * 0.2
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

            if not isinstance(energies, _np.ndarray):
                raise TypeError('Imaging energies must be form an ndarray.')

            images = self.images[-1]

            if energy_bounds:
                for bounds in energy_bounds:
                    if bounds[0] > bounds[1]:
                        raise ValueError('Energy bounds in a tuple must be '
                                         'ordered.')
                    for bound in bounds:
                        if not energies[0] <= bound <= energies[-1]:
                            raise ValueError('Extrapolation would be required.')

                if len(panel_indices) < len(energy_bounds):
                    print('Warning: Fewer panels than energy intervals.')

                integrated = _np.zeros((images.shape[0],
                                        len(energy_bounds),
                                        images.shape[2]), dtype = _np.double)

                intensities = _np.zeros((images.shape[2],
                                         images.shape[1]), dtype = _np.double)

                for i in range(images.shape[0]): # phases
                    for j in range(images.shape[2]): # sky directions
                        intensities[j,:] = images[i,:,j]

                    for k in range(len(energy_bounds)):
                        a, b = energy_bounds[k]
                        _integrated = energy_integrator(threads,
                                                        intensities,
                                                        _np.log10(energies),
                                                        _np.log10(a),
                                                        _np.log10(b))

                        integrated[i,k,:] = _integrated[:]

                images = integrated

            else:
                if len(panel_indices) != len(energies):
                    print('Warning: Fewer panels than energies.')

            if phase_average:
                if phases[0] != 0.0 or phases[-1] != _2pi:
                    raise ValueError('Minimum and maximum phases at which star '
                                     'is imaged must be zero and unity if you '
                                     'are phase averaging.')

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

            if normalise_each_panel: # normalise intensity for each individual panel
                levels = []

                for j in range(images.shape[1]): # at each energy
                    # find extreme intensities over discrete set of image phases
                    # and sky directions
                    MIN = _np.min(images[:,j,:][images[:,j,:] > 0.0])
                    MAX = _np.max(images[:,j,:])
                    levels.append(_np.array([0.0] + list(_np.linspace(MIN,
                                                                      MAX,
                                                                      num_levels))))
            else:
                MIN = _np.min(images[:,:,:][images[:,:,:] > 0.0])
                MAX = _np.max(images[:,:,:])
                levels = _np.array([0.0] + list(_np.linspace(MIN,
                                                             MAX,
                                                             num_levels)))


            # because of default tick formatting and a minus sign,
            # the left and bottom margins need to be different
            left = 0.160 * (fontsize/28.0)
            bottom = 0.125 * (fontsize/28.0)
            right = 0.975
            top = bottom + (right - left)

            ref = self._spacetime

            cmap = colormap or (cm.gray_r if invert else cm.gray)

            fig = Figure(figsize = figsize) #plt.figure(figsize = figsize)
            canvas = FigureCanvas(fig)

            gs = gridspec.GridSpec(panel_layout[0],
                                   panel_layout[1],
                                   left=left, right=right,
                                   bottom=bottom, top=top,
                                   wspace=0.2, hspace=0.2)

            axes = [fig.add_subplot(gs[j]) for j in range(len(panel_indices))]

            for i in range(images.shape[0]):
                for j, idx in enumerate(panel_indices):
                    ax = axes[j]
                    if _np.product(panel_layout) - j - 1 < panel_layout[1]:
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
                    ax.set_facecolor('white' if invert else 'black')

                    lvls = levels if isinstance(levels, _np.ndarray) else levels[idx]

                    ax.tricontourf(X,
                                   Y,
                                   images[i,idx,:],
                                   cmap = cmap,
                                   levels = lvls)

                    # correct the aspect ratio
                    x_view = ax.xaxis.get_view_interval()
                    diff = x_view[1] - x_view[0]
                    ax.xaxis.set_view_interval(x_view[0] - diff * 0.025,
                                               x_view[1] + diff * 0.025)
                    y_view = ax.yaxis.get_view_interval()
                    ax.yaxis.set_view_interval(y_view[1] - diff * 1.025,
                                               y_view[1] + diff * 0.025)

                fig.savefig(file_root + '_%i.png' % i, dpi=dpi)

                for ax in axes:
                    ax.clear()

            for ax in axes:
                ax.cla()
            plt.close(fig)

        return figsize, dpi

    @staticmethod
    def _animate(_file_root, _num_frames, _figsize, _dpi,
                 cycles = 1, repeat = True, repeat_delay = 0.0,
                 ffmpeg_path = None, fps = None, **kwargs):
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

        :param bool repeat:
            Inform *ffmpeg* to enter a loop when video play back commences.

        :param repeat_delay:
            Delay between repeats in milliseconds.

        :param str ffmpeg_path:
            Absolute path to *ffmpeg* executable. If ``None``, defaults
            to matplotlib rcParams settings, but no guarantee that the package
            will be found even if installed on system.

        :param int fps:
            Frames per second. If ``None``, then one cycle (assuming images
            have been precomputed for a complete cycle), consisting of
            so many frames, will exhibit a period of one second.

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

        images = []

        for i in range(num_frames): # load phase-ordered set of images
            filename = file_root + '_%i.png' % i
            img = mgimg.imread(filename)
            imgplot = ax.imshow(img, aspect='auto')
            imgplot.axes.get_xaxis().set_visible(False)
            imgplot.axes.get_yaxis().set_visible(False)
            images.append([imgplot])

        cycles = int(cycles)

        if cycles < 1:
            cycles = 1
        elif cycles > 1:
            for i in range(cycles - 1):
                for j in range(1, num_frames):
                    images.append(images[j])

        ani = animation.ArtistAnimation(fig, images,
                                        repeat = repeat,
                                        repeat_delay = repeat_delay)

        if ffmpeg_path is None:
            print('Warning: no path specified for ffmpeg executable. '
                  'Resorting to rcParams default.')
        else:
            plt.rcParams['animation.ffmpeg_path'] = ffmpeg_path

        if fps is None:
            fps = num_frames # all frames span one second

        # secret keyword argument; not clear whether should be exposed to user
        bitrate = kwargs.get('bitrate', -1)

        filename = file_root + '_animated.mp4'
        print('Writing to disk: %s' % filename)
        ani.save(filename, writer = 'ffmpeg',
                 dpi = dpi, fps = fps, bitrate = bitrate)

        plt.close(fig)

Photosphere._update_doc()

def _veneer(x, y, axes, lw=1.0, length=8):
    """ Make the plots a little more aesthetically pleasing. """
    if x is not None:
        if x[1] is not None:
            axes.xaxis.set_major_locator(MultipleLocator(x[1]))
        if x[0] is not None:
            axes.xaxis.set_minor_locator(MultipleLocator(x[0]))
    else:
        axes.xaxis.set_major_locator(AutoLocator())
        axes.xaxis.set_minor_locator(AutoMinorLocator())

    if y is not None:
        if y[1] is not None:
            axes.yaxis.set_major_locator(MultipleLocator(y[1]))
        if y[0] is not None:
            axes.yaxis.set_minor_locator(MultipleLocator(y[0]))
    else:
        axes.yaxis.set_major_locator(AutoLocator())
        axes.yaxis.set_minor_locator(AutoMinorLocator())

    axes.tick_params(which='major', colors='black', length=length, width=lw)
    axes.tick_params(which='minor', colors='black', length=int(length/2), width=lw)
    plt.setp(axes.spines.values(), linewidth=lw, color='black')
