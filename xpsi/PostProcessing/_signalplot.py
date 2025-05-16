from ._global_imports import *

try:
    import fgivenx
except ImportError:
    _warning('Cannot import fgivenx for conditional posterior contours.')
    fgivenx = None

from .. import Signal

class SignalPlot(object, metaclass=ABCMeta):
    """ Base class for a signal plot.

    Interact with a single :class:`~xpsi.Signal.Signal` instance to plot
    signals.

    :param str fig_dir:
        Directory to write to.

    :param str root_filename:
        Prepended to filename automatically generated figure filename.

    :param str cmap:
         Colormap name from :mod:`matplotlib` to use for :mod:`fgivenx`
         contours. It is advisable to choose a colormap that is darker
         at lower values (in tails of conditional posterior distribution) and
         lighter at higher values (at higher densities in conditional posterior
         distribution).

    """
    __figtype__ = None

    # When you subclass, assign values to the following specials and then
    # do not change the values at runtime unless the methods behave
    # accordingly in response to such changes. This is because the methods
    # might simply assume the values and have not adapt to changes in those
    # values, which might require conditional statements to control the flow
    # of execution. It is best to consider these specials as hard-coded/
    # invariant, and not runtime variables because they describe the nature
    # and distribution of information for a given plot type, instead of being
    # boolean execution switches or settings whose values are independent of the
    # flow of execution.
    __caching_targets__ = None

    __rows__          = None # of large panels (of size cls._panelsize)
    __columns__       = None # of large panels (of size cls._panelsize)
    __ax_rows__       = None # of all axes including colorbar axes
    __ax_columns__    = None # of all axes including colorbar axes
    __height_ratios__ = None # of all axes
    __width_ratios__  = None # of all axes
    __wspace__        = None # of all axes
    __hspace__        = None # of all axes

    @classmethod
    @make_verbose('Declaring plot class settings','Settings declared')
    def declare_settings(cls,
                         panelsize=(8,4),
                         fscale=1.0,
                         rasterized=True,
                         tick_lengths=(8,4),
                         tick_width=1.0,
                         logspace_y=False,
                         tqdm_kwargs={'disable': False},
                         scale_ymin=0.9,
                         scale_ymax=1.1,
                         ny=500,
                         lines_on=False,
                         add_label=True,
                         write=True,
                         extension='.pdf',
                         dpi=300,
                         write_kwargs=None,
                         xticks=None,
                         yticks=None):
        """ Configure subclass with attributes, *before instantiation*.

        These attributes can in principle be safely changed during the runtime.

        :param tuple panelsize:
            The approximate size, ``(width, height)``, in inches *per panel*.

        :param float fscale:
            Scale factor for the fontsize applied to the default.

        :param bool rasterized:
            Beware that if ``False``, white polygon edges may appear in figures.

        :param tuple(int, int) tick_lengths:
            The major and minor tick lengths for all axis objects.

        :param int tick_width:
            The width of every tick.

        :param bool logspace_y:
            Should the :mod:`fgivenx` y-values (for which the posterior mass
            is integrated over points at lower densities than the density at y)
            be logarithimcally or linearly spaced? Appropriate for logarithmic
            y-axes of signal panels.

        :param dict tqdm_kwargs:
            Progress bar keyword arguments.

        :param float scale_ymin:
            The minimum y-value to use, as a fraction of the minimum sample
            value.

        :param float scale_ymax:
            The maximum y-value to use, as a fraction of the minimum sample
            value.

        :param int ny:
            Number of y-values to compute the posterior mass for.

        :param bool lines_on:
            Render contour lines?

        :param bool add_label:
            Add contour colorbar label.

        :param bool write:
            Write figure to disk?

        :param str extension:
            File extension for writing.

        :param int dpi:
            Figure write resolution.

        :param dict write_kwargs:
            Additional writing keyword arguments.

        :param list(list(float)) xticks:
            Determine the x-tick values that will be labelled in each of the subplots
            (including color bars). Automatic ticks and labels are used for the
            specific subplot if providing None instead of a list of floats
            (e.g., xticks=[[0.1,0.5,1.0],None,None]).
            If xticks=None (default), automatic ticks and labels are used for all the
            subplots.

        :param list(list(float)) yticks:
            Determine the y-tick values that will be labelled in each of the subplots
            (including color bars). See param xticks for more details.

        """

        if panelsize[1] >= panelsize[0]:
            yield ('Warning: plots are generally better rendered such '
                   'that panel width is greater than panel height. Use the '
                   '``figsize`` argument to specify figure size per panel.')

        cls._panelsize = panelsize
        cls._fscale = fscale
        cls._fontsize = cls._fscale * 3.25 * min(panelsize[0], panelsize[1])

        cls._rasterized = rasterized
        cls._tick_lengths = tick_lengths
        cls._tick_width = tick_width

        cls._logspace_y = logspace_y
        cls._tqdm_kwargs = tqdm_kwargs
        cls._scale_ymin = scale_ymin
        cls._scale_ymax = scale_ymax
        cls._ny = ny
        cls._lines_on = lines_on
        cls._add_label = add_label

        cls._write = write

        cls._ext = extension
        cls._dpi = int(dpi)
        if not write_kwargs:
            cls._write_kwargs = dict(bbox_inches='tight')

        cls._xticks = xticks
        cls._yticks = yticks
        cls._settings_declared = True

        yield

    def __init__(self,
                 fig_dir = './',
                 root_filename = '',
                 cmap = 'RdPu_r',
                 **kwargs):
        try:
            type(self)._settings_declared
        except AttributeError:
            self.declare_settings() # defaults
        else:
            if not type(self)._settings_declared:
                self.declare_settings() # defaults

        # set fontsize globally for use in subclass initializer for axes
        # creation; note that __exit__ from context finalizes plotting
        # before another plot object is handled and also ensures fontsize
        # is set as the class fontsize in case needed
        rc('font',**{'size': self._fontsize}) # per plot fontsize settings

        self._fig_dir = fig_dir
        self._root_filename = root_filename

        self._cmap = cmap

        # will shadow class attributes if you really do want instance-specific
        # attributes
        for key, value in kwargs.items():
            key = '_' + key if key[0] != '_' else key
            try:
                getattr(self, key)
            except AttributeError:
                pass
            else:
                if key[1] != '_': # try to protect specials
                    setattr(self, key, value)

        self._axes = [] # in subclass init, append axes to automatically veneer

    @property
    def posterior(self):
        try:
            return self._posterior
        except AttributeError:
            return ''

    @posterior.setter
    def posterior(self, ID):
        if not isinstance(ID, _six.string_types):
            raise TypeError('Invalid type for ID.')
        else:
            self._posterior = ID

    @property
    def signal(self):
        return self._signal

    @signal.setter
    def signal(self, obj):
        if not isinstance(obj, Signal):
            raise TypeError('Invalid type for signal object.')
        else:
            self._signal = obj

    def __next__(self):
        """ Redirect for Python 3 compatibility. """
        return next(self)

    def __enter__(self):
        """ Plots are a resource to be managed. """
        rc('font',**{'size': self._fontsize}) # per plot fontsize settings
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        if exc_type is not None:
            raise

        self.finalize()

        try:
            self._axes
        except AttributeError: # assume no axes to veneer
            pass
        else:
            try:
                iter(self._axes)
            except TypeError: # quietly do nothing
                pass
            else:
                iax = 0
                for ax in self._axes:
                    if self._xticks is None and self._yticks is None:
                        self._veneer(ax)
                    elif self._xticks is not None and self._yticks is None:
                        self._veneer(ax,xticks=self._xticks[iax])
                        if self._xticks[iax] is not None:
                            ax.set_xticklabels(self._xticks[iax])
                    elif self._yticks is not None and self._xticks is None:
                        self._veneer(ax,yticks=self._yticks[iax])
                        if self._yticks[iax] is not None:
                            ax.set_yticklabels(self._yticks[iax])
                    else:
                        self._veneer(ax,xticks=self._xticks[iax],yticks=self._yticks[iax])
                        if self._yticks[iax] is not None:
                            ax.set_yticklabels(self._yticks[iax])
                        if self._xticks[iax] is not None:
                            ax.set_xticklabels(self._xticks[iax])
                    iax = iax+1

        if self._write: self.savefig()

        # the following handles possible hanging current figure ref
        # if comment out this line, a figure and axis of default
        # matplotlib size seems to be created somewhere and output
        # in a notebook; need to figure ;) out why that happens
        plt.close()

    @abstractmethod
    def execute(self, thetas, wrapper):
        """ Execute behaviours requiring a loop over samples.

        :param ndarray[n,m] thetas:
            Equally-weighted set of posterior samples.

        :param callable wrapper:
            Wrapper function that interleaves iterative calls `next(self)`
            with code to update the :class:`~.Signal.Signal` instance from a
            cache stored on disk. The wrapper returns a callback function,
            If such a cache is unavailable, the likelihood function is
            evaluated on-the-fly when the callback is called passing a sample
            (parameter vector). In this case, signals to be plotted are
            temporarily cached in the :class:`~.Signal.Signal` instance.
            The signature of the wrapper function is:
            `wrapper(self, delete_me=None, index=None)`,
            where `self` is this a subclass of
            :class:`~.SignalPlot`, `delete_me` is an optional
            attribute name or a list of attribute name in `self.__dict__` to be
            deleted before iterating over samples, and `index` is an optional
            integer to index a tuple returned from `next(self)`.
            The signature of the callback function, if you need to call it
            directly (instead of an external library such as fgivenx calling
            it):
            `callback(None, theta)`,
            where each sample/parameter vector is a row in `thetas`.

        """
        pass

    def __next__(self):
        """ Redirect for Python 3 compatibility. """
        return next(self)

    #@abstractmethod
    #def __next__(self):
    #    """ Execute next iteration. """
    #    pass

    @abstractmethod
    def finalize(self):
        """ Execute instructions to finish plotting.

        Before this method is called, the likelihood object will be updated
        so that the true (injected) signals are cached. In the body of the
        subclass implementation of this method, one can plot the injected
        signals by accessing the attributes of the signal instance.

        """
        pass

    def _veneer(self, ax, xticks=None, yticks=None):
        """ These are globally applicable settings for aesthetics. """

        if not isinstance(ax, _mpl.axes.Axes):
            return

        ax.tick_params(axis='both',
                       which='major',
                       colors='black',
                       length=self._tick_lengths[0],
                       width=self._tick_width,
                       direction='out')
        ax.tick_params(axis='both',
                       which='minor',
                       colors='black',
                       length=self._tick_lengths[1],
                       width=self._tick_width,
                       direction='out')

        for spine in ax.spines.values():
            spine.set_linewidth(self._tick_width)

        if xticks:
            ax.set_xticks(xticks)
        if yticks:
            ax.set_yticks(yticks)

    def _add_contours(self,
                      callback,
                      thetas,
                      x,
                      ax,
                      ax_cb,
                      **kwargs):
        """ Plot contours with :mod:`fgivenx`.

        In order to increase control, we bypass direct use of some functions
        from :mod:`fgivenx.drivers`.

        Avoid :mod:`fgivenx` caching because it caching is implemented in
        :mod:`xpsi.PostProcessing` for derived signals, which are generally
        needed beyond a single application of :mod:`fgivenx`. We thus do not
        currently cache the masses computed by :mod:`fgivenx`, but this is
        acceptable in terms of the speed to replot the contours.

        .. todo::

            Determine whether parallelistion is straighforwardly possible
            within :mod:`fgivenx` given that the callback is complicated
            object. Currently ``OpenMP`` threading is enabled by the
            likelihood object instead.

        """
        if fgivenx is None:
            raise ImportError('Install fgivenx to plot contours.')

        tqdm_kwargs = kwargs.get('tqdm_kwargs', self._tqdm_kwargs)

        fsamps = fgivenx.samples.compute_samples(f=[callback],
                                                 x=x,
                                                 samples=[thetas],
                                                 cache=None,
                                                 parallel=False,
                                                 tqdm_kwargs=tqdm_kwargs)

        ymin = fsamps[~_np.isnan(fsamps)].min(axis=None)
        ymin *= kwargs.get('scale_ymin', self._scale_ymin)
        ymax = fsamps[~_np.isnan(fsamps)].max(axis=None)
        ymax *= kwargs.get('scale_ymax', self._scale_ymax)
        ny = kwargs.get('ny', self._ny)
        logspace = kwargs.get('logspace_y', self._logspace_y)
        if not logspace:
            y = _np.linspace(ymin, ymax, ny)
        else:
            y = _np.logspace(_np.log10(ymin),
                             _np.log10(ymax),
                             ny,
                             base=10.0)

        # posterior mass integration more robust under transformation
        _fsamps = fsamps if not logspace else _np.log10(fsamps)
        _y = y if not logspace else _np.log10(y)

        pmf = fgivenx.mass.compute_pmf(fsamps=_fsamps,
                                       y=_y,
                                       parallel=False,
                                       cache=None, tqdm_kwargs=tqdm_kwargs)

        cb = fgivenx.plot.plot(x=x, y=y, z=pmf,
                           ax=ax,
                           colors=plt.get_cmap(kwargs.get('cmap', self._cmap)),
                           lines=kwargs.get('lines_on', self._lines_on),
                           rasterize_contours=self._rasterized,
                           smooth=False) # no additional smoothing with kernel

        cb = plt.colorbar(cb, cax=ax_cb, ticks=[1,2,3]) # ticks units = sigmas
        cb.ax.set_frame_on(True)
        cb.ax.yaxis.set_minor_locator(AutoMinorLocator())
        if kwargs.get('add_label', self._add_label):
            cb.ax.set_yticklabels([r'$1\sigma$', r'$2\sigma$',
                                   r'$3\sigma$'])

    def _add_signal(self, ax, x, signal, axis=None, **kwargs):
        """ Add signal line as a 1D function. """
        if not kwargs:
            kwargs.update(dict(color='k', ls='-', lw=0.05, alpha=1.0))

        if axis is not None:
            signal = _np.sum(signal, axis=axis)

        ax.plot(x, signal, **kwargs)

    def _get_figure(self):
        """ Create and return figure handle. """
        cls = type(self) # do not permit shadowing instance attributes
        self._fig = plt.figure(figsize = (self._panelsize[0] * cls.__columns__,
                                          self._panelsize[1] * cls.__rows__))

        self._gs = gridspec.GridSpec(cls.__ax_rows__,
                                     cls.__ax_columns__,
                                     height_ratios = cls.__height_ratios__,
                                     width_ratios = cls.__width_ratios__)

        self._gs.update(wspace = cls.__wspace__ * self._fscale,
                        hspace = cls.__hspace__ * self._fscale)

        return self._fig, self._gs

    @property
    def fig(self):
        """ Get the figure object. """
        try:
            return self._fig
        except AttributeError:
            pass # nothing to return

    def _add_subplot(self, i, j):
        """ Add subplot and append a reference to the new subplot axes object.

        To work with a subgrid, you must implement it manually. See, e.g.,
        :class:`~xpsi.PostProcessing._spectrum.SpectrumPlot`.

        """
        subplot = self._fig.add_subplot(self._gs[i,j])
        self._axes.append(subplot)
        return subplot

    @make_verbose('Writing plot to disk', 'Written')
    def savefig(self):
        """ Write figure to file. """

        fileroot = (self._root_filename + '__' if self._root_filename else '')
        fileroot += self.posterior.replace(' ', '_')

        filename = _os.path.join(self._fig_dir,
                                 fileroot + '__' + self.__figtype__ + self._ext)

        yield (self.__class__.__name__ + \
                ' instance plot will be written to path %s' % filename)

        try:
            self._fig.savefig(filename, dpi=self._dpi, **self._write_kwargs)
        except TypeError: # quietly proceed
            self._fig.savefig(filename, dpi=self._dpi)

        yield
