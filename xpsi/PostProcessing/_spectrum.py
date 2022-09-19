from ._global_imports import *

try:
    import fgivenx
except ImportError:
    _warning('Cannot import fgivenx for conditional posterior contours.')
    fgivenx = None

from ..tools.phase_integrator import phase_integrator
from ..tools.energy_interpolator import energy_interpolator
from ..tools.phase_interpolator import phase_interpolator as interp

from ._signalplot import SignalPlot

class SpectrumPlot(SignalPlot):
    """ Plot posterior-averaged channel count-rate spectra.

    The figure contains three panels which share phase as an x-axis:

        * the top panel displays the specific photon flux signal
          from the source, resolved as a function of energy, optionally showing
          both unattenuated and attenuated incident spectra and optionally
          using :mod:`fgivenx`;
        * the center panel displays the posterior expectation of the count-rate
          signal as a function of channel and phase, optionally including an
          expected background signal;
        * the bottom panel displays the phase-integrated (averaged) count-rate
          spectum as a function of channel number, optionally including an
          expected background signal and optionally using :mod:`fgivenx`.

    The top and bottom panels aim to render the conditional posterior
    distribution of the associated signal as a function an energy (proxy)
    variable, ideally with contours to map out the conditional posterior mass.
    These panels have space to optionally display other elements such as:
    the posterior-expected total signal; the posterior-expected component
    signals; the true total and component signals if the ground truth (the
    injected signal correpsonding to some model parameter vector) is known;
    attenuated incident spectra; and the summation of posterior-expected
    total (component-summed) source count-rate signals with posterior-expected
    background count-rate signals.

    The following example is (improved) from :ref:`R19`:

    .. image:: _static/_spectrumplot.png

    :param float rel_num_energies:
        The number of energies desired for interpolation as a fraction of the
        number of energies implemented for the original incident signal
        integration. The energy set will be appropriately spaced and bounded.

    :param int num_phases:
        The number of phases to interpolate the pulse-profile signals at for
        the center panel.

    :param str registered_cmap:
         Colormap name from :mod:`matplotlib` to use for the posterior-expected
         registered signal as a function of channel and phase (center panel).

    :param bool show_components:
        If the :class:`~.Signal.Signal` instance has multiple components (hot
        region signals), display the posterior expectations of those components
        as a function of energy (top panel) and channel (bottom panel).

    :param bool show_attenuated:
        If the source signal is attenuated by the interstellar absorption
        processes, display the posterior-expected attenuated incident specific
        photon flux spectra? This switch also instructs :mod:`fgivenx`, if
        enabled, to generate conditional posterior contours for the attenuated
        spectrum instead of the unattenuated spectrum (top panel). If
        :mod:`fgivenx` is not invoked, this switch instructs the plotting
        of sample-by-sample attenuated total (component-summed) spectra to
        delineate the distribution of conditional posterior mass (top panel).

    :param dict expectation_line_kwargs:
        Keyword arguments for plotting the posterior-expected signal lines (in
        the top and bottom panels).

    :param bool add_background:
        Add an posterior-expected background count-rate signal to the total
        (component-summed) expected source count-rate signals in the center
        and bottom panels?

    :param dict background_line_kwargs:
        Keyword arguments for plotting the posterior-expected spectrum in the
        bottom panel that includes the background signal.

    :param bool use_fgivenx:
        Use :mod:`fgivenx` to plot conditional posterior contours in the
        top and bottom panels?

    :param dict incident_contour_kwargs:
        Keyword arguments for :mod:`fgivenx` incident signal contours (top
        panel) that will take precedence over the corresponding class
        attributes. (See the :class:`~.SignalPlot` class if you choose
        not to modify these attributes on this present subclass.)

    :param dict registered_contour_kwargs:
        Keyword arguments for :mod:`fgivenx` registered signal contours (bottom
        panel) that will take precedence over the corresponding class
        attributes. (See the :class:`~.SignalPlot` class if you choose
        not to modify these attributes on this present subclass.)

    :param plot_truth:
        Plot the ground truth (injected) signal, if known and available, in
        the top and bottom panels.

    :param truth_line_kwargs:
        Keyword arguments for plotting the ground truth signal lines (top and
        bottom panels).

    :param comp_truth_line_kwargs:
        Keyword arguments for plotting the component ground truth signal lines
        (top and bottom panels).

    """
    __figtype__ = 'signalplot_spectrum'

    # do not change at runtime (see base class comments):
    __caching_targets__ = ['shifts',
                           'signals', # count-rate signals
                           # incident specific flux signals
                           'incident_specific_flux_signals']
    __rows__          = 3
    __columns__       = 1
    __ax_rows__       = 2
    __ax_columns__    = 1
    __height_ratios__ = [1,2]
    __width_ratios__  = [1]
    __wspace__        = 0.025
    __hspace__        = 0.175

    @make_verbose('Instantiating a spectrum plotter for posterior checking',
                  'Spectrum plotter instantiated')
    def __init__(self,
                 rel_num_energies=10.0,
                 num_phases=1000,
                 registered_cmap='inferno',
                 show_components=False,
                 show_attenuated=True,
                 expectation_line_kwargs=None,
                 comp_expectation_line_kwargs=None,
                 add_background=False,
                 background_line_kwargs=None,
                 sample_line_kwargs=None,
                 use_fgivenx=False,
                 incident_contour_kwargs=None,
                 registered_contour_kwargs=None,
                 plot_truth=False,
                 truth_line_kwargs=None,
                 comp_truth_line_kwargs=None,
                 **kwargs):

        try:
            _shadow = not self._logspace_y
        except AttributeError:
            _shadow = True

        if _shadow: # shadow class attribute
            kwargs.setdefault('logspace_y', True)
            if not kwargs.get('logspace_y'):
                yield ('Spectrum plots may have conditional-probability '
                       'contour artefacts if logspace_y is not True.')

        super(SpectrumPlot, self).__init__(**kwargs)

        self._rel_num_energies = rel_num_energies
        self._phases = _np.linspace(0.0, 2.0, int(num_phases))

        self._show_attenuated = show_attenuated

        self._add_background = add_background
        if add_background: # count-rate spectrum
            self.__caching_targets__.append('background_signal')

        if use_fgivenx and fgivenx is None:
            raise ImportError('Install fgivenx to plot contours.')

        self._use_fgivenx = use_fgivenx
        if self._use_fgivenx:
            self._incident_contour_kwargs =\
                    incident_contour_kwargs if incident_contour_kwargs else {}
            self._registered_contour_kwargs =\
                registered_contour_kwargs if registered_contour_kwargs else {}

        self._get_figure()
        fig = self._fig
        gs = self._gs
        cls = type(self)

        # for the incident specific flux signal incident on an instrument
        gs_top = gridspec.GridSpecFromSubplotSpec(1, 2,
                                                  subplot_spec=gs[0,0],
                                                  wspace=cls.__wspace__,
                                                  width_ratios=[50,1])

        self._ax_incident = fig.add_subplot(gs_top[0,0])

        # for the count-rate signal registered by an an instrument
        gs_bottom = gridspec.GridSpecFromSubplotSpec(2, 2,
                                                     subplot_spec=gs[1,0],
                                                     wspace=cls.__wspace__,
                                                     hspace=0.125*self._fscale,
                                                     height_ratios=[1,1],
                                                     width_ratios=[50,1])

        self._ax_registered = fig.add_subplot(gs_bottom[0,0])

        self._ax_registered_1d = fig.add_subplot(gs_bottom[1,0])

        self._axes = [self._ax_incident,
                      self._ax_registered,
                      self._ax_registered_1d]

        # incident axis properties
        self._ax_incident.set_xlabel(r'$E$ [keV]')
        self._ax_incident.set_xscale('log')

        self._ax_incident.set_ylabel(r'photons/keV/cm$^{2}$/s')
        self._ax_incident.set_yscale('log')

        # registered axis properties
        for ax in self._axes[1:]:
            ax.set_xscale('log')
        self._ax_registered.tick_params(axis='x', which='both',
                                        labelbottom=False)
        self._ax_registered_1d.set_xlabel('channel')

        self._ax_registered.set_ylabel(r'$\phi$ [cycles]')
        self._ax_registered.yaxis.set_major_locator(MultipleLocator(0.5))
        self._ax_registered.yaxis.set_minor_locator(MultipleLocator(0.1))
        self._ax_registered.set_ylim([0.0,2.0])

        self._ax_registered_1d.set_ylabel('counts/s')
        self._ax_registered_1d.set_yscale('log')

        # colorbars
        if use_fgivenx:
            self._ax_incident_cb = fig.add_subplot(gs_top[0,1])
            self._axes.append(self._ax_incident_cb)
            self._ax_registered_1d_cb = fig.add_subplot(gs_bottom[1,1])
            self._axes.append(self._ax_registered_1d_cb)

        self._ax_registered_cb = fig.add_subplot(gs_bottom[0,1])
        self._axes.append(self._ax_registered_cb)

        self._show_components = show_components
        self._registered_cmap = registered_cmap

        if sample_line_kwargs is not None:
            self._sample_line_kwargs = sample_line_kwargs
        else:
            self._sample_line_kwargs = {}

        self._plot_truth = plot_truth

        if self._plot_truth:
            if truth_line_kwargs is None:
                self._truth_line_kwargs = \
                        dict(color=('b' if self._use_fgivenx else 'darkgreen'),
                             ls='-.',
                             lw=1.0,
                             alpha=1.0)
            else:
                self._truth_line_kwargs = truth_line_kwargs

            if comp_truth_line_kwargs is not None:
                self._comp_truth_line_kwargs = comp_truth_line_kwargs
            else:
                self._comp_truth_line_kwargs = self._truth_line_kwargs

        if expectation_line_kwargs is None:
            color = 'k' if self._use_fgivenx else 'r'
            self._expectation_line_kwargs = dict(color=color,
                                                 ls='-',
                                                 lw=1.0,
                                                 alpha=1.0)
        else:
            self._expectation_line_kwargs = expectation_line_kwargs

        if comp_expectation_line_kwargs is not None:
            self._comp_expectation_line_kwargs = comp_expectation_line_kwargs
        else:
            self._comp_expectation_line_kwargs = self._expectation_line_kwargs

        if background_line_kwargs is None:
            self._background_line_kwargs = dict(color='orange',
                                                 ls='-',
                                                 lw=1.0,
                                                 alpha=1.0)
        else:
            self._background_line_kwargs = background_line_kwargs

        plt.close()

        yield

    @property
    def instruction_set(self):
        return self._instruction_set

    @instruction_set.deleter
    def instruction_set(self):
        try:
            del self._instruction_set
        except AttributeError:
            pass

    @make_verbose('SpectrumPlot object iterating over samples',
                  'SpectrumPlot object finished iterating')
    def execute(self, thetas, wrapper):
        self._num_samples = thetas.shape[0]

        self._energies = self._signal.create_energy_array(self._rel_num_energies)

        if self._use_fgivenx:
            # determine which spectrum to compute contours for
            if self._show_attenuated:
                wrapped = wrapper(self, 'incident_sums')
                self._instruction_set = 1
                # calculate expected unattenuated incident
                for i in range(self._num_samples):
                    _ = wrapped(None, thetas[i,:])

                # rewrap to reset the cache iterator
                wrapped = wrapper(self, 'attenuated_incident_sums')
                self._instruction_set = 0
            else:
                wrapped = wrapper(self, 'incident_sums')
                self._instruction_set = 1

            self._add_incident_contours(wrapped, thetas)

            self._instruction_set = 2
            self._add_registered_contours(wrapper(self, 'registered_sums'),
                                          thetas)

            yield 'Added conditional posterior contours for incident spectrum.'

            if self._add_background:
                self._instruction_set = 3
                wrapped = wrapper(self, 'background_sum')
                for i in range(self._num_samples):
                    wrapped(None, thetas[i,:])

        else:
            del self.instruction_set
            wrapped = wrapper(self, ['attenuated_incident_sums',
                                     'incident_sums',
                                     'registered_sums'])
            for i in range(self._num_samples):
                wrapped(None, thetas[i,:])

        yield

    def __next__(self):
        """ Update posterior expected signals given the updated signal.

        Plots signals if :mod:`fgivenx` is not used, otherwise returns
        callback information for :mod:`fgivenx`.

        .. note::

            You cannot make an iterator from an instance of this class.

        """
        try:
            self._instruction_set
        except AttributeError:
            if self._show_attenuated:
                signal = self._handle_attenuated_incident()
                self._add_signal(self._ax_incident,
                                 self._energies,
                                 signal,
                                 **self._sample_line_kwargs)

            signal = self._handle_incident()

            if not self._show_attenuated:
                self._add_signal(self._ax_incident,
                                 self._energies,
                                 signal,
                                 **self._sample_line_kwargs)

            signal = self._handle_registered()
            self._add_registered_spectrum(self._ax_registered_1d,
                                          signal,
                                          **self._sample_line_kwargs)

            if self._add_background:
                self._handle_background() # nothing to plot here
        else:
            if self._instruction_set == 0:
                return self._handle_attenuated_incident() # end execution here

            if self._instruction_set == 1:
                return self._handle_incident()

            if self._instruction_set == 2:
                return self._handle_registered()

            if self._instruction_set == 3:
                self._handle_background() # nothing to return

        return None # reached if not invoking fgivenx

    def _handle_attenuated_incident(self):
        """ Instructions for handling the attenuated incident spectrum. """
        ref = self._signal
        try:
            self._attenuated_incident_sums
        except AttributeError:
            self._attenuated_incident_sums = [None]*len(ref.incident_specific_flux_signals)

        attenuated_incident = None
        for i, component in enumerate(ref.incident_specific_flux_signals):
            temp = phase_integrator(1.0,
                                    _np.array([0.0,1.0]),
                                    component,
                                    ref.phases[i],
                                    0.0)

            temp = energy_interpolator(1, # threads
                                       temp,
                                       _np.log10(ref.energies),
                                       _np.log10(self._energies)).reshape(-1)

            # requires a user implementation or will raise NotImplementedError
            ref.interstellar(self._energies, temp)

            try:
                attenuated_incident += temp
            except TypeError:
                attenuated_incident = temp

            try:
                self._attenuated_incident_sums[i] += temp
            except TypeError:
                self._attenuated_incident_sums[i] = temp.copy()

        return attenuated_incident

    @property
    def attenuated_incident_sums(self):
        return self._attenuated_incident_sums

    @attenuated_incident_sums.deleter
    def attenuated_incident_sums(self):
        try:
            del self._attenuated_incident_sums
        except AttributeError:
            pass

    @property
    def expected_attenuated_incident(self):
        """ Get the expectations of the incident (component) spectra. """
        return [component/self._num_samples for component \
                                    in self._attenuated_incident_sums]

    def _handle_incident(self):
        """ Instructions for handling the unattenuated incident spectrum. """
        ref = self._signal
        try:
            self._incident_sums
        except AttributeError:
            self._incident_sums = [None] * len(ref.incident_specific_flux_signals)

        incident = None
        for i, component in enumerate(ref.incident_specific_flux_signals):
            temp = phase_integrator(1.0,
                                    _np.array([0.0,1.0]),
                                    component,
                                    ref.phases[i],
                                    0.0)

            temp = energy_interpolator(1, # threads
                                       temp,
                                       _np.log10(ref.energies),
                                       _np.log10(self._energies)).reshape(-1)

            try:
                incident += temp
            except TypeError:
                incident = temp

            try:
                self._incident_sums[i] += temp
            except TypeError:
                self._incident_sums[i] = temp.copy()

        return incident

    @property
    def incident_sums(self):
        return self._incident_sums

    @incident_sums.deleter
    def incident_sums(self):
        try:
            del self._incident_sums
        except AttributeError:
            pass

    @property
    def expected_incident(self):
        """ Get the expectations of the incident (component) spectra. """
        return [component/self._num_samples for component in self._incident_sums]

    def _handle_registered(self):
        """ Instructions for handling the registered spectrum. """
        ref = self._signal
        try:
            self._registered_sums
        except AttributeError:
            self._registered_sums = [None] * len(ref.signals)

        registered = None
        for i, component in enumerate(ref.signals):
            temp = phase_integrator(1.0,
                                    _np.array([0.0,1.0]),
                                    component,
                                    ref.phases[i],
                                    0.0).reshape(-1)

            try:
                registered += temp
            except TypeError:
                registered = temp

            try:
                self._registered_sums[i] += component
            except TypeError:
                self._registered_sums[i] = component

        return registered

    @property
    def registered_sums(self):
        return self._registered_sums

    @registered_sums.deleter
    def registered_sums(self):
        try:
            del self._registered_sums
        except AttributeError:
            pass

    @property
    def expected_registered(self):
        """ Get the expectations of the registered count-rate spectra. """
        return [component/self._num_samples for component in self._registered_sums]

    def _handle_background(self):
        """ Instructions for handling the background spectrum. """
        ref = self._signal
        try:
            self._background_sum
        except AttributeError:
            self._background_sum = None

        background = ref.background_signal

        try:
            self._background_sum += background
        except TypeError:
            self._background_sum = background

        return None

    @property
    def background_sum(self):
        return self._background_sum

    @background_sum.deleter
    def background_sum(self):
        try:
            del self._background_sum
        except AttributeError:
            pass

    @property
    def expected_background(self):
        """ Get the expectation of the background count-rate spectrum. """
        mean = self._background_sum/self._num_samples
        # transform to count rate signal, assuming background signal
        # is in units of counts; subclass and overwrite if you want to change
        return mean/self._signal.data.exposure_time

    @make_verbose('SpectrumPlot object finalizing',
                  'SpectrumPlot object finalized')
    def finalize(self):
        """ Execute final instructions. """
        ref = self._signal
        self._plot_components = self._show_components and ref.num_components > 1

        # add the incident signals
        if self._plot_truth:
            self._add_true_incident_signals()

        self._add_expected_incident_signals()

        # add the registered signals
        if self._plot_truth:
            self._add_true_registered_signals()

        self._add_expected_registered_signals()

    def _add_true_incident_signals(self):
        """ Render ground truth incident (component) signals. """
        ref = self._signal

        total = None
        for component, phases  in zip(ref.incident_specific_flux_signals,
                                      ref.phases):
            temp = phase_integrator(1.0,
                                    _np.array([0.0,1.0]),
                                    component,
                                    phases,
                                    0.0)

            temp = energy_interpolator(1, # threads
                                       temp,
                                       _np.log10(ref.energies),
                                       _np.log10(self._energies)).reshape(-1)

            if self._show_attenuated:
                ref.interstellar(self._energies, temp)

            try:
                total += temp
            except TypeError:
                total = temp

            if self._plot_components:
                self._add_signal(self._ax_incident,
                                 self._energies,
                                 temp,
                                 **self._comp_truth_line_kwargs)

        self._add_signal(self._ax_incident,
                         self._energies,
                         total,
                         **self._truth_line_kwargs)

    def _add_expected_incident_signals(self):
        """ Render posterior-expected incident (component) signals. """
        ax = self._ax_incident
        view_y_bottom = ax.yaxis.get_view_interval()[0]

        ref = self._signal

        total = None
        for component in self.expected_incident:
            try:
                total += component
            except TypeError:
                total = component

            if self._plot_components:
                self._add_signal(ax,
                                 self._energies,
                                 component,
                                 **self._comp_expectation_line_kwargs)

        self._add_signal(ax,
                         self._energies,
                         total,
                         **self._expectation_line_kwargs)

        if self._show_attenuated:
            total = None
            for component in self.expected_attenuated_incident:
                try:
                    total += component
                except TypeError:
                    total = component

                if self._plot_components:
                    self._add_signal(ax,
                                     self._energies,
                                     component,
                                     **self._comp_expectation_line_kwargs)

            self._add_signal(ax,
                             self._energies,
                             total,
                             **self._expectation_line_kwargs)

        ax.set_ylim(bottom = view_y_bottom)

        ax.set_xlim([self._signal.energy_edges[0],
                     self._signal.energy_edges[-1]])
        locmaj = LogLocator(base=10.0, numticks=100)
        ax.yaxis.set_major_locator(locmaj)

        locmin = LogLocator(base=10.0, subs=_np.arange(2,10)*0.1, numticks=100)
        ax.yaxis.set_minor_locator(locmin)
        ax.yaxis.set_minor_formatter(NullFormatter())

    def _add_true_registered_signals(self):
        """ Render ground truth registered (component) signals. """
        ref = self._signal#s

        total = None
        for component, phases in zip(ref.signals, ref.phases):
            temp = phase_integrator(1.0,
                                    _np.array([0.0,1.0]),
                                    component,
                                    phases,
                                    0.0).reshape(-1)
            try:
                total += temp
            except TypeError:
                total = temp

            if self._plot_components:
                self._add_registered_spectrum(self._ax_registered_1d,
                                              temp,
                                              **self._comp_truth_line_kwargs)

        self._add_registered_spectrum(self._ax_registered_1d,
                                      total,
                                      **self._truth_line_kwargs)

        if self._add_background: # another line including the background
            total += ref.background_signal
            self._add_registered_spectrum(self._ax_registered_1d,
                                          total,
                                          **self._background_line_kwargs)

    def _add_expected_registered_signals(self):
        """ Render posterior-expected registered (component) signals. """
        ref = self._signal

        total = None
        for component, phases in zip(self.expected_registered, ref.phases):
            temp = phase_integrator(1.0,
                                    _np.array([0.0,1.0]),
                                    component,
                                    phases,
                                    0.0).reshape(-1)
            try:
                total += temp
            except TypeError:
                total = temp

            if self._plot_components:
                self._add_registered_spectrum(self._ax_registered_1d,
                                              temp,
                                              **self._comp_expectation_line_kwargs)

        # 1D
        self._add_registered_spectrum(self._ax_registered_1d,
                                      total,
                                      **self._expectation_line_kwargs)

        if self._add_background: # another line including the background
            total += self.expected_background
            self._add_registered_spectrum(self._ax_registered_1d,
                                          total,
                                          **self._background_line_kwargs)

        ax = self._ax_registered_1d
        ax.set_xlim([ref.data.channels[0],
                     ref.data.channels[-1]])
        locmaj = LogLocator(base=10.0, numticks=100)
        ax.yaxis.set_major_locator(locmaj)

        locmin = LogLocator(base=10.0, subs=_np.arange(2,10)*0.1, numticks=100)
        ax.yaxis.set_minor_locator(locmin)
        ax.yaxis.set_minor_formatter(NullFormatter())

        # 2D
        total = None
        for component, shift, phases in zip(self.expected_registered,
                                            ref.shifts, ref.phases):
            temp = interp(self._phases,
                          phases,
                          component,
                          shift)
            try:
                total += temp
            except TypeError:
                total = temp

        if self._add_background: # add expectated background
            for i in range(total.shape[1]):
                total[:,i] += self.expected_background

        registered = self._ax_registered.pcolormesh(ref.data.channels,
                                    self._phases,
                                    total.T, # channel number as x-axis
                                    cmap = cm.get_cmap(self._registered_cmap),
                                    linewidth = 0,
                                    rasterized = self._rasterized)

        registered.set_edgecolor('face')
        self._ax_registered.set_xlim([ref.data.channels[0],
                                      ref.data.channels[-1]])

        self._registered_cb = plt.colorbar(registered,
                                           cax=self._ax_registered_cb,
                                           ticks=_get_default_locator(None),
                                           format=_get_default_formatter())
        self._registered_cb.ax.set_frame_on(True)
        self._registered_cb.ax.yaxis.set_minor_locator(AutoMinorLocator())

        self._registered_cb.set_label(label=r'counts/s', labelpad=15)

    def _add_registered_spectrum(self, ax, spectrum, **kwargs):
        """ Add registered spectrum line as a function of channel number. """
        if not kwargs:
            kwargs.update(dict(color='k', linestyle='-', lw=0.05, alpha=1.0))
        elif 'ls' in kwargs:
            kwargs['linestyle'] = kwargs.pop('ls')

        ax.step(self._signal.data.channels,
                spectrum,
                where='mid',
                **kwargs)

    @make_verbose('Adding credible intervals on the incident specific photon '
                  'flux spectrum',
                  'Credible intervals added')
    def _add_incident_contours(self, callback, thetas):
        """ Add contours to 1D incident specific flux spectrum axes objects. """
        self._add_contours(callback, thetas, self._energies,
                           self._ax_incident,
                           self._ax_incident_cb,
                           **self._incident_contour_kwargs)
        label = r'$\pi(\mathrm{photons/keV/cm}^{2}\mathrm{/s};E)$'
        self._ax_incident_cb.set_ylabel(label)

    @make_verbose('Adding credible intervals on the count-rate spectrum',
                  'Credible intervals added')
    def _add_registered_contours(self, callback, thetas):
        """ Add contours to 1D count-rate spectrum axes objects. """
        self._add_contours(callback, thetas, self._signal.data.channels,
                           self._ax_registered_1d,
                           self._ax_registered_1d_cb,
                           **self._registered_contour_kwargs)
        self._ax_registered_1d_cb.set_ylabel(r'$\pi(\mathrm{counts/s};\mathrm{channel})$')
