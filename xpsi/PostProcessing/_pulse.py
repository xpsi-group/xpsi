from ._global_imports import *

try:
    import fgivenx
except ImportError:
    _warning('Cannot import fgivenx for conditional posterior contours.')
    fgivenx = None

from ..tools.phase_interpolator import phase_interpolator as interp

from ._signalplot import SignalPlot

class PulsePlot(SignalPlot):
    """ Plot posterior-averaged channel-summed count-rate pulse profiles.

    The figure contains four panels which share phase as an x-axis:

        * the first (topmost) panel displays the posterior expectation of
          the specific photon flux signal from the source, jointly resolved in
          energy and phase;
        * the second panel displays the energy-integrated photon specific flux
          signal as function of phase for a subset of samples, optionally using
          :mod:`fgivenx`;
        * the third panel displays the posterior expectation of the count-rate
          signal as a function of channel and phase;
        * the last (bottommost) panel displays the channel-summed pulse as
          a function of phase for a subset of samples, optionally using
          :mod:`fgivenx`.

    The second and last panels aim to render the conditional posterior
    distribution of the associated signal as a function phase, ideally with
    contours to map out the conditional posterior mass. These panels have space
    to optionally display other elements such as: the posterior-expected total
    signal; the posterior-expected component signals; and the true total and
    component signals if the ground truth (the injected signal correpsonding
    to some model parameter vector) is known.

    The following example is from :ref:`R19`:

    .. image:: _static/_pulseplot.png

    :param int num_phases:
        The number of phases to interpolate the pulse-profile signals at.

    :param str incident_cmap:
         Colormap name from :mod:`matplotlib` to use for the posterior-expected
         incident signal as a function of energy and phase (top panel).

    :param str registered_cmap:
         Colormap name from :mod:`matplotlib` to use for the posterior-expected
         registered signal as a function of channel and phase (third panel).

    :param bool show_components:
        If the :class:`~.Signal.Signal` instance has multiple components,
        display the posterior expectations of those components as a
        function of phase (second and last panels).

    :param dict expectation_line_kwargs:
        Keyword arguments for plotting the posterior-expected signal lines (in
        the second and last panels).

    :param bool use_fgivenx:
        Use :mod:`fgivenx` to plot conditional posterior contours in the
        second and last panels?

    :param dict incident_contour_kwargs:
        Keyword arguments for :mod:`fgivenx` incident signal contours (second
        panel) that will take precedence over the corresponding class
        attributes. (See the :class:`~.SignalPlot` class if you choose
        not to modify these attributes on this present subclass.)

    :param dict registered_contour_kwargs:
        Keyword arguments for :mod:`fgivenx` registered signal contours (last
        panel) that will take precedence over the corresponding class
        attributes. (See the :class:`~.SignalPlot` class if you choose
        not to modify these attributes on this present subclass.)

    :param plot_truth:
        Plot the ground truth (injected) signal, if known and available, in
        the second and last panels.

    :param truth_line_kwargs:
        Keyword arguments for plotting the ground truth signal lines (second
        and last panels).

    """

    __figtype__ = 'signalplot_pulse'

    # do not change at runtime (see base class comments):
    __caching_targets__ = ['shifts', # hot region phase shifts
                           'signals', # count-rate signals
                           # incident flux signals integrated over energy bins
                           'incident_flux_signals']

    __rows__          = 4
    __columns__       = 1
    __ax_rows__       = 4
    __ax_columns__    = 2
    __height_ratios__ = [1]*4
    __width_ratios__  = [50, 1] # second column for colorbars
    __wspace__        = 0.025
    __hspace__        = 0.125

    @make_verbose('Instantiating a pulse-profile plotter for posterior checking',
                  'Pulse-profile plotter instantiated')
    def __init__(self,
                 num_phases=1000,
                 incident_cmap='inferno',
                 registered_cmap='inferno',
                 show_components=False,
                 expectation_line_kwargs=None,
                 comp_expectation_line_kwargs=None,
                 sample_line_kwargs=None,
                 use_fgivenx=False,
                 incident_contour_kwargs=None,
                 registered_contour_kwargs=None,
                 plot_truth=False,
                 truth_line_kwargs=None,
                 comp_truth_line_kwargs=None,
                 **kwargs):
        super(PulsePlot, self).__init__(**kwargs)

        self._phases = _np.linspace(0.0, 2.0, int(num_phases))

        if use_fgivenx and fgivenx is None:
            raise ImportError('Install fgivenx to plot contours.')

        self._use_fgivenx = use_fgivenx
        if self._use_fgivenx:
            self._incident_contour_kwargs =\
                    incident_contour_kwargs if incident_contour_kwargs else {}
            self._registered_contour_kwargs =\
                registered_contour_kwargs if registered_contour_kwargs else {}

        self._get_figure()

        # for the incident specific flux signal incident on an instrument
        self._ax = self._add_subplot(0,0)
        self._ax_1d = self._add_subplot(1,0)

        # for the count-rate signal registered by an instrument
        self._ax_registered = self._add_subplot(2,0)
        self._ax_registered_1d = self._add_subplot(3,0)

        # properties for phase each axis
        for ax in self._axes:
            ax.xaxis.set_major_locator(MultipleLocator(0.2))
            ax.xaxis.set_minor_locator(MultipleLocator(0.05))
            ax.set_xlim([0.0,2.0])
            if ax is not self._ax_registered_1d:
                ax.tick_params(axis='x', labelbottom=False)
            else:
                ax.set_xlabel('$\phi$ [cycles]')

        self._ax.set_ylabel(r'$E$ [keV]')
        self._ax_1d.set_ylabel(r'photons/cm$^{2}$/s')
        self._ax_registered.set_ylabel('channel')
        self._ax_registered_1d.set_ylabel('counts/s')

        # colorbars
        self._ax_cb = self._add_subplot(0,1)
        self._ax_registered_cb = self._add_subplot(2,1)
        if self._use_fgivenx:
            self._ax_1d_cb = self._add_subplot(1,1)
            self._ax_registered_1d_cb = self._add_subplot(3,1)

        self._show_components = show_components
        self._incident_cmap = incident_cmap
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

        plt.close()

    @property
    def instruction_set(self):
        return self._instruction_set

    @instruction_set.deleter
    def instruction_set(self):
        try:
            del self._instruction_set
        except AttributeError:
            pass

    @make_verbose('PulsePlot object iterating over samples',
                  'PulsePlot object finished iterating')
    def execute(self, thetas, wrapper):
        """ Loop over posterior samples. """
        self._num_samples = thetas.shape[0]
        if self._use_fgivenx:
            self._instruction_set = 0
            self._add_incident_contours(wrapper(self, 'incident_sums'),
                                        thetas)

            yield 'Added conditional posterior contours for incident signal.'

            self._instruction_set = 1
            self._add_registered_contours(wrapper(self, 'registered_sums'),
                                          thetas)

            yield 'Added conditional posterior contours for registered signal.'
        else: # iterate manually instead of driven by fgivenx
            del self.instruction_set
            wrapped = wrapper(self, ['incident_sums', 'registered_sums'])
            for i in range(self._num_samples):
                wrapped(None, thetas[i,:])

        yield

    def __next__(self):
        """ Update posterior expected signals given the updated signal object.

        Plots signals if :mod:`fgivenx` is not used, otherwise returns
        callback information for :mod:`fgivenx`.

        .. note::

            You cannot make an iterator from an instance of this class.

        """
        try:
            self._instruction_set
        except AttributeError:
            incident = self._handle_incident()
            self._add_signal(self._ax_1d,
                             self._phases,
                             incident,
                             **self._sample_line_kwargs)

            registered = self._handle_registered()
            self._add_signal(self._ax_registered_1d,
                             self._phases,
                             registered,
                             **self._sample_line_kwargs)
        else:
            if self._instruction_set == 0:
                return self._handle_incident() # end execution here

            if self._instruction_set == 1:
                return self._handle_registered()

        return None # reached if not invoking fgivenx

    def _handle_incident(self):
        """ Instructions for handling the incident signal. """
        ref = self._signal
        try:
            self._incident_sums
        except AttributeError:
            self._incident_sums = [None]
            self._incident_sums *= len(ref.incident_flux_signals)

        incident = None
        for i, component in enumerate(ref.incident_flux_signals):
            temp = interp(self._phases,
                          ref.phases[i],
                          component,
                          ref.shifts[i])

            try:
                incident += temp
            except TypeError:
                incident = temp

            try:
                self._incident_sums[i] += temp
            except TypeError:
                self._incident_sums[i] = temp.copy()

        return _np.sum(incident, axis=0)

    @property
    def incident_sums(self):
        return self._incident_sums

    @incident_sums.deleter
    def incident_sums(self):
        del self._incident_sums

    @property
    def expected_incident(self):
        """ Get the expectations of the component incident signals. """
        return [component/self._num_samples for component in self._incident_sums]

    def _handle_registered(self):
        """ Instructions for handling the registered signal. """
        ref = self._signal
        try:
            self._registered_sums
        except AttributeError:
            self._registered_sums = [None] * len(ref.signals)

        registered = None
        for i, (signal, shift) in enumerate(zip(ref.signals, ref.shifts)):
            temp = interp(self._phases,
                          ref.phases[i],
                          signal, shift)

            try:
                registered += temp
            except TypeError:
                registered = temp

            try:
                self._registered_sums[i] += temp
            except TypeError:
                self._registered_sums[i] = temp.copy()

        return _np.sum(registered, axis=0)

    @property
    def registered_sums(self):
        return self._registered_sums

    @registered_sums.deleter
    def registered_sums(self):
        del self._registered_sums

    @property
    def expected_registered(self):
        """ Get the expectations of the component registered signals. """
        return [component/self._num_samples for component in self._registered_sums]

    @make_verbose('PulsePlot object finalizing',
                  'PulsePlot object finalized')
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
        for component, shift, phases in zip(ref.incident_flux_signals,
                                            ref.shifts,
                                            ref.phases):
            temp = interp(self._phases,
                          phases,
                          component,
                          shift)

            try:
                total += temp
            except TypeError:
                total = temp

            if self._plot_components:
                self._add_signal(self._ax_1d,
                                 self._phases,
                                 temp,
                                 axis=0,
                                 **self._comp_truth_line_kwargs)

        self._add_signal(self._ax_1d,
                         self._phases,
                         total,
                         axis=0,
                         **self._truth_line_kwargs)

    def _add_expected_incident_signals(self):
        """ Render posterior-expected incident (component) signals. """
        ref = self._signal

        total = None
        for component in self.expected_incident:
            try:
                total += component
            except TypeError:
                total = component

            if self._plot_components:
                self._add_signal(self._ax_1d,
                                 self._phases,
                                 component,
                                 axis=0,
                                 **self._comp_expectation_line_kwargs)
        # 1D
        self._add_signal(self._ax_1d,
                         self._phases,
                         total,
                         axis=0,
                         **self._expectation_line_kwargs)

        self._ax_1d.yaxis.set_major_locator(_get_default_locator(None))
        self._ax_1d.yaxis.set_major_formatter(_get_default_formatter())
        self._ax_1d.yaxis.set_minor_locator(AutoMinorLocator())

        # 2D
        Delta_E = ref.energy_edges[1:] - ref.energy_edges[:-1]
        for i in range(total.shape[1]):
            total[:,i] /= Delta_E # mean specific flux in each interval

        incident = self._ax.pcolormesh(self._phases,
                                       ref.energy_edges,
                                       total,
                                       cmap = cm.get_cmap(self._incident_cmap),
                                       linewidth = 0,
                                       rasterized = self._rasterized)

        incident.set_edgecolor('face')
        self._ax.set_ylim([ref.energy_edges[0],
                           ref.energy_edges[-1]])
        self._ax.set_yscale('log')

        self._incident_cb = plt.colorbar(incident, cax=self._ax_cb,
                                         ticks=_get_default_locator(None),
                                         format=_get_default_formatter())
        self._incident_cb.ax.set_frame_on(True)
        self._incident_cb.ax.yaxis.set_minor_locator(AutoMinorLocator())

        self._incident_cb.set_label(label=r'photons/keV/cm$^{2}$/s',
                                    labelpad=15)

    def _add_true_registered_signals(self):
        """ Render ground truth registered (component) signals. """
        ref = self._signal

        total = None
        for component, shift, phases in zip(ref.signals, ref.shifts, ref.phases):
            temp = interp(self._phases,
                          phases,
                          eomponent,
                          shift)

            try:
                total += temp
            except TypeError:
                total = temp

            if self._plot_components:
                self._add_signal(self._ax_registered_1d,
                                 self._phases,
                                 temp,
                                 axis=0,
                                 **self._truth_line_kwargs)

        self._add_signal(self._ax_registered_1d,
                         self._phases,
                         total,
                         axis=0,
                         **self._truth_line_kwargs)

    def _add_expected_registered_signals(self):
        """ Render posterior-expected registered (component) signals. """
        ref = self._signal

        total = None
        for component in self.expected_registered:
            try:
                total += component
            except TypeError:
                total = component

            if self._plot_components:
                self._add_signal(self._ax_registered_1d,
                                 self._phases,
                                 component,
                                 axis=0,
                                 **self._comp_expectation_line_kwargs)

        # 1D
        self._add_signal(self._ax_registered_1d,
                         self._phases,
                         total,
                         axis=0,
                         **self._expectation_line_kwargs)

        self._ax_registered_1d.yaxis.set_major_locator(_get_default_locator(None))
        self._ax_registered_1d.yaxis.set_major_formatter(_get_default_formatter())
        self._ax_registered_1d.yaxis.set_minor_locator(AutoMinorLocator())

        # 2D
        registered = self._ax_registered.pcolormesh(self._phases,
                                      ref.data.channels,
                                      total,
                                      cmap = cm.get_cmap(self._registered_cmap),
                                      linewidth = 0,
                                      rasterized = self._rasterized)

        registered.set_edgecolor('face')
        self._ax_registered.set_ylim([ref.data.channels[0],
                                      ref.data.channels[-1]])
        self._ax_registered.set_yscale('log')

        self._registered_cb = plt.colorbar(registered,
                                           cax=self._ax_registered_cb,
                                           ticks=_get_default_locator(None),
                                           format=_get_default_formatter())
        self._registered_cb.ax.set_frame_on(True)
        self._registered_cb.ax.yaxis.set_minor_locator(AutoMinorLocator())

        self._registered_cb.set_label(label=r'counts/s', labelpad=15)

    @make_verbose('Adding credible intervals on the incident photon flux '
                  'signal as function of phase',
                  'Credible intervals added')
    def _add_incident_contours(self, callback, thetas):
        """ Add contours to 1D incident photon flux signal axes objects. """
        self._add_contours(callback, thetas, self._phases,
                           self._ax_1d, self._ax_1d_cb,
                           **self._incident_contour_kwargs)
        label = r'$\pi(\mathrm{photons/cm}^{2}\mathrm{/s};\phi)$'
        self._ax_1d_cb.set_ylabel(label)

    @make_verbose('Adding credible intervals on the count-rate '
                  'signal as function of phase',
                  'Credible intervals added')
    def _add_registered_contours(self, callback, thetas):
        """ Add contours to 1D count-rate signal axes objects. """
        self._add_contours(callback, thetas, self._phases,
                           self._ax_registered_1d, self._ax_registered_1d_cb,
                           **self._registered_contour_kwargs)
        self._ax_registered_1d_cb.set_ylabel(r'$\pi(\mathrm{counts/s};\phi)$')
