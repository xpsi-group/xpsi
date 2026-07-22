from pylab import *
from ._global_imports import *
from scipy.ndimage import measurements, gaussian_filter
from scipy import stats as _stats
from ._signalplot import SignalPlot
from matplotlib.patches import Patch as Patch
from matplotlib.lines import Line2D as Line2D

# Optional heavy dependency for Moran's I – graceful fallback if absent
try:
    from libpysal.weights import lat2W as _lat2W
    from esda.moran import Moran as _Moran
    _HAS_ESDA = True
except ImportError:
    _HAS_ESDA = False

# Optional statsmodels for Ljung-Box test
try:
    from statsmodels.stats.diagnostic import acorr_ljungbox as _acorr_ljungbox
    _HAS_STATSMODELS = True
except ImportError:
    _HAS_STATSMODELS = False

class ResidualPlot(SignalPlot):
    """ Plot the count data, the posterior-expected count signal, residuals, clusters and their related distributions.

    The figure contains three upper panels which share phase as an x-axis:

        * the top panel displays the data count numbers in joint channel-phase
          intervals, identically split over two rotational phase cycles;
        * the center panel displays the posterior-expected count signal
          over joint channel-phase intervals;
        * the bottom panel displays the standardised residuals between the
          data and posterior-expected count signal over joint channel-phase
          intervals.

    If requested by plot_clusters, the figure also contains 2 lower panels
        * the first panel displays two figures: the one on the left shows the
          standardised residuals for which their absolute values overpass a chosen
          threshold, and the one on the right present the estimated cluster sizes from
          these residuals.
        * the second panel displays two figures: the one on the left is the residual distribution
          compared with an optimal gaussian one, and the one on the right shows a distribution of
          cluster sizes from residuals which have absolute values above the chosen threshold

    The following example is (improved) from `Riley et al. 2019 <https://ui.adsabs.harvard.edu/abs/2019ApJ...887L..21R/abstract>`_:

    .. image:: _static/_residualplot.png

    :param str data_cmap:
         Colormap name from :mod:`matplotlib` to use for the data count numbers
         over joint channel-phase intervals.

    :param str model_cmap:
         Colormap name from :mod:`matplotlib` to use for the posterior-expected
         count numbers over joint channel-phase intervals.

    :param str residual_cmap:
         Colormap name from :mod:`matplotlib` to use for the residuals between
         the data and posterior-expected count numbers over joint
         channel-phase intervals. A diverging colormap is recommended.

    """

    __figtype__ = 'signalplot_residuals'

    # do not change at runtime (see base class comments):
    __caching_targets__ = ['expected_counts']

    __rows__          = 3
    __columns__       = 1
    __ax_rows__       = 3
    __ax_columns__    = 2
    __height_ratios__ = [1]*3
    __width_ratios__  = [50, 1] # second column for colorbars
    __wspace__        = 0.025
    __hspace__        = 0.125

    @make_verbose('Instantiating a residual plotter for posterior checking',
                  'Residual plotter instantiated')
    def __init__(self,
                 data_cmap='inferno',
                 model_cmap='inferno',
                 residual_cmap='PuOr',
                 parameters_vector=None,
                 plot_pulse=False,
                 legend_pulse_coords=None,
                 n_realisation_per_sample=1,
                 blur_residuals=False,
                 plot_clusters=False,
                 threshlim=2.0,
                 clusters_cmap='PuOr',
                 clustdist_cmap='PuOr',
                 mu=0.0,
                 sigma=1.0,
                 nbins=50,
                 analyse_residual_structure=False,
                 structure_power_n_null=500,
                 structure_power_alpha=0.05,
                 structure_power_low_k_fraction=0.15,
                 structure_power_random_seed=0,
                 structure_variogram_max_lag=20,
                 structure_variogram_n_null=500,
                 structure_variogram_alpha=0.05,
                 structure_variogram_random_seed=0,
                 structure_acf2d_max_lag=8,
                 **kwargs):
        """
        Constructor method for plotting residuals.

        :param str data_cmap:
             Colormap name from :mod:`matplotlib` to use for the data count numbers
             over joint channel-phase intervals.

        :param str model_cmap:
             Colormap name from :mod:`matplotlib` to use for the posterior-expected
             count numbers over joint channel-phase intervals.

        :param str residual_cmap:
             Colormap name from :mod:`matplotlib` to use for the residuals between
             the data and posterior-expected count numbers over joint
             channel-phase intervals. A diverging colormap is recommended.

        :param list parameters_vector:
             List of model parameters to plot the residuals for.

        :param bool plot_pulse:
             Plot the pulse profile?

        :param int n_realisation_per_sample:
             Number of realisations to use to get the errorbars of the posterior predictive.

        :param tuple legend_pulse_coords:
             Coordinates for the legend of the pulse if plotted. If not given, will just use the best location.

        :param bool blur_residuals:
             Blur the residuals?

        :param bool plot_clusters:
             Plot cluster sizes and residual distribution?

        :param float threshlim:
             Threshold above which residuals are classified as clusters.

        :param str clusters_cmap:
             Colormap name from :mod:`matplotlib` to use for cluster sizes.

        :param str clustdist_cmap:
             Colormap name from :mod:`matplotlib` to use for cluster sizes distribution.

        :param float mu:
             Mean for the optimal gaussian distribution to compare with the residuals.

        :param float sigma:
             Standard deviation for the optimal gaussian distribution to compare with the residuals.

        :param bool analyse_residual_structure:
             If ``True``, run five complementary 2D structure tests on the
             standardised residuals after finalising the plot and print a
             summary table to stdout. The tests are:

             1. **2D Power Spectrum** – detects periodic/wave-like patterns
                via the FFT peak-to-median power ratio.
             2. **Moran's I** – global spatial autocorrelation statistic
                (requires ``esda`` and ``libpysal``; skipped if absent).
             3. **Variogram** – measures how correlation decays with
                pixel-distance on the phase–channel grid.
             4. **Runs test** – Wald–Wolfowitz test applied row- and
                column-wise to detect 1-D autocorrelation along each axis.
             5. **Ljung–Box + KS** – serial autocorrelation in the
                flattened residuals and a check for Gaussianity
                (requires ``statsmodels``; KS test always runs).

             The combined verdict (``STRUCTURED`` / ``RANDOM`` /
             ``INCONCLUSIVE``) is stored in
             :attr:`residual_structure_verdict` and the per-test results in
             :attr:`residual_structure_results`.

        :param int structure_power_n_null:
             Number of white-noise simulations used to calibrate the null
             distribution of the 2D power spectrum test.  Default is 500;
             the smallest reportable p-value is ``1/(n_null+1)``.

        :param float structure_power_alpha:
             Significance level of the 2D power spectrum test, before the
             internal Bonferroni factor of 2 over its two statistics.
             Default is 0.05.

        :param float structure_power_low_k_fraction:
             Radial-frequency cut (as a fraction of the Nyquist radius)
             defining the low-k region for the smooth-structure statistic
             of the 2D power spectrum test.  Default is 0.15.

        :param int structure_power_random_seed:
             Seed for the power spectrum null simulations (reproducibility).

        :param int structure_variogram_n_null:
             Number of white-noise simulations used to calibrate the null
             distribution of the variogram test.  Default is 500.

        :param float structure_variogram_alpha:
             Significance level of the variogram test, before the internal
             Bonferroni factor of 2 over its two statistics.  Default 0.05.

        :param int structure_variogram_random_seed:
             Seed for the variogram null simulations (reproducibility).

        :param int structure_variogram_max_lag:
             Maximum pixel-lag used for the variogram test.  Default is 20.

        :param int structure_acf2d_max_lag:
             Maximum pixel-lag used for the acf2d test.  Default is 8.

        :param kwargs:
             Keyword arguments for :class:`SignalPlot`.
        """
        super(ResidualPlot, self).__init__(**kwargs)

        # Setup the class
        self._data_cmap = data_cmap
        self._model_cmap = model_cmap
        self._residual_cmap = residual_cmap
        self._plot_clusters = plot_clusters
        self._plot_pulse = plot_pulse
        self._legend_pulse_coords = legend_pulse_coords
        self._n_realisations_per_model = n_realisation_per_sample
        self._blur_residuals = blur_residuals
        self._analyse_residual_structure = analyse_residual_structure
        self._structure_power_n_null = structure_power_n_null
        self._structure_power_alpha = structure_power_alpha
        self._structure_power_low_k_fraction = structure_power_low_k_fraction
        self._structure_power_random_seed = structure_power_random_seed
        self._structure_variogram_max_lag = structure_variogram_max_lag
        self._structure_variogram_n_null = structure_variogram_n_null
        self._structure_variogram_alpha = structure_variogram_alpha
        self._structure_variogram_random_seed = structure_variogram_random_seed
        self._structure_acf2d_max_lag = structure_acf2d_max_lag
        if not parameters_vector is None:
            self.parameters_vector = parameters_vector

        # Choose the rows for plotting
        self._data_row = 0
        self._model_row = 1 
        self._resid_row = 2

        # Do you want to plot the pulse profile?
        if self._plot_pulse:

            # Move all rows by one
            cls = type(self)
            cls.__rows__            += 1
            cls.__ax_rows__         += 1
            cls.__height_ratios__   = [1] * cls.__ax_rows__
            self._pulse_row = 0
            self._data_row = 1
            self._model_row = 2 
            self._resid_row = 3

        # Do you want to plot clusters ? 
        if self._plot_clusters:
            
            # Add more columns/rows
            cls = type(self)
            cls.__rows__          += 2
            cls.__ax_rows__       += 2
            cls.__ax_columns__    += 3
            cls.__height_ratios__ = [1] * cls.__ax_rows__
            cls.__width_ratios__  = [50, 1, 50, 1, 1] # second column for colorbars
            cls.__wspace__        = 0.1
            cls.__hspace__        = 0.35

            # Add parameters to plots
            self._threshlim = threshlim
            self._mu=mu
            self._sigma=sigma
            self._nbins=nbins
            self._clusters_cmap = clusters_cmap
            self._clustdist_cmap = clustdist_cmap

            # Get the rows
            self._cluster_row = self._resid_row + 1
            self._clustdist_row = self._resid_row + 2 

        # Generate the axes for plotting
        self._get_figure()

        self._ax_data = self._fig.add_subplot(self._gs[self._data_row,:-1])
        self._ax_data_cb = self._fig.add_subplot(self._gs[self._data_row,-1])

        self._ax_model = self._fig.add_subplot(self._gs[self._model_row,:-1])
        self._ax_model_cb = self._fig.add_subplot(self._gs[self._model_row,-1])

        self._ax_resid = self._fig.add_subplot(self._gs[self._resid_row,:-1])
        self._ax_resid_cb = self._fig.add_subplot(self._gs[self._resid_row,-1])

        # Prettify usual plots
        self._ax_resid.set_xlabel(r'$\phi$ [cycles]')
        for ax in (self._ax_data, self._ax_model):
            ax.tick_params(axis='x', labelbottom=False)

        for ax in (self._ax_data, self._ax_model, self._ax_resid):
            ax.set_ylabel('channel')
            ax.xaxis.set_major_locator(MultipleLocator(0.2))
            ax.xaxis.set_minor_locator(MultipleLocator(0.05))
            ax.set_xlim([0.0,2.0])

        # Plot the pulse if required
        if self._plot_pulse:
            self._ax_pulse = self._fig.add_subplot(self._gs[self._pulse_row,:-1])
            self._ax_pulse.set_ylabel(r'Total counts')
            self._ax_pulse.xaxis.set_major_locator(MultipleLocator(0.2))
            self._ax_pulse.xaxis.set_minor_locator(MultipleLocator(0.05))
            self._ax_pulse.tick_params(axis='x', labelbottom=False)
            self._ax_pulse.set_xlim([0.0,2.0])

        # Handle cluster plots if required
        if self._plot_clusters:

            # Generate axes
            self._ax_clres = self._fig.add_subplot(self._gs[self._cluster_row,:2])
            self._ax_clust = self._fig.add_subplot(self._gs[self._cluster_row,2:-1], sharex = self._ax_clres)
            self._ax_clust_cb = self._fig.add_subplot(self._gs[self._cluster_row,-1])
            
            self._ax_rdist = self._fig.add_subplot(self._gs[self._clustdist_row, :2])
            self._ax_cdist = self._fig.add_subplot(self._gs[self._clustdist_row, 2:-1])
            
            # Prettify
            for ax in (self._ax_data, self._ax_model, self._ax_clres, self._ax_clust):
                ax.set_xlabel(r'$\phi$ [cycles]')

            for ax in (self._ax_data, self._ax_model):
                ax.tick_params(axis='x', labelbottom=True)

            self._ax_clres.set_ylabel('channel')
            self._ax_clres.xaxis.set_major_locator(MultipleLocator(0.2))
            self._ax_clres.xaxis.set_minor_locator(MultipleLocator(0.05))
            self._ax_clres.set_xlim([0.0,1.0])
            self._ax_clres.set_title(r'|residuals| > {}'.format(self._threshlim))

            self._ax_clust.set_title('Cluster sizes') 
            self._ax_clust.set_yticklabels([])
            self._ax_clust.set_yticks([])
            
            self._ax_rdist.set_title('Residual distribution')  
            self._ax_rdist.xaxis.set_major_locator(MultipleLocator(1.0))
            self._ax_rdist.xaxis.set_minor_locator(MultipleLocator(0.5))
            self._ax_rdist.set_xlabel('Residuals')

            self._ax_cdist.set_title('Cluster sizes distribution')
            self._ax_cdist.set_xlabel('Cluster sizes')
            self._ax_cdist.tick_params(axis='y', labelright=True, labelleft=False)

        if "yscale" in kwargs:
            self.yscale = kwargs.get("yscale")
        else:
            self.yscale = "log"

        # Restore sizes
        cls = type(self)
        cls.__rows__          = 3
        cls.__columns__       = 1
        cls.__ax_rows__       = 3
        cls.__ax_columns__    = 2
        cls.__height_ratios__ = [1]*3
        cls.__width_ratios__  = [50, 1] # second column for colorbars
        cls.__wspace__        = 0.025
        cls.__hspace__        = 0.125

        plt.close()

    @make_verbose('ResidualPlot object iterating over samples',
                  'ResidualPlot object finished iterating')
    def execute(self, thetas, wrapper):
        """ Loop over posterior samples. """
        self._num_samples = thetas.shape[0]

        wrapped = wrapper(self, 'model_sum')
        for i in range(self._num_samples):
            wrapped(None, thetas[i,:])

    def __next__(self):
        """ Update posterior expected model given the updated signal.

        .. note::

            You cannot make an iterator from an instance of this class.

        """
        try:
            self._model_sum
            self._model_list
        except AttributeError:
            self._model_sum = self._signal.expected_counts.copy()
            self._model_list = [self._signal.expected_counts.copy()]
        else:
            self._model_sum += self._signal.expected_counts
            self._model_list.append( self._signal.expected_counts )

    @property
    def model_sum(self):
        """ Get the current posterior sum of the count numbers. """
        return self._model_sum

    @model_sum.deleter
    def model_sum(self):
        del self._model_sum

    @property
    def model_list(self):
        """ Get the current list of the count numbers. """
        return self._model_list

    @model_list.deleter
    def model_list(self):
        del self._model_list

    @property
    def expected_counts(self):
        """ Get the estimated posterior expectation of the count numbers. """
        return self._model_sum / self._num_samples

    @make_verbose('ResidualPlot object finalizing',
                  'ResidualPlot object finalized')
    def finalize(self):
        self._add_data()
        self._add_expected_counts()
        self._add_residuals()
        if self._plot_clusters:
            self._reveal_clusters()
            self._disp_distributions()
        if self._plot_pulse:
            self.add_pulses()
        if self._analyse_residual_structure:
            self._run_structure_analysis()

    def _set_vminmax(self):
        """ Compute minimum and maximum for data and model colorbars. """
        self._vmin = min(_np.min(self.expected_counts),
                         _np.min(self._signal.data.counts))
        self._vmax = max(_np.max(self.expected_counts),
                         _np.max(self._signal.data.counts))

    def add_pulses( self ):
        """ Display the pulse over the data if requested """

        try:
            self._vmin
        except AttributeError:
            self._set_vminmax()

        # Get the pulse of data and model
        doubles_phases = _np.concatenate( (self._signal.data.phases[:-1] , (self._signal.data.phases[:-1]+1.0)), axis=0 ) 
        pulse_data = self._signal.data.counts.sum( axis = 0 )
        double_pulse_data = _np.concatenate( (pulse_data , pulse_data), axis=0 )
        pulse_model = self.expected_counts.sum( axis = 0 )
        double_pulse_model = _np.concatenate( (pulse_model , pulse_model) )

        # Posterior predictive (PP) : get realizations
        pulse_model_list = [ model.sum( axis = 0 ) for model in self.model_list ]
        pulse_realizations = [ [_np.random.poisson( model ) for _ in range(self._n_realisations_per_model)] for model in pulse_model_list ]
        pulse_realizations_flat = _np.reshape( pulse_realizations , (len( pulse_realizations) * self._n_realisations_per_model,-1) )
        if len( pulse_realizations_flat ) >= 100:
            pulse_PP_q16 = _np.quantile( pulse_realizations_flat, q=0.16, axis=0 )
            pulse_PP_q84 = _np.quantile( pulse_realizations_flat, q=0.84, axis=0 )
            double_pulse_PP_q16 = np.concatenate( (pulse_PP_q16 , pulse_PP_q16), axis=0 )
            double_pulse_PP_q84 = np.concatenate( (pulse_PP_q84 , pulse_PP_q84), axis=0 )
        else:
            print( 'WARNING : Not enough realizations for posterior predictive errorbars. Using Poisson error of the model instead' )
            pulse_PP_errorbars = _np.sqrt( pulse_model )
            double_pulse_PP_errorbars = _np.concatenate( (pulse_PP_errorbars , pulse_PP_errorbars), axis=0 )
            double_pulse_PP_q16 = double_pulse_model - double_pulse_PP_errorbars
            double_pulse_PP_q84 = double_pulse_model + double_pulse_PP_errorbars

        self.bolometric_chi2 = np.sum( (pulse_data - pulse_model)**2 / pulse_model )

        # Plot pulse
        self._ax_pulse.step(x=doubles_phases,
                            y=double_pulse_data,
                            where='mid',
                            color='black')

        self._ax_pulse.step(x=doubles_phases,
                            y=double_pulse_model,
                            where='mid',
                            alpha=0.7,
                            color='blue')

        self._ax_pulse.fill_between(x=doubles_phases,
                                    y1=double_pulse_PP_q16,
                                    y2=double_pulse_PP_q84,
                                    step='mid',
                                    alpha=0.5,
                                    color='steelblue')

        lines = [Line2D([0], [0], color='black', lw=1), (Patch(color='steelblue', alpha=0.5, lw=1), Line2D([0], [0], color='blue', lw=1))]
        legends = ['Data', 'Model']
        if self._legend_pulse_coords is None:
            self._ax_pulse.legend(lines, legends, loc='best', frameon=False)
        else:
            self._ax_pulse.legend(lines, legends, bbox_to_anchor=self._legend_pulse_coords, frameon=False)

    def _add_data(self):
        """ Display data in topmost panel. """

        try:
            self._vmin
        except AttributeError:
            self._set_vminmax()

        #Calculate channel edges by averaging:
        channels = self._signal.data.channels
        channel_edges = _np.zeros((len(self._signal.data.channels)+1))
        channel_edges[1:len(channels)] = (channels[:len(channels)-1]+channels[1:])/2.0
        chandiff1 = (channels[1]-channels[0])/2.0
        chandiff2 = (channels[len(channels)-1]-channels[len(channels)-2])/2.0
        channel_edges[0] = channels[0]-chandiff1
        channel_edges[len(channels)] = channels[len(channels)-1]+chandiff2

        data = self._ax_data.pcolormesh(self._signal.data.phases,
                                        channel_edges,
                                        self._signal.data.counts,
                                        cmap = plt.get_cmap(self._data_cmap),
                                        vmin = self._vmin,
                                        vmax = self._vmax,
                                        linewidth = 0,
                                        rasterized = self._rasterized)
        data.set_edgecolor('face')

        data = self._ax_data.pcolormesh(self._signal.data.phases + 1.0,
                                        channel_edges,
                                        self._signal.data.counts,
                                        cmap = plt.get_cmap(self._data_cmap),
                                        vmin = self._vmin,
                                        vmax = self._vmax,
                                        linewidth = 0,
                                        rasterized = self._rasterized)
        data.set_edgecolor('face')

        self._ax_data.set_ylim([_np.max([channel_edges[0],0.001]),
                                channel_edges[-1]])
        self._ax_data.set_yscale(self.yscale)

        self._data_cb = self._fig.colorbar(data, cax=self._ax_data_cb,
                                     ticks=_get_default_locator(None),
                                     format=_get_default_formatter())
        self._data_cb.ax.set_frame_on(True)
        self._data_cb.ax.yaxis.set_minor_locator(AutoMinorLocator())
        self._data_cb.set_label(label=r'counts/cycle', labelpad=15)

    def _add_expected_counts(self):
        """ Display posterior expectation of model in second panel. """

        try:
            self._vmin
        except AttributeError:
            self._set_vminmax()

        #Calculate channel edges by averaging:
        channels = self._signal.data.channels
        channel_edges = _np.zeros((len(self._signal.data.channels)+1))
        channel_edges[1:len(channels)] = (channels[:len(channels)-1]+channels[1:])/2.0
        chandiff1 = (channels[1]-channels[0])/2.0
        chandiff2 = (channels[len(channels)-1]-channels[len(channels)-2])/2.0
        channel_edges[0] = channels[0]-chandiff1
        channel_edges[len(channels)] = channels[len(channels)-1]+chandiff2

        model = self._ax_model.pcolormesh(self._signal.data.phases,
                                          channel_edges,
                                          self.expected_counts,
                                          cmap = plt.get_cmap(self._model_cmap),
                                          vmin = self._vmin,
                                          vmax = self._vmax,
                                          linewidth = 0,
                                          rasterized = self._rasterized)
        model.set_edgecolor('face')

        model = self._ax_model.pcolormesh(self._signal.data.phases + 1.0,
                                          channel_edges,
                                          self.expected_counts,
                                          cmap = plt.get_cmap(self._model_cmap),
                                          vmin = self._vmin,
                                          vmax = self._vmax,
                                          linewidth = 0,
                                          rasterized = self._rasterized)
        model.set_edgecolor('face')

        self._ax_model.set_ylim([_np.max([channel_edges[0],0.001]),
                                 channel_edges[-1]])
        self._ax_model.set_yscale(self.yscale)

        self._model_cb = self._fig.colorbar(model, cax=self._ax_model_cb,
                                      ticks=_get_default_locator(None),
                                      format=_get_default_formatter())
        self._model_cb.ax.set_frame_on(True)
        self._model_cb.ax.yaxis.set_minor_locator(AutoMinorLocator())

        self._model_cb.set_label(label=r'counts/cycle', labelpad=15)

    def _add_residuals(self):
        """ Display the residuals in the third panel. """

        self._residuals = self.expected_counts - self._signal.data.counts
        self._residuals /= _np.sqrt(self.expected_counts)

        vmax =  _np.max( _np.abs( self._residuals ) )

        #Calculate channel edges by averaging:
        channels = self._signal.data.channels
        channel_edges = _np.zeros((len(self._signal.data.channels)+1))
        channel_edges[1:len(channels)] = (channels[:len(channels)-1]+channels[1:])/2.0
        chandiff1 = (channels[1]-channels[0])/2.0
        chandiff2 = (channels[len(channels)-1]-channels[len(channels)-2])/2.0
        channel_edges[0] = channels[0]-chandiff1
        channel_edges[len(channels)] = channels[len(channels)-1]+chandiff2

        # Plot
        if not self._blur_residuals:
            resid1 = self._ax_resid.pcolormesh(self._signal.data.phases,
                                        channel_edges,
                                        self._residuals,
                                        cmap = plt.get_cmap(self._residual_cmap),
                                        vmin = -vmax,
                                        vmax = vmax,
                                        linewidth = 0,
                                        rasterized = self._rasterized)

            resid2 = self._ax_resid.pcolormesh(self._signal.data.phases + 1.0,
                                        channel_edges,
                                        _np.abs(self._residuals),
                                        cmap = plt.get_cmap(self._residual_cmap),
                                        vmin = -vmax,
                                        vmax = vmax,
                                        linewidth = 0,
                                        rasterized = self._rasterized)

        else:
            k=30
            extended_phases = _np.linspace(self._signal.data.phases[0],self._signal.data.phases[-1:],
                                           num=(self._signal.data.phases.shape[0]-1)*k+1)
            extended_channels = _np.log10( _np.logspace(channel_edges[0],channel_edges[-1:],
                                                        num=(channel_edges.shape[0]-1)*k+1) )

            extended_residuals = _np.repeat( _np.repeat( self._residuals , k , axis=0 ), k , axis=1 )
            blurred_residuals = gaussian_filter( extended_residuals , 
                                                 sigma=k,
                                                 mode=['constant','wrap'] )

            resid1 = self._ax_resid.pcolormesh(_np.squeeze(extended_phases),
                                        _np.squeeze(extended_channels),
                                        blurred_residuals,
                                        cmap = plt.get_cmap(self._residual_cmap),
                                        vmin = -vmax,
                                        vmax = vmax,
                                        linewidth = 0,
                                        rasterized = self._rasterized)

            resid2 = self._ax_resid.pcolormesh(_np.squeeze(extended_phases) + 1.0,
                                        _np.squeeze(extended_channels),
                                        _np.abs(blurred_residuals),
                                        cmap = plt.get_cmap(self._residual_cmap),
                                        vmin = -vmax,
                                        vmax = vmax,
                                        linewidth = 0,
                                        rasterized = self._rasterized)
        
        resid1.set_edgecolor('face')
        resid2.set_edgecolor('face')
        self._ax_resid.axvline(1.0, lw=self._tick_width, color='k')

        self._ax_resid.set_ylim([_np.max([channel_edges[0],0.001]),
                                 channel_edges[-1]])
        self._ax_resid.set_yscale(self.yscale)

        self._resid_cb = self._fig.colorbar(resid2, cax = self._ax_resid_cb,
                                      ticks=AutoLocator())
        self._resid_cb.ax.set_frame_on(True)
        self._resid_cb.ax.yaxis.set_minor_locator(AutoMinorLocator())

        self._resid_cb.set_label(label=r'$(c_{ik}-d_{ik})/\sqrt{c_{ik}}$',
                                 labelpad=15)

    def _reveal_clusters(self):
        """ Display clusters from residuals in the fourth panel. """
        self._residuals = self.expected_counts - self._signal.data.counts
        self._residuals /= _np.sqrt(self.expected_counts)
        self._clusteresid = _np.abs( self._residuals ) >= self._threshlim
        self._lw, self._num = measurements.label(self._clusteresid)
        self._clustarea = measurements.sum(self._clusteresid, self._lw, index=arange(self._lw.max() + 1))
        self._affectedarea = self._clustarea[self._lw]
        vmaxresid =  _np.max( _np.abs( self._residuals ) )
        vmaxarea =  _np.max( self._affectedarea )

        #Calculate channel edges by averaging:
        channels = self._signal.data.channels
        channel_edges = _np.zeros((len(self._signal.data.channels)+1))
        channel_edges[1:len(channels)] = (channels[:len(channels)-1]+channels[1:])/2.0
        chandiff1 = (channels[1]-channels[0])/2.0
        chandiff2 = (channels[len(channels)-1]-channels[len(channels)-2])/2.0
        channel_edges[0] = channels[0]-chandiff1
        channel_edges[len(channels)] = channels[len(channels)-1]+chandiff2
        
        clust1 = self._ax_clres.pcolormesh(self._signal.data.phases,
                                      channel_edges,
                                      _np.where(self._clusteresid, self._residuals, 0),
                                      cmap = plt.get_cmap(self._residual_cmap),
                                      vmin = -vmaxresid,
                                      vmax = vmaxresid,
                                      linewidth = 0,
                                      rasterized = self._rasterized)
        clust1.set_edgecolor('face')

        clust2 = self._ax_clust.pcolormesh(self._signal.data.phases,
                                      channel_edges,
                                      self._affectedarea,
                                      cmap = plt.get_cmap(self._clusters_cmap),
                                      vmin = 0,
                                      vmax = vmaxarea,
                                      linewidth = 0,
                                      rasterized = self._rasterized)
        clust2.set_edgecolor('face')

        self._ax_clres.set_ylim([_np.max([channel_edges[0],0.001]),
                                 channel_edges[-1]])
        self._ax_clres.set_yscale(self.yscale)
        self._ax_clust.set_ylim([_np.max([channel_edges[0],0.001]),
                                 channel_edges[-1]])

        self._clust_cb = self._fig.colorbar(clust2, cax = self._ax_clust_cb,
                                      ticks=AutoLocator())
        self._clust_cb.ax.set_frame_on(True)
        self._clust_cb.ax.yaxis.set_minor_locator(AutoMinorLocator())

        self._clust_cb.set_label(label=r'cluster sizes for |residuals| > {}'.format(self._threshlim),
                                 labelpad=15)
        
    def _disp_distributions(self):
        """ Display residual and cluster distributions in the fifth panel. """
        self._residuals = self.expected_counts - self._signal.data.counts
        self._residuals /= _np.sqrt(self.expected_counts)
        self._clusteresid = _np.abs( self._residuals ) >= self._threshlim
        self._lw, self._num = measurements.label(self._clusteresid)
        self._clustarea = measurements.sum(self._clusteresid, self._lw, index=arange(self._lw.max() + 1))
        self._affectedarea = self._clustarea[self._lw]
        vmaxresid =  _np.max( _np.abs( self._residuals ) )
        vmaxarea =  _np.max( self._affectedarea )
        
        if _np.abs(_np.amin(self._residuals))< _np.abs(_np.amax(self._residuals)):
            minabsresid=(-1.0)*_np.amax(self._residuals)
            maxabsresid=_np.amax(self._residuals)
        else:
            minabsresid=_np.amin(self._residuals)
            maxabsresid=(-1.0)*_np.amin(self._residuals)
        
        residhist, binhist =  _np.histogram(self._residuals, bins=self._nbins, range=[minabsresid, maxabsresid])
        centphase = (binhist[:-1] + binhist[1:]) / 2
        binsize = (maxabsresid-minabsresid)/(self._nbins)
        scale=binsize*8640.0
        f = 1/(self._sigma * _np.sqrt(2 * _np.pi)) * _np.exp( - (centphase - self._mu)**2 / (2 * self._sigma**2) )

        totar = self._clustarea.flatten()
        totar = totar.astype(int)
        count_arr = _np.bincount(totar)
        
        rdist1 = self._ax_rdist.step(centphase, residhist)
        rdist2 = self._ax_rdist.plot(centphase, f*scale, linewidth=2, color='m')
        
        cdist = self._ax_cdist.step(np.linspace(0, vmaxarea.astype(int), len(count_arr)),
                                      count_arr,
                                      where='mid',
                                      )

    # -----------------------------------------------------------------------
    # 2-D residual structure analysis
    # -----------------------------------------------------------------------
    #
    #              1. **2D Power Spectrum** - detects periodic/wave-like patterns
    #                 via the FFT peak-to-median power ratio.
    #              2. **Moran's I** - global spatial autocorrelation statistic
    #                 (requires ``esda`` and ``libpysal``; skipped if absent).
    #              3. **Variogram** - measures how correlation decays with
    #                 pixel-distance on the phase-channel grid.
    #              4. **Runs test** - Wald-Wolfowitz test applied row- and
    #                 column-wise to detect 1-D autocorrelation along each axis.
    #              5. **Directional autocorrelation** - lag-1 correlation along
    #                 the horizontal, vertical, diagonal and anti-diagonal
    #                 directions; sensitive to obliquely clustered residuals
    #                 that the axis-aligned runs test can miss.
    #              6. **Ljung-Box + KS** - serial autocorrelation in the
    #                 flattened residuals and a check for Gaussianity
    #                 (requires ``statsmodels``; KS test always runs).
    # -----------------------------------------------------------------------

    def _run_structure_analysis(self):
        """Run all 2-D structure tests on ``self._residuals`` and print a
        summary.  Results are stored in :attr:`residual_structure_results`
        and :attr:`residual_structure_verdict`.

        Called automatically from :meth:`finalize` when
        ``analyse_residual_structure=True``.
        """
        if not hasattr(self, '_residuals'):
            # Recompute the standardised residuals in case finalize has not
            # stored them yet (should not normally happen, but be defensive).
            self._residuals = (self.expected_counts -
                               self._signal.data.counts) / _np.sqrt(
                self.expected_counts)

        # NOTE: the order of this list defines the reporting order AND the
        # panel order in plot_structure_diagnostics().
        results = [
            self._structure_test_power_spectrum(
                self._residuals,
                n_null=self._structure_power_n_null,
                alpha=self._structure_power_alpha,
                low_k_fraction=self._structure_power_low_k_fraction,
                random_seed=self._structure_power_random_seed),
            self._structure_test_morans_i(self._residuals),
            self._structure_test_variogram(
                self._residuals,
                max_lag=self._structure_variogram_max_lag,
                n_null=self._structure_variogram_n_null,
                alpha=self._structure_variogram_alpha,
                random_seed=self._structure_variogram_random_seed),
            self._structure_test_runs(self._residuals),
            self._structure_test_acf2d(self._residuals, self._structure_acf2d_max_lag),
            self._structure_test_normality_ljungbox(self._residuals),
        ]

        # Aggregate verdict: >=3 of the definitive tests flag structure
        # NOTE: `from pylab import *` shadows the builtin `sum` with
        # numpy's, which rejects generator expressions -- count with a
        # list comprehension instead, immune to the shadowing.
        definitive = [t for t in results if t['structured'] is not None]
        n_structured = len([t for t in definitive if t['structured']])
        n_total = len(definitive)

        if n_structured >= 3:
            verdict = 'STRUCTURED'
        elif n_structured == 0:
            verdict = 'RANDOM'
        else:
            verdict = 'INCONCLUSIVE'

        self.residual_structure_results = results
        self.residual_structure_verdict = verdict

        # -- Pretty-print summary ------------------------------------------
        print('\n' + '=' * 65)
        print('  XPSI 2-D RESIDUAL STRUCTURE ANALYSIS')
        print('  Axes: phase [cycles] x channel')
        print('=' * 65)
        for t in results:
            if t['structured'] is True:
                flag = 'x STRUCTURED'
            elif t['structured'] is False:
                flag = 'v random    '
            else:
                flag = '- skipped   '
            print(f"  [{flag}]  {t['test']}")
            print(f"                {t['interpretation']}")
        print('-' * 65)
        print(f"  VERDICT : {verdict}  "
              f"({n_structured}/{n_total} tests flagged structure)")
        print('=' * 65 + '\n')

    # -- Individual tests ---------------------------------------------------

    @staticmethod
    def _structure_test_power_spectrum(residuals, n_null=500, alpha=0.05,
                                       low_k_fraction=0.15, random_seed=0):
        """2-D power spectrum test, Monte-Carlo calibrated against white noise.

        Two complementary statistics are computed from the 2-D periodogram
        (DC component excluded):

        * ``g_max`` - Fisher's g statistic: maximum power / total power.
          Sensitive to PERIODIC structure, whose power concentrates in a
          few Fourier modes.
        * ``f_lowk`` - fraction of total power at low spatial frequencies
          (radial frequency below ``low_k_fraction`` of the Nyquist
          radius).  Sensitive to SMOOTH, blobby structure, whose power
          concentrates at low k.

        The null distributions of both statistics depend on the grid
        shape, so they are calibrated empirically by simulating
        ``n_null`` white-noise fields of the same shape and standard
        deviation.  One-sided empirical p-values (with an add-one
        correction, so the smallest reportable p is ``1/(n_null+1)``)
        are combined with a Bonferroni factor of 2: structure is flagged
        when ``min(p) < alpha / 2``.

        Unlike a fixed peak/median-ratio threshold, ``alpha`` is a
        genuine significance level (default 0.05), consistent with the
        other tests in the suite, and requires no per-dataset tuning.

        :param residuals: 2-D array of standardised residuals.
        :param int n_null: Number of white-noise simulations for the
            null calibration.  Default 500.
        :param float alpha: Significance level before the Bonferroni
            factor of 2.  Default 0.05.
        :param float low_k_fraction: Radial-frequency cut (as a fraction
            of the Nyquist radius) defining the "low-k" region for the
            smooth-structure statistic.  Default 0.15.
        :param int random_seed: Seed for the null simulations, for
            reproducibility.
        :returns: dict with keys ``test``, ``g_max``, ``g_null_q95``,
            ``p_periodic``, ``f_lowk``, ``f_null_q95``, ``p_smooth``,
            ``min_p``, ``alpha_bonferroni``, ``n_null``, ``power_map``,
            ``structured`` (bool), and ``interpretation`` (str).
        """
        rng = _np.random.default_rng(random_seed)
        ny, nx = residuals.shape

        def _spectrum_stats(field):
            p = _np.abs(_np.fft.fftshift(_np.fft.fft2(field))) ** 2
            p[ny // 2, nx // 2] = 0.0          # remove DC
            total = p.sum()
            ky = _np.fft.fftshift(_np.fft.fftfreq(ny))[:, None]
            kx = _np.fft.fftshift(_np.fft.fftfreq(nx))[None, :]
            kr = _np.sqrt(ky ** 2 + kx ** 2) / _np.sqrt(0.5)  # 0..1
            low = (kr <= low_k_fraction) & (kr > 0)
            return p.max() / total, p[low].sum() / total, p

        g_obs, f_obs, power_map = _spectrum_stats(residuals)

        g_null = _np.empty(n_null)
        f_null = _np.empty(n_null)
        sd = residuals.std()
        for i in range(n_null):
            g_null[i], f_null[i], _ = _spectrum_stats(
                rng.normal(0.0, sd, (ny, nx)))

        # one-sided empirical p-values (add-one correction avoids p = 0)
        p_g = (_np.count_nonzero(g_null >= g_obs) + 1.0) / (n_null + 1.0)
        p_f = (_np.count_nonzero(f_null >= f_obs) + 1.0) / (n_null + 1.0)

        min_p      = min(p_g, p_f)
        alpha_bonf = alpha / 2.0
        structured = bool(min_p < alpha_bonf)

        return {
            'test'             : '2D Power Spectrum (MC-calibrated)',
            'g_max'            : float(g_obs),
            'g_null_q95'       : float(_np.percentile(g_null, 95)),
            'p_periodic'       : float(p_g),
            'f_lowk'           : float(f_obs),
            'f_null_q95'       : float(_np.percentile(f_null, 95)),
            'p_smooth'         : float(p_f),
            'min_p'            : float(min_p),
            'alpha_bonferroni' : alpha_bonf,
            'n_null'           : n_null,
            'power_map'        : power_map,
            'structured'       : structured,
            'interpretation'   : (
                f'periodic: g = {g_obs:.4f} (p = {p_g:.4f});  '
                f'smooth: f_lowk = {f_obs:.3f} (p = {p_f:.4f});  '
                f'alpha/2 = {alpha_bonf:.3f}  '
                f'-> {"structured" if structured else "random"}'
            ),
        }

    @staticmethod
    def _structure_test_morans_i(residuals):
        """Global Moran's I spatial autocorrelation test.

        Moran's I quantifies how similar neighbouring cells are to one
        another.  Under the null hypothesis of complete spatial
        randomness, I ~ 0.  A significantly positive I indicates positive
        autocorrelation (nearby residuals tend to have the same sign).

        Queen contiguity weights (8-connected neighbourhood) are used so
        that row/column AND diagonal neighbours are all considered.

        :param residuals: 2-D array of standardised residuals.
        :returns: dict with keys ``test``, ``I``, ``p_value``,
            ``structured`` (bool or None if skipped), and
            ``interpretation`` (str).

        .. note::
            Requires the optional packages ``esda`` and ``libpysal``
            (``pip install esda libpysal``).  Skipped if unavailable.
        """
        if not _HAS_ESDA:
            return {
                'test': "Moran's I (spatial autocorrelation)",
                'structured': None,
                'interpretation': (
                    'skipped - install esda + libpysal to enable '
                    '(pip install esda libpysal)'
                ),
            }

        ny, nx = residuals.shape
        w = _lat2W(ny, nx, rook=False)  # queen contiguity (8-connected)
        w.transform = 'r'  # row-standardise
        mi = _Moran(residuals.ravel(), w)
        structured = bool(mi.p_sim < 0.05)

        return {
            'test': "Moran's I (spatial autocorrelation)",
            'I': float(mi.I),
            'p_value': float(mi.p_sim),
            'structured': structured,
            'interpretation': (
                f"I = {mi.I:.4f},  p = {mi.p_sim:.4f}  "
                f"-> {'structured' if structured else 'random'}"
            ),
        }

    @staticmethod
    def _structure_test_variogram(residuals, max_lag=20, n_null=500,
                                  alpha=0.05, random_seed=0):
        """Empirical isotropic variogram test, Monte-Carlo calibrated.

        For a random residual field the normalised variogram sits at the
        sill (``gamma(h)/sigma^2 ~ 1``) for every lag - a *pure nugget*
        model.  Spatial structure pulls it away from the sill: smooth
        short-range correlation depresses it at small lags, while
        periodic or banded structure produces a "hole effect" (the
        variogram rises above and/or falls below the sill at larger
        lags).

        Two complementary statistics are therefore computed:

        * ``nugget_deficit`` = ``1 - gamma(1)/sigma^2``.  One-sided:
          a large positive value means strong lag-1 correlation, the
          most common signature of misfit structure.
        * ``max_abs_dev`` = ``max_h |1 - gamma(h)/sigma^2|`` over all
          lags.  Catches hole-effect and longer-range departures that a
          lag-1 (or slope) statistic dilutes or misses entirely.

        Both null distributions depend on the grid shape, so they are
        calibrated by simulating ``n_null`` white-noise fields of the
        same shape.  Empirical p-values (add-one corrected, so the
        smallest reportable p is ``1/(n_null+1)``) are combined with a
        Bonferroni factor of 2: structure is flagged when
        ``min(p) < alpha / 2``.  ``alpha`` is a genuine significance
        level, consistent with the 2D power spectrum test, and replaces
        the earlier fixed slope threshold.

        :param residuals: 2-D array of standardised residuals.
        :param int max_lag: Maximum pixel lag.  Defaults to 20.
        :param int n_null: Number of white-noise simulations for the
            null calibration.  Default 500.
        :param float alpha: Significance level before the Bonferroni
            factor of 2.  Default 0.05.
        :param int random_seed: Seed for the null simulations.
        :returns: dict with keys ``test``, ``lags``, ``semivariance``
            (raw), ``semivariance_norm``, ``sample_variance``,
            ``nugget_deficit``, ``p_nugget``, ``max_abs_dev``,
            ``p_maxdev``, ``min_p``, ``alpha_bonferroni``, ``n_null``,
            ``structured`` (bool), and ``interpretation`` (str).
        """

        def _variogram_curve(field):
            lags = _np.arange(1, max_lag + 1)
            sv = _np.empty(len(lags))
            for k, lag in enumerate(lags):
                h_diffs = (field[:, lag:] - field[:, :-lag]).ravel()
                v_diffs = (field[lag:, :] - field[:-lag, :]).ravel()
                sv[k] = 0.5 * _np.mean(
                    _np.concatenate([h_diffs, v_diffs]) ** 2)
            return lags, sv

        def _stats_pair(field):
            _, sv = _variogram_curve(field)
            svn = sv / (_np.var(field) + 1e-12)
            return (1.0 - svn[0]), _np.max(_np.abs(1.0 - svn))

        lags, semivariance = _variogram_curve(residuals)
        sample_variance    = float(_np.var(residuals))
        semivariance_norm  = semivariance / (sample_variance + 1e-12)

        d_obs, m_obs = _stats_pair(residuals)

        rng    = _np.random.default_rng(random_seed)
        d_null = _np.empty(n_null)
        m_null = _np.empty(n_null)
        sd = residuals.std()
        for i in range(n_null):
            d_null[i], m_null[i] = _stats_pair(
                rng.normal(0.0, sd, residuals.shape))

        # one-sided empirical p-values (add-one correction avoids p = 0)
        p_d = (_np.count_nonzero(d_null >= d_obs) + 1.0) / (n_null + 1.0)
        p_m = (_np.count_nonzero(m_null >= m_obs) + 1.0) / (n_null + 1.0)

        min_p      = min(p_d, p_m)
        alpha_bonf = alpha / 2.0
        structured = bool(min_p < alpha_bonf)

        return {
            'test'              : 'Variogram (MC-calibrated)',
            'lags'              : lags,
            'semivariance'      : semivariance,
            'semivariance_norm' : semivariance_norm,
            'sample_variance'   : sample_variance,
            'nugget_deficit'    : float(d_obs),
            'p_nugget'          : float(p_d),
            'max_abs_dev'       : float(m_obs),
            'p_maxdev'          : float(p_m),
            'min_p'             : float(min_p),
            'alpha_bonferroni'  : alpha_bonf,
            'n_null'            : n_null,
            'structured'        : structured,
            'interpretation'    : (
                f'lag-1 nugget deficit = {d_obs:.3f} (p = {p_d:.4f});  '
                f'max |sill deviation| = {m_obs:.3f} (p = {p_m:.4f});  '
                f'alpha/2 = {alpha_bonf:.3f}  '
                f'-> {"structured" if structured else "random"}'
            ),
        }

    @staticmethod
    def _structure_test_runs(residuals):
        """Wald-Wolfowitz runs test applied per-row and per-column.

        A *run* is a maximal sequence of values all above or all below the
        row/column median.  Too few runs imply positive autocorrelation
        along that axis.  We flag structure when more than 15 % of rows
        *or* columns individually reject randomness at the 5 % level
        (~5 % would be expected by chance).

        .. note::
            This test scans the horizontal and vertical directions only.
            Structure aligned with the grid diagonals is targeted by
            :meth:`_structure_test_directional_autocorr`.

        :param residuals: 2-D array of standardised residuals.
        :returns: dict with keys ``test``, ``frac_rows_significant``,
            ``frac_cols_significant``, ``threshold``, ``structured``
            (bool), and ``interpretation`` (str).
        """

        def _runs_pvalue(x):
            """One-sided p-value for too-few runs (positive autocorrelation)."""
            median = _np.median(x)
            signs = _np.sign(x - median)
            signs = signs[signs != 0]
            n = len(signs)
            if n < 10:
                return 1.0
            n1 = int(_np.sum(signs > 0))
            n2 = int(_np.sum(signs < 0))
            if n1 == 0 or n2 == 0:
                return 0.0
            runs = 1 + int(_np.sum(signs[1:] != signs[:-1]))
            mu_r = 2.0 * n1 * n2 / n + 1.0
            var_r = (2.0 * n1 * n2 * (2.0 * n1 * n2 - n)
                     / (n ** 2 * (n - 1)))
            if var_r <= 0:
                return 1.0
            z = (runs - mu_r) / _np.sqrt(var_r)
            return float(_stats.norm.cdf(z))  # left tail: few runs -> small p

        n_rows, n_cols = residuals.shape
        row_pvals = [_runs_pvalue(residuals[i, :]) for i in range(n_rows)]
        col_pvals = [_runs_pvalue(residuals[:, j]) for j in range(n_cols)]

        frac_rows = float(_np.mean(_np.array(row_pvals) < 0.05))
        frac_cols = float(_np.mean(_np.array(col_pvals) < 0.05))

        threshold = 0.15
        structured = bool((frac_rows > threshold) or (frac_cols > threshold))

        return {
            'test': 'Runs test (per row & column)',
            'frac_rows_significant': frac_rows,
            'frac_cols_significant': frac_cols,
            'threshold': threshold,
            'structured': structured,
            'interpretation': (
                f'rows with structure: {frac_rows:.1%},  '
                f'cols with structure: {frac_cols:.1%}  '
                f'(threshold {threshold:.0%})  '
                f'-> {"structured" if structured else "random"}'
            ),
        }

    @staticmethod
    def _structure_test_acf2d(residuals, max_lag=8):
        """Two-dimensional autocorrelation-function (ACF) scan.

        Generalises a directional lag-1 test to EVERY lag vector
        ``(dy, dx)`` with ``0 <= dy <= max_lag`` and ``|dx| <= max_lag``
        (upper half-plane only, since the ACF is symmetric under
        ``(dy, dx) -> (-dy, -dx)``).  This makes the test sensitive to
        clustered residuals at ANY orientation - horizontal, vertical,
        diagonal, or any oblique angle - and at any separation up to
        ``max_lag`` pixels, so it also catches structure whose
        correlation appears only at larger lags (e.g. regularly spaced
        blobs).

        For each lag vector, the Pearson correlation ``r`` between the
        residual field and its shifted copy is computed.  Under the
        white-noise null hypothesis ``r * sqrt(N)`` is asymptotically
        standard normal, giving a two-sided p-value per lag.  A
        Bonferroni correction over the total number of lags tested is
        applied: structure is flagged when ``min(p) < 0.05 / n_lags``.

        :param residuals: 2-D array of standardised residuals.
        :param int max_lag: Maximum lag in pixels along each axis.
            Default 8, giving ``max_lag * (2 * max_lag + 1) + max_lag``
            = 144 tested lag vectors.
        :returns: dict with keys ``test``, ``acf_map`` (2-D array of
            ``r`` values, NaN where untested), ``max_lag``, ``best_lag``
            (tuple ``(dy, dx)`` of the most significant lag),
            ``best_r``, ``best_angle_deg`` (orientation of the
            strongest correlation, degrees from the +phase axis),
            ``min_p``, ``alpha_bonferroni``, ``n_lags_tested``,
            ``structured`` (bool), and ``interpretation`` (str).
        """
        m, n = residuals.shape
        acf_map = _np.full((max_lag + 1, 2 * max_lag + 1), _np.nan)
        p_values = []
        best = {'z': 0.0, 'dy': 0, 'dx': 1, 'r': 0.0}

        for dy in range(0, max_lag + 1):
            for dx in range(-max_lag, max_lag + 1):
                if dy == 0 and dx <= 0:
                    continue  # skip (0,0) and mirror duplicates
                if dx >= 0:
                    a = residuals[dy:, dx:]
                    b = residuals[:m - dy, :n - dx]
                else:
                    a = residuals[dy:, :n + dx]
                    b = residuals[:m - dy, -dx:]
                a = a.ravel()
                b = b.ravel()
                npairs = a.size
                if npairs < 10:
                    continue
                r = float(_np.corrcoef(a, b)[0, 1])
                z = r * _np.sqrt(npairs)
                p = float(2.0 * _stats.norm.sf(abs(z)))
                acf_map[dy, dx + max_lag] = r
                p_values.append(p)
                if abs(z) > abs(best['z']):
                    best = {'z': z, 'dy': dy, 'dx': dx, 'r': r}

        n_lags = len(p_values)
        alpha_bonf = 0.05 / n_lags
        min_p = min(p_values)
        structured = bool(min_p < alpha_bonf)
        angle = float(_np.degrees(_np.arctan2(best['dy'], best['dx'])))

        return {
            'test': '2D ACF scan (all orientations & lags)',
            'acf_map': acf_map,
            'max_lag': max_lag,
            'best_lag': (best['dy'], best['dx']),
            'best_r': float(best['r']),
            'best_angle_deg': angle,
            'min_p': float(min_p),
            'alpha_bonferroni': alpha_bonf,
            'n_lags_tested': n_lags,
            'structured': structured,
            'interpretation': (
                f"strongest r = {best['r']:+.3f} at lag "
                f"(dchan={best['dy']}, dphase={best['dx']}), "
                f"orientation ~ {angle:.0f} deg;  "
                f"min-p = {min_p:.2e} over {n_lags} lags "
                f"(Bonferroni alpha = {alpha_bonf:.1e})  "
                f"-> {'structured' if structured else 'random'}"
            ),
        }

    @staticmethod
    def _structure_test_normality_ljungbox(residuals):
        """Kolmogorov-Smirnov normality test + Ljung-Box autocorrelation test.

        * **KS test** - compares the marginal distribution of the
          standardised residuals against a fitted normal.
        * **Ljung-Box test** - detects autocorrelation at lags 1-20 in the
          residuals flattened in row-major (phase-fast) order.

        :param residuals: 2-D array of standardised residuals.
        :returns: dict with keys ``test``, ``ks_stat``, ``ks_pvalue``,
            ``ljungbox_min_pvalue`` (``None`` if statsmodels absent),
            ``structured`` (bool), and ``interpretation`` (str).

        .. note::
            The Ljung-Box test requires ``statsmodels``.  If absent, only
            the KS test contributes to the verdict.
        """
        flat = residuals.ravel()
        mu_hat = float(flat.mean())
        sd_hat = float(flat.std())

        ks_stat, ks_p = _stats.kstest(flat, 'norm', args=(mu_hat, sd_hat))
        ks_stat = float(ks_stat)
        ks_p = float(ks_p)

        if _HAS_STATSMODELS:
            lb_result = _acorr_ljungbox(flat, lags=20, return_df=True)
            lb_min_p = float(lb_result['lb_pvalue'].min())
            structured = bool((ks_p < 0.01) or (lb_min_p < 0.05))
            lb_str = f'Ljung-Box min-p = {lb_min_p:.4f}'
        else:
            lb_min_p = None
            structured = bool(ks_p < 0.01)
            lb_str = ('Ljung-Box skipped '
                      '(pip install statsmodels to enable)')

        return {
            'test': 'KS normality + Ljung-Box autocorrelation',
            'ks_stat': ks_stat,
            'ks_pvalue': ks_p,
            'ljungbox_min_pvalue': lb_min_p,
            'structured': structured,
            'interpretation': (
                f'KS p = {ks_p:.4f},  {lb_str}  '
                f'-> {"structured" if structured else "random"}'
            ),
        }

    # -- Diagnostic figure --------------------------------------------------

    def _get_channel_edges(self):
        """Compute channel bin edges by averaging adjacent channel centres.

        Same construction as used by the main channel-phase panels.
        """
        channels = self._signal.data.channels
        channel_edges = _np.zeros(len(channels) + 1)
        channel_edges[1:len(channels)] = (channels[:-1] + channels[1:]) / 2.0
        channel_edges[0] = channels[0] - (channels[1] - channels[0]) / 2.0
        channel_edges[-1] = channels[-1] + (channels[-1] - channels[-2]) / 2.0
        return channel_edges

    def plot_structure_diagnostics(self, figsize=(20, 9), savefile=None):
        """Produce a seven-panel diagnostic figure for the structure tests.

        Must be called *after* :meth:`finalize` with
        ``analyse_residual_structure=True``.

        The panels follow the SAME ORDER as the tests reported by
        :meth:`finalize`:

        * **(0) Input** - standardised residual map over the actual
          phase [cycles] x channel grid (log channel axis), identical
          axes to the main residual panel.
        * **(1) Test 1** - 2-D power spectrum (log scale).
        * **(2) Test 3** - normalised empirical variogram with a dashed
          pure-nugget reference line.
        * **(3) Test 5** - 2D autocorrelation-function map over all lag
          vectors, with the most significant lag marked.
        * **(4) Test 6** - normal Q-Q plot of the flattened residuals.
        * **(5) Summary** - text panel listing all six test verdicts in
          report order (covers Moran's I and the runs test, which have
          no natural plot).

        :param tuple figsize: Figure size in inches.  Default ``(16, 9)``.
        :param str savefile: If given, the figure is saved to this path.
        :raises AttributeError: If called before :meth:`finalize` or with
            ``analyse_residual_structure=False``.
        """
        if not hasattr(self, 'residual_structure_results'):
            raise AttributeError(
                'Call finalize() with analyse_residual_structure=True first.')

        results = self.residual_structure_results
        results_map = {t['test']: t for t in results}
        verdict = self.residual_structure_verdict
        colour_map = {
            'STRUCTURED': '#d62728',
            'RANDOM': '#2ca02c',
            'INCONCLUSIVE': '#ff7f0e',
        }

        from matplotlib import gridspec as _gspec
        fig = plt.figure(figsize=figsize)
        fig.suptitle(
            f'XPSI Residual Structure Analysis - Verdict: {verdict}',
            fontsize=16, fontweight='bold',
            color=colour_map.get(verdict, 'black'))
        gs = _gspec.GridSpec(2, 3, figure=fig, hspace=0.42, wspace=0.32)

        residuals = self._residuals

        # -- Panel (0,0): residual map on the ACTUAL phase-channel grid ----
        ax1 = fig.add_subplot(gs[0, 0])
        vmax = _np.percentile(_np.abs(residuals), 98)
        channel_edges = self._get_channel_edges()
        pcm = ax1.pcolormesh(self._signal.data.phases,
                             channel_edges,
                             residuals,
                             cmap=plt.get_cmap(self._residual_cmap),
                             vmin=-vmax, vmax=vmax,
                             linewidth=0,
                             rasterized=self._rasterized)
        pcm.set_edgecolor('face')
        fig.colorbar(pcm, ax=ax1, fraction=0.046, pad=0.04)
        ax1.set_ylim([_np.max([channel_edges[0], 0.001]), channel_edges[-1]])
        ax1.set_yscale(self.yscale)
        ax1.set_title('Standardised Residuals')
        ax1.set_xlabel(r'$\phi$ [cycles]')
        ax1.set_ylabel('channel')
        ax1.xaxis.set_major_locator(MultipleLocator(0.2))
        ax1.xaxis.set_minor_locator(MultipleLocator(0.05))

        # -- Panel (0,1): Test 1 - 2-D power spectrum ----------------------
        ax2 = fig.add_subplot(gs[0, 1])
        ps_t = results_map.get('2D Power Spectrum (MC-calibrated)', {})
        if 'power_map' in ps_t:
            ax2.imshow(_np.log1p(ps_t['power_map']),
                       aspect='auto', origin='lower', cmap='inferno')
            ax2.set_title(
                f"Test 1: 2D Power Spectrum (log)\n"
                f"p_periodic = {ps_t.get('p_periodic', float('nan')):.4f}, "
                f"p_smooth = {ps_t.get('p_smooth', float('nan')):.4f}")
        else:
            ax2.text(0.5, 0.5, 'No data', ha='center', va='center',
                     transform=ax2.transAxes)
            ax2.set_title('Test 1: 2D Power Spectrum')
        ax2.set_xlabel(r'$k_\phi$')
        ax2.set_ylabel(r'$k_\mathrm{ch}$')

        # -- Panel (0,2): Test 3 - variogram (normalised) ------------------
        ax3 = fig.add_subplot(gs[0, 2])
        vg_t = results_map.get('Variogram (MC-calibrated)', {})
        if 'lags' in vg_t:
            ax3.plot(vg_t['lags'], vg_t['semivariance_norm'],
                     'o-', color='#1f77b4', ms=4, label='empirical')
            ax3.axhline(1.0, ls='--', color='gray', lw=1,
                        label='pure nugget (random)')
            ax3.set_title(
                'Test 3: Variogram (normalised)\n'
                f"p_nugget = {vg_t.get('p_nugget', float('nan')):.4f}, "
                f"p_maxdev = {vg_t.get('p_maxdev', float('nan')):.4f}")
        else:
            ax3.text(0.5, 0.5, 'No data', ha='center', va='center',
                     transform=ax3.transAxes)
            ax3.set_title('Test 3: Variogram (normalised)')
        ax3.set_xlabel('Lag (pixels)')
        ax3.set_ylabel('Normalised semivariance')
        ax3.legend(fontsize=9, loc='best')

        # -- Panel (1,0): Test 5 - 2D ACF scan ------------------------------
        ax4 = fig.add_subplot(gs[1, 0])
        acf_t = results_map.get('2D ACF scan (all orientations & lags)', {})
        if 'acf_map' in acf_t:
            L = acf_t['max_lag']
            amap = acf_t['acf_map']
            amax = _np.nanmax(_np.abs(amap))
            im = ax4.imshow(amap, aspect='auto', origin='lower',
                            cmap='RdBu_r', vmin=-amax, vmax=amax,
                            extent=[-L - 0.5, L + 0.5, -0.5, L + 0.5])
            fig.colorbar(im, ax=ax4, fraction=0.046, pad=0.04,
                         label=r'Pearson $r$')
            bdy, bdx = acf_t['best_lag']
            ax4.plot(bdx, bdy, 'x', color='k', ms=10, mew=2,
                     label=f'strongest lag ({bdy},{bdx})')
            ax4.set_title('Test 5: 2D ACF scan\n'
                          f"min-p = {acf_t.get('min_p', float('nan')):.2e}")
            ax4.set_xlabel(r'$\Delta$ phase bin')
            ax4.set_ylabel(r'$\Delta$ channel bin')
            ax4.legend(fontsize=9, loc='upper right')
        else:
            ax4.text(0.5, 0.5, 'No data', ha='center', va='center',
                     transform=ax4.transAxes)
            ax4.set_title('Test 5: 2D ACF scan')

        # -- Panel (1,1): Test 6 - Q-Q plot --------------------------------
        ax5 = fig.add_subplot(gs[1, 1])
        flat = residuals.ravel()
        (osm, osr), (slope, intercept, _) = _stats.probplot(flat, dist='norm')
        ax5.plot(osm, osr, '.', ms=1, alpha=0.4, color='#1f77b4')
        ax5.plot(osm, slope * _np.array(osm) + intercept,
                 'r-', lw=1.5, label='normal reference')
        ks_t = results_map.get('KS normality + Ljung-Box autocorrelation', {})
        ks_p = ks_t.get('ks_pvalue', float('nan'))
        ax5.set_title(f'Test 6: Normal Q-Q  (KS p = {ks_p:.4f})')
        ax5.set_xlabel('Theoretical quantiles')
        ax5.set_ylabel('Sample quantiles')
        ax5.legend(fontsize=9)

        # -- Panel (1,2): summary of ALL tests in report order -------------
        import textwrap as _textwrap
        ax6 = fig.add_subplot(gs[1, 2])
        ax6.axis('off')
        ax6.set_title('Summary (report order)')
        y = 0.95
        for i, t in enumerate(results, start=1):
            if t['structured'] is True:
                mark, col = 'STRUCTURED', '#d62728'
            elif t['structured'] is False:
                mark, col = 'Random', '#2ca02c'
            else:
                mark, col = 'skipped', 'gray'
            # Short label: first part of the test name, wrapped over up
            # to two lines so long names never collide with the verdict
            label = t['test'].split(' (')[0]
            wrapped = _textwrap.wrap(f'{i}. {label}', width=30)[:2]
            for j, line in enumerate(wrapped):
                ax6.text(0.02, y - j * 0.055, line,
                         transform=ax6.transAxes, fontsize=13, va='top')
            ax6.text(0.98, y, mark, transform=ax6.transAxes,
                     fontsize=13, va='top', ha='right',
                     color=col, fontweight='bold')
            y -= 0.11 + 0.055 * (len(wrapped) - 1)
        ax6.text(0.02, y - 0.02,
                 f'VERDICT: {verdict}',
                 transform=ax6.transAxes, fontsize=13, va='top',
                 fontweight='bold', color=colour_map.get(verdict, 'black'))

        if savefile is not None:
            fig.savefig(savefile, dpi=150, bbox_inches='tight')
            print(f'Diagnostic figure saved -> {savefile}')

        return fig
