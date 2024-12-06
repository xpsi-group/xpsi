from pylab import *
from ._global_imports import *
from scipy.ndimage import measurements, gaussian_filter
from ._signalplot import SignalPlot

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
    __caching_targets__ = [] #  ['expected_counts']

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
                 blur_residuals=False,
                 plot_clusters=False,
                 threshlim=2.0,
                 clusters_cmap='PuOr',
                 clustdist_cmap='PuOr',
                 mu=0.0,
                 sigma=1.0,
                 nbins=50,
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
        self._blur_residuals = blur_residuals
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
                ax.set_xlabel('$\phi$ [cycles]')

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
        except AttributeError:
            self._model_sum = self._signal.expected_counts.copy()
        else:
            self._model_sum += self._signal.expected_counts

    @property
    def model_sum(self):
        """ Get the current posterior sum of the count numbers. """
        return self._model_sum

    @model_sum.deleter
    def model_sum(self):
        del self._model_sum

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

        # Compute the chi squared
        chi2 = self._signal.chi2

        # Plot pulse
        self._ax_pulse.errorbar(x=doubles_phases,
                                y=double_pulse_data,
                                yerr=_np.sqrt( double_pulse_data ),
                                ds='steps-mid',
                                label='Data', color='black')

        self._ax_pulse.errorbar(x=doubles_phases,
                                y=double_pulse_model, 
                                ds='steps-mid',
                                label=r'Model $\chi^2=$'+f'{chi2:.2f}', color='steelblue')

        self._ax_pulse.legend(loc='lower right')

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
                                        cmap = cm.get_cmap(self._data_cmap),
                                        vmin = self._vmin,
                                        vmax = self._vmax,
                                        linewidth = 0,
                                        rasterized = self._rasterized)
        data.set_edgecolor('face')

        data = self._ax_data.pcolormesh(self._signal.data.phases + 1.0,
                                        channel_edges,
                                        self._signal.data.counts,
                                        cmap = cm.get_cmap(self._data_cmap),
                                        vmin = self._vmin,
                                        vmax = self._vmax,
                                        linewidth = 0,
                                        rasterized = self._rasterized)
        data.set_edgecolor('face')

        self._ax_data.set_ylim([_np.max([channel_edges[0],0.001]),
                                channel_edges[-1]])
        self._ax_data.set_yscale(self.yscale)

        self._data_cb = plt.colorbar(data, cax=self._ax_data_cb,
                                     ticks=_get_default_locator(None),
                                     format=_get_default_formatter())
        self._data_cb.ax.set_frame_on(True)
        self._data_cb.ax.yaxis.set_minor_locator(AutoMinorLocator())
        self._data_cb.set_label(label=r'counts', labelpad=15)

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
                                          cmap = cm.get_cmap(self._model_cmap),
                                          vmin = self._vmin,
                                          vmax = self._vmax,
                                          linewidth = 0,
                                          rasterized = self._rasterized)
        model.set_edgecolor('face')

        model = self._ax_model.pcolormesh(self._signal.data.phases + 1.0,
                                          channel_edges,
                                          self.expected_counts,
                                          cmap = cm.get_cmap(self._model_cmap),
                                          vmin = self._vmin,
                                          vmax = self._vmax,
                                          linewidth = 0,
                                          rasterized = self._rasterized)
        model.set_edgecolor('face')

        self._ax_model.set_ylim([_np.max([channel_edges[0],0.001]),
                                 channel_edges[-1]])
        self._ax_model.set_yscale(self.yscale)

        self._model_cb = plt.colorbar(model, cax=self._ax_model_cb,
                                      ticks=_get_default_locator(None),
                                      format=_get_default_formatter())
        self._model_cb.ax.set_frame_on(True)
        self._model_cb.ax.yaxis.set_minor_locator(AutoMinorLocator())

        self._model_cb.set_label(label=r'counts', labelpad=15)

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
                                        cmap = cm.get_cmap(self._residual_cmap),
                                        vmin = -vmax,
                                        vmax = vmax,
                                        linewidth = 0,
                                        rasterized = self._rasterized)

            resid2 = self._ax_resid.pcolormesh(self._signal.data.phases + 1.0,
                                        channel_edges,
                                        _np.abs(self._residuals),
                                        cmap = cm.get_cmap(self._residual_cmap),
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
                                        cmap = cm.get_cmap(self._residual_cmap),
                                        vmin = -vmax,
                                        vmax = vmax,
                                        linewidth = 0,
                                        rasterized = self._rasterized)

            resid2 = self._ax_resid.pcolormesh(_np.squeeze(extended_phases) + 1.0,
                                        _np.squeeze(extended_channels),
                                        _np.abs(blurred_residuals),
                                        cmap = cm.get_cmap(self._residual_cmap),
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

        self._resid_cb = plt.colorbar(resid2, cax = self._ax_resid_cb,
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
                                      cmap = cm.get_cmap(self._residual_cmap),
                                      vmin = -vmaxresid,
                                      vmax = vmaxresid,
                                      linewidth = 0,
                                      rasterized = self._rasterized)
        clust1.set_edgecolor('face')

        clust2 = self._ax_clust.pcolormesh(self._signal.data.phases,
                                      channel_edges,
                                      self._affectedarea,
                                      cmap = cm.get_cmap(self._clusters_cmap),
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

        self._clust_cb = plt.colorbar(clust2, cax = self._ax_clust_cb,
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
