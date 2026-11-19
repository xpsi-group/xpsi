from pylab import *
from ._global_imports import *
from scipy.ndimage import gaussian_filter
from scipy.ndimage import label as _ndi_label
from scipy.ndimage import sum_labels as _ndi_sum_labels
from scipy import stats as _stats
from ._signalplot import SignalPlot
from matplotlib.patches import Patch as Patch
from matplotlib.lines import Line2D as Line2D
import matplotlib as _mpl
from matplotlib.colors import LinearSegmentedColormap as _LinearSegmentedColormap
from matplotlib.colors import BoundaryNorm as _BoundaryNorm
from matplotlib.ticker import NullLocator as _NullLocator
from matplotlib.ticker import NullFormatter as _NullFormatter
from matplotlib.ticker import FuncFormatter as _FuncFormatter
import warnings as _warnings

# Sequential (non-diverging) colormap with the look of matplotlib's 'PuOr':
# light orange -> dark orange -> purple -> dark purple, built from PuOr's own
# colours but skipping its white centre, so lightness decreases monotonically.
# Used by default for the cluster-size panel, whose values are non-negative.
# Registered as 'PuOr_seq' so it can also be passed by name.
PuOr_seq = _LinearSegmentedColormap.from_list(
    'PuOr_seq',
    ['#fdb863', '#e08214', '#b35806', '#8073ac', '#542788', '#2d004b'])
if 'PuOr_seq' not in _mpl.colormaps:
    _mpl.colormaps.register(PuOr_seq)

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

    If ``plot_pulse=True``, a panel with the phase-resolved total counts of
    the data and of the posterior-expected signal (with posterior-predictive
    error bars) is added on top.

    If ``plot_clusters=True``, the figure also contains two lower rows for basic
    diagnostics of residuals (via analysis of clusters):

        * a cluster panel sharing the phase axis of the panels above: over the
          first phase cycle, the standardised residuals with
          ``|residual| >= threshlim``; over the second cycle, the size (in
          pixels) of the cluster each of these residuals belongs to. Clusters
          are groups of adjacent pixels (see ``cluster_connectivity``),
          merged across the periodic phase boundary.
        * on the left, the residual distribution (density) with a fitted
          normal distribution, a reference normal distribution (by default
          N(0,1)) and the p-value of a Kolmogorov-Smirnov test against the
          reference; on the right, the cumulative number of clusters with size
          at least :math:`s`, compared with the range expected from counting
          noise alone (Poisson realisations of the posterior-expected counts),
          and the p-value of the largest cluster.

    Statistical tests of 2D structure in the residuals (power spectrum,
    variogram, autocorrelation, ...) are provided by the separate
    :class:`~xpsi.PostProcessing.ResidualStructurePlot`, which can be passed
    together with this plot.

    The following example is (improved) from `Riley et al. 2019 <https://scixplorer.org/abs/2019ApJ...887L..21R/abstract>`_:

    .. image:: _static/_residualplot.png

    Each instance draws one figure: create a new instance for each call of
    :meth:`~xpsi.PostProcessing.SignalPlotter.plot` (reusing an instance that
    has already been plotted raises a ``RuntimeError``).

    :param str data_cmap:
         Colormap name from :mod:`matplotlib` to use for the data count
         numbers over joint channel-phase intervals.

    :param str model_cmap:
         Colormap name from :mod:`matplotlib` to use for the
         posterior-expected count numbers over joint channel-phase intervals.

    :param str residual_cmap:
         Colormap name from :mod:`matplotlib` to use for the standardised
         residuals (residual panel and, with ``plot_clusters=True``, the
         residuals above ``threshlim``). A diverging colormap is recommended.

    :param list parameters_vector:
         If given, the plot is made for this single vector of model parameter
         values (in the order of the likelihood parameters) instead of for
         posterior samples. :class:`~xpsi.PostProcessing.SignalPlotter` only
         looks at the first plot passed for a posterior, so this then applies
         to all plots passed together.

    :param bool plot_pulse:
         Add a top panel with the pulse profile (counts summed over channels)
         of the data and of the posterior-expected signal, with a
         posterior-predictive band (see ``n_realisation_per_sample``)?

    :param tuple legend_pulse_coords:
         Position of the legend of the pulse panel, passed to
         :mod:`matplotlib` as ``bbox_to_anchor``. If not given, the best
         location is chosen automatically.

    :param int n_realisation_per_sample:
         Number of Poisson realisations drawn from the pulse profile of each
         posterior sample, to build the posterior-predictive band (16th to
         84th percentiles) of the pulse panel (``plot_pulse=True``). If fewer
         than 100 realisations are available in total (number of posterior
         samples times this number), the band is instead the
         posterior-expected pulse plus or minus its square root (Poisson
         error).
         Default 1.

    :param bool blur_residuals:
         Smooth the displayed residual map with a Gaussian filter? This only
         affects the display; clusters and statistics use the unsmoothed
         residuals.

    :param bool plot_clusters:
         Add the cluster panel and the distribution panels (residual
         distribution and cluster-size distribution) described above? The
         parameters ``threshlim``, ``clusters_cmap``, ``mu``, ``sigma``,
         ``nbins``, ``cluster_null_n``, ``cluster_null_random_seed`` and
         ``cluster_connectivity`` are only used if ``True``.

    :param float threshlim:
         Threshold on the absolute value of the standardised residuals:
         pixels with ``|residual| >= threshlim`` are grouped into clusters.
         Default 2.0.

    :param str clusters_cmap:
         Colormap name from :mod:`matplotlib` to use for the cluster sizes.
         Cluster sizes are non-negative, so a sequential colormap is
         recommended. The default, ``'PuOr_seq'``, is a sequential colormap
         with the colours of ``'PuOr'`` (light orange for small clusters to
         dark purple for large ones), registered with :mod:`matplotlib` by
         this module. When all clusters are smaller than 7 pixels, the colour
         bar is discrete, with one colour for each size from 1 to 6. Pixels
         outside any cluster are left blank.

    :param str clustdist_cmap:
         Deprecated and ignored: the cluster-size distribution panel has no
         colormap. Passing it issues a ``FutureWarning``; it will be removed
         in a future release.

    :param float mu:
         Mean of the reference normal distribution (dashed line) drawn over
         the residual distribution. Default 0, as expected for standardised
         residuals. A maximum-likelihood normal fitted to the residuals is
         always drawn as well (solid line); its parameters are stored in
         :attr:`residual_fit_mu` and :attr:`residual_fit_sigma`. A one-sample
         Kolmogorov-Smirnov test of the residuals against the reference
         distribution is reported in the legend; its statistic and p-value are
         stored in :attr:`residual_ref_ks_stat` and
         :attr:`residual_ref_ks_pvalue`.

    :param float sigma:
         Standard deviation of the reference normal distribution (dashed
         line). Default 1, as expected for standardised residuals.

    :param int nbins:
         Number of bins of the residual distribution histogram. The bins span
         a range symmetric about zero. Default 50.

    :param int cluster_null_n:
         Number of Poisson realisations of the posterior-expected counts used
         to build the noise-only reference of the cluster-size distribution
         panel. Each realisation is standardised and clustered exactly like
         the data. Default 1000; the smallest reportable p-value for the
         largest cluster is ``1/(cluster_null_n+1)``. The size of the largest
         cluster and its p-value are stored in :attr:`cluster_max_size` and
         :attr:`cluster_max_pvalue`.

    :param int cluster_null_random_seed:
         Seed for the cluster-size noise realisations (reproducibility).

    :param int cluster_connectivity:
         Pixel connectivity used to group ``|residual| >= threshlim`` pixels
         into clusters: ``8`` (default) joins pixels sharing an edge or a
         corner, so oblique structures stay in one cluster; ``4`` joins
         pixels sharing an edge only. Applied identically to the data and to
         the noise realisations, including across the periodic phase
         boundary.

    :param kwargs:
         Keyword arguments for :class:`SignalPlot` (e.g. ``write``, ``dpi``,
         ``fig_dir``, ``root_filename``), and ``yscale`` (``'log'`` by
         default) for the channel axis of the channel-phase panels.

    """

    __figtype__ = 'signalplot_residuals'

    # do not change at runtime (see base class comments):
    __caching_targets__ = ['expected_counts']

    # When the largest cluster is smaller than this, the cluster-size colour
    # bar is discrete with one colour for each size from 1 to this - 1;
    # otherwise it is continuous
    __discrete_cluster_sizes_below__ = 7

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
                 clusters_cmap='PuOr_seq',
                 clustdist_cmap=None,
                 mu=0.0,
                 sigma=1.0,
                 nbins=50,
                 cluster_null_n=1000,
                 cluster_null_random_seed=0,
                 cluster_connectivity=8,
                 **kwargs):
        """ Constructor method for plotting residuals.

        See the class docstring for the parameters.
        """
        # The 2D structure tests moved to their own plot class; SignalPlot
        # silently ignores unknown keyword arguments, so fail loudly instead
        moved = sorted(k for k in kwargs
                       if k == 'analyse_residual_structure'
                       or k.startswith('structure_'))
        if moved:
            raise TypeError('ResidualPlot: {} moved to the separate plot '
                            'class xpsi.ResidualStructurePlot, which can be '
                            'passed together with ResidualPlot.'
                            .format(', '.join(moved)))

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
        if not parameters_vector is None:
            self.parameters_vector = parameters_vector

        if cluster_connectivity not in (4, 8):
            raise ValueError('cluster_connectivity must be 4 or 8.')
        if clustdist_cmap is not None:
            _warnings.warn('ResidualPlot: clustdist_cmap is deprecated and '
                           'ignored (the cluster-size distribution panel has '
                           'no colormap); it will be removed in a future '
                           'release.', FutureWarning, stacklevel=3)

        # Start from the default layout: the class attributes below are
        # modified for this instance's figure and restored at the end of
        # __init__, but an exception raised in between (e.g. in a notebook)
        # would otherwise leave them modified and break the next instance
        self._reset_layout()

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
            # Two panel rows (clusters, distributions), each preceded by an
            # empty spacer row.  All rows are separated by the same small gap
            # as without clusters (panels without x labels above them); the
            # spacers add room only where it is needed: for the residual
            # panel's x labels and the colour bars above the cluster panel,
            # and for the cluster panel's x labels and the titles of the
            # distribution panels.  Gaps and spacers are in panel heights.
            gap                   = 0.05
            spacer_colourbars     = 0.45
            spacer_titles         = 0.10
            ratios = ([1] * (self._resid_row + 1)
                      + [spacer_colourbars, 1, spacer_titles, 1])
            cls.__rows__          += 2
            cls.__ax_rows__       = len(ratios)
            cls.__ax_columns__    += 3
            cls.__height_ratios__ = ratios
            cls.__width_ratios__  = [50, 1, 50, 1, 1] # second column for colorbars
            cls.__wspace__        = 0.1
            # hspace is a fraction of the AVERAGE row height, which the
            # spacers lower: rescale it so the gap is `gap` panel heights
            cls.__hspace__        = gap * len(ratios) / sum(ratios)

            # Add parameters to plots
            self._threshlim = threshlim
            self._mu=mu
            self._sigma=sigma
            self._nbins=nbins
            self._clusters_cmap = clusters_cmap
            self._cluster_null_n = cluster_null_n
            self._cluster_null_random_seed = cluster_null_random_seed
            self._cluster_connectivity = cluster_connectivity

            # Get the rows
            self._cluster_row = self._resid_row + 2
            self._clustdist_row = self._resid_row + 4

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
            # A single axis spanning the same columns as the panels above:
            # |residuals| > threshlim over phases [0, 1], cluster sizes over
            # [1, 2].  Each half gets a horizontal colour bar above it.
            self._ax_clres = self._fig.add_subplot(self._gs[self._cluster_row,:-1])
            self._ax_clres_cb = self._ax_clres.inset_axes([0.05, 1.04, 0.40, 0.05])
            self._ax_clust_cb = self._ax_clres.inset_axes([0.55, 1.04, 0.40, 0.05])
            
            self._ax_rdist = self._fig.add_subplot(self._gs[self._clustdist_row, :2])
            self._ax_cdist = self._fig.add_subplot(self._gs[self._clustdist_row, 2:-1])
            
            # Prettify (the data and model panels keep sharing the phase
            # axis labelled below the residual panel)
            self._ax_clres.set_xlabel(r'$\phi$ [cycles]')

            self._ax_clres.set_ylabel('channel')
            self._ax_clres.xaxis.set_major_locator(MultipleLocator(0.2))
            self._ax_clres.xaxis.set_minor_locator(MultipleLocator(0.05))
            self._ax_clres.set_xlim([0.0,2.0])
            
            self._ax_rdist.set_title('Residual distribution')  
            # At most ~8 "nice" ticks, whatever the residual range
            self._ax_rdist.xaxis.set_major_locator(MaxNLocator(nbins=8, steps=[1, 2, 2.5, 5, 10]))
            self._ax_rdist.xaxis.set_minor_locator(AutoMinorLocator())
            self._ax_rdist.set_xlabel('Residuals')

            self._ax_cdist.set_title(r'Distribution of cluster size $s$')
            self._ax_cdist.set_xlabel(r'Cluster size $s$ [pixels]')
            self._ax_cdist.set_ylabel(r'Number of clusters with size $\geq s$')
            self._ax_cdist.yaxis.set_label_position('right')
            self._ax_cdist.tick_params(axis='y', which='both',
                                       labelright=True, labelleft=False,
                                       right=True, left=False)

        if "yscale" in kwargs:
            self.yscale = kwargs.get("yscale")
        else:
            self.yscale = "log"

        # Restore sizes
        self._reset_layout()

        plt.close()

    @classmethod
    def _reset_layout(cls):
        """ Reset the figure-layout class attributes to their defaults
        (data, model and residual rows only). """
        cls.__rows__          = 3
        cls.__columns__       = 1
        cls.__ax_rows__       = 3
        cls.__ax_columns__    = 2
        cls.__height_ratios__ = [1]*3
        cls.__width_ratios__  = [50, 1] # second column for colorbars
        cls.__wspace__        = 0.025
        cls.__hspace__        = 0.125

    @make_verbose('ResidualPlot object iterating over samples',
                  'ResidualPlot object finished iterating')
    def execute(self, thetas, wrapper):
        """ Loop over posterior samples. """
        self._check_not_reused()
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

    def _check_not_reused(self):
        """ Raise if this instance has already been plotted.

        The figure is created once, in the constructor, and every run draws
        onto it, so a second run would overlay a second set of curves and
        reuse results stored by the first run.
        """
        if getattr(self, '_finalized', False):
            raise RuntimeError(
                'ResidualPlot instance has already been plotted; create a new '
                'ResidualPlot instance for each plotting call (e.g. instantiate the '
                'plot in the same notebook cell as the plotting call).')

    @make_verbose('ResidualPlot object finalizing',
                  'ResidualPlot object finalized')
    def finalize(self):
        self._check_not_reused()
        self._add_data()
        self._add_expected_counts()
        self._add_residuals()
        if self._plot_clusters:
            self._compute_clusters()
            self._reveal_clusters()
            self._disp_distributions()
        if self._plot_pulse:
            self.add_pulses()
        self._finalized = True

    def _compute_residuals(self):
        """ Compute and store the standardised residuals
        :math:`(c_{ik}-d_{ik})/\\sqrt{c_{ik}}` in ``self._residuals``. """
        self._residuals = ((self.expected_counts - self._signal.data.counts)
                           / _np.sqrt(self.expected_counts))
        return self._residuals

    def _compute_clusters(self):
        """ Identify clusters of residuals with ``|residual| >= threshlim``.

        Stores the threshold mask (``self._clusteresid``), the cluster
        labels and their number (``self._lw``, ``self._num``), the size of
        each labelled cluster (``self._clustarea``, index 0 being the
        background) and the per-pixel size of the cluster each pixel
        belongs to (``self._affectedarea``, 0 outside clusters).

        Pixels are grouped with ``cluster_connectivity`` (8 by default) and,
        rotational phase being periodic, clusters touching across the phase
        boundary are merged into a single cluster (see
        :meth:`_label_clusters`).
        """
        if not hasattr(self, '_residuals'):
            self._compute_residuals()
        self._clusteresid, self._lw, self._num = self._label_clusters(
            self._residuals, self._threshlim, self._cluster_connectivity)
        self._clustarea = _ndi_sum_labels(self._clusteresid, self._lw, index=_np.arange(self._lw.max() + 1))
        self._affectedarea = self._clustarea[self._lw]

    @classmethod
    def _label_clusters(cls, residuals, threshlim, connectivity=8):
        """ Label clusters of pixels with ``|residual| >= threshlim``.

        This is the single definition of a cluster, used for both the data
        and the noise realisations.

        :param residuals: 2-D array (channel x phase) of standardised
            residuals.
        :param float threshlim: Threshold on ``|residual|``.
        :param int connectivity: ``8`` joins pixels sharing an edge or a
            corner; ``4`` joins pixels sharing an edge only.  Clusters
            touching across the periodic phase boundary are merged with the
            same connectivity (:meth:`_wrap_phase_labels`).
        :returns: ``(mask, labels, num)``: the boolean threshold mask, the
            integer label of each pixel (0 for the background) and the
            number of clusters.
        """
        mask = _np.abs(residuals) >= threshlim
        structure = _np.ones((3, 3), dtype=bool) if connectivity == 8 else None
        labels, num = cls._wrap_phase_labels(
            mask, *_ndi_label(mask, structure=structure),
            connectivity=connectivity)
        return mask, labels, num

    def _compute_cluster_size_null(self):
        """ Cluster sizes expected from counting noise alone.

        Draws ``cluster_null_n`` Poisson realisations of the
        posterior-expected counts, standardises each one like the data
        (:meth:`_compute_residuals`) and clusters it like the data
        (:meth:`_label_clusters`).  Stores the list of cluster-size arrays,
        one per realisation, in ``self._cluster_null_sizes``.
        """
        rng = _np.random.default_rng(self._cluster_null_random_seed)
        model = self.expected_counts
        sqrt_model = _np.sqrt(model)
        self._cluster_null_sizes = []
        for _ in range(self._cluster_null_n):
            resid = (model - rng.poisson(model)) / sqrt_model
            _, labels, _ = self._label_clusters(resid, self._threshlim,
                                                self._cluster_connectivity)
            self._cluster_null_sizes.append(_np.bincount(labels.ravel())[1:])

    @staticmethod
    def _wrap_phase_labels(mask, labels, num, connectivity=8):
        """ Merge cluster labels across the periodic phase boundary.

        :param mask: 2-D boolean array (channel x phase) of flagged pixels.
        :param labels: 2-D integer array of cluster labels for ``mask``
            (0 for the background), e.g. from :func:`scipy.ndimage.label`.
        :param int num: Number of clusters in ``labels``.
        :param int connectivity: ``4``: a pixel in the last phase bin
            touches the first phase bin of the same channel only; ``8``:
            also of the adjacent channels (diagonal neighbours).  The
            channel axis is not periodic.
        :returns: ``(labels, num)`` after merging every pair of clusters
            touching across phase = 1 -> 0, relabelled consecutively from 1
            (0 remains the background).
        """
        # Pairs of (channel at the last phase bin, channel at the first
        # phase bin) that are both flagged and touch across the boundary
        last, first = mask[:, -1], mask[:, 0]
        pairs = [(i, i) for i in _np.nonzero(last & first)[0]]
        if connectivity == 8:
            pairs += [(i, i + 1) for i in _np.nonzero(last[:-1] & first[1:])[0]]
            pairs += [(i + 1, i) for i in _np.nonzero(last[1:] & first[:-1])[0]]
        if not pairs:
            return labels, num

        # Union-find over labels; a cluster can wrap on several channels
        # and so be linked to more than one cluster on the other side
        parent = _np.arange(num + 1)

        def _root(k):
            while parent[k] != k:
                parent[k] = parent[parent[k]]
                k = parent[k]
            return k

        for i, j in pairs:
            a, b = _root(labels[i, -1]), _root(labels[j, 0])
            if a != b:
                parent[max(a, b)] = min(a, b)

        roots = _np.array([_root(k) for k in range(num + 1)])
        # np.unique sorts, so the background root 0 keeps index 0
        _, relabel = _np.unique(roots, return_inverse=True)
        return relabel[labels], int(relabel.max())

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

        channel_edges = self._get_channel_edges()

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

        channel_edges = self._get_channel_edges()

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

        self._compute_residuals()

        vmax =  _np.max( _np.abs( self._residuals ) )

        channel_edges = self._get_channel_edges()

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
        """ Display clusters from residuals in the fourth panel.

        Calls :meth:`_compute_clusters` if the clusters have not been
        computed yet. """
        if not hasattr(self, '_affectedarea'):
            self._compute_clusters()
        vmaxresid =  _np.max( _np.abs( self._residuals ) )
        vmaxarea =  _np.max( self._affectedarea )

        channel_edges = self._get_channel_edges()
        
        clust1 = self._ax_clres.pcolormesh(self._signal.data.phases,
                                      channel_edges,
                                      # mask pixels below the threshold so
                                      # they are left blank, like the
                                      # cluster-size half
                                      _np.ma.masked_where(~self._clusteresid,
                                                          self._residuals),
                                      cmap = plt.get_cmap(self._residual_cmap),
                                      vmin = -vmaxresid,
                                      vmax = vmaxresid,
                                      linewidth = 0,
                                      rasterized = self._rasterized)
        clust1.set_edgecolor('face')

        # Mask pixels outside any cluster (size 0) so they are left blank,
        # as in the |residuals| panel, and map sizes from 1 (single pixel)
        # to the largest cluster.  max(..., 1) keeps vmin <= vmax when there
        # is no cluster at all.
        nsizes = int(max(vmaxarea, 1))
        # Small clusters only: a fixed discrete scale with one colour per
        # integer size from 1 to __discrete_cluster_sizes_below__ - 1
        # (whatever the largest size actually present), with bin edges at
        # half-integers so each colour is centred on its size
        clust_discrete = nsizes < self.__discrete_cluster_sizes_below__
        if clust_discrete:
            nsizes = self.__discrete_cluster_sizes_below__ - 1
            clust_cmap = plt.get_cmap(self._clusters_cmap).resampled(nsizes)
            clust_scale = dict(norm = _BoundaryNorm(
                                   _np.arange(0.5, nsizes + 1.0), nsizes))
        else:
            clust_cmap = plt.get_cmap(self._clusters_cmap)
            clust_scale = dict(vmin = 1, vmax = nsizes)
        clust2 = self._ax_clres.pcolormesh(self._signal.data.phases + 1.0,
                                      channel_edges,
                                      _np.ma.masked_equal(self._affectedarea, 0),
                                      cmap = clust_cmap,
                                      linewidth = 0,
                                      rasterized = self._rasterized,
                                      **clust_scale)
        clust2.set_edgecolor('face')

        self._ax_clres.axvline(1.0, lw=self._tick_width, color='k')
        self._ax_clres.set_ylim([_np.max([channel_edges[0],0.001]),
                                 channel_edges[-1]])
        self._ax_clres.set_yscale(self.yscale)

        # Repeat the channel axis (ticks and label) on the right, next to the
        # cluster-size half; a secondary axis follows the parent's scale
        self._ax_clres_right = self._ax_clres.secondary_yaxis('right')
        self._ax_clres_right.set_ylabel('channel')

        # Horizontal colour bars above each half, titled by their label
        self._clres_cb = self._fig.colorbar(clust1, cax = self._ax_clres_cb,
                                      orientation='horizontal',
                                      ticks=AutoLocator())
        self._clust_cb = self._fig.colorbar(clust2, cax = self._ax_clust_cb,
                                      orientation='horizontal',
                                      ticks=(_np.arange(1, nsizes + 1)
                                             if clust_discrete
                                             else AutoLocator()))
        for cb, label in ((self._clres_cb,
                           # mathtext braces clash with str.format: concatenate
                           r'Residuals $\left|(c_{ik}-d_{ik})/\sqrt{c_{ik}}\right|$ > '
                           + f'{self._threshlim}'),
                          (self._clust_cb, r'Cluster size $s$')):
            cb.ax.set_frame_on(True)
            if cb is self._clust_cb and clust_discrete:
                cb.ax.xaxis.set_minor_locator(_NullLocator())
            else:
                cb.ax.xaxis.set_minor_locator(AutoMinorLocator())
            cb.ax.xaxis.set_ticks_position('top')
            cb.ax.xaxis.set_label_position('top')
            cb.set_label(label=label, labelpad=8)
        
    def _disp_distributions(self):
        """ Display residual and cluster distributions in the fifth panel.

        Calls :meth:`_compute_clusters` if the clusters have not been
        computed yet. """
        if not hasattr(self, '_affectedarea'):
            self._compute_clusters()

        # Residual distribution, on a density scale (independent of the
        # number of channel-phase pixels), over a range symmetric about 0
        flat = self._residuals.ravel()
        amax = _np.max(_np.abs(flat))
        rdens, _, _ = self._ax_rdist.hist(flat,
                            bins=self._nbins,
                            range=(-amax, amax),
                            density=True,
                            histtype='stepfilled',
                            alpha=0.5,
                            color='steelblue',
                            label='Residuals')

        # Maximum-likelihood normal fitted to the residuals
        x = _np.linspace(-amax, amax, 400)
        fit_mu, fit_sigma = _stats.norm.fit(flat)
        self.residual_fit_mu = float(fit_mu)
        self.residual_fit_sigma = float(fit_sigma)
        self._ax_rdist.plot(x, _stats.norm.pdf(x, fit_mu, fit_sigma),
                            color='darkblue', linewidth=2,
                            label=(r'Fitted normal ($\mu$=' + f'{fit_mu:.3f}, '
                                   + r'$\sigma$=' + f'{fit_sigma:.3f})'))

        # Reference normal expected for well-modelled standardised residuals,
        # with a one-sample KS test of the (unbinned) residuals against it.
        # A frozen CDF is passed because scipy maps the string 'norm' to a
        # function that does not accept loc/scale arguments.
        ks_stat, ks_p = _stats.kstest(
            flat, _stats.norm(loc=self._mu, scale=self._sigma).cdf)
        self.residual_ref_ks_stat = float(ks_stat)
        self.residual_ref_ks_pvalue = float(ks_p)
        ks_p_str = f'{ks_p:.3f}' if ks_p >= 1e-3 else f'{ks_p:.1e}'
        self._ax_rdist.plot(x, _stats.norm.pdf(x, self._mu, self._sigma),
                            color='m', linewidth=2, linestyle='--',
                            label=(r'Reference $\mathcal{N}$('
                                   + f'{self._mu:g}, {self._sigma:g}), '
                                   + f'KS test p-value = {ks_p_str}'))

        # Leave headroom above the histogram and curves for the legend
        ytop = max(_np.max(rdens),
                   _np.max(_stats.norm.pdf(x, fit_mu, fit_sigma)),
                   _np.max(_stats.norm.pdf(x, self._mu, self._sigma)))
        self._ax_rdist.set_ylim(0.0, 1.5 * ytop)

        self._ax_rdist.set_ylabel('Density')
        self._ax_rdist.legend(fontsize='x-small', frameon=False, loc='upper left')

        # Cluster-size distribution: cumulative number of clusters with size
        # >= s, compared with the band expected from counting noise alone
        if not hasattr(self, '_cluster_null_sizes'):
            self._compute_cluster_size_null()
        null_sizes = self._cluster_null_sizes
        obs_sizes = self._clustarea[1:].astype(int)   # index 0 = background
        obs_max = int(obs_sizes.max()) if obs_sizes.size else 0
        null_max = _np.array([z.max() if z.size else 0 for z in null_sizes])

        s_grid = _np.arange(1, max(obs_max, int(null_max.max()), 1) + 2)
        def _ccdf(sizes):
            return _np.count_nonzero(sizes[None, :] >= s_grid[:, None], axis=1)
        null_ccdf = _np.array([_ccdf(z) for z in null_sizes])
        q = _np.percentile(null_ccdf, [2.5, 16, 50, 84, 97.5], axis=0)

        # One-sided p-value of the largest cluster (add-one corrected)
        n_ge = int(_np.count_nonzero(null_max >= obs_max))
        p_max = (n_ge + 1.0) / (len(null_sizes) + 1.0)
        self.cluster_max_size = obs_max
        self.cluster_max_pvalue = float(p_max)

        ax = self._ax_cdist
        ax.fill_between(s_grid, q[0], q[4], step='post', color='0.85',
                        linewidth=0, label='noise 95%')
        ax.fill_between(s_grid, q[1], q[3], step='post', color='0.65',
                        linewidth=0, label='noise 68%')
        ax.step(s_grid, q[2], where='post', color='0.3', linestyle='--',
                linewidth=1, label='noise median')
        ax.step(s_grid, _ccdf(obs_sizes), where='post', color=PuOr_seq(0.95),
                linewidth=2, label='data')
        if obs_max > 0:
            p_str = (r'$\leq$' if n_ge == 0 else '= ') + f'{p_max:.2g}'
            ax.axvline(obs_max, color=PuOr_seq(0.2), linestyle=':',
                       linewidth=1.5,
                       label=f'largest = {obs_max} (p {p_str})')
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_xlim(0.9, 1.3 * s_grid[-1])
        # Plain-number size ticks; label minor ticks (2, 3, ...) only when
        # the size range is narrow, where there are no other ticks
        _int_ticks = _FuncFormatter(
            lambda v, _: f'{v:g}' if (v >= 1 and float(v).is_integer()) else '')
        ax.xaxis.set_major_formatter(_int_ticks)
        ax.xaxis.set_minor_formatter(
            _int_ticks if s_grid[-1] < 10 else _NullFormatter())
        # Keep the curves in the lower ~55% of the (log) axis, leaving room
        # for the legend above them
        ybottom = 0.5
        ycurves = max(q[4].max(), _ccdf(obs_sizes).max(), 1)
        ax.set_ylim(ybottom, ybottom * (ycurves / ybottom) ** (1.0 / 0.55))
        ax.legend(fontsize='x-small', frameon=False, loc='upper right')

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
