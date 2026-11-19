from ._global_imports import *
from ._signalplot import SignalPlot
from scipy.stats import norm, kstest

class Residual1DPlot(SignalPlot):
    """ Plot the count data, the posterior-expected count signal, and residuals.

    The figure contain one panel :

        * histogram of the standardised residuals between the
          data and posterior-expected count signal over joint channel-phase
          intervals.

    The following example is (improved) from `Riley et al. 2019 <https://scixplorer.org/abs/2019ApJ...887L..21R/abstract>`_:

    .. image:: _static/_1d_residualplot.png

    :param ... Todo

    """

    __figtype__ = 'signalplot_1dresiduals'

    # do not change at runtime (see base class comments):
    __caching_targets__ = ['expected_counts']

    __rows__          = 1
    __columns__       = 1
    __ax_rows__       = 1
    __ax_columns__    = 1
    __height_ratios__ = [1]
    __width_ratios__  = [1] # second column for colorbars
    __wspace__        = 0.025
    __hspace__        = 0.125

    @make_verbose('Instantiating a 1D residual plotter for residual distribution checking',
                  'Residual plotter instantiated')
    def __init__(self,
                 nbins=20,
                 plot_fit=False,
                 mu=0.0,
                 sigma=1.0,
                 parameters_vector=None,
                 **kwargs):
        """
        Constructor method for plotting the 1D residual distribution.

        :param int nbins:
             Number of histogram bins.

        :param bool plot_fit:
             Overlay a maximum-likelihood normal fitted to the residuals?

        :param float mu:
             Mean of the reference normal distribution (dashed line) drawn
             over the histogram.  Default 0, as expected for standardised
             residuals.  A one-sample Kolmogorov-Smirnov test of the residuals
             against the reference is reported in the legend; its statistic
             and p-value are stored in :attr:`residual_ref_ks_stat` and
             :attr:`residual_ref_ks_pvalue`.

        :param float sigma:
             Standard deviation of the reference normal distribution.
             Default 1, as expected for standardised residuals.

        :param list parameters_vector:
             List of model parameters to plot the residuals for.

        :param kwargs:
             Keyword arguments for :class:`SignalPlot`.
        """
        super(Residual1DPlot, self).__init__(**kwargs)

        # Get the input parameters
        self._nbins = int( nbins )
        self._plot_fit = plot_fit
        self._mu = mu
        self._sigma = sigma
        if not parameters_vector is None:
            self.parameters_vector = parameters_vector

         # Get on plotting
        self._get_figure()
        self._ax_1d_residuals = self._add_subplot(0,0)
        self._ax_1d_residuals.set_xlabel(r'$(c_{ik}-d_{ik})/\sqrt{c_{ik}}$')
        self._ax_1d_residuals.set_ylabel('Occurence')

        plt.close()

    @make_verbose('Residual1DPlot object iterating over samples',
                  'Residual1DPlot object finished iterating')
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

    def _check_not_reused(self):
        """ Raise if this instance has already been plotted.

        The figure is created once, in the constructor, and every run draws
        onto it, so a second run would overlay a second set of curves and
        reuse results stored by the first run.
        """
        if getattr(self, '_finalized', False):
            raise RuntimeError(
                'Residual1DPlot instance has already been plotted; create a new '
                'Residual1DPlot instance for each plotting call (e.g. instantiate the '
                'plot in the same notebook cell as the plotting call).')

    @make_verbose('ResidualPlot object finalizing',
                  'ResidualPlot object finalized')
    def finalize(self):
        self._check_not_reused()
        self._add_1d_residuals()
        self._finalized = True

    def _add_1d_residuals(self):
        """ Display the residuals in the third panel. """

        # Get the residuals 
        self._residuals = self.expected_counts - self._signal.data.counts
        self._residuals /= _np.sqrt(self.expected_counts)
        self._residuals = self._residuals.flatten()

        # Plot them
        rdens, _, _ = self._ax_1d_residuals.hist(self._residuals, bins=self._nbins, density=True, alpha=0.5)

        xmin, xmax = _np.min( self._residuals ) , _np.max( self._residuals )
        x = _np.linspace(xmin, xmax, 400)
        ytop = _np.max(rdens)

        # Add the fit if needed
        if self._plot_fit:

            mean,std=norm.fit( self._residuals )
            y = norm.pdf(x, mean, std)
            ytop = max(ytop, _np.max(y))
            self._ax_1d_residuals.plot( x , y , label=r'Fitted Normal $\sigma=$'+f'{std:.3f}', color='darkblue' )
            self._ax_1d_residuals.axvline(mean, color='k', linestyle='dashed', linewidth=2, label=r'$\mu=$'+f'{mean:.5f}')

        # Reference normal expected for well-modelled standardised residuals,
        # with a one-sample KS test of the (unbinned) residuals against it.
        # A frozen CDF is passed because scipy maps the string 'norm' to a
        # function that does not accept loc/scale arguments.
        ks_stat, ks_p = kstest(self._residuals,
                               norm(loc=self._mu, scale=self._sigma).cdf)
        self.residual_ref_ks_stat = float(ks_stat)
        self.residual_ref_ks_pvalue = float(ks_p)
        ks_p_str = f'{ks_p:.3f}' if ks_p >= 1e-3 else f'{ks_p:.1e}'
        y_ref = norm.pdf(x, self._mu, self._sigma)
        ytop = max(ytop, _np.max(y_ref))
        self._ax_1d_residuals.plot(x, y_ref,
                                   color='m', linewidth=2, linestyle='--',
                                   label=(r'Reference $\mathcal{N}$('
                                          + f'{self._mu:g}, {self._sigma:g}), '
                                          + f'KS test p-value = {ks_p_str}'))

        # Leave headroom above the histogram and curves for the legend
        self._ax_1d_residuals.set_ylim(0.0, 1.5 * ytop)
        self._ax_1d_residuals.legend(fontsize='small', frameon=False, loc='upper left')


        # Do the nice plotting
        self._ax_1d_residuals.set_title('Residual histogram')
        self._ax_1d_residuals.set_xlabel(r'$(c_{ik}-d_{ik})/\sqrt{c_{ik}}$')
        self._ax_1d_residuals.set_ylabel('Density')