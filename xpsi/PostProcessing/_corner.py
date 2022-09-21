from __future__ import division, print_function
import numpy as np_

from scipy.special import logsumexp

from ._global_imports import *

from . import _precision
from collections import OrderedDict

from getdist.plots import getSubplotPlotter
from getdist.mcsamples import MCSamples

try:
    from nestcheck.ns_run_utils import get_logw, get_w_rel
    from nestcheck.plots import bs_param_dists
    from nestcheck.error_analysis import run_ci_bootstrap
    from nestcheck.estimators import param_cred, logz
except ImportError:
    _warning('CornerPlotter instances cannot use nestcheck functionality.')
else:
    try:
        from nestcheck.plots import getdist_kde
    except ImportError:
        try:
            from nestcheck.plots import weighted_1d_gaussian_kde
        except ImportError:
            _warning('CornerPlotter instances cannot use nestcheck '
                     'functionality.')
        else:
            _warning('Using native nestcheck KDE instead of GetDist KDE.')

from ._backends import NestedBackend
from ._postprocessor import PostProcessor

class CornerPlotter(PostProcessor):
    """ Plot marginal posterior densities and estimators.

    """


    @fix_random_seed
    @make_verbose('Executing posterior density estimation',
                  'Posterior density estimation complete')
    def plot(self,
             params,
             IDs=None,
             combine=False,
             combine_all=False,
             only_combined=False,
             force_combine=True,
             overwrite_combined=False,
             bootstrap_estimators=True,
             bootstrap_density=False,
             separate_plots=False,
             write=False,
             root_filename='',
             directory='./',
             ext='.pdf',
             dpi=300,
             maxdots=2000,
             **kwargs):
        """ Generate posterior density plots.

        Invokes :mod:`getdist` and :mod:`nestcheck` for nested sampling runs.

        Up to five runs can be plotted natively via :mod:`nestcheck`; beyond
        such a number the plots generally display too much information and
        clarity is lost.

        :param list params:
            List of parameter strings for plotting. Must be shared by all
            posteriors selected with the ``IDs`` argument.

        :param OrderedDict IDs:
            Keys must be string identifiers of :class:`Runs` instances.
            Each dictionary element must be a list of string identifiers,
            each matching objects collected in :class:`Runs` instance
            corresponding to the key. Defaults to ``None``, meaning attempt to
            use as many runs as possible subject to plotting restrictions.

        .. note::

            The order of IDs is used to control the layering of posteriors.
            If there are multiple underlying posteriors (i.e., multiple
            dictionary keys), only one (combined) run per posterior is
            rendered, and the first posterior is rendered on the topmost layer.
            If there is only one underlying posterior (i.e., one dictionary
            keys), then the combined-sample posterior, if available, is plotted
            on the topmost layer, whilst the runs on that posterior are
            rendered on layers underneath in the order specified in the list.
            In either case, estimators are calculated and reported for the
            (combined) run on the topmost layer.

        :param bool combine:
            Additionally combine the runs into a single run for overplotting?
            The overriding setting if there is more than one underlying
            posterior to be plotted, is to attempt to combine runs on each
            posterior if multiple such runs are available, in order to
            reduce information density. If there is a single underlying
            posterior, the user-specified value is respected.

        :param bool combine_all:
            Combine all runs in each :class:`Runs` instance or only those
            for which IDs are provided? Ignored if ``combine`` is ``False``.

        :param bool force_combine:
            Force recombination of elligible run sets, even if a
            combined run is already cached?

        :param bool only_combined:
            Only plot the combined run? Only heeded if a single posterior
            is selected for plotting, and in that case is ignored if
            ``combine`` is ``False``.

        :param bool overwrite_combined:
            Overwrite combined-sample files on disk with the same filename?

        :param bool bootstrap:
            Use :mod:`nestcheck` and :mod:`fgivenx` to bootstrap the runs for
            posterior density error estimation?

        :param bool separate_plots:
            Generate a lower-triangle plot with :mod:`getdist`, and a separate
            error plot with :mod:`nestcheck` (with :mod:`fgivenx` and
            :mod:`getdist`). If ``False`` (default), the diagonal panels of the
            lower-triangle plot are modified by adding the nestcheck output.
            Ignored if ``bootstrap`` is ``False``.

        :param bool write:
            Export the figure?

        :param str root_filename:
            Root filename to prepend to automatically generated name. Can be,
            e.g., a model and/or data set identifier.

        :param str directory:
            If ``None`` defaults to current directory.

        :param str ext:
            File extension for writing. E.g., ``'.png'``.

        :param int dpi:
            Dots-per-square-inch settings for exporting plots.

        :param kwargs:
            Keyword arguments for the :meth:`_plot_density_with_error` and
            :meth:`_plot_triangle` methods. Keyword arguments for line
            properties (width, color, and alpha) for :mod:`getdist` contours and density
            distributions. If ``bootstrap and not separate_plots`` then
            the density distribution linewidth is set to zero if not
            explicitly specified with kwarg ``lw_1d``.
            In addition, keyword arguments for avoiding unnecessary re-drawing of prior samples (``force_draw``, ``prior_samples_fnames`` and ``priors_identical``).

        """
        self.set_subset(IDs, combine, combine_all,
                        force_combine, only_combined,
                        overwrite_combined)
        self.set_params(params)

        if bootstrap_density and not separate_plots:
            if 'lw_1d' not in kwargs: kwargs['lw_1d'] = 0.0
        self._set_line_and_contour_args(**kwargs)

        self._plotter = self._plot_triangle(bootstrap_estimators, **kwargs)

        if bootstrap_density and separate_plots:
            figs = self._plot_density_with_error(**kwargs)
        elif bootstrap_density:
            figs = self._plot_density_with_error(plotter=self._plotter,
                                                 **kwargs)
        if write:
            root_filename = (root_filename + '__' if root_filename else '') + \
                'posteriorDensity__runs_' + \
                '_'.join(str(ID).replace(' ', '') for
                         ID in self.get_attr('ID')) + '__'

            _dpi = dpi
            if maxdots > 0:
                ndots = dpi * len(self.params)
                ndots *= self._plotter.settings.subplot_size_inch
                if ndots > maxdots:
                    dpi = int(maxdots * _dpi / ndots)

            self._plotter.export(fname=root_filename+'triangle'+ext,
                                   adir=directory, dpi=dpi)
            try:
                figs[1].savefig(_os.path.join(directory,
                                              root_filename+'fthetas_1d.pdf'),
                                dpi=_dpi,
                                bbox_inches='tight')
            except IndexError:
                if separate_plots:
                    fname = root_filename + 'params_1d.pdf'
                else:
                    fname = root_filename + 'fthetas_1d.pdf'
                figs[0].savefig(_os.path.join(directory, fname),
                                dpi=_dpi,
                                bbox_inches='tight')
            except (TypeError, NameError):
                pass

        return self._plotter

    @make_verbose('Simulating nested sampling realisations for '
                  'posterior density error visualisation',
                  'Simulated nested sampling realisations and '
                  'plotted posterior density error estimates')
    def _plot_density_with_error(self,
                                 plotter = None,
                                 fthetas = None,
                                 kde_func = None,
                                 kde_kwargs = None,
                                 **kwargs):
        """
        :param plotter:
            A :attr:`getdist.GetDistPlotter` instance if the :mod:`nestcheck`
            output is to be displayed on a lower-triangle plot.

        :param list fthetas:
            Iterable containing functions of parameter vector for which
            density is to be plotted with :mod:`nestcheck` via :mod:`fgivenx`.
            The parameters themselves are handled automatically with ``lambda``
            functions. Additional functions are always plotted using the
            native :mod:`nestcheck` matplotlib figure; the parameter densities
            are be added to a :mod:`getdist` lower-triangle plot is supplied.

        :param func kde_func:
            Function for KDE compatible with :mod:`nestcheck.plots`. Must
            be *weighted* KDE (Higson et al. 2018). If ``None``, uses
            :mod:`getdist` if available, or the native KDE function otherwise.
            If using :mod:`getdist`, the KDE settings are automatically
            retrieved from the first run and applied to :mod:`nestcheck` and
            :mod:`fgivenx` *for all runs*.

        :param kwargs:
            Keyword arguments for :func:`nestcheck.plots.bs_param_dists`.

        TODO:
        ----

            * lims based on credible interval estimate for efficiency?

        """
        nestcheck_bcknds, runs = self._filter_nestcheck_compatible()

        nx = kwargs.pop('nx', 200); ny = kwargs.pop('ny', nx)
        scale_ymax = kwargs.pop('scale_ymax', 1.1)
        n_simulate = kwargs.pop('n_simulate', 200)

        params = self.params.names
        labels = self.params.labels

        # declare how to access parameter samples for each run
        _fthetas = []
        for run in runs:
            func = lambda y: (lambda theta: theta[:,y])
            _fthetas.append([func(run.get_index(param)) for param in params])

        # declare limits for parameters
        bounds = []
        for run in runs:
            bounds.append([list(run.bounds[param]) for param in params])

        _lims = [list(plotter.subplots[i,i].get_xlim()) for i in range(len(params))]
        lims = [_lims for run in runs]
        for _l, _b in zip(lims, bounds): # loop over runs
            for l, b in zip(_l,_b): # loop over parameters
                l[0] = (l[0] if l[0] > b[0] else b[0])
                l[1] = (l[1] if l[1] < b[1] else b[1])

        if kde_func is None:
            try:
                kde_func = getdist_kde
            except NameError:
                kde_func = weighted_1d_gaussian_kde
                kde_kwargs = [None] * len(runs)
            else:
                normalize = kwargs.pop('normalize', False)
                kde_kwargs = []
                for run in runs:
                    kde_kwargs.append(
                            {'settings': run.kde_settings,
                             'ranges': [run.bounds[param] for param in params],
                             'normalize': normalize}
                            )

        lines = kwargs.pop('lines', False)
        parallel = kwargs.pop('parallel', True)
        rasterize_contours = kwargs.pop('rasterize_contours', True)
        tqdm_kwargs = kwargs.pop('tqdm_kwargs', None)

        figsize = kwargs.pop('figsize', _np.array([6.0 * len(params),
                                                   3.0 * len(params)]))

        figs = []
        with verbose(plotter is not None,
                     'Adding density error information to triangle plot',
                     'Added density error information'):
            fig = bs_param_dists(nestcheck_bcknds,
                                 fthetas=_fthetas,
                                 kde_func=kde_func,
                                 kde_kwargs=kde_kwargs,
                                 ftheta_lims=lims,
                                 nx=nx,
                                 ny=ny,
                                 scale_ymax=scale_ymax,
                                 n_simulate=n_simulate,
                                 simulate_weights=True,
                                 getdist_plotter=plotter,
                                 figsize=figsize,
                                 lines=lines,
                                 parallel=parallel,
                                 rasterize_contours=rasterize_contours,
                                 labels=labels,
                                 no_means=True,
                                 tqdm_kwargs=tqdm_kwargs)

        if fig: figs.append(fig)

        if fthetas:
            if not isinstance(fthetas[0], list):
                num_funcs = len(fthetas[0])
            else:
                num_funcs = len(fthetas)
            kde_kwargs['ranges'] = ftheta_lims
            figsize *= float(num_funcs)/len(params)
            if 'ftheta_labels' in kwargs:
                kwargs = {'labels': kwargs['ftheta_labels']}
            else:
                kwargs = {}

            fig = bs_param_dists(nestcheck_bcknds,
                                 fthetas=fthetas,
                                 kde_func=kde_func,
                                 kde_kwargs=kde_kwargs,
                                 ftheta_lims=ftheta_lims,
                                 nx=nx,
                                 ny=ny,
                                 scale_ymax=scale_ymax,
                                 n_simulate=n_simulate,
                                 simulate_weights=True,
                                 figsize=figsize,
                                 lines=lines,
                                 parallel=parallel,
                                 rasterize_contours=rasterize_contours,
                                 **kwargs)
            figs.append(fig)

        return figs if figs else None

    @make_verbose('Constructing lower-triangle posterior density '
                  'plot via Gaussian KDE:',
                  'Constructed lower-triangle posterior density plot')
    def _plot_triangle(self,
                       bootstrap,
                       prior_density = True,
                       KL_divergence = True,
                       KL_base = 'bits',
                       ndraws = int(1e6),
                       param_plot_lims = None,
                       crosshairs = False,
                       filled = False,
                       legend_loc = 'lower left',
                       legend_corner_coords = (0.75,0.75),
                       legend_frameon = False,
                       scale_attrs = None,
                       normalize = True,
                       veneer = False,
                       no_zero = True,
                       no_ylabel = False,
                       label_right = True,
                       no_ytick = False,
                       credible_interval_1d = True,
                       credible_interval_1d_all_show = False,
                       annotate_credible_interval = True,
                       annotate_xy=(0.025,0.915),
                       sixtyeight = True,
                       ninety = False,
                       compute_all_intervals=True,
                       **kwargs):
        """ Call :meth:`getdist.plots.GetDistPlotter.triangle_plot`.

        :param bool prior_density:
            If ``True`` tries to draw samples from the joint prior and plot
            marginal 1D prior densit functions. Silently fails if attempt
            fails due to missing likelihood and prior callables.

        :param bool KL_divergence:
            If `True` and `prior_density` is `True`, estimate and annotate
            credible interval for Kullback-Leibler divergence for each
            parameter in triangle plot.

        :param str KL_base:
            Base for Kullback-Leibler divergence. Options are {'bits', 'nats'}.

        :param int ndraws:
            Number of samples drawn from the joint prior. Ignored if
            ``prior_density is False`` attempt to plot density fails.

        :param dict param_plot_lims:
            Dictionary of viewing ranges for plotting. Keys must be parameter
            names.

        :param bool crosshairs:
            Display parameter truth crosshairs?

        :param bool filled:
            Specify whether the contours are to be filled.

        :param str legend_loc:
            Specify the legend location. Defaults to ``upper right`` because the
            plot is a lower-triangle array.

        :param tuple legend_corner_coords:
            Modifies meaning of ``legend_loc`` to be the coordinates of the
            point on the legend box specified by ``legend_loc``. Pass ``None``
            if not applicable. Defaults to place legend in empty upper region
            of lower-triangle plot.

        :param dict scale_attrs:
            Scale :class:`getdist.plots.GetDistPlotterSettings` attributes
            from the automatic values. E.g., ``{'legend_fontsize': 2.0}``.
            Use string values to set the key attribute to the value attribute.
            Caution: do not rely on ordering of pairs in a dictionary, but
            use an :class:`collections.OrderedDict` instead to heed order.

        :param bool normalize:
            Normalize density distribution in on-diagonal 1D panels?

        :param bool no_zero:
            Hide axes zeros within on-diagonal 1D panels?

        :param bool no_ylabel:
            Hide *probability density* label on diagonal 1D marginal panels?

        :param bool label_right:
            Display *probability density* label on diagonal 1D marginal plots?

        :param bool no_ytick:
            Hide y-axis ticks on diagonal 1D marginal plots?

        :param bool credible_interval_1d:
            Estimate and plot 1D marginal credible intervals? The upper
            and lower quantiles of the interval are estimated via bootstrapping
            with :mod:`nestcheck`, each bounding quantile plotted as a (narrow)
            band bounded by the same quantiles with respect to the bootstrap
            realisations. The interval between the two bands is generally much
            wider and is shaded lighter.

        :param bool credible_interval_1d_all_show:
            Show the 1D marginal credible intervals for all the runs IDs.

        :param bool annotate_credible_interval:
            Annotate each on-diagonal panel with numeric credible interval
            as median +/- distance to symmetric quantiles. Each quantile,
            including the median, is estimated via bootstrapping with
            :mod:`nestcheck`, and the median of each quantile from the bootstrap
            realisations is used for the reported numeric credible interval.

        :param tuple annotate_xy:
            Coordinates as axis fractions for annotation of credible intervals.

        :param bool sixtyeight:
            Should the credible interval, which is symmetric in quantiles about
            the mean, be approximately the 1-\sigma credible interval thus
            containing ~68% of the posterior mass? If ``False`` the interval
            computed is the approximate 2-\sigma interval containing ~95% of
            the posterior mass.

        :param kwargs:

            * additional keyword arguments passed to
              :meth:`getdist.GetDistPlotter.triangle_plot`
            * settings for :mod:`getdist` posterior lower-triangle plotting, applied
              to a :class:`getdist.plots.GetDistPlotSettings` instance
            * setting for not drawing prior samples but reading them from a file instead or saving them in a file if no such already exists ``force_draw = [False,]``
            * setting for choosing a non-default file name to be used if samples are not re-drawn: ``prior_samples_fnames=["name.npy",]``
            * setting for plotting the priors only for one of the runs if they are known to be identical: ``priors_identical=True``

        .. note::

            Using ``subplot_size`` keyword argument (specify in inches) invokes
            automated label fontsizes and tick sizes. If ``width_inch`` is
            used instead, this automation does not occur.

        """

        self.credible_interval_1d_all_show=credible_interval_1d_all_show
        if credible_interval_1d_all_show:
            KL_divergence = False
        try:
            for run in self.subset_to_plot:
                if not isinstance(run, NestedBackend):
                    raise TypeError('Nested sampling backends are required.')
        except AttributeError:
            print('Nested sampling runs are required.')
            raise
        else:
            getdist_bcknds = self.get_attr('getdist_backend')
            getdist_bcknds.reverse()

            line_args = self.get_attr('lines')
            line_args.reverse()

            contour_args = self.get_attr('contours')
            contour_args.reverse()

            if len(getdist_bcknds) == 1:
                legend_labels = None
            elif len(self._subset) > 1:
                legend_labels = self.get_attr('parent_ID')
            else:
                legend_labels = self.get_attr('ID')

            if legend_labels is not None:
                legend_labels.reverse()

        if param_plot_lims is None:
            param_plot_lims = {}

        if param_plot_lims:
            prune = kwargs.get('tick_prune', None)

        # try to set matching :class:`getdist.plots.GetDistPlotSettings` attrs
        plotter = getSubplotPlotter(kwargs.pop('subplot_size', 2),
                                    kwargs.pop('width_inch', None))
        setattr(plotter.settings, 'progress', True)
        setattr(plotter.settings, 'norm_prob_label', 'Probability density')
        setattr(plotter.settings, 'prob_y_ticks', True)
        setattr(plotter.settings, 'thin_long_subplot_ticks', False)
        setattr(plotter.settings, 'tick_prune', None)

        for key in kwargs.copy():
            if hasattr(plotter.settings, key):
                setattr(plotter.settings, key, kwargs[key])
                del kwargs[key]

        if scale_attrs is None:
            scale_attrs = {}

        for key, value in scale_attrs.iteritems():
            if hasattr(plotter.settings, key):
                if isinstance(value, float) or isinstance(value, int):
                    setattr(plotter.settings, key,
                            getattr(plotter.settings, key) * value)
                elif isinstance(value, _six.string_types):
                    if hasattr(plotter.settings, value):
                        setattr(plotter.settings, key,
                                getattr(plotter.settings, value))

        if isinstance(normalize, bool):
            diag1d_kwargs = {'normalized': normalize}
        if isinstance(no_zero, bool):
            diag1d_kwargs['no_zero'] = no_zero
        if isinstance(no_ylabel, bool):
            diag1d_kwargs['no_ylabel'] = no_ylabel
        if isinstance(label_right, bool):
            diag1d_kwargs['label_right'] = label_right
        if isinstance(no_ytick, bool):
            diag1d_kwargs['no_ytick'] = no_ytick
        plotter.triangle_plot(getdist_bcknds,
                               legend_labels = legend_labels,
                               params = self.params.names,
                               filled = filled,
                               legend_loc = legend_loc,
                               line_args = line_args,
                               contour_args = contour_args,
                               diag1d_kwargs = diag1d_kwargs,
                               **kwargs)
        try:
            if legend_corner_coords:
                plotter.legend.set_bbox_to_anchor(legend_corner_coords)
        except AttributeError:
            pass
        else:
            plotter.legend.set_frame_on(legend_frameon)

        # add custom parameter plotting limits and updated autolocation
        with fragile(verbose(param_plot_lims,
                             'Applying bespoke parameter viewing intervals',
                             'Viewing intervals applied')) as condition:
            if not condition: fragile.Break

            params = self.params
            for param, l in param_plot_lims.items():
                j = params.names.index(param)
                for i in range(j,len(params.names)):
                    ax = plotter.subplots[i,j]
                    ax.set_xlim(l)
                    ax.xaxis.set_major_locator(_get_default_locator(prune))
                    ax.xaxis.set_minor_locator(AutoMinorLocator())

                for i in range(j):
                    ax = plotter.subplots[j,i]
                    ax.set_ylim(l)
                    ax.yaxis.set_major_locator(_get_default_locator(prune))
                    ax.yaxis.set_minor_locator(AutoMinorLocator())

            plotter.fig.canvas.draw() # ensure the new locators take effect
            for param in param_plot_lims.keys():
                j = params.names.index(param)

                # deal with x-axes
                axis = plotter.subplots[-1,j].xaxis
                xmin, xmax = axis.get_view_interval()
                width = xmax - xmin
                gap_wanted = width * plotter.settings.tight_gap_fraction
                tick = [x for x in axis.get_major_ticks() if xmin <= x.get_loc() <= xmax]

                if tick[0].get_loc() - xmin < gap_wanted:
                    tick[0].label1.set_visible(False)
                if xmax - tick[-1].get_loc() < gap_wanted:
                    tick[-1].label1.set_visible(False)

                # deal with y-axes
                if j > 0:
                    axis = plotter.subplots[j,0].yaxis
                    xmin, xmax = axis.get_view_interval()
                    width = xmax - xmin
                    gap_wanted = width * plotter.settings.tight_gap_fraction
                    tick = [x for x in axis.get_major_ticks() if xmin <= x.get_loc() <= xmax]

                    if tick[0].get_loc() - xmin < gap_wanted:
                        tick[0].label1.set_visible(False)
                    if xmax - tick[-1].get_loc() < gap_wanted:
                        tick[-1].label1.set_visible(False)

        if prior_density:
            # only report KL divergence for topmost posterior,
            # but plot the priors if available for the other posteriors
            # if priors known to be identical, plot them only once.

            if "priors_identical" in kwargs:
                priors_identical = kwargs.get("priors_identical")
            else:
                priors_identical = False

            if "force_draw" in kwargs:
                force_draw = kwargs.get("force_draw")
            else:
                force_draw = [True for i in range(len(self.subset))]

            for i, posterior in enumerate(self.subset):

                force_draw_i = force_draw[i]

                if "prior_samples_fnames" in kwargs:
                    prior_samples_fname = kwargs.get("prior_samples_fnames")[i]
                else:
                    prior_samples_fname = "prior_samples_"+posterior.ID+".npy"
                self._add_prior_density(plotter, posterior,
                            ndraws, normalize,
                            KL_divergence = KL_divergence if i == 0 else False,
                            KL_base = KL_base,
                            bootstrap = bootstrap,
                            n_simulate = kwargs.get('n_simulate'),
                            force_draw = force_draw_i,
                            prior_samples_fname=prior_samples_fname)

                if (i==0 and priors_identical):
                    break

        if veneer:
            self._veneer_spines_ticks(plotter, **kwargs)
        if crosshairs:
            # only for topmost posterior
            self._add_crosshairs(plotter, self.params.names, self.subset_to_plot[0].truths)


        self.credible_intervals=OrderedDict()

        if credible_interval_1d_all_show:
            for r in range(len(self.subset_to_plot)):

                id=self.get_attr("parent_ID")[r]+"_"+self.get_attr("ID")[r]
                self.r=r
                self.val_cred = []
                self.run = self.subset[0].subset_to_plot[r]

                self._add_credible_interval(plotter,
                                            self.subset[0],
                                            bootstrap=bootstrap,
                                            n_simulate=kwargs.get('n_simulate'),
                                            annotate=annotate_credible_interval,
                                            annotate_xy=annotate_xy,
                                            sixtyeight=sixtyeight,
                                            ninety=ninety,
                                            compute_all_intervals=compute_all_intervals)

                self.credible_intervals[id]=self.val_cred
        else:
                id=self.get_attr("parent_ID")[0]+"_"+self.get_attr("ID")[0]
                self.r=0
                self.val_cred = []
                self.run = self.subset[0].subset_to_plot[0]

                self._add_credible_interval(plotter,
                                            self.subset[0],
                                            bootstrap=bootstrap,
                                            n_simulate=kwargs.get('n_simulate'),
                                            annotate=annotate_credible_interval,
                                            annotate_xy=annotate_xy,
                                            sixtyeight=sixtyeight,
                                            ninety=ninety,
                                            compute_all_intervals=compute_all_intervals)

                self.credible_intervals[id]=self.val_cred


        self._plotter = plotter
        return plotter

    @make_verbose('Adding 1D marginal prior density functions',
                  'Added 1D marginal prior density functions')
    def _add_prior_density(self, plotter, posterior,
                           ndraws, normalize,
                           KL_divergence, KL_base,
                           bootstrap, n_simulate,
                           force_draw,
                           prior_samples_fname):
        """ Crudely estimate the prior density.

        Kullback-Leibler divergence estimated in bits for a combined run or
        the same run for which the credible intervals are calculated.

        """
        run = posterior.subset_to_plot[0]


        #self.samples[posterior.ID]=samples

        yield 'Plotting prior for posterior %s...' % posterior.ID

        l = posterior.likelihood

        if l is None:
            return # quietly do nothing
        elif not hasattr(l, 'prior'):
            return
        elif not hasattr(l.prior, 'draw'):
            return
        elif not callable(l.prior.draw):
            return

        if force_draw:
            samples, _ = l.prior.draw(ndraws, transform=True)
        else:
            samples_npy = prior_samples_fname
            try:
                samples = _np.load(samples_npy)
                print("Not drawing samples from the joint prior. Reading them instead from a pre-computed table:",prior_samples_fname)
            except:
                 samples, _ = l.prior.draw(ndraws, transform=True)
                 _np.save(samples_npy,samples)

        color, lw = (run.contours[key] for key in ('color', 'lw'))

        quantiles = [None] * 3

        with verbose(KL_divergence,
                     'Estimating 1D marginal KL-divergences in %s' % KL_base,
                     'Estimated 1D marginal KL-divergences') as condition:
            for i, ax in enumerate([plotter.subplots[i,i] \
                                for i in range(plotter.subplots.shape[0])]):

                name = self.params.names[i]
                bounds = {name: posterior.bounds[name]}
                settings = {'fine_bins': 1024,
                            'smooth_scale_1D': 0.3,
                            'boundary_correction_order': 1,
                            'mult_bias_correction_order': 1} # adopt from posterior settings or take custom input?

                idx = l.index(name)
                if idx is None: idx = l.prior.index(name)

                bcknd = MCSamples(sampler='uncorrelated',
                                  samples=samples[:,idx],
                                  weights=None,
                                  names=[name],
                                  ranges=bounds,
                                  settings=settings)

                if normalize:
                    bcknd.get1DDensity(name).normalize(by='integral',
                                                       in_place=True)

                x = _np.linspace(ax.xaxis.get_view_interval()[0],
                                 ax.xaxis.get_view_interval()[1],
                                 1000)

                ax.plot(x, bcknd.get1DDensity(name).Prob(x),
                        ls='-.', color=color, lw=lw)

                if not condition: continue # go to next iteration if no KL

                # a prototype Kullback-Leibler divergence callback
                # information in bits
                def KL(ns_run, logw):
                    x = ns_run['theta'][:,posterior.get_index(name)]
                    w_rel = _np.exp(logw - logw.max())
                    where = w_rel > run.kde_settings.get('min_weight_ratio',
                                                         1.0e-30)
                    prior = bcknd.get1DDensity(name).Prob(x[where])
                    p = getdist_kde(x[where], x, w_rel,
                                        ranges=[posterior.bounds[name]],
                                        idx=0,
                                        normalize=normalize,
                                        settings=run.kde_settings)
                    # Due to spline interpolation, very small densities can be
                    # negative, so manually give a small postive value which
                    # does not affect KL integral approximation
                    p[p<=0.0] = p[p>0.0].min()

                    KL = _np.sum(w_rel[where] \
                                   * (_np.log(p) - _np.log(prior))) \
                                   /_np.sum(w_rel[where])

                    if KL_base == 'bits':
                        return KL / _np.log(2.0)
                    elif KL_base == 'nats':
                        return KL
                    else:
                        raise ValueError('Invalid base for KL-divergence.')

                if bootstrap:
                    for j, cred_int in enumerate([0.025, 0.5, 0.975]):
                        quantiles[j] = run_ci_bootstrap(run.nestcheck_backend,
                                                     estimator_list=[KL],
                                                     cred_int=cred_int,
                                                     n_simulate=n_simulate,
                                                     simulate_weights=True,
                                                     flip_skew=True)
                    # KL in bits
                    interval = r'$D_{\mathrm{KL}}=%.2f_{-%.2f}^{+%.2f}$' \
                                                  % (quantiles[1],
                                                     quantiles[1] - quantiles[0],
                                                     quantiles[2] - quantiles[1])

                    yield ('%s KL-divergence = %.4f/-%.4f/+%.4f'
                            % (name,
                               quantiles[1],
                               quantiles[1] - quantiles[0],
                               quantiles[2] - quantiles[1]))

                    if not rcParams['text.usetex']:
                        fontsize = plotter.settings.lab_fontsize - 1
                    else:
                        fontsize = plotter.settings.lab_fontsize

                    ax.set_title(interval, color=color,
                                 fontsize=fontsize)
                else:
                    where = run.samples[:,0] > 0.0

                    ns_run = {'theta': run.samples[where,2:]}
                    divergence = KL(ns_run, _np.log(run.samples[where,0]))

                    yield ('%s KL-divergence = %.4f' % (name, divergence))

                    divergence = (r'$D_{\mathrm{KL}}=%.2f$' % divergence)

                    if not rcParams['text.usetex']:
                        fontsize = plotter.settings.lab_fontsize - 1
                    else:
                        fontsize = plotter.settings.lab_fontsize

                    ax.set_title(divergence, color=color,
                                 fontsize=fontsize)

        yield None

    @make_verbose('Adding 1D marginal credible intervals',
                  'Added 1D marginal credible intervals')
    def _add_credible_interval(self, plotter, posterior, bootstrap, n_simulate,
                               annotate, annotate_xy, sixtyeight,
                               ninety, compute_all_intervals):
        """
        Estimate 1-:math:`\sigma` credible interval in one-dimension on a
        combined run, or if such a run does not exist, on the run with
        the specified ID.

        Calls :func:`nestcheck.estimators.param_cred` for one-tailed weighted
        estimate; two such estimates give a credible interval which is
        symmetric in quantiles with respect to the median. Also calls
        :func:`nestcheck.error_analysis.run_ci_bootstrap` for
        credible interval on quantiles.

        :param bool sixtyeight:
            Plot the 68% credible interval about the median in 1D plots? If
            ``False`` plots 95% credible interval about the median -- i.e.,
            symmetric quantiles about the median.
        """
        diag = [plotter.subplots[i,i] for i in range(plotter.subplots.shape[0])]

        run = self.run #posterior.subset_to_plot[0]

        yield 'Plotting credible regions for posterior %s...' % posterior.ID

        color = run.contours['color']

        # estimator requires closure to be changable
        def get_estimator(quantile, param_ind):
            def estimator(*args, **kwargs):
                return param_cred(*args,
                                  probability=quantile,
                                  param_ind=param_ind,
                                  **kwargs)
            return estimator

        quantiles = [0.159, 0.5, 0.841] if sixtyeight else ([0.05,0.5,0.95] if ninety else [0.025, 0.5, 0.975])

        def format_CI(name, cred, summary, additional=2, sscript=False):

            if len(cred.shape) > 1:
                _qs = (cred[1,1],
                       cred[1,1] - cred[0,1],
                       cred[2,1] - cred[1,1])

            else:
                _qs = (cred[1],
                       cred[1] - cred[0],
                       cred[2] - cred[1])

            _p = max(_precision(_qs[0]), _precision(_qs[1]), _precision(_qs[2]))
            _f = '%.' + str(_p + additional) + 'f'

            if name: name += ' '

            stats = ('%s' % name) + ('CI$_{%i\%%} = ' % summary)

            if sscript:
                stats += (('%s_{-%s}^{+%s}$' % (_f, _f, _f)) % (_qs[0], _qs[1], _qs[2]))
            else:
                stats += (('%s/-%s/+%s$' % (_f, _f, _f)) % (_qs[0], _qs[1], _qs[2]))

            self.val_cred.append([np_.float(_f % _qs[0]),-np_.float(_f % _qs[1]),+np_.float(_f % _qs[2])])

            return stats

        if bootstrap:
            for i, ax in enumerate(diag):
                ind = posterior.get_index(self.params.names[i])
                def calculate_intervals(quantiles):
                    cred = _np.zeros((len(quantiles), len(quantiles)), dtype=_np.double)
                    for j, p in enumerate(quantiles):
                        for k, q in enumerate(quantiles):
                            cred[j,k] = run_ci_bootstrap(run.nestcheck_backend,
                                             estimator_list=[get_estimator(p, ind)],
                                             cred_int=q,
                                             n_simulate=n_simulate,
                                             simulate_weights=True,
                                             flip_skew=True)[0]
                    return cred

                cred = calculate_intervals(quantiles)

                zorder = max([_.zorder for _ in ax.get_children()]) + 1

                ax.axvspan(cred[0,0], cred[0,2], alpha=0.5,
                           facecolor=color,
                           edgecolor=color,
                           linewidth=0,
                           rasterized=True,
                           zorder=zorder)

                ax.axvspan(cred[2,0], cred[2,2], alpha=0.5,
                           facecolor=color,
                           edgecolor=color,
                           linewidth=0,
                           rasterized=True,
                           zorder=zorder)

                ax.axvspan(cred[0,2], cred[2,0], alpha=0.25,
                           facecolor=color,
                           edgecolor=color,
                           linewidth=0,
                           rasterized=True,
                           zorder=zorder)

                if annotate:
                    stats = format_CI('', # parameter name not needed on plot
                                      cred,
                                      68 if sixtyeight else (90 if ninety else 95),
                                      additional=1,
                                      sscript=True)

                    title = ax.get_title()
                    if title:
                        title = stats.center(30) + '\n' + title.center(30)
                    else:
                        title = stats

                    if not rcParams['text.usetex']:
                        fontsize = plotter.settings.lab_fontsize - 1
                    else:
                        fontsize = plotter.settings.lab_fontsize


                    if self.credible_interval_1d_all_show:
                        x_0,_=ax.get_xlim()
                        y_0,y_1=ax.get_ylim()
                        y_pos = y_1+0.11*(self.r+1)*(y_1-y_0)
                        ax.text(x_0, y_pos, title,
                                    color = color,
                                    #ha="center",
                                    alpha =1.,
                                    fontsize=fontsize)

                    else:
                        ax.set_title(title, color=color,
                                     fontsize=fontsize)


                if compute_all_intervals:
                    yield format_CI(self.params.names[i],
                                    cred,
                                    68 if sixtyeight else (90 if ninety else 95))

                    if sixtyeight:
                        yield format_CI(self.params.names[i],
                                        calculate_intervals([0.05, 0.5, 0.95]),
                                        90)
                        yield format_CI(self.params.names[i],
                                        calculate_intervals([0.025, 0.5, 0.975]),
                                        95)
                    elif ninety:
                        yield format_CI(self.params.names[i],
                                        calculate_intervals([0.159, 0.5, 0.841]),
                                        68)
                        yield format_CI(self.params.names[i],
                                        calculate_intervals([0.025, 0.5, 0.975]),
                                        95)
                    else:
                        yield format_CI(self.params.names[i],
                                        calculate_intervals([0.159, 0.5, 0.841]),
                                        68)
                        yield format_CI(self.params.names[i],
                                        calculate_intervals([0.05, 0.5, 0.95]),
                                        90)

        else:
            for i, ax in enumerate(diag):
                ind = posterior.get_index(self.params.names[i])
                def calculate_intervals(quantiles):
                    cred = _np.zeros(len(quantiles), dtype=_np.double)
                    for j, p in enumerate(quantiles):
                        where = run.samples[:,0] > 0.0
                        _t1 = run.samples[where,2:]
                        _t2 = _np.log(run.samples[where,0])
                        cred[j] = get_estimator(p, ind)({'theta': _t1}, _t2)

                    return cred

                cred = calculate_intervals(quantiles)
                zorder = max([_.zorder for _ in ax.get_children()]) + 1

                if self.r==0:
                    ax.axvspan(cred[0], cred[2], alpha=0.25,
                               facecolor=color,
                               edgecolor=color,
                               linewidth=0,
                               rasterized=True,
                               zorder=zorder)


                if annotate:
                    stats = format_CI('', # parameter name not needed on plot
                                      cred,
                                      68 if sixtyeight else (90 if ninety else 95),
                                      additional=1,
                                      sscript=True)

                    title = ax.get_title()
                    if title:
                        title = stats.center(30) + '\n' + title.center(30)
                    else:
                        title = stats

                    if not rcParams['text.usetex']:
                        fontsize = plotter.settings.lab_fontsize - 1
                    else:
                        fontsize = plotter.settings.lab_fontsize

                    if self.credible_interval_1d_all_show:
                        x_0,_=ax.get_xlim()
                        y_0,y_1=ax.get_ylim()
                        y_pos = y_1+0.11*(self.r+1)*(y_1-y_0)
                        ax.text(x_0, y_pos, title,
                                    color = color,
                                    #ha="center",
                                    alpha =1.,
                                    fontsize=fontsize)

                    else:
                        ax.set_title(title, color=color,
                                     fontsize=fontsize)


                if compute_all_intervals:
                    yield format_CI(self.params.names[i],
                                    cred,
                                    68 if sixtyeight else (90 if ninety else 95))

                    if sixtyeight:
                        yield format_CI(self.params.names[i],
                                        calculate_intervals([0.05, 0.5, 0.95]),
                                        90)
                        yield format_CI(self.params.names[i],
                                        calculate_intervals([0.025, 0.5, 0.975]),
                                        95)
                    elif ninety:
                        yield format_CI(self.params.names[i],
                                        calculate_intervals([0.159, 0.5, 0.841]),
                                        68)
                        yield format_CI(self.params.names[i],
                                        calculate_intervals([0.025, 0.5, 0.975]),
                                        95)
                    else:
                        yield format_CI(self.params.names[i],
                                        calculate_intervals([0.159, 0.5, 0.841]),
                                        68)
                        yield format_CI(self.params.names[i],
                                        calculate_intervals([0.05, 0.5, 0.95]),
                                        90)


        self.val_cred=np_.stack(self.val_cred,axis=0)
        yield None

    @staticmethod
    @make_verbose('Adding parameter truth crosshairs',
                  'Added crosshairs')
    def _add_crosshairs(plotter, names, truths):
        """ Add parameter crosshairs to triangle plot. """
        spine = next(plotter.subplots[0,0].spines.itervalues())
        lw = spine.get_linewidth()
        for i, name in enumerate(names):
            true_val = truths[name]
            if true_val is not None:
                for ax in plotter.subplots[:,i]:
                    if ax is not None:
                        ax.axvline(true_val, color='black', ls='-', lw=lw)
                if i > 0:
                    for ax in plotter.subplots[i,:i]:
                        if ax is not None:
                            ax.axhline(true_val, color='black', ls='-', lw=lw)

    @staticmethod
    @make_verbose('Veneering spines and axis ticks',
                  'Veneered')
    def _veneer_spines_ticks(plotter, lengthen=2.0, embolden=2.0,
                             **kwargs):
        """ Embolden spines, and embolden and lengthen ticks. """

        ax = plotter.subplots[0,0]

        major_length = ax.xaxis.majorTicks[0].tick1line.get_markersize()
        major_length *= lengthen
        minor_length = ax.xaxis.minorTicks[0].tick1line.get_markersize()
        minor_length *= lengthen
        lw = ax.xaxis.majorTicks[0].tick1line.get_markeredgewidth()
        lw *= embolden

        for i in range(plotter.subplots.shape[0]):
            for j in range(i+1):
                ax = plotter.subplots[i,j]
                ax.tick_params(which='major', colors='black', length=major_length)
                ax.tick_params(which='minor', colors='black', length=minor_length)
                ax.xaxis.set_tick_params(which='both', width=lw)
                ax.yaxis.set_tick_params(which='both', width=lw)
                for spine in ax.spines.itervalues():
                    spine.set_linewidth(lw)

    def _set_line_and_contour_args(self, lw=1.0, alpha=1.0, **kwargs):
        """ Match the :mod:`nestcheck` color scheme or let the user decide colors using kwargs.

        Always assigns reds (or the first user-defined color) to a combined run if it is found to exist.

        """
        # Color blind friendly
        nestcheck_colors = ['darkred',"darkblue", "darkorange", "darkgreen",'deeppink', 'maroon',
                            'indigo','dimgrey', 'olive']
        if 'line_colors' in kwargs:
            nestcheck_colors = kwargs.get("line_colors")


        for run, color in zip(self.subset_to_plot,
                              nestcheck_colors[:len(self.subset_to_plot)]):
            run.lines =  {'lw': kwargs.get('lw_1d', lw),
                          'color': color,
                          'alpha': alpha}
            run.contours = {'lw': lw, 'color': color, 'alpha': alpha}

    def KL_divergence(self, base='bits', bootstrap=False,
                      quantiles=[0.025, 0.5, 0.975],
                      n_simulate=200, **kwargs):
        """ Kullback-Leibler divergence integral jointly for all parameters.

        E.g., if you want the interval about the median containing divergence
        of 90% of realisations, declare ``quantiles=[0.05,0.5,0.95]``.

        """
        if kwargs:
            self.set_subset(**kwargs)

        nestcheck_bcknds, runs = self._filter_nestcheck_compatible()

        def estimator(ns_run, logw):
            w_rel = _np.exp(logw - logw.max())
            KL = _np.sum(w_rel * ns_run['logl']) / _np.sum(w_rel)
            KL -= logsumexp(logw)

            if base == 'bits':
                return KL / _np.log(2.0)
            elif base == 'nats':
                return KL
            else:
                raise ValueError('Invalid base for KL-divergence.')

        if bootstrap:
            _quantiles = {}
            for bcknd, run in zip(nestcheck_bcknds, runs):
                _quantiles[run.prepend_ID] = [run_ci_bootstrap(bcknd,
                                      estimator_list=[estimator],
                                      cred_int=q,
                                      n_simulate=n_simulate,
                                      simulate_weights=True,
                                      flip_skew=True)[0] for q in quantiles]
            return _quantiles
        else:
            divergence = {}
            for bcknd in nestcheck_bcknds:
                divergence[run.prepend_ID] = estimator(bcknd, get_logw(bcknd))

            return divergence

    def evidence_error(self, quantiles=[0.025,0.5,0.975], n_simulate=200,
                       simulate_weights=True, flip_skew=True, **kwargs):

        """ Estimate evidence error for nestcheck-compatible runs.

        E.g., if you want the interval about the median containing the evidence
        of 90% of realisations, declare ``quantiles=[0.05,0.5,0.95]``.

        """
        if kwargs:
            self.set_subset(**kwargs)

        nestcheck_bcknds, runs = self._filter_nestcheck_compatible()

        _quantiles = {}
        for bcknd, run in zip(nestcheck_bcknds, runs):
            _quantiles[run.prepend_ID] = [run_ci_bootstrap(bcknd,
                                  estimator_list=[logz],
                                  cred_int=q,
                                  n_simulate=n_simulate,
                                  simulate_weights=simulate_weights,
                                  flip_skew=flip_skew)[0] for q in quantiles]
        return _quantiles
