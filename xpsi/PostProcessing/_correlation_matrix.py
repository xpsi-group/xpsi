""" Correlation matrix plots for posterior samples.

This submodule provides the :class:`~.CorrelationMatrix` class, an
alternative to :class:`~xpsi.PostProcessing.CornerPlotter` for
high-dimensional parameter spaces, where lower-triangle corner plots
become too large to read. Instead of plotting all pairwise marginal
densities, a colour-coded matrix of (posterior-weighted) linear or rank
correlation coefficients between all selected parameters is rendered.

The typical usage pattern in a post-processing notebook is identical to
that of the corner plotter::

    matrix = xpsi.PostProcessing.CorrelationMatrix([runs])
    fig = matrix.plot(params=names, IDs=OrderedDict([('run ID', ['run',]),]))

The computed matrices are cached in the
:attr:`CorrelationMatrix.correlation_matrices` attribute for further
programmatic use.

"""
from ._global_imports import *

from ._backends import NestedBackend
from ._postprocessor import PostProcessor

class CorrelationMatrix(PostProcessor):
    """ Plot a correlation matrix of posterior parameter samples.

    An alternative to :class:`~xpsi.PostProcessing.CornerPlotter` for
    complex models: the posterior sample correlation structure of a
    high-dimensional parameter space is summarised in a single
    colour-coded matrix, with the colour scale showing correlation
    (positive) or anti-correlation (negative) between parameter pairs.

    Correlation coefficients are computed from the *weighted* nested
    samples (via the :mod:`getdist` backend sample weights), so the
    matrix faithfully represents the posterior distribution rather than
    the raw distribution of nested sampling points.

    """

    @fix_random_seed
    @make_verbose('Executing posterior correlation estimation',
                  'Posterior correlation estimation complete')
    def plot(self,
             params,
             IDs=None,
             combine=False,
             combine_all=False,
             only_combined=False,
             force_combine=True,
             overwrite_combined=False,
             method='pearson',
             comparison='split',
             triangle='lower',
             cluster=False,
             annotate=True,
             annotate_fontsize=None,
             precision=2,
             threshold=None,
             cmap='RdBu_r',
             cell_size=0.5,
             tick_fontsize=None,
             label_fontsize=None,
             colorbar=True,
             colorbar_label=None,
             grid=True,
             write=False,
             root_filename='',
             directory='./',
             ext='.pdf',
             dpi=300,
             **kwargs):
        r""" Generate a posterior correlation matrix plot.

        Correlations are estimated from the sample backends of the runs
        selected via the ``IDs`` argument, in a manner analogous to
        :meth:`~xpsi.PostProcessing.CornerPlotter.plot`.

        :param list params:
            List of parameter strings for plotting. Must be shared by all
            posteriors selected with the ``IDs`` argument.

        :param OrderedDict IDs:
            Keys must be string identifiers of :class:`Runs` instances.
            Each dictionary element must be a list of string identifiers,
            each matching objects collected in the :class:`Runs` instance
            corresponding to the key. Defaults to ``None``, meaning attempt
            to use as many runs as possible.

        :param bool combine:
            Additionally combine the runs into a single run? (See
            :meth:`~xpsi.PostProcessing.CornerPlotter.plot`.)

        :param bool combine_all:
            Combine all runs in each :class:`Runs` instance or only those
            for which IDs are provided? Ignored if ``combine`` is ``False``.

        :param bool only_combined:
            Only use the combined run? Only heeded if a single posterior
            is selected, and ignored if ``combine`` is ``False``.

        :param bool force_combine:
            Force recombination of eligible run sets, even if a combined
            run is already cached?

        :param bool overwrite_combined:
            Overwrite combined-sample files on disk with the same filename?

        :param str method:
            Correlation estimator. Options are ``'pearson'`` (weighted
            linear correlation; default) and ``'spearman'`` (weighted rank
            correlation, more robust to non-linear monotonic degeneracies).

        :param str comparison:
            Layout when two or more runs are selected. Options:

                * ``'split'`` (default): if *exactly two* runs are
                  selected, render a single matrix whose lower triangle
                  shows the first run and whose upper triangle shows the
                  second run;
                * ``'panels'``: render one matrix panel per run,
                  side-by-side with a shared colour bar.

            If more than two runs are selected, ``'panels'`` is used
            regardless.

        :param str triangle:
            For a single-run plot, render only the ``'lower'`` triangle
            (default; the matrix is symmetric so the upper triangle and
            unit diagonal are redundant), or the ``'full'`` matrix.

        :param bool cluster:
            Reorder the parameters by hierarchical clustering on the
            correlation distance :math:`1 - |r|` (computed for the first
            run), so that blocks of mutually (anti-)correlated parameters
            appear adjacent. Requires :mod:`scipy`. This is often the most
            readable ordering in very high-dimensional spaces.

        :param bool annotate:
            Print the numeric correlation coefficient in each cell?
            For very large matrices annotation can be disabled, or the
            font size adjusted via ``annotate_fontsize``.

        :param float annotate_fontsize:
            Font size for cell annotations. Defaults to an automatic
            size based on the matrix dimension.

        :param int precision:
            Number of decimal digits for cell annotations.

        :param float threshold:
            If not ``None``, cells with :math:`|r| <` ``threshold`` are
            masked (rendered blank and not annotated), visually isolating
            the significant (anti-)correlations.

        :param str cmap:
            A diverging :mod:`matplotlib` colormap name. The colour scale
            is always symmetric about zero on :math:`[-1, 1]`.

        :param float cell_size:
            Approximate size (inches) of each matrix cell, from which the
            figure size is computed.

        :param float tick_fontsize:
            Font size for the parameter tick labels. Defaults to an
            automatic size based on the matrix dimension.

        :param float label_fontsize:
            Font size for panel titles and the colour bar label.

        :param bool colorbar:
            Draw a colour bar?

        :param str colorbar_label:
            Colour bar label. Defaults to an automatic label reporting
            the correlation method.

        :param bool grid:
            Draw thin white separators between cells?

        :param bool write:
            Export the figure?

        :param str root_filename:
            Root filename to prepend to the automatically generated name.
            Can be, e.g., a model and/or data set identifier.

        :param str directory:
            Directory to write to. Defaults to the current directory.

        :param str ext:
            File extension for writing, e.g. ``'.png'``.

        :param int dpi:
            Dots-per-square-inch setting for exporting plots.

        :param kwargs:
            Additional keyword arguments passed to
            :func:`matplotlib.axes.Axes.imshow`.

        :returns: The :class:`matplotlib.figure.Figure` instance.

        .. note::

            The computed matrices and parameter ordering are stored in
            the :attr:`correlation_matrices` attribute, an
            :class:`collections.OrderedDict` with run IDs as keys and
            dictionaries holding the parameter ``names``, ``labels``,
            and the ``matrix`` itself as values.

        """
        self.set_subset(IDs, combine, combine_all,
                        force_combine, only_combined,
                        overwrite_combined)
        self.set_params(params)

        try:
            for run in self.subset_to_plot:
                if not isinstance(run, NestedBackend):
                    raise TypeError('Nested sampling backends are required.')
        except AttributeError:
            print('Nested sampling runs are required.')
            raise

        if method not in ('pearson', 'spearman'):
            raise ValueError("Correlation method must be 'pearson' or "
                             "'spearman'.")

        runs = self.subset_to_plot

        if len(self._subset) > 1:
            run_labels = self.get_attr('parent_ID')
        else:
            run_labels = self.get_attr('ID')

        names = list(self.params.names)
        labels = list(self.params.labels)

        # compute one matrix per run
        matrices = [self._compute_matrix(run, names, method) for run in runs]

        # optionally reorder parameters by hierarchical clustering
        if cluster:
            order = self._cluster_order(matrices[0])
            names = [names[i] for i in order]
            labels = [labels[i] for i in order]
            matrices = [m[_np.ix_(order, order)] for m in matrices]

        self.correlation_matrices = OrderedDict(
            (ID, {'names': names, 'labels': labels, 'matrix': m})
            for ID, m in zip(run_labels, matrices))

        if colorbar_label is None:
            colorbar_label = (r'Spearman rank correlation $r_s$'
                              if method == 'spearman' else
                              r'Pearson correlation $r$')

        fig = self._render(matrices, labels, run_labels,
                           comparison=comparison,
                           triangle=triangle,
                           annotate=annotate,
                           annotate_fontsize=annotate_fontsize,
                           precision=precision,
                           threshold=threshold,
                           cmap=cmap,
                           cell_size=cell_size,
                           tick_fontsize=tick_fontsize,
                           label_fontsize=label_fontsize,
                           colorbar=colorbar,
                           colorbar_label=colorbar_label,
                           grid=grid,
                           **kwargs)
        self._fig = fig

        if write:
            root_filename = (root_filename + '__' if root_filename else '') +\
                'correlationMatrix__runs_' + \
                '_'.join(str(ID).replace(' ', '') for
                         ID in self.get_attr('ID')) + '__'
            fig.savefig(_os.path.join(directory,
                                      root_filename + 'matrix' + ext),
                        dpi=dpi, bbox_inches='tight')

        return fig

    @property
    def fig(self):
        """ Get the last :class:`matplotlib.figure.Figure` instance. """
        return self._fig

    @staticmethod
    def _get_samples_and_weights(run, names):
        """ Extract weighted posterior samples from a run's getdist backend.

        :returns:
            A tuple ``(X, w)`` where ``X`` has shape ``(nsamples,
            nparams)``, with columns ordered as :obj:`names`, and ``w``
            is the sample weight vector.

        """
        gd = run.getdist_backend

        cols = []
        for name in names:
            index = gd.paramNames.numberOfName(name)
            if index < 0:
                raise KeyError('Parameter %s not found in the sample '
                               'backend of run %s.' % (name, run.ID))
            cols.append(gd.samples[:, index])

        return _np.column_stack(cols), _np.asarray(gd.weights, dtype=float)

    @staticmethod
    def _weighted_pearson(X, w):
        """ Weighted linear (Pearson) correlation matrix.

        :param ndarray[n,d] X: Sample matrix.
        :param ndarray[n] w: Sample weights.

        """
        w = w / _np.sum(w)
        mean = w @ X
        diffs = X - mean
        cov = (diffs * w[:, None]).T @ diffs
        sigma = _np.sqrt(_np.diag(cov))
        sigma[sigma == 0.0] = 1.0 # zero-variance guard; correlation -> 0
        corr = cov / _np.outer(sigma, sigma)
        # guard against round-off outside [-1,1]
        return _np.clip(corr, -1.0, 1.0)

    @staticmethod
    def _weighted_ranks(x, w):
        """ Weighted fractional ranks of a sample vector. """
        order = _np.argsort(x, kind='mergesort')
        ranks = _np.empty(x.shape[0], dtype=float)
        cumulative = _np.cumsum(w[order])
        ranks[order] = cumulative - 0.5 * w[order]
        return ranks

    def _compute_matrix(self, run, names, method):
        """ Compute the correlation matrix for a single run. """
        X, w = self._get_samples_and_weights(run, names)

        if method == 'spearman':
            X = _np.column_stack([self._weighted_ranks(X[:, i], w)
                                  for i in range(X.shape[1])])

        return self._weighted_pearson(X, w)

    @staticmethod
    def _cluster_order(matrix):
        """ Hierarchical-clustering ordering on correlation distance. """
        try:
            from scipy.cluster.hierarchy import linkage, leaves_list
            from scipy.spatial.distance import squareform
        except ImportError:
            _warning('Cannot import scipy; skipping parameter clustering.')
            return list(range(matrix.shape[0]))

        distance = 1.0 - _np.abs(matrix)
        _np.fill_diagonal(distance, 0.0)
        distance = 0.5 * (distance + distance.T) # enforce exact symmetry
        Z = linkage(squareform(distance, checks=False), method='average')
        return list(leaves_list(Z))

    def _render(self, matrices, labels, run_labels,
                comparison, triangle,
                annotate, annotate_fontsize, precision, threshold,
                cmap, cell_size, tick_fontsize, label_fontsize,
                colorbar, colorbar_label, grid, **kwargs):
        """ Render the matrix (or matrices) with matplotlib. """
        ndims = len(labels)

        split = (len(matrices) == 2 and comparison == 'split')

        if len(matrices) == 1 or split:
            n_panels = 1
        else:
            n_panels = len(matrices)

        if tick_fontsize is None:
            tick_fontsize = max(5.0, min(11.0, 200.0 / max(ndims, 1) / n_panels))
        if annotate_fontsize is None:
            annotate_fontsize = max(4.0, min(10.0, 160.0 / max(ndims, 1) / n_panels))
        if label_fontsize is None:
            label_fontsize = tick_fontsize + 2.0

        panel = max(3.0, cell_size * ndims)
        figsize = (panel * n_panels + (1.0 if colorbar else 0.0) + 1.2,
                   panel + 1.2)

        fig, axes = plt.subplots(1, n_panels, figsize=figsize, squeeze=False)
        axes = axes[0]
        self._axes = axes

        # assemble per-panel display matrices
        display = []
        titles = []
        if split:
            merged = _np.full((ndims, ndims), _np.nan)
            il = _np.tril_indices(ndims, k=-1)
            iu = _np.triu_indices(ndims, k=1)
            merged[il] = matrices[0][il]
            merged[iu] = matrices[1][iu]
            display.append(merged)
            titles.append('lower: %s  /  upper: %s' % tuple(run_labels[:2]))
        elif len(matrices) == 1:
            m = matrices[0].copy()
            if triangle == 'lower':
                m[_np.triu_indices(ndims, k=0)] = _np.nan
            display.append(m)
            titles.append(str(run_labels[0]))
        else:
            for m, ID in zip(matrices, run_labels):
                _m = m.copy()
                if triangle == 'lower':
                    _m[_np.triu_indices(ndims, k=0)] = _np.nan
                display.append(_m)
                titles.append(str(ID))

        if threshold is not None:
            for m in display:
                m[_np.abs(m) < threshold] = _np.nan

        colormap = _mpl.colormaps[cmap] if hasattr(_mpl, 'colormaps') \
                                        else cm.get_cmap(cmap)
        colormap = colormap.copy()
        colormap.set_bad(color=(0, 0, 0, 0)) # masked cells transparent

        tick_labels = [r'$%s$' % l for l in labels]

        image = None
        for ax, m, title in zip(axes, display, titles):
            image = ax.imshow(m, cmap=colormap, vmin=-1.0, vmax=1.0,
                              origin='upper', interpolation='nearest',
                              **kwargs)

            ax.set_xticks(range(ndims))
            ax.set_yticks(range(ndims))
            ax.set_xticklabels(tick_labels, rotation=90,
                               fontsize=tick_fontsize)
            ax.set_yticklabels(tick_labels, fontsize=tick_fontsize)
            ax.tick_params(length=0)

            if len(matrices) > 1 or split:
                ax.set_title(title, fontsize=label_fontsize)

            if grid:
                ax.set_xticks(_np.arange(ndims + 1) - 0.5, minor=True)
                ax.set_yticks(_np.arange(ndims + 1) - 0.5, minor=True)
                ax.grid(which='minor', color='white', linewidth=0.75)
                ax.tick_params(which='minor', length=0)

            for spine in ax.spines.values():
                spine.set_visible(False)

            if annotate:
                fmt = '%.{}f'.format(int(precision))
                for i in range(ndims):
                    for j in range(ndims):
                        value = m[i, j]
                        if _np.isnan(value):
                            continue
                        color = 'white' if _np.abs(value) > 0.6 else 'black'
                        ax.text(j, i, fmt % value,
                                ha='center', va='center',
                                fontsize=annotate_fontsize,
                                color=color)

        if colorbar and image is not None:
            cbar = fig.colorbar(image, ax=list(axes), fraction=0.046,
                                pad=0.02, shrink=0.85)
            cbar.set_label(colorbar_label, fontsize=label_fontsize)
            cbar.ax.tick_params(labelsize=tick_fontsize)
            cbar.outline.set_visible(False)

        return fig
