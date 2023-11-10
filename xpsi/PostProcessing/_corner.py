from ._global_imports import *
from ._postprocessor import PostProcessor

from ._nestcheck_modifications import getdist_kde

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
            In addition, keyword arguments for avoiding unnecessary re-drawing of prior samples
            (``force_draw``, ``prior_samples_fnames`` and ``priors_identical``).
            Param``precisions`` (a list of integers or Nones) can be used to define the decimal number
            precision for each credible interval plotted. In case of 2 parameters, one can do e.g.
            precisions=[2,None] to use 2 digit decimal precision for the first parameter and use the
            default automatic precision for the second.


        """
        return None
