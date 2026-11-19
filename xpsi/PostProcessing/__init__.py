""" Module for various (posterior) sample post-processing operations.

Users
-----

The typical usage pattern is to use the following classes:

    * :class:`~xpsi.PostProcessing.Runs` to load a set of nested sampling
      runs on one underlying posterior distribution;
    * :class:`~xpsi.PostProcessing.CornerPlotter` to plot 1D and 2D
      marginal posterior density distributions (from a set of
      :class:`~xpsi.PostProcessing.Runs` instances) in lower-triangular
      form, and to compute various posterior estimators;
    * :class:`~xpsi.PostProcessing.CorrelationMatrix` to plot the 
      posterior distributions correlation structures;
    * :class:`~xpsi.PostProcessing.SignalPlotter` to plot thematic
      X-ray signals (or, generally, derived quantities) in posterior-
      expected form or conditional-posterior form, together with ground
      truths (i.e., injected) and in component form when there are two
      or more disjoint surface hot regions;
    * to use the :class:`~xpsi.PostProcessing.SignalPlotter` class, a
      user must supply a set of plot class instances for each posterior
      distribution to be processed, where these native classes are

        + :class:`~xpsi.PostProcessing.ResidualPlot`
        + :class:`~xpsi.PostProcessing.Residual1DPlot`
        + :class:`~xpsi.PostProcessing.ResidualStructurePlot`
        + :class:`~xpsi.PostProcessing.PulsePlot`
        + :class:`~xpsi.PostProcessing.SpectrumPlot`

      and can be configured both on a class and instance level.

Statistical residual diagnostics
--------------------------------

The :class:`~xpsi.PostProcessing.ResidualPlot` class can optionally add
(``plot_clusters=True``) a map of the clusters of large standardised
residuals, the residual distribution compared with a normal distribution
(with a Kolmogorov-Smirnov test), and the distribution of cluster sizes
compared with the range expected from counting noise alone (with a p-value
for the largest cluster). The :class:`~xpsi.PostProcessing.Residual1DPlot`
class shows the residual distribution and its Kolmogorov-Smirnov test on
its own.

The :class:`~xpsi.PostProcessing.ResidualStructurePlot` class quantifies
whether the standardised (data minus posterior-expected) residuals exhibit
2D structure over the joint phase-channel intervals, or are instead
consistent with random noise. It runs a suite of six complementary
statistical tests:

    1. a 2D power spectrum (FFT) test for periodic and smooth patterns;
    2. an empirical variogram test for distance-dependent correlation;
    3. a 2D autocorrelation-function scan over all lag vectors, sensitive
       to clustered residuals at any orientation (including diagonal);
    4. a Kolmogorov-Smirnov normality test combined with a Ljung-Box
       autocorrelation test (optional dependency ``statsmodels``);
    5. Moran's I global spatial autocorrelation (optional dependencies
       ``esda`` and ``libpysal``);
    6. a Wald-Wolfowitz runs test applied per channel and per phase bin.

A combined verdict (``STRUCTURED``, ``RANDOM``, or ``INCONCLUSIVE``) is
printed and stored on the instance, and its figure shows tests 1 to 4 and a
summary of all six tests. It is a separate plot class, which can be passed
to :class:`~xpsi.PostProcessing.SignalPlotter` together with
:class:`~xpsi.PostProcessing.ResidualPlot` for the same posterior, e.g.,
``plots = {'ST': [ResidualPlot(plot_clusters=True), ResidualStructurePlot()]}``;
the posterior-expected counts are then computed only once for both plots.

Customization
-------------

The aforementioned native classes can of course be subclassed if the
aim is to wrap, modify, or add functionality too.

In order to make entirely new plot type, we need to subclass the abstract
base class :class:`~xpsi.PostProcessing._signalplot.SignalPlot`, using the
native subclasses for guidance.

If new plot types are developed, please contribute to the source code via
a pull request. Isolate the new functionality in a proposed submodule,
adding your authorship information, and supply examples of usage in a Jupyter
notebook (either a new tutorial or by extending an existing tutorial as
appropriate). The subclass can then be imported from the submodule namespace
as shown below for the native classes. Then add the submodule to the Sphinx
doc pages. If you can include an example image in the subclass docstring,
please do, as it will be rendered directly in the class documentation.
Finally, if you can point to a journal article that implements the plot type,
please do.

"""
try:
    from xpsi.PostProcessing._global_imports import * # access, e.g., random_seed global
except ImportError:
    from xpsi import _warning
    _warning('Cannot use PostProcessing module.')
    raise
    __all__ = []
else:
    __all__ = ["Runs",
               "SignalPlotter",
               "PulsePlot",
               "SpectrumPlot",
               "ResidualPlot",
               "Residual1DPlot",
               "ResidualStructurePlot"]

    from ._runs import Runs
    from ._signalplotter import SignalPlotter
    from ._residual import ResidualPlot
    from ._1d_residual import Residual1DPlot
    from ._residual_structure import ResidualStructurePlot
    from ._pulse import PulsePlot
    from ._spectrum import SpectrumPlot
    from ._backends import NestedBackend

    def _precision(x):
        """ A solution adapted from Stack Overflow.

        Reference: questions/3018758/determine-precision-and-scale-of-particular-number-in-python

        """
        import math

        max_digits = 14
        int_part = int(abs(x))
        magnitude = 1 if int_part == 0 else int(math.log10(int_part)) + 1
        if magnitude >= max_digits:
            return 0
        frac_part = abs(x) - int_part
        precision = 1 if frac_part == 0.0 else -int(math.log10(frac_part)) + 1
        if precision >= 5:
            if magnitude + 5 >= max_digits:
                return max_digits - magnitude
            else:
                return 5

        return precision

    try:
        from ._corner import CornerPlotter
    except ImportError:
        _warning("Cannot use CornerPlotter module.")
    else:
        __all__.append("CornerPlotter")

    try:
        from ._correlation_matrix import CorrelationMatrix
    except ImportError:
        _warning("Cannot use CorrelationMatrix module.")
    else:
        __all__.append("CorrelationMatrix")

def set_random_seed(value):
    """ Set the NumPy random seed for the post-processing module. """
    from . import _global_imports
    try:
        _global_imports.random_seed = value
    except AttributeError:
        _warning('Could not set random seed for post-processing.')
