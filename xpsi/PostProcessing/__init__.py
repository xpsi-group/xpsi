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
    * :class:`~xpsi.PostProcessing.SignalPlotter` to plot thematic
      X-ray signals (or, generally, derived quantities) in posterior-
      expected form or conditional-posterior form, together with ground
      truths (i.e., injected) and in component form when there are two
      or more disjoint surface hot regions;
    * to use the :class:`~xpsi.PostProcessing.SignalPlotter` class, a
      user must supply a set of plot class instances for each posterior
      distribution to be processed, where these native classes are

        + :class:`~xpsi.PostProcessing.ResidualPlot`
        + :class:`~xpsi.PostProcessing.PulsePlot`
        + :class:`~xpsi.PostProcessing.SpectrumPlot`

      and can be configured both on a class and instance level.

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
    from ._global_imports import * # access, e.g., random_seed global
except ImportError:
    from .. import _warning
    _warning('Cannot use PostProcessing module.')
    raise
    __all__ = []
else:
    __all__ = ["Runs",
               "SignalPlotter",
               "PulsePlot",
               "SpectrumPlot",
               "ResidualPlot"]

    from ._runs import Runs
    from ._signalplotter import SignalPlotter
    from ._residual import ResidualPlot
    from ._pulse import PulsePlot
    from ._spectrum import SpectrumPlot

    try:
        from ._corner import CornerPlotter
    except ImportError:
        _warning("Cannot use CornerPlotter module.")
    else:
        __all__.append("CornerPlotter")

def set_random_seed(value):
    """ Set the NumPy random seed for the post-processing module. """
    from . import _global_imports
    try:
        _global_imports.random_seed = value
    except AttributeError:
        _warning('Could not set random seed for post-processing.')

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
        precision = -int(math.log10(frac_part)) + 1

        return precision


