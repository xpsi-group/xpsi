.. module:: xpsi.PostProcessing

.. _PostProcessing:

PostProcessing
==============

.. autosummary::

    xpsi.PostProcessing.Run
    xpsi.PostProcessing.NSBackend
    xpsi.PostProcessing.Runs
    xpsi.PostProcessing.PostProcessor

.. autoclass:: xpsi.PostProcessing.Run
    :members:
    :special-members: __call__

.. autoclass:: xpsi.PostProcessing.NSBackend
    :members:
    :special-members: __call__
    :show-inheritance:

.. autoclass:: xpsi.PostProcessing.Runs
    :members:
    :special-members: __call__

.. autoclass:: xpsi.PostProcessing.PostProcessor
    :members:
    :special-members: __call__
    :exclude-members: plot_posteriorDensity, plot_pulse_and_spectrum

    .. automethod:: plot_posteriorDensity(self, params, run_IDs=None, combine=False, combine_all=False, only_combined=False, bootstrap_estimators=True, bootstrap_density=False, separate_plots=False, write=False, root_filename='', directory='./', ext='.pdf', dpi=300, maxdots=2000,**kwargs)
    .. automethod:: plot_pulse_and_spectrum(self, params, run_IDs=None, combine=False, combine_all=False, only_combined=False, bootstrap_estimators=True, bootstrap_density=False, separate_plots=False, write=False, root_filename='', directory='./', ext='.pdf', dpi=300, maxdots=2000,**kwargs)
