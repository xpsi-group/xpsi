""" Likelihood function implementations.

Overview
--------

Users and developers can implement the final phase of likelihood function
evaluation as an extension module.

These functions must operate on a channel-phase resolved set of signals
(one per hot region). Each component will have an associated phase shift to
account for. An exposure time must be invoked to scale the signal to yield the
expected count numbers.

A phase-resolved signal means that the signal is evaluated at a sequence
of phase points, and thus the signal is to be integrated over a sequence of
phase intervals if the data (event) space in each channel is discrete. If
one desires an unbinned likelihood function, this can be approximated with
a high number of phase bins using existing extensions, or a new extension
module can be straighforwardly developed to loop over events at greater
computational expense.

Lastly, a background model is required. If a physical background model is
condition upon, an extension can operate with the expected background count
(rate or number) signal. In lieu of such a model, one can invoke the default
background marginalisation protocol; channel-by-channel
background count rate variables will be numerically and rapidly marginalised
over. See the discussion in :ref:`R19` and the X-PSI technical notes to evaluate
whether this is appropriate for coupling with a given model for target source
signal.

"""

from .default_background_marginalisation import precomputation
from .default_background_marginalisation import eval_marginal_likelihood

from ._poisson_likelihood_given_background import poisson_likelihood_given_background
