""" Tests for the background-marginalised likelihood. """

import faulthandler

import numpy as np

from xpsi.likelihoods.default_background_marginalisation import eval_marginal_likelihood
from xpsi.likelihoods.default_background_marginalisation import precomputation

# Maximum time (in seconds) allowed for one likelihood evaluation of a single channel
TIMEOUT = 60.0


def _loglikelihood(exposure_time, counts, signal, phase_shift, support):
    """ Evaluate the background-marginalised log-likelihood. """
    counts = np.array([counts], dtype=np.double)
    signal = np.array([signal], dtype=np.double)
    loglike = eval_marginal_likelihood(exposure_time,
                                       np.linspace(0.0, 1.0, counts.shape[1] + 1),
                                       counts,
                                       (signal,),
                                       (np.linspace(0.0, 1.0, signal.shape[1]),),
                                       np.array([phase_shift]),
                                       precomputation(counts.astype(np.int32)),
                                       np.array([support], dtype=np.double),
                                       1000, 0.0, 1.0e-8, 1.0e-3, 10.0, -1.0e90)[0]
    return loglike


class TestBackgroundMarginalisationIssue518(object):
    """ Regression test for issue #518.

    One energy channel of a synthetic ST data set with no background, at a parameter
    vector where the model signal is zero in some phase bins that contain counts.
    The Newton estimate of the best-fit background used to overshoot to a negative
    value and stop at zero, so that the integrand was normalised at the low end of the
    integration interval. With a tight background support of [0, 1e-11] counts/s, the
    integrand then rose to ~1e211, and GSL cquad never finished (its error estimate
    overflowed). With an unbounded support, the log-likelihood was +inf.

    The reference values were checked against an independent calculation (maximum
    found by bisection, integral computed with scipy.integrate.quad).
    """

    exposure_time = 1.0e6
    phase_shift = 0.51522312
    counts = [0, 0, 6, 131, 445, 1121, 1965, 2967, 3861, 4725, 5207, 5702, 5912, 6022,
              5651, 5293, 4959, 4516, 3855, 3407, 2900, 2377, 2005, 1601, 1195, 902, 599,
              334, 151, 58, 8, 0]
    signal = [0.1312685255677862, 0.1244573842611507, 0.11735800771358744,
              0.110104659697357, 0.10280996148537738, 0.09556550157583055,
              0.0884432158487088, 0.08149750256166845, 0.07476764890001181,
              0.06827980979061828, 0.06204913824314168, 0.05608198755342002,
              0.05037773907206468, 0.04493038960172567, 0.03972988623083532,
              0.03476313858800521, 0.030014885102893008, 0.025468159584118875,
              0.021104738892653817, 0.016903992233303595, 0.012850280618456302,
              0.009153761613574292, 0.0060369913011612125, 0.003588428269002173,
              0.0018201316624806847, 0.0006985786687345404, 0.00013955225824213515,
              2.942156216152547e-06, 1.3666828401156535e-36, 0.0, 0.0, 0.0, 0.0, 0.0,
              0.0, 0.0, 0.0, 0.0, 3.0117271162432395e-10, 5.655308922833262e-05,
              0.0009806827413526417, 0.0038040432854325266, 0.009142175813641658,
              0.017318757504263326, 0.028368022060331424, 0.041997897093318394,
              0.05742927955365583, 0.07312413392268027, 0.08819357752072775,
              0.10231526616772711, 0.11518864065189004, 0.12655365326319978,
              0.13620586164399945, 0.14400508855527985, 0.14988005012774142,
              0.15382729387100896, 0.15590563078819905, 0.15622702770198071,
              0.15494509787599256, 0.15224282805060138, 0.14832037071555773,
              0.14338387244889514, 0.13763594437424936, 0.1312685255677862]

    def _run_with_timeout(self, support, capfd):
        # The integration runs in compiled code that cannot be interrupted from
        # Python. If it hangs, the faulthandler watchdog thread prints the stack
        # and ends the test run, instead of blocking it (e.g. on CI) indefinitely.
        # Output capturing is disabled so that the stack is shown in the log.
        with capfd.disabled():
            faulthandler.dump_traceback_later(TIMEOUT, exit=True)
            try:
                return _loglikelihood(self.exposure_time, self.counts, self.signal,
                                      self.phase_shift, support)
            finally:
                faulthandler.cancel_dump_traceback_later()

    def test_tight_support(self, capfd):
        loglike = self._run_with_timeout([0.0, 1.0e-11], capfd)
        assert np.isfinite(loglike)
        assert np.isclose(loglike, -4521.763801536406, rtol=1.0e-8, atol=0.0)

    def test_unbounded_support(self, capfd):
        loglike = self._run_with_timeout([0.0, -1.0], capfd)
        assert np.isfinite(loglike)
        assert np.isclose(loglike, -2772.718471175409, rtol=1.0e-8, atol=0.0)
