from ._global_imports import *
from scipy import stats as _stats
from ._signalplot import SignalPlot
import textwrap as _textwrap

# Optional heavy dependency for Moran's I – graceful fallback if absent
try:
    from libpysal.weights import lat2W as _lat2W
    from esda.moran import Moran as _Moran
    _HAS_ESDA = True
except ImportError:
    _HAS_ESDA = False

# Optional statsmodels for Ljung-Box test
try:
    from statsmodels.stats.diagnostic import acorr_ljungbox as _acorr_ljungbox
    _HAS_STATSMODELS = True
except ImportError:
    _HAS_STATSMODELS = False


class ResidualStructurePlot(SignalPlot):
    """ Test whether the 2D residuals contain structure, and plot the tests.

    Computes the standardised residuals between the data and the
    posterior-expected count signal over joint channel-phase intervals,
    :math:`(c_{ik}-d_{ik})/\\sqrt{c_{ik}}`, and runs six complementary tests
    of whether they exhibit 2D structure or are consistent with random noise:

        1. **2D power spectrum** – periodic and smooth structure from the
           2D periodogram (calibrated with white-noise simulations).
        2. **Variogram** – how the correlation decays with pixel distance on
           the phase-channel grid (calibrated with white-noise simulations).
        3. **2D autocorrelation-function (ACF) scan** – correlation at every
           lag vector up to ``structure_acf2d_max_lag``, at any orientation.
        4. **KS + Ljung–Box** – Gaussianity of the residuals and serial
           autocorrelation of the flattened residuals (Ljung–Box requires
           ``statsmodels``; KS always runs).
        5. **Moran's I** – global spatial autocorrelation (requires ``esda``
           and ``libpysal``; skipped if absent).
        6. **Runs test** – Wald–Wolfowitz test applied per channel and per
           phase bin.

    A summary of the tests and a combined verdict (``STRUCTURED`` if at
    least three tests flag structure, ``RANDOM`` if none does,
    ``INCONCLUSIVE`` otherwise) is printed when the plot is finalised. The
    per-test results are stored in :attr:`residual_structure_results` and
    the verdict in :attr:`residual_structure_verdict`.

    The figure contains:

        * in the first two columns, one panel per test for tests 1 to 4: the
          2D power spectrum and the normalised variogram (top row), the 2D
          ACF map with the strongest lag marked and a normal Q-Q plot
          (bottom row);
        * in the third column, a summary of all six tests (including tests 5
          and 6, which have no panel): for each test, a short description
          of what it looks for and how to read its panel, and the values
          and threshold that decide its verdict; then the combined verdict.

    The residual map itself is shown by
    :class:`~xpsi.PostProcessing.ResidualPlot`.

    This plot is independent of :class:`~xpsi.PostProcessing.ResidualPlot`.
    Both use the posterior-expected counts, so they can be passed together
    for the same posterior, e.g.::

        plots = {'ST': [xpsi.ResidualPlot(plot_clusters=True),
                        xpsi.ResidualStructurePlot()]}

    They then use the same posterior samples, and with caching (the default
    of :meth:`~xpsi.PostProcessing.SignalPlotter.plot`) the model is
    evaluated only once for both.

    Each instance draws one figure: create a new instance for each call of
    :meth:`~xpsi.PostProcessing.SignalPlotter.plot` (reusing an instance that
    has already been plotted raises a ``RuntimeError``).

    :param list parameters_vector:
         If given, the residuals are computed for this single vector of model
         parameter values (in the order of the likelihood parameters) instead
         of for posterior samples. :class:`~xpsi.PostProcessing.SignalPlotter`
         only looks at the first plot passed for a posterior, so this then
         applies to all plots passed together.

    :param int structure_power_n_null:
         Number of white-noise simulations used to calibrate the null
         distribution of the 2D power spectrum test (test 1).  Default
         1000; the smallest reportable p-value is ``1/(n_null+1)``.

    :param float structure_power_alpha:
         Significance level of the 2D power spectrum test, before the
         internal Bonferroni factor of 2 over its two statistics.
         Default 0.05.

    :param float structure_power_low_k_fraction:
         Radial-frequency cut (as a fraction of the Nyquist radius)
         defining the low-k region for the smooth-structure statistic
         of the 2D power spectrum test.  Default 0.15.

    :param int structure_power_random_seed:
         Seed for the power spectrum null simulations (reproducibility).

    :param int structure_variogram_max_lag:
         Maximum pixel lag used for the variogram test (test 2).
         Default 20.

    :param int structure_variogram_n_null:
         Number of white-noise simulations used to calibrate the null
         distribution of the variogram test.  Default 1000.

    :param float structure_variogram_alpha:
         Significance level of the variogram test, before the internal
         Bonferroni factor of 2 over its two statistics.  Default 0.05.

    :param int structure_variogram_random_seed:
         Seed for the variogram null simulations (reproducibility).

    :param int structure_acf2d_max_lag:
         Maximum pixel lag along each axis for the 2D ACF scan (test 3).
         Default 8.

    :param kwargs:
         Keyword arguments for :class:`SignalPlot` (e.g. ``write``, ``dpi``,
         ``fig_dir``, ``root_filename``).

    """

    __figtype__ = 'signalplot_residual_structure'

    # do not change at runtime (see base class comments):
    __caching_targets__ = ['expected_counts']

    __rows__          = 3
    __columns__       = 2.5
    __ax_rows__       = 2
    # Three columns of panels separated by two empty spacer columns, so the
    # two gaps can be sized independently (wspace is then 0): the second
    # gap holds the y-axis labels of the middle column, which are on the
    # right of its panels.  The third column, for the summary, is wider.
    __ax_columns__    = 5
    __height_ratios__ = [1]*2
    # wider third column for the summary text
    __width_ratios__  = [1, 0.25, 1, 0.25, 1.6]
    __wspace__        = 0.0
    __hspace__        = 0.42

    @make_verbose('Instantiating a residual structure plotter for posterior checking',
                  'Residual structure plotter instantiated')
    def __init__(self,
                 parameters_vector=None,
                 structure_power_n_null=1000,
                 structure_power_alpha=0.05,
                 structure_power_low_k_fraction=0.15,
                 structure_power_random_seed=0,
                 structure_variogram_max_lag=20,
                 structure_variogram_n_null=1000,
                 structure_variogram_alpha=0.05,
                 structure_variogram_random_seed=0,
                 structure_acf2d_max_lag=8,
                 **kwargs):
        """ Constructor method for the residual structure tests and their plot.

        See the class docstring for the parameters.
        """
        super(ResidualStructurePlot, self).__init__(**kwargs)

        self._structure_power_n_null = structure_power_n_null
        self._structure_power_alpha = structure_power_alpha
        self._structure_power_low_k_fraction = structure_power_low_k_fraction
        self._structure_power_random_seed = structure_power_random_seed
        self._structure_variogram_max_lag = structure_variogram_max_lag
        self._structure_variogram_n_null = structure_variogram_n_null
        self._structure_variogram_alpha = structure_variogram_alpha
        self._structure_variogram_random_seed = structure_variogram_random_seed
        self._structure_acf2d_max_lag = structure_acf2d_max_lag
        if not parameters_vector is None:
            self.parameters_vector = parameters_vector

        self._get_figure()

        plt.close()

    @make_verbose('ResidualStructurePlot object iterating over samples',
                  'ResidualStructurePlot object finished iterating')
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
                'ResidualStructurePlot instance has already been plotted; create a new '
                'ResidualStructurePlot instance for each plotting call (e.g. instantiate the '
                'plot in the same notebook cell as the plotting call).')

    @make_verbose('ResidualStructurePlot object finalizing',
                  'ResidualStructurePlot object finalized')
    def finalize(self):
        self._check_not_reused()
        self._compute_residuals()
        self._run_structure_analysis()
        self._draw_diagnostics()
        self._finalized = True

    def _compute_residuals(self):
        """ Compute and store the standardised residuals
        :math:`(c_{ik}-d_{ik})/\\sqrt{c_{ik}}` in ``self._residuals``. """
        self._residuals = ((self.expected_counts - self._signal.data.counts)
                           / _np.sqrt(self.expected_counts))
        return self._residuals

    def _run_structure_analysis(self):
        """Run all 2-D structure tests on ``self._residuals`` and print a
        summary.  Results are stored in :attr:`residual_structure_results`
        and :attr:`residual_structure_verdict`.

        Called automatically from :meth:`finalize`.
        """
        if not hasattr(self, '_residuals'):
            self._compute_residuals()

        # NOTE: the order of this list defines the test numbering used in
        # the printed report AND in the figure: tests with a diagnostic
        # panel come first (in panel order), followed by the tests without
        # a panel.
        results = [
            # -- tests with a diagnostic panel --
            self._structure_test_power_spectrum(
                self._residuals,
                n_null=self._structure_power_n_null,
                alpha=self._structure_power_alpha,
                low_k_fraction=self._structure_power_low_k_fraction,
                random_seed=self._structure_power_random_seed),
            self._structure_test_variogram(
                self._residuals,
                max_lag=self._structure_variogram_max_lag,
                n_null=self._structure_variogram_n_null,
                alpha=self._structure_variogram_alpha,
                random_seed=self._structure_variogram_random_seed),
            self._structure_test_acf2d(self._residuals, self._structure_acf2d_max_lag),
            self._structure_test_normality_ljungbox(self._residuals),
            # -- tests without a diagnostic panel --
            self._structure_test_morans_i(self._residuals),
            self._structure_test_runs(self._residuals),
        ]

        # Aggregate verdict: >=3 of the definitive tests flag structure
        definitive = [t for t in results if t['structured'] is not None]
        n_structured = len([t for t in definitive if t['structured']])
        n_total = len(definitive)

        if n_structured >= 3:
            verdict = 'STRUCTURED'
        elif n_structured == 0:
            verdict = 'RANDOM'
        else:
            verdict = 'INCONCLUSIVE'

        self.residual_structure_results = results
        self.residual_structure_verdict = verdict

        # -- Pretty-print summary ------------------------------------------
        print('\n' + '=' * 65)
        print('  XPSI 2-D RESIDUAL STRUCTURE ANALYSIS')
        print('  Axes: phase [cycles] x channel')
        print('=' * 65)
        for i, t in enumerate(results, start=1):
            if t['structured'] is True:
                flag = 'x STRUCTURED'
            elif t['structured'] is False:
                flag = 'v random    '
            else:
                flag = '- skipped   '
            print(f"  [{flag}]  Test {i}: {t['test']}")
            print(f"                {t['interpretation']}")
        print('-' * 65)
        print(f"  VERDICT : {verdict}  "
              f"({n_structured}/{n_total} tests flagged structure)")
        print('=' * 65 + '\n')

    # -- Individual tests ---------------------------------------------------

    @staticmethod
    def _structure_test_power_spectrum(residuals, n_null=500, alpha=0.05,
                                       low_k_fraction=0.15, random_seed=0):
        """2-D power spectrum test, Monte-Carlo calibrated against white noise.

        Two complementary statistics are computed from the 2-D periodogram
        (DC component excluded):

        * ``g_max`` - Fisher's g statistic: maximum power / total power.
          Sensitive to PERIODIC structure, whose power concentrates in a
          few Fourier modes.
        * ``f_lowk`` - fraction of total power at low spatial frequencies
          (radial frequency below ``low_k_fraction`` of the Nyquist
          radius).  Sensitive to SMOOTH, blobby structure, whose power
          concentrates at low k.

        The null distributions of both statistics depend on the grid
        shape, so they are calibrated empirically by simulating
        ``n_null`` white-noise fields of the same shape and standard
        deviation.  One-sided empirical p-values (with an add-one
        correction, so the smallest reportable p is ``1/(n_null+1)``)
        are combined with a Bonferroni factor of 2: structure is flagged
        when ``min(p) < alpha / 2``.

        Unlike a fixed peak/median-ratio threshold, ``alpha`` is a
        genuine significance level (default 0.05), consistent with the
        other tests in the suite, and requires no per-dataset tuning.

        :param residuals: 2-D array of standardised residuals.
        :param int n_null: Number of white-noise simulations for the
            null calibration.  Default 500.
        :param float alpha: Significance level before the Bonferroni
            factor of 2.  Default 0.05.
        :param float low_k_fraction: Radial-frequency cut (as a fraction
            of the Nyquist radius) defining the "low-k" region for the
            smooth-structure statistic.  Default 0.15.
        :param int random_seed: Seed for the null simulations, for
            reproducibility.
        :returns: dict with keys ``test``, ``g_max``, ``g_null_q95``,
            ``p_periodic``, ``f_lowk``, ``f_null_q95``, ``p_smooth``,
            ``min_p``, ``alpha_bonferroni``, ``n_null``, ``power_map``,
            ``structured`` (bool), and ``interpretation`` (str).
        """
        rng = _np.random.default_rng(random_seed)
        ny, nx = residuals.shape

        def _spectrum_stats(field):
            p = _np.abs(_np.fft.fftshift(_np.fft.fft2(field))) ** 2
            p[ny // 2, nx // 2] = 0.0          # remove DC
            total = p.sum()
            ky = _np.fft.fftshift(_np.fft.fftfreq(ny))[:, None]
            kx = _np.fft.fftshift(_np.fft.fftfreq(nx))[None, :]
            kr = _np.sqrt(ky ** 2 + kx ** 2) / _np.sqrt(0.5)  # 0..1
            low = (kr <= low_k_fraction) & (kr > 0)
            return p.max() / total, p[low].sum() / total, p

        g_obs, f_obs, power_map = _spectrum_stats(residuals)

        g_null = _np.empty(n_null)
        f_null = _np.empty(n_null)
        sd = residuals.std()
        for i in range(n_null):
            g_null[i], f_null[i], _ = _spectrum_stats(
                rng.normal(0.0, sd, (ny, nx)))

        # one-sided empirical p-values (add-one correction avoids p = 0)
        p_g = (_np.count_nonzero(g_null >= g_obs) + 1.0) / (n_null + 1.0)
        p_f = (_np.count_nonzero(f_null >= f_obs) + 1.0) / (n_null + 1.0)

        min_p      = min(p_g, p_f)
        alpha_bonf = alpha / 2.0
        structured = bool(min_p < alpha_bonf)

        return {
            'test'             : '2D Power Spectrum (MC-calibrated)',
            'g_max'            : float(g_obs),
            'g_null_q95'       : float(_np.percentile(g_null, 95)),
            'p_periodic'       : float(p_g),
            'f_lowk'           : float(f_obs),
            'f_null_q95'       : float(_np.percentile(f_null, 95)),
            'p_smooth'         : float(p_f),
            'min_p'            : float(min_p),
            'alpha_bonferroni' : alpha_bonf,
            'n_null'           : n_null,
            'power_map'        : power_map,
            'structured'       : structured,
            'interpretation'   : (
                f'periodic: g = {g_obs:.4f} (p = {p_g:.4f});  '
                f'smooth: f_lowk = {f_obs:.3f} (p = {p_f:.4f});  '
                f'alpha/2 = {alpha_bonf:.3f}  '
                f'-> {"structured" if structured else "random"}'
            ),
        }

    @staticmethod
    def _structure_test_morans_i(residuals):
        """Global Moran's I spatial autocorrelation test.

        Moran's I quantifies how similar neighbouring cells are to one
        another.  Under the null hypothesis of complete spatial
        randomness, I ~ 0.  A significantly positive I indicates positive
        autocorrelation (nearby residuals tend to have the same sign).

        Queen contiguity weights (8-connected neighbourhood) are used so
        that row/column AND diagonal neighbours are all considered.

        :param residuals: 2-D array of standardised residuals.
        :returns: dict with keys ``test``, ``I``, ``p_value``,
            ``structured`` (bool or None if skipped), and
            ``interpretation`` (str).

        .. note::
            Requires the optional packages ``esda`` and ``libpysal``
            (``pip install esda libpysal``).  Skipped if unavailable.
        """
        if not _HAS_ESDA:
            return {
                'test': "Moran's I (spatial autocorrelation)",
                'structured': None,
                'interpretation': (
                    'skipped - install esda + libpysal to enable '
                    '(pip install esda libpysal)'
                ),
            }

        ny, nx = residuals.shape
        w = _lat2W(ny, nx, rook=False)  # queen contiguity (8-connected)
        w.transform = 'r'  # row-standardise
        mi = _Moran(residuals.ravel(), w)
        structured = bool(mi.p_sim < 0.05)

        return {
            'test': "Moran's I (spatial autocorrelation)",
            'I': float(mi.I),
            'p_value': float(mi.p_sim),
            'structured': structured,
            'interpretation': (
                f"I = {mi.I:.4f},  p = {mi.p_sim:.4f}  "
                f"-> {'structured' if structured else 'random'}"
            ),
        }

    @staticmethod
    def _structure_test_variogram(residuals, max_lag=20, n_null=500,
                                  alpha=0.05, random_seed=0):
        """Empirical isotropic variogram test, Monte-Carlo calibrated.

        For a random residual field the normalised variogram sits at the
        sill (``gamma(h)/sigma^2 ~ 1``) for every lag - a *pure nugget*
        model.  Spatial structure pulls it away from the sill: smooth
        short-range correlation depresses it at small lags, while
        periodic or banded structure produces a "hole effect" (the
        variogram rises above and/or falls below the sill at larger
        lags).

        Two complementary statistics are therefore computed:

        * ``nugget_deficit`` = ``1 - gamma(1)/sigma^2``.  One-sided:
          a large positive value means strong lag-1 correlation, the
          most common signature of misfit structure.
        * ``max_abs_dev`` = ``max_h |1 - gamma(h)/sigma^2|`` over all
          lags.  Catches hole-effect and longer-range departures that a
          lag-1 (or slope) statistic dilutes or misses entirely.

        Both null distributions depend on the grid shape, so they are
        calibrated by simulating ``n_null`` white-noise fields of the
        same shape.  Empirical p-values (add-one corrected, so the
        smallest reportable p is ``1/(n_null+1)``) are combined with a
        Bonferroni factor of 2: structure is flagged when
        ``min(p) < alpha / 2``.  ``alpha`` is a genuine significance
        level, consistent with the 2D power spectrum test, and replaces
        the earlier fixed slope threshold.

        :param residuals: 2-D array of standardised residuals.
        :param int max_lag: Maximum pixel lag.  Defaults to 20.
        :param int n_null: Number of white-noise simulations for the
            null calibration.  Default 500.
        :param float alpha: Significance level before the Bonferroni
            factor of 2.  Default 0.05.
        :param int random_seed: Seed for the null simulations.
        :returns: dict with keys ``test``, ``lags``, ``semivariance``
            (raw), ``semivariance_norm``, ``sample_variance``,
            ``nugget_deficit``, ``p_nugget``, ``max_abs_dev``,
            ``p_maxdev``, ``min_p``, ``alpha_bonferroni``, ``n_null``,
            ``structured`` (bool), and ``interpretation`` (str).
        """

        def _variogram_curve(field):
            lags = _np.arange(1, max_lag + 1)
            sv = _np.empty(len(lags))
            for k, lag in enumerate(lags):
                h_diffs = (field[:, lag:] - field[:, :-lag]).ravel()
                v_diffs = (field[lag:, :] - field[:-lag, :]).ravel()
                sv[k] = 0.5 * _np.mean(
                    _np.concatenate([h_diffs, v_diffs]) ** 2)
            return lags, sv

        def _stats_pair(field):
            _, sv = _variogram_curve(field)
            svn = sv / (_np.var(field) + 1e-12)
            return (1.0 - svn[0]), _np.max(_np.abs(1.0 - svn))

        lags, semivariance = _variogram_curve(residuals)
        sample_variance    = float(_np.var(residuals))
        semivariance_norm  = semivariance / (sample_variance + 1e-12)

        d_obs, m_obs = _stats_pair(residuals)

        rng    = _np.random.default_rng(random_seed)
        d_null = _np.empty(n_null)
        m_null = _np.empty(n_null)
        sd = residuals.std()
        for i in range(n_null):
            d_null[i], m_null[i] = _stats_pair(
                rng.normal(0.0, sd, residuals.shape))

        # one-sided empirical p-values (add-one correction avoids p = 0)
        p_d = (_np.count_nonzero(d_null >= d_obs) + 1.0) / (n_null + 1.0)
        p_m = (_np.count_nonzero(m_null >= m_obs) + 1.0) / (n_null + 1.0)

        min_p      = min(p_d, p_m)
        alpha_bonf = alpha / 2.0
        structured = bool(min_p < alpha_bonf)

        return {
            'test'              : 'Variogram (MC-calibrated)',
            'lags'              : lags,
            'semivariance'      : semivariance,
            'semivariance_norm' : semivariance_norm,
            'sample_variance'   : sample_variance,
            'nugget_deficit'    : float(d_obs),
            'p_nugget'          : float(p_d),
            'max_abs_dev'       : float(m_obs),
            'p_maxdev'          : float(p_m),
            'min_p'             : float(min_p),
            'alpha_bonferroni'  : alpha_bonf,
            'n_null'            : n_null,
            'structured'        : structured,
            'interpretation'    : (
                f'lag-1 nugget deficit = {d_obs:.3f} (p = {p_d:.4f});  '
                f'max |sill deviation| = {m_obs:.3f} (p = {p_m:.4f});  '
                f'alpha/2 = {alpha_bonf:.3f}  '
                f'-> {"structured" if structured else "random"}'
            ),
        }

    @staticmethod
    def _structure_test_runs(residuals):
        """Wald-Wolfowitz runs test applied per-row and per-column.

        A *run* is a maximal sequence of values all above or all below the
        row/column median.  Too few runs imply positive autocorrelation
        along that axis.  We flag structure when more than 15 % of rows
        *or* columns individually reject randomness at the 5 % level
        (~5 % would be expected by chance).

        .. note::
            This test scans the horizontal and vertical directions only.
            Structure aligned with the grid diagonals is targeted by
            :meth:`_structure_test_acf2d`.

        :param residuals: 2-D array of standardised residuals.
        :returns: dict with keys ``test``, ``frac_rows_significant``,
            ``frac_cols_significant``, ``threshold``, ``structured``
            (bool), and ``interpretation`` (str).
        """

        def _runs_pvalue(x):
            """One-sided p-value for too-few runs (positive autocorrelation)."""
            median = _np.median(x)
            signs = _np.sign(x - median)
            signs = signs[signs != 0]
            n = len(signs)
            if n < 10:
                return 1.0
            n1 = int(_np.sum(signs > 0))
            n2 = int(_np.sum(signs < 0))
            if n1 == 0 or n2 == 0:
                return 0.0
            runs = 1 + int(_np.sum(signs[1:] != signs[:-1]))
            mu_r = 2.0 * n1 * n2 / n + 1.0
            var_r = (2.0 * n1 * n2 * (2.0 * n1 * n2 - n)
                     / (n ** 2 * (n - 1)))
            if var_r <= 0:
                return 1.0
            z = (runs - mu_r) / _np.sqrt(var_r)
            return float(_stats.norm.cdf(z))  # left tail: few runs -> small p

        n_rows, n_cols = residuals.shape
        row_pvals = [_runs_pvalue(residuals[i, :]) for i in range(n_rows)]
        col_pvals = [_runs_pvalue(residuals[:, j]) for j in range(n_cols)]

        frac_rows = float(_np.mean(_np.array(row_pvals) < 0.05))
        frac_cols = float(_np.mean(_np.array(col_pvals) < 0.05))

        threshold = 0.15
        structured = bool((frac_rows > threshold) or (frac_cols > threshold))

        return {
            'test': 'Runs test (per row & column)',
            'frac_rows_significant': frac_rows,
            'frac_cols_significant': frac_cols,
            'threshold': threshold,
            'structured': structured,
            'interpretation': (
                f'rows with structure: {frac_rows:.1%},  '
                f'cols with structure: {frac_cols:.1%}  '
                f'(threshold {threshold:.0%})  '
                f'-> {"structured" if structured else "random"}'
            ),
        }

    @staticmethod
    def _structure_test_acf2d(residuals, max_lag=8):
        """Two-dimensional autocorrelation-function (ACF) scan.

        Generalises a directional lag-1 test to EVERY lag vector
        ``(dy, dx)`` with ``0 <= dy <= max_lag`` and ``|dx| <= max_lag``
        (upper half-plane only, since the ACF is symmetric under
        ``(dy, dx) -> (-dy, -dx)``).  This makes the test sensitive to
        clustered residuals at ANY orientation - horizontal, vertical,
        diagonal, or any oblique angle - and at any separation up to
        ``max_lag`` pixels, so it also catches structure whose
        correlation appears only at larger lags (e.g. regularly spaced
        blobs).

        For each lag vector, the Pearson correlation ``r`` between the
        residual field and its shifted copy is computed.  Under the
        white-noise null hypothesis ``r * sqrt(N)`` is asymptotically
        standard normal, giving a two-sided p-value per lag.  A
        Bonferroni correction over the total number of lags tested is
        applied: structure is flagged when ``min(p) < 0.05 / n_lags``.

        :param residuals: 2-D array of standardised residuals.
        :param int max_lag: Maximum lag in pixels along each axis.
            Default 8, giving ``max_lag * (2 * max_lag + 1) + max_lag``
            = 144 tested lag vectors.
        :returns: dict with keys ``test``, ``acf_map`` (2-D array of
            ``r`` values, NaN where untested), ``max_lag``, ``best_lag``
            (tuple ``(dy, dx)`` of the most significant lag),
            ``best_r``, ``best_angle_deg`` (orientation of the
            strongest correlation, degrees from the +phase axis),
            ``min_p``, ``alpha_bonferroni``, ``n_lags_tested``,
            ``structured`` (bool), and ``interpretation`` (str).
        """
        m, n = residuals.shape
        acf_map = _np.full((max_lag + 1, 2 * max_lag + 1), _np.nan)
        p_values = []
        best = {'z': 0.0, 'dy': 0, 'dx': 1, 'r': 0.0}

        for dy in range(0, max_lag + 1):
            for dx in range(-max_lag, max_lag + 1):
                if dy == 0 and dx <= 0:
                    continue  # skip (0,0) and mirror duplicates
                if dx >= 0:
                    a = residuals[dy:, dx:]
                    b = residuals[:m - dy, :n - dx]
                else:
                    a = residuals[dy:, :n + dx]
                    b = residuals[:m - dy, -dx:]
                a = a.ravel()
                b = b.ravel()
                npairs = a.size
                if npairs < 10:
                    continue
                r = float(_np.corrcoef(a, b)[0, 1])
                z = r * _np.sqrt(npairs)
                p = float(2.0 * _stats.norm.sf(abs(z)))
                acf_map[dy, dx + max_lag] = r
                p_values.append(p)
                if abs(z) > abs(best['z']):
                    best = {'z': z, 'dy': dy, 'dx': dx, 'r': r}

        n_lags = len(p_values)
        alpha_bonf = 0.05 / n_lags
        min_p = min(p_values)
        structured = bool(min_p < alpha_bonf)
        angle = float(_np.degrees(_np.arctan2(best['dy'], best['dx'])))

        return {
            'test': '2D ACF scan (all orientations & lags)',
            'acf_map': acf_map,
            'max_lag': max_lag,
            'best_lag': (best['dy'], best['dx']),
            'best_r': float(best['r']),
            'best_angle_deg': angle,
            'min_p': float(min_p),
            'alpha_bonferroni': alpha_bonf,
            'n_lags_tested': n_lags,
            'structured': structured,
            'interpretation': (
                f"strongest r = {best['r']:+.3f} at lag "
                f"(dchan={best['dy']}, dphase={best['dx']}), "
                f"orientation ~ {angle:.0f} deg;  "
                f"min-p = {min_p:.2e} over {n_lags} lags "
                f"(Bonferroni alpha = {alpha_bonf:.1e})  "
                f"-> {'structured' if structured else 'random'}"
            ),
        }

    @staticmethod
    def _structure_test_normality_ljungbox(residuals):
        """Kolmogorov-Smirnov normality test + Ljung-Box autocorrelation test.

        * **KS test** - compares the marginal distribution of the
          standardised residuals against a fitted normal.
        * **Ljung-Box test** - detects autocorrelation at lags 1-20 in the
          residuals flattened in row-major (phase-fast) order.

        :param residuals: 2-D array of standardised residuals.
        :returns: dict with keys ``test``, ``ks_stat``, ``ks_pvalue``,
            ``ljungbox_min_pvalue`` (``None`` if statsmodels absent),
            ``structured`` (bool), and ``interpretation`` (str).

        .. note::
            The Ljung-Box test requires ``statsmodels``.  If absent, only
            the KS test contributes to the verdict.
        """
        flat = residuals.ravel()
        mu_hat = float(flat.mean())
        sd_hat = float(flat.std())

        # Pass a frozen CDF rather than ``'norm', args=(mu, sd)``: scipy 1.18
        # maps the string 'norm' to ``special.ndtr``, which takes no
        # loc/scale and raises a TypeError when ``args`` are supplied.
        ks_stat, ks_p = _stats.kstest(flat, _stats.norm(loc=mu_hat, scale=sd_hat).cdf)
        ks_stat = float(ks_stat)
        ks_p = float(ks_p)

        if _HAS_STATSMODELS:
            lb_result = _acorr_ljungbox(flat, lags=20, return_df=True)
            lb_min_p = float(lb_result['lb_pvalue'].min())
            structured = bool((ks_p < 0.01) or (lb_min_p < 0.05))
            lb_str = f'Ljung-Box min-p = {lb_min_p:.4f}'
        else:
            lb_min_p = None
            structured = bool(ks_p < 0.01)
            lb_str = ('Ljung-Box skipped '
                      '(pip install statsmodels to enable)')

        return {
            'test': 'KS normality + Ljung-Box autocorr.',
            'ks_stat': ks_stat,
            'ks_pvalue': ks_p,
            'ljungbox_min_pvalue': lb_min_p,
            'structured': structured,
            'interpretation': (
                f'KS p = {ks_p:.4f},  {lb_str}  '
                f'-> {"structured" if structured else "random"}'
            ),
        }

    # -- Diagnostic figure --------------------------------------------------

    def _draw_diagnostics(self):
        """Draw the diagnostic figure of the structure tests.

        Called automatically from :meth:`finalize`, after the tests.  Tests
        1 to 4 occupy the first two columns (tests 1 and 2 on the top row,
        tests 3 and 4 on the bottom row), and the summary the full third
        column.  The panels follow the numbering of the printed report:

        * **Test 1** - 2-D power spectrum (log scale).
        * **Test 2** - normalised empirical variogram with a dashed
          pure-nugget reference line.
        * **Test 3** - 2D autocorrelation-function map over all lag
          vectors, with the most significant lag marked.
        * **Test 4** - normal Q-Q plot of the flattened residuals
          (KS normality + Ljung-Box test).
        * **Summary** - all six test verdicts in report order: first the
          tests shown in panels (1-4), then the tests without a panel
          (5: Moran's I, 6: runs test), followed by the overall verdict.
        """
        results = self.residual_structure_results
        results_map = {t['test']: t for t in results}
        verdict = self.residual_structure_verdict
        colour_map = {
            'STRUCTURED': '#d62728',
            'RANDOM': '#2ca02c',
            'INCONCLUSIVE': '#ff7f0e',
        }

        fig = self._fig
        fig.suptitle(
            f'XPSI Residual Structure Analysis - Verdict: {verdict}',
            fontsize=20, fontweight='bold',
            color=colour_map.get(verdict, 'black'))
        gs = self._gs

        residuals = self._residuals

        # -- Panel (0,0): Test 1 - 2-D power spectrum ----------------------
        ax2 = fig.add_subplot(gs[0, 0])
        ps_t = results_map.get('2D Power Spectrum (MC-calibrated)', {})
        if 'power_map' in ps_t:
            ax2.imshow(_np.log1p(ps_t['power_map']),
                       aspect='auto', origin='lower', cmap='inferno')
            ax2.set_title(
                f"Test 1: 2D Power Spectrum (log)\n"
                f"p_periodic = {ps_t.get('p_periodic', float('nan')):.4f}, "
                f"p_smooth = {ps_t.get('p_smooth', float('nan')):.4f}")
        else:
            ax2.text(0.5, 0.5, 'No data', ha='center', va='center',
                     transform=ax2.transAxes)
            ax2.set_title('Test 1: 2D Power Spectrum')
        ax2.set_xlabel(r'$k_\phi$')
        ax2.set_ylabel(r'$k_\mathrm{ch}$')

        # -- Panel (0,1): Test 2 - variogram (normalised) ------------------
        ax3 = fig.add_subplot(gs[0, 2])
        vg_t = results_map.get('Variogram (MC-calibrated)', {})
        if 'lags' in vg_t:
            ax3.plot(vg_t['lags'], vg_t['semivariance_norm'],
                     'o-', color='#1f77b4', ms=4, label='empirical')
            ax3.axhline(1.0, ls='--', color='gray', lw=1,
                        label='pure nugget (random)')
            ax3.set_title(
                'Test 2: Variogram (normalised)\n'
                f"p_nugget = {vg_t.get('p_nugget', float('nan')):.4f}, "
                f"p_maxdev = {vg_t.get('p_maxdev', float('nan')):.4f}")
        else:
            ax3.text(0.5, 0.5, 'No data', ha='center', va='center',
                     transform=ax3.transAxes)
            ax3.set_title('Test 2: Variogram (normalised)')
        ax3.set_xlabel('Lag (pixels)')
        ax3.set_ylabel('Normalised semivariance')
        ax3.legend(fontsize=9, loc='best')
        self._y_axis_right(ax3)

        # -- Panel (1,0): Test 3 - 2D ACF scan ------------------------------
        ax4 = fig.add_subplot(gs[1, 0])
        acf_t = results_map.get('2D ACF scan (all orientations & lags)', {})
        if 'acf_map' in acf_t:
            L = acf_t['max_lag']
            amap = acf_t['acf_map']
            amax = _np.nanmax(_np.abs(amap))
            im = ax4.imshow(amap, aspect='auto', origin='lower',
                            cmap='RdBu_r', vmin=-amax, vmax=amax,
                            extent=[-L - 0.5, L + 0.5, -0.5, L + 0.5])
            fig.colorbar(im, ax=ax4, fraction=0.046, pad=0.04,
                         label=r'Pearson $r$')
            bdy, bdx = acf_t['best_lag']
            ax4.plot(bdx, bdy, 'x', color='k', ms=10, mew=2,
                     label=f'strongest lag ({bdy},{bdx})')
            ax4.set_title('Test 3: 2D ACF scan\n'
                          f"min-p = {acf_t.get('min_p', float('nan')):.2e}")
            ax4.set_xlabel(r'$\Delta$ phase bin')
            ax4.set_ylabel(r'$\Delta$ channel bin')
            ax4.legend(fontsize=9, loc='upper right')
        else:
            ax4.text(0.5, 0.5, 'No data', ha='center', va='center',
                     transform=ax4.transAxes)
            ax4.set_title('Test 3: 2D ACF scan')

        # -- Panel (1,1): Test 4 - Q-Q plot --------------------------------
        ax5 = fig.add_subplot(gs[1, 2])
        flat = residuals.ravel()
        (osm, osr), (slope, intercept, _) = _stats.probplot(flat, dist='norm')
        ax5.plot(osm, osr, '.', ms=1, alpha=0.4, color='#1f77b4')
        ax5.plot(osm, slope * _np.array(osm) + intercept,
                 'r-', lw=1.5, label='normal reference')
        ks_t = results_map.get('KS normality + Ljung-Box autocorr.', {})
        ks_p = ks_t.get('ks_pvalue', float('nan'))
        ax5.set_title(f'Test 4: Normal Q-Q  (KS p = {ks_p:.4f})')
        ax5.set_xlabel('Theoretical quantiles')
        ax5.set_ylabel('Sample quantiles')
        ax5.legend(fontsize=9)
        self._y_axis_right(ax5)

        # -- Right column: summary of ALL tests in report order ------------
        ax6 = fig.add_subplot(gs[:, 4])
        ax6.axis('off')
        ax6.set_title('Summary (report order)')
        self._draw_summary(ax6, results, verdict,
                           colour_map.get(verdict, 'black'))

    @staticmethod
    def _y_axis_right(ax):
        """ Put the y ticks, tick labels and axis label on the right. """
        ax.yaxis.tick_right()
        ax.yaxis.set_label_position('right')

    # Plain-language description of each test, shown in grey in the summary
    # column (keyed by the 'test' name returned by each test).  For the
    # tests with a panel, the description says how to read the panel.
    _TEST_DESCRIPTIONS = {
        '2D Power Spectrum (MC-calibrated)': (
            'Splits the residual map into spatial frequencies (k_\u03c6 and '
            'k_ch); bins with zero frequency at the centre; '
            'colour = log of the power. Noise spreads power evenly over the '
            'map. A repeating pattern gives isolated bright spots, and smooth '
            'blobs a bright centre.'),
        'Variogram (MC-calibrated)': (
            'Compares how different two residuals are, as a function of '
            'their distance. Distance between pixels (lag) as a function '
            'of mean squared difference divided by its noise value '
            '(normalised semi-variance). For noise, the curve stays on the '
            'dashed line at 1 at every distance. Structure makes nearby '
            'residuals more alike, so the curve drops below 1 at small lags.'),
        '2D ACF scan (all orientations & lags)': (
            'Correlates the residual map with copies of itself shifted in '
            'phase and channel. Panel: x = shift in phase bins, y = shift in '
            'channel bins; colour = correlation r (red positive, blue '
            'negative); \u00d7 = strongest shift. Noise gives r \u2248 0 '
            'everywhere. Structure gives a red patch around the origin, '
            'elongated in the direction of the structure.'),
        'KS normality + Ljung-Box autocorr.': (
            'KS compares the shape of the residual distribution with a normal '
            'distribution of the same mean and width. Panel (Q-Q plot): '
            'x = values expected for a normal distribution (theoretical '
            'quantiles), y = sorted residuals (sample quantiles). If the shape '
            'is normal, the points follow the red line; tails bending away '
            'from it mean too many or too few large residuals. Ljung-Box (no '
            'panel) checks whether consecutive residuals along phase are '
            'correlated.'),
        "Moran's I (spatial autocorrelation)": (
            'Checks whether each residual resembles its 8 neighbours. '
            'I \u2248 0 for noise; I > 0 when similar residuals sit together.'),
        'Runs test (per row & column)': (
            'Along each channel and each phase bin, counts the runs of '
            'residuals above or below the median. Noise switches side often; '
            'too few runs means neighbouring residuals are correlated.'),
    }

    @staticmethod
    def _format_p(p, n_null=None):
        """ Format a p-value for the summary. With a Monte-Carlo p-value
        at its smallest possible value ``1/(n_null+1)``, show an upper
        bound. """
        if n_null is not None and p <= 1.0 / (n_null + 1.0) * (1 + 1e-9):
            return f'p \u2264 {p:.2g}'
        if p == 0.0:
            return 'p < 1e-300'
        if p < 1e-3:
            return f'p = {p:.1e}'
        return f'p = {p:.3f}'

    @classmethod
    def _summary_quantification(cls, t):
        """ Values of a test and the threshold that decides its verdict,
        as shown in the summary column. """
        name = t['test']
        fp = cls._format_p
        if name.startswith('2D Power Spectrum'):
            n = t['n_null']
            return ('Strongest single frequency: {}. Low-frequency power: '
                    '{:.0%} (95% of noise: \u2264 {:.0%}), {}. Structured if '
                    'either p < {:.3g}.'.format(
                        fp(t['p_periodic'], n), t['f_lowk'], t['f_null_q95'],
                        fp(t['p_smooth'], n), t['alpha_bonferroni']))
        if name.startswith('Variogram'):
            n = t['n_null']
            return ('Value at lag 1: {:.2f} ({}). Largest departure from 1 '
                    'over lags 1\u2013{}: {:.2f} ({}). Structured if either '
                    'p < {:.3g}.'.format(
                        1.0 - t['nugget_deficit'], fp(t['p_nugget'], n),
                        len(t['lags']), t['max_abs_dev'], fp(t['p_maxdev'], n),
                        t['alpha_bonferroni']))
        if name.startswith('2D ACF'):
            dy, dx = t['best_lag']
            bins = lambda n, axis: '{} {} bin{}'.format(
                n, axis, '' if abs(n) == 1 else 's')
            return ('Strongest correlation: r = {:+.2f} for a shift of {} and '
                    '{}, {}. Structured if p < {:.1e} (0.05 divided by the {} '
                    'shifts tested).'.format(
                        t['best_r'], bins(dx, 'phase'), bins(dy, 'channel'),
                        fp(t['min_p']), t['alpha_bonferroni'],
                        t['n_lags_tested']))
        if name.startswith('KS normality'):
            lb = t['ljungbox_min_pvalue']
            lb_str = ('Ljung-Box: smallest {} over lags 1\u201320 (structured '
                      'if < 0.05).'.format(fp(lb)) if lb is not None else
                      'Ljung-Box: skipped (install statsmodels).')
            return 'KS: {} (structured if < 0.01). {}'.format(
                fp(t['ks_pvalue']), lb_str)
        if name.startswith("Moran"):
            if t['structured'] is None:
                return 'Skipped: install esda and libpysal.'
            return ('I = {:.2f}, {}, from 999 random reshufflings of the map. '
                    'Structured if p < 0.05.'.format(
                        t['I'], fp(t['p_value'], 999)))
        if name.startswith('Runs test'):
            return ('{:.0%} of channels and {:.0%} of phase bins have too few '
                    'runs (about 5% expected by chance). Structured if either '
                    'exceeds {:.0%}.'.format(
                        t['frac_rows_significant'], t['frac_cols_significant'],
                        t['threshold']))
        return t.get('interpretation', '')

    def _draw_summary(self, ax, results, verdict, verdict_colour):
        """ Write the summary of the tests in the (axis-less) column ``ax``.

        For each test: its number, name and verdict (font size 11), a
        plain-language description (grey, size 9) and the values that
        decide the verdict (size 9); then the combined verdict.  Text is
        laid out from the top of the column with line heights in physical
        units, so it keeps its size whatever the column height.
        """
        fig = self._fig
        w_in, h_in = fig.get_size_inches()
        pos = ax.get_position()
        ax_w_pt = pos.width * w_in * 72.0
        ax_h_pt = pos.height * h_in * 72.0

        def line_h(fs):                       # line height, axes fraction
            return 1.18 * fs / ax_h_pt

        def n_chars(fs, frac=0.96):           # characters per line
            return max(20, int(frac * ax_w_pt / (0.55 * fs)))

        x0 = 0.02
        y = 1.0 - 0.5 * line_h(11)

        for i, t in enumerate(results, start=1):
            if t['structured'] is True:
                mark, col = 'STRUCTURED', '#d62728'
            elif t['structured'] is False:
                mark, col = 'Random', '#2ca02c'
            else:
                mark, col = 'skipped', 'gray'

            # Name followed by the verdict, the whole title in the colour of
            # the verdict
            label = t['test'].split(' (')[0]
            for line in _textwrap.wrap(f'Test {i}: {label} - {mark}',
                                       width=n_chars(11)):
                ax.text(x0, y, line, transform=ax.transAxes, fontsize=11,
                        va='top', color=col)
                y -= line_h(11)

            for text, colour in ((self._TEST_DESCRIPTIONS.get(t['test'], ''),
                                  '0.45'),
                                 (self._summary_quantification(t), 'black')):
                for line in _textwrap.wrap(text, width=n_chars(9)):
                    ax.text(x0, y, line, transform=ax.transAxes, fontsize=9,
                            va='top', color=colour)
                    y -= line_h(9)
            y -= 0.4 * line_h(11)

        definitive = [t for t in results if t['structured'] is not None]
        n_structured = len([t for t in definitive if t['structured']])
        ax.text(x0, y, f'VERDICT: {verdict}', transform=ax.transAxes,
                fontsize=11, va='top', fontweight='bold', color=verdict_colour)
        y -= line_h(11)
        rule = ('{} of {} tests flagged structure. STRUCTURED if \u2265 3, '
                'RANDOM if 0, otherwise INCONCLUSIVE.'.format(
                    n_structured, len(definitive)))
        for line in _textwrap.wrap(rule, width=n_chars(9)):
            ax.text(x0, y, line, transform=ax.transAxes, fontsize=9,
                    va='top')
            y -= line_h(9)
