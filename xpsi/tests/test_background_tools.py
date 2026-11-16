from types import SimpleNamespace

import numpy as np
import matplotlib.pyplot as plt
import pytest

import xpsi
from xpsi.utilities.BackgroundTools import plotBackgroundSpectrum


def _make_model(data, n_used, seed=0):
    """ Minimal stand-in for a loaded X-PSI model with a single signal.

    Only the attributes accessed by ``plotBackgroundSpectrum`` are provided.
    The model arrays have ``n_used`` channels, i.e. the number of channels
    actually used by the signal (``len(data.channels)``).
    """
    rng = np.random.default_rng(seed)
    component_phases = np.linspace(0.0, 1.0, 5)
    signal = SimpleNamespace(
        data=data,
        num_components=1,
        signals=(rng.uniform(1.0, 2.0, size=(n_used, component_phases.size)),),
        phases=(component_phases,),
        expected_counts=rng.uniform(50.0, 100.0, size=(n_used, 1)),
        background_signal_given_support=rng.uniform(5.0, 10.0, size=n_used),
        instrument=SimpleNamespace(channel_edges=np.array([np.arange(n_used) * 0.1 + 0.3,
                                                           np.arange(n_used) * 0.1 + 0.4])),
        _support=None,
    )
    return SimpleNamespace(signal=signal, likelihood=lambda p, reinitialise=True: None)


def _grouped_trimmed_data():
    """ Grouped data in the state left by ``Signal.__init__`` when grouped channels
    are ignored: ``data.channels`` trimmed, ``data.counts`` not trimmed. """
    n_rows = 10
    counts = np.arange(100.0, 100.0 + n_rows).reshape(n_rows, 1)
    data = xpsi.Data(counts, channels=np.arange(n_rows), phases=np.array([0.0, 1.0]),
                     first=0, last=n_rows - 1, exposure_time=1000.0)
    data.errors = np.full((n_rows, 1), 7.0)
    data._quality = np.zeros(n_rows + 5)
    data._grouping = np.ones(n_rows + 5)
    data._grpdata = True
    # Mimic xpsi/Signal.py grouping branch
    data._channels = data._channels[:-2]
    data._first = data._channels[0]
    data._last = data._channels[-1]
    return data


def _quality_offset_data():
    """ Quality-filtered data whose channel numbers do not start at 0, with
    ``counts`` and ``channels`` of equal length and one empty channel. """
    n_rows = 10
    counts = np.arange(100.0, 100.0 + n_rows).reshape(n_rows, 1)
    counts[3] = 0.0
    data = xpsi.Data(counts, channels=np.arange(30, 30 + n_rows), phases=np.array([0.0, 1.0]),
                     first=0, last=n_rows - 1, exposure_time=1000.0)
    data._quality = np.zeros(30 + n_rows + 5)
    return data


@pytest.mark.parametrize('make_data', [_grouped_trimmed_data, _quality_offset_data],
                         ids=['grouped_trimmed', 'quality_offset'])
def test_background_plot_channel_selection_and_residuals(make_data):
    'Data spectrum, errors and residuals use only the channels used by the signal.'

    data = make_data()
    n_used = len(data.channels)
    model = _make_model(data, n_used)
    signal = model.signal

    fig, (ax, ax_res) = plotBackgroundSpectrum(XPSI_model=model,
                                               samples_path=None,
                                               InstrumentName=None,
                                               plot_params=[0.0],
                                               plot_range=False,
                                               plot_data_errors=True,
                                               plot_residuals=True)

    # Expected data spectrum and uncertainties on the used channels only
    expected_data = data.counts[:n_used].sum(axis=1)
    if getattr(data, 'errors', None) is not None:
        expected_sigma = np.sqrt((data.errors[:n_used]**2).sum(axis=1))
    else:
        expected_sigma = np.sqrt(np.where(expected_data > 0, expected_data, 1.0))
    model_counts = signal.expected_counts.sum(axis=1) + signal.background_signal_given_support
    expected_residuals = (expected_data - model_counts) / expected_sigma

    data_line = ax.containers[0].lines[0]
    assert np.allclose(data_line.get_ydata(), expected_data)

    residuals = ax_res.containers[0].lines[0].get_ydata()
    assert np.all(np.isfinite(residuals))
    assert np.allclose(residuals, expected_residuals)

    plt.close(fig)


def test_background_plot_default_returns_single_axis():
    'Without residuals, the original ``fig, ax`` return signature is kept.'

    data = _quality_offset_data()
    model = _make_model(data, len(data.channels))
    fig, ax = plotBackgroundSpectrum(XPSI_model=model,
                                     samples_path=None,
                                     InstrumentName=None,
                                     plot_params=[0.0],
                                     plot_range=False)
    assert not isinstance(ax, tuple)
    plt.close(fig)
