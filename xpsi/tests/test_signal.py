"""Unit tests for xpsi.Signal (everything except ``__call__``).

Run with:  pytest -v test_signal.py

In real use, ``Data`` and ``Instrument`` are usually loaded from files
(``Data.load``, ``Instrument.from_ogip_fits``). Here no input files are
needed: both classes are built directly from arrays (``make_data``,
``make_instrument``). The writer tests only write to ``tmp_path``.
The other helpers are:

* ``FakeInterstellar``: an ``Interstellar`` with a constant attenuation;
* ``ConcreteSignal``: ``Signal.__call__`` is abstract, so a trivial one is
  provided;
* ``StubSignal``: a bare object with ``energy_edges``, for
  ``construct_energy_array``.

Toy geometry used throughout
----------------------------
* 6 output channels (0..5), 8 input energy intervals (9 edges), 3 phase bins
* response matrix is strictly positive unless a test says otherwise
* ``Signal.register`` tests call the compiled xpsi energy integrators, so they
  need a built xpsi (as does ``import xpsi`` itself).
"""

import numpy as np
import pytest
from astropy.io import fits

from xpsi.Background import Background
from xpsi.Data import Data
from xpsi.Instrument import Instrument, ChannelError, EdgesError
from xpsi.Interstellar import Interstellar
from xpsi.Parameter import Parameter
from xpsi.ParameterSubspace import ParameterSubspace
from xpsi.Signal import Signal, LikelihoodError, construct_energy_array

N_CHAN, N_IN, N_PHASE = 6, 8, 3
EDGES = np.logspace(np.log10(0.2), np.log10(3.0), N_IN + 1)

# --------------------------------------------------------------------------
# Fakes
# --------------------------------------------------------------------------
class FakeInterstellar(Interstellar):
    def __init__(self, factor=0.5, *params):
        self.factor = factor
        super().__init__(*params)

    def attenuation(self, energies):
        return self.factor * np.ones_like(energies)


class FakeBackground(Background):
    """A background of 1 count per input interval and phase bin."""

    def __call__(self, energy_edges, phases):
        self.incident_background = np.ones((len(energy_edges) - 1,
                                            len(phases) - 1))


class ConcreteSignal(Signal):
    """Signal has an abstract ``__call__``; provide a trivial one."""

    def __call__(self, **kwargs):
        pass


def default_matrix():
    return (np.arange(1, N_CHAN * N_IN + 1, dtype=np.double)
            .reshape(N_CHAN, N_IN) / 100.0)


def make_data(**kw):
    """Real ``Data`` built from arrays.

    ``last`` defaults to ``first + number of counts rows - 1``, so
    ``index_range == (first, first + rows)``.
    """
    args = dict(counts=np.ones((N_CHAN, N_PHASE)),
                channels=np.arange(N_CHAN),
                phases=np.linspace(0.0, 1.0, N_PHASE + 1),
                first=0,
                exposure_time=1.0)
    args.update(kw)
    args.setdefault('last', args['first'] + len(args['counts']) - 1)
    return Data(**args)


def make_instrument(matrix=None, edges=None, twod_edges=False):
    if matrix is None:
        matrix = default_matrix()
    if edges is None:
        edges = EDGES.copy()
        if twod_edges:
            edges = np.array([EDGES[:-1], EDGES[1:]])
    return Instrument(matrix, edges, np.arange(N_CHAN), twod_edges=twod_edges)


def make_signal(data=None, instrument=None, **kw):
    data = make_data() if data is None else data
    instrument = make_instrument() if instrument is None else instrument
    return ConcreteSignal(data, instrument, **kw)


@pytest.fixture
def sig():
    return make_signal()


@pytest.fixture
def energies():
    return np.logspace(np.log10(0.1), np.log10(5.0), 40)


def flux(energies, value=1.0):
    """Specific flux component of shape (n_energies, n_phases)."""
    return np.full((len(energies), N_PHASE), value, dtype=np.double)


# --------------------------------------------------------------------------
# Construction: type checks
# --------------------------------------------------------------------------
class TestConstructionTypes:

    def test_is_parameter_subspace(self, sig):
        """A Signal is a ParameterSubspace, so its parameters join the model's.
        """
        assert isinstance(sig, ParameterSubspace)

    def test_bad_data_type(self):
        """Data that is not a ``Data`` instance raises TypeError."""
        with pytest.raises(TypeError, match='Invalid type'):
            ConcreteSignal('not data', make_instrument())

    def test_bad_instrument_type(self):
        """An instrument that is not an ``Instrument`` instance raises
        TypeError.
        """
        with pytest.raises(TypeError, match='Invalid type'):
            ConcreteSignal(make_data(), 'not an instrument')

    def test_bad_background_type(self):
        """A background that is not a ``Background`` instance raises TypeError.
        """
        with pytest.raises(TypeError, match='Invalid type'):
            make_signal(background='nope')

    def test_bad_interstellar_type(self):
        """An interstellar model that is not an ``Interstellar`` instance
        raises TypeError.
        """
        with pytest.raises(TypeError, match='Invalid type'):
            make_signal(interstellar='nope')

    def test_no_background_is_none(self, sig):
        """Without a background, the ``background`` property is None."""
        assert sig.background is None

    def test_no_interstellar_is_none(self, sig):
        """Without an interstellar model, the ``interstellar`` property is
        None.
        """
        assert sig.interstellar is None

    def test_data_and_instrument_accessors(self):
        """The ``data`` and ``instrument`` properties return the objects given
        to the constructor (no trimming requested).
        """
        d, i = make_data(), make_instrument()
        s = ConcreteSignal(d, i)
        assert s.data is d
        assert s.instrument is i

    def test_channel_size_mismatch_raises(self):
        """Data channels that do not match the instrument channels (and no
        channel trimming) raise ChannelError.
        """
        bad = make_data(channels=np.arange(N_CHAN - 1))
        with pytest.raises(ChannelError,
                           match='does not match the declared size'):
            make_signal(data=bad)

    # Instrument names: Data.load sets ``data.instrument`` and
    # Instrument.from_ogip_fits sets ``instrument.name`` (both per instance);
    # they are compared only when both exist.
    @staticmethod
    def _named(data_name=None, instrument_name=None):
        d, i = make_data(), make_instrument()
        if data_name is not None:
            d.instrument = data_name
        if instrument_name is not None:
            i.name = instrument_name
        return d, i

    def test_instrument_name_mismatch_raises(self):
        """Data and instrument naming different instruments raise
        AssertionError.
        """
        d, i = self._named('A', 'B')
        with pytest.raises(AssertionError, match='different instruments'):
            ConcreteSignal(d, i)

    def test_instrument_name_match_accepted(self):
        """Data and instrument naming the same instrument are accepted."""
        d, i = self._named('A', 'A')
        assert ConcreteSignal(d, i).instrument is i

    @pytest.mark.parametrize('data_name, instrument_name',
                             [('A', None), (None, 'B')])
    def test_name_check_skipped_if_one_side_unnamed(self, data_name,
                                                    instrument_name):
        """The instrument name check is skipped when only one of data and
        instrument carries a name.
        """
        d, i = self._named(data_name, instrument_name)
        assert ConcreteSignal(d, i).instrument is i


# --------------------------------------------------------------------------
# Stokes handling
# --------------------------------------------------------------------------
class TestStokes:

    FLAGS = ['isI', 'isQ', 'isU', 'isQn', 'isUn']

    @pytest.mark.parametrize('stokes, flag', [('I', 'isI'), ('Q', 'isQ'),
                                              ('U', 'isU'), ('Qn', 'isQn'),
                                              ('Un', 'isUn')])
    def test_exactly_one_flag_set(self, stokes, flag):
        """Each ``stokes`` value sets its own flag and leaves the four others
        False.
        """
        s = make_signal(stokes=stokes)
        for f in self.FLAGS:
            assert getattr(s, f) == (f == flag)

    def test_default_is_I(self, sig):
        """Without a ``stokes`` argument, only ``isI`` is set."""
        for f in self.FLAGS:
            assert getattr(sig, f) == (f == 'isI')

    @pytest.mark.parametrize('bad', ['V', 'i', '', 'IQ', None, 1])
    def test_invalid_stokes_raises(self, bad):
        """A ``stokes`` value other than I, Q, U, Qn or Un raises TypeError."""
        with pytest.raises(TypeError, match='param stokes'):
            make_signal(stokes=bad)

    @pytest.mark.parametrize('flag', FLAGS)
    def test_flag_setters(self, sig, flag):
        """Each Stokes flag can be set to True and back to False."""
        setattr(sig, flag, True)
        assert getattr(sig, flag) is True
        setattr(sig, flag, False)
        assert getattr(sig, flag) is False


# --------------------------------------------------------------------------
# Phase shift parameter
# --------------------------------------------------------------------------
class TestPhaseShift:

    def test_parameter_exists(self, sig):
        """Every Signal has a ``phase_shift`` parameter."""
        assert 'phase_shift' in sig.names

    def test_default_fixed_at_zero(self, sig):
        """Without bounds or value, ``phase_shift`` is fixed at zero."""
        assert sig.get_param('phase_shift').fixed is True
        assert sig['phase_shift'] == 0.0

    def test_fixed_value_from_values(self):
        """A value in ``values`` fixes ``phase_shift`` at that value."""
        s = make_signal(values={'phase_shift': 0.25})
        assert s['phase_shift'] == pytest.approx(0.25)

    def test_free_with_bounds(self):
        """Bounds make ``phase_shift`` a free parameter with those bounds and
        no value yet.
        """
        s = make_signal(bounds={'phase_shift': (0.0, 1.0)})
        p = s.get_param('phase_shift')
        assert p.fixed is False
        assert p.bounds == (0.0, 1.0)
        assert p.value is None

    def test_bounds_narrower_than_cycle_ok(self):
        """Bounds spanning less than one cycle, including negative phases, are
        accepted.
        """
        s = make_signal(bounds={'phase_shift': (-0.1, 0.1)})
        assert s.get_param('phase_shift').bounds == (-0.1, 0.1)

    @pytest.mark.parametrize('bounds', [(0.0, 1.5), (-1.0, 1.0), (0.5, 0.0)])
    def test_bounds_wider_than_cycle_or_reversed_raises(self, bounds):
        """Bounds spanning more than one cycle, or with lower above upper,
        raise ValueError.
        """
        with pytest.raises(ValueError, match='maximum of one cycle'):
            make_signal(bounds={'phase_shift': bounds})

    @pytest.mark.parametrize('bounds', [(0.0, np.inf), (-np.inf, 0.0),
                                        (np.nan, 0.5)])
    def test_infinite_bounds_raise(self, bounds):
        """Infinite or NaN bounds raise ValueError."""
        with pytest.raises(ValueError, match='must be finite'):
            make_signal(bounds={'phase_shift': bounds})

    def test_none_bounds_raise(self):
        """Bounds of ``(None, None)`` raise ValueError: a free phase shift
        needs both bounds.
        """
        with pytest.raises(ValueError, match='must be specified'):
            make_signal(bounds={'phase_shift': (None, None)})


# --------------------------------------------------------------------------
# Support (background support array)
# --------------------------------------------------------------------------
class TestSupport:

    def test_default_support_shape(self, sig):
        """The default background support has one row per data channel and two
        columns.
        """
        assert sig._support.shape == (N_CHAN, 2)

    def test_default_support_values(self, sig):
        """The default support is ``[0, -1]`` per channel, i.e. no constraint
        on the background.
        """
        assert np.all(sig._support[:, 0] == 0.0)
        assert np.all(sig._support[:, 1] == -1.0)

    def test_explicit_support_kept(self):
        """A support of the data's size is stored as given."""
        support = np.column_stack([np.zeros(N_CHAN), np.ones(N_CHAN)])
        assert make_signal(support=support)._support is support

    def test_wrong_support_shape_raises(self):
        """A support whose size matches neither data nor instrument raises
        TypeError.
        """
        with pytest.raises(TypeError, match='same shape'):
            make_signal(support=np.zeros((N_CHAN + 10, 2)))


# --------------------------------------------------------------------------
# Data channels a subset of the instrument channels
# --------------------------------------------------------------------------
class TestChannelSubset:
    """Data channels 1..4 (4 rows, first=1) with instrument channels 0..5.

    ``first=1`` makes the data ``index_range`` (1, 5) match the instrument
    rows. (The ``Data`` loaders always set ``first=0``; that case is not
    tested here.)
    """

    N_SUB = 4

    def _signal(self, **kw):
        data = make_data(counts=np.ones((self.N_SUB, N_PHASE)),
                         channels=np.arange(1, self.N_SUB + 1), first=1)
        return make_signal(data=data, **kw)

    @staticmethod
    def _support(n_rows):
        # row i is [i, 10 + i], so slicing is visible in the values
        return np.column_stack([np.arange(n_rows, dtype=float),
                                np.arange(n_rows, dtype=float) + 10.0])

    def test_instrument_index_range(self):
        """The instrument rows used are those of the data channels, 1 to 4."""
        assert self._signal()._instrument_index_range_channels == (1, 5)

    def test_default_support_has_data_rows(self):
        """The default support has one row per data channel, not per instrument
        channel.
        """
        assert self._signal()._support.shape == (self.N_SUB, 2)

    def test_support_of_data_size_kept(self):
        """A support with one row per data channel is stored as given."""
        support = self._support(self.N_SUB)
        assert self._signal(support=support)._support is support

    def test_support_of_instrument_size_is_sliced(self):
        """A support with one row per instrument channel is cut to the data
        channels' rows.
        """
        support = self._support(N_CHAN)
        np.testing.assert_array_equal(self._signal(support=support)._support,
                                      support[1:5])


# --------------------------------------------------------------------------
# min_channel / max_channel trimming of data and response
# --------------------------------------------------------------------------
class TestTrimming:
    """Channels 0..5 trimmed by ``Signal(min_channel=..., max_channel=...)``.

    ``min_channel`` alone is not tested: ``Data.trim_data`` then fails on
    its default ``max_channel=-1``.
    """

    def test_trim_both_data(self):
        """``min_channel`` and ``max_channel`` keep data channels 2 to 4 and
        their counts rows.
        """
        s = make_signal(min_channel=2, max_channel=4)
        np.testing.assert_array_equal(s.data.channels, [2, 3, 4])
        assert s.data.counts.shape == (3, N_PHASE)
        assert s.data.index_range == (0, 3)

    def test_trim_both_instrument(self):
        """The response is trimmed to the same channels 2 to 4, rows included.
        """
        s = make_signal(min_channel=2, max_channel=4)
        np.testing.assert_array_equal(s.instrument.channels, [2, 3, 4])
        np.testing.assert_array_equal(s.instrument.matrix,
                                      default_matrix()[2:5])
        assert s._instrument_index_range_channels == (0, 3)

    def test_trim_max_only(self):
        """``max_channel`` alone keeps channels 0 to 3 in both data and
        instrument.
        """
        s = make_signal(max_channel=3)
        np.testing.assert_array_equal(s.data.channels, [0, 1, 2, 3])
        np.testing.assert_array_equal(s.instrument.channels, [0, 1, 2, 3])

    def test_original_data_not_trimmed(self):
        """``original_data`` keeps all channels after trimming."""
        s = make_signal(min_channel=2, max_channel=4)
        np.testing.assert_array_equal(s.original_data.channels,
                                      np.arange(N_CHAN))
        assert s.original_data.counts.shape == (N_CHAN, N_PHASE)

    def test_support_of_original_instrument_size_is_sliced(self):
        """A support sized to the untrimmed instrument is cut to the kept
        channels. Row i is ``[i, 10 + i]``, so rows 2 to 4 must remain.
        """
        support = np.column_stack([np.arange(N_CHAN, dtype=float),
                                   np.arange(N_CHAN, dtype=float) + 10.0])
        s = make_signal(min_channel=2, max_channel=4, support=support)
        np.testing.assert_array_equal(s._support, support[2:5])


# --------------------------------------------------------------------------
# Quality filtering (data._quality, as set by Data.from_pha(checkqual=True))
# --------------------------------------------------------------------------
class TestQualityFiltering:
    """Instrument channels 0..5, channel 1 flagged bad (quality 1).

    The data keeps the 5 good channels, renumbered 0..4 as ``from_pha``
    does; ``_quality`` is indexed by channel number.
    """

    QUALITY = np.array([0, 1, 0, 0, 0, 0])
    CHANNEL_EDGES = np.array([np.arange(N_CHAN, dtype=float),
                              np.arange(N_CHAN, dtype=float) + 1.0])

    def _data(self, quality=True, grpdata=False):
        d = make_data(counts=np.ones((N_CHAN - 1, N_PHASE)),
                      channels=np.arange(N_CHAN - 1))
        if quality:
            d._quality = self.QUALITY.copy()
        d._grpdata = grpdata
        return d

    def _instrument(self, matrix=None):
        matrix = default_matrix() if matrix is None else matrix
        return Instrument(matrix, np.array([EDGES[:-1], EDGES[1:]]),
                          np.arange(N_CHAN),
                          channel_edges=self.CHANNEL_EDGES.copy(),
                          twod_edges=True)

    def test_quality_requires_2d_edges(self):
        """Quality flags with 1-D instrument energy edges raise EdgesError."""
        with pytest.raises(EdgesError, match='good quality'):
            make_signal(data=self._data())

    def test_grouping_requires_2d_edges(self):
        """Grouped data with 1-D instrument energy edges raise EdgesError."""
        with pytest.raises(EdgesError, match='group channels'):
            make_signal(data=self._data(quality=False, grpdata=True))

    def test_bad_channel_row_removed(self):
        """The response row of the bad channel is removed."""
        s = make_signal(data=self._data(), instrument=self._instrument())
        expected = np.delete(default_matrix(), 1, axis=0)
        np.testing.assert_array_equal(s.instrument.matrix, expected)

    def test_instrument_takes_data_channels(self):
        """After filtering, the instrument channels are the data's renumbered
        channels 0 to 4.
        """
        s = make_signal(data=self._data(), instrument=self._instrument())
        np.testing.assert_array_equal(s.instrument.channels, np.arange(5))

    def test_channel_edges_keep_original_channels(self):
        """The channel energy edges are those of the good original channels;
        channel 1 is gone.
        """
        s = make_signal(data=self._data(), instrument=self._instrument())
        np.testing.assert_array_equal(s.instrument.channel_edges[0],
                                      [0, 2, 3, 4, 5])
        np.testing.assert_array_equal(s.instrument.channel_edges[1],
                                      [1, 3, 4, 5, 6])

    def test_empty_input_column_dropped(self):
        """An input energy interval that only reaches the bad channel is
        dropped, with its energy edges.
        """
        m = default_matrix()
        m[:, 0] = 0.0
        m[1, 0] = 1.0
        s = make_signal(data=self._data(), instrument=self._instrument(m))
        assert s.instrument.matrix.shape == (N_CHAN - 1, N_IN - 1)
        np.testing.assert_allclose(s.instrument.energy_edges,
                                   [EDGES[1:-1], EDGES[2:]])


# --------------------------------------------------------------------------
# Interstellar integration in the constructor
# --------------------------------------------------------------------------
class TestInterstellarCheck:

    def test_valid_interstellar_stored(self):
        """An attenuation below 1 is accepted and the model is stored."""
        ism = FakeInterstellar(0.5)
        assert make_signal(interstellar=ism).interstellar is ism

    def test_attenuation_of_one_accepted(self):
        """An attenuation of exactly 1 (no absorption) is accepted."""
        ism = FakeInterstellar(1.0)
        assert make_signal(interstellar=ism).interstellar is ism

    def test_attenuation_above_one_raises(self):
        """An attenuation above 1 would add flux, so it raises ValueError."""
        with pytest.raises(ValueError, match='lower than 1'):
            make_signal(interstellar=FakeInterstellar(1.5))

    def test_unset_parameter_set_to_bounds_midpoint(self):
        """An interstellar parameter without a value is set to its bounds'
        midpoint for the check.
        """
        p = Parameter('dummy', strict_bounds=(0.0, 10.0), bounds=(2.0, 4.0),
                      doc='dummy', symbol=r'$d$', value=None)
        ism = FakeInterstellar(0.5, p)
        make_signal(interstellar=ism)
        assert ism['dummy'] == pytest.approx(3.0)

    def test_interstellar_parameters_in_signal_subspace(self):
        """Interstellar parameters become parameters of the Signal."""
        p = Parameter('dummy', strict_bounds=(0.0, 10.0), bounds=(2.0, 4.0),
                      doc='dummy', symbol=r'$d$', value=None)
        s = make_signal(interstellar=FakeInterstellar(0.5, p))
        assert 'dummy' in s.names

    def test_twod_edges_supported(self):
        """The attenuation check also works with 2-D instrument energy edges.
        """
        s = make_signal(instrument=make_instrument(twod_edges=True),
                        interstellar=FakeInterstellar(0.5))
        assert s.interstellar is not None

    def test_attenuation_above_one_raises_twod_edges(self):
        """With 2-D energy edges, an attenuation above 1 also raises
        ValueError.
        """
        with pytest.raises(ValueError, match='lower than 1'):
            make_signal(instrument=make_instrument(twod_edges=True),
                        interstellar=FakeInterstellar(1.5))


# --------------------------------------------------------------------------
# Photosphere prefix
# --------------------------------------------------------------------------
class TestPhotosphere:

    def test_prefix_stored(self):
        """``photosphere_prefix`` is returned by the ``photosphere`` property.
        """
        assert make_signal(photosphere_prefix='ST').photosphere == 'ST'

    def test_missing_prefix_raises_attribute_error(self, sig):
        """Without a prefix, reading ``photosphere`` raises AttributeError."""
        with pytest.raises(AttributeError):
            sig.photosphere


# --------------------------------------------------------------------------
# Waveband identification
# --------------------------------------------------------------------------
class TestWaveband:

    def test_full_matrix_input_range(self, sig):
        """With a strictly positive response, all input energy intervals are
        used.
        """
        assert sig._input_interval_range == (0, N_IN)

    def test_full_matrix_energy_edges(self, sig):
        """With a strictly positive response, the energy edges are the
        instrument's.
        """
        np.testing.assert_allclose(sig.energy_edges, EDGES)

    def test_energy_mids(self, sig):
        """Energy midpoints are the means of consecutive 1-D edges."""
        np.testing.assert_allclose(sig._energy_mids,
                                   0.5 * (EDGES[:-1] + EDGES[1:]))

    def test_instrument_channel_range(self, sig):
        """With data on all channels, all instrument rows are used."""
        assert sig._instrument_index_range_channels == (0, N_CHAN)

    def test_zero_corners_trim_waveband(self):
        """A zero response in the matrix corners (first channel and first input
        interval, last channel and last interval) removes those two input
        intervals from the waveband.
        """
        m = default_matrix()
        m[0, 0] = 0.0
        m[-1, -1] = 0.0
        s = make_signal(instrument=make_instrument(matrix=m))
        assert s._input_interval_range == (1, N_IN - 1)
        np.testing.assert_allclose(s.energy_edges, EDGES[1:N_IN])
        assert len(s._energy_mids) == len(s.energy_edges) - 1

    def test_twod_edges_shape(self):
        """With 2-D edges, the edges keep their two rows and there is one
        midpoint per interval.
        """
        s = make_signal(instrument=make_instrument(twod_edges=True))
        assert s.energy_edges.shape == (2, N_IN)
        assert s._energy_mids.shape == (N_IN,)

    def test_twod_mids_values(self):
        """With 2-D edges, midpoints are the means of each interval's lower and
        upper edges.
        """
        s = make_signal(instrument=make_instrument(twod_edges=True))
        np.testing.assert_allclose(s._energy_mids,
                                   0.5 * (EDGES[:-1] + EDGES[1:]))


# --------------------------------------------------------------------------
# original_data / data accessors
# --------------------------------------------------------------------------
class TestDataAccessors:

    def test_original_data_is_a_copy(self, sig):
        """``original_data`` is a separate object from ``data``."""
        assert sig.original_data is not sig.data

    def test_original_data_unaffected_by_changes(self, sig):
        """Changing ``data`` counts leaves ``original_data`` unchanged (deep
        copy).
        """
        sig.data.counts[0, 0] = 999.0
        assert sig.original_data.counts[0, 0] == 1.0

    def test_data_setter_accepts_data(self, sig):
        """A ``Data`` instance can replace the data."""
        new = make_data()
        sig.data = new
        assert sig.data is new

    @pytest.mark.parametrize('bad', ['x', 1, None, np.ones(3)])
    def test_data_setter_rejects_other_types(self, sig, bad):
        """Setting ``data`` to anything else raises TypeError."""
        with pytest.raises(TypeError, match='invalid type'):
            sig.data = bad


# --------------------------------------------------------------------------
# cache / store flags
# --------------------------------------------------------------------------
class TestCacheFlag:

    def test_default_off(self, sig):
        """Caching is off by default."""
        assert sig.cache is False

    def test_constructor_cache_true(self):
        """``cache=True`` in the constructor turns caching on."""
        assert make_signal(cache=True).cache is True

    @pytest.mark.parametrize('bad', [1, 0, 'yes', None])
    def test_constructor_non_bool_raises(self, bad):
        """A non-boolean ``cache`` raises TypeError, even 0 and 1."""
        with pytest.raises(TypeError, match='boolean'):
            make_signal(cache=bad)

    # deprecated ``store`` keyword, read from **kwargs
    def test_constructor_store_keyword(self):
        """The deprecated ``store`` keyword turns caching on."""
        assert make_signal(store=True).cache is True

    def test_constructor_store_overrides_cache(self):
        """``store`` takes precedence over ``cache`` when both are given."""
        assert make_signal(cache=False, store=True).cache is True

    @pytest.mark.parametrize('bad', [1, 'yes', None])
    def test_constructor_store_non_bool_raises(self, bad):
        """A non-boolean ``store`` raises TypeError."""
        with pytest.raises(TypeError, match='boolean'):
            make_signal(store=bad)

    def test_cache_setter(self, sig):
        """The ``cache`` setter turns caching on."""
        sig.cache = True
        assert sig.cache is True

    def test_store_alias_reads_cache(self, sig):
        """``store`` reads the ``cache`` flag."""
        sig.cache = True
        assert sig.store is True

    def test_store_setter_updates_cache(self, sig):
        """Setting ``store`` sets ``cache``."""
        sig.store = True
        assert sig.cache is True

    @pytest.mark.parametrize('bad', [1, 'True', None])
    def test_cache_setter_rejects_non_bool(self, sig, bad):
        """Setting ``cache`` to a non-boolean raises ValueError."""
        with pytest.raises(ValueError, match='boolean activation'):
            sig.cache = bad

    @pytest.mark.parametrize('bad', [1, 'True', None])
    def test_store_setter_rejects_non_bool(self, sig, bad):
        """Setting ``store`` to a non-boolean raises ValueError."""
        with pytest.raises(ValueError, match='boolean activation'):
            sig.store = bad


# --------------------------------------------------------------------------
# Simple stored properties
# --------------------------------------------------------------------------
class TestStoredProperties:

    # energies
    def test_energies_roundtrip(self, sig, energies):
        """``energies`` returns the array it was set to."""
        sig.energies = energies
        assert sig.energies is energies

    def test_energies_delete(self, sig, energies):
        """After ``del energies``, reading it raises AttributeError."""
        sig.energies = energies
        del sig.energies
        with pytest.raises(AttributeError):
            sig.energies

    def test_signals_missing_raises(self, sig):
        """Reading ``signals`` before any is registered raises AttributeError.
        """
        with pytest.raises(AttributeError):
            sig.signals

    def test_signals_append(self, sig):
        """Each assignment to ``signals`` appends a component."""
        sig.signals = np.ones((N_CHAN, N_PHASE))
        sig.signals = 2 * np.ones((N_CHAN, N_PHASE))
        assert len(sig.signals) == 2
        assert sig.num_components == 2

    def test_signals_returns_tuple_of_copies(self, sig):
        """``signals`` returns a tuple of copies: changing them leaves the
        stored ones unchanged.
        """
        sig.signals = np.ones((N_CHAN, N_PHASE))
        got = sig.signals
        assert isinstance(got, tuple)
        got[0][...] = -1.0
        assert np.all(sig.signals[0] == 1.0)

    def test_signals_delete(self, sig):
        """After ``del signals``, reading it raises AttributeError."""
        sig.signals = np.ones((N_CHAN, N_PHASE))
        del sig.signals
        with pytest.raises(AttributeError):
            sig.signals

    def test_incident_specific_flux_append_and_copy(self, sig):
        """Assignments to ``incident_specific_flux_signals`` append, and reads
        return copies.
        """
        sig.incident_specific_flux_signals = np.ones(4)
        sig.incident_specific_flux_signals = np.zeros(4)
        got = sig.incident_specific_flux_signals
        assert len(got) == 2
        got[0][...] = 5.0
        assert np.all(sig.incident_specific_flux_signals[0] == 1.0)

    def test_incident_specific_flux_delete(self, sig):
        """After deleting ``incident_specific_flux_signals``, reading it raises
        AttributeError.
        """
        sig.incident_specific_flux_signals = np.ones(4)
        del sig.incident_specific_flux_signals
        with pytest.raises(AttributeError):
            sig.incident_specific_flux_signals

    # incident flux
    def test_incident_flux_append_and_copy(self, sig):
        """``incident_flux_signals`` returns copies: changing them leaves the
        stored ones unchanged.
        """
        sig.incident_flux_signals = np.ones(4)
        got = sig.incident_flux_signals
        got[0][...] = 5.0
        assert np.all(sig.incident_flux_signals[0] == 1.0)

    def test_incident_flux_delete(self, sig):
        """After deleting ``incident_flux_signals``, reading it raises
        AttributeError.
        """
        sig.incident_flux_signals = np.ones(4)
        del sig.incident_flux_signals
        with pytest.raises(AttributeError):
            sig.incident_flux_signals

    # expected counts
    def test_expected_counts_roundtrip(self, sig):
        """``expected_counts`` returns the array it was set to."""
        arr = np.ones((N_CHAN, N_PHASE))
        sig.expected_counts = arr
        assert sig.expected_counts is arr

    def test_expected_counts_delete(self, sig):
        """After ``del expected_counts``, reading it raises AttributeError."""
        sig.expected_counts = np.ones(3)
        del sig.expected_counts
        with pytest.raises(AttributeError):
            sig.expected_counts

    # phases
    def test_phases_wraps_non_list(self, sig):
        """A single phase array is wrapped in a list."""
        sig.phases = np.linspace(0, 1, 5)
        assert len(sig.phases) == 1

    def test_phases_list_kept(self, sig):
        """A list of phase arrays is stored as one entry per hot region."""
        sig.phases = [np.linspace(0, 1, 5), np.linspace(0, 1, 7)]
        assert len(sig.phases) == 2

    def test_phases_returns_copies(self, sig):
        """``phases`` returns copies: changing them leaves the stored arrays
        unchanged.
        """
        sig.phases = [np.linspace(0, 1, 5)]
        sig.phases[0][...] = -1.0
        assert np.all(sig.phases[0] >= 0.0)

    # background signal
    def test_background_signal_ndarray_stored(self, sig):
        """An ndarray assigned to ``background_signal`` is stored as given."""
        arr = np.ones(3)
        sig.background_signal = arr
        assert sig.background_signal is arr

    def test_background_signal_non_ndarray_ignored(self, sig):
        """A non-ndarray assigned to ``background_signal`` is silently ignored.
        """
        sig.background_signal = [1, 2, 3]
        with pytest.raises(AttributeError):
            sig.background_signal

    def test_background_signal_delete(self, sig):
        """After ``del background_signal``, reading it raises AttributeError.
        """
        sig.background_signal = np.ones(3)
        del sig.background_signal
        with pytest.raises(AttributeError):
            sig.background_signal


# --------------------------------------------------------------------------
# Shifts
# --------------------------------------------------------------------------
class TestShifts:

    def test_shifts_include_signal_phase_shift(self):
        """The stored shifts are the hot-region shifts plus the Signal's
        ``phase_shift``.
        """
        s = make_signal(values={'phase_shift': 0.25})
        s.phases = [np.linspace(0, 1, 5), np.linspace(0, 1, 5)]
        s.shifts = np.array([0.1, 0.2])
        np.testing.assert_allclose(s.shifts, [0.35, 0.45])

    def test_shifts_zero_phase_shift(self, sig):
        """With ``phase_shift`` at zero, the shifts are stored unchanged."""
        sig.phases = [np.linspace(0, 1, 5)]
        sig.shifts = np.array([0.3])
        np.testing.assert_allclose(sig.shifts, [0.3])

    def test_shifts_wrong_length_raises(self, sig):
        """A number of shifts different from the number of phase arrays raises
        TypeError.
        """
        sig.phases = [np.linspace(0, 1, 5)]
        with pytest.raises(TypeError, match='1D ndarray'):
            sig.shifts = np.array([0.1, 0.2])

    def test_shifts_non_ndarray_raises(self, sig):
        """Shifts given as a list instead of an ndarray raise TypeError."""
        sig.phases = [np.linspace(0, 1, 5)]
        with pytest.raises(TypeError, match='1D ndarray'):
            sig.shifts = [0.1]

    def test_shifts_delete(self, sig):
        """After ``del shifts``, reading it raises AttributeError."""
        sig.phases = [np.linspace(0, 1, 5)]
        sig.shifts = np.array([0.1])
        del sig.shifts
        with pytest.raises(AttributeError):
            sig.shifts


# --------------------------------------------------------------------------
# Caching targets
# --------------------------------------------------------------------------
class TestCachingTargets:

    def test_not_declared_raises(self, sig):
        """Reading ``caching_targets`` before any are declared raises
        AttributeError.
        """
        with pytest.raises(AttributeError):
            sig.caching_targets

    def test_set_and_get(self, sig):
        """Declared targets are returned as a dict from attribute name to its
        value.
        """
        sig.caching_targets = ['data', 'instrument']
        got = sig.caching_targets
        assert got == {'data': sig.data, 'instrument': sig.instrument}

    def test_target_names(self, sig):
        """``caching_target_names`` lists the declared attribute names."""
        sig.caching_targets = ['data']
        assert sig.caching_target_names == ['data']

    @pytest.mark.parametrize('bad', ['data', ('data',), [1, 2], ['data', 3],
                                     ['no_such_attribute'], None])
    def test_invalid_targets_raise(self, sig, bad):
        """Targets that are not a list of existing attribute names raise
        ValueError.
        """
        with pytest.raises(ValueError, match='Invalid caching targets'):
            sig.caching_targets = bad

    def test_custom_cache_instructions_not_implemented(self, sig):
        """The base class has no custom cache instructions: calling them raises
        NotImplementedError.
        """
        with pytest.raises(NotImplementedError, match='Cache method'):
            sig.execute_custom_cache_instructions()


# --------------------------------------------------------------------------
# Log-likelihood
# --------------------------------------------------------------------------
class TestLogLikelihood:

    def test_finite_stored(self, sig):
        """A finite log-likelihood is stored as given."""
        sig.loglikelihood = -123.4
        assert sig.loglikelihood == -123.4

    def test_zero_stored(self, sig):
        """A log-likelihood of zero is stored (not treated as unset)."""
        sig.loglikelihood = 0.0
        assert sig.loglikelihood == 0.0

    def test_nan_raises(self, sig):
        """A NaN log-likelihood raises LikelihoodError."""
        with pytest.raises(LikelihoodError, match='NaN'):
            sig.loglikelihood = np.nan

    @pytest.mark.parametrize('bad', [np.inf, -np.inf])
    def test_infinite_becomes_minus_inf(self, sig, bad):
        """An infinite log-likelihood, of either sign, is stored as minus
        infinity.
        """
        sig.loglikelihood = bad
        assert sig.loglikelihood == -np.inf

    def test_unset_raises_attribute_error(self, sig):
        """Reading the log-likelihood before it is set raises AttributeError.
        """
        with pytest.raises(AttributeError):
            sig.loglikelihood

    def test_likelihood_error_is_exception(self):
        """LikelihoodError can be caught as an Exception."""
        assert issubclass(LikelihoodError, Exception)


# --------------------------------------------------------------------------
# bolometric chi2
# --------------------------------------------------------------------------
class TestBolometricChi2:

    def test_zero_for_perfect_model(self):
        """Expected counts equal to the data give a chi2 of zero."""
        counts = np.arange(1, N_CHAN * N_PHASE + 1, dtype=float).reshape(
            N_CHAN, N_PHASE)
        s = make_signal(data=make_data(counts=counts))
        s.expected_counts = counts.copy()
        assert s.bolometric_chi2 == pytest.approx(0.0)

    def test_matches_manual_formula(self):
        """The chi2 equals the sum over phases of (data - model)^2 / model,
        with counts summed over channels.
        """
        rng = np.random.default_rng(1)
        counts = rng.integers(5, 50, size=(N_CHAN, N_PHASE)).astype(float)
        model = counts + rng.normal(0, 2, size=counts.shape)
        s = make_signal(data=make_data(counts=counts))
        s.expected_counts = model
        pd, pm = counts.sum(axis=0), model.sum(axis=0)
        assert s.bolometric_chi2 == pytest.approx(np.sum((pd - pm) ** 2 / pm))

    def test_non_negative(self):
        """The chi2 is not negative for a model that differs from the data."""
        s = make_signal(data=make_data(counts=np.full((N_CHAN, N_PHASE), 4.0)))
        s.expected_counts = np.full((N_CHAN, N_PHASE), 6.0)
        assert s.bolometric_chi2 >= 0.0

    def test_missing_expected_counts_raises(self, sig):
        """Without expected counts, the chi2 raises AttributeError."""
        with pytest.raises(AttributeError):
            sig.bolometric_chi2


# --------------------------------------------------------------------------
# create_energy_array
# --------------------------------------------------------------------------
class TestCreateEnergyArray:

    def test_default_length(self, sig, energies):
        """By default, the array has 10 times as many energies as ``energies``.
        """
        sig.energies = energies
        assert len(sig.create_energy_array()) == int(10.0 * len(energies))

    @pytest.mark.parametrize('rel', [0.5, 1.0, 2.0, 7.5])
    def test_length_scales_with_rel_num(self, sig, energies, rel):
        """The number of energies is ``rel_num`` times the length of
        ``energies``.
        """
        sig.energies = energies
        assert len(sig.create_energy_array(rel)) == int(rel * len(energies))

    def test_span_equals_edges_1d(self, sig, energies):
        """With 1-D edges, the array runs from the first to the last energy
        edge.
        """
        sig.energies = energies
        e = sig.create_energy_array()
        assert e[0] == pytest.approx(EDGES[0])
        assert e[-1] == pytest.approx(EDGES[-1])

    def test_span_equals_edges_2d(self, energies):
        """With 2-D edges, the array runs from the lowest to the highest energy
        edge.
        """
        s = make_signal(instrument=make_instrument(twod_edges=True))
        s.energies = energies
        e = s.create_energy_array()
        assert e[0] == pytest.approx(EDGES[0])
        assert e[-1] == pytest.approx(EDGES[-1])

    def test_strictly_increasing(self, sig, energies):
        """Energies are strictly increasing."""
        sig.energies = energies
        assert np.all(np.diff(sig.create_energy_array()) > 0)

    def test_logarithmic_spacing(self, sig, energies):
        """Energies are evenly spaced in log10."""
        sig.energies = energies
        ratios = np.diff(np.log10(sig.create_energy_array()))
        np.testing.assert_allclose(ratios, ratios[0], rtol=1e-8)

    def test_requires_energies(self, sig):
        """Without ``energies`` set, AttributeError is raised."""
        with pytest.raises(AttributeError):
            sig.create_energy_array()


# --------------------------------------------------------------------------
# register (uses compiled xpsi integrators)
# --------------------------------------------------------------------------
class TestRegister:

    def test_single_region_single_component(self, sig, energies):
        """One hot region gives one component of shape (channels, phases)."""
        sig.energies = energies
        sig.register([[flux(energies)]])
        assert sig.num_components == 1
        assert sig.signals[0].shape == (N_CHAN, N_PHASE)

    def test_unit_flux_gives_matrix_times_bin_widths(self, sig, energies):
        """A flux of 1 per keV integrates to the bin width over each input
        interval, so each phase column is the response times the widths.
        """
        sig.energies = energies
        sig.register([[flux(energies, 1.0)]])
        expected = default_matrix() @ np.diff(EDGES)
        for k in range(N_PHASE):
            np.testing.assert_allclose(sig.signals[0][:, k], expected,
                                       rtol=1e-5)

    def test_two_regions_give_two_components(self, sig, energies):
        """Two hot regions give two components."""
        sig.energies = energies
        sig.register([[flux(energies)], [flux(energies)]])
        assert sig.num_components == 2

    def test_reset_true_clears_previous(self, sig, energies):
        """With ``reset=True``, a new registration replaces the previous
        components.
        """
        sig.energies = energies
        sig.register([[flux(energies)]])
        sig.register([[flux(energies)]], reset=True)
        assert sig.num_components == 1

    def test_reset_false_accumulates(self, sig, energies):
        """With ``reset=False``, a new registration adds to the previous
        components.
        """
        sig.energies = energies
        sig.register([[flux(energies)]])
        sig.register([[flux(energies)]], reset=False)
        assert sig.num_components == 2

    def test_components_in_a_region_are_summed(self, sig, energies):
        """Components within a hot region are summed: two fluxes of 1 equal one
        flux of 2.
        """
        sig.energies = energies
        sig.register([[flux(energies, 1.0), flux(energies, 1.0)]])
        two = sig.signals[0]
        sig.register([[flux(energies, 2.0)]])
        np.testing.assert_allclose(two, sig.signals[0], rtol=1e-10)

    def test_signal_scales_linearly_with_flux(self, sig, energies):
        """Tripling the flux triples the registered signal."""
        sig.energies = energies
        sig.register([[flux(energies, 1.0)]])
        one = sig.signals[0]
        sig.register([[flux(energies, 3.0)]])
        np.testing.assert_allclose(sig.signals[0], 3.0 * one, rtol=1e-10)

    def test_threads_give_same_result(self, sig, energies):
        """Two threads give the same signal as one."""
        sig.energies = energies
        sig.register([[flux(energies)]], threads=1)
        one = sig.signals[0]
        sig.register([[flux(energies)]], threads=2)
        assert sig.num_components == 1
        np.testing.assert_allclose(sig.signals[0], one, rtol=1e-12)

    def test_interstellar_applied(self, energies):
        """A constant attenuation of 0.5 halves the registered signal."""
        plain = make_signal()
        plain.energies = energies
        plain.register([[flux(energies)]])

        attenuated = make_signal(interstellar=FakeInterstellar(0.5))
        attenuated.energies = energies
        attenuated.register([[flux(energies)]])

        np.testing.assert_allclose(attenuated.signals[0],
                                   0.5 * plain.signals[0], rtol=1e-10)

    def test_no_caching_by_default(self, sig, energies):
        """Without caching, no incident signals are stored."""
        sig.energies = energies
        sig.register([[flux(energies)]])
        with pytest.raises(AttributeError):
            sig.incident_flux_signals
        with pytest.raises(AttributeError):
            sig.incident_specific_flux_signals

    def test_caching_stores_incident_signals(self, energies):
        """With caching, the summed incident specific flux and the integrated
        flux of each region are stored.
        """
        s = make_signal(cache=True)
        s.energies = energies
        s.register([[flux(energies, 1.0), flux(energies, 2.0)]])
        specific = s.incident_specific_flux_signals
        assert len(specific) == 1
        np.testing.assert_allclose(specific[0], flux(energies, 3.0))
        assert len(s.incident_flux_signals) == 1

    def test_caching_reset_does_not_accumulate(self, energies):
        """With caching, a second registration replaces the cached incident
        signals.
        """
        s = make_signal(cache=True)
        s.energies = energies
        s.register([[flux(energies)]])
        s.register([[flux(energies)]])
        assert len(s.incident_specific_flux_signals) == 1
        assert len(s.incident_flux_signals) == 1

    def test_cached_incident_flux_not_attenuated(self, energies):
        """The cached incident flux is stored before interstellar attenuation.
        """
        s = make_signal(cache=True, interstellar=FakeInterstellar(0.5))
        s.energies = energies
        s.register([[flux(energies)]])
        plain = make_signal(cache=True)
        plain.energies = energies
        plain.register([[flux(energies)]])
        np.testing.assert_allclose(s.incident_flux_signals[0],
                                   plain.incident_flux_signals[0])

    def test_input_flux_not_modified(self, sig, energies):
        """``register`` leaves the input flux arrays unchanged."""
        sig.energies = energies
        comp = flux(energies, 2.0)
        sig.register([[comp]])
        assert np.all(comp == 2.0)

    def test_twod_edges_unit_flux(self, energies):
        """With 2-D edges, a unit flux gives the same response times bin widths
        as with 1-D edges, through the 2-D integrator.
        """
        s = make_signal(instrument=make_instrument(twod_edges=True))
        s.energies = energies
        s.register([[flux(energies, 1.0)]])
        expected = default_matrix() @ np.diff(EDGES)
        for k in range(N_PHASE):
            np.testing.assert_allclose(s.signals[0][:, k], expected,
                                       rtol=1e-5)

    def test_background_registered(self, energies):
        """The background is folded through the response: one count per input
        interval gives the row sums of the matrix.
        """
        bkg = FakeBackground()
        s = make_signal(background=bkg)
        s.energies = energies
        s.register([[flux(energies)]])
        assert bkg.registered_background.shape == (N_CHAN, N_PHASE)
        expected = default_matrix() @ np.ones(N_IN)
        for k in range(N_PHASE):
            np.testing.assert_allclose(bkg.registered_background[:, k],
                                       expected)

    def test_background_type_error_reraised(self, energies, capsys):
        """A TypeError raised by the background is re-raised after a message is
        printed.
        """
        class BrokenBackground(Background):
            def __call__(self, energy_edges, phases):
                raise TypeError('broken')

        s = make_signal(background=BrokenBackground())
        s.energies = energies
        with pytest.raises(TypeError, match='broken'):
            s.register([[flux(energies)]])
        assert 'incident background' in capsys.readouterr().out

    @pytest.mark.parametrize('cache, calls', [(True, 1), (False, 0)])
    def test_custom_cache_instructions_called_when_caching(self, energies,
                                                           cache, calls):
        """Custom cache instructions run once per registration with caching,
        and not at all without.
        """
        class CountingSignal(ConcreteSignal):
            calls = 0

            def execute_custom_cache_instructions(self):
                CountingSignal.calls += 1

        s = CountingSignal(make_data(), make_instrument(), cache=cache)
        s.energies = energies
        s.register([[flux(energies)]])
        assert CountingSignal.calls == calls


# --------------------------------------------------------------------------
# synthesise: argument validation only (no compiled synthesis is run)
# --------------------------------------------------------------------------
class TestSynthesiseValidation:

    def test_creates_directory(self, sig, tmp_path):
        """The output directory is created before the arguments are checked."""
        target = tmp_path / 'newdir'
        with pytest.raises(ValueError,
                           match='exposure_time or expected_source_counts'):
            sig.synthesise(directory=str(target))
        assert target.is_dir()

    def test_fits_requires_instrument_name(self, sig, tmp_path):
        """FITS output without an instrument name raises AssertionError."""
        with pytest.raises(AssertionError, match='Instrument name'):
            sig.synthesise(exposure_time=1.0, format='FITS',
                           directory=str(tmp_path))

    def test_requires_exposure_or_counts(self, sig, tmp_path):
        """Neither ``exposure_time`` nor ``expected_source_counts`` raises
        ValueError.
        """
        with pytest.raises(ValueError,
                           match='exposure_time or expected_source_counts'):
            sig.synthesise(directory=str(tmp_path))

    def test_base_synthesise_not_the_stub(self, sig):
        """``synthesise`` is defined twice in Signal.py; the keyword-based one
        must be the active one.
        """
        import inspect
        params = inspect.signature(sig.synthesise).parameters
        assert 'exposure_time' in params


# --------------------------------------------------------------------------
# synthesise: end to end (uses the compiled synthesis)
# --------------------------------------------------------------------------
N_HOT = 32  # GSL aborts the process if the hot-region phase grid is too short


@pytest.fixture
def ready(energies):
    """Factory of signals ready to synthesise: one flat hot region."""
    def _make(**kw):
        s = make_signal(**kw)
        s.energies = energies
        s.register([[np.ones((len(energies), N_HOT))]])
        s.phases = [np.linspace(0.0, 1.0, N_HOT)]
        s.shifts = np.array([0.0])
        return s
    return _make


class TestSynthesise:

    T = 10.0
    BKG_1D = 3.0 * np.ones(N_CHAN)

    def test_txt_files_written(self, ready, tmp_path):
        """TXT output writes the expected counts and their Poisson realisation.
        """
        ready().synthesise(exposure_time=self.T, directory=str(tmp_path))
        assert (tmp_path / 'synthetic.dat').is_file()
        assert (tmp_path / 'synthetic_realisation.dat').is_file()

    def test_expected_counts_total(self, ready, tmp_path):
        """The total expected counts are the exposure time times the phase-
        averaged signal (the signal is flat in phase).
        """
        s = ready()
        s.synthesise(exposure_time=self.T, directory=str(tmp_path))
        expected = self.T * s.signals[0].sum(axis=0).mean()
        assert s.expected_counts.sum() == pytest.approx(expected, rel=1e-6)

    def test_expected_counts_file_matches_attribute(self, ready, tmp_path):
        """The written expected counts are those stored in ``expected_counts``.
        """
        s = ready()
        s.synthesise(exposure_time=self.T, directory=str(tmp_path))
        np.testing.assert_allclose(np.loadtxt(tmp_path / 'synthetic.dat'),
                                   s.expected_counts, rtol=1e-6)

    def test_seed_reproducible(self, ready, tmp_path):
        """The same seed gives the same realisation."""
        for d in ('a', 'b'):
            ready().synthesise(exposure_time=self.T, seed=42,
                               directory=str(tmp_path / d))
        np.testing.assert_array_equal(
            np.loadtxt(tmp_path / 'a' / 'synthetic_realisation.dat'),
            np.loadtxt(tmp_path / 'b' / 'synthetic_realisation.dat'))

    def test_expected_source_counts(self, ready, tmp_path):
        """``expected_source_counts`` sets the total expected counts."""
        s = ready()
        s.synthesise(expected_source_counts=500, directory=str(tmp_path))
        assert s.expected_counts.sum() == pytest.approx(500, rel=1e-6)

    def test_data_bkg_1d_saved(self, ready, tmp_path):
        """A phase-averaged ``data_BKG`` is saved as given, with its
        realisation.
        """
        ready().synthesise(exposure_time=self.T, data_BKG=self.BKG_1D,
                           save_background=True, directory=str(tmp_path))
        assert (tmp_path / 'synthetic_bkg_realisation.dat').is_file()
        np.testing.assert_allclose(np.loadtxt(tmp_path / 'synthetic_bkg.dat'),
                                   self.BKG_1D)

    def test_background_and_data_bkg_raise(self, ready, tmp_path):
        """A background model and ``data_BKG`` together raise ValueError."""
        s = ready(background=FakeBackground())
        with pytest.raises(ValueError, match='both data_BKG'):
            s.synthesise(exposure_time=self.T, data_BKG=self.BKG_1D,
                         directory=str(tmp_path))

    def test_background_requires_expected_counts(self, ready, tmp_path):
        """A background model without ``expected_background_counts`` raises
        ValueError.
        """
        s = ready(background=FakeBackground())
        with pytest.raises(ValueError, match='expected_background_counts'):
            s.synthesise(exposure_time=self.T, directory=str(tmp_path))

    def test_background_three_phase_bins_raise(self, ready, tmp_path):
        """A background model with data of 3 phase bins (4 edges, the default
        data) raises ValueError.
        """
        s = ready(background=FakeBackground())
        with pytest.raises(ValueError, match='3 or 4 phases'):
            s.synthesise(exposure_time=self.T, expected_background_counts=10,
                         directory=str(tmp_path))

    def test_warning_save_without_background(self, ready, tmp_path,
                                             monkeypatch, capsys):
        """``save_background=True`` without a background warns and writes no
        background file.
        """
        monkeypatch.setattr('xpsi._verbose', True)
        ready().synthesise(exposure_time=self.T, save_background=True,
                           directory=str(tmp_path))
        assert 'no background file will be saved' in capsys.readouterr().out
        assert not (tmp_path / 'synthetic_bkg.dat').exists()

    def test_warning_background_not_saved(self, ready, tmp_path,
                                          monkeypatch, capsys):
        """A background used but not saved gives a warning, and no background
        file.
        """
        monkeypatch.setattr('xpsi._verbose', True)
        ready().synthesise(exposure_time=self.T, data_BKG=self.BKG_1D,
                           directory=str(tmp_path))
        assert 'will not be saved' in capsys.readouterr().out
        assert not (tmp_path / 'synthetic_bkg.dat').exists()

    def test_fits_phase_resolved_writes_evt(self, ready, tmp_path):
        """FITS output of phase-resolved data writes event files."""
        ready().synthesise(exposure_time=self.T, format='FITS',
                           instrument_name='X', directory=str(tmp_path))
        assert (tmp_path / 'synthetic.evt').is_file()
        assert (tmp_path / 'synthetic_realisation.evt').is_file()

    def test_fits_phase_averaged_writes_pha(self, ready, tmp_path):
        """FITS output of data with one phase bin writes spectra, background
        included.
        """
        data = make_data(counts=np.ones((N_CHAN, 1)),
                         phases=np.array([0.0, 1.0]))
        ready(data=data).synthesise(exposure_time=self.T, format='FITS',
                                    instrument_name='X',
                                    data_BKG=self.BKG_1D, save_background=True,
                                    directory=str(tmp_path))
        for stem in ('synthetic', 'synthetic_realisation',
                     'synthetic_bkg', 'synthetic_bkg_realisation'):
            assert (tmp_path / (stem + '.pha')).is_file()


# --------------------------------------------------------------------------
# File writers
# --------------------------------------------------------------------------
class TestWriteTXT:

    def test_2d_counts_roundtrip(self, sig, tmp_path):
        """2-D counts are written as a table and read back unchanged."""
        counts = np.arange(12, dtype=float).reshape(4, 3)
        sig._write_TXT(counts, str(tmp_path / 'out'))
        loaded = np.loadtxt(tmp_path / 'out.dat', ndmin=2)
        np.testing.assert_allclose(loaded, counts)

    def test_1d_counts_become_single_row(self, sig, tmp_path):
        """1-D counts are written as one row."""
        sig._write_TXT(np.arange(5, dtype=float), str(tmp_path / 'out'))
        loaded = np.loadtxt(tmp_path / 'out.dat', ndmin=2)
        assert loaded.shape == (1, 5)

    def test_integer_format(self, sig, tmp_path):
        """``fmt='%u'`` writes the counts as integers."""
        sig._write_TXT(np.array([[1.0, 2.0]]), str(tmp_path / 'out'),
                       fmt='%u')
        assert (tmp_path / 'out.dat').read_text().split() == ['1', '2']

    def test_extra_kwargs_ignored(self, sig, tmp_path):
        """Keywords used by the FITS writers are accepted and ignored."""
        sig._write_TXT(np.ones((2, 2)), str(tmp_path / 'out'),
                       instrument='X', exposure=1.0)
        assert (tmp_path / 'out.dat').exists()


class TestWritePHA:

    KW = dict(instrument='FAKE', exposure=123.0, backscale=2.0)

    def _write(self, sig, tmp_path, counts):
        channels = np.arange(len(counts))
        sig._write_PHA(channels, counts, str(tmp_path / 'spec'), **self.KW)
        return fits.open(tmp_path / 'spec.pha')

    def test_file_created_with_extension(self, sig, tmp_path):
        """The spectrum is written to ``<filename>.pha``."""
        self._write(sig, tmp_path, np.arange(6, dtype=float))
        assert (tmp_path / 'spec.pha').exists()

    def test_columns_roundtrip(self, sig, tmp_path):
        """The CHANNEL and COUNTS columns hold the channels and counts given.
        """
        counts = np.arange(6, dtype=float) + 0.5
        with self._write(sig, tmp_path, counts) as hdul:
            tbl = hdul['SPECTRUM'].data
            np.testing.assert_array_equal(tbl['CHANNEL'], np.arange(6))
            np.testing.assert_allclose(tbl['COUNTS'], counts)

    def test_header_keywords(self, sig, tmp_path):
        """The header holds the instrument, the SPECTRUM class, the exposure
        and the backscale.
        """
        with self._write(sig, tmp_path, np.ones(6)) as hdul:
            h = hdul['SPECTRUM'].header
            assert h['INSTRUME'] == 'FAKE'
            assert h['HDUCLAS1'] == 'SPECTRUM'
            assert h['EXPOSURE'] == 123.0
            assert h['BACKSCAL'] == 2.0

    def test_overwrites_existing(self, sig, tmp_path):
        """Writing again to the same file replaces it."""
        self._write(sig, tmp_path, np.ones(6)).close()
        with self._write(sig, tmp_path, 2 * np.ones(6)) as hdul:
            np.testing.assert_allclose(hdul['SPECTRUM'].data['COUNTS'], 2.0)


class TestWriteEVT:

    COUNTS = np.array([[1, 0, 2, 0],
                       [0, 3, 0, 0],
                       [0, 0, 0, 1]], dtype=float)
    KW = dict(instrument='FAKE', exposure=10.0, backscale=1.0)

    def _write(self, sig, tmp_path):
        sig._write_EVT(np.array([10, 11, 12]), self.COUNTS,
                       str(tmp_path / 'ev'), **self.KW)
        return fits.open(tmp_path / 'ev.evt')

    def test_number_of_events(self, sig, tmp_path):
        """There is one event per count."""
        with self._write(sig, tmp_path) as hdul:
            assert len(hdul['EVENTS'].data) == int(self.COUNTS.sum())

    def test_channels_repeated_per_count(self, sig, tmp_path):
        """Each channel appears as many times as its total counts."""
        with self._write(sig, tmp_path) as hdul:
            pi = hdul['EVENTS'].data['PI']
            assert (pi == 10).sum() == 3
            assert (pi == 11).sum() == 3
            assert (pi == 12).sum() == 1

    def test_phases_are_bin_centres(self, sig, tmp_path):
        """Event phases are the centres of the phase bins."""
        with self._write(sig, tmp_path) as hdul:
            phases = np.unique(hdul['EVENTS'].data['PULSE_PHASE'])
            np.testing.assert_allclose(phases, [0.125, 0.375, 0.625, 0.875])

    @pytest.mark.parametrize('channel, phases', [
        (10, [0.125, 0.625, 0.625]),
        (11, [0.375, 0.375, 0.375]),
        (12, [0.875]),
    ])
    def test_phase_paired_with_channel(self, sig, tmp_path, channel, phases):
        """Each event carries the phase-bin centre of its own counts entry."""
        with self._write(sig, tmp_path) as hdul:
            ev = hdul['EVENTS'].data
            got = np.sort(ev['PULSE_PHASE'][ev['PI'] == channel])
            np.testing.assert_allclose(got, phases)

    def test_header(self, sig, tmp_path):
        """The header holds the EVENTS class, the exposure and the instrument.
        """
        with self._write(sig, tmp_path) as hdul:
            h = hdul['EVENTS'].header
            assert h['HDUCLAS1'] == 'EVENTS'
            assert h['EXPOSURE'] == 10.0
            assert h['INSTRUME'] == 'FAKE'


class TestRoundTrip:
    """Files written by ``Signal`` are read back by ``Data.load``.

    Only the event file for now: ``_write_PHA`` writes no TLMIN1/TLMAX1,
    which ``Data.load`` needs for a spectrum.
    """

    CHANNELS = np.array([10, 11, 12])

    @pytest.mark.parametrize('channels', [CHANNELS, None],
                             ids=['given', 'from_file'])
    def test_evt_roundtrip(self, sig, tmp_path, channels):
        """An event file read by ``Data.load`` gives back the counts, channels,
        phase edges, exposure and instrument, with channels given or read from
        the file.
        """
        sig._write_EVT(self.CHANNELS, TestWriteEVT.COUNTS,
                       str(tmp_path / 'ev'), **TestWriteEVT.KW)
        d = Data.load(str(tmp_path / 'ev.evt'),
                      n_phases=TestWriteEVT.COUNTS.shape[1],
                      channels=channels)
        np.testing.assert_array_equal(d.counts, TestWriteEVT.COUNTS)
        np.testing.assert_array_equal(d.channels, self.CHANNELS)
        np.testing.assert_allclose(d.phases, np.linspace(0.0, 1.0, 5))
        assert d.exposure_time == TestWriteEVT.KW['exposure']
        assert d.instrument == TestWriteEVT.KW['instrument']

class TestWriteDispatch:

    @pytest.fixture
    def recorded(self, sig, monkeypatch):
        calls = []
        for name in ('TXT', 'PHA', 'EVT'):
            monkeypatch.setattr(
                sig, '_write_' + name,
                lambda _n=name, **kw: calls.append((_n, kw)))
        return calls

    KW = dict(channels=np.arange(3), filename='f', instrument='X',
              exposure=1.0, backscale=1.0)

    def _check_forwarded(self, kw, counts):
        """All keyword arguments reach the writer unchanged."""
        assert set(kw) == set(self.KW) | {'counts'}
        np.testing.assert_array_equal(kw['channels'], self.KW['channels'])
        for key in ('filename', 'instrument', 'exposure', 'backscale'):
            assert kw[key] == self.KW[key]
        np.testing.assert_array_equal(kw['counts'], counts)

    def test_txt(self, sig, recorded):
        """TXT format goes to the TXT writer with all arguments unchanged."""
        counts = np.arange(6.0).reshape(3, 2)
        sig._write('TXT', counts=counts, **self.KW)
        assert [c[0] for c in recorded] == ['TXT']
        self._check_forwarded(recorded[0][1], counts)

    def test_fits_phase_resolved_goes_to_evt(self, sig, recorded):
        """FITS counts with several phase columns go to the event writer."""
        counts = np.arange(6.0).reshape(3, 2)
        sig._write('FITS', counts=counts, **self.KW)
        assert [c[0] for c in recorded] == ['EVT']
        self._check_forwarded(recorded[0][1], counts)

    def test_fits_single_phase_column_goes_to_pha_as_1d(self, sig, recorded):
        """FITS counts with one phase column go to the spectrum writer,
        flattened to 1-D.
        """
        counts = np.arange(3.0).reshape(3, 1)
        sig._write('FITS', counts=counts, **self.KW)
        assert [c[0] for c in recorded] == ['PHA']
        kw = recorded[0][1]
        assert kw['counts'].shape == (3,)
        self._check_forwarded(kw, counts[:, 0])

    def test_fits_1d_goes_to_pha(self, sig, recorded):
        """1-D FITS counts go to the spectrum writer."""
        counts = np.arange(3.0)
        sig._write('FITS', counts=counts, **self.KW)
        assert [c[0] for c in recorded] == ['PHA']
        self._check_forwarded(recorded[0][1], counts)

    @pytest.mark.parametrize('bad', ['fits', 'txt', 'CSV', '', None])
    def test_unknown_format_raises(self, sig, recorded, bad):
        """Any format other than TXT or FITS (case-sensitive) raises
        ValueError, and nothing is written.
        """
        with pytest.raises(ValueError, match='Unknown format'):
            sig._write(bad, counts=np.ones(3), **self.KW)
        assert recorded == []


# --------------------------------------------------------------------------
# construct_energy_array (module-level function)
# --------------------------------------------------------------------------
class StubSignal:
    """Minimal stand-in for a Signal. A plain class (identity equality) is
    required because construct_energy_array calls ``list.remove``, which
    compares with ``==``; SimpleNamespace would compare the numpy arrays."""
 
    def __init__(self, energy_edges):
        self.energy_edges = energy_edges
 
 
def stub(lo, hi, twod_edges=False):
    if twod_edges:
        return StubSignal(np.array([[lo], [hi]]))
    return StubSignal(np.array([lo, hi]))

class TestConstructEnergyArray:

    def test_single_signal_length(self):
        """With one signal, the array has the requested number of energies."""
        assert len(construct_energy_array(50, [stub(0.3, 10.0)])) == 50

    def test_single_signal_span(self):
        """With one signal, the array spans its energy edges."""
        e = construct_energy_array(50, [stub(0.3, 10.0)])
        assert e[0] == pytest.approx(0.3)
        assert e[-1] == pytest.approx(10.0)

    def test_output_ascending(self):
        """Energies are strictly increasing."""
        e = construct_energy_array(50, [stub(0.3, 10.0)])
        assert np.all(np.diff(e) > 0)

    def test_log_spaced(self):
        """Energies are evenly spaced in log10."""
        e = construct_energy_array(50, [stub(0.3, 10.0)])
        d = np.diff(np.log10(e))
        np.testing.assert_allclose(d, d[0], rtol=1e-8)

    def test_contiguous_output(self):
        """The array is C-contiguous, as the compiled integrators need."""
        assert construct_energy_array(20, [stub(0.3, 10.0)]).flags['C_CONTIGUOUS']

    def test_2d_edges_supported(self):
        """A signal with 2-D edges spans its lowest to highest edge."""
        e = construct_energy_array(30, [stub(0.3, 10.0, twod_edges=True)])
        assert e[0] == pytest.approx(0.3)
        assert e[-1] == pytest.approx(10.0)

    def test_overlapping_signals_cover_union(self):
        """Overlapping signals give an array spanning the union of their
        ranges.
        """
        e = construct_energy_array(60, [stub(0.3, 5.0), stub(1.0, 10.0)])
        assert e[0] == pytest.approx(0.3)
        assert e[-1] == pytest.approx(10.0)

    def test_order_of_signals_irrelevant(self):
        """The order of the signals does not change the result."""
        a = construct_energy_array(60, [stub(0.3, 5.0), stub(1.0, 10.0)])
        b = construct_energy_array(60, [stub(1.0, 10.0), stub(0.3, 5.0)])
        np.testing.assert_allclose(a, b)

    def test_nested_signal_does_not_extend_range(self):
        """A signal inside another one's range does not change the span."""
        e = construct_energy_array(40, [stub(0.3, 10.0), stub(1.0, 2.0)])
        assert e[0] == pytest.approx(0.3)
        assert e[-1] == pytest.approx(10.0)

    def test_input_list_is_consumed(self):
        """The input list of signals is emptied by the function."""
        signals = [stub(0.3, 5.0), stub(1.0, 10.0)]
        construct_energy_array(40, signals)
        assert signals == []

    def test_max_energy_truncates(self):
        """``max_energy`` inside the range becomes the upper limit."""
        e = construct_energy_array(40, [stub(0.3, 10.0)], max_energy=5.0)
        assert e[-1] == pytest.approx(5.0)
        assert e[0] == pytest.approx(0.3)

    def test_max_energy_above_range_ignored(self):
        """``max_energy`` above the range has no effect."""
        e = construct_energy_array(40, [stub(0.3, 10.0)], max_energy=50.0)
        assert e[-1] == pytest.approx(10.0)

    def test_coverage_gap_left_empty(self):
        """No energies fall in a gap between two signals' ranges."""
        e = construct_energy_array(80, [stub(0.3, 2.0), stub(2.5, 10.0)])
        assert not np.any((e > 2.0) & (e < 2.5))
        assert np.all(np.diff(e) > 0)
        assert e[0] == pytest.approx(0.3)
        assert e[-1] == pytest.approx(10.0)

    def test_mixed_1d_2d_signals(self):
        """Signals with 1-D and 2-D edges combine into one array spanning both.
        """
        e = construct_energy_array(60, [stub(0.3, 5.0),
                                        stub(1.0, 10.0, twod_edges=True)])
        assert len(e) == 60
        assert e[0] == pytest.approx(0.3)
        assert e[-1] == pytest.approx(10.0)

    def test_three_overlapping_signals(self):
        """Three overlapping signals, given out of order, give one increasing
        array over their union.
        """
        e = construct_energy_array(60, [stub(1.0, 3.0), stub(0.3, 2.0),
                                        stub(2.5, 10.0)])
        assert len(e) == 60
        assert e[0] == pytest.approx(0.3)
        assert e[-1] == pytest.approx(10.0)
        assert np.all(np.diff(e) > 0)

    def test_max_energy_inside_gap(self):
        """``max_energy`` inside the gap (2.0, 2.5) drops the upper waveband,
        so the array stops at the gap's lower edge.
        """
        e = construct_energy_array(80, [stub(0.3, 2.0), stub(2.5, 10.0)],
                                   max_energy=2.2)
        assert len(e) == 80
        assert e[0] == pytest.approx(0.3)
        assert e[-1] == pytest.approx(2.0)
        assert np.all(e <= 2.0 * (1 + 1e-12))
