import numpy as np
import pytest
from xpsi.Interstellar import Interstellar

class TestInterstellar:

    class MockInterstellar(Interstellar):
        def attenuation(self, energies):
            return np.exp(-energies)  # Mock implementation for testing

    def setup_method(self):
        self.model = self.MockInterstellar()

    def test_attenuation(self):
        energies = np.array([1., 2., 3.])
        expected = np.exp(-energies)
        result = self.model.attenuation(energies)
        np.testing.assert_allclose(result, expected, rtol=1e-5)

    def test_call_single_dimension(self):
        energies = np.array([1., 2., 3.])
        signal = np.array([10., 20., 30.])
        expected_signal = signal * np.exp(-energies)
        self.model(energies, signal)
        np.testing.assert_array_equal(signal, expected_signal)

    def test_call_two_dimensions(self):
        energies = np.array([1., 2., 3.])
        signal = np.array([[10., 15.], [20., 25.], [30., 35.]])
        expected_signal = signal * np.exp(-energies[:, np.newaxis])
        self.model(energies, signal)
        np.testing.assert_array_equal(signal, expected_signal)

    def test_call_invalid_signal_type_str(self):
        energies = np.array([1., 2., 3.])
        signal = 'invalid'
        with pytest.raises(TypeError, match='Signal must be a numpy.ndarray.'):
            self.model(energies, signal)

    def test_call_invalid_signal_type_list(self):
        energies = np.array([1., 2., 3.])
        signal = [10., 20., 30.]  # List instead of ndarray
        with pytest.raises(TypeError, match='Signal must be a numpy.ndarray.'):
            self.model(energies, signal)

    def test_call_invalid_signal_dimensions(self):
        energies = np.array([1., 2., 3.])
        signal = np.array([[10., 15.]])  # Invalid shape
        with pytest.raises(ValueError, match="non-broadcastable output operand with shape"): # with shape (1,) doesn't match the broadcast shape (3,)"):
            self.model(energies, signal)

    def test_call_invalid_signal_dimensions3D(self):
        energies = np.array([1., 2., 3.])
        signal = np.array([[[10., 15.]]])  # Invalid shape
        with pytest.raises(ValueError, match="Invalid number of signal array dimensions"):
            self.model(energies, signal)
