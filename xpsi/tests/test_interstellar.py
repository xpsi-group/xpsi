import unittest
import numpy as np
from xpsi.Interstellar import Interstellar

class TestInterstellar(Interstellar):
    def attenuation(self, energies):
        # Mock implementation for testing
        return np.ones_like(energies)

class TestInterstellarMethods(unittest.TestCase):
    def setUp(self):
        self.model = TestInterstellar()

    def test_attenuation(self):
        energies = np.array([1, 2, 3])
        expected = np.array([1, 1, 1])
        result = self.model.attenuation(energies)
        np.testing.assert_array_equal(result, expected)

    def test_call_with_1d_signal(self):
        energies = np.array([1, 2, 3])
        signal = np.array([10, 20, 30])
        self.model(energies, signal)
        expected_signal = np.array([10, 20, 30])
        np.testing.assert_array_equal(signal, expected_signal)

    def test_call_with_2d_signal(self):
        energies = np.array([1, 2, 3])
        signal = np.array([[10, 20], [30, 40], [50, 60]])
        self.model(energies, signal)
        expected_signal = np.array([[10, 20], [30, 40], [50, 60]])
        np.testing.assert_array_equal(signal, expected_signal)

    def test_invalid_signal_type(self):
        energies = np.array([1, 2, 3])
        with self.assertRaises(TypeError):
            self.model(energies, list(range(10)))  # Passing a list instead of ndarray

    def test_invalid_signal_dimensions(self):
        energies = np.array([1, 2, 3])
        with self.assertRaises(ValueError):
            self.model(energies, np.array([[1]]))  # Invalid dimensions

if __name__ == '__main__':
    unittest.main()


#
# # File: xpsi/tests/test_Interstellar.py
# import numpy as np
# import pytest
# from xpsi.Interstellar import Interstellar
#
#
# class MockInterstellar(Interstellar):
#     """ Mock class implementing the abstract attenuation method. """
#
#     def attenuation(self, energies):
#         return np.exp(-0.1 * energies)
#
#
# class TestInterstellar:
#     def test_signal_dimension_one(self):
#         energies = np.array([0.1, 0.5, 1.0, 2.0])
#         signal = np.array([10.0, 20.0, 30.0, 40.0])
#         interstellar = MockInterstellar()
#         interstellar(energies, signal)
#         expected_signal = signal * np.exp(-0.1 * energies)
#         print(signal)
#         print(expected_signal)
#         assert np.allclose(signal, expected_signal), "Attenuation failed for 1D signal."
#     #
#     # def test_signal_dimension_two(self):
#     #     energies = np.array([0.1, 0.5, 1.0, 2.0])
#     #     signal = np.array([[10.0, 15.0], [20.0, 25.0], [30.0, 35.0], [40.0, 45.0]])
#     #     interstellar = MockInterstellar()
#     #     interstellar(energies, signal)
#     #     expected_signal = np.multiply(signal, np.exp(-0.1 * energies)[:, np.newaxis])
#     #     assert np.allclose(signal, expected_signal), "Attenuation failed for 2D signal."
#
#     def test_invalid_signal_type(self):
#         energies = np.array([0.1, 0.5, 1.0, 2.0])
#         signal = [10.0, 20.0, 30.0, 40.0]  # Not a numpy array
#         interstellar = MockInterstellar()
#         with pytest.raises(TypeError, match='Signal must be a numpy.ndarray.'):
#             interstellar(energies, signal)
#
#     def test_invalid_signal_dimension(self):
#         energies = np.array([0.1, 0.5, 1.0, 2.0])
#         signal = np.array([[[10.0, 15.0]]])  # 3D array
#         interstellar = MockInterstellar()
#         with pytest.raises(ValueError, match='Invalid number of signal array dimensions.'):
#             interstellar(energies, signal)
#
#     def test_invalid_energy_dimension(self):
#         energies = np.array([[0.1, 0.5], [1.0, 2.0]])  # 2D array
#         signal = np.array([10.0, 20.0, 30.0, 40.0])
#         interstellar = MockInterstellar()
#         with pytest.raises(ValueError, match='Invalid number of energy array dimensions.'):
#             interstellar(energies, signal)