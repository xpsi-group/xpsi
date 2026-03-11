import numpy as np
from xpsi.cellmesh.common_functions import compute_pol_ang_py, disk_block_py

class TestCellmeshCommonFunctions(object):

    def test_compute_pol_ang_basic(self):
        'Testing that the polarization angle calculation gives a finite result.'
        result = compute_pol_ang_py(
            0.1,   # leaves_kdx
            0.2,   # sin_psi
            0.98,  # cos_psi
            0.3,   # sin_alpha
            0.95,  # cos_alpha
            0.4,   # sin_theta_i
            0.92,  # cos_theta_i
            0.5,   # sin_i
            0.86,  # cos_i
            0.6,   # sin_gamma
            0.8,   # cos_gamma
            1.0,   # Grav_z
            0.7,   # mu
            1.0,   # eta
            0.2,   # beta
            1.1,   # Lorentz
            0.9    # cos_xi
        )

        assert np.isfinite(result)

    def reference_chi(self, theta, phi, i):
        """
        Non-relativistic polarization angle from analytic formula.
        """
        sin_chi0 = -np.sin(theta) * np.sin(phi)
        cos_chi0 = np.sin(i)*np.cos(theta) - np.cos(i)*np.sin(theta)*np.cos(phi)
        return np.arctan2(sin_chi0, cos_chi0)

    def test_compute_pol_ang_non_relativistic_limit(self):
        'Testing that the polarization angle calculation gives the correct result in the non-relativistic limit.'

        rng = np.random.default_rng(1234)

        # sample geometry
        for _ in range(10):
            theta = rng.uniform(0.001, np.pi-0.001)
            phi = rng.uniform(0, 2*np.pi)
            i = rng.uniform(0.001, np.pi-0.001)

            print("theta, phi, i: ",theta,phi,i)

            sin_theta = np.sin(theta)
            cos_theta = np.cos(theta)

            sin_i = np.sin(i)
            cos_i = np.cos(i)

            # values that remove relativistic terms
            sin_gamma = 0.0
            cos_gamma = 1.0
            beta = 0.0
            Lorentz = 1.0

            # other parameters (not relevant in this limit)
            sin_psi = 0.3
            cos_psi = np.sqrt(1 - sin_psi**2)
            sin_alpha = 0.4
            cos_alpha = np.sqrt(1 - sin_alpha**2)
            Grav_z = 1.0
            mu = 0.5
            eta = 1.0
            cos_xi = 0.5

            chi_model = compute_pol_ang_py(
                phi,
                sin_psi,
                cos_psi,
                sin_alpha,
                cos_alpha,
                sin_theta,
                cos_theta,
                sin_i,
                cos_i,
                sin_gamma,
                cos_gamma,
                Grav_z,
                mu,
                eta,
                beta,
                Lorentz,
                cos_xi,
            )

            chi_expected = self.reference_chi(theta, phi, i)

            assert np.isclose(chi_model, chi_expected, atol=1e-12)

    def test_disk_blocking_when_Rin_smaller_than_radius(self):
        'Testing that the accretion disk blocks the rays always when the inner disc radius is smaller than the radius of the star.'

        rng = np.random.default_rng(1234)

        for _ in range(10):

            radius = rng.uniform(5.0, 16.0)
            R_in = rng.uniform(0.1, radius * 0.999)  # always smaller than radius
            i = rng.uniform(0.001, np.pi/2 - 0.001)
            psi = rng.uniform(0.001, np.pi - 0.001)
            theta_i = rng.uniform(np.pi/2 + 0.01, np.pi - 0.01)

            cos_i = np.cos(i)
            cos_psi = np.cos(psi)
            cos_theta_i = np.cos(theta_i)

            sin_alpha = rng.uniform(0.1, 0.9)
            r_s_over_r_i = rng.uniform(0.01, 0.2)

            theta_i_over_pi = theta_i / np.pi

            result = disk_block_py(
                R_in,
                cos_i,
                cos_psi,
                cos_theta_i,
                r_s_over_r_i,
                radius,
                sin_alpha,
                theta_i_over_pi
            )

            assert result == 0
