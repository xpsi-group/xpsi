"""Tests of the image-plane ray tracer (``xpsi/pixelmesh``), run through
:meth:`xpsi.PhotospherePlotter.image`.

1. A uniformly emitting, (almost) non-rotating star: the imaged flux must match
   the analytic Schwarzschild result, and must equal the sum of the cached
   sky-map intensities over the image-plane grid.
2. A 600 Hz star with the AlGendy & Morsink (2014) quadrupole: imaging must
   finish. It used to hang because of catastrophic cancellation in the
   quasi-Kerr metric functions F1, F2 at large r. It runs in a subprocess so
   that a regression fails the test rather than stalling the test suite.

These tests assume X-PSI is compiled with the default ``local_variables.pyx``;
they fail if it has been replaced by a custom extension.
"""

import os
import subprocess
import sys
import textwrap

import numpy as np
import pytest

import xpsi
from xpsi.global_imports import _kpc, _km, _GM, _keV, _k_B

pytestmark = pytest.mark.filterwarnings("ignore")

# Physical constants as used by X-PSI's blackbody (surface_radiation_field).
_erg = 1.0e-7
_PLANCK_PHOTON = 5.040366110812353e22  # 2 / (h^3 c^2) in X-PSI's units

MASS, RADIUS, DISTANCE, LOGT = 1.4, 12.0, 0.2, 6.5  # Msun, km, kpc, log10(K)
ENERGIES = np.array([0.5, 1.0, 2.0, 4.0])  # keV


def _star(frequency, cos_inclination):
    spacetime = xpsi.Spacetime(bounds={}, values=dict(frequency=frequency,
                                                      mass=MASS,
                                                      radius=RADIUS,
                                                      distance=DISTANCE,
                                                      cos_inclination=cos_inclination))
    everywhere = xpsi.Everywhere(time_invariant=True, bounds={},
                                 values=dict(temperature=LOGT),
                                 sqrt_num_cells=32, num_rays=128)
    photosphere = xpsi.Photosphere(hot=None, elsewhere=None,
                                   everywhere=everywhere,
                                   values=dict(mode_frequency=frequency))
    xpsi.Star(spacetime, photosphere).update()
    return spacetime, photosphere


class _UniformStarPlotter(xpsi.PhotospherePlotter):
    """ The whole surface emits a blackbody at LOGT. In the compiled default
    ``local_variables`` the primary superseding region is a circle centred
    on the pole; with angular radius pi it covers the star. """

    @property
    def global_variables(self):
        gv = np.zeros(16)
        gv[2] = np.pi    # primary superseding radius
        gv[12] = LOGT    # primary superseding temperature
        return gv


def _area_weights(r_mid, semi_major, semi_minor):
    """ Image-plane area of each ray: a central disc and N rings of N rays,
    with ring boundaries halfway between ring radii (as in integrator.pyx). """
    N = len(r_mid)
    r = np.append(r_mid, 1.0)
    w = np.empty(N * N + 1)
    w[0] = np.pi * (0.5 * r[1]) ** 2
    for i in range(N):
        r1 = 0.5 * (r[i + 1] + r[i])
        r2 = 0.5 * (r[i + 2] + r[i + 1]) if i < N - 1 else 2.0 * r[N] - r1
        w[1 + i * N:1 + (i + 1) * N] = 0.5 * (r2 ** 2 - r1 ** 2) * 2.0 * np.pi / N
    return w * semi_major * semi_minor


@pytest.fixture(scope="module")
def uniform_image():
    spacetime, photosphere = _star(frequency=1.0, cos_inclination=0.5)
    plotter = _UniformStarPlotter(photosphere)
    plotter.image(reimage=True,
                  cache_intensities=0.1,
                  energies=ENERGIES,
                  phases=np.array([0.0]),
                  sqrt_num_rays=100,
                  threads=2,
                  max_steps=100000,
                  epsrel_ray=1.0e-12)
    return spacetime, plotter.images


def test_uniform_star_matches_analytic_flux(uniform_image):
    """ Flux from a uniform blackbody star in Schwarzschild: the star's
    apparent radius is b = R / sqrt(1 - u) with u = r_s / R, and every ray
    carries I(E / sqrt(1 - u)) (1 - u)^(3/2), so the photon flux is
    N(E) = pi b^2 I(E / sqrt(1 - u)) (1 - u)^(3/2) / (D^2 E). At 1 Hz the
    rotational corrections are ~1e-6. The tolerance covers the
    discretisation of the limb by the image-plane grid (+0.5% at 100 rays
    per ring). """
    spacetime, images = uniform_image
    flux = images[0][:, 0]

    R = RADIUS * _km
    u = 2.0 * MASS * _GM / R
    b = R / np.sqrt(1.0 - u)
    D = DISTANCE * _kpc
    kT = _k_B * 10.0 ** LOGT / _keV

    def blackbody(E):  # specific intensity in X-PSI's units, E in keV
        return _erg * _PLANCK_PHOTON * E ** 3 / np.expm1(E / kT)

    # I/E^3 is invariant; the photon flux divides by the observed energy.
    expected = (np.pi * b ** 2 / D ** 2
                * blackbody(ENERGIES / np.sqrt(1.0 - u)) * (1.0 - u) ** 1.5
                / (ENERGIES * _keV))

    assert np.allclose(flux / expected, 1.0, rtol=0.0, atol=0.01)
    # The limb error is the same at all energies for a uniform star; at 1 Hz
    # the energy dependence from rotation is ~1e-6.
    assert np.ptp(flux / expected) < 1.0e-4


def test_flux_equals_sum_of_cached_sky_map(uniform_image):
    """ The integrator's flux must equal the cached intensities summed with
    the image-plane area weights. Catches out-of-bounds area weights and
    misalignment between the cache and the ray map. """
    spacetime, images = uniform_image
    flux = images[0]                              # energy x phase
    weights = _area_weights(images[9], images[10], images[11])
    intensity = np.asarray(images[12], dtype=float)  # phase x energy x ray
    D = DISTANCE * _kpc
    summed = (intensity * weights).sum(axis=-1).T * _erg * _PLANCK_PHOTON / (D ** 2 * _keV)
    assert np.allclose(summed, flux, rtol=1.0e-8, atol=0.0)

    # Every ray that hits the star has a cached intensity, and no other ray.
    hit = images[5] > 0.0
    cached = intensity[0, 0] > 0.0
    assert np.array_equal(hit, cached)


# Code for the subprocess in test_imaging_with_quadrupole_finishes, which a
# fresh Python process runs so that a hanging ray can be killed by a timeout.
# It repeats the contents of _star and _UniformStarPlotter (with the constants
# written out), because it cannot import them from this file.
_QUADRUPOLE_SCRIPT = textwrap.dedent("""
    import numpy as np
    import xpsi
    spacetime = xpsi.Spacetime(bounds={}, values=dict(frequency=600.0, mass=1.4,
                               radius=12.0, distance=0.2, cos_inclination=0.5))
    everywhere = xpsi.Everywhere(time_invariant=True, bounds={},
                                 values=dict(temperature=6.5),
                                 sqrt_num_cells=32, num_rays=128)
    photosphere = xpsi.Photosphere(hot=None, elsewhere=None, everywhere=everywhere,
                                   values=dict(mode_frequency=600.0))
    xpsi.Star(spacetime, photosphere).update()
    assert spacetime.q > 0.3  # AlGendy & Morsink quadrupole, about 0.40 here

    class Plotter(xpsi.PhotospherePlotter):
        @property
        def global_variables(self):
            gv = np.zeros(16)
            gv[2] = np.pi
            gv[12] = 6.5
            return gv

    plotter = Plotter(photosphere)
    plotter.image(reimage=True, energies=np.array([1.0]),
                  phases=np.array([0.0]), sqrt_num_rays=200, threads=2,
                  max_steps=100000, epsrel_ray=1.0e-12)
    flux = plotter.images[0][0, 0]
    assert np.isfinite(flux) and flux > 0.0
    print("FLUX", flux)
""")


def test_imaging_with_quadrupole_finishes(tmp_path):
    """ A 600 Hz star (a = 651 m, q = 0.40) at i = 60 deg with 200 rays per
    ring: before the large-r series in METRIC_qK.pyx, one ray of this image
    never finished. Takes about 10-20 s.

    The subprocess must import the same (compiled) xpsi as this process, not
    the uncompiled source tree that would be found first if it ran in the
    repository root. So it runs in a temporary directory, with the directory
    containing this process's xpsi package put first on its path. """
    timeout = 200  # seconds
    env = dict(os.environ)
    xpsi_parent = os.path.dirname(os.path.dirname(os.path.abspath(xpsi.__file__)))
    env["PYTHONPATH"] = os.pathsep.join(
        [xpsi_parent] + [p for p in env.get("PYTHONPATH", "").split(os.pathsep) if p])
    try:
        proc = subprocess.run([sys.executable, "-c", _QUADRUPOLE_SCRIPT],
                              capture_output=True, text=True, timeout=timeout,
                              env=env, cwd=tmp_path)
    except subprocess.TimeoutExpired:
        pytest.fail("Imaging a 600 Hz star with a quadrupole did not finish "
                    f"within {timeout} s.")
    assert proc.returncode == 0, proc.stderr[-3000:]
    assert "FLUX" in proc.stdout
    assert "Steps exceeded" not in proc.stdout, \
        "Some rays ran out of integration steps (max_steps)."
