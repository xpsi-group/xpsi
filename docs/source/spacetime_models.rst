.. _spacetime_models:

Neutron star space-time models
==============================

X-PSI uses two different treatments of the space-time outside the neutron
star:

* the **Oblate Schwarzschild plus Doppler** (OS+D) approximation, which is used
  to compute pulse profiles, and therefore in all inference runs; and
* **numerical ray tracing in a quasi-Kerr metric**, which is used to image the
  star (for example for sky maps) with :class:`~.PhotospherePlotter.PhotospherePlotter`.

This page describes the physics behind each approach, the papers it is based
on, and the routines in the code that implement it. Equation numbers refer to
the papers listed in :ref:`spacetime_models_references`.

.. contents:: On this page
   :local:
   :depth: 2


Oblate Schwarzschild plus Doppler
---------------------------------

The approximation
~~~~~~~~~~~~~~~~~

The OS+D approximation was introduced by `Morsink et al. (2007) <M07_>`__, hereafter
`M07`_. It extends the older Schwarzschild plus Doppler (S+D) approximation
(`Miller & Lamb 1998 <ML98_>`__; `Poutanen & Gierliński 2003 <PG03_>`__), in which light bending is
computed as if the star were not rotating, and rotation enters only through
special-relativistic Doppler and aberration effects.

Rapid rotation also makes the star oblate. `Cadeau et al. (2007) <C07_>`__ showed, by
ray tracing in the exact numerical metric of rapidly rotating stars, that the
oblate *shape* of the surface has a significant effect on pulse profiles,
whereas frame dragging and the rotation-induced changes to the exterior metric
have only a small effect. The OS+D approximation therefore:

1. places the emitting surface on an **oblate** surface :math:`R(\theta)`;
2. propagates photons from that surface to the observer along **Schwarzschild**
   geodesics, which give the light-bending angles and arrival-time delays; and
3. adds the **Doppler** shift, Doppler boost and aberration caused by the
   motion of the surface, in the same way as in the S+D approximation.

The main effects of oblateness are geometric. The surface normal is tilted
away from the radial direction, which changes both the projected area of a
surface element and which parts of the surface are visible to the observer
(`M07`_, Section 4).

Dimensionless parameters
~~~~~~~~~~~~~~~~~~~~~~~~

The shape of the star and the other "universal" relations below depend on two
dimensionless combinations of the gravitational mass :math:`M`, the
equatorial radius :math:`R_{\rm eq}` and the angular spin frequency
:math:`\Omega = 2\pi f`:

.. math::

    \zeta = \frac{GM}{R_{\rm eq}c^2},
    \qquad
    \epsilon = \frac{\Omega^2 R_{\rm eq}^3}{GM}

(`M07`_, Eqs. 9 and 10; `AlGendy & Morsink 2014 <AM14_>`__, hereafter `AM14`_, Eqs. 1 and 2,
where they are called :math:`x` and :math:`\bar{\Omega}^2`). In X-PSI these
are available as :attr:`~.Spacetime.Spacetime.zeta` and
:attr:`~.Spacetime.Spacetime.epsilon`. :math:`R_{\rm eq}` is the
circumferential radius of the equator, which is the ``radius`` parameter of
:class:`~.Spacetime.Spacetime`.

Shape of the star
~~~~~~~~~~~~~~~~~

By default, X-PSI uses the universal relation for the oblate shape of `AM14`_,

.. math::

    \frac{R(\theta)}{R_{\rm eq}} = 1 + o_2(\zeta, \epsilon)\cos^2\theta,
    \qquad
    o_2 = \epsilon\,(-0.788 + 1.030\,\zeta),

where :math:`\theta` is the colatitude (`AM14`_, Eqs. 20 and 21, with the
coefficients from their Table 1). This replaces the Legendre-polynomial fit
originally given in `M07`_ (their Eq. 8 and Table 1). `AM14`_ fitted the relation
to stars with :math:`\epsilon \lesssim 0.1`, using a more recent set of
equations of state, and imposed :math:`R(\pi/2) = R_{\rm eq}` exactly.

The tilt of the surface is described by the function

.. math::

    f(\theta) = \frac{1}{R}\left(1 - \frac{2GM}{Rc^2}\right)^{-1/2}
                \frac{{\rm d}R}{{\rm d}\theta}

(`M07`_, Eq. 3). From it follow the angle :math:`\gamma` between the surface
normal and the radial direction, :math:`\cos\gamma = (1 + f^2)^{-1/2}` (`M07`_,
Eq. 13), and the proper area of a surface element,
:math:`{\rm d}S = R^2(\theta)\,(1 + f^2)^{1/2}\sin\theta\,{\rm d}\theta\,{\rm d}\phi`
(`M07`_, Eq. 2).

Light bending and time delays
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For a point on the surface at colatitude :math:`\theta` and rotational phase
:math:`\phi`, observed at inclination :math:`i`, the light-bending angle
:math:`\psi` satisfies

.. math::

    \cos\psi = \cos i\cos\theta + \sin i\sin\theta\cos\phi

(`M07`_, Eq. 15). The emission angle :math:`\alpha`, measured from the radial
direction, is related to :math:`\psi` through the impact parameter :math:`b`
(`M07`_, Eqs. 19 and 20), and the relative arrival time of a photon is given by
`M07`_, Eq. 38. Both integrals are evaluated in the Schwarzschild metric at the
local radius :math:`R(\theta)`. `M07`_, Eq. 38 gives the delay relative to a
radial photon from the same radius. Because :math:`R(\theta)` varies over an
oblate surface, X-PSI refers all delays to a radial photon emitted from the
equator, adding the radial light-travel time between :math:`R(\theta)` and
:math:`R_{\rm eq}` (`Salmi et al. 2018 <S18_>`__, Eqs. 9–11, with
:math:`r_{\rm ref} = R_{\rm eq}`).

Because the surface is oblate, some photons that reach the observer are
emitted *towards* the star (:math:`\alpha > \pi/2`). They pass through a
turning point before escaping. Their bending angles and delays are given by
`M07`_, Eqs. 40–45; X-PSI evaluates these integrals numerically rather than
using the series approximations of `M07`_, Eqs. 46 and 47. Higher-order images,
for photons that are bent by more than :math:`\pi`, can also be included.

Doppler effects and aberration
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The surface moves with speed

.. math::

    v = \frac{\Omega R(\theta)\sin\theta}{c}\left(1 - \frac{2GM}{R(\theta)c^2}\right)^{-1/2}

(`M07`_, Eq. 34), measured by a static observer at the surface. The angle
:math:`\xi` between the velocity and the initial photon direction follows from
`M07`_, Eq. 35, and the Doppler factor is

.. math::

    \delta = \frac{\sqrt{1 - v^2}}{1 - v\cos\xi}

(`M07`_, Eq. 33). The angle :math:`\sigma` between the initial photon direction
and the surface *normal*, as seen by a static observer, is

.. math::

    \cos\sigma = \cos\alpha\cos\gamma + \sin\alpha\sin\gamma\cos\beta

(`M07`_, Eq. 17), where :math:`\beta` is the angle of the spherical triangle
shown in `M07`_, Fig. 1. In the frame comoving with the surface, aberration
changes this angle to :math:`\cos\sigma_{\rm em} = \delta\cos\sigma` (`M07`_,
Eq. 36). A surface element is visible when :math:`\cos\sigma > 0`.

The observed photon energy is
:math:`E_{\rm obs} = \delta\,(1 - 2GM/Rc^2)^{1/2}\,E_{\rm em}` (`M07`_, Eq. 32).
The local specific intensity of the atmosphere model is evaluated at
:math:`E_{\rm em}` and :math:`\cos\sigma_{\rm em}`. The observed flux follows
from `M07`_, Eqs. 30, 37 and 39. The solid angle of a surface element (`M07`_,
Eq. 30) contains :math:`\cos\sigma` where a spherical star would have
:math:`\cos\alpha`. This is the "tilt" correction for oblateness.

X-PSI also includes a correction to `M07`_, Eq. 39: the emitting area must
be measured in the frame comoving with the surface, using simultaneity for
comoving observers. In the form of `Bogdanov et al. (2019) <B19_>`__, Eq. 31,
the flux from a surface element is

.. math::

    {\rm d}F(E) = \left(1 - \frac{2GM}{Rc^2}\right)^{1/2}
        \delta^3\, I'(E', \sigma_{\rm em})\, \cos\sigma_{\rm em}\,
        \left|\frac{\partial\cos\alpha}{\partial\cos\psi}\right|
        \frac{{\rm d}S'}{D^2},

where primes denote the comoving frame,
:math:`\cos\sigma_{\rm em} = \delta\cos\sigma` as above, and the comoving
area :math:`{\rm d}S' = \gamma\,{\rm d}S` is larger than the static-frame
area :math:`{\rm d}S` (`M07`_, Eq. 2) by the Lorentz factor
:math:`\gamma = (1 - v^2)^{-1/2}`. Written in static-frame quantities this is
:math:`\delta^4\gamma\cos\sigma\,{\rm d}S`, one factor of :math:`\gamma` more
than `M07`_, Eq. 39. Because :math:`\delta = 1/[\gamma(1 - v\cos\xi)]`, the
:math:`\gamma` from the area cancels the :math:`1/\gamma` inside one Doppler
factor, :math:`\gamma\cos\sigma_{\rm em} = \cos\sigma/(1 - v\cos\xi)`, so the
code evaluates the flux weight as :math:`\delta^3\cos\sigma/(1 - v\cos\xi)`.

The derivation of the comoving area is given by
`Nättilä & Pihajoki (2018) <NP18_>`__, Sect. 2.6, Eqs. 55–59, and by
`Lo et al. (2018) <Lo18_>`__, Sect. 1. `Lo et al. <Lo18_>`__ note that the
factor :math:`\gamma` was missing from all earlier calculations that followed
photons forward in time from the star to the observer, including `M07`_,
`Cadeau et al. (2007) <C07_>`__ and
`Poutanen & Gierliński (2003) <PG03_>`__. Leaving it out makes the flux from
each surface element too small by a factor :math:`1/\gamma`. For a small spot
this changes only the amplitude of the pulse profile. For an extended spot it
also changes the shape. The fractional error is at most :math:`\gamma - 1` at
the fastest-moving point of the spot: for example, about 2% for a spot on the
equator of a star spinning at 600 Hz (`Lo et al. 2018 <Lo18_>`__). The
corrected OS formalism is also given by `Salmi et al. (2018) <S18_>`__,
Sect. 2.1, Eqs. 5–7. X-PSI includes the correction, as discussed in `Bogdanov et al. (2019) <B19_>`__.

For emission that does not change with time and is the same at every
azimuth (for example a uniform surface modelled with
:class:`~.Elsewhere.Elsewhere`, or :class:`~.Everywhere.Everywhere` with
time-invariant emission), X-PSI sums the whole surface at a single instant and
uses a flux weight :math:`\delta^3\cos\sigma` instead of
:math:`\delta^3\cos\sigma/(1 - v\cos\xi)`. The extra factor
:math:`1/(1 - v\cos\xi) = \delta\gamma` for moving, localised emission
arises from binning photons by their arrival time
(`Cadeau et al. 2007 <C07_>`__, Sect. 2.7) together with the comoving area,
and does not apply to a steady, azimuthally symmetric surface.

..
   In integrator.pyx the flux weight is mu * eta**3 / superlum, with
   mu = cos(sigma), eta = delta and superlum = 1 - v cos(xi) (in M07's sign
   convention); integrator_for_time_invariance.pyx uses mu * eta**3.

Effective surface gravity
~~~~~~~~~~~~~~~~~~~~~~~~~

Atmosphere models depend on the local effective gravity, which varies over
the surface of a rotating star. By default, X-PSI uses the universal fit of `AM14`_ for
rapidly rotating stars,

.. math::

    \frac{g(\theta)}{g_0} = 1
        + (c_e\epsilon + d_e\epsilon^2 + f_e\epsilon^3)\sin^2\theta
        + (c_p\epsilon + d_p\epsilon^2 + f_p\epsilon^3 - d_{60}\epsilon^2)\cos^2\theta
        + d_{60}\epsilon^2\cos\theta,

where :math:`g_0 = GM / \bigl(R_{\rm eq}^2\sqrt{1 - 2\zeta}\bigr)` (`AM14`_,
Eqs. 39 and 50, with the coefficients from their Tables 4 and 5).

Choosing the shape
~~~~~~~~~~~~~~~~~~

The shape and effective gravity are chosen together, with the ``name``
argument of :class:`~.Star.Star`. The choice is stored in the spacetime and
can be read from :attr:`~.Spacetime.Spacetime.obl_surfgrav`:

* ``'universal'`` (the default) uses the `AM14`_ oblate shape and effective
  gravity described above.
* ``'sphere'`` uses a spherical surface, :math:`R(\theta) = R_{\rm eq}` and
  :math:`f = 0`, which reduces OS+D to the S+D approximation. The effective
  gravity is then :math:`g_0` everywhere on the surface. This option is meant
  for comparisons, not for modelling rapidly rotating stars.
* ``'psr j1231-1411'`` and ``'psr j1614-2230'`` use relations for the shape
  and effective gravity that are fitted for the spin frequency of one
  particular pulsar, in place of the universal relations. The J1614−2230
  relations are described by `Mauviard et al. (2026) <M26_>`__; the J1231−1411
  relations are from Kini et al. (in preparation).

The shape always has the form :math:`R(\theta)/R_{\rm eq} = 1 +
\epsilon\,o_2(\zeta)\cos^2\theta`; only the function :math:`o_2(\zeta)`
changes. The effective gravity for J1231−1411 keeps the functional form of
`AM14`_, Eq. 50, with refitted coefficients, whereas J1614−2230 uses its own
fitting formula (see the docstrings in
``xpsi/surface_radiation_field/effective_gravity_universal.pyx``). Shorthand
aliases are accepted (for example
``'j1614'``), and are listed in
:attr:`~.Spacetime.Spacetime.STAR_OBLATENESS_SURFGRAV_ALIASES`. Names for
further pulsars are reserved but raise :class:`NotImplementedError` until
their coefficients are added.

..
   Keep this list in step with Spacetime.STAR_OBLATENESS_SURFGRAV_CANONICAL
   and STAR_OBLATENESS_SURFGRAV_NOT_IMPLEMENTED.

Validity
~~~~~~~~

* Frame dragging and the quadrupole moment of the exterior metric are
  neglected. `Cadeau et al. (2007) <C07_>`__ found that these have only a small effect on
  pulse profiles. `M07`_ (Section 5) showed that OS light curves agree closely
  with exact ray tracing for a very oblate test star spinning at 600 Hz.
* The `AM14`_ shape relation was fitted in the slow-rotation regime,
  :math:`\epsilon \lesssim 0.1`, which covers most realistic neutron stars
  spinning below about 600–700 Hz (`AM14`_, Section 2.1). The effective-gravity
  relation was fitted up to the mass-shedding limit.
* The universal relations have a small residual dependence on the equation of
  state, at the per-cent level (`AM14`_).
* For slowly rotating stars (below about 300 Hz) oblateness matters less, and
  the S+D approximation is often sufficient (`M07`_, Section 1).

Where this is implemented
~~~~~~~~~~~~~~~~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 35 30 35

   * - Quantity
     - Reference
     - Code
   * - :math:`\zeta`, :math:`\epsilon`
     - `M07`_ Eqs. 9–10; `AM14`_ Eqs. 1–2
     - :attr:`~.Spacetime.Spacetime.zeta`,
       :attr:`~.Spacetime.Spacetime.epsilon`
   * - Oblate shape :math:`R(\theta)`
     - `AM14`_ Eqs. 20–21, Table 1
     - ``radiusNormalised`` in ``xpsi/cellmesh/mesh_tools.pyx``, with
       :math:`o_2` from ``oblateness_func_o2`` in
       ``xpsi/surface_radiation_field/effective_gravity_universal.pyx``
   * - Tilt function :math:`f(\theta)`
     - `M07`_ Eq. 3
     - ``f_theta`` in ``xpsi/cellmesh/mesh_tools.pyx``
   * - Surface area element
     - `M07`_ Eq. 2
     - integrands in ``xpsi/cellmesh/mesh_tools.pyx``
   * - Normal tilt :math:`\cos\gamma`
     - `M07`_ Eq. 13
     - ``xpsi/cellmesh/global_mesh.pyx`` and the other mesh constructors
   * - Bending angle and time delay, including initially ingoing photons
     - `M07`_ Eqs. 19–20, 38, 40–45; `Salmi et al. (2018) <S18_>`__ Eqs. 9–11 for the
       reference radius
     - ``xpsi/cellmesh/rays.pyx``; lag calibration in
       ``HotRegion.__calibrate_lag`` and ``Everywhere._calibrate_lag``
   * - :math:`\cos\psi`, :math:`v`, :math:`\cos\xi`, :math:`\delta`,
       :math:`\cos\sigma`, :math:`\cos\sigma_{\rm em}`, redshift and flux
       factors
     - `M07`_ Eqs. 15, 17, 30, 32–37, 39; :math:`\gamma` factor from
       `Nättilä & Pihajoki (2018) <NP18_>`__ Eqs. 55–59
     - ``xpsi/cellmesh/integrator.pyx`` and its variants
       (``integrator_for_*.pyx``, ``integratorIQU*.pyx``)
   * - Effective gravity
     - `AM14`_ Eqs. 39, 50, Tables 4–5
     - ``effectiveGravity`` in
       ``xpsi/surface_radiation_field/effective_gravity_universal.pyx``,
       which selects ``grav_algendy_morsink_2014``, ``grav_sphere`` or a
       pulsar-specific function
   * - Choice of relations
     - —
     - ``name`` argument of :class:`~.Star.Star`;
       :attr:`~.Spacetime.Spacetime.obl_surfgrav`

For a description of the OS+D implementation as used for NICER pulse-profile
modelling, see also `Bogdanov et al. (2019) <B19_>`__.


Quasi-Kerr ray tracing (imaging)
--------------------------------

To produce images of the star, such as sky maps and their animations, X-PSI
traces rays *backwards* from a grid of pixels on the image plane of a distant
observer until they hit the stellar surface. This is done by the
``xpsi.pixelmesh`` extension, called from
:meth:`~.PhotospherePlotter.PhotospherePlotter.image`. The method follows
`Psaltis & Johannsen (2012) <PJ12_>`__ and its application to rapidly
rotating neutron stars by `Bauböck et al. (2012) <B12_>`__. Unlike OS+D, the
exterior space-time includes frame dragging and the star's mass quadrupole
moment. The approach is general but much slower than OS+D, and it is not used
for pulse-profile computation during inference. Its uses are visualisation,
and independent checks of OS+D: an earlier version of this image-plane code
was cross-checked against the OS+D pulse profiles in the limit of zero spin
(`Bogdanov et al. 2019 <B19_>`__, Appendix B.5).

In this section, PJ12, B12, GB06, JP10 and C13 refer to
`Psaltis & Johannsen (2012) <PJ12_>`__, `Bauböck et al. (2012) <B12_>`__,
`Glampedakis & Babak (2006) <GB06_>`__,
`Johannsen & Psaltis (2010) <JP10_>`__ and
`Chan, Psaltis & Özel (2013) <C13_>`__.

The quasi-Kerr metric
~~~~~~~~~~~~~~~~~~~~~

The quasi-Kerr metric of `GB06`_ is the Kerr metric in Boyer–Lindquist
coordinates :math:`(t, r, \theta, \phi)` plus a correction that changes the
quadrupole moment,

.. math::

    g^{\mu\nu} = g_{\rm K}^{\mu\nu} + \kappa\, h^{\mu\nu}

(`GB06`_, Eqs. 2.13–2.16; `PJ12`_, Eqs. 1–5; `B12`_, Eqs. 2–6). Here
:math:`a = J/(Mc)` is the spin parameter, which sets the frame dragging, and
:math:`\kappa` is the dimensionless deviation of the quadrupole moment from
that of a Kerr black hole with the same spin:

.. math::

    Q = -M\left(a^2 + \kappa M^2\right)

(`PJ12`_, Eq. 6; `GB06`_, Eq. 2.3, where :math:`\kappa` is called
:math:`\epsilon`). The correction :math:`h^{\mu\nu}` has no
:math:`t\phi` component, so :math:`\kappa` does not change the frame
dragging. Its radial dependence enters through two functions
:math:`\mathcal{F}_1(r)` and :math:`\mathcal{F}_2(r)` (`GB06`_,
Eqs. A.9–A.10). The metric has unphysical pathologies close to the horizon
(`PJ12`_, Sect. 2), but these lie well inside any neutron star.

X-PSI computes the spin from the `AM14`_ universal relation for the moment of
inertia (`AM14`_, Eqs. 13–15), via :attr:`~.Spacetime.Spacetime.a`. It computes
:math:`\kappa` from the `AM14`_ quadrupole relations (`AM14`_,
Eqs. 17–19), via :attr:`~.Spacetime.Spacetime.q`. Note that
:attr:`~.Spacetime.Spacetime.q` is this deviation :math:`\kappa`, not the
quadrupole moment itself, and that it uses the coordinate-invariant
quadrupole with the sign of the :math:`\beta` term corrected relative to the
printed `AM14`_, Eq. 18. The docstring of :attr:`~.Spacetime.Spacetime.q`
gives the derivation and the conventions of the different papers. Because
imaging needs these spin and quadrupole relations,
:meth:`~.PhotospherePlotter.PhotospherePlotter.image` is currently only
available with ``Star(name='universal')``.

Image plane and initial conditions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The image plane lies at a distance :math:`d` from the star, perpendicular to
the line of sight at inclination :math:`i`, and centred on :math:`\phi = 0`
(`JP10`_, Fig. 1; `PJ12`_, Fig. 1). A pixel at image-plane coordinates
:math:`(X, Y)`, which are impact parameters, starts at the position given by
`PJ12`_, Eqs. 15–17, with a wavevector perpendicular to the image plane
(`PJ12`_, Eqs. 18–20; `JP10`_, Eqs. A1–A4). The time component follows from
the null condition (`PJ12`_, Eq. 13), and the photon's conserved impact
parameter :math:`b = L/E` is computed at the same time. X-PSI places the
image plane at :math:`d = 10^4\,b_{\rm max}`, where
:math:`b_{\rm max} = R_{\rm eq}/(1 - 2GM/R_{\rm eq}c^2)^{1/2}` is the
apparent radius of a non-rotating star. Because :math:`d/r \to 1` there, the
initial conditions are computed in extended precision, as `PJ12`_ (below
Eq. 20) advise.

The extent of the image is found first, by bisecting along two lines through
the image plane for the edges of the star, with a 0.5% margin. The pixels are
then laid out on an elliptical polar grid that covers it: :math:`N` rings of
:math:`N` pixels each, plus a central pixel, where :math:`N` is
``sqrt_num_rays``. Ring :math:`j` lies at a fractional radius
:math:`(j/N)^p`, with :math:`p` set by
``image_plane_radial_increment_power``. The default :math:`p = 1/2` gives
rings that enclose equal areas.

Integrating the rays
~~~~~~~~~~~~~~~~~~~~

Each ray is integrated backwards with the state vector
:math:`(t, \phi, r, {\rm d}r/{\rm d}\lambda, \theta,
{\rm d}\theta/{\rm d}\lambda)`. The equations for :math:`t` and :math:`\phi`
are first order, from the conserved energy and angular momentum (`PJ12`_,
Eqs. 7–10). The equations for :math:`r` and :math:`\theta` are the
second-order geodesic equations (`PJ12`_, Eqs. 11–12). X-PSI uses an
adaptive embedded Runge–Kutta–Cash–Karp (4, 5) stepper from the GNU
Scientific Library, with tolerances set by ``epsabs_ray`` and
``epsrel_ray``, where `PJ12`_ used a fixed-order scheme with their own step
rule. After each step, the null condition is monitored as in `PJ12`_
(Eqs. 13–14), and a step that drifts from it is redone with a smaller step
size.

Rays that pass close to the spin axis need special care, because the
:math:`\phi` equation diverges there. When a ray crosses the axis,
:math:`\theta` is reflected back into :math:`[0, \pi]` and :math:`\phi` is
shifted by :math:`\pi` (`C13`_, Eq. 3). If the adaptive stepper cannot take
a step near the axis, a single forward-Euler step carries the ray across it,
following `C13`_ (Sect. 2). A ray that misses the star, runs out of steps,
or cannot be carried across the axis is flagged and contributes no flux.

The stellar surface
~~~~~~~~~~~~~~~~~~~

The surface has the same oblate shape as in OS+D,
:math:`R(\theta)/R_{\rm eq} = 1 + \epsilon\,o_2(\zeta)\cos^2\theta`
(`AM14`_, Eqs. 20–21). In the ray tracer it is imposed in spherical
coordinates :math:`(r_{\rm sph}, \theta_{\rm sph})` built from the
Boyer–Lindquist coordinates,

.. math::

    r_{\rm sph} = \left(r^2 + a^2\sin^2\theta\right)^{1/2},
    \qquad
    r_{\rm sph}\cos\theta_{\rm sph} = r\cos\theta.

`AM14`_ fit their shape in a different radial coordinate, so this is an
approximation. The scale of the surface in :math:`r_{\rm sph}` is chosen so
that the circumferential equatorial radius equals :math:`R_{\rm eq}`
exactly; away from the equator the shape matches `AM14`_ only
approximately. The point where a ray meets the surface is found as in
`B12`_ (Sect. 3): when a step ends inside the star, the ray is interpolated
linearly between the last two steps and the crossing is found by bisection.

Redshift, emission angle and time delays
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The surface rotates rigidly with angular velocity :math:`\Omega`, with
4-velocity :math:`u^\mu \propto (1, 0, 0, \Omega)`. The redshift of a photon
with impact parameter :math:`b` is

.. math::

    1 + z = \frac{1 - b\,\Omega/c}
        {\left(-g_{tt} - 2g_{t\phi}\,\Omega/c - g_{\phi\phi}\,\Omega^2/c^2\right)^{1/2}}

(`Cadeau et al. 2007 <C07_>`__, Eq. 16, written for this metric; `PJ12`_,
Eq. 24). This one factor contains the gravitational redshift, the Doppler
shift and the effect of frame dragging. The emission angle
:math:`\sigma_{\rm em}` is measured from the surface normal in the frame
comoving with the surface (`Cadeau et al. 2007 <C07_>`__, Eqs. 18–20;
`B12`_, Eq. 13). Because the photon is projected directly into the comoving
frame, aberration is included automatically and is not applied as a separate
step as in OS+D.

The coordinate travel time of each ray is measured relative to that of the
central ray. For a non-rotating, spherical star, the exact travel time of a
radial ray from the equator is used instead. The difference, multiplied by
the angular frequency of the emission pattern
(:math:`2\pi` times the ``mode_frequency``), gives a phase lag: a pixel at
observed phase :math:`\phi_{\rm obs}` shows the surface as it was at
rotational phase :math:`\phi_{\rm obs}` minus that lag.

Flux and sky maps
~~~~~~~~~~~~~~~~~

For each pixel that hits the star, the observed specific intensity is

.. math::

    I_{\rm obs}(E) = \frac{I\bigl(E(1 + z), \cos\sigma_{\rm em}\bigr)}{(1 + z)^3},

from the invariance of :math:`I/E^3`, where :math:`I` is the comoving
intensity of the atmosphere model at the emission phase. The flux is the sum
of :math:`I_{\rm obs}` over pixels, each weighted by its solid angle, which
is its image-plane area divided by the square of the distance (`PJ12`_,
Eq. 25). Because the solid angle is measured directly on the observer's sky,
the comoving-area factor :math:`\gamma` and arrival-time binning of OS+D do
not arise here. When intensities are cached, the per-pixel values form the
sky maps that :class:`~.PhotospherePlotter.PhotospherePlotter` plots and
animates.

Accuracy and limitations
~~~~~~~~~~~~~~~~~~~~~~~~

* **Limb.** A pixel counts with its full area if its central ray hits the
  star, so the limb is resolved only to the ring spacing. With the default
  grid (``sqrt_num_rays = 100``, :math:`p = 1/2`), a uniformly bright disc
  comes out about 0.5% too bright. The error decreases with finer grids.
* **Small hot regions** can be missed entirely if the grid is too coarse to
  resolve them; the flux at that phase is then zero.
* **Surface shape.** The shape matches `AM14`_ exactly only at the equator
  (see above).
* **Tests.** ``xpsi/tests/test_pixelmesh.py`` checks that the imaged flux of a
  uniformly emitting, almost non-rotating star matches the analytic
  Schwarzschild result to within the limb error, that the flux equals the sum
  of the cached sky-map intensities, and that imaging a 600 Hz star with the
  `AM14`_ quadrupole completes.

Where the imaging is implemented
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

All files are in ``xpsi/pixelmesh/`` unless stated otherwise.

.. list-table::
   :header-rows: 1
   :widths: 35 30 35

   * - Step
     - Reference
     - Code
   * - Spin :math:`a` and quadrupole deviation :math:`\kappa`
     - `AM14`_ Eqs. 13–19
     - :attr:`~.Spacetime.Spacetime.a`, :attr:`~.Spacetime.Spacetime.q`
   * - Quasi-Kerr metric, :math:`\mathcal{F}_1`, :math:`\mathcal{F}_2`
     - `GB06`_ Eqs. 2.13–2.16, A.9–A.10; `PJ12`_ Eqs. 1–6
     - ``METRIC_qK.pyx``
   * - Image plane and initial conditions
     - `PJ12`_ Eqs. 13–20; `JP10`_ Eqs. A1–A4
     - ``BOUNDARY_CONDITIONS.pyx``
   * - Equations of motion, Christoffel symbols
     - `PJ12`_ Eqs. 7–12
     - ``RODES_qK.pyx``
   * - Ray integration, null check, polar axis
     - `PJ12`_ Eqs. 13–14; `C13`_ Sect. 2, Eq. 3
     - ``RK_IP2S_tracer.pyx``
   * - Extent of the image
     - —
     - ``get_IP_radius.pyx``
   * - Surface shape, intersection, redshift, emission angle
     - `AM14`_ Eqs. 20–21; `B12`_ Sect. 3, Eq. 13; `C07`_ Eqs. 16, 18–20;
       `PJ12`_ Eq. 24
     - ``surfaceBisection.pyx``, ``coordinateTransformation.pyx``
   * - Pixel grid, all rays, time delays
     - —
     - ``globalRayMap.pyx``
   * - Intensities, flux, sky-map cache
     - `PJ12`_ Eq. 25
     - ``integrator.pyx``;
       :meth:`~.PhotospherePlotter.PhotospherePlotter.image`


.. _spacetime_models_references:

References
----------

* `AlGendy, M., & Morsink, S. M. 2014, ApJ, 791, 78 <AM14_>`__
* `Bauböck, M., Psaltis, D., Özel, F., & Johannsen, T. 2012, ApJ, 753, 175 <B12_>`__
* `Bogdanov, S., Lamb, F. K., Mahmoodifar, S., et al. 2019, ApJL, 887, L26 <B19_>`__
* `Cadeau, C., Morsink, S. M., Leahy, D., & Campbell, S. S. 2007, ApJ, 654, 458 <C07_>`__
* `Chan, C.-K., Psaltis, D., & Özel, F. 2013, ApJ, 777, 13 <C13_>`__
* `Glampedakis, K., & Babak, S. 2006, CQGra, 23, 4167 <GB06_>`__
* `Johannsen, T., & Psaltis, D. 2010, ApJ, 718, 446 <JP10_>`__
* `Lo, K. H., Miller, M. C., Bhattacharyya, S., & Lamb, F. K. 2018, ApJ, 854, 187 <Lo18_>`__
  (erratum to `Lo et al. 2013, ApJ, 776, 19 <Lo13_>`__)
* `Mauviard, L., Guillot, S., Suleiman, L., Kini, Y., et al. 2026, arXiv:2609.00172 <M26_>`__
* `Miller, M. C., & Lamb, F. K. 1998, ApJL, 499, L37 <ML98_>`__
* `Morsink, S. M., Leahy, D. A., Cadeau, C., & Braga, J. 2007, ApJ, 663, 1244 <M07_>`__
* `Nättilä, J., & Pihajoki, P. 2018, A&A, 615, A50 <NP18_>`__
* `Poutanen, J., & Gierliński, M. 2003, MNRAS, 343, 1301 <PG03_>`__
* `Psaltis, D., & Johannsen, T. 2012, ApJ, 745, 1 <PJ12_>`__
* `Riley, T. E., Watts, A. L., Bogdanov, S., et al. 2019, ApJL, 887, L21 <R19_>`__
* `Salmi, T., Nättilä, J., & Poutanen, J. 2018, A&A, 618, A161 <S18_>`__

..
   SciX link targets, one per paper. In the text, cite a paper with
   `Author et al. (year) <KEY_>`__ (anonymous link) or `KEY`_.

.. _AM14: https://scixplorer.org/abs/2014ApJ...791...78A/abstract
.. _B12: https://scixplorer.org/abs/2012ApJ...753..175B/abstract
.. _B19: https://scixplorer.org/abs/2019ApJ...887L..26B/abstract
.. _C07: https://scixplorer.org/abs/2007ApJ...654..458C/abstract
.. _C13: https://scixplorer.org/abs/2013ApJ...777...13C/abstract
.. _GB06: https://scixplorer.org/abs/2006CQGra..23.4167G/abstract
.. _JP10: https://scixplorer.org/abs/2010ApJ...718..446J/abstract
.. _Lo13: https://scixplorer.org/abs/2013ApJ...776...19L/abstract
.. _Lo18: https://scixplorer.org/abs/2018ApJ...854..187L/abstract
.. _M07: https://scixplorer.org/abs/2007ApJ...663.1244M/abstract
.. _M26: https://scixplorer.org/abs/arXiv:2609.00172/abstract
.. _ML98: https://scixplorer.org/abs/1998ApJ...499L..37M/abstract
.. _NP18: https://scixplorer.org/abs/2018A%26A...615A..50N/abstract
.. _PG03: https://scixplorer.org/abs/2003MNRAS.343.1301P/abstract
.. _PJ12: https://scixplorer.org/abs/2012ApJ...745....1P/abstract
.. _R19: https://scixplorer.org/abs/2019ApJ...887L..21R/abstract
.. _S18: https://scixplorer.org/abs/2018A%26A...618A.161S/abstract
