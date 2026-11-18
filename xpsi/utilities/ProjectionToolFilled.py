"""
2D projection of X-PSI hot regions, coloured by temperature.

``plot_projection_filled`` renders the stellar surface pixel by pixel and
colours every point that emits by its temperature:

* ``ST``                 -> superseding region at ``super_temperature``.
* ``CST`` / ``EST`` / ``PST``
                         -> superseding region minus the (zero-temperature)
                            omission region. Only the surviving area is
                            coloured, so rings, crescents and arcs appear as
                            they would in X-PSI.
* ``CDT`` / ``EDT`` / ``PDT``
                         -> superseding region at ``super_temperature``;
                            ceding region at ``cede_temperature`` *only where
                            it is not superseded*.

Models with one or two hot regions are supported, using the model strings,
parameter names and points of view described in the "Emitting patterns 2D
projection" tutorial.

No physcis are included (e.g. relativistic and Doppler effects): this is a plain 
orthographic projection of the unit sphere. Emission on the hemisphere facing 
the viewer is drawn solid; emission on the far hemisphere is hatched (or dimmed), 
as if seen through a glass sphere. ``back_view=True`` adds a panel showing the 
opposite side.

``animate_projection_filled`` writes a movie (MP4 or GIF) of the star rotating
as seen by the observer, using the same rendering.

The rotation conventions are those of ``ProjectionTool.py``, which provides
the older outline-only version of this plot (``plot_projection_general``).

This file only needs numpy, scipy and matplotlib, so it can be used without
X-PSI installed.
"""

import numpy as np
from matplotlib.patches import Patch
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.cm import ScalarMappable
from scipy.ndimage import gaussian_filter
import shutil
import subprocess

__all__ = ["plot_projection_filled", "animate_projection_filled",
           "parse_hot_regions", "HOT_GLOW_CMAP",
           "INFERNO_DARK_CMAP"]


# ---------------------------------------------------------------------------
# Colours
# ---------------------------------------------------------------------------

#: Alternative colour map: orange (cooler) -> pink -> magenta -> pale lilac.
HOT_GLOW_CMAP = LinearSegmentedColormap.from_list(
    "xpsi_hot_glow",
    ["#ff7b00", "#ff9a3c", "#ff4f8b", "#ff1fd0", "#ff00ff", "#f2c8ff"])

#: Default colour map for the dark style: inferno without its darkest 20 %,
#: so that the coolest regions do not vanish into the night-sky background.
INFERNO_DARK_CMAP = LinearSegmentedColormap.from_list(
    "inferno_dark", plt.get_cmap("inferno")(np.linspace(0.2, 1.0, 256)))

STYLES = {
    "dark": dict(bg_top="#0e1d3a", bg_bottom="#05070f", grid="#2f6bff",
                 limb="#3f86ff", disk="#0a1630", disk_alpha=0.55,
                 text="#c9d6ff", equator="#5f9bff", pole="#9cc0ff"),
    "light": dict(bg_top="#ffffff", bg_bottom="#ffffff", grid="#7f93c9",
                  limb="#304a8c", disk="#ffffff", disk_alpha=0.6,
                  text="#1a1a1a", equator="#304a8c", pole="#304a8c"),
}

_TAGS = ("p", "s")
_TAGNAME = {"p": "Primary", "s": "Secondary"}


# ---------------------------------------------------------------------------
# Geometry helpers (same conventions as ProjectionTool.transform)
# ---------------------------------------------------------------------------

def _unit(theta, phi_cycles):
    """Unit vector at colatitude theta [rad] and phase phi [cycles]."""
    a = 2.0 * np.pi * phi_cycles
    return np.array([np.sin(theta) * np.cos(a),
                     np.sin(theta) * np.sin(a),
                     np.cos(theta)])


def _rotation(thetaR, phiR, phi0=0.0):
    """3x3 matrix identical to ProjectionTool.transform(thetaR, phiR, V, phi0)."""
    phi0 = 2.0 * np.pi * phi0
    ax = np.array([-np.sin(phi0), np.cos(phi0), 0.0])
    c, s = np.cos(thetaR), np.sin(thetaR)
    K = np.array([[0.0, -ax[2], ax[1]],
                  [ax[2], 0.0, -ax[0]],
                  [-ax[1], ax[0], 0.0]])
    Mtheta = c * np.eye(3) + s * K + (1.0 - c) * np.outer(ax, ax)
    p = 2.0 * np.pi * phiR
    Mphi = np.array([[np.cos(p), -np.sin(p), 0.0],
                     [np.sin(p), np.cos(p), 0.0],
                     [0.0, 0.0, 1.0]])
    return Mphi @ Mtheta


def _view_matrix(theta_sub, phi_sub):
    """Stellar frame -> viewer frame (viewer looks down the +z axis)."""
    return _rotation(-theta_sub, -phi_sub, phi_sub)


def _to_screen(v_obs):
    """Viewer frame -> plot coordinates, as in ProjectionTool (x=y_obs, y=-x_obs)."""
    return v_obs[1], -v_obs[0]


# ---------------------------------------------------------------------------
# Parameter parsing
# ---------------------------------------------------------------------------

def _names(dictVp):
    try:
        return list(dictVp.keys())
    except AttributeError:
        return list(dictVp.names)


def _get(dictVp, names, key):
    if key not in names:
        raise KeyError("Parameter '%s' is required for this model but was "
                       "not found." % key)
    return float(dictVp[key])


def _flavour(hot_name):
    """Return (kind, geometry) e.g. ('ST', 'P') for 'PST', ('DT', 'C') for 'CDT'."""
    if hot_name == "ST":
        return "ST", None
    if len(hot_name) == 3 and hot_name[0] in "CEP" and hot_name[1:] in ("ST", "DT"):
        return hot_name[1:], hot_name[0]
    raise ValueError("Hot spot name '%s' not recognised. Use ST, CST, EST, "
                     "PST, CDT, EDT or PDT." % hot_name)


def _read_region(dictVp, names, tag, hot_name, antiphase=False):
    """ Build one hot region from tagged parameters.

    Mirrors ``ProjectionTool.fillVECTORS`` (and the ``xpsi.HotRegion``
    conventions).

    :param dictVp: Parameter values (or an X-PSI ParameterSubspace).
    :param list names: Available parameter names.
    :param str tag: Hot-region tag, ``'p'`` (primary) or ``'s'`` (secondary).
    :param str hot_name: Hot-spot name, e.g. ``'ST'``, ``'PST'`` or ``'CDT'``.
    :param bool antiphase: If True, offset the secondary by 0.5 cycles.

    :returns:
        A dictionary describing the region, with keys ``name``, ``tag``,
        ``super = (theta, phi, zeta, logT)``,
        ``omit = (theta, phi, zeta)`` or ``None``, and
        ``cede = (theta, phi, zeta, logT)`` or ``None``.
    """
    kind, geom = _flavour(hot_name)
    g = lambda k: _get(dictVp, names, tag + "__" + k)

    phase = g("phase_shift")
    if antiphase and tag == "s":
        phase += 0.5
    s_theta, s_zeta, s_T = g("super_colatitude"), g("super_radius"), g("super_temperature")
    s_phi = phase
    omit = cede = None

    if kind == "ST" and geom is not None:
        o_zeta = g("omit_radius")
        if geom == "C":
            omit = (s_theta, phase, o_zeta)
        else:
            # omission sits at the phase shift; superseding region is offset
            omit = (g("omit_colatitude"), phase, o_zeta)
            s_phi = phase - g("omit_azimuth") / (2.0 * np.pi)
    elif kind == "DT":
        c_zeta, c_T = g("cede_radius"), g("cede_temperature")
        if geom == "C":
            cede = (s_theta, phase, c_zeta, c_T)
        else:
            cede = (g("cede_colatitude"), phase + g("cede_azimuth") / (2.0 * np.pi),
                    c_zeta, c_T)

    return dict(name=hot_name, tag=tag, super=(s_theta, s_phi, s_zeta, s_T),
                omit=omit, cede=cede)


def _antipode(comp):
    if comp is None:
        return None
    theta, phi = comp[0], comp[1]
    phi = phi + 0.5 if phi < 0.5 else phi - 0.5
    return (np.pi - theta, phi) + tuple(comp[2:])


def _antipodal_region(region, tag):
    return dict(name=region["name"], tag=tag,
                super=_antipode(region["super"]),
                omit=_antipode(region["omit"]),
                cede=_antipode(region["cede"]))


def parse_hot_regions(dictVp, model, antiphase=False):
    """
    Translate an X-PSI parameter dictionary + model string into a list of hot
    regions (see ``_read_region`` for the format).

    Model strings follow ProjectionTool: e.g. "ST", "PDT", "ST+PST",
    "ST-S", "PST-Ua", "ST-U". Only one or two hot regions are supported.
    """
    names = _names(dictVp)
    regions = []

    if "+" in model:
        parts = model.split("+")
        if len(parts) > 2:
            raise ValueError("Only models with one or two hot regions are supported.")
        for tag, hot_name in zip(_TAGS, parts):
            regions.append(_read_region(dictVp, names, tag, hot_name, antiphase))

    elif model.endswith("-S"):
        prim = _read_region(dictVp, names, "p", model[:-2])
        regions = [prim, _antipodal_region(prim, "s")]

    elif model.endswith("-Ua"):
        base = model[:-3]
        prim = _read_region(dictVp, names, "p", base)
        sec = _antipodal_region(prim, "s")
        kind, _ = _flavour(base)
        g = lambda k: _get(dictVp, names, "s__" + k)
        sec["super"] = sec["super"][:2] + (g("super_radius"), g("super_temperature"))
        if sec["omit"] is not None:
            sec["omit"] = sec["omit"][:2] + (g("omit_radius"),)
        if sec["cede"] is not None:
            sec["cede"] = sec["cede"][:2] + (g("cede_radius"), g("cede_temperature"))
        regions = [prim, sec]

    elif model.endswith("-U"):
        base = model[:-2]
        regions = [_read_region(dictVp, names, t, base, antiphase) for t in _TAGS]

    elif "-" in model:
        raise ValueError("Model suffix not recognised (use -S, -U or -Ua).")

    else:
        regions = [_read_region(dictVp, names, "p", model)]

    _sanity_check(regions)
    return regions


def _sanity_check(regions):
    for r in regions:
        label = _TAGNAME[r["tag"]]
        for key in ("super", "omit", "cede"):
            comp = r[key]
            if comp is None:
                continue
            theta, phi, zeta = comp[:3]
            if not (0.0 <= theta <= np.pi):
                raise ValueError("%s %s colatitude out of range [0, pi]." % (label, key))
            if not (0.0 <= zeta <= np.pi / 2):
                raise ValueError("%s %s angular radius out of range [0, pi/2]." % (label, key))
            if not (-2.0 <= phi <= 2.0):
                raise ValueError("%s %s phase out of range [-2, 2]." % (label, key))
            if len(comp) == 4 and not (3.0 <= comp[3] <= 15.0):
                raise ValueError("%s %s log10 temperature out of range." % (label, key))
        if r["omit"] is not None and r["name"].startswith("C") and r["omit"][2] >= r["super"][2]:
            print("WARNING: %s hot region: emitting component is completely "
                  "masked by the omitting one (nothing will be coloured)." % label)
        if r["cede"] is not None and r["name"].startswith("C") and r["cede"][2] <= r["super"][2]:
            print("WARNING: %s hot region: ceding component is completely "
                  "masked by the superseding one." % label)


# ---------------------------------------------------------------------------
# Temperature map on the sphere
# ---------------------------------------------------------------------------

def _inside(n, comp):
    """Boolean mask: unit vectors n (3, ...) inside the circular cap comp."""
    c = _unit(comp[0], comp[1])
    return np.tensordot(c, n, axes=1) >= np.cos(comp[2])


def surface_log_temperature(n, regions, warn=True):
    """
    log10(T/K) of the emitting surface at unit vectors n (shape (3, ...)).
    Non-emitting points (cold surface or masked by an omission region) are NaN.
    """
    logT = np.full(n.shape[1:], np.nan)
    overlap = False
    for r in regions:
        Tr = np.full(n.shape[1:], np.nan)
        sup = _inside(n, r["super"])
        if r["cede"] is not None:
            ced = _inside(n, r["cede"]) & ~sup
            Tr[ced] = r["cede"][3]
        if r["omit"] is not None:
            sup &= ~_inside(n, r["omit"])
        Tr[sup] = r["super"][3]

        both = np.isfinite(logT) & np.isfinite(Tr)
        overlap |= bool(both.any())
        logT = np.where(both, np.maximum(logT, Tr), np.where(np.isfinite(Tr), Tr, logT))
    if overlap and warn:
        print("WARNING: different hot regions overlap; the hotter one is shown "
              "where they do (X-PSI itself would simply add their signals).")
    return logT


def _check_overlap(regions, n=200000):
    """Warn once if different hot regions overlap anywhere on the sphere."""
    if len(regions) < 2:
        return
    i = np.arange(n) + 0.5
    z = 1.0 - 2.0 * i / n
    a = np.pi * (1.0 + 5.0 ** 0.5) * i
    rho = np.sqrt(1.0 - z * z)
    surface_log_temperature(np.array([rho * np.cos(a), rho * np.sin(a), z]), regions, warn=True)


# ---------------------------------------------------------------------------
# Point of view
# ---------------------------------------------------------------------------

def _pov_angles(POV, dictVp, regions):
    """Return (theta_sub, phi_sub, label) for the requested point of view."""
    names = _names(dictVp)
    if isinstance(POV, str):
        if POV == "":
            return 0.0, 0.0, "NorthPole"
        if POV == "I":
            if "cos_inclination" in names:
                cosi = float(dictVp["cos_inclination"])
            elif "sin_inclination" in names:
                cosi = np.sqrt(1.0 - float(dictVp["sin_inclination"]) ** 2)
            else:
                raise ValueError("POV 'I' needs cos_inclination or sin_inclination.")
            return np.arccos(cosi), 0.0, "Earth"

        options = {}
        for r in regions:
            T = r["tag"].upper()
            nm = _TAGNAME[r["tag"]]
            if r["omit"] is not None:
                options[T + "E"] = (r["super"], nm + "Emission")
                options[T + "O"] = (r["omit"], nm + "Omission")
            elif r["cede"] is not None:
                options[T + "S"] = (r["super"], nm + "Super")
                options[T + "C"] = (r["cede"], nm + "Cede")
            else:
                options[T] = (r["super"], nm)
        if POV not in options:
            raise ValueError("Point of view '%s' not allowed for this model. "
                             "Possibilities are: %s" % (POV, ["I"] + list(options)))
        comp, label = options[POV]
        return comp[0], comp[1], label

    v = np.asarray(POV, dtype=float)
    if abs(np.linalg.norm(v) - 1.0) > 1e-6:
        print("WARNING: POV vector normalised to unit length (plot_projection_general "
              "does not normalise, so its view will differ for this input).")
    v = v / np.linalg.norm(v)
    theta_sub = np.arccos(np.clip(v[2], -1, 1))
    phi_sub = 0.0
    if abs(np.sin(theta_sub)) > 1e-16:
        a0 = np.clip(v[0] / np.sin(theta_sub), -1, 1)
        phi_sub = (np.arccos(a0) if v[1] >= 0 else np.arccos(a0) + np.pi) / (2 * np.pi)
    return theta_sub, phi_sub, "VECTOR%1.2f_%1.2f_%1.2f" % tuple(v)


# ---------------------------------------------------------------------------
# Drawing helpers
# ---------------------------------------------------------------------------

def _split_near_far(v_obs):
    """Screen coordinates of a curve split into near-side and far-side pieces."""
    x, y = _to_screen(v_obs)
    near = v_obs[2] >= 0
    xn, yn = np.where(near, x, np.nan), np.where(near, y, np.nan)
    xf, yf = np.where(~near, x, np.nan), np.where(~near, y, np.nan)
    return (xn, yn), (xf, yf)


def _circle(theta_c, phi_c, zeta, n=721):
    """Points (3, n) on the small circle of angular radius zeta about (theta_c, phi_c)."""
    a = np.linspace(0, 2 * np.pi, n)
    V = np.array([np.sin(zeta) * np.cos(a), np.sin(zeta) * np.sin(a),
                  np.cos(zeta) * np.ones_like(a)])
    return _rotation(theta_c, phi_c) @ V


def _draw_curve(ax, M, V, color, lw, alpha_near, alpha_far, ls="-", z=3, label=None):
    (xn, yn), (xf, yf) = _split_near_far(M @ V)
    if alpha_far > 0:
        ax.plot(xf, yf, ls, color=color, lw=lw, alpha=alpha_far, zorder=z - 2,
                solid_capstyle="round")
    ax.plot(xn, yn, ls, color=color, lw=lw, alpha=alpha_near, zorder=z,
            solid_capstyle="round", label=label)


# ---------------------------------------------------------------------------
# Main function
# ---------------------------------------------------------------------------

def plot_projection_filled(dictVp, model, POV="", ThetaDisplay="", antiphase=False,
                           SaveFlag=False, dir="", Name="", POVname="",
                           extension=".png", antipodal=False,
                           Tmin=None, Tmax=None, min_span=0.3, cmap=None, style="dark",
                           show_far_side=True, far_alpha=0.45, far_style="hatch", far_legend=True,
                           show_outlines=False, grid_step=15.0, resolution=900,
                           colorbar=True, ax=None, figsize=None, back_view=False,
                           _M=None, _title=None, _warn=True):
    """ Plot a filled, temperature-coloured 2D projection of the hot region(s).

    :param dictVp: Parameter values (or an X-PSI ParameterSubspace).
    :param str model: Model string, e.g. ``"ST"``, ``"PDT"``, ``"ST+PST"``,
        ``"CST-S"``, ``"ST-Ua"`` or ``"ST-U"``.
    :param POV: Point of view: ``""``, ``"I"``, a component code (``"P"``,
        ``"PE"``, ``"PO"``, ``"PS"``, ``"PC"``, ``"S"``, ...) or a 3-vector
        pointing from the stellar centre to the viewer.
    :param str ThetaDisplay: Accepted for compatibility; the latitude/longitude
        grid is always drawn in the stellar frame instead.
    :param bool antiphase: If True, the secondary phase is defined 0.5 cycles
        from the primary's origin.
    :param bool SaveFlag: Whether to save the figure to file.
    :param str dir: Directory in which to save the figure.
    :param str Name: Extra tag appended to the saved file name.
    :param str POVname: Name used in the saved file name when ``POV`` is a
        vector rather than a recognised string.
    :param str extension: Extension of the saved file, e.g. ``".png"``.
    :param bool antipodal: If True, also draw (dashed) the outline of the
        antipode of the primary.
    :param Tmin: Lower log10(T/K) limit of the colour scale (default: span of
        the emitting temperatures in this model). Fix ``Tmin`` and ``Tmax`` to
        compare several figures on the same scale.
    :param Tmax: Upper log10(T/K) limit of the colour scale (see ``Tmin``).
    :param float min_span: Minimum width [dex] of the automatic colour scale,
        so that e.g. 6.10 vs 6.13 is not stretched across the whole map.
    :param cmap: Matplotlib colormap or name (default: inferno; for the dark
        style its darkest 20 % is cut off). ``HOT_GLOW_CMAP`` gives the
        orange/magenta look.
    :param str style: ``"dark"`` (wireframe on a night-sky background) or
        ``"light"``.
    :param bool show_far_side: If True, also show (dimmed) what lies on the far
        hemisphere.
    :param float far_alpha: Opacity of the far-hemisphere layer.
    :param str far_style: ``"hatch"`` (default): far-side emission is drawn with
        diagonal stripes and a dashed outline; ``"dim"``: just dimmed.
    :param bool far_legend: If True, add a small legend explaining near vs far
        side.
    :param bool show_outlines: If True, draw thin outlines of every component,
        including omission regions (dashed) and the masked part of ceding
        regions.
    :param float grid_step: Spacing of the latitude/longitude wireframe [deg].
    :param int resolution: Number of pixels across the stellar disc.
    :param bool colorbar: If True, add a log10(T) colour bar.
    :param ax: Existing axes to draw into (a new figure is made otherwise).
    :param figsize: Figure size (default ``(8.6, 8)``, or ``(16, 8)`` with
        ``back_view``).
    :param bool back_view: If True, add a second panel showing the star from the
        opposite side, as if you walked around it (north stays up). Spots on the
        far side in the main panel appear solid there.

    :returns: ``ax`` (or ``(ax_front, ax_back)`` when ``back_view=True``).
    """
    st = STYLES[style]
    if cmap is None:
        cmap = INFERNO_DARK_CMAP if style == "dark" else plt.get_cmap("inferno")
    elif isinstance(cmap, str):
        cmap = plt.get_cmap(cmap)

    if back_view:
        return _front_back(dictVp, model, POV, antiphase, SaveFlag, dir, Name, POVname,
                           extension, antipodal, Tmin, Tmax, min_span, cmap, style,
                           show_far_side, far_alpha, far_style, far_legend,
                           show_outlines, grid_step, resolution, colorbar,
                           figsize or (16, 8))
    if figsize is None:
        figsize = (8.6, 8)

    regions = parse_hot_regions(dictVp, model, antiphase)
    theta_sub, phi_sub, labPOV = _pov_angles(POV, dictVp, regions)
    M = _view_matrix(theta_sub, phi_sub)      # stellar -> viewer

    if _M is not None:
        M = _M
    norm = _temp_norm(regions, Tmin, Tmax, min_span)
    lo = norm.vmin

    # ---- render the surface ------------------------------------------------
    N = int(resolution)
    pix = (np.arange(N) + 0.5) / N * 2.0 - 1.0
    X, Y = np.meshgrid(pix, pix[::-1])          # screen coordinates
    r2 = X ** 2 + Y ** 2
    on_disk = r2 <= 1.0
    zf = np.sqrt(np.clip(1.0 - r2, 0.0, 1.0))

    def layer(sign, warn):
        v_obs = np.array([-Y, X, sign * zf])     # inverse of _to_screen
        n = np.tensordot(M.T, v_obs, axes=1)     # viewer -> stellar
        logT = surface_log_temperature(n, regions, warn=warn)
        logT[~on_disk] = np.nan
        rgba = cmap(norm(np.nan_to_num(logT, nan=lo)))
        rgba[..., 3] = np.isfinite(logT).astype(float)
        return rgba

    if _warn:
        _check_overlap(regions)
    near = layer(+1.0, warn=False)
    far = layer(-1.0, warn=False) if show_far_side else None

    # Soft anti-aliasing of the disc edge for the filled layers
    edge = np.clip((1.0 - np.sqrt(r2)) * N / 2.0 + 0.5, 0.0, 1.0)
    near[..., 3] *= edge

    # ---- figure ------------------------------------------------------------
    own_fig = ax is None
    if own_fig:
        fig = plt.figure(figsize=figsize)
        ax = fig.add_axes([0.0, 0.0, 0.86 if colorbar else 1.0, 1.0])
    fig = ax.figure
    ext = [-1, 1, -1, 1]
    L = 1.18
    ax.set_xlim(-L, L)
    ax.set_ylim(-L, L)
    ax.set_aspect("equal")
    ax.axis("off")

    # background gradient
    if style == "dark":
        g = np.linspace(0, 1, 256)[:, None] * np.ones((1, 2))
        bg = LinearSegmentedColormap.from_list("bg", [st["bg_bottom"], st["bg_top"]])
        ax.imshow(g, cmap=bg, extent=[-L, L, -L, L], origin="lower",
                  aspect="auto", zorder=-10)
        fig.patch.set_facecolor(st["bg_bottom"])
    else:
        fig.patch.set_facecolor("white")

    # far hemisphere: emission, then wireframe
    if far is not None:
        far_mask = far[..., 3] > 0
        if far_style == "hatch":
            # diagonal stripes in screen space: strong on the stripe, faint between
            period = max(N / 60.0, 4.0)
            stripe = ((X + Y) * N / 2.0 / period) % 1.0 < 0.45
            far[..., 3] *= np.where(stripe, min(1.0, far_alpha * 1.8), far_alpha * 0.35)
        else:
            far[..., 3] *= far_alpha
        ax.imshow(far, extent=ext, origin="upper", interpolation="antialiased", zorder=1)
        if far_style == "hatch" and far_mask.any():
            ax.contour(X, Y, far_mask.astype(float), levels=[0.5], colors=[st["text"]],
                       linewidths=1.0, linestyles="--", alpha=0.8, zorder=1.5)

    # grid in the stellar frame
    step = np.deg2rad(grid_step)
    a = np.linspace(0, 2 * np.pi, 721)
    grid_far = 0.18 if show_far_side else 0.0
    for colat in np.arange(step, np.pi - 1e-6, step):
        V = np.array([np.sin(colat) * np.cos(a), np.sin(colat) * np.sin(a),
                      np.cos(colat) * np.ones_like(a)])
        is_eq = abs(colat - np.pi / 2) < 1e-6
        _draw_curve(ax, M, V, st["equator"] if is_eq else st["grid"],
                    1.6 if is_eq else 0.7, 0.95 if is_eq else 0.75, grid_far, z=4,
                    label="Equator" if is_eq else None)
    t = np.linspace(0, np.pi, 361)
    for lon in np.arange(0, 2 * np.pi - 1e-6, step):
        V = np.array([np.sin(t) * np.cos(lon), np.sin(t) * np.sin(lon), np.cos(t)])
        _draw_curve(ax, M, V, st["grid"], 1.3 if lon == 0 else 0.7,
                    0.95 if lon == 0 else 0.75, grid_far, z=4)

    # translucent stellar disc between the far and near layers
    ax.add_patch(plt.Circle((0, 0), 1.0, color=st["disk"], alpha=st["disk_alpha"],
                            zorder=2, lw=0))

    # near hemisphere emission 
    ax.imshow(near, extent=ext, origin="upper", interpolation="antialiased", zorder=6)

    # limb
    if style == "dark":
        for lw, al in ((8, 0.08), (4, 0.18)):
            ax.add_patch(plt.Circle((0, 0), 1.0, fill=False, color=st["limb"],
                                    lw=lw, alpha=al, zorder=7))
    ax.add_patch(plt.Circle((0, 0), 1.0, fill=False, color=st["limb"], lw=1.6, zorder=7))

    # outlines of components (optional)
    if show_outlines:
        for r in regions:
            for key, ls in (("super", "-"), ("cede", "-"), ("omit", "--")):
                comp = r[key]
                if comp is None:
                    continue
                col = "white" if key == "omit" else cmap(norm(comp[3]))
                if style == "light" and key == "omit":
                    col = "black"
                _draw_curve(ax, M, _circle(*comp[:3]), col, 1.0, 0.9,
                            far_alpha * 0.8 if show_far_side else 0.0, ls=ls, z=8)
    if antipodal:
        prim = regions[0]
        for key in ("super", "cede", "omit"):
            comp = _antipode(prim[key])
            if comp is not None:
                _draw_curve(ax, M, _circle(*comp[:3]), st["text"], 1.2, 0.9,
                            far_alpha * 0.8 if show_far_side else 0.0, ls="-.", z=8,
                            label="Antipode of primary" if key == "super" else None)

    # poles and phase-zero marker
    for vec, lab in (([0, 0, 1], "N"), ([0, 0, -1], "S")):
        vo = M @ np.array(vec, float)
        x, y = _to_screen(vo)
        visible = vo[2] >= 0
        if visible or show_far_side:
            ax.plot(x, y, marker="*", ms=11 if visible else 7, color=st["pole"],
                    alpha=1.0 if visible else far_alpha, zorder=9, ls="none")
            ax.annotate(lab, (x, y), xytext=(6, 6), textcoords="offset points",
                        color=st["pole"], alpha=1.0 if visible else far_alpha,
                        fontsize=12, zorder=9)
    vo = M @ _unit(np.pi / 2, 0.0)
    x, y = _to_screen(vo)
    if vo[2] >= 0 or show_far_side:
        ax.plot(x, y, marker="o", ms=6, mfc="none", mec=st["pole"],
                alpha=1.0 if vo[2] >= 0 else far_alpha, zorder=9, ls="none")
        ax.annotate(r"$\phi=0$", (x, y), xytext=(6, -14), textcoords="offset points",
                    color=st["pole"], alpha=1.0 if vo[2] >= 0 else far_alpha,
                    fontsize=11, zorder=9)

    # colour bar
    if colorbar and own_fig:
        _add_colorbar(fig, fig.add_axes([0.87, 0.2, 0.025, 0.6]), norm, cmap, st, regions)
    elif colorbar:
        cb = fig.colorbar(ScalarMappable(norm=norm, cmap=cmap), ax=ax, fraction=0.04)
        cb.set_label(r"$\log_{10}(T\;[\mathrm{K}])$", color=st["text"])
        cb.ax.tick_params(colors=st["text"])

    ax.set_xlim(-L, L)
    ax.set_ylim(-L, L)
    ax.set_aspect("equal")
    if show_far_side and far_legend and far is not None and far[..., 3].any():
        mid = cmap(0.6)
        handles = [Patch(facecolor=mid, edgecolor="none", label="Near side"),
                   Patch(facecolor=(*mid[:3], far_alpha * 0.5), edgecolor=st["text"],
                         hatch="////" if far_style == "hatch" else None, ls="--",
                         label="Far side (seen through the star)")]
        leg = ax.legend(handles=handles, loc="upper left", fontsize=10, frameon=False, borderaxespad=0.2,
                        labelcolor=st["text"], handlelength=2.2, handleheight=1.2)
        leg.set_zorder(20)
    ax.set_title(_title or "%s   (view: %s)" % (model, labPOV), color=st["text"], fontsize=14)

    if SaveFlag:
        _save(fig, dir, POV, POVname, labPOV, Name, extension)

    return ax


def _temp_norm(regions, Tmin, Tmax, min_span):
    temps = []
    for r in regions:
        temps.append(r["super"][3])
        if r["cede"] is not None:
            temps.append(r["cede"][3])
    lo = min(temps) if Tmin is None else Tmin
    hi = max(temps) if Tmax is None else Tmax
    if hi - lo < min_span:   # avoid exaggerating tiny temperature differences
        mid = 0.5 * (hi + lo)
        lo = mid - 0.5 * min_span if Tmin is None else lo
        hi = mid + 0.5 * min_span if Tmax is None else hi
    return Normalize(lo, hi)


def _add_colorbar(fig, cax, norm, cmap, st, regions):
    cb = fig.colorbar(ScalarMappable(norm=norm, cmap=cmap), cax=cax)
    cb.set_label(r"$\log_{10}(T\;[\mathrm{K}])$", color=st["text"], fontsize=14)
    cb.ax.tick_params(colors=st["text"], labelsize=11)
    cb.outline.set_edgecolor(st["text"])
    for r in regions:
        for key in ("super", "cede"):
            if r[key] is not None:
                cb.ax.axhline(r[key][3], color=st["text"], lw=1.0)
    return cb


def _save(fig, dir, POV, POVname, labPOV, Name, extension):
    if dir and not dir.endswith("/"):
        dir = dir + "/"
    FigName = dir + "ProjectionFilled"
    if isinstance(POV, str):
        FigName += "_from_%s" % labPOV
    elif POVname.strip():
        FigName += "_from_" + POVname
    else:
        raise RuntimeError("Please provide POVname to save file with.")
    FigName += ("_" + Name if Name else "") + extension
    fig.savefig(FigName, dpi=300, bbox_inches="tight", facecolor=fig.get_facecolor())


def _front_back(dictVp, model, POV, antiphase, SaveFlag, dir, Name, POVname, extension,
                antipodal, Tmin, Tmax, min_span, cmap, style, show_far_side, far_alpha,
                far_style, far_legend, show_outlines, grid_step, resolution,
                colorbar, figsize):
    """Main view plus the view from the opposite side, sharing one colour scale."""
    st = STYLES[style]
    if cmap is None:
        cmap = INFERNO_DARK_CMAP if style == "dark" else plt.get_cmap("inferno")
    elif isinstance(cmap, str):
        cmap = plt.get_cmap(cmap)
    regions = parse_hot_regions(dictVp, model, antiphase)
    theta_sub, phi_sub, labPOV = _pov_angles(POV, dictVp, regions)
    M = _view_matrix(theta_sub, phi_sub)
    # walk around the star: rotate the viewer by 180 deg about the screen's
    # vertical axis (viewer x axis), so "up" on the page stays "up"
    M_back = np.diag([1.0, -1.0, -1.0]) @ M
    norm = _temp_norm(regions, Tmin, Tmax, min_span)

    fig = plt.figure(figsize=figsize)
    w = 0.44 if colorbar else 0.5
    ax_f = fig.add_axes([0.0, 0.0, w, 0.94])
    ax_b = fig.add_axes([w, 0.0, w, 0.94])
    common = dict(antiphase=antiphase, antipodal=antipodal, Tmin=norm.vmin, Tmax=norm.vmax,
                  min_span=0.0, cmap=cmap, style=style, far_alpha=far_alpha,
                  far_style=far_style, show_outlines=show_outlines,
                  grid_step=grid_step, resolution=resolution, colorbar=False)
    plot_projection_filled(dictVp, model, POV, show_far_side=show_far_side,
                           far_legend=far_legend, ax=ax_f,
                           _title="Front  (view: %s)" % labPOV, **common)
    plot_projection_filled(dictVp, model, POV, show_far_side=False, ax=ax_b, _M=M_back, _warn=False,
                           _title="Back  (opposite side, north still up)", **common)
    fig.suptitle(model, color=st["text"], fontsize=16)
    if colorbar:
        _add_colorbar(fig, fig.add_axes([0.9, 0.2, 0.014, 0.56]), norm, cmap, st, regions)
    if SaveFlag:
        _save(fig, dir, POV, POVname, labPOV, (Name + "_" if Name else "") + "frontback",
              extension)
    return ax_f, ax_b


# ---------------------------------------------------------------------------
# Animation
# ---------------------------------------------------------------------------

def animate_projection_filled(dictVp, model, filename="rotation.mp4", n_frames=72,
                              fps=24, n_rotations=1, inclination=None, antiphase=False,
                              dpi=100, crf=20, show_phase=True, **kwargs):
    """
    Write a movie of the star rotating, as seen by the observer.

    The observer is fixed at colatitude ``inclination`` (by default
    ``arccos(cos_inclination)`` from ``dictVp``) and the star spins
    prograde about its north pole, which stays at the top of the frame.
    Frame ``k`` is at rotational phase ``t = k * n_rotations / n_frames``
    [cycles], with the X-PSI convention: at phase ``t`` the observer is at
    azimuth ``-t`` in the stellar frame, so a hot region with
    ``phase_shift = s`` crosses the observer's meridian at ``t = -s`` (mod 1),
    as in X-PSI pulse profiles (ignoring the small offsets from e.g. light bending,
    time delays and Doppler boosting, which this projection does not include).

    :param dictVp: Parameter values (or an X-PSI ParameterSubspace), as for
        :func:`plot_projection_filled`.
    :param str model: Model string, as for :func:`plot_projection_filled`.
    :param str filename: Output file. ``.gif`` is written with Pillow (no extra
        software needed, plays inline almost everywhere, but larger files);
        anything else (e.g. ``.mp4``) is encoded as H.264 with ffmpeg, which
        must be installed (e.g. ``conda install ffmpeg``).
    :param int n_frames: Number of frames (per movie, not per rotation).
    :param int fps: Frames per second.
    :param int n_rotations: Number of rotations shown; the movie loops
        seamlessly.
    :param inclination: Observer colatitude [rad]; overrides
        ``cos_inclination`` / ``sin_inclination``.
    :param bool antiphase: As for :func:`plot_projection_filled`.
    :param int dpi: Frame resolution (figure size times dpi gives the pixel
        size).
    :param int crf: H.264 quality for MP4 output (lower = better, larger file).
    :param bool show_phase: If True, put the current rotational phase in the
        title.
    :param kwargs: Passed on to :func:`plot_projection_filled` (e.g. ``Tmin``,
        ``Tmax``, ``cmap``, ``style``, ``show_outlines``, ``far_style``,
        ``resolution``, ``figsize``). ``POV``, ``back_view``, ``ax`` and
        ``SaveFlag`` are not available here.

    :returns: ``filename``.
    """

    for bad in ("POV", "back_view", "ax", "SaveFlag", "_M", "_title", "_warn"):
        if bad in kwargs:
            raise TypeError("'%s' cannot be used with animate_projection_filled." % bad)

    regions = parse_hot_regions(dictVp, model, antiphase)
    _check_overlap(regions)
    if inclination is None:
        theta_i, _, _ = _pov_angles("I", dictVp, regions)
    else:
        theta_i = float(inclination)

    # same colour scale in every frame
    norm = _temp_norm(regions, kwargs.pop("Tmin", None), kwargs.pop("Tmax", None),
                      kwargs.pop("min_span", 0.3))
    kwargs.setdefault("resolution", 600)
    kwargs.setdefault("far_legend", False)

    gif = filename.lower().endswith(".gif")
    if not gif:
        ffmpeg = shutil.which("ffmpeg") or plt.rcParams.get("animation.ffmpeg_path")
        if not ffmpeg or not shutil.which(ffmpeg):
            raise RuntimeError("ffmpeg was not found. Install it (e.g. 'conda install ffmpeg') "
                               "or give a filename ending in '.gif'.")

    frames, proc = [], None
    try:
        for k in range(n_frames):
            t = k * n_rotations / float(n_frames)
            M = _view_matrix(theta_i, -t)
            title = model + ("   (phase %.2f)" % (t % 1.0) if show_phase else "")
            ax = plot_projection_filled(dictVp, model, antiphase=antiphase, Tmin=norm.vmin,
                                        Tmax=norm.vmax, min_span=0.0, _M=M, _title=title, _warn=False,
                                        **kwargs)
            fig = ax.figure
            fig.set_dpi(dpi)
            fig.canvas.draw()
            img = np.asarray(fig.canvas.buffer_rgba())[..., :3]
            plt.close(fig)
            h, w = img.shape[0] - img.shape[0] % 2, img.shape[1] - img.shape[1] % 2
            img = np.ascontiguousarray(img[:h, :w])      # H.264 needs even sizes
            if gif:
                frames.append(img)
                continue
            if proc is None:
                proc = subprocess.Popen(
                    [ffmpeg, "-y", "-loglevel", "error", "-f", "rawvideo",
                     "-pix_fmt", "rgb24", "-s", "%dx%d" % (w, h), "-r", str(fps),
                     "-i", "-", "-c:v", "libx264", "-preset", "veryfast",
                     "-crf", str(crf), "-pix_fmt", "yuv420p",
                     "-movflags", "+faststart", filename],
                    stdin=subprocess.PIPE)
            proc.stdin.write(img.tobytes())
    finally:
        if proc is not None:
            proc.stdin.close()
            if proc.wait() != 0:
                raise RuntimeError("ffmpeg failed to write '%s'." % filename)

    if gif:
        from PIL import Image
        ims = [Image.fromarray(f) for f in frames]
        ims[0].save(filename, save_all=True, append_images=ims[1:],
                    duration=int(round(1000.0 / fps)), loop=0, optimize=True)
    return filename
