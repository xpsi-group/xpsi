import numpy as np
import matplotlib.pyplot as plt
import h5py as h5
from matplotlib.ticker import NullFormatter
import sys
import scipy

from io import StringIO
import os
import math
import time

from matplotlib import rcParams
from matplotlib.ticker import MultipleLocator, AutoLocator, AutoMinorLocator
from matplotlib import gridspec
from matplotlib import cm

import xpsi

"""
This library contains plotting functions that are commonly used in the tutorial notebooks.
This is meant to reduce the length of notebook
"""

rcParams['text.usetex'] = False
rcParams['font.size'] = 14.0


def veneer(x, y, axes, lw=1.0, length=8, yticks=None):
    """ Make the plots a little more aesthetically pleasing. """
    if x is not None:
        if x[1] is not None:
            axes.xaxis.set_major_locator(MultipleLocator(x[1]))
        if x[0] is not None:
            axes.xaxis.set_minor_locator(MultipleLocator(x[0]))
    else:
        axes.xaxis.set_major_locator(AutoLocator())
        axes.xaxis.set_minor_locator(AutoMinorLocator())

    if y is not None:
        if y[1] is not None:
            axes.yaxis.set_major_locator(MultipleLocator(y[1]))
        if y[0] is not None:
            axes.yaxis.set_minor_locator(MultipleLocator(y[0]))
    else:
        axes.yaxis.set_major_locator(AutoLocator())
        axes.yaxis.set_minor_locator(AutoMinorLocator())

    axes.tick_params(which='major', colors='black', length=length, width=lw)
    axes.tick_params(which='minor', colors='black', length=int(length / 2), width=lw)
    plt.setp(axes.spines.values(), linewidth=lw, color='black')
    if yticks:
        axes.set_yticks(yticks)


def plot_1d_pulse(signal=None, photosphere=None, photosphere_phases=None):
    """ Plot hot region signals before and after telescope operation. """
    fig = plt.figure(figsize=(7, 7))
    ax = fig.add_subplot(111)

    ax.set_ylabel('Signal [arbitrary normalisation]')
    ax.set_xlabel('Phase [cycles]')

    if signal is not None:
        temp = np.sum(signal.signals[0], axis=0)
        ax.plot(signal.phases[0], temp / np.max(temp), '-', color='k', lw=0.5, label="Signal (spot 1)")
        temp = np.sum(signal.signals[1], axis=0)
        ax.plot(signal.phases[1], temp / np.max(temp), '-', color='r', lw=0.5, label="Signal (spot 2)")

    if photosphere is not None:
        if photosphere_phases is None:
            print("WARNING: Photosphere_phases is missing. Ignoring 'photosphere' in plot")
        else:
            temp = np.sum(photosphere.signal[0][0], axis=0)
            ax.plot(photosphere_phases[0], temp / np.max(temp),
                    'o-', color='k', lw=0.5, markersize=2,
                    label="Photosphere (spot 1)")
            temp = np.sum(photosphere.signal[1][0], axis=0)
            ax.plot(photosphere_phases[1], temp / np.max(temp), 'o-', color='r', lw=0.5, markersize=2,
                    label="Photosphere (spot 2)")

    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    veneer((0.05, 0.2), (0.05, 0.2), ax)


def plot_2d_pulse(pulse, x, y, label=r'Counts',
                  cmap=cm.magma, vmin=None, vmax=None,
                  rotations=2,
                  ylabel='',
                  ylog=True,
                  yticks=None,
                  cbar_label='',
                  normalized=False,
                  ):
    """ Plot a pulse resolved over a single rotational cycle. """

    if rotations > 1:
        fig = plt.figure(figsize=(14, 7))
    else:
        fig = plt.figure(figsize=(7, 7))

    gs = gridspec.GridSpec(1, 2, width_ratios=[50, 1])
    ax = plt.subplot(gs[0])
    ax_cb = plt.subplot(gs[1])

    # Calculate the centre of phase bins (required by pcolormesh instead of phase edges)
    if np.shape(x)[0] - np.shape(pulse)[-1] == 0:
        _phase_bincenter = x[:]
    elif np.shape(x)[0] - np.shape(pulse)[-1] == 1:
        _phase_bincenter = x[:-1] + 0.5 * (x[1] - x[0])
    else:
        raise Exception(f"Phase array dimension ({len(x)}) does not match pulse data (shape:{np.shape(pulse)})")

    # Replicate and append the phase bin centers for the number of rotations
    tmp = _phase_bincenter
    for r in range(1, rotations):
        _phase_bincenter = np.concatenate((_phase_bincenter, r + tmp))

    # Normalizing
    if normalized:
        pulse = pulse / np.max(pulse)

    # Replicate and append the pulse for the number of rotations
    new_pulse = np.concatenate((pulse, np.tile(pulse, rotations - 1)), axis=1)

    profile = ax.pcolormesh(_phase_bincenter,
                            y,
                            new_pulse,
                            vmin=vmin,
                            vmax=vmax,
                            cmap=cmap,
                            linewidth=0,
                            rasterized=True)

    profile.set_edgecolor('face')

    ax.set_xlim([0.0, rotations])
    if ylog:
        ax.set_yscale('log')
    ax.set_ylabel(ylabel)
    ax.set_xlabel(r'Phase')

    cb = plt.colorbar(profile,
                      cax=ax_cb)

    cb.set_label(label=cbar_label, labelpad=25)
    cb.solids.set_edgecolor('face')

    veneer((0.05, 0.2), (None, None), ax, yticks=yticks)

    plt.subplots_adjust(wspace=0.025)

    if yticks is not None:
        ax.set_yticklabels(yticks)

    cb.outline.set_linewidth(1.0)


def plot_samps(samps, x_idx, y_idx, xlabel, ylabel,
               s=1.0, color='k',
               extras=None, extras_label=None,
               **kwargs,
               ):
    """ Plot samples as 2D scatter plot. With extra lines if requested """
    fig = plt.figure(figsize=(7, 7))
    ax = fig.add_subplot(111)
    ax.scatter(samps[:, x_idx], samps[:, y_idx], s=s, color=color, **kwargs, label='Samples')
    veneer(None, None, ax)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)

    if extras is not None:
        for i, ex in enumerate(extras):
            if isinstance(extras_label, list):
                try:
                    label = extras_label[i]
                except IndexError:
                    label = None
            else:
                label = None
            ax.plot(ex[0], ex[1], label=label)

    ax.legend(loc='center left', bbox_to_anchor=(1.0, 0.5))

    return plt.gca()




def plot_instruments(instruments, xlabel="Energy interval", ylabel="Channel"):
    nb_instru = len(instruments)

    fig, axs = plt.subplots(nb_instru + 1, 1, layout="constrained", figsize=(7, 5 * (nb_instru)))

    # RESPONSE
    for i, inst in enumerate(instruments):
        ax1 = axs[i]
        veneer((25, 100), (50, 100), ax1)
        matrixplot = ax1.imshow(inst.matrix,
                                cmap=cm.viridis,
                                rasterized=True)
        ax1.set_xlabel('Energy interval')
        ax1.set_ylabel(ylabel)

        cbar = fig.colorbar(matrixplot, ax=ax1, shrink=1 / nb_instru, pad=0.01)
        cbar.set_label(r'[cm$^2\,\mathrm{count}\,/\,\mathrm{photon}$]')

    # EFFECTIVE AREA
    ax2 = axs[nb_instru]
    veneer((0.1, 0.5), (50, 250), ax2)

    for inst in instruments:
        ax2.plot((inst.energy_edges[:-1] + inst.energy_edges[1:]) / 2.0,
                 np.sum(inst.matrix, axis=0), label='NICER')  ## TODO:  change NICER by instrument.name
    ax2.legend()
    ax2.set_ylabel(r'Effective area [cm$^{2}$]')
    ax2.set_xlabel('Energy [keV]')


def plot_rmf(matrix,
             x,
             y,
             xlabel="Energy interval",
             ylabel="Channel"):

    fig = plt.figure(figsize=(14, 7))
    ax = fig.add_subplot(111)

    # RESPONSE
    im = ax.pcolormesh(x, y, matrix,
                       cmap=cm.viridis,
                       rasterized=True)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    cbar = fig.colorbar(im, ax=ax, shrink=1, pad=0.01)
    cbar.set_label(r'[cm$^2\,\mathrm{count}\,/\,\mathrm{photon}$]')

    return ax

def plot_arf(instruments,
             ):

    fig = plt.figure(figsize=(7, 7))
    ax = fig.add_subplot(111)

    for i, inst in enumerate(instruments):
        ax.plot((inst.energy_edges[:-1] + inst.energy_edges[1:]) / 2.0,
                np.sum(inst.matrix, axis=0), label=inst.name)  ## TODO:  change NICER by instrument.name

    ax.set_ylabel('Effective area [cm$^{2}$]')
    ax.set_xlabel('Energy [keV]')
    ax.legend(loc='best')
    return ax

def plot_meshes(regions,
                lines=True,
                primary_ticks=(1, 5),
                secondary_ticks=(1, 5)):
    """ Plot representations of the cached meshes.

    Note that the lowest colatitude row of elements is plotted as the
    lowest row, so colatitude increases along the y-axis and azimuth
    increaes along the x-axis. This could be considered as spatially
    inverted if we were looking at the star whilst being oriented
    such that "up" is in the spin direction.

    """

    nb_regions = len(regions)
    nb_temp = 1
    for r in regions:
        if r._HotRegion__cellArea[1] is not None:
            nb_temp = 2
    print(f"Plotting the meshes for {nb_regions} regions, with {nb_temp} temperatures")

    fig = plt.figure(figsize=(nb_regions * 5, 5 * nb_temp))
    width = [50] * nb_regions
    width.append(1)

    gs = gridspec.GridSpec(nb_temp, nb_regions + 1, width_ratios=width, wspace=0.2, hspace=0.2)

    for t in range(nb_temp):
        for r in range(nb_regions):
            ax = plt.subplot(gs[t, r])
            veneer(primary_ticks, primary_ticks, ax)
            z = regions[r]._HotRegion__cellArea[t] / np.max(regions[r]._HotRegion__cellArea[0])

            patches = plt.pcolormesh(z,
                                     vmin=np.min(z),
                                     vmax=np.max(z),
                                     cmap=cm.magma,
                                     linewidth=0.5 if lines else 0.0,
                                     rasterized=True,
                                     edgecolor='black')
            if t == 0:
                ax.set_title(f"Region {r+1}, Temperature 1")
            else:
                ax.set_title(f"Region {r+1}, Temperature 2")

        ax_cb = plt.subplot(gs[t, -1])
        cb = plt.colorbar(patches,
                          cax=ax_cb,
                          ticks=MultipleLocator(0.2))

        cb.set_label(label=r'cell area (normalised by maximum)', labelpad=25)
        cb.solids.set_edgecolor('face')

    # veneer((None, None), (0.05, None), ax_cb)
    cb.outline.set_linewidth(1.0)


def plot_spectrum(all_data,
                  all_labels=None,
                  thickness=None):
    fig = plt.figure(figsize=(10, 10))

    ax = fig.add_subplot(111)
    veneer((5, 25), (None, None), ax)

    for i, d in enumerate(all_data):

        if isinstance(all_labels, list):
            try:
                label = all_labels[i]
            except IndexError:
                label = None
        else:
            label = None

        if isinstance(thickness, list):
            try:
                th = thickness[i]
            except IndexError:
                th = None
        else:
            th = None

        if np.shape(d)[0] != 2:
            ax.step(np.arange(len(d)), d, label=label, lw=th)
        else:
            xvalue = np.arange(np.shape(d)[1])
            ax.fill_between(xvalue, d[0], d[1],
                            alpha=0.5,
                            step='pre',
                            color = 'k',
                            label=label)

    ax.legend()
    ax.set_yscale('log')
    ax.set_ylabel('Counts')
    _ = ax.set_xlabel('Channel')

