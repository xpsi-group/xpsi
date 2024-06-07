#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  7 15:48:20 2023

@author: bas
"""

import numpy as np
from xpsi.global_imports import  _keV, _k_B
k_B_over_keV = _k_B / _keV

def get_T_in_log10_Kelvin(T_keV):
  """
  Converts temperature from keV to log10(K) for a given input (scalar or tuple).

  Args:
      T_keV: The temperature in keV, can be a scalar or a tuple.

  Returns:
      A scalar or tuple containing the temperature in log10(K) for each input element.

  Raises:
      ValueError: If the input is not a scalar or a tuple.
  """

  if isinstance(T_keV, (int, float)):
    # Handle scalar case
    T_log10_Kelvin = np.log10(T_keV / k_B_over_keV)
    return T_log10_Kelvin
  elif isinstance(T_keV, tuple):
    # Handle tuple case
    T_log10_Kelvin_values = []
    for t in T_keV:
      T_log10_Kelvin_values.append(np.log10(t / k_B_over_keV))
    return tuple(T_log10_Kelvin_values)
  else:
    raise ValueError("Input must be a scalar or a tuple.")

def get_keV_from_log10_Kelvin(T_log10_Kelvin):
  """
  Converts temperature from log10(K) to keV for a given input (scalar or tuple).

  Args:
    T_log10_Kelvin: The temperature in log10(K), can be a scalar or a tuple.

  Returns:
    A scalar or tuple containing the temperature in keV for each input element.

  Raises:
    ValueError: If the input is not a scalar or a tuple.
  """

  if isinstance(T_log10_Kelvin, (int, float)):
    # Handle scalar case
    T_keV = np.power(10, T_log10_Kelvin) * k_B_over_keV
    return T_keV
  elif isinstance(T_log10_Kelvin, tuple):
    # Handle tuple case
    T_keV_values = []
    for t in T_log10_Kelvin:
      T_keV_values.append(np.power(10, t) * k_B_over_keV)
    return tuple(T_keV_values)
  else:
    raise ValueError("Input must be a scalar or a tuple.")

def get_mids_from_edges(edges):
    mids_len = len(edges)-1
    mids = np.empty(mids_len)
    for i in range(mids_len):
        mids[i] = (edges[i]+edges[i+1])/2
    return mids

from matplotlib import pyplot as plt

from matplotlib.ticker import MultipleLocator, AutoLocator, AutoMinorLocator
from matplotlib import gridspec
from matplotlib import cm
from xpsi.tools import phase_interpolator

def veneer(x, y, axes, lw=1.0, length=8):
    """
    Adjust the plot aesthetics for better appearance. (written by ChatGPT)
    
    This function adjusts the tick locations, tick lengths, tick widths, and
    spine widths of the given axes object(s) to enhance the appearance of the plot.
    
    Parameters:
    - x (tuple or None): Tuple containing the minor and major tick locator values
                        for the x-axis. If None, default tick locators are used.
    - y (tuple or None): Tuple containing the minor and major tick locator values
                        for the y-axis. If None, default tick locators are used.
    - axes (Axes or iterable of Axes): The axes object(s) to which the adjustments
                                       will be applied.
    - lw (float): Width of the axis spines.
    - length (float): Length of major tick marks in points.
    
    Note:
    - The x and y parameters should be tuples of the form (minor_locator, major_locator),
      where locator is a multiple of the data unit.
    - If only minor_locator is provided in x or y, the major locator will be set to None,
      resulting in default major tick locators.
    - If axes is an iterable, adjustments will be applied to each axes object in the iterable.
    """
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
    axes.tick_params(which='minor', colors='black', length=int(length/2), width=lw)
    plt.setp(axes.spines.values(), linewidth=lw, color='black')

def plot_2D_pulse(z, x, shift, y, ylabel,
                  num_rotations=1.0, res=1000, figsize=(5,3),
                  cm=cm.viridis, normalize=True):
    """ Helper function to plot a phase-energy pulse.

    :param array-like z:
        A pair of *ndarray[m,n]* objects representing the signal at
        *n* phases and *m* values of an energy variable.

    :param ndarray[n] x: Phases the signal is resolved at.

    :param tuple shift: Hot region phase parameters.

    :param ndarray[m] x: Energy values the signal is resolved at.

    """

    fig = plt.figure(figsize = figsize)

    gs = gridspec.GridSpec(1, 2, width_ratios=[50,1], wspace=0.025)
    ax = plt.subplot(gs[0])
    ax_cb = plt.subplot(gs[1])

    new_phase_edges = np.linspace(0.0, num_rotations, res+1)
    new_phase_mids = get_mids_from_edges(new_phase_edges)
    interpolated = phase_interpolator(new_phase_mids,
                                      x,
                                      z[0], shift[0])

    
    if len(z) == 2:
        interpolated += phase_interpolator(new_phase_mids,
                                            x,
                                            z[1], shift[1])
    if normalize:
        profile = ax.pcolormesh(new_phase_mids,
                                 y,
                                 interpolated/np.max(interpolated),
                                 cmap = cm,
                                 linewidth = 0,
                                 rasterized = True)
    elif not normalize:
        profile = ax.pcolormesh(new_phase_mids,
                                 y,
                                 interpolated,
                                 cmap = cm,
                                 linewidth = 0,
                                 rasterized = True,
                                 shading = 'flat')
    
    profile.set_edgecolor('face')

    ax.set_xlim([0.0, num_rotations])
    ax.set_yscale('log')
    ax.set_ylabel(ylabel)
    ax.set_xlabel(r'Phase')
    veneer((0.1, 0.5), (None,None), ax)
    
    if normalize:
        cb = plt.colorbar(profile, cax = ax_cb, ticks = MultipleLocator(0.2))
        cb.set_label(label=r'Signal (arbitrary units)', labelpad=25)
        cb.solids.set_edgecolor('face')
        veneer((None, None), (0.05, None), ax_cb)
        cb.outline.set_linewidth(1.0)
    elif not normalize:
        ticks = np.linspace(0, np.max(z[0]), 10)
        cb = plt.colorbar(profile, cax = ax_cb,  ticks = ticks)
    
    return fig, ax


def custom_subplots(*args, sharex=False, **kwargs):
    '''written by chatgpt'''
    fig, axes = plt.subplots(*args, sharex=sharex, **kwargs)
    if isinstance(axes, np.ndarray):
        for ax in axes.flatten():
            ax.__class__ = CustomAxes
    else:
        axes.__class__ = CustomAxes
    return fig, axes 

from matplotlib.axes import Axes

class CustomAxes(Axes):
    def plot_2D_counts(self, signal, phases, channels, label=r'Counts', cm=cm.jet):
        """ Plot a pulse resolved over a single rotational cycle. """

        if (signal < 0.0).any():
            vmax =  np.max( np.abs( signal ) )
            vmin = -vmax
        else:
            vmax = np.max(signal)
            vmin = np.min(signal)

        profile = self.pcolormesh(phases,
                                   channels,
                                   signal,
                                   cmap = cm,
                                   vmin = vmin,
                                   vmax = vmax,
                                   linewidth = 0,
                                   rasterized = True)

        profile.set_edgecolor('face')
        
        self.set_xlim([0.0, 1.0])
        self.set_yscale('log')
        self.set_ylabel(r'Energy (keV)')
        self.set_xlabel(r'Phase')



        #veneer((0.05, 0.2), (None, None), self)
        #plt.subplots_adjust(wspace = 1)#0.025)
        
        return profile

    def plot_pulse(self, phases_edges, my_data):
    
        # Plot the normalized data
        self.stairs(my_data, phases_edges)
    
        # Label the axes
        self.set_xlabel('Phase')
        # self.set_ylabel('Normalized Counts')
        
        # set lims
        self.set_ylim([np.min(my_data), np.max(my_data)])

    
    def plot_bolometric_pulse(self, phases_edges, my_data):
        """
        Plot a bolometric pulse.
    
        This function plots the bolometric pulse by summing the data along the
        second axis and normalizing it. It then labels the axes accordingly.
    
        Parameters:
        - phases_space (ndarray): Array representing the phases space.
        - my_data (ndarray): Array containing the data.
    
        Returns:
        - fig, ax: The matplotlib figure and axes objects.
        """
    
        # Calculate the mids from edges for phases_space
        phases_mids = get_mids_from_edges(phases_edges)
    
        # Sum the data along the second axis
        summed_data = np.sum(my_data, axis=0)
    
        # Normalize the summed data
        normalized_data = summed_data / np.max(summed_data)
    
        # Plot the normalized data
        self.plot(phases_mids, normalized_data)
    
        # Label the axes
        self.set_xlabel('Phase')
        # self.set_ylabel('Normalized Counts')


    def plot_2D_signal(self, z, x, shift, y, ylabel, num_rotations=1.0, 
                      res=1000, cm=cm.viridis, normalize=True):
        """ Helper function to plot a phase-energy pulse. """

        new_phase_edges = np.linspace(0.0, num_rotations, res+1)
        new_phase_mids = get_mids_from_edges(new_phase_edges)
        interpolated = phase_interpolator(new_phase_mids,
                                          x,
                                          z[0], shift[0])
    
        if len(z) == 2:
            interpolated += phase_interpolator(new_phase_mids,
                                                x,
                                                z[1], shift[1])
        if normalize:
            profile = self.pcolormesh(new_phase_mids,
                                      y,
                                      interpolated/np.max(interpolated),
                                      cmap = cm,
                                      linewidth = 0,
                                      rasterized = True)
        elif not normalize:
            profile = self.pcolormesh(new_phase_mids,
                                      y,
                                      interpolated,
                                      cmap = cm,
                                      linewidth = 0,
                                      rasterized = True,
                                      shading = 'flat')
    
        profile.set_edgecolor('face')
    
        self.set_xlim([0.0, num_rotations])
        self.set_yscale('log')
        self.set_ylabel(ylabel)
        self.set_xlabel(r'Phase')
        veneer((0.1, 0.5), (None,None), self)
        
        return profile
    
import xpsi
class SynthesiseData(xpsi.Data):
    """ Custom data container to enable synthesis. """

    def __init__(self, channels, phases, first, last):

        self.channels = channels
        # print(channels)
        # print(len(channels))
        self._phases = phases

        try:
            self._first = int(first)
            self._last = int(last)
        except TypeError:
            raise TypeError('The first and last channels must be integers.')
        if self._first >= self._last:
            raise ValueError('The first channel number must be lower than the '
                             'the last channel number.')
            
            
def shift_phase_data(phase_edges, bolometric_data, start_phase_index):
    """
    Shifts the phase data to start from a specified phase bin index.

    Parameters:
        phases_edges (array-like): Array of phase edges.
        bolometric_data (array-like): Bolometric data corresponding to the phase bins.
        start_phase_index (int): Index of the phase bin to start from.

    Returns:
        shifted_phases_edges (array): Shifted array of phase edges.
        shifted_bolometric_data (array): Shifted bolometric data.
    """
    num_phases = len(phase_edges)
    start_phase_index %= num_phases  # Ensure start_phase_index is within range

    shifted_phases_edges = np.concatenate([phase_edges[start_phase_index:], phase_edges[:start_phase_index] + 1])
    shifted_bolometric_data = np.concatenate([bolometric_data[start_phase_index:], bolometric_data[:start_phase_index]])

    return shifted_phases_edges, shifted_bolometric_data