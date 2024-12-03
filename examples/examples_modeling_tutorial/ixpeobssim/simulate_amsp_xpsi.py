#!/usr/bin/env python
#
# Copyright (C) 2016--2018, the ixpeobssim team.
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

#To specify the directory where output is saved when running this script:
#export IXPEOBSSIM_DATA=/home/xiaotuo/ixpeobssimdata/ixpeobssimdata_24X

from __future__ import print_function, division
import numpy

import ixpeobssim.config.model_amsp_xpsi as input_model

import ixpeobssim.core.pipeline as pipeline
from ixpeobssim.binning.base import xEventBinningBase
from ixpeobssim.binning.misc import xBinnedPulseProfile
from ixpeobssim.binning.polarization import xBinnedPolarizationCube
from ixpeobssim.utils.misc import pairwise_enum
from ixpeobssim.utils.matplotlib_ import plt, setup_gca, last_line_color
from ixpeobssim.utils.fmtaxis import fmtaxis

DURATION = 600000.
source_name = list(input_model.ROI_MODEL)[0]
PHASE_BINNING = numpy.linspace(0., 1., 11)
ENERGY_BINNING = numpy.array([2., 8.])

def simulate():
    """Run the simulation and fold the events in phase.
    """
    file_list = pipeline.xpobssim(duration=DURATION, startdate=input_model.start_date,
                                  seed=0, deadtime=0.0)#, occult=True, saa=True)
    pipeline.xpphase(*file_list, suffix='folded', **input_model.ephemeris.dict())

def select():
    """Run xpselect in a number of bins in phase.
    """
    file_list = pipeline.file_list('folded')
    for i, (min_, max_) in pairwise_enum(PHASE_BINNING):
        pipeline.xpselect(*file_list, phasemin=min_, phasemax=max_,
                          suffix=pipeline.suffix('phase', i))

def bin_():
    """Create a pulse profile, as well as modulation a modulation for each
    subselection in phase.
    """
    # Pulse profile.
    file_list = pipeline.file_list('folded')
    pipeline.xpbin(*file_list, algorithm='PP')
    # Modulation cubes.
    for i, (min_, max_) in pairwise_enum(PHASE_BINNING):
        file_list = pipeline.file_list('folded', ('phase', i))
        pipeline.xpbin(*file_list, algorithm='PCUBE', ebinalg='LIST',
                       ebinning=ENERGY_BINNING)

def display_pulse_profile():
    """Display the pulse profile.
    """
    shape = (len(ENERGY_BINNING) - 1, len(PHASE_BINNING) - 1)
    emean = numpy.zeros(shape)
    for i, (min_, max_) in pairwise_enum(PHASE_BINNING):
        file_list = pipeline.file_list('folded', ('phase', i), 'pcube')
        pcube = xBinnedPolarizationCube.from_file_list(file_list)
        emean[:,i] = pcube.E_MEAN
    mean_energy = numpy.mean(emean)

    pipeline.figure('pulse profile')
    file_list = pipeline.file_list('folded_pp')
    pp = xBinnedPulseProfile.from_file_list(file_list)
    pp.plot(label='IXPE %d ks' % pp.ontime())
    phase = numpy.linspace(0., 1., 100)
    model_pp = input_model.spec(mean_energy,phase)
    scale = numpy.mean(pp.COUNTS) / numpy.mean(model_pp)
    plt.plot(phase, scale * model_pp, label='Input model @ %.2f keV' % mean_energy)
    setup_gca(ymin=0, legend=True)


def display_pol_degree():
    """Display the polarization degree as a function of the pulse phase.
    """
    phase = numpy.linspace(0., 1., 100)
    phase_bins = xEventBinningBase.bin_centers(PHASE_BINNING)
    shape = (len(ENERGY_BINNING) - 1, len(PHASE_BINNING) - 1)
    pol_deg = numpy.zeros(shape)
    pol_deg_err = numpy.zeros(shape)
    pol_ang = numpy.zeros(shape)
    pol_ang_err = numpy.zeros(shape)
    MDP = numpy.zeros(shape)
    emean = numpy.zeros(shape)
    for i, (min_, max_) in pairwise_enum(PHASE_BINNING):
        file_list = pipeline.file_list('folded', ('phase', i), 'pcube')
        pcube = xBinnedPolarizationCube.from_file_list(file_list)
        pol_deg[:,i] = pcube.PD
        pol_deg_err[:,i] = pcube.PD_ERR
        pol_ang[:,i] = pcube.PA
        pol_ang_err[:,i] = pcube.PA_ERR
        MDP[:,i] = pcube.MDP_99
        emean[:,i] = pcube.E_MEAN

    def data_label(emin, emax):
        return 'IXPE %d ks (%.2f - %.2f keV)' % (pcube.ontime(), emin, emax)

    def model_label(mean_energy):
        return 'Input model @ %.2f keV' % mean_energy

    pipeline.figure('polarization degree')
    for i, (min_, max_) in pairwise_enum(ENERGY_BINNING):
        plt.errorbar(phase_bins, pol_deg[i,:], pol_deg_err[i,:], fmt='o',
                     label=data_label(min_, max_))
        energy = numpy.mean(emean[i,:])
        plt.plot(phase, input_model.pol_deg(energy, phase),
                 color=last_line_color(), label=model_label(energy))
        plt.errorbar(phase_bins, MDP[i,:], 0.0*MDP[i,:], fmt='o',color='red',
                     label="MDP")
    setup_gca(ymin=0, ymax=0.2, legend=True, **fmtaxis.pp_pol_deg)

    pipeline.figure('polarization angle')
    for i, (min_, max_) in pairwise_enum(ENERGY_BINNING):
        plt.errorbar(phase_bins, pol_ang[i,:], pol_ang_err[i,:], fmt='o',
                     label=data_label(min_, max_))
        energy = numpy.mean(emean[i,:])
        plt.plot(phase, numpy.degrees(input_model.pol_ang(energy, phase)),
                 color=last_line_color(), label=model_label(energy))
    setup_gca(ymin=-180.0, ymax=180.0, legend=True, **fmtaxis.pp_pol_ang)


def display():
    """Display all.
    """
    display_pulse_profile()
    display_pol_degree()


def run():
    """Run all.
    """
    simulate()
    select()
    bin_()
    display()


if __name__ == '__main__':
    pipeline.bootstrap_pipeline('model_amsp_xpsi')
