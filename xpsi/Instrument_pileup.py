""" Instrument module for X-PSI modelling. Includes loading of any instrument's response."""

import os
import numpy as np
from astropy.io import fits
from astropy.table import Table

import xpsi

from xpsi import Parameter
from xpsi.utils import make_verbose
from xpsi.Instrument import ResponseError

from PileupModule import XrayPileup

class Instrument_pileup(xpsi.Instrument):

    """ Instrument with pile-up computation """

    def __call__(self, signal, *args):
        """ Overwrite. """
        ## Compute spectrum with pileup
        piled_spectrum = self.pileup.analyze(signal,
                                             alpha=self['grade_migration'],
                                             psf_frac=self['psf_fraction'])

        self._cached_signal = piled_spectrum
        return self._cached_signal

    @classmethod
    @make_verbose('Loading instrument response matrix',
                  'Response matrix loaded')
    def from_ogip_fits(cls,
              Data_path,
              ARF_path,
              RMF_path,
              min_channel=0,
              max_channel=-1,
              min_input=1,
              max_input=-1,
              bounds=dict(),
              values=dict(),
              datafolder=None,
              **kwargs):
        
        """ Load any instrument response matrix. """

        if datafolder:
            ARF_path = os.path.join( datafolder, ARF_path )
            RMF_path = os.path.join( datafolder, RMF_path )
            Data_path = os.path.join( datafolder, Data_path )

        # Open useful values in ARF/RMF    
        with fits.open( ARF_path ) as ARF_hdul:
            ARF_header = ARF_hdul['SPECRESP'].header
        ARF_instr = ARF_header['INSTRUME']
            
        with fits.open( RMF_path ) as RMF_hdul:
            RMF_header = RMF_hdul['MATRIX'].header
        RMF_instr = RMF_header['INSTRUME'] 
        DETCHANS = RMF_header['DETCHANS']
        NUMGRP = RMF_header['NAXIS2']
        TLMIN = RMF_header['TLMIN4']
        TLMAX = RMF_header['TLMAX4']

        # Get the values and change the -1 values if requried
        if max_channel == -1:
            max_channel = DETCHANS -1
        if max_input == -1:
            max_input = NUMGRP
        channels = np.arange( min_channel, max_channel+1 )
        inputs = np.arange( min_input, max_input+1  )

        # Perform routine checks
        assert ARF_instr == RMF_instr
        assert min_channel >= TLMIN and max_channel <= TLMAX
        assert min_input >= 0 and max_input <= NUMGRP

        # If everything in order, get the data
        with fits.open( RMF_path ) as RMF_hdul:
            RMF_MATRIX = RMF_hdul['MATRIX'].data
            RMF_EBOUNDS = RMF_hdul['EBOUNDS'].data

        # Get the channels from the data
        RMF = np.zeros((DETCHANS, NUMGRP))
        for i, (N_GRP, F_CHAN, N_CHAN, RMF_line) in enumerate( zip(RMF_MATRIX['N_GRP'], RMF_MATRIX['F_CHAN'], RMF_MATRIX['N_CHAN'], RMF_MATRIX['MATRIX']) ):

            # Skip if needed
            if N_GRP == 0:
                continue

            # Check the values
            if not isinstance(F_CHAN, np.ndarray ):
                F_CHAN = [F_CHAN]
                N_CHAN = [N_CHAN]

            # Add the values to the RMF
            n_skip = 0 
            for f_chan, n_chan in zip(F_CHAN,N_CHAN):

                if n_chan == 0:
                    continue

                RMF[f_chan:f_chan+n_chan,i] += RMF_line[n_skip:n_skip+n_chan]
                n_skip += n_chan

        # Make the RSP
        ARF = Table.read(ARF_path, 'SPECRESP')
        ARF_area = ARF['SPECRESP']

        # Extract the required matrix
        RSP = RMF * ARF_area
        RSP = RSP[min_channel:max_channel+1,min_input-1:max_input]

        # Find empty columns and lines
        empty_channels = np.all(RSP == 0, axis=1)
        empty_inputs = np.all(RSP == 0, axis=0)
        RSP = RSP[~empty_channels][:,~empty_inputs]
        channels = channels[ ~empty_channels ]
        inputs = inputs[ ~empty_inputs ]
        if empty_inputs.sum() > 0:
            print(f'Triming the response matrix because it contains lines with only 0 values.\n Now min_input={inputs[0]} and max_input={inputs[-1]}')
        if empty_channels.sum() > 0:
            print(f'Triming the response matrix because it contains columns with only 0 values.\n Now min_channel={channels[0]} and max_channel={channels[-1]}')

        # Get the edges of energies for both input and channel
        energy_edges = np.append( ARF['ENERG_LO'][inputs-1], ARF['ENERG_HI'][inputs[-1]-1]).astype(dtype=np.double)
        energies = (ARF['ENERG_LO']+ARF['ENERG_HI'][:])/2

        channel_energy_edges = np.append(RMF_EBOUNDS['E_MIN'][channels],RMF_EBOUNDS['E_MAX'][channels[-1]])
    

        ## Definition of the parameters for the pile-up
        alpha_grade = Parameter('grade_migration',
                    strict_bounds = (0.0,1.0),
                    bounds = bounds.get('grade_migration', None),
                    doc = 'Grade migration factor : probability that the piled event is not rejected as “bad event”',
                    symbol = r'$G_n$',
                    value = values.get('grade_migration', None))
    
        psffrac = Parameter('psf_fraction',
                    strict_bounds = (0.8,1.0),  #min in sherpa 0.85
                    bounds = bounds.get('psf_fraction', None),
                    doc = 'fraction of events in the source extraction region to which pileup will be applied',
                    symbol = r'$PSF_frac',
                    value = values.get('psf_fraction', 0.9))
    
        nregions = Parameter('nregions',
                    strict_bounds = (0.0,10.0),
                    bounds = None,     ##this parameter should always be fixed  - should be 1.0 for point sources
                    doc = 'number of regions to which pileup model will be applied independently',
                    symbol = r'$N_{regions}$',
                    value = values.get('nregions', 1.0))

        g0 = Parameter('g0',
                    strict_bounds = (0.0,1.0),
                    bounds = None,     ##this parameter should always be fixed
                    doc = 'grade correction for single photon detection',
                    symbol = r'$g_0$',
                    value = values.get('g0', 1.0))
    
        npiled = Parameter('npiled',
                    strict_bounds = (0.0,100.0),
                    bounds = None,     ##this parameter should always be fixed 
                    doc = 'number of photons considered for pileup in a single frame',
                    symbol = r'$N_{phot}$',
                    value = values.get('npiled', 5)) #or 30 in sherpa

        with fits.open( Data_path ) as hdul:
            Data_header = hdul['SPECTRUM'].header 

        ## need to add condition that checks that EXPTIME exists !! same for fracexpo
        
        frametime = Data_header['EXPTIME']

        frame = Parameter('frame_time',
                    strict_bounds = (0.0,10.0),
                    bounds = None,     ##this parameter should always be fixed
                    doc = 'good exposure time per frame',
                    symbol = r'$\tau$ (s)',
                    value = frametime )

        frac_expo = ARF_header['FRACEXPO']

        fracexpo = Parameter('frac_expo',
                    strict_bounds = (0.0,1.0),
                    bounds = None,    
                    doc = 'fraction of the frame exposure time to create effective frame exposure time',
                    symbol = r'$f_{expo}$',
                    value = frac_expo) 

        Instrument = cls(RSP,
                        energy_edges,
                        channels,
                        channel_energy_edges,
                        alpha_grade,psffrac,nregions,g0,npiled,frame,fracexpo,
                        **kwargs)
    
        # Add ARF and RMF
        Instrument.RMF = RMF[min_channel:max_channel+1,min_input-1:max_input][~empty_channels][:,~empty_inputs]
        Instrument.ARF = ARF_area[min_input-1:max_input][~empty_inputs]
        Instrument.energies = energies[min_input-1:max_input][~empty_inputs]

        ##Initialization of the pileup module
        pileup = XrayPileup(Instrument)
        Instrument.pileup = pileup

        return Instrument
    
