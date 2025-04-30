import numpy as np
from scipy.fft import fft, ifft
import math
from typing import Tuple
import warnings

class XrayPileup:
    def __init__(self, Instrument):
        """
        Initialize pileup analysis with X-ray data files
        
        Parameters:
        -----------
        Instrument : xpsi.Instrument 
            Instrument object already loaded with the relevant parameters and ARFs/RMFs
        """
    
        # Setup initial parameters
        self.max_num_terms = int(Instrument['npiled'])
        self.num_points = len(Instrument.energies)
        self.frame_time = Instrument['frame_time']

        # Load data files
        self.rmf_data = Instrument.RMF 
        self.arf_data = Instrument.ARF
        self.energies = Instrument.energies
        self.instrument = Instrument


    def perform_pileup(self, spectrum: np.ndarray, alpha: float, psf_frac_init: float) -> np.ndarray:
        """     
            Perform pileup calculation
            Code mainly inspired from the XSPEC equivalent, based on Davis 2001
        """

        num_regions = self.instrument['nregions']
        g0 = self.instrument['g0']
        fracexpo = self.instrument['frac_expo']
        
        # ARF convolution
        arf_s = spectrum * self.arf_data   
        arf_s = np.maximum(arf_s, 0)
        
        # This is an offset we will require in the convolution in case energy_ar(0) is non-zero. 
        # Note that we assume the energy array is linear and starts at a multiple of the energy bin size.
        ioff = -int(self.energies[0] // (self.energies[1] - self.energies[0])) 
        
        # Calculate pileup following original algorithm
        psf_frac = psf_frac_init / num_regions / fracexpo
        arf_s_tmp = arf_s* psf_frac   
        integ_arf_s = np.sum(arf_s_tmp)  
        
        results = arf_s * psf_frac  #term p=1 of the sum

        if integ_arf_s == 0:
            return results
        
        exp_factor = math.exp(- self.frame_time * integ_arf_s / g0)  
        exp_factor *=  num_regions *fracexpo
        if exp_factor == 0:
            return np.zeros(self.num_points)
            
        # Normalize to avoid overflow and perform FFT convolutions
        arf_s_tmp /= integ_arf_s        
        integ_arf_s_n = integ_arf_s        ## for renormalization after
        arf_s_fft = fft(arf_s_tmp)
        factor = 1

        # Compute FFT with offset
        tmpar = [arf_s_tmp[ie + ioff] for ie in range(self.num_points)]
        arf_s_fft_2 = fft(tmpar)

        # Calculate higher order terms
        for i in range(2, self.max_num_terms + 1):

            integ_arf_s_n *= integ_arf_s   #renormalization factor

            # Convolution via FFT
            arf_s_tmp = np.real(ifft(arf_s_fft * arf_s_fft_2 ** (i-1)))
            arf_s_tmp = np.maximum(arf_s_tmp, 0)
            
            # Apply grade migration factor
            factor *= alpha * self.frame_time / i
            results += factor * integ_arf_s_n * arf_s_tmp 

        # Apply final corrections
        results *= exp_factor  

        # Handle non-piled fraction  
        remaining_frac = 1.0 - psf_frac_init
        if remaining_frac > 0:
            results += arf_s * remaining_frac 

        return results 


    def analyze(self, model_spectrum: np.ndarray, alpha: float, psf_frac: float) -> Tuple[np.ndarray, dict]:
        """
        Perform pileup analysis on input spectrum
        
        Parameters:
        -----------
        model_spectrum : np.ndarray
            Input model spectrum
        alpha : float
            Grade migration scale factor :  probability of the piled event to be assigned a "good" grade
        psf_frac : float
            Fraction of PSF in the source extraction region to which pileup will be applied
            
        Returns:
        --------
        piled_spectrum
            Piled up spectrum in counts/sec 
        """
        piled_spectrum = model_spectrum.copy()

        for k in range(len(model_spectrum.shape[1])):

            # Perform pileup calculation for each phase bin
            spectrum_pileup = self.perform_pileup(model_spectrum[:,k], alpha, psf_frac) 
            piled_spectrum[:,k] = spectrum_pileup
        
        # Apply RMF if available
        if self.rmf_data is not None:
            self.rmf_data = np.maximum(self.rmf_data, 0)
            piled_spectrum =  np.dot(self.rmf_data, piled_spectrum) 
        
        return piled_spectrum
