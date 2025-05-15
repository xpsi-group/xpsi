import numpy as np
import pytest
import math
import os
from collections import OrderedDict

import xpsi
from xpsi.global_imports import gravradius
from xpsi.utilities.ProjectionTool import plot_projection_general
from xpsi.utilities.BackgroundTools import readSummary, plotBackgroundSpectrum

class TestPostProcessing(object):

    def setup_class(cls):

        import main as ST
        ST.names=['mass','radius','distance','cos_inclination','p__phase_shift',
                  'p__super_colatitude','p__super_radius','p__super_temperature']
        ST.bounds = {'mass':(1.0,1.6),
                     'radius':(10,13),
                     'distance':(0.5,2.0),
                     'cos_inclination':(0,1),
                     'p__phase_shift':(-0.25, 0.75),
                     'p__super_colatitude':(0.001, math.pi/2 - 0.001),
                     'p__super_radius':(0.001, math.pi/2.0 - 0.001),
                     'p__super_temperature':(6., 7.)}
        ST.labels = {'mass': r"M\;\mathrm{[M}_{\odot}\mathrm{]}",
                      'radius': r"R_{\mathrm{eq}}\;\mathrm{[km]}",
                      'distance': r"D \;\mathrm{[kpc]}",
                      'cos_inclination': r"\cos(i)",
                      'p__phase_shift': r"\phi_{p}\;\mathrm{[cycles]}",
                      'p__super_colatitude': r"\Theta_{\mathrm{spot}}\;\mathrm{[rad]}",
                      'p__super_radius': r"\zeta_{\mathrm{spot}}\;\mathrm{[rad]}",
                      'p__super_temperature': r"\mathrm{log10}(T_{\mathrm{spot}}\;[\mathrm{K}])"}
        ST.names +=['compactness']
        ST.bounds['compactness']=(gravradius(1.0)/16.0, 1.0/3.0)
        ST.labels['compactness']= r"M/R_{\mathrm{eq}}"

        cls.ST = ST

        cls.getdist_kde_settings = {'ignore_rows': 0,
                                 'min_weight_ratio': 1.0e-10,
                                 'contours': [0.683, 0.954, 0.997],
                                 'credible_interval_threshold': 0.001,
                                 'range_ND_contour': 0,
                                 'range_confidence': 0.001,
                                 'fine_bins': 1024,
                                 'smooth_scale_1D': 0.4,
                                 'num_bins': 100,
                                 'boundary_correction_order': 1,
                                 'mult_bias_correction_order': 1,
                                 'smooth_scale_2D': 0.4,
                                 'max_corr_2D': 0.99,
                                 'fine_bins_2D': 512,
                                 'num_bins_2D': 40}

    def test_load_run_works(self):
        'Testing that loading runs does not fail when it should not.'
        ST = self.ST
        getdist_kde_settings = self.getdist_kde_settings
        ST.runs = xpsi.Runs.load_runs(ID='ST',
                                       run_IDs=['run'],
                                       roots=['ST_live_1000_eff_0.3_seed0'],
                                       base_dirs=[os.path.join(os.path.dirname(os.path.abspath(__file__)),'../../examples/examples_fast/Outputs/')],
                                       use_nestcheck=[False],
                                       kde_settings=getdist_kde_settings,
                                       likelihood=ST.likelihood,
                                       names=ST.names,
                                       bounds=ST.bounds,
                                       labels=ST.labels,
                                       implementation='multinest',
                                       overwrite_transformed=True)

    def test_load_run_fails(self):
        'Testing that loading runs fails if lengths of names and bounds do not match.'
        ST = self.ST
        getdist_kde_settings = self.getdist_kde_settings
        with pytest.raises(TypeError):
            _ = xpsi.Runs.load_runs(ID='ST',
                                           run_IDs=['run'],
                                           roots=['ST_live_1000_eff_0.3_seed0'],
                                           base_dirs=[os.path.join(os.path.dirname(os.path.abspath(__file__)),'../../examples/examples_fast/Outputs/')],
                                           use_nestcheck=[False],
                                           kde_settings=getdist_kde_settings,
                                           likelihood=ST.likelihood,
                                           names=ST.names+['extra parameter'],
                                           bounds=ST.bounds,
                                           labels=ST.labels,
                                           implementation='multinest',
                                           overwrite_transformed=True)

    def test_cornerplotter_works(self):
        'Testing that CornerPlotter does not fail when it should not.'
        'Checking also that a credible interval length and a median value are what they should be.'
        ST = self.ST

        pp = xpsi.PostProcessing.CornerPlotter([ST.runs])
        _ = pp.plot(
             params=ST.names,
             IDs=OrderedDict([('ST', ['run',]),]),
             prior_density=True,
             KL_divergence=True,
             ndraws=1e2,
             bootstrap_estimators=False,
             compute_all_intervals=False,
             sixtyeight=True)        

        credible_intervals=pp.credible_intervals
        mass_median = credible_intervals["ST_run"][0][0]
        mass_68interval = credible_intervals["ST_run"][0][2]-credible_intervals["ST_run"][0][1]

        assert np.isclose(1.434, mass_median, rtol=1.0e-3) 
        assert np.isclose(0.09, mass_68interval, rtol=1.0e-3)

    def test_residual_plotting_works(self):
        'Testing that the residual plotter does not fail when it should not.'

        pp2 = xpsi.SignalPlotter([self.ST.runs])
        plots = {'ST': xpsi.Residual1DPlot(
                     parameters_vector= None,
                     nbins=50,
                     plot_fit=True)}
        pp2.plot(IDs=OrderedDict([('ST', ['run']),]),
                combine=False,
                combine_all=False,
                force_combine=False,
                only_combined=False,
                force_cache=True,
                nsamples=2,
                plots = plots)

    def test_pulse_plotting_works(self):
        'Testing that the pulse plotter does not fail when it should not.'

        pp3 = xpsi.SignalPlotter([self.ST.runs])
        pp3.plot(IDs=OrderedDict([('ST', ['run',]),]),  
                nsamples=2,
                plots = {'ST': xpsi.PulsePlot(use_fgivenx=True,
                         num_phases=5,
                         root_filename= 'ST_model_signal',),})

    def test_spectrum_plotting_fails(self):
        'Testing that the spectrum plotter fails if providing negative number of phases.'

        pp4 = xpsi.SignalPlotter([self.ST.runs])
        with pytest.raises(ValueError):
            pp4.plot(IDs=OrderedDict([('ST', ['run',]),]),     
                    nsamples=2,
                    plots = {'ST': xpsi.SpectrumPlot(use_fgivenx=True,
                             rel_num_energies=1.0,
                             num_phases=-1,
                             show_attenuated=False,
                             root_filename= 'ST_spectra',),})

    def test_background_plotting_works(self):
        'Testing that the background plotter does not fail when it should not.'

        samples_path = os.path.join(os.path.dirname(os.path.abspath(__file__)),'../../examples/examples_fast/Outputs/ST_live_1000_eff_0.3_seed0_v2')
        fig, ax = plotBackgroundSpectrum(XPSI_model=self.ST, 
                                samples_path=samples_path, 
                                InstrumentName=None,
                                Nsamples=2,
                                plot_range=False)

    def test_background_plotting_fails(self):
        'Testing that the background plotter fails when trying to plot background uncertainty range'
        'with MultiNest output file having only 1 equally weighted posterior sample.'

        with pytest.raises(IndexError):
            samples_path = os.path.join(os.path.dirname(os.path.abspath(__file__)),'../../examples/examples_fast/Outputs/ST_live_1000_eff_0.3_seed0_v2')
            fig, ax = plotBackgroundSpectrum(XPSI_model=self.ST, 
                                    samples_path=samples_path, 
                                    InstrumentName=None,
                                    Nsamples=2,
                                    plot_range=True)
                                    
    def test_projection_tool_works(self):
        'Testing that the projection tool does not fail when it should not.'

        samples_path = os.path.join(os.path.dirname(os.path.abspath(__file__)),'../../examples/examples_fast/Outputs/ST_live_1000_eff_0.3_seed0_v2')
        AverageP ,SigmaP, BestFitP, MAP_P = readSummary(samples_path,verbose=False)
        self.ST.likelihood(BestFitP,reinitialise=True)
        ax = plot_projection_general(self.ST.likelihood, 'ST', "I","SP", antiphase=False, SaveFlag = False, Name = 'BestFit' )

