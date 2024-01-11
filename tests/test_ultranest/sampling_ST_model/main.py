""""Adjusted from examples_fast using ultranest as sampler and parameter bounds matching create_dataset.py.
"""
import argparse
import numpy as np
import re 
import os 
import shutil

import xpsi
import time

from CustomInstrument import CustomInstrument
from CustomSignal import CustomSignal
from CustomPhotosphere import CustomPhotosphere
from CustomPrior import CustomPrior

def load_data(directory, name, sample_number):
    """Load the synthetic data and the overview file containing the parameters used to create 
    the synthetic data. Also, retrieve the sample number. 
    """
    
    # get information from syndat_overview.txt 
    overview_file_path = os.path.join(directory, "syndat_overview.txt")
    overview_file = np.loadtxt(overview_file_path)

    try: 

        if sample_number is not None:
            data_loaded = np.loadtxt(f"../../synthetic_data/dataset_{sample_number}/syndat_{sample_number}_realisation.dat", dtype=np.double)
        else:    
            file_path = os.path.join(directory, name)
            data_loaded = np.loadtxt(file_path, dtype=np.double)
            sample_number = int(re.findall(r'\d+', name)[0])

    except Exception as e:
        print(e, "Specify the sample number in order to retrieve synthetic data parameters.")
    
    return overview_file, sample_number, data_loaded 

def create_star_model(overview_file, sample_number, data_loaded ): 
    """Create a simple star model (ST) with only mass as free parameter. 
    Return the likelihood and prior.
    """

    # get the synthetic data parameters 
    row = sample_number -1 
    params = overview_file[row] 
    
    radius = params[2]
    distance = params[3]
    cos_inclination = params[4]
    phase_shift = params[5]
    super_colatitude = params[6]
    super_radius = params[7]
    super_temperature = params[8]
    background = params[9]
    exposure_time = params[10]
    # expected_background_counts = params[11]

    # # Data
    data = xpsi.Data(data_loaded,
                        channels=np.arange(10,301),
                        phases=np.linspace(0.0, 1.0, 33),
                        first=0,
                        last=290,
                        exposure_time=exposure_time) 

    # # Instrument settings
    channel_number=np.arange(0,1501)    # The channel nnumber
    energy_low=np.arange(0,15.01, 0.01) # Lower bounds of each channel
    energy_high=energy_low+0.01         # Upper bounds of each channel
    channel_edges=np.array([list(channel_number),list(energy_low),list(energy_high)]).T

    # ARF
    arf_energy_low=[0.1]
    arf_energy_high=[0.105]
    arf_val=[1800]

    counter=1
    while arf_energy_low[-1]<=14.995:
        arf_energy_low.append(arf_energy_low[-1]+0.005)
        arf_energy_high.append(arf_energy_high[-1]+0.005)
        arf_val.append(1800)
        counter +=1

    ARF=np.array([list(arf_energy_low),
                list(arf_energy_high),
                list(arf_val)]).T

    # RMF
    RMF=np.diag(np.full(counter,1))

    Instrument = CustomInstrument.from_response_files(ARF =ARF,
                                                RMF = RMF,
                                                channel_edges =channel_edges,
                                                max_input = 301,
                                                min_input = 10,
                                                channel=[10,301])

    # # Signal
    signal = CustomSignal(data = data,
                        instrument = Instrument,
                        interstellar = None,
                        background = None, 
                        cache = True,
                        workspace_intervals = 1000,
                        epsrel = 1.0e-8,
                        epsilon = 1.0e-3,
                        sigmas = 10.0) 

    # # Space-time
    bounds = dict(mass = (1.0,3.0))

    values = dict(distance = distance,
                radius = radius,
                cos_inclination = cos_inclination,
                frequency = 314.0)
    
    spacetime = xpsi.Spacetime(bounds=bounds,
                            values=values)

    # # Hot-spot
    values = dict(super_colatitude = super_colatitude,
                super_radius = super_radius,
                phase_shift = phase_shift,
                super_temperature = super_temperature)

    hot_spot = xpsi.HotRegion(bounds={},
                            values=values,
                            symmetry=True,
                            omit=False,
                            cede=False,
                            concentric=False,
                            sqrt_num_cells=32,
                            min_sqrt_num_cells=10,
                            max_sqrt_num_cells=64,
                            num_leaves=100,
                            num_rays=200,
                            is_secondary=True,
                            prefix='hot')

    # # Photosphere
    photosphere = CustomPhotosphere(hot = hot_spot, elsewhere = None,
                                    values=dict(mode_frequency = spacetime['frequency']))

    # # Star
    star = xpsi.Star(spacetime = spacetime, photospheres = photosphere)

    # # Prior
    prior = CustomPrior()

    # # Likelihood
    likelihood = xpsi.Likelihood(star = star, signals = signal,
                                num_energies = 128,
                                threads = 1,
                                externally_updated = True,
                                prior = prior)
    
    return likelihood, prior 

def run_sampler(likelihood, prior, directory, name, sample_number, use_ultranest):
    """Run the sampler, either using UltraNest (default) or with MultiNest. 
    """
    # track sampling time 
    start = time.time()

    # create output directory and output filenames 
    if sample_number is not None:
        output_dir = f"../../synthetic_data/dataset_{sample_number}/"
        output_filename = f"syndat_{sample_number}_output.txt"
    else:
        output_dir = directory
        output_filename = f"{name}_output.txt"

    # start sampling with ultranest
    if use_ultranest: 
        sampler = xpsi.Sample.ultranested(likelihood=likelihood, 
                                            prior=prior, 
                                            sampler_params={'log_dir' : output_dir},
                                            runtime_params={'show_status':True, 'update_interval_volume_fraction': 0.1, 'min_num_live_points':5000},
                                            use_stepsampler=False, 
                                            # stepsampler_params={'max_nsteps' : 400}, 
                                            out_filename=output_filename)
        sampler.print_results()
    
    # otherwise sample multinest 
    else: 
        # create dir to temporarily store multinest output 
        os.makedirs(f"{output_dir}multinest/", exist_ok=True)

        # multinest sampler 
        runtime_params = {'resume': False,
                        'importance_nested_sampling': False,
                        'multimodal': False,
                        'n_clustering_params': None,
                        'outputfiles_basename': f"{output_dir}multinest/syndat_{sample_number}_multinest_output",
                        'n_iter_before_update': 50, #100,
                        'n_live_points': 50,#1000,
                        'sampling_efficiency': 0.3, #0.1,
                        'const_efficiency_mode': False,
                        'evidence_tolerance': 0.1,
                        'max_iter': 100, #-1,
                        'seed' : 0, # Fixing the seed
                        'verbose': True}

        xpsi.Sample.nested(likelihood, prior, **runtime_params)

        # keep only the file needed for the PP-plot and remove rest 
        try: 
            os.replace(f"{output_dir}multinest/syndat_{sample_number}_multinest_output.txt", f"{output_dir}syndat_{sample_number}_multinest_output.txt")
            shutil.rmtree(f"{output_dir}multinest/")
        except Exception as e:
            print(e)


    print('Sampling took', (time.time()-start)/60, 'minutes')


if __name__ == '__main__':

    
    parser = argparse.ArgumentParser(description="Sampling synthetic data")

    parser.add_argument("-d", "--directory", action="store", required=False, dest="directory",
                        default="../../synthetic_data/", help="Specify directory with synthetic data input\
                        file (default: ../synthetic_data/)")

    parser.add_argument("-n", "--name", action="store", required=False, dest="name",
                        default="new_synthetic_data_realisation.dat", help="Specify filename for synthetic data input\
                        file (default: new_synthetic_data)")
    
    parser.add_argument("-s", "--sampnumber", action="store", required=False, dest="sample_number",
                        default=None, type=int, help="Specify sample number of the synthetic data input\
                        file (default: 1)")
    
    parser.add_argument("-u", "--ultranest", action=argparse.BooleanOptionalAction, required=False, dest="use_ultranest",
                        default=True, type=bool, help="Specify which sampler to use: uses UltraNest by default, \
                        --no-ultranest to use MultiNest (default: True)")
    
    clargs = parser.parse_args()
    directory = clargs.directory 
    name = clargs.name
    sample_number = clargs.sample_number
    use_ultranest = clargs.use_ultranest

    # load synthetic data and create ST model
    overview_file, sample_number, data_loaded  = load_data(directory, name, sample_number)
    likelihood, prior = create_star_model(overview_file, sample_number, data_loaded )
    
    # start sampling 
    run_sampler(likelihood, prior, directory, name, sample_number, use_ultranest)