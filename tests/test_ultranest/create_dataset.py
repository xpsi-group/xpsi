""""Create x amount of synthetic datasets with random varying parameter values. The datasets 
all have between 10^5 and 10^7 counts, so that they aren't too noisy or take too long to sample.
Currently, no background is assumed. This is done by setting the expected background counts to 0.
"""

import math
import numpy as np
import os
import random 
import subprocess

# create file storing synthetic data parameters 
filepath = "../synthetic_data/syndat_overview.txt"
header = "# sample_number, mass, radius, distance, cos_inclination, phase_shift, super_colatitude, super_radius, super_temperature, background, exposure_time, expected_background_counts"

os.makedirs(os.path.dirname(filepath), exist_ok=True)
with open(filepath, "w+") as f:
    f.write(header)
    f.close()

try: 
    max_sample_size = 10 # amount of synthetic datasets you want to create 
    sample_number = 1 
    while sample_number <= max_sample_size: 

        # pick random parameter values between specified bounds 
        mass = random.uniform(1., 3.)                               # Mass in solar Mass
        radius = random.uniform(10., 15.)                           # Equatorial radius in km
        distance = random.uniform(0.5, 2.)                          # Distance in kpc
        cos_inclination = random.uniform(0., 1.)                    # Cosine of Earth inclination to rotation axis
        phase_shift = random.uniform(-0.25, 0.75)                   # Phase shift
        super_colatitude = random.uniform(0.001, math.pi/2 - 0.001) # Colatitude of the centre of the superseding region
        super_radius = random.uniform(0.001, math.pi/2 - 0.001)     # Angular radius of the (circular) superseding region
        super_temperature = random.uniform(6.5, 7.2)                # Temperature in log 10
        background = random.uniform(1., 3.)                         # Background spectral index : gamma (E^gamma) 
        exposure_time = random.uniform(100, 1000.0)                 # Exposure time in seconds
        expected_background_counts = 0.0                            # No background counts expected 
        
        directory = f"../synthetic_data/dataset_{sample_number}/"
        os.makedirs(os.path.dirname(directory), exist_ok=True)
        name = f"syndat_{sample_number}"
        
        # create synthetic data 
        print(f"Working on dataset {sample_number}")
        completed_process = subprocess.run(["python", "create_synthetic_data.py",
                                            "-d", directory, "-n", name, "-t", str(exposure_time), "-b", str(expected_background_counts),
                                            "-p", str(mass), str(radius), str(distance), str(cos_inclination), str(phase_shift), 
                                            str(super_colatitude), str(super_radius), str(super_temperature), str(background)
                                            ])
        try: 
            # check synthetic data
            new_synthetic_data = np.loadtxt(f"{directory}{name}_realisation.dat")
            total_photon_count = np.sum(new_synthetic_data)

            # discard data that is too noisy (below 10^5 counts)
            # or for which sampling takes too long (above 10^7 counts)
            if 10**5 < total_photon_count < 10**7:
                with open(filepath, "a") as f:
                    f.write(f"""\n{sample_number}\t {mass}\t {radius}\t {distance}\t {cos_inclination}\t {phase_shift}\t {super_colatitude}\t {super_radius}\t {super_temperature}\t {background}\t {exposure_time}\t {expected_background_counts}""")
                    f.close()

                sample_number += 1
        
        except FileNotFoundError:
            "No new synthetic data created within photon count bounds (between 10^5 and 10^7). Trying again."

except KeyboardInterrupt:
    pass

