import numpy as np
import matplotlib.pyplot as plt
import os 

def plotConvergence( samples_path , threshold=0. ):
    """ Plot to asses the convergence of the run.
    
    You should see a mountain and then a rise in the right    
    :param str samples_path: 
        path to the samples file without the extension.

    :param float threshold: 
        Threshold on the probability to keep in the plotting.
    """

    # Extract samples from multimode file if multimode
    multimode = os.path.exists(samples_path+'post_separate.dat')
    if multimode:

        # Open
        with open(samples_path+"post_separate.dat", 'r') as file:
            lines = file.readlines()

        # Convert lines to a numpy array, replacing blank lines with a marker (e.g., None or NaN)
        data = []
        for line in lines:
            stripped = line.strip()
            if stripped:  # Non-blank line
                data.append(np.array(stripped.split(), dtype=float))
            else:  # Blank line
                data.append(None)

        # Separate modes based on None markers
        samples = []
        current_mode = []
        for row in data:
            if row is None:
                if current_mode:
                    samples.append(np.array(current_mode))
                    current_mode = []
            else:
                current_mode.append(row)

        # Append the last mode if not empty
        if current_mode:
            samples.append(np.array(current_mode))

    else:
        samples = [np.loadtxt(samples_path+'.txt')]

    # Order and extract
    Nmodes = len( samples )
    orderLogL = [ np.argsort(-0.5*s[:,1]) for s in samples ]
    LogL_ordered  = [ -0.5*samples[i][ orderLogL[i],1] for i in range(Nmodes) ]
    sample_proba_ordered = [ samples[i][ orderLogL[i],0] for i in range(Nmodes) ]

    if threshold != 0.:
        for i in range(len(samples)):
            LogL_ordered[i] = LogL_ordered[i][sample_proba_ordered[i]>threshold]
            sample_proba_ordered[i] = sample_proba_ordered[i][sample_proba_ordered[i]>threshold]

    # Plot nicely
    fig = plt.figure(figsize=(10,4))
    ax = fig.add_subplot(111)
    ax.grid(linestyle="-.", color='grey',linewidth=0.7)
    for tick in ax.xaxis.get_major_ticks():
        tick.label1.set_fontsize(15)
    for tick in ax.yaxis.get_major_ticks():
        tick.label1.set_fontsize(15)
    ax.set_xlabel('LogL')
    ax.set_ylabel(r'$p(\theta | d )$')
    ax.set_title(f'Ordered sample probability')
    
    # Make the plot
    for i in range(Nmodes):
        ax.plot( LogL_ordered[i], sample_proba_ordered[i], label=f'Mode {i+1}' if multimode else '')
    if multimode:
        ax.legend()

    return fig, ax
