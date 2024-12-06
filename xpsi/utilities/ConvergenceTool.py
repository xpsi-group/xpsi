import numpy as np
import matplotlib.pyplot as plt

def plotConvergence( samples_path , threshold=0. ):
    """ Plot to asses the convergence of the run.
    
    You should see a mountain and then a rise in the right    
    :param str samples_path: 
        path to the samples file without the extension.

    :param float threshold: 
        Threshold on the probability to keep in the plotting.
    """

    # Extract samples
    samples = np.loadtxt(samples_path+'.txt')

    # Order and extract
    orderLogL = np.argsort(-0.5*samples[:,1])
    LogL_ordered  = -0.5*samples[ orderLogL,1]
    sample_proba_ordered = samples[ orderLogL,0]

    if threshold != 0.:
        LogL_ordered = LogL_ordered[sample_proba_ordered>threshold]
        sample_proba_ordered = sample_proba_ordered[sample_proba_ordered>threshold]

    # Plot nicely
    fig = plt.figure(figsize=(10,4))
    ax = fig.add_subplot(111)
    ax.grid(linestyle="-.", color='grey',linewidth=0.7)
    ax.plot( LogL_ordered, sample_proba_ordered )
    for tick in ax.xaxis.get_major_ticks():
        tick.label1.set_fontsize(15)
    for tick in ax.yaxis.get_major_ticks():
        tick.label1.set_fontsize(15)
    ax.set_xlabel('LogL')
    ax.set_ylabel(r'$p(\theta | d )$')
    ax.set_title(f'Ordered sample probability : {samples_path}')
    
    return fig, ax
