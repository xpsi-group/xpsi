import numpy as np
import matplotlib.pyplot as plt


def readModeSummary(samples_path,
             mode_number=0,
             verbose=True):
    """ Get the mode from the summary file. 
    
    :param str samples_path:
        Path to the run files without the extension.
    
    :param int mode_number:
        Number of the mode to recover. If 0, takes everything.
    
    :param bool verbose:
        Should the information be printed?
    """

    # Load the summary
    summary = np.loadtxt(f'{samples_path}summary.txt')
    Nmodes = np.shape(summary)[0] - 1
    assert mode_number <= Nmodes, 'The mode to recover must be between 0 and %i' % Nmodes
    summaryMode = summary[mode_number,:-2]
    Npar = len(summaryMode) // 4    # There are 4 metrics per mode

    # Get the values
    AverageMode = summaryMode[:Npar]
    SigmaMode = summaryMode[Npar:2*Npar]
    BestFitMode = summaryMode[2*Npar:3*Npar]
    MAP_Mode = summaryMode[3*Npar:4*Npar]

    if verbose:
        print( f"Mode Recovered: {mode_number}" )
        print( "AVERAGE: ",AverageMode,len(AverageMode))
        print ("STD DEVIATION: ",SigmaMode,len(SigmaMode))
        print ("MAX LIKELIHOOD: ",BestFitMode,len(BestFitMode))
        print ("MAX POSTERIOR: ",MAP_Mode,len(MAP_Mode),'\n')

    return AverageMode ,SigmaMode, BestFitMode, MAP_Mode


def readSummary(samples_path,
                verbose=True):
    """ Read the summary file of the whole run. 
    
    :param str samples_path:
        Path to the run files without the extension.
    
    :param bool verbose:
        Should the information be printed?
    """
    return readModeSummary( samples_path=samples_path , mode_number=0 , verbose=verbose )


def getSignal( XPSImodel , InstrumentName):
    """ Get the signal depending on the used instrument.
    
    :param obj XPSImodel:
        An instance of an imported model containing the signal.
    
    :param str InstrumentName:
        Name of the instrument used, usually XTI, NICER, PN, MOS1 or MOS2.
    """

    # Check the signal of the instrument
    if not hasattr( XPSImodel, 'signal' ):

        # Check the number of lists
        signals = XPSImodel.signals[0]
        try:
            _ = signals[0].prefix
        except AttributeError:
            signals = XPSImodel.signals
        InstrumentLoaded = [ signals[i].prefix for i in range(len(signals)) ]
        InstrumentExists = InstrumentName in InstrumentLoaded

        # If requested instrument exists, get its associated index
        # Otherwise, try getting the signal directly
        if InstrumentExists:
            idx = InstrumentLoaded.index( InstrumentName )
            signal = signals[idx]
        else:
            print(f'Required instrument could not be found. Please choose from {InstrumentLoaded}')

    else:
        signal = XPSImodel.signal
        if hasattr( signal, 'prefix' ):
            print(f'Loaded the only existing signal {signal.prefix}. Please check that it is from the required instrument')
        else:
            print(f'Loaded the only existing signal. Please check that it is from the required instrument')


    return signal


def extractBKG(p, 
                XPSImodel , 
                InstrumentName,
                given_support = True ):
    """ Extract the background from given parameter values.
    
    :param list p:
        Parameter values.
    
    :param obj XPSImodel:
        An instance of an imported model containing the signal.
    
    :param str InstrumentName:
        Name of the instrument used, usually XTI, NICER, PN, MOS1 or MOS2.
    
    :param bool given_support:
        Should the background be extracted from the given support or the full support
    """
    
    # Compute the likelihood and get the signal then
    XPSImodel.likelihood(p, reinitialise=True)
    signal = getSignal( XPSImodel=XPSImodel, InstrumentName=InstrumentName )
    if given_support:
        try:
            return signal.background_signal_given_support
        except AttributeError:
            print('WARNING: signal.background_signal_given_support does not exist, falling back to background_signal. '
                  'Please check that it is correct')
            return signal.background_signal
    else:
        return signal.background_signal


def plotBackgroundSpectrum( XPSI_model, 
                            samples_path, 
                            InstrumentName,
                            plot_params=None,
                            Nsamples=200,
                            plot_range=True, 
                            yscale='linear',
                            plot_support=False):
    """ Plot the spectrum generated along data for given parameters.
    
    :param obj XPSI_model:
        An instance of an imported model containing the signal.
    
    :param str samples_path:
        Path to the posterior file without the extension.
    
    :param str InstrumentName:
        Name of the instrument used, usually XTI, NICER, PN, MOS1 or MOS2.
    
    :param list | None plot_params:
        Parameter values to plot, plot Best Fit if None.

    :param int Nsamples:
        Number of samples to plot the background from.
    
    :param bool plot_range:
        Should the background uncertainty ranges be plotted.
    
    :param str yscale:
        Scale of the y-axis.

    :param bool plot_support:
        Should the support of the background be plotted ? 
    """
                    
    # Get signal and data
    signal = getSignal(XPSImodel=XPSI_model, InstrumentName=InstrumentName)
    data = signal.data
    Data_Spectrum = data.counts.sum(axis=1)

    # Get the parameters to plot the spectrum and backgound 
    if plot_params is None:
        _ ,_, p, _ = readSummary(samples_path,verbose=False)
    else:
        p = plot_params
    BKG = extractBKG(p=p, XPSImodel=XPSI_model, InstrumentName=InstrumentName, given_support=True) 

    # Get expected counts from both spots
    num_components = signal.num_components
    HotRegion_spectra = np.array( [ np.sum(signal.signals[i], axis=1)/float(len(signal.phases[0]))*data.exposure_time for i in range(num_components)] )
    Expected_Spectrum = signal.expected_counts.sum(axis = 1) 
    print( f"Maximum counts in an energy bin : {np.max( HotRegion_spectra.sum(axis=0) )}")
    
    # Extract channels
    x0 = signal.instrument.channel_edges
    x0 = ( x0[:-1] + x0[1:] ) / 2
    
    # Extract background from samples
    if plot_range:
        posterior = np.loadtxt(samples_path+'post_equal_weights.dat')
        if posterior.ndim == 1:
            raise IndexError("*post_equal_weights.dat has only one sample and thus you cannot plot background uncertainty range. Please set plot_range=False.")
        indexR = np.random.randint(low=0, high=len(posterior)-1, size=Nsamples)
        BKG_A = np.array( [extractBKG( p=posterior[indexR[i]][:-1], XPSImodel=XPSI_model, InstrumentName=InstrumentName, given_support=True) for i in range(Nsamples)] )
        sigma = np.std(BKG_A,axis = 0)
        mean = np.mean(BKG_A,axis = 0)
        print( f"Max deviation : {np.max(sigma)} counts")

    # Do the plotting
    fig, ax = plt.subplots(1,1,figsize=(13,8))
    cm = plt.get_cmap( ['RdPu','BuPu','YlGnBu'][0] ) 
    mycolors = [cm(xcol) for xcol in np.linspace(0,1 , 8)]
    ax.grid(linestyle="-.", color='grey',linewidth=0.7)
    ax.set_ylabel(r"Spectrum [counts]",fontsize=15 )
    ax.set_xlabel("Energy Channel [keV]",fontsize=15 )
    
    # Prepare labels
    lBs = r'$\mathrm{BKG\, MEAN\pm 1\sigma}$'
    lBp = 'BKG@MaxL'
    
    # Plotting outputs
    number_labels = ['Primary@MaxL','Secondary@MaxL','Tertiary@MaxL']
    for i in range(num_components):
        ax.fill_between(x0, HotRegion_spectra[:i].sum(axis=0), HotRegion_spectra[:i+1].sum(axis=0), color=mycolors[i+1], alpha = 0.5, label=number_labels[i] )
    ax.fill_between(x0, HotRegion_spectra.sum(axis=0), HotRegion_spectra.sum(axis=0)+BKG, color=mycolors[num_components+1], alpha = 0.5, label ='BKG@MaxL')
    
    # Plotting actual background
    if plot_range:
        ax.fill_between(x0, np.abs(mean-3*sigma), (mean+3*sigma), color =mycolors[7],alpha = 0.3,label =lBs.replace(r'1\sigma',r'3\sigma'))
        ax.fill_between(x0, np.abs(mean-1*sigma), (mean+1*sigma), color =mycolors[7],alpha = 0.5,label =lBs)
        ax.fill_between(x0, np.abs(mean-2*sigma), (mean+2*sigma), color =mycolors[7],alpha = 0.4,label =lBs.replace(r'1\sigma',r'2\sigma'))
    ax.plot(x0,BKG,color =mycolors[6],label = lBp,lw = 2)

    # Plotting data and expected values
    ax.plot(x0,Expected_Spectrum, color=mycolors[4], lw=3, label='Expected signal')
    ax.plot(x0,Data_Spectrum,'--', color=mycolors[5], lw=2, label='Data light curve')
    
    # Plotting background support
    if plot_support and signal._support is not None and signal._support[signal._support>0].any():
        support = signal._support * data.exposure_time
        ax.fill_between(x0, support[:,0], support[:,1], color='red', alpha = 0.2, label='BKG prior support')

    # Finish the plot
    _ = ax.legend(fontsize=15, loc='upper right' )
    ax.set_yscale( yscale )
    
    return fig, ax
