import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate
from matplotlib.colors import LinearSegmentedColormap
from cmap import Colormap


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
                            clrmap='petroff10',
                            custcols=None, 
                            yscale='linear',
                            plot_support=False,
                            plot_data_errors=False,
                            plot_residuals=False):
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
    
    :param str clrmap:
        Determine the colormap for the plot. Can be customised by setting 'CustomCmap'.
        Otherwise, the available colormaps are 'Accent', 'Dark2', 'tab10_r', the default
        colorblind friendly 'petroff10' or 'petroff10_r'.
    
    :par list | None custcols:
        A list of colors for a customised colormap if required.
        Must be large enough for all elements to plot.
    
    :param str yscale:
        Scale of the y-axis.

    :param bool plot_support:
        Should the support of the background be plotted ?  

    :param bool plot_data_errors:
        Should error bars be plotted on the data spectrum? The errors are
        taken from ``data.errors`` if available (e.g. PHA data loaded with
        ``withuncs=True``), and otherwise assumed to be Poissonian,
        i.e. :math:`\\sqrt{N}` (with :math:`N=0` channels assigned an
        error of 1 count).

    :param bool plot_residuals:
        Should a lower panel with the residuals
        :math:`(\\mathrm{data}-\\mathrm{model})/\\sigma` be added? The model is
        the expected signal plus the background at the plotted parameters
        (``BKG@MaxL``). The uncertainties :math:`\\sigma` are the same as for
        ``plot_data_errors``.

    :returns:
        ``fig, ax`` or, if ``plot_residuals`` is ``True``,
        ``fig, (ax, ax_res)`` where ``ax_res`` is the residuals panel.
    """
                    
    # Get signal and data
    signal = getSignal(XPSImodel=XPSI_model, InstrumentName=InstrumentName)
    data = signal.data

    # Restrict to the channels actually used (grouping in Signal may trim
    # data.channels without trimming data.counts)
    a, b = data.index_range
    Data_Spectrum = data.counts[a:b].sum(axis=1)

    # Uncertainties on the data spectrum, used for error bars and residuals
    if plot_data_errors or plot_residuals:
        if getattr(data, 'errors', None) is not None:
            Data_Errors = np.sqrt( (data.errors[a:b]**2).sum(axis=1) )
        else:
            Data_Errors = np.sqrt( np.where(Data_Spectrum > 0, Data_Spectrum, 1.0) )

    # Get the parameters to plot the spectrum and backgound 
    if plot_params is None:
        _ ,_, p, _ = readSummary(samples_path,verbose=False)
    else:
        p = plot_params
    BKG = extractBKG(p=p, XPSImodel=XPSI_model, InstrumentName=InstrumentName, given_support=True) 

    # Get expected counts from both spots
    num_components = signal.num_components
    HotRegion_spectra = np.array( [ integrate.simpson( y=signal.signals[i], x=signal.phases[i] , axis=1 ) * data.exposure_time for i in range(num_components)] )
    Expected_Spectrum = signal.expected_counts.sum(axis = 1) 
    print( f"Maximum counts in an energy bin : {np.max( HotRegion_spectra.sum(axis=0) )}")
    
    # Extract channels
    x0 = signal.instrument.channel_edges
    if x0.ndim == 2:
        x0 = ( x0[0] + x0[1] ) / 2
    else:
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
    
    # Prepare labels
    number_labels = ['Primary@MaxL','Secondary@MaxL','Tertiary@MaxL']
    lBs = r'$\mathrm{BKG\, MEAN\pm 1\sigma}$'
    lBp = 'BKG@MaxL'
    
    # Do the plotting
    if plot_residuals:
        fig, (ax, ax_res) = plt.subplots(2, 1, figsize=(13,11), sharex=True,
                                         gridspec_kw={'height_ratios': [4,1], 'hspace': 0.03})
    else:
        fig, ax = plt.subplots(1,1,figsize=(13,8))

    if clrmap == 'CustomCmap':
        try:
            assert type(custcols) == list
            assert len(custcols) >= len(number_labels)+5
            clrmap = LinearSegmentedColormap.from_list("CustomCmap", custcols)
            cm=plt.get_cmap(clrmap)
            cm.colors = custcols
        except AssertionError:
            raise AssertionError("For a customised colormap, we need a list of colors large enough for all potential elements to plot: the signal components, the background, the data and the expected signal.")
    elif clrmap in ['petroff10', 'petroff10_r']:
        CMtoplot=Colormap('petroff:'+clrmap)
        cm=CMtoplot.to_matplotlib()
    elif clrmap in ['Accent', 'Dark2', 'tab10_r']:
        cm=plt.get_cmap(clrmap)
    else:
        raise TypeError("If it exists, this colormap doesn't help to distinguish the different components clearly. Please select one of the listed colormaps for the parameter 'clrmap' ('Accent', 'Dark2', 'tab10_r', the default colorblind friendly 'petroff10' or 'petroff10_r') or build your own, setting clrmap='CustomCmap'.")         
    mycolors = [cm(xcol) for xcol in np.linspace(0,1 , len(cm.colors))]
    ax.grid(linestyle="-.", color='grey',linewidth=0.7)
    ax.set_ylabel(r"Spectrum [counts]",fontsize=15 )

    if not plot_residuals:
        ax.set_xlabel("Energy Channel [keV]",fontsize=15 )
    
    # Plotting outputs
    for i in range(num_components):
        ax.fill_between(x0, HotRegion_spectra[:i].sum(axis=0), HotRegion_spectra[:i+1].sum(axis=0), color=mycolors[i], alpha = 0.5, label=number_labels[i] )
    ax.fill_between(x0, HotRegion_spectra.sum(axis=0), HotRegion_spectra.sum(axis=0)+BKG, color=mycolors[num_components], alpha = 0.5, label ='BKG@MaxL')
    
    # Plotting actual background
    if plot_range:
        ax.fill_between(x0, np.abs(mean-3*sigma), (mean+3*sigma), color =mycolors[7],alpha = 0.3,label =lBs.replace(r'1\sigma',r'3\sigma'))
        ax.fill_between(x0, np.abs(mean-1*sigma), (mean+1*sigma), color =mycolors[7],alpha = 0.5,label =lBs)
        ax.fill_between(x0, np.abs(mean-2*sigma), (mean+2*sigma), color =mycolors[7],alpha = 0.4,label =lBs.replace(r'1\sigma',r'2\sigma'))
    ax.plot(x0,BKG,color =mycolors[6],label = lBp,lw = 1)
    

    # Plotting data and expected values
    ax.plot(x0,Expected_Spectrum, color=mycolors[4], lw=2, label='Expected signal')  # lw=3
    if plot_data_errors:
        ax.errorbar(x0, Data_Spectrum, yerr=Data_Errors, fmt='--', color=mycolors[5], lw=1, elinewidth=1, label='Data spectrum')
    else:
        ax.plot(x0,Data_Spectrum,'--', color=mycolors[5], lw=1, label='Data spectrum') # lw=2
    
    # Plotting background support
    if plot_support and signal._support is not None and signal._support[signal._support>0].any():
        support = signal._support * data.exposure_time
        ax.fill_between(x0, support[:,0], support[:,1], color='red', alpha = 0.2, label='BKG prior support')

    # Finish the plot
    _ = ax.legend(fontsize=15, loc='upper right' )
    ax.set_yscale( yscale )

    # Residuals panel
    if plot_residuals:
        # Full model: expected source counts plus the conditional max-likelihood background
        Residuals = (Data_Spectrum - Expected_Spectrum) / Data_Errors
        ax_res.grid(linestyle="-.", color='grey',linewidth=0.7)
        ax_res.axhline(0.0, color='k', lw=1)
        ax_res.errorbar(x0, Residuals, yerr=np.ones_like(Residuals), fmt='.', color=mycolors[5], elinewidth=1, label='Residuals')
        ylim = max(4.0, 1.1 * (np.nanmax(np.abs(Residuals)) + 1.0))
        ax_res.set_ylim(-ylim, ylim)
        ax_res.set_ylabel(r"Residuals [$\sigma$]",fontsize=15 )
        ax_res.set_xlabel("Energy Channel [keV]",fontsize=15 )
        return fig, (ax, ax_res)

    return fig, ax
