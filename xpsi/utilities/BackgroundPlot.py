import numpy as np
import matplotlib.pyplot as plt

# Get the modes
def get_mode(run_files,
             mode_number=0,
             verbose=True):

    # Load the summary
    summary = np.loadtxt(f'{run_files}summary.txt')
    Nmodes = np.shape(summary)[0] - 1
    assert mode_number <= Nmodes, 'The mode to recover must be between 0 and %i' % Nmodes
    summaryMode = summary[mode_number,:-2]
    Npar = len(summaryMode) // 4    # There are 4 metrics per mode

    # Assert that the mode exists
    assert mode_number <= Nmodes
    AverageMode = summaryMode[:Npar]
    SigmaMode = summaryMode[Npar:2*Npar]
    BestFitMode = summaryMode[2*Npar:3*Npar]
    MAP_Mode = summaryMode[3*Npar:4*Npar]

    if verbose:
        print( "AVERAGE: ",AverageMode,len(AverageMode))
        print ("STD DEVIATION: ",SigmaMode,len(SigmaMode))
        print ("MAX LIKELIHOOD: ",BestFitMode,len(BestFitMode))
        print ("MAX POSTERIOR: ",MAP_Mode,len(MAP_Mode),'\n')

    return AverageMode ,SigmaMode, BestFitMode, MAP_Mode


# Get the signal depending on the used instrument
def get_signal( XPSImodel , InstrumentName):

    # Check the instrument
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
        try:
            signal = XPSImodel.signal
            print(f'Loaded the only existing signal {signal.prefix}. Please check that it is from the required instrument')
        except AttributeError:
            raise AttributeError(f'Required instrument {InstrumentName} not existing and no default either')

    return signal


# Extract the background from given parameter values 
def extract_BKG(p, 
                XPSImodel , 
                InstrumentName,
                given_support = True ):
    
    # Compute the likelihood and get the signal then
    XPSImodel.likelihood(p, reinitialise=True)
    signal = get_signal( XPSImodel=XPSImodel, InstrumentName=InstrumentName )
    if given_support:
        return signal.background_signal_given_support
    else:
        return signal.background_signal

# Plot the spectrum generated along data for given parameters
def plot_spectrum_GEN( XPSI_model, 
                        posterior_file, 
                        InstrumentName,
                        Nsamples=200,
                        AmpF=1, 
                        plot_bkg=True, 
                        yscale='linear'):
                    
    # Get signal and data
    signal = get_signal(XPSImodel=XPSI_model, InstrumentName=InstrumentName)
    data = signal.data
    Data_Spectrum = data.counts.sum(axis=1)

    # Get the inferred backgound 
    _ ,_, BestFitPSpectrum, _ = get_mode(posterior_file,mode_number=0,verbose=False)
    BKG = extract_BKG(p=BestFitPSpectrum, XPSImodel=XPSI_model, InstrumentName=InstrumentName) 

    # Get expected counts from both spots
    num_components = signal.num_components
    HotRegion_spectra = np.array( [ np.sum(signal.signals[i], axis=1)/float(len(signal.phases[0]))*data.exposure_time for i in range(num_components)] )
    Expected_Spectrum = signal.expected_counts.sum(axis = 1) 
    print( f"Maximum counts in an energy bin : {np.max( HotRegion_spectra.sum(axis=0) )}")
    
    # Extract channels
    x0 = getattr(XPSI_model,InstrumentName).instrument.channel_edges
    x0 = ( x0[:-1] + x0[1:] ) / 2
    
    # Extract background from samples
    if plot_bkg:
        posterior = np.loadtxt(posterior_file+'post_equal_weights.dat')
        indexR = np.random.randint(low=0, high=len(posterior)-1, size=Nsamples)
        BKG_A = np.array( [extract_BKG( p=posterior[indexR[i]][:-1], XPSImodel=XPSI_model, InstrumentName=InstrumentName) for i in range(Nsamples)] )
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
    lBs = r'$\mathrm{BKG\, MEAN\pm 1\sigma}$' if AmpF ==1 else r'$\mathrm{BKG\, MEAN\pm 1\sigma \times %i}$'%AmpF
    lBp = 'BKG@MaxL'  if AmpF ==1 else r'BKG $\times %i$ @MaxL'%AmpF
    
    # Plotting outputs
    number_labels = ['Primary@MaxL','Secondary@MaxL','Tertiary@MaxL']
    for i in range(num_components):
        ax.fill_between(x0, HotRegion_spectra[:i].sum(axis=0), HotRegion_spectra[:i+1].sum(axis=0), color=mycolors[i+1], alpha = 0.5, label=number_labels[i] )
    ax.fill_between(x0, HotRegion_spectra.sum(axis=0), HotRegion_spectra.sum(axis=0)+BKG, color=mycolors[num_components+1], alpha = 0.5, label ='BKG@MaxL')
    
    # Plotting actual background
    if plot_bkg:
        ax.fill_between(x0, np.abs(mean-1*sigma)*AmpF, (mean+1*sigma)*AmpF, color =mycolors[7],alpha = 0.5,label =lBs)
        ax.fill_between(x0, np.abs(mean-2*sigma)*AmpF, (mean+2*sigma)*AmpF, color =mycolors[7],alpha = 0.4,label =lBs.replace('1\sigma','2\sigma'))
        ax.fill_between(x0, np.abs(mean-3*sigma)*AmpF, (mean+3*sigma)*AmpF, color =mycolors[7],alpha = 0.3,label =lBs.replace('1\sigma','3\sigma'))
    ax.plot(x0,BKG*AmpF,color =mycolors[6],label = lBp,lw = 2)

    # Plotting data and expected values
    ax.plot(x0,Expected_Spectrum, color=mycolors[4], lw=3, label='Expected signal')
    ax.plot(x0,Data_Spectrum,'--', color=mycolors[5], lw=2, label='Data light curve')
    
    # Plotting background support
    if signal.support is not None and signal.support[signal.support>0].any():
        support = signal.support * data.exposure_time
        ax.fill_between(x0, support[:,0]*AmpF, support[:,1]*AmpF, color='red', alpha = 0.2, label='BKG prior support')

    # Finish the plot
    _ = ax.legend(fontsize=15, loc='upper right' )
    ax.set_yscale( yscale )
    
    return fig, ax