import xpsi
import numpy as np
from matplotlib import pyplot as plt
from matplotlib import rcParams
from matplotlib.ticker import MultipleLocator, AutoLocator, AutoMinorLocator
from matplotlib import gridspec
from matplotlib import cm
from scipy.interpolate import interp1d

def shift_phase(phi,shift):
	return (phi + shift) % 1 

def read_pulse_arcmancer():
    datafile = "../../../polcslab/CompSlab/arcman_res/polar_acc_f600_burst_r12_m1.4_d60_i40_x10_obl.csv"
    input = open(datafile, 'r')
    lines = input.readlines()
    input.close()
    Nchain_size = sum(1 for line in open(datafile))
    c_lines = 1
    egrid = 6
    full_chain= [[] for x in range(egrid+1)]

    input = open(datafile, 'r')
    lines = input.readlines()
    input.close()

    for j in range(0,len(full_chain)):
        for i in range(c_lines,Nchain_size): 
            parts = lines[i].split(",")
            full_chain[j].append(float(parts[j]))
	parts = lines[c_lines].split(",")

    full_chain = np.array(full_chain)

    energy_keV = [4.94]
    phase = full_chain[0,:]
    norm_obsF = np.zeros((len(phase), egrid))
    obsF = np.zeros((len(phase),egrid))
    for i in range(1,egrid+1):
        if(i==1):
	    norm_obsF[:,i-1] = full_chain[i,:]/np.max(full_chain[1,:])
        else:
            norm_obsF[:,i-1] = full_chain[i,:]/full_chain[1,:]
        obsF[:,i-1] = full_chain[i,:]
        
    phase_new = shift_phase(np.array(phase),0.04829)        
    return phase_new, norm_obsF[:,0],norm_obsF[:,1],norm_obsF[:,2]

def plot_pulse_resid(photosphere,signals,phasepol=None,qnpol=None,unpol=None,inpol=None,psind0=0,psind=127):

    """ Plot hot region signals before and after telescope operation. """
 
    phases = signals[0][0].phases[0]
 
    plot_arcmancer = True
    if plot_arcmancer:
        phase_arc, pulse_arc, pulse_arcQ, pulse_arcU = read_pulse_arcmancer()                           
   
    
    phot_sig_cut = photosphere.signal[0][0][psind0:psind,:]
    photQ_sig_cut = photosphere.signalQ[0][0][psind0:psind,:]
    photU_sig_cut = photosphere.signalU[0][0][psind0:psind,:]
             
    #print("Energies used in X-PSI:",signals[0][0].energies)

    plt.subplots_adjust(wspace=0, hspace=0)
    fig = plt.figure(figsize=(7,7))
    ax = fig.add_subplot(211)

    ax.set_ylabel('Signal [arbitrary normalisation]')

    
    temp = np.sum(phot_sig_cut, axis=0) 
    I1p = temp
    I1pn = I1p/np.max(I1p)
    ax.plot(signals[0][0].phases[0], I1pn, 'o-', color='k', lw=0.5, markersize=2)

    if (phasepol is not None and inpol is not None):
    	ax.plot(phasepol,inpol,'--',color='red')
    	
    if plot_arcmancer:
        for ipha in range(0,len(phase_arc)-1):
            if(phase_arc[ipha+1] > phase_arc[ipha]):    
                ax.plot(phase_arc[ipha:ipha+2],pulse_arc[ipha:ipha+2],'--',color='blue',dashes=[2,2])    

    ax2 = fig.add_subplot(212)
    #ax2.plot(phasepol,(inpol-I1pn)/inpol,'--',color='red')
    ax2.plot(phasepol,(inpol-I1pn),'--',color='red')    
    #ax2.set_ylabel('Residuals [%]')
    ax2.set_ylabel('Residuals [arbitrary]')    
    ax2.set_xlabel('Phase [cycles]')
    fig.savefig("figs/signalsI_residX.pdf",bbox_inches='tight')
    
    
    plt.subplots_adjust(wspace=0, hspace=0)
    fig = plt.figure(figsize=(7,7))
    ax = fig.add_subplot(211)
    ax.set_ylabel('Signal [arbitrary normalisation]')  
    Q1p = np.sum(photQ_sig_cut, axis=0) 
    Q1pn = np.copy(Q1p)
    for ip in range(len(Q1pn)):
    	if(I1p[ip] > 1e-10):
    		Q1pn[ip] = Q1p[ip]/I1p[ip]
    	else:
    		Q1pn[ip] = 0.0
    ax.plot(signals[0][0].phases[0], Q1pn, 'o-', color='k', lw=0.5, markersize=2)
    if (phasepol is not None and qnpol is not None):
    	ax.plot(phasepol,qnpol,'--',color='red')
    if plot_arcmancer:
        for ipha in range(0,len(phase_arc)-1):
            if(phase_arc[ipha+1] > phase_arc[ipha]):    
                ax.plot(phase_arc[ipha:ipha+2],pulse_arcQ[ipha:ipha+2],'--',color='blue',dashes=[2,2])    	
    ax2 = fig.add_subplot(212)
    ax2.plot(phasepol,(qnpol-Q1pn)/qnpol,'--',color='red')        
    #ax2.plot(phasepol,(qnpol-Q1pn),'--',color='red')
    ax2.set_ylabel('Residuals [%]')        
    #ax2.set_ylabel('Residuals [arbitrary]')        
    ax2.set_xlabel('Phase [cycles]')
    fig.savefig("figs/signalsQ_residX.pdf",bbox_inches='tight')   
    
    plt.subplots_adjust(wspace=0, hspace=0)
    fig = plt.figure(figsize=(7,7))
    ax = fig.add_subplot(211)
    ax.set_ylabel('Signal [arbitrary normalisation]')  
    U1p = np.sum(photU_sig_cut, axis=0) 
    U1pn = np.copy(U1p)
    for ip in range(len(U1pn)):
    	if(I1p[ip] > 1e-10):
    		U1pn[ip] = U1p[ip]/I1p[ip]
    	else:
    		U1pn[ip] = 0.0  
    ax.plot(signals[0][0].phases[0], U1pn, 'o-', color='k', lw=0.5, markersize=2)
    if (phasepol is not None and unpol is not None):
    	ax.plot(phasepol,unpol,'--',color='red')
    if plot_arcmancer:
        for ipha in range(0,len(phase_arc)-1):
            if(phase_arc[ipha+1] > phase_arc[ipha]):    
                ax.plot(phase_arc[ipha:ipha+2],pulse_arcU[ipha:ipha+2],'--',color='blue',dashes=[2,2])         	
    ax2 = fig.add_subplot(212)
    ax2.plot(phasepol,(unpol-U1pn)/unpol,'--',color='red')    
    #ax2.plot(phasepol,(unpol-U1pn),'--',color='red') 
    ax2.set_ylabel('Residuals [%]')               
    #ax2.set_ylabel('Residuals [arbitrary]')    
    ax2.set_xlabel('Phase [cycles]')
    fig.savefig("figs/signalsU_residX.pdf",bbox_inches='tight')      
     






