import numpy as np
from astropy.io import fits

## Reads simulated data (for a single detector unit) made with ixpeobssim pcube method
NPhadat = 10 # number of phases in .fits
def readData_pcube(Filename):

	I = np.zeros((NPhadat))
	Ierr = np.zeros((NPhadat))
	counts_all = np.zeros((NPhadat))
	W2_all = np.zeros((NPhadat))	

	Q = np.zeros((NPhadat))
	Qerr = np.zeros((NPhadat))

	U = np.zeros((NPhadat))
	Uerr = np.zeros((NPhadat))
	
	PDEG = np.zeros((NPhadat))
	PDEG_ERR = np.zeros((NPhadat))
	
	PANG = np.zeros((NPhadat))
	PANG_ERR = np.zeros((NPhadat))

	mu_all = np.zeros(NPhadat)
	emean = np.zeros((NPhadat))

	phase_points = np.linspace(0.0,1.0,NPhadat+1)
	phase = np.zeros((NPhadat))
	for i in range(0,NPhadat):
		phase[i] = (phase_points[i+1]+phase_points[i])/2.0

	for p in range(NPhadat):
		if(p<10):
			hdulist = fits.open(str(Filename) + '_folded_phase000' + str(p) + '_pcube2' + '.fits')
		else:
			hdulist = fits.open(str(Filename) + '_folded_phase00' + str(p) + '_pcube2' + '.fits')
						
		cols1 = hdulist[1].columns
		mu_all[p] = hdulist[1].data["EFFECTIVE_MU"][0]
		Id = hdulist[1].data["I"][0]
		Qd = hdulist[1].data["Q"][0]
		Ud = hdulist[1].data["U"][0]
		PDEG[p] = hdulist[1].data["POL_DEG"][0]
		PDEG_ERR[p] = hdulist[1].data["POL_DEG_ERR"][0]
		PANG[p] = hdulist[1].data["POL_ANG"][0]
		PANG_ERR[p] = hdulist[1].data["POL_ANG_ERR"][0]	
		W2_all[p] = hdulist[1].data["W2"][0]	
		counts_all[p] = hdulist[1].data["COUNTS"][0]
		emean[p] = hdulist[1].data["ENERGY_MEAN"][0]									
			
		I[p] = Id
		Q[p] = Qd/Id
		U[p] = Ud/Id

		Ierr[p] = 1.0 #This is actually not used anywhere at the moment
				
		#See Kislat 2015 formulas A.9a and A.9b.				
		Qerr[p] = np.sqrt((W2_all[p]/(I[p]**2))*(2.0/mu_all[p]**2-Q[p]**2))
		Uerr[p] = np.sqrt((W2_all[p]/(I[p]**2))*(2.0/mu_all[p]**2-U[p]**2))
	
	I = np.nan_to_num(I)
	Q = np.nan_to_num(Q)
	U = np.nan_to_num(U)
	Ierr = np.nan_to_num(Ierr)
	Qerr = np.nan_to_num(Qerr)
	Uerr = np.nan_to_num(Uerr)
	
	keV = np.mean(emean)
	
	return phase, np.array([I]), np.array([Q]), np.array([U]), np.array([Ierr]), np.array([Qerr]), np.array([Uerr]), np.array([keV])

def readData_pcube_combined(Filename):

	NDet = 3 # number of detectors used

	I = np.zeros((NPhadat))
	Ierr = np.zeros((NPhadat))
	counts_all = np.zeros((NPhadat))
	W2_all = np.zeros((NPhadat))	

	Q = np.zeros((NPhadat))
	Qerr = np.zeros((NPhadat))

	U = np.zeros((NPhadat))
	Uerr = np.zeros((NPhadat))
	
	PDEG = np.zeros((NPhadat))
	PDEG_ERR = np.zeros((NPhadat))
	
	PANG = np.zeros((NPhadat))
	PANG_ERR = np.zeros((NPhadat))

	mu_all = np.zeros(NPhadat)
	emean = np.zeros((NPhadat,NDet))

	phase_points = np.linspace(0.0,1.0,NPhadat+1)
	phase = np.zeros((NPhadat))
	for i in range(0,NPhadat):
		phase[i] = (phase_points[i+1]+phase_points[i])/2.0

	for p in range(NPhadat):
		for d in range(1, NDet+1):
			if(p<10):
				hdulist = fits.open(str(Filename) + '_du' + str(d) + '_folded_phase000' + str(p) + '_pcube2' + '.fits')
				cols1 = hdulist[1].columns
				mu = hdulist[1].data["EFFECTIVE_MU"][0]
				#print(p, mu)
				mu_all[p] = mu_all[p] + mu
				Id = hdulist[1].data["I"][0]
				Qd = hdulist[1].data["Q"][0]#/mu[p]
				Ud = hdulist[1].data["U"][0]#/mu[p]
				Pd = hdulist[1].data["POL_DEG"][0]
				Pd_err = hdulist[1].data["POL_DEG_ERR"][0]
				Pad = hdulist[1].data["POL_ANG"][0]
				Pad_err = hdulist[1].data["POL_ANG_ERR"][0]	
				W2 = hdulist[1].data["W2"][0]							
					
				I[p] = I[p] + Id
				Q[p] = Q[p] + Qd
				U[p] = U[p] + Ud
				
				W2_all[p] = W2_all[p]+W2
				
				PDEG[p] = PDEG[p]+Pd
				PANG[p] = PANG[p]+Pad
				PDEG_ERR[p] = PDEG_ERR[p] + Pd_err**2
				PANG_ERR[p] = PANG_ERR[p] + Pad_err**2 
				counts = hdulist[1].data["COUNTS"][0]
				emean[p,d-1] = hdulist[1].data["ENERGY_MEAN"][0]

				counts_all[p] = counts_all[p]+counts

				Ierr[p] = 1.0 #This is actually not used anywhere at the moment			
				
												
			else:
				print("ERROR: We assume no more than 10 phase bins!")
				exit()
		#After sum over detectors finished
		Q[p] = Q[p]/I[p] #/NDet
		U[p] = U[p]/I[p] #/NDet
		PDEG[p] = PDEG[p]/NDet
		PANG[p] = PANG[p]/NDet		
		
		mu_all[p] = mu_all[p]/NDet
		
		#See Kislat 2015 formulas A.9a and A.9b.				
		Qerr[p] = np.sqrt((W2_all[p]/(I[p]**2))*(2.0/mu_all[p]**2-Q[p]**2))
		Uerr[p] = np.sqrt((W2_all[p]/(I[p]**2))*(2.0/mu_all[p]**2-U[p]**2))	
				
		PDEG_ERR[p] = np.sqrt(PDEG_ERR[p])/NDet 
		PANG_ERR[p] = np.sqrt(PANG_ERR[p])/NDet 

	
	I = np.nan_to_num(I)
	Q = np.nan_to_num(Q)
	U = np.nan_to_num(U)
	Ierr = np.nan_to_num(Ierr)
	Qerr = np.nan_to_num(Qerr)
	Uerr = np.nan_to_num(Uerr)
	
	keV = np.mean(emean)
	
	return phase, np.array([I]), np.array([Q]), np.array([U]), np.array([Ierr]), np.array([Qerr]), np.array([Uerr]), np.array([keV])

	
def read_response_IXPE(MRF,RMF,min_input,max_input,min_channel,max_channel):

	hdulist_mrf = fits.open(MRF)
	#cols1 = hdulist_mrf[1].columns
	#print(cols1.info())
	specresp = hdulist_mrf[1].data["SPECRESP"]
	ene_lo = hdulist_mrf[1].data["ENERG_LO"]
	ene_hi = hdulist_mrf[1].data["ENERG_HI"]
	
	hdulist_rmf = fits.open(RMF)	
	matrix = hdulist_rmf[1].data["MATRIX"]
	emin = hdulist_rmf[2].data["E_MIN"]
	emax = hdulist_rmf[2].data["E_MAX"]		
	#cols1 = hdulist_rmf[2].columns

        #matrix = np.ascontiguousarray(RMF[min_input:max_input,20:201].T, dtype=np.double)
        matrix_cut = np.ascontiguousarray(matrix[min_input:max_input,min_channel:max_channel].T, dtype=np.double)
        #print("matrix_cut:")
        #print(matrix_cut)

        #print(ene_lo[min_input:max_input])
        #print(ene_hi[min_input:max_input])
        #print(emin[min_channel:max_channel])
        #print(emax[min_channel:max_channel])

	#print(specresp[min_input:max_input].shape[0]+1)

        edges = np.zeros(specresp[min_input:max_input].shape[0]+1, dtype=np.double)
        edges[0] = ene_lo[min_input]; edges[1:] = ene_hi[min_input:max_input]

        for i in range(matrix_cut.shape[0]):
            matrix_cut[i,:] *= specresp[min_input:max_input]

        #print("matrix_cut:")
        #print(matrix_cut)
        #channels = np.arange(20, 201)
        
	channel_edges = np.zeros(matrix[min_channel:max_channel,0].shape[0]+1, dtype=np.double)
        channel_edges[0] = emin[min_channel]; channel_edges[1:] = emax[min_channel:max_channel]	

	#print(edges)
	#print(channel_edges)
        #exit()

	channels = np.arange(min_channel,max_channel)
	
	#Re-bin channels if necessary for the data product:
	rebin = True
	pcube = True
	if rebin:
		if pcube:
			#Assuming Nchan = 1
			channels = np.array([0])
			channel_edges = np.array([2.0,8.0])#([4.0,8.0])
			matrix_rb = np.zeros((1,len(matrix_cut[0,:])))
			#Calculating just the average here
			for ich in range(0,max_channel-min_channel):
				matrix_rb[0,:] += matrix_cut[ich,:]
			matrix_rb[0,:] = matrix_rb[0,:]/(max_channel-min_channel)
			matrix_cut = matrix_rb	
		else:
			print("Other re-binning options to be implemented.")
			exit()	

	print("matrix_cut=",matrix_cut, len(matrix_cut[0,:]), len(matrix_cut[:,0]))
	return matrix_cut, edges, channels, channel_edges
	
	
	
