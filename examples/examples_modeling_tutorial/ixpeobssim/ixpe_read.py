'''
Helper functions to read IXPE data and response files.
These scripts are still to be checked and likely developed further.
'''


import numpy as np
from astropy.io import fits

NPhadat = 10 # number of phases in .fits

def readData_pcube_ebin(Filename):
      NDet = 3  # number of detectors used
      NEnerg = 1 # number of energy bins in the data

      I = np.zeros((NPhadat,NEnerg))
      Ierr = np.zeros((NPhadat,NEnerg))
      counts_all = np.zeros((NPhadat,NEnerg))
      W2_all = np.zeros((NPhadat,NEnerg))

      Q = np.zeros((NPhadat,NEnerg))
      Qerr = np.zeros((NPhadat,NEnerg))

      U = np.zeros((NPhadat,NEnerg))
      Uerr = np.zeros((NPhadat,NEnerg))

      PDEG = np.zeros((NPhadat,NEnerg))
      PDEG_ERR = np.zeros((NPhadat,NEnerg))

      PANG = np.zeros((NPhadat,NEnerg))
      PANG_ERR = np.zeros((NPhadat,NEnerg))

      MDP_tot = np.zeros((NPhadat))

      mu_all = np.zeros((NPhadat,NEnerg))
      emean = np.zeros((NPhadat, NDet, NEnerg))
      elow = np.zeros((NPhadat, NDet, NEnerg))
      ehigh = np.zeros((NPhadat, NDet, NEnerg))

      phase_points = np.linspace(0.0, 1.0, NPhadat + 1)
      phase = np.zeros((NPhadat))
      for i in range(0, NPhadat):
            phase[i] = (phase_points[i + 1] + phase_points[i]) / 2.0

      for p in range(NPhadat):
            for d in range(1, NDet + 1):
                  if (p < 10):
                        hdulist = fits.open(f'{Filename}_du{d}_folded_phase000{p}_pcube.fits')
                  elif (p < 100):
                        hdulist = fits.open(f'{Filename}_du{d}_folded_phase00{p}_pcube.fits')
                  else:
                        print("ERROR: We assume no more than 100 phase bins!")
                        exit()
                  for ie in range(0,NEnerg):
                        cols1 = hdulist[1].columns
                        #print(cols1)
                        #exit()
                        mu = hdulist[1].data["MU"][ie]
                        mu_all[p] = mu_all[p] + mu
                        Id = hdulist[1].data["I"][ie]
                        Qd = hdulist[1].data["Q"][ie]
                        Ud = hdulist[1].data["U"][ie]
                        Pd = hdulist[1].data["PD"][ie]
                        Pd_err = hdulist[1].data["PD_ERR"][ie]
                        Pad = hdulist[1].data["PA"][ie]
                        Pad_err = hdulist[1].data["PA_ERR"][ie]
                        W2 = hdulist[1].data["W2"][ie]

                        Uerrd = hdulist[1].data["U_ERR"][ie]
                        Qerrd = hdulist[1].data["Q_ERR"][ie]
                        Ierrd = hdulist[1].data["I_ERR"][ie]

                        I[p,ie] = I[p,ie] + Id
                        Q[p,ie] = Q[p,ie] + Qd
                        U[p,ie] = U[p,ie] + Ud

                        Ierr[p,ie] = Ierr[p,ie] + Ierrd**2
                        Qerr[p,ie] = Qerr[p,ie] + Qerrd**2
                        Uerr[p,ie] = Uerr[p,ie] + Uerrd**2

                        W2_all[p,ie] = W2_all[p,ie] + W2

                        PDEG[p,ie] = PDEG[p,ie] + Pd
                        PANG[p,ie] = PANG[p,ie] + Pad
                        PDEG_ERR[p,ie] = PDEG_ERR[p,ie] + Pd_err ** 2
                        PANG_ERR[p,ie] = PANG_ERR[p,ie] + Pad_err ** 2
                        counts = hdulist[1].data["COUNTS"][ie]
                        emean[p, d - 1,ie] = hdulist[1].data["E_MEAN"][ie]
                        elow[p, d - 1, ie] = hdulist[1].data["ENERG_LO"][ie]
                        ehigh[p, d - 1, ie] = hdulist[1].data["ENERG_HI"][ie]

                        counts_all[p,ie] = counts_all[p,ie] + counts

                        #Ierr[p,ie] = Ierr[p,ie] + np.sqrt(Id) ** 2
      # After sum over detectors finished
      Q = Q / I
      U = U / I
      PDEG = PDEG / NDet
      PANG = PANG / NDet
      
      
      mu_all[p] = mu_all[p] / NDet

      # See Kislat 2015 formulas A.9a and A.9b.
      #Qerr = np.sqrt((W2_all / (I ** 2)) * (2.0 / mu_all ** 2 - Q ** 2))
      #Uerr = np.sqrt((W2_all / (I ** 2)) * (2.0 / mu_all ** 2 - U ** 2))

      Ierr = np.sqrt(Ierr)
      Qerr = np.sqrt(Qerr)/I
      Uerr = np.sqrt(Uerr)/I
      
      I = np.nan_to_num(I)
      Q = np.nan_to_num(Q)
      U = np.nan_to_num(U)
      Ierr = np.nan_to_num(Ierr)
      Qerr = np.nan_to_num(Qerr)
      Uerr = np.nan_to_num(Uerr)

      PDEG_ERR = np.sqrt(PDEG_ERR)/NDet
      MDP_tot = 4.29*np.sqrt(W2_all)/(mu_all*I)

      #keV = np.mean(emean, axis=1)
      #keV = np.mean(keV, axis=0)

      elow = np.mean(elow, axis=1)
      elow = np.mean(elow, axis=0)

      ehigh = np.mean(ehigh, axis=1)
      ehigh = np.mean(ehigh, axis=0)

      ebinning_data = np.append(elow[0], ehigh)

      return phase, I, Q, U, Ierr, Qerr, Uerr, PDEG, PDEG_ERR, ebinning_data, MDP_tot


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

        matrix_cut = np.ascontiguousarray(matrix[min_input:max_input,min_channel:max_channel].T, dtype=np.double)

        edges = np.zeros(specresp[min_input:max_input].shape[0]+1, dtype=np.double)
        edges[0] = ene_lo[min_input]; edges[1:] = ene_hi[min_input:max_input]

        for i in range(matrix_cut.shape[0]):
                matrix_cut[i,:] *= specresp[min_input:max_input]

        channel_edges = np.zeros(matrix[min_channel:max_channel,0].shape[0]+1, dtype=np.double)
        channel_edges[0] = emin[min_channel]; channel_edges[1:] = emax[min_channel:max_channel]

        channels = np.arange(min_channel,max_channel)

        #Re-bin channels if necessary for the data product:
        rebin = True
        pcube = True
        if rebin:
                if pcube:
                        #Assuming Nchan = 1
                        channels = np.array([0])
                        channel_edges = np.array([2.0,8.0])
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
