##To plot polarization position angle profile comparisons between our Obl+Schw and arcmancer

#import:
import matplotlib
matplotlib.use('Agg')
from numpy import linspace, logspace, empty, zeros, ones, array, fromfile
from numpy import pi, exp, log, sqrt, sin, cos, arccos, arctan2
from numpy import absolute, sign, floor, ceil
from numpy.polynomial.laguerre import laggauss
from numpy.polynomial.legendre import leggauss
from scipy.interpolate import interp1d
from scipy.special import kn
from matplotlib.pyplot import *
from bisect import bisect
#import find_best_phshift

def shift_phase(phi,shift):
	return (phi + shift) % 1 


#physical constants:
evere=.5109989e6 # electron volts in elecron rest energy 
G=13275412528e1 # G*M_sol in km^3/s^2 
c=299792458e-3 # speed of light in km/s

NEnergy = 128 #281 # 50# 101 # number of energy points (x)
#x_l, x_u = -3.7 , .3 
#IntEnergy = logspace(x_l,x_u,NEnergy), log(1e1)*(x_u-x_l)/(NEnergy-1.)
#x_ene,x_weight=IntEnergy
#energy_keV = x_ene*evere/1e3
energy_keV = [ 1., 1.01975887,  1.03990815,  1.06045555,  1.08140895,  1.10277637,
  1.12456598,  1.14678613,  1.16944532,  1.19255224,  1.21611572,  1.24014479,
  1.26464864,  1.28963666,  1.31511842,  1.34110367,  1.36760236,  1.39462463,
  1.42218084,  1.45028152,  1.47893744,  1.50815956,  1.53795909,  1.56834742,
  1.59933618,  1.63093725,  1.66316273,  1.69602494,  1.72953647,  1.76371015,
  1.79855906,  1.83409655,  1.87033622,  1.90729194,  1.94497787,  1.98340843,
  2.02259833,  2.06256258,  2.10331648,  2.14487563,  2.18725594,  2.23047364,
  2.27454527,  2.31948771,  2.36531816,  2.41205416,  2.45971362,  2.50831477,
  2.55787623,  2.60841697,  2.65995633,  2.71251405,  2.76611025,  2.82076546,
  2.87650059,  2.93333698,  2.99129639,  3.05040102,  3.11067348,  3.17213687,
  3.2348147,   3.29873097,  3.36391015,  3.43037721,  3.49815757,  3.5672772,
  3.63776256,  3.70964062,  3.78293892,  3.8576855 ,  3.93390899,  4.01163858,
  4.09090401,  4.17173564,  4.2541644,   4.33822187,  4.42394022,  4.51135226,
  4.60049147,  4.69139197,  4.78408856,  4.87861672,  4.97501266,  5.07331327,
  5.17355619,  5.2757798,   5.38002323,  5.48632639,  5.59472998,  5.70527551,
  5.81800528,  5.93296247,  6.05019109,  6.16973601,  6.291643, 6.41595873,
  6.54273081,  6.67200775,  6.80383906,  6.93827521,  7.07536767,  7.21516891,
  7.35773247,  7.50311293,  7.65136594,  7.80254826,  7.95671777,  8.11393349,
  8.27425562,  8.43774553,  8.60446582,  8.77448032,  8.9478541,   9.12465356,
  9.30494637,  9.48880157,  9.67628953,  9.86748204, 10.06245231, 10.26127496,
 10.46402612, 10.67078342, 10.88162601, 11.0966346,  11.31589153, 11.53948072,
 11.76748778, 12.  ]



figA = figure(figsize=(14,30), dpi=300)
matplotlib.pyplot.figure(1)
lbfontsz = 25#35#30#25 
lwidth= 2.5#3.5#3.0#2.5#2.0#1.5 
wpad = 10
rc("xtick", labelsize=lbfontsz)
rc("ytick", labelsize=lbfontsz)
rc("axes", linewidth=lwidth)
#figA.clear()
matplotlib.pyplot.rcParams.update({'axes.titlesize': lbfontsz})
matplotlib.pyplot.rcParams.update({'font.size': lbfontsz})
matplotlib.pyplot.rcParams.update({'lines.linewidth': lwidth})
matplotlib.pyplot.rcParams.update({'ytick.major.width': lwidth})
matplotlib.pyplot.rcParams.update({'xtick.major.width': lwidth})
matplotlib.pyplot.rcParams.update({'ytick.major.size': 10.0})
matplotlib.pyplot.rcParams.update({'xtick.major.size': 10.0})
matplotlib.pyplot.rcParams.update({'font.family': 'serif'})
#matplotlib.pyplot.rcParams.update({'font.serif': 'Times'})


plot_only_I = False
plot_all = True
plot_QU = False#True

plot5 = True

matplotlib.pyplot.subplots_adjust(wspace=0, hspace=0)


if(plot_all):
	if not(plot5):
		plotAF=figA.add_subplot(4,1,1,yscale='linear') 
		plotAp=figA.add_subplot(4,1,2)   #
		plotAc=figA.add_subplot(4,1,3)   #
		plotAd=figA.add_subplot(4,1,4)   #
	if(plot5):
		plotAFF=figA.add_subplot(5,1,1,yscale='linear')
		plotAF=figA.add_subplot(5,1,2)
		plotAp=figA.add_subplot(5,1,3)
		plotAc=figA.add_subplot(5,1,4)
		plotAd=figA.add_subplot(5,1,5)   
  
else:	
	plotAc=figA.add_subplot(2,1,1)   #
	plotAd=figA.add_subplot(2,1,2)   #



colors = ["black","blue"]
shapes = np.copy(colors)


for ish in range(1,-1,-1): 

	oblateness=shapes[ish]
	print(ish)
	if(ish == 0):
		#PulsName='pulses/xpatap_rho10f600_Tc_'
		#PulsName='pulses/xpatap_rho10f600_Tc_281_pshift_match_'
		#PulsName='pulses/all_E_281/xpatap_ps21_thom_s21'
		#PulsName='pulses/all_E_128a/xpatap_ps21_burst_s21'
		PulsName='pulses/all_E_128a/xpatap_ps21_thom_s21'		
		NPhase = 150 #121 #150
	if(ish == 1):
	    #pversion = '_test_26052022_thom_corr2'
	    #pversion = 'xpsi_ps21_burstN_s21'
	    pversion = 'xpsi_ps21_thom_s21'
	    #PulsName='pulses/all_E_281/'+pversion		
	    PulsName='pulses/all_E_128a/'+pversion
	    NPhase = 150 #121 #100
	computePulse= True
	plotAtm=not True
	plotPulse=True
	mod=True

	print(PulsName)
	eind = 42 #93 #42 #82 #118 #95

	if(ish==0):
		inFlux = open(PulsName+'FF.bin')
		inphi = open(PulsName+'ff.bin')
		Flux1 = fromfile(inFlux)
		fluxlcurve0 = Flux1[0:len(Flux1):3*NEnergy] #light curve with lowest E
		fluxspec0 = Flux1[0:3*NEnergy:3] #spectrum at phase=0
		ene = eind #82 #118 #82 #118#166#140#166#140 #The chosen energy index
		#print("The chosen energy (keV): ", x_ene[ene]*evere/1e3)
		print("The chosen energy (keV): ", energy_keV[ene])
		fluxlcurve_Iene = Flux1[0+ene*3:len(Flux1):3*NEnergy]
		fluxlcurve_Qene = Flux1[1+ene*3:len(Flux1):3*NEnergy]
		fluxlcurve_Uene = Flux1[2+ene*3:len(Flux1):3*NEnergy]

	else:
		inFlux = open(PulsName+'_F'+str(eind)+'.bin')
		inphi = open(PulsName+'_p'+str(eind)+'.bin')
		inQ = open(PulsName+'_Q'+str(eind)+'.bin')
		inU = open(PulsName+'_U'+str(eind)+'.bin')

		Flux1 = fromfile(inFlux)
		FluxQ = fromfile(inQ)
		FluxU = fromfile(inU)

	phi = fromfile(inphi)

   	   
	#Flux=zeros((NPhase,NEnergy,3))
	#if(ish==0):
	#	print(fluxlcurve_Iene)
	#	print(fluxlcurve_Qene)
	#	print(fluxlcurve_Uene)
	
	labelsize=40#30#20
	fontsize=35#50#35#25
	ticksize=25


	#phase=list(phi/2/pi)+[1.]
	phase=list(phi)#+[1.]
	I=zeros(NPhase)#+1)
	Q=zeros(NPhase)#+1)
	U=zeros(NPhase)#+1)
	for t in range(NPhase):#+1):
		if(ish==0):
			#I[t],Q[t],U[t]=fluxlcurve_Iene[t-1],fluxlcurve_Qene[t-1],fluxlcurve_Uene[t-1]
			#Or when t-1 not used in polpulse (did not test though): 
			I[t],Q[t],U[t]=fluxlcurve_Iene[t],fluxlcurve_Qene[t],fluxlcurve_Uene[t]
		else:
			I[t],Q[t],U[t]=Flux1[t],FluxQ[t],FluxU[t]

	#exit()
	p=sqrt(Q**2+U**2)/I*100
	#PA=arctan2(-U,-Q)*90/pi+90
	PA=arctan2(U,Q)*90/pi+90


	if not(plot_only_I):
		if(plot_QU):
			print("plot_QU option not valid in this version!")
			quit()	
		else:

			plotAc.set_xlim(0,1)
			plotAc.set_ylim(0,180)
			plotAc.tick_params(axis='both', which='major', labelsize=ticksize,direction='in',top=True,right = True)
			plotAc.set_ylabel(r'$\chi\,[\mathrm{deg}]$',fontsize=fontsize)

			if(plot_all):
				plotAp.set_xlim(0,1)
				plotAp.tick_params(axis='both', which='major', labelsize=ticksize,direction='in',top=True,right = True)
				plotAp.set_ylabel(r'$p\,[ \% ]$',fontsize=fontsize)
				plotAp.set_ylabel(r'$\delta F_{\mathrm{Q,U}} / F_{\mathrm{Q,U}}$',fontsize=fontsize)
				plotAF.set_xlim(0,1)
				if(plot5):
					plotAFF.set_xlim(0,1)
					plotAFF.tick_params(axis='both', which='major', labelsize=ticksize,direction='in',top=True,right = True)
					plotAF.set_ylabel(r"$F_{\mathrm{Q,U}}/F_{\mathrm{I}}$",fontsize=fontsize)
					plotAFF.set_ylabel(r"$F_{\mathrm{I}}/F_{\mathrm{I}}^{\mathrm{max}}$",fontsize=fontsize)
				else:
					plotAF.set_ylabel(r"$F_{\mathrm{x}}$",fontsize=fontsize)
				plotAd.tick_params(axis='both', which='major', labelsize=ticksize,direction='in',top=True,right = True)
				plotAd.set_xlim(0,1)
				plotAF.tick_params(axis='both', which='major', labelsize=ticksize,direction='in',top=True,right = True)

	col = colors[ish]

	PA_VP04 = PA
	phase_VP04 = phase
	if(ish == 1): 		
		phshift1 = 0.0 #-0.00213 #-0.04829 #-0.047590550687164 #0.0 #-0.12 #0.0 #-0.047590550687164 #0.0#-0.048315#-0.2517#0.019#0.0#0.019#0.195#-0.07#-0.2517#0.0#0.001#0.008#0.2421#0.2517#0.2535#0.069#0.0#-0.195#-0.18#-0.172#0.0
		phase_new = shift_phase(np.array(phase),phshift1)
		for ipha in range(0,len(phase_new)-1):
			if(phase_new[ipha+1] > phase_new[ipha]):
				#print(PA[ipha],PA[ipha+1])
				if(PA[ipha] > 100.0*PA[ipha+1]):
					plotAc.plot(phase_new[ipha:ipha+2],[PA[ipha],180.0],"-",color="blue")
					#print("pass!")
					##pass
				else:
					plotAc.plot(phase_new[ipha:ipha+2],PA[ipha:ipha+2],"-",color="blue")#,markersize="1.0")
				if(plot_all):
					if plot5:
						plotAFF.plot(phase_new[ipha:ipha+2],I[ipha:ipha+2]/I.max(),color=col)
						plotAF.plot(phase_new[ipha:ipha+2],Q[ipha:ipha+2]/I[ipha:ipha+2],color="red") 
						plotAF.plot(phase_new[ipha:ipha+2],U[ipha:ipha+2]/I[ipha:ipha+2],color="darkorange")
					else:
						plotAF.plot(phase_new[ipha:ipha+2],I[ipha:ipha+2]/I.max(),color=col)
						plotAF.plot(phase_new[ipha:ipha+2],Q[ipha:ipha+2]/I[ipha:ipha+2],color="red")
						plotAF.plot(phase_new[ipha:ipha+2],U[ipha:ipha+2]/I[ipha:ipha+2],color="darkorange")
		PA0_VP04 = PA
		F0_VP04 = I/I.max()
		Q0_VP04 = Q/I
		U0_VP04 = U/I
		phase0_VP04 = phase

	if(ish == 0): 
		phshift2 = 0.0 #0.00217 #0.0#-0.006715 #0.00217 #0.047590550687164#0.0 #-0.047590550687164 #0.0 #-0.04829 #0.0
		phase_new = shift_phase(np.array(phase),phshift2)
		print(len(phase_new))
		for ipha in range(0,len(phase_new)-1):
			if(phase_new[ipha+1] > phase_new[ipha]):
				plotAc.plot(phase_new[ipha:ipha+2],PA[ipha:ipha+2],"--",color=col,dashes=[2,2])
				if(plot_all):
					if plot5:
						plotAFF.plot(phase_new[ipha:ipha+2],I[ipha:ipha+2]/I.max(),"--",color=col,dashes=[2,2])
						plotAF.plot(phase_new[ipha:ipha+2],Q[ipha:ipha+2]/I[ipha:ipha+2],"--",color=col,dashes=[2,2]) 
						plotAF.plot(phase_new[ipha:ipha+2],U[ipha:ipha+2]/I[ipha:ipha+2],"--",color=col,dashes=[2,2])
					else:
						print("This option is not supported!")
						exit()
		PA0_VP04_2 = PA
		PA_acm0 = PA
		phase_acm0 = phase_new
		phase_acm0[len(phase_acm0)-1]=1.0
		F_acm0 = I/I.max()
		Q_acm0 = Q/I
		U_acm0 = U/I

		F0_VP04_2 = I/I.max()
		Q0_VP04_2 = Q/I
		U0_VP04_2 = U/I
		phase0_VP04_2 = phase



compare_to_arcmancer = False #True
if(compare_to_arcmancer): 
	#colors = ["green","blue","black"]
	colors = ["black","black","black"]
	for ic in range(0,1):
		if(ic == 0):
			datafile = "../../../polcslab/CompSlab/arcman_res/polar_acc_f600_burst_r12_m1.4_d60_i40_x10_obl.csv"# (copy).csv"
		if(ic == 1):
			datafile = "../arcmancer/out3/polar_f001_bb_r12_m1.4_d40_i60_x01_sph.csv"
		#input = file(datafile, 'r')
		input = open(datafile, 'r')
		lines = input.readlines()
		input.close()

		Nchain_size = sum(1 for line in open(datafile))
		c_lines = 1
		egrid = 6#5#3
		full_chain= [[] for x in range(egrid+1)]

		input = open(datafile, 'r')
		lines = input.readlines()
		input.close()

		for j in range(0,len(full_chain)):
			for i in range(c_lines,Nchain_size): #not reading comment lines
				parts = lines[i].split(",")
				#print parts
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

		ene = 0
		print("energy_keV = ", energy_keV[ene])

		phshift = 0.0#0.2517#0.25#0.2421#0.2517#0.2535#0.069#0.0#-0.195#-0.18#-0.172#0.0
		phase_new = shift_phase(phase,phshift)
		fluxsp1 = full_chain[1,:]


		#phase_new = shift_phase(phase,-phshift)
		col = colors[ic]
		if not(plot_only_I):
			if(plot_QU):
				print("plot_QU option not valid in this version!")
				quit()	
			else:
				p = full_chain[5,:]*100.0
				PA = full_chain[4,:]*180.0/pi+90.0


				use_PA_avg = False#True
				if use_PA_avg:
					PA = full_chain[6,:]*180/pi+90.0#full_chain[4,:]*180.0/pi+90.0
					for qwe in range(0,len(PA)):
						if(PA[qwe] < 0.0):
							PA[qwe] = PA[qwe]+180.0
						if(PA[qwe] > 180.0):
							PA[qwe] = PA[qwe]-180.0
				for ipha in range(0,len(phase_new)-1):
					if(phase_new[ipha+1] > phase_new[ipha]):
						if(PA[ipha] > 100.0*PA[ipha+1]):
							plotAc.plot(phase_new[ipha:ipha+2],[PA[ipha],180.0],"--",color=col,dashes=[2,2])  
						else:
							plotAc.plot(phase_new[ipha:ipha+2],PA[ipha:ipha+2],"--",color=col,dashes=[2,2])
				if(ic == 0):
					PA_acm0 = PA#norm_obsF[:,ene]#PA
					F_acm0 = norm_obsF[:,ene]#obsF[:,ene]
					Q_acm0 = norm_obsF[:,ene+1]#norm_obsF[:,ene+1]
					U_acm0 = norm_obsF[:,ene+2]#norm_obsF[:,ene+2]
					phase_acm0 = phase_new
					phase_acm0[len(phase_acm0)-1]=1.0

				if(ic == 1):
					PA_acm = PA#norm_obsF[:,ene]#PA
					phase_acm = phase_new
					phase_acm[len(phase_acm)-1]=1.0



		if(plot_all):

			#for i in range(ene,ene+1):
			sind=0
			if(plot5):
				sind=1
				for ipha in range(0,len(phase_new)-1):
					if(phase_new[ipha+1] > phase_new[ipha]):
						plotAFF.plot(phase_new[ipha:ipha+2],norm_obsF[ipha:ipha+2,0],"--",color=colors[0],dashes=[2,2])
			for i in range(sind,ene+3): #sind=0 if plotting all stokes fluxes to same plot
				for ipha in range(0,len(phase_new)-1):
					if(phase_new[ipha+1] > phase_new[ipha]):
						plotAF.plot(phase_new[ipha:ipha+2],norm_obsF[ipha:ipha+2,i],"--",color=colors[i],dashes=[2,2])
  


plot_PA_residuals = True
if(plot_QU):
	plot_PA_residuals = False
if(plot_PA_residuals): 
	col = "black"
	plotAd.set_xlim(0,1)
	plotAd.set_ylim(-0.3,0.3)

	plotAd.tick_params(axis='both', which='major', labelsize=ticksize,direction='in')
	plotAd.tick_params(axis='x', which='major', labelsize=labelsize)
	plotAd.set_ylabel(r'$\chi_{\mathrm{arc}}-\chi_{\mathrm{obl}} \ [\mathrm{deg}]$',fontsize=fontsize)

	PA_VP04_interp = interp1d(phase_VP04,PA_VP04)#,"cubic")
	PA0_VP04_interp = interp1d(shift_phase(np.array(phase0_VP04),phshift1),PA0_VP04,fill_value='extrapolate')# workd with newer scipy


	norm = 1.0#abs(np.max(PA_acm0)-np.min(PA_acm0))

	res_PA = (PA_acm0-PA0_VP04_interp(phase_acm0))/norm
	#print(res_PA)
	plotAd.plot(phase_acm0[abs(res_PA) < 5.0],res_PA[abs(res_PA) < 5.0],"-",markersize=5,color="black")#"red")

	print("sum_PA_err=",sum(res_PA[abs(res_PA)<5.0]))

	plotAp.set_ylim(-0.03,0.03)

	F0_VP04_interp = interp1d(shift_phase(np.array(phase0_VP04),phshift1),F0_VP04,fill_value='extrapolate')# workd with newer scipy
	Q0_VP04_interp = interp1d(shift_phase(np.array(phase0_VP04),phshift1),Q0_VP04,fill_value='extrapolate')# workd with newer scipy
	U0_VP04_interp = interp1d(shift_phase(np.array(phase0_VP04),phshift1),U0_VP04,fill_value='extrapolate')# workd with newer scipy
	res_F = (F_acm0-F0_VP04_interp(phase_acm0))/F_acm0#abs((F_acm0-F0_VP04_interp(phase_acm0))/F_acm0)
	res_Q = (Q_acm0-Q0_VP04_interp(phase_acm0))/Q_acm0
	res_U = (U_acm0-U0_VP04_interp(phase_acm0))/U_acm0
	plotAp.plot(phase_acm0,res_F,"-",markersize=5,color="blue")
	plotAp.plot(phase_acm0[abs(res_Q) < 0.02],res_Q[abs(res_Q) < 0.02],"-",markersize=5,color="red")#"darkblue")
	plotAp.plot(phase_acm0[abs(res_U) < 0.02],res_U[abs(res_U) < 0.02],"-",markersize=5,color="darkorange")#"lightblue")




plotAFF.margins(x=0,y=0)
plotAF.margins(x=0,y=0)
plotAp.margins(x=0,y=0)
plotAc.margins(x=0,y=0)
plotAd.margins(x=0)

plotAF.xaxis.set_major_formatter(matplotlib.pyplot.NullFormatter())
plotAp.xaxis.set_major_formatter(matplotlib.pyplot.NullFormatter())
plotAc.xaxis.set_major_formatter(matplotlib.pyplot.NullFormatter())
plotAd.set_xlabel(r'$\varphi / (2\pi)$',fontsize=fontsize)

#These are already defined elsewhere and top=True,right = True used there (this is used just to adjust labelsizes and pads):
plotAF.tick_params(axis='both', which='major', labelsize=fontsize,pad=wpad)
plotAp.tick_params(axis='both', which='major', labelsize=fontsize,pad=wpad)
plotAc.tick_params(axis='both', which='major', labelsize=fontsize,pad=wpad)
plotAd.tick_params(axis='both', which='major', labelsize=fontsize,pad=wpad)

figA.subplots_adjust(wspace=0, hspace=0)

figA.subplots_adjust(left=0.15)


#align manually instead:
labelx = -0.14#-0.16
plotAFF.yaxis.set_label_coords(labelx, 0.5)
plotAF.yaxis.set_label_coords(labelx, 0.5)
plotAp.yaxis.set_label_coords(labelx, 0.5)
plotAc.yaxis.set_label_coords(labelx, 0.5)
plotAd.yaxis.set_label_coords(labelx, 0.5)


plotAF.set_yticks([-0.04,-0.02,0.0,0.02])
plotAF.set_yticklabels(["-0.04","-0.02","0","0.02"],fontstyle="normal")
plotAF.set_ylim(-0.05,0.03)

plotAp.set_yticks([-0.02,0.0,0.02])
plotAp.set_yticklabels(["-0.02","0","0.02"],fontstyle="normal")

plotAd.set_yticks([-0.2,0.0,0.2])
plotAd.set_yticklabels(["-0.2","0","0.2"],fontstyle="normal")


if(plot5):
	plotAFF.xaxis.set_major_formatter(matplotlib.pyplot.NullFormatter())
	plotAFF.tick_params(axis='both', which='major', labelsize=fontsize, pad=wpad)
	plotAFF.set_ylim(0.0,1.1)
	plotAFF.set_yticks([0.0,0.2,0.4,0.6,0.8,1.0])
	plotAFF.set_yticklabels(["0","0.2","0.4","0.6","0.8","1.0"],fontstyle="normal")

	plotAc.set_yticks([0,50,100,150])
	plotAc.set_yticklabels(["0","50","100","150"],fontstyle="normal")
	figA.subplots_adjust(0.15)#(left=0.175)
else:
	figA.subplots_adjust(left=0.15)

figA.savefig('figs/pulse_compX.pdf',bbox_inches='tight')#.format(e))
figA.clf()




