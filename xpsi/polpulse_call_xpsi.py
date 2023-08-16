#!/usr/bin/python
#This is the main program to calculate polarized pulse profiles (to be run with driver.py or to be called by the sampling method in newmc.py).
#Set SavePulse=True if want to save the output in a file
#This version has an option to use simplThomson atmosphere model as in Salmi+ 2021, or reading precomputed models.
#This code is based on cs_ts2_func.py found in https://github.com/thjsal/CompSlab 
#and originally on cs.py found in https://github.com/belliavesha/CompSlab (developed by Vladislav Loktev).
#Small differences in implemention and computation accuracy may exist between the different versions of the code.


Spectrum={
      2:'Burst',
      1:'simplThomson', 
      0:'FromFile'
}[2]#[1]

oblateness='AlGendy'

#import:
from numpy import linspace, logspace, empty, zeros, ones, array, fromfile
import matplotlib
matplotlib.use('agg')
from numpy import pi, exp, log, sqrt, sin, cos, arccos, arctan2
from numpy import absolute, sign, floor, ceil, argmin
from numpy.polynomial.laguerre import laggauss
from numpy.polynomial.legendre import leggauss
from scipy.interpolate import interp1d#,CubicSpline 
from scipy.interpolate import CubicSpline 
from scipy.interpolate import interp2d
from scipy.special import kn
from matplotlib.pyplot import *
from bisect import bisect


def find_idx(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx

def Planck(x,T):
	constbb= 5.039617e22 
	ex=exp(-x/T)
	evere=.5109989e6 # electron volts in elecron rest energy
	keV = (x*evere)/1e3
	I=constbb*keV**3*ex/(1.0-ex) #in erg cm^-2 s^-1 str^-1 keV^-1

	return I

#Beaming:
def angledep(mu, NZenith):
	a = -0.7
	b = 0
	c= 0.5 /(0.5+1*a/3+1*b/4)
	mudep = zeros(NZenith)
	for i in range(NZenith):
		mudep[i]=c*(1+a*mu[i]+b*mu[i]**2)
	sumi=0
	return mudep

#Beaming function based on the Thomson slab model (read from a file):
def angledep_thom(mu, NZenith, x, e, AtmNamex):
	inI = open(AtmNamex+'I.bin')
	inx = open(AtmNamex+'x.bin')
	inm = open(AtmNamex+'m.bin')
	x2=fromfile(inx)
	mu=fromfile(inm)
	NEnergy=len(x2)
	NZenith=len(mu)
	NMu=int(NZenith/2)
	Intensity=fromfile(inI).reshape((NEnergy,NZenith,2))
	Intensity[:,:,0] = np.nan_to_num(Intensity[:,:,0])
	mudep = zeros(NZenith)
	for d in range(NZenith):
		fI = interp1d(x2, Intensity[:,d,0], kind="linear")
		I = fI(x)
		mudep[d] = I


	sumi=0 #normalize beaming so that \int over \muI\dmu = 0.5 as for isotropic case
	for d in range(NMu,NZenith):
		if(d==NMu):
			dmu = ((mu[d+1]+mu[d])/2.0)-((mu[d]+0.0)/2.0)
		elif(d==NZenith-1):
			dmu = ((1.0+mu[d])/2.0)-((mu[d]+mu[d-1])/2.0)
		else:
			dmu = ((mu[d+1]+mu[d])/2.0)-((mu[d]+mu[d-1])/2.0) 
		sumi = sumi+mu[d]*mudep[d]*dmu
	norm=1.0/(sumi*2.0)
	mudep=mudep*norm

	return mudep
    
def simpl(NEnergy, ear, T):
	seed = zeros(NEnergy)
	x = ear
	ergkev = 1.0/(1.0e3*1.602176565e-12)
     
        
	for i in range(0, NEnergy):
		seed[i] = Planck(x[i], T)

	gamma = 1.8
	switch = 0
	norm = 0.6

	tmparr = zeros(NEnergy)
	enavgam1 = zeros(NEnergy)
	enavgam2 = zeros(NEnergy)
	engam1 = zeros(NEnergy)
	engam2 = zeros(NEnergy)
	photar_sp1 = zeros(NEnergy)
	photar_sp2 = zeros(NEnergy)

	if(gamma==1):
		gamma=1.001

	gamma1=gamma-1
	gamma2=gamma+2

	for i in range(0, NEnergy):
		engam1[i]=x[i]**(-gamma1)
		if(i!=0):
			tmparr[i]=0.
			enavgam1[i]=(0.5*(x[i-1]+x[i]))**gamma1


	if(switch > 0):
		for i in range(1, NEnergy):
			tmparr[i]=tmparr[i]+seed[i]*(1.-enavgam1[i]*engam1[i])
			for j in range(i+1, NEnergy):
				tmparr[j]=tmparr[j]+enavgam1[i]*(engam1[j-1]-engam1[j])*seed[i]

	
	if(switch <= 0):
		gnormUP=(gamma+2.)/(1.+2.*gamma)
		gnormDN=(gamma-1.)/(1.+2.*gamma)
		for k in range(0, NEnergy):
			engam2[k]=x[k]**[gamma2]
			if(k!=0):
				enavgam2[k]=(0.5*(x[k-1]+x[k]))**(-gamma2)

		for i in range(1,NEnergy):
			tmparr[i]=tmparr[i]+seed[i]*((1.-enavgam1[i]*engam1[i])*gnormUP+gnormDN*(1.-engam2[i]*enavgam2[i]))

			for j in range(1,NEnergy):
				if(j<i):
					tmparr[j]=tmparr[j]+enavgam2[i]*(engam2[j]-engam2[j-1])*seed[i]*gnormDN
				if(j>i):
					tmparr[j]=tmparr[j]+enavgam1[i]*(engam1[j-1]-engam1[j])*seed[i]*gnormUP

	for i in range(0, NEnergy):
		photar_sp1[i]=(1.-norm)*seed[i]
		photar_sp2[i]=norm*tmparr[i]

	photar_sp2[0] = photar_sp2[1]*.9
        
	return photar_sp1, photar_sp2

def poldeg(pol, mu):
	return pol*(mu - 1.)/(1. + 3.582*abs(mu))

def poldegfromfile(x, e, d, AtmNamex):
	inI = open(AtmNamex+'I.bin')
	inx = open(AtmNamex+'x.bin')
	inm = open(AtmNamex+'m.bin')
	x2=fromfile(inx)
	mu=fromfile(inm)
	NEnergy=len(x2)
	NZenith=len(mu)
	NMu=NZenith/2
	Intensity=fromfile(inI).reshape((NEnergy,NZenith,2))

	Intensity[:,:,0] = np.nan_to_num(Intensity[:,:,0])
	Intensity[:,:,1] = np.nan_to_num(Intensity[:,:,1])

	fI = interp1d(x2, Intensity[:,d,0], kind="linear")
	I = fI(x)
	fQ = interp1d(x2, Intensity[:,d,1], kind="linear")
	Q = fQ(x)

	
	if(I[e]<1e-10):
		return 0

	p = Q[e]/I[e]
	return p

def compf(mass,eqrad,incl_deg,theta_deg,rho_deg,pol,ekev,ph,spherical=False,antipodal=False,spath="pulse_test_00",savePulse=False,atmos_path="model/atmos_thom/"):
	PulsName=spath

	#physical constants:
	evere=.5109989e6 # electron volts in elecron rest energy 
	G=13275412528e1 # G*M_sol in km^3/s^2 
	c=299792458e-3 # speed of light in km/s

	# Atmosphere parameters: 
	tau_T= 1.0 # Thomson optical depth of thermalization 

	T = 0.002 # 10/evere #  dimensionless photon black body temperature T = k T_bb / m_e c^2 #~  1.0219978 keV

	#precomputations : Done already beforehand

	NMu = 22# 20# 15 # number of propagation zenith angle cosines (\mu) [0,1]
	NZenith = 2*NMu # number of propagation zenith angles (z) [0,pi]
	#IntEnergy = logspace(x_l,x_u,NEnergy), log(1e1)*(x_u-x_l)/(NEnergy-1.) # sample points and weights for integrations over the spectrum computing sorce function
	IntZenith = leggauss(NZenith) #  sample points and weights for integrations over zenith angle in positive and negative directions together   

	mu,mu_weight=IntZenith
	
	mu = np.append(mu,1.0)
	#mu=cos(linspace(-pi/2,pi/2,num=NZenith))
	#print(mu)
	#exit()

	#x_l, x_u = -3.7 , -1.2 # -3.7 , .3 # -3.7, -1.2 # lower and upper bounds of the log_10 energy span
	#NEnergy =  281 # number of energy points (x)	
	#x,x_weight=IntEnergy
	#keV = (x*evere)/1e3
	keV = ekev
	x = (ekev*1e3)/evere
	NEnergy = len(ekev)

	if Spectrum=='simplThomson' : # Initializing Stokes vectors arrays, computing zeroth scattering 
		Intensity=zeros((NEnergy,NZenith,2)) # total intensity of all scattering orders from the slab suface 
		photar_sp1, photar_sp2 = simpl(NEnergy, x, T)
		xNorm = 1

		mudep = (1+2.06*mu*(pol/0.1171))*exp(-tau_T/mu)
		mudepc = np.ones((len(mu))) #exp(-tau_T/mu) #np.ones((len(mu))) #beaming used as comparison case
		sumi=0 #normalize beaming so that \int over \muI\dmu = 0.5 as for isotropic case or for some other preference
		#cbeam=0
		for d in range(NMu,NZenith):
			if(d==NMu):
				dmu = ((mu[d+1]+mu[d])/2.0)-((mu[d]+0.0)/2.0)
			elif(d==NZenith-1):
				dmu = ((1.0+mu[d])/2.0)-((mu[d]+mu[d-1])/2.0)
			else:
				dmu = ((mu[d+1]+mu[d])/2.0)-((mu[d]+mu[d-1])/2.0)
			sumi = sumi+mu[d]*mudep[d]*dmu
			#cbeam = cbeam+mu[d]*mudepc[d]*dmu
		cbeam=0.5 #for isotropic
		Normb0=cbeam/sumi

		mudep_n0 = Normb0*(1.0+2.06*mu*(pol/0.1171))*exp(-tau_T/mu) 


		for e in range(NEnergy):

			if(pol<0.001):
				pdstr="0"
			else:
				pdstr="1171"
			#AtmName_angdep = '../../../CompSlab/res/ixpe_spec_final/tau10_te01/grid_res_pbi_NN_no0/atmos_thom_p'+pdstr   
			#AtmName_angdep = atmos_path+'atmos_thom_n0_p'+pdstr+'_corrPlanck' #angdep from file not including zero scattered photons
			AtmName_angdep = atmos_path+'atmos_thom_n0_p'+pdstr    #angdep from file not including zero scattered photons			
			mudep_thom = angledep_thom(mu,NZenith,x[e],e,AtmName_angdep)
			for d in range(NZenith):
				Intensity[e,d,0]=photar_sp1[e-1]*mudep_n0[d]+photar_sp2[e-1]*mudep_thom[d]
				#Intensity[e,d,1]=Intensity[e,d,0]*poldegfromfile(x,e,d,atmos_path+'atmos_thom_y0_p'+pdstr+'_corrPlanck') #PD from model including also zeroth scattering
				Intensity[e,d,1]=Intensity[e,d,0]*poldegfromfile(x,e,d,atmos_path+'atmos_thom_y0_p'+pdstr) #PD from model including also zeroth scattering

		#AtmName = "/home/tuomo/polcslab/X-PATAP/x-patap/analysis/model/atmos_nsx_like/atmos_thom_csformat_from_callxpsi"
		#outI = open(AtmName+'I.bin','w')
		#outx = open(AtmName+'x.bin','w')
		#outm = open(AtmName+'m.bin','w')
		#Intensity.tofile(outI,format="%e")
		#x.tofile(outx,format="%e")
		#mu.tofile(outm,format="%e")

				


	if Spectrum=='FromFile' :
		#AtmName='../../../CompSlab/res/ixpe_spec_final/tau10_te01/atmos_thom_p11'
		AtmName='/home/tuomo/polcslab/X-PATAP/x-patap/analysis/model/atmos_nsx_like/atmos_thom_csformat'#_from_callxpsi'
		inI = open(AtmName+'I.bin')
		inx = open(AtmName+'x.bin')
		inm = open(AtmName+'m.bin')
		xi=fromfile(inx)
		mui=fromfile(inm)
		NEnergyi=len(xi)
		NZenithi=len(mui) #assumed to be same as NZenith
		Intensityf=fromfile(inI).reshape((NEnergyi,NZenithi,2))
		Intensity=zeros((NEnergy,NZenith,2))
                
		Intensityf[:,:,0] = np.nan_to_num(Intensityf[:,:,0])
		Intensityf[:,:,1] = np.nan_to_num(Intensityf[:,:,1])

		for d in range(0,NZenith):                
			fI = interp1d(xi, Intensityf[:,d,0], kind="linear")
			I = fI(x)
			fQ = interp1d(xi, Intensityf[:,d,1], kind="linear")
			Q = fQ(x)
			Intensity[:,d,0]=I
			Intensity[:,d,1]=Q



	if Spectrum=='Burst' : # Initializing Stokes vectors arrays, computing zeroth scattering 
		Intensity=zeros((NEnergy,NZenith+1,2)) # total intensity of all scattering orders from the slab surface
		for e in range(NEnergy):
			for d in range(NZenith+1):
				#Intensity[e,d,0]=Planck(x[e],T)#*(1 + 2.06*mu[d])#TS TESTING BLACKBODY energy spectrum with burst-beaming only in polarization
				Intensity[e,d,0]=Planck(x[e],T)*(0.421+0.868*mu[d]) #more accurate version
				Intensity[e,d,1]=Intensity[e,d,0]*0.1171*(mu[d] - 1.)/(1. + 3.582*mu[d])
				#if e==42:
				#	print("Q,mu:",Intensity[e,d,1],mu[d])


	NPhi = 120 #500 #120 # Number of equidistant phase points
	NBend= 20 # Number of knots in light bending integrations
	NAlpha= 200#1000 # 10000 # Number of psi/aplha grid points 
	IntBend = leggauss(NBend)
	NZenithBig=5001 #101 #100
	#NZenithBig = NZenith
	bending_exact = False #True
	tdelay_exact = True #False #True #False #True

	phi=linspace(0,2*pi,num=NPhi,endpoint=False,retstep=False)
	#exit()

	#NPhase = 150 #500# 150 # Number of observation phases
	#phase =linspace(0,1,num=NPhase,endpoint=True,retstep=False)

	phase = ph
	NPhase = len(ph)

	phase_obs=zeros(NPhi)
	nu=401.0#1.0#100 #600 # star rotation frequency in Hz
	#M=1.4 # star mass in solar masses
	M=mass #input param
	R_g=M*2.95325 #2.95325024 #2.95325 # gravitational Schwarzschild radius #TS: Made this more accurate
	#R_e=12.0 # equatorial radius of the star in kilometers
	R_e=eqrad #input param

	#Increase the resolution for the following when considering large spots!:
	NRho=1 #20 #1 #4 #4 #20 #4 #40 #4#40#20#4#2#8
	NVarphi=1 #20 #1 #2 #6 #20 #2 #40 #2#40#20#6#4

	# IntVarphi = linspace(0,2*pi,num=NVarphi,endpoint=False,retstep=True)
	IntVarphi = leggauss(NVarphi)
	IntRho = leggauss(NRho)
      
	if oblateness=='AlGendy': # from AlGendy et. al. (2014)
		Omega_bar=2*pi*nu*sqrt(2*R_e**3/R_g)/c
		flattening=(0.788-0.515*R_g/R_e)*Omega_bar**2 
	elif oblateness=='Sphere':
		flattening=0.0
	else:
		flattening=oblateness

	if(spherical):#TS: from input param, 0.0 not working in this code
		flattening = 1e-8

	def Beloborodov(cos_psi):
	    """Beloborodov's approximation for cos_alpha(cos_psi) light bending function
	    takes the cos psi 
	    returns the cos alpha and its derivative
	    """
	    return 1. + (cos_psi - 1.)/redshift**2 ,1./redshift**2
	    
	    
	def Poutanen(u,y):
		return ( 1 - u )*y*( 1 + u*u*y*y/112 - np.e/1e2*u*y*( np.log( 1 - y/2 ) + y/2 ) )

	def Poutanen_der(u,y):
		der_end = (3/2) - (0.5/(1-(y/2)))
		return ( 1 - u )*( 1 + 3*u*u*y*y/112 - np.e/1e2*u*y*(2*np.log( 1 - y/2 ) + y*der_end))		    

	def tlag_PB06(u,y):
		#return y*(1.0+u/8.0*y*(1.0+y*(1.0/3.0-u/14.0)))
		#keeping just the leading term, since accuracy not improved by including more
		return y
		
		 

	def Schwarzschild(R,alpha):
	    """Schwarzschild exact relation between the \psi and \\alpha angles, where
	    \\alpha is the angle between radius vector of the spot and the direction of the outgoing photon near the surface
	    and \psi is the angle between normal and light propagation at the limit of infinite distance.
	    For given distance from the mass center and the emission angle \\alpha 
	    this function returns two numbers: 
	          the corresponding angle \psi 
	          and the time lag over against the fotons emited with zero impact parameter at the radius.
	    """
	    kx,wx=IntBend
	    eps=(1+kx[0])/4e2
	    u=R_g/R 
	    b=sin(alpha)/sqrt(1-u)*R # impact parameter
	    if 2*alpha>pi+eps:
	          cos_3eta=sqrt(27)*R_g/2/b
	          if cos_3eta > 1:
	                return pi+2*eps,0 # the timelag 
	          closest_approach=-2*b/sqrt(3)*cos(arccos(cos_3eta)/3 + 2*pi/3)
	          psi_max, lag_max= Schwarzschild(closest_approach,pi/2.)
	          psi_min, lag_min= Schwarzschild(R,pi-alpha)
	          psi=2*psi_max - psi_min    
	          lag=2*lag_max - lag_min # + 2*(R - closest_approach + R_g*log((R - R_g)/(closest_approach - R_g)))/c 
	          if psi>pi:
	                return pi+eps,lag
	    else:
	          psi=0
	          lag=(R_e - R + R_g*log( (R_e - R_g)/(R - R_g) ) )/c
	          for i in range(NBend):
	                ex=(kx[i]+1)/2
	                q=(2. - ex*ex - u*(1 - ex*ex)**2/(1 - u))*sin(alpha)**2
	                sr=sqrt(cos(alpha)**2+ex*ex*q)
	                if  2*alpha>pi-eps:
	                      dpsi=b/R/sqrt(q)*wx[i] #*2/2
	                else:
	                      dpsi=ex*b/R/sr*wx[i] #*2/2
	                dlag=dpsi*b/c/(1+sr) #*2/2
	                psi+= dpsi
	                lag+= dlag
	    return psi,lag


	#flattening=0
	NRadius=2 + int(flattening*R_e/1e-1)
	NRadius=4+ int(flattening*R_e/1e-1)
	r, dr = linspace(R_e*(1 - flattening),R_e,num=NRadius,retstep=True)
	    
	if(bending_exact or tdelay_exact):    
		alpha, dalpha = linspace(0,arccos(-1/sqrt(2*r[0]/R_g/3)),NAlpha,retstep=True)
		psi=zeros((NRadius,NAlpha))
		dt=zeros((NRadius,NAlpha))
		for d in range(NRadius):
			#print(d)
			for a in range(NAlpha):
				psi[d,a],dt[d,a]=Schwarzschild(r[d],alpha[a])



	Flux=zeros((NPhase,NEnergy,3))
	Flux_obs=zeros((NPhi,NEnergy,3))

	i=pi/180.0*incl_deg#pi*7/18    # line of sight colatitude

	#Integration over spot:
	rho_total=rho_deg*pi/180.0#pi/5 #pi*5/180 # radius of the spot
	theta_center=pi/180.0*theta_deg#pi/4.1

	NSpots=0
	varphi,dvarphi=IntVarphi[0]*pi,IntVarphi[1] *pi
	rho,drho=(IntRho[0]+1)*rho_total/2,(IntRho[1])*rho_total/2

	l=[]
	theta=[]
	dS=[]

	for v in range(NVarphi):
		for rh in range(NRho):
			NSpots+=1   
			cos_theta=cos(theta_center)*cos(rho[rh])+sin(theta_center)*sin(rho[rh])*cos(varphi[v])
			sin_l=sin(rho[rh])*sin(varphi[v])/sqrt(1- cos_theta**2)
			cos_l=sqrt(1- sin_l**2)
			if cos_theta*cos(theta_center)> cos(rho[rh]) : 
				cos_l=-cos_l      
			l.append(arctan2(-sin_l,-cos_l) + pi)
			theta.append(arccos(cos_theta))  
			dS.append(drho[rh]*dvarphi[v]*sin(rho[rh]))
			#print(v,rh,cos_theta,cos_l,sin_l,l[-1],theta[-1],dS[-1])
			if antipodal :
				NSpots+=1
				l.append(arctan2(sin_l,cos_l) + pi)
				theta.append(pi- theta[-1])
				dS.append(dS[-1])
			#print(v,rh,cos_theta,cos_l,sin_l,l[-1],theta[-1],dS[-1])


	sin_i=sin(i)
	cos_i=cos(i)

	BoloFlux=zeros((NPhase,3))
	z=cos(linspace(-pi/2,pi/2,num=NZenithBig))
	logIntensity=zeros((NEnergy,NZenithBig,3))



	#print("mu:",mu[NMu:])
	#exit()

	for e in range(NEnergy):
		IntInt=CubicSpline(mu[NMu:],Intensity[e,NMu:,0],extrapolate=True) # interpolate intensity
		IQ=CubicSpline(mu[NMu:],Intensity[e,NMu:,1],extrapolate=True) #




		for d in range(NZenithBig):
			#print("IQ",IQ(z[d]))
			logIntensity[e,d] = log(max(0,IntInt(z[d]))),log(absolute(IQ(z[d]))),sign(IQ(z[d]))
			#Testing with non-log interpolation:
			#logIntensity[e,d] = max(0,IntInt(z[d])),IQ(z[d]),sign(IQ(z[d]))
			#if e ==42:
			#	print("Q,mu=",logIntensity[e,d][1],z[d])

	mu=z.copy()
	#exit()

	for p in range(NSpots):
		sin_theta=sin(theta[p])
		cos_theta=cos(theta[p])

		R=R_e*(1 - flattening*cos_theta**2) 
		dR=2*R_e*flattening*cos_theta*sin_theta # dR / d\theta

		r1=bisect(r[1:-1],R) 
		r2=r1 + 1
		dr1=(R - r[r1])/dr
		dr2=(r[r2] - R)/dr

		redshift=1.0/sqrt(1.0 - R_g/R) # 1/sqrt(1-R_g/R) = 1+ z = redshift
		f=redshift/R*dR
		sin_gamma=f/sqrt(1 + f**2) # angle gamma is positive towards the north pole 
		cos_gamma=1.0/sqrt(1 + f**2)
		beta=2*pi*nu*R*redshift*sin_theta/c
		Gamma=1.0/sqrt(1.0 - beta**2)
		Gamma1= (1.0-sqrt(1.0 - beta**2) )/ beta

		u=R_g/R

		for t in range(NPhi):
			if True: # find mu
				phi0=phi[t]+l[p]
				sin_phi=sin(phi0)
				cos_phi=cos(phi0)
				cos_psi=cos_i*cos_theta + sin_i*sin_theta*cos_phi
				sin_psi=sqrt(1. - cos_psi**2)

				if(bending_exact or tdelay_exact):
					psi0=arccos(cos_psi) 
					a1=bisect(psi[r1], psi0)
					a2=bisect(psi[r2], psi0)
					if(a1 >= len(psi[r1])):
						a1=len(psi[r1])-1
					if(a2 >= len(psi[r2])):
						a2=len(psi[r2])-1
					psi1=psi[r1,a1]					
					psi2=psi[r2, a2]
					dpsi1=psi1 - psi[r1, a1 - 1]
					dpsi2=psi2 - psi[r2, a2 - 1]
					dpsi=dpsi1*dr2 + dpsi2*dr1
					dalpha1 = dalpha*(psi1 - psi0)/dpsi1
					dalpha2 = dalpha*(psi2 - psi0)/dpsi2
					alpha1=alpha[a1] - dalpha1
					alpha2=alpha[a2] - dalpha2
						
				if(bending_exact):
	
					cos_alpha = cos(alpha2*dr1 + alpha1*dr2) # linear interpolation of alpha(psi)
					sin_alpha = sqrt(1. - cos_alpha**2)
					sin_alpha_over_sin_psi= sin_alpha/sin_psi if sin_psi > 1e-10 else 1./redshift
					#if not sin_psi > 1e-10:
					#    print("sin_psi:",sin_psi)
					dcos_alpha=sin_alpha_over_sin_psi *dalpha/dpsi # d cos\alpha \over d \cos \psi
					
				else:
					# cos_alpha, dcos_alpha=Beloborodov(cos_psi)
					# sin_alpha = sqrt(1. - cos_alpha**2)
					# sin_alpha_over_sin_psi= sin_alpha/sin_psi if sin_psi > 1e-4 else 1./redshift
					
					cos_alpha = 1.0 - Poutanen(u, 1.0 - cos_psi) 
					sin_alpha = sqrt(1. - cos_alpha**2)	
					sin_alpha_over_sin_psi= sin_alpha/sin_psi if sin_psi > 1e-10 else 1./redshift
					#if not sin_psi > 1e-10:
					#    print("HELLO2")
					dcos_alpha = Poutanen_der(u, 1.0 - cos_psi)
				#print("sin_alpha, sin_psi:",sin_alpha,sin_psi)
				#exit()
				#sin_alpha = sin_psi
				#cos_alpha = cos_psi
				#dcos_alpha = 1.0
				#sin_alpha_over_sin_psi = 1.0

				if(tdelay_exact):
					dt1=dt[r1,a1 - 1]*dalpha1/dalpha + dt[r1,a1]*(1. - dalpha1/dalpha)
					dt2=dt[r2,a2 - 1]*dalpha2/dalpha + dt[r2,a2]*(1. - dalpha2/dalpha)
					dphase=(dt1*dr2 + dt2*dr1)*nu # \delta\phi = \phi_{obs} - \phi 
				else:						
					dphase = tlag_PB06(u, 1.0 - cos_psi)*nu*R/c
				
				#dphase = 0
				phase_obs[t]=( phi[t]/2/pi+dphase)%1.

				cos_xi = - sin_alpha_over_sin_psi*sin_i*sin_phi
				delta = 1./Gamma/(1.-beta*cos_xi)
				#if(p==0):
				#	print(phase_obs[t], " ",delta, ", ")
				cos_sigma = cos_gamma*cos_alpha + sin_alpha_over_sin_psi*sin_gamma*(cos_i*sin_theta - sin_i*cos_theta*cos_phi)

				sin_sigma = sqrt(1. - cos_sigma**2)
				mu0=delta*cos_sigma # cos(sigma')
				Omega=dS[p]*mu0*redshift**2*dcos_alpha*Gamma*R*R/cos_gamma


				#print(t,' : \t',mu0,' \t ',dcos_alpha,'\t',cos_alpha,cos_psi,Omega)
				if mu0<0: # this only for speeding up. the backwards intensity is usually zero
					Flux_obs[t]=0
					continue 


				if True: # find chi
					sin_chi_0= - sin_theta*sin_phi # times sin psi
					cos_chi_0=sin_i*cos_theta - sin_theta*cos_i*cos_phi # times sin psi 
					chi_0=arctan2(sin_chi_0,cos_chi_0)

					sin_chi_1=sin_gamma*sin_i*sin_phi*sin_alpha_over_sin_psi #times sin alpha sin sigma 
					cos_chi_1=cos_gamma - cos_alpha*cos_sigma  #times sin alpha sin sigma 
					chi_1=arctan2(sin_chi_1,cos_chi_1)

					sin_lambda=sin_theta*cos_gamma - sin_gamma*cos_theta
					cos_lambda=cos_theta*cos_gamma + sin_theta*sin_gamma
					cos_eps = sin_alpha_over_sin_psi*(cos_i*sin_lambda - sin_i*cos_lambda*cos_phi + cos_psi*sin_gamma) - cos_alpha*sin_gamma
						
					#OR THE ORIGINAL VERSION:
					sin_chi_prime=cos_eps*mu0*Gamma*beta#*the_thing#*delta**3#*delta
					cos_chi_prime=(1. - cos_sigma**2 /(1. - beta*cos_xi))#*the_thing

					chi_prime=arctan2(sin_chi_prime,cos_chi_prime)   

					chi=chi_0 + chi_1 + chi_prime
					#  print(chi,'\t',chi_0/chi,'\t',chi_1/chi ,'\t',  chi_prime/chi )

					sin_2chi=sin(2*chi)
					cos_2chi=cos(2*chi)


			d2=bisect(mu[:-1],mu0)
			d1=d2-1
			mu1,mu2=mu[d1],mu[d2]
			dmu, dmu1, dmu2 = mu2 - mu1, mu0 - mu1, mu2 - mu0
			shift=delta/redshift


			for e in range(NEnergy): 
				x0=x[e]/shift
				e1=bisect(x[1:-1],x0) # not the fastest way? anybody cares? ## seems, that light bending is more time consuming anyways
				e2=e1+1

				x1, x2 = x[e1], x[e2]
				dx, dx1, dx2 = x2 - x1, x0 - x1, x2 - x0
				logIQ = (
				      dx2*dmu2*logIntensity[e1, d1] + 
				      dx2*dmu1*logIntensity[e1, d2] +
				      dx1*dmu2*logIntensity[e2, d1] +
				      dx1*dmu1*logIntensity[e2, d2] # bilinear interpolation of the Stokes parameters
				)/dx/dmu 
				I,Q=exp(logIQ[:2])* shift**3 * Omega
				Q*=logIQ[2]
				#Testing with non-log interpolation:
				#I,Q=logIQ[:2]* shift**3 * Omega
				#Testing with not interpolating at all:
				#I = Planck(x0,T)*(0.421+0.868*mu0)*shift**3 * Omega
				#Q = I*0.1171*(mu0 - 1.)/(1. + 3.582*mu0)
				#Qtest = (0.421+0.868*mu0)*0.1171*(mu0 - 1.)/(1. + 3.582*mu0)

				if I<0: ############
					print('never')
				Flux_obs[t,e]=[I, Q*cos_2chi, Q*sin_2chi]
				if e==42:
				    #print("Q0,Qtest,Q:",Q0,Qtest,Q)
				    print("mu1,mu0,mu2",mu1,mu0,mu2)
				    #print("t, Q, mu1, mu2=",t,Q,mu1,mu2,logIntensity[e2, d1],logIntensity[e2, d2])

		for t in range(NPhase):
			phase0=phase[t]


			#A newer version from Vlad:
			#for t2 in range(NPhi):
			#	t1=t2-1
			#	phase2=phase_obs[t2]
			#	phase1=phase_obs[t1]
			#	dphase10 = (phase1-phase0+0.5)%1-0.5
			#	dphase20 = (phase2-phase0+0.5)%1-0.5
			#	if dphase10<0.0<dphase20 :
			#		break
			#Flux[t]+=(Flux_obs[t1]*dphase20-Flux_obs[t2]*dphase10)/(dphase20-dphase10)

			for t2 in range(NPhi):
				t1=t2-1
				phase2=phase_obs[t2]
				phase1=phase_obs[t1]
				if phase0>phase1 and phase0<phase2:
					break
			else :
				if(phase0>max(phase_obs)):
					phase0=phase0-1
				t2=argmin(phase_obs)
				t1=t2-1
				phase2=phase_obs[t2]
				phase1=phase_obs[t1]-1

			dphase1=phase0-phase1
			dphase2=phase2-phase0
			dphase=phase2-phase1
			Flux[t]+=(Flux_obs[t2]*dphase1+Flux_obs[t1]*dphase2)/dphase 
			#TS: changed here t to t-1 to match the phaseshift with fortran or X-PSI results

	if savePulse:
		outF = open(PulsName + 'FF.bin','w')
		outf = open(PulsName + 'ff.bin','w')
		Flux.tofile(outF,format="%e")
		phase.tofile(outf,format="%e")

	return Flux



