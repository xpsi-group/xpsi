import numpy as np
import matplotlib.pyplot as plt
import h5py as h5
from matplotlib.ticker import NullFormatter
import sys
import scipy

import os
import math
import time

from matplotlib import rcParams
from matplotlib.ticker import MultipleLocator, AutoLocator, AutoMinorLocator
from matplotlib import gridspec
from matplotlib import cm

import xpsi


rc = {"font.family" : "serif",
      "mathtext.fontset" : "stix"}
plt.rcParams.update(rc)
plt.rcParams["font.serif"] = ["Times New Roman"] + plt.rcParams["font.serif"]

ticksize = 15
legsize = 15
labelsize = 15
ticksize = ticksize*1.5
legsize = legsize*1.5
labelsize = labelsize*1.3

rcParams['text.usetex'] = False
rcParams['font.size'] = 14.0


class IpyExit(SystemExit):
    """Exit Exception for IPython.

    Exception temporarily redirects stderr to buffer.
    """
    def __init__(self):
        # print("exiting")  # optionally print some message to stdout, too
        # ... or do other stuff before exit
        sys.stderr = StringIO()

    def __del__(self):
        sys.stderr.close()
        sys.stderr = sys.__stderr__  # restore from backup



def veneer(x, y, axes, lw=1.0, length=8):
    """ Make the plots a little more aesthetically pleasing. """
    if x is not None:
        if x[1] is not None:
            axes.xaxis.set_major_locator(MultipleLocator(x[1]))
        if x[0] is not None:
            axes.xaxis.set_minor_locator(MultipleLocator(x[0]))
    else:
        axes.xaxis.set_major_locator(AutoLocator())
        axes.xaxis.set_minor_locator(AutoMinorLocator())

    if y is not None:
        if y[1] is not None:
            axes.yaxis.set_major_locator(MultipleLocator(y[1]))
        if y[0] is not None:
            axes.yaxis.set_minor_locator(MultipleLocator(y[0]))
    else:
        axes.yaxis.set_major_locator(AutoLocator())
        axes.yaxis.set_minor_locator(AutoMinorLocator())

    for tick in axes.xaxis.get_major_ticks():
        tick.label1.set_fontsize(ticksize)
    for tick in axes.yaxis.get_major_ticks():
        tick.label1.set_fontsize(ticksize)

    axes.tick_params(which='major', colors='black', length=length, width=lw)
    axes.tick_params(which='minor', colors='black', length=int(length/2), width=lw)
    plt.setp(axes.spines.values(), linewidth=lw, color='black')


REFl= [   '__phase_shift',
          '__super_colatitude',
          '__super_radius',
          '__super_temperature',
          '__omit_azimuth',
          '__omit_colatitude',
          '__omit_radius',
          '__cede_azimuth',
          '__cede_colatitude',
          '__cede_radius',
          '__cede_temperature']

ABBl = ['PS','SC','SR','ST','OA','OC','OR','CA','CC','CR','CT']

REF = dict(zip(ABBl,REFl))

def transform(thetaR,phiR,V,phi0=0.0):
        phi0 = phi0*np.pi*2.

        RA = [-np.sin(phi0),np.cos(phi0),0.]
        V = np.matrix(V)
        phiR = phiR*2*np.pi

        Mphi     = np.matrix([[np.cos(phiR),-np.sin(phiR),0.],
                              [np.sin(phiR),np.cos(phiR),0.],
                              [0.,0.,1.]]) #_z
        Mtheta = np.matrix([[np.cos(thetaR)+ RA[0]*RA[0]*(1.-np.cos(thetaR)),
                             RA[0]*RA[1]*(1.-np.cos(thetaR)) - RA[2]*np.sin(thetaR),
                             RA[0]*RA[2]*(1.-np.cos(thetaR)) + RA[1]*np.sin(thetaR),],
                            [RA[0]*RA[1]*(1.-np.cos(thetaR)) + RA[2]*np.sin(thetaR),
                             np.cos(thetaR)+ RA[1]*RA[1]*(1.-np.cos(thetaR)),
                             RA[1]*RA[2]*(1.-np.cos(thetaR)) - RA[0]*np.sin(thetaR),],
                            [RA[0]*RA[2]*(1.-np.cos(thetaR)) - RA[1]*np.sin(thetaR),
                             RA[1]*RA[2]*(1.-np.cos(thetaR)) + RA[0]*np.sin(thetaR),
                            np.cos(thetaR)+ RA[2]*RA[2]*(1.-np.cos(thetaR))]])

        Mtheta[np.abs(Mtheta)<1e-16] = 0.0
        M = Mphi*Mtheta
        if np.shape(V)[0] ==1:
            V = np.transpose(V)
        Vout = M*V
        Vout = np.squeeze(np.array(Vout))
        return Vout

def plot_projection_general(dictVp, model, POV = "", ThetaDisplay = "",antiphase = True, SaveFlag = False, dir = "", Name = "", POVname="", extension = ".png", antipodal = False):
    """
    Vp: vector of parameters
    model: model adopted e.g. ST-U, ST-PST !PLEASE USE "-S" and not just "S" for symmetric models
    POV: Point Of View - location (for string arguments: colatitude) from which the neutron star is observed; if string can be:
        "I" -from the point of view of Earth-: in this case x and y are rotated compared
                                           other configuration to visualise equator horizontally and North on top
        or "P" -primary-,"S", -secondary-, SO" -secondary omission-, "SE" -secondary emission-,  "SC" -secondary cede-, or equivalently "PO", "PE", "PC", depending on the adopted model.
        A vector of three coordinates, whose squadratic sum =1, can also be used to determin the location from which the NS surface is seeing.
    ThetaDisplay: in which pole would you like to see the countours of constant colatitutde? Options are: "SP" for south pole and "NP" for north pole (where i ==0)
    antiphase: is antiphase turned on for the secondary spot?
    SaveFlag: bool: do you want to save the image
    dir: directory path where to save the image
    Name: Name to add to the file to be saved
    POVname: If POV is something outside of the aforementioned options e.g. using a vector of 3 coordinates,
             provide POV name to save file with
    extension: extension of the file to be saved
    antipodal: if you want to show antipodal configuration
    """

    Plab = 'Primary'
    Slab = 'Secondary'
    EMIlab = ' emission'
    OMIlab = ' omission'
    SUPERlab = ' superseding'
    CEDElab = ' ceding'

    cm        = plt.get_cmap('RdBu')#gist_earth; cividis
    nColors   = 7
    mycolors0 = [cm(xcol) for xcol in np.linspace(0,1 , nColors)]


    def fillVECTORS(HStag):
        DICT_VECTOR = []
        try:
            DICT_VECTOR = dictVp.keys()
        except AttributeError:
            DICT_VECTOR = dictVp.names

        if (HStag!='p') and (HStag!='s'):
            print ("ERROR: invalid Key argument (only 'p' and 's' allowed), not ",HStag)
            raise IpyExit

        Mlab = Plab if HStag=='p' else Slab
        labels = [Mlab]

        thetaA =[dictVp[HStag+REF['SC']]]

        print ("Phase",dictVp[HStag+REF['PS']] )
        phi = dictVp[HStag+REF['PS']] +0.5 if ((antiphase) and (HStag=='s')) else dictVp[HStag+REF['PS']]
        print ("Phase",phi, antiphase,HStag)

        phiA = [phi]

        zetaA =[dictVp[HStag+REF['SR']]]
        TA = [dictVp[HStag+REF['ST']]]

        if (HStag+REF['OC']) in DICT_VECTOR:
            thetaA.append(dictVp[HStag+REF['OC']])
            zetaA.append(dictVp[HStag+REF['OR']])
            phi_o = phiA[0]
            phi_o = phi_o
            phi_s = phi_o-dictVp[HStag+REF['OA']]/(2*np.pi)
            phiA = [phi_s,phi_o]
            labels = [Mlab+EMIlab, Mlab + OMIlab]

        elif (HStag+REF['CC']) in DICT_VECTOR:
            thetaA.append(dictVp[HStag+REF['CC']])
            zetaA.append(dictVp[HStag+REF['CR']])
            TA.append(dictVp[HStag+REF['CT']])
            #phi_c = phi_c+0.5 if ((antiphase) and (HStag=='s')) else phi_c
            phi_c = phiA[0]
            phi_c = phi_c+dictVp[HStag+REF['CA']]/(2*np.pi)
            phiA.append(phi_c)
            labels = [Mlab+SUPERlab, Mlab + CEDElab]

        elif (HStag+REF['OR']) in DICT_VECTOR:
            if not('C' in model):
                print ("ERROR: C is not in model name, colatitude omission not present, omission ang. radius present")
                raise IpyExit
            thetaA.append(dictVp[HStag+REF['SC']])
            zetaA.append(dictVp[HStag+REF['OR']])
            phi_s = phiA[0]
            phiA.append(phi_s)
            labels = [Mlab+EMIlab, Mlab + OMIlab]

        elif (HStag+REF['CR']) in DICT_VECTOR:
            if not('C' in model):
                print ("ERROR: C is not in model name, colatitude ceding not present, ceding ang. radius present")
                raise IpyExit

            thetaA.append(dictVp[HStag+REF['SC']])
            zetaA.append(dictVp[HStag+REF['CR']])
            phi_s = phiA[0]
            phiA.append(phi_s)
            TA.append(dictVp[HStag+REF['CT']])
            labels = [Mlab+SUPERlab, Mlab + CEDElab]
        labT_0 = " log(T/K) = %1.2f"%TA[0]
        labels[0] = labels[0] +labT_0
        if len(TA)>1:
            labT_1 = " log(T/K) = %1.2f"%TA[1]
            labels[1] =labels[1] +labT_1
        return phiA,thetaA,zetaA,TA,labels



    def SYMMETRIC(phiA,thetaA,labels):
        phiA_s = np.zeros(len(phiA))
        thetaA_s = np.zeros(len(thetaA))
        labels_s = list(labels)

        for i in range(len(phiA)):
            phiA_s[i] = phiA[i]+0.5 if (phiA[i])<0.5 else phiA[i]-0.5
            thetaA_s[i] = np.pi - thetaA[i]
            labels_s[i] = labels_s[i].replace(Plab,Slab)

        return phiA_s,thetaA_s,labels_s


    phiA_p,thetaA_p,zetaA_p,TA_p,labels_p = fillVECTORS('p')


    phiA_s = []
    thetaA_s = []
    zetaA_s = []
    TA_s = []
    labels_s =[]



    CAall =[]
    if '-S' in model:

        phiA_s,thetaA_s,labels_s = SYMMETRIC(phiA_p,thetaA_p,labels_p)

        zetaA_s = zetaA_p
        TA_s = TA_p

        if len(TA_p)==2:
            iii = 0 if TA_p[0]>TA_p[1] else 6
            CAall=[mycolors0[7-iii],mycolors0[iii],mycolors0[7-iii],mycolors0[iii]]
        else:
            if len(TA_p)==len(labels_p):
                CAall=[mycolors0[6],mycolors0[6]]
            else:
                CAall=[mycolors0[6],'black',mycolors0[6],'black']

    else:
        phiA_s,thetaA_s,zetaA_s,TA_s,labels_s = fillVECTORS('s')
        print ("phiA_s",phiA_s)
        TAall = TA_p+TA_s
        LABall = labels_p +labels_s
        ind = [1,2,4,5] if len(TAall)>2 else [0,6]

        argT = np.argsort(TAall)
        j = 0
        for i in range(len(LABall)):
            if (OMIlab in LABall[i]):
                coli = 'black'
            else:
                ind_i = ind[argT[j]]
                coli = mycolors0[ind_i]
                j = j+1
            CAall.append(coli)
            
     if antipodal:
            
            phiA_anti,thetaA_anti,labels_anti = SYMMETRIC(phiA_p,thetaA_p,labels_p)
            zetaA_anti = zetaA_p

    cosi    = dictVp['cos_inclination']
    phi_inc = 0.0

    cm        = plt.get_cmap('magma')
    nColors   = 8
    mycolors1 = [cm(xcol) for xcol in np.linspace(0,1 , nColors)]


    fig = plt.figure(figsize=(8,8))
    ax = fig.add_subplot(111)
    ax.grid(linestyle="-.", color='grey',linewidth=0.7)


    veneer(None,None, ax)

    alphaOH = 0.35 #OE=Opposite Hemisphere

    nr = 10
    r  = np.linspace(0,1,nr)

    n  = 200
    x2 = np.linspace(0,2.*np.pi,n)

    #drawing circles from pole
    VA = []
    for i in range(len(zetaA_p)):
        Vi  = [np.cos(x2)*np.sin(zetaA_p[i]),np.sin(x2)*np.sin(zetaA_p[i]),np.cos(zetaA_p[i])*np.ones(len(x2))]
        VA.append(Vi)
    for i in range(len(zetaA_s)):
        Vi  = [np.cos(x2)*np.sin(zetaA_s[i]),np.sin(x2)*np.sin(zetaA_s[i]),np.cos(zetaA_s[i])*np.ones(len(x2))]
        VA.append(Vi)



    theta_sub = 0.0
    phi_sub   = 0.0

    NO_POLE_FLAG = False


    allowedPOV = ["I"]
    LABall = labels_p +labels_s

    flag_P = False
    flag_S = False
    for l in LABall:

        if (Plab + OMIlab) in l:
            allowedPOV.append("PO")
            allowedPOV.append("PE")
            flag_P = True
        elif (Slab + OMIlab) in l:
            allowedPOV.append("SO")
            allowedPOV.append("SE")
            flag_S = True
        elif (Plab + CEDElab) in l: # l == Plab + CEDElab
            allowedPOV.append("PS")
            allowedPOV.append("PC")
            flag_P = True
        elif (Slab + CEDElab) in l:
            allowedPOV.append("SS")
            allowedPOV.append("SC")
            flag_S = True


    if not(flag_P):
        allowedPOV.append("P")
    if not(flag_S):
        allowedPOV.append("S")

    if (isinstance(POV, str)) and (POV not in allowedPOV):
        print ("ERROR: point of view not allowed for this model!")
        print ("      POSSIBILITIES are: ",allowedPOV)
        raise IpyExit

    labPOV = ""
    if POV!="":
        NO_POLE_FLAG = True
        if POV=="SO":
            if not('-S' in model):
                theta_sub = dictVp['s'+REF['OC']]
                phi_sub = dictVp['s'+REF['PS']] +0.5 if antiphase else dictVp['s'+REF['PS']]
            else:
                ii = [i for i, s in enumerate(labels_s) if OMIlab in s]
                theta_sub = thetaA_s[ii]
                phi_sub   = phiA_s[ii]
            labPOV ="SecondaryOmission"

        elif POV=="SE":
            if not('-S' in model):
                theta_sub = dictVp['s'+REF['SC']]
                phi_o = dictVp['s'+REF['PS']]+0.5 if antiphase else dictVp['s'+REF['PS']]
                phi_s = phi_o-dictVp['s'+REF['OA']]/(2*np.pi)
                phi_sub = phi_s
            else:
                ii = [i for i, s in enumerate(labels_s) if EMIlab in s]
                theta_sub = thetaA_s[ii]
                phi_sub   = phiA_s[ii]
            labPOV ="SecondaryEmission"

        elif POV=="S":
            if not('-S' in model):
                theta_sub = dictVp['s'+REF['SC']]
                phi_sub = dictVp['s'+REF['PS']]+0.5 if antiphase else dictVp['s'+REF['PS']]
            else:
                theta_sub = thetaA_s[0]
                phi_sub   = phiA_s[0]
            print ("phi_sub",phi_sub)
            labPOV ="Secondary"

        elif POV=="SS":
            if not('-S' in model):
                theta_sub = dictVp['s'+REF['SC']]
                phi_sub = dictVp['s'+REF['PS']] +0.5 if antiphase else dictVp['s'+REF['PS']]
            else:
                ii = [i for i, s in enumerate(labels_s) if SUPERlab in s]
                theta_sub = thetaA_s[ii]
                phi_sub   = phiA_s[ii]
            labPOV ="SecondarySuper"

        elif POV=="SC":
            if not('-S' in model):
                theta_sub = dictVp['s'+REF['CC']]
                phi_s = dictVp['s'+REF['PS']] +0.5 if antiphase else dictVp['s'+REF['PS']]
                phi_c = phi_s+dictVp['p'+REF['CA']]/(2*np.pi)
                phi_sub = phi_c
            else:
                ii = [i for i, s in enumerate(labels_s) if CEDElab in s]
                theta_sub = thetaA_s[ii]
                phi_sub   = phiA_s[ii]
            labPOV ="SecondaryCede"

        elif POV=="P":
            theta_sub = dictVp['p'+REF['SC']]
            phi_sub = dictVp['p'+REF['PS']]
            labPOV ="Primary"
        elif POV=="PO":
            theta_sub = dictVp['p'+REF['OC']]
            phi_sub = dictVp['p'+REF['PS']]
            labPOV ="PrimaryOmission"
        elif POV=="PE":
            theta_sub = dictVp['p'+REF['SC']]
            phi_o = dictVp['p'+REF['PS']]
            phi_s = phi_o-dictVp['p'+REF['OA']]/(2*np.pi)
            phi_sub = phi_s
            labPOV ="PrimaryEmission"
        elif POV=="PS":
            theta_sub = dictVp['p'+REF['SC']]
            phi_sub = dictVp['p'+REF['PS']]
            labPOV ="PrimarySuper"
        elif POV=="PC":
            theta_sub = dictVp['p'+REF['CC']]
            phi_s = dictVp['p'+REF['PS']]
            phi_c = phi_s+dictVp['p'+REF['CA']]/(2*np.pi)
            phi_sub = phi_c
            labPOV ="PrimaryCede"

        elif POV=="I":
            theta_sub = np.arccos(cosi)
            labPOV ="Earth"
        elif not(isinstance(POV, str)):
            a = np.matrix(POV)*np.transpose(np.matrix([0,0.,1.]))
            cosT = a[0,0]
            theta_sub = np.arccos(cosT)

            phi_sub = 0.0
            if np.abs(np.sin(theta_sub))>1e-16:
                a0 = POV[0]/np.sin(theta_sub)
                phi_sub = np.arccos(a0) if POV[1]>=0. else np.arccos(a0)+np.pi
                phi_sub = phi_sub/(2.*np.pi)

            labPOV = "VECTOR%1.2f_%1.2f_%1.2f"

    #change to add weird transformations
    def useT(T):
        return T
    def useP(P):
        return P

    #function to rotate and plot a vector V by theta t[rad] and phi p[cycles]
    signFlag    = True
    label_plotA = []
    label_plot  = ''

    def tranform_plot(t,p,V,MS,signFlag,symb,col,FlagLeg=False,L1="",L2="",cL ="gray"):
        sign       = 'x'
        label_plot = ''
        out        = transform(useT(t),useP(p),V)

        if NO_POLE_FLAG:
            out = transform(useT(-theta_sub),useP(-phi_sub),out,phi_sub)

        alpha_o = 1 if out[2]>0. else alphaOH

        if signFlag:
            sign = 'x' if out[2]>0. else symb;

        if symb and not(signFlag):
            sign = symb

        if FlagLeg and (L1 and L2):
            label_plot = L1 if out[2]>=0. else L2
            if not(label_plot in label_plotA):
                ax.plot([],[],sign, color = cL,markersize = MS, alpha = alpha_o,label = label_plot)
                label_plotA.append(label_plot)
            label_plot = ''

        if FlagLeg and (not(L1) or not(L2)):
            label_plot = L1 if L1 else L2
        ax.plot(out[1],-out[0],sign, color = col,markersize = MS,label = label_plot, alpha = alpha_o)


    def RotateList(tryout):
        if (np.shape(tryout)[1]==1):
            return np.array(tryout)
        tryout_xD    = list(np.abs(tryout[0,1:] -  tryout[0,:-1]))
        tryout_yD    = list(np.abs(tryout[1,1:] -  tryout[1,:-1]))
        tryout_xDord = np.sort(tryout_xD)
        tryout_yDord = np.sort(tryout_yD)

        max_index_x  = tryout_xD.index(tryout_xDord[-1])
        max_index_y  = tryout_yD.index(tryout_yDord[-1])
        indexR       = -1


        if (max_index_x == max_index_y):
            indexR = max_index_x
        elif (tryout_xDord[-1] > 3.*tryout_xDord[-2]):
            indexR = max_index_x
        elif (tryout_yDord[-1] > 3.*tryout_yDord[-2]):
            indexR = max_index_y

        if indexR>=0:
            indexR = indexR+1
            nE = np.shape(tryout)[1]
            Tx = list(tryout)[0][indexR:nE]
            Tx = list(Tx)
            Tx.extend(list(tryout)[0][0:indexR])
            Ty = list(tryout)[1][indexR:nE]
            Ty = list(Ty)
            Ty.extend(list(tryout)[1][0:indexR])
            tryout = np.array(tryout)
            tryout[0,:] = np.array(Tx)
            tryout[1,:] = np.array(Ty)
        return np.array(tryout)



    def transform_plot_lines(TA,PA,VA,color_plot,NO_POLE_FLAGt,style_plot = '-',lw_plot = 1, LG = ""):

        tryout = transform(useT(TA),useP(PA),VA)


        if NO_POLE_FLAGt:
            tryout = transform(useT(-theta_sub),useP(-phi_sub),tryout,phi_sub)

        tryout_p = tryout[:,tryout[2,:] >= 0.]
        tryout_n = tryout[:,tryout[2,:] <  0.]

        if (tryout_p[2,:].any()<1e-15) & (tryout_n[2,:].any()>1e-15):
            tryout_p = tryout[:,tryout[2,:]>np.max(tryout[2,:])]
            tryout_n = np.array(tryout)

        elif (np.shape(tryout_n)[1]!=0.) and (np.shape(tryout_p)[1]!=0.):
            tryout_p = RotateList(tryout_p)
            tryout_n = RotateList(tryout_n)

        ax.plot(tryout_p[1,:],-tryout_p[0,:], style_plot, color = color_plot, lw = lw_plot, label = LG)
        ax.plot(tryout_n[1,:],-tryout_n[0,:], style_plot, color = color_plot, lw = lw_plot, alpha = alphaOH)



    x= np.linspace(0,2.*np.pi,1000)
    yfact = np.linspace(0.,np.pi*0.5,10)
    #plot star edges:
    ax.plot(np.cos(x),np.sin(x),'-', color = 'gray', lw = 1, alpha = 0.5)

    CIRC_A = []
    for i in range(len(yfact)):
        CIRC_i = np.array([np.sin(yfact[i])*np.cos(x),np.sin(yfact[i])*np.sin(x),np.cos(yfact[i])*np.ones(len(x))])
        if ThetaDisplay == "NP":
            transform_plot_lines(0.0,0.0,CIRC_i,'gray',NO_POLE_FLAG,'-.',0.5)
        elif ThetaDisplay == "SP":
            CIRC_i = np.array([np.sin(yfact[i])*np.cos(x),np.sin(yfact[i])*np.sin(x),-np.cos(yfact[i])*np.ones(len(x))])
            transform_plot_lines(0.0,0.0,CIRC_i,'gray',NO_POLE_FLAG,'-.',0.5)
        else:
            transform_plot_lines(0.0,0.0,CIRC_i,'gray',False,'-.',0.5)
        CIRC_A.append(CIRC_i)

    transform_plot_lines(0.0,0.0,CIRC_A[-1],'dimgray',True,'-',2,"Equator")


    #rotating and plotting edges of hot regions
    THETAall = list(thetaA_p)+list(thetaA_s)
    PHIall = list(phiA_p)+list(phiA_s)
    LABall = list(labels_p)+list(labels_s)
    TAall = list(TA_p)+list(TA_s)


    for i in range(len(VA)):
        print ("colatitude: ",THETAall[i],", phase: ",PHIall[i]," of ",LABall[i])
        transform_plot_lines(THETAall[i],PHIall[i],VA[i],CAall[i],NO_POLE_FLAG,'-',3,LABall[i])
        tranform_plot(THETAall[i],PHIall[i],[0,0,1],10,True,'o',CAall[i], True,'Hot Region centers','Hot Region centers -Opposite Hemisphere')

    if antipodal:
        
        for i in range(len(zetaA_anti)):

            ind = np.arange(len(x2))
            x3 = x2[ind%10 ==0]
            Vi  = [np.cos(x3)*np.sin(zetaA_anti[i]),np.sin(x3)*np.sin(zetaA_anti[i]),np.cos(zetaA_anti[i])*np.ones(len(x3))]
            if i ==0:
                transform_plot_lines(thetaA_anti[i],phiA_anti[i],Vi,CAall[i],NO_POLE_FLAG,'-.',2,"Antipode of primary")
            else:
                transform_plot_lines(thetaA_anti[i],phiA_anti[i],Vi,CAall[i],NO_POLE_FLAG,'-.',2)
                
    L1 = r'$\phi=0\,[cycle],\, \theta = \pi/2\,[rad]$'
    tranform_plot(np.pi*0.5,0.0,[0,0,1],5,False,'o', mycolors1[2], True,L1, '')

    L1 = r'$\phi=0.125\,[cycle],\, \theta = \pi/2\,[rad]$'
    tranform_plot(np.pi*0.5,0.125,[0,0,1],5,False,'D', mycolors1[2], True,L1, '')



    #rotating and plotting lines defining the phi = 0 and positive phi
    MultA = [np.zeros(len(r)),r*np.sin(0.0), r*np.cos(0.0)]
    transform_plot_lines(np.pi*0.5,0.0,MultA,mycolors1[2],NO_POLE_FLAG,'--',1)
    transform_plot_lines(np.pi*0.5,0.125,MultA,mycolors1[2],NO_POLE_FLAG,'--',1)

    #rotating and plotting poles
    tranform_plot(0.,0.,[0,0,1],10,False,'*',mycolors1[6], True, "Projected North Pole","Projected North Pole  -Opposite Hemisphere",mycolors1[6])
    tranform_plot(0.,0.,[0,0,-1],25,False,'*',mycolors1[6], True, "Projected South Pole","Projected South Pole  -Opposite Hemisphere",mycolors1[6])

    ax.legend(fontsize = legsize,loc='upper left', bbox_to_anchor=(1, 1.))
    ax.set_xlim([-1.25,1.25])
    ax.set_ylim([-1.25,1.25])

    if SaveFlag:
        if dir:
            dir = dir+'/' if (list(dir)[-1]!= '/') else dir

        FigName = dir+'Projection'
        if isinstance(POV, str):
            FigName = FigName+'_from_%s'%labPOV
            FigName = FigName+'_'+Name+extension if(Name) else FigName+extension
        elif bool(POVname.strip()):
            FigName = FigName+'_from_'+POVname
            FigName = FigName+'_'+Name+extension if(Name) else FigName+extension
        else:
            raise RuntimeError('Please provide POVname to save file with.')
        plt.savefig(FigName, dpi=300, bbox_inches='tight')

    return ax
