"""
The main idea of this tool is to create a 2D image showing the neutron star surface with its hot spots.
The hemisphere facing the observer has solid lines describing the hot spot; the hemisphere opposite to the observer has dimmed colored lines compared to what is shown in the legend.
The centers of a hot spot on the hemisphere facing the observer are tagged with crosses; the ones corresponding to the hot spots on the opposite hemisphere are tagged with circles.
Each hot spot is created with one of the poles as the center and then rotated to the right location.
Further rotation is necessary to account for the point of view.
Colors: omitting regions are marked in black; emitting regions are plotted with RdBu colormap: dark blue -> light blue -> dim red -> dark red, moving from hotter to colder emitting hot spots.
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
from matplotlib.ticker import MultipleLocator, AutoLocator, AutoMinorLocator
import sys
from io import StringIO

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

REF = dict(list(zip(ABBl,REFl)))

class IpyExit(SystemExit):
    """
    Exit Exception for IPython.
    Exception temporarily redirects stderr to buffer.
    """
    def __init__(self):
        print("exiting")  # optionally print some message to stdout, too
        # ... or do other stuff before exit
        sys.stderr = StringIO()

    def __del__(self):
        sys.stderr.close()
        sys.stderr = sys.__stderr__  # restore from backup



def veneer(x, y, axes, lw=1.0, length=8):
    """ Make the plots a little more aesthetically pleasing. Author Thomas E. Riley.

    :param x: Sets an x-axis tick every ``x[0]`` value, and labels the x-axis
        tick values every ``x[1]`` value.
    :param y: Sets a y-axis tick every ``y[0]`` value, and labels the y-axis
        tick values every ``y[1]`` value.
    :param axes: Axes of the plot.
    :param lw: Tick line width in points.
    :param length: Tick length in points.
    """
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
    plt.setp(list(axes.spines.values()), linewidth=lw, color='black')


def transform(thetaR,phiR,V,phi0=0.0):
    """ Rotate vector ``V`` to the right coordinates.

    :param thetaR: Desired colatitude of the vector (could be a hot spot
        centre) [rad].
    :param phiR: Desired phase of the vector (could be a hot spot centre)
        [cycles].
    :param V: Starting vector (could be defining the hot spot circle).
    :param phi0: Phase of the point of view [cycle].

    :returns: The rotated vector describing the hot spot.
    """
    phi0 = phi0*np.pi*2.

    RA = [-np.sin(phi0),np.cos(phi0),0.]
    V = np.asarray(V, dtype=float)
    phiR = phiR*2*np.pi

    Mphi     = np.array([[np.cos(phiR),-np.sin(phiR),0.],
                            [np.sin(phiR),np.cos(phiR),0.],
                            [0.,0.,1.]]) #_z
    Mtheta = np.array([[np.cos(thetaR)+ RA[0]*RA[0]*(1.-np.cos(thetaR)),
                            RA[0]*RA[1]*(1.-np.cos(thetaR)) - RA[2]*np.sin(thetaR),
                            RA[0]*RA[2]*(1.-np.cos(thetaR)) + RA[1]*np.sin(thetaR),],
                        [RA[0]*RA[1]*(1.-np.cos(thetaR)) + RA[2]*np.sin(thetaR),
                            np.cos(thetaR)+ RA[1]*RA[1]*(1.-np.cos(thetaR)),
                            RA[1]*RA[2]*(1.-np.cos(thetaR)) - RA[0]*np.sin(thetaR),],
                        [RA[0]*RA[2]*(1.-np.cos(thetaR)) - RA[1]*np.sin(thetaR),
                            RA[1]*RA[2]*(1.-np.cos(thetaR)) + RA[0]*np.sin(thetaR),
                        np.cos(thetaR)+ RA[2]*RA[2]*(1.-np.cos(thetaR))]])

    Mtheta[np.abs(Mtheta)<1e-16] = 0.0
    M = Mphi @ Mtheta
    if V.ndim == 2 and np.shape(V)[0] == 1:
        V = np.transpose(V)
    Vout = M @ V
    Vout = np.squeeze(Vout)
    return Vout

def plot_projection_general(dictVp, model, POV = "", ThetaDisplay = "",antiphase = False, SaveFlag = False, dir = "", Name = "", POVname="", extension = ".png", antipodal = False):
    """ Print a 2D visualisation of hot spot(s) on a spherical surface.

    :param dictVp: Vector of parameters (or an X-PSI ParameterSubspace).
    :param str model: Model adopted, e.g. ``"ST-U"``, ``"ST+PST"``. Please use
        ``"-S"`` and not just ``"S"`` for symmetric models. Recognised hot spot
        names are ``"ST"``, ``"PST"``, ``"CST"``, ``"EST"``, ``"PDT"``,
        ``"CDT"`` and ``"EDT"``; for a two hot spot model, connect each hot spot
        model with ``"+"`` (first the primary, second the secondary). If the two
        hot spots have the same name, add ``"-S"`` (the secondary is antipodal
        and shares the primary's temperature and radius), ``"-U"`` (the
        properties of the secondary are completely unrelated to the primary) or
        ``"-Ua"`` (the secondary is antipodal to the primary but has an
        independent radius and temperature).
    :param POV: Point of view: the location (for string arguments, the
        colatitude) from which the neutron star is observed. If a string, it can
        be ``"I"`` (from the point of view of Earth; in this case x and y are
        rotated compared to other configurations, to show the equator
        horizontally and North on top), or ``"P"`` (primary), ``"S"``
        (secondary), ``"SO"`` (secondary omission), ``"SE"`` (secondary
        emission), ``"SC"`` (secondary cede), or equivalently ``"PO"``,
        ``"PE"``, ``"PC"``, depending on the adopted model. A vector of three
        coordinates, whose quadratic sum is 1, can also be used to set the
        location from which the surface is seen.
    :param ThetaDisplay: At which pole to show the contours of constant
        colatitude: ``"SP"`` for the south pole and ``"NP"`` for the north pole
        (where i == 0).
    :param antiphase: Whether antiphase is turned on for the secondary spot.
    :param bool SaveFlag: Whether to save the image.
    :param str dir: Directory path where to save the image.
    :param str Name: Name to add to the saved file.
    :param str POVname: If ``POV`` is outside the options above (e.g. a vector
        of 3 coordinates), the name used to save the file.
    :param str extension: Extension of the saved file.
    :param antipodal: Whether to show the antipodal configuration.

    :returns: ``ax`` (``ax = fig.add_subplot(111)``).
    """

    # Set abbreviations for plot labels
    Plab = 'Primary'
    Slab = 'Secondary'
    EMIlab = ' emission'
    OMIlab = ' omission'
    SUPERlab = ' superseding'
    CEDElab = ' ceding'

    cm        = plt.get_cmap('RdBu')#gist_earth; cividis
    nColors   = 7
    mycolors0 = [cm(xcol) for xcol in np.linspace(0,1 , nColors)]

    # Define set of parameters
    DICT_VECTOR = []
    try:
        DICT_VECTOR = list(dictVp.keys())
    except AttributeError:
        DICT_VECTOR = dictVp.names
        
    def hotspot_check(hot_name, params, HStag, asymFlag = False):
        """ Check consistency between the hot spot name and the parameters.

        :param hot_name: Hot spot name according to our naming convention,
            ``(P,E,C)ST`` / ``(P,E,C)DT``.
        :param params: List of parameters.
        :param HStag: Hot spot tag, ``"p"`` for the primary hot spot (also for
            single hot spot models) or ``"s"`` for the secondary hot spot.
        :param asymFlag: If True, allows a "partially derived" hot spot to be
            described only by temperature and radius.
        """
        ### check super
        T_name = HStag+REF['ST']
        R_name = HStag+REF['SR']#'__super_radius'
        P_name = HStag+REF['PS']#'__phase_shift'
        C_name = HStag+REF['SC']#'__super_colatitude'
        if ((not(R_name in params) or not(T_name in params))) or (not(asymFlag) and (not(P_name in params) or not(C_name in params))):
            print ("ERROR! super properties required for \'%s\' hot spot ('p' = primary; 's' = secondary and '' for single hot spot model)"%HStag)
            raise IpyExit
        elif (asymFlag and ((P_name in params) or (C_name in params))):
            print ("WARNING! there are info for a secondary hot spot that will not being used")
            
        if ('DT' in hot_name):
            ### check cede
            T_name = HStag+REF['CT']
            R_name = HStag+REF['CR']
            if ('C' in hot_name) and (not(T_name in params) or not(R_name in params)):
                print ("ERROR! Double temperature (DT) models require the definition of cede properties")
                raise IpyExit
            C_name = HStag+REF['CT']
            A_name = HStag+REF['CA']
            
            if ('C' in hot_name) and ((C_name in params) or (A_name in params)):
                print ("WARNING! there are info for a complex geometry, but they are not being used")
            if (('P' in hot_name) or ('E' in hot_name)) and (not(C_name in params) or not(C_name in params)):
                print ("ERROR! Protruding and Eccentric models require the definition of cede colatitudes and phases")
                raise IpyExit
            if not(('C' in hot_name) ^ ('E' in hot_name) ^ ('P' in hot_name)):
                print ("ERROR! Double temperature models require the setting if concentric (C, hot spot name = CDT), eccentric (E, hot spot name = EDT), protruding (P, hot spot name = PDT), current name is",hot_name)
                raise IpyExit
            if (('C' in hot_name) or ('E' in hot_name)) and dictVp[R_name]<dictVp[HStag+REF['SR']]:
                if ('C' in hot_name):
                    print ("WARNING: ceding component is completely masked by the superseding one")
                else:
                    print ("WARNING: hot spot declared ECCENTRIC, but parameters will render a PROTRUDING configuration")
            
        elif ('ST' in hot_name) and (('C' in hot_name) ^ ('E' in hot_name) ^ ('P' in hot_name)):
            R_name = HStag+REF['OR']
            if ('C' in hot_name) and not((R_name in params)):
                print ("ERROR! models including an omission component require the definition of omit properties")
                raise IpyExit
            C_name = HStag+REF['OC']
            A_name = HStag+REF['OA']
            if (('C' in hot_name) or ('E' in hot_name)) and dictVp[R_name]>dictVp[HStag+REF['SR']]:
                if ('C' in hot_name):
                    print ("WARNING: emitting component is completely masked by the omitting one")
                else:
                    print ("WARNING: hot spot declared ECCENTRIC, but parameters will render a PROTRUDING configuration")
                
            if (not ('P' in hot_name) and not('E' in hot_name)) and ((C_name in params) or (A_name in params)):
                print ("WARNING! there are info for a complex geometry, but they are not being used")
            if (('P' in hot_name) or ('E' in hot_name)) and (not(C_name in params) or not(A_name in params)):
                print ("ERROR! Protruding and Eccentric models require the definition of cede colatitudes and phases")
                raise IpyExit
            
        elif ('ST' == hot_name):
                params_red = [param for param in params if HStag+'_' in param]
                if len(params_red)>4:
                    if not( (('cos_inclination' in params_red) or ('sin_inclination' in params_red)) and (len(params_red)==5)):
                        print ("WARNING! there are info for a complex geometry, but they are not being used")
        
        else:
            print ("ERROR: hot spot name not recognised")
            raise IpyExit
                
        
    def check_model_param_names():
        """ Check compatibility between the specified model and the parameters.

        :returns:
            ``(hot_name, hot_name_s)``, the names of the primary and secondary
            hot spots (according to our naming convention).
        """
        hot_name = ''
        hot_name_s = ''
        if ('-' in model) and (not('-S' in model) and  not('-U' in model) and not('-Ua' in model)):
            print (r"ERROR: model not recognised. If model has derived quantities and does not fall into the categories mentioned below: transform your dictionary to match \"-U\" requirements. Recognised hot spot names are \"ST\", \"\PST\", \"\CST\", \"EST\", \"PDT\", \"CDT\", \"EDT\"; for two hot spot model connect each hot spot model with \"+\" (first for primary; second for secondary). If the two hot spots have the same name, add \"-S\" (if all the secondary hot spot is antipodal and its temperature and radius is the same of the primary); \"-U\" if the properties of the secondary spot is completely unrelated to the primary ones; \"-Ua\" if the secondary hot spot is antipodal to the primary but has independent radius and temperature.")
            raise IpyExit
        if (model.count('+')>1) or ('-U3' in model):
            print ("ERROR: only models with one or two hot spots are supported")
            raise IpyExit

        if ('+' in model) or ('-' in model):
            print ("YOU ARE USING A 2 HOT SPOT MODEL")
            if ('-S' in model):
                ind = model.find('-')
                hot_name = model[:ind]
                hotspot_check(hot_name,DICT_VECTOR,'p')
                params_red = [param for param in DICT_VECTOR if 's_' in param]
                if len(params_red)>0:
                    if not( (('cos_inclination' in params_red) or ('sin_inclination' in params_red)) and (len(params_red)==1)):
                        print ("WARNING! there are info for a complex geometry, but they are not being used")
            
            elif ('-Ua' in model):
                print ("WARNING: the hot spot whose colatitude and phase are derived is assumed to be the secondary; if it is not we suggest the user to redefine the dictionary and the model as an ST-U")
                ind = model.find('-')
                hot_name = model[:ind]
                hotspot_check(hot_name,DICT_VECTOR,'p')
                hotspot_check(hot_name,DICT_VECTOR,'s', asymFlag = True)
                
            elif ('-U' in model):
                ind = model.find('-')
                hot_name = model[:ind]
                hot_name_s = model[:ind]
                hotspot_check(hot_name,DICT_VECTOR,'p')
                hotspot_check(hot_name,DICT_VECTOR,'s')
            
                
            elif ('+' in model):
                ind = model.find('+')
                hot_name = model[:ind]
                hotspot_check(hot_name,DICT_VECTOR,'p')
                hot_name_s = model[ind+1:]
                hotspot_check(hot_name_s,DICT_VECTOR,'s')
                
            else:
                print ("ERROR! MODEL NOT RECOGNISED")
                raise IpyExit
                
        else:
            print ("YOU ARE USING 1 HOT SPOT MODEL")
            hot_name = model
            hotspot_check(hot_name,DICT_VECTOR,'p')

            
        return hot_name,hot_name_s
        
        
    def fillVECTORS(HStag, hot_name, flagUa = False):
        """ Fill vectors describing the geometries and temperatures of the hot
        spots, and labels for the plots.

        :param HStag: Hot spot tag, ``"p"`` for the primary hot spot (also for
            single hot spot models) or ``"s"`` for the secondary hot spot.
        :param hot_name: Name (according to our naming convention) of the hot
            spot.
        :param bool flagUa: True if the model is ``"-Ua"``.

        :returns:
            ``(phiA, thetaA, zetaA, TA, labels)``, where ``phiA`` are the phases
            of the hot spot components [cycle], ``thetaA`` their colatitudes
            [rad], ``zetaA`` their angular radii [rad], ``TA`` their log10
            temperatures [log10(T/K)], and ``labels`` the label strings for each
            component.
        """

        if (HStag!='p') and (HStag!='s'):
            print ("ERROR: invalid Key argument (only 'p' and 's' allowed), not ",HStag)
            raise IpyExit

        if HStag == 'p':
            Mlab = Plab
        elif HStag == 's':
            Mlab = Slab

        labels = [Mlab]
    
        thetaA = []
        if not(flagUa):
            thetaA =[dictVp[HStag+REF['SC']]]
        
        phiA = []
        if not(flagUa):
            phi = dictVp[HStag+REF['PS']] +0.5 if ((antiphase) and (HStag=='s')) else dictVp[HStag+REF['PS']]

            phiA = [phi]

        zetaA =[dictVp[HStag+REF['SR']]]
        TA = [dictVp[HStag+REF['ST']]]

        if (HStag+REF['OR']) in DICT_VECTOR and (('P' in hot_name) ^ ('E' in hot_name)):
            zetaA.append(dictVp[HStag+REF['OR']])
            if not(flagUa):
                thetaA.append(dictVp[HStag+REF['OC']])
                phi_o = phiA[0]
                phi_o = phi_o
                phi_s = phi_o-dictVp[HStag+REF['OA']]/(2*np.pi)
                phiA = [phi_s,phi_o]
            labels = [Mlab+EMIlab, Mlab + OMIlab]
        elif (HStag+REF['OR']) in DICT_VECTOR and ('C' in hot_name):
            zetaA.append(dictVp[HStag+REF['OR']])
            if not(flagUa):
                thetaA.append(dictVp[HStag+REF['SC']])
                phi_s = phiA[0]
                phiA = [phi_s,phi_s]
            labels = [Mlab+EMIlab, Mlab + OMIlab]

        elif (HStag+REF['CR']) in DICT_VECTOR and (('P' in hot_name) ^ ('E' in hot_name)):
            zetaA.append(dictVp[HStag+REF['CR']])
            TA.append(dictVp[HStag+REF['CT']])
            if not(flagUa):
                thetaA.append(dictVp[HStag+REF['CC']])
                phi_c = phiA[0]
                phi_c = phi_c+dictVp[HStag+REF['CA']]/(2*np.pi)
                phiA.append(phi_c)
            labels = [Mlab+SUPERlab, Mlab + CEDElab]
        elif (HStag+REF['CR']) in DICT_VECTOR and ('C' in hot_name):
            zetaA.append(dictVp[HStag+REF['CR']])
            TA.append(dictVp[HStag+REF['CT']])
            if not(flagUa):
                thetaA.append(dictVp[HStag+REF['SC']])
                phi_s = phiA[0]
                phiA.append(phi_s)
            labels = [Mlab+SUPERlab, Mlab + CEDElab]
            
        
        labT_0 = " log(T/K) = %1.2f"%TA[0]
        labels[0] = labels[0] +labT_0
        if len(TA)>1:
            labT_1 = " log(T/K) = %1.2f"%TA[1]
            labels[1] =labels[1] +labT_1
        
        ##### CHECK if VALUES ARE SENSIBLE#####
        if not(flagUa):
            for phi in phiA:
                if (phi >1. or phi < -1.):
                    print ("WARNING: phases out of range (< 1 or > 1)!!!")
                if (phi >2. or phi < -2.):
                    print ("ERROR: phases out of range (< 2 or > 2)")
                    raise IpyExit
        if not(flagUa):
            for theta in thetaA:
                if (theta > np.pi or theta < 0.):
                    print ("ERROR: colatitude out of range (< 0 or > pi)")
                    raise IpyExit
                
        for zeta in zetaA:
            if (zeta > np.pi*0.5 or zeta < 0.):
                print ("ERROR: angular radius out of range (< 0. or > 0.5pi)")
                raise IpyExit
        
        for logT in TA:
            if (logT > 15 or logT < 3):
                print ("ERROR: log10 temperature out of range (> 15 or < 3)")
                raise IpyExit
                
        return phiA,thetaA,zetaA,TA,labels



    def SYMMETRIC(phiA,thetaA,labels):
        """ Set the phases and colatitudes for the antipodal configuration of
        the secondary hot spot.

        :param phiA: Array of phases for the components of the primary hot spot
            [cycle].
        :param thetaA: Array of colatitudes for the components of the primary
            hot spot [rad].
        :param labels: Array of labels for the components of the primary hot
            spot.

        :returns:
            ``(phiA_s, thetaA_s, labels_s)``, the phases [cycle], colatitudes
            [rad] and labels for the components of the secondary (derived) hot
            spot.
        """
        phiA_s = np.zeros(len(phiA))
        thetaA_s = np.zeros(len(thetaA))
        labels_s = list(labels)

        for i in range(len(phiA)):
            phiA_s[i] = phiA[i]+0.5 if (phiA[i])<0.5 else phiA[i]-0.5
            thetaA_s[i] = np.pi - thetaA[i]
            labels_s[i] = labels_s[i].replace(Plab,Slab)

        return phiA_s,thetaA_s,labels_s

    # run model check
    HSname_p, HSname_s =  check_model_param_names()
    
    # Defining phases, colatitudes, angular radii, temperatures and labels for primary hot spot
    phiA_p,thetaA_p,zetaA_p,TA_p,labels_p = fillVECTORS('p',HSname_p)

    # Defining phases, colatitudes, angular radii, temperatures and labels for secondary hot spot (if any)
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
            iii = 1 if TA_p[0]>TA_p[1] else nColors-1
            CAall=[mycolors0[nColors-iii],mycolors0[iii],mycolors0[nColors-iii],mycolors0[iii]]
        else:
            if len(TA_p)==len(labels_p):
                CAall=[mycolors0[nColors-1],mycolors0[nColors-1]]
            else:
                CAall=[mycolors0[nColors-1],'black',mycolors0[nColors-1],'black']

    elif (('+' in model) or ('-' in model)):
        if ('-Ua' in model):
            phiA_s,thetaA_s,_ = SYMMETRIC(phiA_p,thetaA_p,labels_p)
            _,_,zetaA_s,TA_s,labels_s = fillVECTORS('s',HSname_s, flagUa = True)
        else:
            phiA_s,thetaA_s,zetaA_s,TA_s,labels_s = fillVECTORS('s',HSname_s)
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
    
    else:
        zetaA_s =[]
        
        if len(TA_p)==2:
            iii = 1 if TA_p[0]>TA_p[1] else nColors-1
            CAall=[mycolors0[nColors-iii],mycolors0[iii],mycolors0[nColors-iii],mycolors0[iii]]
        else:
            if len(TA_p)==len(labels_p):
                CAall=[mycolors0[nColors-1],mycolors0[nColors-1]]
            else:
                CAall=[mycolors0[nColors-1],'black',mycolors0[nColors-1],'black']
    
    if antipodal:
            
        phiA_anti,thetaA_anti,labels_anti = SYMMETRIC(phiA_p,thetaA_p,labels_p)
        zetaA_anti = zetaA_p

    try:
        cosi    = dictVp['cos_inclination']
    except KeyError:
        cosi    = np.sqrt( 1. - dictVp['sin_inclination']**2 )
    phi_inc = 0.0

    cm        = plt.get_cmap('magma')
    nColors   = 8
    mycolors1 = [cm(xcol) for xcol in np.linspace(0,1 , nColors)]


    alphaOH = 0.35 #OH=Opposite Hemisphere

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

    NO_POLE_FLAG = False # False when origin of 2D projection is the reference point
    
    allowedPOV = ["I"]
    if not('cos_inclination' in DICT_VECTOR) and not('sin_inclination' in DICT_VECTOR):
        print ("WARNING: view from Earth is not allowed")
        allowedPOV = []
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
        elif (Plab + CEDElab) in l:
            allowedPOV.append("PS")
            allowedPOV.append("PC")
            flag_P = True
        elif (Slab + CEDElab) in l:
            allowedPOV.append("SS")
            allowedPOV.append("SC")
            flag_S = True


    if not(flag_P):
        allowedPOV.append("P")
    if not(flag_S) and ("+" in model):
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
                theta_sub = 0.
                phi_sub = dictVp['s'+REF['PS']] +0.5 if antiphase else dictVp['s'+REF['PS']]
                if 'C' in HSname_s:
                    theta_sub = dictVp['s'+REF['SC']]
                else:
                    theta_sub = dictVp['s'+REF['OC']]
            else:
                ii = [i for i, s in enumerate(labels_s) if OMIlab in s]
                theta_sub = thetaA_s[ii[0]]
                phi_sub   = phiA_s[ii[0]]
            labPOV ="SecondaryOmission"

        elif POV=="SE":
            if not('-S' in model):
                theta_sub = dictVp['s'+REF['SC']]
                phi_o = dictVp['s'+REF['PS']]+0.5 if antiphase else dictVp['s'+REF['PS']]
                phi_s = 0.
                phi_sub = 0.
                if 'C' in HSname_s:
                    phi_sub = phi_o
                else:
                    phi_s = phi_o-dictVp['s'+REF['OA']]/(2*np.pi)
                    phi_sub = phi_s
            else:
                ii = [i for i, s in enumerate(labels_s) if EMIlab in s]
                theta_sub = thetaA_s[ii[0]]
                phi_sub   = phiA_s[ii[0]]
            labPOV ="SecondaryEmission"

        elif POV=="S":
            if not(('-S' in model) or ('-Ua' in model)):
                theta_sub = dictVp['s'+REF['SC']]
                phi_sub = dictVp['s'+REF['PS']]+0.5 if antiphase else dictVp['s'+REF['PS']]
            else:
                theta_sub = thetaA_s[0]
                phi_sub   = phiA_s[0]
            labPOV ="Secondary"

        elif POV=="SS":
            if not('-S' in model):
                theta_sub = dictVp['s'+REF['SC']]
                phi_sub = dictVp['s'+REF['PS']] +0.5 if antiphase else dictVp['s'+REF['PS']]
            else:
                ii = [i for i, s in enumerate(labels_s) if SUPERlab in s]
                theta_sub = thetaA_s[ii[0]]
                phi_sub   = phiA_s[ii[0]]
            labPOV ="SecondarySuper"

        elif POV=="SC":
            if not('-S' in model):
                theta_sub = 0.
                phi_s = dictVp['s'+REF['PS']] +0.5 if antiphase else dictVp['s'+REF['PS']]
                phi_c = 0.
                if 'C' in HSname_s:
                    theta_sub = dictVp['s'+REF['SC']]
                    phi_c = phi_s
                else:
                    theta_sub = dictVp['s'+REF['CC']]
                    phi_c = phi_s+dictVp['s'+REF['CA']]/(2*np.pi)
                phi_sub = phi_c
            else:
                ii = [i for i, s in enumerate(labels_s) if CEDElab in s]
                theta_sub = thetaA_s[ii[0]]
                phi_sub   = phiA_s[ii[0]]
            labPOV ="SecondaryCede"

        elif POV=="P":
            theta_sub = dictVp['p'+REF['SC']]
            phi_sub = dictVp['p'+REF['PS']]
            labPOV ="Primary"
        elif POV=="PO":
            theta_sub = dictVp['p'+REF['OC']]
            phi_sub = dictVp['p'+REF['PS']]
            if 'C' in HSname_p:
                theta_sub = dictVp['p' + REF['SC']]
                phi_sub = dictVp['p'+REF['PS']]
            labPOV ="PrimaryOmission"
        elif POV=="PE":
            theta_sub = dictVp['p'+REF['SC']]
            phi_o = dictVp['p'+REF['PS']]
            phi_s = 0.
            if 'C' in HSname_p:
                phi_s = dictVp['p'+REF['PS']]
            else:
                phi_s = phi_o-dictVp['p'+REF['OA']]/(2*np.pi)
            phi_sub = phi_s
            labPOV ="PrimaryEmission"
        elif POV=="PS":
            theta_sub = dictVp['p'+REF['SC']]
            phi_sub = dictVp['p'+REF['PS']]
            labPOV ="PrimarySuper"
        elif POV=="PC":
            theta_sub = 0.
            phi_s = dictVp['p'+REF['PS']]
            phi_c = 0.
            if 'C' in HSname_p:
                theta_sub = dictVp['p'+REF['SC']]
                phi_c = phi_s
            else:
                theta_sub = dictVp['p'+REF['CC']]
                phi_c = phi_s+dictVp['p'+REF['CA']]/(2*np.pi)
            phi_sub = phi_c
            labPOV ="PrimaryCede"

        elif POV=="I":
            theta_sub = np.arccos(cosi)
            labPOV ="Earth"
        elif not(isinstance(POV, str)):
            cosT = np.dot(np.asarray(POV, dtype=float), [0., 0., 1.])
            theta_sub = np.arccos(cosT)

            phi_sub = 0.0
            if np.abs(np.sin(theta_sub))>1e-16:
                a0 = np.clip(POV[0]/np.sin(theta_sub), -1., 1.)
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

    def transform_plot(t,p,V,MS,signFlag,symb,col,FlagLeg=False,L1="",L2="",cL ="gray"):
        """ Plot a specific point on the sphere (mock neutron star) surface.

        :param t: Target colatitude of the point of interest.
        :param p: Target phase of the point of interest.
        :param V: Initial origin (usually the North pole, rarely the South
            pole); shape ``[x, y, z]``.
        :param MS: Marker size used to plot the point of interest.
        :param signFlag: True to plot hot spot centres, False otherwise.
        :param symb: Symbol plotted (unless ``signFlag`` is True and the hot
            spot is visible, in which case the marker is set to ``'x'``).
        :param col: Colour of the marker used to plot the point of interest.
        :param FlagLeg: True to add the legend entry for the marker of interest
            to the plot.
        :param L1: Legend label if the point of interest is on the hemisphere
            facing the observer.
        :param L2: Legend label if the point of interest is on the hemisphere
            opposite to the observer.
        :param cL: Colour of the hot spot centre marker in the legend (only used
            for this type of point, when different legends are used depending on
            whether the hemisphere is facing or opposite to the observer).
        """
        sign       = 'x'
        label_plot = ''
        # unless users modiefied it, useT(x) and useP(x) return x
        # tranform function: rotate vector V (starting origin of the point of interest) to coordinate (t,p) respectively (colatitude, phase)
        out        = transform(useT(t),useP(p),V)

        if NO_POLE_FLAG:
            out = transform(useT(-theta_sub),useP(-phi_sub),out,phi_sub) # changes location of point according to the desired point of view

        alpha_o = 1 if out[2]>0. else alphaOH #out[2] > 0. => point plotted on the observer facing hemisphere

        if signFlag:
            sign = 'x' if out[2]>0. else symb; #out[2] > 0. => point plotted on the observer facing hemisphere
            if sign == symb:
                MS = 3
        if symb and not(signFlag):
            sign = symb
        
        # add legend entry for point, if not already there and if there are different legend depening on the point location (emisphere facing or opposite to the obsever)
        if FlagLeg and (L1 and L2):
            label_plot = L1 if out[2]>=0. else L2
            if not(label_plot in label_plotA):
                ax.plot([],[],sign, color = cL,markersize = MS, alpha = alpha_o,label = label_plot)
                label_plotA.append(label_plot)
            label_plot = ''
        # add legend for the point of interest
        if FlagLeg and (not(L1) or not(L2)):
            label_plot = L1 if L1 else L2
        ax.plot(out[1],-out[0],sign, color = col,markersize = MS,label = label_plot, alpha = alpha_o)


    def RotateList(tryout):
        """ Modify the input vector to avoid plotting lines joining the wrong
        points on the circles.

        In practice it determines the optimal starting point of the vector so
        that the lines of the hot spot are correctly plotted.

        :param tryout: Input list.
        """
        # Check if the input list has only one column
        if (np.shape(tryout)[1]==1):
            return np.array(tryout)
            
        # Calculate the differences between consecutive x and y values
        tryout_xD    = list(np.abs(tryout[0,1:] -  tryout[0,:-1]))
        tryout_yD    = list(np.abs(tryout[1,1:] -  tryout[1,:-1]))
        
        # Sort the differences
        tryout_xDord = np.sort(tryout_xD)
        tryout_yDord = np.sort(tryout_yD)
        
        # Get the index of the largest difference in x and y
        max_index_x  = tryout_xD.index(tryout_xDord[-1])
        max_index_y  = tryout_yD.index(tryout_yDord[-1])
        indexR       = -1

        # Determine the optimal starting index of vector of interest
        if (max_index_x == max_index_y):
            indexR = max_index_x
        elif (tryout_xDord[-1] > 3.*tryout_xDord[-2]):
            indexR = max_index_x
        elif (tryout_yDord[-1] > 3.*tryout_yDord[-2]):
            indexR = max_index_y

        # defining optimal starting point in vector to avoid weird lines appearing connecting non logically consecutive (but consecutive in list) points the hot spot lines
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
        """ Plot the lines defining the hot spots.

        :param TA: Array of colatitudes of the points describing the hot spots.
        :param PA: Array of phases of the points.
        :param VA: Initial vector.
        :param color_plot: Colour of the line.
        :param NO_POLE_FLAGt: Flag, True for a specified point of view.
        :param style_plot: Line style of the hot spot to be plotted.
        :param lw_plot: Line width of the spot to be plotted.
        :param LG: Legend label for the hot spot.
        """
        tryout = transform(useT(TA),useP(PA),VA)


        if NO_POLE_FLAGt:
            tryout = transform(useT(-theta_sub),useP(-phi_sub),tryout,phi_sub)
        
        tryout = np.around(tryout,10)
        tryout_p = tryout[:,tryout[2,:] >= 0.] # find part of the hot spot contour that is located in the hemisphere facing the observer
        tryout_n = tryout[:,tryout[2,:] <  0.] # find part of the hot spot contour that is located in the hemisphere opposite to the observer
        
        # Corrects _p or _n vector in case of very low assigned values present in the array
        if (tryout_p[2,:].any()<1e-15) & (tryout_n[2,:].any()>1e-15):
            tryout_p = tryout[:,tryout[2,:]>np.max(tryout[2,:])]
            tryout_n = np.array(tryout)
        elif (np.shape(tryout_n)[1]!=0.) and (np.shape(tryout_p)[1]!=0.):
            tryout_p = RotateList(tryout_p)
            tryout_n = RotateList(tryout_n)

        # plot lines descring hot spot on hemisphere facing the observer
        ax.plot(tryout_p[1,:],-tryout_p[0,:], style_plot, color = color_plot, lw = lw_plot, label = LG)
        # plot lines descring hot spot on hemisphere opposite to the observer (dimmed lines)
        ax.plot(tryout_n[1,:],-tryout_n[0,:], style_plot, color = color_plot, lw = lw_plot, alpha = alphaOH)

    
    fig = plt.figure(figsize=(8,8))
    ax = fig.add_subplot(111)
    veneer(None,None, ax)

    # coordinates to span any cycle:
    x= np.linspace(0,2.*np.pi,1000)
    yfact = np.linspace(0.,np.pi*0.5,10)
    
    #plot star edges:
    ax.plot(np.cos(x),np.sin(x),'-', color = 'gray', lw = 1, alpha = 0.5)

    # plot circle lines marking different colatitudes from Pole or origine of the 2D plot
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

    # plot equator
    transform_plot_lines(0.0,0.0,CIRC_A[-1],'dimgray',True,'-',2,"Equator")


    #rotating and plotting edges of hot regions
    THETAall = list(thetaA_p)+list(thetaA_s)
    PHIall = list(phiA_p)+list(phiA_s)
    LABall = list(labels_p)+list(labels_s)
    TAall = list(TA_p)+list(TA_s)


    for i in range(len(VA)):
        transform_plot_lines(THETAall[i],PHIall[i],VA[i],CAall[i],NO_POLE_FLAG,'-',3,LABall[i])
        transform_plot(THETAall[i],PHIall[i],[0,0,1],10,True,'o',CAall[i], True,'Hot Region centers','Hot Region centers -Opposite Hemisphere')

    if antipodal:
        
        for i in range(len(zetaA_anti)):

            ind = np.arange(len(x2))
            Vi  = [np.cos(x2)*np.sin(zetaA_anti[i]),np.sin(x2)*np.sin(zetaA_anti[i]),np.cos(zetaA_anti[i])*np.ones(len(x2))]
            if i ==0:
                transform_plot_lines(thetaA_anti[i],phiA_anti[i],Vi,CAall[i],NO_POLE_FLAG,'-.',2,"Antipode of primary")
            else:
                transform_plot_lines(thetaA_anti[i],phiA_anti[i],Vi,CAall[i],NO_POLE_FLAG,'-.',2)
                
    L1 = r'$\phi=0\,[\mathrm{cycle}],\, \theta = \pi/2\,[\mathrm{rad}]$'
    transform_plot(np.pi*0.5,0.0,[0,0,1],5,False,'o', mycolors1[2], True,L1, '')

    L1 = r'$\phi=0.125\,[\mathrm{cycle}],\, \theta = \pi/2\,[\mathrm{rad}]$'
    transform_plot(np.pi*0.5,0.125,[0,0,1],5,False,'D', mycolors1[2], True,L1, '')



    #rotating and plotting lines defining the phi = 0 and positive phi
    MultA = [np.zeros(len(r)),r*np.sin(0.0), r*np.cos(0.0)]
    transform_plot_lines(np.pi*0.5,0.0,MultA,mycolors1[2],NO_POLE_FLAG,'--',1)
    transform_plot_lines(np.pi*0.5,0.125,MultA,mycolors1[2],NO_POLE_FLAG,'--',1)

    #rotating and plotting poles
    transform_plot(0.,0.,[0,0,1],10,False,'*',mycolors1[6], True, "Projected North Pole","Projected North Pole  -Opposite Hemisphere",mycolors1[6])
    transform_plot(0.,0.,[0,0,-1],25,False,'*',mycolors1[6], True, "Projected South Pole","Projected South Pole  -Opposite Hemisphere",mycolors1[6])

    ax.legend(fontsize = legsize,loc='upper left', bbox_to_anchor=(1, 1.))
    ax.set_xlim([-1.25,1.25])
    ax.set_ylim([-1.25,1.25])
    
    plt.axis('off')
    extent = ax.get_window_extent().transformed(fig.dpi_scale_trans.inverted())

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
