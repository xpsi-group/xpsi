cdef size_t atmos_extension = 1

from libc.math cimport log10, pow
from libc.stdio cimport printf

#Blackbody
from xpsi.surface_radiation_field.hot_BB cimport (init_hot_BB,
                                                     free_hot_BB,
                                                     eval_hot_BB,
                                                     eval_hot_norm_BB)
#4D-Numerical
# from xpsi.surface_radiation_field.hot_Num4D cimport (init_hot_Num4D,
#                                                      free_hot_Num4D,
#                                                      eval_hot_Num4D,
#                                                      eval_hot_norm_Num4D)

#Blackbody-burst (with beaming and polarization)
from xpsi.surface_radiation_field.hot_BB_burst cimport (init_hot_BB_burst,
                                                     free_hot_BB_burst,
                                                     eval_hot_BB_burst_I,
                                                     eval_hot_BB_burst_Q,
                                                     eval_hot_norm_BB_burst)
#2D-Numerical (with polarization)
from xpsi.surface_radiation_field.hot_Num2D cimport (init_hot_Num2D,
                                                     free_hot_Num2D,
                                                     eval_hot_Num2D_I,
                                                     eval_hot_Num2D_Q,
                                                     eval_hot_norm_Num2D)


from xpsi.surface_radiation_field.hot_Num2D_split cimport (init_hot_2D,
                                               eval_hot_2D_I,
                                               eval_hot_2D_norm,
                                               free_hot_2D)

#User-defined atmosphere extension (Blackbody by default)
from xpsi.surface_radiation_field.hot_user cimport (init_hot_user,
                                                     free_hot_user,
                                                     eval_hot_user_I,
                                                     eval_hot_user_Q,
                                                     eval_hot_norm_user)

from xpsi.surface_radiation_field.hot_Num5D_split cimport (init_hot_Num5D,
                                               eval_hot_Num5D_I,
                                               eval_hot_norm_Num5D,
                                               free_hot_Num5D,
                                               produce_2D_data_Num5D,
                                               make_atmosphere_2D_Num5D)

from xpsi.surface_radiation_field.hot_Num4D_split cimport (init_hot_Num4D,
                                               eval_hot_Num4D,
                                               eval_hot_norm_Num4D,
                                               free_hot_Num4D,
                                               produce_2D_data_Num4D,
                                               make_atmosphere_2D_Num4D)



from xpsi.global_imports import _keV, _k_B, _h_keV
cdef double k_B = _k_B
cdef double keV = _keV
cdef double h_keV = _h_keV
cdef double k_B_over_keV = k_B / keV

#----------------------------------------------------------------------->>>
cdef void* init_hot(size_t numThreads, const _preloaded *const preloaded, size_t atm_ext) nogil:
    global atmos_extension
    atmos_extension=atm_ext
    if atmos_extension == 1:
        return init_hot_BB(numThreads, preloaded)
    elif atmos_extension == 2:
        return init_hot_Num4D(numThreads, preloaded)
    elif atmos_extension == 3:
        return init_hot_BB_burst(numThreads, preloaded)
    elif atmos_extension == 4:
        return init_hot_Num2D(numThreads, preloaded)
    elif atmos_extension == 5:
        return init_hot_user(numThreads, preloaded)
    elif atmos_extension == 6:
        return init_hot_Num5D(numThreads, preloaded)
    else:
        printf("WARNING: Wrong atmosphere extension provided for hot region(s)."
               "Defaulting to Blackbody (atm_ext=BB).\n")
        return init_hot_BB(numThreads, preloaded)

cdef int free_hot(size_t numThreads, void *const data) nogil:
    if atmos_extension == 1:
        return free_hot_BB(numThreads, data)
    elif atmos_extension == 2:
        return free_hot_Num4D(numThreads, data)
    elif atmos_extension == 3:
        return free_hot_BB_burst(numThreads, data)
    elif atmos_extension == 4:
        return free_hot_Num2D(numThreads, data)
    elif atmos_extension == 5:
        return free_hot_user(numThreads, data)
    elif atmos_extension == 6:
        return free_hot_Num5D(numThreads, data)
    else:
        printf("WARNING: Wrong atmosphere extension provided for hot region(s)."
               "Defaulting to Blackbody (atm_ext=BB).\n")
        return free_hot_BB(numThreads, data)

cdef double eval_hot_I(size_t THREAD,
                     double E,
                     double mu,
                     const double *const VEC,
                     void *const data,
                     size_t beam_opt) nogil:

    cdef:
        double I_hot_beam=0.0
        double I_hot=0.0
        double I_hot_imu=0.0
        double I_fbeam_imu=0.0
        double I_nom=0.0
        double I_denom=0.0
        double I_denom_test=0.0
        double mu_imu=0.0
        double dmu=0.0
        double anorm=0.0
        double VEC_red[2]
        double VEC_special[3]
    cdef size_t imu

    VEC_red[0] = VEC[0]
    VEC_red[1] = VEC[1]

    cdef double E_eff = k_B_over_keV * pow(10.0, VEC[0])

    if atmos_extension == 1:
        I_hot = eval_hot_BB(THREAD,E,mu,VEC_red,data)
    elif atmos_extension == 2:
        E_dataunits = log10(E / E_eff)
        I_dataunits = eval_hot_Num4D(THREAD,E_dataunits,mu,VEC_red,data)
        I_hot = I_dataunits * pow(10.0, 3.0 * VEC[0])
    elif atmos_extension == 3:
        I_hot = eval_hot_BB_burst_I(THREAD,E,mu,VEC_red,data)
    elif atmos_extension == 4:
        E_dataunits = log10(E)
        I_hot = eval_hot_Num2D_I(THREAD,E_dataunits,mu,VEC_red,data)
    elif atmos_extension == 5:
        I_hot = eval_hot_user_I(THREAD,E,mu,VEC_red,data)
    elif atmos_extension == 6:
        VEC_special[0] = VEC[0]
        VEC_special[1] = VEC[1]
        VEC_special[2] = VEC[2]
        E_dataunits=E*0.001956951 #kev to electron rest energy conversion
        I_hot = eval_hot_Num5D_I(THREAD,E_dataunits,mu,VEC_special,data)
    else:
        printf("WARNING: Wrong atmosphere extension provided for hot region(s)."
               "Defaulting to Blackbody (atm_ext=BB).\n")
        I_hot = eval_hot_BB(THREAD,E,mu,VEC_red,data)

    if beam_opt==0:
        return I_hot

    abb = VEC[2]
    bbb = VEC[3]
    cbb = VEC[4]
    dbb = VEC[5]
    nimu = VEC[6]

    if beam_opt==1:
        I_hot_beam = (1.0+abb*((E)**cbb)*mu+bbb*((E)**dbb)*mu**2)*I_hot
    if beam_opt==2:
        anorm = 0.5/(0.5+(1.0/3.0)*abb*E**cbb+(1.0/4.0)*bbb*E**dbb)
        I_hot_beam = anorm*(1.0+abb*((E)**cbb)*mu+bbb*((E)**dbb)*mu**2)*I_hot
    if beam_opt==3:
        #integral using trapezoidal rule (consider changing to Gauss quadrature)
        for imu in range(<size_t> nimu):
            mu_imu = mu_imu+(1.0/nimu)
            if imu==0 or imu==nimu-1:
                    dmu = (0.5/nimu)
            else:
                    dmu = (1.0/nimu)
            if atmos_extension == 1:
                I_hot_imu = eval_hot_BB(THREAD,E,mu_imu,VEC_red,data)
            elif atmos_extension == 2:
                I_hot_imu = eval_hot_Num4D(THREAD,E,mu_imu,VEC_red,data)
            elif atmos_extension == 3:
                I_hot_imu = eval_hot_user_I(THREAD,E,mu_imu,VEC_red,data)
            else:
                I_hot_imu = eval_hot_BB(THREAD,E,mu_imu,VEC_red,data)

            I_fbeam_imu =  (1.0+abb*((E)**cbb)*mu_imu+bbb*((E)**dbb)*mu_imu**2)
            I_denom = I_denom + mu_imu*I_fbeam_imu*I_hot_imu*dmu
            I_nom = I_nom + mu_imu*I_hot_imu*dmu
            #I_denom_test = I_denom_test + mu_imu*dmu
            #printf("%d,%.8e,%.8e,%.8e,%.8e,%.8e\n",imu,mu_imu,dmu,I_denom,I_denom_test,I_denom_test*I_nsx_imu)
        if I_denom == 0.0:
            I_hot_beam = 0.0
        else:
            I_hot_beam = (I_nom/I_denom)*(1.0+abb*((E)**cbb)*mu+bbb*((E)**dbb)*mu**2)*I_hot

    #printf("Normalization constant: %.8e, %.8e, %.8e\n",I_hot, I_nom, I_denom)
    if I_hot_beam < 0.0:
        return 0.0
    return I_hot_beam

cdef double eval_hot_Q(size_t THREAD,
                     double E,
                     double mu,
                     const double *const VEC,
                     void *const data,
                     size_t beam_opt) nogil:
    # Arguments:
    # E = photon energy in keV
    # mu = cosine of ray zenith angle (i.e., angle to surface normal)
    # VEC = variables such as temperature, effective gravity, ...
    # data = numerical model data required for intensity evaluation

    cdef:
        double Q_hot=0.0

    #Q_hot = eval_hot_burst_Q(THREAD,E,mu,VEC_red,data)
    #cdef double I_E
    #I_E = eval_hot_I(THREAD,E,mu,VEC,data,beam_opt)
    #cdef double PD = 0.1171*(mu - 1.)/(1. + 3.582*mu)
    #return PD*I_E

    if atmos_extension == 1:
        Q_hot = 0.0
    elif atmos_extension == 2:
        printf("Warning: Polarimetry not implemented for this atmosphere extension!\n")
        Q_hot = 0.0
    elif atmos_extension == 3:
        Q_hot = eval_hot_BB_burst_Q(THREAD,E,mu,VEC,data)
    elif atmos_extension == 4:
        Q_hot = eval_hot_Num2D_Q(THREAD,E,mu,VEC,data)
    else:
        Q_hot = eval_hot_user_Q(THREAD,E,mu,VEC,data)
    return Q_hot


cdef double eval_hot_norm() nogil:
    if atmos_extension == 1:
        return eval_hot_norm_BB()
    elif atmos_extension == 2:
        return eval_hot_norm_Num4D()
    elif atmos_extension == 3:
        return eval_hot_norm_BB_burst()
    elif atmos_extension == 4:
        return eval_hot_norm_Num2D()
    elif atmos_extension == 5:
        return eval_hot_norm_user()
    elif atmos_extension == 6:
        return eval_hot_norm_Num5D()
    else:
        printf("WARNING: Wrong atmosphere extension provided for hot region(s)."
               "Defaulting to Blackbody (atm_ext=BB).\n")
        return eval_hot_norm_BB()


cdef double* produce_2D_data(size_t THREAD, const double *const VEC, void *const data) nogil:
    if atmos_extension == 6:
        return produce_2D_data_Num5D(THREAD,VEC,data)
    elif atmos_extension == 2:
        return produce_2D_data_Num4D(THREAD,VEC,data)
    else:
        printf("Error: atmosphere extension must either be 6 or 2 when using split atmospheres")

cdef object make_atmosphere_2D(double *I_data, const double *const VEC, void *const data):
    if atmos_extension == 6:
        return make_atmosphere_2D_Num5D(I_data,data)
    elif atmos_extension == 2:
        return make_atmosphere_2D_Num4D(I_data,VEC,data)
    else:
        printf("Error: atmosphere extension must either be 6 or 2 when using split atmospheres")
        
        
cdef void* init_hot_2D_split(size_t numThreads, const _preloaded *const preloaded) nogil:
    return init_hot_2D(numThreads, preloaded)
    
cdef int free_hot_2D_split(size_t numThreads, void *const data) nogil:
    return free_hot_2D(numThreads, data)

cdef double eval_hot_2D_split(size_t THREAD,
                     double E,
                     double mu,
                     const double *const VEC,
                     void *const data) nogil:
    cdef double E_eff = k_B_over_keV * pow(10.0, VEC[0])
    if atmos_extension == 2:
        E_dataunits = log10(E / E_eff)
        I_dataunits = eval_hot_2D_I(THREAD,E_dataunits,mu,data)
        I_hot = I_dataunits * pow(10.0, 3.0 * VEC[0])
    elif atmos_extension == 6:
        E_dataunits=E*0.001956951 #kev to electron rest energy conversion
        I_hot = eval_hot_2D_I(THREAD,E_dataunits,mu,data)
        # printf("I_hot = %.8e\n", I_hot)
    else:
        printf("Error: atmosphere extension must either be 6 or 2 when using split atmospheres")
    return I_hot
        
cdef double eval_hot_2D_norm_split() nogil:
    return eval_hot_2D_norm()
    
        