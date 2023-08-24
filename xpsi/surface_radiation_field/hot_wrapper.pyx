cdef size_t atmos_extension = 1

from libc.stdio cimport printf

#Blackbody
from xpsi.surface_radiation_field.hot_BB cimport (init_hot_BB,
                                                     free_hot_BB,
                                                     eval_hot_BB,
                                                     eval_hot_norm_BB)
#4D-Numerical
from xpsi.surface_radiation_field.hot_Num4D cimport (init_hot_Num4D,
                                                     free_hot_Num4D,
                                                     eval_hot_Num4D,
                                                     eval_hot_norm_Num4D)

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


#User-defined atmosphere extension (Blackbody by default)
from xpsi.surface_radiation_field.hot_user cimport (init_hot_user,
                                                     free_hot_user,
                                                     eval_hot_user_I,
                                                     eval_hot_user_Q,
                                                     eval_hot_norm_user)

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
    else:
        return init_hot_user(numThreads, preloaded)

cdef int free_hot(size_t numThreads, void *const data) nogil:
    if atmos_extension == 1:
        return free_hot_BB(numThreads, data)
    elif atmos_extension == 2:
        return free_hot_Num4D(numThreads, data)
    elif atmos_extension == 3:
        return free_hot_BB_burst(numThreads, data)
    elif atmos_extension == 4:
        return free_hot_Num2D(numThreads, data)
    else:
        return free_hot_user(numThreads, data)

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
    cdef size_t imu

    VEC_red[0] = VEC[0]
    VEC_red[1] = VEC[1]

    if atmos_extension == 1:
        I_hot = eval_hot_BB(THREAD,E,mu,VEC_red,data)
    elif atmos_extension == 2:
        I_hot = eval_hot_Num4D(THREAD,E,mu,VEC_red,data)
    elif atmos_extension == 3:
        I_hot = eval_hot_BB_burst_I(THREAD,E,mu,VEC_red,data)
    elif atmos_extension == 4:
        I_hot = eval_hot_Num2D_I(THREAD,E,mu,VEC_red,data)
    else:
        I_hot = eval_hot_user_I(THREAD,E,mu,VEC_red,data)

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
            else:
                I_hot_imu = eval_hot_user_I(THREAD,E,mu_imu,VEC_red,data)

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
        printf("Warning: Polarimetry not implemented for this atmosphere extension!")
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
    else:
        return eval_hot_norm_user()
