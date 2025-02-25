cdef double F1(double r, double r_s) nogil
    
cdef double F2(double r, double r_s) nogil
    
cdef double dF1dr(double r, double r_s) nogil

cdef double dF2dr(double r, double r_s) nogil
    
cdef long double F1_l(long double r, long double r_s) nogil

cdef long double F2_l(long double r, long double r_s) nogil

cdef long double dF1dr_l(long double r, long double r_s) nogil

cdef long double dF2dr_l(long double r, long double r_s) nogil
    
cdef double g_11(double r,
                 double theta,
                 double r_s,
                 double a,
                 double Sigma,
                 double Delta,
                 double kappa,
                 double sin_theta,
                 double func_theta,
                 double F1_r,
                 double func_r) nogil
    
cdef double g_14(double r,
                 double theta,
                 double r_s,
                 double a,
                 double Sigma,
                 double Delta) nogil
    
cdef double g_22(double r,
                 double theta,
                 double Sigma,
                 double Delta,
                 double kappa,
                 double func_theta,
                 double F1_r,
                 double func_r) nogil
    
cdef double g_33(double r, 
                 double theta,
                 double Sigma,
                 double kappa,
                 double func_theta,
                 double F2_r) nogil
    
cdef double g_44(double r,
                 double theta,
                 double a,
                 double Sigma,
                 double Delta,
                 double kappa,
                 double sin_theta,
                 double func_theta,
                 double F2_r) nogil
    
cdef double det_g(double r,
                  double theta,
                  double r_s,
                  double a,
                  double Sigma,
                  double Delta,
                  double kappa,
                  double func_theta,
                  double F1_r,
                  double F2_r,
                  double func_r,
                  double sin_theta) nogil
    
# Overload/templating in Cython?
cdef long double g_11_l(long double r,
                        long double theta,
                        long double r_s,
                        long double a,
                        long double Sigma,
                        long double Delta,
                        long double kappa,
                        long double sin_theta,
                        long double func_theta,
                        long double F1_r,
                        long double func_r) nogil

cdef long double g_14_l(long double r,
                        long double theta,
                        long double r_s,
                        long double a,
                        long double Sigma,
                        long double Delta) nogil

cdef long double g_22_l(long double r,
                        long double theta,
                        long double Sigma,
                        long double Delta,
                        long double kappa,
                        long double func_theta,
                        long double F1_r,
                        long double func_r) nogil

cdef long double g_33_l(long double r,
                        long double theta,
                        long double Sigma,
                        long double kappa,
                        long double func_theta,
                        long double F2_r) nogil

cdef long double g_44_l(long double r,
                        long double theta,
                        long double a,
                        long double Sigma,
                        long double Delta,
                        long double kappa,
                        long double sin_theta,
                        long double func_theta,
                        long double F2_r) nogil

cdef long double det_g_l(long double r,
                         long double theta,
                         long double r_s,
                         long double a,
                         long double Sigma,
                         long double Delta,
                         long double kappa,
                         long double func_theta,
                         long double F1_r,
                         long double F2_r,
                         long double func_r,
                         long double sin_theta) nogil