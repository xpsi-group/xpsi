from GSL cimport gsl_odeiv2_step, gsl_odeiv2_control, gsl_odeiv2_evolve, gsl_odeiv2_system
from .geometricConfiguration cimport _GEOM

ctypedef struct _RAY:
    double AFFINE
    double BISECT_STATE[6]
    double STATE[7]
    double PREVIOUS[6]
    double PREVIOUS_AFFINE
    double IMPACT
    double XI
    double PREVIOUS_XI
    gsl_odeiv2_step *S
    gsl_odeiv2_control *C
    gsl_odeiv2_evolve *E
    double INIT_STEP
    double STEP
    double PREVIOUS_STEP
    size_t MAXSTEPS
    double MAX_AFFINE
    double X_IP
    double Y_IP
    gsl_odeiv2_system SYS
    double PARAMS[7]
    int NUMSTEPS
    int EVOLVE
    int NUM_SINGULARITY_STEPS
    int NUM_SINGULARITIES
    double SINGULARITY_AFFINE_STEP_SIZE
    double dSTATEdAFFINE[6]

cdef _RAY* alloc_RAY(double epsabs_ray,
                     double epsrel_ray,
                     double INIT_STEP,
                     size_t MAXSTEPS) nogil

cdef void free_RAY(_RAY *const RAY) nogil

cdef int RK(_RAY *const RAY,
            const _GEOM *const GEOM) nogil
