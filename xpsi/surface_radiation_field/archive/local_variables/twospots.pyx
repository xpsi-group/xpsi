#cython: cdivision=True
#cython: boundscheck=False
#cython: nonecheck=False
#cython: wraparound=False

from libc.math cimport sin, cos, acos, log
from libc.stdlib cimport malloc, free
from libc.stdio cimport printf

from .effective_gravity_universal cimport effectiveGravity

cdef int SUCCESS = 0
cdef int ERROR = 1
cdef int HIT = 1
cdef int MISS = 0

cdef int HIT_or_MISS(double theta,
                     double phi,
                     double HYPERSLICE,
                     const double *const global_variables,
                     const storage *const buf) nogil:
    # Use this function to determine whether an arbitrary input ray transports
    # a FINITE quantity of radiation

    cdef:
        double THETA = global_variables[0]
        double PHI = global_variables[1] + HYPERSLICE
        double angRad = global_variables[2]

        double sin_theta = sin(theta)
        double sin_THETA = sin(THETA)
        double psi = acos(cos(THETA) * cos(theta)
                            + sin_THETA * sin_theta * cos(phi) * cos(PHI)
                            + sin_THETA * sin_theta * sin(phi) * sin(PHI))

    if psi <= angRad:
        return HIT # in the primary circular hot region
    else:
        THETA = global_variables[4]
        PHI = global_variables[5] + HYPERSLICE
        angRad = global_variables[6]

        sin_theta = sin(theta)
        sin_THETA = sin(THETA)
        psi = acos(cos(THETA) * cos(theta)
                            + sin_THETA * sin_theta * cos(phi) * cos(PHI)
                            + sin_THETA * sin_theta * sin(phi) * sin(PHI))

        if psi <= angRad:
            return HIT # in the secondary circular hot region
        else:
            return MISS # both hot regions

#----------------------------------------------------------------------->>>
# >>> User modifiable functions.
# >>> Note that the user is entirely free to wrap thread-safe and
# ... non-parallel external C routines from an external library.
# >>> Thus the bodies of the following need not be written explicitly in
# ... the Cython language.
#----------------------------------------------------------------------->>>
cdef storage* init_local_variables(size_t numTHREADS,
                                   const char *const filepath) nogil:

    cdef storage *buf = <storage*> malloc(sizeof(storage))
    buf.local_variables = <double**> malloc(numTHREADS * sizeof(double*))

    cdef size_t THREAD

    for THREAD in range(numTHREADS):
        buf.local_variables[THREAD] = <double*> malloc(2 * sizeof(double))

    # additional memory not required
    buf.ptr = NULL

    return buf

cdef int free_local_variables(size_t numTHREADS, storage *const buf) nogil:

    cdef size_t THREAD

    for THREAD in range(numTHREADS):
        free(buf.local_variables[THREAD])

    free(buf.local_variables)

    # remember to free dynamically allocated memory that buf.ptr points at

    free(buf)

    return SUCCESS

cdef int eval_local_variables(double theta,
                              double phi,
                              double HYPERSLICE,
                              const _GEOM *const GEOM,
                              const double *const global_variables,
                              const storage *const buf,
                              size_t THREAD) nogil:

    cdef double *local_vars = (<double**> buf.local_variables)[THREAD]

    cdef:
        double THETA = global_variables[0]
        double PHI = global_variables[1] + HYPERSLICE
        double angRad = global_variables[2]

        double sin_theta = sin(theta)
        double sin_THETA = sin(THETA)
        double psi = acos(cos(THETA) * cos(theta)
                            + sin_THETA * sin_theta * cos(phi) * cos(PHI)
                            + sin_THETA * sin_theta * sin(phi) * sin(PHI))

    if psi <= angRad:
        # uncomment next line for temperature decay towards boundary
        local_vars[0] = global_variables[3] #- 0.5*psi*psi/(angRad*angRad*log(10.0))
    else:
        THETA = global_variables[4]
        PHI = global_variables[5] + HYPERSLICE
        angRad = global_variables[6]

        sin_theta = sin(theta)
        sin_THETA = sin(THETA)
        psi = acos(cos(THETA) * cos(theta)
                            + sin_THETA * sin_theta * cos(phi) * cos(PHI)
                            + sin_THETA * sin_theta * sin(phi) * sin(PHI))

        if psi <= angRad:
            # uncomment next line for temperature decay towards boundary
            local_vars[0] = global_variables[7] #- 0.5*psi*psi/(angRad*angRad*log(10.0))

    local_vars[1] = effectiveGravity(cos(theta),
                                     GEOM.R_eq,
                                     GEOM.zeta,
                                     GEOM.epsilon)

    return SUCCESS
