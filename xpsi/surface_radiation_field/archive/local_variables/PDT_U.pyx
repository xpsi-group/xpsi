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

cdef double separation(double THETA, double PHI,
                       double theta, double phi, double sin_theta) noexcept nogil:

    cdef:
        double sin_THETA = sin(THETA)
        double psi = acos(cos(THETA) * cos(theta)
                            + sin_THETA * sin_theta * cos(phi) * cos(PHI)
                            + sin_THETA * sin_theta * sin(phi) * sin(PHI))
    return psi

cdef int HIT_or_MISS(double theta,
                     double phi,
                     double HYPERSLICE,
                     const double *const global_variables,
                     const storage *const buf) noexcept nogil:
    # Use this function to determine whether an arbitrary input ray transports
    # a FINITE quantity of radiation

    cdef:
        double super_colat_p = global_variables[0]
        double super_azi_p = global_variables[1] + HYPERSLICE
        double super_radius_p = global_variables[2]
        double cede_colat_p = global_variables[3]
        double cede_azi_p = global_variables[4] + HYPERSLICE
        double cede_radius_p = global_variables[5]

        double super_colat_s = global_variables[6]
        double super_azi_s = global_variables[7] + HYPERSLICE
        double super_radius_s = global_variables[8]
        double cede_colat_s = global_variables[9]
        double cede_azi_s = global_variables[10] + HYPERSLICE
        double cede_radius_s = global_variables[11]

        double super_log10_temperature_p = global_variables[12]
        double super_log10_temperature_s = global_variables[14]

        double sin_theta = sin(theta)

    if separation(cede_colat_p, cede_azi_p,
                  theta, phi, sin_theta) <= cede_radius_p:
        if separation(super_colat_p, super_azi_p,
                      theta, phi, sin_theta) > super_radius_p:
            return HIT # in the primary radiating region
        elif super_log10_temperature_p >= 0.0:
            return HIT
    elif super_log10_temperature_p >= 0.0:
        if separation(super_colat_p, super_azi_p,
                      theta, phi, sin_theta) <= super_radius_p:
            return HIT # in the primary radiating region

    if separation(cede_colat_s, cede_azi_s,
                  theta, phi, sin_theta) <= cede_radius_s:
        if separation(super_colat_s, super_azi_s,
                      theta, phi, sin_theta) > super_radius_s:
            return HIT # in the secondary radiating region
        elif super_log10_temperature_s >= 0.0:
            return HIT
    elif super_log10_temperature_s >= 0.0:
        if separation(super_colat_s, super_azi_s,
                      theta, phi, sin_theta) <= super_radius_s:
            return HIT # in the primary radiating region

    return MISS # both hot regions

#----------------------------------------------------------------------->>>
# >>> User modifiable functions.
# >>> Note that the user is entirely free to wrap thread-safe and
# ... non-parallel external C routines from an external library.
# >>> Thus the bodies of the following need not be written explicitly in
# ... the Cython language.
#----------------------------------------------------------------------->>>
cdef storage* init_local_variables(size_t numTHREADS,
                                   const char *const filepath) noexcept nogil:

    cdef storage *buf = <storage*> malloc(sizeof(storage))
    buf.local_variables = <double**> malloc(numTHREADS * sizeof(double*))

    cdef size_t THREAD

    for THREAD in range(numTHREADS):
        buf.local_variables[THREAD] = <double*> malloc(2 * sizeof(double))

    # additional memory not required
    buf.ptr = NULL

    return buf

cdef int free_local_variables(size_t numTHREADS, storage *const buf) noexcept nogil:

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
                              size_t THREAD) noexcept nogil:

    cdef double *local_vars = (<double**> buf.local_variables)[THREAD]

    cdef:
        double super_colat_p = global_variables[0]
        double super_azi_p = global_variables[1] + HYPERSLICE
        double super_radius_p = global_variables[2]
        double cede_colat_p = global_variables[3]
        double cede_azi_p = global_variables[4] + HYPERSLICE
        double cede_radius_p = global_variables[5]

        double super_colat_s = global_variables[6]
        double super_azi_s = global_variables[7] + HYPERSLICE
        double super_radius_s = global_variables[8]
        double cede_colat_s = global_variables[9]
        double cede_azi_s = global_variables[10] + HYPERSLICE
        double cede_radius_s = global_variables[11]

        double super_log10_temperature_p = global_variables[12]
        double cede_log10_temperature_p = global_variables[13]
        double super_log10_temperature_s = global_variables[14]
        double cede_log10_temperature_s = global_variables[15]

        double sin_theta = sin(theta)

    if separation(cede_colat_p, cede_azi_p,
                  theta, phi, sin_theta) <= cede_radius_p:
        if separation(super_colat_p, super_azi_p,
                      theta, phi, sin_theta) > super_radius_p:
            local_vars[0] = cede_log10_temperature_p
        elif super_log10_temperature_p >= 0.0:
            local_vars[0] = super_log10_temperature_p
    elif super_log10_temperature_p >= 0.0:
        if separation(super_colat_p, super_azi_p,
                      theta, phi, sin_theta) <= super_radius_p:
            local_vars[0] = super_log10_temperature_p

    if separation(cede_colat_s, cede_azi_s,
                  theta, phi, sin_theta) <= cede_radius_s:
        if separation(super_colat_s, super_azi_s,
                      theta, phi, sin_theta) > super_radius_s:
            local_vars[0] = cede_log10_temperature_s
        elif super_log10_temperature_s >= 0.0:
            local_vars[0] = super_log10_temperature_s
    elif super_log10_temperature_s >= 0.0:
        if separation(super_colat_s, super_azi_s,
                      theta, phi, sin_theta) <= super_radius_s:
            local_vars[0] = super_log10_temperature_s

    local_vars[1] = effectiveGravity(cos(theta),
                                     GEOM.R_eq,
                                     GEOM.zeta,
                                     GEOM.epsilon,
                                     GEOM.star_shape_ind)

    return SUCCESS
