#cython: cdivision=True
#cython: boundscheck=False
#cython: nonecheck=False
#cython: wraparound=False

from libc.math cimport sin, cos, acos, log, exp, pow, sqrt
from libc.stdlib cimport malloc, free
from libc.stdio cimport printf

from xpsi.tools.effective_gravity_universal cimport effectiveGravity

cdef int SUCCESS = 0
cdef int ERROR = 1
cdef int HIT = 1
cdef int MISS = 0

cdef int HIT_or_MISS(double theta,
                     double phi,
                     double HYPERSLICE,
                     const double *const global_variables,
                     const storage *const buf) nogil:
    # use this function to determine whether an arbitrary input ray transports
    # a FINITE quantity of radiation
    # NB: HYPERSLICE is a phase parameter with units of radians that the
    # time evolution of the photosphere as it undergoes bulk rotation

    return HIT # independent of spacetime coordinates

#----------------------------------------------------------------------->>>
# >>> User modifiable functions.
# >>> Note that the user is entirely free to wrap thread-safe and
# ... non-parallel external C routines from an external library.
# >>> Thus the bodies of the following need not be written explicitly in
# ... the Cython language.
#----------------------------------------------------------------------->>>
cdef storage* init_local_variables(size_t numTHREADS) nogil:

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

    # r-mode implementation below contributed by F.R.N. Chambers (2020)
    cdef:
        double temperature = global_variables[0]
        double phase = global_variables[1] + HYPERSLICE
        double amplitude = global_variables[2]
        double q = global_variables[3]
        double f, L, mu_sq, alpha, beta, N

    L = -(1.0 + q) / (3.0 * q)

    alpha = q * (L + 1.0)

    beta = q * L

    N = -1.0 * beta / (2.0 * alpha) * exp(1.0 + beta / (2.0 * alpha))

    mu_sq = pow(cos(theta), 2.0)

    f = N * (1.0 - alpha * mu_sq) * exp(-0.5 * beta * mu_sq)

    local_vars[0] = temperature + log((1.0 + amplitude * f * cos(phi - phase))) / log(10.0)

    local_vars[1] = effectiveGravity(cos(theta),
                                     GEOM.R_eq,
                                     GEOM.zeta,
                                     GEOM.epsilon)

    return SUCCESS
