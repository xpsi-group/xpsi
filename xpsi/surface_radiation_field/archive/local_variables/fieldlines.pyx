#cython: cdivision=True
#cython: boundscheck=False
#cython: nonecheck=False
#cython: wraparound=False

from libc.math cimport sin, cos, acos, log, M_PI
from libc.stdlib cimport malloc, free
from libc.stdio cimport printf, fopen, fclose, fscanf, FILE, EOF, rewind

from .effective_gravity_universal cimport effectiveGravity

cdef int SUCCESS = 0
cdef int ERROR = 1
cdef int HIT = 1
cdef int MISS = 0

ctypedef struct fieldlines:
    size_t num_points
    double *colatitude
    double *azimuth
    int *OPEN

cdef int HIT_or_MISS(double theta,
                     double phi,
                     double HYPERSLICE,
                     const double *const global_variables,
                     const storage *const buf) nogil:
    # Use this function to determine whether an arbitrary input ray transports
    # a FINITE quantity of radiation

    cdef size_t i, POINT_INDEX
    cdef double separation
    cdef double min_separation = -1.0
    cdef fieldlines *ptr = <fieldlines*> buf.ptr

    cdef:
        double sin_theta = sin(theta)
        double THETA, sin_THETA, PHI

    for i in range(ptr.num_points):
        THETA = ptr.colatitude[i]
        sin_THETA = sin(THETA)
        PHI = ptr.azimuth[i] + HYPERSLICE

        separation = acos(cos(THETA) * cos(theta)
                            + sin_THETA * sin_theta * cos(phi) * cos(PHI)
                            + sin_THETA * sin_theta * sin(phi) * sin(PHI))

        if separation < min_separation or min_separation < 0.0:
            min_separation = separation
            POINT_INDEX = i

    return ptr.OPEN[POINT_INDEX]

cdef storage* init_local_variables(size_t numTHREADS) nogil:

    cdef storage *buf = <storage*> malloc(sizeof(storage))
    buf.local_variables = <double**> malloc(numTHREADS * sizeof(double*))

    cdef size_t THREAD

    for THREAD in range(numTHREADS):
        buf.local_variables[THREAD] = <double*> malloc(2 * sizeof(double))

    cdef FILE *fp
    # a sequence of points spanning surface with columns:
    # azimuth (degrees), colatitude (degrees), openness (open is True)
    # format discussed with and tested with calculations by Anna Bilous (2020)
    # crude nearest-neighbour interpolation in field-line openness
    fp = fopen("./surface_points.txt", "r") # (insert the correct path)

    load.num_points = 0
    cdef char boolean[10]
    cdef double azimuth, colatitude

    while (fscanf(fp, "%lf %lf %s ",
                  &azimuth,
                  &colatitude,
                  boolean) != EOF):

        load.num_points += 1

    cdef fieldlines *load = <fieldlines*> malloc(sizeof(fieldlines))
    load.colatitude = <double*> malloc(load.num_points * sizeof(double))
    load.azimuth = <double*> malloc(load.num_points * sizeof(double))
    load.OPEN = <int*> malloc(load.num_points * sizeof(int))

    rewind(fp)
    cdef size_t i = 0

    while (fscanf(fp, "%lf %lf %s ",
                 &(load.azimuth[i]),
                 &(load.colatitude[i]),
                 boolean) != EOF):

        load.colatitude[i] *= M_PI/180.0
        load.azimuth[i] *= M_PI/180.0

        if boolean[0] == "T": # True
            load.OPEN[i] = HIT
        else:
            load.OPEN[i] = MISS

        i += 1

    fclose(fp)

    # save a ref
    buf.ptr = load

    return buf

cdef int free_local_variables(size_t numTHREADS, storage *const buf) nogil:

    cdef size_t THREAD

    for THREAD in range(numTHREADS):
        free(buf.local_variables[THREAD])

    free(buf.local_variables)

    cdef fieldlines *ptr = <fieldlines*> buf.ptr

    free(ptr.colatitude)
    free(ptr.azimuth)
    free(ptr.OPEN)
    free(ptr)

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

    local_vars[0] = global_variables[0]

    local_vars[1] = effectiveGravity(cos(theta),
                                     GEOM.R_eq,
                                     GEOM.zeta,
                                     GEOM.epsilon)

    return SUCCESS
