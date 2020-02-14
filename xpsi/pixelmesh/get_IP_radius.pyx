#cython: cdivision=True
#cython: boundscheck=False
#cython: nonecheck=False
#cython: wraparound=False

from libc.math cimport fabs, sqrt, sin, cos
from libc.stdio cimport printf
from xpsi.pixelmesh.RK_IP2S_tracer cimport RK

cdef double compute_imagePlane_radius(const _GEOM *const GEOM,
                                      _RAY *const RAY,
                                      RAY_MAP *const MAP,
                                      int force_circular) nogil:

    # Bisection algorithm to compute an appropriate radial extent for the
    # image plane

    printf("\nCalculating image plane boundary...")

    cdef:
        size_t i, num_bisections
        int status
        double b_max = GEOM.b_max
        double r, r1 = -1.75 * b_max, r2 = -0.5 * b_max
        double L, R, U, D
        size_t terminate = 1000

    num_bisections = 0
    while num_bisections < terminate:
        num_bisections += 1

        L = 0.5 * (r1 + r2)

        RAY.X_IP = L
        RAY.Y_IP = 0.0

        status = RK(RAY, GEOM)

        if RAY.STATE[2] > 0.0:
            r2 = L
        else:
            r1 = L

        # printf("\nr/b_max = %.8f; R(photosphere) = %.8f", r/b_max, RAY.STATE[2])
        # printf("\nNumber of ray steps: %d", RAY.NUMSTEPS)
        # printf("\nb: %.8e", RAY.IMPACT/b_max)
        # printf("\nSpin parameter, a: %.8e", GEOM.a)
        # printf("\nMass quadrupole deviation, kappa: %.8e", GEOM.kappa)

    L *= 1.005

    r1 = 0.5 * b_max
    r2 = 1.75 * b_max

    num_bisections = 0
    while num_bisections < terminate:
        num_bisections += 1

        R = 0.5 * (r1 + r2)

        RAY.X_IP = R
        RAY.Y_IP = 0.0

        status = RK(RAY, GEOM)

        if RAY.STATE[2] > 0.0:
            r1 = R
        else:
            r2 = R

        # printf("\nr/b_max = %.8f; R(photosphere) = %.8f", r/b_max, RAY.STATE[2])
        # printf("\nNumber of ray steps: %d", RAY.NUMSTEPS)
        # printf("\nb: %.8e", RAY.IMPACT/b_max)
        # printf("\nSpin parameter, a: %.8e", GEOM.a)
        # printf("\nMass quadrupole deviation, kappa: %.8e", GEOM.kappa)

    R *= 1.005

    MAP.ORIGIN_X = 0.5 * (L + R)

    if fabs(MAP.ORIGIN_X) < 1.0e-4:
        MAP.ORIGIN_X = 0.0

    #printf("\nx origin: %.8e", MAP.ORIGIN_X)

    MAP.SEMI_MAJOR = fabs(L - MAP.ORIGIN_X)

    if force_circular == 1:
        MAP.SEMI_MINOR = MAP.SEMI_MAJOR
        MAP.ORIGIN_Y = 0.0
        printf("\nImage plane radius: %.2f", MAP.SEMI_MAJOR / b_max)
    else:

        printf("\nImage plane semi-major: %.2f", fabs(MAP.SEMI_MAJOR / b_max))

        r1 = 0.5 * b_max
        r2 = 1.75 * b_max

        num_bisections = 0
        while num_bisections < terminate:
            num_bisections += 1

            U = 0.5 * (r1 + r2)

            RAY.X_IP = MAP.ORIGIN_X
            RAY.Y_IP = U

            status = RK(RAY, GEOM)

            if RAY.STATE[2] > 0.0:
                r1 = U
            else:
                r2 = U

            # printf("\nr/b_max = %.8f; R(photosphere) = %.8f", r/b_max, RAY.STATE[2])
            # printf("\nNumber of ray steps: %d", RAY.NUMSTEPS)
            # printf("\nb: %.8e", RAY.IMPACT/b_max)
            # printf("\nSpin parameter, a: %.8e", GEOM.a)
            # printf("\nMass quadrupole deviation, kappa: %.8e", GEOM.kappa)

        U *= 1.005

        r1 = -1.75 * b_max
        r2 = -0.5 * b_max

        num_bisections = 0
        while num_bisections < terminate:
            num_bisections += 1

            D = 0.5 * (r1 + r2)

            RAY.X_IP = MAP.ORIGIN_X
            RAY.Y_IP = D

            status = RK(RAY, GEOM)

            if RAY.STATE[2] > 0.0:
                r2 = D
            else:
                r1 = D

            # printf("\nr/b_max = %.8f; R(photosphere) = %.8f", r/b_max, RAY.STATE[2])
            # printf("\nNumber of ray steps: %d", RAY.NUMSTEPS)
            # printf("\nb: %.8e", RAY.IMPACT/b_max)
            # printf("\nSpin parameter, a: %.8e", GEOM.a)
            # printf("\nMass quadrupole deviation, kappa: %.8e", GEOM.kappa)

        D *= 1.005

        MAP.ORIGIN_Y = 0.5 * (U + D)

        MAP.SEMI_MINOR = fabs(U - MAP.ORIGIN_Y)

        printf("\nImage plane semi-minor: %.2f", fabs(MAP.SEMI_MINOR / b_max))
