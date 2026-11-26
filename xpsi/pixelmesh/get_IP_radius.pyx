#cython: cdivision=True
#cython: boundscheck=False
#cython: nonecheck=False
#cython: wraparound=False

## Size and centre of the image of the star on the image plane, used by
## globalRayMap.pyx to lay out the pixel grid.
##
## Method: four one-dimensional bisections for the edge of the image, each
## using full ray traces (RK_IP2S_tracer.RK); a ray 'hits' if it ends on the
## surface (STATE[2] > 0):
##   L, R: left and right edges along the X axis (Y = 0),
##   U, D: upper and lower edges along X = ORIGIN_X (elliptical screen only).
## Each edge is pushed outwards by 0.5% as a safety margin. The screen is
## then the ellipse centred on (ORIGIN_X, ORIGIN_Y) with these semi-axes
## (or a circle of radius SEMI_MAJOR if force_circular == 1;
## integrator.pyx calls this with force_circular = 0).
##
## Image-plane coordinates are in metres. b_max (integrator.pyx) is the
## Schwarzschild apparent radius of a sphere of radius R_eq,
## R_eq / sqrt(1 - r_s/R_eq); the search brackets [0.5, 1.75] b_max allow
## for the shift and distortion of the image by rotation.
##
## Checked for 1.4 Msun, 12 km, 600 Hz at inclinations 30, 60 and
## 85 deg: the ellipse contains the whole star, which reaches at most
## 0.9986 of the ellipse radius (at 85 deg). So the 0.5% margin is
## sufficient there, but leaves only ~0.15% to spare; for faster spins or
## more distorted images it may need increasing.
##
## Caveats:
##   - The edges are found along two lines only; the true extreme of the
##     image need not lie on them (hence the margin).
##   - A ray that FAILS (flags -102 to -105) counts as a miss, so failures
##     near the limb shrink the screen; if every ray failed, these
##     bisections would collapse onto the inner bracket 0.5 b_max.
##   - The function is declared double but returns nothing; the caller
##     ignores the return value.
from libc.math cimport fabs, sqrt, sin, cos
from libc.stdio cimport printf
from xpsi.pixelmesh.RK_IP2S_tracer cimport RK

cdef double compute_imagePlane_radius(const _GEOM *const GEOM,
                                      _RAY *const RAY,
                                      RAY_MAP *const MAP,
                                      int force_circular,
                                      int obl_surfgrav_ind) noexcept nogil:

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
        ## Stop each bisection once the bracket is narrower than 1e-12 b_max
        ## (~41 iterations); the bracket reaches double precision after ~55, so
        ## further iterations only cost ray traces. 'terminate' is a safety cap.
        double tolerance = 1.0e-12 * b_max

    num_bisections = 0
    while num_bisections < terminate and r2 - r1 > tolerance:
        num_bisections += 1

        L = 0.5 * (r1 + r2)

        RAY.X_IP = L
        RAY.Y_IP = 0.0

        status = RK(RAY, GEOM, obl_surfgrav_ind)

        if RAY.STATE[2] > 0.0:
            r2 = L
        else:
            r1 = L

        # printf("\nr/b_max = %.8f; R(photosphere) = %.8f", r/b_max, RAY.STATE[2])
        # printf("\nNumber of ray steps: %d", RAY.NUMSTEPS)
        # printf("\nb: %.8e", RAY.IMPACT/b_max)
        # printf("\nSpin parameter, a: %.8e", GEOM.a)
        # printf("\nMass quadrupole deviation, kappa: %.8e", GEOM.kappa)

    ## 0.5% safety margin outwards (L < 0)
    L *= 1.005

    r1 = 0.5 * b_max
    r2 = 1.75 * b_max

    num_bisections = 0
    while num_bisections < terminate and r2 - r1 > tolerance:
        num_bisections += 1

        R = 0.5 * (r1 + r2)

        RAY.X_IP = R
        RAY.Y_IP = 0.0

        status = RK(RAY, GEOM, obl_surfgrav_ind)

        if RAY.STATE[2] > 0.0:
            r1 = R
        else:
            r2 = R

        # printf("\nr/b_max = %.8f; R(photosphere) = %.8f", r/b_max, RAY.STATE[2])
        # printf("\nNumber of ray steps: %d", RAY.NUMSTEPS)
        # printf("\nb: %.8e", RAY.IMPACT/b_max)
        # printf("\nSpin parameter, a: %.8e", GEOM.a)
        # printf("\nMass quadrupole deviation, kappa: %.8e", GEOM.kappa)

    ## 0.5% safety margin outwards
    R *= 1.005

    ## Horizontal centre of the screen. Rotation (frame dragging and Doppler
    ## boosting of the geometry) shifts the image sideways, so this is not 0.
    MAP.ORIGIN_X = 0.5 * (L + R)

    ## (Threshold in metres: effectively only removes round-off.)
    if fabs(MAP.ORIGIN_X) < 1.0e-4:
        MAP.ORIGIN_X = 0.0

    #printf("\nx origin: %.8e", MAP.ORIGIN_X)

    MAP.SEMI_MAJOR = fabs(L - MAP.ORIGIN_X)

    if force_circular == 1:
        ## Circular screen
        MAP.SEMI_MINOR = MAP.SEMI_MAJOR
        MAP.ORIGIN_Y = 0.0
        printf("\nImage plane radius: %.2f", MAP.SEMI_MAJOR / b_max)
    else:

        printf("\nImage plane semi-major: %.2f", fabs(MAP.SEMI_MAJOR / b_max))

        r1 = 0.5 * b_max
        r2 = 1.75 * b_max

        num_bisections = 0
        while num_bisections < terminate and r2 - r1 > tolerance:
            num_bisections += 1

            U = 0.5 * (r1 + r2)

            RAY.X_IP = MAP.ORIGIN_X
            RAY.Y_IP = U

            status = RK(RAY, GEOM, obl_surfgrav_ind)

            if RAY.STATE[2] > 0.0:
                r1 = U
            else:
                r2 = U

            # printf("\nr/b_max = %.8f; R(photosphere) = %.8f", r/b_max, RAY.STATE[2])
            # printf("\nNumber of ray steps: %d", RAY.NUMSTEPS)
            # printf("\nb: %.8e", RAY.IMPACT/b_max)
            # printf("\nSpin parameter, a: %.8e", GEOM.a)
            # printf("\nMass quadrupole deviation, kappa: %.8e", GEOM.kappa)

        ## 0.5% safety margin outwards
        U *= 1.005

        r1 = -1.75 * b_max
        r2 = -0.5 * b_max

        num_bisections = 0
        while num_bisections < terminate and r2 - r1 > tolerance:
            num_bisections += 1

            D = 0.5 * (r1 + r2)

            RAY.X_IP = MAP.ORIGIN_X
            RAY.Y_IP = D

            status = RK(RAY, GEOM, obl_surfgrav_ind)

            if RAY.STATE[2] > 0.0:
                r2 = D
            else:
                r1 = D

            # printf("\nr/b_max = %.8f; R(photosphere) = %.8f", r/b_max, RAY.STATE[2])
            # printf("\nNumber of ray steps: %d", RAY.NUMSTEPS)
            # printf("\nb: %.8e", RAY.IMPACT/b_max)
            # printf("\nSpin parameter, a: %.8e", GEOM.a)
            # printf("\nMass quadrupole deviation, kappa: %.8e", GEOM.kappa)

        ## 0.5% safety margin outwards (D < 0)
        D *= 1.005

        MAP.ORIGIN_Y = 0.5 * (U + D)

        MAP.SEMI_MINOR = fabs(U - MAP.ORIGIN_Y)

        printf("\nImage plane semi-minor: %.2f", fabs(MAP.SEMI_MINOR / b_max))
