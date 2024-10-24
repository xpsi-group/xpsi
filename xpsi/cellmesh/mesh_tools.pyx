#cython: cdivision=True
#cython: boundscheck=False
#cython: nonecheck=False
#cython: wraparound=False

from libc.math cimport M_PI, M_PI_2
from libc.math cimport sqrt, sin, cos, tan, asin, acos, atan, fabs, exp, pow
from libc.math cimport fmin, fmax, floor, ceil

cdef double M_2PI = 2.0 * M_PI
cdef double _c = 2.99792458e8
cdef double _keV = 1.60217662e-16
cdef double _k_B = 1.38064852e-23
cdef double _h = 6.62607004e-34

cdef double radiusNormalised(double mu, 
                             double epsilon, 
                             double zeta, 
                             double R_eq) nogil: # ADDED! 
    """
    Calculate the normalised radius (R(theta)/R_eq) based on the slow-elliptical approximation from Silva et al. (2021) 
    (see equation 16 / R_eq). This approximation is obtained using stars with epsilon <= 0.25.   

    :param mu: Colatitude (cos(theta))
    :type mu: double

    :param epsilon: Dimensionless spin parameter = (omega**2*R_eq**3)/(G*M) defined as sigma in Silva et al. (2021).
    :type epsilon: double

    :param zeta: Compactness = (G*M)/(R_eq*c**2) defined as zeta in AlGendy & Morsink (2007) and as 
        kappa in Silva et al. (2021) 
    :type zeta: double 

    :param R_eq: Equatorial radius of the star 
    :type R_eq: double

    :return: Normalised radius (R(theta)/R_eq)
    :rtype: double 
    """
    cdef : 
        double esq = epsilon # why do this? 
        double esqsq = epsilon * epsilon # why seperately define this? 
    
        # slow-elliptical fit (epsilon<=0.25)
        double e = 1.089 * sqrt(esq) + 0.168 * esq - 0.685 * esq * zeta - 0.802 * esqsq # eccentricity of the star 
        
        # fast-elliptical fit (epsilon>0.25, corresponds to min ~700-800 Hz) 
        # double e = 0.251 + 0.935 * esq + 0.709 * zeta + 0.030 * esq * zeta - 0.472 * esqsq - 2.427 * zeta * zeta

        # slow-elliptical fit (epsilon<=0.25)
        double e = 1.089 * sqrt(esq) + 0.168 * esq - 0.685 * esq * zeta - 0.802 * esqsq
        double a_2 = -1.013 - 0.312 *  esq + 0.930 * esq * zeta - 1.596 * esqsq 
        double a_4 = 0.016 + 0.301 *  esq - 1.261 * esq * zeta + 2.728 * esqsq 

        # fast-elliptical fit (epsilon>0.25) 
        # double e = 0.251 + 0.935 * esq + 0.709 * zeta + 0.030 * esq * zeta - 0.472 * esqsq - 2.427 * zeta * zeta
        # double a2 = -1.265 + 0.220 * esq + 2.651 * zeta + 1.010 * esq * zeta - 1.815 * esqsq - 7.657 * zeta * zeta
        # double a4 = 0.556 - 1.465 * esq - 4.260 * zeta - 2.327 * esq * zeta + 4.921 * esqsq + 12.98 * zeta * zeta

        double g = 1 + a_2 * mu**2 + a_4 * pow(mu, 4) - (1 + a2_ + a_4) * pow(mu, 6)

    return sqrt((1 - e**2) / (1 - e**2 * g)) 

cdef double f_theta(double mu,
                    double radiusNormed, 
                    double epsilon,
                    double zeta) nogil:
    """
    Calculate f(theta) based on the slow-elliptical approximation from Silva et al. (2021) for the derived
    normalised radius (see equation 16 for the un-normalised radius equation). This approximation is 
    obtained using stars with epsilon <= 0.25.   

    :param mu: Colatitude (cos(theta))
    :type mu: double

    :param radiusNormed: Normalised radius (R(theta)/R_eq) (see radiusNormalised)
    :type R_eq: double

    :param epsilon: Dimensionless spin parameter = (omega**2*R_eq**3)/(G*M) defined as sigma in Silva et al. (2021).
    :type epsilon: double

    :param zeta: Compactness = (G*M)/(R_eq*c**2) defined as zeta in AlGendy & Morsink (2007) and as 
        kappa in Silva et al. (2021) 
    :type zeta: double 

    :return: Normalised radius (R(theta)/R_eq)
    :rtype: double 
    """

    cdef double radiusDerivNormed

    double esq = epsilon 
    double esqsq = epsilon * epsilon 

    # slow-elliptical fit (epsilon<=0.25)
    double e = 1.089 * sqrt(esq) + 0.168 * esq - 0.685 * esq * zeta - 0.802 * esqsq
    double a_2 = -1.013 - 0.312 *  esq + 0.930 * esq * zeta - 1.596 * esqsq 
    double a_4 = 0.016 + 0.301 *  esq - 1.261 * esq * zeta + 2.728 * esqsq 

    radiusDerivNormed = (e**2 * sqrt(1 - mu**2)) / (2 * radiusNormed * (1 - e**2 * g)) * (2 * a_2 * mu + 4 * a_4 * pow(mu,3) - 6 * (1 + a_2 + a_4) * pow(mu,5))

    return radiusDerivNormed / (radiusNormed * sqrt(1.0 - 2.0 * zeta / radiusNormed))

cdef double integrand(double theta, void *params) nogil:

    if theta == 0.0:
        return 0.0

    cdef:
        double mu = cos(theta)
        double eval_integrand, radius, f
        double epsilon = (<double*> params)[0]
        double zeta = (<double*> params)[1]
        int av = <int>(<double*> params)[2]

    radius = radiusNormalised(mu, epsilon, zeta)
    f = f_theta(mu, radius, epsilon, zeta)

    if (av == 0):
        eval_integrand = radius * radius * sqrt(1.0 + f * f) * sin(theta)
    else:
        eval_integrand = theta * radius * radius * sqrt(1.0 + f * f) * sin(theta)

    return eval_integrand

cdef double integrateArea(double lower_lim,
                          double upper_lim,
                          double R_eq,
                          double epsilon,
                          double zeta,
                          int av,
                          gsl_integration_cquad_workspace *w) nogil:

    cdef double params[3]
    cdef double integral, error
    cdef size_t nevals

    params[0] = epsilon
    params[1] = zeta
    params[2] = <double> av

    cdef gsl_function F
    F.function = integrand
    F.params = &params

    gsl_integration_cquad(&F,
                          lower_lim,
                          upper_lim,
                          0,
                          1.0e-8,
                          w,
                          &integral,
                          &error,
                          &nevals)

    return R_eq * R_eq * integral

cdef double eval_psi(double theta, double phi, double THETA) nogil:

    cdef double psi = acos(cos(THETA) * cos(theta) + sin(THETA) * sin(theta) * cos(phi))

    return psi

cdef double eval_cedeAzi(double theta, double phi, double psi, double THETA) nogil:

    cdef double azi = sin(theta) * sin(phi) / sin(psi)

    if azi < -1.0 or azi > 1.0:
        return 0.0
    else:
        azi = asin(azi)

    cdef double critical

    if THETA <= 0.0:
        if cos(theta) == 0.0:
            critical = M_PI_2
        else:
            critical = atan(tan(theta) * cos(phi))
            if critical > 0.0:
                critical -= M_PI

        #with gil:
        #    print('critical rotation = %.8e' % critical)
        #    print('phi = %.8e' % phi)

        if phi < 0.0:
            while phi < 0.0:
                phi += M_2PI
        elif phi > M_2PI:
            while phi > M_2PI:
                phi -= M_2PI

        if THETA >= critical:
            if phi <= M_PI_2:
                return azi
            elif M_PI_2 < phi < M_PI:
                return M_PI - azi
            elif M_PI < phi <= 3.0*M_PI_2:
                return M_PI - azi
            elif 3.0*M_PI_2 < phi < M_2PI:
                return azi + 2.0*M_PI
            elif phi == 0.0 or phi == M_2PI:
                return 0.0
            elif phi == M_PI:
                return M_PI
        elif THETA < critical:
            if phi <= M_PI_2:
                return M_PI - azi
            elif M_PI_2 < phi < M_PI:
                return azi
            elif M_PI < phi <= 3.0*M_PI_2:
                return azi + M_2PI
            elif 3.0*M_PI_2 < phi < M_2PI:
                return M_PI - azi
            elif phi == 0.0 or phi == M_2PI:
                return M_PI
            elif phi == M_PI:
                return 0.0
    else:
        if cos(theta) == 0.0:
            critical = M_PI_2
        else:
            critical = atan(tan(theta) * cos(phi))
            if critical < 0.0:
                critical += M_PI

        #with gil:
        #    print('critical rotation = %.8e' % critical)
        #    print('phi = %.8e' % phi)

        if phi < 0.0:
            while phi < 0.0:
                phi += M_2PI
        elif phi > M_2PI:
            while phi > M_2PI:
                phi -= M_2PI

        if THETA <= critical:
            if phi <= M_PI_2:
                return azi
            elif M_PI_2 < phi < M_PI:
                return M_PI - azi
            elif M_PI < phi <= 3.0*M_PI_2:
                return M_PI - azi
            elif 3.0*M_PI_2 < phi < M_2PI:
                return azi + 2.0*M_PI
            elif phi == 0.0 or phi == M_2PI:
                return 0.0
            elif phi == M_PI:
                return M_PI
        elif THETA > critical:
            if phi <= M_PI_2:
                return M_PI - azi
            elif M_PI_2 < phi < M_PI:
                return azi
            elif M_PI < phi <= 3.0*M_PI_2:
                return azi + M_2PI
            elif 3.0*M_PI_2 < phi < M_2PI:
                return M_PI - azi
            elif phi == 0.0 or phi == M_2PI:
                return M_PI
            elif phi == M_PI:
                return 0.0


cdef double eval_phi(double theta, double THETA, double psi) nogil:

    cdef double cos_phi = cos(psi) - cos(THETA) * cos(theta)

    cos_phi /= sin(THETA) * sin(theta)

    if not -1.0 <= cos_phi <= 1.0:
        return -1.0
    else:
        return acos(cos_phi)

cdef double get_interval(double a, double b, double LB, double UB) nogil:

    cdef int a_condition, b_condition

    a_condition = (LB <= a <= UB)
    b_condition = (LB <= b <= UB)

    if a_condition == 0:
        a_condition = 2*<int>(a > UB)
    if b_condition == 0:
        b_condition = 2*<int>(b > UB)

    if a_condition == b_condition:
        if a_condition == 1:
            return b - a
        else:
            return 0.0

    elif b_condition == 1 and a_condition == 0:
        return b - LB

    elif b_condition == 2 and a_condition == 1:
        return UB - a

    elif b_condition == 2 and a_condition == 0:
        return UB - LB

cdef double cell_integrand(double theta, void *params) nogil:

    if theta == 0.0:
        return 0.0

    cdef double epsilon = (<double**> params)[0][0]
    cdef double zeta = (<double**> params)[1][0]
    cdef double colat = (<double**> params)[2][0]
    cdef double cedeRadius = (<double**> params)[3][0]
    cdef double superRadius = (<double**> params)[4][0]
    cdef double superCentreAzi = (<double**> params)[5][0]
    cdef double superCentreColat = (<double**> params)[6][0]
    cdef double phi_a = (<double**>params)[7][0]
    cdef double phi_b = (<double**>params)[8][0]
    cdef int verbose = (<int**>params)[9][0]

    cdef double mu = cos(theta)
    cdef double radius, f

    radius = radiusNormalised(mu, epsilon, zeta)
    f = f_theta(mu, radius, epsilon, zeta)

    cdef double aLB, aUB, cLB, cUB, delta_phi

    aUB = eval_phi(theta, colat, cedeRadius)
    if aUB == -1.0:
        if theta + colat < cedeRadius:
            aUB = M_PI
            aLB = -M_PI
        else:
           return 0.0
    else:
        aLB = -aUB

    if fabs(theta - superCentreColat) > superRadius:
        cLB = cUB = 0.0
    else:
        cUB = eval_phi(theta, superCentreColat, superRadius)
        if cUB == -1.0:
            return 0.0
        else:
            cLB = -cUB
            cUB += superCentreAzi
            cLB += superCentreAzi

        if cUB > M_PI:
            cUB -= M_2PI
        if cLB < -M_PI:
            cLB += M_2PI

    #if verbose:
    #    with gil:
    #        print('theta=%.8e' % theta)
    #        print('cLB=%.8e, cUB=%.8e, aLB=%.8e, aUB=%.8e' % (cLB, cUB, aLB, aUB))
    #        print('phi_a=%.8e, phi_b=%.8e' % (phi_a, phi_b))

    if cUB >= cLB:
        if cLB >= aUB or cUB <= aLB:
            delta_phi = get_interval(phi_a, phi_b, aLB, aUB)
        elif cLB >= aLB and cUB <= aUB:
            delta_phi = get_interval(phi_a, phi_b, aLB, cLB)
            delta_phi += get_interval(phi_a, phi_b, cUB, aUB)
        elif cLB < aLB and aLB < cUB <= aUB:
            delta_phi = get_interval(phi_a, phi_b, cUB, aUB)
        elif aLB <= cLB < aUB and aUB < cUB:
            delta_phi = get_interval(phi_a, phi_b, aLB, cLB)
        elif cLB < aLB and cUB > aUB:
            return 0.0
    else:
        if cUB >= aUB or cLB <= aLB:
            return 0.0
        elif cUB >= aLB and cLB <= aUB:
            delta_phi = get_interval(phi_a, phi_b, cUB, cLB)
        elif cUB < aLB and aLB < cLB <= aUB:
            delta_phi = get_interval(phi_a, phi_b, aLB, cLB)
        elif aLB <= cUB < aUB and aUB < cLB:
            delta_phi = get_interval(phi_a, phi_b, cUB, aUB)
        elif cUB < aLB and cLB > aUB:
            delta_phi = get_interval(phi_a, phi_b, aLB, aUB)

    #if verbose:
    #    with gil:
    #        print('delta_phi=%.8e' % delta_phi)

    return delta_phi * radius * radius * sqrt(1.0 + f * f) * sin(theta)


cdef double theta_given_phi(double phi,
                            double colat,
                            double radius,
                            double interval[]) nogil:

    cdef double x, theta

    x = sin(colat) * sin(phi) / sin(radius)
    if x <= 1.0:
        x = asin(x)
        theta = tan(0.5*(colat - radius)) / sin(0.5*(x - phi))
        theta *= sin(0.5*(x + phi))
        theta = atan(theta)
        if theta < 0.0:
            theta *= -2.0
        else:
            theta *= 2.0

        if not interval[0] <= theta <= interval[1]:
            x = M_PI - x
            theta = tan(0.5*(colat - radius)) / sin(0.5*(x - phi))
            theta *= sin(0.5*(x + phi))
            theta = atan(theta)
            if theta < 0.0:
                theta *= -2.0
            else:
                theta *= 2.0

            if not interval[0] <= theta <= interval[1]:
                return -1.0
            else:
                return theta
        else:
            return theta
    else:
        return -1.0


cdef double integrateCell(double theta_a,
                          double theta_b,
                          double phi_a,
                          double phi_b,
                          double R_eq,
                          double epsilon,
                          double zeta,
                          double colat,
                          double cedeRadius,
                          double superRadius,
                          double superCentreAzi,
                          double superCentreColat,
                          gsl_cq_work *w) nogil:

    cdef void *params[10]
    cdef double integral, error
    cdef int verbose = 0
    cdef size_t nevals, i

    params[0] = &epsilon
    params[1] = &zeta
    params[2] = &colat
    params[3] = &cedeRadius
    params[4] = &superRadius
    params[5] = &superCentreAzi
    params[6] = &superCentreColat
    params[7] = &phi_a
    params[8] = &phi_b
    params[9] = &verbose

    cdef gsl_function F
    F.function = &cell_integrand
    F.params = &params

    gsl_integration_cquad(&F,
                            theta_a,
                            theta_b,
                            0,
                            1.0e-8,
                            w,
                            &integral,
                            &error,
                            &nevals)


    cdef double interval[2]
    interval[0] = theta_a
    interval[1] = theta_b

    cdef double theta_left[2]
    cdef double theta_right[2]
    cdef double phi_top[2]
    cdef double phi_bottom[2]
    cdef double cede_intersect_thetas[2]
    cdef double super_intersect_thetas[2]
    cdef double thetas[4]
    cdef double theta_min, theta_max

    if integral == 0.0:
        #verbose = 1
        # repeat with verbose mode
        #gsl_integration_cquad(&F,
        #                      theta_a,
        #                      theta_b,
        #                      0,
        #                      1.0e-8,
        #                      w,
        #                      &integral,
        #                      &error,
        #                      &nevals)

        # for the cede-region, LHS of element (handedness relative to rotation)
        theta_left[0] = theta_given_phi(phi_a,
                                        colat,
                                        cedeRadius,
                                        interval)
        # for the super-region, LHS of element
        theta_left[1] =  theta_given_phi(phi_a - superCentreAzi,
                                         superCentreColat,
                                         superRadius,
                                         interval)

        # for the cede-region, RHS of element
        theta_right[0] = theta_given_phi(phi_b,
                                         colat,
                                         cedeRadius,
                                         interval)
        # for the super-region, RHS of element
        theta_right[1] = theta_given_phi(phi_b - superCentreAzi,
                                         superCentreColat,
                                         superRadius,
                                         interval)

        # for the cede-region, top of element
        phi_top[0] = eval_phi(theta_a,
                              colat,
                              cedeRadius)
        if phi_top[0] >= 0.0:
            if not phi_a <= phi_top[0] <= phi_b:
                if not phi_a <= -1.0*phi_top[0] <= phi_b:
                    phi_top[0] = -1.0
                elif phi_top[0] < 0.0:
                    phi_top[0] *= -1.0
            elif phi_top[0] < 0.0:
                phi_top[0] *= -1.0

        # for the super-region, top of element
        if fabs(theta_a - superCentreColat) > superRadius:
            phi_top[1] = -1.0
        else:
            phi_top[1] = eval_phi(theta_a,
                                  superCentreColat,
                                  superRadius)
            if phi_top[1] >= 0.0:
                phi_top[1] += superCentreAzi
                if phi_top[1] > M_PI:
                    phi_top[1] -= M_2PI
                if phi_top[1] < -M_PI:
                    phi_top[1] += M_2PI

                if not phi_a <= phi_top[1] <= phi_b:
                    phi_top[1] = -1.0
                elif phi_top[1] < 0.0:
                    phi_top[1] *= -1.0

            if phi_top[1] < 0.0:
                phi_top[1] = eval_phi(theta_a,
                                      superCentreColat,
                                      superRadius)
                if phi_top[1] > 0.0:
                    phi_top[1] *= -1.0
                    phi_top[1] += superCentreAzi
                    if phi_top[1] > M_PI:
                        phi_top[1] -= M_2PI
                    if phi_top[1] < -M_PI:
                        phi_top[1] += M_2PI

                    if not phi_a <= phi_top[1] <= phi_b:
                        phi_top[1] = -1.0
                    elif phi_top[1] < 0.0:
                        phi_top[1] *= -1.0

        # for the cede-region, bottom of element
        phi_bottom[0] = eval_phi(theta_b,
                                 colat,
                                 cedeRadius)
        if phi_bottom[0] >= 0.0:
            if not phi_a <= phi_bottom[0] <= phi_b:
                if not phi_a <= -1.0*phi_bottom[0] <= phi_b:
                    phi_bottom[0] = -1.0
                elif phi_bottom[0] < 0.0:
                    phi_bottom[0] *= -1.0
            elif phi_bottom[0] < 0.0:
                phi_bottom[0] *= -1.0

        # for the super-region, bottom of element
        if fabs(theta_b - superCentreColat) > superRadius:
            phi_bottom[1] = -1.0
        else:
            phi_bottom[1] = eval_phi(theta_b,
                                     superCentreColat,
                                     superRadius)
            if phi_bottom[1] >= 0.0:
                phi_bottom[1] += superCentreAzi
                if phi_bottom[1] > M_PI:
                    phi_bottom[1] -= M_2PI
                if phi_bottom[1] < -M_PI:
                    phi_bottom[1] += M_2PI

                if not phi_a <= phi_bottom[1] <= phi_b:
                    phi_bottom[1] = -1.0
                elif phi_bottom[1] < 0.0:
                    phi_bottom[1] *= -1.0

            if phi_bottom[1] < 0.0:
                phi_bottom[1] = eval_phi(theta_b,
                                         superCentreColat,
                                         superRadius)
                if phi_bottom[1] >= 0.0:
                    phi_bottom[1] *= -1.0
                    phi_bottom[1] += superCentreAzi
                    if phi_bottom[1] > M_PI:
                        phi_bottom[1] -= M_2PI
                    if phi_bottom[1] < -M_PI:
                        phi_bottom[1] += M_2PI

                    if not phi_a <= phi_bottom[1] <= phi_b:
                        phi_bottom[1] = -1.0
                    elif phi_bottom[1] < 0.0:
                        phi_bottom[1] *= -1.0

        # get cede-region intersect element thetas
        if theta_left[0] >= 0.0:
            cede_intersect_thetas[0] = theta_left[0]
            theta_left[0] = -1.0
        elif theta_right[0] >= 0.0:
            cede_intersect_thetas[0] = theta_right[0]
            theta_right[0] = -1.0
        elif phi_top[0] >= 0.0:
            cede_intersect_thetas[0] = theta_a
            phi_top[0] = -1.0
        elif phi_bottom[0] >= 0.0:
            cede_intersect_thetas[0] = theta_b
            phi_bottom[0] = -1.0
        else:
            cede_intersect_thetas[0] = -1.0

        if theta_left[0] >= 0.0:
            cede_intersect_thetas[1] = theta_left[0]
            theta_left[0] = -1.0
        elif theta_right[0] >= 0.0:
            cede_intersect_thetas[1] = theta_right[0]
            theta_right[0] = -1.0
        elif phi_top[0] >= 0.0:
            cede_intersect_thetas[1] = theta_a
            phi_top[0] = -1.0
        elif phi_bottom[0] >= 0.0:
            cede_intersect_thetas[1] = theta_b
            phi_bottom[0] = -1.0
        else:
            cede_intersect_thetas[1] = -1.0

        # get super-region intersect element thetas
        if theta_left[1] >= 0.0:
            super_intersect_thetas[0] = theta_left[1]
            theta_left[1] = -1.0
        elif theta_right[1] >= 0.0:
            super_intersect_thetas[0] = theta_right[1]
            theta_right[1] = -1.0
        elif phi_top[1] >= 0.0:
            super_intersect_thetas[0] = theta_a
            phi_top[1] = -1.0
        elif phi_bottom[1] >= 0.0:
            super_intersect_thetas[0] = theta_b
            phi_bottom[1] = -1.0
        else:
            super_intersect_thetas[0] = -1.0

        if theta_left[1] >= 0.0:
            super_intersect_thetas[1] = theta_left[1]
            theta_left[1] = -1.0
        elif theta_right[1] >= 0.0:
            super_intersect_thetas[1] = theta_right[1]
            theta_right[1] = -1.0
        elif phi_top[1] >= 0.0:
            super_intersect_thetas[1] = theta_a
            phi_top[1] = -1.0
        elif phi_bottom[1] >= 0.0:
            super_intersect_thetas[1] = theta_b
            phi_bottom[1] = -1.0
        else:
            super_intersect_thetas[1] = -1.0

        theta_min = theta_b
        theta_max = theta_a

        thetas[0] = cede_intersect_thetas[0]
        thetas[1] = cede_intersect_thetas[1]
        thetas[2] = super_intersect_thetas[0]
        thetas[3] = super_intersect_thetas[1]

        for i in range(4):
            if thetas[i] >= 0.0:
                if thetas[i] < theta_min:
                    theta_min = thetas[i]
                if thetas[i] > theta_max:
                    theta_max = thetas[i]

        #with gil:
        #    if integral == 0.0:
        #        print('')
        #        print('cell integral is zero, with error %.8e and nevals %i' % (error, nevals))
        #        print('cede-region', cede_intersect_thetas[0], cede_intersect_thetas[1])
        #        print('super-region', super_intersect_thetas[0], super_intersect_thetas[1])
        #        print('temps', thetas[0], thetas[1], thetas[2], thetas[3])
        #        print('theta_min = %.8e, theta_max = %.8e' % (theta_min, theta_max))
        #        print('theta_a = %.8e, theta_b = %.8e' % (theta_a, theta_b))
        #        print('')

        if theta_min < theta_max:
            gsl_integration_cquad(&F,
                                    theta_min,
                                    theta_max,
                                    0,
                                    1.0e-8,
                                    w,
                                    &integral,
                                    &error,
                                    &nevals)

    #with gil:
    #    if integral == 0.0:
    #        print('cell integral is zero, with error %.8e and nevals %i' % (error, nevals))

    return R_eq * R_eq * integral


cdef double spot_integrand(double theta, void *params) nogil:

    if theta == 0.0:
        return 0.0

    cdef double epsilon = (<double*>params)[0]
    cdef double zeta = (<double*>params)[1]
    cdef double colat = (<double*>params)[2]
    cdef double cedeRadius = (<double*>params)[3]
    cdef double superRadius = (<double*>params)[4]
    cdef double superCentreAzi = (<double*>params)[5]
    cdef double superCentreColat = (<double*>params)[6]
    cdef int Lorentz = <int>((<double*>params)[7])
    cdef double R_eq = (<double*>params)[8]
    cdef double Omega = (<double*>params)[9]
    cdef double r_s_over_r

    cdef double mu = cos(theta)
    cdef double radius, f

    radius = radiusNormalised(mu, epsilon, zeta)
    f = f_theta(mu, radius, epsilon, zeta)

    cdef double aLB, aUB, cLB, cUB, delta_phi

    aUB = eval_phi(theta, colat, cedeRadius)
    if aUB == -1.0:
        if theta + colat < cedeRadius:
            aUB = M_PI
            aLB = -M_PI
        else:
           return 0.0
    else:
        aLB = -aUB

    if fabs(theta - superCentreColat) > superRadius:
        cLB = cUB = 0.0
    else:
        cUB = eval_phi(theta, superCentreColat, superRadius)
        if cUB == -1.0:
            return 0.0
        else:
            cLB = -cUB
            cUB += superCentreAzi
            cLB += superCentreAzi

        if cUB > M_PI:
            cUB -= M_2PI
        if cLB < -M_PI:
            cLB += M_2PI

    if cUB >= cLB:
        if cLB >= aUB or cUB <= aLB:
            delta_phi = -aLB + aUB
        elif cLB >= aLB and cUB <= aUB:
            delta_phi = -aLB + cLB
            delta_phi += -cUB + aUB
        elif cLB < aLB and aLB < cUB <= aUB:
            delta_phi = -cUB + aUB
        elif aLB <= cLB < aUB and aUB < cUB:
            delta_phi = -aLB + cLB
        elif cLB < aLB and cUB > aUB:
            return 0.0
    else:
        if cUB >= aUB or cLB <= aLB:
            return 0.0
        elif cUB >= aLB and cLB <= aUB:
            delta_phi = -cUB + cLB
        elif cUB < aLB and aLB < cLB <= aUB:
            delta_phi = -aLB + cLB
        elif aLB <= cUB < aUB and aUB < cLB:
            delta_phi = -cUB + aUB
        elif cUB < aLB and cLB > aUB:
            delta_phi = -aLB + aUB

    if Lorentz == 1:
        r_s_over_r = 2.0 * zeta / radius
        beta = R_eq * radius * Omega * sin(theta) / (_c * sqrt(1.0 - r_s_over_r))
        eval_integrand = 1.0 / sqrt(1.0 - beta*beta)
    else:
        eval_integrand = 1.0

    eval_integrand *= delta_phi * radius * radius * sqrt(1.0 + f * f) * sin(theta)

    return eval_integrand


cdef double integrateSpot(double theta_a,
                          double theta_b,
                          double R_eq,
                          double epsilon,
                          double zeta,
                          double colat,
                          double cedeRadius,
                          double superRadius,
                          double superCentreColat,
                          double superCentreAzi,
                          int Lorentz,
                          double Omega,
                          gsl_cq_work *w) nogil:


    cdef double params[10]
    cdef double integral, error
    cdef size_t nevals

    params[0] = epsilon
    params[1] = zeta
    params[2] = colat
    params[3] = cedeRadius
    params[4] = superRadius
    params[5] = superCentreAzi
    params[6] = superCentreColat
    params[7] = <double> Lorentz
    params[8] = R_eq
    params[9] = Omega

    cdef gsl_function F
    F.function = &spot_integrand
    F.params = &params

    gsl_integration_cquad(&F,
                            theta_a,
                            theta_b,
                            0,
                            1.0e-8,
                            w,
                            &integral,
                            &error,
                            &nevals)

    return R_eq * R_eq * integral

def eval_cedeCentreCoords(double superColatitude,
                          double cedeOffset,
                          double cedeOffsetAzi):

    cdef double cedeColatitude, cedeAzimuth

    if cedeOffset > 0.0:
        cedeColatitude = eval_psi(cedeOffset, cedeOffsetAzi, -superColatitude)
        if cedeColatitude > 0.0:
            cedeAzimuth = eval_cedeAzi(cedeOffset, cedeOffsetAzi, cedeColatitude, -superColatitude)
        else:
            cedeAzimuth = 0.0

        if cedeAzimuth > M_PI:
            while cedeAzimuth > M_PI:
                cedeAzimuth -= M_2PI
        elif cedeAzimuth < -M_PI:
            while cedeAzimuth < -M_PI:
                cedeAzimuth += M_2PI
    else:
        cedeColatitude = superColatitude
        cedeAzimuth = 0.0

    #print('cedeColatitude = %.8e' % cedeColatitude)
    #print('cedeAzimuth = %.8e' % cedeAzimuth)

    return cedeColatitude, cedeAzimuth

def allocate_cells(size_t num_cells,
                   size_t min_sqrt_num_cells,
                   size_t max_sqrt_num_cells,
                   double mass,
                   double r_s,
                   double R_eq,
                   double Omega,
                   double zeta,
                   double epsilon,
                   double superRadius,
                   double cedeRadius,
                   double superColatitude,
                   double cedeColatitude,
                   double superAzimuth,
                   double holeColatitude,
                   double holeRadius,
                   double holeAzimuth,
                   fast_components):

    cdef:
        double super_numCell, super_sqrt_numCell, super_cellArea, super_area
        double super_boundary_phi, superColatitude_cpy
        double super_colat_lims[2]
        double cede_numCell, cede_sqrt_numCell, cede_cellArea, cede_area
        double cede_boundary_phi
        double cede_colat_lims[2]
        double f, y
        double fast_super, fast_cede
        cdef gsl_cq_work *w = gsl_integration_cquad_workspace_alloc(100)

    fast_super = 1.0
    fast_cede = 1.0
    if fast_components is not None and len(fast_components) > 1:
        try:
            fast_super, fast_cede = fast_components
        except TypeError:
            fast_super = 1.0
            fast_cede = 1.0

        if fast_super == 0.0 and fast_cede == 0.0: # invisible (at this res)
            fast_super = 1.0 # just assume equal and base on relative area
            fast_cede = 1.0

    superColatitude_cpy = superColatitude
    # first integrate area of default superseding region
    if superColatitude - superRadius < 0.0 or superColatitude + superRadius > M_PI:
        super_colat_lims[0] = 0.0
        if superColatitude + superRadius > M_PI:
            super_colat_lims[1] = M_PI - superColatitude + superRadius
            superColatitude = M_PI - superColatitude
            holeColatitude = M_PI - holeColatitude
        else:
            super_colat_lims[1] = superColatitude + superRadius

        super_boundary_phi = M_PI
        super_cellArea = M_2PI
    else:
        super_colat_lims[0] = superColatitude - superRadius
        super_colat_lims[1] = superColatitude + superRadius

        # Right-spherical triangle
        super_boundary_phi = asin(sin(superRadius) / sin(superColatitude))
        super_cellArea = 2.0 * super_boundary_phi

    super_cellArea *= integrateArea(super_colat_lims[0],
                                      super_colat_lims[1],
                                      R_eq,
                                      epsilon,
                                      zeta,
                                      0,
                                      w)

    super_area = integrateSpot(super_colat_lims[0],
                               super_colat_lims[1],
                               R_eq,
                               epsilon,
                               zeta,
                               superColatitude,
                               superRadius,
                               holeRadius,
                               holeColatitude,
                               holeAzimuth,
                               0,
                               Omega,
                               w)
    #print(super_area)
    #if super_area == 0.0:
    #    super_area = integrateSpot(super_colat_lims[0],
    #                               super_colat_lims[0] + 0.5*superRadius,#(super_colat_lims[1] + super_colat_lims[0])/2.0,
    #                               R_eq,
    #                               epsilon,
    #                               zeta,
    #                               superColatitude,
    #                               superRadius,
    #                               holeRadius,
    #                               holeColatitude,
    #                               holeAzimuth,
    #                               0,
    #                               Omega,
    #                               w)
    #print(super_area, sqrt(<double>num_cells * 1000.0))
    #print(super_colat_lims[0], super_colat_lims[1])
    #print(R_eq, epsilon, zeta, Omega)
    #print(superColatitude, superRadius)
    #print(holeRadius, holeColatitude, holeAzimuth)

    if cedeRadius == 0.0: # no ceding region
        if super_area == 0.0: # Gaussian integral did not resolve
            super_area = super_cellArea / 1000.0
        super_sqrt_numCell = ceil(sqrt((<double>num_cells) * super_cellArea / super_area))
        #print(super_sqrt_numCell)
        if super_sqrt_numCell < min_sqrt_num_cells:
            super_sqrt_numCell = min_sqrt_num_cells
        elif super_sqrt_numCell > max_sqrt_num_cells:
            super_sqrt_numCell = max_sqrt_num_cells

        if <size_t>(super_sqrt_numCell)%2 != 0:
            super_sqrt_numCell += 1.0

        super_numCell = super_sqrt_numCell * super_sqrt_numCell

        gsl_integration_cquad_workspace_free(w)

        return (super_sqrt_numCell,
                super_numCell,
                0, 0)

    else: # there is a ceding region
        if cedeColatitude - cedeRadius < 0.0 or cedeColatitude + cedeRadius > M_PI:
            cede_colat_lims[0] = 0.0
            if cedeColatitude + cedeRadius > M_PI:
                cede_colat_lims[1] = M_PI - cedeColatitude + cedeRadius
                cedeColatitude = M_PI - cedeColatitude
                superColatitude = M_PI - superColatitude_cpy
            else:
                cede_colat_lims[1] = cedeColatitude + cedeRadius
                superColatitude = superColatitude_cpy

            cede_boundary_phi = M_PI
            cede_cellArea = M_2PI
        else:
            cede_colat_lims[0] = cedeColatitude - cedeRadius
            cede_colat_lims[1] = cedeColatitude + cedeRadius
            superColatitude = superColatitude_cpy

            # Right-spherical triangle
            cede_boundary_phi = asin(sin(cedeRadius) / sin(cedeColatitude))
            cede_cellArea = 2.0 * cede_boundary_phi

        cede_cellArea *= integrateArea(cede_colat_lims[0],
                                          cede_colat_lims[1],
                                          R_eq,
                                          epsilon,
                                          zeta,
                                          0,
                                          w)

        #print(cede_colat_lims[0], cede_colat_lims[1])
        #print(R_eq, epsilon, zeta, Omega)
        #print(cedeColatitude, cedeRadius)
        #print(superRadius, superColatitude, superAzimuth)

        cede_area = integrateSpot(cede_colat_lims[0],
                                     cede_colat_lims[1],
                                     R_eq,
                                     epsilon,
                                     zeta,
                                     cedeColatitude,
                                     cedeRadius,
                                     superRadius,
                                     superColatitude,
                                     superAzimuth,
                                     0,
                                     Omega,
                                     w)

        if super_area == 0.0: # Gausian integral did not resolve
            super_area = super_cellArea / 1000.0
        if cede_area == 0.0: # Gausian integral did not resolve
            cede_area = cede_cellArea / 1000.0

        if fast_super > 0.0:
            f = fast_super / (fast_super + fast_cede)
            #print('Fast super signal/fast cede signal: %.8e' % f)

            y = ((1.0 - f)/f)*(cede_area/super_area) - 1.0

            if y == 0.0:
                super_numCell = 0.5 * <double>num_cells
                cede_numCell = 0.5 * <double>num_cells
            else:
                super_numCell = <double>num_cells * ((sqrt(1.0 + y) - 1.0)/y)
                cede_numCell = <double>num_cells - super_numCell
        else:
            super_numCell = 0.0
            cede_numCell = <double>num_cells

        #print(super_numCell, super_cellArea, super_area)
        super_sqrt_numCell = ceil(sqrt(super_numCell * super_cellArea / super_area))
        #print(cede_numCell, cede_cellArea, cede_area)
        cede_sqrt_numCell = ceil(sqrt(cede_numCell * cede_cellArea / cede_area))

        if super_sqrt_numCell < min_sqrt_num_cells:
            super_sqrt_numCell = min_sqrt_num_cells
        elif super_sqrt_numCell > max_sqrt_num_cells:
            super_sqrt_numCell = max_sqrt_num_cells

        if cede_sqrt_numCell < min_sqrt_num_cells:
            cede_sqrt_numCell = min_sqrt_num_cells
        elif cede_sqrt_numCell > max_sqrt_num_cells:
            cede_sqrt_numCell = max_sqrt_num_cells

        if <size_t>(super_sqrt_numCell)%2 != 0:
            super_sqrt_numCell += 1.0

        if <size_t>(cede_sqrt_numCell)%2 != 0:
            cede_sqrt_numCell += 1.0

        super_numCell = super_sqrt_numCell * super_sqrt_numCell
        cede_numCell = cede_sqrt_numCell * cede_sqrt_numCell

        gsl_integration_cquad_workspace_free(w)

        return (super_sqrt_numCell, super_numCell,
                cede_sqrt_numCell, cede_numCell)
