#cython: cdivision=True
#cython: boundscheck=False
#cython: nonecheck=False
#cython: wraparound=False

from libc.math cimport log

cdef extern from "math.h":
    long double logl(long double x) noexcept nogil

## Expressions are linked to equations in relevant quasi-Kerr papers:
##   Psaltis & Johannsen 2012, ApJ 745, 1 (PJ12)
##   Bauböck et al. 2012, ApJ 753, 175 (B12)
## Both build on Glampedakis & Babak 2006, CQG 23, 4167 (GB06).
## All equation numbers refer to the published journal versions.
##
## The quasi-Kerr metric in Boyer-Lindquist coordinates is given in
## PJ12 Equations (1)-(5), B12 Equations (2)-(6), and GB06 Equations (2.13)-(2.16).
##
## Conventions used throughout this module:
##   - SI units with G = c = 1 absorbed into lengths: r_s = 2GM/c^2 is the
##     Schwarzschild radius, so M in the papers corresponds to r_s/2 here.
##   - a = J/(Mc) is the spin parameter in metres (not dimensionless).
##   - kappa is the dimensionless quadrupole deviation from Kerr: epsilon in
##     GB06 and PJ12, eta*(a/M)^2 in B12. The quadrupole moment is then
##     Q = -M(a^2 + kappa*M^2) (PJ12 Eq. 6; GB06 Eq. 2.3, Q = Q^K - epsilon*M^3).
##   - Index labels 1, 2, 3, 4 denote t, r, theta, phi.
##   - Callers pass in the precomputed helper quantities
##       Sigma      = r^2 + a^2 cos^2(theta)        (PJ12 Eq. 4, B12 Eq. 5)
##       Delta      = r^2 - r_s r + a^2             (PJ12 Eq. 3, B12 Eq. 4)
##       func_theta = 1 - 3 cos^2(theta)            (angular factor in h^uv)
##       func_r     = 1 - r_s/r                     (radial factor in h^uv)
##       F1_r, F2_r = F1(r), F2(r) defined below
##
## The F1, F2 terms contain log(r/(r - r_s)), which diverges at r = r_s. PJ12
## (Section 2) note the quasi-Kerr metric has unphysical pathologies at
## r < ~2.6M for a/M < ~0.4; neutron-star surfaces lie well outside this.

# Large-r series for F1, F2 and their derivatives
## WHY: the closed forms below (GB06 Eqs. A.9, A.10) subtract two terms that
## are each ~ (15/4)(r/r_s) to leave a result ~ (r_s/r)^3. That cancellation
## costs ~4 digits per decade in r, and evaluating log(r/(r - r_s)) for an
## argument this close to 1 costs about one more, so the closed forms
## lose ~5 digits per decade in r: in double precision the relative error is
## ~1e-12 at 5 r_s, ~1e-7 at 50 r_s, and at the image plane (~3.6e4 r_s)
## F1 comes out ~1e8 times too large with the wrong sign. Times kappa, that
## noise makes the adaptive integrator shrink its step to ~1 m and run out of
## steps (or stall) unless 1 - 3cos^2(theta) ~ 0.
##
## FIX: for r > _SERIES_SWITCH * r_s use the Taylor series in u = r_s/r,
##   F(u) = sum_{n=3}^{40} a_n u^n,   r_s dF/dr = -sum_{n=3}^{40} n a_n u^(n+1),
## which has no cancellation. The a_n are exact rationals obtained by
## expanding GB06 Eqs. (A.9)-(A.10) with log(r/(r - r_s)) = sum_k u^k / k;
## the leading terms are F1 = -u^3/8 - ..., F2 = +u^3/8 + ..., i.e.
## F1 -> -(M/r)^3 as required for Q = Q^K - kappa M^3 (GB06 Eq. 2.3).
## Checked against 60-digit arithmetic: maximum relative error ~5e-16 for
## F1, F2 and both derivatives over 4 r_s < r < 1e6 r_s. Below 4 r_s the
## closed forms are kept; their double-precision error there is at most ~1e-12.
## The series is evaluated in Horner form (as flat statements); coefficients
## are generated, not typed by hand.

cdef double _SERIES_SWITCH = 4.0

cdef inline double _F1_series(double u) noexcept nogil:
    ## Horner evaluation, written as flat statements (deeply nested
    ## parentheses exceed the Cython parser's recursion limit on Python 3.14).
    cdef double p = -0.31032229965156793
    p = -0.31021341463414637 + u*p
    p = -0.31009615384615385 + u*p
    p = -0.3099696356275304 + u*p
    p = -0.30983285917496445 + u*p
    p = -0.3096846846846847 + u*p
    p = -0.30952380952380953 + u*p
    p = -0.30934873949579833 + u*p
    p = -0.3091577540106952 + u*p
    p = -0.30894886363636365 + u*p
    p = -0.3087197580645161 + u*p
    p = -0.3084677419354839 + u*p
    p = -0.3081896551724138 + u*p
    p = -0.3078817733990148 + u*p
    p = -0.30753968253968256 + u*p
    p = -0.3071581196581197 + u*p
    p = -0.30673076923076925 + u*p
    p = -0.30625 + u*p
    p = -0.30570652173913043 + u*p
    p = -0.3050889328063241 + u*p
    p = -0.30438311688311687 + u*p
    p = -0.30357142857142855 + u*p
    p = -0.3026315789473684 + u*p
    p = -0.30153508771929827 + u*p
    p = -0.3002450980392157 + u*p
    p = -0.29871323529411764 + u*p
    p = -0.296875 + u*p
    p = -0.29464285714285715 + u*p
    p = -0.2918956043956044 + u*p
    p = -0.28846153846153844 + u*p
    p = -0.2840909090909091 + u*p
    p = -0.2784090909090909 + u*p
    p = -0.2708333333333333 + u*p
    p = -0.2604166666666667 + u*p
    p = -0.24553571428571427 + u*p
    p = -0.22321428571428573 + u*p
    p = -0.1875 + u*p
    p = -0.125 + u*p
    return u*u*u*p

cdef inline double _F2_series(double u) noexcept nogil:
    ## Horner evaluation, written as flat statements (deeply nested
    ## parentheses exceed the Cython parser's recursion limit on Python 3.14).
    cdef double p = 0.04241071428571429
    p = 0.04338649155722327 + u*p
    p = 0.044407894736842105 + u*p
    p = 0.04547817047817048 + u*p
    p = 0.046600877192982455 + u*p
    p = 0.04777992277992278 + u*p
    p = 0.049019607843137254 + u*p
    p = 0.05032467532467533 + u*p
    p = 0.05170036764705882 + u*p
    p = 0.0531524926686217 + u*p
    p = 0.0546875 + u*p
    p = 0.05631256952169077 + u*p
    p = 0.05803571428571429 + u*p
    p = 0.05986590038314176 + u*p
    p = 0.061813186813186816 + u*p
    p = 0.06388888888888888 + u*p
    p = 0.06610576923076923 + u*p
    p = 0.06847826086956521 + u*p
    p = 0.07102272727272728 + u*p
    p = 0.07375776397515528 + u*p
    p = 0.07670454545454546 + u*p
    p = 0.07988721804511278 + u*p
    p = 0.08333333333333333 + u*p
    p = 0.08707430340557276 + u*p
    p = 0.09114583333333333 + u*p
    p = 0.09558823529411764 + u*p
    p = 0.10044642857142858 + u*p
    p = 0.10576923076923077 + u*p
    p = 0.11160714285714286 + u*p
    p = 0.11800699300699301 + u*p
    p = 0.125 + u*p
    p = 0.13257575757575757 + u*p
    p = 0.140625 + u*p
    p = 0.1488095238095238 + u*p
    p = 0.15625 + u*p
    p = 0.16071428571428573 + u*p
    p = 0.15625 + u*p
    p = 0.125 + u*p
    return u*u*u*p

cdef inline double _rs_dF1dr_series(double u) noexcept nogil:
    ## Horner evaluation, written as flat statements (deeply nested
    ## parentheses exceed the Cython parser's recursion limit on Python 3.14).
    cdef double p = -12.412891986062718
    p = -12.098323170731707 + u*p
    p = -11.783653846153847 + u*p
    p = -11.468876518218623 + u*p
    p = -11.153982930298719 + u*p
    p = -10.838963963963964 + u*p
    p = -10.523809523809524 + u*p
    p = -10.208508403361344 + u*p
    p = -9.893048128342246 + u*p
    p = -9.577414772727273 + u*p
    p = -9.261592741935484 + u*p
    p = -8.945564516129032 + u*p
    p = -8.629310344827585 + u*p
    p = -8.312807881773399 + u*p
    p = -7.996031746031746 + u*p
    p = -7.6789529914529915 + u*p
    p = -7.361538461538461 + u*p
    p = -7.04375 + u*p
    p = -6.725543478260869 + u*p
    p = -6.406867588932807 + u*p
    p = -6.087662337662338 + u*p
    p = -5.767857142857143 + u*p
    p = -5.447368421052632 + u*p
    p = -5.12609649122807 + u*p
    p = -4.803921568627451 + u*p
    p = -4.4806985294117645 + u*p
    p = -4.15625 + u*p
    p = -3.830357142857143 + u*p
    p = -3.5027472527472527 + u*p
    p = -3.173076923076923 + u*p
    p = -2.840909090909091 + u*p
    p = -2.5056818181818183 + u*p
    p = -2.1666666666666665 + u*p
    p = -1.8229166666666667 + u*p
    p = -1.4732142857142858 + u*p
    p = -1.1160714285714286 + u*p
    p = -0.75 + u*p
    p = -0.375 + u*p
    return -u*u*u*u*p

cdef inline double _rs_dF2dr_series(double u) noexcept nogil:
    ## Horner evaluation, written as flat statements (deeply nested
    ## parentheses exceed the Cython parser's recursion limit on Python 3.14).
    cdef double p = 1.6964285714285714
    p = 1.6920731707317074 + u*p
    p = 1.6875 + u*p
    p = 1.6826923076923077 + u*p
    p = 1.6776315789473684 + u*p
    p = 1.6722972972972974 + u*p
    p = 1.6666666666666667 + u*p
    p = 1.6607142857142858 + u*p
    p = 1.6544117647058822 + u*p
    p = 1.6477272727272727 + u*p
    p = 1.640625 + u*p
    p = 1.6330645161290323 + u*p
    p = 1.625 + u*p
    p = 1.6163793103448276 + u*p
    p = 1.6071428571428572 + u*p
    p = 1.5972222222222223 + u*p
    p = 1.5865384615384615 + u*p
    p = 1.575 + u*p
    p = 1.5625 + u*p
    p = 1.548913043478261 + u*p
    p = 1.5340909090909092 + u*p
    p = 1.5178571428571428 + u*p
    p = 1.5 + u*p
    p = 1.480263157894737 + u*p
    p = 1.4583333333333333 + u*p
    p = 1.4338235294117647 + u*p
    p = 1.40625 + u*p
    p = 1.375 + u*p
    p = 1.3392857142857142 + u*p
    p = 1.2980769230769231 + u*p
    p = 1.25 + u*p
    p = 1.1931818181818181 + u*p
    p = 1.125 + u*p
    p = 1.0416666666666667 + u*p
    p = 0.9375 + u*p
    p = 0.8035714285714286 + u*p
    p = 0.625 + u*p
    p = 0.375 + u*p
    return -u*u*u*u*p

cdef inline long double _F1_series_l(long double u) noexcept nogil:
    ## Horner evaluation, written as flat statements (deeply nested
    ## parentheses exceed the Cython parser's recursion limit on Python 3.14).
    cdef long double p = (<long double>-1425/<long double>4592)
    p = (<long double>-407/<long double>1312) + u*p
    p = (<long double>-129/<long double>416) + u*p
    p = (<long double>-1225/<long double>3952) + u*p
    p = (<long double>-3485/<long double>11248) + u*p
    p = (<long double>-275/<long double>888) + u*p
    p = (<long double>-13/<long double>42) + u*p
    p = (<long double>-589/<long double>1904) + u*p
    p = (<long double>-925/<long double>2992) + u*p
    p = (<long double>-435/<long double>1408) + u*p
    p = (<long double>-1225/<long double>3968) + u*p
    p = (<long double>-153/<long double>496) + u*p
    p = (<long double>-143/<long double>464) + u*p
    p = (<long double>-125/<long double>406) + u*p
    p = (<long double>-155/<long double>504) + u*p
    p = (<long double>-575/<long double>1872) + u*p
    p = (<long double>-319/<long double>1040) + u*p
    p = (<long double>-49/<long double>160) + u*p
    p = (<long double>-225/<long double>736) + u*p
    p = (<long double>-1235/<long double>4048) + u*p
    p = (<long double>-375/<long double>1232) + u*p
    p = (<long double>-17/<long double>56) + u*p
    p = (<long double>-23/<long double>76) + u*p
    p = (<long double>-275/<long double>912) + u*p
    p = (<long double>-245/<long double>816) + u*p
    p = (<long double>-325/<long double>1088) + u*p
    p = (<long double>-19/<long double>64) + u*p
    p = (<long double>-33/<long double>112) + u*p
    p = (<long double>-425/<long double>1456) + u*p
    p = (<long double>-15/<long double>52) + u*p
    p = (<long double>-25/<long double>88) + u*p
    p = (<long double>-49/<long double>176) + u*p
    p = (<long double>-13/<long double>48) + u*p
    p = (<long double>-25/<long double>96) + u*p
    p = (<long double>-55/<long double>224) + u*p
    p = (<long double>-25/<long double>112) + u*p
    p = (<long double>-3/<long double>16) + u*p
    p = (<long double>-1/<long double>8) + u*p
    return u*u*u*p

cdef inline long double _F2_series_l(long double u) noexcept nogil:
    ## Horner evaluation, written as flat statements (deeply nested
    ## parentheses exceed the Cython parser's recursion limit on Python 3.14).
    cdef long double p = (<long double>19/<long double>448)
    p = (<long double>185/<long double>4264) + u*p
    p = (<long double>27/<long double>608) + u*p
    p = (<long double>175/<long double>3848) + u*p
    p = (<long double>85/<long double>1824) + u*p
    p = (<long double>99/<long double>2072) + u*p
    p = (<long double>5/<long double>102) + u*p
    p = (<long double>31/<long double>616) + u*p
    p = (<long double>225/<long double>4352) + u*p
    p = (<long double>145/<long double>2728) + u*p
    p = (<long double>7/<long double>128) + u*p
    p = (<long double>405/<long double>7192) + u*p
    p = (<long double>13/<long double>224) + u*p
    p = (<long double>125/<long double>2088) + u*p
    p = (<long double>45/<long double>728) + u*p
    p = (<long double>23/<long double>360) + u*p
    p = (<long double>55/<long double>832) + u*p
    p = (<long double>63/<long double>920) + u*p
    p = (<long double>25/<long double>352) + u*p
    p = (<long double>95/<long double>1288) + u*p
    p = (<long double>27/<long double>352) + u*p
    p = (<long double>85/<long double>1064) + u*p
    p = (<long double>1/<long double>12) + u*p
    p = (<long double>225/<long double>2584) + u*p
    p = (<long double>35/<long double>384) + u*p
    p = (<long double>13/<long double>136) + u*p
    p = (<long double>45/<long double>448) + u*p
    p = (<long double>11/<long double>104) + u*p
    p = (<long double>25/<long double>224) + u*p
    p = (<long double>135/<long double>1144) + u*p
    p = (<long double>1/<long double>8) + u*p
    p = (<long double>35/<long double>264) + u*p
    p = (<long double>9/<long double>64) + u*p
    p = (<long double>25/<long double>168) + u*p
    p = (<long double>5/<long double>32) + u*p
    p = (<long double>9/<long double>56) + u*p
    p = (<long double>5/<long double>32) + u*p
    p = (<long double>1/<long double>8) + u*p
    return u*u*u*p

cdef inline long double _rs_dF1dr_series_l(long double u) noexcept nogil:
    ## Horner evaluation, written as flat statements (deeply nested
    ## parentheses exceed the Cython parser's recursion limit on Python 3.14).
    cdef long double p = (<long double>-7125/<long double>574)
    p = (<long double>-15873/<long double>1312) + u*p
    p = (<long double>-2451/<long double>208) + u*p
    p = (<long double>-45325/<long double>3952) + u*p
    p = (<long double>-31365/<long double>2812) + u*p
    p = (<long double>-9625/<long double>888) + u*p
    p = (<long double>-221/<long double>21) + u*p
    p = (<long double>-19437/<long double>1904) + u*p
    p = (<long double>-1850/<long double>187) + u*p
    p = (<long double>-13485/<long double>1408) + u*p
    p = (<long double>-18375/<long double>1984) + u*p
    p = (<long double>-4437/<long double>496) + u*p
    p = (<long double>-1001/<long double>116) + u*p
    p = (<long double>-3375/<long double>406) + u*p
    p = (<long double>-2015/<long double>252) + u*p
    p = (<long double>-14375/<long double>1872) + u*p
    p = (<long double>-957/<long double>130) + u*p
    p = (<long double>-1127/<long double>160) + u*p
    p = (<long double>-2475/<long double>368) + u*p
    p = (<long double>-25935/<long double>4048) + u*p
    p = (<long double>-1875/<long double>308) + u*p
    p = (<long double>-323/<long double>56) + u*p
    p = (<long double>-207/<long double>38) + u*p
    p = (<long double>-4675/<long double>912) + u*p
    p = (<long double>-245/<long double>51) + u*p
    p = (<long double>-4875/<long double>1088) + u*p
    p = (<long double>-133/<long double>32) + u*p
    p = (<long double>-429/<long double>112) + u*p
    p = (<long double>-1275/<long double>364) + u*p
    p = (<long double>-165/<long double>52) + u*p
    p = (<long double>-125/<long double>44) + u*p
    p = (<long double>-441/<long double>176) + u*p
    p = (<long double>-13/<long double>6) + u*p
    p = (<long double>-175/<long double>96) + u*p
    p = (<long double>-165/<long double>112) + u*p
    p = (<long double>-125/<long double>112) + u*p
    p = (<long double>-3/<long double>4) + u*p
    p = (<long double>-3/<long double>8) + u*p
    return -u*u*u*u*p

cdef inline long double _rs_dF2dr_series_l(long double u) noexcept nogil:
    ## Horner evaluation, written as flat statements (deeply nested
    ## parentheses exceed the Cython parser's recursion limit on Python 3.14).
    cdef long double p = (<long double>95/<long double>56)
    p = (<long double>555/<long double>328) + u*p
    p = (<long double>27/<long double>16) + u*p
    p = (<long double>175/<long double>104) + u*p
    p = (<long double>255/<long double>152) + u*p
    p = (<long double>495/<long double>296) + u*p
    p = (<long double>5/<long double>3) + u*p
    p = (<long double>93/<long double>56) + u*p
    p = (<long double>225/<long double>136) + u*p
    p = (<long double>145/<long double>88) + u*p
    p = (<long double>105/<long double>64) + u*p
    p = (<long double>405/<long double>248) + u*p
    p = (<long double>13/<long double>8) + u*p
    p = (<long double>375/<long double>232) + u*p
    p = (<long double>45/<long double>28) + u*p
    p = (<long double>115/<long double>72) + u*p
    p = (<long double>165/<long double>104) + u*p
    p = (<long double>63/<long double>40) + u*p
    p = (<long double>25/<long double>16) + u*p
    p = (<long double>285/<long double>184) + u*p
    p = (<long double>135/<long double>88) + u*p
    p = (<long double>85/<long double>56) + u*p
    p = (<long double>3/<long double>2) + u*p
    p = (<long double>225/<long double>152) + u*p
    p = (<long double>35/<long double>24) + u*p
    p = (<long double>195/<long double>136) + u*p
    p = (<long double>45/<long double>32) + u*p
    p = (<long double>11/<long double>8) + u*p
    p = (<long double>75/<long double>56) + u*p
    p = (<long double>135/<long double>104) + u*p
    p = (<long double>5/<long double>4) + u*p
    p = (<long double>105/<long double>88) + u*p
    p = (<long double>9/<long double>8) + u*p
    p = (<long double>25/<long double>24) + u*p
    p = (<long double>15/<long double>16) + u*p
    p = (<long double>45/<long double>56) + u*p
    p = (<long double>5/<long double>8) + u*p
    p = (<long double>3/<long double>8) + u*p
    return -u*u*u*u*p

# Quadrupole moment correction F(r) functions

## F1(r): GB06 Eq. (A.9), the calligraphic F_1 in h^uv, with M = r_s/2.
## For r > 4 r_s the large-r series above is used instead (see WHY above).
## (Checked symbolically: identical to GB06 after substituting r_s = 2M.)
## Large-r behaviour: F1 -> -(M/r)^3, which is what makes the kappa term
## add -kappa*M^3 to the quadrupole moment.
cdef double F1(double r, double r_s) noexcept nogil:

    if r > _SERIES_SWITCH * r_s:
        return _F1_series(r_s / r)
    return -5.0*(2.0*r - r_s) * (r_s*r_s + 6.0*r_s*r - 6.0*r*r) / (16.0*r_s*r*(r - r_s)) - 15.0*r*(r - r_s)*log(r/(r - r_s)) / (4.0*r_s*r_s)


## F2(r): GB06 Eq. (A.10), the calligraphic F_2 in h^uv, with M = r_s/2.
## For r > 4 r_s the large-r series above is used instead.
## (Checked symbolically: identical to GB06 after substituting r_s = 2M.)
cdef double F2(double r, double r_s) noexcept nogil:

    if r > _SERIES_SWITCH * r_s:
        return _F2_series(r_s / r)
    return 5.0*(r_s*r_s - 3.0*r_s*r - 6.0*r*r) / (8.0*r_s*r) + 15.0*(2.0*r*r - r_s*r_s)*log(r/(r - r_s)) / (8.0*r_s*r_s)


## dF1/dr: analytic derivative of GB06 Eq. (A.9) (not given in the papers).
## Needed for the Christoffel symbols in RODES_qK.pyx (PJ12 Eqs. 11-12).
## (Checked symbolically against d/dr of F1 above.)
cdef double dF1dr(double r, double r_s) noexcept nogil:

    if r > _SERIES_SWITCH * r_s:
        return _rs_dF1dr_series(r_s / r) / r_s

    cdef double temp = (-5.0/(16.0*r_s*r*(r - r_s)))*(4.0*(-9.0*r*r + 9.0*r_s*r - r_s*r_s) - (2.0*r - r_s)*(2.0*r - r_s)*(r_s*r_s + 6.0*r_s*r - 6.0*r*r) / (r*(r - r_s)))

    return temp - (15.0/(4.0*r_s*r_s))*(log(r/(r - r_s))*(2.0*r - r_s) - r_s)


## dF2/dr: analytic derivative of GB06 Eq. (A.10) (not given in the papers).
## Needed for the Christoffel symbols in RODES_qK.pyx (PJ12 Eqs. 11-12).
## (Checked symbolically against d/dr of F2 above.)
cdef double dF2dr(double r, double r_s) noexcept nogil:

    if r > _SERIES_SWITCH * r_s:
        return _rs_dF2dr_series(r_s / r) / r_s
    return (-5.0/(8.0*r_s))*(6.0 + r_s*r_s / (r*r)) + (15.0/(8.0*r_s*r_s))*(4.0*r*log(r/(r - r_s)) - r_s*(2.0*r*r - r_s*r_s) / (r*(r - r_s)))


# Quadrupole moment correction F(r) functions in long double precision.
## Same expressions as above, including the large-r series. The
## long-double versions are used in BOUNDARY_CONDITIONS.pyx, where the image
## plane sits at d = 10^4 b_max and PJ12 (below Eq. 20) warn that d/r_i -> 1
## causes round-off error in the initial photon momentum. Note that long
## double alone does NOT rescue the closed forms at the image plane (their
## error there is still ~1e4 times the true value); the series does.

## F1(r): GB06 Eq. (A.9), long double
cdef long double F1_l(long double r, long double r_s) noexcept nogil:

    if r > _SERIES_SWITCH * r_s:
        return _F1_series_l(r_s / r)
    return -5.0*(2.0*r - r_s) * (r_s*r_s + 6.0*r_s*r - 6.0*r*r) / (16.0*r_s*r*(r - r_s)) - 15.0*r*(r - r_s)*logl(r/(r - r_s)) / (4.0*r_s*r_s)

## F2(r): GB06 Eq. (A.10), long double
cdef long double F2_l(long double r, long double r_s) noexcept nogil:

    if r > _SERIES_SWITCH * r_s:
        return _F2_series_l(r_s / r)
    return 5.0*(r_s*r_s - 3.0*r_s*r - 6.0*r*r) / (8.0*r_s*r) + 15.0*(2.0*r*r - r_s*r_s)*logl(r/(r - r_s)) / (8.0*r_s*r_s)

## dF1/dr: derivative of GB06 Eq. (A.9), long double
cdef long double dF1dr_l(long double r, long double r_s) noexcept nogil:

    if r > _SERIES_SWITCH * r_s:
        return _rs_dF1dr_series_l(r_s / r) / r_s

    cdef long double temp = (-5.0/(16.0*r_s*r*(r - r_s)))*(4.0*(-9.0*r*r + 9.0*r_s*r - r_s*r_s) - (2.0*r - r_s)*(2.0*r - r_s)*(r_s*r_s + 6.0*r_s*r - 6.0*r*r) / (r*(r - r_s)))
    return temp - (15.0/(4.0*r_s*r_s))*(logl(r/(r - r_s))*(2.0*r - r_s) - r_s)

## dF2/dr: derivative of GB06 Eq. (A.10), long double
cdef long double dF2dr_l(long double r, long double r_s) noexcept nogil:

    if r > _SERIES_SWITCH * r_s:
        return _rs_dF2dr_series_l(r_s / r) / r_s
    return (-5.0/(8.0*r_s))*(6.0 + r_s*r_s / (r*r)) + (15.0/(8.0*r_s*r_s))*(4.0*r*logl(r/(r - r_s)) - r_s*(2.0*r*r - r_s*r_s) / (r*(r - r_s)))


# Quasi-Kerr CONTRAVARIANT metric components: g^uv = g_K^uv + kappa h^uv
##
## The papers write the metric in covariant form, g_uv = g^K_uv + kappa h_uv
## (GB06 Eq. 2.13, PJ12 Eq. 1, B12 Eq. 2), but give the correction h^uv only in
## contravariant form (GB06 Eqs. 2.14-2.16, PJ12 Eq. 5, B12 Eq. 6). This module
## works entirely with contravariant components:
##   - g_K^uv is the inverse of the covariant Kerr metric in PJ12 Eq. (2) /
##     B12 Eq. (3). (Checked symbolically: the Kerr parts of g_11 ... g_44
##     below are exactly the inverse of that line element.)
##   - h^uv is added with a PLUS sign, i.e. the listed h^uv are taken to be
##     the contravariant correction itself, not h_uv with indices raised
##     (which would enter with a minus sign at first order). GB06 write their
##     Hartle-Thorne expressions the same way (GB06 Eq. 2.12:
##     g^tt_HT = g^tt_Ka2 + eps(...)).
##     This sign is the one that reproduces GB06 Eq. (2.3): expanding g_tt at
##     large r with F1 -> -(M/r)^3 gives Q = Q^K - kappa*M^3.

## g^tt: Kerr part = -[(r^2+a^2)^2 - a^2 Delta sin^2(theta)] / (Sigma Delta),
##       rewritten using (r^2+a^2)^2 - a^2 Delta sin^2 = (r^2+a^2) Sigma + r_s r a^2 sin^2;
##       correction h^tt = (1 - r_s/r)^-1 (1 - 3cos^2 theta) F1  (GB06 Eq. 2.14, PJ12 Eq. 5)
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
                 double func_r) noexcept nogil:

    cdef double temp = r*r + a*a + r_s * r * a * a * sin_theta * sin_theta / Sigma

    temp *= -1.0 / Delta

    return temp + kappa * func_theta * F1_r / func_r

## g^t phi: Kerr only, -r_s a r / (Sigma Delta); h^t phi = 0 (GB06 Eq. 2.16, PJ12 Eq. 5).
## The quadrupole correction does not change frame dragging (B12, below Eq. 6).
cdef double g_14(double r,
                 double theta,
                 double r_s,
                 double a,
                 double Sigma,
                 double Delta) noexcept nogil:

    return -r_s*a*r / (Sigma*Delta)

## g^rr: Kerr part = Delta / Sigma;
##       correction h^rr = (1 - r_s/r) (1 - 3cos^2 theta) F1  (GB06 Eq. 2.14, PJ12 Eq. 5)
cdef double g_22(double r,
                 double theta,
                 double Sigma,
                 double Delta,
                 double kappa,
                 double func_theta,
                 double F1_r,
                 double func_r) noexcept nogil:

    return Delta / Sigma + kappa*F1_r*func_theta*func_r

## g^theta theta: Kerr part = 1 / Sigma;
##       correction h^theta theta = -(1/r^2) (1 - 3cos^2 theta) F2  (GB06 Eq. 2.15, PJ12 Eq. 5)
cdef double g_33(double r,
                 double theta,
                 double Sigma,
                 double kappa,
                 double func_theta,
                 double F2_r) noexcept nogil:

    return 1.0 / Sigma - kappa * F2_r * func_theta / (r * r)

## g^phi phi: Kerr part = (Delta - a^2 sin^2 theta) / (Sigma Delta sin^2 theta);
##       correction h^phi phi = -(1/(r^2 sin^2 theta)) (1 - 3cos^2 theta) F2  (GB06 Eq. 2.15, PJ12 Eq. 5)
cdef double g_44(double r,
                 double theta,
                 double a,
                 double Sigma,
                 double Delta,
                 double kappa,
                 double sin_theta,
                 double func_theta,
                 double F2_r) noexcept nogil:

    cdef double temp = (Delta - a * a * sin_theta * sin_theta) / (Sigma * Delta)

    temp -= kappa * F2_r * func_theta / (r * r)

    return temp / (sin_theta * sin_theta)

## Determinant of the contravariant (t, phi) block: det = g^tt g^phiphi - (g^tphi)^2.
## Used to invert that 2x2 block and obtain the covariant components:
##   g_tt = g^phiphi / det,   g_tphi = -g^tphi / det,   g_phiphi = g^tt / det.
## (g_rr = 1/g^rr and g_thetatheta = 1/g^thetatheta, since those blocks are diagonal.)
## The covariant components are what enter PJ12 Eqs. (9)-(10) and (13)-(14);
## the denominator g_phiphi g_tt - g_tphi^2 in Eqs. (9)-(10) equals 1/det.
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
                  double sin_theta) noexcept nogil:

    cdef:
        double G_11 = g_11(r, theta, r_s, a, Sigma, Delta, kappa, sin_theta, func_theta, F1_r, func_r)
        double G_14 = g_14(r, theta, r_s, a, Sigma, Delta)
        double G_44 = g_44(r, theta, a, Sigma, Delta, kappa, sin_theta, func_theta, F2_r)

    return G_11*G_44 - G_14*G_14

# Quasi-Kerr CONTRAVARIANT metric components in long double precision:
## same expressions as above (g^uv = g_K^uv + kappa h^uv); used in
## BOUNDARY_CONDITIONS.pyx to set up the photon at the distant image plane.

## g^tt, long double (see g_11)
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
                        long double func_r) noexcept nogil:

    return -1.0*(r*r + a*a + r_s*r*a*a*sin_theta*sin_theta / Sigma) / Delta + kappa*func_theta*F1_r / func_r


## g^t phi, long double (see g_14)
cdef long double g_14_l(long double r,
                        long double theta,
                        long double r_s,
                        long double a,
                        long double Sigma,
                        long double Delta) noexcept nogil:

    return -r_s*a*r / (Sigma*Delta)

## g^rr, long double (see g_22)
cdef long double g_22_l(long double r,
                        long double theta,
                        long double Sigma,
                        long double Delta,
                        long double kappa,
                        long double func_theta,
                        long double F1_r,
                        long double func_r) noexcept nogil:

    return Delta / Sigma + kappa*F1_r*func_theta*func_r

## g^theta theta, long double (see g_33)
cdef long double g_33_l(long double r,
                        long double theta,
                        long double Sigma,
                        long double kappa,
                        long double func_theta,
                        long double F2_r) noexcept nogil:

    return 1.0 / Sigma - kappa*F2_r*func_theta / (r*r)

## g^phi phi, long double (see g_44)
cdef long double g_44_l(long double r,
                        long double theta,
                        long double a,
                        long double Sigma,
                        long double Delta,
                        long double kappa,
                        long double sin_theta,
                        long double func_theta,
                        long double F2_r) noexcept nogil:

    cdef long double temp = Delta - a * a * sin_theta * sin_theta

    temp /= Sigma * Delta * sin_theta * sin_theta

    return  temp - kappa * F2_r * func_theta / (r * r * sin_theta * sin_theta)


## Determinant of the contravariant (t, phi) block, long double (see det_g)
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
                         long double sin_theta) noexcept nogil:

    cdef:
        long double G_11 = g_11_l(r, theta, r_s, a, Sigma, Delta, kappa, sin_theta, func_theta, F1_r, func_r)
        long double G_14 = g_14_l(r, theta, r_s, a, Sigma, Delta)
        long double G_44 = g_44_l(r, theta, a, Sigma, Delta, kappa, sin_theta, func_theta, F2_r)

    return G_11 * G_44 - G_14 * G_14
