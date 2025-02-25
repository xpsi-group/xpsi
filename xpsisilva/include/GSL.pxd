cdef extern from "gsl/gsl_errno.h":

  cdef enum:
    GSL_SUCCESS  = 0 
    GSL_FAILURE  = -1
    GSL_CONTINUE = -2  # iteration has not converged */
    GSL_EDOM     = 1   # input domain error, e.g sqrt(-1) */
    GSL_ERANGE   = 2   # output range error, e.g. exp(1e100) */
    GSL_EFAULT   = 3   # invalid pointer */
    GSL_EINVAL   = 4   # invalid argument supplied by user */
    GSL_EFAILED  = 5   # generic failure */
    GSL_EFACTOR  = 6   # factorization failed */
    GSL_ESANITY  = 7   # sanity check failed - shouldn't happen */
    GSL_ENOMEM   = 8   # malloc failed */
    GSL_EBADFUNC = 9   # problem with user-supplied function */
    GSL_ERUNAWAY = 10  # iterative process is out of control */
    GSL_EMAXITER = 11  # exceeded max number of iterations */
    GSL_EZERODIV = 12  # tried to divide by zero */
    GSL_EBADTOL  = 13  # user specified an invalid tolerance */
    GSL_ETOL     = 14  # failed to reach the specified tolerance */
    GSL_EUNDRFLW = 15  # underflow */
    GSL_EOVRFLW  = 16  # overflow  */
    GSL_ELOSS    = 17  # loss of accuracy */
    GSL_EROUND   = 18  # failed because of roundoff error */
    GSL_EBADLEN  = 19  # matrix, vector lengths are not conformant */
    GSL_ENOTSQR  = 20  # matrix not square */
    GSL_ESING    = 21  # apparent singularity detected */
    GSL_EDIVERGE = 22  # integral or series is divergent */
    GSL_EUNSUP   = 23  # requested feature is not supported by the hardware */
    GSL_EUNIMPL  = 24  # requested feature not (yet) implemented */
    GSL_ECACHE   = 25  # cache limit exceeded */
    GSL_ETABLE   = 26  # table limit exceeded */
    GSL_ENOPROG  = 27  # iteration is not making progress towards solution */
    GSL_ENOPROGJ = 28  # jacobian evaluations are not improving the solution */
    GSL_ETOLF    = 29  # cannot reach the specified tolerance in F */
    GSL_ETOLX    = 30  # cannot reach the specified tolerance in X */
    GSL_ETOLG    = 31  # cannot reach the specified tolerance in gradient */
    GSL_EOF      = 32  # end of file */


  ctypedef void gsl_error_handler_t (const char * reason, const char * file, int line, int gsl_errno) nogil

  gsl_error_handler_t * gsl_set_error_handler_off () nogil

  const char * gsl_strerror (const int gsl_errno) nogil

  gsl_error_handler_t *gsl_set_error_handler(gsl_error_handler_t * new_handler) nogil


cdef extern from "gsl/gsl_math.h":

    double M_E

    double M_LOG2E

    double M_LOG10E

    double M_SQRT2

    double M_SQRT1_2

    double M_SQRT3

    double M_PI

    double M_PI_2

    double M_PI_4

    double M_SQRTPI

    double M_2_SQRTPI

    double M_1_PI

    double M_2_PI

    double M_LN10

    double M_LN2

    double M_LNPI

    double M_EULER

    int  gsl_isnan(double x) nogil

    int  gsl_isinf(double x) nogil

    int  gsl_finite(double x) nogil

    double  gsl_log1p(double x) nogil

    double  gsl_expm1(double x) nogil

    double  gsl_hypot(double x, double y) nogil

    double  gsl_acosh(double x) nogil

    double  gsl_asinh(double x) nogil

    double  gsl_atanh(double x) nogil

    double  gsl_ldexp(double x, int e) nogil

    double  gsl_frexp(double x, int * e) nogil

    double  gsl_pow_int(double x, int n) nogil

    double  gsl_pow_2(double x) nogil

    double  gsl_pow_3(double x) nogil

    double  gsl_pow_4(double x) nogil

    double  gsl_pow_5(double x) nogil

    double  gsl_pow_6(double x) nogil

    double  gsl_pow_7(double x) nogil

    double  gsl_pow_8(double x) nogil

    double  gsl_pow_9(double x) nogil

    int GSL_SIGN(double x) nogil

    int GSL_IS_ODD(int n) nogil

    int GSL_IS_EVEN(int n) nogil

    double GSL_MAX(double a, double  b) nogil

    double GSL_MIN(double a, double  b) nogil

    double  GSL_MAX_DBL(double a, double b) nogil

    double  GSL_MIN_DBL(double a, double b) nogil

    int  GSL_MAX_INT(int a, int b) nogil

    int  GSL_MIN_INT(int a, int b) nogil

    long double  GSL_MAX_LDBL(long double a, long double b) nogil

    long double  GSL_MIN_LDBL(long double a, long double b) nogil

    int  gsl_fcmp(double x, double y, double epsilon) nogil

    # Definition of an arbitrary function with parameters
    ctypedef struct gsl_function:
      double (* function) (double x, void * params) nogil
      void * params

    double GSL_FN_EVAL(gsl_function * F, double x) nogil

    # Definition of an arbitrary function returning two values, r1, r2
    ctypedef struct gsl_function_fdf:
      double (* f) (double x, void * params) nogil
      double (* df) (double x, void * params) nogil
      void (* fdf) (double x, void * params, double * f, double * df) nogil
      void * params

    double GSL_FN_FDF_EVAL_F(gsl_function_fdf * FDF, double x) nogil 

    GSL_FN_FDF_EVAL_DF(gsl_function_fdf * FDF,double x) nogil 

    GSL_FN_FDF_EVAL_F_DF(gsl_function_fdf * FDF,double x, double y,double dy) nogil 


cdef extern from "gsl/gsl_interp.h":
  
  ctypedef struct gsl_interp_accel
  
  ctypedef struct gsl_interp_type

  ctypedef struct gsl_interp:
    const gsl_interp_type * type
    double  xmin
    double  xmax
    size_t  size
    void * state

  
  gsl_interp_type * gsl_interp_linear
  gsl_interp_type * gsl_interp_polynomial
  gsl_interp_type * gsl_interp_cspline
  gsl_interp_type * gsl_interp_cspline_periodic
  gsl_interp_type * gsl_interp_akima
  gsl_interp_type * gsl_interp_akima_periodic
  gsl_interp_type * gsl_interp_steffen
  
  gsl_interp_accel * gsl_interp_accel_alloc() nogil
  
  size_t gsl_interp_accel_find(gsl_interp_accel * a,  double x_array[], size_t size, double x) nogil
  
  int gsl_interp_accel_reset (gsl_interp_accel * a) nogil
  
  void gsl_interp_accel_free(gsl_interp_accel * a) nogil
  
  gsl_interp * gsl_interp_alloc( gsl_interp_type * T, size_t n) nogil
       
  int gsl_interp_init(gsl_interp * obj,  double xa[],  double ya[], size_t size) nogil
  
  char * gsl_interp_name( gsl_interp * interp) nogil
  unsigned int gsl_interp_min_size( gsl_interp * interp) nogil
  
  int gsl_interp_eval_e( gsl_interp * obj,
                     double xa[],  double ya[], double x,
                    gsl_interp_accel * a, double * y) nogil
  
  double gsl_interp_eval( gsl_interp * obj,
                   double xa[],  double ya[], double x,
                  gsl_interp_accel * a) nogil
  
  int gsl_interp_eval_deriv_e( gsl_interp * obj,
                           double xa[],  double ya[], double x,
                          gsl_interp_accel * a,
                          double * d) nogil
  
  double gsl_interp_eval_deriv( gsl_interp * obj,
                         double xa[],  double ya[], double x,
                        gsl_interp_accel * a) nogil
  
  int gsl_interp_eval_deriv2_e( gsl_interp * obj,
                            double xa[],  double ya[], double x,
                           gsl_interp_accel * a,
                           double * d2) nogil
  
  double gsl_interp_eval_deriv2( gsl_interp * obj,
                          double xa[],  double ya[], double x,
                         gsl_interp_accel * a) nogil
  
  int gsl_interp_eval_integ_e( gsl_interp * obj,
                           double xa[],  double ya[],
                          double a, double b,
                          gsl_interp_accel * acc,
                          double * result) nogil
  
  double gsl_interp_eval_integ( gsl_interp * obj,
                         double xa[],  double ya[],
                        double a, double b,
                        gsl_interp_accel * acc) nogil
  
  void gsl_interp_free(gsl_interp * interp) nogil
  
  size_t gsl_interp_bsearch( double x_array[], double x,
                            size_t index_lo, size_t index_hi) nogil
  
cdef extern from "gsl/gsl_spline.h":
  ctypedef struct gsl_spline:
      gsl_interp * interp
      double  * x
      double  * y
      size_t  size

  
  gsl_spline * gsl_spline_alloc( gsl_interp_type * T, size_t size) nogil
       
  int gsl_spline_init(gsl_spline * spline,  double xa[],  double ya[], size_t size) nogil
  
  
  int gsl_spline_eval_e( gsl_spline * spline, double x,
                    gsl_interp_accel * a, double * y) nogil
  
  double gsl_spline_eval( gsl_spline * spline, double x, gsl_interp_accel * a) nogil
  
  int gsl_spline_eval_deriv_e( gsl_spline * spline, double x,
                          gsl_interp_accel * a, double * y) nogil
  
  double gsl_spline_eval_deriv( gsl_spline * spline, double x, gsl_interp_accel * a) nogil
  
  int gsl_spline_eval_deriv2_e( gsl_spline * spline, double x,
                           gsl_interp_accel * a, double * y) nogil
  
  double gsl_spline_eval_deriv2( gsl_spline * spline, double x,
                         gsl_interp_accel * a) nogil
  
  int gsl_spline_eval_integ_e( gsl_spline * spline, double a, double b,
                          gsl_interp_accel * acc, double * y) nogil
  
  double gsl_spline_eval_integ( gsl_spline * spline, double a, double b,
                        gsl_interp_accel * acc) nogil
  
  void gsl_spline_free(gsl_spline * spline) nogil
  

cdef extern from "gsl/gsl_interp2d.h":
  
    ctypedef struct gsl_interp2d_type

    ctypedef struct gsl_interp2d

    gsl_interp2d_type * gsl_interp2d_bilinear
    gsl_interp2d_type * gsl_interp2d_bicubic
  
  
cdef extern from "gsl/gsl_spline2d.h":
  
  ctypedef struct gsl_spline2d
  
  gsl_spline2d * gsl_spline2d_alloc (const gsl_interp2d_type * T, size_t xsize, size_t ysize) nogil
  
  int gsl_spline2d_init (gsl_spline2d * spline, const double xa[], const double ya[], const double za[], size_t xsize, size_t ysize) nogil
  
  void gsl_spline2d_free (gsl_spline2d * spline) nogil
  
  const char * gsl_spline2d_name (const gsl_spline2d * spline) nogil
  
  unsigned int gsl_spline2d_min_size (const gsl_spline2d * spline) nogil
  
  double gsl_spline2d_eval (const gsl_spline2d * spline, const double x, const double y, gsl_interp_accel * xacc, gsl_interp_accel * yacc) nogil

  int gsl_spline2d_set (const gsl_spline2d * interp, double zarr[], const size_t i, const size_t j, const double z) nogil
  
  double gsl_spline2d_get(const gsl_spline2d * interp, const double zarr[], const size_t i, const size_t j) nogil

cdef extern from "gsl/gsl_integration.h":

    ctypedef struct gsl_integration_workspace
    ctypedef struct gsl_integration_qaws_table
    ctypedef struct  gsl_integration_qawo_table
    cdef enum:
      GSL_INTEG_GAUSS15 = 1
      GSL_INTEG_GAUSS21 = 2
      GSL_INTEG_GAUSS31 = 3
      GSL_INTEG_GAUSS41 = 4
      GSL_INTEG_GAUSS51 = 5
      GSL_INTEG_GAUSS61 = 6
    cdef enum gsl_integration_qawo_enum:
      GSL_INTEG_COSINE, GSL_INTEG_SINE

    int  gsl_integration_qng(gsl_function *f, double a, double b, double epsabs, double epsrel, double * result, double * abserr, size_t * neval) nogil

    gsl_integration_workspace *  gsl_integration_workspace_alloc(size_t n) nogil

    void  gsl_integration_workspace_free(gsl_integration_workspace * w) nogil

    int  gsl_integration_qag(gsl_function *f, double a, double b, double epsabs, double epsrel, size_t limit, int key, gsl_integration_workspace * workspace, double * result, double * abserr) nogil

    int  gsl_integration_qags(gsl_function * f, double a, double b, double epsabs, double epsrel, size_t limit, gsl_integration_workspace * workspace, double *result, double *abserr) nogil

    int  gsl_integration_qagp(gsl_function * f, double *pts, size_t npts, double epsabs, double epsrel, size_t limit, gsl_integration_workspace * workspace, double *result, double *abserr) nogil

    int  gsl_integration_qagi(gsl_function * f, double epsabs, double epsrel, size_t limit, gsl_integration_workspace * workspace, double *result, double *abserr) nogil

    int  gsl_integration_qagiu(gsl_function * f, double a, double epsabs, double epsrel, size_t limit, gsl_integration_workspace * workspace, double *result, double *abserr) nogil

    int  gsl_integration_qagil(gsl_function * f, double b, double epsabs, double epsrel, size_t limit, gsl_integration_workspace * workspace, double *result, double *abserr) nogil

    int  gsl_integration_qawc(gsl_function *f, double a, double b, double c, double epsabs, double epsrel, size_t limit, gsl_integration_workspace * workspace, double * result, double * abserr) nogil

    gsl_integration_qaws_table *  gsl_integration_qaws_table_alloc(double alpha, double beta, int mu, int nu) nogil

    int  gsl_integration_qaws_table_set(gsl_integration_qaws_table * t, double alpha, double beta, int mu, int nu) nogil

    void  gsl_integration_qaws_table_free(gsl_integration_qaws_table * t) nogil

    int  gsl_integration_qaws(gsl_function * f, double a, double b, gsl_integration_qaws_table * t, double epsabs, double epsrel, size_t limit, gsl_integration_workspace * workspace, double *result, double *abserr) nogil

    gsl_integration_qawo_table *  gsl_integration_qawo_table_alloc(double omega, double L,  gsl_integration_qawo_enum sine, size_t n) nogil

    int  gsl_integration_qawo_table_set(gsl_integration_qawo_table * t, double omega, double L,  gsl_integration_qawo_enum sine) nogil

    int  gsl_integration_qawo_table_set_length(gsl_integration_qawo_table * t, double L) nogil

    void  gsl_integration_qawo_table_free(gsl_integration_qawo_table * t) nogil

    int  gsl_integration_qawo(gsl_function * f, double a, double epsabs, double epsrel, size_t limit, gsl_integration_workspace * workspace, gsl_integration_qawo_table * wf, double *result, double *abserr) nogil

    int  gsl_integration_qawf(gsl_function * f, double a, double epsabs, size_t limit, gsl_integration_workspace * workspace, gsl_integration_workspace * cycle_workspace, gsl_integration_qawo_table * wf, double *result, double *abserr) nogil

    ctypedef struct gsl_integration_cquad_ival

    ctypedef struct gsl_integration_cquad_workspace

    gsl_integration_cquad_workspace * gsl_integration_cquad_workspace_alloc (const size_t n) nogil

    void gsl_integration_cquad_workspace_free (gsl_integration_cquad_workspace * w) nogil

    int gsl_integration_cquad (const gsl_function * f, double a, double b, double epsabs, double epsrel, gsl_integration_cquad_workspace * ws, double *result, double *abserr, size_t * nevals) nogil

cdef extern from "gsl/gsl_sf_legendre.h":

  double  gsl_sf_legendre_P1(double x) nogil

  double  gsl_sf_legendre_P2(double x) nogil

  double  gsl_sf_legendre_P3(double x) nogil

  int  gsl_sf_legendre_P1_e(double x, gsl_sf_result * result) nogil

  int  gsl_sf_legendre_P2_e(double x, gsl_sf_result * result) nogil

  int  gsl_sf_legendre_P3_e(double x, gsl_sf_result * result) nogil

  double  gsl_sf_legendre_Pl(int l, double x) nogil

  int  gsl_sf_legendre_Pl_e(int l, double x, gsl_sf_result * result) nogil

  int  gsl_sf_legendre_Pl_array(int lmax, double x, double result_array[]) nogil

  double  gsl_sf_legendre_Q0(double x) nogil

  int  gsl_sf_legendre_Q0_e(double x, gsl_sf_result * result) nogil

  double  gsl_sf_legendre_Q1(double x) nogil

  int  gsl_sf_legendre_Q1_e(double x, gsl_sf_result * result) nogil

  double  gsl_sf_legendre_Ql(int l, double x) nogil

  int  gsl_sf_legendre_Ql_e(int l, double x, gsl_sf_result * result) nogil

  double  gsl_sf_legendre_Plm(int l, int m, double x) nogil

  int  gsl_sf_legendre_Plm_e(int l, int m, double x, gsl_sf_result * result) nogil

  int  gsl_sf_legendre_Plm_array(int lmax, int m, double x, double result_array[]) nogil

  double  gsl_sf_legendre_sphPlm(int l, int m, double x) nogil

  int  gsl_sf_legendre_sphPlm_e(int l, int m, double x, gsl_sf_result * result) nogil

  int  gsl_sf_legendre_sphPlm_array(int lmax, int m, double x, double result_array[]) nogil

  int  gsl_sf_legendre_array_size(int lmax, int m) nogil

  ctypedef enum gsl_sf_legendre_t:
    GSL_SF_LEGENDRE_SCHMIDT,
    GSL_SF_LEGENDRE_SPHARM,
    GSL_SF_LEGENDRE_FULL,
    GSL_SF_LEGENDRE_NONE

  int gsl_sf_legendre_array(const gsl_sf_legendre_t norm,
                          const size_t lmax, const double x,
                          double result_array[]) nogil

  int gsl_sf_legendre_array_e(const gsl_sf_legendre_t norm,
                            const size_t lmax, const double x,
                            const double csphase,
                            double result_array[]) nogil

  size_t gsl_sf_legendre_array_n(const size_t lmax) nogil

  size_t gsl_sf_legendre_array_index(const size_t l, const size_t m) nogil

  double  gsl_sf_conicalP_half(double lambd, double x) nogil

  int  gsl_sf_conicalP_half_e(double lambd, double x, gsl_sf_result * result) nogil

  double  gsl_sf_conicalP_mhalf(double lambd, double x) nogil

  int  gsl_sf_conicalP_mhalf_e(double lambd, double x, gsl_sf_result * result) nogil

  double  gsl_sf_conicalP_0(double lambd, double x) nogil

  int  gsl_sf_conicalP_0_e(double lambd, double x, gsl_sf_result * result) nogil

  double  gsl_sf_conicalP_1(double lambd, double x) nogil

  int  gsl_sf_conicalP_1_e(double lambd, double x, gsl_sf_result * result) nogil

  double  gsl_sf_conicalP_sph_reg(int l, double lambd, double x) nogil

  int  gsl_sf_conicalP_sph_reg_e(int l, double lambd, double x, gsl_sf_result * result) nogil

  double  gsl_sf_conicalP_cyl_reg(int m, double lambd, double x) nogil

  int  gsl_sf_conicalP_cyl_reg_e(int m, double lambd, double x, gsl_sf_result * result) nogil

  double  gsl_sf_legendre_H3d_0(double lambd, double eta) nogil

  int  gsl_sf_legendre_H3d_0_e(double lambd, double eta, gsl_sf_result * result) nogil

  double  gsl_sf_legendre_H3d_1(double lambd, double eta) nogil

  int  gsl_sf_legendre_H3d_1_e(double lambd, double eta, gsl_sf_result * result) nogil

  double  gsl_sf_legendre_H3d(int l, double lambd, double eta) nogil

  int  gsl_sf_legendre_H3d_e(int l, double lambd, double eta, gsl_sf_result * result) nogil

  int  gsl_sf_legendre_H3d_array(int lmax, double lambd, double eta, double result_array[]) nogil

cdef extern from "gsl/gsl_sf_result.h":
    ctypedef struct gsl_sf_result:
        double val
        double err

    ctypedef struct gsl_sf_result_e10:
        double val
        double err
        int    e10
    
cdef extern from "stdio.h":
    ctypedef struct FILE
    FILE *fopen(char *path, char *mode) nogil
    int fclose(FILE *strea) nogil
    cdef FILE *stdout
    int scanf(char *format, ...) nogil

cdef extern from "gsl/gsl_block_double.h":

    ctypedef struct gsl_block:
        size_t size
        double * data

    gsl_block *  gsl_block_alloc(size_t n) nogil

    gsl_block *  gsl_block_calloc(size_t n) nogil

    void  gsl_block_free(gsl_block * b) nogil

    int  gsl_block_fread(FILE * stream, gsl_block * b) nogil

    int  gsl_block_fwrite(FILE * stream, gsl_block * b) nogil

    int  gsl_block_fscanf(FILE * stream, gsl_block * b) nogil

    int  gsl_block_fprintf(FILE * stream, gsl_block * b, char * format) nogil

    size_t gsl_block_size (gsl_block * b) nogil
    double * gsl_block_data (gsl_block * b) nogil
    
cdef extern from "gsl/gsl_vector.h":

    ctypedef struct gsl_vector:
        size_t size
        size_t stride
        double *data
        gsl_block *block
        int owner

    ctypedef struct gsl_vector_view:
        gsl_vector vector

    ctypedef struct gsl_vector_const_view:
        gsl_vector vector


    # Allocation
    gsl_vector *  gsl_vector_alloc(size_t n) nogil

    gsl_vector *  gsl_vector_calloc(size_t n) nogil

    gsl_vector_alloc_from_block(gsl_block * b, size_t offset,
                              size_t n, size_t stride) nogil

    gsl_vector *gsl_vector_alloc_from_vector(gsl_vector * v,
                         size_t offset, size_t n, size_t stride) nogil

    void  gsl_vector_free(gsl_vector * v) nogil

    # Views
    gsl_vector_view  gsl_vector_view_array(double *base, size_t n) nogil

    gsl_vector_view  gsl_vector_subvector(gsl_vector *v, size_t offset, size_t n) nogil

    gsl_vector_view  gsl_vector_view_array_with_stride(double * base, size_t stride, size_t n) nogil

    gsl_vector_const_view  gsl_vector_const_view_array(double *base, size_t n) nogil

    gsl_vector_const_view  gsl_vector_const_view_array_with_stride(double * base, size_t stride, size_t n) nogil

    gsl_vector_const_view  gsl_vector_const_subvector(gsl_vector * v, size_t offset, size_t n) nogil

    gsl_vector_view  gsl_vector_subvector_with_stride(gsl_vector *v, size_t offset, size_t stride, size_t n) nogil

    gsl_vector_const_view  gsl_vector_const_subvector_with_stride(gsl_vector * v, size_t offset, size_t stride, size_t n) nogil


    # Operations
    double  gsl_vector_get(gsl_vector * v, size_t i) nogil

    void  gsl_vector_set(gsl_vector * v, size_t i, double x) nogil

    double *  gsl_vector_ptr(gsl_vector * v, size_t i) nogil

    double *  gsl_vector_const_ptr(gsl_vector * v, size_t i) nogil

    void  gsl_vector_set_zero(gsl_vector * v) nogil

    void  gsl_vector_set_all(gsl_vector * v, double x) nogil

    int  gsl_vector_set_basis(gsl_vector * v, size_t i) nogil

    # Reading and writing vectors
    int  gsl_vector_fread(FILE * stream, gsl_vector * v) nogil

    int  gsl_vector_fwrite(FILE * stream, gsl_vector * v) nogil

    int  gsl_vector_fscanf(FILE * stream, gsl_vector * v) nogil

    int  gsl_vector_fprintf(FILE * stream, gsl_vector * v, char * format) nogil

    # Copying or exchanging elements
    int  gsl_vector_memcpy(gsl_vector * dest, gsl_vector * src) nogil

    int  gsl_vector_reverse(gsl_vector * v) nogil

    int  gsl_vector_swap(gsl_vector * v, gsl_vector * w) nogil

    int  gsl_vector_swap_elements(gsl_vector * v, size_t i, size_t j) nogil

    # Finding maximum and minimum elements of vectors

    double  gsl_vector_max(gsl_vector * v) nogil

    double  gsl_vector_min(gsl_vector * v) nogil

    void  gsl_vector_minmax(gsl_vector * v, double * min_out, double * max_out) nogil

    size_t  gsl_vector_max_index(gsl_vector * v) nogil

    size_t  gsl_vector_min_index(gsl_vector * v) nogil

    void  gsl_vector_minmax_index(gsl_vector * v, size_t * imin, size_t * imax) nogil

    # Vector operations
    int  gsl_vector_add(gsl_vector * a, gsl_vector * b) nogil

    int  gsl_vector_sub(gsl_vector * a, gsl_vector * b) nogil

    int  gsl_vector_mul(gsl_vector * a, gsl_vector * b) nogil

    int  gsl_vector_div(gsl_vector * a, gsl_vector * b) nogil

    int  gsl_vector_scale(gsl_vector * a, double x) nogil

    int  gsl_vector_add_constant(gsl_vector * a, double x) nogil

    int  gsl_vector_isnull(gsl_vector * v) nogil
    
    
cdef extern from "gsl/gsl_matrix_double.h":

    ctypedef struct gsl_matrix:
        size_t size1
        size_t size2
        size_t tda
        double * data
        gsl_block * block
        int owner

    ctypedef struct gsl_matrix_view:
        gsl_matrix matrix

    ctypedef struct gsl_matrix_const_view

    # Allocation
    gsl_matrix *  gsl_matrix_alloc(size_t n1, size_t n2) nogil

    gsl_matrix *  gsl_matrix_calloc(size_t n1, size_t n2) nogil

    gsl_matrix *  gsl_matrix_alloc_from_block(gsl_block * b, size_t offset, size_t n1, size_t n2, size_t d2) nogil

    gsl_matrix * gsl_matrix_alloc_from_matrix (gsl_matrix * m,  size_t k1,  size_t k2,  size_t n1,  size_t n2) nogil

    gsl_vector * gsl_vector_alloc_row_from_matrix (gsl_matrix * m,  size_t i) nogil

    gsl_vector * gsl_vector_alloc_col_from_matrix (gsl_matrix * m,  size_t j) nogil

    void  gsl_matrix_free(gsl_matrix * m) nogil

    # Views
    gsl_matrix_view  gsl_matrix_submatrix(gsl_matrix * m, size_t k1, size_t k2, size_t n1, size_t n2) nogil

    gsl_vector_view  gsl_matrix_row(gsl_matrix * m, size_t i) nogil

    gsl_vector_view  gsl_matrix_column(gsl_matrix * m, size_t j) nogil

    gsl_vector_view  gsl_matrix_diagonal(gsl_matrix * m) nogil

    gsl_vector_view  gsl_matrix_subdiagonal(gsl_matrix * m, size_t k) nogil

    gsl_vector_view  gsl_matrix_superdiagonal(gsl_matrix * m, size_t k) nogil

    gsl_matrix_view  gsl_matrix_view_array(double * base, size_t n1, size_t n2) nogil

    gsl_matrix_view  gsl_matrix_view_array_with_tda(double * base, size_t n1, size_t n2, size_t tda) nogil

    gsl_matrix_view  gsl_matrix_view_vector(gsl_vector * v, size_t n1, size_t n2) nogil

    gsl_matrix_view  gsl_matrix_view_vector_with_tda(gsl_vector * v, size_t n1, size_t n2, size_t tda) nogil

    gsl_matrix_const_view  gsl_matrix_const_submatrix(gsl_matrix * m, size_t k1, size_t k2, size_t n1, size_t n2) nogil

    gsl_vector_const_view  gsl_matrix_const_row(gsl_matrix * m, size_t i) nogil

    gsl_vector_const_view  gsl_matrix_const_column(gsl_matrix * m, size_t j) nogil

    gsl_vector_const_view  gsl_matrix_const_diagonal(gsl_matrix * m) nogil

    gsl_vector_const_view  gsl_matrix_const_subdiagonal(gsl_matrix * m, size_t k) nogil

    gsl_vector_const_view  gsl_matrix_const_superdiagonal(gsl_matrix * m, size_t k) nogil

    gsl_matrix_const_view  gsl_matrix_const_view_array(double * base, size_t n1, size_t n2) nogil

    gsl_matrix_const_view  gsl_matrix_const_view_array_with_tda(double * base, size_t n1, size_t n2, size_t tda) nogil

    gsl_matrix_const_view  gsl_matrix_const_view_vector(gsl_vector * v, size_t n1, size_t n2) nogil

    gsl_matrix_const_view  gsl_matrix_const_view_vector_with_tda(gsl_vector * v, size_t n1, size_t n2, size_t tda) nogil


    # Operations
    double  gsl_matrix_get(gsl_matrix * m, size_t i, size_t j) nogil

    void  gsl_matrix_set(gsl_matrix * m, size_t i, size_t j, double x) nogil

    double *  gsl_matrix_ptr(gsl_matrix * m, size_t i, size_t j) nogil

    double *  gsl_matrix_const_ptr(gsl_matrix * m, size_t i, size_t j) nogil

    void  gsl_matrix_set_zero(gsl_matrix * m) nogil

    void  gsl_matrix_set_identity(gsl_matrix * m) nogil

    void  gsl_matrix_set_all(gsl_matrix * m, double x) nogil

    # Reading and writing matrices
    int  gsl_matrix_fread(FILE * stream, gsl_matrix * m) nogil

    int  gsl_matrix_fwrite(FILE * stream, gsl_matrix * m) nogil

    int  gsl_matrix_fscanf(FILE * stream, gsl_matrix * m) nogil

    int  gsl_matrix_fprintf(FILE * stream, gsl_matrix * m, char * format) nogil

    # Copying or exchanging elements
    int  gsl_matrix_memcpy(gsl_matrix * dest, gsl_matrix * src) nogil

    int  gsl_matrix_swap(gsl_matrix * m1, gsl_matrix * m2) nogil

    int  gsl_matrix_swap_rows(gsl_matrix * m, size_t i, size_t j) nogil

    int  gsl_matrix_swap_columns(gsl_matrix * m, size_t i, size_t j) nogil

    int  gsl_matrix_swap_rowcol(gsl_matrix * m, size_t i, size_t j) nogil

    int  gsl_matrix_transpose(gsl_matrix * m) nogil

    int  gsl_matrix_transpose_memcpy(gsl_matrix * dest, gsl_matrix * src) nogil

    # Finding maximum and minimum elements of matrices
    double  gsl_matrix_max(gsl_matrix * m) nogil

    double  gsl_matrix_min(gsl_matrix * m) nogil

    void  gsl_matrix_minmax(gsl_matrix * m, double * min_out, double * max_out) nogil

    void  gsl_matrix_max_index(gsl_matrix * m, size_t * imax, size_t * jmax) nogil

    void  gsl_matrix_min_index(gsl_matrix * m, size_t * imax, size_t * jmax) nogil

    void  gsl_matrix_minmax_index(gsl_matrix * m, size_t * imin, size_t * jmin, size_t * imax, size_t * jmax) nogil

    int  gsl_matrix_isnull(gsl_matrix * m) nogil

    # Matrix operations
    int  gsl_matrix_add(gsl_matrix * a, gsl_matrix * b) nogil

    int  gsl_matrix_sub(gsl_matrix * a, gsl_matrix * b) nogil

    int  gsl_matrix_mul_elements(gsl_matrix * a, gsl_matrix * b) nogil

    int  gsl_matrix_div_elements(gsl_matrix * a, gsl_matrix * b) nogil

    int  gsl_matrix_scale(gsl_matrix * a, double x) nogil

    int  gsl_matrix_add_constant(gsl_matrix * a, double x) nogil

    int gsl_matrix_add_diagonal (gsl_matrix * a,  double x) nogil

    # The functions below are obsolete
    int  gsl_matrix_get_row(gsl_vector * v, gsl_matrix * m, size_t i) nogil

    int  gsl_matrix_get_col(gsl_vector * v, gsl_matrix * m, size_t j) nogil

    int  gsl_matrix_set_row(gsl_matrix * m, size_t i, gsl_vector * v) nogil

    int  gsl_matrix_set_col(gsl_matrix * m, size_t j, gsl_vector * v) nogil
  

cdef extern from "gsl/gsl_odeiv2.h":
    
    #define GSL_ODEIV_FN_EVAL(S,t,y,f)  (*((S)->function))(t,y,f,(S)->params)
    
    ctypedef struct gsl_odeiv2_system:
        int (* function) (double t,  double y[], double dydt[], void * params) nogil
        int (* jacobian) (double t,  double y[], double * dfdy, double dfdt[], void * params) nogil
        size_t dimension
        void * params

    ctypedef struct gsl_odeiv2_step
    ctypedef struct gsl_odeiv2_control
    ctypedef struct gsl_odeiv2_evolve
    ctypedef struct gsl_odeiv2_driver

    ctypedef struct gsl_odeiv2_step_type

    gsl_odeiv2_step_type *gsl_odeiv2_step_rk2
    gsl_odeiv2_step_type *gsl_odeiv2_step_rk4
    gsl_odeiv2_step_type *gsl_odeiv2_step_rkf45
    gsl_odeiv2_step_type *gsl_odeiv2_step_rkck
    gsl_odeiv2_step_type *gsl_odeiv2_step_rk8pd
    gsl_odeiv2_step_type *gsl_odeiv2_step_rk2imp
    gsl_odeiv2_step_type *gsl_odeiv2_step_rk4imp
    gsl_odeiv2_step_type *gsl_odeiv2_step_bsimp
    gsl_odeiv2_step_type *gsl_odeiv2_step_rk1imp
    gsl_odeiv2_step_type *gsl_odeiv2_step_msadams
    gsl_odeiv2_step_type *gsl_odeiv2_step_msbdf

    gsl_odeiv2_step * gsl_odeiv2_step_alloc(gsl_odeiv2_step_type * T, size_t dim) nogil
    int  gsl_odeiv2_step_reset(gsl_odeiv2_step * s) nogil
    void gsl_odeiv2_step_free(gsl_odeiv2_step * s) nogil

    char * gsl_odeiv2_step_name(gsl_odeiv2_step * s) nogil
    unsigned int gsl_odeiv2_step_order(gsl_odeiv2_step * s) nogil

    int  gsl_odeiv2_step_apply(gsl_odeiv2_step * s, double t, double h, double y[], double yerr[],  double dydt_in[], double dydt_out[],  gsl_odeiv2_system * dydt) nogil

    int  gsl_odeiv2_step_set_driver(gsl_odeiv2_step * s, gsl_odeiv2_driver * d) nogil

    #ctypedef struct gsl_odeiv2_control

    ctypedef struct gsl_odeiv2_control_type

    cdef enum:
        GSL_ODEIV_HADJ_DEC = -1
        GSL_ODEIV_HADJ_NIL = 0  
        GSL_ODEIV_HADJ_INC = 1

    gsl_odeiv2_control * gsl_odeiv2_control_alloc(gsl_odeiv2_control_type * T) nogil
    int gsl_odeiv2_control_init(gsl_odeiv2_control * c, double eps_abs, double eps_rel, double a_y, double a_dydt) nogil
    void gsl_odeiv2_control_free(gsl_odeiv2_control * c) nogil
    int gsl_odeiv2_control_hadjust (gsl_odeiv2_control * c, gsl_odeiv2_step * s,  double y0[],  double yerr[],  double dydt[], double * h) nogil
    char * gsl_odeiv2_control_name(gsl_odeiv2_control * c) nogil
    int gsl_odeiv2_control_errlevel(gsl_odeiv2_control * c, double y, double dydt, double h, size_t ind, double *errlev) nogil
    int gsl_odeiv2_control_set_driver(gsl_odeiv2_control * c, gsl_odeiv2_driver * d) nogil

    gsl_odeiv2_control * gsl_odeiv2_control_standard_new(double eps_abs, double eps_rel, double a_y, double a_dydt) nogil
    gsl_odeiv2_control * gsl_odeiv2_control_y_new(double eps_abs, double eps_rel) nogil
    gsl_odeiv2_control * gsl_odeiv2_control_yp_new(double eps_abs, double eps_rel) nogil

    gsl_odeiv2_control * gsl_odeiv2_control_scaled_new(double eps_abs, double eps_rel, double a_y, double a_dydt,  double scale_abs[], size_t dim) nogil

    #ctypedef struct gsl_odeiv2_evolve

    gsl_odeiv2_evolve * gsl_odeiv2_evolve_alloc(size_t dim) nogil
    int gsl_odeiv2_evolve_apply(gsl_odeiv2_evolve * e, gsl_odeiv2_control * con, gsl_odeiv2_step * step, gsl_odeiv2_system * dydt, double * t, double t1, double * h, double y[]) nogil
    int gsl_odeiv2_evolve_apply_fixed_step(gsl_odeiv2_evolve * e, gsl_odeiv2_control * con, gsl_odeiv2_step * step, gsl_odeiv2_system * dydt, double * t, double h0, double y[]) nogil
    int gsl_odeiv2_evolve_reset(gsl_odeiv2_evolve * e) nogil
    void gsl_odeiv2_evolve_free(gsl_odeiv2_evolve * e) nogil
    int gsl_odeiv2_evolve_set_driver(gsl_odeiv2_evolve * e, gsl_odeiv2_driver * d) nogil

    #ctypedef struct gsl_odeiv2_driver

    gsl_odeiv2_driver *gsl_odeiv2_driver_alloc_y_new(gsl_odeiv2_system *sys, gsl_odeiv2_step_type * T, double hstart, double epsabs, double epsrel) nogil
    gsl_odeiv2_driver *gsl_odeiv2_driver_alloc_yp_new(gsl_odeiv2_system *sys, gsl_odeiv2_step_type * T, double hstart, double epsabs, double epsrel) nogil
    gsl_odeiv2_driver *gsl_odeiv2_driver_alloc_scaled_new(gsl_odeiv2_system *sys, gsl_odeiv2_step_type * T, double hstart, double epsabs, double epsrel, double a_y, double a_dydt, double scale_abs[]) nogil
    gsl_odeiv2_driver *gsl_odeiv2_driver_alloc_standard_new(gsl_odeiv2_system *sys, gsl_odeiv2_step_type * T, double hstart, double epsabs, double epsrel, double a_y, double a_dydt) nogil
    int gsl_odeiv2_driver_set_hmin(gsl_odeiv2_driver *d, double hmin) nogil
    int gsl_odeiv2_driver_set_hmax(gsl_odeiv2_driver *d, double hmax) nogil
    int gsl_odeiv2_driver_set_nmax(gsl_odeiv2_driver *d, unsigned long int nmax) nogil
    int gsl_odeiv2_driver_apply(gsl_odeiv2_driver *d, double *t, double t1, double y[]) nogil
    int gsl_odeiv2_driver_apply_fixed_step(gsl_odeiv2_driver *d, double *t, double h, unsigned long int n, double y[]) nogil
    int gsl_odeiv2_driver_reset(gsl_odeiv2_driver *d) nogil
    int gsl_odeiv2_driver_free(gsl_odeiv2_driver *d) nogil


cdef extern from "gsl/gsl_rng.h":

    ctypedef struct gsl_rng

    gsl_rng * gsl_rng_alloc(const gsl_rng_type *T) nogil
    void gsl_rng_set(const gsl_rng *r, unsigned long int s) nogil
    void gsl_rng_free(gsl_rng *r) nogil

    ctypedef struct gsl_rng_type

    gsl_rng_type *gsl_rng_default

    const gsl_rng_type * gsl_rng_env_setup() nogil

cdef extern from "gsl/gsl_randist.h":

    unsigned int gsl_ran_poisson(const gsl_rng *r, double mu) nogil
    void gsl_ran_poisson_array(const gsl_rng *r, size_t n, unsigned int array[],
                               double mu) nogil
    double gsl_ran_poisson_pdf(const unsigned int k, const double mu) nogil
