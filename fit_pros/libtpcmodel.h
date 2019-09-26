/// @file libtpcmodel.h
/// @brief Header file for libtpcmodel.
/// @author Vesa Oikonen
///
#ifndef _LIBTPCMODEL_H
#define _LIBTPCMODEL_H
/*****************************************************************************/

/*****************************************************************************/
#include "tpcclibConfig.h"
/*****************************************************************************/
//#include <omp.h>
/*****************************************************************************/
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <ctype.h>
#include <math.h>
#include <float.h>
#include <time.h>
/*****************************************************************************/
#include "libtpcmisc.h"
/*****************************************************************************/

/*****************************************************************************/
#ifndef MAX_PARAMETERS
/** Max nr of parameters */
#define MAX_PARAMETERS 50
#endif
#ifndef MAX_PARAMS
/** Max nr of parameters */
#define MAX_PARAMS MAX_PARAMETERS
#endif
/*****************************************************************************/

/*****************************************************************************/
/** Mersenne Twister state vector length */
#define TPCCLIB_MERTWI_NN 312
/** Mersenne Twister required constant */
#define TPCCLIB_MERTWI_A UINT64_C(0xB5026F5AA96619E9)
/** Struct needed in functions for Mersenne Twister MT19937.
 *  Contents are initiated in seed functions. */
typedef struct MERTWI {
  /** Constant N, set to by seed function. */
  unsigned int n;
  /** Constant M, set by seed function. */
  unsigned int m;
  /** Constant MATRIX_A, set by seed function. */
  uint64_t a;
  /** Constant UM, most significant 33 bits, set by seed function. */
  uint64_t um;
  /** Constant LM, least significant 31 bits, set by seed function. */
  uint64_t lm;
  /** The array for the state vector */
  uint64_t mt[TPCCLIB_MERTWI_NN]; 
  /** Index of mt array; mti==NN+1 means mt[NN] is not initialized */
  uint64_t mti;
} MERTWI;
/* mertwi */
uint32_t mertwiSeed32(void);
uint64_t mertwiSeed64(void);
void mertwiInit(MERTWI *mt);
void mertwiInitWithSeed64(MERTWI *mt, uint64_t seed);
void mertwiInitByArray64(MERTWI *mt, uint64_t init_key[], uint64_t key_length);
uint64_t mertwiRandomInt64(MERTWI *mt);
int64_t mertwiRandomInt63(MERTWI *mt);
double mertwiRandomDouble1(MERTWI *mt);
double mertwiRandomDouble2(MERTWI *mt);
double mertwiRandomDouble3(MERTWI *mt);
/*****************************************************************************/

/*****************************************************************************/
/** BOBYQA return codes */
typedef enum {
   BOBYQA_INVALID_ARGS = -1,
   BOBYQA_OUT_OF_MEMORY = -2,
   BOBYQA_ROUNDOFF_LIMITED = -3,
   BOBYQA_FAIL = -4, /* generic fail code */
   BOBYQA_SUCCESS = 0, /* generic success code */
   BOBYQA_MINF_MAX_REACHED = 1,
   BOBYQA_FTOL_REACHED = 2,
   BOBYQA_XTOL_REACHED = 3,
   BOBYQA_MAXEVAL_REACHED = 4,
   BOBYQA_RELFTOL_REACHED = 5,
   BOBYQA_ABSFTOL_REACHED = 6
} bobyqa_result;

/** BOBYQA objective function */
typedef double (*bobyqa_func)(int n, const double *x, void *func_data);

/** bobyqa data in a struct by VO */
typedef struct {
  /** N is the number of fitted variables and must be at least two */
  int n;
  /** NPT is the number of interpolation conditions. Its value must be in
      the interval [N+2,(N+1)(N+2)/2]. Choices that exceed 2*N+1 are not
      recommended. */
  int npt;
  /** Initial values of the variables to be fitted must be set in
   *  X(1),X(2),...,X(N).
   *  They will be changed to the values that give the least calculated F. */
  double *x;
  /** Size of array X */
  int x_size;
  /** Scale factors for fitted function parameters */
  double *xscale;
  /** Size of array XSCALE */
  int xscale_size;
  /** Nr of all parameters, including the fixed parameters */
  int nfull;
  /** Same as above, but this list contains all parameters required by the
   *  minimized function, i.e. both fitted and fixed parameters;
   *  Array size is nfull. */
  double *xfull;
  /** Index list of fitted parameters, that is, at which positions in xfull[]
   *  the parameters in x[] must be placed. */
  int *xplace; 
  /** For I=1,2,...,N, XL(I) and XU(I) must provide the lower and upper
      bounds, respectively, on X(I). The construction of quadratic models
      requires XL(I) to be strictly less than XU(I) for each I. Further,
      the contribution to a model from changes to the I-th variable is
      damaged severely by rounding errors if XU(I)-XL(I) is too small. */
  double *xl;
  /** Size of array XL */
  int xl_size;
  /** see xl. */
  double *xu;
  /** Size of array XU */
  int xu_size;
  /** RHOBEG and RHOEND must be set to the initial and final values of a trust
      region radius, so both must be positive with RHOEND no greater than
      RHOBEG. Typically, RHOBEG should be about one tenth of the greatest
      expected change to a variable, while RHOEND should indicate the
      accuracy that is required in the final values of the variables. An
      error return occurs if any of the differences XU(I)-XL(I), I=1,...,N,
      is less than 2*RHOBEG. */
  double rhobeg;
  /** see rhobeg */
  double rhoend; 
  /** Stopping rule: maximum allowed function value */
  double minf_max;
  /** Stopping rule: relative tolerance to function value */
  double ftol_rel;
  /** Stopping rule: absolute tolerance to function value */
  double ftol_abs;
  /** Stopping rule: max nr of function evaluations */ 
  int maxeval;
  /** Number of function evaluations */
  int nevals;

  /** SUBROUTINE CALFUN (N,X,F) has to be provided by the user. It must set
      F to the value of the objective function for the current values of the
      variables X(1),X(2),...,X(N), which are generated automatically in a
      way that satisfies the bounds given in XL and XU. */
  bobyqa_func objf;
  /** Pointer to data passed to objective function */
  void *objf_data;
  /** Minimum value of function */
  double minf;

  /** Pointer to total working memory; just for testing etc */
  double *wmptr;
  /** Pointer to working memory allocated by bobyqa_set_memory();
   *  bobyqa_free_memory() will free this memory and set it to NULL.
   *  If memory was not allocated by bobyqa_set_memory(), this is NULL. */
  double *lwmptr;
  /** Pointer to int working memory allocated by bobyqa_set_memory();
   *  bobyqa_free_memory() will free this memory and set it to NULL. */
  int *liwmptr;

  /** XBASE holds a shift of origin that should reduce the contributions
      from rounding errors to values of the model and Lagrange functions.
      XOPT is set to the displacement from XBASE of the trust region centre. */
  double *xbase; 
  /** Size of array XBASE */
  int xbase_size;
  /** XPT is a two-dimensional array that holds the coordinates of the
      interpolation points relative to XBASE. */
  double *xpt;
  /** Size of array XPT */
  int xpt_size;
  /** FVAL holds the values of F at the interpolation points. */
  double *fval;
  /** Size of array FVAL */
  int fval_size;
  /** All the components of every XOPT are going to satisfy the bounds
      SL(I) .LEQ. XOPT(I) .LEQ. SU(I), with appropriate equalities when
      XOPT is on a constraint boundary. */
  double *xopt;
  /** Size of array XOPT */
  int xopt_size;
  /** GOPT holds the gradient of the quadratic model at XBASE+XOPT. */
  double *gopt;
  /** Size of array GOPT */
  int gopt_size;
  /** HQ holds the explicit second derivatives of the quadratic model. */
  double *hq;
  /** Size of array HQ */
  int hq_size;
  /** PQ contains the parameters of the implicit second derivatives of the
      quadratic model. */
  double *pq;
  /** Size of array PQ */
  int pq_size;
  /** BMAT holds the last N columns of H. */
  double *bmat;
  /** Size of array BMAT */
  int bmat_size;
  /** ZMAT holds the factorization of the leading NPT by NPT submatrix of H,
      this factorization being ZMAT times ZMAT^T, which provides both the
      correct rank and positive semi-definiteness. */
  double *zmat; 
  /** Size of array ZMAT */
  int zmat_size; 
  /** NDIM is the first dimension of BMAT and has the value NPT+N. */
  int ndim;
  /** SL and SU hold the differences XL-XBASE and XU-XBASE, respectively. */
  double *sl;
  /** Size of array SL */
  int sl_size;
  /** SL and SU hold the differences XL-XBASE and XU-XBASE, respectively. */
  double *su;
  /** Size of array SU */
  int su_size;
  /** XNEW is chosen by SUBROUTINE TRSBOX or ALTMOV. Usually XBASE+XNEW is the
      vector of variables for the next call of CALFUN. XNEW also satisfies
      the SL and SU constraints in the way that has just been mentioned. */
  double *xnew; 
  /** Size of array XNEW */
  int xnew_size; 
  /** XALT is an alternative to XNEW, chosen by ALTMOV, that may replace XNEW
      in order to increase the denominator in the updating of UPDATE. */
  double *xalt;
  /** Size of array XALT */
  int xalt_size;
  /** DTRIAL is reserved for a trial step from XOPT, usually XNEW-XOPT. */
  double *dtrial;
  /** Size of array DTRIAL */
  int dtrial_size;
  /** VLAG contains the values of the Lagrange functions at a new point X.
   *  They are part of a product that requires VLAG to be of length NDIM. */
  double *vlag;
  /** Size of array VLAG */
  int vlag_size;

  /** W2NPT is working memory of size 2*NPT */
  double *w2npt;
  /** Size of array W2NPT */
  int w2npt_size;
  /** Working memory of size NDIM */
  double *wndim;
  /** Size of array WNDIM */
  int wndim_size;
  /** WN is working memory of size N */
  double *wn;
  /** Size of array WN */
  int wn_size;
  /** GNEW holds the gradient of the quadratic model at XOPT+D. It is updated
      when DTRIAL is updated. Size N. */
  double *gnew;
  /** Size of array GNEW */
  int gnew_size;
  /** XBDI is a working space vector of size N.
   *  For I=1,2,...,N, the element XBDI(I) is
      set to -1.0, 0.0, or 1.0, the value being nonzero if and only if the
      I-th variable has become fixed at a bound, the bound being SL(I) or
      SU(I) in the case XBDI(I)=-1.0 or XBDI(I)=1.0, respectively. This
      information is accumulated during the construction of XNEW. */
  double *xbdi;
  /** Size of array XBDI */
  int xbdi_size;
  /** The arrays S, HS and HRED are also used for working space, all of size N.
      They hold the current search direction, and the changes in the gradient
      of Q along S and the reduced D, respectively, where the reduced D is the
      same as D, except that the components of the fixed variables are zero. */
  double *s;
  /** Size of array S */
  int s_size;
  /** The arrays S, HS and HRED are also used for working space, all of size N.
      They hold the current search direction, and the changes in the gradient
      of Q along S and the reduced D, respectively, where the reduced D is the
      same as D, except that the components of the fixed variables are zero. */
  double *hs;
  /** Size of array HS */
  int hs_size;
  /** The arrays S, HS and HRED are also used for working space, all of size N.
      They hold the current search direction, and the changes in the gradient
      of Q along S and the reduced D, respectively, where the reduced D is the
      same as D, except that the components of the fixed variables are zero. */
  double *hred;
  /** Size of array HRED */
  int hred_size;
  /** GLAG is a vector of length N for the gradient of the KNEW-th Lagrange
   *  function at XOPT used in ALTMOV. */
  double *glag;
  /** Size of array GLAG */
  int glag_size;
  /** HCOL is a vector of length NPT for the second derivative coefficients
   *  of the KNEW-th Lagrange function used in ALTMOV. */
  double *hcol;
  /** Size of array HCOL */
  int hcol_size;
  /** W is a working space vector of length 2N that is going to hold the
      constrained Cauchy step from XOPT of the Lagrange function, followed
      by the downhill version of XALT when the uphill step is calculated
      in ALTMOV. */
  double *ccstep;
  /** Size of array CCSTEP */
  int ccstep_size;

  /** Verbose level; if zero, then nothing is printed to stderr or stdout */
  int verbose;

  /** Variables for exchange inside BOBYQA */

/// @cond
  /** CRVMIN is set to zero in TRSBOX if D reaches the trust region boundary.
   *  Otherwise it is set to the least curvature of H that occurs in
   *  the conjugate gradient searches that are not restricted by any
   *  constraints. The value CRVMIN=-1.0 is set, however, if all of these
   *  searches are constrained. */
  double _crvmin;
/// @endcond
  /** Extra variable defined in bobyqb */
  int ntrits;
  /** Extra variable defined in bobyqb */
  double rho;
  /** Extra variable defined in bobyqb */
  int nresc;
  /** Extra variable defined in bobyqb */
  double delta;
  /** Extra variable defined in bobyqb */
  double diffa;
  /** Extra variable defined in bobyqb */
  double diffb;
  /** Extra variable defined in bobyqb */
  double diffc;
  /** Extra variable defined in bobyqb */
  double ratio;
  /** Extra variable defined in bobyqb */
  int itest;
  /** Extra variable defined in bobyqb */
  int nfsav;
  /** Extra variable defined in bobyqb */
  int kopt;
  /** Extra variable defined in bobyqb */
  double fsave;
  /** Extra variable defined in bobyqb */
  double vquad;
  /** Extra variable defined in bobyqb */
  double fopt;
  /** Extra variable defined in bobyqb */
  double dsq;
  /** Extra variable defined in bobyqb */
  double xoptsq;
  /** Extra variable defined in bobyqb */
  int nptm;
  /** Extra variable defined in bobyqb */
  double alpha;
  /** Extra variable defined in bobyqb */
  double beta;
  /** Extra variable defined in bobyqb */
  double dnorm;
  /** Return code */
  int rc;
  /** Extra variable defined in bobyqb */
  double newf;
  /** Extra variable defined in bobyqb */
  int knew;
  /** Extra variable defined in bobyqb */
  int kbase;
  /** Extra variable defined in bobyqb */
  double denom;
  /** Extra variable defined in bobyqb */
  double delsq;
  /** Extra variable defined in bobyqb */
  double scaden;
  /** Extra variable defined in bobyqb */
  double biglsq;
  /** Extra variable defined in bobyqb */
  double distsq;
  /** Extra variable defined in bobyqb */
  double cauchy;
  /** Extra variable defined in bobyqb */
  double adelt;
  /** Nr of subfunction calls */
  int prelim_nr;
  /** Nr of subfunction calls */
  int rescue_nr;
  /** Nr of subfunction calls */
  int altmov_nr;
  /** Nr of subfunction calls */
  int trsbox_nr;
  /** Nr of subfunction calls */
  int update_nr;
} bobyqa_data;
/*****************************************************************************/

/*****************************************************************************/
/* aic */
double aicSS(double ss, const int n, const int k);
int parFreeNr(const int n, double *pLower, double *pUpper);
int aicWeights(double *aic, double *w, int n);
double aicWeightedAvg(double *w, double *p, int n);
double aicModel(double *w, int n);
/*****************************************************************************/

/*****************************************************************************/
/* bobyqa */
bobyqa_result bobyqb(
  bobyqa_data *bdata
);
bobyqa_result bobyqa(
  int n,
  int npt,
  double *x,   
  const double *xl,
  const double *xu, 
  const double *dx,
  const double rhoend,
  double xtol_rel,
  double minf_max,
  double ftol_rel,
  double ftol_abs,
  int maxeval,
  int *nevals,
  double *minf,
  double (*f)(int n, double *x, void *objf_data),
  void *objf_data,
  double *working_space,
  int verbose
);
int bobyqa_minimize_single_parameter(
  bobyqa_data *bdata
);
char *bobyqa_rc(
  bobyqa_result rc
);
int fixed_params(
  int n, const double *lower, const double *upper, const double *delta
);
int bobyqa_working_memory_size(
  int n, int fitted_n, int npt, bobyqa_data *bdata
);
bobyqa_result bobyqa_set_memory(
  int n, int fitted_n, int npt, bobyqa_data *bdata, double *wm
);
bobyqa_result bobyqa_free_memory(
  bobyqa_data *bdata
);
bobyqa_result bobyqa_reset_memory(
  bobyqa_data *bdata
);
void bobyqa_print(
  bobyqa_data *bdata, int sw, FILE *fp
);
double bobyqa_x_funcval(
  bobyqa_data *bdata, double *x
);
void bobyqa_xfull(
  bobyqa_data *bdata
);
bobyqa_result bobyqa_set_optimization(
  int full_n,
  double *x,
  const double *dx,
  const double *xl,
  const double *xu,
  const double rhoend,
  double xtol_rel,
  double minf_max,
  double ftol_rel,
  double ftol_abs,
  int maxeval,
  double (*f)(int n, double *x, void *objf_data),
  void *objf_data,
  int verbose,
  bobyqa_data *bdata
);
/*****************************************************************************/

/*****************************************************************************/
/* bootstrap */
int bootstrap(
  int iterNr,
  double *cLim1, double *cLim2, double *SD, double *parameter,
  double *lowlim, double *uplim, int frameNr, double *origTac,
  double *fitTac, double *bsTAC, int parNr, double *weight,
  double (*objf)(int, double*, void*), char *status, int verbose
);
int bootstrapr(
  int iterNr,
  double *cLim1, double *cLim2, double *SD, double *parameter,
  double *lowlim, double *uplim, int frameNr, double *origTac,
  double *fitTac, double *bsTAC, int parNr, double *weight,
  double (*objf)(int, double*, void*), char *status, int verbose, double *matrix
);
/*****************************************************************************/

/*****************************************************************************/
/* bvls */
int bvls(
  int key, const /*unsigned*/ int m, const /*unsigned*/ int n,
  double *a, double *b, double *bl, double *bu, double *x,
  double *w, double *act, double *zz, int *istate, int *iter,
  int verbose
);
int llsqWght(int N, int M, double **A, double *a, double *b, double *weight);
int llsqWghtSquared(int N, int M, double **A, double *a, double *b, double *weight);
/*****************************************************************************/

/*****************************************************************************/
/* constraints */
int modelCheckParameters(
  int par_nr, double *lower_p, double *upper_p, double *test_p,
  double *accept_p, double *penalty
);
int modelCheckLimits(
  int par_nr, double *lower_p, double *upper_p, double *test_p
);
int fitExpDecayNNLS(
  double *x, double *y, int n, double fittime, double kmin, double kmax,
  int pnr, double *a, double *k, int *fnr, int verbose
);
/*****************************************************************************/

/*****************************************************************************/
/* gaussdev */

/** Seed for random number generator */
long int GAUSSDEV_SEED;

unsigned int drandSeed(short int seed);
double gaussdev();
double gaussdev2();
void init_gaussdev();
double drand();
int rand_range(
  int nr, double *d, double low, double up, int type
);
/*****************************************************************************/

/*****************************************************************************/
/* hholder */
double householder_transform(double *v, int N);
int householder_hm(
  double tau, double *vector, double **matrix, int rowNr, int columnNr
);
int householder_hv(
  double tau, int size, double *v, double *w
);
double householder_norm(double *v, int size);
/*****************************************************************************/

/*****************************************************************************/
/* integr */

/** Verbose prints from integr.c */
int INTEGR_TEST;

int interpolate(
  double *x, double *y, int nr,
  double *newx, double *newy, double *newyi, double *newyii, int newnr
);
int finterpolate(
  float *x, float *y, int nr,
  float *newx, float *newy, float *newyi, float *newyii, int newnr
);
int integrate(
  double *x, double *y, int nr,
  double *yi
);
int fintegrate(
  float *x, float *y, int nr,
  float *yi
);
int petintegrate(
  double *x1, double *x2, double *y, int nr, 
  double *newyi, double *newyii
);
int fpetintegrate(
  float *x1, float *x2, float *y, int nr, 
  float *newyi, float *newyii
);
int interpolate4pet(
  double *x, double *y, int nr,
  double *newx1, double *newx2, 
  double *newy, double *newyi, double *newyii, int newnr
);
int finterpolate4pet(
  float *x, float *y, int nr,
  float *newx1, float *newx2, 
  float *newy, float *newyi, float *newyii, int newnr
);
int petintegral(
  double *x1, double *x2, double *y, int nr,   
  double *ie, double *iie
);
int fpetintegral(
  float *x1, float *x2, float *y, int nr,   
  float *ie, float *iie
);
int petintegrate2fe(
  double *x1, double *x2, double *y, int nr,    
  double *e, double *ie, double *iie
);
int fpetintegrate2fe(
  float *x1, float *x2, float *y, int nr,    
  float *e, float *ie, float *iie
);
/*****************************************************************************/

/*****************************************************************************/
/* llsqwt */

/** Verbose prints from LLSQWT */
int LLSQWT_TEST;

int llsqwt(
  double *x, double *y, int n, double *wx, double *wy, double tol, double *w,
  double *ic, double *slope, double *nwss, double *sic, double *sslope,
  double *cx, double *cy
);
int best_llsqwt(
  double *x, double *y, double *wx, double *wy, int nr, int min_nr, int mode,
  double *slope, double *ic, double *nwss, double *sslope, double *sic,
  double *cx, double *cy, int *bnr
);
int llsqperp(double *x, double *y, int nr,
  double *slope, double *ic, double *ssd
);
int llsqperp3(double *x, double *y, int nr,
  double *slope, double *ic, double *ssd
);
int quadratic(double a, double b, double c, double *m1, double *m2);
int medianline(double *x, double *y, int nr, double *slope, double *ic);
/*****************************************************************************/

/*****************************************************************************/
/* lms */
double least_median_of_squares(
  double *data, int n
);
/*****************************************************************************/

/*****************************************************************************/
/* lts */
/// @cond
#ifndef CHI2INV_1
#define CHI2INV_1 0.45493642311957
#endif
/// @endcond

int least_trimmed_square(
  double data[], long int n, double *mean, double *variance
);
/*****************************************************************************/

/*****************************************************************************/
/* median */
double d_kth_smallest(double *data, int n, int k);
double dmedian(double *data, int n);
double dmean(double *data, int n, double *sd);
double dmean_nan(double *data, int n, double *sd, int *vn);
/*****************************************************************************/

/*****************************************************************************/
/* mestim */
double mEstim(
  double *data, int nr, int iterNr, double cutoff
);
double huber(
  double x, double b
);
/*****************************************************************************/

/*****************************************************************************/
/* mtga */

/** Min nr of points in MTGA line fit */
#define MTGA_BEST_MIN_NR 5

/** MTGA line fit method */
typedef enum {
  /** Traditional line fit, same as Pearson's correlation coefficient */
  PEARSON, 
  /** Simple non-iterative perpendicular line fitting (Varga & Szabo, 2002) */
  PERP,
  /** Iterative method for linear least-squares fit with errors in both
      coordinates (York 1966, Lybanon 1984, Reed 1992). */
  LLSQWT,
  /** Median-based distribution-free estimation of slope and intercept
   *  (Siegel, 1982); note that this is not LMS */
  MEDIAN
} linefit_method;

/** MTGA line fit range */
typedef enum {
  PRESET, EXCLUDE_BEGIN, EXCLUDE_END
} linefit_range;

int patlak_data(
  int data_nr, double *i, double *ii, double *c, double *x, double *y
);
int logan_data(
  int data_nr, double *i, double *ii, double *c, double *ci,
  double k2, double *x, double *y
);
int mtga_best_perp(
  double *x, double *y, int nr,
  double *slope, double *ic, double *ssd, int *fnr
);
/*****************************************************************************/

/*****************************************************************************/
/* nnls */
int nnls(
  double **a, int m, int n, double *b, double *x,
  double *rnorm, double *w, double *zz, int *index
);
int nnlsWght(
  int N, int M, double **A, double *b, double *weight
);
int nnlsWghtSquared(
  int N, int M, double **A, double *b, double *sweight
);
/*****************************************************************************/

/*****************************************************************************/
/* normaldistr */
double ndtr(double a);
double normal_pvalue_2(double x);
double normal_pvalue_1(double x);
/*****************************************************************************/

/*****************************************************************************/
/* pearson */
/** Verbose prints from pearson */
int PEARSON_TEST;

int pearson(
  double *x, double *y, int nr,
  double *k, double *kSD, double *b, double *bSD, double *r, double *ySD
);
int pearson2(
  double *x, double *y, char *is, int nr,
  double *k, double *kSD, double *b, double *bSD, double *r, double *ySD
);
int pearson3(
  double *x, double *y, int nr,
  double *k, double *kSD, double *b, double *bSD, double *r, double *ySD
);
int pearson4(
  double *x, double *y, int nr, double start, double end,
  double *k, double *kSD, double *b, double *bSD, double *r, double *ySD
);
int best_pearson(
  double *x, double *y, int nr, int min_nr, int *first, int *last,
  double *k, double *kSD, double *b, double *bSD, double *r, double *ySD
);
int mean(
  double *x, double *y, int nr,
  double *xmean, double *xsd, double *ymean, double *ysd
);
int regr_line(
  double *x, double *y, int n, double *m, double *c
);
int highest_slope(
  double *x, double *y, int n, int slope_n,
  double *m, double *c, double *xi, double *xh
);
int highest_slope_after(
  double *x, double *y, int n, int slope_n,
  double x_start, double *m, double *c, double *xi, double *xh
);
/*****************************************************************************/

/*****************************************************************************/
/* powell */
/// @cond
extern int POWELL_LINMIN_MAXIT;
/// @endcond

int powell(double *p, double *delta, int parNr, double ftol, int *iterNr,
  double *fret, double (*_fun)(int, double*, void*), void *fundata, int verbose
);
/*****************************************************************************/

/*****************************************************************************/
/* qr */
int qr(
  double **A, int m, int n, double *B, double *X, double *rnorm,
  double *tau, double *res, double **wws, double *ws
);
int qr_decomp(
  double **a, int M, int N, double *tau, double **cchain, double *chain
);
int qr_solve(
  double **QR, int M, int N, double *tau, double *b, double *x,
  double *residual, double *resNorm, double **cchain, double *chain
);
int qr_weight(
  int N, int M, double **A, double *b, double *weight, double *ws
);
int qrLH(
  const unsigned int m, const unsigned int n,
  double *a, double *b, double *x, double *r2
);
/*****************************************************************************/

/*****************************************************************************/
/* runs_test */
int runs_test(
  double *data1, double *data2, int N, double alpha, double *p
);
void residuals(
  double *d1, double *d2, int Nr, int *Rnr, int *Nminus, int *Nplus
);
int mrl_between_tacs(
  double y1[], double y2[], int n
);
/*****************************************************************************/

/*****************************************************************************/
/* shuffle */
void random_shuffle(int *array, int n);
void randperm(int *array, int n, int a);
/*****************************************************************************/

/*****************************************************************************/
/* simplex */
/** Verbose prints from simplex */
int SIMPLEX_TEST;

double simplex(
  double (*_fun)(double*),
  int parNr,
  double *par,
  double *delta,
  double maxerr,
  int maxiter
);
/*****************************************************************************/

/*****************************************************************************/
/* simulate */
/** Verbose prints from simulation functions */
int SIMULATE_TEST;
/*****************************************************************************/
int simC3s(
  double *t, double *ca, int nr, double k1, double k2,
  double k3, double k4, double k5, double k6,
  double *ct, double *cta, double *ctb, double *ctc
);
int simC3p(
  double *t, double *ca, int nr, double k1, double k2,
  double k3, double k4, double k5, double k6,
  double *ct, double *cta, double *ctb, double *ctc
);
int simC3vs(
  double *t, double *ca, double *cb, int nr,
  double k1, double k2, double k3, double k4, double k5, double k6,
  double f, double vb, double fa,
  double *cpet, double *cta, double *ctb, double *ctc,
  double *ctab, double *ctvb
);
int simC3vp(
  double *t, double *ca, double *cb, int nr,
  double k1, double k2, double k3, double k4, double k5, double k6,
  double f, double vb, double fa,
  double *cpet, double *cta, double *ctb, double *ctc,
  double *ctab, double *ctvb
);
int simC2l(
  double *t, double *ca, int nr, double k1, double k2,
  double k3, double kLoss, double *ct, double *cta, double *ctb
);
int simC2vl(
  double *t, double *ca, double *cb, int nr,
  double k1, double k2, double k3, double kL,
  double f, double vb, double fa,
  double *cpet, double *cta, double *ctb,
  double *ctab, double *ctvb
);
int simC3vpKLoss(
  double *t, double *ca, double *cb, int nr,
  double k1, double k2, double k3, double k4, double k5, double k6,
  double kLoss, double f, double vb, double fa,
  double *cpet, double *cta, double *ctb, double *ctc,
  double *ctab, double *ctvb
);
int simRTCM(
  double *t, double *cr, int nr, double R1, double k2,
  double k3, double k4, double *ct, double *cta, double *ctb
);
int simSRTM(
  double *t, double *cr, int nr, double R1, double k2,
  double BP, double *ct
);
int simTRTM(
  double *t, double *cr, int nr, double R1, double k2,
  double k3, double *ct
);
int simHuangmet(
  double *t, double *ctot, int nr,
  double k01, double k12, double k21, double k03, double k34, double k43,
  double *c0, double *c1, double *c3 
);
int simTPCMOD0009c(
  double *t, double *ctot, int nr, double km, double k1m,
  double k2m, double k3m, double k4m, double *ca, double *cm
);
int simMBF(
  double *t, double *ci, int nr, double k1, double k2, double Vfit, double *ct
);
int simC1(
  double *t, double *ca, int nr, double k1, double k2, double *ct
);
int simC3DIvs(
  double *t, double *ca1, double *ca2, double *cb, int nr,
  double k1, double k2, double k3, double k4, double k5, double k6,
  double k1b, double k2b, double f, double vb, double fa,
  double *scpet, double *sct1, double *sct2, double *sct3, double *sct1b,
  double *sctab, double *sctvb
);
int simC4DIvp(
  double *t, double *ca1, double *ca2, double *cb, int nr,
  double k1, double k2, double k3, double k4, double k5, double k6, double k7,
  double km, double k1b, double k2b, double f, double vb, double fa,
  double *scpet, double *sct1, double *sct2, double *sct3,
  double *sct1b, double *sctab, double *sctvb,
  int verbose
);
int simC4DIvs(
  double *t, double *ca1, double *ca2, double *cb, int nr,
  double k1, double k2, double k3, double k4, double k5, double k6, double k7,
  double km, double k1b, double k2b, double f, double vb, double fa,
  double *scpet, double *sct1, double *sct2, double *sct3,
  double *sct1b, double *sctab, double *sctvb,
  int verbose
);
int simDispersion(
  double *x, double *y, int n,
  double tau1, double tau2, double *tmp
);
int simOxygen(
  double *t, double *ca1, double *ca2, double *ca1i, double *ca2i, const int n,
  const double k1a, const double k2a, const double km, 
  const double k1b, const double k2b, 
  const double vb, const double fa,
  double *scpet, double *sct1, double *sct2, double *sctab, 
  double *sctvb1, double *sctvb2, double *scvb1, double *scvb2,
  const int verbose
);


/// @cond
/* Deprecated functions. Please don't use these any more */
#define autointerpolateDFT dftAutointerpolate
#define c3sSIM simC3s
#define c3pSIM simC3p
#define c3vsSIM simC3vs
#define c3vpSIM simC3vp
#define c2lSIM simC2l
#define c2vlSIM simC2vl
#define rtcmSIM simRTCM
#define srtmSIM simSRTM
#define trtmSIM simTRTM
#define sim_dispersion simDispersion
/// @endcond
/*****************************************************************************/

/*****************************************************************************/
/* o2 */

/** Arterial oxygen saturation fraction */
#define DEFAULT_SAO2 0.97
/** Half-saturation pressure p50 (kPa) for hemoglobin */
#define DEFAULT_P50HB 3.6
/** Half-saturation pressure p50 (kPa) p50 for myoglobin */
#define DEFAULT_P50MB 0.319
/** Hill coefficient n for hemoglobin */
#define DEFAULT_NHB 2.7
/** Hemoglobin concentration in blood (mg/g) */
#define DEFAULT_CHB 150.0
/** Myoglobin concentration in muscle (mg/g) */
#define DEFAULT_CMB 4.7

double mo2k1k2(
  const double OER, const double SaO2, const double p50Hb, const double p50Mb,
  const double nHb, const double cHb, const double cMb,
  const int verbose
);
double mo2pO2(
  const double OER, const double K1k2, const double SaO2, const double p50Mb,
  const double cHb, const double cMb,
  const int verbose
);
/*****************************************************************************/

/*****************************************************************************/
/* tgo */
/** Biased (1) or even (0) parameter distribution */
extern int TGO_SQUARED_TRANSF;
/** Local optimization outside (0) or inside (1) iTGO */
extern int TGO_LOCAL_INSIDE;
/** Local optimization method is Powell-Brent (0) or Bobyqa (1) */
extern int TGO_LOCAL_OPT;

/** TGO point */
typedef struct {
  /** index of min */
  int topomin;
  /** min values */
  double fvalue;
  /** paramaters at min */
  double par[MAX_PARAMS];
  /** parameter deltas */
  double delta[MAX_PARAMS]; // added by VO 2011-11-25
  /** parameter ranges */
  double fvalrange; // added by VO 2011-11-25
} TGO_POINT;

int tgo(
  double *lowlim, double *uplim, double (*objf)(int, double*, void*),
  void *objfData, int dim, int neighNr, double *fmin, double *gmin,
  int samNr, int tgoNr, int verbose
);    
void tgoRandomParameters(
  TGO_POINT *p, int parNr, int sNr,
  double *low, double *up
);  
void tgoRandomParametersST(
  TGO_POINT *p, int parNr, int sNr,
  double *low, double *up
);
/*****************************************************************************/

/*****************************************************************************/
/* nlopt1d */
int nlopt1D(
  double (*_fun)(double, void*), void *_fundata,
  double x, double xl, double xu, double delta, double tol, const int maxeval,
  double *nx, double *nf, int verbose
);
/*****************************************************************************/

/*****************************************************************************/
#endif // _LIBTPCMODEL_H
