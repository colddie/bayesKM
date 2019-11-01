#ifdef __cplusplus
extern "C"
{
#endif

// C header here

#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <string.h>
#include <dlfcn.h>
#include <time.h>
// #include <iostream>

#ifdef __cplusplus
}
#endif


#ifndef MAX_PARAMETERS
/** Max nr of parameters */
#define MAX_PARAMETERS 50
#endif
#ifndef MAX_PARAMS
/** Max nr of parameters */
#define MAX_PARAMS MAX_PARAMETERS
#endif

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

extern "C" int tgo(
  double *lowlim, double *uplim, double (*objf)(int, double*, void*),
  void *objfData, int dim, int neighNr, double *fmin, double *gmin,
  int samNr, int tgoNr, int verbose
);    
extern "C" void tgoRandomParameters(
  TGO_POINT *p, int parNr, int sNr,
  double *low, double *up
);  
extern "C" void tgoRandomParametersST(
  TGO_POINT *p, int parNr, int sNr,
  double *low, double *up
);


/* gaussdev */

/** Seed for random number generator */
// long int GAUSSDEV_SEED;

extern "C" unsigned int drandSeed(short int seed);
extern "C" double gaussdev();
extern "C" double gaussdev2();
extern "C" void init_gaussdev();
extern "C" double drand();
extern "C" int rand_range(
  int nr, double *d, double low, double up, int type
);


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

/* bobyqa */
extern "C"  bobyqa_result bobyqb(
  bobyqa_data *bdata
);
extern "C"  bobyqa_result bobyqa(
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
extern "C"  int bobyqa_minimize_single_parameter(
  bobyqa_data *bdata
);
extern "C"  char *bobyqa_rc(
  bobyqa_result rc
);
extern "C"  int fixed_params(
  int n, const double *lower, const double *upper, const double *delta
);
extern "C"  int bobyqa_working_memory_size(
  int n, int fitted_n, int npt, bobyqa_data *bdata
);
extern "C"  bobyqa_result bobyqa_set_memory(
  int n, int fitted_n, int npt, bobyqa_data *bdata, double *wm
);
extern "C"  bobyqa_result bobyqa_free_memory(
  bobyqa_data *bdata
);
extern "C"  bobyqa_result bobyqa_reset_memory(
  bobyqa_data *bdata
);
extern "C"  void bobyqa_print(
  bobyqa_data *bdata, int sw, FILE *fp
);
extern "C"  double bobyqa_x_funcval(
  bobyqa_data *bdata, double *x
);
extern "C"  void bobyqa_xfull(
  bobyqa_data *bdata
);
extern "C"  bobyqa_result bobyqa_set_optimization(
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



/* powell */
/// @cond
extern int POWELL_LINMIN_MAXIT;
/// @endcond

extern "C" int powell(double *p, double *delta, int parNr, double ftol, int *iterNr,
  double *fret, double (*_fun)(int, double*, void*), void *fundata, int verbose
);



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



/* bootstrap */
int bootstrapr(
  int iterNr,
  double *cLim1, double *cLim2, double *SD, double *parameter,
  double *lowlim, double *uplim, int frameNr, double *origTac,
  double *fitTac, double *bsTAC, int parNr, double *weight,
  double (*objf)(int, double*, void*), char *status, int verbose, double *matrix
);

int temp_roundf(float e);