/// @file powell.c
/// @author Vesa Oikonen
/// @brief Powell function minimization routines.
///
///  Based on Numerical recipes in C (Press et al.).
///
/******************************************************************************/
#include "libtpcmodel.h"
/******************************************************************************/
/** Max iterations for linear minimization inside powell() */
int POWELL_LINMIN_MAXIT=100;
/******************************************************************************/
/* Local variables */
/// @cond
int _powell_ncom;
double _powell_pcom[MAX_PARAMETERS], _powell_xicom[MAX_PARAMETERS];
double (*_powellFunc)(int, double*, void*);
void *_powellFuncData;
int _powell_func_calls;
/******************************************************************************/
/* Local functions */
void _powell_linmin(double *p, double *xi, int n, double *fret, int *itnr);
double _powell_brent(double ax, double bx, double cx, double tol, double *xmin,
      int *itnr, int dim);
double _powell_f1dim(double x, int dim);
void _powell_mnbrak(double *ax, double *bx, double *cx, double *fa, double *fb,
      double *fc, int dim);
/* Local "inline" functions */
double _powell_sqr(double x) {return (x*x);}
void _powell_shft(double *a, double *b, double *c, double *d) {*a=*b; *b=*c; *c=*d;}
double _powell_fmax(double a, double b) {return((a>b) ? a:b);}
/// @endcond
/******************************************************************************/

/******************************************************************************/
/** Powell function minimization routine.
    @return Returns 0, if succesful, 1 if required tolerance was not reached, 
    2 if initial guess does not give finite function value,
    3 if final function value is NaN or infinite,
    and >3 in case of another error.
    @sa simplex, tgo, bobyqa, nlopt1D
 */
int powell(
  /** Initial guess and final set of parameters */
  double *p,
  /** Initial changes for parameters, ==0 if fixed */
  double *delta,
  /** Nr of parameters */
  int parNr,
  /** Fractional tolerance (for WSS); 0<ftol<1 */
  double ftol,
  /** Max nr of iterations, and nr of required iters */
  int *iterNr,
  /** Function return value (WSS) at minimum */
  double *fret,
  /** Function to minimize (must return the WSS) */
  double (*_fun)(int, double*, void*),
  /** Pointer to data which is passed on to the function; NULL if not needed */
  void *fundata,
  /** Verbose level; if zero, then nothing is printed into stdout or stderr */
  int verbose
) {
  int pbIterNr;
  int i, j, ibig, iterMax, fixed[MAX_PARAMETERS];
  double xi[MAX_PARAMETERS][MAX_PARAMETERS]; /* Matrix for directions */
  double pt[MAX_PARAMETERS], ptt[MAX_PARAMETERS], xit[MAX_PARAMETERS];
  double del, fp, fptt, t;
  double origp[MAX_PARAMETERS];
  int ftol_reached=0;


  if(verbose>0) printf("in powell(,,%d,%g,%d,,)\n", parNr, ftol, *iterNr);
  if(verbose>1) {
    printf("Initial parameter guesses and deltas:\n");
    for(i=0; i<parNr; i++) printf("  %g  %g\n", p[i], delta[i]);
  }
  *fret=nan("");
  if(p==NULL) return(11);
  if(delta==NULL) return(12);
  if(parNr<1) return(21);
  if(ftol<=0.0) return(22);
  if(ftol>=1.0) return(23);
  if((*iterNr)<1) return(24);

  /* SetUp */
  _powellFunc=_fun;
  _powellFuncData=fundata;
  iterMax=*iterNr; /* save the max nr of iterations */
  _powell_ncom=parNr;
  /* Function value at initial point */
  _powell_func_calls=1;
  *fret=(*_powellFunc)(parNr, p, _powellFuncData);
  if(verbose>10) printf("initial point fret=%g\n", *fret);
  if(!isfinite(*fret)) {
    if(verbose>0) printf("in powell(): objf failed at initial point.\n");
    *fret=nan(""); return(2);
  }
  /* Save the initial point (pt[] will be changed later) */
  for(j=0; j<parNr; j++) origp[j]=pt[j]=p[j];
  /* Check which parameters are fixed */
  for(i=0; i<parNr; i++) if(fabs(delta[i])<1.0e-20) fixed[i]=1; else fixed[i]=0;

  /* Initiate matrix for directions */
  for(i=0; i<parNr; i++)
    for(j=0; j<parNr; j++)
      if(i==j) xi[i][j]=delta[i]; else xi[i][j]=0.0;

  /* Iterate */
  for(*iterNr=1; ; (*iterNr)++) {
    if(verbose>2) printf("  iteration %d\n", *iterNr);
    fp=*fret; ibig=0; del=0.0; /* largest function decrease */

    /* In each iteration, loop over all directions in the set */
    for(i=0; i<parNr; i++) {
      if(fixed[i]) continue; /* do nothing with fixed parameters */
      for(j=0; j<parNr; j++) if(fixed[j]) xit[j]=0.0; else xit[j]=xi[j][i];
      fptt=*fret;
      /* minimize along direction xit */
      pbIterNr=POWELL_LINMIN_MAXIT;
      _powell_linmin(p, xit, parNr, fret, &pbIterNr);
      if(verbose>3) printf("iterNr in _powell_linmin() with p%d: %d\n",
                           i, pbIterNr);
      if(fabs(fptt-(*fret))>del) {del=fabs(fptt-(*fret)); ibig=i;}
    }
    if(verbose>20) printf("fret=%g  fp=%g\n", *fret, fp);

    /* Check if done */
#if(0)
    if(2.0*fabs(fp-(*fret)) <= ftol*(fabs(fp)+fabs(*fret))) break;
#else
    if(2.0*fabs(fp-(*fret)) <= ftol*(fabs(fp)+fabs(*fret))) {
      if(ftol_reached>0 || (*iterNr)>=iterMax) break; else ftol_reached++;
    } else ftol_reached=0;
#endif
    if((*iterNr)>=iterMax) {
      if(verbose>0) printf("max iterations nr exceeded in powell().\n");
      break;
    }
    /* Construct the extrapolated point and the average direction moved */
    for(j=0; j<parNr; j++) {
      ptt[j]=2.0*p[j]-pt[j]; xit[j]=p[j]-pt[j];
      pt[j]=p[j]; /* save the old starting point */
    }
    fptt=(*_powellFunc)(parNr, ptt, _powellFuncData); _powell_func_calls++;
    if(fptt<fp) {
      t=2.0*(fp-2.0*(*fret)+fptt)*_powell_sqr(fp-(*fret)-del)-del*_powell_sqr(fp-fptt);
      if(t<0.0) {
        pbIterNr=POWELL_LINMIN_MAXIT;
        _powell_linmin(p, xit, parNr, fret, &pbIterNr);
        if(verbose>3) printf("iterNr in _powell_linmin(): %d\n", pbIterNr);
        for(j=0; j<parNr; j++) {
          xi[j][ibig]=xi[j][parNr-1]; xi[j][parNr-1]=xit[j];}
      }
    }
  } /* next iteration */
  if(verbose>1) {
    printf("iterNr := %d\n", *iterNr);
    printf("nr of function calls := %d\n", _powell_func_calls);
  }

  if(isnan(*fret) || !isfinite(*fret)) {
    if(verbose>10) printf("powell() fails and returns the initial point.\n");
    if(verbose>11) {
      if(isnan(*fret)) printf("  fret := NaN\n");
      if(isfinite(*fret)) printf("  fret := overflow/underflow\n");
    }
    // if failed, then return initial guess
    for(j=0; j<parNr; j++) p[j]=origp[j];
    // and call function again so that any data saved there is correct
    *fret=(*_powellFunc)(parNr, p, _powellFuncData);
    return(3);
  }
  // and call function again so that any data saved there is correct
  *fret=(*_powellFunc)(parNr, p, _powellFuncData);
  if((*iterNr)>=iterMax) return(1);
  if(verbose>0) printf("out of powell() in good order.\n");
  return(0);
}
/******************************************************************************/

/******************************************************************************/
/// @cond
void _powell_linmin(double *p, double *xi, int n, double *fret, int *itnr)
{
  int i;
  double xx, xmin, fx, fb, fa, bx, ax;

  _powell_ncom=n;
  for(i=0; i<n; i++) {_powell_pcom[i]=p[i]; _powell_xicom[i]=xi[i];}
  ax=0.0; xx=1.0;
  _powell_mnbrak(&ax, &xx, &bx, &fa, &fx, &fb, n);
  *fret=_powell_brent(ax, xx, bx, 2.0e-4, &xmin, itnr, n);
  for(i=0; i<n; i++) {xi[i]*=xmin; p[i]+=xi[i];}
}
/******************************************************************************/
double _powell_brent(
  double ax, double bx, double cx, double tol, double *xmin, int *itnr, int dim
) {
  //const int ITMAX = 100;
  const double CGOLD = 0.3819660;
  const double ZEPS = 1.0E-10;
  int iterMax;
  double a, b, d=0.0, etemp, fu, fv, fw, fx, p, q, r;
  double e=0.0, tol1, tol2, u, v, w, x, xm;

  a=(ax<cx ? ax:cx); b=(ax>cx ? ax:cx); x=w=v=bx;
  fw=fv=fx=_powell_f1dim(x, dim);
  iterMax=*itnr;
  for(*itnr=0; *itnr<iterMax; (*itnr)++) {
    xm=0.5*(a+b); tol2=2.0*(tol1=tol*fabs(x)+ZEPS);
    if(fabs(x-xm)<=(tol2-0.5*(b-a))) {*xmin=x; return(fx);}
    if(fabs(e)>tol1) {
      r=(x-w)*(fx-fv); q=(x-v)*(fx-fw); p=(x-v)*q-(x-w)*r;
      q=2.0*(q-r); if(q>0.0) p=-p; q=fabs(q);
      etemp=e; e=d;
      if(fabs(p)>=fabs(0.5*q*etemp) || p<=q*(a-x) || p>=q*(b-x))
        d=CGOLD*(e=(x>=xm ? a-x:b-x));
      else {
        d=p/q; u=x+d;
        if(u-a<tol2 || b-u<tol2) d=copysign(tol1, xm-x);
      }
    } else {d=CGOLD*(e=(x>=xm ? a-x:b-x));}
    u=(fabs(d)>=tol1 ? x+d : x+copysign(tol1, d));
    fu=_powell_f1dim(u, dim);
    if(fu<=fx) {
      if(u>=x) a=x; else b=x;
      _powell_shft(&v, &w, &x, &u); _powell_shft(&fv, &fw, &fx, &fu);
    } else {
      if(u<x) a=u; else b=u;
      if(fu<=fw || w==x) {v=w; w=u; fv=fw; fw=fu;}
      else if(fu<=fv || v==x || v==w) {v=u; fv=fu;}
    }
  }
  *xmin=x;
  return(fx);
}
/******************************************************************************/
double _powell_f1dim(double x, int dim)
{
  int i;
  double f, xt[MAX_PARAMETERS];

  for(i=0; i<_powell_ncom; i++) xt[i]=_powell_pcom[i]+x*_powell_xicom[i];
  f=(*_powellFunc)(dim, xt, _powellFuncData); _powell_func_calls++;
  return(f);
}
/******************************************************************************/
void _powell_mnbrak(
  double *ax, double *bx, double *cx, double *fa, double *fb,
  double *fc, int dim
) {
  const double GOLD = 1.618034;
  const double GLIMIT = 100.0;
  const double TINY = 1.0e-20;
  double ulim, u, r, q, fu, dum=0.0;

  *fa=_powell_f1dim(*ax, dim); *fb=_powell_f1dim(*bx, dim);
  if(*fb>*fa) {_powell_shft(&dum, ax, bx, &dum); _powell_shft(&dum, fb, fa, &dum);}
  *cx=(*bx)+GOLD*(*bx-*ax); *fc=_powell_f1dim(*cx, dim);
  while((*fb)>(*fc)) {
    r=(*bx-*ax)*(*fb-*fc); q=(*bx-*cx)*(*fb-*fa);
    u=(*bx)-((*bx-*cx)*q-(*bx-*ax)*r)/(2.0*copysign(_powell_fmax(fabs(q-r),TINY),q-r));
    ulim=(*bx)+GLIMIT*(*cx-*bx);
    if(((*bx)-u)*(u-(*cx)) > 0.0) {
      fu=_powell_f1dim(u, dim);
      if(fu < *fc) {*ax=(*bx); *bx=u; *fa=(*fb); *fb=fu; return;}
      else if(fu > *fb) {*cx=u; *fc=fu; return;}
      u=(*cx)+GOLD*(*cx-*bx);
      fu=_powell_f1dim(u, dim);
    } else if((*cx-u)*(u-ulim) > 0.0) {
      fu=_powell_f1dim(u, dim);
      if(fu < *fc) {
        q=*cx+GOLD*(*cx-*bx); r=_powell_f1dim(u, dim);
        _powell_shft(bx, cx, &u, &q); _powell_shft(fb, fc, &fu, &r);
      }
    } else if((u-ulim)*(ulim-*cx) >= 0.0) {
      u=ulim; fu=_powell_f1dim(u, dim);
    } else {
      u=(*cx)+GOLD*(*cx-*bx);
      fu=_powell_f1dim(u, dim);
    }
    _powell_shft(ax, bx, cx, &u); _powell_shft(fa, fb, fc, &fu);
  }
}
/******************************************************************************/
/// @endcond

/******************************************************************************/
