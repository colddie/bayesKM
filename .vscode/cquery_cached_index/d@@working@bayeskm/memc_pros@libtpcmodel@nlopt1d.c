/// @file nlopt1d.c
/// @author Vesa Oikonen
/// @brief Nonlinear one-dimensional optimization.
///
/******************************************************************************/
#include "libtpcmodel.h"
/******************************************************************************/

/******************************************************************************/
/** Local one-dimensional minimization by bracketing.
    @return 0 in case of no errors.
    @sa simplex, tgo, bobyqa, powell, pearson
 */
int nlopt1D(
  /** Pointer to the function to be minimized. It must be defined in main program as: 
      double func(double p, void *data); where p is the parameter, and data points to
      user-defined data structure needed by the function. */
  double (*_fun)(double, void*),
  /** Pointer to data that will be passed on to the _fun(). */
  void *_fundata,
  /** Parameter initial value. */
  double x,
  /** Parameter lower limit. */
  double xl,
  /** Parameter upper limit. */
  double xu,
  /** Initial parameter delta. */
  double delta,
  /** Required tolerance. */
  double tol,
  /** Maximum number of function evaluations. */
  const int maxeval,
  /** Pointer for optimized parameter value. */
  double *nx,
  /** Pointer for function value at optimum; NULL if not needed. */
  double *nf,
  /** Verbose level */
  int verbose
) {
  if(verbose>0) {
    printf("%s(f, fdata, %g, %g, %g, %g, %g, %d, ...)\n", __func__, x, xl, xu, delta, tol, maxeval); 
    fflush(stdout);
  }


  /* Check the input */
  if(_fun==NULL) return(1);
  if(xl>xu || x<xl || x>xu) return(1);
  if(!(delta>0.0) || !(tol>0.0)) return(1);
  if(tol>=delta) return(1);
  if(maxeval<5) return(1);
  if(nx==NULL) return(1);

  double begin=xl;
  double end=xu;
  int nevals=0;

  double p1=0, p2=0, p3=0, f1=0, f2=0, f3=0;
  /* Calculate function value with initial parameter value */
  double minf=f2=_fun(x, _fundata); nevals++;
  /* If parameter is fixed, then this was all that we can do */
  if(xl>=xu) {
    *nx=x; if(nf!=NULL) *nf=minf;
    return(0);
  }

  /* Find three bracketing points such that f1 > f2 < f3.
     Do this by generating a sequence of points expanding away from 0.
     Also note that, in the following code, it is always the
     case that p1 < p2 < p3. */
  p2=x;

  /* Start by setting a starting set of 3 points that are inside the bounds */
  p1=p2-delta; if(p1>begin) p1=begin;
  p3=p2+delta; if(p3<end) p3=end;
  /* Compute their function values */
  f1=_fun(p1, _fundata); nevals++;
  f3=_fun(p3, _fundata); nevals++;
  if(p2==p1 || p2==p3) {
    p2=0.5*(p1+p3);
    f2=minf=_fun(p2, _fundata); nevals++;
  }

  /* Now we have 3 points on the function. 
     Start looking for a bracketing set such that f1 > f2 < f3 is the case. */
  double jump_size=delta;
  while(!(f1>f2 && f2<f3)) {
    /* check for hitting max_iter */
    if(verbose>5) printf("  bracketing: nevals=%d\n", nevals);
    if(nevals >= maxeval) {
      *nx=p2; minf=f2; if(nf!=NULL) *nf=minf; 
      return(0);
    }
    /* check if required tolerance was reached */
    if((p3-p1)<tol) { //if (p3-p1 < eps)
      if(verbose>1) printf("  max tolerance was reached during bracketing\n");
      if(f1<f2 && f1<f3) {
        *nx=p1; minf=f1; if(nf!=NULL) *nf=minf;
        return(0);
      }
      if(f2<f1 && f2<f3) {
        *nx=p2; minf=f2; if(nf!=NULL) *nf=minf;
        return(0);
      }
      *nx=p3; minf=f3; if(nf!=NULL) *nf=minf;
      return(0);
    }
    if(verbose>6) printf("    jump_size=%g\n", jump_size);
    /* if f1 is small then take a step to the left */
    if(f1<f3) { 
      /* check if the minimum is colliding against the bounds. If so then pick
         a point between p1 and p2 in the hopes that shrinking the interval will
         be a good thing to do.  Or if p1 and p2 aren't differentiated then try
         and get them to obtain different values. */
      if(p1==begin || (f1==f2 && (end-begin)<jump_size )) {
        p3=p2; f3=f2; p2=0.5*(p1+p2);
        f2=minf=_fun(p2, _fundata); nevals++;
      } else {
        /* pick a new point to the left of our current bracket */
        p3=p2; f3=f2; p2=p1; f2=f1;
        p1-=jump_size; if(p1<begin) p1=begin;
        f1=_fun(p1, _fundata); nevals++;
        jump_size*=2.0;
      }
    } else { // otherwise f3 is small and we should take a step to the right
      /* check if the minimum is colliding against the bounds. If so then pick
         a point between p2 and p3 in the hopes that shrinking the interval will
         be a good thing to do.  Or if p2 and p3 aren't differentiated then
         try and get them to obtain different values. */
      if(p3==end || (f2==f3 && (end-begin)<jump_size)) {
        p1=p2; f1=f2; p2=0.5*(p3+p2);
        f2=minf=_fun(p2, _fundata); nevals++;
      } else {
        /* pick a new point to the right of our current bracket */
        p1=p2; f1=f2; p2=p3; f2=f3;
        p3+=jump_size; if(p3>end) p3=end;
        f3=minf=_fun(p3, _fundata); nevals++;
        jump_size*=2.0;
      }
    }
  }
  if(verbose>4) printf("  brackets ready\n");

  /* Loop until we have done the max allowable number of iterations or
     the bracketing window is smaller than eps.
     Within this loop we maintain the invariant that: f1 > f2 < f3 and 
     p1 < p2 < p3. */
  double d, d2;
  const double tau=0.1;
  double p_min, f_min;
  while((nevals<maxeval) && (p3-p1>tol)) {
    if(verbose>5) printf("  main loop: nevals=%d\n", nevals);

    //p_min = lagrange_poly_min_extrap(p1,p2,p3, f1,f2,f3);
    d=f1*(p3*p3-p2*p2) + f2*(p1*p1-p3*p3) + f3*(p2*p2-p1*p1);
    d2=2.0*(f1*(p3-p2) + f2*(p1-p3) + f3*(p2-p1));
    if(d2==0.0 || !isfinite(d/=d2)) { // d=d/d2
      p_min=p2;
    } else {
      if(p1<=d && d<=p3) {p_min=d;}
      else {p_min=d; if(p1>p_min) p_min=p1; if(p3<p_min) p_min=p3;}
    }

    /* make sure p_min isn't too close to the three points we already have */
    if(p_min<p2) {
      d=(p2-p1)*tau;
      if(fabs(p1-p_min)<d) p_min=p1+d;
      else if(fabs(p2-p_min)<d) p_min=p2-d;
    } else {
      d=(p3-p2)*tau;
      if(fabs(p2-p_min)<d) p_min=p2+d;
      else if(fabs(p3-p_min)<d) p_min=p3-d;
    }

    /* make sure one side of the bracket isn't super huge compared to the other
       side.  If it is then contract it. */
    double bracket_ratio=fabs(p1-p2)/fabs(p2-p3);
    if(!(bracket_ratio<100.0 && bracket_ratio>0.01)) {
      /* Force p_min to be on a reasonable side. */
      if(bracket_ratio>1.0 && p_min>p2) p_min=0.5*(p1+p2);
      else if(p_min<p2) p_min=0.5*(p2+p3);
    }

    /* Compute function value at p_min */
    *nx=p_min;
    f_min=minf=_fun(p_min, _fundata); nevals++;

    /* Remove one of the endpoints of our bracket depending on where the new point falls */
    if(p_min<p2) {
      if(f1>f_min && f_min<f2) {p3=p2; f3=f2; p2=p_min; f2=f_min;}
      else {p1=p_min; f1=f_min;}
    } else {
      if(f2>f_min && f_min<f3) {p1=p2; f1=f2; p2=p_min; f2=f_min;}
      else {p3=p_min; f3=f_min;}
    }
  }
  *nx=p2; minf=f2; if(nf!=NULL) *nf=minf;
  return(0);
}
/******************************************************************************/

/******************************************************************************/
