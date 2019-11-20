/// @file constraints.c
/// @author Vesa Oikonen
/// @brief Setting and checking fit parameter constraints and limits.
///
/****************************************************************************/
#include "libtpcmodel.h"
/****************************************************************************/

/****************************************************************************/
/** Check that model parameters are within given limits.
 *  If not, then compute a penalty factor.
\return Return the number of parameters that are inside or at the limits;
    0 in case of error or if all are outside of the limits.
 */
int modelCheckParameters(
  /** Nr of parameters */
  int par_nr,
  /** Lower limits */
  double *lower_p,
  /** Upper limits */
  double *upper_p,
  /** Parameters to test */
  double *test_p,
  /** Pointer to corrected parameters (NULL if not needed) */
  double *accept_p,
  /** Pointer to variable in which the possible penalty factor will be written;
      1 if no penalty, or >1. Set to NULL if not needed. */
  double *penalty
) {
  int pi, accept_nr=0;
  double range;
  
  if(penalty!=NULL) *penalty=1.0;
  for(pi=0; pi<par_nr; pi++) {
    range=upper_p[pi]-lower_p[pi]; if(range<=0.0) range=1.0;
    if(test_p[pi]<lower_p[pi]) {
      if(accept_p!=NULL) accept_p[pi]=lower_p[pi];
      if(penalty!=NULL) *penalty += (lower_p[pi]-test_p[pi])/range;
    } else if(test_p[pi]>upper_p[pi]) {
      if(accept_p!=NULL) accept_p[pi]=upper_p[pi];
      if(penalty!=NULL) *penalty += (test_p[pi]-upper_p[pi])/range;
    } else {
      if(accept_p!=NULL) accept_p[pi]=test_p[pi];
      accept_nr++;
    }
  }
  return accept_nr;
}
/****************************************************************************/

/****************************************************************************/
/** Check if model parameters have collided with given limits.

    If parameter is fixed (equal lower and upper limit) then it is not counted
    as collision.
   @return Return the number of parameters that too close to the limits;
   0 in case all are well inside the limits.
 */
int modelCheckLimits(
  /** Nr of parameters */
  int par_nr,
  /** Lower limits */
  double *lower_p,
  /** Upper limits */
  double *upper_p,
  /** Parameters to test */
  double *test_p
) {
  int pi, collision_nr=0;
  double range, range_factor=0.00001;
  
  for(pi=0; pi<par_nr; pi++) {
    // fixed parameters are not checked
    range=upper_p[pi]-lower_p[pi]; if(range<=0.0) continue;
    // if parameter exceeds the limit, then that is a collision for sure
    if(test_p[pi]<=lower_p[pi]) {
      collision_nr++; continue;
    } else if(test_p[pi]>upper_p[pi]) {
      collision_nr++; continue;
    }
    // but in addition, accepted range is a bit smaller than constrained range
    // except if limit is 0, then it is ok that parameter is close to zero
    range*=range_factor;
    if(test_p[pi]<lower_p[pi]+range) {
      if(fabs(lower_p[pi])>=1.0E-06) collision_nr++;
    } else if(test_p[pi]>upper_p[pi]-range) {
      if(fabs(upper_p[pi])>=1.0E-06) collision_nr++;
    }
  }
  return collision_nr;
}
/****************************************************************************/

/****************************************************************************/
/** @brief Estimate initial values for sum of exponentials to be fitted on 
    decaying x,y-data.
  
    @return 0 when successful, otherwise <>0.
 */
int fitExpDecayNNLS(
  /** Pointer to array of x values; must not contain NaNs; data is not changed;
      samples must be sorted by increasing x. */
  double *x,
  /** Pointer to array of y values; must not contain NaNs; data is not changed. */
  double *y,
  /** Nr of x and y samples (x and y array lengths); data is not changed. */
  int n,
  /** Fittime is usually set to <=0 or a very high value to start the search
      for the best line fit from the last x,y sample towards the first sample;
      However, to exclude the end phase you may want to set fittime to include
      only certain time range from the beginning. */
  double fittime,
  /** Minimum eigenvalue, must be >0; for example, 1.0E-06. */
  double kmin,
  /** Maximum eigenvalue; for example, 1.0E+03. */
  double kmax,
  /** Max number of exp functions to save; length of a[] and k[] arrays. */
  int pnr,
  /** Pointer to an array of length pnr where exp function coefficients will be 
      written; enter NULL if not needed. */
  double *a,
  /** Pointer to an array of length pnr where exp function eigenvalues will be 
      written; enter NULL if not needed. */
  double *k,
  /** The number of fitted exponentials will be written in here; 
      note that this number may be higher than pnr; enter NULL if not needed. */
  int *fnr,
  /** Verbose level; if zero, then nothing is printed to stderr or stdout */
  int verbose
) {
  if(verbose>0) 
    printf("fitExpDecayNNLS(x, y, %d, %g, %g, %g, ...)\n", n, fittime, kmin, kmax);
  if(n<3) return(1);
  if(fittime>0.0) {
    while(x[n-1]>fittime && n>0) n--;
    if(verbose>1) printf("  n := %d\n", n);
  } 
  if(kmin<1.0E-100) kmin=1.0E-100;
  if(kmax<=kmin) return(1);

  /* Allocate memory for NNLS */
  int      NNLS_N=100;
  int      nn, mm, nnls_n, nnls_m, nnls_index[NNLS_N];
  double  *nnls_a[NNLS_N], *nnls_b, *nnls_zz, nnls_x[NNLS_N], *nnls_mat,
           nnls_wp[NNLS_N], *dptr, nnls_rnorm;

  nnls_n=NNLS_N;
  nnls_m=n;
  nnls_mat=(double*)malloc(((nnls_n+2)*nnls_m)*sizeof(double));
  if(nnls_mat==NULL) return(2);
  for(nn=0, dptr=nnls_mat; nn<nnls_n; nn++) {nnls_a[nn]=dptr; dptr+=nnls_m;}
  nnls_b=dptr; dptr+=nnls_m; nnls_zz=dptr;

  /* Set exponent function decay parameters */
  double epar[NNLS_N];
  {
    double elnmin, elnmax;
    elnmin=log(kmin); elnmax=log(kmax);
    //double elnmin=-12.5, elnmax=6.7;
    double d, r;
    r=elnmax-elnmin;
    d=r/(double)(NNLS_N-1);
    for(nn=0; nn<nnls_n; nn++) epar[nn]=-exp(elnmin+(double)nn*d);
  }

  /* Fill NNLS matrix */
  for(mm=0; mm<nnls_m; mm++) {
    nnls_b[mm]=y[mm];
  }
  for(nn=0; nn<nnls_n; nn++)
    for(mm=0; mm<nnls_m; mm++)
      nnls_a[nn][mm]=exp(epar[nn]*x[mm]);

  /* NNLS */
  int ret;
  ret=nnls(nnls_a, nnls_m, nnls_n, nnls_b, nnls_x, &nnls_rnorm,
           nnls_wp, nnls_zz, nnls_index);
  if(ret>1) {
    if(verbose>0) fprintf(stderr, "Error: NNLS solution not possible.\n");
    free(nnls_mat);
    return(3);
  } else if(ret==1) {
    if(verbose>0) fprintf(stderr, "Warning: max iteration count exceeded in NNLS.\n");
  }
  if(verbose>2) {
    printf("NNLS results:\n");
    for(nn=0; nn<nnls_n; nn++) 
      printf("\t%e\t%g\t%g\n", epar[nn], nnls_x[nn], nnls_wp[nn]);
  }
  if(verbose>1) {
    printf("Reasonable NNLS results:\n");
    for(nn=0; nn<nnls_n; nn++) 
      if(nnls_wp[nn]==0.0) printf("\t%e\t%g\n", epar[nn], nnls_x[nn]);
  }

  /* Reset exponent function decay parameters with the clusters found above */
  {
    if(verbose>1) printf("Cluster means:\n");
    int i, j, nr;
    double ev;
    i=j=0;    
    while(1) {
      /* jump over functions with wp<0 */
      while(i<nnls_n && nnls_wp[i]<0.0) i++;
      if(i==nnls_n) break;
      /* calculate mean of this cluster */
      nr=0; ev=0.0;
      while(i<nnls_n && nnls_wp[i]==0.0) {nr++; ev+=epar[i]; i++;}
      ev/=(double)nr;
      if(verbose>1) printf("mean_e := %e\n", ev);
      epar[j++]=ev;
    }
    nnls_n=j;
  }

  /* Fill NNLS matrix with these functions */
  for(mm=0; mm<nnls_m; mm++) {
    nnls_b[mm]=y[mm];
  }
  for(nn=0; nn<nnls_n; nn++)
    for(mm=0; mm<nnls_m; mm++)
      nnls_a[nn][mm]=exp(epar[nn]*x[mm]);

  /* NNLS */
  ret=nnls(nnls_a, nnls_m, nnls_n, nnls_b, nnls_x, &nnls_rnorm,
           nnls_wp, nnls_zz, nnls_index);
  if(ret>1) {
    if(verbose>0) fprintf(stderr, "Error: NNLS solution not possible.\n");
    free(nnls_mat);
    return(3);
  } else if(ret==1) {
    if(verbose>0) fprintf(stderr, "Warning: max iteration count exceeded in NNLS.\n");
  }
  if(verbose>1) {
    printf("NNLS results:\n");
    for(nn=0; nn<nnls_n; nn++) 
      printf("\t%e\t%g\t%g\n", epar[nn], nnls_x[nn], nnls_wp[nn]);
  }
  
  /* Return the results */
  int nr=0;
  for(nn=0; nn<nnls_n; nn++) {
    if(nnls_wp[nn]<0.0) continue;
    if(nr>=pnr) {nr++; continue;}
    if(a!=NULL) a[nr]=nnls_x[nn];
    if(k!=NULL) k[nr]=epar[nn];
    nr++;
  }
  if(fnr!=NULL) *fnr=nr;

  free(nnls_mat);
  return(0);
}
/*****************************************************************************/

/*****************************************************************************/
