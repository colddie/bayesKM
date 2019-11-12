/// @file pearson.c
/// @author Vesa Oikonen
/// @brief Pearson's correlation coefficient and regression line.
///
#include "libtpcmodel.h"
/******************************************************************************/

/******************************************************************************/
/** Calculate slope and intercept of a line and Pearson's correlation coefficient.
    @sa regr_line, best_pearson, highest_slope, pearson2, pearson3, pearson4, mean,
        imgsegmPearson, nlopt1D
    @return If successful, function returns value 0.
*/

int PEARSON_TEST;

int pearson(
  /** data x values */
  double *x,
  /** data y values */
  double *y,
  /** number of data sample values */
  int nr,
  /** slope */
  double *k,
  /** S.D. of slope */
  double *kSD,
  /** y axis intercept */
  double *b,
  /** S.D. of y axis intercept */
  double *bSD,
  /** Pearson's correlation coefficient, r */
  double *r,
  /** Residual variance of y values */
  double *ySD
) {
  int i;
  double e, f, g, meanx, meany, kk, bb, ke, be, rr, sumsdcy=0.0;
  double sumx=0.0, sumy=0.0, sumsx=0.0, sumsdx=0.0, sumsdy=0.0, sumdxdy=0.0;
  double sumxy=0.0, sumsy=0.0;


  if(PEARSON_TEST) {
    fprintf(stdout, "pearson(x[], y[], %d, *k, *kSD, *b, *bSD, *r, *ySD)\n", nr);
    fflush(stdout);
  }
  /* Check that there is some data */
  if(x==NULL || y==NULL || nr<2) return(1);

  /* If only 2 datapoints */
  if(nr==2) {
    f=x[1]-x[0]; if(fabs(f)<1.0E-50) return(1);
    *k=(y[1]-y[0])/f; *b=y[0]-(*k)*x[0];
    *kSD=*bSD=*ySD=0.; *r=1.;
    return(0);
  }

  /* Calculate (x,y) sums and means */
  for(i=0; i<nr; i++) {
    sumx+=x[i]; sumy+=y[i];
    sumsx+=x[i]*x[i]; sumsy+=y[i]*y[i];
    sumxy+=x[i]*y[i];
  }
  meanx=sumx/(double)nr;
  meany=sumy/(double)nr;
  /* and then based on means */
  for(i=0; i<nr; i++) {
    f=x[i]-meanx; sumsdx+=f*f;
    g=y[i]-meany; sumsdy+=g*g;
    sumdxdy+=f*g;
  }
  if(sumsdx<1.0e-50 || sumsdy<1.0e-50) return(3);
  /* Regression coefficient */
  kk=sumdxdy/sumsdx; *k=kk;
  /* Intercept with y axis */
  bb=(sumsdx*sumy - sumx*sumdxdy)/((double)nr*sumsdx); *b=bb;
  /* Errors */
  for(i=0; i<nr; i++) {
    f=kk*x[i]+bb-y[i];
    sumsdcy+=f*f;
  }
  /* Deviation of y values */
  if(sumsdcy<=1.0e-12) e=0.0; else e=sqrt(sumsdcy/(double)(nr-2)); *ySD=e;
  /* SD of slope and intercept */
  ke=e/sqrt(sumsdx); be=e/sqrt((double)nr-sumx*sumx/sumsx);
  *kSD=ke; *bSD=be;
  be=sqrt(sumsdcy/(double)(nr-2))/sqrt((double)nr-(sumx*sumx)/sumsx);
  /* Pearson's correlation coefficient */
  rr=(sumxy-((sumx*sumy)/(double)nr)) /
     sqrt((sumsx-sumx*sumx/(double)nr)*(sumsy-sumy*sumy/(double)nr));
  /* Correct for small sample size */
  if(nr>4) rr*=1.0+(1.0-rr*rr)/(double)(2*(nr-4));
  *r=rr;
  if(PEARSON_TEST) {
    fprintf(stdout, "k=%14.7e +- %14.7e\n", kk, ke);
    fprintf(stdout, "b=%14.7e +- %14.7e\n", bb, be);
    fprintf(stdout, "r=%14.7e ySD=%14.7e\n", rr, e);
    fflush(stdout);
  }

  return(0);
}
/******************************************************************************/

/******************************************************************************/
/** Calculate slope and intercept of a line and Pearson's correlation coefficient.

    Array char is[] specifies whether single (x,y) points are used in the fit.
    @sa regr_line, best_pearson, pearson, highest_slope
    @return If successful, function returns value 0.
*/
int pearson2(
  /** data x values */
  double *x,
  /** data y values */
  double *y,
  /** Switch values: 0=do not use this point */
  char *is,
  /** number of data sample values */
  int nr,
  /** slope */
  double *k,
  /** S.D. of slope */
  double *kSD,
  /** y axis intercept */
  double *b,
  /** S.D. of y axis intercept */
  double *bSD,
  /** Pearson's correlation coefficient, r */
  double *r,
  /** Residual variance of y values */
  double *ySD
) {
  int i, j;
  double *nx, *ny;

  /* Allocate memory for new data pointers  */
  nx=(double*)calloc(nr, sizeof(double)); if(nx==NULL) return 1;
  ny=(double*)calloc(nr, sizeof(double)); if(ny==NULL) {free((char*)nx); return 1;}

  /* Copy data to pointers */
  for(i=0, j=0; i<nr; i++) if(is[i]) {nx[j]=x[i]; ny[j]=y[i]; j++;}

  /* Use pearson() */
  i=pearson(nx, ny, j, k, kSD, b, bSD, r, ySD); free((char*)nx); free((char*)ny);
  return (i);
}
/******************************************************************************/

/******************************************************************************/
/** Calculate slope and intercept of a line and Pearson's correlation coefficient.

    Data points may contain NaN's.
    @sa regr_line, pearson, best_pearson, highest_slope
    @return If successful, function returns value 0.
*/
int pearson3(
  /** data x values */
  double *x,
  /** data y values */
  double *y,
  /** number of data sample values */
  int nr,
  /** slope */
  double *k,
  /** S.D. of slope */
  double *kSD,
  /** y axis intercept */
  double *b,
  /** S.D. of y axis intercept */
  double *bSD,
  /** Pearson's correlation coefficient, r */
  double *r,
  /** Residual variance of y values */
  double *ySD
) {
  int i, j;
  double *nx, *ny;

  /* Allocate memory for new data pointers  */
  nx=(double*)calloc(nr, sizeof(double)); if(nx==NULL) return 1;
  ny=(double*)calloc(nr, sizeof(double)); if(ny==NULL) {free((char*)nx); return 1;}

  /* Copy data to pointers  */
  for(i=0, j=0; i<nr; i++)
    if(!isnan(x[i]) && !isnan(y[i])) {nx[j]=x[i]; ny[j]=y[i]; j++;}

  /* Use pearson()  */
  i=pearson(nx, ny, j, k, kSD, b, bSD, r, ySD); free((char*)nx); free((char*)ny);
  return (i);
}
/******************************************************************************/

/******************************************************************************/
/** @brief Calculate slope and intercept of a line and Pearson's correlation coefficient.
 
    Data points may contain NaN's. Fit start and end times are specified.
    @sa regr_line, best_pearson, highest_slope
    @return If successful, function returns value 0.
 */
int pearson4(
  /** data x values */
  double *x,
  /** data y values */
  double *y,
  /** number of data sample values */
  int nr,
  /** fit start time */
  double start,
  /** fit end time */
  double end,
  /** slope */
  double *k,
  /** S.D. of slope */
  double *kSD,
  /** y axis intercept */
  double *b,
  /** S.D. of y axis intercept */
  double *bSD,
  /** Pearson's correlation coefficient, r */
  double *r,
  /** Residual variance of y values */
  double *ySD
) {
  int i, j;
  double *nx, *ny;

  if(PEARSON_TEST) {
    fprintf(stdout, "pearson4(x[], y[], %d, %g, %g, *k, *kSD, *b, *bSD, *r, *ySD)\n",
                    nr, start, end);
    fflush(stdout);
  }
  /* Allocate memory for new data pointers */
  if((nx=(double*)calloc(nr, sizeof(double)))==NULL) return -3;
  if((ny=(double*)calloc(nr, sizeof(double)))==NULL) {free((char*)nx); return -3;}

  /* Copy data to pointers */
  for(i=0, j=0; i<nr; i++)
    if(x[i]>=start && x[i]<=end && !isnan(x[i]) && !isnan(y[i])) {
      nx[j]=x[i]; ny[j]=y[i]; j++;}

  /* Use pearson() */
  i=pearson(nx, ny, j, k, kSD, b, bSD, r, ySD);
  free((char*)nx); free((char*)ny);
  return i;
}
/******************************************************************************/

/******************************************************************************/
/** Find the best linear fit to double data (x[], y[]) with nr points.

    Data may contain NaN's, which are not used.
    @sa highest_slope, regr_line, pearson
    @return Returns the nr of points actually used, or 0, in case of error.
 */
int best_pearson(
  /** Data x values */
  double *x,
  /** Data y values */
  double *y,
  /** Number of data sample values */
  int nr,
  /** Minimum nr of data points to use */
  int min_nr,
  /** Index [0..last-2] of the first point to use initially, and after fitting */
  int *first,
  /** Index [first+1..nr-1] of the last point to use initially, and after fitting */
  int *last,
  /** slope */
  double *k,
  /** S.D. of slope */
  double *kSD,
  /** y axis intercept */
  double *b,
  /** S.D. of y axis intercept */
  double *bSD,
  /** Pearson's correlation coefficient, r */
  double *r,
  /** Residual variance of y values */
  double *ySD
) {
  int i, n, m, x1, x2, b1, b2, bm;
  double *nx, *ny;
  double tk, tkSD, tb, tbSD, tr, tySD;


  if(PEARSON_TEST) {
    fprintf(stdout, "best_pearson(x, y, %d, %d, %d, %d, k, kSD, b, bSD, r, ySD)\n",
      nr, min_nr, *first, *last);
    fflush(stdout);
  }
  /* Remove NaN's and those outside range first-last */
  /* Allocate memory for new data pointers          */
  nx=(double*)calloc(nr, sizeof(double)); if(nx==NULL) return 0;
  ny=(double*)calloc(nr, sizeof(double)); if(ny==NULL) {free((char*)nx); return 0;}
  /* Copy data to pointers */
  if(*last>nr-1) *last=nr-1;
  if(*first<0) *first=0;
  for(i=*first, n=0; i<=*last; i++)
    if(!isnan(x[i]) && !isnan(y[i])) {nx[n]=x[i]; ny[n]=y[i]; n++;}

  /* Check that we have enough points */
  if(n<2 || n<min_nr) {free((char*)nx); free((char*)ny); return 0;}
  if(n==min_nr) {
    i=pearson(nx, ny, n, k, kSD, b, bSD, r, ySD);
    free((char*)nx); free((char*)ny);
    if(i) return 0; else return n;
  }

  /* OK, let's start  */
  *k=*kSD=*b=*bSD=*r=*ySD=0.; b1=b2=bm=0;
  for(x1=0, m=n; n-x1>min_nr; x1++) {
    /*printf("I x1=%i x2=%i m=%i n=%i\n", x1, x2, m, n);*/
    for(x2=n-1, m=x2-x1+1; m>=min_nr; x2--, m--) {
      if(pearson(nx+x1, ny+x1, m, &tk, &tkSD, &tb, &tbSD, &tr, &tySD))
        continue;
      if((tr>*r) ||
          (tr==*r && m>bm) ||
          (tr==*r && m==bm && x1>b1) ||
          (tr==*r && m==bm && x1==b1 && tk>*k) ) {
        bm=m; b1=x1; b2=x2;
        *k=tk; *kSD=tkSD; *b=tb; *bSD=tbSD; *r=tr; *ySD=tySD;
      }
      /*printf("Fit range: %2i-%2i  ;  best fit: %2i-%2i\n", x1, x2, b1, b2);*/
    }
  }

  /* Set first&last  */
  for(i=*first; i<=*last; i++)
    if(x[i]==nx[b1] && y[i]==ny[b1]) {*first=i; break;}
  for(i=*last; i>=*first; i--)
    if(x[i]==nx[b2] && y[i]==ny[b2]) {*last=i; break;}
  /*printf("FIRST=%i  LAST=%i\n", *first, *last);*/

  free((char*)nx); free((char*)ny);
  return bm;
}
/******************************************************************************/

/******************************************************************************/
/** Calculates the mean and SD of data. Data (y data) may contain NaN's.
    @sa dmean, dmedian, dmean_nan, regr_line, pearson
    @return Returns !=0 in case of an error.
 */
int mean(
  /** Data x values */
  double *x,
  /** Data y values */
  double *y,
  /** Number of data sample values */
  int nr,
  /** Calculated x mean */
  double *xmean,
  /** Calculated SD of x mean */
  double *xsd,
  /** Calculated y mean */
  double *ymean,
  /** Calculated SD of y mean */
  double *ysd
) {
  int   i, n;
  double sumsqr, sqrsum;

  for(i=0, sumsqr=sqrsum=0.0, n=0; i<nr; i++) if(!isnan(x[i]) && !isnan(y[i])) {
    sumsqr+=x[i]*x[i]; sqrsum+=x[i]; n++;}
  if(n<=0) return -1;
  *xmean=sqrsum/(double)n;
  if(n==1) *xsd=0.0; else {
    sqrsum*=sqrsum;
    *xsd=sqrt( (sumsqr - sqrsum/(double)n) / (double)(n-1) );
  }
  for(i=0, sumsqr=sqrsum=0.0, n=0; i<nr; i++) if(!isnan(x[i]) && !isnan(y[i])) {
    sumsqr+=y[i]*y[i]; sqrsum+=y[i]; n++;}
  if(n<=0) return -1;
  *ymean=sqrsum/(double)n;
  if(n==1) *ysd=0.0; else {
    sqrsum*=sqrsum;
    *ysd=sqrt( (sumsqr - sqrsum/(double)n) / (double)(n-1) );
  }
  return 0;
}
/******************************************************************************/

/******************************************************************************/
/** Calculates regression line slope (m) and y axis intercept.

    Data (x and y data) may contain NaN's.
    @sa highest_slope, pearson, imgsegmPearson
    @return Returns 0 if ok.
 */
int regr_line(
  /** An array of x axis values. */
  double *x,
  /** An array of y axis values. */
  double *y,
  /** The number of values in x and y arrays. */
  int n,
  /** Pointer where calculated slope is written. */
  double *m,
  /** Pointer where calculated y axis intercept is written. */
  double *c
) {
  double xsum=0.0, ysum=0.0, x2sum=0.0, xysum=0.0, delta;
  int i, nn=0;

  /* Check the data */
  if(x==NULL || y==NULL) return(1);
  if(n<2) return(2);

  /* Compute */
  for(i=0, nn=0; i<n; i++) {
    if(isnan(x[i]) || isnan(y[i])) continue;
    xsum+=x[i]; ysum+=y[i]; x2sum+=x[i]*x[i]; xysum+=x[i]*y[i]; nn++;
  }
  if(nn<2) return(2);
  delta=(double)nn*x2sum - xsum*xsum; if(delta==0.0) return(3);
  *m=((double)nn*xysum - xsum*ysum)/delta;
  *c=(x2sum*ysum - xsum*xysum)/delta;
  return(0);
}
/******************************************************************************/

/******************************************************************************/
/** Finds the regression line with the highest slope for x,y data.
    @sa highest_slope_after, regr_line, best_pearson
    @return Return 0 if ok.
 */
int highest_slope(
  /** An array of x axis values */
  double *x,
  /** An array of y axis values */
  double *y,
  /** The number of values in x and y arrays */
  int n,
  /** The number of samples used to fit the line */
  int slope_n,
  /** Pointer where calculated slope is written; NULL if not needed */
  double *m,
  /** Pointer where calculated y axis intercept is written; NULL if not needed */
  double *c,
  /** Pointer where calculated x axis intercept is written; NULL if not needed */
  double *xi,
  /** Pointer where the place (x) of the highest slope is written;
      NULL if not needed */
  double *xh
) {
  int i, ret, i_at_max=0;
  double slope, ic, max_slope, ic_at_max=0.0;

  /* Check the data */
  if(x==NULL || y==NULL) return(1);
  if(n<2) return(2);
  if(slope_n>n) return(3);

  /* Compute */
  max_slope=-1.0E200;
  for(i=0; i<=n-slope_n; i++) {
    ret=regr_line(x+i, y+i, slope_n, &slope, &ic); if(ret) continue;
    if(slope>max_slope) {max_slope=slope; ic_at_max=ic; i_at_max=i;}
  } /* next slope */
  if(max_slope==-1.0E200) return(10);
  if(m!=NULL) *m=max_slope;
  if(c!=NULL) *c=ic_at_max;
  if(xi!=NULL) {if(max_slope!=0.0) *xi=-ic_at_max/max_slope; else *xi=0.0;}
  if(xh!=NULL) {
    *xh=0.0; for(i=i_at_max; i<i_at_max+slope_n; i++) *xh+=x[i];
    *xh/=(double)slope_n;
  }

  return(0);
}
/******************************************************************************/

/******************************************************************************/
/** Finds the regression line with the highest slope for x,y data after specified x.
    @sa highest_slope, regr_line
    @return Return 0 if ok.
 */
int highest_slope_after(
  /** An array of x axis values. */
  double *x,
  /** An array of y axis values. */
  double *y,
  /** The number of values in x and y arrays. */
  int n,
  /** The number of samples used to fit the line. */
  int slope_n,
  /** Estimation start x value, samples with smaller x are ignored;
      can usually be set to zero. */
  double x_start,
  /** Pointer where calculated slope is written; NULL if not needed. */
  double *m,
  /** Pointer where calculated y axis intercept is written; NULL if not needed. */
  double *c,
  /** Pointer where calculated x axis intercept is written; NULL if not needed. */
  double *xi,
  /** Pointer where the place (x) of the highest slope is written; NULL if not needed. */
  double *xh
) {
  int i, ret, i_at_max=0;
  double slope, ic, max_slope, ic_at_max=0.0;

  /* Check the data */
  if(x==NULL || y==NULL) return(1);
  if(n<2) return(2);
  if(slope_n>n) return(3);

  /* Compute */
  max_slope=-1.0E200;
  for(i=0; i<=n-slope_n; i++) if(x[i]>=x_start) {
    ret=regr_line(x+i, y+i, slope_n, &slope, &ic); if(ret) continue;
    if(slope>max_slope) {max_slope=slope; ic_at_max=ic; i_at_max=i;}
  } /* next slope */
  if(max_slope==-1.0E200) return(10);
  if(m!=NULL) *m=max_slope;
  if(c!=NULL) *c=ic_at_max;
  if(xi!=NULL) {if(max_slope!=0.0) *xi=-ic_at_max/max_slope; else *xi=0.0;}
  if(xh!=NULL) {
    *xh=0.0; for(i=i_at_max; i<i_at_max+slope_n; i++) *xh+=x[i];
    *xh/=(double)slope_n;
  }

  return(0);
}
/******************************************************************************/

/******************************************************************************/
