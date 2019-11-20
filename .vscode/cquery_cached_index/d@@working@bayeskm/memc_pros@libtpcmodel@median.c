/// @file median.c
/// @author Vesa Oikonen
/// @brief Calculation of median value.
///
/******************************************************************************/
#include "libtpcmodel.h"
/******************************************************************************/

/******************************************************************************/
/** Returns the kth smallest value in data[0..n-1]. Array is partially sorted.
    Algorithm is based on the book Wirth N. Algorithms + data structures = programs. 
    Englewood Cliffs, Prentice-Hall, 1976.
    @return Returns the kth smallest value in data[0..n-1].
*/
double d_kth_smallest(
  /** Pointer to data; array is partially sorted. */
  double *data,
  /** Length of data array. */
  int n,
  /** kth smallest value will be returned. */
  int k
) {
  int i, j, l, m;
  double x, s;

  l=0; m=n-1;
  while(l<m) {
    x=data[k]; i=l; j=m;
    do {
      while(data[i]<x) i++;
      while(x<data[j]) j--;
      if(i<=j) {s=data[i]; data[i]=data[j]; data[j]=s; i++; j--;}
    } while(i<=j);
    if(j<k) l=i;
    if(k<i) m=j;
  }
  return(data[k]);
}
/******************************************************************************/

/******************************************************************************/
/** Returns the median in array data[0..n-1]. Array is partially sorted.
    Algorithm is based on the book Wirth N. Algorithms + data structures = programs. 
    Englewood Cliffs, Prentice-Hall, 1976.
    @sa dmean, mean, fmedian
    @return Returns the median in array data[0..n-1].
*/
double dmedian(
  /** Pointer to data; array is partially sorted. */
  double *data,
  /** Length of data array. */
  int n
) {
  int k;
  double d1, d2;

  if(n<1) return(0.0);
  if(n%2) {
    k=(n-1)/2; return(d_kth_smallest(data, n, k));
  } else {
    k=n/2; d1=d_kth_smallest(data, n, k-1); d2=d_kth_smallest(data, n, k);
    return(0.5*(d1+d2));
  }
}
/******************************************************************************/

/******************************************************************************/
/** Returns the mean in array data[0..n-1], and optionally calculates also
    the (sample) standard deviation of the mean.
    @sa dmedian, mean, dmean_nan, fmean
    @return Returns the mean in array data[0..n-1].
*/
double dmean(
  /** Pointer to data; data is not changed in any way. */
  double *data,
  /** Length of data array. */
  int n,
  /** Pointer to variable where SD will be written; enter NULL if not needed. */
  double *sd
) {
  int i;
  double sumsqr=0.0, sqrsum=0.0, avg;

  if(n<1 || data==NULL) {if(sd!=NULL) *sd=0.0; return(0.0);}

  for(i=0; i<n; i++) {sumsqr+=data[i]*data[i]; sqrsum+=data[i];}
  avg=sqrsum/(double)n; if(sd==NULL) return(avg);
  if(n==1) {
    *sd=0.0;
  } else {
    sqrsum*=sqrsum;
    *sd=sqrt( (sumsqr - sqrsum/(double)n) / (double)(n-1) );
  }
  return(avg);
}
/******************************************************************************/

/******************************************************************************/
/** Returns the mean in array data[0..n-1], and optionally calculates also
    the standard deviation of the mean. Data may contain missing samples marked as NaNs.
    @sa dmedian, mean, dmean
    @return Returns the mean in array data[0..n-1].
*/
double dmean_nan(
  /** Pointer to data; data is not changed in any way. */
  double *data,
  /** Length of data array. */
  int n,
  /** Pointer to variable where SD will be written; enter NULL if not needed. */
  double *sd,
  /** Pointer to variable where number of valid (not NaN) samples will be written; 
      enter NULL if not needed. */
  int *vn
) {
  int i, m;
  double sumsqr=0.0, sqrsum=0.0, avg;

  if(n<1 || data==NULL) {
    if(sd!=NULL) *sd=0.0; 
    if(vn!=NULL) *vn=0; 
    return(0.0);
  }
  for(i=m=0; i<n; i++) if(!isnan(data[i])) m++;
  if(vn!=NULL) *vn=m;
  if(m<1) {
    if(sd!=NULL) *sd=nan(""); 
    return(nan(""));
  }
  

  for(i=0; i<n; i++) if(!isnan(data[i])) {
    sumsqr+=data[i]*data[i]; sqrsum+=data[i];
  }
  avg=sqrsum/(double)m; if(sd==NULL) return(avg); // SD not requested
  if(m==1) {
    *sd=0.0;
  } else {
    sqrsum*=sqrsum;
    *sd=sqrt( (sumsqr - sqrsum/(double)m) / (double)(m-1) );
  }
  return(avg);
}
/******************************************************************************/

/******************************************************************************/
