/** @file doubleutil.c
 *  @brief Working with doubles.
 *  @author Vesa Oikonen
 */
/*****************************************************************************/

/*****************************************************************************/
#include "tpcclibConfig.h"
/*****************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <ctype.h>
#include <unistd.h>
#include <math.h>
#include "libtpcmisc.h"
/*****************************************************************************/

/*****************************************************************************/
/** Verifies that given two doubles have the same value inside given limits.
 *  Values are considered to match also if both are NaNs.
 *  @return 1 if match is found, and 0 if not.
 *  @author Vesa Oikonen
 *  @sa doubleMatchRel, doubleMachEps
 */
int doubleMatch(
  /** First value */
  const double v1,
  /** Second value */
  const double v2,
  /** Limit for absolute difference (if <0 then test will fail every time) */
  const double lim
) {
  if(isnan(v1) && isnan(v2)) return 1;
  if(isnan(v1) || isnan(v2)) return 0;
  if(v1==v2) return 1;
  if(isnan(lim) || lim<0.0) return 0;
  if(fabs(v1-v2)<=lim) return 1;
  return 0;
}
/*****************************************************************************/

/*****************************************************************************/
/** Verifies that given two doubles have the same value inside given relative
 *  limit |2*(v1-v2)/(v1+v2)|.
 *  Values are considered to match also if both are NaNs, but not if mean
 *  is zero (unless both are exactly zero).
 *  @return 1 if match is found, and 0 if not.
 *  @author Vesa Oikonen
 *  @sa doubleMatch, doubleMachEps
 */
int doubleMatchRel(
  /** First value */
  const double v1,
  /** Second value */
  const double v2,
  /** Limit for relative difference (if <0 then test will fail every time) */
  const double lim
) {
  if(isnan(v1) && isnan(v2)) return 1;
  if(isnan(v1) || isnan(v2)) return 0;
  if(v1==v2) return 1;
  if(isnan(lim)) return 0;
  double mean;
  mean=0.5*(v1+v2); if(!isnormal(mean)) return 0;
  if(fabs((v1-v2)/mean)<=lim) return 1;
  return 0;
}
/*****************************************************************************/

/*****************************************************************************/
/** Estimates the machine epsilon, the upper bound on the relative error due 
 *  to rounding in floating point arithmetic, within one order of magnitude
 *  of the true machine epsilon.
 *  Standard C library should also have DBL_EPSILON in float.h.
 *  @return Estimate of machine epsilon.
 *  @author Vesa Oikonen
 *  @sa doubleMatchRel, doubleMatch
 */
double doubleMachEps()
{
  double macheps=1.0;
  do {macheps/=2.0;} while((1.0+macheps/2.0)!=1.0);
  return(macheps);
}
/*****************************************************************************/

/*****************************************************************************/
/** Copy double values from the 2nd array to the first.
 */
void doubleCopy(
  /** Target array */
  double *t,
  /** Source array */
  double *s,
  /** Length of arrays */
  const unsigned int n
) {
  unsigned int i;
  if(t==NULL || s==NULL || n<1) return;
  for(i=0; i<n; i++) t[i]=s[i];
}
/*****************************************************************************/

/*****************************************************************************/
/** Find the maximum value in given double array.
    @return The index [0..n-1] of maximum value in array; 0 is returned also
     in case of errors.
 */
unsigned int doubleMaxIndex(
  /** Pointer to double array. */
  double *a,
  /** Length of array. */
  const unsigned int n
) {
  if(a==NULL || n<1) return(0);
  unsigned int i, mi=0;
  double mv=nan("");
  for(i=0; i<n; i++) if(isnan(mv) || a[i]>mv) {mv=a[i]; mi=i;}
  return(mi);
}
/*****************************************************************************/

/*****************************************************************************/
/** Calculate the sum of values in given double array.
    @return The sum of array values.
    @sa doubleMean
 */
double doubleSum(
  /** Pointer to double array. */
  double *a,
  /** Length of array. */
  const unsigned int n
) {
  double s=0.0;
  if(a==NULL || n<1) return(s);
  for(unsigned int i=0; i<n; i++) if(!isnan(a[i])) s+=a[i];
  return(s);
}
/*****************************************************************************/

/*****************************************************************************/
/** Calculate the mean of values in given double array.
    @return The mean of array values, or NaN in case of an error.
    @sa doubleSum
 */
double doubleMean(
  /** Pointer to double array. */
  double *a,
  /** Length of array. */
  const unsigned int n
) {
  if(a==NULL || n<1) return(nan(""));
  double s=0.0;
  unsigned int i, ci=0;
  for(i=0; i<n; i++) if(!isnan(a[i])) {ci++; s+=a[i];}
  if(ci<1) return(nan(""));
  return(s/(double)ci);
}
/*****************************************************************************/

/*****************************************************************************/
/** Returns the length of array consisting of only positive (>0 and not NaN)
    values.
   @return Returns the index of first nonpositive value in the given array.
   @sa doubleCSpanPositives
 */
int doubleSpanPositives(
  /** Pointer to the array. */
  double *a,
  /** Length of the array. */
  const int n
) {
  if(a==NULL || n<1) return(0);
  int i=0;
  for(i=0; i<n; i++) if(!(a[i]>0.0)) break;
  return(i);
}
/*****************************************************************************/

/*****************************************************************************/
/** Returns the length of array consisting of non-positive (<=0 and NaN) values.
   @return Returns the index of first positive value in the given array.
   @sa doubleSpanPositives
 */
int doubleCSpanPositives(
  /** Pointer to the array. */
  double *a,
  /** Length of the array. */
  const int n
) {
  if(a==NULL || n<1) return(0);
  int i=0;
  for(i=0; i<n; i++) if(a[i]>0.0) break;
  return(i);
}
/*****************************************************************************/

/*****************************************************************************/
/// @cond
/* Local functions */
static int statDoubleCompAsc(const void *i, const void *j)
{
  const double *di = (const double *)i;
  const double *dj = (const double *)j;
  return(*di > *dj) - (*di < *dj);
}
static int statDoubleCompDesc(const void *i, const void *j)
{
  const double *di = (const double *)i;
  const double *dj = (const double *)j;
  return(*di < *dj) - (*di > *dj);
}
/// @endcond

/** Sort the given double array into ascending or descending order.
 *  @author Vesa Oikonen
 */
void statSortDouble(
  /** Pointer to data array of size n */
  double *data,
  /** Length of data array */
  unsigned int n,
  /** Ascending (0) or descending (<>0) order */
  int order
) {
  if(n<2 || data==NULL) return;
  if(order==0) qsort(data, n, sizeof(double), statDoubleCompAsc);
  else qsort(data, n, sizeof(float), statDoubleCompDesc);
}
/*****************************************************************************/

/*****************************************************************************/
/// @cond
/* Local functions */
static int statFloatCompAsc(const void *i, const void *j)
{
  const float *di = (const float *)i;
  const float *dj = (const float *)j;
  return(*di > *dj) - (*di < *dj);
}
static int statFloatCompDesc(const void *i, const void *j)
{
  const float *di = (const float *)i;
  const float *dj = (const float *)j;
  return(*di < *dj) - (*di > *dj);
}
/// @endcond

/** Sort the given float array into ascending or descending order.
 *  @author Vesa Oikonen
 */
void statSortFloat(
  /** Pointer to data array of size n */
  float *data,
  /** Length of data array */
  unsigned int n,
  /** Ascending (0) or descending (<>0) order */
  int order
) {
  if(n<2 || data==NULL) return;
  if(order==0) qsort(data, n, sizeof(float), statFloatCompAsc);
  else qsort(data, n, sizeof(float), statFloatCompDesc);
}
/*****************************************************************************/

/*****************************************************************************/
