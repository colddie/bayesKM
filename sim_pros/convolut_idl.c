/** @file convolut.c
 *  @brief Linear convolution for discrete data.
 */
/*****************************************************************************/
#include "tpcclibConfig.h"
/*****************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
/*****************************************************************************/
#include "tpccm.h"
/*****************************************************************************/

/*****************************************************************************/
/** @brief Calculates the convolution sum of a discrete real data set data[0..n-1] and a 
    discretized response function kernel[0..m].
    @details Convolution is not aware of the step size (default is 1); if step size is not 1
    (it usually isn't), the step size must be taken into account either when computing the kernel
    or by scaling the convolution sum.
    @remark This function does not use FFT, which would be faster with very large dataset,
    but slightly less precise.
    @sa tacInterpolateToEqualLengthFrames
    @return Function returns 0 when successful, or 1 if input data is not valid.
*/
int convolut_idl(int argc, float * argv[])
{
// (
//   /** Data array of length n-1 to be convolved, including any user-defined zero-padding. */
//   double *data,
//   /** Nr of data values. */
//   const int n,
//   /** Response function values in an array of length m. */
//   double *kernel,
//   /** Length of kernel array. */
//   const int m,
//   /** The convolved sum data is returned in out[0..n-1]; this must not overlap the input data. */
//   double *out
// ) {

  double *data;
  int     n;
  double *kernel;
  int     m;
  double *out;
printf("this is a test!");
  /* read in parameters */	
  data            =  (double *)  argv[0];
  n               =  *(int *)    argv[1];
  kernel          =  (double*)   argv[2];
  m               =  *(int *)    argv[3];
  out             =  (double*)   argv[4];

printf("this is a test!%d %d",n,m);
  if(n<1 || m<1 || data==NULL || kernel==NULL || out==NULL || out==data) return 1;

  for(int di=0; di<n; di++) out[di]=0.0;

  for(int di=m-1; di<n; di++) {
    int dj=di;
    for(int k=0; k<m; k++)
      out[di] += data[dj--] * kernel[k];
  }

  for(int di=0; di<n && di<m-1; di++) {
    int k=0;
    for(int dj=di; dj>=0; dj--)
      out[di] += data[dj] * kernel[k++];
  }
printf("this is a test!%f %f",out[0],out[200]);
  return 0;
}
/*****************************************************************************/

// /*****************************************************************************/
// /** Check whether given values have steady intervals, with the first interval starting at zero. 
//     Optionally return the interval, or in case of uneven intervals, the shortest interval.
//    @return Returns 1 if values are equidistant, and 0 if intervals are uneven,
//     or cannot be determined.
//  */
// int simIsSteadyInterval(
//   /** Array of values to check; must be sorted in ascending order. */
//   double *x,
//   /** Size of the data array; at least two. */
//   const int n,
//   /** Pointer for the interval, either the common or the shortest interval;
//       enter NULL, if not needed. */
//   double *f
// ) {
//   if(f!=NULL) *f=nan("");
//   if(x==NULL || n<2) return(0);
//   //if(!(fabs(x[0])>1.0E-12)) return(0); // on purpose, to catch NaN

//   double d, freq=nan("");
//   int isdif=0;
//   for(int i=1; i<n; i++) {
//     d=x[i]-x[i-1];
//     if(!(d>1.0E-12)) return(0); // interval must be > 0
//     if(isnan(freq)) {freq=d; continue;} // first interval
//     if(fabs(freq-d)<1.0E-12) continue; // same interval as before
//     isdif++; // different
//     if(d<freq) freq=d; // save smaller interval
//   }
//   if(f!=NULL) *f=freq;
//   if(isdif>0) return(0); // variable intervals
//   /* Check that first sample time is at the midpoint of interval starting from 0 */
//   d=0.5*freq; //printf("d=%g x[0]=%g\n", d, x[0]);
//   if(fabs(x[0]-d)>1.0E-06) return(0);
//   return(1);
// }
// /*****************************************************************************/

// /*****************************************************************************/
