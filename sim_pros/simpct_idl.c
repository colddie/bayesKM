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


int simpct_idl(int argc, float * argv[])
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

  double *ts;
  double *ctt;
  int    frameNr;
  double cbf;
  double mtt;
  double *tac;

  /* read in parameters */	
  ts            =  (double *)    argv[0];
  ctt           =  (double *)    argv[1];
  frameNr       =  *(int *)      argv[2];
  cbf           =  *(double *)   argv[3];
  mtt           =  *(double*)    argv[4];
  tac           =  (double *)    argv[5];

  cbf         = cbf/6000.0;              //% ml/100ml/min = 1/(100*60s)
  double  data[frameNr];
  int     n = frameNr;
  int     m = frameNr;
  for (int i=0;i<n;i++) { 
    data[i]=cbf*exp( -(ts[i]-mtt)); 
    if (ts[i] < mtt) { data[i] = cbf; }
    }

  if(n<1 || m<1 || data==NULL || ctt==NULL || tac==NULL || tac==data) { return 1; }

  for(int di=0; di<n; di++) tac[di]=0.0;

  for(int di=m-1; di<n; di++) {
    int dj=di;
    for(int k=0; k<m; k++)
      tac[di] += data[dj--] * ctt[k];
  }

  for(int di=0; di<n && di<m-1; di++) {
    int k=0;
    for(int dj=di; dj>=0; dj--)
      tac[di] += data[dj] * ctt[k++];
  }

 cbf         = cbf*6000.0;    // need to scale back to original!  
// printf("test here! %d %f %f",frameNr,cbf,mtt);
// printf("test here! %f %f %f %f %f %f %f %f %f",tac[0],tac[50],tac[100],tac[150],tac[200],tac[250],tac[300],tac[400],tac[450]);
  return 0;
}
/*****************************************************************************/

