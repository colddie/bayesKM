/// @file lts.c
/// @author Jussi Tohka, Kaisa Sederholm, Vesa Oikonen
/// @brief Least trimmed squares estimates for univariate location and variance.
///
///  The algorithm (exact) is described in  P.J. Rousseeuw
///  and A.M. Leroy: Robust Regression and Outlier Detection.
///  John Wiley & Sons 1987.
///  Algorithm from N. Wirth's book, implementation by N. Devillard.
///  This code in public domain.
///
/*****************************************************************************/
#include "libtpcmodel.h"
/*****************************************************************************/
/* local function definitions */
int ltsQSort(const void *par1, const void *par2);
/*****************************************************************************/

/*****************************************************************************/
/** Least trimmed squares estimates for univariate location and variance.
    Data samples are expected to be truly real valued (i.e too many samples
    having the same value might lead to problems. 
    The algorithm (exact) is described in
    P.J. Rousseeuw and A.M. Leroy: Robust Regression and Outlier Detection
    John Wiley & Sons 1987.
\return Returns 0, if successful.
 */
int least_trimmed_square(
  /** Input vector of n sample values; data samples are expected to be truly
   *  real valued (i.e too many samples having the same value might lead
   *  to problems.  */
  double data[],
  /** Number of samples */
  long int n,
  /** Output: Mean of sample values */
  double *mean,
  /** Output: Variance of sample values */
  double *variance
) {
  int i, j, h, h2;
  double score, best_score, loc, best_loc, old_sum, new_sum, medd;
  double old_power_sum, new_power_sum;
  double *scaled_data;

  h = n - n/2;
  h2 = n/2;

  qsort(data, n, sizeof(double),ltsQSort); 
  
  old_sum=old_power_sum=0.0;
  for(i=0; i<h; i++) {
    old_sum = old_sum + data[i];
    old_power_sum = old_power_sum + data[i]*data[i];
  }

  loc = old_sum/h;  
  /* For better understanding of the algorithm: 
    O(N^2) implementation of the algorithm would compute score as:
    score = 0.0;
    for(i = 0;i < h;i++) {
      score = score + (data[i] - loc)*(data[i] - loc);
    } 
    But there is a faster way to this: */
  score = old_power_sum - old_sum*loc; 
  best_score = score;
  best_loc = loc;
  for(j=1; j<h2+1; j++) {
    new_sum = old_sum - data[j-1] + data[h-1+j];
    old_sum = new_sum;
    loc = old_sum/h;
    new_power_sum = old_power_sum - data[j-1]*data[j-1] 
                  + data[h-1+j]*data[h-1+j];
    old_power_sum = new_power_sum;
    score = old_power_sum - old_sum*loc; 
    if(score < best_score) {
      best_score = score;
      best_loc = loc;
    }
  }  
  *mean = best_loc;

  /* For the variance, it is needed to calculate the ellipsoid covering
     one half of samples. This is not implemented optimally here because
     data has already been sorted. */
  scaled_data = malloc(n*sizeof(double));
  if(scaled_data == NULL) return(1);
  for(i=0; i<n; i++)
    scaled_data[i] = (data[i]-best_loc)*((h-1)/best_score)*(data[i]-best_loc);
  medd = dmedian(scaled_data, n);
  free(scaled_data);
  *variance = (best_score/(h-1))*(medd/CHI2INV_1);
  return(0);
}
/*****************************************************************************/

/*****************************************************************************/
/** Compares two numbers.
\return Returns the -1 if value1<value2, 1 if value1>value2 and 0 otherwise
*/
int ltsQSort(
  /** value nr 1*/
  const void *par1, 
  /** value nr 2*/
  const void *par2
) {
  if( *((double*)par1) < *((double*)par2)) return(-1);
  else if( *((double*)par1) > *((double*)par2)) return(1);
  else return(0);
}
/*****************************************************************************/

/*****************************************************************************/

