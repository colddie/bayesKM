/// @file lms.c
/// @author Kaisa Sederholm, Vesa Oikonen
/// @brief Least median of squares estimate for single data.
///
/*****************************************************************************/
#include "libtpcmodel.h"
/*****************************************************************************/
/* local function definitions */
/// @cond
int lmsQSort(const void *par1, const void *par2);
/// @endcond
/*****************************************************************************/
/** @brief Fit a constant (horisontal straight line) to the data by minimising
    the median of squared residuals.

    The algorithm is described in
    P.J. Rousseeuw: Least Median of Squares Regression, Journal of the 
    American Statistical Association, Vol. 79, No. 388 (1984), 871-880.
\return Returns the LMS estimate.
 */
double least_median_of_squares(
  /** Data array */
  double *data,
  /** Number of data values */
  int n
) {
  int i, odd=1, half=floor(n/2), smallnr;
  double small, *help;
  if(fmod((double)n, 2.0)<1e-99) odd=0;
  double *halfs_data;

  halfs_data=(double*)malloc((half+odd) * sizeof(double)); 

  help=data;
  /* sort data in ascending order*/
  qsort(help, n, sizeof(double), lmsQSort); 

  /* if n is even number */
  for(i=0; i<half+odd; i++) {
    halfs_data[i]=data[half+i]-data[i];
  }

  i=smallnr=0;
  for(i=1, small=halfs_data[0]; i<half+odd; i++) {
    if(halfs_data[i]<small) {
      small=halfs_data[i];
      smallnr=i;
    }
  }

  return (data[half+smallnr]+data[smallnr])/2.0;
}
/*****************************************************************************/

/*****************************************************************************/
/** Compares two numbers.
\return Returns -1 if value1<value2, 1 if value1>value2 and 0 otherwise.
*/
int lmsQSort(
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
