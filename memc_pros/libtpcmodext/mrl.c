/// @file mrl.c
/// @brief Calculation of maximum run length.
/// @author Vesa Oikonen
///
/*****************************************************************************/

/*****************************************************************************/
#include "libtpcmodext.h"
/*****************************************************************************/

/*****************************************************************************/
/** Return the maximum run length between given n length arrays of data */
int mrl_between_tacs(
  /** Array of data1; may contain NaNs */
  double *y1,
  /** Array of data2; may contain NaNs */
  double *y2,
  /** Nr of samples in array 1 and 2 */
  int n
) {
  int i, mrl=0, rl=0;
  char last_sign=0, sign;

  if(n<1 || y1==NULL || y2==NULL) return(0);
  for(i=0; i<n; i++) {
    if(isnan(y1[i]) || isnan(y2[i])) continue;
    if(y1[i]>y2[i]) sign=1; else if(y1[i]<y2[i]) sign=-1; else sign=0;
    if(sign!=last_sign) {
      rl=0; last_sign=sign;
    } else {
      if(sign!=0) {rl++; if(rl>mrl) mrl=rl;}
    }
  }
  return(mrl);
}
/*****************************************************************************/

/*****************************************************************************/
