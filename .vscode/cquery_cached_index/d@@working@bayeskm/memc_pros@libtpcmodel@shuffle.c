/// @file shuffle.c
/// @author Vesa Oikonen
/// @brief Random shuffle and related functions.
///
/******************************************************************************/
#include "libtpcmodel.h"
/******************************************************************************/

/******************************************************************************/
/** Random shuffle: arrange the n elements of array in random order.
 *  Only effective if N is much smaller than RAND_MAX.
 */
void random_shuffle(int *array, int n)
{
  if(n<=1 || array==NULL) return;
  int i, j, tmp;
  for(i=0; i<n-1; i++) {
    j=i+rand()/(RAND_MAX/(n-i)+1);
    tmp=array[j]; array[j]=array[i]; array[i]=tmp;
  }
}
/******************************************************************************/

/******************************************************************************/
/** Random permutation: given (allocated) array of length n is filled with
 *  random permutation of numbers in the range [a:n+a-1]; that is, each number
 *  once in random order.
 */
void randperm(int *array, int n, int a)
{
  if(n<1 || array==NULL) return;
  int i;
  for(i=0; i<n; i++) array[i]=i+a;
  random_shuffle(array, n);
}
/******************************************************************************/

/******************************************************************************/
