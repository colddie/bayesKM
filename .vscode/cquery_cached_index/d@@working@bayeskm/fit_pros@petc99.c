/// @file petc99.c
/// @author Vesa Oikonen, Calle Laakkonen
///
/// @brief This file contains the ISO C99 functions that are not yet available in
///  all C compilers but that are required by the PET library.
/// @note The functions are named as temp_functionname() to prevent problems with
///  compilers that already have these functions.
///
/*****************************************************************************/
#include "libtpcmisc.h"
/*****************************************************************************/

/*****************************************************************************/
/*!
 * int roundf(float e)  - Rounds up float e to nearest int
 *
 * @param e float value
 * @return rounded integer
 */
int temp_roundf(float e)
{
#if defined(__STDC_VERSION__) && __STD_VERSION__>=199901L
  return(roundf(e));
#else
  if(e<0.0) {
      return (int)(e-0.5);
  } else {
      return (int)(e+0.5);
  }
#endif
}
/****************************************************************************/

/****************************************************************************/

