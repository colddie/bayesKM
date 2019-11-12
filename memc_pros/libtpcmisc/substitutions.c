/// @file substitutions.c
/// @author Harri Merisaari, Vesa Oikonen
/// @brief Functions for ANSI substitutions for better compatibility.
///
/*****************************************************************************/
#include "libtpcmisc.h"
//#include <stddef.h>
/*****************************************************************************/

/*****************************************************************************/
#ifndef HAVE_STRDUP
/*!
 * Allocates memory and copies into
 * it the string addressed by s, including the terminating
 * character. User should free the allocated memory.
 *
 * @param s input string
 * @return pointer to allocated string
 */
char* strdup(
  const char *s)
{
  void *r;
  size_t length;
  length=strlen(s)+1;
  r=malloc(length);
  if(r==NULL) return NULL;
  return (char*)memcpy(r, s, length);
}
#endif
/*****************************************************************************/

/*****************************************************************************/
#ifndef HAVE_STRCASESTR
/** Case-insensitive version of strstr().
 *  @return a pointer to the beginning of the first occurrence, or NULL
 *  if not found.
 */
char *strcasestr(
  /** Pointer to string in which substring needle is seached */
  const char *haystack,
  /** Pointer to substring which is searced for in source string haystack */
  const char *needle
) {
  if(!haystack || !*haystack || !needle || !*needle) return 0;

  const char *s=haystack, *p=needle;
  do {
    if(!*p) return(char*)haystack;
    if((*p==*s) || (tolower(*p)==tolower(*s))) {
      p++; s++;
    } else {
      p=needle; if(!*s) return(NULL);
      s=++haystack;
    }
  } while(1);
  return *p ? NULL : (char*)haystack;
}
#endif // HAVE_STRCASESTR
/*****************************************************************************/

/*****************************************************************************/

