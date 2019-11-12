/// @file strext.c
/// @author Vesa Oikonen
/// @brief Utility functions for processing strings.
///
/*****************************************************************************/
#include "libtpcmisc.h"
#include <string.h>
/*****************************************************************************/

/*****************************************************************************/
/** The strTokenNr() function returns the number of tokens in the string pointed to by str1. 
    The characters making up the string pointed to by str2 are the delimiters that 
    determine the token.
    @sa strTokenNCpy, strChrCount
    @return Returns the nr of tokens.
 */  
int strTokenNr(
  /** String from where tokens are calculated; not modified in any way. */
  const char *str1,
  /** String containing character delimiters. */
  const char *str2
) {
  int i=0, n=0;
  char *cptr;
  if(str1==NULL || str2==NULL || strlen(str1)==0 || strlen(str2)==0) return(0);

  cptr=(char*)str1;
  do {
    // pass delimiter characters
    i=strspn(cptr, str2); cptr+=i;
    // pass characters between delimiters
    i=strcspn(cptr, str2); cptr+=i; if(i>0) n++; 
  } while(i>0);
  return(n);
}
/*****************************************************************************/

/*****************************************************************************/
/** The strTokenNCpy() function copies the i'th token in the string pointed to
    by str1 into string pointed to by str3. The characters making up the string
    pointed to by str2 are the delimiters that determine the token.
    @sa strTokenNr, strTokenDup
    @return Returns the length of token, 0 if no token(s) found.
 */  
int strTokenNCpy(
  /** String from where tokens are searched; not modified in any way. */
  const char *str1,
  /** String containing character delimiters. */
  const char *str2,
  /** Token number to copy (1..nr of tokens). */
  int i,
  /** String array into where the token is copied; string will be null terminated. */
  char *str3,
  /** Length of str3, including terminal null */
  int count
) {
  int j=0, n=0;
  char *cptr;
  if(str1==NULL || str2==NULL || strlen(str1)==0 || strlen(str2)==0) return(0);
  if(i<1 || str3==NULL || count<2) return(0);

  cptr=(char*)str1;
  do {
    // pass delimiter characters
    j=strspn(cptr, str2); cptr+=j;
    // pass characters between delimiters
    j=strcspn(cptr, str2); if(j>0) n++;
    // if this is the required token nr, then stop here
    if(n==i) {
      if(j>count-1) j=count-1;
      strncpy(str3, cptr, j); str3[j]=(char)0;
      break;
    }
    cptr+=j;
  } while(j>0);
  if(n>i) {str3[0]=(char)0; return(0);}
  return(j);
}
/*****************************************************************************/

/*****************************************************************************/
/** Search the string s1 for the first token. The characters making
    up the string s2 are the delimiters that determine the tokens.
    @return Returns pointer to a copy of the token string, or NULL in case of 
    an error or if no token found.
    @post Remember to free the memory from the returned pointer after last use.
    @author Vesa Oikonen
 */  
char *strTokenDup(
  /** String from where tokens are searched; not modified in any way. */
  const char *s1,
  /** String containing character delimiters. */
  const char *s2,
  /** Index of s1 where the token ended; set to NULL, if not needed. */
  int *next
) {
  if(next!=NULL) *next=0;
  if(s1==NULL) return NULL;

  char *s3=NULL, *cptr;
  size_t j;
  
  /* If no delimiters, then return copy of s1 */
  if(s2==NULL || strlen(s2)<1) {
    s3=strdup(s1); if(next!=NULL) *next=strlen(s1);
    return s3;
  } 
  /* Pass initial delimiter characters */
  cptr=(char*)s1; j=strspn(cptr, s2); cptr+=j; if(next!=NULL) *next=j;
  /* calculate characters between delimiters */
  j=strcspn(cptr, s2); if(j==0) {return NULL;}
  if(next!=NULL) *next+=j;
  /* Allocate space for token */
  s3=calloc(j+1, sizeof(char)); if(s3==NULL) return NULL;
  strlcpy(s3, cptr, j+1);
  return s3;
}
/*****************************************************************************/

/*****************************************************************************/
/** Count how many times specified characters are found in a string.
    Search is case-sensitive.
    @sa strTokenNr, strReplaceChar
    @return Returns the nr of characters of str2 found in str1. 
 */
int strChrCount(
  /** String to search for characters; not modified. */
  const char *str1,
  /** String containing characters which are searched for; not modified. */
  const char *str2
) {
  unsigned int n=0, i, j;
  if(str1==NULL || str2==NULL || strlen(str1)==0 || strlen(str2)==0) return n;
  for(i=0; i<strlen(str1); i++)
    for(j=0; j<strlen(str2); j++)
      if(str1[i]==str2[j]) n++;
  return n;
}
/*****************************************************************************/

/*****************************************************************************/
/** Replace certain characters in string with another character. */
void strReplaceChar(
  /** Pointer to string in which the character is replaced. */
  char *str,
  /** Character to be replaced. */
  char c1,
  /** Character to use instead. If NULL, then only the first character is replaced. */
  char c2
) {
  char *cptr;
  if(strlen(str)==0) return;
  while((cptr=strchr(str, c1))!=NULL) *cptr=c2;
  return;
}
/*****************************************************************************/

/*****************************************************************************/
#ifndef HAVE_STRNLEN
/** Safer version of strlen, in case the argument s is not NUL terminated string. 
    Computes the length of string s, but never scans beyond the n first bytes of the string. 
    @remark Included in POSIX and GCC, so this implementation may not be needed.
    @return same as strlen() or n, whichever is smaller. 
 */
size_t strnlen(
  /** Pointer to string, or character array, that may not be NULL terminated. */
  const char *s,
  /** The actual length of buffer allocated for the string;
      for example, string could have been allocated as char s[n]; */ 
  size_t n
) {
  if(s==NULL) return(0);
  char *ps=(char*)s;
  size_t i=0;
  while(i<n && *ps!='\0') {i++; ps++;}
  return(i);
}
#endif // HAVE_STRNLEN
/*****************************************************************************/

/*****************************************************************************/
#ifndef HAVE_STRLCAT
/** Safer version of strncat. 
    At most dstsize-1 characters are appended from the source string to destination string.
    Destination string will be NUL terminated, unless dstsize <= strlen(dst). 
    @remark Included in POSIX but not in GCC.
    @return the size of the buffer that would have been needed for the destination string; 
    if >=dstsize, then truncation occurred. 
 */
size_t strlcat(
  /** Destination string. */
  char *dst,
  /** Source string. */
  const char *src,
  /** The actual length of buffer allocated for the destination string;
      for example, destination string has been allocated as char dst[dstsize]; */
  size_t dstsize
) {
  char *d;
  const char *s=src;
  size_t dlen, n;

  /* Find the current length of dst */
  dlen=strnlen(dst, dstsize);
  if(s==NULL) return(dlen);
  n=dstsize-dlen;
  if(n==0) return(dlen+strlen(s));
  d=dst+dlen;
  while(*s!='\0') {
    if(n!=1) {*d=*s; d++; n--;}
    s++;
  }
  *d='\0';
  return(dlen+(s-src));
}
#endif // HAVE_STRLCAT
/*****************************************************************************/

/*****************************************************************************/
#ifndef HAVE_STRLCPY
/** Safer version of strncpy or strcpy.

    At most dstsize-1 characters are copied from the source string to the destination string. 
    Destination string will be NUL terminated. 
    @remark Included in POSIX but not in GCC.
    @return the size of the buffer that would have been needed for the destination string; 
    if >=dstsize, then truncation occurred. 
 */
size_t strlcpy(
  /** Destination string. */
  char *dst,
  /** Source string. */
  const char *src,
  /** The actual length of buffer allocated for the destination string;
      for example, destination string has been allocated as char dst[dstsize]; */
  size_t dstsize
) {
  if(src==NULL) return(0);

  char *d=dst;
  const char *s=src;
  size_t n;

  /* Copy as many chars as allowed */
  n=dstsize;
  if(n!=0) while(--n!=0) {*d=*s; if(*d=='\0') {d++; s++; break;} d++; s++;}
  if(n==0) { // not enough space, add NUL, and check how much space were needed
    if(dstsize!=0) *d='\0';
    while(*s++) {}
  }
  return(s-src-1);
}
#endif // HAVE_STRLCPY
/*****************************************************************************/

/*****************************************************************************/
/** Version of strncpy() which as usual copies s2 to s1, but without any
    space characters or line end characters that may be around the string s2.
    @return the length of the new string s1.
    @author Vesa Oikonen
 */
int strncpyCleanSpaces(
  /** Pointer to pre-allocated result string with length of at least maxlen
   *  characters, including NULL character. */
  char *s1,
  /** Pointer to the original string. */
  const char *s2,
  /** Max length of s1, including the trailing zero. */
  int maxlen
) {
  if(s1==NULL) return(0);
  s1[0]=(char)0; if(maxlen<1) return(0);
  if(maxlen<2) {strcpy(s1, ""); return(0);}
  if(s2==NULL || strlen(s2)<1) return(0);

  char *cptr;
  int i;

  cptr=(char*)s2; cptr+=strspn(s2, "\t\n\r ");
  strlcpy(s1, cptr, maxlen); i=strlen(s1); if(i<1) return(0);
  cptr=s1+(i-1);
  while(i>0) {
    if(*cptr!='\t' && *cptr!='\n' && *cptr!='\r' && *cptr!=' ') break;
    i--; s1[i]=(char)0; cptr--;
  }
  return(i);
}
/*****************************************************************************/

/*****************************************************************************/
/** Removes any initial and trailing space characters from specified string s.

    Space characters in the middle of the string are not removed.
    @return 0 when successful, otherwise >0.
    @author Vesa Oikonen
 */
int strCleanSpaces(
  /** Pointer to the string. */
  char *s
) {
  if(s==NULL) return 0;
  int len=strlen(s); if(len<0) return 0;
  char *s2; s2=strdup(s); if(s2==NULL) return(1);
  int n=strncpyCleanSpaces(s2, s, len+1);
  if(n<1) strcpy(s, ""); else strcpy(s, s2);
  free(s2);
  return 0;
}
/*****************************************************************************/

/*****************************************************************************/
