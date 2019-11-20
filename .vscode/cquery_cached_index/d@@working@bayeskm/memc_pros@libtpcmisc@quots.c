/// @file quots.c
/// @author Vesa Oikonen
/// @brief Functions for processing strings with quotation marks.
///
/****************************************************************************/
#include "libtpcmisc.h"
/*****************************************************************************/

/*****************************************************************************/
/** The strstr_noquotation() function returns a pointer to the first occurrence
 *  in the string pointed to by str1, excluding parts that are inside quotation
 *  marks "" or '', of the string pointed to by str2.
 *
\return Returns NULL pointer if no match is found, or if found, then pointer
    to the first occurrence.
 */ 
char *strstr_noquotation(
  /** Pointer to string to be searched */
  const char *str1,
  /** Pointer to string with quotation marks */
  const char *str2
) {
  unsigned int i, test_len;
  unsigned int single_quotation=0;
  unsigned int double_quotation=0;
  char *cptr;
  
  if(str1==NULL) return((char*)NULL);
  if(str2==NULL) return((char*)str1);
  test_len=strlen(str2); if(test_len<1) return((char*)str1);
  for(i=0, cptr=(char*)str1; i<strlen(str1); i++, cptr++) {
    if(*cptr=='\'') {
      if(single_quotation==0) single_quotation=1; else single_quotation=0;
      continue;
    }
    if(*cptr=='\"') {
      if(double_quotation==0) double_quotation=1; else double_quotation=0;
      continue;
    }
    if(single_quotation==1 || double_quotation==1) continue;
    if(strncmp(cptr, str2, test_len)==0) return(cptr);
  }
  return((char*)NULL);
}
/*****************************************************************************/

/*****************************************************************************/
/** Copy str2 to str1, removing any quotation marks around the string,
 *  and making sure that string fits to str2.
\return Returns the length of the new string str2.
 */
int strnCopyClean(
  /** Pointer to pre-allocated result string with length of at least maxlen
   *  characters, including NULL character */
  char *str1,
  /** Pointer to the original string; not changed in this function */
  const char *str2,
  /** Max length of str1, including the end NULL */
  int maxlen
) {
  char *cptr;
  int i;

  if(str1==NULL || maxlen<1) return(0);
  str1[0]=(char)0; if(str2==NULL || strlen(str2)<1) return(0);
  cptr=(char*)str2; cptr+=strspn(str2, "\"\'\t\n\r ");
  strncpy(str1, cptr, maxlen-1); str1[maxlen-1]=(char)0;
  i=strlen(str1); if(i<1) return(0);
  cptr=str1+(i-1);
  while(i>0) {
    if(*cptr!='\"' && *cptr!='\'' && *cptr!='\t' && *cptr!='\n' &&
       *cptr!='\r' && *cptr!=' ')
           break;
    i--; str1[i]=(char)0; cptr--;
  }
  return(i);
}
/*****************************************************************************/

/*****************************************************************************/

