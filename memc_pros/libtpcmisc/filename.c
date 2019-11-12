/// @file filename.c
/// @author Vesa Oikonen
/// @brief Functions for editing file names.
///
/*****************************************************************************/
#include "libtpcmisc.h"
/*****************************************************************************/

/*****************************************************************************/
/** Remove path from file name.
 *  @author Vesa Oikonen
 * */
void filenameRmPath(
  /** Pointer to string */
  char *s
) {
  if(s==NULL || strlen(s)<1) return;
  char *cptr;
  cptr=strrchr(s, '/'); if(cptr==NULL) cptr=strrchr(s, '\\');
  if(cptr==NULL) return;
  cptr++;
  int i, n=strlen(cptr);
  for(i=0; i<n; i++, cptr++) s[i]=*cptr;
  s[i]=(char)0;
  return;
}
/*****************************************************************************/

/*****************************************************************************/
/** Remove the last extension from file name.
 *  @author Vesa Oikonen
 *  @return 1 if extension was found (and removed), 0 if not.
 * */
int filenameRmExtension(
  /** Pointer to string */
  char *s
) {
  if(s==NULL || strlen(s)<1) return(0);
  char *cptr;
  cptr=strrchr(s, '.'); if(cptr==NULL) return(0);
  if(cptr[1]=='/' || cptr[1]=='\\') return(0);
  *cptr=(char)0;
  return(1);
}
/*****************************************************************************/

/*****************************************************************************/
/** Remove all extensions from file name.
 *  @author Vesa Oikonen
 * */
void filenameRmExtensions(
  /** Pointer to string */
  char *s
) {
  if(s==NULL || strlen(s)<1) return;
  while(filenameRmExtension(s)) {}
  return;
}
/*****************************************************************************/

/*****************************************************************************/
/** Check if string fname matches string key, which may contain wildcards
 *  ? and *.
 * @return 1 if strings do match and 0 if not.
 */
int fnmatch(
  /** filename that is evaluated */
  const char *fname,
  /** key string which may contain wildcards '?' and '*' */
  const char *key
) {
  if(fname==NULL || key==NULL) return(0);
  char *key_ptr=NULL, *fname_ptr=NULL;

  while((*key)&&(*fname)) {
    if((*key=='?')||(*key==*fname)) {
      key++; fname++;
    } else if(*key=='*') {
      if(*(key+1)==*fname) {key_ptr=(char*)key++; fname_ptr=(char*)fname+1;
      } else {
        fname++;
        if(*(key+1)=='?') {key_ptr=(char*)key++; fname_ptr=(char*)fname;}
      }
    } else if((key_ptr!=NULL) && (*fname_ptr)) {
      return(fnmatch(key_ptr, fname_ptr));
    } else {
      return(0);
    }
  }
  if((*fname)&&(key_ptr!=NULL)) {return(fnmatch(key_ptr, fname_ptr));}
  else {if(*key=='*') key++; return(*key==*fname);}
}
/*****************************************************************************/

/*****************************************************************************/
/** Case-independent check whether string fname matches string key,
 *  which may contain wildcards ? and *.
 *  @return 1 if strings do match and 0 if not.
 */
int fncasematch(
  /** filename that is evaluated */
  const char *fname,
  /** key string which may contain wildcards '?' and '*' */
  const char *key
) {
  if(fname==NULL || key==NULL) return(0);
  char *key_ptr=NULL, *fname_ptr=NULL;

  while((*key)&&(*fname)) {
    if((*key=='?')||(toupper((int)*key)==toupper((int)*fname))) {
      key++; fname++;
    } else if(*key=='*') {
      if(toupper((int)*(key+1))==toupper((int)*fname)) {
        key_ptr=(char*)key++; fname_ptr=(char*)fname+1;
      } else {
        fname++;
        if(*(key+1)=='?') {key_ptr=(char*)key++; fname_ptr=(char*)fname;}
      }
    } else if((key_ptr!=NULL) && (*fname_ptr)) {
      return(fnmatch(key_ptr, fname_ptr));
    } else {
      return(0);
    }
  }
  if((*fname)&&(key_ptr!=NULL)) {return(fnmatch(key_ptr, fname_ptr));}
  else {if(*key=='*') key++; return(toupper((int)*key)==toupper((int)*fname));}
}
/*****************************************************************************/

/*****************************************************************************/
/** @brief Get the last extension of a filename.
 *
 *  Extension(s) in pathname are not searched for. 
 *  Note that pointer points to the original string.
 *  
 *  @return Pointer to the filename extension starting with '.', or NULL if no 
 *  extension is found.
 */
char *filenameGetExtension(
  /** Pointer to string; string is not edited here. */
  char *s
) {
  if(s==NULL || strlen(s)<1) return((char*)NULL);
  /* Identify path */
  char *pptr=strrchr(s, '/'); if(pptr==NULL) pptr=strrchr(s, '\\');
  if(pptr==NULL) pptr=s; else pptr++;
  /* Search for the last '.' after the path */
  char *cptr=strrchr(pptr, '.'); if(cptr==NULL) return(cptr);
  /* If filename starts with '.', that is not counted as extension */
  if(strlen(pptr)==strlen(cptr)) return((char*)NULL);
  return(cptr);
}
/*****************************************************************************/

/*****************************************************************************/
/** @brief Get all extensions of a filename.
 *
 *  Extension(s) in pathname are not searched for. 
 *  Note that pointer points to the original string.
 *  
 *  @return Pointer to the filename extension starting with '.', or NULL if no 
 *  extension is found.
 */
char *filenameGetExtensions(
  /** Pointer to string; string is not edited here. */
  char *s
) {
  if(s==NULL || strlen(s)<1) return((char*)NULL);
  /* Identify path */
  char *pptr=strrchr(s, '/'); if(pptr==NULL) pptr=strrchr(s, '\\');
  if(pptr==NULL) pptr=s; else pptr++;
  /* Search for the last '.' after the path, ignoring first character
     in case filename starts with '.' */
  char *cptr=strchr(pptr+1, '.');
  return(cptr);
}
/*****************************************************************************/

/*****************************************************************************/
