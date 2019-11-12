/// @file iftsrch.c
/// @author Vesa Oikonen
/// @brief Search functions for IFT contents.
///
/******************************************************************************/
#include "libtpcmisc.h"
/******************************************************************************/

/******************************************************************************/
/*!
 * Find the key in the IFT and return the index [0..keyNr-1].
 * Key is case insensitive.
 * 
 * @param ift Pointer to existing IFT
 * @param key Pointer to the key string; contents are replaced by
 * the correct key string
 * @return -1 if key was not found, or other negative value in case of an error
 */
int iftGet(IFT *ift, char *key) {
  int li;

  if(IFT_TEST) printf("iftGet(*ift, \"%s\")\n", key);
  if(ift==NULL) {return(-10);}
  if(key==NULL || strlen(key)<1) {iftSetStatus(ift, IFT_FAULT); return(-11);}
  
  /* Search the list */
  for(li=0; li<ift->keyNr; li++) {
    if(strcasecmp(ift->item[li].key, key)==0) {
      strcpy(key, ift->item[li].key);
      iftSetStatus(ift, IFT_OK); return(li);
    }
  }
  iftSetStatus(ift, IFT_KEYNOTFOUND);
  return(-1);
}
/******************************************************************************/

/******************************************************************************/
/*!
 * Find the Nth key with similar name in the IFT and return the index
 * [0..keyNr-1]. Key is case insensitive.
 *
 * @param ift Pointer to existing IFT
 * @param key Pointer to the key string; contents are replaced by
 * the correct key string
 * @param n Nth (1..) incidence of key is searched.
 * @return -1 if key was not found, or other negative value in case of an error
 */
int iftGetNth(IFT *ift, char *key, int n) {
  int li, found_nr=0;

  if(IFT_TEST) printf("iftGetNth(*ift, \"%s\", %d)\n", key, n);
  if(ift==NULL) {return(-10);}
  if(key==NULL || strlen(key)<1) {iftSetStatus(ift, IFT_FAULT); return(-11);}
  if(n<1) {iftSetStatus(ift, IFT_FAULT); return(-11);}
  
  /* Search the list */
  for(li=0; li<ift->keyNr; li++) {
    if(strcasecmp(ift->item[li].key, key)==0) {
      strcpy(key, ift->item[li].key); found_nr++;
      if(n==found_nr) {iftSetStatus(ift, IFT_OK); return(li);}
    }
  }
  iftSetStatus(ift, IFT_KEYNOTFOUND);
  return(-1);
}
/******************************************************************************/

/******************************************************************************/
/*!
 * Find the Nth item of IFT where the specified string is found in the key.
 * Comparison is case sensitive.
 * 
 * @param ift Pointer to existing IFT
 * @param str Pointer to the case-sensitive (partial) key string
 * @param n Nth (1..keyNr-1) incidence of value is searched.
 * @return -1 if key was not found, or other negative value in
 case of an error, and the index [0..keyNr-1] if matching key is found.
 */
int iftFindNthKey(IFT *ift, char *str, int n) {
  int li, found_nr=0;

  if(IFT_TEST) printf("iftFindNthKey(*ift, \"%s\", %d)\n", str, n);
  if(ift==NULL) {return(-10);}
  if(str==NULL || strlen(str)<1) {iftSetStatus(ift, IFT_FAULT); return(-11);}
  if(n<1) {iftSetStatus(ift, IFT_FAULT); return(-11);}
  
  /* Search the list */
  for(li=0; li<ift->keyNr; li++) {
    if(strstr(ift->item[li].key, str)!=NULL) {
      found_nr++;
      if(n==found_nr) {iftSetStatus(ift, IFT_OK); return(li);}
    }
  }
  iftSetStatus(ift, IFT_KEYNOTFOUND);
  return(-1);
}
/******************************************************************************/

/******************************************************************************/
/*!
 * Find the Nth item of IFT where the specified string is found in the value.
 *   Comparison is case sensitive.
 *
 * @param ift Pointer to existing IFT
 * @param str Pointer to the case-sensitive (partial) value string
 * @param n Nth (1..keyNr-1) incidence of value is searched.
 * @return -1 if key was not found, or other negative value in
 * case of an error, and the index [0..keyNr-1] if matching value is found.
 */
int iftFindNthValue(IFT *ift, char *str, int n) {
  int li, found_nr=0;

  if(IFT_TEST) printf("iftFindNthValue(*ift, \"%s\", %d)\n", str, n);
  if(ift==NULL) {return(-10);}
  if(str==NULL || strlen(str)<1) {iftSetStatus(ift, IFT_FAULT); return(-11);}
  if(n<1) {iftSetStatus(ift, IFT_FAULT); return(-11);}
  
  /* Search the list */
  for(li=0; li<ift->keyNr; li++) {
    if(strstr(ift->item[li].value, str)!=NULL) {
      found_nr++;
      if(n==found_nr) {iftSetStatus(ift, IFT_OK); return(li);}
    }
  }
  iftSetStatus(ift, IFT_VALUENOTFOUND);
  return(-1);
}
/******************************************************************************/

/******************************************************************************/
/** Finds the specified key in the IFT, starting from specified index.
 *  Key is case insensitive.
\return Returns the index of key/value, -1 if key or value was not found,
    and <-1 in case of an error.
 */
int iftGetFrom(
  /** Pointer to existing IFT */
  IFT *ift,
  /** Index [0..keyNr-1] from which the search is started */
  int si,
  /** Pointer to the key string; search is case-insensitive */
  const char *key
) {
  int li;

  if(IFT_TEST) printf("iftGetFrom(*ift, %d, \"%s\")\n", si, key);
  if(ift==NULL) {return(-10);}
  if(key==NULL || strlen(key)<1) {iftSetStatus(ift, IFT_FAULT); return(-11);}
  if(si<0) {iftSetStatus(ift, IFT_FAULT); return(-12);}
  
  /* Search the list */
  for(li=si; li<ift->keyNr; li++) {
    if(strcasecmp(ift->item[li].key, key)==0) {
      iftSetStatus(ift, IFT_OK); return(li);
    }
  }
  iftSetStatus(ift, IFT_KEYNOTFOUND);
  return(-1);
}
/******************************************************************************/

/******************************************************************************/
/** Finds the index with specified key and value in the IFT, starting from
 *  specified index. Key and value are case insensitive.
\return Returns the index of key/value, -1 if key or value was not found,
    and <-1 in case of an error.
 */
int iftGetFullmatchFrom(
  /** Pointer to existing IFT */
  IFT *ift,
  /** Index [0..keyNr-1] from which the search is started */
  int si,
  /** Pointer to the key string; search is case-insensitive */
  const char *key,
  /** Pointer to the value string; search is case-insensitive */
  const char *value
) {
  int li;

  if(IFT_TEST)
    printf("iftGetFullmatchFrom(*ift, %d, \"%s\", \"%s\")\n", si, key, value);
  if(ift==NULL) {return(-10);}
  if(key==NULL) {iftSetStatus(ift, IFT_FAULT); return(-11);}
  if(value==NULL) {iftSetStatus(ift, IFT_FAULT); return(-12);}
  if(si<0) {iftSetStatus(ift, IFT_FAULT); return(-13);}
  
  /* Search the list */
  for(li=si; li<ift->keyNr; li++) {
    if(strcasecmp(ift->item[li].key, key)!=0) continue;
    if(strcasecmp(ift->item[li].value, value)!=0) continue;
    iftSetStatus(ift, IFT_OK); return(li);
  }
  iftSetStatus(ift, IFT_KEYNOTFOUND);
  return(-1);
}
/******************************************************************************/

/******************************************************************************/
/** Finds the specified key string from IFT struct, and reads the corresponding
 *  value as float.
\return Returns the index of key/value, -1 if key or value was not found,
    and <-1 in case of an error.
 */
int iftGetFloatValue(
  /** Pointer to existing IFT */
  IFT *ift,
  /** Index [0..keyNr-1] from which the search is started */
  int si,
  /** Pointer to the key string; search is case-insensitive */
  const char *key,
  /** Pointer to float variable where value is written; NaN is written in case
   *  of an error. */
  float *value
) {
  int li;

  if(IFT_TEST) printf("iftGetFloatValue(*ift, \"%s\", *value)\n", key);
  if(ift==NULL) {return(-10);}
  if(value==NULL) {iftSetStatus(ift, IFT_FAULT); return(-10);}
  *value=nanf("");
  if(key==NULL || strlen(key)<1) {iftSetStatus(ift, IFT_FAULT); return(-11);}
  if(si<0) {iftSetStatus(ift, IFT_FAULT); return(-12);}
  iftSetStatus(ift, IFT_VALUENOTFOUND);
  /* Search the list */
  for(li=si; li<ift->keyNr; li++) {
    if(strcasecmp(ift->item[li].key, key)==0) {
      if(ift->item[li].value==NULL || strlen(ift->item[li].value)<1) return -1;
      (void)sscanf(ift->item[li].value, "%f", value);
      if(isnan(*value)) return -1;
      iftSetStatus(ift, IFT_OK); return(li);
    }
  }
  iftSetStatus(ift, IFT_KEYNOTFOUND);
  return(-1);
}
/******************************************************************************/

/******************************************************************************/
/** Finds the specified key string from IFT struct, and reads the corresponding
 *  value as double.
\return Returns the index of key/value, -1 if key or value was not found,
    and <-1 in case of an error.
 */
int iftGetDoubleValue(
  /** Pointer to existing IFT */
  IFT *ift,
  /** Index [0..keyNr-1] from which the search is started */
  int si,
  /** Pointer to the key string; search is case-insensitive */
  const char *key,
  /** Pointer to double variable where value is written; NaN is written in case
   *  of an error. */
  double *value
) {
  int li;

  if(IFT_TEST) printf("iftGetDoubleValue(*ift, \"%s\", *value)\n", key);
  if(ift==NULL) {return(-10);}
  if(value==NULL) {iftSetStatus(ift, IFT_FAULT); return(-10);}
  *value=nan("");
  if(key==NULL || strlen(key)<1) {iftSetStatus(ift, IFT_FAULT); return(-11);}
  if(si<0) {iftSetStatus(ift, IFT_FAULT); return(-12);}
  iftSetStatus(ift, IFT_VALUENOTFOUND);
  /* Search the list */
  for(li=si; li<ift->keyNr; li++) {
    if(strcasecmp(ift->item[li].key, key)==0) {
      if(ift->item[li].value==NULL || strlen(ift->item[li].value)<1) return -1;
      (void)sscanf(ift->item[li].value, "%lf", value);
      if(isnan(*value)) return -1;
      iftSetStatus(ift, IFT_OK); return(li);
    }
  }
  iftSetStatus(ift, IFT_KEYNOTFOUND);
  return(-1);
}
/******************************************************************************/

/******************************************************************************/
/** Finds the specified key string from IFT struct, and reads the corresponding
 *  value as int.
\return Returns the index of key/value, -1 if key or value was not found,
    and <-1 in case of an error.
 */
int iftGetIntValue(
  /** Pointer to existing IFT */
  IFT *ift,
  /** Index [0..keyNr-1] from which the search is started */
  int si,
  /** Pointer to the key string; search is case-insensitive */
  const char *key,
  /** Pointer to int variable where value is written; -9999 is written in case
   *  of an error. */
  int *value
) {
  int li;

  if(IFT_TEST) printf("iftGetFloatValue(*ift, \"%s\", *value)\n", key);
  if(ift==NULL) {return(-10);}
  if(value==NULL) {iftSetStatus(ift, IFT_FAULT); return(-10);}
  *value=-9999;
  if(key==NULL || strlen(key)<1) {iftSetStatus(ift, IFT_FAULT); return(-11);}
  if(si<0) {iftSetStatus(ift, IFT_FAULT); return(-12);}
  /* Search the list */
  iftSetStatus(ift, IFT_VALUENOTFOUND);
  for(li=si; li<ift->keyNr; li++) {
    if(strcasecmp(ift->item[li].key, key)==0) {
      if(ift->item[li].value==NULL || strlen(ift->item[li].value)<1) return -1;
      (void)sscanf(ift->item[li].value, "%d", value);
      if(*value==-9999) return -1;
      iftSetStatus(ift, IFT_OK); return(li);
    }
  }
  iftSetStatus(ift, IFT_KEYNOTFOUND);
  return(-1);
}
/******************************************************************************/

/******************************************************************************/
/** Find the nr of occurrences of the specified key in the IFT.
 *  Key is case insensitive.
 *
\return Returns the nr of occurrences of the key.
 */
int iftGetKeyNr(
  /** Pointer to existing IFT */
  IFT *ift,
  /** Pointer to the key string */
  const char *key
) {
  int li, found_nr=0;

  if(IFT_TEST) printf("iftGetKeyNr(*ift, \"%s\")\n", key);
  if(ift==NULL) {return(0);}
  if(key==NULL || strlen(key)<1) {iftSetStatus(ift, IFT_FAULT); return(0);}
  iftSetStatus(ift, IFT_KEYNOTFOUND);
  
  /* Search the list */
  for(li=0; li<ift->keyNr; li++)
    if(strcasecmp(ift->item[li].key, key)==0) found_nr++;
  if(found_nr>0) iftSetStatus(ift, IFT_OK);
  else iftSetStatus(ift, IFT_KEYNOTFOUND);
  return(found_nr);
}
/*****************************************************************************/

/*****************************************************************************/
