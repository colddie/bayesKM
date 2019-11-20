/// @file intex.c
/// @author Vesa Oikonen, Calle Laakkonen
/// @brief Expansion of positive integers specified in a string.
///
/*****************************************************************************/
#include "libtpcmisc.h"
/*****************************************************************************/

/*****************************************************************************/

/** Call this (once) before using INT_LIST struct for the first time. */
void intInit(
  /** Pointer to INT_LIST struct */
  INT_list *l
) {
  if(l==NULL) return;
  l->nr=0;
  l->i=NULL;
  return;
}

/** Free the memory allocated in the INT_LIST struct. */
void intEmpty(
  /** Pointer to INT_LIST struct */
  INT_list *l
) {
  if(l==NULL) return;
  if(l->nr>0) free(l->i);
  l->nr=0;
  l->i=NULL;
  return;
}

/*!
 * Existing list is freed and all data is cleared. Deprecated.
 * Expanded integers are listed in list.i[] in increasing order.
 *
 * @param text Integer expressions to be expanded, e.g. 0-8,12,34-28
 * @param list Pointer for int list data
 * @return 0 if ok and at least one integer is listed.
 * @sa integerListAddFromString, integerListExpandFromString.
 */
int intExpand(char *text, INT_list *list) {
  int j;
  char *p, *t;
  int first, last, swap, intMax=65536;

  /* Check the arguments */
  if(strlen(text)<1) return(1);
  intEmpty(list);

  /* Expand */
  p=strtok(text, " ,;.&\t\n\r\0");
  while(p!=NULL) {
    t=p; first=last=-1;
    while((*t!='-') && (!isdigit((int)*t)) && (*t)) t++;
    if(*t=='-') {
      while((!isdigit((int)*t)) && (*t)) t++;
      if(isdigit((int)*t)) {first=0; last=atoi(t);}
    } else if(isdigit((int)*t)) {
      first=atoi(t); /*if (first==0) first=1;*/
      while((isdigit((int)*t)) && (*t)) t++;
      if(*t == '-') {
        t++; while((!isdigit((int)*t)) && (*t)) t++;
        if(isdigit((int)*t)) last=atoi(t); else last=intMax;
      }
    }
    if((first>=0) && (last>=0)) {
      if(first>last) {swap=first; first=last; last=swap;}
      if(last>intMax) {if(first<=intMax) last=intMax; else last=0;}
      for(j=first; j<=last && list->nr<intMax; j++)
        if(_intexadd(list, j)<0) return(2);
    } else if(first>=0) {
      if(first<=intMax) if(_intexadd(list, first)<0) return(3);
    } else if(last>=0) {
      if(last>=intMax) if(_intexadd(list, last)<0) return(4);
    }
    p=strtok(NULL, " ,;.&\t\n\r\0");
  }
  if(list->nr<1) return(1);
  return(0);
}
/*****************************************************************************/

/*****************************************************************************/
/*!
 * int _intexadd(int *list, int a) ; local function.  Deprecated.
 *
 * @param list
 * @param a
 */
int _intexadd(INT_list *list, int a) {
  int i, j, n;

  /* Check if list is yet empty */
  if(list->nr==0) {
    list->i=(int*)malloc(sizeof(int)); if(list->i==NULL) return(-1);
    /* Put the first integer to list and return */
    list->nr=1; list->i[0]=a; return(1);
  }
  n=list->nr;
  /* Check through the existing list */
  for(i=0; i<n; i++) {
    /* if it already is listed, just return */
    if(list->i[i]==a) return(0);
    /* make room for this integer */
    list->i=(int*)realloc(list->i, (n+1)*sizeof(int)); if(list==NULL) return(-1);
    if(list->i[i]>a) {for(j=n-1; j>=i; j--) list->i[j+1]=list->i[j]; break;}
  }
  list->i[i]=a; list->nr=n+1;
  return(1);
}
/*****************************************************************************/

/*****************************************************************************/
/*!
 * Merges two lists and returns the result. (the originals are not touched)
 * Duplicate entries are removed. Deprecated.
 * 
 * @param list1 The first list
 * @param list2 The second list
 * @return pointer to the new combined list.
 */
INT_list intMerge(INT_list *list1, INT_list *list2) {
  int r,count=0,l1=0,l2=0;
  INT_list newlist;
  int *tmplist;
  int found;
  tmplist=(int*)malloc(sizeof(int)*(list1->nr+list2->nr));
  while(l1<list1->nr || l2<list2->nr) {
    if(l1<list1->nr) {
      tmplist[count]=list1->i[l1];
      count++;
      l1++;
    }
    found=0;
    if(l2<list2->nr) {
      for(r=0;r<count;r++) {
        if(tmplist[r]==list2->i[l2]) found++;
      }
      if(found<1) {
        tmplist[count]=list2->i[l2]; count++;
      }
      l2++;
    }
  }
  newlist.i=(int*)malloc(sizeof(int)*count);
  memcpy(newlist.i,tmplist,count*sizeof(int));
  newlist.nr=count;
  free(tmplist);
  return newlist;
}

/*****************************************************************************/

/*****************************************************************************/
/** Call this (once) before using INTEGER_LIST struct for the first time.
 *  @sa integerListEmpty, integerListSort, integerListAddFromString, integerListExpandFromString.
 *  @return Returns <>0 in case of an error.
 * */
int integerListInit(
  /** Pointer to INTEGER_LIST struct */
  INTEGER_LIST *l
) {
  if(l==NULL) return -1;
  l->nr=l->_allocNr=0;
  l->list=NULL;
  return 0;
}

/** Free the memory allocated in the INTEGER_LIST struct.
 *  @sa integerListInit, integerListSort, integerListAdd.
 *  @return Returns <>0 in case of an error.
 * */
int integerListEmpty(
  /** Pointer to INTEGER_LIST struct */
  INTEGER_LIST *l
) {
  if(l==NULL) return -1;
  if(l->_allocNr>0) free(l->list);
  l->nr=l->_allocNr=0;
  l->list=NULL;
  return 0;
}

/** Add one integer to INTEGER_LIST.
   @sa integerListInit, integerListSort, integerListAddFromString, integerListExpandFromString.
   @return Returns the number of added integers (0 or 1), or <0 in case of an error.
 */
int integerListAdd(
  /** Pointer to initiated list */
  INTEGER_LIST *l,
  /** Integer value to add */
  int v,
  /** Add integer to the list only if it is new (0=no, 1=yes) */
  int ifnew
) {
  int i;
  if(l==NULL) return(-1);

  /* If only new value is to be added, then check that this is new */
  if(ifnew!=0) {for(i=0; i<l->nr; i++) if(l->list[i]==v) return(0);}
  /* Add value to list if there is space left, and quit */
  if(l->_allocNr>l->nr) {l->list[l->nr++]=v; return(1);}
  /* Allocate more space */
  if(l->_allocNr==0) l->list=(int*)malloc(10*sizeof(int));
  else l->list=(int*)realloc(l->list, (10+l->_allocNr)*sizeof(int));
  if(l->list==NULL) {l->nr=l->_allocNr=0; return(-2);}
  l->_allocNr+=10;
  /* Add value to list */
  l->list[l->nr++]=v;
  return(1);
}

/** Sort INTEGER_LIST 
 *  @sa integerListInit, integerListAdd, integerListAddFromString, integerListExpandFromString.
 *  @return Returns <>0 in case of an error.
 */
int integerListSort(
  /** Pointer to INTEGER_LIST struct */
  INTEGER_LIST *l
) {
  int i, j, v;

  if(l==NULL) return -1;
  if(l->nr<2) return 0;
  for(i=0; i<l->nr; i++) for(j=i+1; j<l->nr; j++) {
    if(l->list[i]>l->list[j]) {
      v=l->list[i]; l->list[i]=l->list[j]; l->list[j]=v;
    }
  }
  return 0;
}
/*****************************************************************************/

/*****************************************************************************/
/** Read a list of integer values from given string with given delimiters.
 *  @sa integerListInit, integerListSort, integerListAdd, integerListExpandFromString.
 *  @return The number of added integer values, or <0 in case of an error.
 *  @author Vesa Oikonen
 */
int integerListAddFromString(
  /** Pointer to string from which the integers are read, for example
   *  "2,3,6,8". */
  const char *s1,
  /** String containing character delimiters, for example ", ". */
  const char *s2,
  /** Pointer to INTEGER_LIST struct; previous contents are preserved. */
  INTEGER_LIST *l,
  /** Add integer to the list only if it is new (0=no, 1=yes) */
  const int ifnew
) {
  if(l==NULL) return(-1);
  if(s1==NULL || s2==NULL) return(0);
  
  /* Get the nr of tokens */
  int i, j, m, n, v;
  n=strTokenNr((char*)s1, s2); if(n<1) return(0);
  /* Read the values */
  char tmp[128];
  for(i=j=0; i<n; i++) {
    if(strTokenNCpy(s1, s2, 1+i, tmp, 128)<1) return(-2);
    if(atoi_with_check(tmp, &v)) return(-3);
    m=integerListAdd(l, v, ifnew); if(m<0) return(-4);
    j+=m;
  }
  return(j);
}
/*****************************************************************************/

/*****************************************************************************/
/** Read ranges and individual integer values from given string with given 
 *  delimiters.
 *  @sa integerListInit, integerListSort, integerListAdd, integerListAddFromString.
 *  @return The number of added integer values, or <0 in case of an error.
 *  @author Vesa Oikonen
 */
int integerListExpandFromString(
  /** Pointer to string from which the integers are read, for example
   *  "0-8,12,32-28" or "0..8, 12, 28..34". */
  const char *s1,
  /** String containing character delimiters, for example ", ". */
  const char *s2,
  /** Pointer to INTEGER_LIST struct; previous contents are preserved. */
  INTEGER_LIST *l,
  /** Add integer to the list only if it is new (0=no, 1=yes) */
  const int ifnew
) {
  if(l==NULL) return(-1);
  if(s1==NULL || s2==NULL) return(0);
  
  /* Get the nr of tokens */
  int n=strTokenNr((char*)s1, s2); if(n<1) return(0);
  /* Read the values */
  char tmp[128], tmp2[128], *t, *tail;
  int i, j, m, first, last, sw;
  for(i=j=0; i<n; i++) {
    if(strTokenNCpy(s1, s2, 1+i, tmp, 128)<1) return(-2);
    t=tmp; errno=0;
    first=strtol(t, &tail, 10); if(errno) return(-3);
    if(*tail) {
      strcpy(tmp2, tail); t=tmp2;
      if(*t=='-') t++;
      else if(*t=='.') {t++; if(*t=='.') t++; else return(-4);}
      else return(-5);
      if(!*t) return(-6);
      last=strtol(t, &tail, 10); if(errno) return(-7);
      if(*tail) return(-8);
    } else {
      last=first;
    }

    if(first>last) {sw=first; first=last; last=sw;}
    for(int v=first; v<=last; v++) {
      m=integerListAdd(l, v, ifnew); if(m<0) return(-10);
      j+=m;
    }
  }
  return(j);
}
/*****************************************************************************/

/*****************************************************************************/
