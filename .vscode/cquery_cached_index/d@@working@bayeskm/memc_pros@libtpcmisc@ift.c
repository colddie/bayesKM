/// @file ift.c
/// @author Vesa Oikonen
/// @brief Functions for basic processing of IFT data structure.
///
/*****************************************************************************/
#include "libtpcmisc.h"
/*****************************************************************************/

int IFT_TEST;
/*****************************************************************************/
/** IFT struct status */
static const char *ift_status[] = {
  /*  0 */ "ok",
  /*  1 */ "fault in calling routine",
  /*  2 */ "out of memory",
  /*  3 */ "cannot open file",
  /*  4 */ "cannot write file",
  /*  5 */ "unsupported file type",
  /*  6 */ "key not found",
  /*  7 */ "file contains no data",
  /*  8 */ "value not found",
//  /*  9 */ "key not found",
  0
  };
/*****************************************************************************/

/*****************************************************************************/
/*!
 * Sets IFT status
 *
 * @param ift target IFT structure
 * @param status new status value
 */
void iftSetStatus(IFT *ift, int status) {
  if(ift==NULL) return;
  if(status<0 || status>8) ift->status=ift_status[IFT_FAULT];
  else ift->status=ift_status[status];
}
/*****************************************************************************/

/*****************************************************************************/
/*!
 * Initiate IFT structure. This should be called once before first use.
 *
 * @param ift target ift structure
 */
void iftInit(IFT *ift) {
  if(IFT_TEST) printf("iftInit()\n");
  if(ift==NULL) return;
  memset(ift, 0, sizeof(IFT));
  ift->_memNr=ift->keyNr=0; ift->item=NULL;
  ift->data=NULL; ift->datasize=0;
}
/*****************************************************************************/

/*****************************************************************************/
/*!
 * Free memory allocated for IFT. All contents are destroyed.
 *
 * @param ift target IFT structure
 */
void iftEmpty(IFT *ift) {
  int i;

  if(IFT_TEST) printf("iftEmpty()\n");
  if(ift==NULL) return;
  for(i=0; i<ift->_memNr; i++) {
    if(ift->item[i].key!=NULL) free(ift->item[i].key);
    if(ift->item[i].value!=NULL) free(ift->item[i].value);
  }
  free(ift->item); ift->item=NULL; ift->_memNr=ift->keyNr=0;
  ift->data=NULL; ift->datasize=0;
}
/*****************************************************************************/

/*****************************************************************************/
/** Add specified key and its value to the IFT.
 *  Also comment type (first character pointed to) can be added.
 *  Either key or value can be empty, but not both of them.
\return Returns 0 if ok. Sets ift->status.
 */
int iftPut(
  /** Pointer to initiated IFT; previous contents are not changed */
  IFT *ift,
  /** Key string; can be empty ("" or NULL) */
  char *key,
  /** Value string; can be empty ("" or NULL) */
  char *value,
  /** Pointer to comment character, e.g. '#' or ';' or '!';
   *  can be empty ("" or NULL) */
  char *cmt_type
) {
  int i, add_nr=1000;

  if(ift==NULL) {return(1);}
  if((key==NULL || strlen(key)<1) && (value==NULL || strlen(value)<1)) {
    iftSetStatus(ift, IFT_FAULT); return(2);}
  if(IFT_TEST) {
    printf("iftPut(ift, ");
    if(key!=NULL) printf("\"%s\", ", key); else printf("NULL, ");
    if(value!=NULL) printf("\"%s\", ", value); else printf("NULL, ");
    if(cmt_type!=NULL) printf("\"%1.1s\")\n", cmt_type); else printf("NULL)\n");
  }

  /* If necessary, allocate more memory for items */
  if(ift->_memNr<=ift->keyNr) {
    ift->item=realloc(ift->item, (ift->_memNr+add_nr)*sizeof(IFT_KEY_AND_VALUE));
    if(ift->item==NULL) {iftSetStatus(ift, IFT_NOMEMORY); return(5);}
    ift->_memNr+=add_nr;
    for(i=ift->keyNr; i<ift->_memNr; i++) {
      ift->item[i].type=(char)0;
      ift->item[i].sw=(short int)0;
      ift->item[i].key=NULL;
      ift->item[i].value=NULL;
    }
  }
  
  /* Set the contents */
  /* type */
  if(cmt_type!=NULL) ift->item[ift->keyNr].type=cmt_type[0];
  /* key */
  if(key!=NULL) ift->item[ift->keyNr].key=strdup(key);
  else ift->item[ift->keyNr].key=strdup("");
  if(ift->item[ift->keyNr].key==NULL) {iftSetStatus(ift, IFT_NOMEMORY); return(7);}
  /* value */
  if(value!=NULL) ift->item[ift->keyNr].value=strdup(value);
  else ift->item[ift->keyNr].value=strdup("");
  if(ift->item[ift->keyNr].value==NULL) {iftSetStatus(ift, IFT_NOMEMORY); return(8);}

  ift->keyNr++;
  iftSetStatus(ift, IFT_OK);
  return(0);
}
/*****************************************************************************/

/*****************************************************************************/
/** Add specified key and its floating point (double) value to the IFT.
 *  Also comment type (first character pointed to) can be added.
 *  Key can be empty.
\return Returns 0 if ok. Sets ift->status.
 */
int iftPutDouble(
  /** Pointer to initiated IFT; previous contents are not changed */
  IFT *ift,
  /** Key string; can be empty ("" or NULL) */
  char *key,
  /** Value as double */
  double value,
  /** Pointer to comment character, e.g. '#' or ';' or '!';
   *  can be empty ("" or NULL) */
  char *cmt_type
) {
  char dstr[128];
  sprintf(dstr, "%g", value);
  return(iftPut(ift, key, dstr, cmt_type));
}
/*****************************************************************************/

/*****************************************************************************/
/*!
 * Remove the specified item from IFT.
 *
 * @param ift Pointer to existing IFT
 * @param item Index [0..keyNr-1] of key and value to delete
 * @return 0 if ok.
 */
int iftDeleteItem(IFT *ift, int item) {
  int i;

  if(IFT_TEST) printf("iftDeleteItem(*ift, %d)\n", item);
  if(ift==NULL) return(1);
  iftSetStatus(ift, IFT_FAULT);
  if(ift==NULL) {return(1);}
  if(item<0 || item>=ift->keyNr) {return(2);}

  if(ift->item[item].key!=NULL) free(ift->item[item].key);
  if(ift->item[item].value!=NULL) free(ift->item[item].value);
  ift->item[item].key=ift->item[item].value=NULL;
  for(i=item+1; i<ift->keyNr; i++) {
    ift->item[i-1].type=ift->item[i].type;
    ift->item[i-1].sw=ift->item[i].sw;
    ift->item[i-1].key=ift->item[i].key;
    ift->item[i-1].value=ift->item[i].value;
    ift->item[i].key=ift->item[i].value=NULL;
  }
  ift->keyNr--;
  iftSetStatus(ift, IFT_OK);
  return(0);
}
/*****************************************************************************/

/*****************************************************************************/
/*!
 * Replaces specified value in IFT with a new value.
 * 
 * @param ift Pointer to initiated IFT
 * @param item Index [0..keyNr-1] of key and value
 * @param value Value string; can be empty ("" or NULL)
 * @return 0 if ok.
 */
int iftReplaceNthValue(
  IFT *ift, int item, char *value
) {
  if(ift==NULL) {return(1);}
  if(item>=ift->keyNr) {iftSetStatus(ift, IFT_FAULT); return(2);}
  if(IFT_TEST) printf("iftReplaceNthValue(ift, %d, %s)\n", item, value);
  /* Delete old value */
  if(ift->item[item].value!=NULL) free(ift->item[item].value);
  /* Set new value */
  if(value!=NULL) ift->item[item].value=strdup(value);
  else ift->item[item].value=strdup("");
  if(ift->item[item].value==NULL) {iftSetStatus(ift, IFT_NOMEMORY); return(8);}
  iftSetStatus(ift, IFT_OK);
  return(0);
}
/*****************************************************************************/

/*****************************************************************************/
/** Make a copy (duplicate) of IFT struct.
\return Returns 0 (IFT_OK) when successful, otherwise an appropriate 
    ift status code.
 */
int iftdup(
  /** Pointer to IFT struct to be copied */
  IFT *ift1,
  /** Pointer to initiated IFT struct; any previous contents are deleted. */
  IFT *ift2
) {
  int ret, li;

  if(IFT_TEST) printf("iftdup(*ift1, *ift2)\n");
  /* Check the input */
  if(ift1==NULL || ift2==NULL) return IFT_FAULT;

  /* Empty the new IFT */
  iftEmpty(ift2);

  /* Copy the contents */
  ift2->type=ift1->type;
  ift2->status=ift1->status;
  for(li=0; li<ift1->keyNr; li++) {
    ret=iftPut(ift2, ift1->item[li].key, ift1->item[li].value, NULL);
    if(ret!=IFT_OK) {iftEmpty(ift2); return(ret);} 
    ift2->item[li].type=ift1->item[li].type;
    ift2->item[li].sw=ift1->item[li].sw;
  }
  // keyNr was set by iftPut()
  return IFT_OK;
}
/*****************************************************************************/

/*****************************************************************************/
