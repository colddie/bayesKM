/// @file xeleris.c
/// @author Vesa Oikonen
/// @brief I/O functions for GEMS Xeleris TAC files.
///
/*****************************************************************************/
#include "libtpccurveio.h"
/*****************************************************************************/

/*****************************************************************************/
/** Read Xeleris TAC file into DFT data structure.
    Any previous content of DFT is deleted.
\return Returns nonzero in case an error is encountered and sets dfterrmsg.
 */
int xelRead(
  /** Name of Xeleris TAC file to be read */
  char *filename,
  /** Pointer to DFT data where to regional TAC data is read */
  DFT *dft
) {
  int ii, si, ri, fi, n, ret;
  char tmp[1024], tmp2[512], *cptr;
  IFT ift;
  float f[6];
  

  /* Check the arguments */
  if(filename==NULL || dft==NULL || strlen(filename)<1) {
    strcpy(dfterrmsg, "program error"); return(1);
  }
  
  /* Read the file contents */
  iftInit(&ift); ret=iftRead(&ift, filename, 0);
  if(ret) {strcpy(dfterrmsg, ift.status); iftEmpty(&ift); return(2);}
  
  /* Check that this actually is a Xeleris TAC file */
  strcpy(tmp, "Image Position"); ii=iftGet(&ift, tmp);
  if(ii<0) {strcpy(dfterrmsg, "wrong format"); iftEmpty(&ift); return(3);}
  strcpy(tmp, "XAxis"); ii=iftGet(&ift, tmp);
  if(ii<0) {strcpy(dfterrmsg, "wrong format"); iftEmpty(&ift); return(3);}
  strcpy(tmp, "YAxis"); ii=iftGet(&ift, tmp);
  if(ii<0) {strcpy(dfterrmsg, "wrong format"); iftEmpty(&ift); return(3);}
  
  /* Find the data title, and assume that lines after that contain the TACs */
  strcpy(tmp, "Curve     X    Dur    Max    Min   Mean StdDev");
  ii=iftFindNthValue(&ift, tmp, 1);
  if(ii<0) {strcpy(dfterrmsg, "wrong format"); iftEmpty(&ift); return(4);}
  si=ii+1;
  
  /* Determine the number of TACs and time frames */
  ri=fi=0; strcpy(tmp2, ""); n=1;
  for(ii=si; ii<ift.keyNr; ii++) {
    if(sscanf(ift.item[ii].value, "%s", tmp)!=1) {
      strcpy(dfterrmsg, "wrong format"); iftEmpty(&ift); return(5);}
    if(strcmp(tmp, tmp2)==0) {
      n++; if(n>fi) fi=n;
    } else {
      ri++; strcpy(tmp2, tmp); if(n>fi) fi=n; n=1;
    }
    //printf("tmp=%s tmp2=%s n=%d, ri=%d fi=%d\n", tmp, tmp2, n, ri, fi);
  }
  if(ri<1 || fi<1) {strcpy(dfterrmsg, "wrong format"); iftEmpty(&ift); return(6);}
  
  /* Allocate memory for DFT data */
  if(dftSetmem(dft, fi, ri)) {
    strcpy(dfterrmsg, "out of memory"); iftEmpty(&ift); return(7);}
  dft->frameNr=fi; dft->voiNr=ri;
  
  /* Read TAC data and fill DFT */
  ri=fi=0; strcpy(tmp2, ""); n=1;
  for(ii=si; ii<ift.keyNr; ii++) {
    if(sscanf(ift.item[ii].value, "%s %f %f %f %f %f %f", 
       tmp, f, f+1, f+2, f+3, f+4, f+5)!=7) {
      strcpy(dfterrmsg, "wrong format"); iftEmpty(&ift); return(5);}
    if(strcmp(tmp, tmp2)==0) n++;
    else {ri++; strcpy(tmp2, tmp); if(n>fi) fi=n; n=1;}
    //printf("tmp=%s tmp2=%s n=%d, ri=%d fi=%d\n", tmp, tmp2, n, ri, fi);
    /* TAC name */
    strlcpy(dft->voi[ri-1].voiname, tmp, MAX_REGIONSUBNAME_LEN+1);
    //dft->voi[ri-1].voiname[MAX_REGIONSUBNAME_LEN]='\0';
    /* Frame time */
    if(ri==1) {
      dft->x1[n-1]=f[0]; dft->x2[n-1]=f[0]+f[1];
      dft->x[n-1]=0.5*(dft->x1[n-1]+dft->x2[n-1]);
    } else {
      if(fabs(dft->x1[n-1]-f[0])>1.0E-12 || fabs(dft->x2[n-1]-f[0]-f[1])>1.0E-12) {
        strcpy(dfterrmsg, "wrong format"); iftEmpty(&ift); return(8);}
    }
    /* Concentrations */
    dft->voi[ri-1].y[n-1]=f[4];
  }
  
  /* Add image position to the TAC name, if possible */
  strcpy(tmp, "Image Position"); ii=iftGet(&ift, tmp);
  if(ii>=0) {
    f[0]=atof(ift.item[ii].value); sprintf(tmp, "%-6.0f", f[0]); tmp[6]='\0';
    for(ri=0; ri<dft->voiNr; ri++) {
      strcpy(dft->voi[ri].place, tmp);
      sprintf(dft->voi[ri].name, "%s . %s", dft->voi[ri].voiname,
              dft->voi[ri].place);
    }
  }

  /* Determine the concentration (y axis) unit */
  strcpy(tmp, "YAxis"); ii=iftGet(&ift, tmp);
  if(ii>=0 && strncasecmp(ift.item[ii].value, "Uptake (Bqml)", 8)==0) {
    cptr=ift.item[ii].value+8;
    if(     strncasecmp(cptr, "Bqml", 4)==0) strcpy(dft->unit, "Bq/ml");
    else if(strncasecmp(cptr, "kBqml", 5)==0) strcpy(dft->unit, "kBq/ml");
    else if(strncasecmp(cptr, "MBqml", 5)==0) strcpy(dft->unit, "MBq/ml");
    else if(strncasecmp(cptr, "Bqcc", 4)==0) strcpy(dft->unit, "Bq/ml");
    else if(strncasecmp(cptr, "kBqcc", 5)==0) strcpy(dft->unit, "kBq/ml");
    else if(strncasecmp(cptr, "MBqcc", 5)==0) strcpy(dft->unit, "MBq/ml");
  }

  /* Determine the time (x axis) unit */
  dft->timeunit=TUNIT_UNKNOWN;
  strcpy(tmp, "XAxis"); ii=iftGet(&ift, tmp);
  if(ii>=0) {
    if(strcasecmp(ift.item[ii].value, "sec")==0) dft->timeunit=TUNIT_SEC;
    else if(strcasecmp(ift.item[ii].value, "min")==0) dft->timeunit=TUNIT_MIN;
  }
  
  iftEmpty(&ift);

  /* Set the study number based on filename */
  (void)studynr_from_fname(filename, dft->studynr);

  /* Set the rest of DFT "header" */
  dft->_type=1;
  dft->timetype=3;
  dft->isweight=0;

  return(0);
}
/*****************************************************************************/

/*****************************************************************************/
