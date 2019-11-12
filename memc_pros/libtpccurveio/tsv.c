/// @file tsv.c
/// @author Vesa Oikonen
/// @brief I/O functions for Amide *.tsv TAC files.
///
/*****************************************************************************/
#include "libtpccurveio.h"
/*****************************************************************************/

/*****************************************************************************/
/** Read Amide TAC file (*.tsv) into DFT data structure.
    Any previous content of DFT is deleted.
\return Returns nonzero in case an error is encountered and sets dfterrmsg.
 */
int tsvRead(
  /** Name of Amide TAC file (*.tsv) to be read */
  char *filename,
  /** Pointer to DFT data where to regional TAC data is read */
  DFT *dft
) {
  int ii, si, ri, fi, n, ret;
  char tmp[1024], tmp2[512], *cptr;
  IFT ift;
  float f[13];
  

  /* Check the arguments */
  if(filename==NULL || dft==NULL || strlen(filename)<1) {
    strcpy(dfterrmsg, "program error"); return(1);
  }
  
  /* Read the file contents */
  iftInit(&ift); ret=iftRead(&ift, filename, 0);
  if(ret) {strcpy(dfterrmsg, ift.status); iftEmpty(&ift); return(2);}

  /* Check that this actually is a Amide TAC file */
  strcpy(tmp, "amide"); ii=iftGet(&ift, tmp);
  if(ii<0) {strcpy(dfterrmsg, "wrong format"); iftEmpty(&ift); return(3);}
  if(strncasecmp(ift.item[ii].value, "ROI Analysis File - ", 16)!=0) {
    strcpy(dfterrmsg, "wrong format"); iftEmpty(&ift); return(3);}
  sprintf(dft->comments, "# Amide %s\n# original_filename := %s\n",
    ift.item[ii].value, filename);

  /* Get the number of ROIs */
  strcpy(tmp, "ROI");
  ri=0; while((ii=iftGetNth(&ift, tmp, ri+1))>=0) ri++;
  //printf("ri=%d\n", ri);
  if(ri<1) {strcpy(dfterrmsg, "wrong format"); iftEmpty(&ift); return(4);}
  
  /* Get the max frame number */
  for(ii=0, fi=0; ii<ift.keyNr; ii++) {
    if(ift.item[ii].type=='#') continue;
    n=sscanf(ift.item[ii].value, "%f %f %f %f %f %f %f %f %f %f %f %f %f",
      f, f+1, f+2, f+3, f+4, f+5, f+6, f+7, f+8, f+9, f+10, f+11, f+12);
    if(n<13) continue;
    n=temp_roundf(f[0])+1; // Amide frames start from 0
    if(n>fi) fi=n;
  }
  //printf("fi=%d\n", fi);
  if(fi<1) {strcpy(dfterrmsg, "wrong format"); iftEmpty(&ift); return(5);}
  
  /* Allocate memory for DFT data */
  if(dftSetmem(dft, fi, ri)) {
    strcpy(dfterrmsg, "out of memory"); iftEmpty(&ift); return(7);}
  dft->frameNr=fi; dft->voiNr=ri;
  
  /*
   *  Read one ROI at a time
   */   
  strcpy(tmp, "ROI"); ri=0;
  while((ii=iftGetNth(&ift, tmp, ri+1))>=0) {
    /* Get ROI name */
    strcpy(tmp2, ift.item[ii].value); cptr=strtok(tmp2, " \t\n\r");
    if(cptr!=NULL) {
      strncpy(dft->voi[ri].voiname, cptr, MAX_REGIONSUBNAME_LEN);
      dft->voi[ri].voiname[MAX_REGIONSUBNAME_LEN]='\0';
    } else {
      char buf[128]; snprintf(buf, 128, "%03d", ri+1);
      char *p=buf+strlen(buf)-3;
      snprintf(dft->voi[ri].voiname, MAX_REGIONSUBNAME_LEN, "VOI%s", p);
    }
    /* Get the next data set (output filename will be based on this) */
    for(si=ii+1; si<ift.keyNr; si++) {
      if(ift.item[si].type!='#') continue;
      if(strcasecmp(ift.item[si].key, "Data Set")==0) break;
    }
    if(si<0) {strcpy(dfterrmsg, "wrong format"); iftEmpty(&ift); return(8);}
    strcpy(tmp2, ift.item[si].value); cptr=strtok(tmp2, " \t\n\r");
    if(cptr!=NULL) {
      strncpy(dft->voi[ri].name, cptr, MAX_REGIONNAME_LEN);
      dft->voi[ri].name[MAX_REGIONNAME_LEN]='\0';
    }
    /* Read the frame data until the next data set */
    for(ii=si+1, fi=0; ii<ift.keyNr; ii++)
      if(ift.item[ii].type!='#') break;
    for(; ii<ift.keyNr; ii++) {
      if(ift.item[ii].type=='#') break;
      n=sscanf(ift.item[ii].value, "%f %f %f %f %f %f %f %f %f %f %f %f %f",
        f, f+1, f+2, f+3, f+4, f+5, f+6, f+7, f+8, f+9, f+10, f+11, f+12);
      if(n<13) continue;
      /* Get mean activity concentration */
      dft->voi[ri].y[fi]=f[5];
      /* Get VOI size */
      if(fi==0)
        dft->voi[ri].size=f[10];
      /* Get frame times, if this is the first ROI */
      if(ri==0) {
        dft->x[fi]=f[2];
        dft->x1[fi]=f[2]-0.5*f[1];
        dft->x2[fi]=f[2]+0.5*f[1];
      }
      fi++;
    } // next frame
    ri++;
  }
  iftEmpty(&ift);

  /* Set the study number based on Data Set */
  (void)studynr_from_fname(dft->voi[0].name, dft->studynr);

  /* Set the rest of DFT "header" */
  dft->_type=1;
  dft->isweight=0;
  dftTimeunitToDFT(dft, petTunit(TUNIT_SEC)); // time units are in sec
  dftUnitToDFT(dft, CUNIT_UNKNOWN); // conc units are not known

  return(0);
}
/*****************************************************************************/

/*****************************************************************************/
