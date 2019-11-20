/// @file idwc.c
/// @author Vesa Oikonen
/// @brief IO functions for IDWC TAC data.
///
/******************************************************************************/
/** Max line length in IDWC file */
#define MAX_IDWC_LINE_LEN 512
/******************************************************************************/
#include "libtpccurveio.h"
#include <unistd.h>
/******************************************************************************/

/******************************************************************************/
/** Write DFT data into IDWC file format.
    If file exists, a backup file (+BACKUP_EXTENSION) is written.
\return Returns nonzero in case an error is encountered and sets dfterrmsg.
 */
int idwcWrite(
  /** Pointer to DFT data that is written in IDWC format */
  DFT *dft,
  /** Name of IDWC file to be written; also "stdout" is accepted */
  char *filename
) {
  int ri, fi, n;
  char tmp[1024], is_stdout=0;
  FILE *fp;


  /* Check that there is some data to write */
  if(dft->voiNr<1 || dft->frameNr<1 || filename==NULL) {
    strcpy(dfterrmsg, "no data"); return(1);}

  /* Check if writing to stdout */
  if(!strcasecmp(filename, "stdout")) is_stdout=1;

  /* Check if file exists; backup, if necessary */
  if(!is_stdout && access(filename, 0) != -1) {
    strcpy(tmp, filename); strcat(tmp, BACKUP_EXTENSION);
    if(access(tmp, 0) != -1) remove(tmp);
    rename(filename, tmp);
  }

  /* Open output file */
  if(is_stdout) fp=(FILE*)stdout;
  else if((fp = fopen(filename, "w")) == NULL) {
    strcpy(dfterrmsg, "cannot open file"); return(2);}

  /* Write sample number */
  n=fprintf(fp, "%d\n", dft->frameNr);
  if(n<2) {
    strcpy(dfterrmsg, "cannot write file");
    fclose(fp); return(3);
  }
  /* write data lines */
  for(ri=0; ri<dft->voiNr; ri++) {
    for(fi=0; fi<dft->frameNr; fi++) {
      n=fprintf(fp, "%6.1f %18.14f %18.14f %3d\n",
        dft->x[fi], dft->voi[ri].y[fi], dft->w[fi], ri+1);
      if(n<8) {
        strcpy(dfterrmsg, "cannot write file");
        fclose(fp); return(4);
      }
    }
  }
  /* close file */
  fclose(fp);

  return(0);
}
/******************************************************************************/

/******************************************************************************/
/** Read IDWC file into DFT data structure.
    Any previous content of DFT is deleted.
\return Returns nonzero in case an error is encountered and sets dfterrmsg.
 */
int idwcRead(
  /** Name of IDWC file to be read */
  char *filename,
  /** Pointer to DFT data where to regional TAC data is read */
  DFT *dft
) {
  int ri, fi, n;
  char tmp[MAX_IDWC_LINE_LEN], *lptr, *cptr;
  FILE *fp;

  /* Check the arguments */
  if(filename==NULL || dft==NULL || strlen(filename)<1) {
    strcpy(dfterrmsg, "program error"); return(1);
  }

  /* Open file */
  fp=fopen(filename, "r");
  if(fp==NULL) {
    strcpy(dfterrmsg, "cannot open file"); return(2);
  }

  /* Read the line telling the sample number */
  fi=0; do {
    if(fgets(tmp, MAX_IDWC_LINE_LEN, fp)==NULL) {
      strcpy(dfterrmsg, "wrong format"); fclose(fp); return(3);}
    lptr=tmp;
    /* Read first token, and check for empty lines as well */
    cptr=strtok(lptr, "; \t\n\r"); if(cptr==NULL) continue;
    /* Check for comment line */
    if(cptr[0]=='#' || cptr[0]==';') continue;
    /* Read the sample number */
    fi=atoi(cptr); break;
  } while(1);
  if(fi<1) {strcpy(dfterrmsg, "wrong format"); fclose(fp); return(3);}

  /* Read the line number */
  n=0; while(fgets(tmp, MAX_IDWC_LINE_LEN, fp)!=NULL) {
    lptr=tmp;
    /* Read first token, and check for empty lines as well */
    cptr=strtok(lptr, "; \t\n\r"); if(cptr==NULL) continue;
    /* Check for comment line */
    if(cptr[0]=='#' || cptr[0]==';') continue;
    n++;
  }
  /* Calculate the number of TACs */
  ri=n/fi;  //printf("frameNr=%d voiNr=%d\n", fi, ri);
  if(ri<1) {strcpy(dfterrmsg, "wrong format"); fclose(fp); return(4);}


  /* Allocate memory for data */
  if(dftSetmem(dft, fi, ri)) {
    strcpy(dfterrmsg, "out of memory"); fclose(fp); return(11);}
  dft->frameNr=fi; dft->voiNr=ri;

  /* Read the data */
  rewind(fp); n=ri=fi=0;
  do {
    if(fgets(tmp, MAX_IDWC_LINE_LEN, fp)==NULL) {
      strcpy(dfterrmsg, "wrong format");
      fclose(fp); dftEmpty(dft); return(3);
    }
    lptr=tmp;
    /* Read first token, and check for empty lines as well */
    cptr=strtok(lptr, "; \t\n\r"); if(cptr==NULL) continue;
    /* Check for comment line */
    if(cptr[0]=='#' || cptr[0]==';') continue;
    n++; //printf("%d (fi=%d, ri=%d): %s\n", n, fi, ri, tmp);
    /* This time forget the sample number */
    if(n==1) {continue;}
    /* read the sample time */
    dft->x[fi]=atof(cptr);
    /* read sample value */
    cptr=strtok(NULL, "; \t\n\r"); if(cptr==NULL) continue;
    dft->voi[ri].y[fi]=atof(cptr);
    /* read sample weight */
    cptr=strtok(NULL, "; \t\n\r"); if(cptr==NULL) continue;
    dft->w[fi]+=atof(cptr);
    /* read TAC number */
    cptr=strtok(NULL, "; \t\n\r"); if(cptr==NULL) continue;
    if(fi==0) sprintf(dft->voi[ri].voiname, "%-*.*s",
       MAX_REGIONSUBNAME_LEN, MAX_REGIONSUBNAME_LEN, cptr);
    else if(strncasecmp(dft->voi[ri].voiname, cptr, strlen(cptr))!=0) {
      strcpy(dfterrmsg, "wrong format");
      //printf("'%s' '%s'\n", dft->voi[ri].voiname, cptr);
      fclose(fp); dftEmpty(dft); return(4);
    }
    fi++;
    if(fi==dft->frameNr) {fi=0; ri++;}
    if(ri==dft->voiNr) break;
  } while(1);
  dft->voiNr=ri;
  for(fi=0; fi<dft->frameNr; fi++) dft->w[fi]/=(double)dft->voiNr;

  /* close file */
  fclose(fp);

  /* Set DFT "header" */
  dft->_type=1;
  dft->timetype=0; dftFrametimes(dft);
  dft->timeunit=TUNIT_SEC; /* sec */
  dft->isweight=1;

  return(0);
}
/******************************************************************************/

/******************************************************************************/
