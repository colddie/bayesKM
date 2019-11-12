/// @file if.c
/// @author Vesa Oikonen
/// @brief IO functions for IF TAC data.
///
/******************************************************************************/
/** Max line length in IF file */
#define MAX_IF_LINE_LEN 512
/******************************************************************************/
#include "libtpccurveio.h"
#include <unistd.h>
/******************************************************************************/

/******************************************************************************/
/** Write metabolite corrected plasma TAC and blood TAC into IF file format.
    If file exists, a backup file (+BACKUP_EXTENSION) is written.
\return Returns nonzero in case an error is encountered and sets dfterrmsg.
 */
int ifWrite(
  /** Pointer to DFT data that will be written in IF format: first TAC
      must be the metabolite corrected plasma, and the 2nd TAC must be
      the whole blood TAC. */
  DFT *dft,
  /** Name of IF file to be written; also "stdout" is accepted */
  char *filename
) {
  int fi, n;
  char tmp[1024], is_stdout=0;
  FILE *fp;


  /* Check that there is some data to write */
  if(dft->voiNr<2 || dft->frameNr<1 || filename==NULL) {
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
  /* Write data lines */
  for(fi=0; fi<dft->frameNr; fi++) {
    n=fprintf(fp, "%f\t%f\t%f\n",
      dft->x[fi], dft->voi[0].y[fi], dft->voi[1].y[fi]);
    if(n<6) {
      strcpy(dfterrmsg, "cannot write file");
      fclose(fp); return(4);
    }
  }
  /* close file */
  fclose(fp);

  return(0);
}
/******************************************************************************/

/******************************************************************************/
/** Read IF file into DFT data structure, metabolite corrected plasma as TAC #1
    and whole blood as TAC #2.
    Any previous content of DFT is deleted.
\return Returns nonzero in case an error is encountered and sets dfterrmsg.
 */
int ifRead(
  /** Name of IDWC file to be read */
  char *filename,
  /** Pointer to DFT data where to regional TAC data is read */
  DFT *dft
) {
  int fi, n;
  char tmp[MAX_IF_LINE_LEN], *lptr, *cptr;
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
  fi=n=0; do {
    if(fgets(tmp, MAX_IF_LINE_LEN, fp)==NULL) {
      strcpy(dfterrmsg, "wrong format"); fclose(fp); return(3);}
    lptr=tmp; n++;
    /* Read first token, and check for empty lines as well */
    cptr=strtok(lptr, "; \t\n\r"); if(cptr==NULL) continue;
    /* Check for comment line */
    if(cptr[0]=='#' || cptr[0]==';') continue;
    /* Read the sample number */
    fi=atoi(cptr); break;
  } while(1);
  if(fi<1) {strcpy(dfterrmsg, "wrong format"); fclose(fp); return(3);}

  /* Allocate memory for data */
  if(dftSetmem(dft, fi, 2)) {
    strcpy(dfterrmsg, "out of memory"); fclose(fp); return(11);}
  dft->frameNr=fi; dft->voiNr=2;

  /* Read the data */
  n=fi=0;
  do {
    if(fgets(tmp, MAX_IF_LINE_LEN, fp)==NULL) {
      strcpy(dfterrmsg, "wrong format");
      fclose(fp); dftEmpty(dft); return(3);
    }
    lptr=tmp;
    /* Read first token, and check for empty lines as well */
    cptr=strtok(lptr, "; \t\n\r"); if(cptr==NULL) continue;
    /* Check for comment line */
    if(cptr[0]=='#' || cptr[0]==';') continue;
    n++; //printf("%d (fi=%d): %s\n", n, fi, tmp);
    /* read the sample time */
    dft->x[fi]=atof(cptr);
    /* read metabolite corrected plasma */
    cptr=strtok(NULL, "; \t\n\r"); if(cptr==NULL) continue;
    dft->voi[0].y[fi]=atof(cptr);
    /* read whole blood */
    cptr=strtok(NULL, "; \t\n\r"); if(cptr==NULL) continue;
    dft->voi[1].y[fi]=atof(cptr);
    fi++;
    if(fi==dft->frameNr) break;
  } while(1);

  /* close file */
  fclose(fp);

  /* Set DFT "header" */
  dft->_type=1;
  dft->timetype=0;
  dft->timeunit=TUNIT_SEC; /* sec */
  dft->isweight=0;
  strcpy(dft->voi[0].voiname, "Plasma");
  strcpy(dft->voi[1].voiname, "Blood");

  return(0);
}
/******************************************************************************/

/******************************************************************************/
