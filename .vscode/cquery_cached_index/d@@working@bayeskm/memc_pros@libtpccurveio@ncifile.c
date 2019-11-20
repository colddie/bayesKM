/// @file ncifile.c
/// @author Vesa Oikonen
/// @brief IO for old TPC TAC formats *.roi.nci and *.roi.kbq format.
/// @note Needed onlyjust for compatibility with the old TAC files
///       from Turku PET Centre.
///
/*****************************************************************************/
#include "libtpccurveio.h"
#include <unistd.h>
/*****************************************************************************/
/** Max line length in *.roi.kbq files */
#define ROIKBQ_MAX_LINE_LEN 8192
/*****************************************************************************/
/* Local functions */
/// @cond
void roikbq_strip_spaces(char *s);
void roikbq_move_to_next_line(FILE *fp);
int roikbq_fgets(FILE *fp, int nr, char *s);
/// @endcond
/*****************************************************************************/

/*****************************************************************************/
/** Write DFT contents in *.roi.kbq format.
\return Returns 0 if ok.
*/
int roikbqWrite(
  /** Pointer to DFT */
  DFT *dft, 
  /** Filename */
  char *fname
) {
  int ri, fi;
  char tmp[FILENAME_MAX];
  FILE *fp;

  /* Check that there is some data to write */
  if(dft->voiNr<1 || dft->frameNr<1) {
    strcpy(dfterrmsg, "no data"); return(1);}

  /* Check if file exists; backup, if necessary */
  if(access(fname, 0) != -1) {
    strcpy(tmp, fname); strcat(tmp, BACKUP_EXTENSION);
    if(access(tmp, 0) != -1) remove(tmp);
    rename(fname, tmp);
  }

  /* Open output file */
  if((fp = fopen(fname, "w")) == NULL) {
    strcpy(dfterrmsg, "cannot open file"); fclose(fp); return(2);}

  /* Write title lines */

  /* Write 1st title line (Program name and ROI names) */
  fprintf(fp, "%-15.15s", "cpt2nci 3");
  for(ri=0; ri<dft->voiNr; ri++) {
    fprintf(fp, " %-6.6s", dft->voi[ri].voiname);
    if(strcmp(dft->voi[ri].hemisphere, ".")==0) fprintf(fp, " %-6.6s", " ");
    else fprintf(fp, " %-6.6s", dft->voi[ri].hemisphere);
  }
  fprintf(fp, "\n");

  /* Write 2nd title line (Study nr and Planes) */
  fprintf(fp, "%-15.15s", dft->studynr);
  for(ri=0; ri<dft->voiNr; ri++)
    fprintf(fp, " %-13.13s", dft->voi[ri].place);
  fprintf(fp, "\n");

  /* Write 3rd title line (ROI volumes) */
  fprintf(fp, "%-15.15s", "Time (min)");
  for(ri=0; ri<dft->voiNr; ri++)
    fprintf(fp, " %-13.1f", dft->voi[ri].size);
  fprintf(fp, "\n");

  /* Write times and data of each frame */
  for(fi=0; fi<dft->frameNr; fi++) {
    fprintf(fp, "%8.3f       ", dft->x[fi]);
    for(ri=0; ri<dft->voiNr; ri++)
      if(!isnan(dft->voi[ri].y[fi]))
        fprintf(fp, " %-13.6e", dft->voi[ri].y[fi]);
      else
        fprintf(fp, "       .      ");
    fprintf(fp, "\n");
  }

  /* Close file */
  fclose(fp);

  return(0);
}
/*****************************************************************************/

/*****************************************************************************/
/** Read an old TAC file in *.roi.kbq / *.roi.nci format.
\return Returns 0 if ok.
*/
int roikbqRead(
  /** Pointer to filename */
  char *fname,
  /** Pointer to an empty but initiated DFT structure where data is written */
  DFT *dft
) {
  FILE *fp;
  int i, j, k, l, roi_nr=0, frame_nr=0, ret, nr=0;
  char c, line[ROIKBQ_MAX_LINE_LEN], tmp[25], s1[25], s2[25];
  double a, b;


  /* Check the arguments */
  if(dft==NULL || fname==NULL) {
    strcpy(dfterrmsg, "program fault"); return(1);}
  /* Empty data */
  dftEmpty(dft);

  /* Open file */
  fp=fopen(fname, "r");
  if(fp==NULL) {strcpy(dfterrmsg, "cannot open file"); return(2);}

  /* Check file type */
  rewind(fp); i=0;
  while(fgets(line, ROIKBQ_MAX_LINE_LEN-1, fp) != NULL) {
    roikbq_strip_spaces(line); if(!strlen(line) || line[0]=='#') continue;
    if(!isalpha((int)line[0])) continue;
    if(strncasecmp(line, "cpt2nci 3", 9)==0) {i=1; break;}
  }
  if(i==0) {strcpy(dfterrmsg, "unsupported file format"); fclose(fp); return(3);}

  /* Get data size */
  rewind(fp); roi_nr=frame_nr=0;
  /* read first title line */
  c=fgetc(fp);
  if(c=='#') roikbq_move_to_next_line(fp); else ungetc(c, fp);
  i=roikbq_fgets(fp, 16, tmp);
  if(i<16) {strcpy(dfterrmsg, "unsupported file format"); fclose(fp); return(3);}
  while(i>13) {
    i=roikbq_fgets(fp, 14, tmp); if(i>11 && isalnum((int)tmp[0])) roi_nr++;
  }
  /* read two title lines more */
  roikbq_move_to_next_line(fp); roikbq_move_to_next_line(fp);
  /* now start reading data lines */
  do {
    c=fgetc(fp); if(c==EOF) break; if(c!='\n') frame_nr++;
    roikbq_move_to_next_line(fp);
  } while(c!=EOF);
  if(roi_nr<1 || frame_nr<1) {
    strcpy(dfterrmsg, "unsupported file format"); fclose(fp); return(3);}

  /* Allocate memory for DFT */
  ret=dftSetmem(dft, frame_nr, roi_nr);
  if(ret) {strcpy(dfterrmsg, "out of memory"); fclose(fp); return(4);}
  dft->_type=1;

  /* Read data */
  /* Set file pointer to the beginning of first title line */
  rewind(fp); c=fgetc(fp);
  if(c=='#') roikbq_move_to_next_line(fp); else ungetc(c, fp);
  /* Read the names of ROIs and Hemispheres */
  j=roikbq_fgets(fp, 16, tmp);
  if(j<16) {strcpy(dfterrmsg, "unsupported file format"); fclose(fp); return(3);}
  i=0; while (j>13) {
    j=roikbq_fgets(fp, 14, tmp); if(i>=roi_nr) continue;
    k=sscanf(tmp, "%s %s", s1, s2); if(k<1) continue;
    if(k==1) s2[0]='\0';
    strcpy(dft->voi[dft->voiNr+i].voiname, s1);
    strcpy(dft->voi[dft->voiNr+i].hemisphere, s2);
    if(s2[0]!='\0')
      snprintf(dft->voi[dft->voiNr+i].name, MAX_REGIONNAME_LEN, "%s %s", s1, s2);
    else
      snprintf(dft->voi[dft->voiNr+i].name, MAX_REGIONNAME_LEN, "%s .", s1);
    i++;
  }
  nr=i;
  /* Read the name of place (planes) in 2nd title line */
  c=fgetc(fp); if(c=='#') roikbq_move_to_next_line(fp); else ungetc(c, fp);
  j=roikbq_fgets(fp, 16, tmp);
  if(j<16) {strcpy(dfterrmsg, "unsupported file format"); fclose(fp); return(3);}
  strncpy(dft->studynr, tmp, 6); dft->studynr[6]=(char)0;
  roikbq_strip_spaces(dft->studynr);
  i=0; while(j>13) {
    j=roikbq_fgets(fp, 14, tmp); if(i>=nr) continue;
    k=sscanf(tmp, "%s", s1); if(k<1) s1[0]='\0';
    strcpy(dft->voi[dft->voiNr+i].place, s1);
    if(s1[0]!='\0') {
      strcat(dft->voi[dft->voiNr+i].name, " ");
      strcat(dft->voi[dft->voiNr+i].name, s1);
    } else strcat(dft->voi[dft->voiNr+i].name, " .");
    i++;
  }
  for(j=i; j<nr; j++) strcpy(dft->voi[dft->voiNr+j].place, "");
  /* Read the ROI volumes in 3rd title line */
  c=fgetc(fp); if(c=='#') roikbq_move_to_next_line(fp); else ungetc(c, fp);
  j=roikbq_fgets(fp, 16, tmp);
  if(j<16) {strcpy(dfterrmsg, "unsupported file format"); fclose(fp); return(3);}
  if(strstr(tmp, "Times")) dft->timetype=3;
  else if(strstr(tmp, "Start")) dft->timetype=1;
  else if(strstr(tmp, "End")) dft->timetype=2;
  else dft->timetype=0;
  if(strstr(tmp, "sec")) dft->timeunit=TUNIT_SEC; else dft->timeunit=TUNIT_MIN;
  i=0; while(j>13) {
    j=roikbq_fgets(fp, 14, tmp); if(i>=nr) continue;
    if(sscanf(tmp, "%lf", &dft->voi[dft->voiNr+i].size)<1)
      dft->voi[dft->voiNr+i].size=0.;
    i++;
  }
  for(j=i; j<nr; j++) dft->voi[dft->voiNr+j].size=0.;
  /* Read frame data */
  k=0; c=fgetc(fp); if(c=='#') roikbq_move_to_next_line(fp); else ungetc(c, fp);
  while((j=roikbq_fgets(fp, 16, tmp))>15) {
    if(k>=frame_nr) break;
    roikbq_strip_spaces(tmp); if(strlen(tmp)==0) continue;
    if(dft->timetype==0) {
      dft->x[k]=atof(tmp);
      dft->x1[k]=dft->x2[k]=nan("");
    } else if(dft->timetype==3) {
      sscanf(tmp, "%lf %lf", &a, &b); dft->x1[k]=a; dft->x2[k]=b;
      dft->x[k]=0.5*(dft->x1[k] + dft->x2[k]);
    } else if(dft->timetype==1) {
      dft->x1[k]=dft->x[k]=dft->x2[k]=atof(tmp);
    } else if(dft->timetype==2) {
      dft->x2[k]=dft->x[k]=dft->x1[k]=atof(tmp);
    }
    if(dft->x[k]<-3.e38) dft->x[k]=nan("");
    i=0; while(j>13) {
      j=roikbq_fgets(fp, 14, tmp); if(i>=nr || j<13) continue;
      roikbq_strip_spaces(tmp);
      if(!strlen(tmp) || !strcmp(tmp, ".")) dft->voi[dft->voiNr+i].y[k]=nan("");
      else {
        dft->voi[dft->voiNr+i].y[k]=atof(tmp);
        if(dft->voi[dft->voiNr+i].y[k]<-3.e38)
          dft->voi[dft->voiNr+i].y[k]=nan("");
      }
      i++;
    }
    for(l=i; l<nr; l++) dft->voi[dft->voiNr+l].y[k]=nan("");
    k++;
    c=fgetc(fp); if(c=='#') roikbq_move_to_next_line(fp); else ungetc(c, fp);
  }
  for(i=k; i<frame_nr; i++) {
    dft->x[i]=0.;
    for(j=dft->voiNr; j<dft->voiNr+nr; j++) dft->voi[j].y[i]=nan("");
  }
  /* Set voiNr and frameNr and close the file */
  dft->voiNr=nr; dft->frameNr=frame_nr; fclose(fp);
  /* Set data unit based on filename */
  if(strstr(fname, ".nci")==NULL && strstr(fname, ".NCI")==NULL)
    strcpy(dft->unit, "kBq/ml");
  else strcpy(dft->unit, "nCi/ml");
  /* Set weights in DFT to 1.0 */
  dft->isweight=0;
  for(i=0; i<dft->frameNr; i++) dft->w[i]=1.0;
  /* quit */  
  return(0);
}
/*****************************************************************************/

/*****************************************************************************/
/*
 *  Local functions
 */
/// @cond

/* Removes spaces (and tabs and newlines etc) from string */
void roikbq_strip_spaces(char *s)
{
  int i, len;

  len=strlen(s);
  for(i=len-1; i>=0; i--) if(isspace((int)s[i])) s[i]=(char)0; else break;
  len=i;
  for(i=0; i<len; i++) if(!isspace((int)s[i])) {if(i>0) strcpy(s, &s[i]); break;}
}

/* Moves file pointer to the start of next line, which is not a comment line
   or empty line. */
void roikbq_move_to_next_line(FILE *fp)
{
  char ch;

  do {
    do {ch=fgetc(fp);} while(ch!=EOF && ch!='\n' && ch!='\r');
    if(ch!=EOF) ch=fgetc(fp);
  } while(ch!=EOF && (ch=='\n' || ch=='\r' || ch=='#'));
  ungetc(ch, fp);
}

/* Copy max nr characters from text file fp;
   stops if linefeed or EOF is found, but does not include these,
   as original fgets does.
   Returns the actual number of chars read.
*/
int roikbq_fgets(FILE *fp, int nr, char *s)
{
  int i=0;
  char ch;

  do {
    ch=fgetc(fp);
    if(i<nr && ch!=EOF && ch!='\n' && ch!='\r') {*s++=ch; i++;} else break;
  } while(i<nr && ch!=EOF && ch!='\n' && ch!='\r');
  *s='\0';
  return(i);
}
/// @endcond
/*****************************************************************************/

/*****************************************************************************/
