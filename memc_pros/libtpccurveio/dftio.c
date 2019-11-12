/// @file dftio.c
/// @author Vesa Oikonen
/// @brief Contains I/O functions for DFT files.
/// @bug Line length cannot be larger than max int value, but this is
///       not taken into consideration. This may limit the curve number
///       in some systems.
///
/*****************************************************************************/
#include "libtpccurveio.h"
#include <unistd.h>
/*****************************************************************************/
/** Nr of decimals for concentration values */
int DFT_NR_OF_DECIMALS = 3;
/*****************************************************************************/

/*****************************************************************************/
/** Read TAC file contents into specified DFT data structure. Reads standard DFT files, 
    plain DFT files, and some other formats.
    @sa dftFormat, dftWrite
    @return Returns 0 when successful, in case of error sets dfterrmsg.
 */
int dftRead(
  /** Name of file to be read. */
  char *filename,
  /** Pointer to initiated DFT struct where data will be written; any old content is deleted. */
  DFT *data
) {
  const int verbose=0;
  if(verbose>0) {printf("%s(%s, dft)\n", __func__, filename); fflush(stdout);}

  FILE *fp;
  char *cptr, *line, *lptr, temp[128];
  int ret, i, j, c, type=0, voiNr=0, frameNr=0, longest=0;
  double f;
  IFT ift;


  /* Empty data */
  if(data==NULL) {strcpy(dfterrmsg, "invalid data"); return 1;}
  dftEmpty(data);

  /* Check if file can be opened for read; then close it for now */
  if(verbose>1) {printf("opening file\n"); fflush(stdout);}
  fp=fopen(filename, "r");
  if(fp==NULL) {strcpy(dfterrmsg, "cannot open file"); return 1;}
  fclose(fp);

  /* Try to identify the file format */
  type=dftFormat(filename);
  if(type==DFT_FORMAT_UNKNOWN) {
    strcpy(dfterrmsg, "unknown file format"); return 1;
  } else if(type==DFT_FORMAT_FIT) {
    strcpy(dfterrmsg, "cannot read fit file"); return 1;
  }

  /* Read some special formats */
  if(type==DFT_FORMAT_NCI) { // TPC format before DFT
    if(verbose>1) {printf("calling roikbqRead()\n"); fflush(stdout);}
    ret=roikbqRead(filename, data); if(ret!=0) return ret;
    dftFrametimes(data);
    return(0);
  } else if(type==DFT_FORMAT_IDWC) {
    if(verbose>1) {printf("calling idwcRead()\n"); fflush(stdout);}
    ret=idwcRead(filename, data); if(ret!=0) return ret;
    if(strlen(data->studynr)==0) studynr_from_fname(filename, data->studynr);
    return(0);
  } else if(type==DFT_FORMAT_IF) {
    if(verbose>1) {printf("calling iftRead()\n"); fflush(stdout);}
    ret=ifRead(filename, data); if(ret!=0) return ret;
    if(strlen(data->studynr)==0) studynr_from_fname(filename, data->studynr);
    return(0);
  } else if(type==DFT_FORMAT_CSV_INT || type==DFT_FORMAT_CSV_UK) {
    CSV csv;
    csvInit(&csv); //CSV_TEST=100;
    if(verbose>1) {printf("calling csvRead()\n"); fflush(stdout);}
    ret=csvRead(&csv, filename);
    if(ret==0) ret=csv2dft(&csv, data);
    csvEmpty(&csv);
    if(ret==0) {
      /* Try to read information from an interfile-type header */
      IFT ift; iftInit(&ift);
      if(iftRead(&ift, filename, 1)==0 && ift.keyNr>0) dft_fill_hdr_from_IFT(data, &ift);
      iftEmpty(&ift);
      /* study number, too */
      if(strlen(data->studynr)==0) studynr_from_fname(filename, data->studynr);
      /* Set DFT comments */
      dftSetComments(data);
      return(0);
    }
    // if not ok, then just drop through with plain format
    type=DFT_FORMAT_PLAIN;
#if(0)
    if(ret!=0) {strcpy(dfterrmsg, "cannot read file"); return ret;}
    ret=csv2dft(&csv, data); csvEmpty(&csv);
    if(ret!=0) {strcpy(dfterrmsg, "invalid CSV format"); return ret;}
    if(strlen(data->studynr)==0) studynr_from_fname(filename, data->studynr);
    return(0);
#endif
  }

  /* Try to read supported TAC formats, not the others */
  if(type!=DFT_FORMAT_PLAIN && type!=DFT_FORMAT_STANDARD &&
     type!=DFT_FORMAT_IFT && type!=DFT_FORMAT_PMOD
  ) {
    strcpy(dfterrmsg, "unsupported file format"); return 1;
  }

  /* Try to read information from an interfile-type header */
  if(verbose>1) {printf("calling iftRead()\n"); fflush(stdout);}
  iftInit(&ift);
  if(iftRead(&ift, filename, 1)!=0) iftEmpty(&ift);

  /* Open file */
  if(verbose>1) {printf("re-opening file\n"); fflush(stdout);}
  fp=fopen(filename, "r");
  if(fp==NULL) {
    strcpy(dfterrmsg, "cannot open file"); iftEmpty(&ift); return 1;
  }

  /* Get the length of the longest line */
  i=0; while((c=fgetc(fp))!=EOF) {
    if(c==10 || c==13) {if(i>longest) longest=i; i=0;} else i++;}
  if(i>longest) longest=i;
  rewind(fp); longest+=2;
  if(verbose>1) {printf("  longest := %d\n", longest); fflush(stdout);}

  /* and allocate memory for string of that length */
  line=(char*)malloc((longest+1)*sizeof(char));
  if(line==NULL) {strcpy(dfterrmsg, "out of memory"); iftEmpty(&ift); fclose(fp); return 2;}

  /* Get the frame number */
  i=0; while(fgets(line, longest, fp)!=NULL && *line) {
    if(line[0]=='#') continue;
    cptr=strchr(line, '#'); if(cptr!=NULL) *cptr=(char)0;
    for(j=0; j<(int)strlen(line); j++)
      if(isalnum((int)line[j]) || line[j]=='.') {i++; break;}
  }
  rewind(fp); frameNr=i;
  if(type==DFT_FORMAT_STANDARD) frameNr-=4; /* title lines */
  if(type==DFT_FORMAT_PMOD) frameNr-=1;
  if(verbose>1) {printf("frameNr := %d\n", frameNr); fflush(stdout);}
  if(frameNr<1) {
    strcpy(dfterrmsg, "contains no data");
    iftEmpty(&ift); free(line); fclose(fp); return 1;
  }

  /* Get the number of curves */
  /* find first line that is not an empty or a comment line */
  while(fgets(line, longest, fp)!=NULL) {
    if(line[0]=='#') continue;
    /* If PMOD file, then read it from the title */
    if(type==DFT_FORMAT_PMOD) {voiNr=dftGetPmodTitle(NULL, line); break;}
    /* In plain data file that is first frame, which starts with time */
    /* In normal DFT file that is 1st title line, which starts with 'DFT' */
    char *nline; nline=strdup(line); lptr=nline;
    cptr=strtok(lptr, " \t\n\r"); if(cptr==NULL) {free(nline); continue;}
    voiNr=0; while((cptr=strtok(NULL, " \t\n\r"))!=NULL) voiNr++;
    free(nline);
    break;
  }
  rewind(fp);
  if(verbose>1) {printf("voiNr := %d\n", voiNr); fflush(stdout);}
  if(voiNr<1) {
    strcpy(dfterrmsg, "contains no curves");
    iftEmpty(&ift); free(line); fclose(fp); return 1;
  }

  /* Allocate memory for data */
  if(verbose>1) {printf("allocating memory\n"); fflush(stdout);}
  if(dftSetmem(data, frameNr, voiNr)) {
    strcpy(dfterrmsg, "out of memory");
    iftEmpty(&ift); free(line); fclose(fp); return 2;
  }

  /* For plain data files, set defaults for title data */
  if(type==DFT_FORMAT_PLAIN || type==DFT_FORMAT_IFT) {
    if(verbose>1) {printf("setting title defaults for plain data\n"); fflush(stdout);}
    int u, n; u=voiNr; n=1; while((u/=10)>=1) n++;
    if(n>MAX_REGIONSUBNAME_LEN) n=MAX_REGIONSUBNAME_LEN;
    for(i=0; i<voiNr; i++) {
      snprintf(data->voi[i].voiname, MAX_REGIONSUBNAME_LEN+1, "%0*d", n, i+1);
      strcpy(data->voi[i].name, data->voi[i].voiname);
    }
    /* Default: times in minutes, and frame mid time */
    data->timeunit=TUNIT_UNKNOWN; data->timetype=DFT_TIME_MIDDLE;
  }

  /*
   *  Try to read information from an interfile-type header
   */
  if(ift.keyNr>0) {
    if(verbose>1) {printf("calling dft_fill_hdr_from_IFT\n"); fflush(stdout);}
    dft_fill_hdr_from_IFT(data, &ift);
  }
  iftEmpty(&ift);

  /* Try to read information from a single title line in PMOD files */
  if(type==DFT_FORMAT_PMOD) {
    if(fgets(line, longest, fp)==NULL) {
      strcpy(dfterrmsg, "wrong format"); free(line); fclose(fp); return 101;}
    if(verbose>1) {printf("calling dftGetPmodTitle\n"); fflush(stdout);}
    dftGetPmodTitle(data, line);
    // do not rewind!
  }

  /* 
   *  Read DFT title lines, if they exist
   */
  i=0; 
  if(type==DFT_FORMAT_STANDARD) {
    if(verbose>1) {printf("reading DFT title lines\n"); fflush(stdout);}
    do {
      if(fgets(line, longest, fp)==NULL) {
        strcpy(dfterrmsg, "wrong format"); free(line); fclose(fp); return 102;}
      lptr=line;
      /* Check for comment line */
      if(line[0]=='#') { /* Read comment only if there is space left */
#if(1)
        strlcat(data->comments, line, _DFT_COMMENT_LEN);
#else
        if(strlen(data->comments)+strlen(line)<_DFT_COMMENT_LEN) strcat(data->comments, line);
#endif
        continue;
      }
      cptr=strchr(line, '#'); if(cptr!=NULL) *cptr=(char)0;
      /* Read first token, and check for empty lines as well */
      cptr=strtok(lptr, " \t\n\r"); if(cptr==NULL) continue; else i++;
      if(i==1) { /* VOI names */
        for(j=0; j<voiNr; j++) {
          if((cptr=strtok(NULL, " \t\n\r"))==NULL) {
            strcpy(dfterrmsg, "wrong format"); free(line); fclose(fp); return 103;}
          if(strlen(cptr)>1 || *cptr!='.') {
            strncpy(data->voi[j].name, cptr, MAX_REGIONNAME_LEN);
            data->voi[j].name[MAX_REGIONNAME_LEN]=(char)0;
            strncpy(data->voi[j].voiname, cptr, MAX_REGIONSUBNAME_LEN);
            data->voi[j].voiname[MAX_REGIONSUBNAME_LEN]=(char)0;
            /* if name is long, then divide it into name subfields */
            if(strlen(cptr)>MAX_REGIONSUBNAME_LEN) { 
              strncpy(data->voi[j].hemisphere, cptr+MAX_REGIONSUBNAME_LEN, MAX_REGIONSUBNAME_LEN);
              data->voi[j].hemisphere[MAX_REGIONSUBNAME_LEN]=(char)0;
              if(strlen(cptr)>2*MAX_REGIONSUBNAME_LEN) {
                strncpy(data->voi[j].place, cptr+2*MAX_REGIONSUBNAME_LEN, MAX_REGIONSUBNAME_LEN);
                data->voi[j].place[MAX_REGIONSUBNAME_LEN]=(char)0;
              }
            }
          } else {
            int u, n; u=voiNr; n=1; while((u/=10)>=1) n++; if(n>6) n=6;
            snprintf(data->voi[j].voiname, 7, "%0*d", n, j+1);
            strcpy(data->voi[j].name, data->voi[j].voiname);
          }
        }
      } else if(i==2) { /* 2nd VOI names (hemispheres) */
        if(strcmp(cptr, ".")!=0) {
          strlcpy(data->studynr, cptr, MAX_STUDYNR_LEN);
        } else strcpy(data->studynr, "");
        for(j=0; j<voiNr; j++) {
          if((cptr=strtok(NULL, " \t\n\r"))==NULL) {
            strcpy(dfterrmsg, "missing field on 2nd line");
            free(line); fclose(fp); return 104;
          }
          if(strlen(cptr)>1 || *cptr!='.') {
            strncpy(data->voi[j].hemisphere, cptr, MAX_REGIONSUBNAME_LEN);
            data->voi[j].hemisphere[MAX_REGIONSUBNAME_LEN]=(char)0;
            strcat(data->voi[j].name, " ");
            strlcat(data->voi[j].name, data->voi[j].hemisphere, MAX_REGIONNAME_LEN);
          } else {
            strlcat(data->voi[j].name, " .", MAX_REGIONNAME_LEN+1);
          }
        }
        /* check that there are no more contents */
        if((cptr=strtok(NULL, " \t\n\r"))!=NULL) {
          strcpy(dfterrmsg, "wrong format"); free(line); fclose(fp); return 105;}
      } else if(i==3) { /* unit and VOI place (planes) */
        /* first, copy the unit */
        strlcpy(data->unit, cptr, 13); //data->unit[12]=(char)0;
        /* then, copy the names */
        int len, ii=0;
        j=0;
        while(j<voiNr) {
          if((cptr=strtok(NULL, " \t\n\r"))==NULL) {
            strcpy(dfterrmsg, "missing field on 3rd line");
            free(line); fclose(fp); return 106;
          }
          len=strlen(cptr);
          /* check if cunit consists of two parts */
          if(ii==0 && cptr[0]=='(' && cptr[len-1]==')') {ii++; continue;}
          /* copy the name */
          if(len>1 || *cptr!='.') {
            strlcpy(data->voi[j].place, cptr, MAX_REGIONSUBNAME_LEN+1);
            //data->voi[j].place[MAX_REGIONSUBNAME_LEN]=(char)0;
            strcat(data->voi[j].name, " ");
            strlcat(data->voi[j].name, data->voi[j].place, MAX_REGIONNAME_LEN+1);
          } else {
            strlcat(data->voi[j].name, " .", MAX_REGIONNAME_LEN+1);
          }
          j++; ii++;
        }
      } else if(i==4) { /* time type, time unit, and VOI sizes */
        /* time type */
        if(strcasecmp(cptr, "Time")==0) data->timetype=DFT_TIME_MIDDLE;
        else if(strcasecmp(cptr, "Times")==0) data->timetype=DFT_TIME_STARTEND;
        else if(strcasecmp(cptr, "Start")==0) data->timetype=DFT_TIME_START;
        else if(strcasecmp(cptr, "End")==0) data->timetype=DFT_TIME_END;
        else if(strcasecmp(cptr, "Distance")==0) data->timetype=DFT_TIME_MIDDLE;
        else if(strcasecmp(cptr, "Distances")==3) data->timetype=DFT_TIME_STARTEND;
        else {strcpy(dfterrmsg, "wrong format"); free(line); fclose(fp); return 108;}
        /* time unit */
        if((cptr=strtok(NULL, " \t\n\r"))==NULL) {
          strcpy(dfterrmsg, "wrong format"); free(line); fclose(fp); return 109;}
        strcpy(temp, ""); j=sscanf(cptr, "(%127s)", temp);
        j=strlen(temp)-1; if(j>=0 && temp[j]==')') temp[j]=(char)0;
        data->timeunit=petTunitId(temp);
        if(data->timeunit<0) {
          strcpy(dfterrmsg, "wrong format"); free(line); fclose(fp); return 110;}
        /* volumes */
        for(j=0; j<voiNr; j++) {
          if((cptr=strtok(NULL, " \t\n\r"))==NULL) {
            strcpy(dfterrmsg, "wrong format"); free(line); fclose(fp); return 111;}
          if(strlen(cptr)==1 && *cptr=='.') data->voi[j].size=0.0;
          else data->voi[j].size=atof_dpi(cptr);
        }
      }
    } while(i<4); // Standard DFT title lines are now read
  }

  /* Read data lines */
  if(verbose>1) {printf("reading data lines\n"); fflush(stdout);}
  i=0; while(fgets(line, longest, fp)!=NULL) {

    /* Check for comment line */
    if(line[0]=='#') { /* Read comment only if there is space left */
      strlcat(data->comments, line, _DFT_COMMENT_LEN);
      continue;
    }
    // remove end-of-line comment    
    cptr=strchr(line, '#'); if(cptr!=NULL) *cptr=(char)0;

    /* Read first token, and check for empty lines as well */
    char *nline; nline=strdup(line); lptr=nline;
    cptr=strtok(lptr, " \t\n\r"); if(cptr==NULL) {free(nline); continue;}
    if(i<frameNr) {
      /* Read time(s) */
      if(atof_with_check(cptr, &f)!=0) {
        strcpy(dfterrmsg, "wrong format"); 
        free(line); free(nline); fclose(fp); return 130;
      }
      if(data->timetype==DFT_TIME_STARTEND) {
        data->x1[i]=f; cptr=strtok(NULL, " \t\n\r");
        if(atof_with_check(cptr, &data->x2[i])!=0) {
          strcpy(dfterrmsg, "wrong format"); 
          free(line); free(nline); fclose(fp); return 131;
        }
        data->x[i]=0.5*(data->x1[i]+data->x2[i]);
      } else data->x[i]=f;
      /* curve data */
      for(j=0; j<voiNr; j++) {
        if((cptr=strtok(NULL, " \t\n\r"))==NULL) {
          strcpy(dfterrmsg, "wrong format"); 
          free(line); free(nline); fclose(fp); return 132;
        }
        if(strlen(cptr)==1 && *cptr=='.') 
          data->voi[j].y[i]=nan("");
        else {
          if(atof_with_check(cptr, &data->voi[j].y[i])!=0) {
            strcpy(dfterrmsg, "wrong format"); 
            free(line); free(nline); fclose(fp); return 133;
          }
          c=dec_nr(cptr);
          if(c>DFT_NR_OF_DECIMALS && c<11) DFT_NR_OF_DECIMALS=c;
        }
      }
    }
    free(nline);
    i++;
  }
  if(i!=frameNr) {
    strcpy(dfterrmsg, "wrong format"); fclose(fp); free((char*)line); return 134;}

  /* Close file */
  fclose(fp); free((char*)line);

  /* Set voiNr and frameNr, and type */
  data->voiNr=voiNr; data->frameNr=frameNr;
  if(type==DFT_FORMAT_IFT) {
    data->_type=DFT_FORMAT_PLAIN;
  } else if(type==DFT_FORMAT_PMOD) {
    data->_type=type; // keeping PMOD format
    // unless TAC names could not be read
    for(i=0; i<data->voiNr; i++) if(strlen(data->voi[i].name)<1) {
      data->_type=DFT_FORMAT_PLAIN; break;} 
  } else {
    data->_type=type;
  }

  /* Set study number from filename, if necessary */
  if(strlen(data->studynr)==0) studynr_from_fname(filename, data->studynr);

  /* If weight 'voi' was read, move it to its own place */
  for(i=0; i<data->voiNr; i++) if(!strcasecmp(data->voi[i].voiname, "weight")) {
    data->isweight=1;
    for(j=0; j<data->frameNr; j++) data->w[j]=data->voi[i].y[j];
    /* Move the following VOIs one step backwards */
    for(c=i+1; c<data->voiNr; c++) if(dftCopyvoi(data, c, c-1)) {
      strcpy(dfterrmsg, "cannot read weight"); return 4;}
    data->voiNr--;
    break;
  }
  /* If no weights were found, set uniform weight */
  if(!data->isweight) for(i=0; i<data->frameNr; i++) data->w[i]=1.0;

  /* Calculate frame mid times, or frame start and end times;
     this works only if time units are known */
  dftFrametimes(data);

  if(CSV_TEST>100) dftPrint(data);

  return(0);
}
/*****************************************************************************/

/*****************************************************************************/
/** Determine the type of DFT file. This will replace dftType().
 *  Note that only some of formats are currently identified, and
 *  identification does not mean that dftRead() supports the format.
 *  @return Returns DFT_FORMAT_UNKNOWN or other format id defined in dft.h.
 *  @sa dftRead, dftWrite
 */
int dftFormat(
  /** Pointer to file name; this string is not modified. */
  char *fname
) {
  if(CSV_TEST>0) {printf("dftFormat('%s')\n", fname); fflush(stdout);}

  FILE *fp;
  char tmp[256];
  int c;


  /* Open file */
  fp=fopen(fname, "r");
  if(fp==NULL) return DFT_FORMAT_UNKNOWN;

  /* Binary data? */
  while((c=fgetc(fp))!=EOF) 
    if(!isalnum(c) && !isspace(c) && !isgraph(c) && c!=169) {
      fclose(fp); return DFT_FORMAT_UNKNOWN;
    }
  rewind(fp);

  /* File with one title line without comment mark? */
  while(fgets(tmp, 10, fp) != NULL) {if(strlen(tmp)) break;}
  if(strncasecmp(tmp, "Time[", 5)==0) {fclose(fp); return DFT_FORMAT_PMOD;}
  if(strncasecmp(tmp, "Start[", 6)==0) {fclose(fp); return DFT_FORMAT_PMOD;}
  if(strcasestr(tmp, "Start\tEnd\t")!=NULL) {fclose(fp); return DFT_FORMAT_PMOD;}
  rewind(fp);

  /* Read the first line that is not empty or comment line */
  while(fgets(tmp, 32, fp) != NULL) {
    if(strlen(tmp)==0) continue;
    if(tmp[0]=='#') continue;
    if(strncmp(tmp, "//", 2)==0) continue;
    break;
  }
  rewind(fp);

  /* Check for identification strings */
  if(strncasecmp(tmp, "DFT", 3)==0) {fclose(fp); return DFT_FORMAT_STANDARD;}
  else if(strncasecmp(tmp, "FIT1", 3)==0) {fclose(fp); return DFT_FORMAT_FIT;}
  else if(strncasecmp(tmp, "cpt", 3)==0) {fclose(fp); return DFT_FORMAT_NCI;}

  /* Identify certain filename extensions */
  if(fncasematch(fname, "*.idwc")==1 || fncasematch(fname, "*.idw")==1) {
    fclose(fp); return DFT_FORMAT_IDWC;}
  if(fncasematch(fname, "*.if")==1) {fclose(fp); return DFT_FORMAT_IF;}

  fclose(fp);

  /* Try to read as CSV */
  CSV csv; csvInit(&csv);
  if(csvRead(&csv, fname)==0) {
    int format=DFT_FORMAT_UNKNOWN;
    if(csv.separator==';') format=DFT_FORMAT_CSV_INT;
    else if(csv.separator==',') format=DFT_FORMAT_CSV_UK;
    else if(csv.separator=='\t') {
      int i, commas=0, dots=0;
      for(i=0; i<csv.nr; i++) {
        if(strchr(csv.c[i].content, ',')!=NULL) commas++;
        else if(strchr(csv.c[i].content, '.')!=NULL) dots++;
      }
      if(dots>commas) format=DFT_FORMAT_CSV_UK; else format=DFT_FORMAT_CSV_INT;
    }
    if(CSV_TEST>1) printf("  format=%d\n", format);
    csvEmpty(&csv);
    if(format!=DFT_FORMAT_UNKNOWN) return(format);
  }

#if(0)
  if(fncasematch(fname, "*.csv")==1) {
    /* If ';' is found, then international format, otherwise UK format */
    while((c=fgetc(fp))!=EOF)
      if(c==';') {fclose(fp); return DFT_FORMAT_CSV_INT;}
    fclose(fp); return DFT_FORMAT_CSV_UK;
  }
  if(fncasematch(fname, "*.tsv")==1) {
    /* If ',' is found, then international format, otherwise UK format */
    while((c=fgetc(fp))!=EOF)
      if(c==',') {fclose(fp); return DFT_FORMAT_CSV_INT;}
    fclose(fp); return DFT_FORMAT_CSV_UK;
  }
  fclose(fp);
#endif

  return DFT_FORMAT_PLAIN;
}
/*****************************************************************************/

/*****************************************************************************/
/** Determine the type of DFT file.
   @deprecated
   Replace calls to dftType() by dftFormat(), but note that return values have changed.
  @sa dftFormat
  @return 0=unknown; 1=normal DFT; 2=plain data; 3=fit file; 4=nci file; 5=KI file
 */
int dftType(FILE *fp)
{
  char tmp[256];
  int c;

  /* Find first line that is not empty or comment line */
  rewind(fp);
  while(fgets(tmp, 4, fp) != NULL) {if(strlen(tmp) && tmp[0]!='#') break;}
  rewind(fp);

  /* Check for identification strings */
  if(strncasecmp(tmp, "DFT", 3)==0) return 1;
  else if(strncasecmp(tmp, "FIT1", 3)==0) return 3;
  else if(strncasecmp(tmp, "cpt", 3)==0) return 4;

  /* Binary data? */
  while((c=fgetc(fp))!=EOF) 
    if(!isalnum(c) && !isspace(c) && !isgraph(c) && c!=169) {
      rewind(fp); return 0;
    }

  /* File with one title line without comment mark? */
  rewind(fp);
  while(fgets(tmp, 10, fp) != NULL) {if(strlen(tmp)) break;}
  if(strncasecmp(tmp, "Time", 4)==0) {rewind(fp); return 5;}

  rewind(fp); return 2;
}
/*****************************************************************************/

/*****************************************************************************/
/** Prints to stdout the contents of DFT data structure.
    Mainly for testing purposes.
 */
void dftPrint(DFT *data)
{
  int voi, frame;

  printf("Number of curves: %d     Number of data points: %d\n",
    data->voiNr, data->frameNr);
  printf("Study: '%s'  Unit: '%s'\n", data->studynr, data->unit);
  printf("Time unit and type: %d %d\n", data->timeunit, data->timetype);
  if(strlen(data->radiopharmaceutical))
    printf("Radiopharmaceutical: %s\n", data->radiopharmaceutical);
  if(strlen(data->isotope)) printf("Isotope: %s\n", data->isotope);
  if(strlen(data->scanStartTime))
    printf("Scan start time: %s\n", data->scanStartTime);
  if(strlen(data->injectionTime))
    printf("Injection time: %s\n", data->injectionTime);
  if(data->decayCorrected==DFT_DECAY_CORRECTED)
    printf("Corrected for physical decay: yes\n");
  else if(data->decayCorrected==DFT_DECAY_NOTCORRECTED)
    printf("Corrected for physical decay: no\n");
  printf("_datasize = %d\n", data->_dataSize);
  for(voi=0; voi<data->voiNr; voi++) {
    if(strlen(data->voi[voi].name)>0)
      printf("\nROI name: '%s' Size: %g\n",
        data->voi[voi].name, data->voi[voi].size);
    else
      printf("\nROI name: '%s' '%s' '%s'  Size: %g\n",
        data->voi[voi].voiname, data->voi[voi].hemisphere, data->voi[voi].place,
        data->voi[voi].size);
    for(frame=0; frame<data->frameNr; frame++) {
      printf("%03d:  %11.3e %11.3e %11.3e    %11.3e %11.3e %11.3e\n",
        frame+1, data->x[frame], data->x1[frame], data->x2[frame],
        data->voi[voi].y[frame], data->voi[voi].y2[frame],
        data->voi[voi].y3[frame]);
    }
  }
  printf("Comments:\n");
  if(strlen(data->comments)>0) printf("%s\n", data->comments);
  printf("Weights:\n");
  if(data->isweight)
    for(frame=0; frame<data->frameNr; frame++)
      printf(" %03d  %11.3e %11.3e  %11.3e\n", frame+1,
        data->x1[frame], data->x2[frame], data->w[frame]);
  else
    printf(" contains no weights.\n");
  return;
}
/*****************************************************************************/

/*****************************************************************************/
/** Write DFT data, usually containing regional time-activity curves (TACs) into specified file.

    The file format specified in data is applied.
    Number of decimals can be determined by changing global variable DFT_NR_OF_DECIMALS.
    @sa dftRead
    @return Returns 0 when successful, in case of error sets dfterrmsg.
 */
int dftWrite(
  /** Pointer to DFT struct which contents are to be written. */
  DFT *data,
  /** File name where DFT contents are written. If file exists, original file
      is renamed to a backup file. */
  char *filename
) {
  int i, j, n, sw, mfw, prec;
  char tmp[1024], tmp2[128], is_stdout=0;
  FILE *fp;


  /* Check that there is some data to write */
  if(data->voiNr<1 || data->frameNr<1) {
    strcpy(dfterrmsg, "no data"); return 1;}

  /* If format if DFT_FORMAT_HTML or extension is *.HTM(L),
     write in HTML format */
  if(data->_type==DFT_FORMAT_HTML || fncasematch(filename, "*.htm")==1 ||
     fncasematch(filename, "*.html")==1)
  {
    return(dftWriteHTML(data, filename, 1));
  }

  /* Check if writing to stdout */
  if(!strcasecmp(filename, "stdout")) is_stdout=1;

  /* Set minimum field width and precision for concentrations */
  mfw=11; if(DFT_NR_OF_DECIMALS>3) mfw+=DFT_NR_OF_DECIMALS-3;
  prec=0; if(DFT_NR_OF_DECIMALS>prec) prec=DFT_NR_OF_DECIMALS;

  /* Check if file exists; backup, if necessary */
  if(!is_stdout) (void)backupExistingFile(filename, NULL, NULL);

  /* Open output file */
  if(is_stdout) fp=(FILE*)stdout;
  else if((fp = fopen(filename, "w")) == NULL) {strcpy(dfterrmsg, "cannot open file"); return 2;}

  /* Write title lines */
  if(data->_type==DFT_FORMAT_STANDARD) {
    /* 1st line with filetype identification string and region names */
    n=fprintf(fp, "%s", DFT_VER);
    if(n==0) {strcpy(dfterrmsg, "disk full"); if(!is_stdout) fclose(fp); return 3;}
    for(i=0; i<data->voiNr; i++)
      fprintf(fp,"\t%s", data->voi[i].voiname);
    if(data->isweight) fprintf(fp,"\t%s", "weight");
    fprintf(fp, "\n");
    /* 2nd line with study identification and 2nd name (hemisphere) */
    if(strlen(data->studynr)) 
      strlcpy(tmp, data->studynr, 11);
    else 
      strcpy(tmp, ".");
    fprintf(fp, "%s", tmp);
    for(i=0; i<data->voiNr; i++) {
      if(strlen(data->voi[i].hemisphere)) strcpy(tmp, data->voi[i].hemisphere);
      else strcpy(tmp, ".");
      fprintf(fp, "\t%s", tmp);
    }
    if(data->isweight) fprintf(fp,"\t%s", ".");
    fprintf(fp, "\n");
    /* 3rd line with unit and plane names */
    if(strlen(data->unit)) strcpy(tmp, data->unit); else strcpy(tmp, ".");
    fprintf(fp, "%s", tmp);
    for(i=0; i<data->voiNr; i++) {
      if(strlen(data->voi[i].place)) strcpy(tmp, data->voi[i].place);
      else strcpy(tmp, ".");
      fprintf(fp, "\t%s", tmp);
    }
    if(data->isweight) fprintf(fp,"\t%s", ".");
    fprintf(fp, "\n");
    /* 4th line with time type & unit and region volumes */
    switch(data->timetype) {
      case DFT_TIME_MIDDLE:
        if(data->timeunit==TUNIT_MM || data->timeunit==TUNIT_UM || data->timeunit==TUNIT_CM)
          strcpy(tmp, "Distance ");
        else
          strcpy(tmp, "Time ");
        break;
      case DFT_TIME_START: strcpy(tmp, "Start "); break;
      case DFT_TIME_END: strcpy(tmp, "End "); break;
      case DFT_TIME_STARTEND:
        if(data->timeunit==TUNIT_MM || data->timeunit==TUNIT_UM || data->timeunit==TUNIT_CM)
          strcpy(tmp, "Distances ");
        else
          strcpy(tmp, "Times ");
        break;
    }
    snprintf(tmp2, 128, "(%s)", petTunit(data->timeunit)); strcat(tmp, tmp2);
    fprintf(fp, "%s", tmp);
    for(i=0; i<data->voiNr; i++) {
      if(data->voi[i].size>=0.0) sprintf(tmp, "%.*e", prec, data->voi[i].size);
      else strcpy(tmp, ".");
      fprintf(fp, "\t%s", tmp);
    }
    if(data->isweight) fprintf(fp,"\t%s", ".");
    fprintf(fp, "\n");
  } else if(data->_type==DFT_FORMAT_PMOD) { // write PMOD title line
    char *cptr, tunit[128], cunit[128]; // Set units to format accepted by PMOD
    if(data->timeunit==TUNIT_SEC) strcpy(tunit, "seconds");
    else if(data->timeunit==TUNIT_MIN) strcpy(tunit, "minutes");
    else strcpy(tunit, petTunit(data->timeunit));
    strcpy(cunit, data->unit); cptr=strstr(cunit, "mL");
    if(cptr!=NULL) {*cptr='c'; cptr++; *cptr='c';}
    if(data->timetype==DFT_TIME_STARTEND)
      fprintf(fp, "start[%s]\tend[%s]", tunit, cunit);
    else
      fprintf(fp, "time[%s]", tunit);
    for(i=0; i<data->voiNr; i++) {
      /* write roi names, must not include spaces */
      if(strchr(data->voi[i].name, ' ')==NULL) {
        fprintf(fp, "\t%s", data->voi[i].name);
      } else {
        fprintf(fp, "\t%s", data->voi[i].voiname);
        if(strlen(data->voi[i].hemisphere)>0)
          fprintf(fp, "-%s", data->voi[i].hemisphere);
        if(strlen(data->voi[i].place)>0)
          fprintf(fp, "-%s", data->voi[i].place);
      }
      /* write calibration unit when necessary */
      if(i==0 && data->timetype!=DFT_TIME_STARTEND)
        fprintf(fp, "[%s]", cunit);
    }
    if(data->isweight) fprintf(fp,"\t%s", "weight");
    fprintf(fp, "\n");
  }

  /* Make sure that only frame mid times are saved without titles */
  if(data->_type==DFT_FORMAT_PLAIN && data->timetype!=DFT_TIME_MIDDLE) sw=0;
  else sw=data->timetype;

  /* Write data */
  for(j=0; j<data->frameNr; j++) {
    /* Time(s) */
    switch(sw) {
      case 0: fprintf(fp, "%.5f", data->x[j]); break;
      case 1: fprintf(fp, "%.5f", data->x1[j]); break;
      case 2: fprintf(fp, "%.5f", data->x2[j]); break;
      case 3: fprintf(fp, "%.5f\t%.5f", data->x1[j], data->x2[j]); break;
    }
    /* y values */
    for(i=0; i<data->voiNr; i++) {
      if(!isnan(data->voi[i].y[j])) sprintf(tmp, "%.*e", prec, data->voi[i].y[j]);
      else strcpy(tmp, ".");
      fprintf(fp, "\t%s", tmp);
    }
    if(data->isweight) {
      if(data->_type==DFT_FORMAT_STANDARD || data->_type==DFT_FORMAT_PMOD) {
        if(!isnan(data->w[j])) sprintf(tmp, "%.*e", prec, data->w[j]);
        else strcpy(tmp, "."); 
        fprintf(fp,"\t%s", tmp);
      }
    }
    n=fprintf(fp, "\n");
    if(n==0) {strcpy(dfterrmsg, "disk full"); if(!is_stdout) fclose(fp); return 3;}
  }

  /* Write comments */
  if(1) { //data->_type!=DFT_FORMAT_PMOD) {
    n=strnlen(data->comments, _DFT_COMMENT_LEN);
    if(n) for(i=0; i<n; i++) {
      if(i>0 && data->comments[i]=='#' && (data->comments[i-1]!='\n' && data->comments[i-1]!='\r'))
        fputc('\n', fp);
      fputc(data->comments[i], fp);
    }
  }

  /* Close file */
  if(!is_stdout) {fflush(fp); fclose(fp);}

  return 0;
}
/*****************************************************************************/

/*****************************************************************************/
/** Write DFT contents in HTML table format
    If file exists, a backup file (+BACKUP_EXTENSION) is written also.
    If "stdout" is given as filename, output is directed to stdout
    In case of an error, description is written in dfterrmsg.
\return Returns 0 if ok.
*/
int dftWriteHTML(
  /** Input DFT */
  DFT *dft,
  /** HTML filename */
  char *fname,
  /** Table orientation: 1=original, 2=transposed */
  int orientation
) {
  int ri, fi, is_stdout=0, ret;
  char tmp[FILENAME_MAX];
  FILE *fp;


  /* Check input */
  strcpy(dfterrmsg, "invalid input to dftWriteHTML()");
  if(dft==NULL || dft->frameNr<1 || dft->voiNr<1) return(1);
  if(fname==NULL || strlen(fname)<1) return(1);
  strcpy(dfterrmsg, "");
  /* Check if writing to stdout */
  if(!strcasecmp(fname, "stdout")) is_stdout=1;

  /* Check if file exists; backup, if necessary */
  if(!is_stdout && access(fname, 0) != -1) {
    strcpy(tmp, fname); strcat(tmp, BACKUP_EXTENSION); rename(fname, tmp);}
  strcpy(dfterrmsg, "cannot write file");

  /* Open output file */
  if(is_stdout) fp=(FILE*)stdout;
  else if((fp=fopen(fname, "w"))==NULL) {strcpy(dfterrmsg, "cannot open file"); return(2);}

  /* Write XHTML 1.1 doctype and head */
  if(dftWriteXHTML11_doctype(fp) || dftWriteXHTML11_head(fp, "")) {
    strcpy(dfterrmsg, "disk full"); if(!is_stdout) fclose(fp); return(3);
  }

  /* Start writing the body of the HTML file */
  ret=fprintf(fp, "<body>\n");
  if(ret==0) {strcpy(dfterrmsg, "disk full"); if(!is_stdout) fclose(fp); return(3);}

  /* Start the div for table, and the table */
  fprintf(fp, "\n<div id=\"table\">\n");
  fprintf(fp, "<table>\n");

  /* Write the title lines, if there are titles, into the head of the table */
  if(dft->_type!=DFT_FORMAT_PLAIN) {
    fprintf(fp, "<thead>\n");
    if(orientation!=2) {
      /* Empty cell and Region names */
      fprintf(fp, "<tr><td> </td>\n");
      for(ri=0; ri<dft->voiNr; ri++)
        fprintf(fp, "<th>%s</th>\n", dft->voi[ri].voiname);
      fprintf(fp, "</tr>\n");
      /* Study number and Hemispheres */
      fprintf(fp, "<tr><th>%s</th>\n", dft->studynr);
      for(ri=0; ri<dft->voiNr; ri++)
        fprintf(fp, "<th>%s</th>\n", dft->voi[ri].hemisphere);
      fprintf(fp, "</tr>\n");
      /* Unit and Places */
      fprintf(fp, "<tr><th>%s</th>\n", dft->unit);
      for(ri=0; ri<dft->voiNr; ri++)
        fprintf(fp, "<th>%s</th>\n", dft->voi[ri].place);
      fprintf(fp, "</tr>\n");
      /* Time unit and Volumes */
      fprintf(fp, "<tr><th>%s</th>\n", petTunit(dft->timeunit));
      for(ri=0; ri<dft->voiNr; ri++)
        fprintf(fp, "<th>%g</th>\n", dft->voi[ri].size);
      fprintf(fp, "</tr>\n");
    } else if(orientation==2) {
      /* Study number and Unit */
      fprintf(fp, "<tr><th>%s</th>\n", dft->studynr);
      fprintf(fp, "<th>%s</th></tr>\n", dft->unit);  
    }
    /* End the head of table */
    fprintf(fp, "</thead>\n");
  }

  /* Write the DFT data into the body of the table */
  fprintf(fp, "<tbody>\n");

  /* If transposed, write titles */
  if(orientation==2) {
    fprintf(fp, "<tr><th>Region</th><th>Hemisphere</th><th>Plane</th>\n");
    fprintf(fp, "<th>Volume</th></tr>\n");
  }
  if(orientation!=2) { /* Not transposed */
    for(fi=0; fi<dft->frameNr; fi++) {
      if(fi%2) strcpy(tmp, "evenframe"); else strcpy(tmp, "oddframe");
      /* Time */
      fprintf(fp, "<tr class=\"%s\"><th>%g</th>\n", tmp, dft->x[fi]);
      /* Values */
      for(ri=0; ri<dft->voiNr; ri++)
        fprintf(fp, "<th>%g</th>", dft->voi[ri].y[fi]);
      fprintf(fp, "</tr>\n");
    }
  } else { /* transposed */
    for(ri=0; ri<dft->voiNr; ri++) {
      if(ri%2) strcpy(tmp, "evenframe"); else strcpy(tmp, "oddframe");
      /* Region names and volume */
      fprintf(fp, "<tr class=\"%s\"><th>%s</th><th>%s</th><th>%s</th>\n", tmp,
        dft->voi[ri].voiname, dft->voi[ri].hemisphere, dft->voi[ri].place);
      fprintf(fp, "<td>%g</td>\n", dft->voi[ri].size);
      /* Values */
      for(fi=0; fi<dft->frameNr; fi++) {
        fprintf(fp, "<td>%g</td>", dft->voi[ri].y[fi]);
      }
      fprintf(fp, "</tr>\n");
    }
  }
  /* End the table body and table */
  fprintf(fp, "</tbody></table>\n");

  /* Stop writing the body of the HTML file, and end the file */
  fprintf(fp, "</div>\n");
  ret=fprintf(fp, "</body></html>\n\n");
  if(ret==0) {
    strcpy(dfterrmsg, "disk full");
    if(!is_stdout) fclose(fp);
    return(3);
  }

  /* Close file */
  if(!is_stdout) {fflush(fp); fclose(fp);}
  strcpy(dfterrmsg, "");

  return(0);
}
/*****************************************************************************/

/*****************************************************************************/
/** Write XHTML 1.1 doctype into an opened file pointer.
\return Returns 0 if successful, non-zero in case of a failure.
 */
int dftWriteXHTML11_doctype(
  FILE *fp
) {
  int n;
  if(fp==NULL) return(1);
  n=fprintf(fp, "<!DOCTYPE html PUBLIC \"-//W3C//DTD XHTML 1.1//EN\" ");
  n+=fprintf(fp, "\"http://www.w3.org/TR/xhtml11/DTD/xhtml11.dtd\">\n");
  if(n<20) return(2);
  n=fprintf(fp, "<html xmlns=\"http://www.w3.org/1999/xhtml\" xml:lang=\"en\">\n\n");
  if(n<20) return(2);
  return(0);
}
/*****************************************************************************/

/*****************************************************************************/
/** Write XHTML 1.1 head for DFT file into an opened file pointer.
\return Returns 0 if successful, non-zero in case of a failure.
 */
int dftWriteXHTML11_head(
  /** File pointer where to write */
  FILE *fp,
  /** Author name, for example software name */
  char *author_name
) {
  int n;
  if(fp==NULL) return(1);

  n=fprintf(fp, "<head>\n"); if(n<6) return(2);
  fprintf(fp, "  <title>PET data</title>\n");
  fprintf(fp, "  <meta http-equiv=\"content-type\" content=\"text/html; ");
  fprintf(fp, "charset=iso-8859-1\" />\n");
  fprintf(fp, "  <meta http-equiv=\"content-language\" content=\"en-gb\" />\n");
  fprintf(fp, "  <meta name=\"author\" content=\"%s\" />\n", author_name);
  fprintf(fp, "  <meta name=\"ProgId\" content=\"Excel.Sheet\" />\n");
  fprintf(fp, "  <link rel=\"icon\" href=\"http://www.turkupetcentre.net/favi");
  fprintf(fp, "con.ico\" type=\"image/x-icon\" />\n");
  fprintf(fp, "  <link rel=\"shortcut icon\" href=\"http://www.turkupetcentre");
  fprintf(fp, ".net/favicon.ico\" type=\"image/x-icon\" />\n");
  /* write local CSS with basic settings in case that external CSS is not available */
  fprintf(fp, "  <style type=\"text/css\">\n");
  fprintf(fp, "    thead {background-color:#999999; color:black;}\n");
  fprintf(fp, "    table {text-align:left; width:100%%; border-collapse:coll");
  fprintf(fp, "apse; empty-cells:show;}\n");
  fprintf(fp, "    td {border:1px solid black;}\n");
  fprintf(fp, "    .oddframe {background-color:#FFFFFF; color:black;}\n");
  fprintf(fp, "    .evenframe {background-color:#CCCCCC; color:black;}\n");
  fprintf(fp, "    #regcontainer ul {margin-left:0; padding-left:0;}\n");
  fprintf(fp, "    #regcontainer ul li {display:inline; list-style-type:none;}\n");
  fprintf(fp, "    #regcontainer a {padding:2px 4px;}\n");
  fprintf(fp, "    <!--table\n");
  fprintf(fp, "    	{mso-displayed-decimal-separator:\"\\.\";\n");
  fprintf(fp, "    	mso-displayed-thousand-separator:\" \";}\n");
  fprintf(fp, "    -->\n");
  fprintf(fp, "  </style>\n");
  /* load external CSS with more fancy settings */
  fprintf(fp, "  <link rel=\"stylesheet\" type=\"text/css\" ");
  fprintf(fp, "href=\"http://www.turkupetcentre.net/dft.css\" />\n");
  fprintf(fp, "</head>\n");
  if(n<7) return(2);
  return(0);
}
/*****************************************************************************/

/*****************************************************************************/
/** Read certain keys from IFT and set DFT fields accordingly.
\return Returns the nr of identified keys.
 */
int dft_fill_hdr_from_IFT(
  /** Pointer to allocated DFT struct where information will be written */
  DFT *dft,
  /** Pointer to IFT struct from where information is retrieved */
  IFT *ift
) {
  int ki, ri, ok_nr=0, ret;
  char keystr[256], *cptr;

  //printf("dft_fill_hdr_from_IFT()\n");

  /* Check for study number */
  strcpy(keystr, "studynr"); ki=iftGet(ift, keystr);
  if(ki==-1) {strcpy(keystr, "study number"); ki=iftGet(ift, keystr);}
  if(ki==-1) {strcpy(keystr, "study_number"); ki=iftGet(ift, keystr);}
  if(ki>=0) {
    strlcpy(dft->studynr, ift->item[ki].value, MAX_STUDYNR_LEN);
    ok_nr++;
  }

  /* Check for time unit */
  strcpy(keystr, "timeunit"); ki=iftGet(ift, keystr);
  if(ki==-1) {strcpy(keystr, "time unit"); ki=iftGet(ift, keystr);}
  if(ki==-1) {strcpy(keystr, "time_unit"); ki=iftGet(ift, keystr);}
  if(ki==-1) {strcpy(keystr, "Time units"); ki=iftGet(ift, keystr);}
  if(ki>=0) {
    ret=petTunitId(ift->item[ki].value);
    if(ret>=0 && ret!=TUNIT_UNKNOWN) {dft->timeunit=ret; ok_nr++;}
  }

  /* Check for sample unit */
  strcpy(keystr, "unit"); ki=iftGet(ift, keystr);
  if(ki==-1) {strcpy(keystr, "Activity units"); ki=iftGet(ift, keystr);}
  if(ki>=0) {
    strncpy(dft->unit, ift->item[ki].value, 12); dft->unit[12]=(char)0;
    ok_nr++;
  }

  /* Check for region names */
  ri=0;
  do {
    strcpy(keystr, "voiname"); ki=iftGetNth(ift, keystr, ri+1);
    if(ki>=0) {
      strncpy(dft->voi[ri].name, ift->item[ki].value, MAX_REGIONNAME_LEN);
      dft->voi[ri].name[MAX_REGIONNAME_LEN]=(char)0;
      rnameSplit(ift->item[ki].value, dft->voi[ri].voiname,
        dft->voi[ri].hemisphere, dft->voi[ri].place, MAX_REGIONSUBNAME_LEN);
    }
    ri++; ok_nr++;
  } while(ki>=0 && ri<dft->_voidataNr);

  /* Check for region volumes */
  strcpy(keystr, "sizes"); ki=iftGet(ift, keystr);
  if(ki==-1) {strcpy(keystr, "volumes"); ki=iftGet(ift, keystr);}
  if(ki>=0 && strlen(ift->item[ki].value)) {
    ri=0; cptr=strtok(ift->item[ki].value, " ;\t");
    while(cptr!=NULL && ri<dft->_voidataNr) {
      if(strcmp(cptr, ".")==0) {ri++; continue;}
      dft->voi[ri++].size=atof_dpi(cptr);
      cptr=strtok(NULL, " ;\t");
    }
    ok_nr++;
  }

  /* Check for the name of radiopharmaceutical */
  strcpy(keystr, "radiopharmaceutical"); ki=iftGet(ift, keystr);
  if(ki>=0 && strlen(ift->item[ki].value)) {
    strncpy(dft->radiopharmaceutical, ift->item[ki].value, 31);
    dft->radiopharmaceutical[31]=(char)0;
    ok_nr++;
  }

  /* Check for isotope name */
  strcpy(keystr, "isotope"); ki=iftGet(ift, keystr);
  if(ki>=0 && strlen(ift->item[ki].value)) {
    strncpy(dft->isotope, ift->item[ki].value, 6); dft->isotope[6]=(char)0;
    ok_nr++;
  }

  /* Check if there is anything about correction for physical decay */
  strcpy(keystr, "decay_correction"); ki=iftGet(ift, keystr);
  if(ki==-1) {strcpy(keystr, "decay correction"); ki=iftGet(ift, keystr);}
  if(ki>=0 && strlen(ift->item[ki].value)) {
    if(strncasecmp(ift->item[ki].value, "Yes", 1)==0) {
      dft->decayCorrected=DFT_DECAY_CORRECTED; ok_nr++;
    } else if(strncasecmp(ift->item[ki].value, "No", 1)==0) {
      dft->decayCorrected=DFT_DECAY_NOTCORRECTED; ok_nr++;
    }
  }

  /* Check for injection time */
  strcpy(keystr, "injection time"); ki=iftGet(ift, keystr);
  if(ki==-1) {strcpy(keystr, "injection_time"); ki=iftGet(ift, keystr);}
  if(ki>=0 && strlen(ift->item[ki].value)) {
    /* check if time is in correct format; otherwise try to correct it */
    if(ift->item[ki].value[4]=='-' && ift->item[ki].value[7]=='-' &&
       ift->item[ki].value[13]==':' && ift->item[ki].value[16]==':') {
      strncpy(dft->injectionTime, ift->item[ki].value, 19);
      dft->injectionTime[19]=(char)0; ok_nr++;
    } else if(ift->item[ki].value[2]=='.' && ift->item[ki].value[5]=='.' &&
              ift->item[ki].value[13]==':' && ift->item[ki].value[16]==':') {
      sprintf(dft->injectionTime, "%4.4s-%2.2s-%2.2s %8.8s",
        ift->item[ki].value+6, ift->item[ki].value+3, ift->item[ki].value,
        ift->item[ki].value+11);
      ok_nr++;
    }
  }

  /* Check for scan start time */
  strcpy(keystr, "scan start time"); ki=iftGet(ift, keystr);
  if(ki==-1) {strcpy(keystr, "scan_start_time"); ki=iftGet(ift, keystr);}
  if(ki==-1) {strcpy(keystr, "scan_start"); ki=iftGet(ift, keystr);}
  if(ki==-1) {strcpy(keystr, "scan start"); ki=iftGet(ift, keystr);}
  if(ki>=0 && strlen(ift->item[ki].value)) {
    /* check if time is in correct format; otherwise try to correct it */
    if(ift->item[ki].value[4]=='-' && ift->item[ki].value[7]=='-' &&
       ift->item[ki].value[13]==':' && ift->item[ki].value[16]==':') {
      strncpy(dft->scanStartTime, ift->item[ki].value, 19);
      dft->scanStartTime[19]=(char)0; ok_nr++;
    } else if(ift->item[ki].value[2]=='.' && ift->item[ki].value[5]=='.' &&
              ift->item[ki].value[13]==':' && ift->item[ki].value[16]==':') {
      sprintf(dft->scanStartTime, "%4.4s-%2.2s-%2.2s %8.8s",
        ift->item[ki].value+6, ift->item[ki].value+3, ift->item[ki].value,
        ift->item[ki].value+11);
      ok_nr++;
    }
  }

  return(ok_nr);
}
/*****************************************************************************/

/*****************************************************************************/
/** Read single title line from PMOD files and set DFT fields accordingly.
 *  Alternatively, reads the number of regions in PMOD title line.
\return Returns 0 if title contents were successfully saved in DFT struct;
    if pointer to DFT struct is not specified, then the number of regions
    is returned.
 */
int dftGetPmodTitle(
  /** Pointer to allocated DFT struct where information will be written;
   *  Enter NULL, if only the nr of regions is to be returned. */
  DFT *dft,
  /** Pointer to string containing the title line; string is not modified */
  char *title_line
) {
  //printf("dftGetPmodTitle()\n");
  //printf(" '%s'\n", title_line);
  //if(dft==NULL) printf(" dft=NULL\n");

  int ti=0, ri, ret, tabs=0, timetype=DFT_TIME_MIDDLE;
  char *cptr, *cptr2, rname[MAX_REGIONNAME_LEN+1];
  char unit[MAX_UNITS_LEN+1], separs[8];

  /* Check the input */
  if(strlen(title_line)<1) return 1;

  /* Count the nr of tabs in title line */
  cptr=title_line; while(*cptr) {if(*cptr=='\t') tabs++; cptr++;}
  /* If tabs found, then do not use space as field separator */
  if(tabs>0) strcpy(separs, "\t\n\r"); else strcpy(separs, " \t\n\r");

  /* Make a copy of title line */
  char ntl[1+strlen(title_line)], *tl;
  strcpy(ntl, title_line); tl=ntl;
  cptr=strtok(tl, separs); if(cptr==NULL) return 2;
  ti=ri=0;
  while(cptr!=NULL) {
    if(ti==0) {
      if(strlen(cptr)>6 && strncasecmp(cptr, "Time[", 5)==0) {
        strncpy(unit, cptr+5, MAX_UNITS_LEN); unit[MAX_UNITS_LEN]=(char)0;
        cptr2=strchr(unit, ']'); if(cptr2!=NULL) *cptr2=(char)0;
        if(dft!=NULL) dft->timeunit=petTunitId(unit);
        timetype=DFT_TIME_MIDDLE;
      } else if(strlen(cptr)>6 && strncasecmp(cptr, "start[", 6)==0) {
        strncpy(unit, cptr+6, MAX_UNITS_LEN); unit[MAX_UNITS_LEN]=(char)0;
        cptr2=strchr(unit, ']'); if(cptr2!=NULL) *cptr2=(char)0;
        if(dft!=NULL) dft->timeunit=petTunitId(unit);
        timetype=DFT_TIME_STARTEND;
      }
    } else if(ti==1 && timetype==DFT_TIME_STARTEND) {
      if(strlen(cptr)>5 && strncasecmp(cptr, "end[", 4)==0) {
        strncpy(unit, cptr+4, MAX_UNITS_LEN); unit[MAX_UNITS_LEN]=(char)0;
        cptr2=strchr(unit, ']'); if(cptr2!=NULL) *cptr2=(char)0;
        ret=petCunitId(unit);
        if(ret!=CUNIT_UNKNOWN && dft!=NULL) strcpy(dft->unit, petCunit(ret));
      }
    } else {
      // get tac name
      strncpy(rname, cptr, MAX_REGIONNAME_LEN);
      rname[MAX_REGIONNAME_LEN]=(char)0;
      cptr2=strchr(rname, '['); if(cptr2!=NULL) *cptr2=(char)0;
      while((cptr2=strchr(rname, '_'))!=NULL) *cptr2=' ';
      if(dft!=NULL && ri<dft->_voidataNr) strcpy(dft->voi[ri].name, rname);
      // get unit
      cptr2=strchr(cptr, '['); 
      if(cptr2!=NULL) {
        strncpy(unit, cptr2+1, MAX_UNITS_LEN); unit[MAX_UNITS_LEN]=(char)0;
        cptr2=strchr(unit, ']'); if(cptr2!=NULL) *cptr2=(char)0;
        ret=petCunitId(unit);
        if(ret!=CUNIT_UNKNOWN && dft!=NULL) strcpy(dft->unit, petCunit(ret));
      }
      ri++;
    }
    ti++;
    cptr=strtok(NULL, separs);
  }
  if(dft==NULL) {
    return(ri);
  }

  dft->timetype=timetype;
  
  //for(int i=0; i<dft->_voidataNr; i++) printf("    '%s'\n", dft->voi[i].name);  
  
  /* Split ROI names */
  ti=ri;
  for(ri=0; ri<ti; ri++)
    rnameSplit(dft->voi[ri].name, dft->voi[ri].voiname,
          dft->voi[ri].hemisphere, dft->voi[ri].place, MAX_REGIONSUBNAME_LEN);

  return(0);
}
/*****************************************************************************/

/*****************************************************************************/

