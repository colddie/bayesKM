/// @file result.c
/// @author Vesa Oikonen
/// @brief IO for result files and handling RES struct data.
///
/*****************************************************************************/
#include "libtpccurveio.h"
#include <unistd.h>
/*****************************************************************************/
/* Local functions */
/// @cond
int resQSortComp(const void *f1, const void *f2);
int resQSortName(const void *voi1, const void *voi2);
/// @endcond
/*****************************************************************************/

/*****************************************************************************/
/** Free memory allocated for results. All data are cleared. */
void resEmpty(RES *res)
{
  if(res==NULL) return;
  if(res->_voidataNr>0) {
    free((char*)(res->voi));
    res->_voidataNr=0;
  }
  res->voiNr=0;
  res->parNr=0;
  res->studynr[0]=(char)0;
  for(int pi=0; pi<MAX_RESPARAMS; pi++) {
    strcpy(res->parname[pi], ""); strcpy(res->parunit[pi], "");}
  res->titleline[0]=res->unitline[0]=(char)0;
  res->program[0]=(char)0; res->refroi[0]=(char)0; res->datarange[0]=(char)0;
  res->datanr=0; res->fitmethod[0]=(char)0;
  res->datafile[0]=res->reffile[0]=(char)0;
  res->plasmafile[0]=res->plasmafile2[0]=res->bloodfile[0]=(char)0;
  res->density=res->lc=res->concentration=res->beta=0.0;
  res->Vb=-1.0;
  res->fA=-1.0;
  res->E=-1.0;
}
/*****************************************************************************/

/*****************************************************************************/
/** Initiate RES structure. This should be called once before first use. */
void resInit(RES *res)
{
  if(res==NULL) return;
  memset(res, 0, sizeof(RES));
  res->_voidataNr=0; res->voiNr=0; res->parNr=0;
  //res->Vb=-1.0;
  //res->fA=-1.0;
  //res->E=-1.0;
  resEmpty(res);
}
/*****************************************************************************/

/*****************************************************************************/
/** Allocate memory for result data. Old data is destroyed. */
int resSetmem(
  /** Pointer to initiated and possibly allocated result data */
  RES *res,
  /** Nr of regional results */
  int voiNr
) {
  int ri, pi;

  /* Check that there is something to do */
  if(res==NULL || voiNr<1) return(1);

  /* Clear previous data, but only if necessary */
  if(res->_voidataNr>0 || res->voiNr>0) resEmpty(res);

  /* Allocate memory for regional curves */
  res->voi=(ResVOI*)calloc(voiNr, sizeof(ResVOI));
  if(res->voi==NULL) return(2);
  res->_voidataNr=voiNr;

  /* Set SDs and CLs to NA */
  for(ri=0; ri<res->_voidataNr; ri++) for(pi=0; pi<MAX_RESPARAMS; pi++)
    res->voi[ri].sd[pi]=res->voi[ri].cl1[pi]=res->voi[ri].cl2[pi]=nan("");

  return(0);
}
/*****************************************************************************/

/*****************************************************************************/
/** Fix result parameter names and units, so that both representations are
 *  filled correctly, that is, the new string lists *parname[] and *parunit[],
 *  and the deprecated titleline[] and unitline[].
 *  New representation, if filled, always overwrites the deprecated one.
 *  Units are assumed to follow parameter name representation.
 */
void resFixParnames(
  /** Pointer to RES struct */
  RES *res
) {
  int i, len, n;
  char *cptr, tmp[1024], *lptr;

  if(RESULT_TEST>0) printf("resFixParnames(*res)\n");
  if(res==NULL) return;
  if(res->parNr<1) return;
  if(res->parNr>MAX_RESPARAMS) res->parNr=MAX_RESPARAMS;

  /* If new string lists are filled, then copy those to old representation */
  for(i=n=0; i<res->parNr; i++) {
    len=strlen(res->parname[i]); if(len<1) continue;
    if(strcmp(res->parname[i], ".")==0) continue;
    n++;
  }
  if(n>0) { // copy names and units
    strcpy(res->titleline, "");
    for(i=0; i<res->parNr; i++) {
      if(1023<(1+strlen(res->titleline)+strlen(res->parname[i]))) break;
      if(i>0) strcat(res->titleline, " ");
      len=strlen(res->parname[i]);
      if(len<1) strcat(res->titleline, ".");
      else strcat(res->titleline, res->parname[i]);
    }
    strcpy(res->unitline, "");
    for(i=0; i<res->parNr; i++) {
      if(1023<(1+strlen(res->unitline)+strlen(res->parunit[i]))) break;
      if(i>0) strcat(res->unitline, " ");
      len=strlen(res->parunit[i]);
      if(len<1) strcat(res->unitline, ".");
      else strcat(res->unitline, res->parunit[i]);
    }
    if(RESULT_TEST>1) {
      printf("Parameter names and units:\n");
      for(i=0; i<res->parNr; i++)
        printf("  %d: '%s' '%s'\n", i+1, res->parname[i], res->parunit[i]);
      printf("Created titleline: %s\n", res->titleline);
      printf("Created unitline: %s\n", res->unitline);
    }
    return;
  }

  /* If new string lists are not filled, then get them from deprecated strings */
  for(i=0; i<res->parNr; i++) {
    strcpy(res->parname[i], "");
    strcpy(res->parunit[i], "");
  }
  strcpy(tmp, res->titleline); lptr=tmp; cptr=strtok(lptr, " \t\n\r");
  i=0; while(cptr!=NULL && i<res->parNr) {
    if(strcmp(cptr, ".")==0) {i++; continue;}
    strncpy(res->parname[i], cptr, MAX_RESPARNAME_LEN);
    res->parname[i][MAX_RESPARNAME_LEN]=(char)0;
    cptr=strtok(NULL, " \t\n\r"); i++;
  }
  strcpy(tmp, res->unitline); lptr=tmp; cptr=strtok(lptr, " \t\n\r");
  i=0; while(cptr!=NULL && i<res->parNr) {
    if(strcmp(cptr, ".")==0) {i++; continue;}
    strncpy(res->parunit[i], cptr, MAX_RESPARNAME_LEN);
    res->parunit[i][MAX_RESPARNAME_LEN]=(char)0;
    cptr=strtok(NULL, " \t\n\r"); i++;
  }
  if(RESULT_TEST>1) {
    printf("Original titleline: %s\n", res->titleline);
    printf("Original unitline: %s\n", res->unitline);
    printf("Resolved parameter names and units:\n");
    for(i=0; i<res->parNr; i++)
      printf("  %d: '%s' '%s'\n", i+1, res->parname[i], res->parunit[i]);
  }
  return;
}
/*****************************************************************************/

/*****************************************************************************/
/** Print to stdout the contents of RES data structure. */
void resPrint(
  /** Pointer to result data */
  RES *res
) {
  resWrite(res, "stdout", 0);
}
/*****************************************************************************/

/*****************************************************************************/
/** Read RES file contents to the specified data structure.
\return In case of an error, >0 is returned, and a description is written in
    reserrmsg.
 */
int resRead(
  /** Result filename */
  char *filename,
  /** Pointer to initiated RES structure; any previous contents are deleted */
  RES *res,
  /** Verbose level; if <=0, then nothing is printed into stdout */
  int verbose
) {
  FILE *fp;
  char *cptr, line[1024], *lptr, tmp[1024];
  int i, n;
  fpos_t file_loc;


  if(verbose>0) printf("resRead(%s, *res);\n", filename);
  if(res==NULL) return(1);
  /* Empty data */
  resEmpty(res);
  res->isweight=-1; /* unknown */

  /* Open file; note that 'b' is required for fgetpos() and fsetpos() to work
     correctly */
  fp=fopen(filename, "rb");
  if(fp==NULL) {strcpy(reserrmsg, "cannot open file"); return(1);}

  /*
   * Read data, each result set separately, saving only the first one
   */
  strcpy(reserrmsg, "wrong format");

  /* Read program name */
  if(verbose>1) printf("reading program name\n");
  while(fgets(line, 1024, fp)!=NULL) {
    /* Ignore empty and comment lines */
    if(strlen(line)<4 || line[0]=='#') continue; else break;
  }
  /* Check for string (c) or (C) */
  if(strstr(line, "(c)") || strstr(line, "(C)")) {
    i=strlen(line)-1; while(i>0 && isspace(line[i])) i--; line[i+1]=(char)0;
    if(i<4) {fclose(fp); return(2);}
    strcpy(res->program, line);
  } else {fclose(fp); return(2);}

  /* Read calculation date and time */
  if(verbose>1) printf("reading date and time\n");
  while(fgets(line, 1024, fp)!=NULL) if(strlen(line)>2 && line[0]!='#') break;
  if(strncasecmp(line, "Date:", 5)) {fclose(fp); return(4);}
  cptr=&line[5]; while(isblank(cptr[0])) cptr++;
  if(cptr!=NULL) {
    if(verbose>3) printf("date_str := %s", cptr);
    struct tm st;
    if(get_datetime(cptr, &st, verbose-3)==0) res->time=timegm(&st);
    else if(get_date(cptr, &st)==0) res->time=timegm(&st);
    else res->time=(time_t)0;
  }

  /* Read studynr, datafiles, ref region, data range, etc */
  if(verbose>1) printf("reading headers\n");
  do {
    /* Omit comment and too short lines */
    while((lptr=fgets(line, 1024, fp))!=NULL) {
      if(strlen(line)>2 && line[0]!='#') break;
      if(lptr!=NULL && verbose>6)
        printf("omitted line[%d] := %s", (int)strlen(line), line);
    }
    if(lptr==NULL) break;
    strcpy(tmp, line);
    n=0; if(verbose>3) printf("line[%d] := %s", (int)strlen(line), line); 
    if(strncasecmp(line, "Study", 5)==0) {
      lptr=&tmp[6]; cptr=strtok(lptr, " \t\n\r"); n=1;
      if(cptr!=NULL && strlen(cptr)<1024) {
        strlcpy(res->studynr, cptr, MAX_STUDYNR_LEN+1);
      }  
    } else if(strncasecmp(line, "Data file", 9)==0) {
      lptr=&tmp[10]; cptr=strtok(lptr, " \t\n\r"); n=1;
      if(cptr!=NULL && strlen(cptr)<1024) strcpy(res->datafile, cptr);
    } else if(strncasecmp(line, "ROI file", 8)==0) {
      lptr=&tmp[9]; cptr=strtok(lptr, " \t\n\r"); n=1;
      if(cptr!=NULL && strlen(cptr)<1024) strcpy(res->datafile, cptr);
    } else if(strncasecmp(line, "Plasma file", 11)==0) {
      cptr=strchr(tmp, ':')+1; while(isspace((int)*cptr)) cptr++; n=1;
      if(cptr!=NULL && strlen(cptr)<1024) strcpy(res->plasmafile, cptr);
      cptr=res->plasmafile+strlen(res->plasmafile)-1;
      while(isspace((int)*cptr)) {*cptr=(char)0; cptr--;}
    } else if(strncasecmp(line, "2nd Plasma file", 15)==0) {
      cptr=strchr(tmp, ':')+1; while(isspace((int)*cptr)) cptr++; n=1;
      if(cptr!=NULL && strlen(cptr)<1024) strcpy(res->plasmafile2, cptr);
      cptr=res->plasmafile2+strlen(res->plasmafile2)-1;
      while(isspace((int)*cptr)) {*cptr=(char)0; cptr--;}
    } else if(strncasecmp(line, "Blood file", 10)==0) {
      lptr=&tmp[11]; cptr=strtok(lptr, " \t\n\r"); n=1;
      if(cptr!=NULL && strlen(cptr)<1024) strcpy(res->bloodfile, cptr);
    } else if(strncasecmp(line, "Reference file", 14)==0) {
      lptr=&tmp[15]; cptr=strtok(lptr, " \t\n\r"); n=1;
      if(cptr!=NULL && strlen(cptr)<1024) strcpy(res->reffile, cptr);
    } else if(strncasecmp(line, "Reference region", 16)==0) {
      cptr=strchr(tmp, ':')+1; while(isspace((int)*cptr)) cptr++; n=1;
      if(cptr!=NULL && strlen(cptr)<64) strcpy(res->refroi, cptr);
      cptr=res->refroi+strlen(res->refroi)-1;
      while(isspace((int)*cptr)) {*cptr=(char)0; cptr--;}
    } else if(strncasecmp(line, "Fit time", 8)==0 ||
              strncasecmp(line, "Data range", 10)==0) {
      cptr=strchr(tmp, ':')+1; while(isspace((int)*cptr)) cptr++; n=1;
      if(cptr!=NULL && strlen(cptr)<128) strcpy(res->datarange, cptr);
      cptr=res->datarange+strlen(res->datarange)-1;
      while(isspace((int)*cptr)) {*cptr=(char)0; cptr--;}
    } else if(strncasecmp(line, "Data nr", 7)==0) {
      cptr=strchr(tmp, ':')+1; while(isspace((int)*cptr)) cptr++; n=1;
      if(cptr!=NULL) res->datanr=atoi(cptr);
    } else if(strncasecmp(line, "Fit method", 10)==0) {
      cptr=strchr(tmp, ':')+1; while(isspace((int)*cptr)) cptr++; n=1;
      if(cptr!=NULL && strlen(cptr)<128) strcpy(res->fitmethod, cptr);
      cptr=res->fitmethod+strlen(res->fitmethod)-1;
      while(isspace((int)*cptr)) {*cptr=(char)0; cptr--;}
    } else if(strncasecmp(line, "Tissue density", 14)==0) {
      cptr=strchr(tmp, ':')+1; while(isspace((int)*cptr)) cptr++; n=1;
      if(cptr!=NULL) res->density=atof_dpi(cptr);
    } else if(strncasecmp(line, "Lumped constant", 15)==0) {
      cptr=strchr(tmp, ':')+1; while(isspace((int)*cptr)) cptr++; n=1;
      if(cptr!=NULL) res->lc=atof_dpi(cptr);
    } else if(strncasecmp(line, "Concentration", 13)==0) {
      cptr=strchr(tmp, ':')+1; while(isspace((int)*cptr)) cptr++; n=1;
      if(cptr!=NULL) res->concentration=atof_dpi(cptr);
    } else if(strncasecmp(line, "Beta", 4)==0) {
      cptr=strchr(tmp, ':')+1; while(isspace((int)*cptr)) cptr++; n=1;
      if(cptr!=NULL) res->beta=atof_dpi(cptr);
    } else if(strncasecmp(line, "Vb", 2)==0) {
      cptr=strchr(tmp, ':')+1; while(isspace((int)*cptr)) cptr++; n=1;
      if(cptr!=NULL) res->Vb=atof_dpi(cptr);
    } else if(strncasecmp(line, "fA", 2)==0) {
      cptr=strchr(tmp, ':')+1; while(isspace((int)*cptr)) cptr++; n=1;
      if(cptr!=NULL) res->fA=atof_dpi(cptr);
    } else if(strncasecmp(line, "Extraction", 10)==0) {
      cptr=strchr(tmp, ':')+1; while(isspace((int)*cptr)) cptr++; n=1;
      if(cptr!=NULL) res->E=atof_dpi(cptr);
    } else if(strncasecmp(line, "Weighting", 9)==0) {
      cptr=strchr(tmp, ':')+1; while(isspace((int)*cptr)) cptr++; n=1;
      if(strncasecmp(cptr, "yes", 1)==0) res->isweight=1;
      else if(strncasecmp(cptr, "no", 1)==0) res->isweight=0;
      else res->isweight=-1;
    } else if(strncasecmp(line, "Data was weighted", 17)==0) {
      res->isweight=1; n=1;
    } else if(strncasecmp(line, "Data was not weighted", 21)==0) {
      res->isweight=0; n=1;
    } else if(strncasecmp(line, "Region", 6)==0) {
      /* Header end, stop here */
      n=0;
    } else { /* Ignore all other header lines */
      n=1;
    }
  } while(n);
  if(verbose>6) printf("quit reading headers\n");

  /* Read the result parameter title line */
  if(verbose>1) printf("reading parameter titles\n");
  if(verbose>6) printf("using previously read line[%d] := %s",
                           (int)strlen(line), line);
  if(strncasecmp(line, "Region", 6)) {fclose(fp); return(10);}
  strcpy(tmp, line);
  lptr=strpbrk(tmp+6, " \n\r\t"); if(lptr==NULL) {fclose(fp); return(10);}
  cptr=strtok(lptr, " \n\r\t");
  n=0;
  if(cptr!=NULL) {
    if(strcmp(cptr, ".")!=0) {
      strncpy(res->parname[n], cptr, MAX_RESPARNAME_LEN);
      res->parname[n][MAX_RESPARNAME_LEN]=(char)0;
    } else {strcpy(res->parname[n], "");}
    if(verbose>5) printf("  parname[%d] := '%s'\n", n, res->parname[n]);
    n++;
  } else {strcpy(res->parname[n], "");}
  while(cptr!=NULL) {
    cptr=strtok(NULL, " \n\r\t");
    if(cptr!=NULL && n<MAX_RESPARAMS) {
      if(strcmp(cptr, ".")!=0) {
        strncpy(res->parname[n], cptr, MAX_RESPARNAME_LEN);
        res->parname[n][MAX_RESPARNAME_LEN]=(char)0;
      } else {strcpy(res->parname[n], "");}
      if(verbose>5) printf("  parname[%d] := '%s'\n", n, res->parname[n]);
      n++;
    }
  }
  res->parNr=n; if(verbose>1) printf("parNr := %d\n", res->parNr);

  /* Read the result parameter unit line */
  if(verbose>2) printf("seeking unit line...\n");
  if(fgetpos(fp, &file_loc)!=0) {fclose(fp); return(20);}
  while(fgets(line, 1024, fp)!=NULL && strlen(line)<3) {
    if(verbose>6) printf("omitted line[%d] := %s", (int)strlen(line), line);
    /* save the file position where unit line or results really start */
    if(fgetpos(fp, &file_loc)!=0) {fclose(fp); return(20);}
  }
  if(verbose>5) printf("possible unit line[%d]: %s", (int)strlen(line), line);
  if(strncasecmp(line, "# Units :", 7)==0 ||
     strncasecmp(line, "Units : ", 7)==0 ||
     strncasecmp(line, "# Units: ", 7)==0 ||
     strncasecmp(line, "Units: ", 6)==0)
  {
    if(verbose>1) printf("reading parameter units\n");
    strcpy(tmp, line); lptr=strchr(tmp+5, ':');
    if(lptr!=NULL) lptr++; else lptr=strpbrk(tmp+6, " \n\r\t");
    cptr=strtok(lptr, " \n\r\t"); n=0;
    if(cptr!=NULL) {
      if(strcmp(cptr, ".")!=0) {
        strncpy(res->parunit[n], cptr, MAX_RESPARNAME_LEN);
        res->parunit[n][MAX_RESPARNAME_LEN]=(char)0;
      } else {strcpy(res->parunit[n], "");}
      if(verbose>5) printf("  parunit[%d] := '%s'\n", n, res->parunit[n]);
      n++;
    } else {strcpy(res->parunit[n], "");}
    while(cptr!=NULL) {
      cptr=strtok(NULL, " \n\r\t");
      if(cptr!=NULL && n<MAX_RESPARAMS) {
        if(strcmp(cptr, ".")!=0) {
          strncpy(res->parunit[n], cptr, MAX_RESPARNAME_LEN);
          res->parunit[n][MAX_RESPARNAME_LEN]=(char)0;
        } else strcpy(res->parunit[n], "");
        if(verbose>5) printf("  parunit[%d] := '%s'\n", n, res->parunit[n]);
        n++;
      }
    }
  } else {
    if(verbose>5) printf("  ... not identified as unit line.\n");
    /* move file pointer to the previous place */
    if(fsetpos(fp, &file_loc)!=0) {fclose(fp); return(20);}
  }

  /* Read the nr of result lines */
  if(verbose>1) printf("reading nr of results\n");
  /* Set bookmark to the start of result lines */
  if(fgetpos(fp, &file_loc)!=0) {fclose(fp); return(21);}
  n=0;
  while(fgets(line, 1024, fp)!=NULL) {
    if(verbose>6) printf("line[%d] := %s", (int)strlen(line), line);
    if(line[0]=='#' || line[0]==';') continue;
    i=strlen(line); if(i<2) continue; if(i<3) break;
    n++;
  }
  /* Return file pointer to the start of result lines */
  if(fsetpos(fp, &file_loc)!=0) {fclose(fp); return(22);}
  if(verbose>1) printf("nr of result lines is %d\n", n);
  if(n<1) {
    strcpy(reserrmsg, "invalid result lines");
    fclose(fp); return(23);
  }

  /* Allocate memory for regional results */
  if(verbose>2) printf("allocating memory\n");
  if(resSetmem(res, n)) {
    strcpy(reserrmsg, "cannot allocate memory");
    fclose(fp); return(25);
  }

  /* Read regional results */
  if(verbose>1) printf("reading results to memory\n");
  int separtab=0; // 0=space as separator, 1=tab as separator
  char separstr[12];
  res->voiNr=0;
  while(fgets(line, 1024, fp)!=NULL) {
    if(verbose>6) printf("line[%d] := %s", (int)strlen(line), line);
    if(line[0]=='#' || line[0]==';') continue;
    i=strlen(line); if(i<2) continue; if(i<3) break;
    strcpy(tmp, line);
    if(verbose>2) printf("reading result %d\n", 1+res->voiNr);

    /* Read region names */
    if(strchr(tmp, '\t')==NULL) separtab=0; else separtab=1;
    if(separtab) strcpy(separstr, "\t\n\r"); else strcpy(separstr, " \t\n\r");
    int tokenNr=strTokenNr(tmp, separstr);
    if(verbose>20) printf("  tokenNr := %d\n", tokenNr);
    if(!separtab) { // old format with space as separator
      if(verbose>20) printf("separator: space\n");
      n=strTokenNCpy(tmp, separstr, 1, res->voi[res->voiNr].voiname, MAX_REGIONSUBNAME_LEN+1);
      if(n==0) {fclose(fp); return(31);}
      n=strTokenNCpy(tmp, separstr, 2, res->voi[res->voiNr].hemisphere, MAX_REGIONSUBNAME_LEN+1);
      if(n==0) {fclose(fp); return(31);}
      n=strTokenNCpy(tmp, separstr, 3, res->voi[res->voiNr].place, MAX_REGIONSUBNAME_LEN+1);
      if(n==0) {fclose(fp); return(31);}
      rnameCatenate(res->voi[res->voiNr].name, MAX_REGIONNAME_LEN, res->voi[res->voiNr].voiname,
                    res->voi[res->voiNr].hemisphere, res->voi[res->voiNr].place, ' ');
    } else { // space can exist only in TAC name
      if(verbose>20) printf("separator: tab\n");
      n=strTokenNCpy(tmp, separstr, 1, res->voi[res->voiNr].name, MAX_REGIONNAME_LEN+1);
      if(n==0) {fclose(fp); return(31);}
      rnameSplit(res->voi[res->voiNr].name, res->voi[res->voiNr].voiname,
                 res->voi[res->voiNr].hemisphere, res->voi[res->voiNr].place, MAX_REGIONSUBNAME_LEN);
    }
    if(strcmp(res->voi[res->voiNr].voiname, ".")==0) res->voi[res->voiNr].voiname[0]=(char)0;
    if(strcmp(res->voi[res->voiNr].hemisphere, ".")==0) res->voi[res->voiNr].hemisphere[0]=(char)0;
    if(strcmp(res->voi[res->voiNr].place, ".")==0) res->voi[res->voiNr].place[0]=(char)0;
    if(verbose>18) {
      printf("  voiname := '%s'\n", res->voi[res->voiNr].voiname);
      printf("  hemisphere := '%s'\n", res->voi[res->voiNr].hemisphere);
      printf("  place := '%s'\n", res->voi[res->voiNr].place);
    }

    /* Read results, continuing from the pointer after region names */
    int tokeni=2; if(!separtab) tokeni=4;
    char buf[128];
    for(i=0; i<MAX_RESPARAMS && tokeni<=tokenNr; i++, tokeni++) {
      n=strTokenNCpy(tmp, separstr, tokeni, buf, 128);
      if(n==0) {fclose(fp); return(32);}
      if(strlen(buf)==1 && buf[0]=='.') res->voi[res->voiNr].parameter[i]=nan("");
      else res->voi[res->voiNr].parameter[i]=atof_dpi(buf);
    }
    if(verbose>5) printf("  for '%s' parNr:=%d\n", res->voi[res->voiNr].name, i);
    //if(res->voiNr==0) {res->parNr=i;} else if(i<res->parNr) res->parNr=i;
    if(i<res->parNr) {
      if(verbose>0) printf("Warning: smaller parNr %d on region '%s'\n", i, 
        res->voi[res->voiNr].name);
      res->parNr=i;
    }

    /* If 'region' name implies that this was confidence limit or sd, then */
    /* move the values into correct places, and do not increase the voiNr */
    if(res->voiNr==0) {res->voiNr++; continue;}
    if(strcasecmp(res->voi[res->voiNr].voiname, "CL")==0) {
      if(strcmp(res->voi[res->voiNr].hemisphere, "95%")==0) {
        if(strcasecmp(res->voi[res->voiNr].place, "Lower")==0)
          for(i=0; i<res->parNr; i++)
            res->voi[res->voiNr-1].cl1[i]=res->voi[res->voiNr].parameter[i];
        else if(strcasecmp(res->voi[res->voiNr].place, "Upper")==0)
          for(i=0; i<res->parNr; i++)
            res->voi[res->voiNr-1].cl2[i]=res->voi[res->voiNr].parameter[i];
        continue;
      }
    } else if(strcasecmp(res->voi[res->voiNr].voiname, "SD")==0) {
      for(i=0; i<res->parNr; i++)
        res->voi[res->voiNr-1].sd[i]=res->voi[res->voiNr].parameter[i];
      continue;
    }
    res->voiNr++;
  }
  if(res->parNr==0) {fclose(fp); return(33);}
  if(verbose>0) printf("nr of results: %d ; nr of parameters: %d\n",
    res->voiNr, res->parNr);

  /* Seek for other results in the same file */
  while(fgets(line, 1024, fp)!=NULL) {
    /* Ignore empty and comment lines */
    if(strlen(line)<3 || line[0]=='#') continue; else break;
  }
  /* Check again for string (c) or (C) */
  if(strstr(line, "(c)") || strstr(line, "(C)")) {
    fprintf(stderr, 
    "Warning: %s contains more than one set of results; only the 1st one is used.\n",
      filename);
  }

  /* Close file */
  fclose(fp);
  strcpy(reserrmsg, "");

  /* Fill studynr if it was not found in file */
  if(!res->studynr[0]) studynr_from_fname(filename, res->studynr);
  /* Set also deprecated parameter name and unit representations, for now */
  resFixParnames(res);

  return(0);
}
/*****************************************************************************/

/*****************************************************************************/
/** Write calculation results into specified file.
  If file exists, a backup file (+BACKUP_EXTENSION) is written also.
  If "stdout" is given as filename, output is directed to stdout.
  If filename extension is *.htm(l), file is saved in HTML format.
\return In case of an error, >0 is returned, and a description is written in reserrmsg.
*/
int resWrite(
  /** Pointer to result data */
  RES *res,
  /** Output filename */
  char *filename,
  /** Verbose level; if <=0, then nothing is printed into stdout */
  int verbose
) {
  int i, j, n;
  char tmp[1024], is_stdout=0, *cptr;
  FILE *fp;
  int partype[MAX_RESPARAMS]; /* 0=int, 1=double, 2=exp */
  int colwidth[MAX_RESPARAMS];
  double *p;


  if(verbose>1) printf("resWrite(*res, %s, %d)\n", filename, verbose);
  /* Check that there is some data to write */
  if(res==NULL) {strcpy(reserrmsg, "error in result data"); return(1);}
  if(res->voiNr<1) {strcpy(reserrmsg, "no result data"); return(1);}

  /* Write results in HTML format, if necessary */
  cptr=strrchr(filename, '.');
  if(cptr!=NULL && (!strncasecmp(cptr, ".htm", 4)))
    return(resWriteHTML(res, filename, verbose));

  /* Check if writing to stdout */
  if(!strcasecmp(filename, "stdout")) is_stdout=1;

  /* Check if file exists; backup, if necessary */
  if(!is_stdout) (void)backupExistingFile(filename, NULL, NULL);

  resFixParnames(res);

  /* Open output file */
  if(is_stdout) fp=(FILE*)stdout;
  else if((fp = fopen(filename, "w")) == NULL) {
    strcpy(reserrmsg, "cannot open file"); return(2);}

  /* Program name */
  n=fprintf(fp, "%s\n\n", res->program);
  if(n==0) {
    strcpy(reserrmsg, "disk full");
    if(!is_stdout) fclose(fp);
    return(3);
  }

  /* Write calculation date and time */
  if(!ctime_r_int(&res->time, tmp)) strcpy(tmp, "1900-01-01 00:00:00"); 
  fprintf(fp, "Date:\t%s\n", tmp);

  /* Write the studynr */
  if(res->studynr[0]) fprintf(fp, "Study:\t%s\n", res->studynr);

  /* Write the names of the original datafiles */
  if(res->datafile[0])   fprintf(fp, "Data file:\t%s\n", res->datafile);
  if(res->plasmafile[0]) fprintf(fp, "Plasma file:\t%s\n", res->plasmafile);
  if(res->plasmafile2[0]) fprintf(fp, "2nd Plasma file:\t%s\n", res->plasmafile2);
  if(res->bloodfile[0])  fprintf(fp, "Blood file:\t%s\n", res->bloodfile);
  if(res->reffile[0])    fprintf(fp, "Reference file:\t%s\n", res->reffile);
  if(res->refroi[0])     fprintf(fp, "Reference region:\t%s\n", res->refroi);

  /* Write data range etc */
  if(res->datarange[0])  fprintf(fp, "Data range:\t%s\n", res->datarange);
  if(res->datanr>0)      fprintf(fp, "Data nr:\t%d\n", res->datanr);
  if(res->fitmethod[0])  fprintf(fp, "Fit method:\t%s\n", res->fitmethod);

  /* Write constants */
  if(res->density>0.0)   fprintf(fp, "Tissue density:\t%g\n", res->density);
  if(res->lc>0.0)        fprintf(fp, "Lumped constant:\t%g\n", res->lc);
  if(res->concentration>0.0)
                         fprintf(fp, "Concentration:\t%g\n", res->concentration);
  if(res->beta>0.0)      fprintf(fp, "Beta:\t%g\n", res->beta);
  if(res->Vb>=0.0)       fprintf(fp, "Vb:\t%g %%\n", res->Vb);
  if(res->fA>=0.0)       fprintf(fp, "fA:\t%g %%\n", res->fA);
  if(res->E>=0.0)        fprintf(fp, "Extraction:\t%g\n", res->E);

  /* Weighting */
  if(res->isweight>0) strcpy(tmp, "yes");
  else if(res->isweight==0) strcpy(tmp, "no");
  else strcpy(tmp, "unknown");
  fprintf(fp, "Weighting:\t%s\n", tmp);

  /* Determine column widths: */
  for(j=0; j<res->parNr; j++) colwidth[j]=1;
  /* should column be printed as integers (0), floats (1) or exponentials (2)? */
  for(j=0; j<res->parNr; j++) {
    partype[j]=resParameterPrintType(res, j);
  }
  /* min width required by titles */
  for(j=0; j<res->parNr; j++) colwidth[j]=strlen(res->parname[j]);
  if(verbose>2) {
    printf("col widths after titles were checked:\n");
    for(i=0; i<res->parNr; i++)
      printf("  par%d : partype=%d colwidth=%d\n", i+1, partype[i], colwidth[i]);
  }
  /* min width required by units */
  for(j=0; j<res->parNr; j++) {
    n=strlen(res->parunit[j]); if(n>colwidth[j]) colwidth[j]=n;
  }
  if(verbose>2) {
    printf("col widths after units were checked:\n");
    for(i=0; i<res->parNr; i++)
      printf("  par%d : partype=%d colwidth=%d\n", i+1, partype[i], colwidth[i]);
  }
  /* widths required by result values */
  for(i=0; i<res->voiNr; i++) {
    p=res->voi[i].parameter;
    for(j=0; j<res->parNr; j++) {
      if(isnan(p[j])) continue;
      if(p[j]>=0) n=4; else n=3;
      switch(partype[j]) {
        case 0: sprintf(tmp, "%.0f", p[j]); break;
        case 1: sprintf(tmp, "%.*f", n, p[j]); break;
        case 2:
        default: sprintf(tmp, "%.*e", n, p[j]); break;
      }
      n=strlen(tmp); if(n>colwidth[j]) colwidth[j]=n;
    }
  }
  if(verbose>2) {
    printf("col widths after result values were checked:\n");
    for(i=0; i<res->parNr; i++)
      printf("  par%d : partype=%d colwidth=%d\n", i+1, partype[i], colwidth[i]);
  }

  /* Title line */
  if(verbose>4) {
    printf("  writing title line with %d parameter(s)\n", res->parNr); 
    fflush(stdout);
  }
  fprintf(fp, "\n%s\t", "Region");
  for(i=0; i<res->parNr; i++) {
    if(strlen(res->parname[i])<1) strcpy(tmp, ".");
    else strcpy(tmp, res->parname[i]);
    fprintf(fp, "\t%s", tmp);
  }
  fprintf(fp, "\n");

  /* Write units, if they exist, currently as comment line */
  for(i=j=0; i<res->parNr; i++) if(strlen(res->parunit[i])>0) j++;
  if(j>0) {
    if(verbose>4) {
      printf("  writing units line with %d parameter(s)\n", j); 
      fflush(stdout);
    }
    fprintf(fp, "%s", "# Units:");
    //fprintf(fp, "%s ", "# Units :");
    for(i=0; i<res->parNr; i++) {
      if(strlen(res->parunit[i])<1) strcpy(tmp, ".");
      else strcpy(tmp, res->parunit[i]);
      fprintf(fp, "\t%s", tmp);
    }
    fprintf(fp, "\n");
  }
  fflush(fp);

  /* Write regional results */
  if(verbose>4) {
    printf("  writing %d regional results\n", res->voiNr); fflush(stdout);
  }
  for(i=0; i<res->voiNr; i++) {
    if(verbose>6) {printf("    writing region %d\n", 1+i); fflush(stdout);}
#if(1) // new version
    if(res->voi[i].name[0]) {
      fprintf(fp, "%s", res->voi[i].name);
    } else {
      if(res->voi[i].voiname[0]) strcpy(tmp, res->voi[i].voiname); else strcpy(tmp, ".");
      fprintf(fp, "%.*s ", MAX_REGIONSUBNAME_LEN, tmp);
      if(res->voi[i].hemisphere[0]) strcpy(tmp, res->voi[i].hemisphere); else strcpy(tmp, ".");
      fprintf(fp, "%.*s ", MAX_REGIONSUBNAME_LEN, tmp);
      if(res->voi[i].place[0]) strcpy(tmp, res->voi[i].place); else strcpy(tmp, ".");
      fprintf(fp, "%.*s", MAX_REGIONSUBNAME_LEN, tmp);
    }
#else // previous version
    if(res->voi[i].voiname[0])
      strcpy(tmp, res->voi[i].voiname); else strcpy(tmp, ".");
    fprintf(fp, "%.*s ", MAX_REGIONSUBNAME_LEN, tmp);
    if(res->voi[i].hemisphere[0])
      strcpy(tmp, res->voi[i].hemisphere); else strcpy(tmp, ".");
    fprintf(fp, "%.*s ", MAX_REGIONSUBNAME_LEN, tmp);
    if(res->voi[i].place[0])
      strcpy(tmp, res->voi[i].place); else strcpy(tmp, ".");
    fprintf(fp, "%.*s", MAX_REGIONSUBNAME_LEN, tmp);
#endif
    p=res->voi[i].parameter;
    for(j=0; j<res->parNr; j++) {
      if(verbose>15) {printf("      writing par %d\n", 1+j); fflush(stdout);}
      if(isnan(p[j])) {fprintf(fp, "\t."); continue;}
      switch(partype[j]) {
        case 0: fprintf(fp, "\t%.0f", p[j]); break;
        case 1:
          if(p[j]>=0) n=4; else n=3;
          fprintf(fp, "\t%.*f", n, p[j]);
          break;
        default:
          if(p[j]>=0) n=4; else n=3;
          fprintf(fp, "\t%.*e", n, p[j]);
        break;
      }
    }
    fprintf(fp, "\n"); fflush(fp);
    /* Write SD's, if they exist */
    if(verbose>25) {printf("    sd?\n"); fflush(stdout);}
    //if(res->voi[i].sd==NULL) {printf("NULL\n"); fflush(stdout);}
    if(verbose>25) {printf("    parNr=%d\n", res->parNr); fflush(stdout);}
    n=0; for(int j=0; j<res->parNr; j++) {if(!isnan(res->voi[i].sd[j])) n++;}
    if(verbose>25) {printf("      n=%d\n", n); fflush(stdout);}
    if(n>0) {
      fprintf(fp, "SD . .");
      for(j=0; j<res->parNr; j++) {
        if(verbose>15) {printf("      writing sd %d\n", 1+j); fflush(stdout);}
        if(!isnan(res->voi[i].sd[j])) {
          switch(partype[j]) {
            case 0: fprintf(fp, "\t%.0f", res->voi[i].sd[j]); break;
            case 1:
              if(res->voi[i].sd[j]>=0) n=4; else n=3;
              fprintf(fp, "\t%.*f", n, res->voi[i].sd[j]);
              break;
            default:
              if(res->voi[i].sd[j]>=0) n=4; else n=3;
              fprintf(fp, "\t%.*e", n, res->voi[i].sd[j]);
              break;
          }
        } else {
          fprintf(fp, "\t.");
        }
      }
      fprintf(fp, "\n"); fflush(fp);
    }
    /* Write lower confidence limits, if they exist */
    if(verbose>25) {printf("    cl1?\n"); fflush(stdout);}
    n=0; for(int j=0; j<res->parNr; j++) if(!isnan(res->voi[i].cl1[j])) n++;
    if(verbose>25) {printf("      n=%d\n", n); fflush(stdout);}
    if(n>0) {
      fprintf(fp, "CL 95%% Lower");
      for(j=0; j<res->parNr; j++) {
        if(verbose>15) {printf("      writing CL1 %d\n", 1+j); fflush(stdout);}
        if(!isnan(res->voi[i].cl1[j])) {
          switch(partype[j]) {
            case 0: fprintf(fp, "\t%.0f", res->voi[i].cl1[j]); break;
            case 1:
              if(res->voi[i].cl1[j]>=0) n=4; else n=3;
              fprintf(fp, "\t%.*f", n, res->voi[i].cl1[j]);
              break;
            default:
              if(res->voi[i].cl1[j]>=0) n=4; else n=3;
              fprintf(fp, "\t%.*e", n, res->voi[i].cl1[j]);
              break;
          }
        } else {
          fprintf(fp, "\t.");
        }
      }
      fprintf(fp, "\n");
    }
    /* Write upper confidence limits, if they exist */
    if(verbose>25) {printf("    cl2?\n"); fflush(stdout);}
    n=0; for(int j=0; j<res->parNr; j++) if(!isnan(res->voi[i].cl2[j])) n++;
    if(verbose>25) {printf("      n=%d\n", n); fflush(stdout);}
    if(n>0) {
      fprintf(fp, "CL 95%% Upper");
      for(j=0; j<res->parNr; j++) {
        if(verbose>15) {printf("      writing CL2 %d\n", 1+j); fflush(stdout);}
        if(!isnan(res->voi[i].cl2[j])) {
          switch(partype[j]) {
            case 0: fprintf(fp, "\t%.0f", res->voi[i].cl2[j]); break;
            case 1:
              if(res->voi[i].cl2[j]>=0) n=4; else n=3;
              fprintf(fp, "\t%.*f", n, res->voi[i].cl2[j]);
              break;
            default:
              if(res->voi[i].cl2[j]>=0) n=4; else n=3;
              fprintf(fp, "\t%.*e", n, res->voi[i].cl2[j]);
              break;
          }
        } else {
          fprintf(fp, "\t.");
        }
      }
      fprintf(fp, "\n"); fflush(fp);
    }
  } /* next region */

  /* Close file */
  if(!is_stdout) {fflush(fp); fclose(fp);}
  strcpy(reserrmsg, "");
  if(verbose>1) {printf("resWrite() done.\n"); fflush(stdout);}

  return(0);
}
/*****************************************************************************/

/*****************************************************************************/
/** Write calculation results into specied XHTML 1.1 file.
    If file exists, a backup file (+BACKUP_EXTENSION) is written also.
    If "stdout" is given as filename, output is directed to stdout
\return In case of an error, >0 is returned, and a description is written in
    reserrmsg.
 */
int resWriteHTML(
  /** Pointer to result data */
  RES *res,
  /** Output file name */
  char *fname,
  /** Verbose level; if <=0, then nothing is printed into stdout */
  int verbose
) {
  int n, is_stdout=0;
  char tmp[1024];
  FILE *fp;


  if(verbose>0) printf("resWriteHTML(*res, %s, %d)\n", fname, verbose);
  /* Check that there is some data to write */
  if(res==NULL) {strcpy(reserrmsg, "error in result data"); return(1);}
  if(res->voiNr<1) {strcpy(reserrmsg, "no result data"); return(1);}
  /* Check if writing to stdout */
  if(!strcasecmp(fname, "stdout")) is_stdout=1;

  resFixParnames(res);

  /* Check if file exists; backup, if necessary */
  if(!is_stdout && access(fname, 0) != -1) {
    strcpy(tmp, fname); strcat(tmp, BACKUP_EXTENSION);
    if(access(tmp, 0)!=-1) remove(tmp);
    rename(fname, tmp);
  }
  strcpy(reserrmsg, "cannot write file");

  /* Open output file */
  if(is_stdout) fp=(FILE*)stdout;
  else if((fp=fopen(fname, "w"))==NULL) {
    strcpy(reserrmsg, "cannot open file"); return(2);}

  /* Write XHTML 1.1 doctype and head */
  if(resWriteXHTML11_doctype(fp) || resWriteXHTML11_head(fp, res->program)) {
    strcpy(reserrmsg, "disk full");
    if(!is_stdout) fclose(fp);
    return(3);
  }

  /* Start writing the body of the HTML file */
  fprintf(fp, "\n<body>\n");

  /* Start the div for tables */
  fprintf(fp, "\n<div id=\"tables\">\n");

  /* Write results into a table */
  if(resWriteHTML_table(res, fp)!=0) {
    strcpy(reserrmsg, "disk full");
    if(!is_stdout) fclose(fp);
    return(3);
  }

  /* Stop writing the body of the HTML file, and end the file */
  fprintf(fp, "</div>\n");
  n=fprintf(fp, "</body></html>\n\n");
  if(n==0) {
    strcpy(reserrmsg, "disk full");
    if(!is_stdout) fclose(fp);
    return(3);
  }

  /* Close file */
  if(!is_stdout) {fflush(fp); fclose(fp);}
  strcpy(reserrmsg, "");

  return(0);
}
/*****************************************************************************/

/*****************************************************************************/
/** Write XHTML 1.1 doctype into an opened file pointer.
\return Returns 0 if successful, non-zero in case of a failure.
 */
int resWriteXHTML11_doctype(
  /** Output file pointer */
  FILE *fp
) {
  int n;
  if(fp==NULL) return(1);
  n=fprintf(fp, "<!DOCTYPE html PUBLIC \"-//W3C//DTD XHTML 1.1//EN\" ");
  n+=fprintf(fp, "\"http://www.w3.org/TR/xhtml11/DTD/xhtml11.dtd\">\n");
  if(n<20) return(2);
  n=fprintf(fp,"<html xmlns=\"http://www.w3.org/1999/xhtml\" xml:lang=\"en\">\n\n");
  if(n<20) return(2);
  return(0);
}
/*****************************************************************************/

/*****************************************************************************/
/** Write XHTML 1.1 head for PET results file into an opened file pointer.
\return Returns 0 if successful, non-zero in case of a failure.
 */
int resWriteXHTML11_head(
  /** File pointer where to write */
  FILE *fp,
  /** Author name, for example software name */
  char *author_name
) {
  int n;
  if(fp==NULL) return(1);

  n=fprintf(fp, "<head>\n"); if(n<6) return(2);
  fprintf(fp, "  <title>Tabulated PET results</title>\n");
  fprintf(fp, "  <meta http-equiv=\"content-type\" content=\"text/html; ");
  fprintf(fp, "charset=iso-8859-1\" />\n");
  fprintf(fp, "  <meta http-equiv=\"content-language\" content=\"en-gb\" />\n");
  fprintf(fp, "  <meta name=\"description\" content=\"Regional PET results\" />\n");
  fprintf(fp, "  <meta name=\"author\" content=\"%s\" />\n", author_name);
  fprintf(fp, "  <meta name=\"ProgId\" content=\"Excel.Sheet\" />\n");
  fprintf(fp, "  <link rel=\"icon\" href=\"http://www.turkupetcentre.net/favi");
  fprintf(fp, "con.ico\" type=\"image/x-icon\" />\n");
  fprintf(fp, "  <link rel=\"shortcut icon\" href=\"http://www.turkupetcentre");
  fprintf(fp, ".net/favicon.ico\" type=\"image/x-icon\" />\n");
  /* write local CSS with basic settings in case external CSS is not available */
  fprintf(fp, "  <style type=\"text/css\">\n");
  fprintf(fp, "    thead {background-color:#999999; color:black;}\n");
//fprintf(fp, "    table {text-align:left; width:100%%; border-collapse:colla");
  fprintf(fp, "    table {text-align:left; border-collapse:collapse;");
  fprintf(fp, " empty-cells:show;}\n");
//fprintf(fp, "    td {border:1px solid black;}\n");
  fprintf(fp, "    .oddroi {background-color:#FFFFFF; color:black;}\n");
  fprintf(fp, "    .evenroi {background-color:#CCCCCC; color:black;}\n");
  fprintf(fp, "    .sd {background-color:#999999; color:blue;}\n");
  fprintf(fp, "    .cl1 {background-color:#999999; color:green;}\n");
  fprintf(fp, "    .cl2 {background-color:#999999; color:green;}\n");
  fprintf(fp, "    .oddstudy {background-color:#FFFFFF; color:black;}\n");
  fprintf(fp, "    .evenstudy {background-color:#CCCCCC; color:black;}\n");
  fprintf(fp, "    .oddsum {background-color:#999999; color:black;}\n");
  fprintf(fp, "    .evensum {background-color:#CCCCCC; color:black;}\n");
  fprintf(fp, "    #regcontainer ul {margin-left:0; padding-left:0;}\n");
  fprintf(fp, "    #regcontainer ul li {display:inline; list-style-type:none;}\n");
  fprintf(fp, "    #regcontainer a {padding:2px 4px;}\n");
  fprintf(fp, "    <!--table\n");
  fprintf(fp, "    	{mso-displayed-decimal-separator:\"\\.\";\n");
  fprintf(fp, "    	mso-displayed-thousand-separator:\" \";}\n");
  fprintf(fp, "    -->\n");
  fprintf(fp, "  </style>\n");
  /* load external CSS with more fancy settings */
  fprintf(fp, "  <link rel=\"stylesheet\" type=\"text/css\" href=\"http://www");
  fprintf(fp, ".turkupetcentre.net/result.css\" />\n");
  fprintf(fp, "</head>\n");
  if(n<7) return(2);
  return(0);
}
/*****************************************************************************/

/*****************************************************************************/
/** Write calculation results as one HTML table into an opened (X)HTML file.
\return Returns 0 if successful, non-zero in case of a failure.
 */
int resWriteHTML_table(
  /** Pointer to result data */
  RES *res,
  /** Output file pointer */
  FILE *fp
) {
  int i, j, n;
  char *cptr, tmp[1024];
  int partype[MAX_RESPARAMS]; /* 0=int, 1=double, 2=exp */
  double *p;

  if(RESULT_TEST>0) printf("resWriteHTML_table(*res, fp)\n");
  if(fp==NULL || res==NULL) return(1);

  /* Write the title lines to the head of table */
  fprintf(fp, "<table>\n");
  fprintf(fp, "<thead>\n");
  /* Program name */
  strcpy(tmp, res->program);
  cptr=strcasestr(tmp, "(c)");
  if(cptr!=NULL) {
    *cptr=(char)0;
     i=strlen(tmp)-1; while(i>0 && isspace(tmp[i])) i--; tmp[i+1]=(char)0;
  }
  if(tmp[0])
    fprintf(fp, "<tr align=left><th>Program:</th><td>%s</td></tr>\n", tmp);
  /* Write calculation date and time */
  if(!ctime_r_int(&res->time, tmp)) strcpy(tmp, "1900-01-01 00:00:00");
  fprintf(fp, "<tr align=left><th>Date:</th><td>%s</td></tr>\n", tmp);
  /* Write the studynr */
  if(res->studynr[0])
    fprintf(fp, "<tr align=left><th>Study:</th><td>%s</td></tr>\n",
            res->studynr);
  /* Write the names of the original datafiles */
  if(res->datafile[0])
    fprintf(fp, "<tr align=left><th>Data file:</th><td>%s</td></tr>\n",
            res->datafile);
  if(res->plasmafile[0])
    fprintf(fp, "<tr align=left><th>Plasma file:</th><td>%s</td></tr>\n",
            res->plasmafile);
  if(res->plasmafile2[0])
    fprintf(fp, "<tr align=left><th>2nd Plasma file:</th><td>%s</td></tr>\n",
            res->plasmafile2);
  if(res->bloodfile[0])
    fprintf(fp, "<tr align=left><th>Blood file:</th><td>%s</td></tr>\n",
            res->bloodfile);
  if(res->reffile[0])
    fprintf(fp, "<tr align=left><th>Reference file:</th><td>%s</td></tr>\n",
            res->reffile);
  if(res->refroi[0])
    fprintf(fp, "<tr align=left><th>Reference region:</th><td>%s</td></tr>\n",
            res->refroi);
  /* Write data range etc */
  if(res->datarange[0])
    fprintf(fp, "<tr align=left><th>Data range:</th><td>%s</td></tr>\n",
            res->datarange);
  if(res->datanr>0)
    fprintf(fp, "<tr align=left><th>Data nr:</th><td>%d</td></tr>\n",
            res->datanr);
  if(res->fitmethod[0])
    fprintf(fp, "<tr align=left><th>Fit method:</th><td>%s</td></tr>\n",
            res->fitmethod);
  /* Write constants */
  if(res->density>0.0)
    fprintf(fp, "<tr align=left><th>Tissue density:</th><td>%g</td></tr>\n",
            res->density);
  if(res->lc>0.0)
    fprintf(fp, "<tr align=left><th>Lumped constant:</th><td>%g</td></tr>\n",
            res->lc);
  if(res->concentration>0.0)
    fprintf(fp, "<tr align=left><th>Concentration:</th><td>%g</td></tr>\n",
            res->concentration);
  if(res->beta>0.0)
    fprintf(fp, "<tr align=left><th>Beta:</th><td>%g</td></tr>\n", res->beta);
  if(res->Vb>=0.0)
    fprintf(fp, "<tr align=left><th>Vb:</th><td>%g %%</td></tr>\n", res->Vb);
  if(res->fA>=0.0)
    fprintf(fp, "<tr align=left><th>fA:</th><td>%g %%</td></tr>\n", res->fA);
  if(res->E>0.0)
    fprintf(fp, "<tr align=left><th>Extraction:</th><td>%g</td></tr>\n", res->E);
  /* Weighting */
  if(res->isweight>0) strcpy(tmp, "yes");
  else if(res->isweight==0) strcpy(tmp, "no");
  else strcpy(tmp, "unknown");
  fprintf(fp, "<tr align=left><th>Weighting:</th><td>%s</td></tr>\n", tmp);
  /* End the head of table */
  fprintf(fp, "</thead>\n");

  /* Collect info on column types */
  for(j=0; j<res->parNr; j++) {
    /* should column be printed as integers, floats or exponentials? */
    partype[j]=resParameterPrintType(res, j);
  }

#if(0)
  int rnameSubfields=resRNameSubfieldExists(res);
#endif

  /* Write the result data to the table body */
  fprintf(fp, "<tbody>\n");
  /* Write the result title line */
#if(1)
  fprintf(fp,"<tr align=left><th>Region</th>\n");
#else
  if(rnameSubfields==0)
    fprintf(fp,"<tr align=left><th>Region</th>\n");
  else
    fprintf(fp,"<tr align=left><th>Region</th><th>Hemisphere</th><th>Plane</th>\n");
#endif
  for(i=0; i<res->parNr; i++) fprintf(fp, "<th>%s</th>", res->parname[i]);
  fprintf(fp, "</tr>\n");
  /* Write regional results */
  for(i=0; i<res->voiNr; i++) {
    if(i%2) strcpy(tmp, "evenroi"); else strcpy(tmp, "oddroi");
    fprintf(fp, "<tr class=\"%s\">", tmp);
    fprintf(fp, "<th>%s</th>", res->voi[i].name);
#if(0)
    if(rnameSubfields!=0) {
      fprintf(fp, "<th>%s</th>", res->voi[i].hemisphere);
      fprintf(fp, "<th>%s</th>", res->voi[i].place);
    }
#endif
    p=res->voi[i].parameter;
    for(j=0; j<res->parNr; j++) switch(partype[j]) {
      case 0: fprintf(fp, "<td>%.0f</td>", p[j]); break;
      case 1: fprintf(fp, "<td>%.4f</td>", p[j]); break;
      default: fprintf(fp, "<td>%.4E</td>", p[j]); break;
    }
    fprintf(fp, "</tr>\n");
    /* Write SD's, if they exist */
    p=res->voi[i].sd; for(j=n=0; j<res->parNr; j++) if(!isnan(p[j])) n++;
    if(n>0) {
      fprintf(fp, "<tr class=\"sd\">");
      fprintf(fp, "<th>%s</th>", "SD");
#if(0)
      if(rnameSubfields!=0) fprintf(fp, "<th></th><th></th>");
#endif
      for(j=0; j<res->parNr; j++) {
        if(!isnan(p[j])) switch(partype[j]) {
          case 0: fprintf(fp, "<td>%.0f</td>", p[j]); break;
          case 1: fprintf(fp, "<td>%.4f</td>", p[j]); break;
          default: fprintf(fp, "<td>%.4E</td>", p[j]); break;
        } else fprintf(fp, "<td></td>");
      }
      fprintf(fp, "</tr>\n");
    }
    /* Write lower confidence limits, if they exist */
    p=res->voi[i].cl1; for(j=n=0; j<res->parNr; j++) if(!isnan(p[j])) n++;
    if(n>0) {
      fprintf(fp, "<tr class=\"cl1\">");
#if(1)
      fprintf(fp, "<th>CL95%%L</th>");
#else
      if(rnameSubfields!=0)
        fprintf(fp, "<th>%s</th><th>%s</th><th>%s</th>", "CL", "95%", "Lower");
      else
        fprintf(fp, "<th>CL95%%L</th>");
#endif
      for(j=0; j<res->parNr; j++) {
        if(!isnan(p[j])) switch(partype[j]) {
          case 0: fprintf(fp, "<td>%.0f</td>", p[j]); break;
          case 1: fprintf(fp, "<td>%.4f</td>", p[j]); break;
          default: fprintf(fp, "<td>%.4E</td>", p[j]); break;
        } else fprintf(fp, "<td></td>");
      }
      fprintf(fp, "</tr>\n");
    }
    /* Write upper confidence limits, if they exist */
    p=res->voi[i].cl2; for(j=n=0; j<res->parNr; j++) if(!isnan(p[j])) n++;
    if(n>0) {
      fprintf(fp, "<tr class=\"cl2\">");
#if(1)
      fprintf(fp, "<th>CL95%%U</th>");
#else
      if(rnameSubfields!=0)
        fprintf(fp, "<th>%s</th><th>%s</th><th>%s</th>", "CL", "95%", "Upper");
      else
        fprintf(fp, "<th>CL95%%U</th>");
#endif
      for(j=0; j<res->parNr; j++) {
        if(!isnan(p[j])) switch(partype[j]) {
          case 0: fprintf(fp, "<td>%.0f</td>", p[j]); break;
          case 1: fprintf(fp, "<td>%.4f</td>", p[j]); break;
          default: fprintf(fp, "<td>%.4E</td>", p[j]); break;
        } else fprintf(fp, "<td></td>");
      }
      fprintf(fp, "</tr>\n");
    }
  }
  /* End the body of the table and table */
  fprintf(fp, "</tbody></table>\n");

  return(0);
}
/*****************************************************************************/

/*****************************************************************************/
/** Set study number based on filename.
 *  @sa studynr_from_fname()
 *  @return Non-zero in case of an error.
 */
int resFName2study(
  char *fname, 
  char *studyNumber
) {
  return(studynr_from_fname(fname, studyNumber));
}
/*****************************************************************************/

/*****************************************************************************/
/** Calculate the median and the lowest and highest value in the
    specified double array data of length nr.
    Note that array is sorted in this function.
    NULL pointer may be specified to function in place of an unwanted
    return parameter.
\return Returns 0 if succesfull.
 */
int resMedian(
  /** Array of data */
  double *data,
  /** Length of data array */
  int nr,
  /** Pointer where median is written */
  double *median,
  /** Pointer where min is written */
  double *min,
  /** Pointer where max is written */
  double *max
) {
  /* Check the arguments */
  if(data==NULL) return(1);
  if(nr<1) return(2);
  /* Sort data in increasing order */
  qsort(data, nr, sizeof(double), resQSortComp);
  /* Get minimum and maximum */
  if(min!=NULL) *min=data[0];
  if(max!=NULL) *max=data[nr-1];
  /* Calculate median */
  if(median!=NULL) {
    if(nr%2) *median=data[(nr-1)/2];
    else *median=0.5*(data[(nr/2)-1]+data[nr/2]);
  }
  return(0);
}
/*****************************************************************************/

/*****************************************************************************/
/** Calculate the mean and sd in the specified double array data of length nr.
    NULL pointer may be specified to function in place of an unwanted
    return parameter.
\return Returns 0 if succesfull.
 */
int resMean(
  /** Array of data */
  double *data,
  /** Length of data array */
  int nr,
  /** Pointer where mean is written */
  double *mean,
  /** Pointer where S.D. is written */
  double *sd
) {
  int i;
  double sum, ssum, v;

  /* Check the arguments */
  if(data==NULL) return(1);
  if(nr<1) return(2);
  /* Calculate avg and sd */
  for(i=0, sum=ssum=0.0; i<nr; i++) {
    sum+=data[i]; ssum+=data[i]*data[i];
  }
  if(mean!=NULL) *mean=sum/(double)nr;
  if(sd!=NULL) {
    if(nr>1) v=(ssum-sum*sum/(double)nr)/(double)(nr-1); else v=0.0;
    if(v>0.0) *sd=sqrt(v); else *sd=0.0;
  }
  return(0);
}
/*****************************************************************************/

/*****************************************************************************/
/// @cond
int resQSortComp(const void *f1, const void *f2)
{
  if(*(double*)f1<*(double*)f2) return(-1);
  else if(*(double*)f1>*(double*)f2) return(1);
  else return(0);
}
/// @endcond
/*****************************************************************************/

/*****************************************************************************/
/** Sort RES regions by region name. */
void resSortByName(
  /* Pointer to RES struct */
  RES *res
) {
  if(res==NULL || res->voiNr<=1) return;
  qsort(res->voi, res->voiNr, sizeof(ResVOI), resQSortName);
  return;
}
/*****************************************************************************/

/*****************************************************************************/
/// @cond
int resQSortName(const void *voi1, const void *voi2)
{
  int res;

  res=strcasecmp(
       ((const ResVOI*)voi1)->voiname, ((const ResVOI*)voi2)->voiname );
  if(res!=0) return(res);
  res=strcasecmp(
       ((const ResVOI*)voi1)->hemisphere, ((const ResVOI*)voi2)->hemisphere );
  if(res!=0) return(res);
  res=strcasecmp( ((const ResVOI*)voi1)->place, ((const ResVOI*)voi2)->place );
  return(res);
}
/// @endcond
/*****************************************************************************/

/*****************************************************************************/
/** Copy result main header information to another result structure. */
int resCopyMHeader(
  /* Pointer to RES struct to copy header from */
  RES *res1,
  /* Pointer to RES struct to copy header to */
  RES *res2
) {
  int i;

  if(res1==NULL || res2==NULL) return(1);
  strcpy(res2->program, res1->program);
  res2->time=res1->time;
  res2->parNr=res1->parNr;
  strcpy(res2->studynr, res1->studynr);
  strcpy(res2->datafile, res1->datafile);
  strcpy(res2->reffile, res1->reffile);
  strcpy(res2->plasmafile, res1->plasmafile);
  strcpy(res2->plasmafile2, res1->plasmafile2);
  strcpy(res2->bloodfile, res1->bloodfile);
  strcpy(res2->refroi, res1->refroi);
  strcpy(res2->datarange, res1->datarange);
  res2->datanr=res1->datanr;
  strcpy(res2->fitmethod, res1->fitmethod);
  res2->density=res1->density;
  res2->lc=res1->lc;
  res2->concentration=res1->concentration;
  res2->beta=res1->beta;
  res2->Vb=res1->Vb;
  res2->fA=res1->fA;
  res2->E=res1->E;
  res2->isweight=res1->isweight;
  for(i=0; i<res1->parNr; i++) {
    strcpy(res2->parname[i], res1->parname[i]);
    strcpy(res2->parunit[i], res1->parunit[i]);
  }
  strcpy(res2->titleline, res1->titleline);
  strcpy(res2->unitline, res1->unitline);
  return(0);
}
/*****************************************************************************/

/*****************************************************************************/
/** Delete specified region (0..voiNr-1) from the structure.
\return Returns 0 if ok.
 */
int resDelete(
  /** Pointer to the result data */
  RES *res,
  /** TAC index (0..voiNr-1) */
  int voi
) {
  int i;

  /* Check that region exists */
  if(res==NULL || voi>res->voiNr-1 || voi<0) return(1);
  /* If it is the last one, then just decrease the voiNr */
  if(voi==res->voiNr) {res->voiNr--; return(0);}
  /* Otherwise we have to move the following regions in its place */
  for(i=voi+1; i<res->voiNr; i++)
    memcpy(&res->voi[i-1], &res->voi[i], sizeof(ResVOI));
  res->voiNr--;
  return(0);
}
/*****************************************************************************/

/*****************************************************************************/
/** Select VOIs (sets sw=1), whose names are matching specified string.
    If no string is specified, then all VOIs are selected.
\return Returns the number of matches, or <0, if an error occurred.
 */
int resSelect(
  /** Pointer to the result data */
  RES *data,
  /** Region name string */
  char *name
) {
  unsigned int i, j, n;
  char *p, n1[128], n2[128], n3[128], tmp[128], sname[1024], *lptr;

  if(data==NULL) return(-1);
  /* Select all, if no string was specified */
  if(name==NULL || strlen(name)==0) {
    for(i=0; i<(unsigned int)data->voiNr; i++) data->voi[i].sw=1; 
    return(data->voiNr);
  }
  /* Make a copy of 'name' and use it */
  strcpy(sname, name); lptr=sname;
  /* Check if string contains several substrings (hemisphere and place) */
  n1[0]=n2[0]=n3[0]=(char)0;
  p=strtok(lptr, " ,;\n\t|"); if(p!=NULL) strcpy(n1, p); else return(-1);
  p=strtok(NULL, " ,;\n\t|"); if(p!=NULL) {
    strcpy(n2, p); p=strtok(NULL, " ,;\n\t|"); if(p!=NULL) strcpy(n3, p);}
  /* Convert strings to lowercase */
  for(i=0; i<strlen(n1); i++) n1[i]=tolower(n1[i]);
  for(i=0; i<strlen(n2); i++) n2[i]=tolower(n2[i]);
  for(i=0; i<strlen(n3); i++) n3[i]=tolower(n3[i]);
  /* Search through the data */
  for(i=0, n=0; i<(unsigned int)data->voiNr; i++) {
    data->voi[i].sw=0;
    snprintf(tmp, 128, "%s%s%s", data->voi[i].voiname, data->voi[i].hemisphere,
                           data->voi[i].place);
    for(j=0; j<strlen(tmp); j++) tmp[j]=tolower(tmp[j]);
    if(strstr(tmp, n1)==NULL) continue;
    if(n2[0] && strstr(tmp, n2)==NULL) continue;
    if(n3[0] && strstr(tmp, n3)==NULL) continue;
    data->voi[i].sw=1; n++;
  }
  return(n);
}
/*****************************************************************************/

/*****************************************************************************/
/** Select the VOIs that have matching region name or number.
    Sets sw=1 or sw=0. This will replace resSelect().
\return Returns the number of selected VOIs, or <0 in case of an error.
 */
int resSelectRegions(
  /** Pointer to RES data where VOIs are selected */
  RES *res,
  /** Name or VOI number which is searched */
  char *region_name,
  /** 1=Non-matching VOIs are deselected, 0=Old selections are preserved */
  int reset
) {
  int ri, match_nr=0;

  /* Check the input */
  if(res==NULL || res->voiNr<1 || strlen(region_name)<1) return(-1);
  /* Reset all selections if required */
  if(reset!=0) for(ri=0; ri<res->voiNr; ri++) res->voi[ri].sw=0;
  /* Check each VOI */
  for(ri=0; ri<res->voiNr; ri++) {
    if(rnameMatch(res->voi[ri].name, ri+1, region_name)!=0) {
      res->voi[ri].sw=1; match_nr++;
    }
  }
  return(match_nr);
}
/*****************************************************************************/

/*****************************************************************************/
/** Determine whether the result parameter should be printed as
    integer (0), float (1), or exponential (2).
\return Returns the code, or <0 in case of an error.
 */
int resParameterPrintType(
  /** Pointer to the result struct */
  RES *res,
  /** Index of the parameter to test */
  int parIndex
) {
  int vi, partype=0;
  double x, m=0.0, pint;

  if(res==NULL || res->voiNr<1 || parIndex<0 || parIndex>=res->parNr)
    return(-1);
  for(vi=0; vi<res->voiNr; vi++) {
    x=res->voi[vi].parameter[parIndex];
    if(modf(x, &pint)!=0.0) partype=1;
    x=fabs(x); if(x>m) m=x;
  }
  if(partype==1 && (m>=10.0 || m<0.1)) partype=2;
  return(partype);
}
/*****************************************************************************/

/*****************************************************************************/
/** Check if result structure contains duplicate region names.
\return Returns 0 if not, 1 if duplicates are found, and <0 in case of an error.
 */  
int resIsDuplicateNames(
  /** Pointer to the result data */
  RES *res
) {
  int ri, rj;
  
  if(res==NULL) return(-1);
  if(res->voiNr<2) return(0);
  for(ri=0; ri<res->voiNr-1; ri++) for(rj=ri+1; rj<res->voiNr; rj++)
    if(strcasecmp(res->voi[ri].name, res->voi[rj].name)==0) return(1);
  return(0);
}
/*****************************************************************************/

/*****************************************************************************/
/** Check whether result header field values are the same.
    Fields that are not checked: program, time, titleline.
\return Returns 0 in case of match, and <>0 if not matching.
 */
int resMatchHeader(
  /** Pointers to the result data that are tested */
  RES *res1,
  /** Pointers to the result data that are tested */
  RES *res2
) {
  if(res1==NULL || res2==NULL) return(1);
  if(res1->voiNr!=res2->voiNr) return(3);
  if(res1->parNr!=res2->parNr) return(4);
  if(strcasecmp(res1->datafile, res2->datafile)!=0) return(6);
  if(strcasecmp(res1->reffile, res2->reffile)!=0) return(7);
  if(strcasecmp(res1->plasmafile, res2->plasmafile)!=0) return(8);
  if(strcasecmp(res1->plasmafile2, res2->plasmafile2)!=0) return(9);
  if(strcasecmp(res1->bloodfile, res2->bloodfile)!=0) return(10);
  if(strcasecmp(res1->refroi, res2->refroi)!=0) return(11);
  if(strcasecmp(res1->datarange, res2->datarange)!=0) return(12);
  if(res1->isweight!=res2->isweight) return(13);
  if(res1->density!=res2->density) return(14);
  if(res1->lc!=res2->lc) return(15);
  if(res1->beta!=res2->beta) return(16);
  if(res1->concentration!=res2->concentration) return(17);
  if(res1->Vb!=res2->Vb) return(18);
  if(res1->datanr!=res2->datanr) return(19);
  if(strcasecmp(res1->fitmethod, res2->fitmethod)!=0) return(20);
  if(res1->fA!=res2->fA) return(21);
  if(res1->E!=res2->E) return(22);
  /* Less important */
  if(strcasecmp(res1->studynr, res2->studynr)!=0) return(5);
  return(0);
}
/*****************************************************************************/

/*****************************************************************************/
/** Check whether result region names are the same.
\return Returns 0 in case of match, and 1 if not matching.
 */
int resMatchRegions(
  /** Pointers to the result data that are tested */
  RES *res1,
  /** Pointers to the result data that are tested */
  RES *res2
) {
  int ri, m=0;

  if(res1==NULL || res2==NULL || res1->voiNr!=res2->voiNr) return(1);
  for(ri=0; ri<res1->voiNr; ri++) {
    m=0;
//  if(strcmp(res1->voi[ri].name, res2->voi[ri].name)!=0) return(1);
    if(strcmp(res1->voi[ri].voiname, res2->voi[ri].voiname)!=0) m++;
    if(strcmp(res1->voi[ri].hemisphere, res2->voi[ri].hemisphere)!=0) m++;;
    if(strcmp(res1->voi[ri].place, res2->voi[ri].place)!=0) m++;
    if(m>0) {
      if(RESULT_TEST>5) {
        printf("  '%s' vs '%s'\n", res1->voi[ri].voiname, res2->voi[ri].voiname);
        printf("  '%s' vs '%s'\n", res1->voi[ri].hemisphere, res2->voi[ri].hemisphere);
        printf("  '%s' vs '%s'\n", res1->voi[ri].place, res2->voi[ri].place);
      }
      return(1);
    }
  }
  return(0);
}
/*****************************************************************************/

/*****************************************************************************/
/** Check whether result parameter names are the same.
\return Returns 0 in case of match, and 1 if not matching.
 */
int resMatchParameternames(
  /** Pointers to the result data that are tested */
  RES *res1,
  /** Pointers to the result data that are tested */
  RES *res2
) {
  int i;

  if(res1==NULL || res2==NULL || res1->parNr!=res2->parNr) return(1);

  resFixParnames(res1);
  resFixParnames(res2);

  for(i=0; i<res1->parNr; i++) {
    if(strcasecmp(res1->parname[i], res2->parname[i])!=0) {
      if(RESULT_TEST>1)
        printf("  Parameter names '%s' and '%s' do not match\n",
               res1->parname[i], res2->parname[i]);
      return(1);
    }
    if(strcasecmp(res1->parunit[i], res2->parunit[i])!=0) {
      if(RESULT_TEST>1)
        printf("  Parameter units '%s' and '%s' do not match\n",
               res1->parunit[i], res2->parunit[i]);
      return(1);
    }
  }
  return(0);
}
/*****************************************************************************/

/*****************************************************************************/
/** Check whether result parameter values are the same.
\return Returns 0 in case of match, and <>0 if not matching.
 */
int resMatchParameters(
  /** Pointers to the result data that are tested */
  RES *res1,
  /** Pointers to the result data that are tested */
  RES *res2,
  /** Parameter index (0..parNr-1) that is verified; <0, if all */
  int test_par,
  /** Test limit (how exact match is required) */
  double test_limit,
  /** Test (1) or do not test (0) SD and Confidence limits */
  int test_sd
) {
  int ri, pi;
  double s, v1, v2;

  if(res1==NULL || res2==NULL || res1->voiNr!=res2->voiNr) return(1);
  if(res1->parNr!=res2->parNr) {
    if(test_par<0) return(1);
    if(test_par+1>res1->parNr || test_par+1>res2->parNr) return(1);
  }
  for(ri=0; ri<res1->voiNr; ri++) {
    for(pi=0; pi<res1->parNr; pi++) {
      if(test_par>=0 && test_par!=pi) continue;
      v1=res1->voi[ri].parameter[pi]; v2=res2->voi[ri].parameter[pi];
      if(isnan(v1) && isnan(v2)) continue; // ok if both are NaN
      /* Parameter values */
      if(RESULT_TEST>5) printf("pi=%d ri=%d\n", pi, ri);
      if(RESULT_TEST>5) printf(" %g vs %g\n", v1, v2);
      s=fabs(v1+v2); if(isnan(s)) return(2); // Either one is NaN
      if(test_limit<=0.0) { // values are requested to match exactly
        if(v1!=v2) return(2);
      } if(s==0.0 || v1==0.0 || v2==0.0) { // relative matching cant be done
        if(fabs(v1-v2)>test_limit) return(2);
      } else { // test relative difference against given limit
        if(fabs((v1-v2)/s)>test_limit) return(2);
      }
      if(test_sd!=0) {
        /* SD */
        v1=res1->voi[ri].sd[pi]; v2=res2->voi[ri].sd[pi];
        if(RESULT_TEST>5) printf(" SD: %g vs %g\n", v1, v2);
        if(isnan(v1) && isnan(v2)) {
          // both are NA, that is ok
        } else if(isnan(v1) || isnan(v2)) {
          // one is NA, that is not ok
          return(3);
        } else if(test_limit<=0.0) { // values are requested to match exactly
          if(v1!=v2) return(3);
        } if(s==0.0) { // relative matching cant be done
          if(fabs(v1-v2)>test_limit) return(3);
        } else { // test relative difference agains given limit
          if(fabs((v1-v2)/s)>test_limit) return(3);
        }
        /* CL1 */
        v1=res1->voi[ri].cl1[pi]; v2=res2->voi[ri].cl1[pi];
        if(RESULT_TEST>5) printf(" CL1: %g vs %g\n", v1, v2);
        if(isnan(v1) && isnan(v2)) {
          // both are NA, that is ok
        } else if(isnan(v1) || isnan(v2)) {
          // one is NA, that is not ok
          return(4);
        } else if(test_limit<=0.0) { // values are requested to match exactly
          if(v1!=v2) return(4);
        } if(s==0.0) { // relative matching cant be done
          if(fabs(v1-v2)>test_limit) return(4);
        } else { // test relative difference agains given limit
          if(fabs((v1-v2)/s)>test_limit) return(4);
        }
        /* CL2 */
        v1=res1->voi[ri].cl2[pi]; v2=res2->voi[ri].cl2[pi];
        if(RESULT_TEST>5) printf(" CL2: %g vs %g\n", v1, v2);
        if(isnan(v1) && isnan(v2)) {
          // both are NA, that is ok
        } else if(isnan(v1) || isnan(v2)) {
          // one is NA, that is not ok
          return(5);
        } else if(test_limit<=0.0) { // values are requested to match exactly
          if(v1!=v2) return(5);
        } if(s==0.0) { // relative matching cant be done
          if(fabs(v1-v2)>test_limit) return(5);
        } else { // test relative difference agains given limit
          if(fabs((v1-v2)/s)>test_limit) return(5);
        }
      }
    }
  }
  return(0);
}
/******************************************************************************/

/******************************************************************************/
/** Check whether the two sets of result parameter values are similar within
 *  a given absolute range.  
\return Returns 0 in case of match, and <>0 if not matching.
 */
int resMatchParametersAbs(
  /** Pointers to the result data that are tested */
  RES *res1,
  /** Pointers to the result data that are tested */
  RES *res2,
  /** Parameter index (0..parNr-1) that is verified; <0, if all */
  int test_par,
  /** Test limit; positive value, below which the absolute difference must be.*/
  double test_limit,
  /** Test (1) or do not test (0) SD and Condifence limits */
  int test_sd
) {
  int ri, pi;
  double s, v1, v2;

  if(res1==NULL || res2==NULL || res1->voiNr!=res2->voiNr) return(1);
  if(res1->parNr!=res2->parNr) {
    if(test_par<0) return(1);
    if(test_par+1>res1->parNr || test_par+1>res2->parNr) return(1);
  }
  if(test_limit<0.0) return(1);
  for(ri=0; ri<res1->voiNr; ri++) {
    for(pi=0; pi<res1->parNr; pi++) {
      if(test_par>=0 && test_par!=pi) continue;
      v1=res1->voi[ri].parameter[pi]; v2=res2->voi[ri].parameter[pi];
      if(isnan(v1) && isnan(v2)) continue;
      /* Parameter values */
      if(RESULT_TEST>5) printf("pi=%d ri=%d\n", pi, ri);
      if(RESULT_TEST>5) printf(" %g vs %g\n", v1, v2);
      s=fabs(v1-v2);
      if(isnan(s) || s>test_limit) return(2);
      if(test_sd!=0) {
        /* SD */
        v1=res1->voi[ri].sd[pi]; v2=res2->voi[ri].sd[pi];
        if(isnan(v1) && !isnan(v2)) return(3);
        if(!isnan(v1) && isnan(v2)) return(3);
        if(!isnan(v1) && !isnan(v2)) {
          if(RESULT_TEST>8) printf(" SD: %g vs %g\n", v1, v2);
          s=fabs(v1-v2); if(s>test_limit) return(3);
        }
        /* CL1 */
        v1=res1->voi[ri].cl1[pi]; v2=res2->voi[ri].cl1[pi];
        if(isnan(v1) && !isnan(v2)) return(4);
        if(!isnan(v1) && isnan(v2)) return(4);
        if(!isnan(v1) && !isnan(v2)) {
          if(RESULT_TEST>8) printf(" CL1: %g vs %g\n", v1, v2);
          s=fabs(v1-v2); if(s>test_limit) return(4);
        }
        /* CL2 */
        v1=res1->voi[ri].cl2[pi]; v2=res2->voi[ri].cl2[pi];
        if(isnan(v1) && !isnan(v2)) return(5);
        if(!isnan(v1) && isnan(v2)) return(5);
        if(!isnan(v1) && !isnan(v2)) {
          if(RESULT_TEST>8) printf(" CL2: %g vs %g\n", v1, v2);
          s=fabs(v1-v2); if(s>test_limit) return(5);
        }
      }
    }
  }
  return(0);
}
/*****************************************************************************/

/*****************************************************************************/
/** Check whether region name subfields exist in any region.
\return Returns 1 if hemisphere exists, 2 if place exists, 3 if both exist,
    and 0 if neither exists, and <0 in case of an error. 
 */
int resRNameSubfieldExists(
  /** Pointer to RES struct */
  RES *res
) {
  int ri, m=0, n=0;

  if(res==NULL) return(-1);
  if(res->voiNr<1) return(-1);
  for(ri=0; ri<res->voiNr; ri++) {
    if(strlen(res->voi[ri].hemisphere)>0 &&
       strcmp(res->voi[ri].hemisphere, ".")!=0)
      m++;
    if(strlen(res->voi[ri].place)>0 &&
       strcmp(res->voi[ri].place, ".")!=0)
      n++;
  }
  ri=0; if(m>0) ri+=1; if(n>0) ri+=2;
  return(ri);
}
/*****************************************************************************/

/*****************************************************************************/
