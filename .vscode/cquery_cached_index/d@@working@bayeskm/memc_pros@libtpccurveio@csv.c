/// @file csv.c
/// @author Vesa Oikonen
/// @brief I/O functions for CSV files (comma-separated values).
///
/*****************************************************************************/
#include "libtpccurveio.h"
/*****************************************************************************/

/*****************************************************************************/
/** Initiate CSV struct */
void csvInit(
  /** Pointer to CSV struct */
  CSV *csv
) {
  if(csv==NULL) return;
  csv->c=NULL;
  csv->nr=csv->row_nr=csv->col_nr=0;
  csv->separator=(char)0;
}

/** Free allocated memory in CSV struct */
void csvEmpty(
  /** Pointer to CSV struct */
  CSV *csv
) {
  int i;
  if(csv==NULL) return;
  for(i=0; i<csv->nr; i++)
    if(csv->c[i].content!=NULL) free(csv->c[i].content);
  free(csv->c); csv->c=NULL;
  csv->nr=csv->row_nr=csv->col_nr=0;
  csv->separator=(char)0;
}
/*****************************************************************************/

/*****************************************************************************/
/** Read CSV file.
 *  @return Non-zero in case of an error.
 */
int csvRead(
  /** Pointer to CSV struct */
  CSV *csv,
  /** Filename */
  char *fname
) {
  if(CSV_TEST>2) {printf("csvRead('%s')\n", fname); fflush(stdout);}

  FILE *fp;
  int i, nr, ret, nonprintable=0, inside_quotes=0, previous, col_nr=0;
  const int MAX_CSV_FIELD_LENGTH=1024;
  char buf[MAX_CSV_FIELD_LENGTH+1];
  int tabnr, spacenr, commanr;
  
  
  if(csv==NULL || fname==NULL) return CSV_ERROR;
  /* Open file */
  fp=fopen(fname, "r"); if(fp==NULL) return CSV_CANNOTOPEN;

  /* Check the file size */
  nr=nonprintable=0; while((ret=fgetc(fp))!=EOF) {
    if(iscntrl(ret) && ret!=13 && ret!=10 && ret!=9) {nonprintable=1; break;}
    nr++;
  }
  if(CSV_TEST>0) printf("filesize := %d\n", nr);
  if(nr<2) {fclose(fp); return CSV_INVALIDFORMAT;}
  if(nr>5000000) {fclose(fp); return CSV_TOOBIG;}
  rewind(fp);

  /* Determine the field separator (unless set outside) */
  if(csv->separator==(char)0) {
    /* Check if ; or tab or space or comma character is found outside double quotes */
    inside_quotes=0; nr=0; tabnr=0; spacenr=0; commanr=0;
    while((ret=fgetc(fp))!=EOF) {
      if(ret=='"') {
        if(inside_quotes==0) inside_quotes=1; else inside_quotes=0;
	continue;
      }
      if(inside_quotes==1) continue;
      if(ret==';') nr++;
      else if(ret=='\t') tabnr++;
      else if(ret==',') commanr++;
      else if(ret==' ') spacenr++;
    }
    if(CSV_TEST>0) {
      printf("semicolon_nr := %d\n", nr);
      printf("tab_nr := %d\n", tabnr);
      printf("comma_nr := %d\n", commanr);
      printf("space_nr := %d\n", spacenr);
    }
    /* If at least one, then assume that ; is the separator, otherwise , or tab */
    if(nr>0) csv->separator=';';
    else if(tabnr>0) csv->separator='\t'; 
    else if(commanr>spacenr) csv->separator=',';
    else csv->separator=' ';
    rewind(fp);
  }
  if(CSV_TEST>0) printf("separator := '%c'\n", csv->separator);
  /* We will not accept space as separator for CSV */
  if(csv->separator==' ') {fclose(fp); return CSV_INVALIDFORMAT;}

  /* Determine the number of fields in CSV file */
  inside_quotes=0; nr=0; previous=0;
  while((ret=fgetc(fp))!=EOF) {
    if((previous==13 || previous==10) && (ret==13 || ret==10)) {previous=ret; continue;}
    if(ret=='"') {
      if(inside_quotes==0) inside_quotes=1; else inside_quotes=0;
      previous=ret; continue;
    }
    if(inside_quotes==0) {
      if(ret==csv->separator) {
        nr++; previous=ret;
	continue;
      }
      if( (ret==13 || ret==10) && previous!=13 && previous!=10) {
        nr++; previous=ret;
	continue;
      }
    } //printf("%c", (char)ret);
    previous=ret;
  }
  rewind(fp); if(CSV_TEST>0) printf("field_nr := %d\n", nr);

  /* Allocate memory for fields */
  csv->c=(CSV_item*)calloc(nr, sizeof(CSV_item));
  if(csv->c==NULL) {fclose(fp); return CSV_OUTOFMEMORY;}
  csv->nr=nr;

  /* Copy field contents from CSV file */
  if(CSV_TEST>0) printf("  copying contents...\n");
  inside_quotes=0; nr=0; previous=0; i=0; col_nr=0;
  while((ret=fgetc(fp))!=EOF) {
    if((previous==13 || previous==10) && (ret==13 || ret==10)) {previous=ret; continue;}
    if(ret=='"') {
      if(inside_quotes==0) inside_quotes=1; else inside_quotes=0;
      previous=ret; continue;
    }
    if(inside_quotes==0) {
      if(ret==csv->separator) {
        buf[i]=(char)0; strCleanSpaces(buf); if(CSV_TEST>10) printf("'%s'\n", buf);
        if(nr>=csv->nr) {printf("  index overflow\n"); break;}
        csv->c[nr].content=strdup(buf);
/*
        if(i>0) {
          if(nr>=csv->nr) {printf("  index overflow\n"); break;}
          csv->c[nr].content=(char*)malloc(i+1);
          if(csv->c[nr].content!=NULL) strcpy(csv->c[nr].content, buf);
        } //if(CSV_TEST>10) printf("---\n");
*/
        csv->c[nr].row=1+csv->row_nr; csv->c[nr].col=1+col_nr;
	i=0; nr++; col_nr++; previous=ret;
	continue;
      }
      if( (ret==13 || ret==10) && previous!=13 && previous!=10) {
        buf[i]=(char)0; strCleanSpaces(buf); if(CSV_TEST>10) printf("'%s'\n", buf); 
        if(nr>=csv->nr) {printf("   index overflow\n"); break;}
        csv->c[nr].content=strdup(buf);
        if(CSV_TEST>10) printf("===\n");
/*
        if(i>0) {
          if(nr>=csv->nr) {printf("   index overflow\n"); break;}
          csv->c[nr].content=(char*)malloc(i+1);
          if(csv->c[nr].content!=NULL) strcpy(csv->c[nr].content, buf);
        } if(CSV_TEST>10) printf("===\n");
*/
        col_nr++; if(col_nr>csv->col_nr) csv->col_nr=col_nr;
        csv->c[nr].row=1+csv->row_nr; csv->c[nr].col=col_nr;
        i=0; nr++; col_nr=0; previous=ret; csv->row_nr++;
	continue;
      }
    }
    if(i<MAX_CSV_FIELD_LENGTH) buf[i]=(char)ret;
    i++;
    previous=ret;
  }
  if(CSV_TEST>0) printf("  ... copied: nr=%d\n", nr);

  fclose(fp);
  return CSV_OK;
}
/*****************************************************************************/

/*****************************************************************************/
/** Reads different CSV formats into DFT struct
\return Returns 0 when successful, otherwise an error code.
 */
int csv2dft(
  /** Pointer to CSV data to be converted */
  CSV *csv,
  /** Pointer to empty DFT struct which will be allocated and filled here */
  DFT *dft
) {
  int ret, i, u, n, m;

  if(CSV_TEST>2) {printf("csv2dft()\n"); fflush(stdout);}
  if(csv==NULL || dft==NULL) return CSV_ERROR;
  if(csv->row_nr<1 || csv->col_nr<1) return CSV_INVALIDFORMAT;

  /* Is this LinkSet format? */
  if(strcasecmp(csv->c[0].content, "LinkSet")==0) {
    ret=csv2dft_linkset(csv, dft);
    if(ret!=CSV_OK) {
      if(CSV_TEST>2) printf("reading LinkSet CSV format failed.\n");
    }
    return(ret);
  }

  /* Maybe mat file? */
  ret=csv2dft_mat(csv, dft);
  if(ret==CSV_OK) {
    if(CSV_TEST>2) printf("reading Mat CSV format successful.\n");
    return(ret);
  }

  /* Many CSV formats are impossible to identify, therefore we will just
     try to convert different formats until one succeeds or all are failed */
  if(CSV_TEST>2) printf("trying to read 1st CSV format\n");
  ret=csv2dft_a(csv, dft);
  if(ret!=CSV_OK) {
    if(CSV_TEST>2) printf("reading 1st CSV format failed; trying 2nd format\n");
    ret=csv2dft_b(csv, dft);
  }
  if(ret!=CSV_OK) {
    if(CSV_TEST>2) printf("2nd CSV format failed\n");
  }
  if(ret!=CSV_OK) {dft->_type=DFT_FORMAT_PLAIN; return ret;}
  /* Make sure that TAC names are filled */
  u=dft->voiNr; n=1; while((u/=10)>=1) n++;
  if(n>MAX_REGIONSUBNAME_LEN) n=MAX_REGIONSUBNAME_LEN;
  for(i=0, m=0; i<dft->voiNr; i++) {
    if(strlen(dft->voi[i].voiname)<1 || strcmp(dft->voi[i].voiname, ".")==0) {
      snprintf(dft->voi[i].voiname, 7, "%0*d", n, i+1);
      strcpy(dft->voi[i].name, dft->voi[i].voiname);
      m++;
    }
  }
  /* If none of TACs had a name, then set DFT plain format */
  if(m==dft->voiNr) dft->_type=DFT_FORMAT_PLAIN;

  if(CSV_TEST>3) dftPrint(dft);
  return CSV_OK;
}
/*****************************************************************************/

/*****************************************************************************/
/** Reads simple and Inveon type 1 data into DFT struct
\return Returns 0 when successful, otherwise an error code.
 */
int csv2dft_a(
  /** Pointer to CSV data to be converted */
  CSV *csv,
  /** Pointer to empty DFT struct which will be allocated and filled here */
  DFT *dft
) {
  int ri, ci, sci, ii, ret;
  char *cptr;

  if(CSV_TEST>2) {printf("csv2dft_a()\n"); fflush(stdout);}
  if(csv==NULL || dft==NULL) return CSV_ERROR;
  if(csv->row_nr<1 || csv->col_nr<1) return CSV_INVALIDFORMAT;


  if(CSV_TEST>2)
    for(int i=0; i<csv->nr; i++)
      printf("row=%d col=%d content='%s'\n", csv->c[i].row, csv->c[i].col, csv->c[i].content);


  /* Allocate memory for DFT */
  dftEmpty(dft);
  if(CSV_TEST>2) {
    printf("frame_nr=%d voi_nr=%d\n", csv->row_nr, csv->col_nr-1); 
    fflush(stdout);
  }
  ret=dftSetmem(dft, csv->row_nr, csv->col_nr-1);
  if(ret!=0) return CSV_OUTOFMEMORY;
  /* Set DFT defaults */
  dft->timetype=DFT_TIME_MIDDLE;
  dft->_type=DFT_FORMAT_STANDARD;
  dft->isweight=0;
  dftUnitToDFT(dft, CUNIT_UNKNOWN);
  dft->timeunit=TUNIT_UNKNOWN;
  for(ri=0; ri<csv->row_nr; ri++) dft->w[ri]=1.0;
  for(ci=0; ci<csv->col_nr-1; ci++) dft->voi[ci].sw=0;

  /* Fill DFT */
  ri=0;
  for(ii=0; ii<csv->nr;) {
    // goto start of row
    for(; ii<csv->nr && csv->c[ii].col!=1; ii++) {}
    if(ii==csv->nr) break;
    if(CSV_TEST>10) {
      printf("\nline start at %d\n", ii); 
      printf("  ri=%d\n", ri); 
      fflush(stdout);
    }
    // ignore line with empty first column
    if(csv->c[ii].content==NULL) {
      if(CSV_TEST>11) {printf("  empty first column\n"); fflush(stdout);}
      ii++; continue;
    }
    // ignore comment line
    if(csv->c[ii].content[0]=='#') {
      if(CSV_TEST>11) {printf("  comment line\n"); fflush(stdout);}
      ii++; continue;
    }
    // ignore line that does not start with number
    // unless it contains TAC titles
    if(!isdigit(csv->c[ii].content[0]) && csv->c[ii].content[0]!='-') {
      if(strstr(csv->c[ii].content, "Time")==NULL &&
         strstr(csv->c[ii].content, "TIME")==NULL &&
         strstr(csv->c[ii].content, "time")==NULL) {
        if(CSV_TEST>11) {
          printf("  not a numerical value or title\n"); fflush(stdout);}
        ii++; continue;
      }
      dft->_type=DFT_FORMAT_STANDARD;
      if(strncasecmp(csv->c[ii].content, "Start time", 10)==0 &&
         strncasecmp(csv->c[ii+1].content, "End time", 8)==0) {
        dft->timetype=DFT_TIME_STARTEND;
        if(CSV_TEST>6) printf("timetype := %d\n", dft->timetype);
      }
      if(CSV_TEST>7) {
        printf("first title field := '%s'\n", csv->c[ii].content);
        fflush(stdout);
      }
      if(strstr(csv->c[ii].content, "min")!=NULL) dft->timeunit=TUNIT_MIN;
      else if(strstr(csv->c[ii].content, "sec")!=NULL) dft->timeunit=TUNIT_SEC;
      else dft->timeunit=TUNIT_UNKNOWN;
      ii++;
      
      if(dft->timetype==DFT_TIME_MIDDLE) sci=2; else {sci=3; ii++;}
      for(ci=sci; ci<=csv->col_nr && ii<csv->nr; ci++, ii++) {
        if(CSV_TEST>2) {
          printf("col=%d row=%d\n", csv->c[ii].col, csv->c[ii].row);
          if(CSV_TEST>3) printf("ci=%d ii=%d\n", ci, ii);
          fflush(stdout);
        }
        if(csv->c[ii].col!=ci) {dftEmpty(dft); return CSV_NOTABLE;}
        if(csv->c[ii].content!=NULL) {
          /* Check if this column should be omitted from DFT */
          if(strstr(csv->c[ii].content, " - Time")!=NULL) {
            if(CSV_TEST>2) printf("  ignored time column.\n");
            dft->voi[ci-sci].sw=1; continue;
          }
          if(strstr(csv->c[ii].content, "(upper bound)")!=NULL) {
            if(CSV_TEST>2) printf("  ignored upper bound column.\n");
            dft->voi[ci-sci].sw=2; continue;
          }
          if(strstr(csv->c[ii].content, "(lower bound)")!=NULL) {
            if(CSV_TEST>2) printf("  ignored lower bound column.\n");
            dft->voi[ci-sci].sw=3; continue;
          }
          if(strstr(csv->c[ii].content, "(standard deviation)")!=NULL) {
            if(CSV_TEST>2) printf("  ignored s.d. column.\n");
            dft->voi[ci-sci].sw=4; continue;
          }
          /* Search calibration unit from name */
          if(petCunitId(dft->unit)==CUNIT_UNKNOWN) {
            if(strstr(csv->c[ii].content, "(Bq/ml)")!=NULL)
              dftUnitToDFT(dft, CUNIT_BQ_PER_ML);
            else if(strstr(csv->c[ii].content, "(kBq/ml)")!=NULL)
              dftUnitToDFT(dft, CUNIT_KBQ_PER_ML);
            else if(strstr(csv->c[ii].content, "(MBq/ml)")!=NULL)
              dftUnitToDFT(dft, CUNIT_MBQ_PER_ML);
            else if(strstr(csv->c[ii].content, "(% ID/g)")!=NULL)
              dftUnitToDFT(dft, CUNIT_PIDM);
            else if(strstr(csv->c[ii].content, "Bq/ml")!=NULL)
              dftUnitToDFT(dft, CUNIT_BQ_PER_ML);
            else if(strstr(csv->c[ii].content, "kBq/ml")!=NULL)
              dftUnitToDFT(dft, CUNIT_KBQ_PER_ML);
            else if(strstr(csv->c[ii].content, "MBq/ml")!=NULL)
              dftUnitToDFT(dft, CUNIT_MBQ_PER_ML);
            else if(strstr(csv->c[ii].content, "% ID/g")!=NULL)
              dftUnitToDFT(dft, CUNIT_PIDM);
          }
        }
        if(csv->c[ii].content==NULL) {
          sprintf(dft->voi[ci-sci].name, "%d", ci-sci+1);
          strcpy(dft->voi[ci-sci].voiname, dft->voi[ci-sci].name);
        } else {
          cptr=strstr(csv->c[ii].content, " - "); if(cptr!=NULL) *cptr=(char)0;
          strncpy(dft->voi[ci-sci].name, csv->c[ii].content, MAX_REGIONNAME_LEN);
          dft->voi[ci-sci].name[MAX_REGIONNAME_LEN]=(char)0;
          rnameSplit(csv->c[ii].content, dft->voi[ci-sci].voiname,
  	  dft->voi[ci-sci].hemisphere, dft->voi[ci-sci].place,
          MAX_REGIONSUBNAME_LEN);
        }
        if(CSV_TEST>8) printf("name[%d]=%s\n", ci-sci, dft->voi[ci-sci].name);
      }
    }
    
    // check that allocated DFT frame nr is not exceeded
    if(ri>=csv->row_nr) {dftEmpty(dft); return CSV_INVALIDFORMAT;}
    // read the sample time
    if(dft->timetype==DFT_TIME_MIDDLE) {
      dft->x[ri]=atof_dpi(csv->c[ii++].content);
      if(CSV_TEST>3) printf("x[%d]=%g\n", ri, dft->x[ri]);
    } else {
      dft->x1[ri]=atof_dpi(csv->c[ii++].content);
      dft->x2[ri]=atof_dpi(csv->c[ii++].content);
      dft->x[ri]=0.5*(dft->x1[ri]+dft->x2[ri]);
      if(CSV_TEST>3)
        printf("x1[%d]=%g x2[%d]=%g\n", ri, dft->x1[ri], ri, dft->x2[ri]);
    }
    // read the sample values
    if(dft->timetype==DFT_TIME_MIDDLE) sci=2; else sci=3;
    for(ci=sci; ci<=csv->col_nr && ii<csv->nr; ci++, ii++) {
      if(CSV_TEST>2) {
        printf("  col=%d row=%d\n", csv->c[ii].col, csv->c[ii].row);
        if(CSV_TEST>3) printf("  ci=%d ii=%d\n", ci, ii);
        fflush(stdout);
      }
      if(csv->c[ii].col!=ci) {dftEmpty(dft); return CSV_NOTABLE;}
      if(strlen(csv->c[ii].content)==0 || strcmp(csv->c[ii].content, ".")==0)
        dft->voi[ci-sci].y[ri]=nan("");
      else
       dft->voi[ci-sci].y[ri]=atof_dpi(csv->c[ii].content);
      if(CSV_TEST>4)
        printf("  y[%d][%d]=%g\n", ri, ci-sci, dft->voi[ci-sci].y[ri]);
    }
    ri++; //if(ri>csv->row_nr) {dftEmpty(dft); return CSV_INVALIDFORMAT;}
  }
  if(CSV_TEST>1) printf("  %d frame(s) read from CSV\n", ri);
  if(ri<1) {dftEmpty(dft); return CSV_INVALIDFORMAT;}
  dft->frameNr=ri;
  dft->voiNr=csv->col_nr-1; if(dft->timetype==DFT_TIME_STARTEND) dft->voiNr--;

  /* Remove those DFT VOIs which where above set to be deleted (where sw!=0) */
  for(ci=dft->voiNr-1, ret=0; ci>=0; ci--) if(dft->voi[ci].sw!=0) {
    ret=dftDelete(dft, ci); if(ret!=0) break;}
  if(ret!=0) {dftEmpty(dft); return CSV_INVALIDFORMAT;}
  if(dft->voiNr<1) {dftEmpty(dft); return CSV_INVALIDFORMAT;}

  /* Calculate frame mid times, or frame start and end times;
     this works only if time units are known */
  dftFrametimes(dft);

  return CSV_OK;
}
/*****************************************************************************/

/*****************************************************************************/
/** Reads Inveon type 2 data into DFT struct
\return Returns 0 when successful, otherwise an error code.
 */
int csv2dft_b(
  /** Pointer to CSV data to be converted */
  CSV *csv,
  /** Pointer to empty DFT struct which will be allocated and filled here */
  DFT *dft
) {
  int ri, fi, fip, ii, ret;
  char *cptr, *cptr2, tmp[256];
  double v1, v2;

  if(CSV_TEST>2) {printf("csv2dft_b()\n"); fflush(stdout);}
  if(csv==NULL || dft==NULL) return CSV_ERROR;
  //printf("row_nr=%d col_nr=%d\n", csv->row_nr, csv->col_nr);
  if(csv->row_nr<4 || csv->col_nr!=9) return CSV_INVALIDFORMAT;
  dftEmpty(dft);

  /* Check the format; first line (containing titles) */
  if(strcasecmp(csv->c[0].content, "#Subject ID")!=0) return CSV_INVALIDFORMAT;
  if(strcasecmp(csv->c[1].content, "Subject Weight")!=0) return CSV_INVALIDFORMAT;
  if(strcasecmp(csv->c[2].content, "Subject Sex")!=0) return CSV_INVALIDFORMAT;
  if(strcasecmp(csv->c[3].content, "Unique Series ID")!=0) return CSV_INVALIDFORMAT;
  if(strcasecmp(csv->c[4].content, "Series Date")!=0) return CSV_INVALIDFORMAT;
  if(strcasecmp(csv->c[5].content, "Series Description")!=0) return CSV_INVALIDFORMAT;

  /* Check the format; third line (containing titles) */
  if(strcasecmp(csv->c[12].content, "#Name")!=0) return CSV_INVALIDFORMAT;
  if(strcasecmp(csv->c[13].content, "Volume (mm^3)")!=0) return CSV_INVALIDFORMAT;
  if(strcasecmp(csv->c[14].content, "Mean")!=0) return CSV_INVALIDFORMAT;
  if(strcasecmp(csv->c[15].content, "SD")!=0) return CSV_INVALIDFORMAT;
  if(strcasecmp(csv->c[16].content, "Min")!=0) return CSV_INVALIDFORMAT;
  if(strcasecmp(csv->c[17].content, "Max")!=0) return CSV_INVALIDFORMAT;
  if(strcasecmp(csv->c[18].content, "Frame Index")!=0) return CSV_INVALIDFORMAT;
  if(strncasecmp(csv->c[19].content, "Mid time (sec)", 10)!=0) return CSV_INVALIDFORMAT;
  if(strncasecmp(csv->c[20].content, "Duration (sec)", 10)!=0) return CSV_INVALIDFORMAT;

  /* Calculate the number of ROIs and time frames */
  ri=1; fi=0; fip=-1; ii=21; cptr=csv->c[ii].content;
  //printf("cell[%d] := '%s'\n", ii, cptr); fflush(stdout);
  for(; ii<csv->nr; ii+=9) {
    cptr2=csv->c[ii].content;
    //printf("cell[%d] := '%s'\n", ii, cptr2); fflush(stdout);
    if(strcmp(cptr, cptr2)==0) fi++;
    else {
      ri++; cptr=cptr2;
      if(fip<0) fip=fi; else if(fi!=fip) return CSV_INVALIDFORMAT;
      fi=1;
    }
    //printf("ri=%d fi=%d fip=%d\n", ri, fi, fip);
  }
  //printf("ri=%d fi=%d\n", ri, fi);

  /* Allocate memory for DFT */
  if(CSV_TEST>2) {printf("frame_nr=%d voi_nr=%d\n", fi, ri); fflush(stdout);}
  ret=dftSetmem(dft, fi, ri);
  if(ret!=0) return CSV_OUTOFMEMORY;
  dft->voiNr=ri; dft->frameNr=fi;
  dft->timetype=DFT_TIME_STARTEND;
  dft->_type=DFT_FORMAT_STANDARD;
  dft->isweight=0;
  dftUnitToDFT(dft, CUNIT_UNKNOWN); dft->timeunit=TUNIT_UNKNOWN;
  for(fi=0; fi<dft->frameNr; fi++) dft->w[fi]=1.0;

  /* Fill DFT */
  /* time unit */
  ii=19;
  if(strstr(csv->c[ii].content, "min")!=NULL) dft->timeunit=TUNIT_MIN;
  else if(strstr(csv->c[ii].content, "sec")!=NULL) dft->timeunit=TUNIT_SEC;
  else dft->timeunit=TUNIT_UNKNOWN;
  /* study number */
  ii=6;
  cptr=csv->c[ii].content; strlcpy(dft->studynr, cptr, MAX_STUDYNR_LEN);
  cptr=strchr(dft->studynr, '.'); if(cptr!=NULL) *cptr=(char)0;
  cptr=strchr(dft->studynr, ','); if(cptr!=NULL) *cptr=(char)0;
  cptr=strchr(dft->studynr, ' '); if(cptr!=NULL) *cptr=(char)0;
  /* subject weight */
  ii=7; if(ii>=csv->nr) {dftEmpty(dft); return CSV_INVALIDFORMAT;}
  v1=atof_dpi(csv->c[ii].content);
  if(v1>0.0) sprintf(dft->comments, "# weight := %g\n", v1);
  /* scan start time */
  ii=10; if(ii>=csv->nr) {dftEmpty(dft); return CSV_INVALIDFORMAT;}
  if(strlen(csv->c[ii].content)>9) {
    sprintf(tmp, "# scan_start_time := %s\n", csv->c[ii].content);
    strcat(dft->comments, tmp);
  }
  /* frame times */
  for(fi=0; fi<dft->frameNr; fi++) {
    ii= 21 + fi*9 + 7; if(ii>csv->nr-2) {dftEmpty(dft); return CSV_INVALIDFORMAT;}
    //printf("cell[%d] := '%s'\n", ii, csv->c[ii].content); fflush(stdout);
    v1=atof_dpi(csv->c[ii].content); v2=atof_dpi(csv->c[ii+1].content);
    dft->x[fi]=v1; dft->x1[fi]=v1-0.5*v2; dft->x2[fi]=v1+0.5*v2;
  }
  /* region names, volumes, and concentrations */
  for(ri=0; ri<dft->voiNr; ri++) {
    /* ROI name */
    ii= 21 + ri*dft->frameNr*9;
    if(ii>=csv->nr) {dftEmpty(dft); return CSV_INVALIDFORMAT;}
    //printf("ri=%d cell[%d] := '%s'\n", ri, ii, csv->c[ii].content); fflush(stdout);
    strncpy(dft->voi[ri].name, csv->c[ii].content, MAX_REGIONNAME_LEN);
    dft->voi[ri].name[MAX_REGIONNAME_LEN]=(char)0;
    rnameSplit(csv->c[ii].content, dft->voi[ri].voiname,
  	       dft->voi[ri].hemisphere, dft->voi[ri].place,
               MAX_REGIONSUBNAME_LEN);
    /* Volume */
    ii++; if(ii>=csv->nr) {dftEmpty(dft); return CSV_INVALIDFORMAT;}
    dft->voi[ri].size=atof_dpi(csv->c[ii].content);
    /* Frame concentrations */
    ii++; for(fi=0; fi<dft->frameNr; fi++) {
      if((ii+6)>=csv->nr) {dftEmpty(dft); return CSV_INVALIDFORMAT;}
      /* Get concentration */
      dft->voi[ri].y[fi]=atof_dpi(csv->c[ii].content);
      /* Get concentration SD (probably not needed) */
      dft->voi[ri].y2[fi]=atof_dpi(csv->c[ii+1].content);
      /* check that frame times are correct */
      v1=atof_dpi(csv->c[ii+5].content);
      if(dft->x[fi]!=v1) {dftEmpty(dft); return CSV_INVALIDFORMAT;}
      ii+=9;
    } // next frame
  } // next region

  return CSV_OK;
}
/*****************************************************************************/

/*****************************************************************************/
/** Reads LinkSet data into DFT struct.
\return Returns 0 when successful, otherwise an error code.
 */
int csv2dft_linkset(
  /** Pointer to CSV data to be converted */
  CSV *csv,
  /** Pointer to empty DFT struct which will be allocated and filled here */
  DFT *dft
) {
  int ri, ci, fi, sci, ii, ret;
  double v;
  char *cptr;

  /* Check input data */
  if(CSV_TEST>2) {printf("csv2dft_linkset()\n"); fflush(stdout);}
  if(csv==NULL || dft==NULL) return CSV_ERROR;
  if(csv->nr<2 || csv->row_nr<1 || csv->col_nr<1) return CSV_INVALIDFORMAT;
  if(strcasecmp(csv->c[0].content, "LinkSet")!=0) return CSV_INVALIDFORMAT;
  /* Get the ROI nr */
  for(ii=1, ri=0; ii<csv->nr; ii++) if(csv->c[ii].col==1) {
    if(csv->c[ii].content==NULL) continue;
    if(strncmp(csv->c[ii].content, "VOI:", 4)==0) ri++; 
  }
  if(CSV_TEST>2) {
    printf("frame_nr=%d voi_nr=%d\n", csv->col_nr-2, ri); 
    fflush(stdout);
  }
  if(ri<1 || csv->col_nr<3) return CSV_INVALIDFORMAT;
  /* Delete previous DFT contents */
  dftEmpty(dft);
  /* Allocate memory for DFT */
  ret=dftSetmem(dft, csv->col_nr-2, ri); if(ret!=0) return CSV_OUTOFMEMORY;
  dft->voiNr=ri;
  dft->frameNr=csv->col_nr-2;

  /* Set DFT defaults */
  dft->timetype=DFT_TIME_MIDDLE;
  dft->_type=DFT_FORMAT_STANDARD;
  dft->isweight=0;
  dftUnitToDFT(dft, CUNIT_UNKNOWN);
  dft->timeunit=TUNIT_UNKNOWN;
  for(fi=0; fi<dft->frameNr; fi++) dft->w[fi]=1.0;
  for(ri=0; ri<dft->voiNr; ri++) dft->voi[ri].sw=0;
  
  /* Fill DFT contents */
  for(ri=0; ri<dft->voiNr; ri++) {
    if(CSV_TEST>3) printf("reading VOI %d\n", 1+ri);
    /* Search the (ri+1)th VOI */
    sci=0;
    for(ii=1; ii<csv->nr; ii++) if(csv->c[ii].col==1) {
      if(csv->c[ii].content==NULL) continue;
      if(strncmp(csv->c[ii].content, "VOI:", 4)==0) sci++;
      if(sci==ri+1) break; 
    }
    if(sci!=ri+1) break; // not found
    // ii is not the index of the start of data for one VOI
    if(CSV_TEST>5) printf("  ri=%d ii=%d row=%d col=%d\n",
                          ri, ii, csv->c[ii].row, csv->c[ii].col);
    /* Take the ROI number from there in case better name is not found later */
    strncpy(dft->voi[ri].name, csv->c[ii].content+4, MAX_REGIONNAME_LEN);
    dft->voi[ri].name[MAX_REGIONNAME_LEN]=(char)0;
    /* Read time unit from the next field */
    if(CSV_TEST>4 && ri==0) printf("reading time unit\n");
    cptr=strchr(csv->c[ii+1].content, '(');
    if(cptr!=NULL) {
      ci=petTunitId(cptr+1);
      if(ri==0) {
        dft->timeunit=ci;
      } else if(dft->timeunit!=ci) {
        if(CSV_TEST>0) printf("different time units.\n");
        return CSV_INVALIDFORMAT;
      }
    }
    if(CSV_TEST>4 && ri==0) printf("time unit: %s\n", petTunit(dft->timeunit));
    /* Read the times */
    for(fi=0; fi<dft->frameNr; fi++) {
      ci=ii+fi+2; if(ci>=csv->nr) return CSV_INVALIDFORMAT;
      cptr=csv->c[ci].content;
      if(atof_with_check(cptr, &v)!=0) return CSV_INVALIDFORMAT;
      if(ri==0) {
        dft->x[fi]=v;
      } else {
        if(fabs(v-dft->x[fi])>1.0E-003) return CSV_INVALIDFORMAT;
      }
    } // next frame time
    /* First column on the 2nd row should contain the region name */
    if(CSV_TEST>4) printf("reading VOI name\n");
    ii+=csv->col_nr; if(ii>=csv->nr) return CSV_INVALIDFORMAT;
    cptr=csv->c[ii].content;
    if(cptr!=NULL) {
      strncpy(dft->voi[ri].name, cptr, MAX_REGIONNAME_LEN);
      dft->voi[ri].name[MAX_REGIONNAME_LEN]=(char)0;
    } // default name was set before
    rnameSplit(dft->voi[ri].name, dft->voi[ri].voiname,
  	       dft->voi[ri].hemisphere, dft->voi[ri].place,
               MAX_REGIONSUBNAME_LEN);
    /* VOI average values should be two rows below this */
    if(CSV_TEST>4) printf("reading VOI values\n");
    ii+=2*csv->col_nr; if(ii>=csv->nr) return CSV_INVALIDFORMAT;
    /* Check that from 2nd column, and get from there the unit, too */
    ii++; if(ii>=csv->nr) return CSV_INVALIDFORMAT;
    cptr=csv->c[ii].content; if(cptr==NULL) return CSV_INVALIDFORMAT;
    if(strncasecmp(cptr, "Average", 7)!=0) return CSV_INVALIDFORMAT;
    cptr=strchr(csv->c[ii].content, '(');
    if(cptr!=NULL) {
      if(CSV_TEST>4 && ri==0)
        printf("reading activity unit from string: '%s'\n", cptr+1);
      ci=petCunitId(cptr+1);
      if(ri==0) {
        strcpy(dft->unit, petCunit(ci));
      } else if(ci!=petCunitId(dft->unit)) {
        if(CSV_TEST>0) printf("different concentration units.\n");
        return CSV_INVALIDFORMAT;
      }
      if(CSV_TEST>5 && ri==0)
        printf("unit := %s (%d)\n", petCunit(ci), ci);
    }
    /* then read the concentrations */
    for(fi=0; fi<dft->frameNr; fi++) {
      ci=ii+fi+1; if(ci>=csv->nr) return CSV_INVALIDFORMAT;
      cptr=csv->c[ci].content;
      if(atof_with_check(cptr, &v)!=0) return CSV_INVALIDFORMAT;
      dft->voi[ri].y[fi]=v;
    } // next frame time 
  }
  if(ri<dft->voiNr) return CSV_INVALIDFORMAT;

  return CSV_OK;
}
/*****************************************************************************/

/*****************************************************************************/
/** Reads Mat data into DFT struct.
\return Returns 0 when successful, otherwise an error code.
 */
int csv2dft_mat(
  /** Pointer to CSV data to be converted */
  CSV *csv,
  /** Pointer to empty DFT struct which will be allocated and filled here */
  DFT *dft
) {
  int ri, fi, ret;
  char *cptr;

  /* Check input data */
  if(CSV_TEST>2) {printf("csv2dft_mat()\n"); fflush(stdout);}
  if(csv==NULL || dft==NULL) return CSV_ERROR;
  if(csv->nr<4 || csv->row_nr<2 || csv->col_nr<2) return CSV_INVALIDFORMAT;
  if(!csvIsRegular(csv)) return CSV_INVALIDFORMAT;

  /* Delete previous DFT contents */
  dftEmpty(dft);
  /* Allocate memory for DFT */
  ret=dftSetmem(dft, csv->col_nr-1, csv->row_nr-1); if(ret!=0) return CSV_OUTOFMEMORY;
  dft->voiNr=csv->row_nr-1;
  dft->frameNr=csv->col_nr-1;
  if(CSV_TEST>2) {printf("frame_nr=%d voi_nr=%d\n", dft->frameNr, dft->voiNr); fflush(stdout);}

  /* Set DFT defaults */
  dft->timetype=DFT_TIME_STARTEND;
  dft->_type=DFT_FORMAT_STANDARD;
  dft->isweight=0;
  dftUnitToDFT(dft, CUNIT_UNKNOWN);
  dft->timeunit=TUNIT_UNKNOWN;
  for(fi=0; fi<dft->frameNr; fi++) dft->w[fi]=1.0;
  for(ri=0; ri<dft->voiNr; ri++) dft->voi[ri].sw=0;

  if(CSV_TEST>200)
    for(int i=0; i<csv->nr; i++)
      printf("row=%d col=%d content='%s'\n", csv->c[i].row, csv->c[i].col, csv->c[i].content);

  /* First cell may contain study number */
  cptr=csvCell(csv, 1, 1);
  if(cptr!=NULL && cptr[0] && strnlen(cptr, 20)<20) {
    strlcpy(dft->studynr, cptr, MAX_STUDYNR_LEN);
    if(CSV_TEST>3) printf("studynr := %s\n", dft->studynr);
  }
  
  /* Fill DFT contents */
  ret=0;
  for(ri=0; ri<dft->voiNr; ri++) {
    if(CSV_TEST>3) printf("reading VOI %d\n", 1+ri);
    /* Get ROI name */
    cptr=csvCell(csv, ri+2, 1); //printf("cptr='%s'\n", cptr);
    if(cptr==NULL) {ret++; break;}
    int n=strlen(cptr); if(n<3) {ret++; break;}
    if(cptr[n-1]=='\'') {cptr[n-1]=(char)0; n--;} else {ret++; break;}
    if(cptr[0]=='\'') {cptr++; n--;} else {ret++; break;}
    strlcpy(dft->voi[ri].name, cptr, MAX_REGIONNAME_LEN);
    rnameSplit(dft->voi[ri].name, dft->voi[ri].voiname, dft->voi[ri].hemisphere, dft->voi[ri].place,
               MAX_REGIONSUBNAME_LEN);
    /* Get concentrations */
    for(fi=0; fi<dft->frameNr; fi++) {
      cptr=csvCell(csv, ri+2, fi+2); // printf("cptr='%s'\n", cptr);
      ret=atof_with_check(cptr, &dft->voi[ri].y[fi]);
      if(ret) break;
    }
    if(ret!=0) break;
  }
  if(ret!=0) {dftEmpty(dft); return CSV_INVALIDFORMAT;}

  /* Get frame times */
  if(CSV_TEST>3) printf("reading frames\n");
  ret=0;
  for(fi=0; fi<dft->frameNr; fi++) {
    cptr=csvCell(csv, 1, fi+2); if(cptr==NULL) {ret++; break;}
    int n=strlen(cptr); if(n<3) {ret++; break;}
    char *cptr2=strchr(cptr+1, '-'); if(cptr2!=NULL) {*cptr2=(char)0; cptr2++;}
    ret=atof_with_check(cptr, &dft->x1[fi]);
    if(ret==0) ret=atof_with_check(cptr2, &dft->x2[fi]);
    if(ret) break;
    dft->x[fi]=0.5*(dft->x1[fi]+dft->x2[fi]);
  }
  if(ret!=0) {dftEmpty(dft); return CSV_INVALIDFORMAT;}

  return CSV_OK;
}
/*****************************************************************************/

/*****************************************************************************/
/** Check whether CSV is regular, that is, each row contain the same number of columns.
 *  @author Vesa Oikonen
 *  @sa csvRead
 *  @return Returns 1 if CSV is regular, 0 if not.
 */
int csvIsRegular(
  /** Pointer to CSV */
  CSV *csv
) {
  if(csv==NULL) return(0);
  if(csv->nr<2) return(1);
  int i, r, n=0, m=0;
  i=0; r=csv->c[i].row; m++;
  for(i=1; i<csv->nr; i++) {
    if(r==csv->c[i].row) {m++; continue;}
    if(n>0 && m!=n) return(0);
    r=csv->c[i].row; n=m; m=1;
  }
  if(n>0 && m!=n) return(0);
  return(1);
}
/*****************************************************************************/

/*****************************************************************************/
/** Get the CVS field contents in specified row and column.
 *  @return Returns pointer to the content string, or NULL if not found.
 *  @author Vesa Oikonen
 *  @sa csvRead
 */
char* csvCell(
  /** Pointer to CSV */
  CSV *csv,
  /** CSV row index */
  int row,
  /** CSV column index */
  int col
) {
  if(csv==NULL) return((char*)NULL);
  for(int i=0; i<csv->nr; i++) 
    if(csv->c[i].row==row && csv->c[i].col==col)
      return(csv->c[i].content);
  return((char*)NULL); 
}
/*****************************************************************************/

/*****************************************************************************/
