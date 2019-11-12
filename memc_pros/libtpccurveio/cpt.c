/// @file cpt.c
/// @author Vesa Oikonen
/// @brief Reading and writing CPT files.
///
/*****************************************************************************/
#include "libtpccurveio.h"
/*****************************************************************************/

/*****************************************************************************/
/** Print the compilation date and time to specified FILE pointer  */
void libcpt_printdate(FILE *fp)
{
  fprintf(fp, "libcpt compiled on %s %s\n", __DATE__, __TIME__);
}
/*****************************************************************************/

/*****************************************************************************/
/** Split region name into 1-3 subparts of given max length.
\return Returns the number of subparts.
 */
int cptrnameSplit(
  /** Region name to split (string is not edited) */
  char *rname,
  /** Pointer to 1st subname (anatomical region) */
  char *name1,
  /** Pointer to 2nd subname (usually hemisphere) */
  char *name2,
  /** Pointer to 3rd subname (usually image plane) */
  char *name3,
  /** Max lenght of subnames, excluding terminal null */
  int max_name_len
) {
  char temp[MAX_REGIONNAME_LEN+1], *cptr, *cptr2;
  int nr=0;

  if(rname==NULL || name1==NULL || name2==NULL || name3==NULL) return(0);
  if(max_name_len<1) return(0);
  name1[0]=name2[0]=name3[0]=(char)0;
  strncpy(temp, rname, MAX_REGIONNAME_LEN); temp[MAX_REGIONNAME_LEN]=(char)0;
  cptr=strtok(temp, " _\t\n\r"); if(cptr==NULL) return(nr);
  cptr2=strstr(cptr, "dx"); if(cptr2==NULL) cptr2=strstr(cptr, "sin");
  if(cptr2!=NULL) {strcpy(name2, cptr2); *cptr2=(char)0; nr++;}
  else strcpy(name2, ".");
  strncpy(name1, cptr, max_name_len); name1[max_name_len]=(char)0; nr++;
  cptr=strtok(NULL, " _\t\n\r"); if(cptr==NULL) return(nr);
  strncpy(name3, cptr, max_name_len); name3[max_name_len]=(char)0; nr++;
  return(nr);
}
/*****************************************************************************/

/*****************************************************************************/
/** Read TACs from CPT file into a DFT.
\return Returns 0 if successful, sets cpterrmsg in case of an error.
 */
int cptReadOne(
  /** CPT filename */
  char *cptfile,
  /** Pointer to DFT where TACs are read; must be empty and initialized */
  DFT *dft
) {
  int fi, fj, ri, ii, li, rn, fn, n, ret, colNr=10, lineNr=0;
  char *cptr, tmp[512];
  IFT ift;
  int title_line, roi_col, roi_nr=0, frame_nr=0;


  /* Check input */
  if(cptfile==NULL || strlen(cptfile)<1 || dft==NULL) {
    strcpy(cpterrmsg, "program error"); return(1);
  }
  dftEmpty(dft);

  /*
   *  Read CPT file as IFT
   */
  iftInit(&ift);
  ret=iftRead(&ift, cptfile, 0);
  if(ret) {
    strcpy(cpterrmsg, ift.status);
    iftEmpty(&ift); return(4);
  }
  if(CPT_TEST>0) iftWrite(&ift, "stdout");
  
  /*
   *  Find the line of data titles
   */
  title_line=iftFindNthValue(&ift, " ROI Avg ", 1);
  if(title_line<0) {
    sprintf(cpterrmsg, "unsupported filetype");
    iftEmpty(&ift); return(6);
  }
  /* From title line, check if data consists of 10 or 11 columns */
  strcpy(tmp, ift.item[title_line].value);
  cptr=strtok(tmp, " \t\n\r");
  if(strcasecmp(cptr, "Frame")!=0) {
    sprintf(cpterrmsg, "unsupported filetype");
    iftEmpty(&ift); return(7);
  }
  cptr=strtok(NULL, " \t\n\r");
  colNr=10; if(strcasecmp(cptr, "ROI")!=0) colNr++;
  if(colNr==10) roi_col=1; else roi_col=2;
  if(CPT_TEST>0) printf("title_line=%d colNr=%d\n", title_line, colNr);

  /*
   *  Mark data lines with sw=1
   */
  for(ii=lineNr=0; ii<ift.keyNr; ii++) {
    /* Data can not start before title line and the unit line below it */
    if(ii<title_line+2) {ift.item[ii].sw=0; continue;}
    /* Comment lines cannot contain data, or at least we don't read it */
    if(ift.item[ii].type!=' ' && ift.item[ii].type!='\0') {
      ift.item[ii].sw=0; continue;}
    /* and there cannot be any key */
    if(ift.item[ii].key!=NULL && strlen(ift.item[ii].key)>0) {
      ift.item[ii].sw=0; continue;}
    /* Check that column number matches */
    n=0; strncpy(tmp, ift.item[ii].value, 511); tmp[511]=(char)0;
    cptr=strtok(tmp, " \t\n\r");
    while(cptr!=NULL) {n++; cptr=strtok(NULL, " \t\n\r");}
    if(n!=colNr) {ift.item[ii].sw=0; continue;}
    ift.item[ii].sw=1; lineNr++;
  }
  if(lineNr<1)  {
    sprintf(cpterrmsg, "unsupported filetype");
    iftEmpty(&ift); return(8);
  }
  if(CPT_TEST>0) printf("lineNr=%d\n", lineNr);

  /*
   *  Read the data into a new, local table
   */
  double tactable[colNr][lineNr];
  for(ii=0, n=0; ii<ift.keyNr; ii++) if(ift.item[ii].sw==1) {
    ret=sscanf(ift.item[ii].value, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
      &tactable[0][n], &tactable[1][n], &tactable[2][n], &tactable[3][n],
      &tactable[4][n], &tactable[5][n], &tactable[6][n], &tactable[7][n],
      &tactable[8][n], &tactable[9][n], &tactable[10][n] );
    n++; /*printf("ret=%d\n", ret);*/
  }
  if(CPT_TEST>0) { /* print the table */
    printf("\n");
    for(fi=0; fi<lineNr; fi++) {
      for(ri=0; ri<colNr; ri++) printf(" %g", tactable[ri][fi]); 
      printf("\n");
    }
    printf("\n");
  }

  /*
   *  Compute the number of ROIs and Frames from the data lines
   */
  for(fi=1, roi_nr=frame_nr=1; fi<lineNr; fi++) {
    /* have these been seen before? */
    for(fj=0, fn=rn=0; fj<fi; fj++) {
      if(tactable[roi_col][fj]==tactable[roi_col][fi]) rn++;
      if(tactable[0][fj]==tactable[0][fi]) fn++;
    }
    if(rn==0) roi_nr++;
    if(fn==0) frame_nr++;
  }
  if(CPT_TEST>0) printf("roi_nr=%d frame_nr=%d\n", roi_nr, frame_nr);
  /* Also, test that the highest frame number equals number of frames */
  for(fi=0, n=0; fi<lineNr; fi++) if(tactable[0][fi]>n) n=tactable[0][fi];
  if(n!=frame_nr) {
    sprintf(cpterrmsg, "frames are not consequential");
    iftEmpty(&ift); return(9);
  }
  /* And also, test that roi_nr*frame_nr == nr of data lines */
  if(roi_nr*frame_nr!=lineNr) {
    sprintf(cpterrmsg, "missing or extra samples");
    iftEmpty(&ift); return(10);
  }

  /*
   *  Copy data into DFT
   */

  /* Allocate memory for DFT */
  ret=dftSetmem(dft, frame_nr, roi_nr);
  if(ret) {
    sprintf(cpterrmsg, "cannot allocate memory");
    iftEmpty(&ift); dftEmpty(dft); return(11);
  }
  dft->frameNr=frame_nr; dft->voiNr=roi_nr; dft->_type=1;

  /* Make a list of ROI ID numbers */
  int roi_id[roi_nr];
  n=0; roi_id[n++]=tactable[roi_col][0];
  for(fi=1; fi<lineNr; fi++) {
    for(fj=0; fj<fi; fj++)
      if(tactable[roi_col][fj]==tactable[roi_col][fi]) continue;
    roi_id[n++]=(int)(0.5+tactable[roi_col][fi]);
  }
  if(CPT_TEST>0) {
    printf("List of ROI ID numbers:\n");
    for(ri=0; ri<roi_nr; ri++) printf("   %d : %d\n", ri+1, roi_id[ri]);
  }

  /* Extract data for one ROI at a time */
  for(ri=0; ri<roi_nr; ri++) {
    fi=0; /* all frames of this ROI */
    for(li=0; li<lineNr; li++) if(roi_id[ri]==(int)(0.5+tactable[roi_col][li])) {
      /* Activity average value */
      dft->voi[ri].y[fi]=tactable[colNr-8][li];
      /* Only once for each ROI */
      if(fi==0) {
        /* ROI volume */
        dft->voi[ri].size=tactable[colNr-1][li];
        /* ROI name (may be changed later) */
        sprintf(dft->voi[ri].voiname, "ROI%03d", roi_id[ri]);
        if(colNr>10) sprintf(dft->voi[ri].place, "Pl%04.0f", tactable[1][ri]);
        snprintf(dft->voi[ri].name, MAX_REGIONNAME_LEN, "%s . %s",
	  dft->voi[ri].voiname, dft->voi[ri].place);
      }
      /* Get frame times from the first ROI */
      if(ri==0) {
        dft->x1[fi]=tactable[colNr-4][li];
        dft->x2[fi]=dft->x1[fi]+tactable[colNr-3][li];
        dft->x[fi]=0.5*(dft->x1[fi]+dft->x2[fi]);
      }
      fi++;
    } /* next frame */
  } /* next ROI */
  /* Convert frame times into minutes */
  dftSec2min(dft); dft->timetype=3;

  /*
   *  Get studynumber from filename
   */
  studynr_from_fname(cptfile, dft->studynr);

  /*
   *  Try to find the data unit
   */
  ii=iftFindNthValue(&ift, "In units of ", 1);
  if(ii>=0) {
    strncpy(tmp, ift.item[ii].value+12, 511); tmp[511]=(char)0;
    cptr=strtok(tmp, " \t\n\r,;");
    strncpy(dft->unit, cptr, MAX_UNITS_LEN);
    dft->unit[MAX_UNITS_LEN]=(char)0;
  } else {
    strcpy(tmp, "Units"); ii=iftGet(&ift, tmp);
    if(ii>=0) {
      strncpy(dft->unit, ift.item[ii].value, MAX_UNITS_LEN);
      dft->unit[MAX_UNITS_LEN]=(char)0;
    }   
  }

  /*
   *  Try to extract the region name (if only one ROI)
   */
  if(dft->voiNr==1) {
    ii=iftFindNthKey(&ift, "Using ROI ", 1);
    if(ii>=0) {
      cptr=strchr(ift.item[ii].key, '\"');
      if(cptr!=NULL) {
        strncpy(dft->voi[0].name, cptr+1, MAX_REGIONNAME_LEN);
        dft->voi[0].name[MAX_REGIONNAME_LEN]=(char)0;
        cptr=strchr(dft->voi[0].name, '\"');
        if(cptr!=NULL) *cptr=(char)0;
      }
    }
    cptrnameSplit(dft->voi[0].name,
      dft->voi[0].voiname, dft->voi[0].hemisphere, dft->voi[0].place,
      MAX_REGIONSUBNAME_LEN);
  }

  /*
   *  Try to find the plane number
   */
  ii=iftFindNthKey(&ift, "Plane ", 1);
  if(ii>=0) {
    n=atoi(ift.item[ii].key+6);
    if(CPT_TEST>0) printf("Plane %d\n", n);
    if(n>0) for(ri=0; ri<dft->voiNr; ri++)
      snprintf(dft->voi[ri].place, 7, "Pl%04d", n);
  } else {
    ii=iftFindNthValue(&ift, "Plane ", 1);
    if(ii>=0) {
      n=atoi(ift.item[ii].value+6);
      if(CPT_TEST>0) printf("Plane %d\n", n);
      if(n>0) for(ri=0; ri<dft->voiNr; ri++)
        snprintf(dft->voi[ri].place, 7, "Pl%04d", n);
    }
  }

  /* Set weights in DFT to 1.0 */
  dft->isweight=0;
  for(fi=0; fi<dft->frameNr; fi++) dft->w[fi]=1.0;


  iftEmpty(&ift);
  return(0);
}
/*****************************************************************************/

/*****************************************************************************/
/** Write TAC data in CPT (Imagetool) format. If TACs are from different
    planes, then each plane will be saved in its own file.
\return Returns 0 if successful; in case of an error, sets string cpterrmsg.
*/
int cptWrite(
  /** TAC data to write */
  DFT *dft,
  /** CPT path and filename without extension,
      because this function may need to add plane number before .cpt */
  char *filename,
  /** Specific CPT format: 0=default, others not yet supported */
  int cpt_format
) {
  int ri, fi, ret, n, roi_id;
  FILE *fp;
  char cptfile[FILENAME_MAX], tmp[256];

  /* Check input */
  if(dft==NULL || dft->voiNr<1 || dft->frameNr<1 || strlen(filename)<1) {
    strcpy(cpterrmsg, "program error"); return(1);
  }
  if(cpt_format!=0) strcpy(cpterrmsg, "cpt format not supported yet");
  /* Sort data by plane */
  ret=dftSortPlane(dft);
  if(ret) {
    strcpy(cpterrmsg, "error in data file"); return(2);
  }
  /* Convert times to seconds if necessary */
  if(dft->timeunit==TUNIT_MIN) dftMin2sec(dft);

  
  for(ri=0; ri<dft->voiNr; ri++) {
    /* Construct CPT filename */
    strcpy(cptfile, filename);
    if(strlen(dft->voi[ri].place)>0 && strcmp(dft->voi[ri].place, ".")!=0) {
      snprintf(cptfile, FILENAME_MAX, "%s_%s.cpt", filename, dft->voi[ri].place);
    } else {
      snprintf(cptfile, FILENAME_MAX, "%s.cpt", filename);
    }
    /* Open file */
    fp=fopen(cptfile, "w");
    if(fp==NULL) {
      sprintf(cpterrmsg, "cannot open file for write"); return(5);
    }
    /* Write 1st title */
    if(strlen(dft->unit)==0) strcpy(tmp, "unknown"); else strcpy(tmp, dft->unit);
    //fprintf(fp, "# In units of %s per pixel per second\n", tmp);
    fprintf(fp, "# In units of %s\n", tmp);
    /* Write 2nd title */
    if(strncasecmp(dft->voi[ri].place, "Pl", 2)==0)
      sprintf(tmp, "%d", atoi(dft->voi[ri].place+2));
    else if(strlen(dft->voi[ri].place)==0 || strcmp(dft->voi[ri].place, ".")==0)
      strcpy(tmp, "1");
    else
      strcpy(tmp, dft->voi[ri].place);
    fprintf(fp, "Plane %-6.6s Scan Start Date (d m y): 1 1 1980     Scan Start Time (h m s): 0 0 0\n\n", tmp);
    /* Write 3rd title */
    fprintf(fp, "Frame  ROI ID        ROI Avg    #pixels    ROI Total   %%Stdev    Offset   Duration   ROI Surf.     ROI Vol.\n");
    fprintf(fp, "                                (screen)                          (sec)     (sec)     mmxmm       mmxmmxmm\n");
    /* Write all frames */
    for(fi=0, n=ri; fi<dft->frameNr; fi++) {
      n=ri;
      /* current region */
      strcpy(tmp, dft->voi[ri].voiname);
      if(strlen(dft->voi[ri].hemisphere)>0 && strcmp(dft->voi[ri].hemisphere, ".")!=0) {
        strcat(tmp, " "); strcat(tmp, dft->voi[ri].hemisphere);}
      strcpy(tmp, "");
      if(strncasecmp(dft->voi[ri].voiname, "ROI", 3)==0 && atoi(dft->voi[ri].voiname+3)>0)
        roi_id=atoi(dft->voi[ri].voiname+3); else roi_id=1;
      fprintf(fp, "%-6d %-3d %-8.8s %11.4e %5d     %10.4e %7.1f  %9.1f %9.1f    %10.4e    %10.4e\n", 
        fi+1, roi_id, tmp, dft->voi[ri].y[fi], 0, 0.0, 0.0,
        dft->x1[fi], dft->x2[fi]-dft->x1[fi], 0.0, dft->voi[ri].size );
      /* the rest of regions on the same plane */
      while(n<dft->voiNr-1 && strcasecmp(dft->voi[ri].place, dft->voi[n+1].place)==0) {
        n++; strcpy(tmp, dft->voi[n].voiname);
        if(strlen(dft->voi[n].hemisphere)>0 && strcmp(dft->voi[n].hemisphere, ".")!=0) {
          strcat(tmp, " "); strcat(tmp, dft->voi[n].hemisphere);}
        strcpy(tmp, "");
        if(strncasecmp(dft->voi[n].voiname, "ROI", 3)==0 && atoi(dft->voi[n].voiname+3)>0)
          roi_id=atoi(dft->voi[n].voiname+3); else roi_id=n-ri+1;
        fprintf(fp, "%-6d %-3d %-8.8s %11.4e %5d     %10.4e %7.1f  %9.1f %9.1f    %10.4e    %10.4e\n", 
          fi+1, roi_id, tmp, dft->voi[n].y[fi], 0, 0.0, 0.0,
          dft->x1[fi], dft->x2[fi]-dft->x1[fi], 0.0, dft->voi[n].size );
      }
    }
    ri=n;
#if(0)
    /* DON'T DO THIS; at least "asetaatti" program does not work with this */
    /* In the end, write the number of frames */
    fprintf(fp, "\n# %d Frame(s) analyzed.\n", dft->frameNr);
#endif
    fclose(fp);
  }

  return(0);
}
/*****************************************************************************/

/*****************************************************************************/
