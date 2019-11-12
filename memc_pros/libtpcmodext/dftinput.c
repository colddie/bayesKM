/// @file dftinput.c
/// @brief Procedures for handling model input data.
/// @author Vesa Oikonen
///
/*****************************************************************************/

/*****************************************************************************/
#include "libtpcmodext.h"
/*****************************************************************************/

/*****************************************************************************/

/*****************************************************************************/
/** Read DFT format input TAC data to match the time frames in the
    specified tissue DFT file.
    Instead of input filename, a reference region name can be given: then all
    the matching tacs (based on roi names) are moved from the roi data to the
    input data, with the best match as first tac.
    Input data is interpolated (if necessary) and also integrated to y2[].
    Input time and concentration units are tried to be converted to be the same
    as in tissue data.
   @return Returns 0 if succesful, and >0 in case of an error, and specifically
    101 in case input TAC is not valid.
 */
int dftReadinput(
  /** Input data, previous contents are cleared */
  DFT *input,
  /** PET data containing frames and possible input regions */
  DFT *tissue,
  /** Input data filename, or region name in PET data */
  char *filename,
  /** Type of input, as found out here; 1 and 2 =tac, 3=fit, 5=region name;
      enter NULL, if not needed */
  int *filetype,
  /** First time of input data before interpolation (in tissue time units);
      use this to check that required time range was measured; 
      NULL if not needed */
  double *ti1,
  /** Last time of input data before interpolation (in tissue time units);
      use this to check that required time range was measured; 
      NULL if not needed */
  double *ti2,
  /** Set to <>0 to verify that input TAC does not start too late and has
   *  decent peak to provide reliable AUC0-T. */
  int verifypeak,
  /** Pointer to a string (allocated for at least 64 chars) where error message
      or other execution status will be written; enter NULL, if not needed */   
  char *status,
  /** Verbose level; if zero, then nothing is printed to stderr or stdout */
  int verbose
) {
  int ri, n, ret, ftype=0;
  DFT temp;
  FILE *fp;

  if(verbose>0)
    printf("dftReadinput(inp, tis, %s, type, ti1, ti2, %d, status, %d)\n",
           filename, verifypeak, verbose);
  if(filetype!=NULL) *filetype=0;

  /* Check that PET data is ok */
  if(input==NULL || tissue==NULL || filename==NULL) {
    if(status!=NULL) strcpy(status, "program error");
    return 1;
  }
  if(tissue->frameNr<1 || tissue->voiNr<1) {
    if(status!=NULL) strcpy(status, "no pet data");
    return 2;
  }

  /* Delete any previous input data and initiate temp data */
  dftEmpty(input); dftInit(&temp);

  /* Try to figure out what is the 'filename' */
  ftype=0;
  /* Can we open it as a file? */
  fp=fopen(filename, "r");
  if(fp!=NULL) {
    /* File can be opened for read; close it for now */
    if(verbose>1) printf("  file can be openened for reading.\n");
    if(filetype!=NULL) *filetype=1;
    fclose(fp);

    /* Try to identify the file format */
    ftype=dftFormat(filename);
    if(ftype==DFT_FORMAT_UNKNOWN) {
      if(filetype!=NULL) *filetype=0;
      if(status!=NULL) strcpy(status, "unknown file format");
      return 3;
    } else if(ftype==DFT_FORMAT_FIT) {
      if(filetype!=NULL) *filetype=3;
      if(status!=NULL) strcpy(status, "cannot read fit file"); 
      return 3;
    }
    if(verbose>2) printf("  fileformat=%d\n", ftype);

    /* Try to read it */
    if((ret=dftRead(filename, &temp))!=0) {
      if(filetype!=NULL) *filetype=0;
      if(status!=NULL) sprintf(status, "cannot read file (%d)", ret);
      return(2);
    }
    /* Convert input time units to the same as in tissue data */
    (void)dftTimeunitConversion(&temp, tissue->timeunit);
    /* Check the tissue and plasma TAC concentration units */
    ret=dftUnitConversion(&temp, petCunitId(tissue->unit));
    if(ret!=0) {
      sprintf(status, "check the units of input and tissue data");
    }
    /* Tell user what was the original input time range */
    if(temp.timetype==DFT_TIME_STARTEND) {
      if(ti1!=NULL) *ti1=temp.x1[0];
      if(ti2!=NULL) *ti2=temp.x2[temp.frameNr-1];
    } else {
      if(ti1!=NULL) *ti1=temp.x[0];
      if(ti2!=NULL) *ti2=temp.x[temp.frameNr-1];
    }
    /* Verify the peak if requested */
    if(verifypeak!=0) {
      ret=dftVerifyPeak(&temp, 0, verbose-2, status);
      //if(ret!=0) sprintf(status, "input TAC should start at time zero");
      if(ret>0) {dftEmpty(&temp); return 101;}
    }
    /* Interpolate and integrate data to pet times */
    ret=dftInterpolate(&temp, tissue, input, status, verbose);
    dftEmpty(&temp);
    if(ret!=0) return 4;

  } else {

    /* it's not a file, at least an accessible file, but is it a region name? */
    if(filetype!=NULL) *filetype=5;
    /* Select ROIs that match the specified input name */
    n=dftSelectRegions(tissue, filename, 1);
    if(n<=0) {
      if(status!=NULL) sprintf(status, "cannot find region");
      return 7;
    }
    if(n==tissue->voiNr) {
      if(status!=NULL) sprintf(status, "all regions do match");
      return 8;
    }
    /* one or more regions was found ; move them to input data */
    ret=dftdup(tissue, input); if(ret==0) {
      ri=0; ret=0; while(ri<input->voiNr && ret==0) {
        if(input->voi[ri].sw==0) ret=dftDelete(input, ri); else ri++;
      }
      if(ret==0) {
        ri=0; ret=0; while(ri<tissue->voiNr && ret==0) {
          if(tissue->voi[ri].sw!=0) ret=dftDelete(tissue, ri); else ri++;
        }
      }
    }
    if(ret!=0) {
      if(status!=NULL) sprintf(status, "cannot separate input regions");
      dftEmpty(input); return 9;
    }
    /* Try to select the best reference ROI */
    ri=dftSelectBestReference(input); if(ri<0) {
      if(status!=NULL) sprintf(status, "cannot separate input regions");
      dftEmpty(input); return 10;
    }
    /* And move it to the first place */
    if(ri>0) {dftMovevoi(input, ri, 0); ri=0;}
    if(verbose>1) printf("selected ref region := %s\n", input->voi[0].name);
    /* Verify the peak if requested */
    if(verifypeak!=0) {
      ret=dftVerifyPeak(input, 0, verbose-2, status);
      //if(ret!=0) sprintf(status, "input TAC should start at time zero");
      if(ret>0) {dftEmpty(input); return 101;}
    }

    /* Calculate integrals */
    for(ri=0; ri<input->voiNr; ri++) {
      if(input->timetype==DFT_TIME_STARTEND) {
        ret=petintegral(input->x1, input->x2, input->voi[ri].y, input->frameNr,
              input->voi[ri].y2, input->voi[ri].y3);
        if(ti1!=NULL) *ti1=input->x1[0];
        if(ti2!=NULL) *ti2=input->x2[input->frameNr-1];
      } else {
        ret=interpolate(input->x, input->voi[ri].y, input->frameNr,
               input->x, NULL, input->voi[ri].y2,
               input->voi[ri].y3, input->frameNr);
        if(ti1!=NULL) *ti1=input->x[0];
        if(ti2!=NULL) *ti2=input->x[input->frameNr-1];
      }
      if(ret) {
        if(status!=NULL) sprintf(status, "cannot integrate input");
        dftEmpty(input); return(11);
      }
    }
  }

  return(0);
}
/*****************************************************************************/

/*****************************************************************************/
/** Read DFT format reference region TAC data and add it into DFT already
 *  containing other tissue TACs, if reference region TAC(s) are be given
 *  in a separate file. Alternatively the reference region name can be given,
 *  which will then be selected from existing tissue TACs.
 *  Reference TACs are marked in DFT with sw=2 (best name match) or sw=1;
 *  if reference region file contains several TACs then the one which contains 
 *  name 'mean' or 'avg' or has shortest total name length is given value sw=2. 
 *  When necessary, reference data is interpolated and units converted to
 *  match the existing tissue DFT.
 *  Reference TAC is also integrated into y2[].
\return Returns the number of reference TACs, and <=0 in case of an error.
 */
int dftReadReference(
  /** PET data containing existing tissue TACs, possibly also ref regions */
  DFT *tissue,
  /** Reference TAC filename, or region name in previous data */
  char *filename,
  /** Type of input, as found out here; 1 and 2 =tac, 3=fit, 5=region name;
      enter NULL, if not needed */
  int *filetype,
  /** Index of the best reference region; enter NULL if not needed */
  int *ref_index,
  /** Pointer to a string (allocated for at least 64 chars) where error message
      or other execution status will be written; enter NULL, if not needed */     
  char *status,
  /** Verbose level; if zero, then nothing is printed to stderr or stdout */
  int verbose
) {
  int ri, n, ret;
  DFT temp;
  FILE *fp;
  int ftype=0;

  if(verbose>0) printf("dftReadReference(tis, %s, type, i, status, %d)\n", filename, verbose);

  /* Check that input is ok */
  if(tissue==NULL || filename==NULL || strlen(filename)<1) {
    if(status!=NULL) strcpy(status, "program error");
    return -1;
  }
  if(tissue->frameNr<1 || tissue->voiNr<1) {
    if(status!=NULL) strcpy(status, "no pet data");
    return -2;
  }

  /* Check if we can identify the reference as an unsupported file */
  ret=dftFormat(filename);
  if(ret==DFT_FORMAT_FIT) {
    if(filetype!=NULL) *filetype=3;
    if(status!=NULL) strcpy(status, "cannot read fit file");
    return -3;
  }


  /* Can we open it as a file? */
  fp=fopen(filename, "r");
  if(fp!=NULL) {
    fclose(fp);
    /* Try to read the reference as a TAC file */
    dftInit(&temp);
    if((ret=dftRead(filename, &temp))!=0) {
      if(status!=NULL) sprintf(status, "cannot read file (%d)", ret);
      return -4;
    }
    if(filetype!=NULL) *filetype=1;

    /* Convert ref time units to the same as in tissue data */
    (void)dftTimeunitConversion(&temp, tissue->timeunit);
    /* Check the tissue and ref TAC concentration units */
    if((ret=dftUnitConversion(&temp, petCunitId(tissue->unit)))!=0) {
      sprintf(status, "check the units of reference and tissue data");
    }
    /* Interpolate and integrate data to pet times */
    /* this also verifies that ref data does not need extensive extrapolation */
    ret=dftInterpolateInto(&temp, tissue, status, verbose);
    if(ret!=0) {dftEmpty(&temp); return -5;}
    /* Set switches to define which are reference regions and which are not */
    for(ri=0; ri<tissue->voiNr-temp.voiNr; ri++) tissue->voi[ri].sw=0;
    for(; ri<tissue->voiNr; ri++) tissue->voi[ri].sw=1;
    /* Find the best reference region */
    n=temp.voiNr; dftEmpty(&temp);
    if(n==1) {
      ri=tissue->voiNr-n;
    } else {
      if((ri=dftSelectBestReference(tissue))<0) {
        sprintf(status, "cannot select the best reference region");
        dftEmpty(&temp); return -6;
      }
    }
    tissue->voi[ri].sw=2; if(ref_index!=NULL) *ref_index=ri;
    if(verbose>1) printf("selected ref region := %s\n", tissue->voi[ri].name);
    if(status!=NULL) sprintf(status, "%d reference curve(s) read", n);

    return n;

  } // reference file was read and processed 


  /* So it's not a file, at least an accessible file, but is it region name? */
  ftype=5; if(filetype!=NULL) *filetype=ftype;
  /* Select ROIs that match the specified input name */
  n=dftSelectRegions(tissue, filename, 1);
  if(verbose>1) printf("nr of ref regions := %d/%d\n", n, tissue->voiNr);
  if(n<=0) {
    if(status!=NULL) sprintf(status, "cannot find region");
    return -7;
  }
  if(n==tissue->voiNr && tissue->voiNr>1) {
    if(status!=NULL) sprintf(status, "all regions do match");
    return -8;
  }

  /* Try to select the best reference ROI */
  ri=dftSelectBestReference(tissue); if(ri<0) {
    if(status!=NULL) sprintf(status, "cannot select the best reference region");
    return -9;
  }
  tissue->voi[ri].sw=2; if(ref_index!=NULL) *ref_index=ri;
  if(verbose>1) printf("selected ref region := %s\n", tissue->voi[ri].name);

  /* Calculate integrals */
  for(ri=0; ri<tissue->voiNr; ri++) if(tissue->voi[ri].sw>0) {
    if(tissue->timetype==DFT_TIME_STARTEND) {
      ret=petintegrate(tissue->x1, tissue->x2, tissue->voi[ri].y,
              tissue->frameNr, tissue->voi[ri].y2, NULL);
    } else {
      ret=integrate(tissue->x, tissue->voi[ri].y, tissue->frameNr,
              tissue->voi[ri].y2);
    }
    if(ret) {
      if(status!=NULL) sprintf(status, "cannot integrate input");
      return(-11);
    }
  }

  if(status!=NULL) sprintf(status, "%d reference curve(s) read", n);
  return(n);
}
/*****************************************************************************/

/*****************************************************************************/
/** Read tissue and input data for modeling. Time units are converted to min
 *  and input calibration units to units of tissue data.
 *  Input data is NOT interpolated to tissue times.
 *  If input data extends much further than fit duration, the last part is 
 *  removed to save computation time in simulations.
 *  Input data is verified not to end too early, and not to start too late.

\return Returns 0 when successful, otherwise a non-zero value. 
 */ 
int dftReadModelingData(
  /** Tissue data filename; one time sample is sufficient here */
  char *tissuefile,
  /** 1st input data filename */
  char *inputfile1,
  /** 2nd input data filename (or NULL if only not needed) */
  char *inputfile2,
  /** 3rd input data filename (or NULL if only not needed) */
  char *inputfile3,
  /** Fit duration (in minutes); shortened if longer than tissue data */
  double *fitdur,
  /** Nr of time frames (samples) in tissue data that are inside fitdur
   *  will be written here */
  int *fitframeNr,
  /** Pointer to initiated DFT into which tissue data will be written */
  DFT *tis,
  /** Pointer to initiated DFT into which input data (plasma and/or blood) 
   *  will be written */
  DFT *inp,
  /** Give file pointer (for example stdout) where log information is printed;
   *  NULL if not needed; warnings will be printed in stderr anyway. */
  FILE *loginfo,
  /** Verbose level; if zero, then only warnings are printed into stderr */
  int verbose,
  /** Pointer to a string (allocated for at least 64 chars) where error message
      or other execution status will be written; enter NULL, if not needed */     
  char *status
) {
  int ret, fi, n, first, last, ii, input_nr=0;
  DFT tmpdft;
  double f, starttime, endtime;
  FILE *wfp;
  char *fname;


  if(loginfo!=NULL && verbose>0) {
    fprintf(loginfo, "dftReadModelingData(\n");
    fprintf(loginfo, "  %s,\n", tissuefile);
    fprintf(loginfo, "  %s,\n", inputfile1);
    fprintf(loginfo, "  %s,\n", inputfile2);
    fprintf(loginfo, "  %s,\n", inputfile3);
    fprintf(loginfo, "  %g,\n", *fitdur);
    fprintf(loginfo, "  *fitframeNr, *tis, *inp, *loginfo, %d, *status\n",
            verbose);
    fprintf(loginfo, ")\n");
  }
  /* Check the input */
  if(status!=NULL) sprintf(status, "program error");
  if(tis==NULL || inp==NULL) return -1;
  if(tissuefile==NULL || strlen(tissuefile)<1) return -2;
  ret=0;
  if(inputfile1==NULL || strlen(inputfile1)<1) return -3; else input_nr++;
  if(inputfile2!=NULL && strlen(inputfile2)>0) input_nr++;
  if(inputfile3!=NULL && strlen(inputfile3)>0) {
    if(input_nr<2) return -4;
    input_nr++;
  }
  if(status!=NULL) sprintf(status, "arguments validated");
  /* Warnings will be printed into stderr */
  wfp=stderr;

  /* Delete any previous data and initiate temp data */
  dftEmpty(inp); dftEmpty(tis); dftInit(&tmpdft);

  /*
   *  Read tissue data
   */
  if(loginfo!=NULL && verbose>0)
    fprintf(loginfo, "reading tissue data in %s\n", tissuefile);
  if(dftRead(tissuefile, tis)) {
    if(status!=NULL)
      sprintf(status, "cannot read '%s': %s", tissuefile, dfterrmsg);
    return(2);
  }
  /* Check for NaN's */
  ret=dft_nr_of_NA(tis); if(ret>0) {
    if(status!=NULL) sprintf(status, "missing sample(s) in %s", tissuefile);
    dftEmpty(tis); return(2);
  }
  /* Sort the data by increasing sample times */
  dftSortByFrame(tis);
  /* Do not check frame number; static scan is fine here */
  /* Make sure that there is no overlap in frame times */
  if(tis->timetype==DFT_TIME_STARTEND) {
    if(loginfo!=NULL && verbose>0)
      fprintf(loginfo, "checking frame overlap in %s\n", tissuefile);
    ret=dftDeleteFrameOverlap(tis);
    if(ret) {
      if(status!=NULL)
        sprintf(status, "%s has overlapping frame times", tissuefile);
      dftEmpty(tis); return(2);
    }
  }
  if(tis->timeunit==TUNIT_UNKNOWN) {
    fprintf(wfp, "Warning: tissue sample time units not known.\n");
  }

  /*
   *  Read first input data
   */
  if(loginfo!=NULL && verbose>0)
    fprintf(loginfo, "reading input data in %s\n", inputfile1);
  if(dftRead(inputfile1, inp)) {
    if(status!=NULL)
      sprintf(status, "cannot read '%s': %s", inputfile1, dfterrmsg);
    dftEmpty(tis); return(3);
  }
  /* Check and correct the sample time unit */
  if(tis->timeunit==TUNIT_UNKNOWN) tis->timeunit=inp->timeunit;
  else if(inp->timeunit==TUNIT_UNKNOWN) inp->timeunit=tis->timeunit;
  if(inp->timeunit==TUNIT_UNKNOWN) {
    fprintf(wfp, "Warning: input sample time units not known.\n");
  }
  if(tis->timeunit!=inp->timeunit &&
     dftTimeunitConversion(inp, tis->timeunit)!=0) {
    if(status!=NULL)
      sprintf(status, "tissue and plasma do have different time units");
    dftEmpty(tis); dftEmpty(inp);
    return(3);
  }
  if(inp->voiNr>1) {
    fprintf(wfp, "Warning: using only first TAC in %s\n", inputfile1);
    inp->voiNr=1;
  }
  if(inp->frameNr<4) {
    if(status!=NULL) sprintf(status, "%s contains too few samples", inputfile1);
    dftEmpty(tis); dftEmpty(inp); return(3);
  }
  /* Sort the data by increasing sample times */
  dftSortByFrame(inp);

  /*
   *  Read following input files, if required
   */
  for(ii=2; ii<=input_nr; ii++) {
    if(ii==2) fname=inputfile2;
    else fname=inputfile3;

    /* Make room for one more curve in input data */
    ret=dftAddmem(inp, 1);
    if(ret) {
      if(status!=NULL) sprintf(status, "cannot allocate more memory");
      dftEmpty(tis); dftEmpty(inp); return(4);
    }

    /* Read blood data */
    dftInit(&tmpdft);
    if(loginfo!=NULL && verbose>0)
      fprintf(loginfo, "reading input data in %s\n", fname);
    if(dftRead(fname, &tmpdft)) {
      if(status!=NULL) sprintf(status, "cannot read '%s': %s", fname, dfterrmsg);
      dftEmpty(tis); dftEmpty(inp); return(4);
    }
    if(tmpdft.frameNr<4) {
      if(status!=NULL) sprintf(status, "%s contains too few samples", fname);
      dftEmpty(tis); dftEmpty(inp); dftEmpty(&tmpdft); return(4);
    }

    /* Check and correct the sample time unit */
    if(tis->timeunit==TUNIT_UNKNOWN) tis->timeunit=tmpdft.timeunit;
    else if(tmpdft.timeunit==TUNIT_UNKNOWN) tmpdft.timeunit=tis->timeunit;
    if(tmpdft.timeunit==TUNIT_UNKNOWN) {
      fprintf(wfp, "Warning: blood sample time units not known.\n");
    }
    if(inp->timeunit!=tmpdft.timeunit &&
       dftTimeunitConversion(&tmpdft, inp->timeunit)!=0) {
      if(status!=NULL)
        sprintf(status, "two input data are in different time units");
      dftEmpty(tis); dftEmpty(inp); dftEmpty(&tmpdft);
      return(4);
    }
    /* Sort the data by increasing sample times */
    dftSortByFrame(&tmpdft);

    /* Copy to input data */
    if(tmpdft.voiNr>1) {
      fprintf(wfp, "Warning: using only first TAC in %s\n", fname);
      tmpdft.voiNr=1;
    }
    if(loginfo!=NULL && verbose>1)
      fprintf(loginfo, "interpolating %d samples into %d samples.\n",
        tmpdft.frameNr, inp->frameNr);
    ret=dftInterpolateInto(&tmpdft, inp, status, verbose);
    if(ret) {
      if(loginfo!=NULL && verbose>0)
        fprintf(loginfo, "dftInterpolateInto() := %d\n", ret);
      dftEmpty(tis); dftEmpty(inp); dftEmpty(&tmpdft); return(4);
    }
    /* Remove the originally read blood data */
    dftEmpty(&tmpdft);

  } // next input file

  /* Set time unit to min */
  if(loginfo!=NULL && verbose>1)
    fprintf(loginfo, "setting time units to min.\n");
  ret=dftTimeunitConversion(tis, TUNIT_MIN); if(ret)
    fprintf(wfp, "Warning: check that regional data times are in minutes.\n");
  ret=dftTimeunitConversion(inp, TUNIT_MIN); if(ret)
    fprintf(wfp, "Warning: check that input data times are in minutes.\n");

  /* 
   *  Check the input data
   */
  if(loginfo!=NULL && verbose>0)
    fprintf(loginfo, "checking input data\n");
  if(inp->frameNr<4) {
    if(status!=NULL) sprintf(status, "%s contains too few samples", inputfile1);
    dftEmpty(tis); dftEmpty(inp); return(4);
  }
  /* Check for NaN's */
  ret=dft_nr_of_NA(inp); if(ret>0) {
    if(status!=NULL) sprintf(status, "missing sample(s) in data");
    dftEmpty(tis); dftEmpty(inp); return(4);
  }
  /* Check that input and tissue time ranges are about the same */
  if(tis->x[tis->frameNr-1]>10.0*inp->x[inp->frameNr-1] ||
     tis->x[tis->frameNr-1]<0.10*inp->x[inp->frameNr-1]) {
    fprintf(wfp, "Warning: you might need to check the sample time units.\n");
  }
  /* Check and give a warning, if the first value is the highest */
  for(fi=1, n=0, f=inp->voi[0].y[0]; fi<inp->frameNr; fi++)
    if(inp->voi[0].y[fi]>=f) {f=inp->voi[0].y[fi]; n=fi;}
  if(n<2) fprintf(wfp, "Warning: check the first input sample values.\n");

  /* Check the tissue and blood TAC concentration units */
  ret=dftUnitConversion(inp, petCunitId(tis->unit));
  if(ret!=0) {
    fprintf(wfp, "Note: check the units of input and tissue data.\n");
  }

  /*
   *  Check and set fit time length
   */
  if(loginfo!=NULL && verbose>0)
    fprintf(loginfo, "checking and setting fit time length\n");
  /* Set fit duration */
  starttime=0; endtime=*fitdur;
  *fitframeNr=fittime_from_dft(tis, &starttime, &endtime,
                               &first, &last, verbose);
  if(loginfo!=NULL && verbose>1) {
    fprintf(loginfo, "tis.frameNr := %d\n", tis->frameNr);
    fprintf(loginfo, "starttime := %g\n", starttime);
    fprintf(loginfo, "endtime := %g\n", endtime);
    fprintf(loginfo, "first := %d\n", first);
    fprintf(loginfo, "last := %d\n", last);
    fprintf(loginfo, "fitframeNr := %d\n", *fitframeNr);
  }
  *fitdur=endtime;
  /* Check that input data does not end much_before fitdur */
  if(inp->timetype==DFT_TIME_STARTEND) f=inp->x2[inp->frameNr-1]; 
  else f=inp->x[inp->frameNr-1];
  if(*fitdur>1.2*f) {
    if(status!=NULL) sprintf(status, "input TAC is too short");
    dftEmpty(inp); dftEmpty(tis); return(5);
  }
  /* Cut off too many input samples to make calculation faster */
  f=*fitdur; 
  if(loginfo!=NULL && verbose>0) printf("Input TAC cutoff at %g min\n", f); 
  for(fi=0; fi<inp->frameNr; fi++) if(inp->x[fi]>f) break;
  if(fi<inp->frameNr) fi++; 
  inp->frameNr=fi;
  if(inp->frameNr<4) {
    if(status!=NULL)
      sprintf(status, "too few samples in specified fit duration.\n");
    dftEmpty(inp); dftEmpty(tis); return(5);
  }
  if(loginfo!=NULL && verbose>1) {
    fprintf(loginfo, "dft.frameNr := %d\ninp.frameNr := %d\nfitdur := %g\n",
    tis->frameNr, inp->frameNr, *fitdur);
    fprintf(loginfo, "fitframeNr := %d\n", *fitframeNr);
  }

  if(status!=NULL) sprintf(status, "ok");
  return 0;
}
/*****************************************************************************/

/*****************************************************************************/
/** Robust search of the min and max values of DFT TAC data.
 *  Data may contain NaN's, and individual outliers are not taken as min or max.
 * @sa dftMinMaxTAC, dftMinMax
 * @return Returns 0 if successful.
 */
int dftRobustMinMaxTAC(
  /** Pointer to the DFT TAC data to search */
  DFT *dft,
  /** Index of the only TAC which is searched for min and max; <0 if all */
  int tacindex,
  /** Pointer to X at TAC min; set to NULL if not needed */
  double *minx,
  /** Pointer to X at TAC max; set to NULL if not needed */
  double *maxx,
  /** Pointer to min Y; set to NULL if not needed */
  double *miny,
  /** Pointer to max Y; set to NULL if not needed */
  double *maxy,
  /** Index of min TAC; set to NULL if not needed */
  int *mini,
  /** Index of max TAC; set to NULL if not needed */
  int *maxi,
  /** Index of min sample; set to NULL if not needed */
  int *mins,
  /** Index of max sample; set to NULL if not needed */
  int *maxs,
  /** Verbose level; if zero, then nothing is printed to stderr or stdout */
  int verbose
) {
  int ri, fi, i1, i2, s1, s2;
  double x, x1, x2, y1, y2;
  int n, runnh, maxrunnh, runns, maxrunns;
  double ym;
  int run1h, run2h, maxrun1h, maxrun2h, run1s, run2s, maxrun1s, maxrun2s;

  if(verbose>0)
    printf("dftRobustMinMaxTAC(dft, %d, minx, maxx, miny, maxy, mini, maxi, mins, maxs, %d)\n",
           tacindex, verbose);

  if(dft==NULL) return(1);
  if(tacindex>=dft->voiNr) return(2);
  if(dft->voiNr<1 || dft->frameNr<1) return(3);

  /* One TAC at a time */
  double list[dft->frameNr];
  x1=x2=y1=y2=nan(""); i1=i2=s1=s2=0;
  for(ri=0; ri<dft->voiNr; ri++) {
    if(tacindex>=0 && ri!=tacindex) continue;
    //fprintf(stderr, "Region: %s\n", dft->voi[ri].name);
    /* Calculate median of y values */
    for(fi=n=0; fi<dft->frameNr; fi++) {
      if(isnan(dft->voi[ri].y[fi])) continue;
      if(dft->timetype==DFT_TIME_STARTEND) {
        if(isnan(dft->x1[fi]) || isnan(dft->x2[fi])) continue;}
      else {if(isnan(dft->x[fi])) continue;}
      list[n++]=dft->voi[ri].y[fi];
    }
    if(n<10) {
      maxrunnh=maxrunns=dft->frameNr;
      maxrun1h=maxrun1s=0;
      maxrun2h=maxrun2s=dft->frameNr-1;
    } else {
      ym=dmedian(list, n);
      //fprintf(stderr, "median := %g\n", ym);
      /* Search the longest run of values which are larger/smaller than
         the median */
      runnh=maxrunnh=run1h=run2h=maxrun1h=maxrun2h=0;
      runns=maxrunns=run1s=run2s=maxrun1s=maxrun2s=0;
      for(fi=0; fi<dft->frameNr; fi++) {
        if(isnan(dft->voi[ri].y[fi])) continue;
        if(dft->timetype==DFT_TIME_STARTEND) {
          if(isnan(dft->x1[fi]) || isnan(dft->x2[fi])) continue;}
        else {if(isnan(dft->x[fi])) continue;}
        /* max */
        if(dft->voi[ri].y[fi]>ym) {
          runnh++; if(runnh==1) {run1h=run2h=fi;} else {run2h=fi;}
        } else {
          if(runnh>maxrunnh) {maxrunnh=runnh; maxrun1h=run1h; maxrun2h=run2h;}
          runnh=0;
        }
        /* min */
        if(dft->voi[ri].y[fi]<ym) {
          runns++; if(runns==1) {run1s=run2s=fi;} else {run2s=fi;}
        } else {
          if(runns>maxrunns) {maxrunns=runns; maxrun1s=run1s; maxrun2s=run2s;}
          runns=0;
        }
      }
      /* Get range also from the end part */
      if(runnh>maxrunnh) {maxrunnh=runnh; maxrun1h=run1h; maxrun2h=run2h;}
      if(runns>maxrunns) {maxrunns=runns; maxrun1s=run1s; maxrun2s=run2s;}
      /* If longest range includes just one sample, then take the whole range */
      if(maxrunnh<2) {maxrunnh=dft->frameNr; maxrun1h=0; maxrun2h=dft->frameNr-1;}
      if(maxrunns<2) {maxrunns=dft->frameNr; maxrun1s=0; maxrun2s=dft->frameNr-1;}
    }
    if(verbose>12) {
      fprintf(stderr, "longest run for max: %g - %g\n",
        dft->x[maxrun1h], dft->x[maxrun2h]);
      fprintf(stderr, "longest run for min: %g - %g\n",
        dft->x[maxrun1s], dft->x[maxrun2s]);
    }
    /* Inside the range, search for max */
    for(fi=maxrun1h; fi<=maxrun2h; fi++) {
      if(isnan(dft->voi[ri].y[fi])) continue;
      if(dft->timetype==DFT_TIME_STARTEND) {
        if(isnan(dft->x1[fi]) || isnan(dft->x2[fi])) continue;
        x=0.5*(dft->x1[fi]+dft->x2[fi]);
      } else {
        if(isnan(dft->x[fi])) continue;
        x=dft->x[fi];
      }
      if(isnan(y2) || y2<dft->voi[ri].y[fi]) {
        y2=dft->voi[ri].y[fi]; i2=ri; x2=x; s2=fi;}
    }
    /* Inside the range, search for min */
    for(fi=maxrun1s; fi<=maxrun2s; fi++) {
      if(isnan(dft->voi[ri].y[fi])) continue;
      if(dft->timetype==DFT_TIME_STARTEND) {
        if(isnan(dft->x1[fi]) || isnan(dft->x2[fi])) continue;
        x=0.5*(dft->x1[fi]+dft->x2[fi]);
      } else {
        if(isnan(dft->x[fi])) continue;
        x=dft->x[fi];
      }
      if(isnan(y1) || y1>dft->voi[ri].y[fi]) {
        y1=dft->voi[ri].y[fi]; i1=ri; x1=x; s1=fi;}
    }

  } /* next TAC */

  if(minx!=NULL) {if(isnan(x1)) return(11); else *minx=x1;}
  if(maxx!=NULL) {if(isnan(x2)) return(12); else *maxx=x2;}
  if(miny!=NULL) {if(isnan(y1)) return(13); else *miny=y1;}
  if(maxy!=NULL) {if(isnan(y2)) return(14); else *maxy=y2;}
  if(mini!=NULL) {if(isnan(y1)) return(13); else *mini=i1;}
  if(maxi!=NULL) {if(isnan(y2)) return(14); else *maxi=i2;}
  if(mins!=NULL) {if(isnan(y1)) return(13); else *mins=s1;}
  if(maxs!=NULL) {if(isnan(y2)) return(14); else *maxs=s2;}
  return(0);
}
/*****************************************************************************/

/*****************************************************************************/
/** Verify that specified (input) TAC contains peak and that it does not
    start too late to get reliable estimate of AUC.
    @return Returns 0 in case data is suitable, < 0 if there may be a problem,
    1 in case of an error, and 2 if peak seems to be missed. 
 */
int dftVerifyPeak(
  /** Pointer to TAC data which is verified */
  DFT *dft,
  /** Index of TAC to be verified; <0 in case all are (separately) verified */
  int index,
  /** Verbose level; if zero, then nothing is printed to stderr or stdout */
  int verbose,
  /** Pointer to a string (allocated for at least 64 chars) where error message
      or other execution status will be written; enter NULL, if not needed */     
  char *status
) {
  int ri, fi, ret, mini, maxi, starti, warn=0;
  double minx, maxx, miny, maxy, startx, dif, lowesty;

  if(verbose>0) printf("dftVerifyPeak(dft, %d, %d)\n", index, verbose);
  /* Check the arguments */
  if(status!=NULL) strcpy(status, "program error");
  if(dft==NULL || dft->frameNr<1 || dft->voiNr<1) return 1;
  if(index>=dft->voiNr) return 1;

  /* If less than 3 samples, then this can not be suitable as input TAC */
  if(dft->frameNr<3) {
    if(status!=NULL) strcpy(status, "too few samples");
    return 2;
  }

  /* Make sure data is sorted by increasing sample time */
  dftSortByFrame(dft);

  /* Go through all TACs that are requested */
  for(ri=0; ri<dft->voiNr; ri++) {
    if(index>=0 && index!=ri) continue;
    if(verbose>1) printf("checking region %d: %s\n", 1+ri, dft->voi[ri].name);

    /* Get the TAC min and max */
    ret=dftMinMaxTAC(dft, ri, &minx, &maxx, &miny, &maxy, NULL, NULL,
                     &mini, &maxi);
    if(ret!=0) {
      if(verbose>0) printf("Error %d in dftMinMaxTAC()\n", ret);
      if(status!=NULL) strcpy(status, "invalid TAC");
      return 1;
    }

    /* Check that there are no big negative values */
    if(maxy<=0.0) {
      if(verbose>0) printf("TAC has no positive values.\n");
      if(status!=NULL) strcpy(status, "no positive TAC values");
      return 2;
    }
    if(miny<0.0) {
      if((-miny)>0.40*maxy) {
        if(verbose>0) printf("TAC has high negative value(s).\n");
        if(status!=NULL) strcpy(status, "too high negative TAC values");
        return 2;
      }
      if((-miny)>0.02*maxy) {
        if(verbose>1) printf("TAC has negative value(s).\n");
        warn++;
      }
    }

    /* Get the first sample time, considering missing values */
    startx=1.0E+10; starti=dft->frameNr-1;
    for(fi=0; fi<dft->frameNr; fi++) {
      if(isnan(dft->voi[ri].y[fi])) continue;
      if(dft->timetype==DFT_TIME_STARTEND) {
        if(isnan(dft->x1[fi])) continue;
        if(isnan(dft->x2[fi])) continue;
        startx=dft->x1[fi]; starti=fi; break;
      } else {
        if(isnan(dft->x[fi])) continue;
        startx=dft->x[fi]; starti=fi; break;
      }
    }
    if(verbose>2) printf("first time sample at %g\n", startx);

    /* If the first sample is the peak, then check that sample time is not too
       late */
    if(maxi==starti) {
      if(verbose>2) printf("Peak at the first sample.\n");
      if(status!=NULL) strcpy(status, "input TAC should start at time zero");
      /* Calculate peak time to initial sample frequency ratio */
      if(dft->timetype==DFT_TIME_STARTEND) {
        dif=dft->x1[maxi]/(dft->x2[maxi]-dft->x1[maxi]);
      } else {
        //dif=dft->x[maxi]/(dft->x[maxi+1]-dft->x[maxi]);
        for(fi=maxi+1; fi<dft->frameNr; fi++)
          if(!isnan(dft->x[fi]) && dft->x[fi]>dft->x[maxi]) break;
        if(fi==dft->frameNr) dif=1.0E+10;
        else dif=dft->x[maxi]/(dft->x[fi]-dft->x[maxi]);
      }
      if(dif>0.3) {
        if(verbose>0) printf("Peak at the first sample which is not at zero.\n");
        if(verbose>1) printf("dif := %g\n", dif);
      }
      if(dif>5.0) {
        return 2;
      } else if(dif>1.0) {
        /* Not bad, if peak is still much higher than tale */
        if(maxy>20.*miny) warn++; else return 2;
        if(verbose>1) printf("good peak/tail -ratio\n");
      } else if(dif>0.3) warn++;
    }

    /* Search the lowest value before the peak */
    lowesty=dft->voi[ri].y[starti];
    for(fi=starti+1; fi<maxi; fi++)
      if(dft->voi[ri].y[fi]<lowesty) lowesty=dft->voi[ri].y[fi];
    if(verbose>2) printf("lowest value before peak: %g\n", lowesty);

    /* If first sample is closer to peak and relatively high, then that is
       or may be a problem */
    if(maxi>starti && startx>0.001 && startx>0.75*maxx) {
      if(dft->voi[ri].y[starti]>0.66*maxy && lowesty>0.05*maxy && mini>maxi) {
        if(verbose>0) printf("The first sample is relatively late and high.\n");
        if(status!=NULL) strcpy(status, "input TAC should start at time zero");
        if(verbose>2) {
          printf("startx=%g\n", startx);
          printf("starty=%g\n", dft->voi[ri].y[starti]);
          printf("maxx=%g\n", maxx);
          printf("maxy=%g\n", maxy);
        }
        return 2;
      }
    }
    if(maxi>starti && startx>0.001 && startx>0.5*maxx) {
      if(dft->voi[ri].y[starti]>0.5*maxy && lowesty>0.05*maxy && mini>maxi) {
        if(verbose>1) printf("The first sample is relatively late and high.\n");
        warn++;
      }
    }
    if(verbose>5) {
      printf("startx=%g\n", startx);
      printf("starty=%g\n", dft->voi[ri].y[starti]);
      printf("maxx=%g\n", maxx);
      printf("maxy=%g\n", maxy);
    }

    /* If peak is not much higher than lowest value, that may indicate
       a problem */
    if(maxy<1.5*miny) {
      if(verbose>0) printf("TAC does not have a clear peak.\n");
      if(status!=NULL) strcpy(status, "input TAC peak missing");
      return 2;
    }
    if(maxy<5.0*miny) {
      if(verbose>1) printf("TAC does not have a clear peak.\n");
      warn++;
    }

  } // next TAC

  if(verbose>0 && warn>0) printf("%d warning(s)\n", warn);
  if(warn>0) {
    if(status!=NULL) strcpy(status, "input TAC is not optimal");
    return -1;
  }
  if(status!=NULL) strcpy(status, "input TAC ok");
  return 0;
}
/*****************************************************************************/

/*****************************************************************************/
