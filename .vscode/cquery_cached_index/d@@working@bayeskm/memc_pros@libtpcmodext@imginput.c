/// @file imginput.c
/// @brief Procedures for handling model input data.
/// @author Vesa Oikonen
///
/*****************************************************************************/

/*****************************************************************************/
#include "libtpcmodext.h"
/*****************************************************************************/

/*****************************************************************************/
/** Read PET image and input TAC for modelling.

    Input calibration units are converted to units of tissue data.
    Input data is optionally interpolated to tissue times.
    If input data extends much further than fit duration, the last part is 
    removed to save computation time in simulations.
    Input data is optionally verified not to end too early, and not to
    start too late.

    @return Returns 0 when successful, otherwise a non-zero value. 
 */ 
int imgReadModelingData(
  /** PET dynamic image filename; one time frame is sufficient here */
  char *petfile,
  /** Optional SIF. SIF is automatically read with NIfTI and Analyze formats,
   *  if found in the same folder and files are consistently named.
   *  SIF filename needs to be specified here only if that is not the case,
   *  or if other image formats do not contain information on frame times or
   *  tracer isotope, or if weights from SIF are needed.
   *  Enter NULL if not needed. */
  char *siffile,
  /** 1st input data filename */
  char *inputfile1,
  /** 2nd input data filename (or NULL if only not needed) */
  char *inputfile2,
  /** 3rd input data filename (or NULL if only not needed) */
  char *inputfile3,
  /** Fit duration (in minutes); shortened if longer than PET data;
   *  enter <=0 to use all data in image; input is verified to not be much
   *  shorter than fitfur. */
  double *fitdur,
  /** Nr of time frames (samples) in PET data that are inside fitdur
   *  will be written here; pointer MUST be provided. */
  int *fitframeNr,
  /** Pointer to initiated IMG struct into which PET data will be written */
  IMG *img,
  /** Pointer to initiated DFT into which input data (plasma, blood, and/or
   *  reference tissue TAC) will be written without interpolation;
   *  enter NULL if not needed. Times will be in seconds, like in image. */
  DFT *inp,
  /** Pointer to initiated DFT into which input data (plasma, blood, and/or
   *  reference tissue TAC) will be written, interpolated to PET frames;
   *  enter NULL if not needed.  Times will be in seconds, like in image. */
  DFT *iinp,
  /** Set to <>0 to verify that input TAC does not start too late and has
   *  decent peak to provide reliable AUC0-T. */
  int verifypeak,
  /** Give file pointer (for example stdout) where log information is printed;
   *  NULL if not needed; warnings will be printed in stderr anyway. */
  FILE *loginfo,
  /** Verbose level; if zero, then only warnings are printed into stderr */
  int verbose,
  /** Pointer to a string (allocated for at least 256 chars) where error message
      or other execution status will be written; enter NULL, if not needed */     
  char *status
) {
  int ret, fi, n, first, last, ii, input_nr=0;
  DFT tmpdft, dft;
  SIF sif;
  double f, starttime, endtime;
  FILE *wfp;
  char *fname, tmp[512];


  if(loginfo!=NULL && verbose>0) {
    fprintf(loginfo, "imgReadModelingData(\n");
    fprintf(loginfo, "  %s,\n", petfile);
    fprintf(loginfo, "  %s,\n", inputfile1);
    fprintf(loginfo, "  %s,\n", inputfile2);
    fprintf(loginfo, "  %s,\n", inputfile3);
    fprintf(loginfo, "  %g,\n", *fitdur);
    fprintf(loginfo, "  *fitframeNr, *tis, *inp, *iinp *loginfo, %d, *status\n", verbose);
    fprintf(loginfo, ")\n");
  }
  /* Check the function input */
  if(status!=NULL) sprintf(status, "program error");
  if(img==NULL || (inp==NULL && iinp==NULL)) return -1;
  if(petfile==NULL || strlen(petfile)<1) return -2;
  ret=0;
  if(inputfile1==NULL || strlen(inputfile1)<1) return -3; else input_nr++;
  if(inputfile2!=NULL && strlen(inputfile2)>0) input_nr++;
  if(inputfile3!=NULL && strlen(inputfile3)>0) {
    if(input_nr<2) return -4;
    input_nr++;
  }
  if(fitframeNr==NULL || fitdur==NULL) return -5;
  if(status!=NULL) sprintf(status, "arguments validated");
  /* Warnings will be printed into stderr */
  wfp=stderr;

  /* Check that input file(s) exist, because image is read first, and it
     might take a while */
  if(loginfo!=NULL && verbose>1) fprintf(loginfo, "checking access to input files\n");
  if(access(inputfile1, 0) == -1) {
    sprintf(tmp, "cannot read '%s'", inputfile1);
    if(status!=NULL) strcpy(status, tmp);
    else if(loginfo!=NULL) fprintf(loginfo, "Error: %s\n", tmp);
    else fprintf(wfp, "Error: %s\n", tmp);
    return(1);
  }
  for(ii=2; ii<=input_nr; ii++) {
    if(ii==2) fname=inputfile2; else fname=inputfile3;
    if(access(fname, 0) == -1) {
      sprintf(tmp, "cannot read '%s'", fname);
      if(status!=NULL) strcpy(status, tmp);
      else if(loginfo!=NULL) fprintf(loginfo, "Error: %s\n", tmp);
      else fprintf(wfp, "Error: %s\n", tmp);
      return(1);
    }
  }

  /* Delete any previous data and initiate temp data */
  dftEmpty(inp); dftEmpty(iinp); imgEmpty(img);


  /*
   *  Read PET image
   */
  if(loginfo!=NULL && verbose>0) fprintf(loginfo, "reading image %s\n", petfile);
  imgInit(img);
  ret=imgRead(petfile, img); fflush(stderr); fflush(stdout);
  if(ret) {
    sprintf(tmp, "cannot read '%s': %s", petfile, img->statmsg);
    if(status!=NULL) strcpy(status, tmp);
    else if(loginfo!=NULL) fprintf(loginfo, "Error: %s\n", tmp);
    else fprintf(wfp, "Error: %s\n", tmp);
    if(loginfo!=NULL) fprintf(loginfo, "imgRead(%s, img) := %d\n", petfile, ret);
    return(2);
  }
  /* Check if PET data is raw or image */
  if(img->type!=IMG_TYPE_IMAGE) {
    sprintf(tmp, "%s is not an image.\n", petfile);
    if(status!=NULL) strcpy(status, tmp);
    else if(loginfo!=NULL) fprintf(loginfo, "Error: %s\n", tmp);
    else fprintf(wfp, "Error: %s\n", tmp);
    imgEmpty(img); return(3);
  }
  if(loginfo!=NULL && verbose>2)
    fprintf(loginfo, "image contains %d frames and %d planes.\n", img->dimt, img->dimz);

  /* Read SIF, if specified */
  if(siffile!=NULL && siffile[0]) {
    if(loginfo!=NULL && verbose>0) fprintf(loginfo, "reading SIF %s\n", siffile);
    sifInit(&sif);
    ret=sifRead(siffile, &sif);
    if(ret) {
      sprintf(tmp, "cannot read '%s': %s", siffile, siferrmsg);
      if(status!=NULL) strcpy(status, tmp);
      else if(loginfo!=NULL) fprintf(loginfo, "Error: %s\n", tmp);
      else fprintf(wfp, "Error: %s\n", tmp);
      imgEmpty(img); return(4);
    }
    /* Get isotope from SIF, if available */
    {
      double hl;
      hl=hlFromIsotope(sif.isotope_name);
      if(hl>0.0) {
        img->isotopeHalflife=60.0*hl;
        if(loginfo!=NULL && verbose>0) fprintf(loginfo, "isotope code read from %s\n", siffile);
      }
    }
    /* Get frame times from SIF if not available in image */
    if(!imgExistentTimes(img) && img->dimt<=sif.frameNr) {
      for(fi=0; fi<img->dimt; fi++) {
        img->start[fi]=sif.x1[fi]; img->end[fi]=sif.x2[fi];
        img->mid[fi]=0.5*(sif.x1[fi]+sif.x2[fi]);
      }
      if(loginfo!=NULL && verbose>0) fprintf(loginfo, "image frame times read from %s\n", siffile);
    }
    sifEmpty(&sif);
  }

  /* Check that image frame times are available */
  if(loginfo!=NULL && verbose>1) fprintf(loginfo, "checking image contents\n");
  if(!imgExistentTimes(img)) {
    strcpy(tmp, "image frame times not available");
    if(status!=NULL) strcpy(status, tmp);
    else if(loginfo!=NULL) fprintf(loginfo, "Error: %s\n", tmp);
    else fprintf(wfp, "Error: %s\n", tmp);
    imgEmpty(img); return(5);
  }
  /* Make sure that there is no overlap in image frames */ 
  if(loginfo!=NULL && verbose>1) fprintf(loginfo, "checking frame overlap in %s\n", petfile);
  ret=imgDeleteFrameOverlap(img);
  if(ret) {
    strcpy(tmp, "image has overlapping frame times");
    if(status!=NULL) strcpy(status, tmp);
    else if(loginfo!=NULL) fprintf(loginfo, "Error: %s\n", tmp);
    else fprintf(wfp, "Error: %s\n", tmp);
    imgEmpty(img); return(5);
  }
  /* Check fit end time with the PET image */
  if(*fitdur<=0.0) *fitdur=1.0E+99;
  else if(*fitdur<1.0E+10) *fitdur*=60.0; // fit time must be in sec
  *fitframeNr=fittime_from_img(img, fitdur, verbose-2);
  if(*fitframeNr<1 || *fitdur<=3.5) {
    strcpy(tmp, "image has no data in fit time range");
    if(status!=NULL) strcpy(status, tmp);
    else if(loginfo!=NULL) fprintf(loginfo, "Error: %s\n", tmp);
    else fprintf(wfp, "Error: %s\n", tmp);
    imgEmpty(img); return(5);
  }
  *fitdur/=60.0; // back to minutes
  if(loginfo!=NULL && verbose>3) {
    fprintf(loginfo, "fittimeFinal := %g min\n", *fitdur);
    fprintf(loginfo, "fitframeNr := %d\n", *fitframeNr);
  }


  /*
   *  Read first input data
   */
  if(loginfo!=NULL && verbose>0) fprintf(loginfo, "reading input data in %s\n", inputfile1);
  dftInit(&dft);
  if(dftRead(inputfile1, &dft)) {
    sprintf(tmp, "cannot read '%s': %s", inputfile1, dfterrmsg);
    if(status!=NULL) strcpy(status, tmp);
    else if(loginfo!=NULL) fprintf(loginfo, "Error: %s\n", tmp);
    else fprintf(wfp, "Error: %s\n", tmp);
    imgEmpty(img); return(11);
  }
  /* Check TAC nr */
  if(dft.voiNr>1) {
    if(verbose>0) {
      strcpy(tmp, "only first TAC is used as input");
      if(loginfo!=NULL) fprintf(loginfo, "Warning: %s.\n", tmp);
      else fprintf(wfp, "Warning: %s.\n", tmp);
    }
    dft.voiNr=1;
  }
  /* check if file contains NAs (missing values) */
  if(dft_nr_of_NA(&dft)>0) {
    strcpy(tmp, "missing values in input file");
    if(status!=NULL) strcpy(status, tmp);
    else if(loginfo!=NULL) fprintf(loginfo, "Error: %s.\n", tmp);
    else fprintf(wfp, "Error: %s.\n", tmp);
    imgEmpty(img); dftEmpty(&dft); return(12);
  }
  /* Check the concentration units */
  ret=cunit_check_dft_vs_img(&dft, img, tmp, verbose-2);
  if(ret==0) {
    if(loginfo!=NULL && verbose>3) fprintf(loginfo, "%s\n", tmp);
  } else if(ret<0) {
    fprintf(wfp, "Warning: %s\n", tmp);
  } else {
    if(status!=NULL) strcpy(status, tmp);
    else if(loginfo!=NULL) fprintf(loginfo, "Error: %s.\n", tmp);
    else fprintf(wfp, "Error: %s.\n", tmp);
    imgEmpty(img); dftEmpty(&dft); return(13);
  }
  if(verbose>3)
    printf("input time range := %g - %g %s\n",
      dft.x[0], dft.x[dft.frameNr-1], petTunit(dft.timeunit));
  /* Set time unit to sec like in IMG */
  if(dft.timeunit==TUNIT_UNKNOWN) {
    if(dftEndtime(&dft)<0.2*imgEndtime(img)) dft.timeunit=TUNIT_MIN;
    else dft.timeunit=TUNIT_SEC;
    fprintf(wfp, "Warning: assuming that times are in %s in %s\n",
            petTunit(dft.timeunit), inputfile1);
  }
  /* Set time unit to sec like in IMG */
  ret=dftTimeunitConversion(&dft, TUNIT_SEC);
  if(verbose>1)
    printf("input time range := %g - %g %s\n",
           dft.x[0], dft.x[dft.frameNr-1], petTunit(dft.timeunit));
  
  /* Verify the peak if requested; only for the first input TAC, because
     for example metabolite TAC may not show bolus shape */
  if(verifypeak!=0) {
    if(loginfo!=NULL && verbose>1) fprintf(loginfo, "verifying input peak\n");
    ret=dftVerifyPeak(&dft, 0, verbose-5, tmp);
    if(ret>0) {
      if(status!=NULL) strcpy(status, tmp);
      else if(loginfo!=NULL) fprintf(loginfo, "Error: %s.\n", tmp);
      else fprintf(wfp, "Error: %s.\n", tmp);
      imgEmpty(img); dftEmpty(&dft); return(14);
    }
  }
  if(verbose>5)
    printf("input time range := %g - %g %s\n",
           dft.x[0], dft.x[dft.frameNr-1], petTunit(dft.timeunit));


  /*
   *  Read following input files, if required
   */
  dftInit(&tmpdft);
  for(ii=2; ii<=input_nr; ii++) {
    if(ii==2) fname=inputfile2; else fname=inputfile3;

    /* Make room for one more curve in input data */
    ret=dftAddmem(&dft, 1);
    if(ret) {
      strcpy(tmp, "cannot allocate more memory");
      if(status!=NULL) strcpy(status, tmp);
      else if(loginfo!=NULL) fprintf(loginfo, "Error: %s\n", tmp);
      else fprintf(wfp, "Error: %s\n", tmp);
      imgEmpty(img); dftEmpty(&dft); return(15);
    }
    /* Read input TAC */
    if(loginfo!=NULL && verbose>0) fprintf(loginfo, "reading input data in %s\n", fname);
    if(dftRead(fname, &tmpdft)) {
      sprintf(tmp, "cannot read '%s': %s", fname, dfterrmsg);
      if(status!=NULL) strcpy(status, tmp);
      else if(loginfo!=NULL) fprintf(loginfo, "Error: %s\n", tmp);
      else fprintf(wfp, "Error: %s\n", tmp);
      imgEmpty(img); dftEmpty(&dft); return(16);
    }
    /* Check TAC nr */
    if(tmpdft.voiNr>1) {
      if(verbose>0) {
        strcpy(tmp, "only first TAC is used as input");
        if(loginfo!=NULL) fprintf(loginfo, "Warning: %s.\n", tmp);
        else fprintf(wfp, "Warning: %s.\n", tmp);
      }
      dft.voiNr=1;
    }
    /* check if file contains NAs (missing values) */
    if(dft_nr_of_NA(&tmpdft)>0) {
      strcpy(tmp, "missing values in input file");
      if(status!=NULL) strcpy(status, tmp);
      else if(loginfo!=NULL) fprintf(loginfo, "Error: %s\n", tmp);
      else fprintf(wfp, "Error: %s\n", tmp);
      imgEmpty(img); dftEmpty(&dft); dftEmpty(&tmpdft); return(17);
    }
    /* Check and correct the sample time unit */
    if(tmpdft.timeunit==TUNIT_UNKNOWN) {
      if(dftEndtime(&tmpdft)<0.2*imgEndtime(img)) tmpdft.timeunit=TUNIT_MIN;
      else tmpdft.timeunit=TUNIT_SEC;
      fprintf(wfp, "Warning: assuming that times are in %s in %s\n",
              petTunit(tmpdft.timeunit), fname);
    }
    dftTimeunitConversion(&tmpdft, dft.timeunit);
    /* Check the concentration units */
    ret=cunit_check_dft_vs_img(&tmpdft, img, tmp, verbose-2);
    if(ret==0) {
      if(loginfo!=NULL && verbose>3) fprintf(loginfo, "%s\n", tmp);
    } else if(ret<0) {
      fprintf(wfp, "Warning: %s\n", tmp);
    } else {
      if(status!=NULL) strcpy(status, tmp);
      else if(loginfo!=NULL) fprintf(loginfo, "Error: %s\n", tmp);
      else fprintf(wfp, "Error: %s\n", tmp);
      imgEmpty(img); dftEmpty(&dft); dftEmpty(&tmpdft); return(18);
    }
    /* Copy to input data */
    if(loginfo!=NULL && verbose>1)
      fprintf(loginfo, "interpolating %d samples into %d samples.\n",
        tmpdft.frameNr, dft.frameNr);
    ret=dftInterpolateInto(&tmpdft, &dft, tmp, verbose);
    if(ret) {
      if(loginfo!=NULL) fprintf(loginfo, "dftInterpolateInto() := %d\n", ret);
      if(status!=NULL) strcpy(status, tmp);
      else if(loginfo!=NULL) fprintf(loginfo, "Error: %s\n", tmp);
      else fprintf(wfp, "Error: %s\n", tmp);
      imgEmpty(img); dftEmpty(&dft); dftEmpty(&tmpdft); return(19);
    }
    /* Remove the originally read input data */
    dftEmpty(&tmpdft);

  } // next input file
  if(verbose>10)
    printf("input time range := %g - %g %s\n",
           dft.x[0], dft.x[dft.frameNr-1], petTunit(dft.timeunit));


  /*
   *  Check and set input data range vs fit time length
   */
  if(loginfo!=NULL && verbose>0) fprintf(loginfo, "checking and setting input sample time range\n");
  /* Check that input data does not end _much_ before fitdur; 
     dftInterpolateForIMG() below tests this again. */
  starttime=0; endtime=*fitdur; // start and end times always in min
  n=fittime_from_dft(&dft, &starttime, &endtime, &first, &last, verbose-1);
  if(loginfo!=NULL && verbose>2) {
    fprintf(loginfo, "starttime := %g min\n", starttime);
    fprintf(loginfo, "endtime := %g min\n", endtime);
    fprintf(loginfo, "sample_nr := %d\n", n);
  }
  if(*fitdur>1.5*endtime && (*fitdur-endtime)>0.15) {
    strcpy(tmp, "input TAC is too short");
    if(status!=NULL) strcpy(status, tmp);
    else if(verbose>0) fprintf(wfp, "Error: %s\n", tmp);
    imgEmpty(img); dftEmpty(&dft); return(21);
  }
#if(0)
  /* Do NOT test this here! It would prevent using this function in
     reading steady-state studies etc */ 
  if(n<4) {
    strcpy(tmp, "too few input samples in specified fit duration");
    if(status!=NULL) strcpy(status, tmp);
    else if(loginfo!=NULL) fprintf(loginfo, "Error: %s\n", tmp);
    else fprintf(wfp, "Error: %s\n", tmp);
    imgEmpty(img); dftEmpty(&dft); return(22);
  }
#endif
  /* Cut off too many input samples to make calculation faster */
  if(verbose>10) {
    printf("input time range := %g - %g %s\n",
           dft.x[0], dft.x[dft.frameNr-1], petTunit(dft.timeunit));
    printf("fitdur := %g min\n", *fitdur);
  }
  f=*fitdur*60.0;
  if(dftEndtime(&dft)>f) {
    if(loginfo!=NULL && verbose>0) fprintf(loginfo, "Input TAC cutoff at %g sec\n", f); 
    for(fi=0; fi<dft.frameNr; fi++) if(dft.x[fi]>f) break;
    if(fi<dft.frameNr) fi++; 
    dft.frameNr=fi;
  }
  if(verbose>10)
    printf("input time range := %g - %g %s\n",
           dft.x[0], dft.x[dft.frameNr-1], petTunit(dft.timeunit));

  /* Copy input data to given pointer, if required */
  if(inp!=NULL) {
    if(loginfo!=NULL && verbose>0) fprintf(loginfo, "copying input TAC\n");
    ret=dftdup(&dft, inp);
    if(ret!=0) {
      strcpy(tmp, "cannot copy TAC contents");
      if(status!=NULL) strcpy(status, tmp);
      else if(loginfo!=NULL) fprintf(loginfo, "Error: %s.\n", tmp);
      else fprintf(wfp, "Error: %s\n", tmp);
      imgEmpty(img); dftEmpty(&dft); return(31);
    }
  }
    
  /* Interpolate input data to given pointer, if required */
  if(iinp!=NULL) {
    if(loginfo!=NULL && verbose>0) fprintf(loginfo, "interpolating input TAC\n");
    ret=dftInterpolateForIMG(&dft, img, *fitframeNr, iinp, NULL, NULL, verbose, tmp);
    if(ret!=0) {
      if(status!=NULL) strcpy(status, tmp);
      else if(loginfo!=NULL) fprintf(loginfo, "Error: %s.\n", tmp);
      else fprintf(wfp, "Error: %s\n", tmp);
      imgEmpty(img); dftEmpty(&dft); return(33);
    }
  }

  dftEmpty(&dft);
  if(status!=NULL) sprintf(status, "ok");
  return(0);
}
/*****************************************************************************/

/*****************************************************************************/
