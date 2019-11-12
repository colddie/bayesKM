/// @file misc_model.c
/// @brief Miscellaneous functions for PET modelling.
/// @author Vesa Oikonen
///
/*****************************************************************************/

/*****************************************************************************/
#include "libtpcmodext.h"
/*****************************************************************************/

/*****************************************************************************/
/** Interpolate (and integrate) TAC data to sample times that are given with IMG data.
    Input frame lengths are taken into account if they are available in DFT
    or if the framing is the same as in IMG data.
   @return Returns 0 if successful, and <>0 in case of an error.
 */ 
int dftInterpolateForIMG(
  /** Data which is interpolated. */
  DFT *input,
  /** Data to which (sample times) the interpolation is done. */
  IMG *img,
  /** Number of IMG frames that are needed; can be set to 0 if all frames can be included. */
  int frame_nr,
  /** Pointer to initiated DFT into which interpolated values and integrals will be written at 
      input sample times and units. */
  DFT *output,
  /** First time of input data before interpolation (in input time units);
      use this to check that required time range was measured; NULL if not needed. */
  double *ti1,
  /** Last time of input data before interpolation (in input time units);
      use this to check that required time range was measured; NULL if not needed. */
  double *ti2,
  /** Verbose level; set to zero to not to print any comments. */
  int verbose,
  /** Pointer to a string (allocated for at least 64 chars) where error message
      or other execution status will be written; enter NULL, if not needed. */
  char *status
) {
  /* Check input */
  if(verbose>0) {
    printf("%s(*inp, *img, %d, *out, *ti1, *ti2, %d, status)\n", __func__, frame_nr, verbose);
    fflush(stdout);
  }
  if(input==NULL || img==NULL || output==NULL) {
    if(status!=NULL) strcpy(status, "program error");
    return 1;
  }
  if(input->frameNr<1 || input->voiNr<1 || img->dimt<1) {
    if(status!=NULL) strcpy(status, "no pet data");
    return 2;
  }
  if(input->timeunit!=TUNIT_MIN && input->timeunit!=TUNIT_SEC) {
    if(status!=NULL) strcpy(status, "unknown time units");
    return 3;
  }
  if(frame_nr<1 || frame_nr>img->dimt) frame_nr=img->dimt;
  /* Set ti1 and ti2 */
  if(ti1!=NULL) {
    if(input->timetype==DFT_TIME_STARTEND)
      *ti1=input->x1[0]; else *ti1=input->x[0];
  }
  if(ti2!=NULL) {
    if(input->timetype==DFT_TIME_STARTEND) *ti2=input->x2[input->frameNr-1];
    else *ti2=input->x[input->frameNr-1];
  }

  /* Delete any previous data */
  dftEmpty(output);

  /* Allocate memory for interpolated data */
  if(verbose>10) printf("allocating memory for interpolated data\n");
  if(dftAllocateWithHeader(output, frame_nr, input->voiNr, input)) {
    if(status!=NULL) strcpy(status, "memory allocation error");
    return 11;
  }
  output->voiNr=input->voiNr; output->frameNr=frame_nr;

  /* Set output times */
  if(copy_times_from_img_to_dft(img, output, verbose)!=0) {
    if(status!=NULL) strcpy(status, "frame time error");
    dftEmpty(output); return 12;
  }
  if(verbose>10) {
    printf("time range := %g - %g %s\n", output->x[0],
           output->x[output->frameNr-1], petTunit(output->timeunit));
    printf("  timetype := %d\n", output->timetype);
  }

  // Check if input and output data do have the same frame times already
  if(check_times_dft_vs_dft(input, output, verbose)==1 && input->frameNr>=output->frameNr) {
    // copy the values directly and integrate in place
    if(verbose>10) printf("frame times are assumed to be the same\n");
    int ret=0;
    for(int ri=0; ri<output->voiNr && ret==0; ri++) {
      for(int fi=0; fi<output->frameNr; fi++)
        output->voi[ri].y[fi]=input->voi[ri].y[fi];
      ret=petintegral(output->x1, output->x2, output->voi[ri].y,
             output->frameNr, output->voi[ri].y2, output->voi[ri].y3);
    }
    if(ret) {
      if(status!=NULL) sprintf(status, "cannot interpolate (%d)", ret);
      dftEmpty(output); return 15;
    }
    return 0; // that's it then
  }
  if(verbose>10) printf("frame times are not the same\n");

  // Check that there is no need for extrapolation in the start
  if(dftInterpolateCheckStart(input, output, status, verbose)<0) {
    dftEmpty(output); return 16;
  }
  // Check that there is no need for excess extrapolation in the end either
  if(dftInterpolateCheckEnd(input, output, status, verbose)<0) { 
    dftEmpty(output); return 17;
  }

  // Interpolate and integrate input data to tissue sample times,
  // taking into account tissue frame lengths
  {
    int ret=0;
    for(int ri=0; ri<output->voiNr && ret==0; ri++) {
      ret=interpolate4pet(input->x, input->voi[ri].y, input->frameNr,
                  output->x1, output->x2, output->voi[ri].y, output->voi[ri].y2,
                  output->voi[ri].y3, output->frameNr);
    }
    if(ret) {
      if(status!=NULL) sprintf(status, "cannot interpolate (%d)", ret);
      dftEmpty(output); return 18;
    }
  }

  return 0;
}
/*****************************************************************************/

/*****************************************************************************/
/** Integration of dynamic image from time1 to time2, storing integrals in iimg,
    which is allocated here.
    Frames do not have to be continuous in time. Time unit in integral is sec.
    Raw data (sinogram) must be divided by frame durations before calling this.
    You may want to test the time range before applying this routine, because
    this function accepts relatively large extrapolation. 
   @sa imgFrameIntegral, dftInterpolateForIMG, imgExistentTimes
   @return Returns errstatus, which is STATUS_OK (0) when call was successful,
    and >0 in case of an error.
*/
int imgTimeIntegral(
  /** IMG data; preferably dynamic, if static, then the specified time range
      must match with the frame time. */
  IMG *img,
  /** Time where to start integration (sec). */
  float t1,
  /** Time to stop integration (sec); must be higher than t1. */
  float t2,
  /** Pointer to initiated but empty AUC IMG data. */
  IMG *iimg,
  /** Calculate integral or average: 0=integral, 1=average. */
  int calc_mode,
  /** Pointer to a string (allocated for at least 128 chars) where error message
      or other execution status will be written; enter NULL, if not needed. */
  char *status,
  /** Verbose level; if zero, then nothing is printed to stderr or stdout. */
  int verbose
) {
  int ret;
  double accept_tdif=1.0; // in seconds

  if(verbose>0) {
    printf("%s(img, %g, %g, iimg, %d, status, %d)\n", __func__, t1, t2, calc_mode, verbose);
    fflush(stdout);
  }

  /* Check the arguments */
  if(status!=NULL) sprintf(status, "program error");
  if(img==NULL || iimg==NULL || t1<0.0 || t2<0.0) return STATUS_FAULT;
  if(img->status!=IMG_STATUS_OCCUPIED) return STATUS_FAULT;
  float fdur=t2-t1; if(fdur<=0.0) return STATUS_FAULT;
  if(img->dimt<1) return STATUS_FAULT;

  /* Check that time range matches with image frame */
  if(img->dimt==1) { // static image (one frame only)
    /* check that specified times match with image frame times, but accept 1 s differences */
    if(fabs(img->end[0]-t2)>accept_tdif || fabs(img->start[0]-t1)>accept_tdif) {
      if(status!=NULL) sprintf(status, 
        "for static image the integration time range must be exactly as long as the scan");
      return STATUS_FAULT;
    }
    /* Set the times, in case there was a small difference */
    img->end[0]=t2; img->start[0]=t1; img->mid[0]=0.5*(t1+t2);
  } else { // Dynamic image (more than one frame)
    /* First PET frame must start before 1/3 of the integration time */
    /* Last PET frame must end after 2/3 of integration time */
    if(img->start[0]>(0.66*t1+0.34*t2) || img->end[img->dimt-1]<(0.34*t1+0.66*t2)) {
      if(status!=NULL) sprintf(status, "integration time range oversteps data range");
      return STATUS_FAULT;
    }
  }
  if(verbose>10) printf("t1=%g t2=%g fdur=%g\n", t1, t2, fdur);

  /* Get the first and last frame index that resides inside integration time */
  int f1, f2;
  if(img->dimt==1) {
    f1=f2=0;
  } else {
    f1=f2=-1;
    for(int fi=0; fi<img->dimt; fi++) {
      if(f1<0) {if(img->start[fi]>=t1 && img->end[fi]<=t2) f1=f2=fi;}
      if(f1>=0) {if(t2>=img->end[fi]) f2=fi;}
    }
  }
  if(verbose>10) printf("f1=%d f2=%d\n", f1, f2);

  float aucroi, t[2], auc[2];
  if(f1>=0 && f2>=0) {

    /* Integrate over the frames that are included in time range as a whole */
    ret=imgFrameIntegral(img, f1, f2, iimg, verbose-1);
    if(ret) {
      if(status!=NULL) sprintf(status, "cannot integrate (%d)", ret);
      return STATUS_FAULT;
    }

    /* If necessary, add the partial integrals */
    if(img->start[f1]>t1) {
      t[0]=t1; t[1]=img->start[f1];
      if(verbose>20) printf("t[0]=%g t[1]=%g\n", t[0], t[1]); 
      for(int zi=0; zi<img->dimz; zi++)
        for(int yi=0; yi<img->dimy; yi++)
          for(int xi=0; xi<img->dimx; xi++) {
            ret=finterpolate(img->mid, img->m[zi][yi][xi], img->dimt, t, NULL, auc, NULL, 2);
            if(ret) aucroi=0.0; else aucroi=auc[1]-auc[0];
            iimg->m[zi][yi][xi][0]+=aucroi;
          }
    }
    if(t2>img->end[f2]) {
      t[0]=img->end[f2]; t[1]=t2;
      if(verbose>20) printf("t[0]=%g t[1]=%g\n", t[0], t[1]);
      for(int zi=0; zi<img->dimz; zi++)
        for(int yi=0; yi<img->dimy; yi++)
          for(int xi=0; xi<img->dimx; xi++) {
            ret=finterpolate(img->mid, img->m[zi][yi][xi], img->dimt, t, NULL, auc, NULL, 2);
            if(ret) aucroi=0.0; else aucroi=auc[1]-auc[0];
            iimg->m[zi][yi][xi][0]+=aucroi;
          }
    }

  } else { // no full frames inside integration range

    /* Create the AUC image */
    ret=imgAllocateWithHeader(iimg, img->dimz, img->dimy, img->dimx, 1, img);
    if(ret!=0) {
      if(status!=NULL) sprintf(status, "cannot setup integral image");
      return STATUS_FAULT;
    }

    t[0]=t1; t[1]=t2;
    for(int zi=0; zi<img->dimz; zi++)
      for(int yi=0; yi<img->dimy; yi++)
        for(int xi=0; xi<img->dimx; xi++) {
          ret=finterpolate(img->mid, img->m[zi][yi][xi], img->dimt, t, NULL, auc, NULL, 2);
          if(ret) aucroi=0.0; else aucroi=auc[1]-auc[0];
          iimg->m[zi][yi][xi][0]+=aucroi;
        }

  }

  /* Set output image time frame */
  iimg->end[0]=t2; iimg->start[0]=t1; iimg->mid[0]=0.5*(t1+t2);

  /* If required, then calculate average by dividing integral with time. */
  /* Set units accordingly */
  if(calc_mode!=0) { // average
    ret=imgArithmConst(iimg, fdur, ':', 1.0E+10, verbose-1);
    if(ret!=0) {
      if(status!=NULL) sprintf(status, "cannot divide integral image");
      return STATUS_FAULT;
    }
    iimg->unit=img->unit;
    if(status!=NULL) sprintf(status, "average image [%g,%g] calculated",t1,t2);
  } else { // integral
    if(img->unit==CUNIT_KBQ_PER_ML) iimg->unit=CUNIT_SEC_KBQ_PER_ML;
    else iimg->unit=CUNIT_UNKNOWN;
    if(status!=NULL) sprintf(status, "integral image [%g,%g] calculated",t1,t2);
  }

  imgSetStatus(iimg, STATUS_OK);
  return(0);
}
/*****************************************************************************/

/*****************************************************************************/
/** Allocate memory for DFT based on information in IMG.
   @sa sifAllocateWithIMG
   @return Returns 0 if successful, otherwise <>0.
 */
int dftAllocateWithIMG(
  /** Pointer to initiated DFT struct which will be allocated here. */
  DFT *dft,
  /** Nr of TACs to be allocated in DFT; if zero, then TACs for all IMG
      pixels is allocated, but it may lead to out of memory error. */
  int tacNr,
  /** Pointer to IMG struct from where necessary information is read. */
  IMG *img
) {
  int ri, fi, ret;

  // Check the input data
  if(dft==NULL || img==NULL) return 1;
  if(img->status!=IMG_STATUS_OCCUPIED) return 2;
  if(img->dimt<1) return 3;

  // Get tacNr from image dimensions if necessary
  if(tacNr<1) {
    tacNr=img->dimz*img->dimx*img->dimy;
    if(tacNr<1) return 4;
  }

  // Allocate memory
  ret=dftSetmem(dft, img->dimt, tacNr);
  if(ret) return(100+ret);
  dft->voiNr=tacNr;
  dft->frameNr=img->dimt;

  // Set header contents
  dft->timetype=DFT_TIME_STARTEND;
  for(fi=0; fi<dft->frameNr; fi++) {
    dft->x[fi]=img->mid[fi];
    dft->x1[fi]=img->start[fi];
    dft->x2[fi]=img->end[fi];
  }
  dft->isweight=0;
  strncpy(dft->unit, imgUnit(img->unit), MAX_UNITS_LEN);
  dft->timeunit=TUNIT_SEC;
  dft->_type=DFT_FORMAT_STANDARD;
  for(ri=0; ri<dft->voiNr; ri++) {
    snprintf(dft->voi[ri].voiname, 6, "%06d", ri+1);
    strcpy(dft->voi[ri].name, dft->voi[ri].voiname);
  }
  strlcpy(dft->studynr, img->studyNr, MAX_STUDYNR_LEN);

  return 0;
}
/*****************************************************************************/

/*****************************************************************************/
/** Convert SIF data to DFT data.
   @sa sifAllocateWithIMG, dftAllocateWithIMG
   @return Returns 0 if successful, otherwise >0.
 */
int sif2dft(
  /** Pointer to SIF struct containing data to be converted. */
  SIF *sif,
  /** Pointer to initiated but empty DFT struct. */
  DFT *dft
) {
  int fi, ri, tacNr, ret;

  /* Check input */
  if(dft==NULL || sif==NULL) return 1;
  if(sif->frameNr<1) return 2;
  tacNr=sif->colNr-2; if(tacNr<=0) tacNr=1;

  /* Allocate memory */
  ret=dftSetmem(dft, sif->frameNr, tacNr);
  if(ret) return(100+ret);
  dft->voiNr=tacNr;
  dft->frameNr=sif->frameNr;

  // Set header contents
  dft->_type=DFT_FORMAT_STANDARD;
  dft->timetype=DFT_TIME_STARTEND;
  for(fi=0; fi<dft->frameNr; fi++) {
    dft->x[fi]=0.5*(sif->x1[fi]+sif->x2[fi]);
    dft->x1[fi]=sif->x1[fi];
    dft->x2[fi]=sif->x2[fi];
  }
  dft->timeunit=TUNIT_SEC;
  dft->isweight=0;
  strcpy(dft->unit, imgUnit(CUNIT_COUNTS));
  for(ri=0; ri<dft->voiNr; ri++) {
    /* Set ROI name */
    if(ri==0) sprintf(dft->voi[ri].voiname, "Prompt");
    else if(ri==1) sprintf(dft->voi[ri].voiname, "Random");
    else if(ri==2) sprintf(dft->voi[ri].voiname, "True");
    else if(ri==3) sprintf(dft->voi[ri].voiname, "Weight");
    else sprintf(dft->voi[ri].voiname, "%06d", ri+1);
    strcpy(dft->voi[ri].name, dft->voi[ri].voiname);
    /* Copy TAC values */
    if(ri==0 && ri<sif->colNr-2) {
      for(fi=0; fi<dft->frameNr; fi++) dft->voi[ri].y[fi]=sif->prompts[fi];
    } else if(ri==1 && ri<sif->colNr-2) {
      for(fi=0; fi<dft->frameNr; fi++) dft->voi[ri].y[fi]=sif->randoms[fi];
    } else if(ri==2 && ri<sif->colNr-2) {
      for(fi=0; fi<dft->frameNr; fi++) dft->voi[ri].y[fi]=sif->trues[fi];
    } else if(ri==3 && ri<sif->colNr-2) {
      for(fi=0; fi<dft->frameNr; fi++) dft->voi[ri].y[fi]=sif->weights[fi];
    } else {
      for(fi=0; fi<dft->frameNr; fi++) dft->voi[ri].y[fi]=0.0;
    }
  }
  strcpy(dft->studynr, sif->studynr);
  strcpy(dft->isotope, hlCorrectIsotopeCode(sif->isotope_name));

  return 0;
}
/*****************************************************************************/

/*****************************************************************************/
/** Allocate memory for SIF based on information in IMG.

    Frame times are set, 'counts' are optionally set to mean of
    voxel values times frame lengths, scaled to 0-1E7.
    Weights are not calculated.

   @sa sif2dft, dftAllocateWithIMG
   @return Returns 0 if successful, otherwise <>0.
 */
int sifAllocateWithIMG(
  /** Pointer to initiated SIF struct which will be allocated here. */
  SIF *sif,
  /** Pointer to IMG struct from where necessary information is read. */
  IMG *img,
  /** Set (1) or do not set (0) counts, based on sum of all image voxel values and frame lengths. */
  int doCounts,
  /** Verbose level; if zero, then only warnings are printed into stderr. */
  int verbose
) {
  int fi, ret;
  float *cs;
  double f, mf;

  if(verbose>0) {printf("%s(*sif, *img, %d, ...)\n", __func__, doCounts); fflush(stdout);}
  /* Check the input data */
  if(sif==NULL || img==NULL) return 1;
  if(img->status!=IMG_STATUS_OCCUPIED) return 2;
  if(img->dimt<1) return 3;
  if(doCounts<0 || doCounts>1) return 4;
  
  /* Delete any previous SIF contents */
  sifEmpty(sif);
  
  /* Allocate memory for SIF data */
  ret=sifSetmem(sif, img->dimt); if(ret!=0) return(10+ret);

  /* Set SIF information */
  sif->version=1;
  sif->colNr=4;
  strcpy(sif->isotope_name, imgIsotope(img));
  strcpy(sif->studynr, img->studyNr);
  sif->scantime=img->scanStart;
  for(fi=0; fi<img->dimt; fi++) {sif->x1[fi]=img->start[fi]; sif->x2[fi]=img->end[fi];}
  
  if(doCounts==0) return 0; // that's it then
  
  /* Calculate average curve of all pixels */
  if(verbose>1) printf("calculate image average curve.\n");
  cs=(float*)malloc(img->dimt*sizeof(float));
  if(cs==NULL) {sifEmpty(sif); return(20);}
  ret=imgAverageTAC(img, cs); if(ret!=0) {sifEmpty(sif); return(20+ret);}
  /* Multiply average curve with frame durations, and get the max */
  for(fi=0, mf=0.0; fi<sif->frameNr; fi++) {
    f=sif->x2[fi]-sif->x1[fi]; if(f<=0.1) f=0.1;
    cs[fi]*=f; if(cs[fi]>mf) mf=cs[fi];
  }
  /* Put scaled counts into SIF */
  if(mf>0.0) f=1.0E+007/mf; else f=1.0;
  for(fi=0, mf=0.0; fi<sif->frameNr; fi++) {
    sif->prompts[fi]=sif->trues[fi]=cs[fi]*f;
    sif->randoms[fi]=0.0;
  }
  free(cs);
  if(verbose>2) sifPrint(sif);
  return 0;
}
/*****************************************************************************/

/*****************************************************************************/
