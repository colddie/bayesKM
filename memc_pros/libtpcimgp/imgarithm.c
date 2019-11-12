/// @file imgarithm.c
/// @author Vesa Oikonen
/// @brief Misc arithmetical routines for processing 4D IMG data.
///
/*****************************************************************************/
#include "libtpcimgp.h"
/*****************************************************************************/

/*****************************************************************************/
/** Simple arithmetics between matching IMG planes and frames.
    Specify the operation as one of characters +, -, /, :, *, ., x.
    Results that are higher than ulimit are set to ulimit.
    @sa imgArithmConst, imgArithmFrame, imgRawCountsPerTime, imgConvertUnit
    @return Returns 0, if ok.
 */
int imgArithm(
  /** The first IMG data. */
  IMG *img1,
  /** The second IMG data. */
  IMG *img2,
  /** Operation, one of the characters +, -, /, :, *, ., x. */
  char operation,
  /** Results that are higher than ulimit are set to ulimit; set to <=0 if not needed. */
  float ulimit,
  /** Verbose level; if <=0, then nothing is printed into stderr. */
  int verbose
) {
  int pi, yi, xi, fi, ret;

  if(verbose>0) printf("imgArithm(img1, img2, '%c', %g, %d)\n", operation, ulimit, verbose);

  /* Check the arguments */
  if(img1->status!=IMG_STATUS_OCCUPIED || img2->status!=IMG_STATUS_OCCUPIED) {
    if(verbose>0) fprintf(stderr, "Invalid image status.\n");
    return(1);
  }
  ret=0;
  if(img1->dimx!=img2->dimx || img1->dimy!=img2->dimy) ret=1;
  /* Check the plane numbers */
  if(img1->dimz!=img2->dimz) ret=2;
  if(ret==0) for(pi=0; pi<img1->dimz; pi++)
    if(img1->planeNumber[pi]!=img2->planeNumber[pi]) ret=2;
  /* Check the frame numbers */
  if(img1->dimt!=img2->dimt) ret=3;
  if(ret>0) {
    if(verbose>0) fprintf(stderr, "Image dimensions do not match.\n");
    return(ret);
  }

  /* Operate */
  switch(operation) {
    case '+':
      for(pi=0; pi<img1->dimz; pi++) for(yi=0; yi<img1->dimy; yi++)
        for(xi=0; xi<img1->dimx; xi++) for(fi=0; fi<img1->dimt; fi++)
          img1->m[pi][yi][xi][fi]+=img2->m[pi][yi][xi][fi];
      break;
    case '-':
      for(pi=0; pi<img1->dimz; pi++) for(yi=0; yi<img1->dimy; yi++)
        for(xi=0; xi<img1->dimx; xi++) for(fi=0; fi<img1->dimt; fi++)
          img1->m[pi][yi][xi][fi]-=img2->m[pi][yi][xi][fi];
      break;
    case '/':
    case ':':
      for(pi=0; pi<img1->dimz; pi++) for(yi=0; yi<img1->dimy; yi++)
        for(xi=0; xi<img1->dimx; xi++) for(fi=0; fi<img1->dimt; fi++) {
          if(fabs(img2->m[pi][yi][xi][fi])>1.0E-5)
            img1->m[pi][yi][xi][fi]/=img2->m[pi][yi][xi][fi];
          else img1->m[pi][yi][xi][fi]=0.0;
        }
      break;
    case '*':
    case 'x':
    case '.':
      for(pi=0; pi<img1->dimz; pi++) for(yi=0; yi<img1->dimy; yi++)
        for(xi=0; xi<img1->dimx; xi++) for(fi=0; fi<img1->dimt; fi++)
          img1->m[pi][yi][xi][fi]*=img2->m[pi][yi][xi][fi];
      break;
    default:
      if(verbose>0) fprintf(stderr, "Invalid operation.\n");
      return(10);
  }

  /* Check for limits */
  if(ulimit>0.0) {
    for(pi=0; pi<img1->dimz; pi++) for(yi=0; yi<img1->dimy; yi++)
      for(xi=0; xi<img1->dimx; xi++) for(fi=0; fi<img1->dimt; fi++)
        if(img1->m[pi][yi][xi][fi]>ulimit) img1->m[pi][yi][xi][fi]=ulimit;
  }

  return(0);
}
/*****************************************************************************/

/*****************************************************************************/
/** Simple arithmetics between IMG and specified constant.

    Specify the operation as one of characters +, -, /, :, *, ., x.
    Results that are higher than ulimit are set to ulimit.
    @sa imgArithm, imgScale, imgInv, imgArithmFrame, imgRawCountsPerTime, imgConvertUnit
    @return Returns 0, if ok.
 */
int imgArithmConst(
  /** IMG data which is modified. */
  IMG *img,
  /** Constant value which is used to modify the data. */
  float operand,
  /** Operation, one of the characters +, -, /, :, *, ., x. */
  char operation,
  /** Results that are higher than ulimit are set to ulimit;
      set to <=0 if not needed. */
  float ulimit,
  /** Verbose level; if <=0, then nothing is printed into stderr. */
  int verbose
) {
  int pi, yi, xi, fi;

  if(verbose>0) printf("imgArithConst(img, %g, '%c', %g, %d)\n", operand, operation, ulimit, verbose);

  /* Check the arguments */
  if(img->status!=IMG_STATUS_OCCUPIED) {
    if(verbose>0) fprintf(stderr, "Invalid image status.\n");
    return(1);
  }
  switch(operation) {
    case '+':
      for(pi=0; pi<img->dimz; pi++) for(yi=0; yi<img->dimy; yi++)
        for(xi=0; xi<img->dimx; xi++) for(fi=0; fi<img->dimt; fi++)
          img->m[pi][yi][xi][fi]+=operand;
      break;
    case '-':
      for(pi=0; pi<img->dimz; pi++) for(yi=0; yi<img->dimy; yi++)
        for(xi=0; xi<img->dimx; xi++) for(fi=0; fi<img->dimt; fi++)
          img->m[pi][yi][xi][fi]-=operand;
      break;
    case '/':
    case ':':
      if(fabs(operand)<1.0e-100) return(2);
      for(pi=0; pi<img->dimz; pi++) for(yi=0; yi<img->dimy; yi++)
        for(xi=0; xi<img->dimx; xi++) for(fi=0; fi<img->dimt; fi++)
          img->m[pi][yi][xi][fi]/=operand;
      break;
    case '*':
    case 'x':
    case '.':
      for(pi=0; pi<img->dimz; pi++) for(yi=0; yi<img->dimy; yi++)
        for(xi=0; xi<img->dimx; xi++) for(fi=0; fi<img->dimt; fi++)
          img->m[pi][yi][xi][fi]*=operand;
      break;
    default:
      if(verbose>0) fprintf(stderr, "Invalid operation.\n");
      return(10);
  }

  /* Check for limits */
  if(ulimit>0.0) {
    for(pi=0; pi<img->dimz; pi++) for(yi=0; yi<img->dimy; yi++)
      for(xi=0; xi<img->dimx; xi++) for(fi=0; fi<img->dimt; fi++)
        if(img->m[pi][yi][xi][fi]>ulimit) img->m[pi][yi][xi][fi]=ulimit;
  }

  return(0);
}
/*****************************************************************************/

/*****************************************************************************/
/** Simple arithmetics between matching IMG planes and the first frame of img2.

    Specify the operation as one of characters +, -, /, :, *, ., x.
    Results that are higher than ulimit are set to ulimit.
    @sa imgArithmConst, imgArithm, imgFrameIntegral
    @return Returns 0, if ok.
 */
int imgArithmFrame(
  /** IMG data which is modified. */
  IMG *img1,
  /** Operand IMG data, only the first frame is used. */
  IMG *img2,
  /** Operation, one of the characters +, -, /, :, *, ., x. */
  char operation,
  /** Results that are higher than ulimit are set to ulimit; set to <=0 if not needed. */
  float ulimit,
  /** Verbose level; if <=0, then nothing is printed into stderr. */
  int verbose
) {
  int pi, yi, xi, fi, ret;

  if(verbose>0) printf("imgArithFrame(img1, img2, '%c', %g, %d)\n", operation, ulimit, verbose);

  /* Check the arguments */
  if(img1->status!=IMG_STATUS_OCCUPIED || img2->status!=IMG_STATUS_OCCUPIED) {
    if(verbose>0) fprintf(stderr, "Invalid image status.\n");
    return(1);
  }
  ret=0;
  if(img1->dimx!=img2->dimx || img1->dimy!=img2->dimy) ret=1;
  /* Check the plane numbers */
  if(img1->dimz!=img2->dimz) ret=2;
  if(ret==0) for(pi=0; pi<img1->dimz; pi++)
    if(img1->planeNumber[pi]!=img2->planeNumber[pi]) ret=2;
  if(ret>0) {
    if(verbose>0) fprintf(stderr, "Image dimensions do not match.\n");
    return(ret);
  }

  /* Operate */
  switch(operation) {
    case '+':
      for(pi=0; pi<img1->dimz; pi++) for(yi=0; yi<img1->dimy; yi++)
        for(xi=0; xi<img1->dimx; xi++) for(fi=0; fi<img1->dimt; fi++)
          img1->m[pi][yi][xi][fi]+=img2->m[pi][yi][xi][0];
      break;
    case '-':
      for(pi=0; pi<img1->dimz; pi++) for(yi=0; yi<img1->dimy; yi++)
        for(xi=0; xi<img1->dimx; xi++) for(fi=0; fi<img1->dimt; fi++)
          img1->m[pi][yi][xi][fi]-=img2->m[pi][yi][xi][0];
      break;
    case '/':
    case ':':
      for(pi=0; pi<img1->dimz; pi++) for(yi=0; yi<img1->dimy; yi++)
        for(xi=0; xi<img1->dimx; xi++) {
          if(fabs(img2->m[pi][yi][xi][0])>1.0E-8)
            for(fi=0; fi<img1->dimt; fi++)
              img1->m[pi][yi][xi][fi]/=img2->m[pi][yi][xi][0];
          else for(fi=0; fi<img1->dimt; fi++) img1->m[pi][yi][xi][fi]=0.0;
        }
      break;
    case '*':
    case 'x':
    case '.':
      for(pi=0; pi<img1->dimz; pi++) for(yi=0; yi<img1->dimy; yi++)
        for(xi=0; xi<img1->dimx; xi++) for(fi=0; fi<img1->dimt; fi++)
          img1->m[pi][yi][xi][fi]*=img2->m[pi][yi][xi][0];
      break;
    default:
      if(verbose>0) fprintf(stderr, "Invalid operation.\n");
      return(10);
  }

  /* Check for limits */
  if(ulimit>0.0) {
    for(pi=0; pi<img1->dimz; pi++) for(yi=0; yi<img1->dimy; yi++)
      for(xi=0; xi<img1->dimx; xi++) for(fi=0; fi<img1->dimt; fi++)
        if(img1->m[pi][yi][xi][fi]>ulimit) img1->m[pi][yi][xi][fi]=ulimit;
  }

  return(0);
}
/*****************************************************************************/

/*****************************************************************************/
/** Replace IMG data values by their natural logarithms.
    @sa imgLog10, imgAbs, imgArithm, imgInv
    @return Returns 0, if ok.
 */
int imgLn(
  /** Pointer to IMG data. */
  IMG *img
) {
  if(img->status<IMG_STATUS_OCCUPIED) return(1);
  if(img->dimt<1 || img->dimz<1 || img->dimy<1 || img->dimx<1) return(2);
  for(int pi=0; pi<img->dimz; pi++)
    for(int yi=0; yi<img->dimy; yi++)
      for(int xi=0; xi<img->dimx; xi++)
        for(int fi=0; fi<img->dimt; fi++) {
          if(img->m[pi][yi][xi][fi]<=0.0) img->m[pi][yi][xi][fi]=0.0;
          else img->m[pi][yi][xi][fi]=log(img->m[pi][yi][xi][fi]);
        }
  return(0);
}
/*****************************************************************************/

/*****************************************************************************/
/** Replace IMG data values by their log10 values.
    @sa imgLn, imgAbs, imgArithm
    @return Returns 0, if ok.
 */
int imgLog10(
  /** Pointer to IMG data. */
  IMG *img
) {
  if(img->status<IMG_STATUS_OCCUPIED) return(1);
  if(img->dimt<1 || img->dimz<1 || img->dimy<1 || img->dimx<1) return(2);
  for(int pi=0; pi<img->dimz; pi++)
    for(int yi=0; yi<img->dimy; yi++)
      for(int xi=0; xi<img->dimx; xi++)
        for(int fi=0; fi<img->dimt; fi++) {
          if(img->m[pi][yi][xi][fi]<=0.0) img->m[pi][yi][xi][fi]=0.0;
          else img->m[pi][yi][xi][fi]=log10(img->m[pi][yi][xi][fi]);
        }
  return(0);
}
/*****************************************************************************/

/*****************************************************************************/
/** Replace IMG data values by their absolute values.
    @sa imgLog10, imgLn, imgArithm
    @return Returns 0, if ok.
 */
int imgAbs(
  /** Pointer to IMG data. */
  IMG *img
) {
  if(img->status<IMG_STATUS_OCCUPIED) return(1);
  if(img->dimt<1 || img->dimz<1 || img->dimy<1 || img->dimx<1) return(2);
  for(int pi=0; pi<img->dimz; pi++)
    for(int yi=0; yi<img->dimy; yi++)
      for(int xi=0; xi<img->dimx; xi++)
        for(int fi=0; fi<img->dimt; fi++) 
          img->m[pi][yi][xi][fi]=fabsf(img->m[pi][yi][xi][fi]);
  return(0);
}
/*****************************************************************************/

/*****************************************************************************/
/** Replace IMG data values by their inverse value (1/activity).

    Pixels with <=0 are set to zero.
    @sa imgLn, imgLog10, imgAbs, imgArithm
    @return Returns 0, if ok.
 */
int imgInv(
  /** Pointer to IMG data. */
  IMG *img
) {
  if(img->status<IMG_STATUS_OCCUPIED) return(1);
  if(img->dimt<1 || img->dimz<1 || img->dimy<1 || img->dimx<1) return(2);
  for(int pi=0; pi<img->dimz; pi++)
    for(int yi=0; yi<img->dimy; yi++)
      for(int xi=0; xi<img->dimx; xi++)
        for(int fi=0; fi<img->dimt; fi++) {
          if(img->m[pi][yi][xi][fi]<=0.0) img->m[pi][yi][xi][fi]=0.0;
          else {
            img->m[pi][yi][xi][fi]=1.0/img->m[pi][yi][xi][fi];
            if(!isfinite(img->m[pi][yi][xi][fi])) img->m[pi][yi][xi][fi]=0.0;
          }
        }
  return(0);
}
/*****************************************************************************/

/*****************************************************************************/
/** Integration from first frame (0..last) to last frame (first..dimt) to iimg,
    which is allocated here.

    Frames do not have to be continuous in time. Time unit in integral is sec.
    Raw data (sinogram) must be divided by frame durations before calling this.
    If dynamic image data does not contain frame times (e.g. Analyze image)
    then just the sum is calculated.

    @sa imgArithm, imgRawCountsPerTime, imgTimeIntegral
    @return Returns 0 if ok, otherwise >0.
*/
int imgFrameIntegral(
  /** (Dynamic) IMG data. */
  IMG *img,
  /** First frame to include in AUC; 0..dimt-1. */
  int first,
  /** Last frame to include in AUC; 0..dimt-1. */
  int last,
  /** Pointer to initiated and empty AUC IMG data. */
  IMG *iimg,
  /** Verbose level; if <=0, then nothing is printed into stderr. */
  int verbose
) {
  int zi, yi, xi, fi, ret, times_exist;
  float fstart, fend, dur, x, y, k;

  if(verbose>0) printf("imgFrameIntegral(img, %d, %d, iimg, %d)\n", first, last, verbose);

  /* Check the arguments */
  if(img==NULL || iimg==NULL || first<0 || first>last) return(1);
  if(img->status!=IMG_STATUS_OCCUPIED) return(2);
  times_exist=imgExistentTimes(img);
  fstart=img->start[first]; fend=img->end[last];
  if(verbose>1) printf("  time_range := %g - %g\n", fstart, fend);

  /* Allocate memory for the integral */
  imgEmpty(iimg);
  ret=imgAllocateWithHeader(iimg, img->dimz, img->dimy, img->dimx, 1, img);
  if(ret) {
    if(verbose>0) fprintf(stderr, "Error %d in allocating memory for sum image.\n", ret);
    imgEmpty(iimg); return(3);
  }

  /* Set header information */
  iimg->start[0]=fstart; iimg->end[0]=fend;
  iimg->mid[0]=0.5*(iimg->start[0]+iimg->end[0]);
  if(img->type==IMG_TYPE_RAW) iimg->decayCorrection=IMG_DC_NONCORRECTED;
  else iimg->decayCorrection=IMG_DC_CORRECTED;
  iimg->decayCorrFactor[0]=0.0;
  /* Change also the unit, if possible */
  if(img->type==IMG_TYPE_IMAGE) {
    if(img->unit==CUNIT_KBQ_PER_ML) iimg->unit=CUNIT_SEC_KBQ_PER_ML;
  }

  /* Integrate the first frame */
  fi=first;
  if(times_exist) {
    dur=img->end[fi]-img->start[fi]; if(dur<0.0) {imgEmpty(iimg); return(4);}
  } else dur=1.0;
  for(zi=0; zi<img->dimz; zi++)
    for(yi=0; yi<img->dimy; yi++)
      for(xi=0; xi<img->dimx; xi++)
        iimg->m[zi][yi][xi][0]=dur*img->m[zi][yi][xi][fi];
  /* Add integrals of the following frames */
  for(fi=first+1; fi<=last; fi++) {
    if(times_exist) {
      dur=img->end[fi]-img->start[fi]; if(dur<0.0) {imgEmpty(iimg); return(4);}
    } else dur=1.0;
    /* Add the integral of this frame to the previous integral */
    for(zi=0; zi<img->dimz; zi++)
      for(yi=0; yi<img->dimy; yi++)
        for(xi=0; xi<img->dimx; xi++)
          iimg->m[zi][yi][xi][0]+=dur*img->m[zi][yi][xi][fi];
    if(times_exist) {
      /* Check whether frames are contiguous */
      dur=img->start[fi]-img->end[fi-1]; if(dur<=1.0E-10) continue;
      /* When not, calculate the integral between frames */
      x=0.5*(img->start[fi]+img->end[fi-1]);
      for(zi=0; zi<img->dimz; zi++)
        for(yi=0; yi<img->dimy; yi++)
          for(xi=0; xi<img->dimx; xi++) {
            k=(img->m[zi][yi][xi][fi]-img->m[zi][yi][xi][fi-1])
              /(img->mid[fi]-img->mid[fi-1]);
            y=img->m[zi][yi][xi][fi-1]+k*(x-img->mid[fi-1]);
            iimg->m[zi][yi][xi][0]+=dur*y;
          }
    }    
  } /* next frame */

  /* Set frame times */
  iimg->start[0]=img->start[first];
  iimg->end[0]=img->end[last];
  iimg->mid[0]=0.5*(iimg->start[0]+iimg->end[0]);

  return(0);
}
/*****************************************************************************/

/*****************************************************************************/
/** Divide or multiply raw data (sinogram) counts by frame duration.

    If IMG is not raw data, division is quietly not done.
    Error is returned if sinogram does not contain frame times, except if
    sinogram contains only one time frame, also then division is not done.
    @sa imgArithmConst, imgArithm, imgFrameIntegral
    @return Returns 0 if ok.
*/
int imgRawCountsPerTime(
  /** IMG containing raw data. */
  IMG *img,
  /** 1=division, 0=multiply. */
  int operation
) {
  int fi, zi, xi, yi, notime=0;
  float f;

  /* Check the input data */
  if(operation!=0 && operation!=1) return(1);
  if(img==NULL) return(2);
  if(img->status!=IMG_STATUS_OCCUPIED) return(3);
  if(img->type!=IMG_TYPE_RAW) return(0);
  for(fi=0; fi<img->dimt; fi++) {
    f=img->end[fi]-img->start[fi]; if(f<=1.0E-12) notime++;
  }
  if(notime>0) {
    if(img->dimt>1) return(4);
    else return(0); /* just one frame; might be parametric or transmission sinogram,
                       therefore no error. */
  }
  for(fi=0; fi<img->dimt; fi++) {
    f=img->end[fi]-img->start[fi]; if(operation==1) f=1.0/f;
    for(zi=0; zi<img->dimz; zi++)
      for(yi=0; yi<img->dimy; yi++)
        for(xi=0; xi<img->dimx; xi++)
          img->m[zi][yi][xi][fi]*=f;
  }
  return(0);
}
/*****************************************************************************/

/*****************************************************************************/
/** Converts the unit of pixel values in IMG based to specified unit string.
    @sa imgArithmConst, imgArithm, imgFrameIntegral
    @return Returns 0 if successful.
 */
int imgConvertUnit(
  /** Pointer to IMG struct. */
  IMG *img,
  /** String containing the new unit. */
  char *unit
) {
  int ret, new_unit=0;
  float conversion_factor=1.0;

  if(img==NULL || unit==NULL) return(1);
  if(img->unit==CUNIT_UNKNOWN) return(1);
  /* Try to identify requested unit */
  new_unit=imgUnitId(unit); if(new_unit<0) return(1-new_unit);
  /* Check if unit needs no conversion */
  if(img->unit==new_unit) return(0);
  /* Get conversion factor */
  if(img->unit==CUNIT_KBQ_PER_ML && new_unit==CUNIT_BQ_PER_ML)
    conversion_factor=1000.0;
  else if(img->unit==CUNIT_BQ_PER_ML && new_unit==CUNIT_KBQ_PER_ML)
    conversion_factor=0.001;
  else if(img->unit==CUNIT_KBQ_PER_ML && new_unit==CUNIT_NCI_PER_ML)
    conversion_factor=27.027;
  else if(img->unit==CUNIT_NCI_PER_ML && new_unit==CUNIT_KBQ_PER_ML)
    conversion_factor=0.037;
  else if(img->unit==CUNIT_NCI_PER_ML && new_unit==CUNIT_BQ_PER_ML)
    conversion_factor=37.0;
  else if(img->unit==CUNIT_KBQ_PER_ML && new_unit==CUNIT_MBQ_PER_ML)
    conversion_factor=0.001;
  else if(img->unit==CUNIT_MBQ_PER_ML && new_unit==CUNIT_KBQ_PER_ML)
    conversion_factor=1000.0;
  else if(img->unit==CUNIT_PER_SEC && new_unit==CUNIT_PER_MIN)
    conversion_factor=60.0;
  else if(img->unit==CUNIT_PER_MIN && new_unit==CUNIT_PER_SEC)
    conversion_factor=1.0/60.0;
  else if(img->unit==CUNIT_ML_PER_ML && new_unit==CUNIT_ML_PER_DL)
    conversion_factor=0.01;
  else if(img->unit==CUNIT_ML_PER_DL && new_unit==CUNIT_ML_PER_ML)
    conversion_factor=100.0;
  else
    return(10);

  /* Convert pixel values */
  ret=imgArithmConst(img, conversion_factor, '*', FLT_MAX, 0);
  if(ret) return(10+ret);

  /* Set the unit id */
  img->unit=new_unit;

  return(0);
}
/*****************************************************************************/

/*****************************************************************************/
