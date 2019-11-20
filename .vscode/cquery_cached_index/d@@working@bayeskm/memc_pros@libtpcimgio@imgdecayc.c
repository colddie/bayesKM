/// @file imgdecayc.c
/// @author Vesa Oikonen
/// @brief Physical decay and isotopes in IMG.
///
/*****************************************************************************/
#include "libtpcimgio.h"
/*****************************************************************************/

/*****************************************************************************/
/*!
 *  Corrects (mode=1) or removes correction (mode=0) for physical decay.
 *  Removal is based on existing decay correction factors, when possible.
 *
 * @param image pointer to IMG data
 * @param mode 0=Remove decay correction; 1=Correct for decay
 * @return 0 if ok, 1 image status is not 'occupied',
 * 2 decay already corrected/not corrected, 3 image frame times missing
 */
int imgDecayCorrection(IMG *image, int mode) {
  int pi, fi, i, j;
  float lambda;
  float cf, dur;

  /* Check for arguments */
  if(image->status!=IMG_STATUS_OCCUPIED) return(1);
  if(image->isotopeHalflife<=0.0) return(1);
  /* Existing/nonexisting decay correction is an error */
  if(mode==1 && image->decayCorrection!=IMG_DC_NONCORRECTED) return(2);
  if(mode==0 && image->decayCorrection!=IMG_DC_CORRECTED) return(2);

  /* All time frames */
  for(fi=0; fi<image->dimt; fi++) {
    dur=image->end[fi]-image->start[fi];
    if(image->end[fi]>0.0) {
      if(mode==0 && image->decayCorrFactor[fi]>1.000001) {
        /* if decay correction is to be removed, and factor is known,
           then use it */
        cf=1.0/image->decayCorrFactor[fi];
      } else {
        lambda=hl2lambda(image->isotopeHalflife); if(lambda<0.0) return(1);
        /* remove decay correction by giving negative lambda */
        if(mode==0) lambda=-lambda;
        if(fi==image->dimt-1 && image->end[fi]<=0.0) return(3);
        cf=hlLambda2factor_float(lambda, image->start[fi], dur);
      }
      if(IMG_TEST) printf("applied_dc_factor[%d] := %g\n", fi+1, cf);
      /* Set decay correction factor inside IMG for future */
      if(mode==0) {
        image->decayCorrFactor[fi]=1.0;
      } else {
        image->decayCorrFactor[fi]=cf;
      }
      /* All planes, all matrices */
      for(pi=0; pi<image->dimz; pi++)
        for(i=0; i<image->dimy; i++)
          for(j=0; j<image->dimx; j++)
            image->m[pi][i][j][fi]*=cf;
      if(mode==0) image->decayCorrection=IMG_DC_NONCORRECTED;
      else image->decayCorrection=IMG_DC_CORRECTED;
      /* in some cases left unchanged! */
    }
  } /* next frame */
  return(0);
}
/*****************************************************************************/

/*****************************************************************************/
/*!
 * Returns pointer to string describing the isotope in image data
 *
 * @param img image stucture
 * @return pointer to string
 */
char *imgIsotope(IMG *img) {
  return(hlIsotopeCode(hlIsotopeFromHalflife(img->isotopeHalflife/60.0)));
}
/*****************************************************************************/

/*****************************************************************************/
/*!
 * Sets (mode=1) or removes (mode=0) decay correction factors in IMG.
 * IMG pixel data is not changed.
 *
 * @param image pointer to IMG data
 * @param mode factors are calculated for decay correction (1) or
 * for removing decay correction (0)
 * @return 0 if ok, 1 image status is not 'occupied', 2 invalid exponent value,
 * 3 image frame times are missing
 */
int imgSetDecayCorrFactors(IMG *image, int mode) {
  int fi;
  float lambda, cf, dur;

  /* Check for arguments */
  if(image->status!=IMG_STATUS_OCCUPIED) return(1);
  if(image->isotopeHalflife<=0.0) return(1);

  /* Check that image contains frame times */
  if(mode!=0 && image->end[image->dimt-1]<=0.0) return(3);

  /* All time frames */
  for(fi=0; fi<image->dimt; fi++) {
    if(mode==0) {
      image->decayCorrFactor[fi]=1.0;
    } else {
      dur=image->end[fi]-image->start[fi];
      if(image->end[fi]>0.0) {
        lambda=hl2lambda(image->isotopeHalflife); if(lambda<0.0) return(2);
        cf=hlLambda2factor_float(lambda, image->start[fi], dur);
        image->decayCorrFactor[fi]=cf;
      }
    }
  } /* next frame */
  if(mode==0) image->decayCorrection=IMG_DC_NONCORRECTED;
  else image->decayCorrection=IMG_DC_CORRECTED;
  return(0);
}
/*****************************************************************************/

/*****************************************************************************/
/** Corrects image data for branching fraction (mode=1) or removes correction
 *  (mode=0). Removal is primarily based on branching factor stored in 
 *  IMG struct, secondarily on isotope; after removal, branching factor is
 *  set to 1, and pixel values and calibration factor are multiplied with it.
 *  Correction is based on branching fractions in branch.h; pixel values and
 *  calibration factor are divided by it, and its value is stored in IMG struct.
 *  
 *  Note that this function can not know if branching fraction correction is
 *  included in the data (as it usually is) or not.
 *
\return Returns 0 if ok.
 */
int imgBranchingCorrection(
  /** Pointer to IMG data */
  IMG *image,
  /** Branching fraction correction (1) or removal of correction (0) */
  int mode,
  /** Verbose level; if zero, then nothing is printed into stdout or stderr */
  int verbose,
  /** Pointer to allocated string where error message will be written;
   *  NULL, if not needed. */
  char *status
) {
  int pi, fi, i, j, isotope;
  float bf;
  float cf;

  if(verbose>0) printf("imgBranchingCorrection(*img, %d, %d, *status)\n",
    mode, verbose);
  /* Check for arguments */
  if(status!=NULL) strcpy(status, "invalid input");
  if(image->status!=IMG_STATUS_OCCUPIED) return(1);
  if(image->isotopeHalflife<=0.0) {
    if(verbose>0) printf("Error: unknown isotope.\n");
    if(status!=NULL) strcpy(status, "unknown isotope");
    return(2);
  }

  /* If branching factor not found in IMG, then get it based on half-life */
  bf=image->branchingFraction;
  if(bf<=0.0 || bf>=1.0) {
    isotope=hlIsotopeFromHalflife(image->isotopeHalflife/60.0);
    if(verbose>2) printf("isotope := %d\n", isotope);
    bf=branchingFraction(isotope);
  }
  if(bf<=0.0 || bf>=1.0) {
    if(verbose>0) printf("Error: branching fraction unknown for the isotope.\n");
    if(status!=NULL) strcpy(status, "branching fraction unknown for the isotope");
    return(3);
  }

  /* Multiply with BF to remove correction, and divide to correct */ 
  if(mode==0) cf=bf; else cf=1.0/bf;
  /* Process pixel values */
  if(verbose>1) printf("multiplying data by %g\n", cf);
  for(fi=0; fi<image->dimt; fi++) {
    for(pi=0; pi<image->dimz; pi++)
      for(i=0; i<image->dimy; i++)
        for(j=0; j<image->dimx; j++)
          image->m[pi][i][j][fi]*=cf;
  }

  /* Fix header contents, too */
  if(image->calibrationFactor>0.0) image->calibrationFactor*=cf;
  if(mode==0) image->branchingFraction=0.0;
  else image->branchingFraction=bf;

  if(status!=NULL) strcpy(status, "ok");
  return(0);
}
/*****************************************************************************/

/*****************************************************************************/

