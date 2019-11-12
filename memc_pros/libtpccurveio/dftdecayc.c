/// @file dftdecay.c
/// @author Vesa Oikonen
/// @brief Physical decay and isotopes in DFT.
///
/*****************************************************************************/
#include "libtpccurveio.h"
/*****************************************************************************/

/*****************************************************************************/
/** Corrects TAC data for physical decay, or removed the correction.

    @details Weights are not modified.

    @return Returns 0 if conversion was successful, and <> 0 if failed.
 */   
int dftDecayCorrection(
  /** Pointer to existing DFT data; status of decay correction in DFT
   *  is not verified, but set in this function;
   *  DFT must contain valid sample time unit. */
  DFT *dft,
  /** Half-life of isotope in minutes; enter <=0, if correct isotope code
   *  is given in DFT */
  double hl,
  /** 0=Remove decay correction; 1=Correct for decay */
  int mode,
  /** Apply (1) or do not apply (0) correction to y[] data in DFT */
  int y,
  /** Apply (1) or do not apply (0) correction to y2[] data in DFT */
  int y2,
  /** Apply (1) or do not apply (0) correction to y3[] data in DFT */
  int y3,
  /** Pointer to a string (allocated for at least 64 chars) where error message
      or other execution status will be written; enter NULL, if not needed */   
  char *status,
  /** Verbose level; if zero, then nothing is printed to stderr or stdout */
  int verbose
) {
  int ri, fi, isotopeID=-1;
  double lambda, dc;

  if(verbose>0) printf("dftDecayCorrection(dft, %g, %d, %d, %d, %d, ...)\n",
    hl, mode, y, y2, y3);

  /* Check the input */
  if(status!=NULL) strcpy(status, "invalid input");
  if(dft==NULL) return(1);
  if(dft->voiNr<1 || dft->frameNr<1) return(1);
  
  /* Check that DFT contains time unit */
  if(dft->timeunit!=TUNIT_SEC && dft->timeunit!=TUNIT_MIN && 
     dft->timeunit!=TUNIT_HOUR) 
  {
    if(verbose>0) printf("dft->timeunit := %d\n", dft->timeunit);
    if(status!=NULL) strcpy(status, "sample time unit is not specified");
    return(2);
  }
  
  /* Check and set isotope code and half-life */
  if(hl>1.0E-10) { 
    /* Try to identify the isotope from given half-life */
    isotopeID=hlIsotopeFromHalflife(hl);
    if(isotopeID>=0) {
      /* If identified, then set the code in DFT */
      strcpy(dft->isotope, hlIsotopeCode(isotopeID));
      if(verbose>1) printf("  isotope := %s\n", dft->isotope);
    } else {
      /* If not identified, that may just be because isotope is not yet listed
         in TPC library, but give a warning */
      fprintf(stderr, "Warning: halflife %g min is not identified.\n", hl);
    }
  } else {
    /* Check that valid isotope code is found in DFT */
    hl=hlFromIsotope(dft->isotope); if(hl<=0.0) {
      if(verbose>0) printf("dft->isotope := %s\n", dft->isotope);
      if(status!=NULL) strcpy(status, "valid isotope is not specified");
      return(11);
    }
    if(verbose>1) printf("  half-life := %g min\n", hl);
  }
  
  /*
   *  Calculate the lambda
   */
  if(dft->timeunit==TUNIT_SEC) hl*=60.0;
  else if(dft->timeunit==TUNIT_HOUR) hl/=60.0;
  lambda=hl2lambda(hl); if(mode==0) lambda=-lambda;
  if(verbose>1) printf("lambda := %e\n", lambda);

  /*
   *  Decay correction / removal
   */
  if(verbose>2) {
    if(mode!=0) printf("decay correction\n");
    else printf("removing decay correction\n");
  }
  for(fi=0; fi<dft->frameNr; fi++) {
    /* Calculate decay factor */
    if(dft->timetype==DFT_TIME_STARTEND) {
      if(isnan(dft->x1[fi]) || isnan(dft->x2[fi])) continue;
      dc=hlLambda2factor(lambda, dft->x1[fi], dft->x2[fi]-dft->x1[fi]);
    } else {
      if(isnan(dft->x[fi])) continue;
      dc=hlLambda2factor(lambda, dft->x[fi], 0.0);
    }
    if(verbose>4) printf("  %10.4f  ->  %e\n", dft->x[fi], dc);
    /* Correct all regions */
    for(ri=0; ri<dft->voiNr; ri++) {
      if(y!=0 && !isnan(dft->voi[ri].y[fi])) dft->voi[ri].y[fi]*=dc;
      if(y2!=0 && !isnan(dft->voi[ri].y2[fi])) dft->voi[ri].y2[fi]*=dc;
      if(y3!=0 && !isnan(dft->voi[ri].y3[fi])) dft->voi[ri].y3[fi]*=dc;
    }
  }
  if(mode!=0) {
    dft->decayCorrected=DFT_DECAY_CORRECTED;
    if(status!=NULL) strcpy(status, "decay corrected");
  } else {
    dft->decayCorrected=DFT_DECAY_NOTCORRECTED;
    if(status!=NULL) strcpy(status, "decay correction removed");
  }
  return(0);
} 
/*****************************************************************************/

/*****************************************************************************/
