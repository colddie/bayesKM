/// @file weight_model.c
/// @brief Weights for PET data modelling.
/// @author Vesa Oikonen
///
/*****************************************************************************/

/*****************************************************************************/
#include "libtpcmodext.h"
/*****************************************************************************/

/*****************************************************************************/
/** Add weights based on sample frequency or frame length.

    Existing weights are overwritten.

    @return Returns 0 when successful, otherwise <>0.
    @sa dftRead, dftSortByFrame, sifWeight, sifWeightByFrames, 
        dftDecayCorrection, imgSetWeights, dftWSampleNr
 */
int dftWeightByFreq(
  /** Samples/frames must be sorted by sample time, but duplicate samples
      are allowed */
  DFT *dft
) {
  int fi, fi1, fi2;
  double f, t, t1, t2, sumw=0.0;

  if(dft==NULL || dft->frameNr<1) return 1;
  if(dft->frameNr==1) {
    dft->w[0]=1.0;
    dft->isweight=1;
    return 0;
  }

  if(dft->timetype==DFT_TIME_STARTEND) { /* weights by frame lengths */

    for(fi=0; fi<dft->frameNr; fi++) {
      dft->w[fi]=dft->x2[fi]-dft->x1[fi]; 
      sumw+=dft->w[fi];
    }
    /* scale weights so that sum of weights equals sample number */
    sumw/=(double)dft->frameNr;
    for(fi=0; fi<dft->frameNr; fi++) dft->w[fi]/=sumw;

  } else if(dft->timetype==DFT_TIME_MIDDLE) { /* weights by sample distance */

    for(fi=0; fi<dft->frameNr; fi++) {
      t=t1=t2=dft->x[fi]; //printf("fi=%d t=%g\n", fi, t);
      /* Find the closest sample time before this one */
      for(fi1=fi; fi1>=0; fi1--) {t1=dft->x[fi1]; if(t1<t) break;} 
      /* Find the closest sample time after this one */ 
      for(fi2=fi; fi2<dft->frameNr; fi2++) {t2=dft->x[fi2]; if(t2>t) break;} 
      //printf("  t1=%g t2=%g\n", t1, t2);
      /* Mean sample distance */
      f=0.0;
      if(t1<t) f+=t-t1; else f+=t2-t;
      if(t2>t) f+=t2-t; else f+=t-t1;
      f*=0.5; if(f<=0.0) f=1.0;
      /* Set initial weight */
      //printf("    f=%g\n", f);
      dft->w[fi]=f; sumw+=f;
    }
    /* Scale weights */
    sumw/=(double)dft->frameNr;
    for(fi=0; fi<dft->frameNr; fi++) dft->w[fi]/=sumw;

  } else
    return 2;

  dft->isweight=1;
  return(0);
}
/*****************************************************************************/

/*****************************************************************************/
/** Add weights based on average voxel values and frame lengths, just frame 
    lengths, or set all weights to 1 (no weighting).

    Existing weights are overwritten.
    Decay correction is not considered in calculation of weights.

    @return Returns 0 when successful, otherwise <>0.
    @sa dftWeightByFreq, sifWeight, sifWeightByFrames 
 */
int imgSetWeights(
  /** Pointer to IMG struct; image data is used to calculate relational
      frame weights, which are also written into IMG struct. */
  IMG *img,
  /** Weighting method: 0=based on mean voxel values times frame lengths;
      1=based on frame lengths, 2=set weights to 1 (no weighting). */
  int wmet,
  /** Verbose level; if zero, then only warnings are printed into stderr */
  int verbose
) {
  if(verbose>0) printf("imgSetWeights(*img, %d, ...)\n", wmet);
  /* Check the function input */
  if(img==NULL) return(1);
  if(img->status!=IMG_STATUS_OCCUPIED) return(2);
  if(img->dimt<1 || img->dimx<1 || img->dimy<1 || img->dimz<1) return(3);
  if(wmet<0 || wmet>2) return(4);
  img->isWeight=0;

  /* If one time frame, then weight is 1 */
  if(img->dimt==1) {
    img->weight[0]=1.0; 
    if(wmet==0 || wmet==1) img->isWeight=1; 
    return(0);
  }
  
  int fi, ret;
  double f, fm;

  if(wmet==2) { /* If weights are to be set to 1 (removed) */
    for(fi=0; fi<img->dimt; fi++) img->weight[fi]=1.0;
    return(0);
  }

  if(wmet==0) { /* weight based on frame mean values and durations */
    SIF sif;
    sifInit(&sif);
    ret=sifAllocateWithIMG(&sif, img, 1, verbose); if(ret) return(100+ret);
    sifModerateTrues(&sif, 100.0);
    sifWeight(&sif, 0.0); if(verbose>2) sifPrint(&sif);
    for(fi=0; fi<img->dimt; fi++) img->weight[fi]=sif.weights[fi];
    sifEmpty(&sif);
  } else if(wmet==1) { /* weight based on frame durations */
    for(fi=0, fm=0.0; fi<img->dimt; fi++) {
      f=img->end[fi]-img->start[fi]; fm+=f;
      img->weight[fi]=f;
    }
    /* scale weights so that sum of weights equals sample number */
    if(fm>0.0) f=(double)img->dimt/fm; else f=1.0;
    for(fi=0; fi<img->dimt; fi++) img->weight[fi]*=f;
  }
  img->isWeight=1;
  return(0);
}
/*****************************************************************************/

/*****************************************************************************/
/** Get the number of samples in DFT that have weight > 0.
    Missing (NaN) sample values are included as long as weight is not missing.
    If weights are not set, then nr of all samples is returned.
   @return The number of samples with weight above zero.
   @sa dftWeightByFreq, aicSS, parFreeNr
*/
int dftWSampleNr(
  /** Pointer to the DFT struct. */
  DFT *tac
) {
  if(tac==NULL || tac->voiNr<1 || tac->frameNr<1) return(0);
  if(tac->isweight==0) return(tac->frameNr);
  int n=0;
  for(int i=0; i<tac->frameNr; i++) if(tac->w[i]>0.0) n++;
  return(n);
}
/*****************************************************************************/

/*****************************************************************************/
