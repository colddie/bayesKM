/// @file weight.c
/// @author Vesa Oikonen
/// @brief Functions for setting weight factors based on SIF.
///
/*****************************************************************************/
#include "libtpcimgio.h"
/*****************************************************************************/

/*****************************************************************************/
/** @brief Calculate weights for frames in SIF data based on true counts. 
     Weights are normalized to have an average of 1.0.

  @details Weights are calculated from formula
   weight=(frame duration)^2 / (trues in a frame).
   Before calling this routine, trues must be calculated as total counts -
   randoms.
   Counts in SIF are not corrected for physical decay. Therefore, isotope
   halflife must be known, if weights are to be calculated for decay
   corrected TACs. Isotope halflife must be set to 0, if weights are used
   for TACs that are not corrected for decay.
  
   Reference: Mazoyer BM, Huesman RH, Budinger TF, Knittel BL. 
   Dynamic PET data analysis. J Comput Assist Tomogr 1986; 10:645-653.

  @sa sifInit, sifRead, sifWrite, sifEmpty, sifModerateTrues, sifModerateWeights,
      sifExistentCounts, sifWeightByFrames, dftWeightByFreq, dftDecayCorrection
 */
void sifWeight(
  /** Pointer to SIF data. */
  SIF *data,
  /** Halflife (in seconds) of the isotope;
      If halflife is 0, the weights are calculated for non-decay corrected data.
      If halflife is >0, the weights are calculated using decay corrected
      trues, but trues data in SIF are not changed. */
  double halflife
) {
  int i;
  double f, d;

  if(SIF_TEST) printf("sifWeight(*sif, %g)\n", halflife);
  /* Calculate weights */
  for(i=0; i<data->frameNr; i++) {
    if(data->trues[i]<1.0) data->trues[i]=1.0;
    f=data->x2[i]-data->x1[i]; if(f<=0.0) f=1.0;
    if(halflife<=1.0E-8) 
      d=1.0;
    else
      d=exp( ((data->x1[i]+data->x2[i])/2.0)*0.693147/halflife );
    data->weights[i]=(f*f)/(d*data->trues[i]);
    /*printf("%3d %g %g\n", i, data->trues[i], data->weights[i]);*/
  }

  /* Scale weights so that average weight is 1.0 */
  sifWeightNorm(data);

  return;
}
/*****************************************************************************/

/*****************************************************************************/
/** @brief Calculate weights for frames in SIF data based on frame lengths.
     Weights are normalized to have an average of 1.0.

  @details Weights are calculated from the simple formula
   weight=(frame duration), if isotope halflife is set to 0, and using formula
   weight=(frame duration)*exp(-t*ln(2)/halflife), if isotope halflife is
   given (less weight for late frames, which may be more suitable for 
   PET data that is corrected for physical decay).

  @sa sifInit, sifRead, sifWrite, sifEmpty, sifWeight, dftWeightByFreq
 */
void sifWeightByFrames(
  /** Pointer to SIF data. */
  SIF *data,
  /** Halflife (in seconds) of the isotope;
      If halflife is >0, the late frames are given less weight. */
  double halflife
) {
  int i;
  double f, d;

  if(SIF_TEST) printf("sifWeightByFrames(*sif, %g)\n", halflife);
  /* Calculate weights */
  for(i=0; i<data->frameNr; i++) {
    f=data->x2[i]-data->x1[i]; if(f<=0.0) f=1.0;
    if(halflife<=1.0E-8) 
      d=1.0;
    else
      d=exp( -((data->x1[i]+data->x2[i])/2.0)*0.693147/halflife );
    data->weights[i]=f*d;
  }

  /* Scale weights so that average weight is 1.0 */
  sifWeightNorm(data);

  return;
}
/*****************************************************************************/

/*****************************************************************************/
/** Scale weights in SIF data so that average weight is 1.0, and the sum of
    weights equals the nr of frames.
  @sa sifInit, sifRead, sifWrite, sifEmpty, sifWeight, sifWeightByFrames
*/
void sifWeightNorm(
  /** Pointer to SIF data. */
  SIF *d
) {
  if(SIF_TEST) printf("sifWeightNorm(*sif)\n");
  int i;
  double f=0.0;
  for(i=0; i<d->frameNr; i++) f+=d->weights[i];
  f/=(double)d->frameNr;
  for(i=0; i<d->frameNr; i++) d->weights[i]/=f;
  return;
}
/*****************************************************************************/

/*****************************************************************************/
/** Moderate the trues in SIF. 

    True values in SIF are used to calculate weight factors for time frames.
    If trues are very low in certain frames, the weight factors in other frames
    may become very low.
    This function finds the maximum trues, and adds max/limit to each trues
    value, if min trues < max trues / limit.
    Negative trues are always eliminated.

    @sa sifWeight, sifModerateWeights, sifRead, sifWrite, sifWeightNorm
 */
void sifModerateTrues(
  /** Pointer to SIF in which the trues are moderated */
  SIF *sif,
  /** Max trues / limit is added to all trues values; 100.0 might be good */ 
  double limit
) {
  if(SIF_TEST) printf("sifModerateTrues(*sif, %g)\n", limit);

  if(sif==NULL || sif->frameNr<2) return;
  if(limit<=1.0) return;

  int fi;
  double w, f;
  for(w=f=sif->trues[0], fi=1; fi<sif->frameNr; fi++) {
    if(sif->trues[fi]>w) w=sif->trues[fi];
    else if(sif->trues[fi]<f) f=sif->trues[fi];
  }
  if(f*limit<w) {
    for(w/=limit, fi=0; fi<sif->frameNr; fi++)
      if(sif->trues[fi]>0.0) sif->trues[fi]+=w; else sif->trues[fi]=w;
  } else {
    for(fi=0; fi<sif->frameNr; fi++)
      if(sif->trues[fi]<0.0) sif->trues[fi]=0.0;
  }
}
/*****************************************************************************/

/*****************************************************************************/
/** Moderate the weights in SIF. 

    This function finds the maximum weight, and adds max/limit to each weight
    value (except if weight is 0), if min weight < max/limit.
    Negative weights are set to zero.

    @sa sifWeight, sifModerateTrues, sifRead, sifWrite, sifWeightNorm
 */
void sifModerateWeights(
  /** Pointer to SIF in which the weights are moderated */
  SIF *sif,
  /** Max weight / limit is added to all weights; 100.0 might be good */ 
  double limit
) {
  if(SIF_TEST) printf("sifModerateWeights(*sif, %g)\n", limit);

  if(sif==NULL || sif->frameNr<2) return;
  if(limit<=1.0) return;

  int fi;
  double w, f;
  w=f=nan("");
  for(fi=0; fi<sif->frameNr; fi++) {
    if(isnan(sif->weights[fi])) continue;
    if(sif->weights[fi]<=0.0) {sif->weights[fi]=0.0; continue;}
    if(isnan(w)) w=sif->weights[fi];
    if(isnan(f)) f=sif->weights[fi];
    if(sif->weights[fi]>w) w=sif->weights[fi];
    if(sif->weights[fi]<f) f=sif->weights[fi];
  }
  if(isnan(w) || isnan(f)) return;

  if(f*limit<w) {
    for(w/=limit, fi=0; fi<sif->frameNr; fi++)
      if(sif->weights[fi]>0.0) sif->weights[fi]+=w;
  }
}
/*****************************************************************************/

/*****************************************************************************/
/** Verify that SIF contains prompts and randoms.

  @return Returns 0 if neither prompts or randoms can be found, 
    1 or 2 if prompts or randoms can be found, and 
    3 if both prompts and randoms are there.

  @sa sifInit, sifRead, sifWeight
 */
int sifExistentCounts(
  /** Pointer to SIF struct */
  SIF *sif
) {
  int fi, p=0, r=0;
  double v1, v2;
  if(sif==NULL || sif->frameNr<1) return 0;
  /* If just one frame, then value > 0 is fine */
  if(sif->frameNr==1) {
    if(sif->prompts[0]>0.00000001) p=1;
    if(sif->randoms[0]>0.00000001) r=2;
    return(p+r);
  }
  /* Otherwise, check also that frames have different count level */
  for(fi=1; fi<sif->frameNr; fi++) {
    v1=sif->prompts[fi]-sif->prompts[fi-1]; if(fabs(v1)>0.001) p=1;
    v2=sif->randoms[fi]-sif->randoms[fi-1]; if(fabs(v2)>0.001) r=2;
    if((p+r)>2) break;
  }
  return(p+r);
}
/*****************************************************************************/

/*****************************************************************************/
