/// @file imgframe.c
/// @author Vesa Oikonen
/// @brief Functions for setting image frame times.
///
/*****************************************************************************/
#include "libtpcimgp.h"
/*****************************************************************************/

/*****************************************************************************/
/** Check for gaps or overlap between frame times. Gap before the first frame is ignored.
    @sa imgExistentTimes, imgDeleteFrameOverlap, imgTimeIntegral, imgFrameGapFill
    @return 0, if no overlaps or gaps are found, 1 if overlaps are found, 
     2 if gaps are found, and 3 if both overlaps and gaps are found.
 */
int imgFramesCheck(
  /** Pointer to IMG struct containing the 4D image data. Data is not modified. */
  IMG *img,
  /** Verbose level; if zero, then nothing is printed to stderr or stdout. */
  int verbose
) {
  if(verbose>0) printf("imgFramesCheck(*img)\n");
  if(img==NULL || img->status!=IMG_STATUS_OCCUPIED || img->dimt<2) return(0);

  int gapNr=0, overlapNr=0;
  for(int fi=1; fi<img->dimt; fi++) {
    float gap=img->start[fi] - img->end[fi-1];
    if(verbose>2 && fabs(gap)>1.0E-06)
      printf("gap between frames %d and %d: %g\n", fi, fi+1, gap);
    if(gap>1.0E-06) gapNr++;
    else if(gap<-1.0E-06) overlapNr++;
  }
  if(verbose>1) {
    printf("  %d overlap(s)\n", overlapNr);
    printf("  %d gap(s)\n", gapNr);
    fflush(stdout);
  }
  int ret=0;
  if(overlapNr>0) ret+=1;
  if(gapNr>0) ret+=2;
  return(ret);
}
/*****************************************************************************/

/*****************************************************************************/
/** Fill gaps between time frames by extending adjacent frames over the gap.
    Overlaps, and gap before the first frame is ignored.
    @sa imgFramesCheck, imgDeleteFrameOverlap, imgExistentTimes, imgTimeIntegral
    @return 0 if successful, >0 in case of an error.
 */
int imgFrameGapFill(
  /** Pointer to IMG struct containing the 4D image data. */
  IMG *img,
  /** Verbose level; if zero, then nothing is printed to stderr or stdout. */
  int verbose
) {
  if(verbose>0) printf("imgFrameGapFill(*img)\n");
  if(img==NULL || img->status!=IMG_STATUS_OCCUPIED || img->dimt<2) return(0);

  for(int fi=1; fi<img->dimt; fi++) {
    float gap=img->start[fi] - img->end[fi-1];
    if(gap<1.0E-07) continue;
    if(verbose>2) printf("gap between frames %d and %d: %g\n", fi, fi+1, gap);
    img->start[fi] -= 0.5*gap; img->mid[fi]=0.5*(img->start[fi]+img->end[fi]);
    img->end[fi-1]=img->start[fi]; img->mid[fi-1]=0.5*(img->start[fi-1]+img->end[fi-1]);
  }
  return(0);
}
/*****************************************************************************/

/*****************************************************************************/
/** Correct frame times if frames are slightly overlapping or
    have small gaps in between.
    Large gap is not corrected and it does not lead to an error.
    @sa imgFramesCheck, imgFrameGapFill, imgTimeIntegral
    @return If overlap is considerable (>1 s), or other error is encountered,
    function returns a non-zero value. Otherwise 0 is returned.
 */
int imgDeleteFrameOverlap(
  /** Pointer to IMG struct containing the 4D image data. */
  IMG *img
) {
  int fi;
  float overlap, overlap_limit=1.8, flen1, flen2;

  if(IMG_TEST) fprintf(stdout, "imgDeleteFrameOverlap()\n");
  if(img->status!=IMG_STATUS_OCCUPIED || img->dimt<1) return(1);
  for(fi=0; fi<img->dimt-1; fi++) {
    overlap=img->end[fi] - img->start[fi+1];
    if(overlap==0.0) continue; // no gap or overlap
    else if(overlap<-overlap_limit) continue; // gap is large, then do nothing
    else if(overlap>overlap_limit) return(2); // overlap is large: error
    /* Correct the small gap/overlap by making frame durations more similar */
    flen1=img->end[fi]-img->start[fi]; flen2=img->end[fi+1]-img->start[fi+1];
    if(overlap>0.0) { // overlap
      if(flen1>flen2) img->end[fi]=img->start[fi+1]; else img->start[fi+1]=img->end[fi];
    } else { // gap
      if(flen1>flen2) img->start[fi+1]=img->end[fi]; else img->end[fi]=img->start[fi+1];
    }
  }
  return(0);
}
/*****************************************************************************/

/*****************************************************************************/
/** Correct frame times so that frames are not overlapping.
    @sa imgDeleteFrameOverlap, imgTimeIntegral
    @return If overlap is considerable (>1 s), or other error is encountered,
    function returns a non-zero value. Otherwise 0 is returned.
 */
int imgDeleteFrameOverlap_old(
  /** Pointer to IMG struct containing the 4D image data. */
  IMG *img
) {
  int fi;
  float overlap;

  if(IMG_TEST) fprintf(stdout, "imgDeleteFrameOverlap()\n");
  if(img->status!=IMG_STATUS_OCCUPIED || img->dimt<1) return(1);
  for(fi=0; fi<img->dimt-1; fi++) {
    overlap=img->end[fi] - img->start[fi+1];
    if(overlap==0.0) continue;
    else if(overlap>1.0) return(2);
    img->end[fi]=img->start[fi+1];
  }
  return(0);
}
/*****************************************************************************/

/*****************************************************************************/
/** Smooth dynamic image data over specified number of time frames.
    @details Average is weighted by frame durations. Gaps or overlaps
    in frame times are not taken into account.
    Do not use this for quantitative analysis, but only for robust peak search etc.
    @return Non-zero value, if error is encountered, otherwise 0 is returned.
 */
int imgSmoothOverFrames(
  /** Pointer to IMG struct containing the 4D image data. */
  IMG *img, 
  /** Nr of frames to average; n must be an odd number and at least 3. */
  int n
) {
  int zi, yi, xi, fi, fj, m, f1, f2;

  if(IMG_TEST) fprintf(stdout, "imgSmoothOverFrames(img, %d)\n", n);
  if(img->status!=IMG_STATUS_OCCUPIED) return(1);
  if(n<3) n=3; else if((n%2)==0) return(1);
  if(img->dimt<n) return(0); // too few frames for smoothing
  m=n/2;
  double orig[img->dimt], vsum, fsum, fdur;
  for(zi=0; zi<img->dimz; zi++) {
    for(yi=0; yi<img->dimy; yi++) for(xi=0; xi<img->dimx; xi++) {
      /* preserve original data for now */
      for(fi=0; fi<img->dimt; fi++) orig[fi]=img->m[zi][yi][xi][fi];
      /* frame-by-frame */
      for(fi=0; fi<img->dimt; fi++) {
        /* set frame range */
        f1=fi-m; if(f1<0) f1=0; 
        f2=fi+m; if(f2>img->dimt-1) f2=img->dimt-1;
        /* mean */
        fsum=vsum=0.0;
        for(fj=f1; fj<=f2; fj++) {
          fdur=img->end[fj]-img->start[fj];
          vsum+=fdur*orig[fj];
          fsum+=fdur;
        }
        if(fsum<1.0E-010) return(2);
        img->m[zi][yi][xi][fi]=vsum/fsum;
      }
    }
  }
  return(0);
}
/*****************************************************************************/

/*****************************************************************************/

