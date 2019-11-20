/// @file imgeval.c
/// @brief Functions for extracting TACs from image data.
/// @author Vesa Oikonen
///
/*****************************************************************************/
#include "libtpcimgp.h"
/*****************************************************************************/

/*****************************************************************************/
/** Calculates an average time-activity curve of all pixels or bins in the specified IMG data.

    @sa imgAverageMaskTAC, imgAverageAUC, imgRangeMinMax
    @return Returns 0, if ok.
 */
int imgAverageTAC(
  /** (Dynamic) IMG data from which TAC is computed. */
  IMG *img,
  /** Allocated float array for the TAC. */
  float *tac
) {
  return(imgAverageMaskTAC(img, NULL, tac));
}
/*****************************************************************************/

/*****************************************************************************/
/** Calculates an average time-activity curve of pixels or bins in the specified IMG data.

    Mask image specifies the pixels that are included in the average.
    If all pixels are to be averaged, then NULL can be given instead of mask image.

    @sa imgMaskTAC, imgThresholdingLowHigh, imgThresholdMaskCount, imgAverageAUC, imgAverageTAC
    @return Returns 0, if ok.
 */
int imgAverageMaskTAC(
  /** (Dynamic) IMG data from which cluster TAC is computed. */
  IMG *img,
  /** Mask: 0=excluded, otherwise included.
      Enter NULL to include all pixels in the average. */
  IMG *timg,
  /** Allocated float array for the TAC. */
  float *tac
) {
  int zi, yi, xi, fi, pxlNr;
  
  if(img==NULL || img->status<IMG_STATUS_OCCUPIED) return(1);
  if(tac==NULL) return(2);
  if(timg!=NULL) {
    if(timg->status<IMG_STATUS_OCCUPIED) return(3);
    /* check that image dimensions are the same */
    if(timg->dimz!=img->dimz || timg->dimz<1) return(5);
    if(timg->dimy!=img->dimy || timg->dimy<1) return(5);
    if(timg->dimx!=img->dimx || timg->dimx<1) return(5);
  }
  /* compute the nr of pixels in the mask */
  if(timg==NULL) {
    pxlNr=img->dimz*img->dimy*img->dimx; 
  } else {
    pxlNr=0;
    for(zi=0; zi<timg->dimz; zi++)
      for(yi=0; yi<timg->dimy; yi++) 
        for(xi=0; xi<timg->dimx; xi++)
          if(timg->m[zi][yi][xi][0]!=0.0) pxlNr++;
  }
  if(pxlNr<1) return(4);

  /* compute the TAC */
  for(fi=0; fi<img->dimt; fi++) {
    tac[fi]=0.0;
    for(zi=0; zi<img->dimz; zi++) {
      for(yi=0; yi<img->dimy; yi++) {
        for(xi=0; xi<img->dimx; xi++) {
          if(timg==NULL) tac[fi]+=img->m[zi][yi][xi][fi];
          else if(timg->m[zi][yi][xi][0]!=0.0) tac[fi]+=img->m[zi][yi][xi][fi];
        }
      }
    }
    tac[fi]/=(float)pxlNr;
    /*printf("pxlNr=%d/%d\n", pxlNr, img->dimz*img->dimy*img->dimx);*/
  }
  return(0);
}
/*****************************************************************************/

/*****************************************************************************/
/** Calculates the Area-Under-Curve of an average time-activity curve of all
    pixels or bins in the specified IMG data.

    @sa imgAverageMaskTAC, imgAverageTAC, imgAvg
    @return Returns 0, if ok.
 */
int imgAverageAUC(
  /** (Dynamic) IMG data. */
  IMG *img,
  /** Pointer to a float were the average TAC AUC will be written. */
  float *avgauc
) {
  int pi, yi, xi, fi, pxlNr;
  float fv, fl;
  
  if(img->status<IMG_STATUS_OCCUPIED) return(1);
  if(avgauc==NULL) return(2);
  *avgauc=0.0;
  pxlNr=img->dimz*img->dimy*img->dimx; if(pxlNr<1) return(3);
  for(fi=0; fi<img->dimt; fi++) {
    fv=0.0; fl=img->end[fi]-img->start[fi]; if(fl<=0.0) return(4);
    for(pi=0; pi<img->dimz; pi++)
      for(yi=0; yi<img->dimy; yi++)
        for(xi=0; xi<img->dimx; xi++)
          fv+=img->m[pi][yi][xi][fi];
    fv/=(float)pxlNr;
    *avgauc+=fv*fl;
  }
  return(0);
}
/*****************************************************************************/

/*****************************************************************************/
/** Calculate TAC as weighted average of voxels in specified image data
    with relative weights given in a mask image.

    @sa imgsegmThresholdMask, imgAverageMaskTAC, imgThresholdingLowHigh,
        imgThresholdMask

    @return Returns 0 if successful.
 */
int imgMaskTAC(
  /** Pointer to allocated image from which the weighted TAC is calculated. */
  IMG *img,
  /** Pointer to mask image; x, y, and z dimensions must be the same as in
      the image to which the mask is applied. */
  IMG *mask,
  /** Pointer to an array where weighted pixel averages are written;
      it must be allocated for size >= dimt. */
  double *tac,
  /** Verbose level; set to <=0 to prevent all prints to stdout. */
  int verbose
) {
  if(verbose>0) printf("imgMaskTAC()\n");

  int xi, yi, zi, fi;
  double w;

  if(img->status<IMG_STATUS_OCCUPIED) return(1);
  if(mask->status<IMG_STATUS_OCCUPIED) return(2);
  if(mask->dimz!=img->dimz) return(3);
  if(mask->dimy!=img->dimy) return(4);
  if(mask->dimx!=img->dimx) return(5);
  if(img->dimt<1 || mask->dimt<1) return(6);
  if(tac==NULL) return(7);

  /* Initiate TAC */
  for(fi=0; fi<img->dimt; fi++) tac[fi]=0.0;
  /* Add weighted sum to TAC */
  w=0.0;
  for(zi=0; zi<mask->dimz; zi++)
    for(yi=0; yi<mask->dimy; yi++)
      for(xi=0; xi<mask->dimx; xi++) if(mask->m[zi][yi][xi][0]>0.0) {
        for(fi=0; fi<img->dimt; fi++)
          tac[fi]+=mask->m[zi][yi][xi][0]*img->m[zi][yi][xi][fi];
        w+=mask->m[zi][yi][xi][0];
      }
  /* Divide sum TAC by sum weights */
  for(fi=0; fi<img->dimt; fi++) tac[fi]/=w;
  if(verbose>1) printf("mask_sum := %g\n", w);
  if(w<=0.0 && verbose>0)
    fprintf(stderr, "Warning: zero mask applied to image.\n");

  return(0);
}
/*****************************************************************************/

/*****************************************************************************/
/** Get the list of ROIs in the mask image.

    @sa imgMaskCount, imgThresholdMaskCount, imgThresholdMask

    @return Returns 0 if successful.
*/
int imgMaskRoiNr(
  /** Pointer to mask image. Pixel values <=0 represent pixels outside any ROI, and 
      each rounded positive integer value represents one ROI. */
  IMG *img,
  /** Initiated list of integers. */
  INTEGER_LIST *list
) {
  int xi, yi, zi, v, ret;

  /* Check the input */
  if(img==NULL || list==NULL) return 1;

  /* Clear any previous list contents */
  if(list->nr>0) {integerListEmpty(list);}

  /* Go through the mask image */
  for(zi=0; zi<img->dimz; zi++)
    for(yi=0; yi<img->dimy; yi++)
      for(xi=0; xi<img->dimx; xi++) {
        v=temp_roundf(img->m[zi][yi][xi][0]);
        if(v>0) {ret=integerListAdd(list, v, 1); if(ret<0) return(100+ret);}
      }
  /* Sort the list */
  integerListSort(list);

  return 0;
}
/*****************************************************************************/

/*****************************************************************************/
/** Calculate TAC as average of voxels in image data, including only voxels
    which have specified value in the mask image.

    @sa imgsegmThresholdMask, imgMaskTAC, imgMaskDilate, imgMaskErode

    @return Returns the nr of voxels included, or <0 in case of an error.
 */
int imgVoiMaskTAC(
  /** Pointer to allocated image from which the mean TAC is calculated. */
  IMG *img,
  /** Pointer to mask image; x, y, and z dimensions must be the same as in
      the image to which the mask is applied. */
  IMG *mask,
  /** Voxels with this value in mask image are included in the average. */
  int mv,
  /** Pointer to an array where weighted pixel averages are written;
      it must be allocated for size >= dimt. */
  double *tac,
  /** Verbose level; set to <=0 to prevent all prints to stdout. */
  int verbose
) {
  if(verbose>0) printf("imgVoiMaskTAC(img, mask, %d, tac, %d)\n", mv, verbose);

  int xi, yi, zi, fi;
  double mf;

  if(img->status<IMG_STATUS_OCCUPIED) return(-1);
  if(mask->status<IMG_STATUS_OCCUPIED) return(-2);
  if(mask->dimz!=img->dimz) return(-3);
  if(mask->dimy!=img->dimy) return(-4);
  if(mask->dimx!=img->dimx) return(-5);
  if(img->dimt<1 || mask->dimt<1) return(-6);
  if(tac==NULL) return(-7);

  /* Initiate TAC */
  for(fi=0; fi<img->dimt; fi++) tac[fi]=0.0;

  /* Add the sum to TAC */
  int pxlNr[img->dimt], badNr[img->dimt];
  for(fi=0; fi<img->dimt; fi++) pxlNr[fi]=badNr[fi]=0;
  for(zi=0; zi<mask->dimz; zi++) {
    for(yi=0; yi<mask->dimy; yi++) {
      for(xi=0; xi<mask->dimx; xi++) {
        mf=(int)temp_roundf(mask->m[zi][yi][xi][0]);
        if(mf!=mv) continue;
        for(fi=0; fi<img->dimt; fi++) {
          if(isnan(img->m[zi][yi][xi][fi])) {badNr[fi]++; continue;}
          tac[fi]+=img->m[zi][yi][xi][fi]; pxlNr[fi]++;
        }
      }
    }
  }

  /* Divide sum TAC by sum weights */
  int minNr=pxlNr[0];
  for(fi=0; fi<img->dimt; fi++) {
    if(badNr[fi]>0 && verbose>1) {
      fprintf(stderr, "Warning: %d missing pixel values in frame %d.\n", badNr[fi], 1+fi); 
      fflush(stderr);
    }
    if(pxlNr[fi]==0) {
      if(verbose>0) {
        fprintf(stderr, "Warning: zero valid pixels in frame %d.\n", 1+fi); fflush(stderr);
      }
      tac[fi]=nan("");
    }
    if(verbose>1) {
      printf("%d valid pixels inside mask in frame %d\n", pxlNr[fi], 1+fi); fflush(stdout);
    }
    tac[fi]/=(double)pxlNr[fi];
    if(pxlNr[fi]<minNr) minNr=pxlNr[fi];
  }

  return(minNr);
}
/*****************************************************************************/

/*****************************************************************************/
