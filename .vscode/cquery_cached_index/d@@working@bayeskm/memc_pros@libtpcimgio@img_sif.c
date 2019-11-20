/// @file img_sif.c
/// @author Vesa Oikonen
/// @brief Copying data from IMG to SIF and vice versa.
///
/******************************************************************************/
#include "libtpcimgio.h"
/******************************************************************************/

/******************************************************************************/
/** Set IMG contents based on data in SIF.
   @return Returns 0 if successful.
 */ 
int sif2img(
  /** Pointer to SIF struct from which content is copied to IMG. */
  SIF *sif,
  /** Pointer to IMG. Must be initiated with imgInit(), and allocated
      if frame time or count information is to be copied. */
  IMG *img,
  /** Select whether SIF header contents are copied (1) or not (0) */
  int copy_header,
  /** Select whether SIF frame times are copied (1) or not (0). */
  int copy_frames,
  /** Select whether SIF count contents are copied (1) or not (0). */
  int copy_counts,
  /** Verbose level; if zero, then nothing is printed to stderr or stdout. */
  int verbose
) {
  int fi;

  if(verbose>0) printf("sif2img(sif, img, %d, %d, %d, ...)\n", copy_header,
                        copy_frames, copy_counts);
  if(img==NULL || sif==NULL) return 1;

  if(copy_header) {
    if(verbose>1) printf("  copying header.\n");
    img->scanStart=sif->scantime;
    img->isotopeHalflife=60.0*hlFromIsotope(sif->isotope_name);
    if(strlen(sif->studynr)>0 && strcmp(sif->studynr, ".")!=0)
      strlcpy(img->studyNr, sif->studynr, MAX_STUDYNR_LEN+1);
    else
      strcpy(img->studyNr, "");
  }

  if(copy_frames) {
    if(verbose>1) printf("  copying frame times.\n");
    if(sif->frameNr!=img->dimt) return(3);
    for(fi=0; fi<img->dimt; fi++) {
      img->start[fi]=sif->x1[fi];
      img->end[fi]=sif->x2[fi];
      img->mid[fi]=0.5*(img->start[fi]+img->end[fi]);
    }
  }

  if(copy_counts) {
    if(verbose>1) printf("  copying count data.\n");
    if(sif->frameNr!=img->dimt) return(3);
    if(sif->colNr<4) return(4);
    for(fi=0; fi<img->dimt; fi++) {
      img->prompts[fi]=sif->prompts[fi];
      img->randoms[fi]=sif->randoms[fi];
      img->weight[fi]=sif->weights[fi];
    }
  }

  return 0;
}
/*****************************************************************************/

/*****************************************************************************/
/** Set SIF contents based on data in IMG.
    @return Returns 0 if successful.
 */ 
int img2sif(
  /** Pointer to IMG struct from which content is copied to SIF. */
  IMG *img,
  /** Pointer to SIF. Must be initiated with sifInit(), but SIF is allocated here if necessary. */
  SIF *sif,
  /** Select whether header contents are copied (1) or not copied (0) to SIF. */
  int copy_header,
  /** Select whether frame times are copied (1) or not copied (0) to SIF. */
  int copy_frames,
  /** Select whether counts are copied (1) or not copied (0) to SIF,
      or created if they do not exist (2). */
  int copy_counts,
  /** Verbose level; if zero, then nothing is printed to stderr or stdout. */
  int verbose
) {
  int fi;
  if(verbose>0) printf("img2sif(img, sif, %d, %d, %d, ...)\n", copy_header,
                        copy_frames, copy_counts);

  /* Check the arguments */
  if(img==NULL || sif==NULL) return 1;
  if(img->dimt<1) return(1);

  /* Verify that IMG contains frame times */
  if(!imgExistentTimes(img)) {
    if(verbose>0) printf("  image does not contain frame times.\n");
    /* If not, then frame times cannot be copied */
    copy_frames=0;
    /* and counts can not be created */
    if(copy_counts==2) copy_counts=1;
  }

  /* Check if count data needs to be created */
  if(copy_counts==2 && imgExistentCounts(img)!=0)
    copy_counts=1;

  /* Verify that IMG contains isotope information */
  if(img->isotopeHalflife<=0.0) {
    if(verbose>0) printf("  image does not contain isotope halflife.\n");
    /* not, then count data can not be created */
    if(copy_counts==2) copy_counts=1;
  }

  /* Allocate memory for SIF if necessary */
  if((copy_frames || copy_counts) && sif->frameNr!=img->dimt) {
    if(sifSetmem(sif, img->dimt)!=0) return(3);
  }

  /* copy SIF header */
  if(copy_header) {
    if(verbose>1) printf("  copying header fields.\n");
    sif->scantime=img->scanStart;
    sif->colNr=4;
    sif->version=1;
    strcpy(sif->isotope_name, imgIsotope(img));
    if(strlen(img->studyNr)>0 && strcmp(img->studyNr, ".")!=0)
      strlcpy(sif->studynr, img->studyNr, MAX_STUDYNR_LEN+1);
    else
      strcpy(sif->studynr, "");
  }

  /* copy frame times */
  if(copy_frames) {
    if(verbose>1) printf("  copying frame times.\n");
    for(fi=0; fi<img->dimt; fi++) {
      sif->x1[fi]=img->start[fi]; sif->x2[fi]=img->end[fi];
    }
  }

  /* copy or create counts, if required */
  if(copy_counts==2) {
    int i, j, k, pxlNr;
    double v;
    if(verbose>1) printf("  creating count data.\n");
    /* Calculate average of pixel values */
    pxlNr=img->dimz*img->dimx*img->dimy;
    for(fi=0; fi<img->dimt; fi++) {
      v=0.0;
      for(k=0; k<img->dimz; k++)
        for(j=0; j<img->dimy; j++)
	  for(i=0; i<img->dimx; i++){
	    v+=img->m[k][j][i][fi];
	  }
      sif->trues[fi]=v/(double)pxlNr;
    }
    /* Multiply by frame durations, unless we have raw data */
    if(img->type!=IMG_TYPE_RAW)
      for(fi=0; fi<img->dimt; fi++)
        sif->trues[fi]*=(img->end[fi]-img->start[fi]);
    /* Remove decay correction, unless we have raw data */
    if(img->type!=IMG_TYPE_RAW && (img->decayCorrection==IMG_DC_UNKNOWN
         || img->decayCorrection==IMG_DC_CORRECTED))
    {
      double lambda, cf, dur;
      lambda=-hl2lambda(img->isotopeHalflife);
      for(fi=0; fi<img->dimt; fi++) {
        dur=img->end[fi]-img->start[fi];
        cf=hlLambda2factor(lambda, img->start[fi], dur);
        if(cf>0.0) sif->trues[fi]*=cf;
      }
    }
    /* Scale values to a realistic level */
    v=sif->trues[0];
    for(fi=1; fi<sif->frameNr; fi++) if(sif->trues[fi]>v) v=sif->trues[fi];
    v=2.0E+07/v;
    for(fi=0; fi<sif->frameNr; fi++) sif->trues[fi]*=v;
    /* Prompts are set to Trues and Randoms are set to zero */
    for(fi=0; fi<sif->frameNr; fi++) {
      sif->prompts[fi]=sif->trues[fi];
      sif->randoms[fi]=0.0;
      if(sif->trues[fi]<1.0) sif->trues[fi]=1.0;
      sif->weights[fi]=img->weight[fi];
    }
  } else if(copy_counts!=0) {
    if(verbose>1) printf("  copying count data.\n");
    for(fi=0; fi<img->dimt; fi++) {
      sif->prompts[fi]=img->prompts[fi];
      sif->randoms[fi]=img->randoms[fi];
      sif->trues[fi]=sif->prompts[fi]-sif->randoms[fi];
      if(sif->trues[fi]<1.0) sif->trues[fi]=1.0;
      sif->weights[fi]=img->weight[fi];
    }
  }

  return 0;
}
/******************************************************************************/

/******************************************************************************/
