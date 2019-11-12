/// @file img_mtga.c
/// @brief Functions for computing pixel-by-pixel the MTGA
///        (Gjedde-Patlak and Logan plot).
/// @author Vesa Oikonen
///
/*****************************************************************************/

/*****************************************************************************/
#include "libtpcmodext.h"
/*****************************************************************************/

/*****************************************************************************/
/** Computing pixel-by-pixel the graphical analysis for irreversible PET tracers (Gjedde-Patlak plot).
    @sa img_logan, llsqperp, patlak_data, mtga_best_perp
    @return Returns 0 if successful, and >0 in case of an error.
 */
int img_patlak(
  /** Pointer to the TAC data to be used as model input. Sample times in minutes.
      Curve is interpolated to PET frame times, if necessary. */
  DFT *input,
  /** Pointer to dynamic PET image data.
      Image and input data must be in the same calibration units. */
  IMG *dyn_img,
  /** The range of frames where line is fitted, given as the frame start here
      and next the end index, i.e. [0..frame_nr-1]. */ 
  int start,
  /** The range of frames where line is fitted, given as the frame start above
      and here the end index, i.e. [0..frame_nr-1]. */ 
  int end,
  /** Use the whole range or based on data leave out points from the beginning;
      PRESET or EXCLUDE_BEGIN. */
  linefit_range fit_range,
  /** Threshold as fraction of input AUC. */
  float thrs,
  /** Pointer to initiated IMG structure where Ki values will be placed. */
  IMG *ki_img,
  /** Pointer to initiated IMG structure where plot y axis intercept values 
      will be placed; enter NULL, if not needed. */
  IMG *ic_img,
  /** Pointer to initiated IMG structure where the number of plot data points
      actually used in the fit is written; enter NULL, when not needed. */
  IMG *nr_img,
  /** Pointer to a string (allocated for at least 64 chars) where error message
      or other execution status will be written; enter NULL, if not needed. */
  char *status,
  /** Verbose level; if zero, then nothing is printed to stderr or stdout. */
  int verbose
) {
  int zi, yi, xi, fi, nr, pn, ret=0;
  DFT tac;
  double *plotData, *xaxis, *yaxis, slope, ic, f;


  if(verbose>0) {
    printf("%s(input, dyn_img, %d, %d, range, %g, ki_img, ", __func__, start, end, thrs);
    if(ic_img==NULL) printf("NULL, "); else printf("ic_img, ");
    if(nr_img==NULL) printf("NULL, "); else printf("nr_img, ");
    if(status==NULL) printf("NULL)\n"); else printf("status)\n");
  }
  /* Initial check for the arguments */
  if(status!=NULL) sprintf(status, "invalid data");
  if(dyn_img->status!=IMG_STATUS_OCCUPIED || dyn_img->dimt<1) return(1);
  if(input==NULL || input->frameNr<1) return(2);
  nr=1+end-start; if(nr<2) return(3);
  if(end>dyn_img->dimt-1 || start<0) return(4);
  if(ki_img==NULL) return(5);
  /* Convert input time units to min */
  if(input->timeunit==TUNIT_SEC) dftTimeunitConversion(input, TUNIT_MIN);

  /* Check that input contains samples until at least 80% of line fit range */
  if(verbose>1) {
    printf("input_last_sample_time := %g\n", input->x[input->frameNr-1]);
    printf("patlak_start_time := %g\n", dyn_img->start[start]/60.0);
    printf("patlak_end_time := %g\n", dyn_img->end[end]/60.0);
  }
  if(input->x[input->frameNr-1] < (0.2*dyn_img->mid[start]+0.8*dyn_img->mid[end])/60.0) {
    sprintf(status, "too few input samples"); return(6);
  }

  /* Allocate memory for interpolated and integrated input and 
     tissue pixel TACs */
  dftInit(&tac); if(dftSetmem(&tac, nr, 2)!=0) {sprintf(status, "out of memory"); return(11);}
  strcpy(tac.voi[0].voiname, "input"); strcpy(tac.voi[1].voiname, "tissue");
  tac.voiNr=2; tac.frameNr=nr;
  for(fi=0; fi<tac.frameNr; fi++) {
    tac.x1[fi]=dyn_img->start[start+fi]/60.;
    tac.x2[fi]=dyn_img->end[start+fi]/60.;
    tac.x[fi]=dyn_img->mid[start+fi]/60.;
  }

  /* If input sample times are much different from PET frame times,
     as they usually are if arterial plasma is used as input, then use 
     interpolate4pet(), otherwise with image-derived input use petintegral()
  */
  int equal_times;
  equal_times=check_times_dft_vs_img(dyn_img, input, verbose-1);
  if(equal_times==1) {
    if(verbose>1) printf("copying input curve and using petintegral()\n");
    /* Copy frame start and end times from image to input data */
    ret=copy_times_from_img_to_dft(dyn_img, input, verbose-1);
    if(ret) {
      sprintf(status, "cannot set input sample times");
      dftEmpty(&tac); return(12);
    }
    /* integral must be calculated from the zero time */
    ret=petintegral(input->x1, input->x2, input->voi[0].y, input->frameNr,
      input->voi[0].y2, input->voi[0].y3);
    /* because times are the same, the values can be copied directly */
    if(ret==0) for(fi=0; fi<tac.frameNr; fi++) {
      tac.voi[0].y[fi]=input->voi[0].y[start+fi];
      tac.voi[0].y2[fi]=input->voi[0].y2[start+fi];
      tac.voi[0].y3[fi]=input->voi[0].y3[start+fi];
    }
  } else {
    if(verbose>1) printf("using interpolate4pet() for input curve\n");
    ret=interpolate4pet(
      input->x, input->voi[0].y, input->frameNr,
      tac.x1, tac.x2, tac.voi[0].y, tac.voi[0].y2, tac.voi[0].y3, tac.frameNr
    );
  }
  if(ret) {
    dftEmpty(&tac);
    sprintf(status, "cannot interpolate input data"); return(12);
  }
  if(verbose>3) dftPrint(&tac);

  /*
   *  Allocate result images and fill the header info
   */
  /* Ki image */
  imgEmpty(ki_img);
  ret=imgAllocateWithHeader(ki_img, dyn_img->dimz, dyn_img->dimy, dyn_img->dimx, 1, dyn_img);
  if(ret) {
    sprintf(status, "cannot setup memory for Ki image");
    imgEmpty(ki_img); dftEmpty(&tac); return(21);
  }
  ki_img->unit=CUNIT_ML_PER_ML_PER_MIN; /* mL/(mL*min) */
  ki_img->decayCorrection=IMG_DC_NONCORRECTED; ki_img->isWeight=0;
  ki_img->start[0]=dyn_img->start[start]; ki_img->end[0]=dyn_img->end[end];
  if(verbose>9) imgInfo(ki_img);
  /* Ic image */
  if(ic_img!=NULL) {
    imgEmpty(ic_img);
    ret=imgAllocateWithHeader(ic_img, dyn_img->dimz, dyn_img->dimy, dyn_img->dimx, 1, dyn_img);
    if(ret) {
      sprintf(status, "cannot setup memory for Ic image");
      imgEmpty(ic_img); imgEmpty(ki_img); dftEmpty(&tac); return(23);
    }
    ic_img->unit=CUNIT_ML_PER_ML; /* mL/mL */
    ic_img->decayCorrection=IMG_DC_NONCORRECTED; ic_img->isWeight=0;
    ic_img->start[0]=dyn_img->start[start]; ic_img->end[0]=dyn_img->end[end];
    if(verbose>9) imgInfo(ic_img);
  }
  /* Nr image */
  if(nr_img!=NULL) {
    imgEmpty(nr_img);
    ret=imgAllocateWithHeader(nr_img, dyn_img->dimz, dyn_img->dimy, dyn_img->dimx, 1, dyn_img);
    if(ret) {
      sprintf(status, "cannot setup memory for nr image");
      imgEmpty(nr_img); imgEmpty(ic_img); imgEmpty(ki_img); dftEmpty(&tac);
      return(25);
    }
    nr_img->unit=CUNIT_UNITLESS;
    nr_img->decayCorrection=IMG_DC_NONCORRECTED; nr_img->isWeight=0;
    nr_img->start[0]=dyn_img->start[start]; nr_img->end[0]=dyn_img->end[end];
    if(verbose>9) imgInfo(ic_img);
  }

  /* Allocate memory for graphical analysis plot data */
  plotData=malloc(2*nr*sizeof(double));
  if(plotData==NULL) {
    sprintf(status, "cannot allocate memory for plots");
    imgEmpty(ic_img); imgEmpty(ki_img); dftEmpty(&tac); return(25);
  }
  xaxis=plotData; yaxis=plotData+nr;

  /* Calculate threshold */
  thrs*=tac.voi[0].y2[tac.frameNr-1];
  if(verbose>1) printf("  threshold-AUC := %g\n", thrs);

  /*
   *  Compute pixel-by-pixel
   */
  if(verbose>1) printf("computing MTGA pixel-by-pixel\n");
  float pxlauc[dyn_img->dimt];
  int best_nr;
  for(zi=0; zi<dyn_img->dimz; zi++) {
    for(yi=0; yi<dyn_img->dimy; yi++) {
      for(xi=0; xi<dyn_img->dimx; xi++) {
        /* Initiate pixel output values */
        ki_img->m[zi][yi][xi][0]=0.0;
        if(ic_img!=NULL) ic_img->m[zi][yi][xi][0]=0.0;
        /* Check for threshold */
        ret=fpetintegral(dyn_img->start, dyn_img->end, dyn_img->m[zi][yi][xi], 
                         dyn_img->dimt, pxlauc, NULL);
        if(ret) continue;
        if((pxlauc[dyn_img->dimt-1]/60.0) < thrs) continue;
        /* Calculate Patlak plot data */
        for(fi=0; fi<tac.frameNr; fi++) tac.voi[1].y[fi]=dyn_img->m[zi][yi][xi][start+fi];
        pn=patlak_data(nr, tac.voi[0].y, tac.voi[0].y2, tac.voi[1].y, xaxis, yaxis);
        /* Line fit */
        if(fit_range==PRESET || pn<MTGA_BEST_MIN_NR) {
          ret=llsqperp(xaxis, yaxis, pn, &slope, &ic, &f);
          best_nr=pn;
        } else {
          ret=mtga_best_perp(xaxis, yaxis, pn, &slope, &ic, &f, &best_nr);
        }
        if(ret==0) {
          ki_img->m[zi][yi][xi][0]=slope;
          if(ic_img!=NULL) ic_img->m[zi][yi][xi][0]=ic;
          if(nr_img!=NULL) nr_img->m[zi][yi][xi][0]=best_nr;
        }
      } /* next column */
    } /* next row */
  } /* next plane */

  free(plotData); dftEmpty(&tac);
  return(0);
}
/*****************************************************************************/

/*****************************************************************************/
/** Computing pixel-by-pixel the graphical analysis for reversible PET tracers (Logan plot).
    @sa img_patlak, llsqperp, logan_data, mtga_best_perp
    @return Returns 0 if successful, and >0 in case of an error.
 */
int img_logan(
  /** Pointer to the TAC data to be used as model input. Sample times in minutes.
      Curve is interpolated to PET frame times, if necessary. */
  DFT *input,
  /** Pointer to dynamic PET image data.
      Image and input data must be in the same calibration units. */
  IMG *dyn_img,
  /** The range of frames where line is fitted, given as the frame start here
      and next the end index, i.e. [0..frame_nr-1]. */
  int start,
  /** The range of frames where line is fitted, given as the frame start above
      and here the end index, i.e. [0..frame_nr-1]. */
  int end,
  /** Use the whole range or based on data leave out points from the beginning;
      PRESET or EXCLUDE_BEGIN. */
  linefit_range fit_range,
  /** Threshold as fraction of input AUC. */
  float thrs,
  /** Reference region k2; set to <=0 if not needed. */
  double k2,
  /** Pointer to initiated IMG structure where Vt (or DVR) values will be placed. */
  IMG *vt_img,
  /** Pointer to initiated IMG structure where plot y axis intercept values 
      times -1 will be placed; enter NULL, if not needed. */
  IMG *ic_img,
  /** Pointer to initiated IMG structure where the number of plot data points
      actually used in the fit is written; enter NULL, when not needed. */
  IMG *nr_img,
  /** Pointer to a string (allocated for at least 64 chars) where error message
      or other execution status will be written; enter NULL, if not needed. */
  char *status,
  /** Verbose level; if zero, then nothing is printed to stderr or stdout. */
  int verbose
) {
  int zi, yi, xi, fi, nr, pn, ret=0;
  DFT tac;
  double *plotData, *xaxis, *yaxis, slope, ic, f;
  double aucrat;


  if(verbose>0) {
    printf("img_logan(%s, dyn_img, %d, %d, %g, %g, ki_img, ", __func__, start, end, thrs, k2);
    if(ic_img==NULL) printf("NULL, "); else printf("ic_img, ");
    if(nr_img==NULL) printf("NULL, "); else printf("nr_img, ");
    if(status==NULL) printf("NULL)\n"); else printf("status)\n");
  }
  /* Initial check for the arguments */
  if(status!=NULL) sprintf(status, "invalid data");
  if(dyn_img->status!=IMG_STATUS_OCCUPIED || dyn_img->dimt<1) return(1);
  if(input==NULL || input->frameNr<1) return(2);
  nr=1+end-start; if(nr<2) return(3);
  if(end>dyn_img->dimt-1 || start<0) return(4);
  if(vt_img==NULL) return(5);
  /* Convert input time units to min */
  if(input->timeunit==TUNIT_SEC) dftTimeunitConversion(input, TUNIT_MIN);

  /* Check that input contains samples until at least 80% of line fit range */
  if(verbose>1) {
    printf("input_last_sample_time := %g\n", input->x[input->frameNr-1]);
    printf("logan_start_time := %g\n", dyn_img->start[start]/60.0);
    printf("logan_end_time := %g\n", dyn_img->end[end]/60.0);
  }
  if(input->x[input->frameNr-1] < (0.2*dyn_img->mid[start]+0.8*dyn_img->mid[end])/60.0) {
    sprintf(status, "too few input samples"); return(6);
  }

  /* Allocate memory for interpolated and integrated input and 
     tissue pixel TACs */
  dftInit(&tac); if(dftSetmem(&tac, nr, 2)!=0) {sprintf(status, "out of memory"); return(11);}
  strcpy(tac.voi[0].voiname, "input"); strcpy(tac.voi[1].voiname, "tissue");
  tac.voiNr=2; tac.frameNr=nr;
  for(fi=0; fi<tac.frameNr; fi++) {
    tac.x1[fi]=dyn_img->start[start+fi]/60.;
    tac.x2[fi]=dyn_img->end[start+fi]/60.;
    tac.x[fi]=dyn_img->mid[start+fi]/60.;
  }

  /* If input sample times are much different from PET frame times,
     as they usually are if arterial plasma is used as input, then use 
     interpolate4pet(), otherwise with image-derived input use petintegral()
  */
  int equal_times;
  equal_times=check_times_dft_vs_img(dyn_img, input, verbose-1);
  if(equal_times==1) {
    if(verbose>1) printf("copying input curve and using petintegral()\n");
    /* Copy frame start and end times from image to input data */
    ret=copy_times_from_img_to_dft(dyn_img, input, verbose-1);
    if(ret) {
      sprintf(status, "cannot set input sample times");
      dftEmpty(&tac); return(12);
    }
    /* integral must be calculated from the zero time */
    ret=petintegral(input->x1, input->x2, input->voi[0].y, input->frameNr,
      input->voi[0].y2, input->voi[0].y3);
    /* because times are the same, the values can be copied directly */
    if(ret==0) for(fi=0; fi<tac.frameNr; fi++) {
      tac.voi[0].y[fi]=input->voi[0].y[start+fi];
      tac.voi[0].y2[fi]=input->voi[0].y2[start+fi];
      tac.voi[0].y3[fi]=input->voi[0].y3[start+fi];
    }
  } else {
    if(verbose>1) printf("using interpolate4pet() for input curve\n");
    ret=interpolate4pet(
      input->x, input->voi[0].y, input->frameNr,
      tac.x1, tac.x2, tac.voi[0].y, tac.voi[0].y2, tac.voi[0].y3, tac.frameNr
    );
  }
  if(ret) {
    dftEmpty(&tac);
    sprintf(status, "cannot interpolate input data"); return(12);
  }
  if(verbose>3) dftPrint(&tac);

  /*
   *  Allocate result images and fill the header info
   */
  /* Vt image */
  imgEmpty(vt_img);
  ret=imgAllocateWithHeader(vt_img, dyn_img->dimz, dyn_img->dimy, dyn_img->dimx, 1, dyn_img);
  if(ret) {
    sprintf(status, "cannot setup memory for Vt image");
    imgEmpty(vt_img); dftEmpty(&tac); return(21);
  }
  vt_img->unit=CUNIT_ML_PER_ML;
  vt_img->decayCorrection=IMG_DC_NONCORRECTED; vt_img->isWeight=0;
  vt_img->start[0]=dyn_img->start[start]; vt_img->end[0]=dyn_img->end[end];
  if(verbose>9) imgInfo(vt_img);
  /* Ic image */
  if(ic_img!=NULL) {
    imgEmpty(ic_img);
    ret=imgAllocateWithHeader(ic_img, dyn_img->dimz, dyn_img->dimy, dyn_img->dimx, 1, dyn_img);
    if(ret) {
      sprintf(status, "cannot setup memory for Ic image");
      imgEmpty(ic_img); imgEmpty(vt_img); dftEmpty(&tac); return(23);
    }
    ic_img->unit=CUNIT_UNITLESS; /* mL/mL */
    ic_img->decayCorrection=IMG_DC_NONCORRECTED; ic_img->isWeight=0;
    ic_img->start[0]=dyn_img->start[start]; ic_img->end[0]=dyn_img->end[end];
    if(verbose>9) imgInfo(ic_img);
  }
  /* Nr image */
  if(nr_img!=NULL) {
    imgEmpty(nr_img);
    ret=imgAllocateWithHeader(nr_img, dyn_img->dimz, dyn_img->dimy, dyn_img->dimx, 1, dyn_img);
    if(ret) {
      sprintf(status, "cannot setup memory for nr image");
      imgEmpty(nr_img); imgEmpty(ic_img); imgEmpty(vt_img); dftEmpty(&tac);
      return(25);
    }
    nr_img->unit=CUNIT_UNITLESS;
    nr_img->decayCorrection=IMG_DC_NONCORRECTED; nr_img->isWeight=0;
    nr_img->start[0]=dyn_img->start[start]; nr_img->end[0]=dyn_img->end[end];
    if(verbose>9) imgInfo(ic_img);
  }

  /* Allocate memory for graphical analysis plot data */
  plotData=malloc(2*nr*sizeof(double));
  if(plotData==NULL) {
    sprintf(status, "cannot allocate memory for plots");
    imgEmpty(ic_img); imgEmpty(vt_img); dftEmpty(&tac); 
    return(25);
  }
  xaxis=plotData; yaxis=plotData+nr;
  float pxlauc[dyn_img->dimt];

  /* Calculate threshold */
  thrs*=tac.voi[0].y2[tac.frameNr-1];
  if(verbose>1) printf("  threshold-AUC := %g\n", thrs);

  /*
   *  Compute pixel-by-pixel
   */
  if(verbose>1) printf("computing MTGA pixel-by-pixel\n");
  int best_nr;
  for(zi=0; zi<dyn_img->dimz; zi++) {
    for(yi=0; yi<dyn_img->dimy; yi++) {
      for(xi=0; xi<dyn_img->dimx; xi++) {
        /* Initiate pixel output values */
        vt_img->m[zi][yi][xi][0]=0.0;
        if(ic_img!=NULL) ic_img->m[zi][yi][xi][0]=0.0;
        if(nr_img!=NULL) nr_img->m[zi][yi][xi][0]=0.0;
        /* Copy TTAC values */
        for(fi=0; fi<tac.frameNr; fi++) tac.voi[1].y[fi]=dyn_img->m[zi][yi][xi][start+fi];
        /* Compute TTAC AUC(0-t) and check for threshold */
        ret=fpetintegral(dyn_img->start, dyn_img->end, dyn_img->m[zi][yi][xi], 
                         dyn_img->dimt, pxlauc, NULL);
        if(ret) continue;
        if((pxlauc[dyn_img->dimt-1]/60.0) < thrs) continue;
        for(fi=0; fi<tac.frameNr; fi++)
          tac.voi[1].y2[fi]=pxlauc[start+fi]/60.0; // conc*sec -> conc*min
        /* Calculate Logan plot data */
        pn=logan_data(nr, tac.voi[0].y, tac.voi[0].y2, tac.voi[1].y, tac.voi[1].y2, k2, xaxis, yaxis);
        /* Line fit */
        if(fit_range==PRESET || pn<MTGA_BEST_MIN_NR) {
          ret=llsqperp(xaxis, yaxis, pn, &slope, &ic, &f);
          best_nr=pn;
        } else {
          ret=mtga_best_perp(xaxis, yaxis, pn, &slope, &ic, &f, &best_nr);
        }
        if(ret!=0) continue; // line fit failed
        /* Use 10xAUCratio as upper limit to prevent image where only 
           noise-induced hot spots can be seen */
        aucrat=tac.voi[1].y2[tac.frameNr-1]/tac.voi[0].y2[tac.frameNr-1];
        if(slope>10.0*aucrat) {
          if(verbose>50) printf("%g > 10 x %g\n", slope, aucrat);
          slope=10.0*aucrat;
        }
        /* copy result to pixels */
        vt_img->m[zi][yi][xi][0]=slope;
        if(ic_img!=NULL) ic_img->m[zi][yi][xi][0]=-ic;
        if(nr_img!=NULL) nr_img->m[zi][yi][xi][0]=best_nr;
      } /* next column */
    } /* next row */
  } /* next plane */

  free(plotData); dftEmpty(&tac);
  return(0);
}
/*****************************************************************************/

/*****************************************************************************/
