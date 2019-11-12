/// @file img_k1.c
/// @brief Functions for computing K1 pixel-by-pixel.
/// @author Vesa Oikonen
/// @todo Needs to be rewritten, including processing of units.
/*****************************************************************************/

/*****************************************************************************/
#include "libtpcmodext.h"
/*****************************************************************************/

/*****************************************************************************/
/** Computing pixel-by-pixel the K1 for irreversible PET tracers
    using previously determined Ki (K1*k3/(k2+k3) and bilinear regression.
\return Returns 0 if succesful, and >1 in case of an error.
 */
int img_k1_using_ki(
  /** Pointer to the TAC data to be used as model input. Sample times in minutes.
      Curve is interpolated to PET frame times, if necessary */
  DFT *input,
  /** Pointer to dynamic PET image data.
      Image and input data must be in the same calibration units */
  IMG *dyn_img,
  /** Nr of frames that will be included in the fit [3-frame_nr] */ 
  int frame_nr,
  /** Pointer to previously calculated Ki image */
  IMG *ki_img,
  /** Pointer to initiated IMG structure where K1 values will be placed */
  IMG *k1_img,
  /** Pointer to initiated IMG structure where (k2+k3) values 
      will be placed; enter NULL, if not needed */
  IMG *k2k3_img,
  /** Pointer to a string (allocated for at least 64 chars) where error message
      or other execution status will be written; enter NULL, if not needed */     
  char *status,
  /** Verbose level; if zero, then nothing is printed to stderr or stdout */
  int verbose
) {
  int zi, yi, xi, fi, ret=0;
  int nnls_n=2, nnls_m, m;
  DFT tac;
  clock_t fitStart, fitFinish;


  if(verbose>0) {
    printf("%s(input, dyn_img, %d, ki_img, k1_img, ", __func__, frame_nr);
    if(k2k3_img==NULL) printf("NULL, "); else printf("k2k3_img, ");
    if(status==NULL) printf("NULL)\n"); else printf("status)\n");
  }
  /* Initial check for the arguments */
  if(status!=NULL) sprintf(status, "invalid data");
  if(dyn_img->status!=IMG_STATUS_OCCUPIED || dyn_img->dimt<1) return(1);
  if(ki_img->status!=IMG_STATUS_OCCUPIED || ki_img->dimt<1) return(1);
  if(input==NULL || input->frameNr<1) return(1);
  if(frame_nr>dyn_img->dimt) return(1);
  if(k1_img==NULL) return(1);
  if(frame_nr<3) {
    if(status!=NULL) sprintf(status, "invalid fit time"); 
    return(1);
  }
  
  /* Check that input contains samples until at least 80% of line fit range */
  if(verbose>1) {
    printf("input_last_sample_time := %g\n", input->x[input->frameNr-1]);
    printf("k1_start_time := %g\n", dyn_img->start[0]/60.0);
    printf("k1_end_time := %g\n", dyn_img->end[frame_nr-1]/60.0);
  }
  if(input->x[input->frameNr-1] <
      (0.2*dyn_img->mid[0]+0.8*dyn_img->mid[frame_nr-1])/60.0)
  {
    sprintf(status, "too few input samples"); return(2);
  }

  fitStart=clock();

  /* Allocate memory for interpolated and integrated input and 
     tissue pixel TACs */
  dftInit(&tac); if(dftSetmem(&tac, frame_nr, 2)!=0) {
    sprintf(status, "out of memory"); return(3);
  }
  strcpy(tac.voi[0].voiname, "input"); strcpy(tac.voi[1].voiname, "tissue");
  tac.voiNr=2; tac.frameNr=frame_nr;
  for(fi=0; fi<tac.frameNr; fi++) {
    tac.x1[fi]=dyn_img->start[fi]/60.;
    tac.x2[fi]=dyn_img->end[fi]/60.;
    tac.x[fi]=dyn_img->mid[fi]/60.;
  }

  /* If input sample times are much different from PET frame times,
     as they usually are if arterial plasma is used as input, then use 
     interpolate4pet(), otherwise with image-derived input use petintegral()
  */
  ret=0;
  if(input->frameNr<dyn_img->dimt) ret=1;
  for(fi=0; fi<tac.frameNr; fi++) {
    if(verbose>8) printf("  %g  vs   %g\n", input->x1[fi], tac.x1[fi]);
    if(input->x1[fi]>tac.x1[fi]+0.034 || input->x1[fi]<tac.x1[fi]-0.034) {
      ret++; continue;}
    if(input->x2[fi]>tac.x2[fi]+0.034 || input->x2[fi]<tac.x2[fi]-0.034) ret++;
  }
  if(ret>0) {
    if(verbose>1) printf("using interpolate4pet() for input curve\n");
    ret=interpolate4pet(
      input->x, input->voi[0].y, input->frameNr,
      tac.x1, tac.x2, tac.voi[0].y, tac.voi[0].y2, tac.voi[0].y3, tac.frameNr
    );
  } else {
    if(verbose>1) printf("copying input curve and using petintegral()\n");
    /* because times are the same, the values can be copied directly */
    for(fi=0; fi<tac.frameNr; fi++) tac.voi[0].y[fi]=input->voi[0].y[fi];
    ret=petintegral(tac.x1, tac.x2, tac.voi[0].y, tac.frameNr,
      tac.voi[0].y2, tac.voi[0].y3);
  }
  if(ret) {
    dftEmpty(&tac); sprintf(status, "cannot interpolate input data");
    if(verbose>0) printf("  ret := %d\n", ret);
    return(4);
  }
  if(verbose>3) dftPrint(&tac);

  /*
   *  Allocate result images and fill the header info
   */
  /* K1 image */
  imgEmpty(k1_img);
  ret=imgAllocate(k1_img, dyn_img->dimz, dyn_img->dimy, dyn_img->dimx, 1);
  if(ret) {
    sprintf(status, "cannot allocate memory for K1 image");
    imgEmpty(k1_img); dftEmpty(&tac); return(5);
  }
  ret=imgCopyhdr(dyn_img, k1_img);
  if(ret) {
    sprintf(status, "cannot copy header info for result image");
    imgEmpty(k1_img); dftEmpty(&tac); return(5);
  }
  k1_img->unit=CUNIT_ML_PER_ML_PER_MIN;
  k1_img->decayCorrection=IMG_DC_NONCORRECTED; k1_img->isWeight=0;
  k1_img->start[0]=dyn_img->start[0]; k1_img->end[0]=dyn_img->end[frame_nr-1];
  if(verbose>9) imgInfo(k1_img);

  /* (k2+k3) image */
  if(k2k3_img!=NULL) {
    imgEmpty(k2k3_img);
    ret=imgAllocate(k2k3_img, dyn_img->dimz, dyn_img->dimy, dyn_img->dimx, 1);
    if(ret) {
      sprintf(status, "cannot allocate memory for (k2+k3) image");
      imgEmpty(k2k3_img); imgEmpty(k1_img); dftEmpty(&tac); return(6);
    }
    ret=imgCopyhdr(dyn_img, k2k3_img);
    if(ret) {
      sprintf(status, "cannot copy header info for result image");
      imgEmpty(k2k3_img); imgEmpty(k1_img); dftEmpty(&tac); return(6);
    }
    k2k3_img->unit=CUNIT_PER_MIN;
    k2k3_img->decayCorrection=IMG_DC_NONCORRECTED; k2k3_img->isWeight=0;
    k2k3_img->start[0]=dyn_img->start[0];
    k2k3_img->end[0]=dyn_img->end[frame_nr-1];
    if(verbose>9) imgInfo(k2k3_img);
  }

  /*
   *  Allocate memory required by NNLS; C99 !!!
   */
  if(verbose>1) fprintf(stdout, "allocating memory for NNLS\n");
  nnls_m=frame_nr;
  double *nnls_a[nnls_n], nnls_mat[nnls_n][nnls_m], nnls_b[nnls_m],
          nnls_zz[nnls_m];
  double nnls_x[nnls_n], nnls_wp[nnls_n], nnls_rnorm;
  int nnls_index[nnls_n];
  nnls_a[0]=nnls_mat[0];
  nnls_a[1]=nnls_mat[1];

  /*
   *  Compute pixel-by-pixel
   */
  if(verbose>1) printf("computing K1 pixel-by-pixel\n");
  for(zi=0; zi<dyn_img->dimz; zi++) {
    for(yi=0; yi<dyn_img->dimy; yi++) {
      for(xi=0; xi<dyn_img->dimx; xi++) {
        /* Initiate pixel output values */
        k1_img->m[zi][yi][xi][0]=0.0;
        if(k2k3_img!=NULL) k2k3_img->m[zi][yi][xi][0]=0.0;

        /* Copy and integrate pixel curve */
        for(m=0; m<nnls_m; m++) tac.voi[1].y[m]=dyn_img->m[zi][yi][xi][m];
        ret=petintegral(tac.x1, tac.x2, tac.voi[1].y, tac.frameNr,
                        tac.voi[1].y2, NULL);
        if(ret) continue;
        /* if AUC at the end is <= 0, then do nothing */
        if(tac.voi[1].y2[nnls_m-1]<=0.0) continue;       

        /* Fill the NNLS data matrix */
        for(m=0; m<nnls_m; m++) {
          nnls_a[0][m]=tac.voi[0].y2[m];
          nnls_a[1][m]=ki_img->m[zi][yi][xi][0]*tac.voi[0].y3[m]-tac.voi[1].y2[m];
          nnls_b[m]=tac.voi[1].y[m];
        }
        
        /* NNLS */
        ret=nnls(nnls_a, nnls_m, nnls_n, nnls_b, nnls_x, &nnls_rnorm,
                 nnls_wp, nnls_zz, nnls_index);
        if(ret>1) { /* no solution is possible */
          continue;
        } else if(ret==1) { /* max iteration count exceeded */ }
        k1_img->m[zi][yi][xi][0]=nnls_x[0];
        if(k2k3_img!=NULL) k2k3_img->m[zi][yi][xi][0]=nnls_x[1];
        
      } /* next column */
    } /* next row */
  } /* next plane */
  dftEmpty(&tac);

  fitFinish=clock(); //printf("CLOCKS_PER_SEC=%ld\n", CLOCKS_PER_SEC);
  //printf("%ld - %ld\n", fitFinish, fitStart);
  if(verbose>0 && fitFinish!=(clock_t)(-1)) printf("done in %g seconds.\n",
    (double)(fitFinish - fitStart) / (double)CLOCKS_PER_SEC );
  
  return(0);
}
/*****************************************************************************/

/*****************************************************************************/

