/// @file extrapolate.c
/// @brief Procedures for extrapolating PET TAC data.
/// @author Vesa Oikonen
///
/*****************************************************************************/

/*****************************************************************************/
#include "libtpcmodext.h"
/*****************************************************************************/

/*****************************************************************************/
/** Extrapolation of exponentially decreasing tail of PET radiotracer plasma
 *  curves. This is accomplished by fitting line to the end-part of the plot
 *  of the natural logarithm of tracer concentration against time.
 *  By default, the included data points are determined based on maximum
 *  adjusted R^2 from at least three points; fit to the larger 
 *  number of points is used in case difference in adjusted R^2 values is not
 *  larger than 0.0001. 
 *  This function is used and tested with program extrapol.
\return Returns 0 when successful, otherwise <>0.
 */
int extrapolate_monoexp(
  /** Pointer to original TAC data */
  DFT *dft,
  /** By default, the search for the best line fit is started from the last
   *  sample towards the first sample; set fit time to -1 to use this default.
   *  However, if the end phase is unreliable or very noisy, you may want to
   *  set fittime to include only certain time range from the beginning.
   *  Function will write here the fittime that was actually used. */
  double *fittime,
  /** The minimum number of samples used in searching the best fit; at least 2,
   *  but 3 is recommended. If data is very noisy, then this number may need
   *  to be increased.
   *  Function will write here the nr of samples that was actually used.
   *  This can be used as an alternative to mintime or in addition to it. */
  int *min_nr,
  /** The maximum number of samples used in searching the best fit; 
   *  must be higher than min_nr, or set to -1 to not to limit the number. */
  int max_nr,
  /** Minimum time range used in searching the best fit. If data is very noisy,
   *  then this may need to be set, otherwise setting mintime to -1 will use
   *  the default.
   *  This can be used as an alternative to min_nr or in addition to it. */
  double mintime,
  /** Last extrapolated sample time in same time units than in the data */
  double extr_to,
  /** Pointer to data for extrapolated TACs. Struct must be initiated.
   *  Any existing data is deleted. */
  DFT *ext,
  /** Give file pointer (for example stdout) where log information is printed;
   *  NULL if not needed */
  FILE *loginfo,
  /** Pointer to a string (allocated for at least 64 chars) where error message
      or other execution status will be written; enter NULL, if not needed */     
  char *status
) {
  int fi, ri, ret, n, to, from, to_min, from_min, orig_max_nr;
  int fitNr, snr;
  double kel, c0=0.0, *cx, *cy;
  double slope, slope_sd, ic, ic_sd, r, f, adj_r2, adj_r2_max;
  double extr_sampl;


  /* Check the input */
  if(status!=NULL) sprintf(status, "program error");
  if(dft==NULL || ext==NULL) return -1;
  if(*min_nr<2) *min_nr=3;  // Setting erroneous value to a recommended value
  if(max_nr>-1 && max_nr<*min_nr) return -2;
  orig_max_nr=max_nr;

  /* Set the last sample to be fitted */
  fitNr=dft->frameNr;
  if(*fittime>0.0) {
    while(dft->x[fitNr-1]>*fittime && fitNr>0) fitNr--;
    if(loginfo!=NULL) fprintf(loginfo, "fitNr := %d\nTime range := %g - %g\n", 
      fitNr, dft->x[0], dft->x[fitNr-1]);
  }
  *fittime=dft->x[fitNr-1];
  /* Check that there are at least 3 samples */
  if(fitNr<3) { /* min_nr is checked later */
    if(status!=NULL) sprintf(status, "too few samples for extrapolation");
    return(2);
  }
  
  /* If mintime was set, then get min_nr based on that */
  if(mintime>0.0) {
    fi=0; while(dft->x[fi]<(*fittime-mintime) && fi<fitNr) fi++;
    if(--fi<=0) {
      if(status!=NULL) sprintf(status, "required minimum fit range too large");
      return(2);
    }
    *min_nr=((fitNr-fi)>*min_nr?fitNr-fi:*min_nr);
  }
  if(loginfo!=NULL) fprintf(loginfo, "final_min_nr := %d\n", *min_nr);


  /*
   *  Initiate data for extrapolated data
   */
  /* Determine the sampling time for extrapolated range */
  extr_sampl=0.5*(dft->x[dft->frameNr-1]-dft->x[dft->frameNr-3]);
  if(loginfo!=NULL) fprintf(loginfo, "extr_sampl=%g\n", extr_sampl);
  if(extr_sampl<1.0E-8) {
    if(status!=NULL) sprintf(status, "check sample times");
    return(2);
  }
  /* Determine the nr of extrapolated samples */
  f=extr_to-dft->x[dft->frameNr-1];
  if(f<=0.0) {
    if(status!=NULL) sprintf(status, "no extrapolation is needed");
    return(2);
  }
  n=ceil(f/extr_sampl);
  if(loginfo!=NULL) fprintf(loginfo, "  extr_sampl=%g n=%d\n", extr_sampl, n);
  /* Allocate memory */
  dftEmpty(ext);
  ret=dftSetmem(ext, dft->frameNr+n, dft->voiNr);
  if(ret) {
    if(status!=NULL) sprintf(status, "error in memory allocation.\n");
    return(4);
  }
  ext->frameNr=dft->frameNr+n; ext->voiNr=dft->voiNr;
  /* Set sample times */
  for(fi=0; fi<dft->frameNr; fi++) {
    ext->x[fi]=dft->x[fi]; ext->x1[fi]=dft->x1[fi]; ext->x2[fi]=dft->x2[fi];}
  if(dft->timetype==DFT_TIME_MIDDLE) {
    for(; fi<ext->frameNr; fi++) {
      ext->x[fi]=ext->x[fi-1]+extr_sampl;
      ext->x1[fi]=ext->x2[fi-1]; ext->x2[fi]=ext->x[fi]+0.5*extr_sampl;
    }
  } else {
    for(; fi<ext->frameNr; fi++) {
      ext->x1[fi]=ext->x2[fi-1]; ext->x2[fi]=ext->x1[fi]+extr_sampl;
      ext->x[fi]=0.5*(ext->x1[fi]+ext->x2[fi]);
    }
  }
  /* Copy "header" information */
  dftCopymainhdr(dft, ext);
  for(ri=0; ri<dft->voiNr; ri++) dftCopyvoihdr(dft, ri, ext, ri);
  /* Copy existing values */
  for(ri=0; ri<dft->voiNr; ri++) for(fi=0; fi<dft->frameNr; fi++)
    ext->voi[ri].y[fi]=dft->voi[ri].y[fi];


  /*
   *  Make ln transformation for TACs
   */
  if(loginfo!=NULL) fprintf(loginfo, "ln transformation\n");
  for(ri=0; ri<dft->voiNr; ri++) for(fi=0; fi<dft->frameNr; fi++) {
    if(!isnan(dft->voi[ri].y[fi]) && dft->voi[ri].y[fi]>0.0)
      dft->voi[ri].y2[fi]=log(dft->voi[ri].y[fi]);
    else
      dft->voi[ri].y2[fi]=nan("");
  }


  /*
   *  Compute best linear fit to the end of ln transformed TACs
   */
  if(loginfo!=NULL) fprintf(loginfo, "linear fitting\n");
  cx=(double*)malloc(2*fitNr*sizeof(double)); if(cx==NULL) {
    if(status!=NULL) sprintf(status, "out of memory");
    return(6);
  }
  cy=cx+fitNr;
  for(ri=0; ri<dft->voiNr; ri++) {
    kel=0.0; max_nr=orig_max_nr;
    /* Print TAC name, if more than one was found */
    if(dft->voiNr>1 && loginfo!=NULL)
      fprintf(loginfo, "%s :\n", dft->voi[ri].name);

    /* Copy appropriate TAC data */
    for(fi=n=0; fi<fitNr; fi++)
      if(dft->x[fi]>0.0 && !isnan(dft->voi[ri].y2[fi])) {
        cx[n]=dft->x[fi]; cy[n++]=dft->voi[ri].y2[fi];}
    if(n<*min_nr) {
      if(status!=NULL) sprintf(status, "check the datafile (%d<%d)", n, *min_nr);
      free(cx); return(7);
    }
    if(max_nr<=0) max_nr=n;
    /* Search the plot range that gives the max adjusted R^2 */
    from_min=to_min=-1; adj_r2_max=-9.99E+99;
    for(from=n-max_nr, to=n-1; from<1+n-*min_nr; from++) {
      snr=(to-from)+1;
      /* Calculation of linear regression using pearson() */
      ret=pearson(
        cx+from, cy+from, snr, &slope, &slope_sd, &ic, &ic_sd, &r, &f
      );
      if(ret==0) {
        adj_r2= 1.0 - ((1.0-r*r)*(double)(snr-1)) / (double)(snr-2); 
      } else {
        adj_r2=-9.99E+99;
      }
      //if(cv<cv_min) {
      if(adj_r2>adj_r2_max+0.0001) {
        adj_r2_max=adj_r2; from_min=from; to_min=to; kel=-slope; c0=exp(ic);
      }
      if(loginfo!=NULL)
        fprintf(loginfo, "  adj_r2=%g from=%d (%g)\n", adj_r2, from, *(cx+from) );
    }
    if(from_min<0) {
      if(status!=NULL) sprintf(status, "check the datafile");
      free(cx); return(7);
    }
    /* Check the sign of slope */
    if(kel>=0.0) {
      /* negative slope i.e. positive kel is ok */
      if(loginfo!=NULL)
         fprintf(loginfo, "  k(el)=%g adj_r2=%g C(0)=%g; %d data points.\n",
           kel, adj_r2_max, c0, (to_min-from_min)+1 );
      /* Calculate the extrapolated values */
      for(fi=fitNr; fi<ext->frameNr; fi++) {
        ext->voi[ri].y[fi]=c0*exp(-kel*ext->x[fi]);
      }
    } else {
      /* If slope is positive, then extrapolate with line f(x)=b+0*t */
      for(fi=fitNr-1, f=0.0, n=0; fi>=0 && n<3; fi--)
        if(!isnan(dft->voi[ri].y[fi])) {f+=dft->voi[ri].y[fi]; n++;}
      if(n>0) f/=(double)n;
      if(status!=NULL)
        sprintf(status, "curve end is not descending; extrapolation with horizontal line determined as avg of %d samples", n);
      /* Calculate the extrapolated values */
      for(fi=fitNr; fi<ext->frameNr; fi++) {
        ext->voi[ri].y[fi]=f;
      }
    }
  } /* next curve */
  free(cx);

  if(status!=NULL) sprintf(status, "ok");
  return 0;
}
/*****************************************************************************/

/*****************************************************************************/
/** Interpolates TACs to automatically determined sample times with
    smaller intervals in the beginning.
    Only data in y arrays are interpolated; data in y2 and y3 are not used.
\return Function returns 0 when succesful, else a value >= 1.
*/
int dftAutointerpolate(
  /** Data to be interpolated is read from this struct */
  DFT *dft,
  /** Interpolated data is written in this struct; must be initiated;
   *  any previous content is deleted */
  DFT *dft2,
  /** The length of interpolated/extrapolated data */
  double endtime,
  /** Verbose level; if zero, then nothing is printed to stderr or stdout */
  int verbose
) {
  int i, newnr, maxNr=10000;
  double t, dt;


  if(verbose>0) printf("dftAutointerpolate(dft1, dft2, %g)\n", endtime);
  /* Check the arguments */
  if(dft->frameNr<1 || dft->voiNr<1) return 1;
  if(endtime<1.0 || endtime>1.0E12) return 2;
  if(dft->timetype!=DFT_TIME_STARTEND && dft->timetype!=DFT_TIME_MIDDLE)
    return 10;

  /* Calculate the number of interpolated data points */
  t=0.0; dt=0.02; newnr=1;
  while(t+dt<endtime && newnr<maxNr-1) {
    t+=dt; dt*=1.05; newnr++;
  }
  dt=endtime-t; t=endtime; newnr++;
  if(verbose>1) printf("newnr := %d\n", newnr);

  /* Allocate memory for interpolated data */
  dftEmpty(dft2); if(dftSetmem(dft2, newnr, dft->voiNr)) return 3;
  /* Copy header info */
  (void)dftCopymainhdr(dft, dft2); dft2->voiNr=dft->voiNr;
  for(i=0; i<dft->voiNr; i++) if(dftCopyvoihdr(dft, i, dft2, i)) return 4;

  /* Set times */
  if(dft->timetype==DFT_TIME_STARTEND) {

    dft2->timetype=DFT_TIME_STARTEND;
    t=0.0; dt=0.02; i=0;
    if(verbose>1) printf("%05d: %12.5f  %10.5f\n", i, t, dt);
    dft2->x1[i]=t; dft2->x2[i]=t+dt; dft2->x[i]=0.5*(dft2->x1[i]+dft2->x2[i]);
    i++;
    while((t+2.5*dt)<endtime && newnr<maxNr-2) {
      t+=dt; dt*=1.05;
      if(verbose>1) printf("%05d: %12.5f  %10.5f\n", i, t, dt);
      dft2->x1[i]=t; dft2->x2[i]=t+dt; dft2->x[i]=0.5*(dft2->x1[i]+dft2->x2[i]);
      i++;
    }
    t+=dt; dt=endtime-t;
    if(verbose>1) printf("%05d: %12.5f  %10.5f\n", i, t, dt);
    dft2->x1[i]=t; dft2->x2[i]=t+dt; dft2->x[i]=0.5*(dft2->x1[i]+dft2->x2[i]);
    i++; dft2->frameNr=i;

    /* Interpolate */
    for(i=0; i<dft->voiNr; i++)
      if(interpolate4pet(dft->x, dft->voi[i].y, dft->frameNr,
         dft2->x1, dft2->x2, dft2->voi[i].y, NULL, NULL, dft2->frameNr)) {
        dftEmpty(dft2); return 5;
      }

  } else if(dft->timetype==DFT_TIME_MIDDLE) {

    dft2->timetype=DFT_TIME_MIDDLE;
    t=0.0; dt=0.02; i=0;
    if(verbose>1) printf("%05d: %12.5f  %10.5f\n", i, t, dt);
    dft2->x1[i]=t; dft2->x2[i]=t+dt; dft2->x[i]=0.5*(dft2->x1[i]+dft2->x2[i]);
    i++;
    while((t+2.5*dt)<endtime && newnr<maxNr-2) {
      t+=dt; dt*=1.05;
      if(verbose>1) printf("%05d: %12.5f  %10.5f\n", i, t, dt);
      dft2->x1[i]=t; dft2->x2[i]=t+dt; dft2->x[i]=0.5*(dft2->x1[i]+dft2->x2[i]);
      i++;
    }
    t+=dt; dt=endtime-t;
    if(verbose>1) printf("%05d: %12.5f  %10.5f\n", i, t, dt);
    dft2->x1[i]=t; dft2->x2[i]=t+2.0*dt;
    dft2->x[i]=0.5*(dft2->x1[i]+dft2->x2[i]);
    i++; dft2->frameNr=i;

    /* Interpolate */
    for(i=0; i<dft->voiNr; i++)
      if(interpolate4pet(dft->x, dft->voi[i].y, dft->frameNr,
         dft2->x1, dft2->x2, dft2->voi[i].y, NULL, NULL, dft2->frameNr)) {
        dftEmpty(dft2); return 5;
      }

  }

  return 0;
}
/*****************************************************************************/

/*****************************************************************************/
/** Doubles the TAC sample number by making each sample/frame into two
 *  by linear interpolation.
 *  Only data in y arrays are interpolated; data in y2 and y3 are not used.
\return Function returns 0 when succesful, else a value >= 1.
*/
int dftDoubleFrames(
  /** Data to be interpolated is read from this struct */
  DFT *dft,
  /** Interpolated data is written in this struct; must be initiated;
   *  any previous content is deleted */
  DFT *dft2
) {
  int ri, fi, fj, ret;
  double f;


  /* Check the arguments */
  if(dft==NULL || dft2==NULL) return 1;
  if(dft->frameNr<1 || dft->voiNr<1) return 2;
  if(dft->timetype!=DFT_TIME_STARTEND && dft->x[0]<0.0) return 3;

  /* Allocate memory for interpolated data */
  dftEmpty(dft2); if(dftSetmem(dft2, 2*dft->frameNr, dft->voiNr)) return 11;
  /* Copy header info */
  (void)dftCopymainhdr(dft, dft2);
  dft2->voiNr=dft->voiNr; dft2->frameNr=2*dft->frameNr;
  for(ri=0; ri<dft->voiNr; ri++) if(dftCopyvoihdr(dft, ri, dft2, ri)) return 12;

  /* Set new time frames and interpolate */
  if(dft->timetype==DFT_TIME_STARTEND) {
    for(fi=fj=0; fi<dft->frameNr; fi++, fj+=2) {
      f=0.5*(dft->x1[fi]+dft->x2[fi]);
      dft2->x1[fj]=dft->x1[fi]; dft2->x2[fj]=f;
      dft2->x[fj]=0.5*(dft2->x1[fj]+dft2->x2[fj]);
      dft2->x1[fj+1]=f; dft2->x2[fj+1]=dft->x2[fi];
      dft2->x[fj+1]=0.5*(dft2->x1[fj+1]+dft2->x2[fj+1]);
      for(ri=0; ri<dft->voiNr; ri++)
        dft2->voi[ri].y[fj]=dft2->voi[ri].y[fj+1]=dft->voi[ri].y[fi];
    }
    ret=0;
  } else {
    for(fi=fj=0; fi<dft->frameNr; fi++) {
      if(dft->x[fi]<=0.0) {
        /* just copy, no doubles */
        dft2->x[fj++]=dft->x[fi]; continue;
      }
      if(fi==0) f=0.5*dft->x[fi];
      else f=0.5*(dft->x[fi-1]+dft->x[fi]);
      dft2->x[fj++]=f; dft2->x[fj++]=dft->x[fi];
    }
    dft2->frameNr=fj;
    for(ri=0, ret=0; ri<dft->voiNr; ri++) {
      ret=interpolate(dft->x, dft->voi[ri].y, dft->frameNr,
        dft2->x, dft2->voi[ri].y, NULL, NULL, dft2->frameNr);
      if(ret) break;
    }
  }
  if(ret) return 20+ret;
  return 0;
}
/*****************************************************************************/

/*****************************************************************************/
/** Interpolates TACs to automatically determined sample times with
    smaller intervals in the beginning.
    Only data in y arrays are interpolated; data in y2 and y3 are not used.
\return Function returns 0 when succesful, else a value >= 1.
*/
int dftDivideFrames(
  /** Data to be interpolated into more time frames is read from this struct */
  DFT *dft,
  /** Region index [0..voiNr-1] that is interpolated; <0 means all VOIs */
  int voi_index, 
  /** Nr of extra time frames that are created from each original frame;
   *  valid numers are 1-100; 1 doubles the frame number */
  int add_nr,
  /** Interpolated data is written in this struct; must be initiated,
   *  may be allocated but do not need to be; any previous content is deleted
   */
  DFT *dft2,
  /** Verbose level; if zero, then nothing is printed to stderr or stdout */
  int verbose
) {
  int ret, ri, fi, fj, i, new_frameNr=0, new_voiNr=0;
  double fdur;

  if(verbose>0) {
    printf("dftDivideFrames(*dft, %d, %d, *dft2)\n", voi_index, add_nr);
    fflush(stdout);
  }
  /* Check the arguments */
  if(dft==NULL || dft2==NULL) return 1;
  if(dft->frameNr<1 || dft->voiNr<1) return 2;
  if(add_nr<1 || add_nr>100) return 3;
  if(voi_index>=dft->voiNr) return 4;

  /* Calculate how many frames and TACs will be needed */
  if(dft->timetype==DFT_TIME_STARTEND) new_frameNr=(add_nr+1)*dft->frameNr;
  else new_frameNr=(add_nr+1)*dft->frameNr-1;
  if(voi_index<0) new_voiNr=dft->voiNr; else new_voiNr=1;
  if(verbose>1) {
    printf("new_frameNr := %d\n", new_frameNr);
    printf("new_voiNr := %d\n", new_voiNr);
    fflush(stdout);
  }

  /* Allocate memory for interpolated data if necessary */
  if(new_frameNr>dft2->frameNr || new_voiNr>dft2->_voidataNr) {
    if(verbose>1) {
      printf("deleting old data and allocating new\n"); fflush(stdout);}
    dftEmpty(dft2);
    if(dftSetmem(dft2, new_frameNr, dft->voiNr)) return 11;
  }

  /* Copy header info */
  ret=dftCopymainhdr(dft, dft2); if(ret!=0) return 12;
  dft2->voiNr=new_voiNr; dft2->frameNr=new_frameNr;
  if(voi_index>=0) ret=dftCopyvoihdr(dft, voi_index, dft2, 0);
  else {
    for(ri=0; ri<dft->voiNr; ri++) {
      ret=dftCopyvoihdr(dft, ri, dft2, ri); if(ret!=0) break;}
  }
  if(ret!=0) return 13;

  /* Set new frame times and interpolate */
  ret=0;
  if(dft->timetype==DFT_TIME_STARTEND) {
    for(fi=0, fj=0; fi<dft->frameNr; fi++) {
      fdur=(dft->x2[fi]-dft->x1[fi])/(double)(add_nr+1);
      for(i=0; i<add_nr+1; i++) {
        fj=fi*(add_nr+1)+i;
        dft2->x1[fj]=dft->x1[fi]+fdur*(double)i;
        dft2->x2[fj]=dft2->x1[fj]+fdur;
        dft2->x[fj]=0.5*(dft2->x1[fj]+dft2->x2[fj]);
        if(voi_index>=0)
          dft2->voi[0].y[fj]=dft->voi[voi_index].y[fi];
        else
          for(ri=0; ri<dft->voiNr; ri++) dft2->voi[ri].y[fj]=dft->voi[ri].y[fi];
      }
    }
    dft2->frameNr=fj+1;
  } else {
    dft2->x[0]=dft->x[0];
    for(fi=1, fj=0; fi<dft->frameNr; fi++) {
      fdur=(dft->x[fi]-dft->x[fi-1])/(double)(add_nr+1);
      for(i=0; i<add_nr+1; i++) {
        fj=(fi-1)*(add_nr+1) + i + 1;
        dft2->x[fj]=dft->x[fi-1]+fdur*(double)(i+1);
      }
    }
    dft2->frameNr=fj+1;
    if(voi_index>=0) {
      ret=interpolate(dft->x, dft->voi[voi_index].y, dft->frameNr,
        dft2->x, dft2->voi[0].y, NULL, NULL, dft2->frameNr);
    } else {
      for(ri=0; ri<dft->voiNr; ri++) {
        ret=interpolate(dft->x, dft->voi[ri].y, dft->frameNr,
                        dft2->x, dft2->voi[ri].y, NULL, NULL, dft2->frameNr);
        if(ret) break;
      }
    }
  }
  if(ret!=0) return 21;

  return 0;
}
/*****************************************************************************/

/*****************************************************************************/
/** Fits line to the end-part of TACs.
 *  By default, the included data points are determined based on maximum
 *  adjusted R^2 from at least three points; regression to the larger 
 *  number of points is used in case difference in adjusted R^2 values is not
 *  larger than 0.0001. 
\return Returns 0 when successful, otherwise <>0.
 */
int dft_end_line(
  /** Pointer to original TAC data */
  DFT *dft,
  /** By default, the search for the best line fit is started from the last
   *  sample towards the first sample; set fit time to -1 to use this default.
   *  However, if the end phase is unreliable or very noisy, you may want to
   *  set fittime to include only certain time range from the beginning.
   *  Function will write here the fittime that was actually used. */
  double *fittime,
  /** The minimum number of samples used in searching the best fit; at least 2,
   *  but 3 is recommended. If data is very noisy, then this number may need
   *  to be increased.
   *  Function will write here the nr of samples that was actually used.
   *  This can be used as an alternative to mintime or in addition to it. */
  int *min_nr,
  /** The maximum number of samples used in searching the best fit; 
   *  must be higher than min_nr, or set to -1 to not to limit the number. */
  int max_nr,
  /** Minimum time range used in searching the best fit. If data is very noisy,
   *  then this may need to be set, otherwise setting mintime to -1 will use
   *  the default.
   *  This can be used as an alternative to min_nr or in addition to it. */
  double mintime,
  /** Linear fitting can be applied to all data subsets in the fit time range
   *  (check_impr=0), or fitting is stopped when increasing n does not improve
   *  the adjusted R^2 (check_impr=1); the latter mode is for compatibility
   *  with WinNonlin. */ 
  int check_impr,
  /** Pointer to data for fitted parameters. Struct must be initiated.
   *  Any existing data is deleted. Additionally, adjusted R^2 is written
   *  as 3rd (non-documented) parameter. */
  FIT *fit,
  /** Give file pointer (for example stdout) where log information is printed;
   *  NULL if not needed */
  FILE *loginfo,
  /** Pointer to a string (allocated for at least 64 chars) where error message
      or other execution status will be written; enter NULL, if not needed */     
  char *status
) {
  int fi, ri, ret, n, to, from, to_min, from_min, orig_max_nr;
  int fitNr, snr;
  double *cx, *cy;
  double adj_r2, slope, slope_sd, ic, ic_sd, r, f, ySD_min;
  double adj_r2_max, adj_r2_prev=-10.0, ic_min, slope_min;


  /* Check the input */
  if(status!=NULL) sprintf(status, "program error");
  if(dft==NULL || fit==NULL) return -1;
  if(*min_nr<2) *min_nr=3;  // Setting erroneous value to a recommended value
  if(max_nr>-1 && max_nr<*min_nr) return -2;
  orig_max_nr=max_nr;

  /* Set the last sample to be fitted */
  fitNr=dft->frameNr;
  if(*fittime>0.0) {
    while(dft->x[fitNr-1]>*fittime && fitNr>0) fitNr--;
    if(loginfo!=NULL) fprintf(loginfo, "fitNr := %d\nTime range := %g - %g\n", 
      fitNr, dft->x[0], dft->x[fitNr-1]);
  }
  *fittime=dft->x[fitNr-1];
  /* Check that there are at least 3 samples */
  if(fitNr<3) { /* min_nr is checked later */
    if(status!=NULL) sprintf(status, "too few samples for linear fit");
    return(2);
  }
  
  /* If mintime was set, then get min_nr based on that */
  if(mintime>0.0) {
    fi=0; while(dft->x[fi]<(*fittime-mintime) && fi<fitNr) fi++;
    if(--fi<=0) {
      if(status!=NULL) sprintf(status, "required minimum fit range too large");
      return(2);
    }
    *min_nr=((fitNr-fi)>*min_nr?fitNr-fi:*min_nr);
  }
  if(loginfo!=NULL) fprintf(loginfo, "final_min_nr := %d\n", *min_nr);

  /* Allocate data for fits */
  if((ret=fit_allocate_with_dft(fit, dft))!=0) {
    if(loginfo!=NULL) fprintf(loginfo, "Error %d: cannot allocate memory for fits.\n", ret);
    if(status!=NULL) sprintf(status, "cannot allocate memory for fits");
    return 4;
  }
  /* and for x,y data to be fitted */  
  cx=(double*)malloc(2*fitNr*sizeof(double)); if(cx==NULL) {
    if(status!=NULL) sprintf(status, "out of memory");
    return(6);
  }
  cy=cx+fitNr;

  /* Fit each TAC */
  if(loginfo!=NULL) fprintf(loginfo, "linear fitting\n");
  for(ri=0; ri<dft->voiNr; ri++) {
    max_nr=orig_max_nr;
    /* Print TAC name, if more than one was found */
    if(dft->voiNr>1 && loginfo!=NULL)
      fprintf(loginfo, "%s :\n", dft->voi[ri].name);
    /* Set header */
    fit->voi[ri].parNr=2;
    fit->voi[ri].type=101;
    /* Copy appropriate TAC data */
    for(fi=n=0; fi<fitNr; fi++)
      if(dft->x[fi]>0.0 && !isnan(dft->voi[ri].y[fi])) {
        cx[n]=dft->x[fi]; cy[n++]=dft->voi[ri].y[fi];}
    if(n<*min_nr) {
      if(status!=NULL) sprintf(status, "check the datafile (%d<%d)", n, *min_nr);
      free(cx); return(7);
    }
    if(max_nr<=0) max_nr=n;
    /* Search the plot range that gives the max adjusted R^2 */
    from_min=to_min=-1; adj_r2_max=-9.99E+99; ic_min=slope_min=0.0; ySD_min=0.0;
    for(from=n-*min_nr, to=n-1; from>=n-max_nr; from--) {
      snr=(to-from)+1;
      /* Calculation of linear regression using pearson() */
      ret=pearson(
        cx+from, cy+from, snr, &slope, &slope_sd, &ic, &ic_sd, &r, &f
      );
      if(ret==0) {
        adj_r2= 1.0 - ((1.0-r*r)*(double)(snr-1)) / (double)(snr-2); 
        if(adj_r2<0.0 && loginfo!=NULL)
          fprintf(loginfo, "  r=%g; snr=%d\n", r, snr);
      } else {
        adj_r2=-9.99E+99;
      }
      if(adj_r2>adj_r2_max-0.0001) {
        adj_r2_max=adj_r2; from_min=from; to_min=to;
        ic_min=ic; slope_min=slope; ySD_min=f;
      }
      if(loginfo!=NULL)
        fprintf(loginfo, "  adj_r2=%g from=%d (%g)\n",
          adj_r2, from, *(cx+from) );
      /* check for improvement, if required */
      if(check_impr!=0 && adj_r2_prev>-1.0 && adj_r2>0.0 && adj_r2<adj_r2_prev)
        break;      

      adj_r2_prev=adj_r2;
    }
    if(from_min<0) {
      if(status!=NULL) sprintf(status, "check the datafile");
      free(cx); return(7);
    }
    if(loginfo!=NULL) fprintf(loginfo, "  adj_r2_max=%g.\n", adj_r2_max );
    /* Set fit line parameters */
    fit->voi[ri].p[0]=ic_min;
    fit->voi[ri].p[1]=slope_min;
    fit->voi[ri].p[2]=adj_r2_max;
    fit->voi[ri].wss=ySD_min;
    fit->voi[ri].start=cx[from_min]; fit->voi[ri].end=cx[to_min];
    fit->voi[ri].dataNr=(to_min-from_min)+1;
  } // next curve


  free(cx);
  if(status!=NULL) sprintf(status, "ok");
  return 0;
}
/******************************************************************************/

/******************************************************************************/
/** Natural logarithm (ln) transformation for TAC concentrations.
\return Returns 0 when successful, otherwise <>0.
 */
int dft_ln(
  /** Pointer to original TAC data */
  DFT *dft1,
  /** Pointer to allocated memory for ln transformed TAC data; enter NULL,
      if original data is to be overwritten by ln transformed values. */
  DFT *dft2
) {
  int ri, fi, ok_nr=0;
  DFT *out;

  if(dft1==NULL || dft1->voiNr<1 || dft1->frameNr<1) return 1;
  if(dft2!=NULL) out=dft2; else out=dft1;

  for(ri=0; ri<dft1->voiNr; ri++) for(fi=0; fi<dft1->frameNr; fi++) {
    if(!isnan(dft1->voi[ri].y[fi]) && dft1->voi[ri].y[fi]>0.0) {
      out->voi[ri].y[fi]=log(dft1->voi[ri].y[fi]); ok_nr++;
    } else {
      out->voi[ri].y[fi]=nan("");
    }
  }
  if(ok_nr>0) return 0; else return 2;
}
/******************************************************************************/

/******************************************************************************/

