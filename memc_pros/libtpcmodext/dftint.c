/// line 202
/// @file dftint.c
/// @brief Functions for interpolating and integrating DFT data.
/// @author Vesa Oikonen
///
/*****************************************************************************/

/*****************************************************************************/
#include "libtpcmodext.h"
/*****************************************************************************/

/*****************************************************************************/
/** Verify that data to-be-interpolated does not need too much extrapolation in the beginning.
   @sa dftInterpolateCheckEnd, dftInterpolate, dftInterpolateInto, dftTimeIntegral
   @return Returns 0 if data is fine, 1 if it starts late but extrapolation can be
    done reliably, and -1 if extrapolation in the beginning would be too risky.
 */
int dftInterpolateCheckStart(
  /** Data that will be verified for reliable interpolation */
  DFT *input,
  /** Data containing sample times for interpolated data */
  DFT *output,
  /** Pointer to a string (allocated for at least 64 chars) where error message
      or other execution status will be written; enter NULL, if not needed */
  char *status,
  /** Verbose level; if zero, then nothing is printed to stderr or stdout */
  int verbose
) {
  int ri, itunit;
  double t1, t2, f;

  if(verbose>0) printf("dftInterpolateCheckStart()\n");

  // Check the input
  if(status!=NULL) sprintf(status, "program error");
  if(input==NULL || output==NULL) return -1;
  if(input->frameNr<1 || output->frameNr<1) return -1;
  if(input->frameNr==1 && output->frameNr>1) {
    if(status!=NULL) sprintf(status, "too short data for interpolation");
    return -1;
  }
  if(status!=NULL) sprintf(status, "ok");

  /* Try to make sure that time units are the same */
  dftMatchTimeunits(output, input, &itunit, verbose);

  /* Check the start */
  if(input->timetype==DFT_TIME_STARTEND && output->timetype==DFT_TIME_STARTEND) {
    t1=input->x1[0]; t2=output->x1[0];
  } else {
    t1=input->x[0];  t2=output->x[0];
  }
  if(0.95*t1>t2) {
    if(verbose>1) printf("t1 := %g\nt2 := %g\n", t1, t2);
    /* Check that first value is relatively low, so that there will not be any
       error when the initial part is assumed to be a straight line from
       (0,0) to the first sample point */
    f=0.25*dft_kBqMax(input);
    for(ri=0; ri<input->voiNr; ri++) if(input->voi[ri].y[0] > f) {
      if(status!=NULL) strcpy(status, "data starts too late");
      dftTimeunitConversion(input, itunit); // back to original time units
      return -1;
    }
    if(status!=NULL) strcpy(status, "data starts late");
    dftTimeunitConversion(input, itunit); // back to original time units
    return 1;
  }
  dftTimeunitConversion(input, itunit); // back to original time units
  return 0;
}
/*****************************************************************************/

/*****************************************************************************/
/** Verify that data to-be-interpolated does not need too much extrapolation in the end, 
    and that samples are not too sparse.

    Time units in DFTs can be different, if specified. 

    @sa dftInterpolateCheckStart, dftInterpolateInto, dftInterpolate
    @return Returns 0 if data is fine, 1 if extrapolation can be done, but there may
    be too few samples, and -1 if extrapolation in the end is impossible.
 */
int dftInterpolateCheckEnd(
  /** Data that will be verified for reliable interpolation */
  DFT *input,
  /** Data containing sample times for interpolated data */
  DFT *output,
  /** Pointer to a string (allocated for at least 64 chars) where error message
      or other execution status will be written; enter NULL, if not needed */     
  char *status,
  /** Verbose level; if zero, then nothing is printed to stderr or stdout */
  int verbose
) {
  int fi, n, itunit;
  double t1, t2;

  if(verbose>0) printf("dftInterpolateCheckEnd()\n");

  // Check the input
  if(status!=NULL) sprintf(status, "program error");
  if(input==NULL || output==NULL) return -1;
  if(input->frameNr<1 || output->frameNr<1) return -1;
  if(input->frameNr==1 && output->frameNr>1) {
    if(status!=NULL) sprintf(status, "too short data for interpolation");
    return -1;
  }
  if(status!=NULL) sprintf(status, "ok");

  /* Try to make sure that time units are the same */
  dftMatchTimeunits(output, input, &itunit, verbose);

  /* Check that input data extends almost to the end of output range */
  if(input->timetype==DFT_TIME_STARTEND && output->timetype==DFT_TIME_STARTEND) {
    t1=input->x2[input->frameNr-1]; t2=output->x2[output->frameNr-1]; 
  } else {
    t1=input->x[input->frameNr-1];  t2=output->x[output->frameNr-1];
  }
  if(t1<0.95*t2) {
    if(status!=NULL) strcpy(status, "too short data for interpolation");
    if(verbose>1) printf("t1 := %g\nt2 := %g\n", t1, t2);
    dftTimeunitConversion(input, itunit); // back to original time units
    return -1;
  }

  /* Check that input sample frequency is not too low in the end */
  if(output->frameNr>3) {
    if(output->timetype==DFT_TIME_STARTEND) {
      t1=output->x1[output->frameNr-3]; 
      t2=output->x2[output->frameNr-1];
    } else { 
      t1=output->x[output->frameNr-3]; 
      t2=output->x[output->frameNr-1];
    }
    for(fi=n=0; fi<input->frameNr; fi++)
      if(input->x[fi]>=t1 && input->x[fi]<=t2) n++;
    if(n<1 || (n<2 && t2>input->x[input->frameNr-1])) {
      if(status!=NULL) strcpy(status, "too sparse sampling for interpolation");
      if(verbose>1) printf("n=%d t1=%g t2=%g\n", n, t1, t2);
      dftTimeunitConversion(input, itunit); // back to original time units
      return -1; // Error
    }
    if(n<2 || (n<3 && t2>input->x[input->frameNr-1])) {
      if(status!=NULL) strcpy(status, "too sparse sampling for interpolation");
      if(verbose>1) printf("n=%d t1=%g t2=%g\n", n, t1, t2);
      dftTimeunitConversion(input, itunit); // back to original time units
      return 1; // Warning
    }
  } 
  dftTimeunitConversion(input, itunit); // back to original time units
  return 0;
}
/*****************************************************************************/

/*****************************************************************************/
/** Interpolate (and integrate) TAC data to sample times that are given with another TAC data.
    PET frame lengths are taken into account if available in tissue DFT struct.

    Input frame lengths are taken into account if the framing is same as with tissue data.

   @sa dftInterpolateInto, dftTimeIntegral
   @return Returns 0 if successful, and <>0 in case of an error.
 */ 
int dftInterpolate(
  /** Data which is interpolated; make sure that time unit is the same as in
   *  tissue data; time range is checked to cover the tissue data */
  DFT *input,
  /** Data to which (sample times) the interpolation is done */
  DFT *tissue,
  /** Pointer to initiated DFT into which interpolated values and integrals
   *  will be written */
  DFT *output,
  /** Pointer to a string (allocated for at least 64 chars) where error message
      or other execution status will be written; enter NULL, if not needed */
  char *status,
  /** Verbose level; if zero, then nothing is printed to stderr or stdout */
  int verbose
) {
  int fi, ri, ret;

  if(verbose>0) printf("dftInterpolate()");

  // Check the input
  if(input==NULL || tissue==NULL || output==NULL) {
    if(status!=NULL) sprintf(status, "program error");
    return 1;
  }

  // If input and tissue data have the same frame times already, then
  // copy frame times and timetype into input
  if(tissue->timetype==DFT_TIME_STARTEND && 
     check_times_dft_vs_dft(tissue, input, verbose)==1 &&
     input->frameNr<=tissue->frameNr)
  {
    for(fi=0; fi<input->frameNr; fi++) {
      input->x[fi]=tissue->x[fi];
      input->x1[fi]=tissue->x1[fi]; input->x2[fi]=tissue->x2[fi];
    }
    input->timetype=tissue->timetype;
  }

  // Check that there is no need for excess extrapolation
  ret=dftInterpolateCheckEnd(input, tissue, status, verbose);
  // if(ret<0) return 2;           /////////////// force to interpolate and integrate
  ret=dftInterpolateCheckStart(input, tissue, status, verbose);
  // if(ret<0) return 2;

  // Delete any previous output data
  dftEmpty(output);

  // Allocate memory for interpolated data
  if(dftSetmem(output, tissue->frameNr, input->voiNr)) {
    if(status!=NULL) strcpy(status, "memory allocation error");
    return 3;
  }
  output->voiNr=input->voiNr; output->frameNr=tissue->frameNr;

  // Copy header information
  dftCopymainhdr(input, output);
  for(ri=0; ri<input->voiNr; ri++) dftCopyvoihdr(input, ri, output, ri);

  // Copy frame information
  output->isweight=tissue->isweight;
  for(fi=0; fi<tissue->frameNr; fi++) {
    output->x[fi]=tissue->x[fi];
    output->x1[fi]=tissue->x1[fi]; output->x2[fi]=tissue->x2[fi];
    output->w[fi]=tissue->w[fi];
  }
  output->timetype=tissue->timetype;

  // Check if input and tissue data do have the same frame times already
  if(check_times_dft_vs_dft(tissue, input, verbose)==1 &&
     input->frameNr>=tissue->frameNr)
  {
    // copy the values directly and integrate in place
    for(ri=0, ret=0; ri<output->voiNr && ret==0; ri++) {
      for(fi=0; fi<tissue->frameNr; fi++)
        output->voi[ri].y[fi]=input->voi[ri].y[fi];
      if(output->timetype==3)
        ret=petintegral(output->x1, output->x2, output->voi[ri].y,
               output->frameNr, output->voi[ri].y2, output->voi[ri].y3);
      else
        ret=interpolate(output->x, output->voi[ri].y, output->frameNr,
               output->x, output->voi[ri].y, output->voi[ri].y2,
               output->voi[ri].y3, output->frameNr);
    }
    if(ret) {
      if(status!=NULL) sprintf(status, "cannot interpolate (%d)", ret);
      dftEmpty(output); return 5;
    }
    return 0; // that's it then
  }

  // Interpolate and integrate input data to tissue sample times,
  // taking into account tissue frame lengths if available
  for(ri=0, ret=0; ri<output->voiNr && ret==0; ri++) {
    if(output->timetype==DFT_TIME_STARTEND)
      ret=interpolate4pet(input->x, input->voi[ri].y, input->frameNr,
                output->x1, output->x2, output->voi[ri].y, output->voi[ri].y2,
                output->voi[ri].y3, output->frameNr);
    else
      ret=interpolate(input->x, input->voi[ri].y, input->frameNr,
                output->x, output->voi[ri].y, output->voi[ri].y2,
                output->voi[ri].y3, output->frameNr);
  }
  if(ret) {
    if(status!=NULL) sprintf(status, "cannot interpolate (%d)", ret);
    dftEmpty(output); return 6;
  }

  return 0;
}
/*****************************************************************************/

/*****************************************************************************/
/** Interpolate (and integrate) TAC data to sample times that are given with another TAC data. 
    New TAC is written into existing TAC data.
    PET frame lengths are taken into account if available in tissue DFT struct.
    Input frame lengths are taken into account if the framing is same as with tissue data.
   @sa dftTimeIntegral, dftInterpolate, dftInterpolateCheckStart, dftInterpolateCheckEnd
   @return Returns 0 if successful, and <>0 in case of an error.
 */ 
int dftInterpolateInto(
  /** Data which is interpolated; make sure that time unit is the same as in
   *  tissue data; time range is checked to cover the tissue data */
  DFT *inp,
  /** Data into which the interpolation is done */
  DFT *tis,
  /** Pointer to a string (allocated for at least 64 chars) where error message
      or other execution status will be written; enter NULL, if not needed */     
  char *status,
  /** Verbose level; if zero, then nothing is printed to stderr or stdout */
  int verbose
) {
  int fi, ri, ret;

  if(verbose>0) printf("dftInterpolateInto()\n");

  // Check the input
  if(inp==NULL || tis==NULL) {
    if(status!=NULL) sprintf(status, "program error");
    return 1;
  }
  ret=dft_nr_of_NA(inp); if(ret>0) {
    if(status!=NULL) sprintf(status, "missing sample(s)");
    return 2;
  }

  // Check that there is no need for excess extrapolation
  ret=dftInterpolateCheckEnd(inp, tis, status, verbose);
  if(ret<0) return 3;
  ret=dftInterpolateCheckStart(inp, tis, status, verbose);
  if(ret<0) return 4;

  // Allocate memory for interpolated data
  if(dftAddmem(tis, inp->voiNr)) {
    if(status!=NULL) strcpy(status, "memory allocation error");
    return 5;
  }

  // Copy header information
  for(ri=0; ri<inp->voiNr; ri++)
    dftCopyvoihdr(inp, ri, tis, tis->voiNr+ri);

  // Check if input and tissue data do have the same frame times already
  if(check_times_dft_vs_dft(inp, tis, verbose)==1 &&
     inp->frameNr>=tis->frameNr)
  {
    // copy the values directly and integrate in place
    for(ri=0, ret=0; ri<inp->voiNr && ret==0; ri++) {
      for(fi=0; fi<tis->frameNr; fi++)
        tis->voi[tis->voiNr+ri].y[fi]=inp->voi[ri].y[fi];
      if(tis->timetype==DFT_TIME_STARTEND)
        ret=petintegral(tis->x1, tis->x2, tis->voi[tis->voiNr+ri].y,
               tis->frameNr, tis->voi[tis->voiNr+ri].y2, 
               tis->voi[tis->voiNr+ri].y3);
      else
        ret=interpolate(tis->x, tis->voi[tis->voiNr+ri].y, tis->frameNr,
               tis->x, tis->voi[tis->voiNr+ri].y, 
               tis->voi[tis->voiNr+ri].y2,
               tis->voi[tis->voiNr+ri].y3, tis->frameNr);
    }
    if(ret) {
      if(status!=NULL) sprintf(status, "cannot interpolate (%d)", ret);
      dftEmpty(tis); return 7;
    }
    tis->voiNr+=inp->voiNr;
    return 0; // that's it then
  }

  // Interpolate and integrate input data to tissue sample times,
  // taking into account tissue frame lengths if available
  for(ri=0, ret=0; ri<inp->voiNr && ret==0; ri++) {
    if(tis->timetype==DFT_TIME_STARTEND)
      ret=interpolate4pet(inp->x, inp->voi[ri].y, inp->frameNr,
                tis->x1, tis->x2, tis->voi[tis->voiNr+ri].y, 
                tis->voi[tis->voiNr+ri].y2,
                tis->voi[tis->voiNr+ri].y3, tis->frameNr);
    else
      ret=interpolate(inp->x, inp->voi[ri].y, inp->frameNr,
                tis->x, tis->voi[tis->voiNr+ri].y, 
                tis->voi[tis->voiNr+ri].y2,
                tis->voi[tis->voiNr+ri].y3, tis->frameNr);
  }
  if(ret) {
    if(status!=NULL) sprintf(status, "cannot interpolate (%d)", ret);
    return 9;
  }
  tis->voiNr+=inp->voiNr;

  return 0;
}
/*****************************************************************************/

/*****************************************************************************/
/** Integration of regional TAC data from time1 to time2, i.e. AUC(t1,t2).

    You may want to test the time range before applying this routine, because
    this function accepts relatively large extrapolation large. 

   @sa dftInterpolateInto, dftInterpolate, dftInterpolateCheckStart, dftInterpolateCheckEnd
   @return Returns 0 when call was successful, and >0 in case of an error.
*/
int dftTimeIntegral(
  /** Regional TAC data in DFT struct. Number of samples must be at least one.
   *  If only one sample, then the integration time range must match with 
   *  the possible frame start and times.
   *  Frames do not have to be continuous in time. */
  DFT *dft,
  /** Time where to start integration (same unit as in TACs) */
  double t1,
  /** Time to stop integration (same unit as in TACs); must be higher than t1,
   *  except that t1=t2 is acceptable when calc_mode=average */
  double t2,
  /** Pointer to initiated but empty AUC DFT data (output) */
  DFT *idft,
  /** Calculate integral or average: 0=integral, 1=average */
  int calc_mode,
  /** Pointer to a string (allocated for at least 128 chars) where error message
      or other execution status will be written; enter NULL, if not needed */     
  char *status,
  /** Verbose level; if zero, then nothing is printed to stderr or stdout */
  int verbose
) {
  int ri, fi, ret, f1, f2;
  double fdur, aucroi, t[2], auc[2];
  double accept_tdif=1.0; // in seconds

  if(verbose>0)
    printf("dftTimeIntegral(dft, %g, %g, idft, %d, status, %d)\n", t1, t2, calc_mode, verbose);

  /* Check the arguments */
  if(status!=NULL) sprintf(status, "program error");
  if(dft==NULL || idft==NULL || t1<0.0 || t2<0.0) return STATUS_FAULT;
  fdur=t2-t1; if(fdur<0.0) return STATUS_FAULT;
  if(fdur==0.0 && (calc_mode!=1 || dft->frameNr>1)) return STATUS_FAULT;
  if(dft->frameNr<1 || dft->voiNr<1) return STATUS_FAULT;

  /* Set time unit based things */
  if(dft->timeunit==TUNIT_MIN) accept_tdif/=60.0;

  /* Create the DFT struct for AUC values */
  ret=dftAllocateWithHeader(idft, 1, dft->voiNr, dft);
  if(ret!=0) {
    if(status!=NULL) sprintf(status, "cannot setup AUC data (%d)", ret);
    return STATUS_FAULT;
  }
  /* 'frame' start and end time are taken from integration time */
  idft->timetype=DFT_TIME_STARTEND; 
  if(calc_mode==1 && fdur==0.0) idft->timetype=DFT_TIME_MIDDLE;
  idft->_type=DFT_FORMAT_STANDARD;
  /* Initially integral is zero */
  for(ri=0; ri<idft->voiNr; ri++) idft->voi[ri].y[0]=0.0;
  //dftPrint(dft);

  /* Different paths if frame start and end times are taken into account or not */
  if(dft->timetype==DFT_TIME_STARTEND) {

    /* Check that time range matches with pet frame */
    if(dft->frameNr==1) { // static data (one frame only)
      /* check that specified times match with frame times, but
         accept 1 s differences */
      if(fabs(dft->x2[0]-t2)>accept_tdif || fabs(dft->x1[0]-t1)>accept_tdif) {
        if(status!=NULL) {
           strcpy(status, "for static data the integration time range must");
           strcat(status, " be exactly as long as the scan");
        }
        return STATUS_FAULT;
      }
      /* Set the times, in case there was a small difference */
      dft->x2[0]=t2; dft->x1[0]=t1;
    } else { // Dynamic data (more than one frame)
      /* First frame must start before 1/3 of the integration time */
      /* Last frame must end after 2/3 of integration time */
      if(dft->x1[0]>(0.66*t1+0.34*t2) ||
         dft->x2[dft->frameNr-1]<(0.34*t1+0.66*t2)) {
        if(status!=NULL)
          strcpy(status, "integration time range oversteps data range");
        return STATUS_FAULT;
      }
    }

    /* Get the first and last frame index that resides inside integration time*/
    if(dft->frameNr==1) {
      f1=f2=0;
    } else {
      for(fi=0, f1=f2=-1; fi<dft->frameNr; fi++) {
        if(f1<0) {if(dft->x1[fi]>=t1 && dft->x2[fi]<=t2) f1=f2=fi;}
        if(f1>=0) {if(t2>=dft->x2[fi]) f2=fi;}
      }
    }
    if(verbose>1) printf("f1=%d f2=%d\n", f1, f2);
    if(f1>=0 && f2>=0) {
      /* Integrate over the frames that are included in time range as a whole */
      fi=f1;
      for(ri=0; ri<dft->voiNr; ri++) {
        idft->voi[ri].y[0]+=(dft->x2[fi]-dft->x1[fi])*dft->voi[ri].y[fi];
      }
      for(fi=f1+1; fi<=f2; fi++) {
        for(ri=0; ri<dft->voiNr; ri++) {
          idft->voi[ri].y[0]+=(dft->x2[fi]-dft->x1[fi])*dft->voi[ri].y[fi];
        }
        /* Check whether frames are contiguous */
        if(dft->x1[fi]==dft->x2[fi-1]) continue;
        /* When not, calculate the area of an imaginary frame */
        double x, a;
        x=(dft->x1[fi]+dft->x2[fi-1])/2.0;
        for(ri=0; ri<dft->voiNr; ri++) {
          a=(dft->x1[fi]-dft->x2[fi-1])*
            (dft->voi[ri].y[fi]-(dft->voi[ri].y[fi]-dft->voi[ri].y[fi-1])
             *(dft->x2[fi]+dft->x1[fi]-2.0*x)
             /(dft->x2[fi]+dft->x1[fi]-dft->x2[fi-1]-dft->x1[fi-1]));
          idft->voi[ri].y[0]+=a;
        }
      }

      /* If necessary, add the partial integrals */
      if(dft->x1[f1]>t1) {
        t[0]=t1; t[1]=dft->x1[f1];
        if(verbose>5) printf("t[0]=%g t[1]=%g\n", t[0], t[1]); 
        for(ri=0; ri<dft->voiNr; ri++) {
          ret=interpolate(dft->x, dft->voi[ri].y, dft->frameNr, t, NULL, auc, NULL, 2);
          if(ret) aucroi=0.0; else aucroi=auc[1]-auc[0];
          idft->voi[ri].y[0]+=aucroi;
        }
      }
      if(t2>dft->x2[f2]) {
        t[0]=dft->x2[f2]; t[1]=t2;
        if(verbose>5) printf("t[0]=%g t[1]=%g\n", t[0], t[1]);
        for(ri=0; ri<dft->voiNr; ri++) {
          ret=interpolate(dft->x, dft->voi[ri].y, dft->frameNr, t, NULL, auc, NULL, 2);
          if(ret) aucroi=0.0; else aucroi=auc[1]-auc[0];
          idft->voi[ri].y[0]+=aucroi;
        }
      }

    } else { // no full frames inside integration range

      t[0]=t1; t[1]=t2; 
      for(ri=0; ri<dft->voiNr; ri++) {
        ret=interpolate(dft->x, dft->voi[ri].y, dft->frameNr, t, NULL, auc, NULL, 2);
        if(ret) aucroi=0.0; else aucroi=auc[1]-auc[0];
        idft->voi[ri].y[0]+=aucroi;
      }

    }

    /* Set output image time frame */
    idft->x2[0]=t2; idft->x1[0]=t1; idft->x[0]=0.5*(t1+t2);

  } else if(dft->timetype==DFT_TIME_MIDDLE) { // do not consider frame lengths

    /* If average calculation was required, and sample times match the required
       time range, then directly take the values; otherwise, calculate AUC */
    if(calc_mode==1 && dft->x[0]==t1 && dft->x[0]==t2 && dft->frameNr==1) {
      for(ri=0; ri<dft->voiNr; ri++)
        idft->voi[ri].y[0]=dft->voi[ri].y[0];
    } else {

      /* First sample time must be before 1/3 of the integration time */
      /* Last sample time must end after 2/3 of integration time */
      if(dft->x[0]>(0.66*t1+0.34*t2) ||
         dft->x[dft->frameNr-1]<(0.34*t1+0.66*t2)) {
        if(status!=NULL) sprintf(status, "integration time range oversteps data range");
        return STATUS_FAULT;
      }

      t[0]=t1; t[1]=t2; 
      if(verbose>5) printf("t[0]=%g t[1]=%g\n", t[0], t[1]); 
      for(ri=0; ri<dft->voiNr; ri++) {
        ret=interpolate(dft->x, dft->voi[ri].y, dft->frameNr, t, NULL, auc, NULL, 2);
        if(ret) idft->voi[ri].y[0]=0.0; else idft->voi[ri].y[0]=auc[1]-auc[0];
      }
    }
    /* Set output image time frame */
    idft->x2[0]=t2; idft->x1[0]=t1; idft->x[0]=0.5*(t1+t2);

  } else {
    if(status!=NULL)
      sprintf(status, "frame mid times or start and end times required");
    return STATUS_FAULT;
  }

  /* If required, then calculate average by dividing integral with time */
    /* Set units accordingly */
  if(calc_mode!=0) { // average
    if(fdur>0.0)
      for(ri=fi=0; ri<idft->voiNr; ri++) 
        idft->voi[ri].y[fi]/=fdur;
    if(status!=NULL) sprintf(status, "average TAC [%g,%g] calculated", t1, t2);
  } else { // integral
    int unit;
    unit=petCunitId(idft->unit);
    strcpy(idft->unit, dftUnit(CUNIT_UNKNOWN));
    if(unit==CUNIT_KBQ_PER_ML) {
      if(idft->timeunit==TUNIT_MIN)
        strcpy(idft->unit, dftUnit(CUNIT_MIN_KBQ_PER_ML));
      else if(idft->timeunit==TUNIT_SEC)
        strcpy(idft->unit, dftUnit(CUNIT_SEC_KBQ_PER_ML));
    }
    if(status!=NULL) sprintf(status, "TAC integral [%g,%g] calculated", t1, t2);
  }
  
  return(0);
}
/*****************************************************************************/

/*****************************************************************************/
/** Calculate simple derivatives from regional PET TACs. Old version.
 *  Requires that frame start and end times are known.
\return Returns 0 if successul.
 */
int dftDerivative_old(
  /** DFT struct containing the regional tissue TACs */
  DFT *dft,
  /** Allocated DFT struct for derivatives */
  DFT *deriv,
  /** Pointer to allocated string where the status or error message is written;
   *  set to NULL if not needed */
  char *status
) {
  int ri, fi;
  double fdur;

  /* check input */
  if(status!=NULL) sprintf(status, "invalid input for dftDerivative()");
  if(dft==NULL || dft->frameNr<1 || dft->voiNr<1) return 1;
  if(deriv==NULL || deriv->frameNr<dft->frameNr || deriv->voiNr<dft->voiNr)
    return 2;
  if(dft->timetype!=3) {
    if(status!=NULL) sprintf(status, "frame start and end times are required");
    return 3;
  }

  /* calculate derivative */
  for(fi=0; fi<dft->frameNr; fi++) {
    fdur=dft->x2[fi]-dft->x1[fi];
    if(fdur<=1.0E-10) { 
      for(ri=0; ri<dft->voiNr; ri++) deriv->voi[ri].y[fi]=0.0;
      continue;
    }
    for(ri=0; ri<dft->voiNr; ri++) {
      deriv->voi[ri].y[fi]=dft->voi[ri].y[fi];
      if(fi>0) deriv->voi[ri].y[fi]-=dft->voi[ri].y[fi-1];
      deriv->voi[ri].y[fi]/=fdur;
    }
  }
  return 0;
}
/*****************************************************************************/

/*****************************************************************************/
/** Calculate simple derivatives from regional PET TACs.
 *  This must not be used for any quantitative purpose.
\return Returns 0 if successul.
 */
int dftDerivative(
  /** DFT struct containing the regional tissue TACs */
  DFT *dft,
  /** Allocated DFT struct for derivatives */
  DFT *deriv,
  /** Pointer to allocated string where the status or error message is written;
   *  set to NULL if not needed */
  char *status
) {
  int ri, fi;
  double fdur;

  /* check input */
  if(status!=NULL) sprintf(status, "invalid input for dftDerivative()");
  if(dft==NULL || dft->frameNr<1 || dft->voiNr<1) return 1;
  if(deriv==NULL || deriv->frameNr<dft->frameNr || deriv->voiNr<dft->voiNr)
    return 2;
  if(dft->timetype!=DFT_TIME_MIDDLE && dft->timetype!=DFT_TIME_STARTEND) {
    if(status!=NULL)
      sprintf(status, "frame start and end times or mid times are required");
    return 3;
  }

  /* calculate frame mid times if necessary */
  if(dft->timetype==DFT_TIME_STARTEND)
    for(fi=0; fi<dft->frameNr; fi++)
      dft->x[fi]=0.5*(dft->x1[fi]+dft->x2[fi]);

  /* calculate derivative */
  if(dft->x[0]<=1.0E-20) for(ri=0; ri<dft->voiNr; ri++)
    deriv->voi[ri].y[0]=0.0;
  else for(ri=0; ri<dft->voiNr; ri++)
    deriv->voi[ri].y[0]=dft->voi[ri].y[0]/dft->x[0];
  for(fi=1; fi<dft->frameNr; fi++) {
    fdur=dft->x[fi]-dft->x[fi-1];
    if(fdur<=1.0E-20) for(ri=0; ri<dft->voiNr; ri++)
      deriv->voi[ri].y[fi]=0.0;
    else for(ri=0; ri<dft->voiNr; ri++)
      deriv->voi[ri].y[fi]=(dft->voi[ri].y[fi]-dft->voi[ri].y[fi-1])/fdur;
  }  
  for(ri=0; ri<dft->voiNr; ri++) for(fi=0; fi<dft->frameNr-1; fi++) {
    deriv->voi[ri].y[fi]+=deriv->voi[ri].y[fi+1];
    deriv->voi[ri].y[fi]*=0.5;
  }

  return 0;
}
/*****************************************************************************/

/*****************************************************************************/
