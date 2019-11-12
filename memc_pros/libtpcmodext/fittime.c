/// @file fittime.c
/// @brief Check and set fit duration from data.
/// @author Vesa Oikonen
///
/*****************************************************************************/

/*****************************************************************************/
#include "libtpcmodext.h"
/*****************************************************************************/

/*****************************************************************************/
/** Reset user-defined fit time range to comply with DFT data.
   @sa fittime_from_img, getActualSamplenr, dftEndtime
   @return Returns the number of samples included in the fit range, or <0 in case of an error.
 */
int fittime_from_dft(
  /** Pointer to DFT containing (regional tissue) data;
      times can be in minutes or seconds, as long as units are defined. */
  DFT *dft,
  /** Pointer containing originally the requested fit start time (min).
      This is changed to contain the time of the first included frame.
      Unit must be minutes.
      Initially, set to <0 to start from the beginning of the data. */
  double *startTime,
  /** Pointer containing originally the requested fit end time (min).
      This is changed to contain the time of the last included frame.
      Unit must be minutes.
      Initially, set to <0 or to a very large value to reach to the end of data.
   */
  double *endTime,
  /** Function writes the index of the first included sample (frame) here. */
  int *first,
  /** Function writes the index of the last included sample (frame) here. */
  int *last,
  /** Verbose level; if zero, then nothing is printed to stderr or stdout. */
  int verbose
) {
  int fi, dataNr;
  
  if(verbose>0) printf("%s(*dft, %g, %g, first, last)\n", __func__, *startTime, *endTime);
  if(dft==NULL || startTime==NULL || endTime==NULL) return(-1);
  *first=*last=0;
  if(dft->frameNr<=0) {*startTime=*endTime=0.0; return 0;}
  /* Change start and end times to seconds if necessary */
  if(dft->timeunit==TUNIT_SEC) {*startTime*=60.; *endTime*=60.;}
  /* Check that data range is not outside required range */
  if(dft->x[dft->frameNr-1]<*startTime || dft->x[0]>*endTime) {
    *startTime=*endTime=0.0; return 0;}
  /* Get first and last data point */
  for(fi=0, *first=dft->frameNr; fi<dft->frameNr; fi++)
    if(dft->x[fi]>=*startTime) {*first=fi; break;}
  for(*last=fi=*first; fi<dft->frameNr; fi++)
    if(dft->x[fi]<=*endTime) *last=fi; else break;
  if(*first>=dft->frameNr) {*startTime=*endTime=0.0; return 0;}
  /* Correct fit range to frame start and end times */
  *startTime=(dft->timetype==DFT_TIME_STARTEND?dft->x1[*first]:dft->x[*first]);
  *endTime=(dft->timetype==DFT_TIME_STARTEND?dft->x2[*last]:dft->x[*last]);
  /* Calculate the number of data points in the fit range */
  dataNr=(*last)-(*first)+1;
  /* Change start and end times back to minutes if necessary */
  if(dft->timeunit==TUNIT_SEC) {*startTime/=60.; *endTime/=60.;}
  return dataNr;
}
/*****************************************************************************/

/*****************************************************************************/
/** Get the IMG frame end time of the last frame that is inside (mid time before) 
    the specified maximum fittime.
   @sa fittime_from_dft, imgEndtime
   @return Returns the number of IMG frames included in the fittime, or <0 in case of an error.
 */
int fittime_from_img(
  /** Pointer to IMG */
  IMG *img,
  /** Pointer containing originally the fit time maximum, after this the last included IMG frame 
      end time. Unit must be seconds.
      Initially, set to <0 or to a very large value to include all frames. */
  double *fittime,
  /** Verbose level; if zero, then nothing is printed to stderr or stdout. */
  int verbose
) {
  int fi, fitdimt;
  
  if(verbose>0) printf("%s(*img, %g)\n", __func__, *fittime);
  if(img==NULL || fittime==NULL) return(-1);
  if(img->dimt<=0) {*fittime=0.0; return 0;}
  if(*fittime<0) {*fittime=img->end[img->dimt-1]; return img->dimt;}
  for(fi=0, fitdimt=0; fi<img->dimt; fi++)
    if(img->mid[fi]>*fittime) break; else fitdimt++;
  *fittime=img->end[fitdimt-1];
  if(verbose>1) {
    printf("  fitdimt := %d\n", fitdimt);
    printf("  fittime := %g\n", *fittime);
  }
  return fitdimt;
}
/*****************************************************************************/

/*****************************************************************************/
/** Check whether DFT sample times are the same (or very close to) as the frame times in IMG. 
    This would suggest that DFT data originates from the same or similar PET scan. 
    Specified data sets are not edited.
    If frame nr is different, then only the commong frames are compared.
   @sa check_times_dft_vs_dft, copy_times_from_img_to_dft, dftEndtime, imgEndtime
   @return Returns 0 in case of no match, 1 if times do match, and <1 in case of an error.
 */
int check_times_dft_vs_img(
  /** Pointer to IMG data; times must be in sec as usual. */
  IMG *img,
  /** Pointer to DFT data; times can be both seconds or minutes, if unit is correctly set; 
      sample nr may be different than in IMG. */
  DFT *dft,
  /** Verbose level; if zero, then nothing is printed to stderr or stdout. */
  int verbose
) {
  int fi, n, smaller_dimt=0;
  double f, ts;
  double accepted_timedif=2.2; // sec

  if(verbose>0) printf("%s(*img, *dft)\n", __func__);
  if(img==NULL || dft==NULL) return -1;
  if(img->dimt<=0 || dft->frameNr<=0) return 0;
  
  /* Get the smaller frame nr */
  if(img->dimt<dft->frameNr) smaller_dimt=img->dimt; else smaller_dimt=dft->frameNr;

  /* With short study the accepted time difference must be shorter */
  f=0.01*img->end[img->dimt-1]; if(accepted_timedif>f) accepted_timedif=f;
  if(verbose>1) printf("accepted_timedif := %g [s]\n", accepted_timedif);

  /* Convert DFT times to sec if necessary */
  if(dft->timeunit==TUNIT_MIN) ts=60.0; else ts=1.0;

  /* Compare sample times frame-by-frame */
  for(fi=0, n=0; fi<smaller_dimt; fi++) {
    if(dft->timetype==DFT_TIME_MIDDLE) { // check frame mid times
      f=fabs(img->mid[fi]-dft->x[fi]*ts);
      if(verbose>10) printf("timedif[%d] := %g\n", fi, f);
      if(f>accepted_timedif) n++;
      continue;
    }
    if(dft->timetype==DFT_TIME_START || dft->timetype==DFT_TIME_STARTEND) {
      f=fabs(img->start[fi]-dft->x1[fi]*ts);
      if(verbose>10) printf("timedif[%d] := %g\n", fi, f);
      if(f>accepted_timedif) {n++; continue;}
    }
    if(dft->timetype==DFT_TIME_END || dft->timetype==DFT_TIME_STARTEND) {
      f=fabs(img->end[fi]-dft->x2[fi]*ts);
      if(verbose>10) printf("timedif[%d] := %g\n", fi, f);
      if(f>accepted_timedif) n++;
    }
  }
  if(verbose>2) printf("nr of different frame times := %d\n", n);

  if(n==0) return 1; else return 0; 
} 
/*****************************************************************************/

/*****************************************************************************/
/** Check whether sample times are the same (or very close to) in two DFT structs. 
    Data sets are not edited. If DFT structs contain different sample number, then only common nr
    of samples are compared.
   @sa check_times_dft_vs_img, copy_times_from_img_to_dft
   @return Returns 0 in case of no match, 1 if times do match, and <1 in case of an error.
 */
int check_times_dft_vs_dft(
  /** Pointer to first DFT data; times can be both seconds or minutes, if unit is correctly set. */
  DFT *dft1,
  /** Pointer to second DFT data; times can be both seconds or minutes, if unit is correctly set; 
      sample nr may be different than in dft1. */
  DFT *dft2,
  /** Verbose level; if zero, then nothing is printed to stderr or stdout. */
  int verbose
) {
  int fi, n, smaller_frameNr=0;
  double f, ts1, ts2;
  double accepted_timedif=2.2; // sec

  if(verbose>0) printf("%s(*img, *dft)\n", __func__);
  if(dft1==NULL || dft2==NULL) return -1;
  if(dft1->frameNr<=0 || dft2->frameNr<=0) return 0;
  
  /* Which has less samples? */
  if(dft1->frameNr<dft2->frameNr) smaller_frameNr=dft1->frameNr;
  else smaller_frameNr=dft2->frameNr;

  /* Convert times to sec if necessary (and possible) */
  ts1=ts2=1.0;
  if(dft1->timeunit!=TUNIT_UNKNOWN && dft1->timeunit!=TUNIT_UNKNOWN) {
    if(dft1->timeunit==TUNIT_MIN) ts1=60.0;
    if(dft2->timeunit==TUNIT_MIN) ts2=60.0;
  }
  if(verbose>1) {
    printf("dft1->timetype := %d\n", dft1->timetype);
    printf("dft2->timetype := %d\n", dft2->timetype);
    if(verbose>2) {
      printf("time range 1 := %g - %g %s\n", dft1->x[0],
             dft1->x[dft1->frameNr-1], petTunit(dft1->timeunit));
      printf("time range 2 := %g - %g %s\n", dft2->x[0],
             dft2->x[dft2->frameNr-1], petTunit(dft2->timeunit));
    }
  }

  /* With short study the accepted time difference must be shorter */
  f=0.01*dft1->x2[dft1->frameNr-1]*ts1;
  if(accepted_timedif>f) accepted_timedif=f;
  if(verbose>1) printf("accepted_timedif := %g [s]\n", accepted_timedif);

  /* Compare sample times frame-by-frame */
  for(fi=0, n=0; fi<smaller_frameNr; fi++) {
    if(dft1->timetype==DFT_TIME_MIDDLE) { // check frame mid times
      f=fabs(dft1->x[fi]*ts1-dft2->x[fi]*ts2);
      if(verbose>10) printf("timedif[%d] := %g\n", fi, f);
      if(verbose>12) printf("  %g vs %g\n", ts1*dft1->x[fi], ts2*dft2->x[fi]);
      if(f>accepted_timedif) n++;
      continue;
    }
    if(dft1->timetype==DFT_TIME_START || dft1->timetype==DFT_TIME_STARTEND) {
      f=fabs(dft1->x1[fi]*ts1-dft2->x1[fi]*ts2);
      if(verbose>10) printf("timedif[%d] := %g\n", fi, f);
      if(verbose>12) printf("  %g vs %g\n", ts1*dft1->x1[fi], ts2*dft2->x1[fi]);
      if(f>accepted_timedif) {n++; continue;}
    }
    if(dft1->timetype==DFT_TIME_END || dft1->timetype==DFT_TIME_STARTEND) {
      f=fabs(dft1->x2[fi]*ts1-dft2->x2[fi]*ts2);
      if(verbose>10) printf("timedif[%d] := %g\n", fi, f);
      if(verbose>12) printf("  %g vs %g\n", ts1*dft1->x2[fi], ts2*dft2->x2[fi]);
      if(f>accepted_timedif) n++;
    }
  }
  if(verbose>2) printf("nr of different frame times := %d\n", n);

  if(n==0) return 1; else return 0; 
} 
/*****************************************************************************/

/*****************************************************************************/
/** Copies frame times (especially start and end times, but also mid times) from an IMG data into 
    DFT data, and sets DFT 'header' to indicate that frame start and end times are present.
   @sa check_times_dft_vs_img, dftEndtime, imgEndtime, dftMatchTimeunits
   @return Returns 0 if successful and <>0 in case of an error.
 */
int copy_times_from_img_to_dft(
  /** Pointer to IMG data; times must be in sec as usual. */
  IMG *img,
  /** Pointer to DFT data; times can be both seconds or minutes, if unit is correctly set; 
      sample nr may be smaller than in IMG. */
  DFT *dft,
  /** Verbose level; if zero, then nothing is printed to stderr or stdout. */
  int verbose
) {
  int fi, times_changed=0;

  if(verbose>0) printf("%s(*img, *dft)\n", __func__);
  if(img==NULL || dft==NULL) return 1;
  if(img->dimt<=0 || dft->frameNr<=0) return 0;
  
  /* If more DFT samples, then that is an error */
  if(img->dimt<dft->frameNr) return 2;

  /* Convert DFT times to sec if necessary */
  if(dft->timeunit==TUNIT_MIN) {dftMin2sec(dft); times_changed=1;}

  /* Copy the frame times frame-by-frame */
  for(fi=0; fi<dft->frameNr; fi++) {
    dft->x1[fi]=img->start[fi];
    dft->x2[fi]=img->end[fi];
    dft->x[fi]=img->mid[fi];
  }
  dft->timetype=DFT_TIME_STARTEND;

  /* Convert DFT times back to min if necessary */
  if(times_changed!=0) dftSec2min(dft);

  return 0; 
} 
/*****************************************************************************/

/*****************************************************************************/
/** Returns the actual TAC sample number, not including NaNs, samples with negative x, 
    duplicate samples, or samples with zero weights (if data is weighted).
   @sa dftEndtime, fittime_from_dft
   @return Returns sample number.
 */
int getActualSamplenr(
  /** Pointer to TAC data in DFT struct; must be sorted by increasing x. */
  DFT *dft,
  /** Region index [0..voiNr-1]. */
  int ri
) {
  int n, fi;
  double last_x=-1.0E+99;

  if(dft==NULL || dft->frameNr<1 || ri<0 || ri>=dft->voiNr) return 0;
  for(fi=n=0; fi<dft->frameNr; fi++) {
    if(isnan(dft->x[fi])) continue;
    if(isnan(dft->voi[ri].y[fi])) continue;
    if(dft->x[fi]<0.0) continue;
    if(dft->isweight!=0 && dft->w[fi]<=0.0) continue;
    if(dft->x[fi]==last_x) continue;
    n++; last_x=dft->x[fi];
  }
  return n;
}
/*****************************************************************************/

/*****************************************************************************/
/** Get TAC end time. Sample times are assumed to be sorted to increasing order.
   @sa imgEndtime, fittime_from_dft, getActualSamplenr
   @return Returns the TAC end time, not converting the time units.
 */
double dftEndtime(
  /** Pointer to DFT TAC structure. */
  DFT *dft
) {
  int fi;
  if(dft==NULL || dft->frameNr<1) return(0.0);
  for(fi=dft->frameNr-1; fi>=0; fi--) {
    if(dft->timetype==DFT_TIME_MIDDLE) {
      if(!isnan(dft->x[fi])) return(dft->x[fi]);
    } else if(dft->timetype==DFT_TIME_STARTEND) {
      if(!isnan(dft->x2[fi])) return(dft->x2[fi]);
    } else if(dft->timetype==DFT_TIME_START) {
      if(!isnan(dft->x[fi])) return(dft->x[fi]);
    } else if(dft->timetype==DFT_TIME_END) {
      if(!isnan(dft->x[fi])) return(dft->x[fi]);
    } else
      return(0.0);
  }
  return(0.0);
}
/*****************************************************************************/
/** Get IMG end time. Frame times are assumed to be sorted to increasing order.
   @sa dftEndtime, fittime_from_img
   @return Returns the last frame end time. By default, times are in sec.
 */
double imgEndtime(
  /** Pointer to IMG structure. */
  IMG *img
) {
  int fi;
  if(img==NULL || img->dimt<1) return(0.0);
  for(fi=img->dimt-1; fi>=0; fi--)
    if(!isnan(img->end[fi])) return(img->end[fi]);
  return(0.0);
}
/*****************************************************************************/

/*****************************************************************************/
/** Make sure that time units in two DFT structs are the same, converting units when necessary, 
    optionally saving original units so that units can be converted back to what they were using 
    dftTimeunitConversion().
   @sa dftTimeunitConversion, dftEndtime, fittime_from_dft
   @return Returns 0 if the was no need for time unit conversion, 1 if time units in dft2 were 
    converted, and <0 if units were not identified or in case of an error.
 */ 
int dftMatchTimeunits(
  /** Pointer to DFT struct 1. */
  DFT *dft1,
  /** Pointer to DFT struct 2; sample time units are changed to match the data in dft1.  */
  DFT *dft2,
  /** Pointer for original time unit in DFT 2; enter NULL if not needed. */
  int *tunit2,
  /** Verbose level; if zero, then nothing is printed to stderr or stdout. */
  int verbose
) {
  int ret;

  if(verbose>0) printf("%s()\n", __func__);
  if(dft1==NULL || dft2==NULL) return(-1);
  /* Save original units if required */
  if(tunit2!=NULL) *tunit2=dft2->timeunit;
  /* Check that time units are available */
  if(dft1->timeunit!=TUNIT_MIN && dft1->timeunit!=TUNIT_SEC) {
    if(verbose>0) printf("  unknown time units in dft1\n");
    return(-2);
  }
  if(dft2->timeunit!=TUNIT_MIN && dft2->timeunit!=TUNIT_SEC) {
    if(verbose>0) printf("  unknown time units in dft2\n");
    return(-3);
  }
  /* Check if they are the same */
  if(dft1->timeunit==dft2->timeunit) {
    if(verbose>1) printf("  time units are the same in dft1 and dft2.\n");
    return(0);
  }
  /* Conversion to dft2 */
  if(verbose>1) printf("  time units in dft2 converted from %s to %s\n",
    petTunit(dft2->timeunit), petTunit(dft1->timeunit));
  ret=dftTimeunitConversion(dft2, dft1->timeunit);
  if(ret!=0) {
    if(verbose>0) printf("  time unit conversion failed\n");
    return(-10-ret);
  }
  return(1);
}
/*****************************************************************************/

/*****************************************************************************/
