/// @file noise.c
/// @brief Noise simulation for PET modelling.
/// @author Vesa Oikonen
///
/*****************************************************************************/

/*****************************************************************************/
#include "libtpcmodext.h"
/*****************************************************************************/

/*****************************************************************************/
/** Calculate SD for PET radioactivity concentration data to be used
    to simulate noise.

    Note that SD is dependent on the time units.

    Reference: Varga & Szabo. J Cereb Blood Flow Metab 2002;22(2):240-244.

    @return Returns 0 when successful, otherwise <>0.
 */
int noiseSD4Simulation(
  /** Sample radioactivity concentration (decay corrected to zero time) */
  double y,
  /** Radioactivity measurement (frame) start time
      (in same units as the halflife) */
  double t1,
  /** Radioactivity measurement (frame) duration
      (in same units as the halflife) */
  double dt,
  /** Isotope halflife (in same units as the sample time); enter 0 if not
      to be considered */
  double hl,
  /** Proportionality factor. Note that it is inside square root. */
  double a,
  /** Pointer in which SD is written */
  double *sd,
  /** Pointer to a string (allocated for at least 64 chars) where error message
      or other execution status will be written; enter NULL, if not needed */   
  char *status,
  /** Verbose level; if zero, then nothing is printed to stderr or stdout */
  int verbose
) {
  double d, lambda, f;

  if(verbose>0) printf("noiseSD4Simulation(%g, %g, %g, %g, %g, *sd, ...)\n",
                       y, t1, dt, hl, a);
  if(status!=NULL) strcpy(status, "invalid data");
  if(sd==NULL) return 1;
  if(t1<0.0) return 2;
  if(dt<=0.0) return 3;

  /* Decay factor (<=1) */
  if(hl<=0.0) {
    d=1.0; lambda=0.0;
  } else {
    if(status!=NULL) strcpy(status, "invalid half-life");
    lambda=hl2lambda(hl); if(lambda<0.0) return 4;
    d=hlLambda2factor(-lambda, t1, dt); if(d<0.0) return 5;
  }
  /* SD */
  f=y*d*dt;
  if(f<=0.0) *sd=0.0;
  else *sd=y*sqrt(a/f);

  if(status!=NULL) strcpy(status, "ok");
  return(0);
}
/*****************************************************************************/

/*****************************************************************************/
/** Calculate SD for noise simulation from TAC data.

    Sample times will be converted to minutes if necessary.

    Reference: Varga & Szabo. J Cereb Blood Flow Metab 2002;22(2):240-244.

    @return Returns 0 when successful, otherwise <>0.
 */
int noiseSD4SimulationFromDFT(
  /** Pointer to TAC data in DFT struct, based on which the SD for each frame
      is calculated. Contents are not changed.
      Struct must contain correct values for the isotope, time unit,
      and frame times (start and end). Status of decay correction is used,
      and if not known, then data is assumed to be decay corrected. */
  DFT *dft,
  /** TAC index [0..voiNr-1] which is used to calculate the SD in case DFT
      contains more than one TAC. Non-effective if DFT contains only one TAC.
      Set to <0 if mean of all TACs is to be used. */
  int index,
  /** Proportionality factor */
  double pc,
  /** Pointer to allocated array in which SD is written */
  double *sd,
  /** Pointer to a string (allocated for at least 64 chars) where error message
      or other execution status will be written; enter NULL, if not needed */   
  char *status,
  /** Verbose level; if zero, then nothing is printed to stderr or stdout */
  int verbose
) {
  int ret, fi;
  double hl=0.0, t1, deltat;
  double *y;
  DFT mean;

  if(verbose>0) printf("noiseSD4SimulationFromDFT(DFT, %d, %g, sd[], ...)\n",
                       index, pc);
  if(status!=NULL) strcpy(status, "invalid data");
  if(dft==NULL || sd==NULL) return 1;
  if(dft->voiNr<1 || dft->frameNr<1) return 2;
  if(dft->voiNr>1 && index>=dft->voiNr) return 3;

  /* Check that valid frame times are available */
  if(status!=NULL) strcpy(status, "invalid frame times");
  for(fi=0; fi<dft->frameNr; fi++)
    if(dft->x2[fi]<=dft->x1[fi] || dft->x1[fi]<0.0) return 4;
  if(status!=NULL) strcpy(status, "missing time unit");
  if(dft->timeunit!=TUNIT_SEC && dft->timeunit!=TUNIT_MIN) return 5;

  /* Check if isotope is available, if it is needed */
  if(dft->decayCorrected==DFT_DECAY_CORRECTED || 
     dft->decayCorrected==DFT_DECAY_UNKNOWN)
  {
    hl=hlFromIsotope(dft->isotope);
    if(status!=NULL) strcpy(status, "missing isotope halflife");
    if(hl<=0.0) return 6;
  }
  if(verbose>1) printf("halflife := %g\n", hl);

  /* Set pointer to activity concentration data */
  dftInit(&mean);
  if(dft->voiNr==1) {
    y=dft->voi[0].y;
  } else if(index>=0) {
    y=dft->voi[index].y;
  } else {
    /* Compute TAC mean */
    ret=dftMeanTAC(dft, &mean); if(ret!=0) {
      if(status!=NULL) strcpy(status, "cannot calculate mean TAC");
      return 8;
    }
    /* and use that */
    y=mean.voi[0].y;
  }

  /* Compute SD for each sample time */
  if(pc<=0.0) pc=1.0;
  for(fi=0; fi<dft->frameNr; fi++) {
    t1=dft->x1[fi]; deltat=dft->x2[fi]-dft->x1[fi];
    if(dft->timeunit==TUNIT_SEC) {t1/=60.0; deltat/=60.0;}
    ret=noiseSD4Simulation(y[fi], t1, deltat, hl, pc, sd+fi, status, verbose);
    if(ret!=0) break;
  }
  dftEmpty(&mean);
  if(ret!=0) {
    if(status!=NULL)
      sprintf(status, "cannot calculate SD for noise simulation (%d)", ret);
    return(100+ret);
  }

  if(status!=NULL) strcpy(status, "ok");
  return 0;
}
/*****************************************************************************/

/*****************************************************************************/
