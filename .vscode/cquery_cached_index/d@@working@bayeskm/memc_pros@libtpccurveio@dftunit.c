/// @file dftunit.c
/// @author Vesa Oikonen
/// @brief Setting DFT calibration unit.
///
/*****************************************************************************/
#include "libtpccurveio.h"
/*****************************************************************************/

/*****************************************************************************/
/** Set DFT calibration unit string */
void dftUnitToDFT(DFT *dft, int dunit)
{
  strcpy(dft->unit, dftUnit(dunit));
}
/*****************************************************************************/

/*****************************************************************************/
/** Conversion of the DFT calibration unit.
 *  Changes both data values and unit string. 
 *  Currently available conversions are:
 *    MBq/cc <-> kBq/cc <-> Bq/cc <-> nCi/cc <-> uCi/cc
 *    Bq <-> kBq <-> MBq <-> GBq <-> nCi <-> uCi <-> mCi
\return Returns 0 if conversion was successful, and <> 0 if failed.
 */   
int dftUnitConversion(
  /** Pointer to existing DFT data, whose calibration unit is to be changed */
  DFT *dft,
  /** New unit code */
  int dunit
) {
  int current_dunit=CUNIT_UNKNOWN;
  double f=1.0;
  int ri, fi;
  int unit_type=0; /* 0=unknown; 1=concentration; 2=dose */

  /* Check the input */
  if(dft==NULL || dunit<0) return(1);
  
  /* Identify the current unit */
  current_dunit=petCunitId(dft->unit);
  if(current_dunit==CUNIT_UNKNOWN) return(2);
  /* If unit is already the requested one, then return */
  if(current_dunit==dunit) return(0);

  /* Determine the conversion factor to kBq/cc or MBq */
  switch(current_dunit) {
    case CUNIT_BQ_PER_ML:  f*=0.001; unit_type=1; break;
    case CUNIT_KBQ_PER_ML: f*=1.0; unit_type=1; break;
    case CUNIT_MBQ_PER_ML: f*=1000.0; unit_type=1; break;
    case CUNIT_NCI_PER_ML: f*=0.037; unit_type=1; break;
    case CUNIT_UCI_PER_ML: f*=37.0; unit_type=1; break;
    case CUNIT_BQ:         f*=0.000001; unit_type=2; break;
    case CUNIT_KBQ:        f*=0.001; unit_type=2; break;
    case CUNIT_MBQ:        f*=1.0; unit_type=2; break;
    case CUNIT_GBQ:        f*=1000.0; unit_type=2; break;
    case CUNIT_NCI:        f*=0.000037; unit_type=2; break;
    case CUNIT_UCI:        f*=0.037; unit_type=2; break;
    case CUNIT_MCI:        f*=37.0; unit_type=2; break;
    default: return(2);
  }

  /* Determine the conversion factor from kBq/cc or MBq to the required unit */
  if(unit_type==1) switch(dunit) { // concentrations
    case CUNIT_BQ_PER_ML:  f*=1000.0; break;
    case CUNIT_KBQ_PER_ML: f*=1.0; break;
    case CUNIT_MBQ_PER_ML: f*=0.001; break;
    case CUNIT_NCI_PER_ML: f/=0.037; break;
    case CUNIT_UCI_PER_ML: f/=37.0; break;
    default: return(3);
  } else if(unit_type==2) switch(dunit) { // doses
    case CUNIT_BQ:         f*=1000000.0; break;
    case CUNIT_KBQ:        f*=1000.0; break;
    case CUNIT_MBQ:        f*=1.0; break;
    case CUNIT_GBQ:        f*=0.001; break;
    case CUNIT_NCI:        f/=0.000037; break;
    case CUNIT_UCI:        f/=0.037; break;
    case CUNIT_MCI:        f/=37.0; break;
    default: return(3);
  } else
    return(3);

  /* Convert the data values */
  if(f==1.0) return(0);
  for(fi=0; fi<dft->frameNr; fi++) {
    for(ri=0; ri<dft->voiNr; ri++) {
      if(!isnan(dft->voi[ri].y[fi])) dft->voi[ri].y[fi]*=f;
      if(!isnan(dft->voi[ri].y2[fi])) dft->voi[ri].y2[fi]*=f;
      if(!isnan(dft->voi[ri].y3[fi])) dft->voi[ri].y3[fi]*=f;
    }
  }
  
  /* Set the new unit string */
  dftUnitToDFT(dft, dunit);

  return(0);
} 
/*****************************************************************************/

/*****************************************************************************/
/** Set DFT timeunit; does not change the sample times.
\return Returns 0 if successful, and <> 0 if invalid timeunit.
 */
int dftTimeunitToDFT(DFT *dft, const char *timeunit)
{
  if(dft==NULL || timeunit==NULL) return(1);
  int tunit=petTunitId(timeunit);
  if(tunit<0) return(2); else dft->timeunit=tunit;
  return(0);
}
/*****************************************************************************/

/*****************************************************************************/
/** Conversion of the DFT timeunit.
 *  Changes both data values and timeunit code. 
 *  Currently available conversions are:
 *    min <-> sec
\return Returns 0 if conversion was successful, and <> 0 if failed.
 */   
int dftTimeunitConversion(
  /** Pointer to existing DFT data, whose timeunit is to be changed */
  DFT *dft,
  /** New timeunit code */
  int tunit
) {
  /* Check the input */
  if(dft==NULL || tunit<0) return(1);
  
  /* Check if unit already is as required */
  if(dft->timeunit==tunit) return(0);
  
  /* Do the conversion, if supported */
  if(dft->timeunit==TUNIT_MIN && tunit==TUNIT_SEC)
    dftMin2sec(dft);
  else if(dft->timeunit==TUNIT_SEC && tunit==TUNIT_MIN)
    dftSec2min(dft);
  else
    return(2);

  return(0);
} 
/*****************************************************************************/

/*****************************************************************************/
/** Change time unit from min to sec, without checking original unit. */
void dftMin2sec(DFT *dft)
{
  int i;

  for(i=0; i<dft->frameNr; i++) {
    if(!isnan(dft->x[i]))  dft->x[i]*=60.;
    if(!isnan(dft->x1[i])) dft->x1[i]*=60.;
    if(!isnan(dft->x1[i])) dft->x2[i]*=60.;
  }
  dft->timeunit=TUNIT_SEC;
}
/*****************************************************************************/

/*****************************************************************************/
/** Change time unit from sec to min, without checking original unit. */
void dftSec2min(DFT *dft)
{
  int i;

  for(i=0; i<dft->frameNr; i++) {
    if(!isnan(dft->x[i]))  dft->x[i]/=60.; 
    if(!isnan(dft->x1[i])) dft->x1[i]/=60.; 
    if(!isnan(dft->x1[i])) dft->x2[i]/=60.;
  }
  dft->timeunit=TUNIT_MIN;
}
/*****************************************************************************/

/*****************************************************************************/
