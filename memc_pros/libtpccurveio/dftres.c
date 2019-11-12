/// @file dftres.c
/// @author Vesa Oikonen
/// @brief Utilities for setting up results structure based on DFT.
///
/*****************************************************************************/
#include "libtpccurveio.h"
/*****************************************************************************/

/*****************************************************************************/
/** Allocate memory for regional results based on information in DFT.
    @sa fit_allocate_with_dft, dftToResult
    @return Returns 0 if successful, otherwise <>0.
 */
int res_allocate_with_dft(
  /** Pointer to initiated RES struct which will be allocated here and filled with ROI names etc. */
  RES *res,
  /** Regional data from where necessary information is read. */
  DFT *dft
) {
  int ri;

  //printf("res_allocate_with_dft()\n"); fflush(stdout);
  // Check the input data
  if(res==NULL || dft==NULL || dft->voiNr<1) return 1;
  // Allocate memory
  if(resSetmem(res, dft->voiNr)!=0) return 2;
  res->voiNr=dft->voiNr;
  // Copy header information
  strcpy(res->studynr, dft->studynr);
  res->Vb=-1.0;
  res->fA=-1.0;
  res->E=-1.0;
  res->time=time(NULL); // Set current time to results
  res->isweight=dft->isweight;
  /* Copy region names, etc */
  for(ri=0; ri<dft->voiNr; ri++) {
    strcpy(res->voi[ri].name, dft->voi[ri].name);
    strcpy(res->voi[ri].voiname, dft->voi[ri].voiname);
    strcpy(res->voi[ri].hemisphere, dft->voi[ri].hemisphere);
    strcpy(res->voi[ri].place, dft->voi[ri].place);
  }
  /* Set data range */
  if(dft->timetype==DFT_TIME_STARTEND)
    sprintf(res->datarange, "%g - %g %s",
      dft->x1[0], dft->x2[dft->frameNr-1], petTunit(dft->timeunit) );
  else
    sprintf(res->datarange, "%g - %g %s",
      dft->x[0], dft->x[dft->frameNr-1], petTunit(dft->timeunit) );
  res->datanr=dft->frameNr;

  return 0;
}
/*****************************************************************************/

/*****************************************************************************/
/** Copy the contents (both header and data) of DFT struct into RES struct.
    @sa fitToResult, res_allocate_with_dft
    @return Returns 0 if successful, in case of an error >0, and <0 if warning is suggested.
 */
int dftToResult(
  /** Regional data from where necessary information is read. */
  DFT *dft,
  /** Pointer to initiated RES struct which will be allocated here. */
  RES *res,
  /** Pointer to a string (allocated for at least 64 chars) where error message
      or other execution status will be written; enter NULL, if not needed. */
  char *status
) {
  int ri, fi;

  // Check the input data
  if(status!=NULL) sprintf(status, "program error");
  if(res==NULL || dft==NULL || dft->voiNr<1 || dft->frameNr<1) return 1;
  // Allocate memory and copy most of headers
  if(res_allocate_with_dft(res, dft) != 0) {
    if(status!=NULL) sprintf(status, "cannot setup results data");
    return 2;
  }
  res->parNr=dft->frameNr; if(res->parNr>MAX_RESPARAMS) {
    sprintf(status, "only %d frames can be copied to results", MAX_RESPARAMS);
    res->parNr=MAX_RESPARAMS;
  }
  // Set parameter titles and units
  for(fi=0; fi<res->parNr; fi++) {
    sprintf(res->parname[fi], "%d", fi+1);
    strcpy(res->parunit[fi], dft->unit);
  }
  // Copy regional values
  for(ri=0; ri<dft->voiNr; ri++) for(fi=0; fi<res->parNr; fi++)
    res->voi[ri].parameter[fi]=dft->voi[ri].y[fi];

  if(dft->frameNr>MAX_RESPARAMS) return -1;
  /* Set also deprecated parameter name and unit representations, for now */
  resFixParnames(res);
  return 0;
}
/*****************************************************************************/

/*****************************************************************************/
