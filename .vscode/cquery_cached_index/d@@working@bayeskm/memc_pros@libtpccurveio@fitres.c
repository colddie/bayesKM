/// @file fitres.c
/// @author Vesa Oikonen
/// @brief Utility functions for working with FIT struct.
/// @sa mathfunc.c
///
/*****************************************************************************/
#include "libtpccurveio.h"
/*****************************************************************************/

/*****************************************************************************/
/** Allocate memory for regional function fits based on information in DFT.
    @return Returns 0 if successful, otherwise <>0.
 */
int fit_allocate_with_dft(
  /** Pointer to initiated FIT struct which will be allocated here and filled with ROI names etc. */
  FIT *fit,
  /** Regional data from where necessary information is read. */
  DFT *dft
) {
  int ri;

  //printf("fit_allocate_with_dft()\n"); fflush(stdout);
  // Check the input data
  if(fit==NULL || dft==NULL || dft->voiNr<1) return 1;
  // Allocate memory
  if(fitSetmem(fit, dft->voiNr)!=0) return 2;
  fit->voiNr=dft->voiNr;
  /* Set header contents */
  fit->time=time(NULL); // Set current time to results
  strcpy(fit->unit, dft->unit);
  fit->timeunit=dft->timeunit;
  /* Copy region names, etc */
  for(ri=0; ri<dft->voiNr; ri++) {
    strcpy(fit->voi[ri].name, dft->voi[ri].name);
    strcpy(fit->voi[ri].voiname, dft->voi[ri].voiname);
    strcpy(fit->voi[ri].hemisphere, dft->voi[ri].hemisphere);
    strcpy(fit->voi[ri].place, dft->voi[ri].place);
    fit->voi[ri].dataNr=dft->frameNr;
    if(dft->timetype==DFT_TIME_STARTEND) {
      fit->voi[ri].start=dft->x1[0]; fit->voi[ri].end=dft->x2[dft->frameNr-1];
    } else {
      fit->voi[ri].start=dft->x[0]; fit->voi[ri].end=dft->x[dft->frameNr-1];
    }
  }

  return 0;
}
/*****************************************************************************/

/*****************************************************************************/
/** Conversion of FIT contents to RES.
    @sa res_allocate_with_dft, dftToResult
    @return Returns 0 when successful.
 */
int fitToResult(
  /** Pointer to FIT structure, contents of which are written to RES struct. */
  FIT *fit,
  /** Pointer to initiated RES struct where FIT contents are written;
      any previous contents are removed. */
  RES *res,
  /** Pointer to a string (allocated for at least 64 chars) where error message
      or other execution status will be written; enter NULL, if not needed. */
  char *status
) {
  int pi, ri, maxParNr, ret;


  //printf("fitToResult()\n"); fflush(stdout);
  /* Check input */
  if(status!=NULL) sprintf(status, "invalid data");
  if(fit==NULL || res==NULL) return 1;
  if(fit->voiNr<1) return 2;

  /* Determine max parameter number in fits */
  for(ri=0, maxParNr=0; ri<fit->voiNr; ri++)
    if(fit->voi[ri].parNr>maxParNr) maxParNr=fit->voi[ri].parNr;
  //printf("maxParNr := %d\n", maxParNr); fflush(stdout);

  /* Allocate memory for results */
  resEmpty(res); ret=resSetmem(res, fit->voiNr); 
  if(ret) {
    if(status!=NULL) sprintf(status, "cannot allocate memory");
    return 11;
  }
  /* Copy titles & filenames */
  //printf("copy header contents\n"); fflush(stdout);
  if(strlen(fit->program)>0 && strlen(fit->program)<512) 
    snprintf(res->program, 1024, "%.512s (c) 2014", fit->program);
  else 
    strcpy(res->program, "fitToResult (c) 2014");
  strcpy(res->datafile, fit->datafile);
  strcpy(res->studynr, fit->studynr);
  res->time=fit->time;
  /* Copy region names, etc */
  //printf("copy region names\n"); fflush(stdout);
  res->voiNr=fit->voiNr;
  for(ri=0; ri<fit->voiNr; ri++) {
    strcpy(res->voi[ri].name, fit->voi[ri].name);
    strcpy(res->voi[ri].voiname, fit->voi[ri].voiname);
    strcpy(res->voi[ri].hemisphere, fit->voi[ri].hemisphere);
    strcpy(res->voi[ri].place, fit->voi[ri].place);
  }
  /* Copy sample number, if equal in all TACs */
  for(ri=1, ret=0; ri<fit->voiNr; ri++)
    if(fit->voi[ri].dataNr!=fit->voi[0].dataNr) ret++;
  if(ret==0) res->datanr=fit->voi[0].dataNr;
  /* Set parameter names */
  //printf("set parameter names\n"); fflush(stdout);
  res->parNr=maxParNr+2; // function id and wss too
  strcpy(res->parname[0], "Func");
  for(pi=0; pi<maxParNr; pi++) sprintf(res->parname[pi+1], "p%d", pi+1);
  strcpy(res->parname[pi+1], "WSS");
  /* Copy parameter values */
  //printf("copy parameter values\n"); fflush(stdout);
  for(ri=0; ri<fit->voiNr; ri++) {
    res->voi[ri].parameter[0]=fit->voi[ri].type; // function id
    for(pi=0; pi<maxParNr; pi++) { // function parameters
      if(pi>=fit->voi[ri].parNr) res->voi[ri].parameter[pi+1]=0.0;
      else res->voi[ri].parameter[pi+1]=fit->voi[ri].p[pi];
    }
    res->voi[ri].parameter[pi+1]=fit->voi[ri].wss; // wss
  }
  /* Set also deprecated parameter name and unit representations, for now */
  //printf("set deprecated info\n"); fflush(stdout);
  resFixParnames(res);
  return(0);
}
/*****************************************************************************/

/*****************************************************************************/
