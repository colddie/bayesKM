/// @file resift.c
/// @author Vesa Oikonen
/// @brief Utility functions for working with RES and IFT struct.
///
/*****************************************************************************/
#include "libtpccurveio.h"
/*****************************************************************************/

/*****************************************************************************/
/** Copy results in RES struct into IFT struct.
 *  @return Returns 0 when successful. 
 */ 
int res2ift(
  /** Pointer to RES struct */
  RES *res,
  /** Pointer to initiated IFT struct */
  IFT *ift,
  /** Verbose level; if zero, then nothing is printed into stdout or stderr */
  int verbose
) {
  int ri, pi, ret;
  char tmp[1024], tmp2[1024], tmp3[1024];
  struct tm st;

  if(verbose>0) printf("res2ift()\n");

  /* Check the input */
  if(res==NULL || ift==NULL) {
    if(verbose>0) fprintf(stderr, "Error: invalid input\n");
    return(1);
  }

  /* Delete any previous content */
  iftEmpty(ift);
  
  /* Set header info */
  if(strlen(res->program)>0) {
    ret=iftPut(ift, "program", res->program, NULL);
    if(ret!=0) {
      if(verbose>0) fprintf(stderr, "Error: cannot set IFT content\n");
      return(3);
    }
  }
  if(gmtime_r(&res->time, &st)!=NULL) {
    strftime(tmp, 256, "%Y-%m-%d %H:%M:%S", &st);
    iftPut(ift, "date", tmp, NULL);
  }
  if(res->studynr[0]) iftPut(ift, "studynr", res->studynr, NULL);
  if(res->datafile[0]) iftPut(ift, "datafile", res->datafile, NULL);
  if(res->plasmafile[0]) iftPut(ift, "plasmafile", res->plasmafile, NULL);
  if(res->plasmafile2[0]) iftPut(ift, "plasmafile2", res->plasmafile2, NULL);
  if(res->bloodfile[0]) iftPut(ift, "bloodfile", res->bloodfile, NULL);
  if(res->reffile[0]) iftPut(ift, "reffile", res->reffile, NULL);
  if(res->refroi[0]) iftPut(ift, "refroi", res->refroi, NULL);

  if(res->datarange[0]) iftPut(ift, "datarange", res->datarange, NULL);
  if(res->datanr>0) {
    snprintf(tmp, 1024, "%d", res->datanr); iftPut(ift, "datanr", tmp, NULL);}
  if(res->fitmethod[0]) iftPut(ift, "fitmethod", res->fitmethod, NULL);

  if(res->density>0) {
    snprintf(tmp, 1024, "%g", res->density); iftPut(ift, "density", tmp, NULL);}
  if(res->lc>0) {
    snprintf(tmp, 1024, "%g", res->lc); iftPut(ift, "lc", tmp, NULL);}
  if(res->concentration>0) {
    snprintf(tmp, 1024, "%g", res->concentration);
    iftPut(ift, "concentration", tmp, NULL);
  }
  if(res->beta>0) {
    snprintf(tmp, 1024, "%g", res->beta); iftPut(ift, "beta", tmp, NULL);}
  if(res->Vb>0) {
    snprintf(tmp, 1024, "%g", res->Vb); iftPut(ift, "Vb", tmp, NULL);}
  if(res->fA>0) {
    snprintf(tmp, 1024, "%g", res->fA); iftPut(ift, "fA", tmp, NULL);}
  if(res->E>0) {
    snprintf(tmp, 1024, "%g", res->E); iftPut(ift, "E", tmp, NULL);}
  if(res->isweight>0) strcpy(tmp, "yes");
  else if(res->isweight==0) strcpy(tmp, "no");
  else strcpy(tmp, "unknown");
  iftPut(ift, "weighting", tmp, NULL);

  /* Set each ROI and parameter value */
  for(ri=0; ri<res->voiNr; ri++) {
    if(res->voiNr>1) {
      strcpy(tmp, res->voi[ri].name); rnameRmDots(tmp, NULL);
      strcat(tmp, "_");
    } else {
      strcpy(tmp, "");
    }
    for(pi=0; pi<res->parNr; pi++) {
      snprintf(tmp2, 1024, "%s%s", tmp, res->parname[pi]);
      if(isnan(res->voi[ri].parameter[pi])) strcpy(tmp3, "");
      else sprintf(tmp3, "%g", res->voi[ri].parameter[pi]);
      if(strlen(res->parunit[pi])) {
        strcat(tmp3, " "); strcat(tmp3, res->parunit[pi]);}
      iftPut(ift, tmp2, tmp3, NULL);
      if(!isnan(res->voi[ri].sd[pi])) {
        snprintf(tmp2, 1024, "%s%s_%s", tmp, res->parname[pi], "SD");
        snprintf(tmp3, 1024, "%g", res->voi[ri].sd[pi]);
        if(strlen(res->parunit[pi])) {
          strcat(tmp3, " "); strcat(tmp3, res->parunit[pi]);}
        iftPut(ift, tmp2, tmp3, NULL);
      }
      if(!isnan(res->voi[ri].cl1[pi])) {
        snprintf(tmp2, 1024, "%s%s_%s", tmp, res->parname[pi], "CL1");
        snprintf(tmp3, 1024, "%g", res->voi[ri].cl1[pi]);
        if(strlen(res->parunit[pi])) {
          strcat(tmp3, " "); strcat(tmp3, res->parunit[pi]);}
        iftPut(ift, tmp2, tmp3, NULL);
      }
      if(!isnan(res->voi[ri].cl2[pi])) {
        snprintf(tmp2, 1024, "%s%s_%s", tmp, res->parname[pi], "CL2");
        snprintf(tmp3, 1024, "%g", res->voi[ri].cl2[pi]);
        if(strlen(res->parunit[pi])) {
          strcat(tmp3, " "); strcat(tmp3, res->parunit[pi]);}
        iftPut(ift, tmp2, tmp3, NULL);
      }
    }
  }

  return(0);
}
/*****************************************************************************/

/*****************************************************************************/
