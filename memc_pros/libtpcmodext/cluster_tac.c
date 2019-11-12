/// @file cluster_tac.c
/// @brief Clustering and segmentation for PET modeling.
/// @author Vesa Oikonen
///
/*****************************************************************************/

/*****************************************************************************/
#include "libtpcmodext.h"
/*****************************************************************************/

/*****************************************************************************/
/** Allocates memory and calculates the values for average TACs for clusters.
\return Returns 0 if ok, and >0 in case of an error.
 */
int clusterTACs(
  /** Dynamic image */
  IMG *dimg,
  /** Cluster image */
  IMG *cimg,
  /** Highest cluster ID */
  int nr,   
  /** Pointer to initiated but empty DFT data */
  DFT *tac,
  /** Verbose level; if zero, then only warnings are printed into stderr */
  int verbose
) {
  int fi, ret, clusterID, n;

  if(verbose>0)
    printf("clusterTACs(dimg, cimg, %d, tac, %d)\n", nr, verbose);
  /* Check the arguments */
  if(dimg==NULL || cimg==NULL || nr<1 || tac==NULL) return(1);
  if(dimg->dimt<1 || cimg->dimt<1) return(1);
  if(cimg->dimx!=dimg->dimx || cimg->dimy!=dimg->dimy || cimg->dimz!=dimg->dimz)
    return(2);

  /* Allocate memory for the TACs */
  dftEmpty(tac);
  ret=dftSetmem(tac, dimg->dimt, nr+1); if(ret) return(3);
  
  /* Set TAC info */
  tac->voiNr=0; tac->frameNr=dimg->dimt; tac->_type=1;
  for(fi=0; fi<tac->frameNr; fi++) {
    tac->x1[fi]=dimg->start[fi];
    tac->x2[fi]=dimg->end[fi];
    tac->x[fi] =dimg->mid[fi];
  }
  tac->timetype=DFT_TIME_STARTEND;
  tac->timeunit=TUNIT_SEC;
  strcpy(tac->unit, imgUnit(dimg->unit));

  /* Calculate one cluster at a time */
  float y[dimg->dimt];
  for(clusterID=1; clusterID<=nr; clusterID++) {
    char buf[128]; snprintf(buf, 128, "%06d", clusterID);
    char *p=buf+strlen(buf)-6;
    snprintf(tac->voi[clusterID-1].voiname, MAX_REGIONSUBNAME_LEN+1, "%s", p);
    n=imgsegmClusterMean(dimg, cimg, clusterID, y, verbose);
    if(verbose>1) printf("  clusterID%d -> %d pixels\n", clusterID, n);
    if(n<0) return(5); else if(n==0) return(6);
    for(fi=0; fi<tac->frameNr; fi++) tac->voi[clusterID-1].y[fi]=(double)y[fi];
    tac->voi[clusterID-1].size=n*dimg->sizex*dimg->sizey*dimg->sizez;
    tac->voiNr++;
  }
  /* and once more for cluster 0, i.e. the thresholded pixels;  */
  /* note that it is possible that there is no cluster 0 at all */
  clusterID=0;
  sprintf(tac->voi[tac->voiNr].voiname, "%06d", clusterID);
  strcpy(tac->voi[tac->voiNr].name, tac->voi[tac->voiNr].voiname);
  ret=imgsegmClusterMean(dimg, cimg, clusterID, y, verbose);
  if(ret<0) return(7);
  if(ret>0) {
    for(fi=0; fi<tac->frameNr; fi++) tac->voi[tac->voiNr].y[fi]=(double)y[fi];
    tac->voi[tac->voiNr].size=ret*dimg->sizex*dimg->sizey*dimg->sizez;
    tac->voiNr++;
  }

  return(0);
}
/*****************************************************************************/

/*****************************************************************************/
