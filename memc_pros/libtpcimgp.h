/// @file libtpcimgp.h
/// @brief Header file for libtpcimgp.
/// @author Vesa Oikonen
///
#ifndef _LIBTPCIMGP_H
#define _LIBTPCIMGP_H
/*****************************************************************************/

/*****************************************************************************/
#include "tpcclibConfig.h"
/*****************************************************************************/
//#include <omp.h>
/*****************************************************************************/
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <strings.h>
#include <ctype.h>
#include <float.h>
#include <unistd.h>
#include <locale.h>
#include <time.h>
/*****************************************************************************/
#include "libtpcmisc.h"
#include "libtpcimgio.h"
/*****************************************************************************/

/*****************************************************************************/
/** PET image color scale */
#define	PET_GRAYSCALE     1
/** PET image color scale */
#define	PET_GRAYSCALE_INV 2
/** PET image color scale */
#define PET_RAINBOW       3
/** PET image color scale */
#define PET_RAINBOW_WB    4
/*****************************************************************************/
/** 3D image point */
typedef struct point{
  /** x coordinate */
  float x;
  /** y coordinate */
  float y;
  /** z coordinate */
  float z;
} point;
/*****************************************************************************/

/*****************************************************************************/
/* imgarithm */
int imgArithm(IMG *img1, IMG *img2, char operation, float ulimit, int verbose);
int imgArithmConst(IMG *img, float operand, char operation, float ulimit, int verbose);
int imgArithmFrame(IMG *img1, IMG *img2, char operation, float ulimit, int verbose);
int imgLn(IMG *img);
int imgLog10(IMG *img);
int imgAbs(IMG *img);
int imgInv(IMG *img);
int imgFrameIntegral(IMG *img, int first, int last, IMG *iimg, int verbose);
int imgRawCountsPerTime(IMG *img, int operation);
int imgConvertUnit(IMG *img, char *unit);
/*****************************************************************************/

/*****************************************************************************/
/* imgeval.c */
int imgAverageTAC(IMG *img, float *tac);
int imgAverageMaskTAC(IMG *img, IMG *timg, float *tac);
int imgAverageAUC(IMG *img, float *avgauc);
int imgMaskTAC(IMG *img, IMG *mask, double *tac, int verbose);
int imgMaskRoiNr(IMG *img, INTEGER_LIST *list);
int imgVoiMaskTAC(IMG *img, IMG *mask, int mv, double *tac, int verbose);
/*****************************************************************************/

/*****************************************************************************/
/* imgfilter */
int imgFillGaussKernel(
  float **kernel, float stdev, int size
);
int imgConvolute2D(
  float ***data, float **buffer, int frame, int width, int height,
  float **kernel, int size, int border, int verbose, char *errmsg
); 
int imgGaussianFilter(
  IMG *img, int plane, int frame, float gauss_sd, int size, int border,
  int verbose, char *errmsg 
);
int imgFast2DGaussianFilter(
  IMG *img, int plane, int frame, float gauss_sd, int step_nr,
  int verbose, char *errmsg 
);
int imgFast3DGaussianFilter(
  IMG *img, int frame, float gauss_sd, int step_nr, int verbose, char *errmsg
);
int imgFast1DGaussianFilter(
  IMG *img, float gauss_sd, int step_nr, int verbose, char *errmsg
);
/// @cond
float **mallocMatrix(float w, float h);
/// @endcond
float **imgGaussKernel(int size);
void imgFreeKernel(float **kernel, int size);
void imgConvoluteData(float ***data, float **buffer, int frame, 
  int width, int height, float **kernel, int size
);
int imgConvolute(IMG *img, int frame, int plane, float **kernel, int size);
/*****************************************************************************/

/*****************************************************************************/
/* imgflips */
void imgFlipHorizontal(IMG *img);
void imgFlipVertical(IMG *img);
void imgFlipPlanes(IMG *img);
int imgFlipRight(IMG *img);
int imgFlipAbove(IMG *img);
/*****************************************************************************/

/*****************************************************************************/
/* imgframe */
int imgFramesCheck(IMG *img, int verbose);
int imgFrameGapFill(IMG *img, int verbose);
int imgDeleteFrameOverlap_old(IMG *img);
int imgDeleteFrameOverlap(IMG *img);
int imgSmoothOverFrames(IMG *img, int n);
/* Deprecated functions. Please don't use these anymore */
/// @cond
#define NEW_imgDeleteFrameOverlap imgDeleteFrameOverlap
/// @endcond
/*****************************************************************************/

/*****************************************************************************/
/* imgscanner */
int imgSetScanner(IMG *img, int scanner_type);
/*****************************************************************************/

/*****************************************************************************/
/* imgsegm */
int imgsegmThresholdMask(
  IMG *img, float minValue, float maxValue, IMG *timg
);
int imgsegmThresholdByMask(
  IMG *img, IMG *templatexx, float minValue, float maxValue     // conflict with C++ template..
);
int imgsegmThreshold(
  IMG *img, float minValue, float maxValue
);
int imgsegmMaskToCluster(
  IMG *img
);
int imgsegmFindMaxOutsideClusters(
  IMG *sumimg, IMG *cluster, float *max, int *plane, int *row, int *col
);
int imgsegmClusterExpand(
  IMG *cluster, IMG *sum, IMG *dynamic, int clusterID,
  int pi, int ri, int ci, int pj, int rj, int cj, float CVlim, float CClim,
  int verbose
);
float imgsegmPearson(
  float *x, float *y, int nr
);
int imgsegmClusterMean(
  IMG *dimg, IMG *cimg, int clusterID, float *avg, int verbose
);
int imgsegmCheckNeighbours(
  IMG *cimg, int pi, int ri, int ci
);
int imgsegmFindBestNeighbour(
  IMG *dimg, IMG *cimg, int pi, int ri, int ci
);
int imgsegmSimilar(
  IMG *input, int smoothDim, int smoothNr, IMG *output, int verbose
);
int imgsegmCalcMRL(
  float y1[], float y2[], int n
);
/*****************************************************************************/

/*****************************************************************************/
/* imgthreshold */
int imgThresholding(IMG *img, float threshold_level, int *thr_nr);
int imgThresholdingLowHigh(
  IMG *img, float lower_threshold_level, float upper_threshold_level,
  IMG *timg, int *lower_thr_nr, int *upper_thr_nr
);
int imgOutlierFilter(IMG *img, float limit);
int imgThresholdMaskCount(
  IMG *img, float minValue, float maxValue, IMG *timg, int *count
);
int imgThresholdMask(
  IMG *img, float minValue, float maxValue, IMG *timg
);
int imgThresholdByMask(IMG *img, IMG *templt, float thrValue);
void imgCutoff(IMG *image, float cutoff, int mode);
/*****************************************************************************/

/*****************************************************************************/
/* imgtiff */
int tiffWriteImg(
  IMG *img, int plane, int frame, float *maxvalue, int colorscale,
  char *fname, int matXdim, int matYdim,
  int verbose, char *status
);
/*****************************************************************************/

/*****************************************************************************/
/* imgtransform */
int img2cube(IMG *img1, int dim, IMG *img2);
void imgScale(IMG *src, IMG *targ, float zoom, int method);
void integerScale(int frame, float ***src, float **targ, int width, int height, int zoom);
/*****************************************************************************/

/*****************************************************************************/
/* mask.c */
unsigned int imgMaskCount(IMG *img);
int imgMaskErode(IMG *img, IMG *se);
int imgMaskDilate(IMG *img, IMG *se);
int imgStructuringElement(IMG *img, const int structuring_element, int verbose);
void imgMaskInvert(IMG *img);
int imgMaskConjunction(IMG *mask1, IMG *mask2);
int imgMaskRegionLabeling(IMG *mask1, IMG *mask2, int *n, int verbose);
int imgMaskFloodFill(IMG *m, int sz, int sy, int sx, int label, int *n, int verbose);
/*****************************************************************************/

/*****************************************************************************/
/* point */
int pRound(float);
float getDistance(point, point);
float getAngle(point, point);
/*****************************************************************************/

/*****************************************************************************/
#endif // _LIBTPCIMGP_H
