/// @file libtpcmodext.h
/// @brief Header file for libtpcmodext.
/// @author Vesa Oikonen
///
/*****************************************************************************/
#ifndef _LIBTPCMODEXT_H
#define _LIBTPCMODEXT_H
/*****************************************************************************/
#include "tpcclibConfig.h"
/*****************************************************************************/
//#include <omp.h>
/*****************************************************************************/
#include "libtpcmisc.h"
#include "libtpcmodel.h"
#include "libtpcsvg.h"
#include "libtpccurveio.h"
#include "libtpcimgio.h"
#include "libtpcimgp.h"
/*****************************************************************************/

/*****************************************************************************/
// bf_model.c
int bf_srtm(
  double *t, double *cr, int n, int bfNr, double t3min, double t3max, DFT *bf
);
int bfRadiowater(
  DFT *input, DFT *tissue, DFT *bf, int bfNr, double k2min, double k2max,
  char *status, int verbose
);
int bfIrr2TCM(
  DFT *input, DFT *tissue, DFT *bf, int bfNr, double thetamin, double thetamax,
  char *status, int verbose
);
/*****************************************************************************/

/*****************************************************************************/
// dftinput.c
int dftReadinput(
  DFT *input, DFT *tissue, char *filename, int *filetype,
  double *ti1, double *ti2, int verifypeak, char *status, int verbose
);
int dftReadReference(
  DFT *tissue, char *filename, int *filetype, int *ref_index, char *status,
  int verbose
);
int dftReadModelingData(
  char *tissuefile, char *inputfile1, char *inputfile2, char *inputfile3,
  double *fitdur, int *fitframeNr, DFT *tis, DFT *inp,
  FILE *loginfo, int verbose, char *status
);
int dftRobustMinMaxTAC(
  DFT *dft, int tacindex, double *minx, double *maxx, double *miny, double *maxy,
  int *mini, int *maxi, int *mins, int *maxs, int verbose
);
int dftVerifyPeak(
  DFT *dft, int index, int verbose, char *status
);
/*****************************************************************************/

/*****************************************************************************/
// imginput.c
int imgReadModelingData(
  char *petfile, char *siffile, char *inputfile1, char *inputfile2,
  char *inputfile3, double *fitdur, int *fitframeNr, IMG *img,
  DFT *inp, DFT *iinp, int verifypeak, FILE *loginfo,
  int verbose, char *status
);
/*****************************************************************************/

/*****************************************************************************/
// extrapolate.c
int extrapolate_monoexp(
  DFT *dft, double *fittime, int *min_nr, int max_nr, double mintime,
  double extr_to, DFT *ext, FILE *loginfo, char *status
);
int dftAutointerpolate(
  DFT *dft, DFT *dft2, double endtime, int verbose
);
int dftDoubleFrames(
  DFT *dft, DFT *dft2
);
int dftDivideFrames(
  DFT *dft, int voi_index, int add_nr, DFT *dft2, int verbose
);
int dft_end_line(DFT *dft, double *fittime, int *min_nr, int max_nr,
  double mintime, int check_impr, FIT *fit, FILE *loginfo, char *status
);
int dft_ln(DFT *dft1, DFT *dft2);
/*****************************************************************************/

/*****************************************************************************/
// fittime.c
int fittime_from_dft(
  DFT *dft, double *startTime, double *endTime, int *first, int *last, int verbose
);
int fittime_from_img(IMG *img, double *fittime, int verbose);
int check_times_dft_vs_img(IMG *img, DFT *dft, int verbose);
int check_times_dft_vs_dft(DFT *dft1, DFT *dft2, int verbose);
int copy_times_from_img_to_dft(IMG *img, DFT *dft, int verbose);
int getActualSamplenr(DFT *dft, int ri);
double dftEndtime(DFT *dft);
double imgEndtime(IMG *img);
int dftMatchTimeunits(DFT *dft1, DFT *dft2, int *tunit2, int verbose);
/*****************************************************************************/

/*****************************************************************************/
// img_mtga.c
int img_patlak(
  DFT *input, IMG *dyn_img, int start, int end, linefit_range fit_range, 
  IMG *ki_img, IMG *ic_img, IMG *nr_img,
  char *status, int verbose
);
int img_logan(
  DFT *input, IMG *dyn_img, int start, int end, linefit_range fit_range,
  double k2, IMG *vt_img, IMG *ic_img, IMG *nr_img,
  char *status, int verbose
);
/*****************************************************************************/

/*****************************************************************************/
// img_k1.c
int img_k1_using_ki(
  DFT *input, IMG *dyn_img, int frame_nr, IMG *ki_img, IMG *k1_img, IMG *k2k3_img,
  char *status, int verbose
);
/*****************************************************************************/
// dftint.c
int dftInterpolateCheckStart(
  DFT *input, DFT *output, char *status, int verbose
);
int dftInterpolateCheckEnd(
  DFT *input, DFT *output, char *status, int verbose
);
int dftInterpolate(
  DFT *input, DFT *tissue, DFT *output, char *status, int verbose
);
int dftInterpolateInto(
  DFT *inp, DFT *tis, char *status, int verbose
);
int dftTimeIntegral(
  DFT *dft, double t1, double t2, DFT *idft, int calc_mode, char *status, int verbose
);
int dftDerivative(
  DFT *dft, DFT *deriv, char *status
);
int dftDerivative_old(
  DFT *dft, DFT *deriv, char *status
);
/*****************************************************************************/

/*****************************************************************************/
// misc_model.c
int dftInterpolateForIMG(
  DFT *input, IMG *img, int frame_nr, DFT *output, double *ti1, double *ti2,
  int verbose, char *status
);
int imgTimeIntegral(
  IMG *img, float t1, float t2, IMG *iimg, int calc_mode, char *status, int verbose
);
int dftAllocateWithIMG(
  DFT *dft, int tacNr, IMG *img
);
int sif2dft(
  SIF *sif, DFT *dft
);
int sifAllocateWithIMG(
  SIF *sif, IMG *img, int doCounts, int verbose
);
/*****************************************************************************/

/*****************************************************************************/
// mrl.c
int mrl_between_tacs(double *y1, double *y2, int n);
/*****************************************************************************/

/*****************************************************************************/
// noise.c
int noiseSD4Simulation(
  double y, double t1, double dt, double hl, double a,
  double *sd, char *status, int verbose
);
int noiseSD4SimulationFromDFT(
  DFT *dft, int index, double pc, double *sd, char *status, int verbose
);
/*****************************************************************************/

/*****************************************************************************/
// plotdata.c
int plot_svg(
  DFT *dft, RES *res, int first, int last, char *main_title,
  char *x_title, char *y_title, char *fname, int verbose
);
int plotdata(
  DFT *dft, RES *res, int first, int last,
  char *mtitle, char *xtitle, char *ytitle, char *fname
);
int plotdata_as_dft(
  DFT *dft, char *fname
);
/*****************************************************************************/

/*****************************************************************************/
// plotfit.c
int plot_fit_svg(
  DFT *dft1, DFT *dft2, char *main_title, char *fname, int verbose
);
int plot_fitrange_svg(
  DFT *dft1, DFT *dft2, char *main_title,
  double x1, double x2, double y1, double y2, char *fname, int verbose
);
/*****************************************************************************/

/*****************************************************************************/
// units_check.c
int cunit_check_dft_vs_img(
  DFT *dft, IMG *img, char *errmsg, int verbose
);
/*****************************************************************************/

/*****************************************************************************/
// weight_model.c
int dftWeightByFreq(DFT *dft);
int imgSetWeights(IMG *img, int wmet, int verbose);
int dftWSampleNr(DFT *tac);
/*****************************************************************************/

/*****************************************************************************/
// cluster_tac.c
int clusterTACs(
  IMG *dimg, IMG *cimg, int nr, DFT *tac, int verbose
);
/*****************************************************************************/

/*****************************************************************************/
#endif
