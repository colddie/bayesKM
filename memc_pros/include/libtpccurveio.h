// @file libtpccurveio.h
/// @brief Header file for libtpccurveio.
/// @author Vesa Oikonen
///

#ifdef __cplusplus
extern "C" {
#endif


#ifndef _LIBTPCCURVEIO_H_
#define _LIBTPCCURVEIO_H_
/*****************************************************************************/

/*****************************************************************************/
#include "tpcclibConfig.h"
/*****************************************************************************/
#include "libtpcmisc.h"
/*****************************************************************************/
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <ctype.h>
#include <math.h>
#include <time.h>
/*****************************************************************************/

/*****************************************************************************/
#ifndef BACKUP_EXTENSION
/** Backup file extension */
#define BACKUP_EXTENSION ".bak"
#endif 

/*****************************************************************************/

/*****************************************************************************/
/// @cond
#ifndef DFT_VER
#define DFT_VER "DFT1"
#endif
#ifndef _DFT_COMMENT_LEN
#define _DFT_COMMENT_LEN 16384
#endif
/// @endcond
/*****************************************************************************/

/** Definitions for one DFT curve */
typedef struct {
  /* Name of curve */
  /** Name of region, including hemisphere, plane etc;
      this will replace the voiname, hemisphere and place. */
  char          name[MAX_REGIONNAME_LEN+1];
  /** Name of region */
  char          voiname[MAX_REGIONSUBNAME_LEN+1];
  /** dx/sin/avg */
  char          hemisphere[MAX_REGIONSUBNAME_LEN+1];
  /** Image plane or other description */
  char          place[MAX_REGIONSUBNAME_LEN+1];
  /** Volume of region (mm x mm x mm by default) */
  double        size;
  /* Pointer to curve data and two modified curves */
  /** Pointer to original TAC */
  double       *y;
  /** Pointer to 1st modified TAC (for user) */
  double       *y2;
  /** Pointer to 2nd modified TAC (for user) */
  double       *y3;
  /** Temporary switch for outside procedures */
  char          sw;
  /** Temporary switch for outside procedures */
  char          sw2;
  /** Temporary switch for outside procedures */
  char          sw3;
  } Voi;

/** Definitions for DFT, a set of curves */
typedef struct {
  /* Number of data frames (points), and number of planes (curves) */
  /** Nr of samples (frames) in each TAC */
  int           frameNr;
  /** Nr of regional TACs */
  int           voiNr;
  /** Name of study (study number) */
  char          studynr[MAX_STUDYNR_LEN+1];
  /* Units */
  /** Unit of y values */
  char          unit[MAX_UNITS_LEN+1];
  /** Unit of x values: TUNIT_UNKNOWN, TUNIT_SEC, TUNIT_MIN, ... */
  int           timeunit;
  /* Study information */
  /** Name of radiopharmaceutical */
  char          radiopharmaceutical[32];
  /** Isotope (for example C-11, Cu-62, unknown) */
  char          isotope[8];
  /** Decay correction: DFT_DECAY_UNKNOWN, DFT_DECAY_CORRECTED,
   *                    DFT_DECAY_NOTCORRECTED */
  char          decayCorrected;
  /** Scan start date and time: YYYY-MM-DD hh:mm:ss */
  char          scanStartTime[20];
  /** Tracer injection date and time: YYYY-MM-DD hh:mm:ss */
  char          injectionTime[20];
  /* Specified frame time */
  /** Available frame times: DFT_TIME_MIDDLE, DFT_TIME_START, DFT_TIME_END,
   *                         DFT_TIME_STARTEND */
  int           timetype;
  /* Pointers to time values */
  /** Middle frame time */
  double       *x;
  /** Frame start time  */
  double       *x1;
  /** Frame end time    */
  double       *x2;
  /** Pointers to curves */
  Voi          *voi;
  /** Pointer to weight factors */
  double       *w;
  /** Variable indicating whether weights are present (0 or 1) */
  int           isweight;
  /** String for comments */
  char          comments[_DFT_COMMENT_LEN+1];
  /* Internal variables for DFT procedures */
/// @cond
  /** Internal variable: Size of allocated memory (doubles) */
  int           _dataSize;
  /** Internal variable: Pointer to memory */
  double       *_data;
  /** Internal variable: Number of allocated curves (VOIs)*/
  int           _voidataNr;
/// @endcond
  /** Internal variable: 0=plain datafile, 1=standard DFT, ...;
                         use defines DFT_FORMAT_X */
  int           _type;
  } DFT;

/*****************************************************************************/

/*****************************************************************************/
#ifndef MAX_RESPARAMS
/** Max nr of parameters */
#define MAX_RESPARAMS 100
#endif
#ifndef MAX_RESPARNAME_LEN
/** Max length of parameter names and units */
#define MAX_RESPARNAME_LEN 15
#endif

/** Definitions for one RES region */
typedef struct {
  /* Name of the curve */
  /** Name of region, including hemisphere, plane etc;
      this will some day replace the voiname, hemisphere and place. */
  char      name[MAX_REGIONNAME_LEN+1];
  /** Name of region */
  char      voiname[MAX_REGIONSUBNAME_LEN+1];
  /** dx/sin/avg */
  char      hemisphere[MAX_REGIONSUBNAME_LEN+1];
  /** Image plane or other description */
  char      place[MAX_REGIONSUBNAME_LEN+1];
  /* Parameters and their SD's and CL's*/
  /** Array of result values */
  double    parameter[MAX_RESPARAMS];
  /** Array of result SD's   */
  double    sd[MAX_RESPARAMS];
  /** Lower 95% confidence interval */
  double    cl1[MAX_RESPARAMS];
  /** Upper 95% confidence interval */
  double    cl2[MAX_RESPARAMS];
  /** Temporary switch for user */
  int       sw;
  /** Temporary switch for user */
  int       sw2;
} ResVOI;

/** Definitions for a set of RES regions */
typedef struct {
  /** Program that produced the results */
  char      program[1024];
  /** Calculation date and time */
  time_t    time;
  /** Number of regions */
  int       voiNr;
  /** Number of parameters, <=MAX_RESPARAMS */
  int       parNr;
  /** Name of study (study number) */
  char      studynr[MAX_STUDYNR_LEN+1];
  /* Names of original datafiles */
  /** Name of original tissue datafile */
  char      datafile[FILENAME_MAX];
  /** Name of original ref datafile */
  char      reffile[FILENAME_MAX];
  /** Name of original plasmafile */
  char      plasmafile[FILENAME_MAX];
  /** Name of second original plasmafile */
  char      plasmafile2[FILENAME_MAX];
  /** Name of original bloodfile */
  char      bloodfile[FILENAME_MAX];
  /** Name of reference region */
  char      refroi[64];
  /** Free field describing fit time range */
  char      datarange[128];
  /** Number of data values used in modelling */
  int       datanr;
  /** Free text field describing fit method */
  char      fitmethod[128];
  /* Parameters affecting the results */
  /** 0=Data was not weighted, 1=Data was weighted, -1=not known */
  int       isweight;
  /** Tissue density (g/ml) */
  double    density;
  /** Lumped Constant (unitless) */
  double    lc;
  /** Beta */
  double    beta;
  /** Plasma concentration of native substrate, e.g. glucose */
  double    concentration;
  /** Vb percentage */
  double    Vb;
  /** fA percentage (arterial volume of Vb) */
  double    fA;
  /** Extraction fraction */
  double    E;
  /** List of parameter names */
  char      parname[MAX_RESPARAMS][MAX_RESPARNAME_LEN+1];
  /** List of parameter units */
  char      parunit[MAX_RESPARAMS][MAX_RESPARNAME_LEN+1];
  /** Parameter names separated by space(s); deprecated */
  char      titleline[1024];
  /** Parameter units separated by space(s); deprecated */
  char      unitline[1024];
  /** Pointers to regional curves */
  ResVOI   *voi;
  /* Internal variables for RESULT procedures */
  /** Internal variable: Number of allocated curves (VOIs); do not change */
  int      _voidataNr;
} RES;

/*****************************************************************************/

/*****************************************************************************/
/// @cond
#define FIT_VER "FIT1"
/// @endcond

#ifndef MAX_FITPARAMS
/** Max nr of parameters in FIT */
#define MAX_FITPARAMS 100
#endif

/** FIT functions */
enum mathfuncs {
  MF_LEVEL=100, MF_LINE, MF_POL2, MF_POL3, MF_POL4, MF_POL5, MF_POL6, MF_POL7,
  MF_POL8, MF_POL9,
  MF_RATF11=211, MF_RATF21=221, MF_RATF22=222, MF_RATF32=232, MF_RATF33=233,
  MF_EXP1=301, MF_EXP2, MF_EXP3, MF_EXP4, MF_EXP5,
  MF_LUNDQVIST=321, MF_LUNDQVIST2, MF_LUNDQVIST3,
  MF_EXPBOLUSINF=331, MF_EXPBOLUSINF_RW=332, MF_MF_EXPBOLUSINF_AZ=334,
  MF_PK11195=351,
  MF_HILL=841, MF_1MHILL=842, MF_1MHILL_ADE=843, MF_HILL_B=844, MF_AMHILL=845,
  MF_EHILL_PAR=846, MF_EHILL_MET=847,
  MF_EHILL2_PAR=848, MF_EHILL2_MET=849,
  MF_MAMEDE=851, MF_1MMAMEDE,
  MF_MAYER_PAR=861, MF_MAYER_MET, MF_EMAYER_PAR, MF_EMAYER_MET,
  MF_HILL3M_PAR=871, MF_HILL3M_M1, MF_HILL3M_M2, MF_HILL3M_M3, 
  MF_PF3M_PAR=881, MF_PF3M_M1, MF_PF3M_M2, MF_PF3M_M3, 
  MF_RATF33D=1232,
  MF_FENGM2=1313, MF_FENGM2E=1314,
  MF_GAMMAV=1401, MF_GAMMAVB=1402, MF_GAMMAVR=1403,
  MF_WEIBULLCDF_D=1421, MF_WEIBULLCDF_DD=1423,
  MF_SURGE=1431, MF_SURGE_TRAD=1432, MF_SURGE_RECIRC=1433, MF_P2B_SRC=1434,
  MF_HILL_D=1801, MF_HILL_DD=1811, MF_HILL_SDD=1821,
  MF_IMGPROFILE=2111,
  MF_GRAHAM_INP=9501, MF_GRAHAM_EINP, MF_GRAHAM_INPM,
  MF_HUANG_MET=9601, MF_CARSON_EMET, MF_NEW_MET,
  MF_MLMCM=9701,
};

/** Definitions for one curve */
typedef struct {
  /** Name of the curve */
  /** Name of region, including hemisphere, plane etc;
      this will replace the voiname, hemisphere and place. */
  char      name[MAX_REGIONNAME_LEN+1];
  /** Name of region */
  char      voiname[MAX_REGIONSUBNAME_LEN+1];
  /** dx/sin/avg */
  char      hemisphere[MAX_REGIONSUBNAME_LEN+1];
  /** Image plane or other description */
  char      place[MAX_REGIONSUBNAME_LEN+1];
  /** Number (type) of function */
  int       type;
  /** The number of parameters */
  int       parNr;
  /** Fit start time */
  double    start;
  /** Fit end time */
  double    end;
  /** Number of data points in the fit */
  int       dataNr;
  /** Fitted parameters */
  double    p[MAX_FITPARAMS];
  /** (weighted) sum-of-squares */
  double    wss;
  /** Temporary switch for outside procedures */
  char      sw;
  /** Temporary switch for outside procedures */
  char      sw2;
  /** Temporary switch for outside procedures */
  char      sw3;
} FitVOI;

/** Definitions for a set of curves */
typedef struct {
  /** Number of regions */
  int       voiNr;
  /** Name of original datafile */
  char      datafile[FILENAME_MAX];
  /** Name of study (study number) */
  char      studynr[MAX_STUDYNR_LEN+1];
  /** Unit of concentration */
  char      unit[MAX_UNITS_LEN+1];
  /** Time unit: TUNIT_UNKNOWN, TUNIT_SEC, TUNIT_MIN, ... */
  int       timeunit;
  /** Pointers to regional curves */
  FitVOI   *voi;
  /** Fit date and time */
  time_t    time;
  /** Program name */
  char      program[1024];
  /** Internal variables: Number of allocated curves (VOIs) */
  int      _voidataNr; /*  */
} FIT;

/*****************************************************************************/

/*****************************************************************************/
/* cpt */

/** Verbose prints from CPT functions */
int CPT_TEST;

/** Error message from CPT functions */
char cpterrmsg[128];

//void libcpt_printdate(FILE *fp);
int cptrnameSplit(
  char *rname, char *name1, char *name2, char *name3, int max_name_len);
int cptReadOne(char *cptfile, DFT *dft);
int cptWrite(DFT *dft, char *filename, int cpt_format);
/*****************************************************************************/

/*****************************************************************************/
/* csv */

/** CSV struct status */
enum {CSV_OK, CSV_ERROR, CSV_CANNOTOPEN, CSV_INVALIDFORMAT, CSV_TOOBIG,
    CSV_OUTOFMEMORY, CSV_NOTABLE
    };

/** Verbose prints from CSV functions */
int CSV_TEST;

/** CSV item */
typedef struct {
  /** CSV row, 1..row_nr */
  int row;
  /** CSV column, 1..col_nr */
  int col;
  /** CSV item content as string */
  char *content;
} CSV_item;

/** CSV struct */
typedef struct {
  /** Pointer to CSV item */
  CSV_item *c;
  /** Nr of CSV items */
  int nr;
  /** Nr of rows in CSV */
  int row_nr;
  /** max column number per row */
  int col_nr;
  /** Column separator character, usually ; or , */
  char separator;
} CSV;
/*****************************************************************************/
void csvInit(CSV *csv);
void csvEmpty(CSV *csv);
int csvRead(CSV *csv, char *fname);
int csv2dft(CSV *csv, DFT *dft);
int csv2dft_a(CSV *csv, DFT *dft);
int csv2dft_b(CSV *csv, DFT *dft);
int csv2dft_linkset(CSV *csv, DFT *dft);
int csv2dft_mat(CSV *csv, DFT *dft);
int csvIsRegular(CSV *csv);
char* csvCell(CSV *csv, int row, int col);
/*****************************************************************************/

/*****************************************************************************/
/* dft */

/** Error message from DFT functions */
char dfterrmsg[64];

/** Nr of decimals for concentration values */
extern int DFT_NR_OF_DECIMALS;

/** TAC file format (DFT type) */
#define DFT_FORMAT_UNKNOWN  -1
/** TAC file format (DFT type) */
#define DFT_FORMAT_PLAIN     0
/** TAC file format (DFT type) */
#define DFT_FORMAT_STANDARD  1
/** TAC file format (DFT type) */
#define DFT_FORMAT_IFT       2
/** TAC file format (DFT type) */
#define DFT_FORMAT_FIT       3
/** TAC file format (DFT type) */
#define DFT_FORMAT_NCI       4
/** TAC file format (DFT type) */
#define DFT_FORMAT_PMOD      5
/** TAC file format (DFT type) */
#define DFT_FORMAT_CSV_INT   6
/** TAC file format (DFT type) */
#define DFT_FORMAT_CSV_UK    7
/** TAC file format (DFT type) */
#define DFT_FORMAT_CPT       8
/** TAC file format (DFT type) */
#define DFT_FORMAT_IDWC      9
/** TAC file format (DFT type) */
#define DFT_FORMAT_IF        10
/** TAC file format (DFT type) */
#define DFT_FORMAT_XML       11
/** TAC file format (DFT type) */
#define DFT_FORMAT_HTML      12
/** TAC file format (DFT type) */
#define DFT_FORMAT_XELERIS   13

/** Definition for DFT (frame) time type */
#define DFT_TIME_MIDDLE      0
/** Definition for DFT (frame) time type */
#define DFT_TIME_START       1
/** Definition for DFT (frame) time type */
#define DFT_TIME_END         2
/** Definition for DFT (frame) time type */
#define DFT_TIME_STARTEND    3

/** Definition for Decay correction status in DFT struct */
#define DFT_DECAY_UNKNOWN      0
/** Definition for Decay correction status in DFT struct */
#define DFT_DECAY_CORRECTED    1
/** Definition for Decay correction status in DFT struct */
#define DFT_DECAY_NOTCORRECTED 2
/*****************************************************************************/
void dftInit(DFT *data);
void dftEmpty(DFT *data);
int dftSetmem(DFT *data, int frameNr, int voiNr);
int dftAddmem(DFT *data, int voiNr);
int dftAdd(DFT *data1, DFT *data2, int voi);
int dftSelect(DFT *data, char *name);
int dftSelectRegions(DFT *dft, char *region_name, int reset);
int dftSelectBestReference(DFT *dft);
void dftFrametimes(DFT *data);
int dftOverflow(DFT *data);
int dftCopyvoi(DFT *data, int from, int to);
int dftMovevoi(DFT *dft, int from, int to);
int dftDelete(DFT *dft, int voi);
int dftCopymainhdr(DFT *dft1, DFT *dft2);
int dftCopymainhdr2(DFT *dft1, DFT *dft2, int ow);
int dftCopyvoihdr(DFT *dft1, int from, DFT *dft2, int to);
int dftdup(DFT *dft1, DFT *dft2);
int dftAllocateWithHeader(DFT *dft, int frameNr, int voiNr, DFT *dft_from);
int dftAddnullframe(DFT *data);
int dftSort(DFT *data);
int dftSortPlane(DFT *data);
int dft_nr_of_NA(DFT *dft);
int dftNAfill(DFT *dft);
int dftMinMax(DFT *dft, double *minx, double *maxx, double *miny, double *maxy);
int dftMinMaxTAC(
  DFT *dft, int tacindex, double *minx, double *maxx,
  double *miny, double *maxy, int *mini, int *maxi, int *mins, int *maxs
);
int dftMaxY(DFT *dft, double t1, double t2, double *miny, double *maxy);
double dft_kBqMin(DFT *data);
double dft_kBqMax(DFT *data);
int dftSortByFrame(DFT *dft);
int dftDeleteFrameOverlap(DFT *dft);
int dftDeleteFrameOverlap_old(DFT *dft);
int dftRemoveTimeRange(DFT *dft, double startT, double endT);
void dftSetComments(DFT *dft);
int dftFillInitialGap(DFT *dft);
int dftAddSpaceForFrames(DFT *dft, int nr_to_add);
void dftRNameSimplify(DFT *dft, int hemisphere, int place);
int dftMeanTAC(DFT *dft, DFT *mean);
int dftValidNr(DFT *dft, double tstart, double tstop, int index);

/* Deprecated functions. Please don't use these anymore */
/// @cond
#define initDFT dftInit
#define emptyDFT dftEmpty
#define setmemDFT dftSetmem
#define addmemDFT dftAddmem
#define addDFT dftAdd
#define selectDFT dftSelect
#define frametimesDFT dftFrametimes
#define overflowDFT dftOverflow
#define copyvoiDFT dftCopyvoi
#define copymainhdrDFT dftCopymainhdr
#define copyvoihdrDFT dftCopyvoihdr
#define addnullframeDFT dftAddnullframe
#define sortDFT dftSort
#define na_fillDFT dftNAfill
/// @endcond
/*****************************************************************************/

/*****************************************************************************/
/* dftdecayc */
int dftDecayCorrection(
  DFT *dft, double hl, int mode, int y, int y2, int y3,
  char *status, int verbose
);
/*****************************************************************************/

/*****************************************************************************/
/* dftio */
int dftRead(char *filename, DFT *data);
int dftWrite(DFT *data, char *filename);
void dftPrint(DFT *data);
int dftFormat(char *fname);
int dftType(FILE *fp);
int dftWriteHTML(DFT *dft, char *fname, int orientation);
int dftWriteXHTML11_doctype(FILE *fp);
int dftWriteXHTML11_head(FILE *fp, char *author_name);
int dft_fill_hdr_from_IFT(DFT *dft, IFT *ift);
int dftGetPmodTitle(DFT *dft, char *title_line);
/*****************************************************************************/

/*****************************************************************************/
/* dftres */
int res_allocate_with_dft(RES *res, DFT *dft);
int dftToResult(DFT *dft, RES *res, char *status);
/*****************************************************************************/

/*****************************************************************************/
/* dftunit */

void dftUnitToDFT(DFT *dft, int dunit);
int dftUnitConversion(DFT *dft, int dunit);
int dftTimeunitToDFT(DFT *dft, const char *timeunit);
int dftTimeunitConversion(DFT *dft, int tunit);
void dftMin2sec(DFT *data);
void dftSec2min(DFT *data);

/// @cond
/* Deprecated functions. Please don't use these anymore */
#define dftUnitId petCunitId
#define dftUnit petCunit
#define dftTimeunitId petTunitId
#define dftTimeunit petTunit
/* DFTUNITs are deprecated unit definitions. Please don't use these anymore */
#define DFTUNIT_UNKNOWN CUNIT_UNKNOWN
#define DFTUNIT_CPS CUNIT_CPS
#define DFTUNIT_COUNTS CUNIT_COUNTS
#define DFTUNIT_KBQ_PER_ML CUNIT_KBQ_PER_ML
#define DFTUNIT_SEC_KBQ_PER_ML CUNIT_SEC_KBQ_PER_ML
#define DFTUNIT_PER_SEC CUNIT_PER_SEC
#define DFTUNIT_PER_MIN CUNIT_PER_MIN
#define DFTUNIT_ML_PER_ML CUNIT_ML_PER_ML
#define DFTUNIT_ML_PER_DL CUNIT_ML_PER_DL
#define DFTUNIT_ML_PER_ML_PER_MIN CUNIT_ML_PER_ML_PER_MIN
#define DFTUNIT_ML_PER_DL_PER_MIN CUNIT_ML_PER_DL_PER_MIN
#define DFTUNIT_UNITLESS CUNIT_UNITLESS
#define DFTUNIT_NCI_PER_ML CUNIT_NCI_PER_ML
#define DFTUNIT_MBQ_PER_ML CUNIT_MBQ_PER_ML
#define DFTUNIT_BQ_PER_ML CUNIT_BQ_PER_ML
#define DFTUNIT_UCI_PER_ML CUNIT_UCI_PER_ML
#define DFTUNIT_UMOL_PER_MIN_PER_100G CUNIT_UMOL_PER_MIN_PER_100G
#define DFTUNIT_MG_PER_MIN_PER_100G CUNIT_MG_PER_MIN_PER_100G
/* DFTTIMEs are deprecated unit definitions. Please don't use these anymore */
#define DFTTIME_UNKNOWN TUNIT_UNKNOWN
#define DFTTIME_SEC TUNIT_SEC
#define DFTTIME_MIN TUNIT_MIN
#define DFTTIME_UM TUNIT_UM
#define DFTTIME_MM TUNIT_MM
/// @endcond
/*****************************************************************************/

/*****************************************************************************/
/* fitres */
int fit_allocate_with_dft(FIT *fit, DFT *dft);
int fitToResult(FIT *fit, RES *res, char *status);
/*****************************************************************************/

/*****************************************************************************/
/* idwc */
int idwcWrite(DFT *dft, char *filename);
int idwcRead(char *filename, DFT *dft);
/*****************************************************************************/

/*****************************************************************************/
/* if */
int ifWrite(DFT *dft, char *filename);
int ifRead(char *filename, DFT *dft);
/*****************************************************************************/

/*****************************************************************************/
/* mathfunc */

/** Verbose prints from FIT functions */
int MATHFUNC_TEST;

/** Error message from FIT functions */
char fiterrmsg[64];
/*****************************************************************************/
void fitEmpty(FIT *fit);
void fitInit(FIT *fit);
int fitSetmem(FIT *fit, int voiNr);
void fitPrint(FIT *fit);
int fitWrite(FIT *fit, char *filename);
int fitRead(char *filename, FIT *fit, int verbose);
int fitEval(FitVOI *r, double x, double *y);
int fitEvaltac(FitVOI *r, double *x, double *y, int dataNr);
int fitFunctionname(int type, char *str);
int fitFunctionformat(int type, char *str);
int fitIntegralEval(FitVOI *r, double x, double *yi);
int fitIntegralEvaltac(FitVOI *r, double *x, double *yi, int dataNr);
int fitDerivEval(FitVOI *r, double x, double *yd);
int fitDerivEvaltac(FitVOI *r, double *x, double *yd, int dataNr);

/// @cond
/* Deprecated functions. Please don't use these anymore */
#define emptyFIT fitEmpty
#define initFIT fitInit
#define setmemFIT fitSetmem
#define printFIT fitPrint
#define writeFIT fitWrite
#define readFIT fitRead
#define evalFIT fitEval
#define evaltacFIT fitEvaltac
#define functionnameFIT fitFunctionname
#define functionformatFIT fitFunctionformat
#define ievalFIT fitIntegralEval
#define ievaltacFIT fitIntegralEvaltac
#define devalFIT fitDerivEval
#define devaltacFIT fitDerivEvaltac
/// @endcond
/*****************************************************************************/

/*****************************************************************************/
/* ncifile */
int roikbqWrite(DFT *dft, char *fname);
int roikbqRead(char *fname, DFT *dft);
/*****************************************************************************/

/*****************************************************************************/
/* resift */
int res2ift(RES *res, IFT *ift, int verbose);
/*****************************************************************************/

/*****************************************************************************/
/* result */

/** Verbose prints from RES functions */
int RESULT_TEST;

/** Error message from RES functions */
char reserrmsg[64];

void resEmpty(RES *res);
void resInit(RES *res);
int resSetmem(RES *res, int voiNr);
void resFixParnames(RES *res);
void resPrint(RES *res);
int resRead(char *filename, RES *res, int verbose);
int resWrite(RES *res, char *filename, int verbose);
int resWriteHTML(RES *res, char *fname, int verbose);
int resFName2study(char *fname, char *studyNumber);
int resMedian(double *data, int nr, double *median, double *min, double *max);
int resMean(double *data, int nr, double *mean, double *sd);
void resSortByName(RES *res);
int resCopyMHeader(RES *res1, RES *res2);
int resDelete(RES *res, int voi);
int resSelect(RES *data, char *name);
int resSelectRegions(RES *res, char *region_name, int reset);
int resParameterPrintType(RES *res, int parIndex);
int resWriteXHTML11_doctype(FILE *fp);
int resWriteXHTML11_head(FILE *fp, char *author_name);
int resWriteHTML_table(RES *res, FILE *fp);
int resIsDuplicateNames(RES *res);
int resMatchHeader(RES *res1, RES *res2);
int resMatchRegions(RES *res1, RES *res2);
int resMatchParameternames(RES *res1, RES *res2);
int resMatchParameters(
  RES *res1, RES *res2, int test_par, double test_limit, int test_sd);
int resMatchParametersAbs(
  RES *res1, RES *res2, int test_par, double test_limit, int test_sd);
int resRNameSubfieldExists(RES *res);
/*****************************************************************************/

/*****************************************************************************/
/* tsv.c */
int tsvRead(char *filename, DFT *dft);
/*****************************************************************************/

/*****************************************************************************/
/* xeleris.c */
int xelRead(char *filename, DFT *dft);
/*****************************************************************************/

/*****************************************************************************/
#endif /* LIBTPCCURVEIO */


#ifdef __cplusplus
}
#endif
