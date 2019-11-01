/// @file libtpcmisc.h
/// @brief Header file for libtpcmisc.
/// @author Vesa Oikonen
///
#ifndef _LIBTPCMISC_H_
#define _LIBTPCMISC_H_
/*****************************************************************************/

/*****************************************************************************/
#include "tpcclibConfig.h"
/*****************************************************************************/
#include <errno.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <ctype.h>
#include <stdint.h>
#include <unistd.h>
/*****************************************************************************/

/*****************************************************************************/
/** Isotope branching ratio */
#define BRANCHING_O 0.999
/** Isotope branching ratio */
#define BRANCHING_C 0.998
/** Isotope branching ratio */
#define BRANCHING_Cu64 0.174
/** Isotope branching ratio */
#define BRANCHING_N 0.998
/** Isotope branching ratio */
#define BRANCHING_F 0.967
/** Isotope branching ratio */
#define BRANCHING_Ge 0.891
/** Isotope branching ratio */
#define BRANCHING_Ga 0.891
/** Isotope branching ratio */
#define BRANCHING_Rb 0.950
/*****************************************************************************/
/** Isotope halflife in minutes */
#define HL_O15 2.05 /* 123 s */
/** Isotope halflife in minutes */
#define HL_N13 10.0
/** Isotope halflife in minutes */
#define HL_C11 20.4
/** Isotope halflife in minutes */
#define HL_F18 109.8
/** Isotope halflife in minutes */
#define HL_Ge68 396000.0 /* 275 d */
/** Isotope halflife in minutes */
#define HL_Ga68 68.0
/*****************************************************************************/
/* The following halflifes have not been checked from the reference;         */
/* they are thus meant to be used only during program development period     */
/** Isotope halflife in minutes; not verified from the reference */
#define HL_Br75 98.0
/** Isotope halflife in minutes; not verified from the reference */
#define HL_Br76 978.33 /* 58700 s */
/** Isotope halflife in minutes; not verified from the reference */
#define HL_Cu62 9.7    /* 582 s */
/** Isotope halflife in minutes; not verified from the reference */
#define HL_Cu64 762.018 /* 12.7003 h */
/** Isotope halflife in minutes; not verified from the reference */
#define HL_Fe52 4980.0
/** Isotope halflife in minutes; not verified from the reference */
#define HL_Na22 1368000.0
/** Isotope halflife in minutes; not verified from the reference */
#define HL_O14 1.1818
/** Isotope halflife in minutes; not verified from the reference */
#define HL_Rb82 1.25   /* 75 s */
/** Isotope halflife in minutes; not verified from the reference */
#define HL_Zn62 558.0
/** Isotope halflife in minutes; not verified from the reference */
#define HL_I124 6013.44 /* 4.176 */
/*****************************************************************************/
/** isotope_code */
enum {
  TPCISOT_BR_75, TPCISOT_BR_76, TPCISOT_CU_62, TPCISOT_CU_64, 
  TPCISOT_FE_52, TPCISOT_GA_68, TPCISOT_GE_68, TPCISOT_NA_22,
  TPCISOT_RB_82, TPCISOT_ZN_62, TPCISOT_F_18, TPCISOT_C_11,
  TPCISOT_N_13, TPCISOT_O_15, TPCISOT_O_14, TPCISOT_I_124,
  TPCISOT_UNKNOWN
};
/*****************************************************************************/
/* Add ln(2) if it is not defined */
#ifndef M_LN2
/** ln(2) */
#define M_LN2       0.69314718055994530942
#endif
/*****************************************************************************/
#ifndef MAX_UNITS_LEN
/** Max length of units string (+1), based on ECAT7 format */
#define MAX_UNITS_LEN 31
#endif
/*****************************************************************************/
/** Data y units */
enum {
  /*  0 */ CUNIT_UNKNOWN,
  /*  1 */ CUNIT_CPS,
  /*  2 */ CUNIT_COUNTS,
  /*  3 */ CUNIT_KBQ_PER_ML,
  /*  4 */ CUNIT_SEC_KBQ_PER_ML,
  /*  5 */ CUNIT_PER_SEC,
  /*  6 */ CUNIT_PER_MIN,
  /*  7 */ CUNIT_ML_PER_ML,
  /*  8 */ CUNIT_ML_PER_DL,
  /*  9 */ CUNIT_ML_PER_ML_PER_MIN,
  /* 10 */ CUNIT_ML_PER_DL_PER_MIN,
  /* 11 */ CUNIT_UNITLESS,
  /* 12 */ CUNIT_NCI_PER_ML,
  /* 13 */ CUNIT_MBQ_PER_ML,
  /* 14 */ CUNIT_BQ_PER_ML,
  /* 15 */ CUNIT_UCI_PER_ML,
  /* 16 */ CUNIT_UMOL_PER_MIN_PER_100G,
  /* 17 */ CUNIT_MG_PER_MIN_PER_100G,
  /* 18 */ CUNIT_UMOL_PER_MIN_PER_DL,
  /* 19 */ CUNIT_MG_PER_MIN_PER_DL,
  /* 20 */ CUNIT_PERCENTAGE,
  /* 21 */ CUNIT_KCPS,
  /* 22 */ CUNIT_MIN_KBQ_PER_ML,
  /* 23 */ CUNIT_BQ,
  /* 24 */ CUNIT_KBQ,
  /* 25 */ CUNIT_MBQ,
  /* 26 */ CUNIT_GBQ,
  /* 27 */ CUNIT_NCI,
  /* 28 */ CUNIT_UCI,
  /* 29 */ CUNIT_MCI,
  /* 30 */ CUNIT_PID,
  /* 31 */ CUNIT_PIDM,
  /* 32 */ CUNIT_PIDV,
  /* 33 */ CUNIT_G_PER_ML, // SUV unit
  /* 34 */ CUNIT_ML_PER_G // SUV unit
};
/** Data x units */
enum {
  /*  0 */ TUNIT_UNKNOWN,
  /*  1 */ TUNIT_SEC,
  /*  2 */ TUNIT_MIN,
  /*  3 */ TUNIT_UM,
  /*  4 */ TUNIT_MM,
  /*  5 */ TUNIT_CM,
  /*  6 */ TUNIT_M,
  /*  7 */ TUNIT_HOUR,
  /*  8 */ TUNIT_MONTH,
  /*  9 */ TUNIT_YEAR,
  /* 10 */ TUNIT_MSEC
};
/*****************************************************************************/
#ifndef MAX_REGIONNAME_LEN
/** Max length of Region name (+1) */
#define MAX_REGIONNAME_LEN 20
#endif
#ifndef MAX_REGIONSUBNAME_LEN
/** Max length of Region name subfield (+1) */
#define MAX_REGIONSUBNAME_LEN 6
#endif
/*****************************************************************************/
#ifndef MAX_STUDYNR_LEN
/** Max length of Study number (+1) */
#define MAX_STUDYNR_LEN 255
#endif
/*****************************************************************************/

/*****************************************************************************/
/* backup */
int backupExistingFile(char *filename, char *backup_ext, char *status);
int fileCopy(char *filename1, char *filename2, char *status);
/*****************************************************************************/

/*****************************************************************************/
/* branch */
float branchingFraction(int isotope);
/*****************************************************************************/

/*****************************************************************************/
/* datetime */
#ifndef HAVE_GMTIME_R
struct tm* gmtime_r(const time_t* t, struct tm* tm);
#endif
#ifndef HAVE_LOCALTIME_R
struct tm* localtime_r(const time_t* t, struct tm* tm);
#endif
#ifndef HAVE_TIMEGM
time_t timegm(struct tm *tm);
#else
extern time_t timegm(struct tm *tm); // needed at least in OSX
#endif
char* ctime_r_int(const time_t *t, char *buf);
int isdate(char *str);
int isdate2(char *str, char *intdate);
int isdate3(char *str, char *intdate);
int isdate4(int dateint, int *year, int *month, int *day);
int istime(char *str);
int isdatetime(char *str, char *intdate);
int get_datetime(char *str, struct tm *date, int verbose);
int get_date(char *str, struct tm *date);
long int math_div(long int a, long int b);
int isleapyear(long int year);
long int leaps_between(long int year1, long int year2);
void time_to_tm(time_t totalsecs, int offset, struct tm *result);
double tmDifference(struct tm *tm1, struct tm *tm0);
void tmAdd(int s, struct tm *d);
/*****************************************************************************/

/*****************************************************************************/
/* decpoint */
int dec_comma_is(char *str);
int dec_separator(char *str);
void dec_separator_change(char *str, int decsep);
double atof_dpi(char *str);
int dec_nr(char *str);
int atof_with_check(char *double_as_string, double *result_value);
char *strPtrToNextValue(char *str, char **nxtp);
int atoi_with_check(const char *int_as_string, int *result_value);
/*****************************************************************************/

/*****************************************************************************/
/* filename */
void filenameRmPath(char *s);
int filenameRmExtension(char *s);
void filenameRmExtensions(char *s);
int fnmatch(const char *fname, const char *key);
int fncasematch(const char *fname, const char *key);
char *filenameGetExtension(char *s);
char *filenameGetExtensions(char *s);
/*****************************************************************************/

/*****************************************************************************/
/* halflife */
char *hlIsotopeCode(int isotope);
double hlFromIsotope(char *isocode);
double hl2lambda(double halflife);
double hlLambda2factor(double lambda, double frametime, double framedur);
float hlLambda2factor_float(float lambda, float frametime, float framedur);
char *hlCorrectIsotopeCode(char *isocode);
int hlIsotopeFromHalflife(double halflife);
/// @cond
/* Deprecated function names. Please don't use these anymore */
#define lambda2factor hlLambda2factor
#define lambda2factor_float hlLambda2factor_float
/// @endcond
/*****************************************************************************/

/*****************************************************************************/
/* idcrypt (deprecated) */
const char *libpet_idcrypt_version(void);
int id_crypt(const char *string,const char *key,char *out,int decrypt);
/*****************************************************************************/

/*****************************************************************************/
/* ift */
/** Definitions for ift status message */
enum {IFT_OK, IFT_FAULT, IFT_NOMEMORY, IFT_CANNOTREAD, IFT_CANNOTWRITE,
      IFT_UNKNOWNFORMAT, IFT_KEYNOTFOUND, IFT_NODATA, IFT_VALUENOTFOUND};
/*****************************************************************************/
/** IFT item */
typedef struct {
  /** Key (comment) type character: space, #, ; */
  char type;
  /** Undefined short integer for the use of programmer */
  short int sw;
  /** Pointer to the NULL terminated key string; NULL if not allocated */
  char *key;
  /** Pointer to the NULL terminated key value string; NULL if not allocated */
  char *value;
} IFT_KEY_AND_VALUE;
/** IFT struct */
typedef struct {
/// @cond
  /** Number of allocated places for keys and values */
  int _memNr;
/// @endcond
  /** Number of stored keys and their values */
  int keyNr;
  /** Type of the parameter file:
      0=unknown, 1=interfile ':=' , 2=setup '=', 3=result ':', 4=space ' ' 
      5=tab, 6=',', 7=';'
  */
  int type;
  /** Pointer to a status message */
  const char *status;
  /** List of key-value -pairs */
  IFT_KEY_AND_VALUE *item;
  /** Size of binary data (in bytes); not yet supported */
  size_t datasize;
  /** Pointer to binary data; not yet supported */
  unsigned char *data;
} IFT;
/*****************************************************************************/
/** Verbose prints from IFT functions */
int IFT_TEST;
/*****************************************************************************/
//void libift_printdate(FILE *fp);
void iftSetStatus(IFT *ift, int status);
void iftInit(IFT *ift);
void iftEmpty(IFT *ift);
int iftPut(IFT *ift, char *key, char *value, char *cmt_type);
int iftPutDouble(IFT *ift, char *key, double value, char *cmt_type);
int iftRead(IFT *ift, char *filename, int is_key_required);
char *iftReadValue(char *filename, char *keystr);
int iftWriteItem(IFT *ift, int item, FILE *fp);
int iftWrite(IFT *ift, char *filename);
int defRead(IFT *ift, char *filename);
int iftGet(IFT *ift, char *key);
int iftGetNth(IFT *ift, char *key, int n);
int iftFindNthKey(IFT *ift, char *str, int n);
int iftFindNthValue(IFT *ift, char *str, int n);
int iftGetKeyNr(IFT *ift, const char *key);
int iftGetFrom(IFT *ift, int si, const char *key);
int iftGetFullmatchFrom(IFT *ift, int si, const char *key, const char *value);
int iftGetFloatValue(IFT *ift, int si, const char *key, float *value);
int iftGetDoubleValue(IFT *ift, int si, const char *key, double *value);
int iftGetIntValue(IFT *ift, int si, const char *key, int *value);
int iftDeleteItem(IFT *ift, int item);
int iftReplaceNthValue(IFT *ift, int item, char *value);
int iftdup(IFT *ift1, IFT *ift2);
/*****************************************************************************/

/*****************************************************************************/
/* substitutions */
#ifndef HAVE_STRDUP
char* strdup(const char* s);
#endif
///// #ifndef HAVE_STRCASESTR
// extern char *strcasestr(const char *haystack, const char *needle);    /// turn off because the conflict with armadillo!
///// #endif
/*****************************************************************************/

/*****************************************************************************/
/* intex */
/** Integer list for functions intExpand() and intMerge(). 
 *  Deprecated, should be replaced by INTEGER_LIST. */
typedef struct {
  /** Nr of integers */
  int nr;
  /** List of integers */
  int *i;
} INT_list;

/** Integer list for functions integerListInit(), integerListEmpty(), 
 *  integerListAdd(), and integerListSort(). */
typedef struct {
  /** Nr of integers */
  int nr;
  /** Allocated list size */
  int _allocNr;
  /** List of integers */
  int *list;
} INTEGER_LIST;
/*****************************************************************************/
/* Deprecated functions. */
void intInit(INT_list *list);
void intEmpty(INT_list *list);
int intExpand(char *text, INT_list *list);
INT_list intMerge(INT_list *list1,INT_list *list2);
int _intexadd(INT_list *list, int a);
/* Recommended functions */
int integerListInit(INTEGER_LIST *l);
int integerListEmpty(INTEGER_LIST *l);
int integerListAdd(INTEGER_LIST *l, int v, int ifnew);
int integerListSort(INTEGER_LIST *l);
int integerListAddFromString(
  const char *s1, const char *s2, INTEGER_LIST *l, const int ifnew
);
int integerListExpandFromString(
  const char *s1, const char *s2, INTEGER_LIST *l, const int ifnew
);
/*****************************************************************************/

/*****************************************************************************/
/* petc99 */
int temp_roundf(float e);
/*****************************************************************************/

/*****************************************************************************/
/* petunits */
int petCunitId(const char *unit);
int petTunitId(const char *timeunit);
char *petCunit(int cunit);
char *petTunit(int tunit);
int cunitFromFilename(char *fname);
/*****************************************************************************/

/*****************************************************************************/
/* proginfo */
int tpcProcessStdOptions(
  const char *s, int *print_usage, int *print_version, int *verbose_level);
void tpcProgramName(const char *program, int version, int copyright, 
  char *prname, int n);
void tpcPrintUsage(const char *program, char *text[], FILE *fp);
int tpcHtmlUsage(const char *program, char *text[], const char *path);
void tpcPrintBuild(const char *program, FILE *fp);
/*****************************************************************************/

/*****************************************************************************/
/* quots */
char *strstr_noquotation(const char *str1, const char *str2);
int strnCopyClean(char *str1, const char *str2, int maxlen);
/*****************************************************************************/

/*****************************************************************************/
/* readfile */
/** STR_TOKEN_LIST struct for functions in readfile.c */
typedef struct {
  /** Number of available string tokens */
  int token_nr;
  /** Number of allocated list items */
  int list_size;
  /** List of string tokens */
  char **tok;  
} STR_TOKEN_LIST;

void str_token_list_init(STR_TOKEN_LIST *lst);
void str_token_list_empty(STR_TOKEN_LIST *lst);
int str_token_list_add(STR_TOKEN_LIST *lst, char *new_item);
int str_token_list_del(STR_TOKEN_LIST *lst, int item);
int str_token_list_read(const char *filename, STR_TOKEN_LIST *lst);
int textfileReadLines(const char *filename, STR_TOKEN_LIST *lst);
int readStrtokens(const char *filename, char ***toklist);
int asciiCommentLine(const char *line, int *cont);
/*****************************************************************************/

/*****************************************************************************/
/* rname */
int rnameSplit(
  char *rname, char *name1, char *name2, char *name3, int max_name_len
);
int rnameMatch(char *rname, int rnr, char *test_str);
int rnameRmDots(char *rname1, char *rname2);
int rnameCatenate(
  char *rname, int max_rname_len, char *name1, char *name2, char *name3, 
  char space
);
int roinameExists(char *roiname);
/*****************************************************************************/

/*****************************************************************************/
/* strext */
int strTokenNr(const char *str1, const char *str2);
int strTokenNCpy(
  const char *str1, const char *str2, int i, char *str3, int count
);
char *strTokenDup(const char *s1, const char *s2, int *next);
int strChrCount(const char *str1, const char *str2);
void strReplaceChar(char *str, char c1, char c2);
#ifndef HAVE_STRNLEN
size_t strnlen(const char *s, size_t n);
#else
extern size_t strnlen(const char *s, size_t n);
#endif
#ifndef HAVE_STRLCAT
size_t strlcat(char *dst, const char *src, size_t dstsize);
#else
extern size_t strlcat(char *dst, const char *src, size_t dstsize);
#endif
#ifndef HAVE_STRLCPY
size_t strlcpy(char *dst, const char *src, size_t dstsize);
#else
extern size_t strlcpy(char *dst, const char *src, size_t dstsize);
#endif
int strncpyCleanSpaces(char *s1, const char *s2, int maxlen);
int strCleanSpaces(char *s);
/*****************************************************************************/

/*****************************************************************************/
/* studynr */
int studynr_in_fname(char *fname, char *studynr);
int studynr_from_fname(char *fname, char *studynr);
int studynr_from_fname2(char *fname, char *studynr, int force);
int studynr_match(char *studynr1, char *studynr2);
int studynr_validity_check2(char *studynr, int zero_ok);
int studynr_validity_check(char *studynr);
int studynr_rm_zeroes(char *studynr);
int studynr_to_lowercase(char *studynr);
/*****************************************************************************/

/*****************************************************************************/
/* swap */
int little_endian();
void swap(void *orig, void *new1, int size);    // was new, to avoid C++ keyword
void swabip(void *buf, int size);
void swawbip(void *buf, int size);
void swawip(void *buf, int size);
void printf32bits(void *buf);
/*****************************************************************************/

/*****************************************************************************/
/* doubleutil */
int doubleMatch(const double v1, const double v2, const double lim);
int doubleMatchRel(const double v1, const double v2, const double lim);
double doubleMachEps();
void doubleCopy(double *t, double *s, const unsigned int n);
unsigned int doubleMaxIndex(double *a, const unsigned int n);
double doubleSum(double *a, const unsigned int n);
double doubleMean(double *a, const unsigned int n);
int doubleSpanPositives(double *a, const int n);
int doubleCSpanPositives(double *a, const int n);
void statSortDouble(double *data, unsigned int n, int order);
void statSortFloat(float *data, unsigned int n, int order);
/*****************************************************************************/

/*****************************************************************************/
#endif /* LIBTPCMISC */
