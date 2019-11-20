/** @file fitmbf.c
    @brief Fits Hidehiro Iida's myocardial [O-15]H2O one-tissue compartment 
           model to regional PET TAC data.
    @copyright (c) Turku PET Centre
    @author Vesa Oikonen
 */
/// @cond
/*****************************************************************************/
#include "tpcclibConfig.h"
/*****************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <string.h>
#include <time.h>
/*****************************************************************************/
#include "libtpcmisc.h"
#include "libtpcmodel.h"
#include "libtpccurveio.h"
#include "libtpcsvg.h"
#include "libtpcmodext.h"
/*****************************************************************************/

/*****************************************************************************/
int parNr=3;
DFT input;
double *petmeas, *petsim, *weight;
double pc=0.9464, Beta=0.91;
double pmin[MAX_PARAMETERS], pmax[MAX_PARAMETERS];
int fitframeNr;
double wss_wo_penalty=0.0;
/*****************************************************************************/
/** Definitions for model parameter indices */
enum parameter {
  CM_FLOW, CM_PTF, CM_VA, CM_RMBF, CM_WSS
};
/*****************************************************************************/
/* Local functions */
double mbfFunc(int parNr, double *p, void*);
double mbfFunc2(int parNr, double *p, void*);
/*****************************************************************************/

/*****************************************************************************/
static char *info[] = {
  "Non-linear fitting of Iida's MBF model (1, 2) as represented in (3) to",
  "regional dynamic PET [O-15]H2O study data.",
  "The model parameters are myocardial blood flow in perfusable tissue (ptMBF),",
  "perfusable tissue fraction (PTF), and arterial blood volume and spillover",
  "(Va); in addition, mean blood flow in the myocardial region (rMBF), and",
  "weighted sum-of-squares (WSS) are reported.",
  " ",
  "The same method is applied in Carimas, and for clinical work use of Carimas",
  "is recommended; however, it is possible to save regional TACs in Carimas or",
  "other software and use those with this program.",
  " ",
  "User must provide the regional TAC file (tacfile, in DFT or PMOD format),",
  "and the names or numbers of LV cavity (lvcav) and whole myocardial (myoc)",
  "TAC inside the TAC file, and filename for the results.",
  "LV cavity and whole myocardial ROI TACs are used to estimate a spill-in",
  "corrected arterial blood TAC, which is then used as model input for",
  "the smaller myocardial regions; to omit this step and use the LV cavity TAC",
  "directly as input, enter 'none' in place of the myocardial ROI name.",
  " ",
  "Usage: @P [Options] tacfile lvcav myoc resultfile",
  " ",
  "Options:",
  " -lim[=<filename>]",
  "     Specify the constraints for model parameters;",
  "     This file with default values can be created by giving this",
  "     option as the only command-line argument to this program.",
  "     Without filename the default values are printed on screen.",
  "     Parameter can be fixed to a certain value by setting its",
  "     lower and upper limit to that value.",
  " -beta=<Beta value>",
  "     Enter the Beta value (from [O-15]CO study); by default 0.91.",
  " -pH2O=<Partition coefficient for water>",
  "     Enter the partition coefficient of water; 0.9464 by default.",
  " -end=<Fit end time (sec)>",
  "     By default line is fitted to the end of data. Use this option to enter",
  "     the fit end time.",
  " -SD[=<y|N>]",
  "     Standard deviations are calculated and saved in results (Y, default),",
  "     or not calculated (n).",
  "     Program runs a lot faster if SD and CL are not calculated.",
  " -CL[=<y|N>]",
  "     95% Confidence limits are calculated and saved in results (y), or",
  "     not calculated (N, default).",
  " -input=<Filename>",
  "     Save arterial concentration curves, estimated from LV cavity and whole",
  "     myocardial TACs, into specified TAC file.",
  " -fit=<Filename>",
  "     Fitted regional TACs are written in DFT format.",
  "     Input TAC sample times are corrected by the median of fitted time",
  "     delay values and saved; resulting input file can be used with imgflow,",
  "     or as input to this program to have common time delay for all regions.",
  " -svg=<Filename>",
  "     Fitted and measured TACs are plotted in specified SVG file.",
  " -stdoptions", // List standard options like --help, -v, etc
  " ",
  "Example:",
  "     @P -beta=0.91 s2345.tac 'lv Pl06' 'whole' s2345mbf.res",
  " ",
  "References:",
  "1. Iida H, Rhodes CG, de Silva R, Yamamoto Y, Araujo LI, Maseri A, Jones T.",
  "   Myocardial tissue fraction - correction for partial volume effects and",
  "   measure of tissue viability. J Nucl Med 1991; 32:2169-2175.",
  "2. Iida H, Rhodes CG, de Silva R, Araujo LI, Bloomfield P, Lammertsma AA,",
  "   Jones T. Use of the left ventricular time-activity curve as a noninvasive",
  "   input function in dynamic oxygen-15-water positron emission tomography.",
  "   J Nucl Med 1992; 33:1669-1677.",
  "3. Oikonen V. Model equations for myocardial perfusion studies with [15O]H2O",
  "   PET. http://www.turkupetcentre.net/reports/tpcmod0005.pdf",
  " ",
  "See also: sim_mbf, b2t_h2o, simimyoc, fit_h2o, dftweigh, rescoll",
  " ",
  "Keywords: TAC, modelling, myocardium, perfusion, radiowater, 1TCM",
  0};
/*****************************************************************************/

/*****************************************************************************/
/* Turn on the globbing of the command line, since it is disabled by default in
   mingw-w64 (_dowildcard=0); in MinGW32 define _CRT_glob instead, if necessary;
   In Unix&Linux wildcard command line processing is enabled by default. */
/*
#undef _CRT_glob
#define _CRT_glob -1
*/
int _dowildcard = -1;
/*****************************************************************************/

/*****************************************************************************/
/**
 *  Main
 */
int main(int argc, char **argv)
{
  int        ai, help=0, version=0, verbose=1;
  double     def_pmin[MAX_PARAMETERS], def_pmax[MAX_PARAMETERS];
  char       tacfile[FILENAME_MAX], resfile[FILENAME_MAX], 
             fitfile[FILENAME_MAX], svgfile[FILENAME_MAX],
             lvcavname[FILENAME_MAX], myocname[FILENAME_MAX],
             limfile[FILENAME_MAX], inputfile[FILENAME_MAX];
  int        doBootstrap=0, doSD=0, doCL=0; // 0=no, 1=yes
  double     fittime=-1.0;
  int        fittedparNr=0;
  int        ret, originallyMinutes=0;


#ifdef MINGW
  // Use Unix/Linux default of two-digit exponents in MinGW on Windows
  _set_output_format(_TWO_DIGIT_EXPONENT);
#endif

  /* Set parameter initial values and constraints */
  /* ptMBF */ def_pmin[0]=0.00;   def_pmax[0]=10.0; // mL/(mL*min)
  /* PTF   */ def_pmin[1]=0.05;   def_pmax[1]=1.0;  // mL/mL
  /* Va    */ def_pmin[2]=0.05;   def_pmax[2]=0.99; // mL/mL

  /*
   *  Get arguments
   */
  if(argc==1) {tpcPrintUsage(argv[0], info, stderr); return(1);}
  tacfile[0]=resfile[0]=limfile[0]=inputfile[0]=(char)0;
  svgfile[0]=fitfile[0]=lvcavname[0]=myocname[0]=(char)0;
  /* Get options first, because may affect what arguments are read */
  for(ai=1; ai<argc; ai++) if(*argv[ai]=='-') {
    char *cptr=argv[ai]+1; if(*cptr=='-') cptr++; if(cptr==NULL) continue;
    if(tpcProcessStdOptions(argv[ai], &help, &version, &verbose)==0) continue;
    if(strncasecmp(cptr, "CL", 2)==0) {
      if(strlen(cptr)==2) {doCL=1; continue;}
      cptr+=2; if(*cptr=='=') {
        cptr++;
        if(*cptr=='Y' || *cptr=='y') {doCL=1; continue;}
        if(*cptr=='N' || *cptr=='n') {doCL=0; continue;}
      }
    } else if(strncasecmp(cptr, "SD", 2)==0) {
      if(strlen(cptr)==2) {doSD=1; continue;}
      cptr+=2; if(*cptr=='=') {
        cptr++;
        if(*cptr=='Y' || *cptr=='y') {doSD=1; continue;}
        if(*cptr=='N' || *cptr=='n') {doSD=0; continue;}
      }
    } else if(strncasecmp(cptr, "LIM=", 4)==0 && strlen(cptr)>4) {
      strlcpy(limfile, cptr+4, FILENAME_MAX); continue;
    } else if(strcasecmp(cptr, "LIM")==0) {
      strcpy(limfile, "stdout"); continue;
    } else if(strncasecmp(cptr, "BETA=", 5)==0 && strlen(cptr)>5) {
      if(!atof_with_check(cptr+5, &Beta) && Beta>0.0 && Beta<=1.0) continue;
    } else if(strncasecmp(cptr, "PH2O=", 5)==0) {
      if(!atof_with_check(cptr+5, &pc) && pc>0.0 && pc<=1.0) continue;
    } if(strncasecmp(cptr, "INPUT=", 6)==0) {
      strlcpy(inputfile, cptr+6, FILENAME_MAX); 
      if(strlen(inputfile)>0) continue;
    } else if(strncasecmp(cptr, "SVG=", 4)==0) {
      strlcpy(svgfile, cptr+4, FILENAME_MAX); 
      if(strlen(svgfile)>0) continue;
    } else if(strncasecmp(cptr, "FIT=", 4)==0) {
      strlcpy(fitfile, cptr+4, FILENAME_MAX); 
      if(strlen(fitfile)>0) continue;
    } else if(strncasecmp(cptr, "END=", 4)==0) {
      if(!atof_with_check(cptr+4, &fittime) && fittime>10.0) continue;
    }
    fprintf(stderr, "Error: invalid option '%s'.\n", argv[ai]);
    return(1);
  } else break;
  
  /* Print help or version? */
  if(help==2) {tpcHtmlUsage(argv[0], info, ""); return(0);}
  if(help) {tpcPrintUsage(argv[0], info, stdout); return(0);}
  if(version) {tpcPrintBuild(argv[0], stdout); return(0);}

  /* Process other arguments, starting from the first non-option */
  if(ai<argc) {strlcpy(tacfile, argv[ai++], FILENAME_MAX);}
  if(ai<argc) {strlcpy(lvcavname, argv[ai++], FILENAME_MAX);}
  if(ai<argc) {
    strlcpy(myocname, argv[ai++], FILENAME_MAX);
    if(!strcasecmp(myocname, "NONE") || !strcasecmp(myocname, "'NONE'") ||
       !strcasecmp(myocname, "NO") || !strcasecmp(myocname, "0")) 
      myocname[0]=(char)0;
  }
  if(ai<argc) {strlcpy(resfile, argv[ai++], FILENAME_MAX);}
  if(ai<argc) {
    /* we should never get this far */
    fprintf(stderr, "Error: too many arguments: '%s'.\n", argv[ai]);
    return(1);
  }
  if(doSD || doCL) doBootstrap=1; else doBootstrap=0;

  /* If only filename for initial values was given, then write one
     with default contents, and exit */
  if(limfile[0] && !tacfile[0]) {
    IFT ift; iftInit(&ift);
    /* Check that initial value file does not exist */
    if(strcasecmp(limfile, "stdout")!=0 && access(limfile, 0) != -1) {
      fprintf(stderr, "Error: parameter constraint file %s exists.\n", limfile);
      return(9);
    }
    if(verbose>1) printf("writing parameter constraints file\n");
    /* Create parameter file */
    iftPutDouble(&ift, "ptMBF_lower", def_pmin[0], NULL);
    iftPutDouble(&ift, "ptMBF_upper", def_pmax[0], NULL);
    iftPutDouble(&ift, "PTF_lower", def_pmin[1], NULL);
    iftPutDouble(&ift, "PTF_upper", def_pmax[1], NULL);
    iftPutDouble(&ift, "Va_lower", def_pmin[2], NULL);
    iftPutDouble(&ift, "Va_upper", def_pmax[2], NULL);
    if(iftWrite(&ift, limfile)) {
      fprintf(stderr, "Error in writing '%s': %s\n", limfile, ift.status);
      iftEmpty(&ift); return(9);
    }
    if(strcasecmp(limfile, "stdout")!=0)
      fprintf(stdout, "Parameter file %s with initial values written.\n",
              limfile);
    iftEmpty(&ift); return(0);
  }

  /* In verbose mode print arguments and options */
  if(verbose>1) {
    printf("limfile := %s\n", limfile);
    printf("tacfile := %s\n", tacfile);
    printf("lvcavname := %s\n", lvcavname);
    printf("myocname := %s\n", myocname);
    printf("resfile := %s\n", resfile);
    printf("fitfile := %s\n", fitfile);
    printf("svgfile := %s\n", svgfile);
    printf("inputfile := %s\n", inputfile);
    printf("beta := %g\n", Beta);
    printf("pH2O := %g\n", pc);
    printf("doBootstrap := %d\n", doBootstrap);
    printf("doSD := %d\n", doSD);
    printf("doCL := %d\n", doCL);
    if(fittime>0.0) printf("requested_fittime := %g\n", fittime);
    fflush(stdout);
  }

  /* Did we get all the information that we need? */
  if(!resfile[0]) {
    fprintf(stderr, "Error: missing command-line argument; use option --help\n");
    return(1);
  }


  /*
   *  Read model parameter initial values and upper and lower limits
   *  if file for that was given
   */
  if(limfile[0]) {
    IFT ift; iftInit(&ift);
    double v;
    if(verbose>1) printf("reading %s\n", limfile);
    if(iftRead(&ift, limfile, 1)) {
      fprintf(stderr, "Error in reading '%s': %s\n", limfile, ift.status);
      return(9);
    }
    if(verbose>10) iftWrite(&ift, "stdout");
    int n=0;
    /* ptMBF */
    if(iftGetDoubleValue(&ift, 0, "ptMBF_lower", &v)>=0) {def_pmin[0]=v; n++;}
    if(iftGetDoubleValue(&ift, 0, "ptMBF_upper", &v)>=0) {def_pmax[0]=v; n++;}
    if(iftGetDoubleValue(&ift, 0, "MBF_lower", &v)>=0) {def_pmin[0]=v; n++;}
    if(iftGetDoubleValue(&ift, 0, "MBF_upper", &v)>=0) {def_pmax[0]=v; n++;}
    /* PTF */
    if(iftGetDoubleValue(&ift, 0, "PTF_lower", &v)>=0) {def_pmin[1]=v; n++;}
    if(iftGetDoubleValue(&ift, 0, "PTF_upper", &v)>=0) {def_pmax[1]=v; n++;}
    /* Va */
    if(iftGetDoubleValue(&ift, 0, "Va_lower", &v)>=0) {def_pmin[2]=v; n++;}
    if(iftGetDoubleValue(&ift, 0, "Va_upper", &v)>=0) {def_pmax[2]=v; n++;}
    iftEmpty(&ift);
    if(n==0) {fprintf(stderr, "Error: invalid parameter file.\n"); return(9);}
  }
  /* Check that these limits are ok */
  {
    int pi, ret=0;
    fittedparNr=0;
    for(pi=0; pi<parNr; pi++) {
      if(def_pmin[pi]<0.0) ret++;
      if(def_pmax[pi]<def_pmin[pi]) ret++;
      if(def_pmax[pi]>def_pmin[pi]) fittedparNr++;
    }
    if(ret) {
      fprintf(stderr, "Error: invalid parameter constraints.\n");
      return(9);
    }
    if(fittedparNr==0) {
      fprintf(stderr, "Error: no model parameters left free for fitting.\n");
      return(9);
    }
  }
  if(verbose>1) {
    fflush(stdout); printf("Parameter constraints:\n");
    for(int pi=0; pi<parNr; pi++) {
      printf("def_pmin[%d] := %g\n", pi+1, def_pmin[pi]);
      printf("def_pmax[%d] := %g\n", pi+1, def_pmax[pi]);
    }
    printf("fittedParNr := %d\n", fittedparNr);
    fflush(stdout);
  }
  /* Convert MBF constraints to per second */
  def_pmin[0]/=60.0; def_pmax[0]/=60.0;


  /*
   *  Read data file
   */
  if(verbose>1) printf("reading '%s'.\n", tacfile);
  DFT dft; dftInit(&dft);
  ret=dftRead(tacfile, &dft);
  if(ret) {
    fprintf(stderr, "Error in reading '%s': %s\n", tacfile, dfterrmsg);
    if(verbose>1) printf("  ret :=%d\n", ret);
    return(2);
  }
  /* Check for NaN's */
  ret=dft_nr_of_NA(&dft); if(ret>0) {
    fprintf(stderr, "Error: missing sample(s) in %s.\n", tacfile);
    dftEmpty(&dft); return(2);
  }
  /* Sort the data by increasing sample times */
  dftSortByFrame(&dft);
  /* Guess time units, if necessary */
  if(dft.timeunit==TUNIT_UNKNOWN) {
    if(dft.x[dft.frameNr-1]>20.0) {
      if(verbose>1) printf("Note: assuming that times are in seconds.\n");
      dft.timeunit=TUNIT_SEC;
    } else {
      if(verbose>1) printf("Note: assuming that times are in minutes.\n");
      dft.timeunit=TUNIT_MIN;
    }
  }
  /* Set weights to one for data that did not contain weights */
  if(dft.isweight==0) for(int i=0; i<dft.frameNr; i++) dft.w[i]=1.0;
  if(verbose>3) {
    fprintf(stdout, "common_data_weights := %g", dft.w[0]);
    for(int i=1; i<dft.frameNr; i++) fprintf(stdout, ", %g", dft.w[i]);
    fprintf(stdout, "\n");
  }
  /* Convert time to sec if necessary */
  if(dft.timeunit==TUNIT_MIN) {
    dftMin2sec(&dft); originallyMinutes=1;
  }
  /* Check region number */
  if(dft.voiNr<2) {
    fprintf(stderr, "Error: check the contents of datafile.\n");
    dftEmpty(&dft); return(2);
  }
  /* Make sure that there is no overlap in frame times */
  if(dft.timetype==DFT_TIME_STARTEND) {
    if(dftDeleteFrameOverlap(&dft)) {
      fprintf(stderr, "Error: file has overlapping frame times.\n");
      dftEmpty(&dft); return(2);
    }
  }
  /* Set fit duration and get the number of fitted sample points */
  double starttime=0;
  double endtime=1.0E+90; if(fittime>0.0) endtime=fittime/60.0;
  int first, last;
  fitframeNr=fittime_from_dft(&dft, &starttime, &endtime,
                               &first, &last, verbose-1);
  if(verbose>2) {
    printf("frameNr := %d\n", dft.frameNr);
    printf("starttime := %g\n", starttime);
    printf("endtime := %g\n", endtime);
    printf("first := %d\n", first);
    printf("last := %d\n", last);
    printf("fitframeNr := %d\n", fitframeNr);
    fflush(stdout);
  }
  fittime=60.0*endtime;
  /* Check frame number */
  if(dftValidNr(&dft, 0.0, fittime, -1)<4) {
    fprintf(stderr, "Error: check the contents of datafile.\n");
    dftEmpty(&dft); return(2);
  }


  /*
   *  Find myocardial TAC and set wmroi = its index
   *  Find LV TAC and set lvroi = its index
   */
  int wmroi=-1;
  if(myocname[0]) {
    if(verbose>1) printf("searching for (whole) myocardium ROI.\n");
    int n, i; 
    n=dftSelectRegions(&dft, myocname, 1);
    if(verbose>1) printf("nr of myoc regions := %d/%d\n", n, dft.voiNr);
    if(n<=0) {
      fprintf(stderr, "Error: cannot find myoc region.\n");
      dftEmpty(&dft); return(2);
    }
    if(n==dft.voiNr) {
      fprintf(stderr, "Error: all regions match myoc name.\n");
      dftEmpty(&dft); return(2);
    }
    /* Try to select the best match */
    i=dftSelectBestReference(&dft); if(i<0) {
      fprintf(stderr, "Error: cannot select the best myoc region.\n");
      dftEmpty(&dft); return(2);
    }
    wmroi=i;
  }
 
  int lvroi=-1;
  {
    if(verbose>1) printf("searching for LV cavity ROI.\n");
    int n, i; 
    n=dftSelectRegions(&dft, lvcavname, 1);
    if(verbose>1) printf("nr of lvcav regions := %d/%d\n", n, dft.voiNr);
    if(n<=0) {
      fprintf(stderr, "Error: cannot find lvcav region.\n");
      dftEmpty(&dft); return(2);
    }
    if(n==dft.voiNr) {
      fprintf(stderr, "Error: all regions match lvcav name.\n");
      dftEmpty(&dft); return(2);
    }
    /* Try to select the best match */
    i=dftSelectBestReference(&dft); if(i<0) {
      fprintf(stderr, "Error: cannot select the best lvcav region.\n");
      dftEmpty(&dft); return(2);
    }
    lvroi=i;
  }
  
  if(wmroi==lvroi) {
    fprintf(stderr, "Error: cannot determine lvcav or myoc TAC.\n");
    dftEmpty(&dft); return(2);
  }
  if(verbose>1) {
    printf("selected lvcav region := %s\n", dft.voi[lvroi].name);
    if(wmroi>=0) printf("selected myoc region := %s\n", dft.voi[wmroi].name);
  }

  /* Allocate an extra TAC for the bootstrap */
  int bs_index=0;
  if(doBootstrap) {
    ret=dftAddmem(&dft, 1); if(ret) {
      fprintf(stderr, "Error: cannot allocate more memory.\n");
      dftEmpty(&dft); return(3);
    }
    bs_index=dft.voiNr;
    strcpy(dft.voi[bs_index].voiname, "BS");
    strcpy(dft.voi[bs_index].name, "BS");
  }



  /*
   *  Prepare the room for results
   */
  if(verbose>1) printf("initializing result data\n");
  RES res; resInit(&res);
  if(res_allocate_with_dft(&res, &dft)!=0) {
    fprintf(stderr, "Error: cannot setup memory for results.\n");
    dftEmpty(&dft); return(4);
  }
  /* Copy titles & filenames */
  tpcProgramName(argv[0], 1, 1, res.program, 256);
  strcpy(res.datafile, tacfile);
  if(wmroi>=0) strlcpy(res.refroi, dft.voi[wmroi].name, 64);
  strcpy(res.fitmethod, "TGO");
  /* Constants */
  res.beta=Beta; res.Vb=-1.0;
  res.isweight=dft.isweight;
  /* Set data range */
  sprintf(res.datarange, "%g - %g min", starttime, endtime);
  res.datanr=fitframeNr;
  /* Set current time to results */
  res.time=time(NULL);


  /* Set parameter number, including also the extra "parameters"
     and the parameter names and units */
  res.parNr=5;
  {
    int pi;
    pi=0; strcpy(res.parname[pi], "ptMBF"); 
    strcpy(res.parunit[pi], "mL/(min*mL)");
    pi++; strcpy(res.parname[pi], "PTF"); strcpy(res.parunit[pi], "mL/mL");
    pi++; strcpy(res.parname[pi], "Va");
    strcpy(res.parunit[pi], "mL/mL");
    pi++; strcpy(res.parname[pi], "rMBF"); 
    strcpy(res.parunit[pi], "mL/(min*mL)");
    pi++; strcpy(res.parname[pi], "WSS"); strcpy(res.parunit[pi], "");
  }

  /*
   *  Allocate memory for input curve (arterial blood concentration)
   */
  dftInit(&input);
  ret=dftSetmem(&input, dft.frameNr, 1);
  if(ret) {
    fprintf(stderr, "Error: cannot allocate memory for input TAC.\n");
    dftEmpty(&dft); resEmpty(&res); return(4);
  }
  input.voiNr=1; input.frameNr=dft.frameNr;
  (void)dftCopymainhdr(&dft, &input);
  (void)dftCopyvoihdr(&dft, lvroi, &input, 0);
  for(int i=0; i<input.frameNr; i++) {
    input.x[i]=dft.x[i]; input.x1[i]=dft.x1[i]; input.x2[i]=dft.x2[i];
  }

  /*
   *  Allocate memory for fitted curves
   */
  DFT fit; dftInit(&fit);
  ret=dftdup(&dft, &fit);
  if(ret) {
    fprintf(stderr, "Error: cannot allocate memory for fitted curves.\n");
    if(verbose>1) printf("  ret :=%d\n", ret);
    dftEmpty(&dft); dftEmpty(&input); resEmpty(&res);
    return(4);
  }


  /*
   *  Determine the arterial input TAC, either by directly using LV cavity,
   *  or by using whole myocardial TAC and MBF model
   */
  for(int i=0; i<input.frameNr; i++) {
    input.voi[0].y[i]=dft.voi[lvroi].y[i];
  }
  if(wmroi<0) {
    if(verbose>1) printf("Note: using LV cavity directly as the input.\n");
  } else {
    if(verbose>1) printf("starting myoc fitting\n");
    int tgoNr, neighNr, iterNr;
    double *sd, *cl1, *cl2;
    double wss;

    /* Set parameter constraints */
    for(int pi=0; pi<parNr; pi++) {
      pmin[pi]=def_pmin[pi];  pmax[pi]=def_pmax[pi];
    }

    /* set data pointers */
    petmeas=dft.voi[wmroi].y; petsim=fit.voi[wmroi].y;
    weight=dft.w;

    /* Fit */
    if(verbose>2) printf("  fitting\n");
    TGO_LOCAL_INSIDE = 1;
    TGO_SQUARED_TRANSF = 0;
    tgoNr=200;
    neighNr=20;
    iterNr=0;
    ret=tgo(
      pmin, pmax, mbfFunc, NULL, parNr, neighNr,
      &wss, res.voi[wmroi].parameter, tgoNr, iterNr, 
      verbose-8);
    if(ret>0) {
      fprintf(stderr, "\nError in optimization (%d).\n", ret);
      dftEmpty(&dft); dftEmpty(&input); dftEmpty(&fit); resEmpty(&res);
      return(5);
    }
    /* Correct fitted parameters to match constraints like inside function */
    (void)modelCheckParameters(parNr, pmin, pmax, res.voi[wmroi].parameter,
                               res.voi[wmroi].parameter, NULL);
    wss=wss_wo_penalty;
    res.voi[wmroi].parameter[res.parNr-1]=wss;

    /* Print measured and fitted muscle TAC */
    if(verbose>5) {
      printf("     Measured  Fitted    Weight:\n");
      for(int fi=0; fi<fitframeNr; fi++)
        printf("  %2d  %8.2e  %8.2e  %8.2e\n", fi+1, petmeas[fi], petsim[fi],
               weight[fi]);
    }

    /* Bootstrap */
    if(doBootstrap) {
      if(verbose>1) printf("  bootstrapping\n");
      char buf[64];
      /* bootstrap needs and changes petmeas[] and petsim[], */
      /* therefore reset these pointers */
      petmeas=dft.voi[bs_index].y;
      petsim=dft.voi[bs_index].y2;
      /* set pointer for SD and CL arrays */
      if(doSD) sd=res.voi[wmroi].sd; else sd=NULL;
      if(doCL) {cl1=res.voi[wmroi].cl1; cl2=res.voi[wmroi].cl2;} else cl1=cl2=NULL;
      ret=bootstrap(0, cl1, cl2, sd, res.voi[wmroi].parameter, pmin, pmax,
                    fitframeNr, dft.voi[wmroi].y, fit.voi[wmroi].y, petmeas,
                    parNr, dft.w, mbfFunc, buf, verbose-6);
      if(ret) {
        fprintf(stderr, "\nError in bootstrap: %s\n", buf);
        for(int pi=0; pi<parNr; pi++) {
          if(doSD) sd[pi]=nan(""); 
          if(doCL) cl1[pi]=cl2[pi]=nan("");
        }
      }
      // return data pointers back to what they were
      petmeas=dft.voi[wmroi].y; petsim=fit.voi[wmroi].y;
    }

    /*
     *  Calculate arterial blood curve to be used as input with other regions
     */
    {
      double alpha, Va;
      alpha=res.voi[wmroi].parameter[1];
      Va=res.voi[wmroi].parameter[2];
      for(int fi=0; fi<dft.frameNr; fi++) {
        /*
        input.voi[0].y[fi]=
          ((1.0-Beta)*petmeas[fi]-alpha*input[fi])/(Va*(1.0-Beta)-alpha*Beta);
        */
        /* This is the way the input was calculated in fitmbf 2.0 */
        input.voi[0].y[fi]=
          ((1.0-Beta)*petsim[fi]-alpha*dft.voi[lvroi].y[fi])
          / (Va*(1.0-Beta)-alpha*Beta);
      }
    }
  }


  /*
   *  Fit myocardial regions with arterial blood input
   *  originating from fit of whole myocardial ROI with LV input
   */
   
  /* One region at a time */
  for(int ri=0; ri<dft.voiNr; ri++) {
  
    /* Do not fit LV TAC */
    if(ri==lvroi) continue;
    /* Do not fit whole myocardium again */
    if(ri==wmroi) continue;

    if(verbose>1) printf("starting %s fitting\n", dft.voi[ri].name);
    int tgoNr, neighNr, iterNr;
    double *sd, *cl1, *cl2;
    double wss;

    /* Set parameter constraints */
    for(int pi=0; pi<parNr; pi++) {
      pmin[pi]=def_pmin[pi];  pmax[pi]=def_pmax[pi];
    }

    /* set data pointers */
    petmeas=dft.voi[ri].y; petsim=fit.voi[ri].y;
    weight=dft.w;

    /* Fit */
    if(verbose>2) printf("  fitting\n");
    TGO_LOCAL_INSIDE = 1;
    TGO_SQUARED_TRANSF = 0;
    tgoNr=200;
    neighNr=20;
    iterNr=0;
    ret=tgo(
      pmin, pmax, mbfFunc2, NULL, parNr, neighNr,
      &wss, res.voi[ri].parameter, tgoNr, iterNr, 
      verbose-8);
    if(ret>0) {
      fprintf(stderr, "\nError in optimization (%d).\n", ret);
      dftEmpty(&dft); dftEmpty(&input); dftEmpty(&fit); resEmpty(&res);
      return(5);
    }
    /* Correct fitted parameters to match constraints like inside function */
    (void)modelCheckParameters(parNr, pmin, pmax, res.voi[ri].parameter,
                               res.voi[ri].parameter, NULL);
    wss=wss_wo_penalty;
    res.voi[ri].parameter[res.parNr-1]=wss;

    /* Print measured and fitted muscle TAC */
    if(verbose>5) {
      printf("     Measured  Fitted    Weight:\n");
      for(int fi=0; fi<fitframeNr; fi++)
        printf("  %2d  %8.2e  %8.2e  %8.2e\n", fi+1, petmeas[fi], petsim[fi],
               weight[fi]);
    }

    /* Bootstrap */
    if(doBootstrap) {
      if(verbose>1) printf("  bootstrapping\n");
      char buf[64];
      /* bootstrap needs and changes petmeas[] and petsim[], */
      /* therefore reset these pointers */
      petmeas=dft.voi[bs_index].y;
      petsim=dft.voi[bs_index].y2;
      /* set pointer for SD and CL arrays */
      if(doSD) sd=res.voi[ri].sd; else sd=NULL;
      if(doCL) {cl1=res.voi[ri].cl1; cl2=res.voi[ri].cl2;} else cl1=cl2=NULL;
      ret=bootstrap(0, cl1, cl2, sd, res.voi[ri].parameter, pmin, pmax,
                    fitframeNr, dft.voi[ri].y, fit.voi[ri].y, petmeas,
                    parNr, dft.w, mbfFunc2, buf, verbose-6);
      if(ret) {
        fprintf(stderr, "\nError in bootstrap: %s\n", buf);
        for(int pi=0; pi<parNr; pi++) {
          if(doSD) sd[pi]=nan(""); 
          if(doCL) cl1[pi]=cl2[pi]=nan("");
        }
      }
      // return data pointers back to what they were
      petmeas=dft.voi[ri].y; petsim=fit.voi[ri].y;
    }

  } // next TAC




  /* Delete LV region from results */
  (void)resDelete(&res, lvroi);

  /* Convert MBF estimates from 1/sec to 1/min values */
  for(int ri=0; ri<res.voiNr; ri++) {
    res.voi[ri].parameter[0]*=60.;
    res.voi[ri].sd[0]*=60.;
    res.voi[ri].cl1[0]*=60.;
    res.voi[ri].cl2[0]*=60.;
  }
  /* Calculate rMBF from ptMBF and PTF */
  for(int ri=0; ri<res.voiNr; ri++) {
    res.voi[ri].parameter[3]=res.voi[ri].parameter[0]*res.voi[ri].parameter[1];
  }



  /*
   *  Print results on screen
   */
  if(verbose>0) {resPrint(&res); fprintf(stdout, "\n");}


  /*
   *  Save results
   */
  if(verbose>1) printf("saving results\n");
  if(resWrite(&res, resfile, verbose-3)!=0) {
    fprintf(stderr, "Error in writing '%s': %s\n", resfile, reserrmsg);
    dftEmpty(&dft); dftEmpty(&input); dftEmpty(&fit); resEmpty(&res);
    return(11);
  }
  if(verbose>1) fprintf(stdout, "Model parameters written in %s\n", resfile);



  /*
   *  Convert TAC time units back to minutes, if necessary
   */
  if(originallyMinutes) {
    dftSec2min(&dft); dftSec2min(&input); dftSec2min(&fit);
  }


  /*
   *  Saving and/or plotting of fitted TACs
   */
  if(svgfile[0] || fitfile[0]) {

    /* Save fitted TACs */
    if(fitfile[0]) {
      if(verbose>1) printf("saving fitted curves\n");
      if(dftWrite(&fit, fitfile)) {
        fprintf(stderr, "Error in writing '%s': %s\n", fitfile, dfterrmsg);
      } else if(verbose>0) printf("fitted TACs written in %s\n", fitfile);
    }

    /* Save SVG plot of fitted and original data */
    if(svgfile[0]) {
      if(verbose>1) printf("saving SVG plot\n");
      char tmp[64];
      sprintf(tmp, "MBF fit ");
      if(strlen(dft.studynr)>0) strlcat(tmp, dft.studynr, 64);
      ret=plot_fitrange_svg(&dft, &fit, tmp, 0.0, 1.03*dft.x[fitframeNr-1],
                            0.0, nan(""), svgfile, verbose-8);
      if(ret) {
        fprintf(stderr, "Error (%d) in writing '%s'.\n", ret, svgfile);
      } else if(verbose>0) printf("plots written in %s\n", svgfile);
    }

  }


  /*
   *  Save Ca in a file
   */
  if(inputfile[0]) {
    if(verbose>1) printf("saving arterial blood data in %s\n", inputfile);
    if(dftWrite(&input, inputfile)) {
      fprintf(stderr, "Error in writing %s: %s\n", inputfile, dfterrmsg);
      dftEmpty(&dft); dftEmpty(&input); dftEmpty(&fit); resEmpty(&res);
      return(11);
    }
    if(verbose>0) printf("Estimated arterial blood TAC saved in %s\n", inputfile);
  }


  dftEmpty(&dft); dftEmpty(&input); dftEmpty(&fit); resEmpty(&res);

  return(0);
}
/*****************************************************************************/

/*****************************************************************************
 *
 *  Functions to be minimized
 *
 *****************************************************************************/
/* Function for whole myocardial muscle and LV */
double mbfFunc(int parNr, double *p, void *fdata)
{
  int fi, ret;
  double wss=0.0, d, K1, k2, Vfit, flow, Va, alpha;
  double pa[MAX_PARAMETERS], penalty=1.0;

  /* Check parameters against the constraints */
  ret=modelCheckParameters(parNr, pmin, pmax, p, pa, &penalty);
  if(fdata) {}
  
  /* Calculate K1, k2 and Vfit */
  flow=pa[0]; alpha=pa[1]; Va=pa[2];
  K1=(flow/Beta)*(alpha+Va/pc);
  k2=flow*(1.0/pc+(1.0-Beta)/Beta);
  Vfit=Va/Beta;

  /* Simulate whole muscle curve and compute weighted SS */
  ret=simMBF(input.x, input.voi[0].y, input.frameNr, K1, k2, Vfit, petsim);
  if(ret) {
    fprintf(stderr, "error %d in simulation\n", ret);
    return(nan(""));
  }

  for(fi=0; fi<fitframeNr; fi++) if(weight[fi]>0.0) {
    d=petmeas[fi]-petsim[fi]; wss+=weight[fi]*d*d;
  }
  wss_wo_penalty=wss;
  wss*=penalty;

  return(wss);
}
/* function for smaller myocardial regions and Ca */
double mbfFunc2(int parNr, double *p, void *fdata)
{
  int fi, ret;
  double wss=0.0, d, K1, k2, Vfit, flow, Va, alpha;
  double pa[MAX_PARAMETERS], penalty=1.0;

  /* Check parameters against the constraints */
  ret=modelCheckParameters(parNr, pmin, pmax, p, pa, &penalty);
  if(fdata) {}  

  /* Calculate K1, k2 and Vfit */
  flow=pa[0]; alpha=pa[1]; Va=pa[2];
  K1=flow*(alpha+Va/pc);
  k2=flow/pc;
  Vfit=Va;

  /* Simulate muscle curve and compute weighted SS */
  ret=simMBF(input.x, input.voi[0].y, input.frameNr, K1, k2, Vfit, petsim);
  if(ret) {
    fprintf(stderr, "error %d in simulation\n", ret);
    return(nan(""));
  }

  for(fi=0; fi<fitframeNr; fi++) if(weight[fi]>0.0) {
    d=petmeas[fi]-petsim[fi]; wss+=weight[fi]*d*d;
  }
  wss_wo_penalty=wss;
  wss*=penalty;

  return(wss);
}
/*****************************************************************************/

/*****************************************************************************/
/// @endcond
