/** @file fitk5.c
 *  @brief Estimates the parameters of irreversible 3-tissue compartment model.
 *  @copyright (c) Turku PET Centre
 *  @author Vesa Oikonen
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
const int parNr=6;
DFT input, dft;
double *petmeas, *petsim;
double fVb=-1.0;
double pmin[MAX_PARAMETERS], pmax[MAX_PARAMETERS];
double fk1k2=-1.0, fk3k4=-1.0;
int fitframeNr;
double wss_wo_penalty=0.0;
/** Definitions for model parameter indices */
enum parameter {
  CM_K1, CM_K1K2, CM_K3, CM_K3K4, CM_K5, CM_VB, CM_KI, CM_WSS, CM_AIC
};
/*****************************************************************************/
/* Local functions */
double func_sk5(int parNr, double *p, void*);
double func_pk5(int parNr, double *p, void*);
int resK1k2Median(char *filename, double *median, double *mean, char *msg);
/*****************************************************************************/

/*****************************************************************************/
static char *info[] = {
  "Non-linear fitting of irreversible three-tissue compartment model to plasma",
  "input, blood, and tissue time-activity curves (PTAC, BTAC, and TTAC) to",
  "estimate parameters K1, K1/k2, k3, k3/k4, k5, and Vb.",
  " ",
  "Model with two of tissue compartments in parallel:",
  "  _____         ____         ____  ",
  " |     |   K1  |    |   k3  |    | ",
  " |  Ca | ----> |    | ----> | C2 | ",
  " |_____|       |    | <---- |____| ",
  "  _____        | C1 |   k4   ____  ",
  " |     |       |    |       |    | ",
  " |  Cv | <---- |    | ----> | C3 | ",
  " |_____|   k2  |____|   k5  |____| ",
  " ",                               
  "Model with compartments in series:",
  "  _____         ____         ____         ____  ",
  " |     |   K1  |    |   k3  |    |  k5   |    | ",
  " |  Ca | ----> |    | ----> | C2 | ----> | C3 | ",
  " |_____|       |    | <---- |____|       |____| ",
  "  _____        | C1 |   k4                      ",  
  " |     |       |    | ",
  " |  Cv | <---- |    | ",
  " |_____|   k2  |____| ",
  " ",              
  "Sample times must be in minutes.",
  " ",
  "Usage: @P [Options] ptacfile btacfile ttacfile endtime resultfile",
  " ",
  "Options:",
  " -lim[=<filename>]",
  "     Specify the constraints for model parameters;",
  "     This file with default values can be created by giving this",
  "     option as the only command-line argument to this program.",
  "     Without file name the default values are printed on screen.",
  " -SD[=<y|N>]",
  "     Standard deviations are estimated and saved in results (y),",
  "     or not calculated (N, default).",
  "     Program runs a lot faster if SD and CL are not calculated.",
  " -CL[=<y|N>]",
  "     95% Confidence limits are estimated and saved in results (y), or",
  "     not calculated (N, default).",
  " -model=<parallel|series>",
  "     Specify the model that is fitted to data; compartments in parallel",
  "     (default) or in series.",
  " -Vb=<Vb(%)>",
  "     Enter a fixed Vb; fitted by default.",
  " -fk1k2=<<value> || <result filename>>",
  "     K1/k2 is constrained to the given value in all regions; if result",
  "     file name is entered, then K1/k2 is constrained to the median of",
  "     regional K1/k2 values in the result file.",
/* Not tested, thus commented out
   Also, it might be more useful to use this option to fix k3/k4, at least
   in parallel model...
  " -r=<Reference region id or filename>",
  "     Optional reference region is used to constrain K1/k2 in other regions;",
  "     with option -rmod=3 also k3 and k4 are fitted to reference region data,",
  "     thus any large brain region (for example cortex) could be used here.",
  " -rmod=<2|3>",
  "     Specify the model (2- or 3-compartments) that is fitted to reference",
  "     region data; by default 2-CM (1-tissue compartment model).",
*/
  " -fit=<Filename>",
  "     Fitted regional TACs are written in DFT format.",
  " -svg=<Filename>",
  "     Fitted and measured TACs are plotted in specified SVG file.",
  " -stdoptions", // List standard options like --help, -v, etc
  " ",
  "Example 1: estimate K1, K1/k2, k3, k3/k4, k5 and Vb",
  "     @P -svg=ua929fit.svg ua929ap.bld ua929ab.bld ua929.tac 90 ua929.res",
  " ",
  "Example 2: estimate K1, k3, k3/k4, and k5; Vb is constrained to 1.5%, and",
  "K1/k2 is constrained to the regional median",
  "     @P -Vb=1.5 ua929ap.kbq ua929ab.kbq ua929.tac 90 tmp.res",
  "     @P -Vb=1.5 -fk1k2=tmp.res ua929ap.kbq ua929ab.kbq ua929.tac 90 ua929.res",
  " ",
  "See also: logan, fitk2, fitk3, fitk4, p2t_v3c, dftweigh, rescoll",
  " ",
  "Keywords: TAC, modelling, Ki, k5, irreversible uptake, 3TCM",
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

/*****************************************************************************
 *
 *  Main
 *
 *****************************************************************************/
int main(int argc, char *argv[])
{
  int          ai, help=0, version=0, verbose=1;
  int          ri, fi, pi, m, n, ret;
  int          model=1; // 0=series, 1=parallel
  int          ref=-1, refAdded=0, inputtype;
  char         dfile[FILENAME_MAX], pfile[FILENAME_MAX], bfile[FILENAME_MAX],
               rfile[FILENAME_MAX], ffile[FILENAME_MAX], svgfile[FILENAME_MAX],
               limfile[FILENAME_MAX];
  char        *cptr, tmp[FILENAME_MAX], refname[FILENAME_MAX];
  double       fitdur, K1, k2, k3, k4, k5, f;
  RES          res;
  IFT          ift;
  int          doBootstrap=0, doSD=0, doCL=0;
  double      *sd, *cl1, *cl2;
  double       def_pmin[MAX_PARAMETERS], def_pmax[MAX_PARAMETERS];
  int          tgoNr=0, neighNr=0, iterNr=0, fittedparNr=0;
  double     (*func)(int, double*, void*);

#ifdef MINGW
  // Use Unix/Linux default of two-digit exponents in MinGW on Windows
  _set_output_format(_TWO_DIGIT_EXPONENT);
#endif

  //drandSeed(195); // tgo() does it

  /* Set parameter initial values and constraints */
  /* K1    */ def_pmin[CM_K1]=0.0;       def_pmax[CM_K1]=5.0;
  /* K1/k2 */ def_pmin[CM_K1K2]=0.00001; def_pmax[CM_K1K2]=10.0;
  /* k3    */ def_pmin[CM_K3]=0.0;       def_pmax[CM_K3]=2.0;
  /* k3/k4 */ def_pmin[CM_K3K4]=0.00001; def_pmax[CM_K3K4]=2.0;
  /* k5    */ def_pmin[CM_K5]=0.0;       def_pmax[CM_K5]=2.0;
  /* Vb    */ def_pmin[CM_VB]=0.0;       def_pmax[CM_VB]=0.50;


  /*
   *  Get arguments
   */
  if(argc==1) {tpcPrintUsage(argv[0], info, stderr); return(1);}
  dfile[0]=pfile[0]=bfile[0]=rfile[0]=ffile[0]=refname[0]=svgfile[0]=(char)0;
  limfile[0]=(char)0;
  fitdur=fVb=-1.0; 
  iftInit(&ift); resInit(&res); 
  dftInit(&dft); dftInit(&input);
  /* Get options first, because it affects what arguments are read */
  for(ai=1; ai<argc; ai++) if(*argv[ai]=='-') {
    cptr=argv[ai]+1; if(*cptr=='-') cptr++; if(cptr==NULL) continue;
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
    } else if(strncasecmp(cptr, "R=", 2)==0 && strlen(cptr)>2) {
      strlcpy(refname, cptr+2, FILENAME_MAX); continue;
    } else if(strncasecmp(cptr, "MODEL=", 6)==0) {
      cptr+=6;
      if(strncasecmp(cptr, "PARALLEL", 1)==0) {model=1; continue;}
      if(strncasecmp(cptr, "SERIES", 1)==0) {model=0; continue;}
    } else if(strncasecmp(cptr, "FK1K2=", 6)==0 && strlen(cptr)>6) {
      /* Try to read specified K1/k2 from option */
      fk1k2=atof_dpi(cptr+6); def_pmin[CM_K1K2]=def_pmax[CM_K1K2]=fk1k2;
      if(fk1k2>0.0) continue;
      /* If not successful, then try to read median K1/k2 from specified
         result file */
      if(resK1k2Median(cptr+6, &fk1k2, NULL, tmp)==0) {
        def_pmin[CM_K1K2]=def_pmax[CM_K1K2]=fk1k2; continue;
      }
      if(verbose>1) fprintf(stderr, "Error: %s\n", tmp);  
    } else if(strncasecmp(cptr, "FK3K4=", 6)==0 && strlen(cptr)>6) {
      /* Try to read specified K3/k4 from option */
      fk3k4=atof_dpi(cptr+6); def_pmin[CM_K3K4]=def_pmax[CM_K3K4]=fk3k4;
      if(fk3k4>0.0) continue;
    } else if(strncasecmp(cptr, "Vb=", 3)==0 && strlen(cptr)>3) {
      fVb=0.01*atof_dpi(cptr+3);
      if(fVb>=0.0 && fVb<1.0) {
        if(fVb<0.01) fprintf(stderr, "Warning: Vb was set to %g%%\n", 100.*fVb);
        def_pmin[CM_VB]=def_pmax[CM_VB]=fVb;
        continue;
      }
      fVb=-1.0;
    } else if(strncasecmp(cptr, "FIT=", 4)==0) {
      strlcpy(ffile, cptr+4, FILENAME_MAX); if(strlen(ffile)>0) continue;
    } else if(strncasecmp(cptr, "SVG=", 4)==0) {
      strlcpy(svgfile, cptr+4, FILENAME_MAX); if(strlen(svgfile)>0) continue;
    }
    fprintf(stderr, "Error: invalid option '%s'.\n", argv[ai]);
    return(1);
  } else break;
  
  /* Print help or version? */
  if(help==2) {tpcHtmlUsage(argv[0], info, ""); return(0);}
  if(help) {tpcPrintUsage(argv[0], info, stdout); return(0);}
  if(version) {tpcPrintBuild(argv[0], stdout); return(0);}

  /* Process other arguments, starting from the first non-option */
  for(; ai<argc; ai++) {
    if(!pfile[0]) {
      strlcpy(pfile, argv[ai], FILENAME_MAX); continue;
    } else if(!bfile[0]) {
      strlcpy(bfile, argv[ai], FILENAME_MAX); continue;
    } else if(!dfile[0]) {
      strlcpy(dfile, argv[ai], FILENAME_MAX); continue;
    } else if(fitdur<0) {
      if(!atof_with_check(argv[ai], &fitdur) && fitdur>=0.0) continue;
      fprintf(stderr, "Error: invalid fit time '%s'.\n", argv[ai]);
      return(1);
    } else if(!rfile[0]) {
      strlcpy(rfile, argv[ai], FILENAME_MAX); continue;
    }
    /* we should never get this far */
    fprintf(stderr, "Error: too many arguments: '%s'.\n", argv[ai]);
    return(1);
  }
  if(doSD || doCL) doBootstrap=1; else doBootstrap=0;

  /* In verbose mode print arguments and options */
  if(verbose>1) {
    printf("pfile := %s\n", pfile);
    //printf("bfile := %s\n", bfile); Later!
    printf("dfile := %s\n", dfile);
    printf("rfile := %s\n", rfile);
    printf("ffile := %s\n", ffile);
    printf("svgfile := %s\n", svgfile);
    printf("limfile := %s\n", limfile);
    printf("refname := %s\n", refname);
    printf("model := %d\n", model);
    printf("fitdur := %g\n", fitdur);
    //if(fVb>=0.0) printf("fVb := %g\n", fVb); Later!
    //if(fk1k2>0.0) printf("fk1k2 := %g\n", fk1k2); Later!
    //if(fk3k4>0.0) printf("fk3k4 := %g\n", fk3k4); Later!
    printf("doBootstrap := %d\n", doBootstrap);
    printf("doSD := %d\n", doSD);
    printf("doCL := %d\n", doCL);
  }


  /* If only file name for parameter constraints was given, then write one
     with default contents, and exit */
  if(limfile[0] && !pfile[0]) {
    /* Check that initial value file does not exist */
    if(strcasecmp(limfile, "stdout")!=0 && access(limfile, 0) != -1) {
      fprintf(stderr, "Error: parameter constraint file %s exists.\n", limfile);
      return(9);
    }
    if(verbose>1 && strcasecmp(limfile, "stdout")!=0) 
      printf("writing parameter constraints file\n");
    /* Create parameter file */
    iftPutDouble(&ift, "K1_lower", def_pmin[CM_K1], NULL);
    iftPutDouble(&ift, "K1_upper", def_pmax[CM_K1], NULL);
    iftPutDouble(&ift, "K1k2_lower", def_pmin[CM_K1K2], NULL);
    iftPutDouble(&ift, "K1k2_upper", def_pmax[CM_K1K2], NULL);
    iftPutDouble(&ift, "k3_lower", def_pmin[CM_K3], NULL);
    iftPutDouble(&ift, "k3_upper", def_pmax[CM_K3], NULL);
    iftPutDouble(&ift, "k3k4_lower", def_pmin[CM_K3K4], NULL);
    iftPutDouble(&ift, "k3k4_upper", def_pmax[CM_K3K4], NULL);
    iftPutDouble(&ift, "k5_lower", def_pmin[CM_K5], NULL);
    iftPutDouble(&ift, "k5_upper", def_pmax[CM_K5], NULL);
    iftPutDouble(&ift, "Vb_lower", def_pmin[CM_VB], NULL);
    iftPutDouble(&ift, "Vb_upper", def_pmax[CM_VB], NULL);
    ret=iftWrite(&ift, limfile);
    if(ret) {
      fprintf(stderr, "Error in writing '%s': %s\n", limfile, ift.status);
      iftEmpty(&ift); return(9);
    }
    if(strcasecmp(limfile, "stdout")!=0)
      fprintf(stdout, "Parameter file %s with initial values written.\n", 
              limfile);
    iftEmpty(&ift); return(0);
  }


  /* Did we get all the information that we need? */
  if(fitdur==0) fitdur=1.0E+100; 
  else if(fitdur<0) {tpcPrintUsage(argv[0], info, stderr); return(1);}
  if(!rfile[0]) {
    fprintf(stderr, "Error: missing command-line argument; use option --help\n");
    return(1);
  }


  /*
   *  Read model parameter upper and lower limits
   *  if file for that was given
   */
  if(limfile[0]) {
    if(verbose>1) printf("reading %s\n", limfile);
    ret=iftRead(&ift, limfile, 1);
    if(ret) {
      fprintf(stderr, "Error in reading '%s': %s\n", limfile, ift.status);
      return(9);
    }
    if(verbose>2) iftWrite(&ift, "stdout");
    double v;
    int n=0;
    /* K1 */
    if(iftGetDoubleValue(&ift, 0, "K1_lower", &v)>=0) {def_pmin[CM_K1]=v; n++;}
    if(iftGetDoubleValue(&ift, 0, "K1_upper", &v)>=0) {def_pmax[CM_K1]=v; n++;}
    /* K1/k2 */
    if(iftGetDoubleValue(&ift, 0, "K1k2_lower", &v)>=0) {def_pmin[CM_K1K2]=v; n++;}
    if(iftGetDoubleValue(&ift, 0, "K1k2_upper", &v)>=0) {def_pmax[CM_K1K2]=v; n++;}
    /* k3 */
    if(iftGetDoubleValue(&ift, 0, "k3_lower", &v)>=0) {def_pmin[CM_K3]=v; n++;}
    if(iftGetDoubleValue(&ift, 0, "k3_upper", &v)>=0) {def_pmax[CM_K3]=v; n++;}
    /* k3/k4 */
    if(iftGetDoubleValue(&ift, 0, "k3k4_lower", &v)>=0) {def_pmin[CM_K3K4]=v; n++;}
    if(iftGetDoubleValue(&ift, 0, "k3k4_upper", &v)>=0) {def_pmax[CM_K3K4]=v; n++;}
    /* k5 */
    if(iftGetDoubleValue(&ift, 0, "k5_lower", &v)>=0) {def_pmin[CM_K5]=v; n++;}
    if(iftGetDoubleValue(&ift, 0, "k5_upper", &v)>=0) {def_pmax[CM_K5]=v; n++;}
    /* Vb */
    if(iftGetDoubleValue(&ift, 0, "Vb_lower", &v)>=0) {def_pmin[CM_VB]=v; n++;}
    if(iftGetDoubleValue(&ift, 0, "Vb_upper", &v)>=0) {def_pmax[CM_VB]=v; n++;}
    iftEmpty(&ift);
    if(n==0) {fprintf(stderr, "Error: invalid parameter file.\n"); return(9);}
  }
  /* Check that these are ok */
  for(pi=n=0, ret=0; pi<parNr; pi++) {
    if(def_pmin[pi]<0.0) ret++;
    if(def_pmax[pi]<def_pmin[pi]) ret++;
    if(def_pmax[pi]>def_pmin[pi]) n++;
  }
  if(ret) {
    fprintf(stderr, "Error: invalid parameter constraints.\n");
    return(9);
  }
  if(n==0) {
    fprintf(stderr, "Error: no model parameters left free for fitting.\n");
    return(9);
  }

  /* Fixed/fitted Vb */
  if(fVb>=0.0) def_pmin[CM_VB]=def_pmax[CM_VB]=fVb;
  if(def_pmin[4]==def_pmax[CM_VB]) fVb=def_pmin[CM_VB];
  if(fVb==0.0) strcpy(bfile, "");
  if(verbose>1) {
    printf("bfile := %s\n", bfile);
    if(fVb>=0.0) printf("fVb := %g\n", fVb);
  }
  /* Fixed/fitted K1/k2 */
  if(fk1k2>0.0) def_pmin[CM_K1K2]=def_pmax[CM_K1K2]=fk1k2;
  else if(def_pmin[CM_K1K2]==def_pmax[CM_K1K2]) fk1k2=def_pmin[CM_K1K2];
  if(verbose>1) {
    if(fk1k2>0.0) printf("fk1k2 := %g\n", fk1k2);
  }
  /* Fixed/fitted K3/k4 */
  if(fk3k4>0.0) def_pmin[CM_K3K4]=def_pmax[CM_K3K4]=fk3k4;
  else if(def_pmin[CM_K3K4]==def_pmax[CM_K3K4]) fk3k4=def_pmin[CM_K3K4];
  if(verbose>1) {
    if(fk3k4>0.0) printf("fk3k4 := %g\n", fk3k4);
  }


  /*
   *  Read tissue and input data
   */
  if(verbose>1) printf("reading tissue and input data\n");
  ret=dftReadModelingData(dfile, pfile, bfile, NULL, &fitdur, 
                          &fitframeNr, &dft, &input, stdout, verbose-2, tmp);
  if(ret!=0) {
    fprintf(stderr, "Error: %s\n", tmp);
    return(2);
  }
  if(fitframeNr<6 || input.frameNr<6) {
    fprintf(stderr, "Error: too few samples in specified fit duration.\n");
    dftEmpty(&input); dftEmpty(&dft); return(2);
  }
  /* If there is no blood TAC, then create a zero blood TAC */
  if(input.voiNr<2) {
    if(verbose>2) printf("setting blood tac to zero\n");
    ret=dftAddmem(&input, 1);
    if(ret) {
      fprintf(stderr, "Error: cannot allocate more memory.\n");
      dftEmpty(&dft); dftEmpty(&input); return(3);
    }
    strcpy(input.voi[1].voiname, "blood");
    strcpy(input.voi[1].name, input.voi[1].voiname);
    for(fi=0; fi<input.frameNr; fi++) input.voi[1].y[fi]=0.0;
    input.voiNr=2;
  }
  if(verbose>10) dftPrint(&dft);
  if(verbose>10) dftPrint(&input);
  /* Print the weights */
  if(verbose>2) {
    fprintf(stdout, "common_data_weights := %g", dft.w[0]);
    for(fi=1; fi<dft.frameNr; fi++) fprintf(stdout, ", %g", dft.w[fi]);
    fprintf(stdout, "\n");
  }


  /*
   *  Read reference TAC
   */
  /* Check if user even wants any reference region */
  if(!refname[0]) {
    if(verbose>1) printf("no reference region data\n");
    ref=-1;
  } else {
    if(verbose>1) printf("reading reference region data\n");
    if((n=dftReadReference(&dft, refname, &inputtype, &ref, tmp, verbose-3))<1) {
      fprintf(stderr, "Error in reading '%s': %s\n", refname, tmp);
      if(verbose>2) printf("dftReadReference()=%d\n", n);
      dftEmpty(&dft); dftEmpty(&input); return(6);
    }
    if(verbose>30) dftPrint(&dft); 
    if(n>1) {
      fprintf(stderr, "Warning: %s selected of %d reference regions.\n",
        dft.voi[ref].name, n);
      if(verbose>2)
        fprintf(stdout, "selected reference region := %s\n", dft.voi[ref].name);
    }
    if(inputtype==5) { // Reference region name was given
      refAdded=0; strcpy(refname, "");
    } else { // reference file was given; ref TACs may have to be deleted later
      refAdded=1;
    }
    if(verbose>15) dftPrint(&dft);
    if(verbose>1) printf("Reference region: %s\n", dft.voi[ref].name );
  }

  /* Allocate an extra TAC for the bootstrap */
  if(doBootstrap) {
    ret=dftAddmem(&dft, 1); if(ret) {
      fprintf(stderr, "Error: cannot allocate more memory.\n");
      dftEmpty(&dft); dftEmpty(&input); return(9);
    }
    strcpy(dft.voi[dft.voiNr].voiname, "BS");
    strcpy(dft.voi[dft.voiNr].name, "BS");
  }
  if(verbose>10) dftPrint(&dft);  


  /*
   *  Prepare the room for the results
   */
  if(verbose>1) printf("initializing result data\n");
  ret=res_allocate_with_dft(&res, &dft); if(ret!=0) {
    fprintf(stderr, "Error: cannot setup memory for results.\n");
    dftEmpty(&input); dftEmpty(&dft); return(7);
  }
  /* Copy titles & file names */
  tpcProgramName(argv[0], 1, 1, res.program, 256);
  strcpy(res.datafile, dfile);
  strcpy(res.plasmafile, pfile);
  strcpy(res.bloodfile, bfile);
  if(ref>=0) sprintf(res.refroi, "%s", dft.voi[ref].name);
  if(refname[0]) strcpy(res.reffile, refname);
  strcpy(res.fitmethod, "TGO");
  /* Constants */
  res.isweight=dft.isweight;
  if(fVb>=0.0) res.Vb=100.0*fVb;
  /* Set data range */
  sprintf(res.datarange, "%g - %g %s", 0.0, fitdur, dftTimeunit(dft.timeunit));
  res.datanr=fitframeNr;
  /* Set current time to results */
  res.time=time(NULL);
  /* Set parameter number, including also the extra "parameters"
     and the parameter names and units */
  res.parNr=9; 
  pi=0; strcpy(res.parname[pi], "K1"); strcpy(res.parunit[pi], "ml/(min*ml)");
  pi++; strcpy(res.parname[pi], "K1/k2"); strcpy(res.parunit[pi], "");
  pi++; strcpy(res.parname[pi], "k3"); strcpy(res.parunit[pi], "1/min");
  pi++; strcpy(res.parname[pi], "k3/k4"); strcpy(res.parunit[pi], "");
  pi++; strcpy(res.parname[pi], "k5"); strcpy(res.parunit[pi], "1/min");
  pi++; strcpy(res.parname[pi], "Vb"); strcpy(res.parunit[pi], "%");
  pi++; strcpy(res.parname[pi], "Ki"); strcpy(res.parunit[pi], "ml/(min*ml)");
  pi++; strcpy(res.parname[pi], "WSS"); strcpy(res.parunit[pi], "");
  pi++; strcpy(res.parname[pi], "AIC"); strcpy(res.parunit[pi], "");

  /* Set model function */
  if(model==0) func=func_sk5; else func=func_pk5;


  /*
   *  Fit at first the reference ROI(s), if required
   */
  if(ref>=0) for(ri=0; ri<dft.voiNr; ri++) if(dft.voi[ri].sw>0) {
    if(verbose>0)
      fprintf(stdout, "fitting %s as reference region\n", dft.voi[ri].name);
    /* Initiate values */
    petmeas=dft.voi[ri].y; petsim=dft.voi[ri].y2;
    /* Set constraints */
    for(pi=0; pi<parNr; pi++) { pmin[pi]=def_pmin[pi]; pmax[pi]=def_pmax[pi]; }
    for(pi=fittedparNr=0; pi<parNr; pi++) if(pmax[pi]>pmin[pi]) fittedparNr++;
    if(verbose>3) {
      printf("  ref_constraints :=");
      for(pi=0; pi<parNr; pi++) printf(" [%g,%g]", pmin[pi], pmax[pi]);
      printf("\n");
      printf("fittedparNr := %d\n", fittedparNr);
    }
    /* Fit */
    TGO_LOCAL_INSIDE=0;
    TGO_SQUARED_TRANSF=1;
    // n=fittedparNr; if(n>4) n=4; tgoNr=50+50*n; neighNr=8;
    tgoNr=50+30*fittedparNr;
    neighNr=6*fittedparNr;
    iterNr=0;
    ret=tgo(
      pmin, pmax, func, NULL, parNr, neighNr, &res.voi[ri].parameter[CM_WSS],
      res.voi[ri].parameter, tgoNr, iterNr, verbose-8);
    if(ret>0) {
      fprintf(stderr, "Error in optimization (%d).\n", ret);
      dftEmpty(&input); dftEmpty(&dft); resEmpty(&res); return(8);
    }
    /* Correct fitted parameters to match constraints like inside function */
    (void)modelCheckParameters(parNr, pmin, pmax, res.voi[ri].parameter,
                               res.voi[ri].parameter, NULL);
    res.voi[ri].parameter[CM_WSS]=wss_wo_penalty;
    /* If k3 or k3/k4 is zero, then both k3 and k3/k4, and k5, should be zero */
    if(res.voi[ri].parameter[CM_K3]<1.0E-20 || 
       res.voi[ri].parameter[CM_K3K4]<1.0E-20)
    { 
      res.voi[ri].parameter[CM_K3]=res.voi[ri].parameter[CM_K3K4]=0.0;
      res.voi[ri].parameter[CM_K5]=0.0;
    }
    if(verbose>1) {
      fprintf(stdout, "  K1/k2 := %g\n", res.voi[ri].parameter[CM_K1K2]);
      //fprintf(stdout, "  k3/k4 := %g\n", res.voi[ri].parameter[CM_K3K4]);
    }
    if(verbose>4) {
      printf("Original and fitted TACs:\n");
      for(fi=0; fi<fitframeNr; fi++)
        printf("  %8.3f  %9.3f   %9.3f\n", dft.x[fi], petmeas[fi], petsim[fi]);
      printf("Original and fitted TACs:\n");
    }
    /* Fix K1/k2 to non-reference regions, but only after all reference
       regions are fitted */
    if(ri==ref) {
      fk1k2=res.voi[ri].parameter[CM_K1K2];
      if(verbose>2) {fprintf(stdout, "  fixed K1/k2 := %g\n", fk1k2);}
    }
    /* Bootstrap */
    if(doBootstrap) {
      if(verbose>2) printf("  bootstrapping\n");
      /* bootstrap changes measured and simulated data, therefore use copies */
      petmeas=dft.voi[dft.voiNr].y2; petsim=dft.voi[dft.voiNr].y3;
      if(doSD) sd=res.voi[ri].sd; else sd=NULL;
      if(doCL) {cl1=res.voi[ri].cl1; cl2=res.voi[ri].cl2;} else cl1=cl2=NULL;
      ret=bootstrap(
        0, cl1, cl2, sd,
        res.voi[ri].parameter, pmin, pmax, fitframeNr,
        // measured original TAC, not modified
        dft.voi[ri].y,
        // fitted TAC, not modified
	      dft.voi[ri].y2,
        // tissue TAC noisy data is written to be used by objf
	      petmeas, 
	      parNr, dft.w, func, tmp, verbose-4
      );
      if(ret) {
        fprintf(stderr, "Error in bootstrap: %s\n", tmp);
        for(pi=0; pi<parNr; pi++) {
          if(doSD) sd[pi]=nan(""); 
          if(doCL) cl1[pi]=cl2[pi]=nan("");
        }
      }
    }
    /* Calculate AIC based on nr of parameters that actually are fitted */
    for(pi=n=0; pi<parNr; pi++) if(pmax[pi]>pmin[pi]) n++;
    if(verbose>2) printf("nr_of_fitted_parameters := %d\n", n);
    for(fi=m=0; fi<fitframeNr; fi++) if(dft.w[fi]>0.0) m++;
    if(verbose>2) printf("nr_of_fitted_samples := %d\n", m);
    res.voi[ri].parameter[CM_AIC]=aicSS(res.voi[ri].parameter[CM_WSS], m, n);

    /* Calculate Ki */
    K1=res.voi[ri].parameter[CM_K1];
    k2=K1/res.voi[ri].parameter[CM_K1K2];
    k3=res.voi[ri].parameter[CM_K3];
    k4=k3/res.voi[ri].parameter[CM_K3K4];
    k5=res.voi[ri].parameter[CM_K5];
    if(model==0) {
      f=k2*k4+k2*k5+k3*k5;
      if(f<=1.0E-100) res.voi[ri].parameter[CM_KI]=0.0;
      else res.voi[ri].parameter[CM_KI]=K1*k3*k5/f;
    } else {
      f=k2+k5;
      if(f<=1.0E-100) res.voi[ri].parameter[CM_KI]=0.0;
      else res.voi[ri].parameter[CM_KI]=K1*k5/f;
    }

    /* Convert Vb fraction to percentage */
    res.voi[ri].parameter[CM_VB]*=100.;
    if(doSD) res.voi[ri].sd[CM_VB]*=100.;
    if(doCL) {res.voi[ri].cl1[CM_VB]*=100.; res.voi[ri].cl2[CM_VB]*=100.;}

  } // next reference region


  /*
   *  Fit other than reference regions
   */
  if(verbose>0) {fprintf(stdout, "fitting regional TACs: "); fflush(stdout);}
  if(verbose>1) printf("\n");
  /* Set K1/k2 constraint based on above fitted reference region, if required */
  if(fk1k2>0.0) def_pmin[CM_K1K2]=def_pmax[CM_K1K2]=fk1k2;
  /* One ROI at a time */
  for(ri=0; ri<dft.voiNr; ri++) if(dft.voi[ri].sw==0) {
    if(verbose>2) printf("\n  %d %s:\n", ri, dft.voi[ri].name);

    /* Initiate values */
    petmeas=dft.voi[ri].y; petsim=dft.voi[ri].y2;

    /* Set constraints */
    for(pi=0; pi<parNr; pi++) { pmin[pi]=def_pmin[pi]; pmax[pi]=def_pmax[pi]; }
    for(pi=fittedparNr=0; pi<parNr; pi++) if(pmax[pi]>pmin[pi]) fittedparNr++;
    if(verbose>3) {
      printf("  constraints :=");
      for(pi=0; pi<parNr; pi++) printf(" [%g,%g]", pmin[pi], pmax[pi]);
      printf("\n");
      printf("fittedparNr := %d\n", fittedparNr);
    }

    /* Fit */
    TGO_LOCAL_INSIDE=0;
    TGO_SQUARED_TRANSF=1;
    // n=fittedparNr; if(n>4) n=4; tgoNr=50*50*n; neighNr=8; iterNr=0;
    tgoNr=50+30*fittedparNr;
    neighNr=6*fittedparNr;
    iterNr=0;
    ret=tgo(
      pmin, pmax, func, NULL, parNr, neighNr, &res.voi[ri].parameter[CM_WSS],
      res.voi[ri].parameter, tgoNr, iterNr, verbose-8);
    if(ret>0) {
      fprintf(stderr, "\nError in optimization (%d).\n", ret);
      dftEmpty(&input); dftEmpty(&dft); resEmpty(&res); return(8);
    }
    /* Correct fitted parameters to match constraints like inside function */
    (void)modelCheckParameters(parNr, pmin, pmax, res.voi[ri].parameter,
                               res.voi[ri].parameter, NULL);
    res.voi[ri].parameter[CM_WSS]=wss_wo_penalty;
    /* If k3 or k3/k4 is zero, then both k3 and k3/k4, and k5, should be zero */
    if(res.voi[ri].parameter[CM_K3]<1.0E-20 || 
       res.voi[ri].parameter[CM_K3K4]<1.0E-20) 
    {
      res.voi[ri].parameter[CM_K3]=res.voi[ri].parameter[CM_K3K4]=0.0;
      res.voi[ri].parameter[CM_K5]=0.0;
    }

    /* Bootstrap */
    if(doBootstrap) {
      if(verbose>2) printf("  bootstrapping\n");
      /* bootstrap changes measured and simulated data, therefore use copies */
      petmeas=dft.voi[dft.voiNr].y2; petsim=dft.voi[dft.voiNr].y3;
      if(doSD) sd=res.voi[ri].sd; else sd=NULL;
      if(doCL) {cl1=res.voi[ri].cl1; cl2=res.voi[ri].cl2;} else cl1=cl2=NULL;
      ret=bootstrap(
        0, cl1, cl2, sd,
        res.voi[ri].parameter, pmin, pmax, fitframeNr,
        // measured original TAC, not modified
        dft.voi[ri].y,
        // fitted TAC, not modified
        dft.voi[ri].y2,
        // tissue TAC noisy data is written to be used by objf
        petmeas, 
        parNr, dft.w, func, tmp, verbose-4
      );
      if(ret) {
        fprintf(stderr, "Error in bootstrap: %s\n", tmp);
        for(pi=0; pi<parNr; pi++) {
          if(doSD) sd[pi]=nan(""); 
          if(doCL) cl1[pi]=cl2[pi]=nan("");
        }
      }
    }

    /* Calculate Ki */
    K1=res.voi[ri].parameter[CM_K1];
    k2=K1/res.voi[ri].parameter[CM_K1K2];
    k3=res.voi[ri].parameter[CM_K3];
    k4=k3/res.voi[ri].parameter[CM_K3K4];
    k5=res.voi[ri].parameter[CM_K5];
    if(model==0) {
      f=k2*k4+k2*k5+k3*k5;
      if(f<=1.0E-100) res.voi[ri].parameter[CM_KI]=0.0;
      else res.voi[ri].parameter[CM_KI]=K1*k3*k5/f;
    } else {
      f=k2+k5;
      if(f<=1.0E-100) res.voi[ri].parameter[CM_KI]=0.0;
      else res.voi[ri].parameter[CM_KI]=K1*k5/f;
    }

    /* Calculate AIC, based on nr of parameters that actually are fitted */
    for(pi=n=0; pi<parNr; pi++) if(pmax[pi]>pmin[pi]) n++;
    if(verbose>2) printf("nr_of_fitted_parameters := %d\n", n);
    for(fi=m=0; fi<fitframeNr; fi++) if(dft.w[fi]>0.0) m++;
    if(verbose>2) printf("nr_of_fitted_samples := %d\n", m);
    res.voi[ri].parameter[CM_AIC]=aicSS(res.voi[ri].parameter[CM_WSS], m, n);

    /* Convert Vb fraction to percentage */
    res.voi[ri].parameter[CM_VB]*=100.;
    if(doSD) res.voi[ri].sd[CM_VB]*=100.;
    if(doCL) {res.voi[ri].cl1[CM_VB]*=100.; res.voi[ri].cl2[CM_VB]*=100.;}

    /* done with this region */
    if(dft.voiNr>2 && verbose==1) {fprintf(stdout, "."); fflush(stdout);}

  } /* next region */
  if(verbose>0) {fprintf(stdout, "\n"); fflush(stdout);}


  /*
   *  Print results on screen
   */
  if(verbose>0) {resPrint(&res); fprintf(stdout, "\n");}


  /*
   *  Save results
   */
  if(verbose>1) printf("saving results\n");
  ret=resWrite(&res, rfile, verbose-3);
  if(ret) {
    fprintf(stderr, "Error in writing '%s': %s\n", rfile, reserrmsg);
    dftEmpty(&dft); dftEmpty(&input); resEmpty(&res);
    return(11);
  }
  if(verbose>0) fprintf(stdout, "Model parameters written in %s\n", rfile);


  /*
   *  Saving and/or plotting of fitted TACs
   */
  if(svgfile[0] || ffile[0]) {

    /* Create a DFT containing fitted TACs */
    char tmp[64];
    DFT dft2;
    dftInit(&dft2); ret=dftdup(&dft, &dft2);
    if(ret) {
      fprintf(stderr, "Error: cannot save fitted curves.\n");
      dftEmpty(&dft); dftEmpty(&input); resEmpty(&res);
      return(21);
    }
    for(ri=0; ri<dft.voiNr; ri++) for(fi=0; fi<fitframeNr; fi++)
      dft2.voi[ri].y[fi]=dft2.voi[ri].y2[fi];
    dft2.frameNr=fitframeNr;

    /* Save SVG plot of fitted and original data */
    if(svgfile[0]) {
      if(verbose>1) printf("saving SVG plot\n");
      sprintf(tmp, "K1-k5 fit: ");
      if(strlen(dft.studynr)>0) strcat(tmp, dft.studynr);
      ret=plot_fitrange_svg(&dft, &dft2, tmp, 0.0, 1.02*dft.x[fitframeNr-1],
                            0.0, nan(""), svgfile, verbose-8);
      if(ret) {
        fprintf(stderr, "Error (%d) in writing '%s'.\n", ret, svgfile);
        dftEmpty(&dft2); dftEmpty(&dft); dftEmpty(&input); resEmpty(&res);
        return(30+ret);
      }
      if(verbose>0) printf("Plots written in %s\n", svgfile);
    }

    /* Delete reference region(s) from the data */
    if(refAdded!=0) {
      for(ri=dft2.voiNr-1; ri>=0; ri--) if(dft2.voi[ri].sw!=0)
        dftDelete(&dft2, ri);
    }

    /* Save fitted TACs */
    if(ffile[0]) {
      if(verbose>1) printf("saving fitted curves\n");
      tpcProgramName(argv[0], 1, 0, tmp, 128);
      sprintf(dft2.comments, "# program := %s\n", tmp);    
      if(dftWrite(&dft2, ffile)) {
        fprintf(stderr, "Error in writing '%s': %s\n", ffile, dfterrmsg);
        dftEmpty(&dft2); dftEmpty(&dft); dftEmpty(&input); resEmpty(&res);
        return(22);
      }
      if(verbose>0) printf("Fitted TACs written in %s\n", ffile);
    }

    dftEmpty(&dft2);
  }

  dftEmpty(&dft); dftEmpty(&input); resEmpty(&res);
  return(0);
}
/*****************************************************************************/

/*****************************************************************************
 *
 *  Functions to be minimized
 *
 *****************************************************************************/

/* Compartments in series */
double func_sk5(int parNr, double *p, void *fdata)
{
  int fi, ret;
  double Vb, k2, k4, d, wss=0.0;
  double pa[MAX_PARAMETERS], penalty=1.0;

  /* Check parameters against the constraints */
  ret=modelCheckParameters(parNr, pmin, pmax, p, pa, &penalty);
  if(fdata) {}
  /* Calculate individual parameters */
  k2=pa[CM_K1]/pa[CM_K1K2];
  k4=pa[CM_K3]/pa[CM_K3K4];
  Vb=pa[CM_VB];

  /* Simulate the tissue PET TAC */
  ret=simC3vs(
    input.x, input.voi[0].y, input.voi[1].y, input.frameNr,
    pa[CM_K1], k2, pa[CM_K3], k4, pa[CM_K5], 0.0, 0.0, Vb, 1.0,
    input.voi[0].y2, NULL, NULL, NULL, NULL, NULL);
  if(ret) {
    printf("error %d in simulation\n", ret);
    return(nan(""));
  }

  /* Interpolate & integrate to measured PET frames */
  if(dft.timetype==DFT_TIME_STARTEND)
    ret=interpolate4pet(
      input.x, input.voi[0].y2, input.frameNr,
      dft.x1, dft.x2, petsim, NULL, NULL, fitframeNr);
  else
    ret=interpolate(
      input.x, input.voi[0].y2, input.frameNr,
      dft.x, petsim, NULL, NULL, fitframeNr);
  if(ret) {
    printf("error %d in interpolation\n", ret);
    return(nan(""));
  }

  /* Calculate error */
  for(fi=0, wss=0.0; fi<fitframeNr; fi++) if(dft.w[fi]>0.0) {
    d=petmeas[fi]-petsim[fi]; 
    wss+=dft.w[fi]*d*d;
  }
  wss_wo_penalty=wss;
  wss*=penalty;

  return(wss);
}


/* Compartments in parallel */
double func_pk5(int parNr, double *p, void *fdata)
{
  int fi, ret;
  double Vb, k2, k4, d, wss=0.0;
  double pa[MAX_PARAMETERS], penalty=1.0;


  /* Check parameters against the constraints */
  ret=modelCheckParameters(parNr, pmin, pmax, p, pa, &penalty);
  if(fdata) {}
  /* Calculate individual parameters */
  k2=pa[CM_K1]/pa[CM_K1K2];
  k4=pa[CM_K3]/pa[CM_K3K4];
  Vb=pa[CM_VB];

  /* Simulate the tissue PET TAC */
  ret=simC3vp(
    input.x, input.voi[0].y, input.voi[1].y, input.frameNr,
    pa[CM_K1], k2, pa[CM_K3], k4, pa[CM_K5], 0.0, 0.0, Vb, 1.0,
    input.voi[0].y2, NULL, NULL, NULL, NULL, NULL);
  if(ret) {
    printf("error %d in simulation\n", ret);
    return(nan(""));
  }

  /* Interpolate & integrate to measured PET frames */
  if(dft.timetype==DFT_TIME_STARTEND)
    ret=interpolate4pet(
      input.x, input.voi[0].y2, input.frameNr,
      dft.x1, dft.x2, petsim, NULL, NULL, fitframeNr);
  else
    ret=interpolate(
      input.x, input.voi[0].y2, input.frameNr,
      dft.x, petsim, NULL, NULL, fitframeNr);
  if(ret) {
    printf("error %d in interpolation\n", ret);
    return(nan(""));
  }

  /* Calculate error */
  for(fi=0, wss=0.0; fi<fitframeNr; fi++) if(dft.w[fi]>0.0) {
    d=petmeas[fi]-petsim[fi]; 
    wss+=dft.w[fi]*d*d;
  }
  wss_wo_penalty=wss;
  wss*=penalty;

  return(wss);
}
/*****************************************************************************/

/*****************************************************************************/
/** Read specified result file, search K1/k2 parameter and compute the
 *  mean and/or median of regional K1/k2 values
 *  @return Returns 0 if successful, <>0 otherwise.
 */ 
int resK1k2Median(
  /** Name of result file from which mean and median K1/k2 is calculated */
  char *filename,
  /** Pointer to median K1/k2; NULL if not needed */
  double *median,
  /** Pointer to mean K1/k2; NULL if not needed */
  double *mean,
  /** Pointer to string for error message; NULL if not needed */
  char *msg
) {
  RES res;
  int ri, col, n;
  char *cptr;
  double *lista;

  /* Read result file */
  resInit(&res);
  if(resRead(filename, &res, 0)!=0) {
    if(msg!=NULL) sprintf(msg, "cannot read %s: %s", filename, reserrmsg);
    return(1);
  }
  /* Which column contains K1/k2 */
  col=-1; ri=0; cptr=strtok(res.titleline, " \t");
  while(cptr!=NULL && col<0) {
    if(strcasecmp(cptr, "K1/k2")==0) col=ri;
    else if(strcasecmp(cptr, "K1k2")==0) col=ri;
    ri++; cptr=strtok(NULL, " \t");
  }
  if(ri<0) {
    if(msg!=NULL) sprintf(msg, "K1/k2 not found in %s", filename);
    resEmpty(&res); return(3);
  }
  n=res.voiNr;
  /* Calculate mean and median */
  lista=(double*)malloc(n*sizeof(double));
  for(ri=0; ri<n; ri++) lista[ri]=res.voi[ri].parameter[col];
  resEmpty(&res);
  if(mean!=NULL) {
    *mean=dmean(lista, n, NULL);
    if(*mean<=0.0) {
      if(msg!=NULL) sprintf(msg, "invalid K1/k2 mean in %s", filename);
      free(lista); return(6);
    }
  }
  if(median!=NULL) {
    *median=dmedian(lista, n);
    if(*median<=0.0) {
      if(msg!=NULL) sprintf(msg, "invalid K1/k2 median in %s", filename);
      free(lista); return(5);
    }
  }
  free(lista);
  return 0;
}
/*****************************************************************************/

/*****************************************************************************/
/// @endcond
