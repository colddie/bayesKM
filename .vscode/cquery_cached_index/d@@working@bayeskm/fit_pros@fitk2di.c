/** @file fitk2di.c
 *  @brief Estimates the parameters of 1-tissue compartmental model
           with dual input.
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
static char *info[] = {
  "Non-linear fitting of dual input compartment model, with one tissue",
  "compartment for each input (parent tracer and its labeled metabolite):",
  " ",
  "  _____   K1p   _____   ",
  " | Cap | ----> | Ctp |  ",
  " |_____| <---- |_____|  ",
  "          k2p     |     ",
  "                km|     ",
  "                  v     ",
  "  _____   K1m   _____   ",
  " | Cam | ----> | Ctm |  ",
  " |_____| <---- |_____|  ",
  "          k2m           ",
  " ",
  "Sample times must be in minutes.",
  " ",
  "Usage: @P [Options] ptacfile mtacfile btacfile ttacfile endtime resultfile",
  " ",
  "Options:",
  " -lim[=<filename>]",
  "     Specify the constraints for model parameters;",
  "     This file with default values can be created by giving this",
  "     option as the only command-line argument to this program.",
  "     Without filename the default values are printed on screen.",
  " -SD=<y|N>",
  "     Standard deviations are calculated and saved in results (y),",
  "     or not calculated (N, default).",
  "     Program runs a lot faster if SD and CL are not calculated.",
  " -CL=<y|N>",
  "     95% Confidence limits are calculated and saved in results (y), or",
  "     not calculated (N, default).",
  " -Vb=<Vb(%)>",
  "     Enter a fixed Vb; fitted by default.",
  "     If Vb (vascular blood volume) is pre-corrected or to be ignored, set",
  "     it to 0; btacfile can be set to 'none'.",
  " -ref=<Reference region name or filename>",
  "     Specified reference region is fitted using different set of model",
  "     parameter constraints; not necessary if reference region is given",
  "     with one of the following options -BPnd, -BPp, or -DVR.",
  " -<BPnd|BPp|DVR>=<Reference region name or filename>",
  "     Optional reference region is used to calculate BPnd, BPp, or DVR;",
  "     BPnd=DVroi/DVref-1, BPp=DVroi-DVref, and DVR=DVroi/DVref",
  " -refVfm=refVfp",
  "     In reference region Vfm is set to equal Vfp=1.",
  " -mc=<Filename>",
  "     Fit-based metabolite corrected regional TACs are written in the file.",
  " -fit=<Filename>",
  "     Fitted regional TACs are written in the file.",
  " -svg=<Filename>",
  "     Fitted and measured TACs are plotted in specified SVG file.",
  " -stdoptions", // List standard options like --help, -v, etc
  " ",
  "Example 1: fitting with default settings",
  "     @P ia919apc.kbq ia919apm.kbq ia919ab.kbq ia919.dft 60 a919k2di.res",
  " ",
  "Example 2: Vb is constrained to 0%; DVRs are calculated by dividing DVs",
  "by the DV of region 'cer'",
  "     @P -Vb=0 -R=cer p25apc.bld p25apm.bld none p25.tac 60 p25k2di.res",
  " ",
  "See also: fitk2, logan, fitk4, p2t_di, dftweigh, dftcbv",
  " ",
  "Keywords: TAC, modelling, distribution volume, reversible uptake, dual-input",
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
const int parNr=6, parVb=5;
DFT input, dft, fit;
double *petmeas, *petsim;
double fVb=-1.0;
int is_this_ref=0;
int fixed_ref_Vfm_eq_Vfp=0; // If <>0 then Ref region Vfm = Ref region Vfp
double pmin[MAX_PARAMETERS], pmax[MAX_PARAMETERS];
int fitframeNr;
// Parameter names in the limit file and in actual fitting
static char *parName[] = {"K1p", "Vfp", "R1m", "Vfm", "km", "Vb", 0};
/*****************************************************************************/
/* Local functions */
double func1TCMdi(int parNr, double *p, void*);
/*****************************************************************************/

/*****************************************************************************/
/**
 *  Main
 */
int main(int argc, char **argv)
{
  int        ai, help=0, version=0, verbose=1;
  int        ri, fi, pi, n, ret;
  int        ref=-1, refAdded=0, inputtype;
  int        bp_type=0; // 0=no, 1=DVR, 2=BPnd, 3=BPp
  char       dfile[FILENAME_MAX], bfile[FILENAME_MAX],
             pfile[FILENAME_MAX], mfile[FILENAME_MAX],
             rfile[FILENAME_MAX], ffile[FILENAME_MAX],
             mcfile[FILENAME_MAX],
             limfile[FILENAME_MAX], svgfile[FILENAME_MAX];
  char       refname[FILENAME_MAX];
  char      *cptr, tmp[FILENAME_MAX];
  double     f, fitdur, wss;
  RES        res;
  int        doBootstrap=0, doSD=0, doCL=0; // 0=no, 1=yes
  double    *sd, *cl1, *cl2;
  double     def_pmin[MAX_PARAMETERS], def_pmax[MAX_PARAMETERS];
  double     def_pmin_ref[MAX_PARAMETERS], def_pmax_ref[MAX_PARAMETERS];
  int        tgoNr=0, neighNr=0, iterNr=0, fittedparNr=0;


#ifdef MINGW
  // Use Unix/Linux default of two-digit exponents in MinGW on Windows
  _set_output_format(_TWO_DIGIT_EXPONENT);
#endif

  /* Set parameter initial values and constraints */
  /* K1p         */  def_pmin[0]=0.0;  def_pmax[0]=10.0;
  /* Vfp=K1p/k2p */  def_pmin[1]=0.0;  def_pmax[1]=500.0;
  /* R1m=K1m/K1p */  def_pmin[2]=0.0;  def_pmax[2]=10.0;
  /* Vfm=K1m/k2m */  def_pmin[3]=0.0;  def_pmax[3]=10.0;
  /* km          */  def_pmin[4]=0.0;  def_pmax[4]=0.0;
  /* Vb          */  def_pmin[5]=0.0;  def_pmax[5]=0.10;
  /* Same for possible reference region */
  /* K1p         */  def_pmin_ref[0]=0.0;  def_pmax_ref[0]=10.0;
  /* Vfp=K1p/k2p */  def_pmin_ref[1]=0.0;  def_pmax_ref[1]=2.0;
  /* R1m=K1m/K1p */  def_pmin_ref[2]=0.0;  def_pmax_ref[2]=10.0;
  /* Vfm=K1m/k2m */  def_pmin_ref[3]=0.0;  def_pmax_ref[3]=2.0;
  /* km          */  def_pmin_ref[4]=0.0;  def_pmax_ref[4]=0.0;
  /* Vb          */  def_pmin_ref[5]=0.0;  def_pmax_ref[5]=0.10;


  /*
   *  Get arguments
   */
  if(argc==1) {tpcPrintUsage(argv[0], info, stderr); return(1);}
  dfile[0]=pfile[0]=mfile[0]=bfile[0]=refname[0]=(char)0;
  rfile[0]=ffile[0]=limfile[0]=svgfile[0]=mcfile[0]=(char)0;
  fitdur=fVb=-1.0; 
  resInit(&res); dftInit(&input); dftInit(&dft); dftInit(&fit);
  /* Get options first, because it may affect what arguments are read */
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
    } else if(strncasecmp(cptr, "Vb=", 3)==0 && strlen(cptr)>3) {
      fVb=0.01*atof_dpi(cptr+3);
      if(fVb>=0.0 && fVb<1.0) {
        if(fVb<0.01) fprintf(stderr, "Warning: Vb was set to %g%%\n", 100.*fVb);
        def_pmin[parVb]=def_pmax[parVb]=fVb;
        def_pmin_ref[parVb]=def_pmax_ref[parVb]=fVb;
        continue;
      }
      fVb=-1.0;
    } else if(strncasecmp(cptr, "REF=", 4)==0) {
      strlcpy(refname, cptr+4, FILENAME_MAX); if(strlen(refname)>0.0) continue;
    } else if(strncasecmp(cptr, "DVR=", 4)==0) {
      bp_type=1; strlcpy(refname, cptr+4, FILENAME_MAX); 
      if(strlen(refname)>0.0) continue;
    } else if(strncasecmp(cptr, "BPnd=", 5)==0) {
      bp_type=2; strlcpy(refname, cptr+5, FILENAME_MAX); 
      if(strlen(refname)>0.0) continue;
    } else if(strncasecmp(cptr, "BPp=", 4)==0) {
      bp_type=3; strlcpy(refname, cptr+4, FILENAME_MAX); 
      if(strlen(refname)>0.0) continue;
    } else if(strcasecmp(cptr, "refVfm=refVfp")==0) {
      fixed_ref_Vfm_eq_Vfp=1; continue;
    } else if(strncasecmp(cptr, "MC=", 3)==0) {
      strlcpy(mcfile, cptr+3, FILENAME_MAX); if(strlen(mcfile)>0) continue;
    } else if(strncasecmp(cptr, "FIT=", 4)==0) {
      strlcpy(ffile, cptr+4, FILENAME_MAX); if(strlen(ffile)>0) continue;
    } else if(strncasecmp(cptr, "SVG=", 4)==0) {
      strlcpy(svgfile, cptr+4, FILENAME_MAX); if(strlen(svgfile)>0) continue;
    }
    fprintf(stderr, "Error: invalid option '%s'.\n", argv[ai]);
    return(1);
  }

  
  /* Print help or version? */
  if(help==2) {tpcHtmlUsage(argv[0], info, ""); return(0);}
  if(help) {tpcPrintUsage(argv[0], info, stdout); return(0);}
  if(version) {tpcPrintBuild(argv[0], stdout); return(0);}

  /* Process other arguments, starting from the first non-option */
  for(ai=1; ai<argc; ai++) if(*argv[ai]!='-') {
    if(!pfile[0]) {strlcpy(pfile, argv[ai], FILENAME_MAX); continue;}
    else if(!mfile[0]) {strlcpy(mfile, argv[ai], FILENAME_MAX); continue;}
    else if(!bfile[0]) {strlcpy(bfile, argv[ai], FILENAME_MAX); continue;}
    else if(!dfile[0]) {strlcpy(dfile, argv[ai], FILENAME_MAX); continue;}
    else if(fitdur<0) {fitdur=atof_dpi(argv[ai]); if(fitdur>0.0) continue;}
    else if(!rfile[0]) {strlcpy(rfile, argv[ai], FILENAME_MAX); continue;}
    /* we should never get this far */
    fprintf(stderr, "Error: too many arguments: '%s'.\n", argv[ai]);
    return(1);
  }
  if(doSD || doCL) doBootstrap=1; else doBootstrap=0;
  /* User may not have blood file, and has entered 'none' instead */
  if(strcasecmp(bfile, "NONE")==0) {strcpy(bfile, ""); fVb=0.0;}

  /* If only filename for initial values was given, then write one
     with default contents and exit */
  if(limfile[0] && !pfile[0]) {
    /* Check that initial value file does not exist */
    if(strcasecmp(limfile, "stdout")!=0 && access(limfile, 0) != -1) {
      fprintf(stderr, "Error: parameter constraint file %s exists.\n", limfile);
      return(9);
    }
    if(verbose>1) printf("writing parameter constraints file\n");
    IFT ift; iftInit(&ift);
    /* Create parameter file */
    for(pi=0; pi<parNr; pi++) {
      sprintf(tmp, "%s_lower", parName[pi]);
      iftPutDouble(&ift, tmp, def_pmin[pi], NULL);
      sprintf(tmp, "%s_upper", parName[pi]); 
      iftPutDouble(&ift, tmp, def_pmax[pi], NULL);
    }
    for(pi=0; pi<parNr; pi++) {
      sprintf(tmp, "ref_%s_lower", parName[pi]);
      iftPutDouble(&ift, tmp, def_pmin_ref[pi], NULL);
      sprintf(tmp, "ref_%s_upper", parName[pi]);
      iftPutDouble(&ift, tmp, def_pmax_ref[pi], NULL);
    }
    if((ret=iftWrite(&ift, limfile))!=0) {
      fprintf(stderr, "Error in writing '%s': %s\n", limfile, ift.status);
      iftEmpty(&ift); return(9);
    }
    if(strcasecmp(limfile, "stdout")!=0) fprintf(stdout, 
                   "Parameter file %s with initial values written.\n", limfile);
    iftEmpty(&ift); return(0);
  }

  /* Did we get all the information from user that we need? */
  if(fitdur==0) fitdur=1.0E+100; 
  else if(fitdur<0) {tpcPrintUsage(argv[0], info, stderr); return(1);}
  if(!rfile[0]) {
    fprintf(stderr, "Error: missing command-line argument; use option --help\n");
    return(1);
  }

  /* In verbose mode print arguments and options */
  if(verbose>1) {
    printf("pfile := %s\n", pfile);
    printf("mfile := %s\n", mfile);
    printf("dfile :=%s\n", dfile);
    printf("rfile := %s\n", rfile);
    printf("mcfile := %s\n", mcfile);
    printf("ffile := %s\n", ffile);
    printf("svgfile := %s\n", svgfile);
    printf("limfile := %s\n", limfile);
    printf("bp_type := %d\n", bp_type);
    printf("refname := %s\n", refname);
    printf("fitdur := %g\n", fitdur);
    printf("doBootstrap := %d\n", doBootstrap);
    printf("doSD := %d\n", doSD);
    printf("doCL := %d\n", doCL);
    printf("fixed_ref_Vfm_eq_Vfp := %d\n", fixed_ref_Vfm_eq_Vfp);
  }


  /*
   *  Read model parameter initial values and upper and lower limits
   *  if file for that was given
   */
  if(limfile[0]) {
    IFT ift; iftInit(&ift);
    double v;
    if(verbose>1) printf("reading %s\n", limfile);
    ret=iftRead(&ift, limfile, 1);
    if(ret) {
      fprintf(stderr, "Error in reading '%s': %s\n", limfile, ift.status);
      return(9);
    }
    if(verbose>10) iftWrite(&ift, "stdout");
    int n=0;
    for(pi=0; pi<parNr; pi++) {
      sprintf(tmp, "%s_lower", parName[pi]);
      if(iftGetDoubleValue(&ift, 0, tmp, &v)>=0) {def_pmin[pi]=v; n++;}
      sprintf(tmp, "%s_upper", parName[pi]);
      if(iftGetDoubleValue(&ift, 0, tmp, &v)>=0) {def_pmax[pi]=v; n++;}
    }
    for(pi=0; pi<parNr; pi++) {
      sprintf(tmp, "ref_%s_lower", parName[pi]);
      if(iftGetDoubleValue(&ift, 0, tmp, &v)>=0) {def_pmin_ref[pi]=v; n++;}
      sprintf(tmp, "ref_%s_upper", parName[pi]);
      if(iftGetDoubleValue(&ift, 0, tmp, &v)>=0) {def_pmax_ref[pi]=v; n++;}
    }
    iftEmpty(&ift);
    if(n==0) {fprintf(stderr, "Error: invalid parameter file.\n"); return(9);}
  }
  /* Refine constraints based on command-line options */
  if(fixed_ref_Vfm_eq_Vfp!=0) {
    def_pmax_ref[3]=def_pmin_ref[3]=0.0;
  }
  /* Check that these are ok */
  for(pi=n=0, ret=0; pi<parNr; pi++) {
    if(def_pmin[pi]<0.0) ret++;
    if(def_pmax[pi]<def_pmin[pi]) ret++;
    if(def_pmax[pi]>def_pmin[pi]) n++;
  }
  if(ret>0) {
    fprintf(stderr, "Error: invalid parameter constraints.\n");
    return(9);
  }
  if(n==0) {
    fprintf(stderr, "Error: no model parameters left free for fitting.\n");
    return(9);
  }
  /* The same for reference region */
  for(pi=n=0, ret=0; pi<parNr; pi++) {
    if(def_pmin_ref[pi]<0.0) ret++;
    if(def_pmax_ref[pi]<def_pmin_ref[pi]) ret++;
    if(def_pmax_ref[pi]>def_pmin_ref[pi]) n++;
  }
  if(ret>0) {
    fprintf(stderr, "Error: invalid reference region parameter constraints.\n");
    return(9);
  }
  if(n==0) {
    fprintf(stderr, "Error: no ref model parameters left free for fitting.\n");
    return(9);
  }

  /* Fixed/fitted Vb */
  if(fVb>=0.0) 
    def_pmin[parVb]=def_pmax[parVb]=def_pmin_ref[parVb]=def_pmax_ref[parVb]=fVb;
  if(def_pmin[parVb]==def_pmax[parVb]) fVb=def_pmin[parVb];
  if(fVb==0.0) strcpy(bfile, "");
  if(verbose>1) {
    printf("bfile := %s\n", bfile);
    //if(fVb>=0.0) printf("fVb := %g\n", fVb);
  }


  /*
   *  Read tissue and input data
   */
  if(verbose>1) printf("reading tissue and input data\n");
  ret=dftReadModelingData(dfile, pfile, mfile, bfile, &fitdur, 
                        &fitframeNr, &dft, &input, stdout, verbose-2, tmp);
  if(ret!=0) {
    fprintf(stderr, "Error: %s\n", tmp);
    dftEmpty(&dft); dftEmpty(&input); return(2);
  }
  if(fitframeNr<parNr+1 || input.frameNr<parNr+1) {
    fprintf(stderr, "Error: too few samples in specified fit duration.\n");
    dftEmpty(&input); dftEmpty(&dft); return(2);
  }
  if(input.voiNr<2) {
    fprintf(stderr, "Error: valid plasma TACs must be provided.\n");
    dftEmpty(&input); dftEmpty(&dft); return(2);
  }
  /* If there is no blood TAC, then create a zero blood TAC */
  if(input.voiNr<3) {
    if(verbose>2) printf("setting blood tac to zero\n");
    ret=dftAddmem(&input, 1);
    if(ret) {
      fprintf(stderr, "Error: cannot allocate more memory.\n");
      dftEmpty(&dft); dftEmpty(&input); return(3);
    }
    strcpy(input.voi[2].voiname, "blood");
    for(fi=0; fi<input.frameNr; fi++) input.voi[2].y[fi]=0.0;
    input.voiNr=3;
    /* and make sure that Vb=0 */
    def_pmin[parVb]=def_pmax[parVb]=fVb=0.0;
    def_pmin_ref[parVb]=def_pmax_ref[parVb]=0.0;
  }
  if(verbose>1) {
    if(fVb>=0.0) printf("fVb := %g\n", fVb);
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
    if((n=dftReadReference(&dft, refname, &inputtype, &ref, tmp, verbose-3))<1) 
    {
      fprintf(stderr, "Error in reading '%s': %s\n", refname, tmp);
      if(verbose>2) printf("dftReadReference()=%d\n", n);
      dftEmpty(&dft); dftEmpty(&input); return(6);
    }
    if(verbose>30) dftPrint(&dft); 
    if(n>1)
      fprintf(stderr, "Warning: %s selected of %d reference regions.\n",
        dft.voi[ref].name, n);
    if(verbose>1)
      fprintf(stdout, "selected reference region := %s\n", dft.voi[ref].name);
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
  strcpy(res.plasmafile2, mfile);
  strcpy(res.bloodfile, bfile);
  if(ref>=0) sprintf(res.refroi, "%s", dft.voi[ref].name);
  if(refname[0]) strcpy(res.reffile, refname);
  strcpy(res.fitmethod, "TGO");
  /* Constants */
  if(fVb>=0.0) res.Vb=100.0*fVb; else res.Vb=-1.0;
  res.isweight=dft.isweight;
  /* Set data range */
  sprintf(res.datarange, "%g - %g %s", 0.0, fitdur, dftTimeunit(dft.timeunit));
  res.datanr=fitframeNr;
  /* Set current time to results */
  res.time=time(NULL);
  /* Set parameter number, including also the extra "parameters"
     and the parameter names and units */
  res.parNr=parNr+1; if(bp_type>0) res.parNr++;
  pi=0; strcpy(res.parname[pi], "K1p"); strcpy(res.parunit[pi], "ml/(min*ml)");
  pi++; strcpy(res.parname[pi], "K1p/k2p"); strcpy(res.parunit[pi], "ml/ml");
  pi++; strcpy(res.parname[pi], "K1m/K1p"); strcpy(res.parunit[pi], "");
  pi++; strcpy(res.parname[pi], "K1m/k2m"); strcpy(res.parunit[pi], "ml/ml");
  pi++; strcpy(res.parname[pi], "km"); strcpy(res.parunit[pi], "1/min");
  pi++; strcpy(res.parname[pi], "Vb"); strcpy(res.parunit[pi], "%");
  if(bp_type>0) {
    pi++;
    if(bp_type==1) {
      strcpy(res.parname[pi], "DVR"); strcpy(res.parunit[pi], "ml/ml");
    } else if(bp_type==2) {
      strcpy(res.parname[pi], "BPnd"); strcpy(res.parunit[pi], "");
    } else {strcpy(res.parname[pi], "BPp"); strcpy(res.parunit[pi], "");}
  }
  pi++; strcpy(res.parname[pi], "WSS"); strcpy(res.parunit[pi], "");

  /*
   *  Allocate memory for fitted curves
   */
  ret=dftdup(&dft, &fit);
  if(ret) {
    fprintf(stderr, "Error %d in memory allocation for fitted curves.\n", ret);
    dftEmpty(&input); dftEmpty(&dft); resEmpty(&res);
    return(8);
  }


  /*
   *  Fit ROIs
   */
  if(verbose>0) {
    fprintf(stdout, "fitting regional TACs: ");
    if(verbose>1) fprintf(stdout, "\n");
  }
  fflush(stdout);
  for(ri=0; ri<dft.voiNr; ri++) {
    if(verbose>2) printf("\n  %d %s\n", ri, dft.voi[ri].name);

    /* Initiate values */
    petmeas=dft.voi[ri].y; petsim=fit.voi[ri].y;

    /* Set constraints */
    if(ri!=ref) {
      is_this_ref=0;
      for(pi=0; pi<parNr; pi++) {pmin[pi]=def_pmin[pi]; pmax[pi]=def_pmax[pi];}
    } else {
      is_this_ref=1; if(verbose>2) printf("\n  this is reference region\n");
      for(pi=0; pi<parNr; pi++) {
        pmin[pi]=def_pmin_ref[pi]; pmax[pi]=def_pmax_ref[pi];}
    }
    for(pi=fittedparNr=0; pi<parNr; pi++) if(pmax[pi]>pmin[pi]) fittedparNr++;
    if(ri==0 && verbose>1) {
      printf("  constraints :=");
      for(pi=0; pi<parNr; pi++) printf(" [%g,%g]", pmin[pi], pmax[pi]);
      printf("\n");
      printf("fittedparNr := %d\n", fittedparNr);
    }

    /* Fit */
    TGO_LOCAL_INSIDE = 0;
    TGO_SQUARED_TRANSF = 1;
    tgoNr=60+30*fittedparNr; /* 0, 100 */
    neighNr=6*fittedparNr; /* 6, 10 */
    iterNr=0;
    ret=tgo(
      pmin, pmax, func1TCMdi, NULL, parNr, neighNr,
      &wss, res.voi[ri].parameter, tgoNr, iterNr, verbose-8);
    if(ret>0) {
      fprintf(stderr, "\nError in optimization (%d).\n", ret);
      dftEmpty(&input); dftEmpty(&dft); dftEmpty(&fit); resEmpty(&res); 
      return(9);
    }
    /* Correct fitted parameters to match constraints like inside function */
    (void)modelCheckParameters(parNr, pmin, pmax, res.voi[ri].parameter,
                               res.voi[ri].parameter, &f);
    wss/=f; // remove any penalties from WSS

    /* Bootstrap */
    if(doBootstrap) {
      if(verbose>2) printf("\n  bootstrapping\n");
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
        fit.voi[ri].y,
        // tissue TAC noisy data is written to be used by objf
        petmeas, 
        parNr, dft.w, func1TCMdi, tmp, verbose-4
      );
      if(ret) {
        fprintf(stderr, "Error in bootstrap: %s\n", tmp);
        for(pi=0; pi<parNr; pi++) {
          if(doSD) sd[pi]=nan(""); 
          if(doCL) cl1[pi]=cl2[pi]=nan("");
        }
      }
    }
    
    /* Set fixed parameter values */
    if(is_this_ref && fixed_ref_Vfm_eq_Vfp!=0) {
      res.voi[ri].parameter[3]=res.voi[ri].parameter[1];
    }

    /* Set results wss */
    res.voi[ri].parameter[res.parNr-1]=wss;

    /* done with this region */
    if(dft.voiNr>2 && verbose==1) {
      fprintf(stdout, "."); fflush(stdout);}

  } /* next region */
  if(verbose>0) {fprintf(stdout, "\n"); fflush(stdout);}

  /* Convert Vb fractions to percent's */
  for(ri=0; ri<res.voiNr; ri++){ 
    res.voi[ri].parameter[parVb]*=100.0;
    if(!isnan(res.voi[ri].cl1[parVb])) res.voi[ri].cl1[parVb]*=100.;
    if(!isnan(res.voi[ri].cl2[parVb])) res.voi[ri].cl2[parVb]*=100.;
    if(!isnan(res.voi[ri].sd[parVb])) res.voi[ri].sd[parVb]*=100.;
  }
  /* Calculate DVR, BPnd or BPp */
  if(bp_type==1 || bp_type==2) {
    if(fabs(res.voi[ref].parameter[1])>1.0E-10)
      for(ri=0; ri<res.voiNr; ri++) {
        res.voi[ri].parameter[res.parNr-2]=
          res.voi[ri].parameter[1]/res.voi[ref].parameter[1];
        if(bp_type==2) res.voi[ri].parameter[res.parNr-2]-=1.0;
      }
    else
      for(ri=0; ri<res.voiNr; ri++)
        res.voi[ri].parameter[res.parNr-2]=0.0;
  }
  if(bp_type==3) {
    for(ri=0; ri<res.voiNr; ri++)
      res.voi[ri].parameter[res.parNr-2]=
        res.voi[ri].parameter[1]-res.voi[ref].parameter[1];
  }
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
    dftEmpty(&dft); dftEmpty(&fit); dftEmpty(&input); resEmpty(&res);
    return(11);
  }
  if(verbose>0) fprintf(stdout, "Model parameters written in %s\n", rfile);


  /*
   *  Saving regional metabolite-corrected TACs
   */
  if(mcfile[0]) {
    if(verbose>1) printf("calculating mc curves\n");

    /* Create a duplicate DFT, to not mess up the original or fitted data */
    DFT dft2;
    dftInit(&dft2); ret=dftdup(&dft, &dft2);
    if(ret) {
      fprintf(stderr, "Error: cannot make mc curves.\n");
      dftEmpty(&dft); dftEmpty(&input); resEmpty(&res);
      return(31);
    }
    dft2.frameNr=fitframeNr;

    /* Simulate total tissue TAC, excluding metabolite in tissue */
    double K1p, k2p, K1m, k2m, km, Vb;
    for(ri=0; ri<dft2.voiNr; ri++) {
      /* Set necessary parameters */
      Vb=0.01*res.voi[ri].parameter[parVb];
      K1p=res.voi[ri].parameter[0];
      if(res.voi[ri].parameter[1]>0.0) k2p=K1p/res.voi[ri].parameter[1];
      else k2p=0.0;
      K1m=res.voi[ri].parameter[2]*res.voi[ri].parameter[0];
      if(res.voi[ri].parameter[3]>0.0) k2m=K1m/res.voi[ri].parameter[3];
      else k2m=0.0;
      km=res.voi[ri].parameter[4];
      /* Simulate */
      ret=simC4DIvp(
        input.x, input.voi[0].y, input.voi[1].y, input.voi[2].y, input.frameNr,
        K1p, k2p, 0, 0, 0, 0, 0, km, K1m, k2m,
        0.0, Vb, 1.0, input.voi[0].y2, NULL, NULL, NULL, NULL, NULL, NULL,
	verbose-20); 
      if(ret) {
        if(verbose>1) printf("error %d in simulation\n", ret);
        fprintf(stderr, "Error: cannot calculate metabolite-free curve.\n");
        dftEmpty(&dft2); dftEmpty(&dft); dftEmpty(&input); resEmpty(&res);
        return(32);
      }
      /* Interpolate to measured PET frames */
      if(dft2.timetype==DFT_TIME_STARTEND)
        ret=interpolate4pet(
          input.x, input.voi[0].y2, input.frameNr,
          dft2.x1, dft2.x2, dft2.voi[ri].y, NULL, NULL, dft2.frameNr);
      else
        ret=interpolate(
          input.x, input.voi[0].y2, input.frameNr,
          dft2.x, dft2.voi[ri].y, NULL, NULL, dft2.frameNr);
      if(ret) {
        if(verbose>1) printf("error %d in interpolation\n", ret);
        fprintf(stderr, "Error: cannot interpolate metabolite-free curve.\n");
        dftEmpty(&dft2); dftEmpty(&dft); dftEmpty(&input); resEmpty(&res);
        return(33);
      }
    } // next region

    /* Save metabolite corrected TACs */
    if(verbose>1) printf("saving mc curves\n");
    char tmp[64];
    tpcProgramName(argv[0], 1, 0, tmp, 64);
    sprintf(dft2.comments, "# program := %s\n", tmp);    
    if(dftWrite(&dft2, mcfile)) {
      fprintf(stderr, "Error in writing '%s': %s\n", mcfile, dfterrmsg);
      dftEmpty(&dft2); dftEmpty(&dft); dftEmpty(&input); resEmpty(&res);
      return(34);
    }
    if(verbose>0) printf("MC TACs written in %s\n", mcfile);

    dftEmpty(&dft2);
  }


  /*
   *  Saving and/or plotting of fitted TACs
   */
  if(svgfile[0] || ffile[0]) {

    char tmp[64];
    fit.frameNr=fitframeNr;

    /* Save SVG plot of fitted and original data */
    if(svgfile[0]) {
      if(verbose>1) printf("saving SVG plot\n");
      sprintf(tmp, "1TCM fit with dual input: ");
      if(strlen(dft.studynr)>0) strcat(tmp, dft.studynr);
      ret=plot_fitrange_svg(&dft, &fit, tmp, 0.0, 1.02*dft.x[fitframeNr-1],
                            0.0, nan(""), svgfile, verbose-5);
      if(ret) {
        fprintf(stderr, "Error (%d) in writing '%s'.\n", ret, svgfile);
        dftEmpty(&fit); dftEmpty(&dft); dftEmpty(&input); resEmpty(&res);
        return(30+ret);
      }
      if(verbose>0) printf("Plots written in %s\n", svgfile);
    }

    /* Delete reference region(s) from the data */
    if(refAdded!=0) {
      for(ri=fit.voiNr-1; ri>=0; ri--) if(fit.voi[ri].sw!=0)
        dftDelete(&fit, ri);
    }

    /* Save fitted TACs */
    if(ffile[0]) {
      if(verbose>1) printf("saving fitted curves\n");
      tpcProgramName(argv[0], 1, 0, tmp, 64);
      sprintf(fit.comments, "# program := %s\n", tmp);    
      if(dftWrite(&fit, ffile)) {
        fprintf(stderr, "Error in writing '%s': %s\n", ffile, dfterrmsg);
        dftEmpty(&fit); dftEmpty(&dft); dftEmpty(&input); resEmpty(&res);
        return(22);
      }
      if(verbose>0) printf("Fitted TACs written in %s\n", ffile);
    }

  }

  dftEmpty(&dft); dftEmpty(&fit); dftEmpty(&input); resEmpty(&res);
 
  return(0);
}
/*****************************************************************************/

/*****************************************************************************
 *
 *  Functions to be minimized
 *
 *****************************************************************************/
double func1TCMdi(int parNr, double *p, void *fdata)
{
  int fi, ret;
  double K1p, k2p, K1m, k2m, km, Vb, d, wss=0.0;
  double pa[MAX_PARAMETERS], penalty=1.0;


  /* Check parameters against the constraints */
  ret=modelCheckParameters(parNr, pmin, pmax, p, pa, &penalty);
  if(fdata) {}

  /* Calculate actual model parameters */
  K1p=pa[0]; K1m=pa[2]*K1p; Vb=pa[parVb];
  if(pa[1]>0.0) k2p=K1p/pa[1]; else k2p=0.0;
  if(is_this_ref==0 || fixed_ref_Vfm_eq_Vfp==0) {
    if(pa[3]>0.0) k2m=K1m/(pa[3]); else k2m=0.0;
  } else { // Optionally, Vfm=Vfp in reference region
    if(pa[1]>0.0) k2m=K1m/(pa[1]); else k2m=0.0;
  }
  km=pa[4];

  /* Simulate the tissue PET TAC */
  ret=simC4DIvp(
    input.x, input.voi[0].y, input.voi[1].y, input.voi[2].y, input.frameNr,
    K1p, k2p, 0, 0, 0, 0, 0, km, K1m, k2m,
    0.0, Vb, 1.0, input.voi[0].y2, NULL, NULL, NULL, NULL, NULL, NULL, 0); 
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
  for(fi=0; fi<fitframeNr; fi++) if(dft.w[fi]>0.0) {
    d=petmeas[fi]-petsim[fi];
    wss+=dft.w[fi]*d*d;
  }
  wss*=penalty;
  if(0) printf("%g  %g  %g  %g  %g  => %g\n", K1p, k2p, K1m, k2m,Vb,wss);

  return(wss);
}
/*****************************************************************************/

/*****************************************************************************/
/// @endcond
