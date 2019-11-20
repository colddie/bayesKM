/** @file fit_frtm.c
 *  @brief NLLSQ fitting of the parameters of (full) reference tissue compartmental model 
 *   to PET TTACs.
 *  @remark This TPCCLIB version is based on previous fit_frtm version 2.6.0 / 2013-06-26. 
 *  @todo Change simulation to get pre-calculated input integral, which takes into account
 *   the frame lengths.
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
const int parNr=4;
int fitframeNr=0;
double *t, *cr, *ct, *tis, *w; /* These are pointers, not allocated */
double pmin[MAX_PARAMS], pmax[MAX_PARAMS];
double wss_wo_penalty=0.0;
/* Local functions */
double frtmFunc(int parNr, double *p, void*);
/*****************************************************************************/

/*****************************************************************************/
static char *info[] = {
  "NLLSQ estimation of R1 (=K1/K1'), k2, k3, and BPnd (binding potential)",
  "using the (full) reference tissue compartment model, FRTM/RTCM (1,3).",
  "Assumption is that K1/k2 is the same in all brain regions, but 1TCM with",
  "plasma input does not need to fit the tissue curves satisfactorily",
  "as is assumed in SRTM (2), and Cref(t) is not assumed to be the same as",
  "Cfree(t) as is assumed in the ratio methods.",
  " ",
  "Usage: @P [Options] ttacfile reference endtime resultfile",
  " ",
  "TTAC file can be in DFT or PMOD format. Sample times must be in minutes.",
  "If TTAC file contains weights, those are used in the NLLSQ fitting.",
  "Reference region TAC can be given separate TAC file or as the name or number",
  "of the reference region in TTAC file.",
  " ",
  "Options:",
  " -DVR",
  "     Instead of BPnd, program saves the DVR (=BPnd+1) values.",
  " -lim=<filename>",
  "     Specify the constraints for model parameters;",
  "     This file with default values can be created by giving this option",
  "     as the only command-line argument to this program.",
  " -SD[=<y|N>]",
  "     Standard deviations are calculated and saved in results (y), or",
  "     not calculated (n).",
  " -CL[=<y|N>]",
  "     95% Confidence limits are calculated and saved in results (y), or",
  "     not calculated (n).",
  " -w1",
  "     All weights are set to 1.0 (no weighting); by default, weights in",
  "     TTAC file are used, if available.",
  " -wf",
  "     Weight by sampling interval.",
  " -fit=<Filename>",
  "     Fitted regional TACs are written in file.",
  " -svg=<Filename>",
  "     Fitted and measured TACs are plotted in specified SVG file.",
  " -stdoptions", // List standard options like --help, -v, etc
  " ",
  " ",
  "Values of R1, k2, k3, and BPnd are written in the specified result file.",
  "Fitted curves are written in DFT format, if file name is given.",
  " ",
  "Example 1: file a789.tac contains regions-of-interest and reference region,",
  "with name 'cereb all'. The whole time range is used in the fit.",
  "     @P a789.tac 'cereb all' 999 a789.res",
  " ",
  "Example 2: Reference region TAC is in a separate file, a789ref.tac;",
  "standard deviations and confidence limits are also estimated.",
  "     @P -SD=y -CL=y a789.tac a789ref.tac 999 a789.res",
  " ",
  "References:",
  "1. Cunningham VJ, Hume SP, Price GR, Ahier RG, Cremer JE, Jones AKP.",
  "   Compartmental analysis of diprenorphine binding to opiate receptors",
  "   in the rat in vivo and its comparison with equilibrium data in vitro.",
  "   J Cereb Blood Flow Metab 1991;11:1-9.",
  "2. Lammertsma AA, Hume SP. Simplified reference tissue model for PET",
  "   receptor studies. Neuroimage 1996;4:153-158.",
  "3. Oikonen V, Sederholm K. TPCMOD0002: Model equations for reference tissue",
  "   compartmental models. http://www.turkupetcentre.net/reports/tpcmod0002.pdf",
  " ",
  "See also: bfmsrtm, dftweigh, rescoll, logan, fit_srtm, sim_rtcm",
  " ",
  "Keywords: TAC, modelling, binding potential, RTCM, reference input",
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
  int     ai, help=0, version=0, verbose=1;
  char    rtacfile[FILENAME_MAX], ttacfile[FILENAME_MAX], resfile[FILENAME_MAX],
          fitfile[FILENAME_MAX], svgfile[FILENAME_MAX], limfile[FILENAME_MAX];
  char   *cptr, tmp[256];
  double  fitdur=nan("");
  int     weights=0; // 0=default, 1=no weighting, 2=frequency
  int     doDVR=0;
  int     doBootstrap=0, doSD=0, doCL=0;
  double *sd, *cl1, *cl2;
  double  def_pmin[MAX_PARAMETERS], def_pmax[MAX_PARAMETERS];
  int     ret, n;

#ifdef MINGW
  // Use Unix/Linux default of two-digit exponents in MinGW on Windows
  _set_output_format(_TWO_DIGIT_EXPONENT);
#endif

  /* Set parameter initial values and constraints */
  /* R1  */ def_pmin[0]=0.001;     def_pmax[0]=10.0;
  /* k2  */ def_pmin[1]=0.000001;  def_pmax[1]=1.0;
  /* k3  */ def_pmin[2]=0.0;       def_pmax[2]=1.0;
  /* BP  */ def_pmin[3]=0.0;       def_pmax[3]=60.0; // may be reset later

  /*
   *  Get arguments
   */
  if(argc==1) {tpcPrintUsage(argv[0], info, stderr); return(1);}
  rtacfile[0]=ttacfile[0]=resfile[0]=fitfile[0]=svgfile[0]=limfile[0]=(char)0;
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
    } else if(strcasecmp(cptr, "DVR")==0) {
      doDVR=1; continue;
    } else if(strcasecmp(cptr, "W1")==0) {
      weights=1; continue;
    } else if(strcasecmp(cptr, "WF")==0) {
      weights=2; continue;
    } else if(strncasecmp(cptr, "FIT=", 4)==0) {
      strlcpy(fitfile, cptr+4, FILENAME_MAX); if(strlen(fitfile)>0) continue;
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
  if(ai<argc) strlcpy(ttacfile, argv[ai++], FILENAME_MAX);
  if(ai<argc) strlcpy(rtacfile, argv[ai++], FILENAME_MAX);
  if(ai<argc) {
    if(atof_with_check(argv[ai], &fitdur)!=0 || fitdur<0.0) {
      fprintf(stderr, "Error: invalid fit time: '%s'.\n", argv[ai]);
      return(1);
    }
    if(fitdur==0) fitdur=1.0E+10;
    ai++;
  }
  if(ai<argc) strlcpy(resfile, argv[ai++], FILENAME_MAX);
  if(ai<argc) {
    fprintf(stderr, "Error: invalid argument '%s'.\n", argv[ai]);
    return(1);
  }
  if(doSD || doCL) doBootstrap=1; else doBootstrap=0;

  /* In verbose mode print arguments and options */
  if(verbose>1) {
    printf("ttacfile := %s\n", ttacfile);
    printf("reference := %s\n", rtacfile);
    printf("resfile := %s\n", resfile);
    printf("fitfile := %s\n", fitfile);
    printf("svgfile := %s\n", svgfile);
    printf("limfile := %s\n", limfile);
    printf("required_fittime := %g min\n", fitdur);
    printf("weights := %d\n", weights);
    printf("doDVR := %d\n", doDVR);
    printf("doBootstrap := %d\n", doBootstrap);
    printf("doSD := %d\n", doSD);
    printf("doCL := %d\n", doCL);
  }

  /* If only file name for parameter constraints was given, then write one
     with default contents, and exit */
  if(limfile[0] && !ttacfile[0]) {
    /* Check that initial value file does not exist */
    if(strcasecmp(limfile, "stdout")!=0 && access(limfile, 0) != -1) {
      fprintf(stderr, "Error: parameter constraint file %s exists.\n", limfile);
      return(9);
    }
    if(verbose>1 && strcasecmp(limfile, "stdout")!=0) 
      printf("writing parameter constraints file\n");
    /* Create parameter file */
    IFT ift; iftInit(&ift);
    iftPutDouble(&ift, "R1_lower", def_pmin[0], NULL);
    iftPutDouble(&ift, "R1_upper", def_pmax[0], NULL);
    iftPutDouble(&ift, "k2_lower", def_pmin[1], NULL);
    iftPutDouble(&ift, "k2_upper", def_pmax[1], NULL);
    iftPutDouble(&ift, "k3_lower", def_pmin[2], NULL);
    iftPutDouble(&ift, "k3_upper", def_pmax[2], NULL);
    iftPutDouble(&ift, "BP_lower", def_pmin[3], NULL);
    iftPutDouble(&ift, "BP_upper", def_pmax[3], NULL);
    ret=iftWrite(&ift, limfile);
    if(ret) {
      fprintf(stderr, "Error in writing '%s': %s\n", limfile, ift.status);
      iftEmpty(&ift); return(9);
    }
    if(strcasecmp(limfile, "stdout")!=0)
      fprintf(stdout, "Parameter file %s with initial values written.\n", limfile);
    iftEmpty(&ift); return(0);
  }


  /* Did we get all the information that we need? */
  if(!resfile[0]) {
    fprintf(stderr, "Error: missing command-line argument; use option --help\n");
    return(1);
  }

  /*
   *  Read model parameter upper and lower limits
   *  if file for that was given
   */
  if(limfile[0]) {
    if(verbose>1) printf("reading %s\n", limfile);
    IFT ift; iftInit(&ift);
    double v;
    ret=iftRead(&ift, limfile, 1);
    if(ret) {
      fprintf(stderr, "Error in reading '%s': %s\n", limfile, ift.status);
      return(9);
    }
    if(verbose>10) iftWrite(&ift, "stdout");
    int n=0;
    if(iftGetDoubleValue(&ift, 0, "R1_lower", &v)>=0) {def_pmin[0]=v; n++;}
    if(iftGetDoubleValue(&ift, 0, "R1_upper", &v)>=0) {def_pmax[0]=v; n++;}
    if(iftGetDoubleValue(&ift, 0, "k2_lower", &v)>=0) {def_pmin[1]=v; n++;}
    if(iftGetDoubleValue(&ift, 0, "k2_upper", &v)>=0) {def_pmax[1]=v; n++;}
    if(iftGetDoubleValue(&ift, 0, "k3_lower", &v)>=0) {def_pmin[2]=v; n++;}
    if(iftGetDoubleValue(&ift, 0, "k3_upper", &v)>=0) {def_pmax[2]=v; n++;}
    if(iftGetDoubleValue(&ift, 0, "BP_lower", &v)>=0) {def_pmin[3]=v; n++;}
    if(iftGetDoubleValue(&ift, 0, "BP_upper", &v)>=0) {def_pmax[3]=v; n++;}
    iftEmpty(&ift);
    if(n==0) {fprintf(stderr, "Error: invalid parameter file.\n"); return(9);}
  }
  /* Check that these are ok */
  n=0; ret=0;
  for(int pi=0; pi<parNr; pi++) {
    if(verbose>3) printf(" %d %g %g\n", pi+1, def_pmin[pi], def_pmax[pi]);
    if(def_pmin[pi]<0.0 && pi!=3) ret++; // Lower limit for BP can be negative
    if(def_pmax[pi]<def_pmin[pi]) ret++;
    if(def_pmax[pi]>def_pmin[pi]) n++;
    if(verbose>3 && ret>0) printf("   -> invalid\n");
  }
  if(ret!=0) {
    fprintf(stderr, "Error: invalid parameter constraints.\n");
    return(9);
  }
  if(n==0) {
    fprintf(stderr, "Error: no model parameters left free for fitting.\n");
    return(9);
  }
  if(verbose>1) {
    printf("Parameter constraints:\n");
    for(int pi=0; pi<parNr; pi++) {
      printf("def_pmin[%d] := %g\n", pi+1, def_pmin[pi]);
      printf("def_pmax[%d] := %g\n", pi+1, def_pmax[pi]);
    }
  }



  /*
   *  Read TTAC file
   */
  if(verbose>1) printf("reading %s\n", ttacfile);
  DFT dft; dftInit(&dft);
  if(dftRead(ttacfile, &dft)) {
    fprintf(stderr, "Error in reading '%s': %s\n", ttacfile, dfterrmsg);
    return(2);
  }
  /* Check for NA's */
  if(dft_nr_of_NA(&dft) > 0) {
    fprintf(stderr, "Error: missing sample(s) in %s\n", ttacfile);
    dftEmpty(&dft); return(2);
  }
  /* Sort the data by increasing sample times */
  dftSortByFrame(&dft);
  /* Set time unit to min */
  ret=dftTimeunitConversion(&dft, TUNIT_MIN);
  if(ret) fprintf(stderr, "Warning: check that regional data times are in minutes.\n");
  /* Remove frame overlaps and gaps */
  if(dft.timetype==DFT_TIME_STARTEND) {
    if(verbose>2) fprintf(stdout, "checking frame overlap in %s\n", ttacfile);
    ret=dftDeleteFrameOverlap(&dft);
    if(ret) {
      fprintf(stderr, "Error: %s has overlapping frame times.\n", ttacfile);
      dftEmpty(&dft); return(2);
    }
  }
  /* Set fit duration */
  int first, last;
  double starttime, endtime;
  starttime=0.0; endtime=fitdur;
  fitframeNr=fittime_from_dft(&dft, &starttime, &endtime, &first, &last, verbose-2);
  if(fitframeNr<5) {
    fprintf(stderr, "Error: too few data points for a decent fit.\n");
    dftEmpty(&dft); return(2);
  }
  if(verbose>2) {
    printf("dft.frameNr := %d\n", dft.frameNr);
    printf("starttime := %g\n", starttime);
    printf("endtime := %g\n", endtime);
    printf("first := %d\n", first);
    printf("last := %d\n", last);
    printf("fitframeNr := %d\n", fitframeNr);
  }
  fitdur=endtime;
  /* Check that there is not any significant delay in the beginning of the data */
  if(dft.timetype==DFT_TIME_STARTEND) {
    if(dft.x1[0]>0.45) {
      fprintf(stderr, "Error: TACs must start at time zero.\n");
      dftEmpty(&dft); return(2);
    }
    if(dft.x1[0]>0.0833333) {
      fprintf(stderr, "Warning: TACs should start at time zero.\n");
    }
  }
  if(verbose>2) printf("Tissue calibration unit := %s\n", dft.unit);

  /* Add data weights, if requested */
  if(weights==1) {
    dft.isweight=0; 
    for(int i=0; i<dft.frameNr; i++) dft.w[i]=1.0;
  } else if(weights==2) {
    if(dftWeightByFreq(&dft)!=0) {
      fprintf(stderr, "Error: cannot set data weights.\n");
      dftEmpty(&dft); return(2);
    }
  } else if(dft.isweight==0) {
    fprintf(stderr, "Warning: data is not weighted.\n");
  }
  /* Print the weights */
  if(verbose>2) {
    fprintf(stdout, "common_data_weights := %g", dft.w[0]);
    for(int i=1; i<dft.frameNr; i++) fprintf(stdout, ", %g", dft.w[i]);
    fprintf(stdout, "\n");
  }


  /*
   *  Read reference TAC
   */
  if(verbose>1) printf("\nreading reference\n");
  int inputtype=-1, ref=-1;
  ret=dftReadReference(&dft, rtacfile, &inputtype, &ref, tmp, verbose-1);
  if(ret<=0) {
    fprintf(stderr, "Error in reading reference input: %s\n", tmp);
    dftEmpty(&dft); if(verbose>1) printf("ret := %d\n", ret);
    return(3);
  }
  if(ret>1) {
    fprintf(stderr, "Warning: several reference regions found: %s selected.\n", dft.voi[ref].name);
  } else if(verbose>1) printf("reference_region := %s\n", dft.voi[ref].name);
  if(verbose>2) printf("inputtype := %d\n", inputtype);

  /* Calculate tissue integrals to determine TAC-based constraints for BP value */
  if(!limfile[0]) {
    if(dft.timetype==DFT_TIME_STARTEND)
      for(int ri=0; ri<dft.voiNr; ri++)
        petintegrate(dft.x1, dft.x2, dft.voi[ri].y, fitframeNr, dft.voi[ri].y3, NULL);
    else 
      for(int ri=0; ri<dft.voiNr; ri++)
        integrate(dft.x, dft.voi[ri].y, fitframeNr, dft.voi[ri].y3);
  }

  /* Allocate an extra TAC for the bootstrap */
  int bsi=-1;
  if(doBootstrap) {
    ret=dftAddmem(&dft, 1); if(ret) {
      fprintf(stderr, "Error: cannot allocate more memory.\n");
      dftEmpty(&dft); return(4);
    }
    bsi=dft.voiNr;
    strcpy(dft.voi[bsi].voiname, "BS");
    strcpy(dft.voi[bsi].name, "BS");
  }



  /*
   *  Prepare the room for results
   */
  if(verbose>1) printf("initializing result data\n");
  RES res; resInit(&res);
  ret=res_allocate_with_dft(&res, &dft); if(ret!=0) {
    fprintf(stderr, "Error: cannot set-up memory for results.\n");
    dftEmpty(&dft); return(4);
  }
  /* Copy titles & filenames */
  tpcProgramName(argv[0], 1, 1, res.program, 256);
  strcpy(res.datafile, ttacfile);
  if(inputtype!=5 && rtacfile[0]) strcpy(res.reffile, rtacfile);
  if(ref>=0) strcpy(res.refroi, dft.voi[ref].name);
  strcpy(res.fitmethod, "TGO");
  /* Constants */
  res.isweight=dft.isweight;
  /* Set data range */
  sprintf(res.datarange, "%g - %g %s", 0.0, fitdur, petTunit(dft.timeunit));
  res.datanr=fitframeNr;
  /* Set current time to results */
  res.time=time(NULL);
  /* Set parameter number, including also the extra "parameters"
     and the parameter names and units */
  res.parNr=5;
  {
    int pi;
    pi=0; strcpy(res.parname[pi], "R1"); strcpy(res.parunit[pi], "");
    pi++; strcpy(res.parname[pi], "k2"); strcpy(res.parunit[pi], "1/min");
    pi++; strcpy(res.parname[pi], "k3"); strcpy(res.parunit[pi], "1/min");
    pi++; if(doDVR==0) strcpy(res.parname[pi], "BP"); else strcpy(res.parname[pi], "DVR");
          strcpy(res.parunit[pi], "");
    pi++; strcpy(res.parname[pi], "WSS"); strcpy(res.parunit[pi], "");
  }




  /*
   *  Fit one VOI at a time
   */
  if(verbose>0) printf("\nfitting...\n");
  int tgoNr=0, neighNr=0, iterNr=0;
  double wss;
  /* Set common data pointers */
  t=dft.x; cr=dft.voi[ref].y;
  double refIntegral=dft.voi[ref].y3[fitframeNr-1];
  /* Fit model to one TAC at a time */
  for(int ri=0; ri<dft.voiNr; ri++) if(ri!=ref) {

    if(verbose>1) printf("Region %d %s\n", ri+1, dft.voi[ri].name);
    /* Set data pointers */
    tis=dft.voi[ri].y; ct=dft.voi[ri].y2; w=dft.w;
    double *p=res.voi[ri].parameter;

    /* Set common parameter constraints */
    for(int pi=0; pi<parNr; pi++) {pmin[pi]=def_pmin[pi]; pmax[pi]=def_pmax[pi];}
    /* Set BP limits based on the integrals of roi and ref TACs
       if those were not given in limfile */
    if(!limfile[0] && refIntegral>0.0) {
      double a, b, c;
      a=(dft.voi[ri].y3[fitframeNr-1]/refIntegral);
      if(a<1.0) a=1.0; 
      b=0.0*a; c=5.0*a; 
      pmin[3]=b; pmax[3]=c;
    }
    if(verbose>3) {
      printf("Parameter constraints:\n");
      for(int pi=0; pi<parNr; pi++) printf("  %10.3E - %10.3E\n", pmin[pi], pmax[pi]);
    }

    /* Fit */
    if(verbose>2) printf("  fitting curve...\n");
    TGO_LOCAL_INSIDE=0;
    TGO_SQUARED_TRANSF=0;
    tgoNr=260;
    neighNr=20;
    ret=tgo(pmin, pmax, frtmFunc, NULL, parNr, neighNr, &wss, p, tgoNr, iterNr, verbose-8);
    if(ret>0) {
      fprintf(stderr, "Error in optimization (%d).\n", ret);
      dftEmpty(&dft); resEmpty(&res); return(6);
    }
    if(verbose>3) {
      for(int pi=0; pi<parNr; pi++) printf(" %g", p[pi]);
      printf(" -> WSS=%g\n", p[parNr]); fflush(stdout);
    }
    /* Correct fitted parameters to match constraints like inside the function */
    (void)modelCheckParameters(parNr, pmin, pmax, p, p, NULL);
    p[parNr]=wss=wss_wo_penalty;
    if(verbose>2) printf("wss := %g\nfitframeNr := %d\n", wss, fitframeNr);

    /* Bootstrap */
    if(doBootstrap) {
      if(verbose>2) printf("  bootstrapping...\n");
      /* bootstrap changes measured and simulated data, therefore use copies */
      tis=dft.voi[bsi].y; ct=dft.voi[bsi].y2;
      if(doSD) sd=res.voi[ri].sd; else sd=NULL;
      if(doCL) {cl1=res.voi[ri].cl1; cl2=res.voi[ri].cl2;} else cl1=cl2=NULL;
      ret=bootstrap(
        0, cl1, cl2, sd, p, pmin, pmax, fitframeNr,
        // measured and fitted original TAC, not modified
        dft.voi[ri].y, dft.voi[ri].y2,
        // tissue TAC noisy data is written to be used by objf
        tis, 
        parNr, w, frtmFunc, tmp, verbose-5
      );
      if(ret) {
        fprintf(stderr, "Error in bootstrap: %s\n", tmp);
        for(int pi=0; pi<parNr; pi++) {
          if(doSD) sd[pi]=nan(""); 
          if(doCL) cl1[pi]=cl2[pi]=nan("");
        }
      }
      // back to what pointers were
      tis=dft.voi[ri].y; tis=dft.voi[ri].y2;
    }

  } /* Next VOI */
  if(verbose>0) {fprintf(stdout, "\n"); fflush(stdout);}

  /* Delete reference region(s) from the results unless it already existed in data */
  if(inputtype==5) {
    resDelete(&res, ref);
  } else {
    for(int i=dft.voiNr-1; i>=0; i--) if(dft.voi[i].sw!=0) {
      resDelete(&res, i);
    }
    ref=-1;
  }

  /*
   *  Convert BP to DVR, if necessary
   */
  if(doDVR!=0) {
    if(verbose>1) printf("converting BP to DVR\n");
    for(int ri=0; ri<res.voiNr; ri++) {
      res.voi[ri].parameter[3]+=1.0;
      if(doCL && !isnan(res.voi[ri].cl1[3])) res.voi[ri].cl1[3]+=1.0;
      if(doCL && !isnan(res.voi[ri].cl2[3])) res.voi[ri].cl2[3]+=1.0;
    }
  }

  /*
   *  Print results on screen
   */
  if(verbose>0) {resPrint(&res); fprintf(stdout, "\n");}


  /*
   *  Save results
   */
  if(verbose>1) printf("saving results in %s\n", resfile);
  ret=resWrite(&res, resfile, verbose-5);
  if(ret) {
    fprintf(stderr, "Error in writing '%s': %s\n", resfile, reserrmsg);
    resEmpty(&res); dftEmpty(&dft);
    return(11);
  }
  if(verbose>0) fprintf(stdout, "Model parameters written in %s\n", resfile);



  /*
   *  Saving and/or plotting of fitted TACs
   */
  if(svgfile[0] || fitfile[0]) {

    /* Create a DFT containing fitted TACs */
    char tmp[64];
    DFT dft2;
    dftInit(&dft2); ret=dftdup(&dft, &dft2);
    if(ret) {
      fprintf(stderr, "Error: cannot save fitted curves.\n");
      dftEmpty(&dft); resEmpty(&res);
      return(21);
    }
    for(int ri=0; ri<dft.voiNr; ri++) 
      if(ri!=ref) 
        for(int fi=0; fi<fitframeNr; fi++)
          dft2.voi[ri].y[fi]=dft2.voi[ri].y2[fi];
    dft2.frameNr=fitframeNr;

    /* Save SVG plot of fitted and original data */
    if(svgfile[0]) {
      if(verbose>1) printf("saving SVG plot\n");
      sprintf(tmp, "FRTM fit ");
      if(strlen(dft.studynr)>0) strcat(tmp, dft.studynr);
      ret=plot_fitrange_svg(&dft, &dft2, tmp, 0.0, 1.02*dft.x[fitframeNr-1],
                            0.0, nan(""), svgfile, verbose-8);
      if(ret) {
        fprintf(stderr, "Error (%d) in writing '%s'.\n", ret, svgfile);
        dftEmpty(&dft2); dftEmpty(&dft); resEmpty(&res);
        return(30+ret);
      }
      if(verbose>0) printf("Plots written in %s\n", svgfile);
    }

    /* Delete reference region(s) from the data, unless it already existed in data */
    if(inputtype!=5) {
      for(int i=dft2.voiNr-1; i>=0; i--) if(dft2.voi[i].sw!=0) dftDelete(&dft2, i);
    }

    /* Save fitted TACs */
    if(fitfile[0]) {
      if(verbose>1) printf("saving fitted curves\n");
      tpcProgramName(argv[0], 1, 0, tmp, 128);
      sprintf(dft2.comments, "# program := %s\n", tmp);
      if(dftWrite(&dft2, fitfile)) {
        fprintf(stderr, "Error in writing '%s': %s\n", fitfile, dfterrmsg);
        dftEmpty(&dft2); dftEmpty(&dft); resEmpty(&res);
        return(22);
      }
      if(verbose>0) printf("Fitted TACs written in %s\n", fitfile);
    }

    dftEmpty(&dft2);
  }


  resEmpty(&res);
  dftEmpty(&dft);

  return(0);
}
/*****************************************************************************/

/*****************************************************************************
 *
 *  Functions to be minimized
 *
 *****************************************************************************/
double frtmFunc(int parNr, double *p, void *fdata)
{
  int ret;
  double R1, k2, k3, k4, BP, d, wss=0.0;
  double pa[MAX_PARAMETERS], penalty=1.0;


  /* Check parameters against the constraints */
  ret=modelCheckParameters(parNr, pmin, pmax, p, pa, &penalty);
  if(fdata) {}
  /* Get parameters */
  R1=pa[0]; k2=pa[1]; k3=pa[2]; BP=pa[3];
  if(BP>0.0) k4=k3/BP; else k4=0.0;

  /* Simulate the tissue PET TAC */
  ret=simRTCM(t, cr, fitframeNr, R1, k2, k3, k4, ct, NULL, NULL);
  if(ret) {
    fprintf(stderr, "  error %d in simulation\n", ret);
    return(nan(""));
  }

  /* Calculate error */
  for(int i=0; i<fitframeNr; i++) if(w[i]>0.0) {
    d=ct[i]-tis[i]; 
    wss+=w[i]*d*d;
  }
  wss_wo_penalty=wss;
  wss*=penalty;
  if(0) printf("R1=%g  k2=%g  k3=%g  BP=%g  => %g\n", R1, k2, k3, BP, wss);

  return(wss);
}
/*****************************************************************************/

/*****************************************************************************/
/// @endcond
