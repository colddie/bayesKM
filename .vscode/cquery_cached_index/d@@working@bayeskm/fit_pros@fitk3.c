/** @file fitk3.c
 *  @brief Estimates the parameters of irreversible 2-tissue compartment model.
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
DFT input, dft;
double *petmeas, *petsim;
double fVb=-1.0;
double pmin[MAX_PARAMETERS], pmax[MAX_PARAMETERS];
double fk1k2;
int fitframeNr;
double wss_wo_penalty=0.0;
/*****************************************************************************/
/* Local functions */
double cm3Func(int parNr, double *p, void*);
/*****************************************************************************/

/*****************************************************************************/
static char *info[] = {
  "Non-linear fitting of two-tissue compartment model to plasma input, blood,",
  "and tissue time-activity curves (PTAC, BTAC, and TTAC) to estimate",
  "parameters K1, k2, k3, and Vb. Sample times must be in minutes.",
  " ",
  "    ______       ___________________           ",
  "   |      |  K1 |        k3         |          ",
  "   |      | --> |      ------->     |          ",
  "   |  Ca  | <-- |   C1 <------- C2  |          ",
  "   |      |  k2 |        k4=0       |          ",
  "   |______|     |___________________|          ",
  " ",
  "Usage: @P [Options] ptacfile btacfile ttacfile endtime resultfile",
  " ",
  "Options:",
  " -lim[=<filename>]",
  "     Specify the constraints for model parameters;",
  "     This file with default values can be created by giving this",
  "     option as the only command-line argument to this program.",
  "     Without filename the default values are printed on screen.",
  " -SD[=<y|N>]",
  "     Standard deviations are calculated and saved in results (y),",
  "     or not calculated (N, default).",
  "     Program runs a lot faster if SD and CL are not calculated.",
  " -CL[=<y|N>]",
  "     95% Confidence limits are calculated and saved in results (y), or",
  "     not calculated (N, default).",
  " -Vb=<Vb(%)>",
  "     Enter a fixed Vb; fitted by default.",
  " -r=<Reference region id or filename>",
  "     Optional reference region is used to constrain K1/k2 in other regions;",
  "     Also k3 is fitted to reference region data, thus any large region",
  "     (for example cortex) can be used here.",
  " -fit=<Filename>",
  "     Fitted regional TACs are written in DFT format.",
  " -svg=<Filename>",
  "     Fitted and measured TACs are plotted in specified SVG file.",
  " -stdoptions", // List standard options like --help, -v, etc
  " ",
  "Example 1: estimate K1, K1/k2, k3 and Vb",
  "     @P ua919ap.bld ua919ab.bld ua919.tac 60 ua919k3.res",
  " ",
  "Example 2: estimate K1 and k3; Vb is constrained to 4.5% and K1/k2 is",
  "constrained to K1/k2 estimated from region 'occip'",
  "     @P -Vb=4.5 -r=occip ua919ap.bld ua919ab.bld ua919.tac 60 ua919.res",
  " ",
  "See also: patlak, lhsol, p2t_v3c, dftweigh, dftcbv, rescoll",
  " ",
  "Keywords: TAC, modelling, irreversible uptake, k3, Ki, 2TCM",
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
  int          ai, help=0, version=0, verbose=1;
  int          ri, fi, pi, m, n, ret;
  int          lambda_k3=0;
  int          ref=-1, refAdded=0, inputtype;
  char         dfile[FILENAME_MAX], pfile[FILENAME_MAX], bfile[FILENAME_MAX],
               rfile[FILENAME_MAX], ffile[FILENAME_MAX], limfile[FILENAME_MAX];
  char         svgfile[FILENAME_MAX];
  char        *cptr, refname[FILENAME_MAX], tmp[FILENAME_MAX];
  double       fitdur, wss, aic, K1, k2, k3;
  RES          res;
  IFT          ift;
  int          doBootstrap=0, doSD=0, doCL=0;
  double      *sd, *cl1, *cl2;
  double       def_pmin[MAX_PARAMETERS], def_pmax[MAX_PARAMETERS];
  int          tgoNr=0, neighNr=0, iterNr=0, fittedparNr=0;

#ifdef MINGW
  // Use Unix/Linux default of two-digit exponents in MinGW on Windows
  _set_output_format(_TWO_DIGIT_EXPONENT);
#endif

  /* Set parameter initial values and constraints */
  /* K1    */ def_pmin[0]=0.0;       def_pmax[0]=5.0;
  /* K1/k2 */ def_pmin[1]=0.00001;   def_pmax[1]=10.0;
  /* k3    */ def_pmin[2]=0.0;       def_pmax[2]=2.0;
  /* Vb    */ def_pmin[3]=0.0;       def_pmax[3]=0.08;


  /*
   *  Get arguments
   */
  if(argc==1) {tpcPrintUsage(argv[0], info, stderr); return(1);}
  dfile[0]=pfile[0]=bfile[0]=rfile[0]=ffile[0]=refname[0]=limfile[0]=(char)0;
  svgfile[0]=(char)0;
  fitdur=fVb=-1.0; 
  iftInit(&ift); resInit(&res); dftInit(&dft); dftInit(&input);
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
    } else if(strcasecmp(cptr, "LK3")==0) {
      lambda_k3=1; continue;
    } else if(strncasecmp(cptr, "R=", 2)==0 && strlen(cptr)>2) {
      strlcpy(refname, cptr+2, FILENAME_MAX); continue;
    } else if(strncasecmp(cptr, "Vb=", 3)==0 && strlen(cptr)>3) {
      fVb=0.01*atof_dpi(cptr+3);
      if(fVb>=0.0 && fVb<1.0) {
        if(fVb<0.01) fprintf(stderr, "Warning: Vb was set to %g%%\n", 100.*fVb);
        def_pmin[3]=def_pmax[3]=fVb;
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
    printf("dfile := %s\n", dfile);
    printf("rfile := %s\n", rfile);
    printf("ffile := %s\n", ffile);
    printf("svgfile := %s\n", svgfile);
    printf("limfile := %s\n", limfile);
    printf("refname := %s\n", refname);
    printf("lambda_k3 := %d\n", lambda_k3);
    printf("fitdur := %g\n", fitdur);
    printf("doBootstrap := %d\n", doBootstrap);
    printf("doSD := %d\n", doSD);
    printf("doCL := %d\n", doCL);
  }

  /* If only filename for parameter constraints was given, then write one
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
    iftPutDouble(&ift, "K1_lower", def_pmin[0], NULL);
    iftPutDouble(&ift, "K1_upper", def_pmax[0], NULL);
    iftPutDouble(&ift, "K1k2_lower", def_pmin[1], NULL);
    iftPutDouble(&ift, "K1k2_upper", def_pmax[1], NULL);
    iftPutDouble(&ift, "k3_lower", def_pmin[2], NULL);
    iftPutDouble(&ift, "k3_upper", def_pmax[2], NULL);
    iftPutDouble(&ift, "Vb_lower", def_pmin[3], NULL);
    iftPutDouble(&ift, "Vb_upper", def_pmax[3], NULL);
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

  /* Did we get all the information from user that we need? */
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
    double v;
    if(verbose>1) printf("reading %s\n", limfile);
    ret=iftRead(&ift, limfile, 1);
    if(ret) {
      fprintf(stderr, "Error in reading '%s': %s\n", limfile, ift.status);
      return(9);
    }
    if(verbose>10) iftWrite(&ift, "stdout");
    int n=0;
    /* K1 */
    if(iftGetDoubleValue(&ift, 0, "K1_lower", &v)>=0) {def_pmin[0]=v; n++;}
    if(iftGetDoubleValue(&ift, 0, "K1_upper", &v)>=0) {def_pmax[0]=v; n++;}
    /* K1/k2 */
    if(iftGetDoubleValue(&ift, 0, "K1k2_lower", &v)>=0) {def_pmin[1]=v; n++;}
    if(iftGetDoubleValue(&ift, 0, "K1k2_upper", &v)>=0) {def_pmax[1]=v; n++;}
    /* k3 */
    if(iftGetDoubleValue(&ift, 0, "k3_lower", &v)>=0) {def_pmin[2]=v; n++;}
    if(iftGetDoubleValue(&ift, 0, "k3_upper", &v)>=0) {def_pmax[2]=v; n++;}
    /* Vb */
    if(iftGetDoubleValue(&ift, 0, "Vb_lower", &v)>=0) {def_pmin[3]=v; n++;}
    if(iftGetDoubleValue(&ift, 0, "Vb_upper", &v)>=0) {def_pmax[3]=v; n++;}
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
  if(fVb>=0.0) def_pmin[3]=def_pmax[3]=fVb;
  if(def_pmin[3]==def_pmax[3]) fVb=def_pmin[3];
  if(fVb==0.0) strcpy(bfile, "");
  if(verbose>1) {
    printf("bfile := %s\n", bfile);
    if(fVb>=0.0) printf("fVb := %g\n", fVb);
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
  if(fitframeNr<4 || input.frameNr<4) {
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
  /* Copy titles & filenames */
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
  pi++; strcpy(res.parname[pi], "Vb"); strcpy(res.parunit[pi], "%");
  pi++; strcpy(res.parname[pi], "Ki"); strcpy(res.parunit[pi], "ml/(min*ml)");
  pi++; strcpy(res.parname[pi], "k3*K1/k2"); strcpy(res.parunit[pi], "1/min");
  pi++; strcpy(res.parname[pi], "k3/(k2+k3)"); strcpy(res.parunit[pi], "");
  pi++; strcpy(res.parname[pi], "WSS"); strcpy(res.parunit[pi], "");
  pi++; strcpy(res.parname[pi], "AIC"); strcpy(res.parunit[pi], "");


  /*
   *  Fit at first the reference ROI(s), if required
   */
  if(ref>=0) for(ri=0; ri<dft.voiNr; ri++) if(dft.voi[ri].sw>0) {
    if(verbose>0)
      fprintf(stdout, "fitting %s as reference region\n", dft.voi[ri].name);
    /* Initiate values */
    petmeas=dft.voi[ri].y; petsim=dft.voi[ri].y2;
    /* Set constraints */
    pmin[0]=def_pmin[0];    pmax[0]=def_pmax[0];   /* K1    */
    pmin[1]=def_pmin[1];    pmax[1]=def_pmax[1];   /* K1/k2 */
    pmin[2]=def_pmin[2];    pmax[2]=def_pmax[2];   /* k3    */
    pmin[3]=def_pmin[3];    pmax[3]=def_pmax[3];   /* Vb    */
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
    // neighNr=5; tgoNr=300; iterNr=0;
    tgoNr=50+25*fittedparNr;
    neighNr=6*fittedparNr;
    ret=tgo(
      pmin, pmax, cm3Func, NULL, parNr, neighNr,
      &wss, res.voi[ri].parameter, tgoNr, iterNr, verbose-8);
    if(ret>0) {
      fprintf(stderr, "Error in optimization (%d).\n", ret);
      dftEmpty(&input); dftEmpty(&dft); resEmpty(&res); return(8);
    }
    /* Correct fitted parameters to match constraints like inside the function */
    (void)modelCheckParameters(parNr, pmin, pmax, res.voi[ri].parameter,
                               res.voi[ri].parameter, NULL);
    wss=wss_wo_penalty;
    /* Fix K1/k2 to non-reference regions, but only after all reference
       regions are fitted */
    if(ri==ref) {
      fk1k2=res.voi[ri].parameter[1];
      if(verbose>2) {fprintf(stdout, "  K1/k2 := %g\n", fk1k2);}
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
        parNr, dft.w, cm3Func, tmp, verbose-4
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
    aic=aicSS(wss, m, n);

    /* Set results wss and aic */
    res.voi[ri].parameter[res.parNr-2]=wss;
    res.voi[ri].parameter[res.parNr-1]=aic;

    /* Convert Vb fraction to percentage */
    res.voi[ri].parameter[3]*=100.;
    if(doSD) res.voi[ri].sd[3]*=100.;
    if(doCL) {res.voi[ri].cl1[3]*=100.; res.voi[ri].cl2[3]*=100.;}
    /* Calculate Ki, lambda*k3, and k3/(k2+k3) */
    K1=res.voi[ri].parameter[0];
    k2=K1/res.voi[ri].parameter[1];
    k3=res.voi[ri].parameter[2];
    /* Ki */ 
    res.voi[ri].parameter[4]=K1*k3/(k2+k3);
    /* lambda*k3 */ 
    res.voi[ri].parameter[5]=res.voi[ri].parameter[1]*res.voi[ri].parameter[2];
    /* k3/(k2+k3) */ 
    res.voi[ri].parameter[6]=k3/(k2+k3);

  } // next reference region


  /*
   *  Fit other than reference regions
   */
  if(verbose>0) {fprintf(stdout, "fitting regional TACs: "); fflush(stdout);}
  if(verbose>1) printf("\n");
  for(ri=0; ri<dft.voiNr; ri++) if(dft.voi[ri].sw==0) {
    if(verbose>2) printf("\n  %d %s:\n", ri, dft.voi[ri].name);

    /* Initiate values */
    petmeas=dft.voi[ri].y; petsim=dft.voi[ri].y2;

    /* Set constraints */
    pmin[0]=def_pmin[0];    pmax[0]=def_pmax[0];   /* K1    */
    pmin[1]=def_pmin[1];    pmax[1]=def_pmax[1];   /* K1/k2 */
    if(ref>=0) {
      pmin[1]=pmax[1]=fk1k2;                       /* K1/k2 */
    }
    pmin[2]=def_pmin[2];    pmax[2]=def_pmax[2];   /* k3    */
    pmin[3]=def_pmin[3];    pmax[3]=def_pmax[3];   /* Vb    */
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
    // tgoNr=300; iterNr=0; neighNr=5;
    tgoNr=50+25*fittedparNr;
    neighNr=6*fittedparNr;
    iterNr=0;
    ret=tgo(
      pmin, pmax, cm3Func, NULL, parNr, 5,
      &wss, res.voi[ri].parameter, 300, 0, verbose-8);
    if(ret>0) {
      fprintf(stderr, "\nError in optimization (%d).\n", ret);
      dftEmpty(&input); dftEmpty(&dft); resEmpty(&res); return(8);
    }
    /* Correct fitted parameters to match constraints like inside the function */
    (void)modelCheckParameters(parNr, pmin, pmax, res.voi[ri].parameter,
                               res.voi[ri].parameter, NULL);
    wss=wss_wo_penalty;

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
        parNr, dft.w, cm3Func, tmp, verbose-4
      );
      if(ret) {
        fprintf(stderr, "Error in bootstrap: %s\n", tmp);
        for(pi=0; pi<parNr; pi++) {
          if(doSD) sd[pi]=nan(""); 
          if(doCL) cl1[pi]=cl2[pi]=nan("");
        }
      }
    }

    /* Calculate AIC, based on nr of parameters that actually are fitted */
    for(pi=n=0; pi<parNr; pi++) if(pmax[pi]>pmin[pi]) n++;
    if(verbose>2) printf("nr_of_fitted_parameters := %d\n", n);
    for(fi=m=0; fi<fitframeNr; fi++) if(dft.w[fi]>0.0) m++;
    if(verbose>2) printf("nr_of_fitted_samples := %d\n", m);
    aic=aicSS(wss, m, n);

    /* Set results wss and aic */
    res.voi[ri].parameter[res.parNr-2]=wss;
    res.voi[ri].parameter[res.parNr-1]=aic;

    /* Convert Vb fraction to percentage */
    res.voi[ri].parameter[3]*=100.;
    if(doSD) res.voi[ri].sd[3]*=100.;
    if(doCL) {res.voi[ri].cl1[3]*=100.; res.voi[ri].cl2[3]*=100.;}
    /* Calculate Ki, lambda*k3, and k3/(k2+k3) */
    K1=res.voi[ri].parameter[0];
    k2=K1/res.voi[ri].parameter[1];
    k3=res.voi[ri].parameter[2];
    /* Ki */ 
    res.voi[ri].parameter[4]=K1*k3/(k2+k3);
    /* lambda*k3 */ 
    res.voi[ri].parameter[5]=res.voi[ri].parameter[1]*res.voi[ri].parameter[2];
    /* k3/(k2+k3) */ 
    res.voi[ri].parameter[6]=k3/(k2+k3);

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
      sprintf(tmp, "K1-k3 fit: ");
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
double cm3Func(int parNr, double *p, void *fdata)
{
  int fi, ret;
  double Vb, k2, k3, d, wss=0.0;
  double pa[MAX_PARAMETERS], penalty=1.0;


  /* Check parameters against the constraints */
  ret=modelCheckParameters(parNr, pmin, pmax, p, pa, &penalty);
  if(fdata) {}
  /* Calculate k2 and k3 */
  k2=pa[0]/pa[1]; k3=pa[2]; if(fVb>=0.0) Vb=fVb; else Vb=pa[3];

  /* Simulate the tissue PET TAC */
  ret=simC3vs(
    input.x, input.voi[0].y, input.voi[1].y, input.frameNr,
    pa[0], k2, k3, 0.0, 0.0, 0.0, 0.0, Vb, 1.0,
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
  if(0) printf("K1=%g  k2=%g  k3=%g  Vb=%g  => %g\n",
    pa[0], k2, pa[2], Vb, wss);

  return(wss);
}
/*****************************************************************************/

/*****************************************************************************/
/// @endcond
