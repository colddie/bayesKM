/** @file fitk3.c
 *  @brief Estimates the parameters of irreversible 2-tissue compartment model.
 *  @copyright (c) Turku PET Centre
 *  @author Vesa Oikonen
 */

/*****************************************************************************/
// #include "tpcclibConfig.h"
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
static const int parNr=4;
DFT input, data;
double *petmeas, *petsim;
static double fVb=-1.0;
double pmin[MAX_PARAMETERS], pmax[MAX_PARAMETERS];
double fk1k2;
int fitframeNr;
static double wss_wo_penalty=0.0;

/*****************************************************************************/
/* Local functions */
double cm3Func(int parNr, double *p, void*);
/*****************************************************************************/

// struct COM {
//    int parNr = 4;
//    DFT input, data;
//    double *petmeas, *petsim;
//    double fvb;
//    double pmin[MAX_PARAMETERS], pmax[MAX_PARAMETERS];
//    double fk1k2;
//    int fitframeNr;
//    double wss_wo_penalty = 0.0;
// };


/**
 *  Main
 */
int tcm2_idl(int argc, char **argv)
{

  int          ai, help=0, version=0, verbose=1;
  int          fi, pi, m, n, ret;
  int          lambda_k3=0;
  int          ref=-1, refAdded=0;   // currently reference region is not enabled;
  char        *cptr, refname[FILENAME_MAX], tmp[FILENAME_MAX];
  double       fitdur=1.0E+10;
  double       wss, aic, K1, k2, k3, Vb,Ki;
  RES          res;
  IFT          ift;
  int          doBootstrap=0, doSD=0, doCL=0;
  double      *sd, *cl1, *cl2;
  double      *def_pmin, *def_pmax;
  int          tgoNr=0, neighNr=0, iterNr=0, fittedparNr=0;

  int          dataNr=0, first, last;
  double      *t0, *t1, *tac, *ctt, *output, *weights, *bmatrix; //, *matrix;
  int          voiNr = 1;
  unsigned int    frameNr, isweight = 0, logan_mode = 0, 
                  bootstrapIter, directbp=0,ri =0, inputtype=0;
  // int          fVb = -1.0;

  const char *debugfile = "debug.txt";

  dftInit(&data); dftInit(&input); resInit(&res);

#ifdef MINGW
  // Use Unix/Linux default of two-digit exponents in MinGW on Windows
  _set_output_format(_TWO_DIGIT_EXPONENT);
#endif

  /* 
   make sure interpolation and integration has been done before hand!
  */
  frameNr= *(unsigned int*) argv[0];
  t0     =  (double*) argv[1];
  t1     =  (double*) argv[2];
  tac    =  (double*) argv[3];
  ctt    =  (double*) argv[4];
  output =  (double*) argv[5];
  verbose= *(unsigned int*) argv[6]; 
  isweight =  *(unsigned int*) argv[7];
  weights = (double*) argv[8];
  def_pmin = (double*) argv[9];
  def_pmax = (double*) argv[10];
  fVb      = *(double*) argv[11];
  doSD     = *(unsigned int*) argv[12]; 
  doCL     = *(unsigned int*) argv[13]; 
  bootstrapIter = *(unsigned int*) argv[14]; 
  bmatrix   = (double*) argv[15];
  if(doSD || doCL) doBootstrap=1; else doBootstrap=0;
//   /* Set parameter initial values and constraints */
//   /* K1    */ def_pmin[0]=0.0;       def_pmax[0]=5.0;
//   /* K1/k2 */ def_pmin[1]=0.00001;   def_pmax[1]=10.0;
//   /* k3    */ def_pmin[2]=0.0;       def_pmax[2]=2.0;
//   /* Vb    */ def_pmin[3]=0.0;       def_pmax[3]=0.08;


 /* Check that these are ok */
  n=0; ret=0;
  for(int pi=0; pi<parNr; pi++) {
    if(verbose>3) printf(" %d %g %g\n", pi+1, def_pmin[pi], def_pmax[pi]);
    // Let user set negative lower limits if (s)he so wishes, but upper limits must be positive
    //if(def_pmin[pi]<0.0 && pi!=2) ret++; // Lower limit for BP can be negative
    if(def_pmax[pi]<=0.0) ret++; // Upper limit must be > 0
    if(def_pmax[pi]<def_pmin[pi]) ret++;
    if(def_pmax[pi]>def_pmin[pi]) n++;
    if(verbose>3 && ret>0) printf("   -> invalid\n");
  }
  if(ret!=0) {
    printf( "Error: invalid parameter constraints.\n");
    return(9);
  }
  if(n==0) {
    printf("Error: no model parameters left free for fitting.\n");
    return(9);
  }
  if(verbose>1) {
    printf("Parameter constraints:\n");
    for(int pi=0; pi<parNr; pi++) {
      printf("def_pmin[%d] := %g\n", pi+1, def_pmin[pi]);
      printf("def_pmax[%d] := %g\n", pi+1, def_pmax[pi]);
    }
  }


  /* Fixed/fitted Vb */
  if(fVb>=0.0) def_pmin[3]=def_pmax[3]=fVb;
  if(def_pmin[3]==def_pmax[3]) fVb=def_pmin[3];
  // if(fVb==0.0) strcpy(bfile, "");
  if(verbose>1) {
    // printf("bfile := %s\n", bfile);
    if(fVb>=0.0) printf("fVb := %g\n", fVb);
  }


  ret=dftSetmem(&data, frameNr, voiNr);  
  // ret=dftSetmem(&temp, frameNr, voiNr);
  ret=dftSetmem(&input, frameNr, voiNr);

  /* Set voiNr and frameNr, and type */
  data.voiNr=voiNr; data.frameNr=frameNr;
  data.isweight=isweight;
  data._type=DFT_FORMAT_PLAIN;
  char cnr[1] = "1";
  strcpy(cnr, data.studynr);
  char  cunit[6]="kBq/mL";
  strcpy(cunit, data.unit);
  data.timeunit=2;
  data.timetype=3;
  input.voiNr=voiNr; input.frameNr=frameNr;
  input.isweight=isweight;
  input._type=DFT_FORMAT_PLAIN;
  strcpy(cnr, input.studynr);
  strcpy(cunit, input.unit);
  input.timeunit=2;
  input.timetype=3;

  for (int i=0; i<frameNr; i++) { 
    data.x1[i] = *(t0+i);
    data.x2[i] = t1[i];
    data.x[i]=0.5*(data.x1[i]+data.x2[i]);
    data.voi[ri].y[i]= tac[i];

    input.x1[i] = *(t0+i);
    input.x2[i] = t1[i];
    input.x[i]=0.5*(input.x1[i]+input.x2[i]);
    input.voi[ri].y[i]= ctt[i];

    if(data.isweight) { data.w[i]=weights[i]; input.w[i]=weights[i];    }
    if(!data.isweight) { data.w[i]=1.0;  input.w[i]=1.0; }
}

printf("tissue data...\n");
dftPrint(&data); 
printf("input data...\n");
dftPrint(&input); 

/* Sort the data by increasing sample times */
  dftSortByFrame(&data);
  /* Set time unit to min */
  ret=dftTimeunitConversion(&data, TUNIT_MIN);
  if(ret) printf( "Warning: check that regional data times are in minutes.\n");
  /* Remove frame overlaps and gaps */
  if(data.timetype==DFT_TIME_STARTEND) {
    if(verbose>2) printf( "checking frame overlap in \n");
    ret=dftDeleteFrameOverlap(&data);
    if(ret) {
      printf("Error: %s has overlapping frame times.\n");
      dftEmpty(&data); return(2);
    }
  }

  /* Set fit duration */
  //   int first, last;
  double starttime, endtime;
  starttime=0.0; endtime=fitdur;
  fitframeNr=fittime_from_dft(&data, &starttime, &endtime, &first, &last, verbose-2);
  if(fitframeNr<4) {
    printf( "Error: too few data points for a decent fit. %d\n",fitframeNr);
    dftEmpty(&data); return(2);
  }
  if(verbose>2) {
    printf("dft.frameNr := %d\n", data.frameNr);
    printf("starttime := %g\n", starttime);
    printf("endtime := %g\n", endtime);
    printf("first := %d\n", first);
    printf("last := %d\n", last);
    printf("fitframeNr := %d\n", fitframeNr);
  }
  fitdur=endtime;
  /* Check that there is not any significant delay in the beginning of the data */
  if(data.timetype==DFT_TIME_STARTEND) {
    if(data.x1[0]>0.45) {
      printf("Error: TACs must start at time zero.\n");
      dftEmpty(&data); return(2);
    }
    if(data.x1[0]>0.0833333) {
      printf("Warning: TACs should start at time zero.\n");
    }
  }
  if(verbose>2) printf("Tissue calibration unit := %s\n", data.unit);

  /* Print the weights */
  if(verbose>2) {
    printf( "common_data_weights := %g", data.w[0]);
    for(int i=1; i<data.frameNr; i++) printf( ", %g", data.w[i]);
    printf("\n");
  }





 /* Allocate an extra TAC for the bootstrap */
  if(doBootstrap) {
    ret=dftAddmem(&data, 1); if(ret) {
      printf("Error: cannot allocate more memory.\n");
      dftEmpty(&data); dftEmpty(&input); return(9);
    }
    strcpy(data.voi[data.voiNr].voiname, "BS");
    strcpy(data.voi[data.voiNr].name, "BS");
  }
  if(verbose>10) dftPrint(&data);  



/*
   *  Prepare the room for the results
   */
  if(verbose>1) printf("initializing result data\n");
  ret=res_allocate_with_dft(&res, &data); if(ret!=0) {
    printf( "Error: cannot setup memory for results.\n");
    dftEmpty(&input); dftEmpty(&data); return(7);
  }
  /* Copy titles & filenames */
  tpcProgramName(argv[0], 1, 1, res.program, 256);
  // strcpy(res.datafile, dfile);
  // strcpy(res.plasmafile, pfile);
  // strcpy(res.bloodfile, bfile);
  if(ref>=0) sprintf(res.refroi, "%s", data.voi[ref].name);
  // if(refname[0]) strcpy(res.reffile, refname);
  strcpy(res.fitmethod, "TGO");
  /* Constants */
  res.isweight=data.isweight;
  if(fVb>=0.0) res.Vb=100.0*fVb;
  /* Set data range */
  sprintf(res.datarange, "%g - %g %s", 0.0, fitdur, dftTimeunit(data.timeunit));
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
   *  Fit other than reference regions
   */
  if(verbose>0) {printf( "fitting regional TACs: ");}
  if(verbose>1) printf("\n");
  for(ri=0; ri<data.voiNr; ri++) if(data.voi[ri].sw==0) {
    if(verbose>2) printf("\n  %d %s:\n", ri, data.voi[ri].name);

    /* Initiate values */
    petmeas=data.voi[ri].y; petsim=data.voi[ri].y2;

    /* Set constraints */
    pmin[0]=def_pmin[0];    pmax[0]=def_pmax[0];   /* K1    */
    pmin[1]=def_pmin[1];    pmax[1]=def_pmax[1];   /* K1/k2 */
    // if(ref>=0) {
    //   pmin[1]=pmax[1]=fk1k2;                       /* K1/k2 */
    // }
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
      printf( "\nError in optimization (%d).\n", ret);
      dftEmpty(&input); dftEmpty(&data); resEmpty(&res); return(8);
    }
    /* Correct fitted parameters to match constraints like inside the function */
    (void)modelCheckParameters(parNr, pmin, pmax, res.voi[ri].parameter,
                               res.voi[ri].parameter, NULL);
    wss=wss_wo_penalty;


  // COM com;
  // com.wss_wo_penalty = wss_wo_penalty;
  // com.fitframeNr =fitframeNr;
  // com.fk1k2 = fk1k2;
  // com.pmin = pmin;
  // com.pmax = pmax;
  // com.fvb = fvb;
  // com.petmeas = petmeas;
  // com.petsim = petsim;
  // com.input = input;
  // com.data = data;
  // com. parNr = parNr;

    /* Bootstrap */
    if(doBootstrap) {
      if(verbose>2) printf("  bootstrapping\n");
      /* bootstrap changes measured and simulated data, therefore use copies */
      petmeas=data.voi[data.voiNr].y2; petsim=data.voi[data.voiNr].y3;
      if(doSD) sd=res.voi[ri].sd; else sd=NULL;
      if(doCL) {cl1=res.voi[ri].cl1; cl2=res.voi[ri].cl2;} else cl1=cl2=NULL;
      // ret=bootstrap(
      //   bootstrapIter, cl1, cl2, sd,
      //   res.voi[ri].parameter, pmin, pmax, fitframeNr,
      //   // measured original TAC, not modified
      //   data.voi[ri].y,
      //   // fitted TAC, not modified
      //   data.voi[ri].y2,
      //   // tissue TAC noisy data is written to be used by objf
      //   petmeas, 
      //   parNr, data.w, cm3Func, tmp, verbose-4
      // );

      ret=bootstrapr(
        bootstrapIter, cl1, cl2, sd,
        res.voi[ri].parameter, pmin, pmax, fitframeNr,
        // measured original TAC, not modified
        data.voi[ri].y,
        // fitted TAC, not modified
        data.voi[ri].y2,
        // tissue TAC noisy data is written to be used by objf
        petmeas, 
        parNr, data.w, cm3Func, tmp, verbose-4,bmatrix
      );

      if(ret) {
        printf( "Error in bootstrap: %s\n", tmp);
        for(pi=0; pi<parNr; pi++) {
          if(doSD) sd[pi]=nan(""); 
          if(doCL) cl1[pi]=cl2[pi]=nan("");
        }
      }
    }

    /* Calculate AIC, based on nr of parameters that actually are fitted */
    for(pi=n=0; pi<parNr; pi++) if(pmax[pi]>pmin[pi]) n++;
    if(verbose>2) printf("nr_of_fitted_parameters := %d\n", n);
    for(fi=m=0; fi<fitframeNr; fi++) if(data.w[fi]>0.0) m++;
    if(verbose>2) printf("nr_of_fitted_samples := %d\n", m);

    aic=aicSS(wss, m, n);

    /* Set results wss and aic */
    res.voi[ri].parameter[res.parNr-2]=wss;
    res.voi[ri].parameter[res.parNr-1]=aic;

    /* Convert Vb fraction to percentage */
    // res.voi[ri].parameter[3]*=100.;
    // if(doSD) res.voi[ri].sd[3]*=100.;
    // if(doCL) {res.voi[ri].cl1[3]*=100.; res.voi[ri].cl2[3]*=100.;}  
    /* Calculate Ki, lambda*k3, and k3/(k2+k3) */
    K1=res.voi[ri].parameter[0];
    k2=K1/res.voi[ri].parameter[1];
    k3=res.voi[ri].parameter[2];
    Vb=res.voi[ri].parameter[3];
    /* Ki */ 
    Ki=K1*k3/(k2+k3);    
    /* lambda*k3 */ 
    res.voi[ri].parameter[5]=res.voi[ri].parameter[1]*res.voi[ri].parameter[2];
    /* k3/(k2+k3) */ 
    res.voi[ri].parameter[6]=k3/(k2+k3);

    /* done with this region */
    // if(data.voiNr>2 && verbose==1) {fprintf(stdout, ".");}

  } /* next region */


if(verbose>0) {printf( "\n");   resPrint(&res);}
  output[0] = K1;     // see below
  output[1] = k2;
  output[2] = k3;
  output[3] = Vb;
  output[4] = Ki;
  output[5] = wss;
  output[6] = aic;
  printf("test me!\n");
  if(doSD) {output[7] = res.voi[0].sd[0];output[8] = res.voi[0].sd[1];output[9] = res.voi[0].sd[2];output[10] = res.voi[0].sd[3]; }
  printf("test me!\n");
    //  if(doSD) { res->voi[i].cl1[j]; res->voi[i].cl2[j]; }
  // /* Delete reference region(s) from the results unless it already existed in data */
  // if(inputtype==5) {
  //   resDelete(&res, ref);
  // } else {
  //   for(int i=data.voiNr-1; i>=0; i--) if(data.voi[i].sw!=0) {
  //     resDelete(&res, i);
  //   }
  //   ref=-1;
  // }

  resEmpty(&res);
  dftEmpty(&data);
  dftEmpty(&input);
  // dftEmpty(&temp);
  return(0);
}

//   pi=0; strcpy(res.parname[pi], "K1"); strcpy(res.parunit[pi], "ml/(min*ml)");
//   pi++; strcpy(res.parname[pi], "K1/k2"); strcpy(res.parunit[pi], "");
//   pi++; strcpy(res.parname[pi], "k3"); strcpy(res.parunit[pi], "1/min");
//   pi++; strcpy(res.parname[pi], "Vb"); strcpy(res.parunit[pi], "%");
//   pi++; strcpy(res.parname[pi], "Ki"); strcpy(res.parunit[pi], "ml/(min*ml)");
//   pi++; strcpy(res.parname[pi], "k3*K1/k2"); strcpy(res.parunit[pi], "1/min");
//   pi++; strcpy(res.parname[pi], "k3/(k2+k3)"); strcpy(res.parunit[pi], "");
//   pi++; strcpy(res.parname[pi], "WSS"); strcpy(res.parunit[pi], "");
//   pi++; strcpy(res.parname[pi], "AIC"); strcpy(res.parunit[pi], "");



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
  // ret=simC3vs(
  //   input.x, input.voi[0].y, input.voi[1].y, input.frameNr,
  //   pa[0], k2, k3, 0.0, 0.0, 0.0, 0.0, Vb, 1.0,
  //   input.voi[0].y2, NULL, NULL, NULL, NULL, NULL);
  ret=simC2( input.x,input.voi[0].y,input.frameNr,pa[0],k2,k3,0.0,
      petsim,NULL,NULL);
  if(ret) {
    printf("error %d in simulation\n", ret);
    return(nan(""));
  }

  // /* Interpolate & integrate to measured PET frames */
  // if(data.timetype==DFT_TIME_STARTEND)
  //     ret=petintegral(data.x1, data.x2, input.voi[0].y2, data.frameNr,
  //                     petsim, NULL);
  //   else
  //     ret=integrate(data.x, input.voi[0].y2, data.frameNr, petsim);
  //   ret=interpolate4pet(
  //     input.x, input.voi[0].y2, input.frameNr,
  //     data.x1, data.x2, petsim, NULL, NULL, fitframeNr);
  // else
  //   ret=interpolate(
  //     input.x, input.voi[0].y2, input.frameNr,
  //     data.x, petsim, NULL, NULL, fitframeNr);
  if(ret) {
    printf("error %d in interpolation\n", ret);
    return(nan(""));
  }

  /* Calculate error */
  for(fi=0, wss=0.0; fi<fitframeNr; fi++) if(data.w[fi]>0.0) {
    d=petmeas[fi]-petsim[fi]; 
    wss+=data.w[fi]*d*d;
  }
  wss_wo_penalty=wss;
  wss*=penalty;
  if(0) printf("K1=%g  k2=%g  k3=%g  Vb=%g  => %g\n",
    pa[0], k2, pa[2], Vb, wss);

  return(wss);
}
/*****************************************************************************/