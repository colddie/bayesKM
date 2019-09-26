/*****************************************************************************/
#include "srtm_idl.h"
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
// const int parNr=3;
// int fitframeNr=0;
// double *t, *cr, *ct, *tis, *w; /* These are pointers, not allocated */
// double pmin[MAX_PARAMS], pmax[MAX_PARAMS];
// double wss_wo_penalty=0.0;
// /* Local functions */
// double srtmFunc(int parNr, double *p, void*);
/*****************************************************************************/


/**
 *  Main
 */
extern "C" int SRTM_IDL::srtm_idl(int argc, char **argv)
{
  int     ai, help=0, version=0, verbose=1;
  char   *cptr, tmp[256];
  double  fitdur=1.0E+10;
  int     doDVR=0;
  int     doBootstrap=0, doSD=0, doCL=0;
  double *sd, *cl1, *cl2;
  double  *def_pmin, *def_pmax;
  int     ret, n, bootstrapIter;

  int        dataNr=0, first, last;
  double    *t0, *t1, *tac, *ctt, *output, *weights, *bmatrix; //, *matrix;
  int       voiNr = 2;
  unsigned int    frameNr, isweight = 0, logan_mode = 0, directbp=0,ri =0,ref=1, inputtype=0;

  const char *debugfile = "debug.txt";

  DFT data, input, temp; // res;
  dftInit(&data); dftInit(&temp); dftInit(&input); //resInit(&res);

    fitframeNr=0;
    wss_wo_penalty=0.0;

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
  doSD     = *(unsigned int*) argv[11]; 
  doCL     = *(unsigned int*) argv[12]; 
  bootstrapIter = *(unsigned int*) argv[13]; 
  bmatrix   = (double*) argv[14];
  if(doSD || doCL) doBootstrap=1; else doBootstrap=0;
//   /* Set parameter initial values and constraints */
//   /* R1  */ def_pmin[0]=0.001;     def_pmax[0]=10.0;
//   /* k2  */ def_pmin[1]=0.000001;  def_pmax[1]=10.0;
//   /* BP  */ def_pmin[2]=0.0;       def_pmax[2]=60.0; // may be reset later


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

  ret=dftSetmem(&data, frameNr, voiNr);  
  ret=dftSetmem(&temp, frameNr, voiNr);
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

  for (int i=0; i<frameNr; i++) { 
    data.x1[i] = *(t0+i);
    data.x2[i] = t1[i];
    data.x[i]=0.5*(data.x1[i]+data.x2[i]);
    data.voi[ri].y[i]= tac[i];

    data.voi[ref].y[i]= ctt[i];
    if(data.isweight) data.w[i]=weights[i];
    if(!data.isweight) data.w[i]=1.0;
}


dftPrint(&data);

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
//   if(fitframeNr<4) {
//     printf( "Error: too few data points for a decent fit.\n");
//     dftEmpty(&data); return(2);
//   }
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

//   /* Add data weights, if requested */
//   if(weights==1) {
//     dft.isweight=0; 
//     for(int i=0; i<dft.frameNr; i++) dft.w[i]=1.0;
//   } else if(weights==2) {
//     if(dftWeightByFreq(&dft)!=0) {
//       fprintf(stderr, "Error: cannot set data weights.\n");
//       dftEmpty(&dft); return(2);
//     }
//   } else if(dft.isweight==0) {
//     fprintf(stderr, "Warning: data is not weighted.\n");
//   }
  /* Print the weights */
  if(verbose>2) {
    printf( "common_data_weights := %g", data.w[0]);
    for(int i=1; i<data.frameNr; i++) printf( ", %g", data.w[i]);
    printf("\n");
  }


//   /* Check that original input data started early enough, otherwise AUC(0-T)
//      might be wrong */ 
//   if(istart>0.3) {
//     printf("Warning: input TAC should start at time zero.\n");
//   }


 /* Integrate tissue data */
  if(verbose>1) printf("integrating tissue data\n");
  for(ri=0; ri<data.voiNr; ri++) {             // include both tissue and input
    if(data.timetype==DFT_TIME_STARTEND)
      for(int ri=0; ri<data.voiNr; ri++)
        petintegrate(data.x1, data.x2, data.voi[ri].y, fitframeNr, data.voi[ri].y3, NULL);
    else 
      for(int ri=0; ri<data.voiNr; ri++)
        integrate(data.x, data.voi[ri].y, fitframeNr, data.voi[ri].y3);
    if(ret) {
      printf( "Error in integration of tissue data. %d \n", ret);
      dftEmpty(&data); dftEmpty(&temp); dftEmpty(&input); return(2);
    }
  }
  

    data.voi[ref].y2 = data.voi[ref].y3;            //?
    dftPrint(&data);


  /* Allocate an extra TAC for the bootstrap */
  int bsi=-1;
  if(doBootstrap) {
    ret=dftAddmem(&data, 1); if(ret) {
      printf( "Error: cannot allocate more memory.\n");
      dftEmpty(&data); return(4);
    }
    bsi=data.voiNr;
    strcpy(data.voi[bsi].voiname, "BS");
    strcpy(data.voi[bsi].name, "BS");
  }


  /*
   *  Prepare the room for results
   */
  if(verbose>1) printf("initializing result data\n");
  RES res; resInit(&res);
  ret=res_allocate_with_dft(&res, &data); if(ret!=0) {
    printf( "Error: cannot set-up memory for results.\n");
    dftEmpty(&data); return(4);
  }
  /* Copy titles & filenames */
  tpcProgramName(argv[0], 1, 1, res.program, 256);
//   strcpy(res.datafile, ttacfile);
//   if(inputtype!=5 && rtacfile[0]) strcpy(res.reffile, rtacfile);
  if(ref>=0) strcpy(res.refroi, data.voi[ref].name);
  strcpy(res.fitmethod, "TGO");
  /* Constants */
  res.isweight=data.isweight;
  /* Set data range */
  sprintf(res.datarange, "%g - %g %s", 0.0, fitdur, petTunit(data.timeunit));
  res.datanr=fitframeNr;
  /* Set current time to results */
  res.time=time(NULL);
  /* Set parameter number, including also the extra "parameters"
     and the parameter names and units */
  res.parNr=4;
  {
    int pi;
    pi=0; strcpy(res.parname[pi], "R1"); strcpy(res.parunit[pi], "");
    pi++; strcpy(res.parname[pi], "k2"); strcpy(res.parunit[pi], "1/min");
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
  t=data.x; cr=data.voi[ref].y;
  double refIntegral=data.voi[ref].y3[fitframeNr-1];
  /* Fit model to one TAC at a time */
  for(int ri=0; ri<data.voiNr; ri++) if(ri!=ref) {

    if(verbose>1) printf("Region %d %s\n", ri+1, data.voi[ri].name);
    /* Set data pointers */
    tis=data.voi[ri].y; ct=data.voi[ri].y2; w=data.w;
    double *p=res.voi[ri].parameter;

    /* Set common parameter constraints */
    for(int pi=0; pi<parNr; pi++) {pmin[pi]=def_pmin[pi]; pmax[pi]=def_pmax[pi];}
    /* Set BP limits based on the integrals of roi and ref TACs
       if those were not given in limfile */
    if(refIntegral>0.0) {
      double a, b, c;
      a=(data.voi[ri].y3[fitframeNr-1]/refIntegral);
      if(a<1.0) a=1.0; 
      b=0.0*a; c=5.0*a; 
      pmin[2]=b; pmax[2]=c;
    }
    if(verbose>3) {
      printf("Parameter constraints:\n");
      for(int pi=0; pi<parNr; pi++) printf("  %10.3E - %10.3E\n", pmin[pi], pmax[pi]);
    }

    /* Fit */
    if(verbose>2) printf("  fitting curve...\n");
    TGO_LOCAL_INSIDE=0;
    TGO_SQUARED_TRANSF=0;
    tgoNr=220;
    neighNr=20;
    ret=tgo(pmin, pmax, srtmFunc, NULL, parNr, neighNr, &wss, p, tgoNr, iterNr, verbose-8);
    if(ret>0) {
      printf( "Error in optimization (%d).\n", ret);
      dftEmpty(&data); resEmpty(&res); return(6);
    }
    if(verbose>3) {
      for(int pi=0; pi<parNr; pi++) printf(" %g", p[pi]);
      printf(" -> WSS=%g\n", p[parNr]); 
    }
    /* Correct fitted parameters to match constraints like inside the function */
    (void)modelCheckParameters(parNr, pmin, pmax, p, p, NULL);
    p[parNr]=wss=wss_wo_penalty;
    if(verbose>2) printf("wss := %g\nfitframeNr := %d\n", wss, fitframeNr);



    /* Bootstrap */
    if(doBootstrap) {
      if(verbose>2) printf("  bootstrapping...\n");
      /* bootstrap changes measured and simulated data, therefore use copies */
      tis=data.voi[bsi].y; ct=data.voi[bsi].y2;
      if(doSD) sd=res.voi[ri].sd; else sd=NULL;
      if(doCL) {cl1=res.voi[ri].cl1; cl2=res.voi[ri].cl2;} else cl1=cl2=NULL;

    // FILE *pfile = fopen(debugfile, "a+");
    // fprintf(pfile, "plasma %d",bsi);
    // fclose(pfile);
    // // fprintf(pfile, "plasma %f %f %f %f\n", p[0], data.voi[ri].y[0],data.voi[ri].y2[0],tis[0]);

    //   matrix=(double*)malloc(parNr*bootstrapIter*sizeof(double));
      ret=bootstrapr(
        bootstrapIter, cl1, cl2, sd, p, pmin, pmax, fitframeNr,   // was 0
        // measured and fitted original TAC, not modified
        data.voi[ri].y, data.voi[ri].y2,
        // tissue TAC noisy data is written to be used by objf
        tis, 
        parNr, w, srtmFunc, tmp, verbose-5, bmatrix
      );
// printf("return full sampling matrix... %f %f\n", matrix[0], matrix[100] );
// for(int i=0; i<parNr*bootstrapIter; i++) { bmatrix[i]=matrix[i]; }

      if(ret) {
        printf( "Error in bootstrap: %s\n", tmp);
        for(int pi=0; pi<parNr; pi++) {
          if(doSD) sd[pi]=nan(""); 
          if(doCL) cl1[pi]=cl2[pi]=nan("");
        }
      }
      // back to what pointers were
      tis=data.voi[ri].y; tis=data.voi[ri].y2;
    }

  } /* Next VOI */
  
  if(verbose>0) {printf( "\n");   resPrint(&res);}
  output[0] = res.voi[0].parameter[0];
  output[1] = res.voi[0].parameter[1];
  output[2] = res.voi[0].parameter[2];
  if(doSD) {output[3] = res.voi[0].sd[0];output[4] = res.voi[0].sd[1];output[5] = res.voi[0].sd[2];}
    //  if(doSD) { res->voi[i].cl1[j]; res->voi[i].cl2[j]; }
  /* Delete reference region(s) from the results unless it already existed in data */
  if(inputtype==5) {
    resDelete(&res, ref);
  } else {
    for(int i=data.voiNr-1; i>=0; i--) if(data.voi[i].sw!=0) {
      resDelete(&res, i);
    }
    ref=-1;
  }


  
  resEmpty(&res);
  dftEmpty(&data);
  dftEmpty(&temp);
//   if(doBootstrap) { free(matrix); }
  return(0);
}



/*****************************************************************************
 *
 *  Functions to be minimized
 *
 *****************************************************************************/
double SRTM_IDL::srtmFunc(int parNr, double *p, void *fdata)
{
  int ret;
  double R1, k2, BP, d, wss=0.0;
  double pa[MAX_PARAMETERS], penalty=1.0;


  /* Check parameters against the constraints */
  ret=modelCheckParameters(parNr, pmin, pmax, p, pa, &penalty);
  if(fdata) {}
  /* Get parameters */
  R1=pa[0]; k2=pa[1]; BP=pa[2];

  /* Simulate the tissue PET TAC */
  ret=simSRTM(t, cr, fitframeNr, R1, k2, BP, ct);
  if(ret) {
    printf( "  error %d in simulation\n", ret);
    return(nan(""));
  }

  /* Calculate error */
  for(int i=0; i<fitframeNr; i++) if(w[i]>0.0) {
    d=ct[i]-tis[i]; 
    wss+=w[i]*d*d;
  }
  wss_wo_penalty=wss;
  wss*=penalty;
  if(0) printf("R1=%g  k2=%g  BP=%g  => %g\n", R1, k2, BP, wss);

  return(wss);
}
/*****************************************************************************/