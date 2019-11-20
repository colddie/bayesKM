/** @file regfur.c
 *  @brief Estimation of FUR from regional PET TAC data.
 *  @details Fraction Uptake Rate (FUR) is related to the tracer net influx 
 *           rate (Ki) calculated using Patlak multiple-time graphical analysis. 
 *  @copyright (c) Turku PET Centre
 *  @author Vesa Oikonen
 */
/// @cond
/******************************************************************************/
#include "tpcclibConfig.h"
/******************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
/******************************************************************************/
#include "libtpcmisc.h"
#include "libtpcmodel.h"
#include "libtpccurveio.h"
#include "libtpcmodext.h"
/*****************************************************************************
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
#ifndef DEFAULT_LC
#define DEFAULT_LC 1.00
#endif
#ifndef DEFAULT_DENSITY
#define DEFAULT_DENSITY 1.00
#endif
/*****************************************************************************/

/*****************************************************************************/
/**
 *  Main
 */
int regfur_idl(int argc, char **argv)
{
  int        ai, help=0, version=0, verbose=1;
  int        ret, ri=0;
  char       inpfile[FILENAME_MAX], petfile[FILENAME_MAX];
  char       curfile[FILENAME_MAX], outfile[FILENAME_MAX], tmp[1024], *cptr;
  DFT        input, auc, data, avg;
  RES        res;
  double     startTime=-1.0, endTime=-1.0, aucTime=-1.0;
  double     LC=-1.0, Ca=-1.0, density=-1.0;
  int        fur_mode=0; // 0=traditional Ct/iCp, 1=derivative dCt/Cp
  int        voiNr = 1;
  double     *t0, *t1, *tac, *ctt, *output;
  // int        tstart, tstop;
  unsigned int    frameNr;
  const char *debugfile = "debug.txt";

  dftInit(&data); dftInit(&input); dftInit(&auc); dftInit(&avg); //resInit(&res);

  /* 
   make sure interpolation and integration has been done before hand!
  */
  frameNr= *(unsigned int*) argv[0];
  t0     =  (double*) argv[1];
  //t1     =  (double*) argv[2];
  tac    =  (double*) argv[2];
  ctt    =  (double*) argv[3];
  startTime = *(double*) argv[4];
  endTime  = *(double*) argv[5]; 
  output =  (double*) argv[6];
  verbose= *(unsigned int*) argv[7];
  fur_mode = *(unsigned int*) argv[8]; 


  /* Allocate memory for data */
  if(verbose>1) {printf("allocating memory\n"); }
  if(dftSetmem(&data, frameNr, voiNr)) {
    printf( "out of memory\n");}
  if(dftSetmem(&auc, frameNr, voiNr)) {
    printf("out of memory\n");}
  if(dftSetmem(&avg, frameNr, voiNr)) {
    printf("out of memory\n");}
  if(dftSetmem(&input, frameNr, voiNr)) {
    printf("out of memory\n");}

  if(dft_nr_of_NA(&data)>0) {  // check if file contains NAs (missing values)
    printf( "Error: missing values in ");
    dftEmpty(&data); return(2);
  }
  if(dft_nr_of_NA(&auc)>0) {  // check if file contains NAs (missing values)
    printf( "Error: missing values in ");
    dftEmpty(&auc); return(2);
  }
  if(dft_nr_of_NA(&avg)>0) {  // check if file contains NAs (missing values)
    printf( "Error: missing values in ");
    dftEmpty(&avg); return(2);
  }
  if(dft_nr_of_NA(&input)>0) {  // check if file contains NAs (missing values)
    printf("Error: missing values in ");
    dftEmpty(&input); return(2);
  }


  /* Set voiNr and frameNr, and type */
  data.voiNr=voiNr; data.frameNr=frameNr;
  data._type=DFT_FORMAT_PLAIN;
  char cnr[1] = "1";
  strcpy(cnr, data.studynr);
  char  cunit[6]="kBq/mL";
  strcpy(cunit, data.unit);
  data.timeunit=2;
  data.timetype=3;
  input.voiNr=voiNr; input.frameNr=frameNr;
  input._type=DFT_FORMAT_PLAIN;
  strcpy(cnr, data.studynr);
  strcpy(cunit, input.unit); 
  input.timeunit=2;
  input.timetype=3;


  for (int i=0; i<frameNr; i++) { 
    data.x1[i] = t0[i];
    data.x2[i] = t0[i];
    data.x[i]=0.5*(data.x1[i]+data.x2[i]);
    data.voi[ri].y[i]= tac[i];

    input.x1[i] = t0[i];
    input.x2[i] = t0[i];
    input.x[i]=0.5*(input.x1[i]+input.x2[i]);
    input.voi[ri].y[i]= ctt[i];
    }



  // /* Sort the samples by time in case data is catenated from several curves */
  // (void)dftSortByFrame(&data);

  // /* Make sure that there is no overlap in frame times */
  // if(data.timetype==DFT_TIME_STARTEND) {
  //   if(verbose>2) printf( "checking frame overlap in %s\n", petfile);
  //   ret=dftDeleteFrameOverlap(&data);
  //   if(ret) {
  //     printf( "Error: %s has overlapping frame times.\n", petfile);
  //     dftEmpty(&data); return(2);
  //   }
  // }

  /* Set time unit to min */
  ret=dftTimeunitConversion(&data, TUNIT_MIN);
  if(ret) printf( "Warning: check that regional data times are in minutes.\n");
  /* If user did not specify start and end times, then get those from data */
  if(endTime<=1.0E-02) {
    if(data.timetype==DFT_TIME_STARTEND) {
      startTime=data.x1[0]; endTime=data.x2[data.frameNr-1];
    } else {
      startTime=data.x[0]; endTime=data.x[data.frameNr-1];
    }
    if(verbose>1) {
      printf("startTime := %g min\n", startTime);
      printf("endTime := %g min\n", endTime);
    }
  }



  // /*
  //  *  Read input data
  //  */
  // if(verbose>1) printf("Reading input file %s\n", inpfile);
  // if(dftRead(inpfile, &input)) {
  //   printf("Error in reading '%s': %s\n", inpfile, dfterrmsg);
  //   dftEmpty(&data); dftEmpty(&avg); return(4);
  // }
  // if(input.voiNr>1) {
  //   printf( "Warning: only first TAC is used as input.\n");
  //   input.voiNr=1;
  // }
  // if(dft_nr_of_NA(&input)>0) {  // check if file contains NAs (missing values)
  //   printf( "Error: missing values in %s\n", inpfile);
  //   dftEmpty(&data); dftEmpty(&avg); dftEmpty(&input); return(4);
  // }

  // /* Sort the samples by time in case data is catenated from several curves */
  // (void)dftSortByFrame(&data);

  // /* Set time unit to min */
  // ret=dftTimeunitConversion(&input, TUNIT_MIN);
  // if(ret) printf( "Warning: check that input times are in minutes.\n");

  // /*
  //  *  Check the regional and plasma TAC concentration units
  //  */
  // ret=dftUnitConversion(&input, dftUnitId(data.unit));
  // if(ret!=0) {
  //   printf( "Warning: check the units of input and regional data.\n");
  // }


  if(verbose>9) {
    printf("\nInput data:\n");
    dftPrint(&input);
    printf("\nTissue data:\n");
    dftPrint(&data);
  }

  /*
   *  Calculate the average (or slope) over specified time range
   */
  if(verbose>1) printf("calculating average\n");
  ret=dftTimeIntegral(&data, startTime, endTime, &avg, 1, tmp, verbose-3);
  if(ret!=0) {
    printf("Error: %s.\n", tmp);
    if(verbose>2)
      printf("dftTimeIntegral(data, %g, %g, avg, 1, tmp) := %d\n", startTime, endTime, ret);
    dftEmpty(&data); dftEmpty(&avg); return(3);
  }
  if(verbose>1) printf("%s.\n", tmp);
  if(fur_mode==1) {
    // int ri;
    double k, ksd, b, bsd, r, ysd;
    if(verbose>0) {
      printf("calculating slope\n"); fflush(stdout);}
    // for(ri=0; ri<data.voiNr; ri++) {
      ret=pearson4(data.x, data.voi[ri].y, data.frameNr, startTime, endTime,
                  &k, &ksd, &b, &bsd, &r, &ysd);
      // if(ret!=0) break;
      avg.voi[ri].y[0]=k;
    // }
    if(ret!=0) {
      printf( "Error: tissue slope calculation not successful.\n");
      if(verbose>1) printf("      ret := %d\n", ret);
      dftEmpty(&data); dftEmpty(&avg); return(3);
    }
  }
  if(verbose>2) {
    printf("Regional tissue value or derivative\n");
    for(int ri=0; ri<avg.voiNr; ri++)
      printf("%s : %g\n", avg.voi[ri].name, avg.voi[ri].y[0]);
  }



  /*
   *  Calculate and save FUR curve, if required;
   *  this is not used in calculation of regional FUR
   */
  // if(curfile[0]) {
    if(verbose>1) {printf("calculating FUR curve\n"); fflush(stdout);}
    DFT fur;
    int fi;   //, ri;
    dftInit(&fur); ret=dftdup(&data, &fur);
    if(ret!=0) {
      printf( "Warning: cannot allocate memory for FUR curves\n");
      // goto failed;
    }
    fur.frameNr=0;
    for(fi=0; fi<data.frameNr; fi++) {
      if(data.x[fi]<startTime || data.x[fi]>endTime) continue; 
      /* AUC 0-t for each frame */
      if(data.x[fi]>0.0) {
        ret=dftTimeIntegral(&input, 0.0, data.x[fi], &auc, 0, tmp, verbose-4);
        if(ret!=0) {printf("Warning (%d): %s\n", ret, tmp); break;}
        if(auc.voi[0].y[0]<1.0E-006) continue;
      } else continue;
      fur.x1[fur.frameNr]=data.x1[fi]; fur.x2[fur.frameNr]=data.x2[fi];
      fur.x[fur.frameNr]=data.x[fi]; fur.w[fur.frameNr]=data.w[fi];
      /* Divide each region by AUC */
      // for(ri=0; ri<fur.voiNr; ri++)
        fur.voi[ri].y[fur.frameNr]= data.voi[ri].y[fi] / auc.voi[0].y[0];
      fur.frameNr++;
    }
    // if(ret!=0) goto failed;
    // /* Write FUR curve */
  // }


  /*
   *  Calculate input integral from 0 to PET middle time point,
   *  or input value at PET middle time point
   */
  if(aucTime<=0.0) aucTime=0.5*(startTime+endTime);
  if(fur_mode==0)
    ret=dftTimeIntegral(&input, 0.0, aucTime, &auc, 0, tmp, verbose-3);
  else
    ret=dftTimeIntegral(&input, startTime, endTime, &auc, 1, tmp, verbose-3);
  if(ret!=0) {
    printf("Error (%d): %s\n", ret, tmp);
    dftEmpty(&avg); dftEmpty(&input); dftEmpty(&auc); dftEmpty(&data); 
    return(6);
  }
  if(verbose>1) {
    if(fur_mode==0)
      printf("AUC[%g-%g] := %g\n", auc.x1[0], auc.x2[0], auc.voi[0].y[0]);
    else
      printf("Input[%g-%g] := %g\n", auc.x1[0], auc.x2[0], auc.voi[0].y[0]);
  }


  /*
   *  Divide average regional data by input AUC, or
   *  regional concentration derivative by input
   */
  // for(int ri=0; ri<avg.voiNr; ri++) 
  avg.voi[ri].y[0]/=auc.voi[0].y[0];
  dftUnitToDFT(&avg, CUNIT_PER_MIN);

  if(verbose>9) {
    printf("\nAVG:\n");
    dftPrint(&avg);
    printf("\nAUC:\n");
    dftPrint(&auc);
    printf("\nFUR:\n");
    dftPrint(&fur);
  }

  // read output 
  output[0] = avg.voi[ri].y[0];
  for (int i=1;i<fur.frameNr+1;i++) {
    output[i] = fur.voi[ri].y[i-1];
  }





  // /*
  //  *  Calculate metabolic rate, if necessary
  //  */
  // if(Ca>0.0) {
  //   double MRf;
  //   MRf=100.*Ca/(density*LC);
  //   if(verbose>1)
  //     printf( "converting FUR to metabolic rate with factor %g\n", MRf);
  //   for(int ri=0; ri<avg.voiNr; ri++) avg.voi[ri].y[0]*=MRf;
  //   dftUnitToDFT(&avg, CUNIT_UMOL_PER_MIN_PER_100G);
  // }


  // /*
  //  *  Save FUR/MR data
  //  */
  // ret=backupExistingFile(outfile, NULL, tmp); if(ret!=0) {
  //   printf("Error: %s\n", tmp); 
  //   dftEmpty(&avg); dftEmpty(&data); return(11);
  // }
  // if(verbose>2) printf("%s\n", tmp);
  // /* If filename extension is .dft, .csv, or .dat, then save in DFT, CSV-UK, or 
  //    simple format, respectively, otherwise as RES */ 
  // cptr=strrchr(outfile, '.'); 
  // if(cptr!=NULL && strcasecmp(cptr, ".DFT")==0) {
  //   dftSetComments(&avg);
  //   if(dftWrite(&avg, outfile)) {
  //     printf( "Error in writing %s: %s\n", outfile, dfterrmsg);
  //     dftEmpty(&avg); dftEmpty(&data); return(12);
  //   }
  // } else if(cptr!=NULL && strcasecmp(cptr, ".CSV")==0) {
  //   avg._type=DFT_FORMAT_CSV_UK;
  //   avg.timetype=DFT_TIME_MIDDLE;
  //   if(dftWrite(&avg, outfile)) {
  //     printf( "Error in writing %s: %s\n", outfile, dfterrmsg);
  //     dftEmpty(&avg); dftEmpty(&data); return(12);
  //   }
  // } else if(cptr!=NULL && strcasecmp(cptr, ".DAT")==0) {
  //   avg._type=DFT_FORMAT_PLAIN;
  //   avg.timetype=DFT_TIME_MIDDLE;
  //   if(dftWrite(&avg, outfile)) {
  //     printf( "Error in writing %s: %s\n", outfile, dfterrmsg);
  //     dftEmpty(&avg); dftEmpty(&data); return(12);
  //   }
  // } else {
  //   /* Copy DFT data into RES */
  //   ret=dftToResult(&avg, &res, tmp);
  //   if(ret!=0) {
  //     printf("Error: %s.\n", tmp);
  //     dftEmpty(&avg); resEmpty(&res); dftEmpty(&data); return(13);
  //   }
  //   /* Set more header information */
  //   tpcProgramName(argv[0], 1, 1, res.program, 256);
  //   cptr=strrchr(petfile, '/'); if(cptr==NULL) cptr=strrchr(petfile, '\\');
  //   if(cptr==NULL) cptr=petfile; else cptr++; strcpy(res.datafile, cptr);
  //   cptr=strrchr(inpfile, '/'); if(cptr==NULL) cptr=strrchr(inpfile, '\\');
  //   if(cptr==NULL) cptr=inpfile; else cptr++; strcpy(res.plasmafile, cptr);
  //   if(Ca>0.0) {
  //     res.concentration=Ca;
  //     res.lc=LC;
  //     res.density=density;
  //   }
  //   if(Ca>0.0) strcpy(res.parname[0], "MR"); else strcpy(res.parname[0], "FUR");
  //   /* Save RES */
  //   ret=resWrite(&res, outfile, verbose-3); 
  //   if(ret) {
  //     printf( "  Error (%d) in writing file %s\n", ret, outfile);
  //     dftEmpty(&avg); resEmpty(&res); dftEmpty(&data); return(11);
  //   }
  //   resEmpty(&res);
  // }
  // if(verbose>0) {
  //   if(Ca<=0.0) printf( "FUR(s) saved in %s.\n", outfile);
  //   else printf( "MRs saved in %s.\n", outfile);
  // }

  /* Free memory */
  /* Input AUC is not needed anymore */
  /* Original input curve is not needed anymore */
  dftEmpty(&input);
  dftEmpty(&auc);
  dftEmpty(&avg);
  dftEmpty(&data);

  return(0);
}
/*****************************************************************************/

/*****************************************************************************/
/// @endcond
