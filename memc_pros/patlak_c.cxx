/** @file patlak.c
 *  @brief Regional Patlak plot.
 *  @details Estimation of the tracer net influx rate from regional PET data
             using Gjedde-Patlak multiple-time graphical analysis. 
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
#include "libtpcsvg.h"
#include "libtpcmodext.h"
#include "meKineticRigid.h"
/*****************************************************************************/

/*****************************************************************************/
#ifndef DEFAULT_LC
#define DEFAULT_LC 1.00
#endif
#ifndef DEFAULT_DENSITY
#define DEFAULT_DENSITY 1.00
#endif
#ifndef BAD_FIT
#define BAD_FIT 9.999E19
#endif
#ifndef MAX_PARAMETERS
#define MAX_PARAMETERS 6
#endif
/*****************************************************************************/



/*****************************************************************************/
/**
 *  Main
 */
extern "C" int patlak_c(unsigned int frameNr, double *t0, double *t1,double *tac, double *ctt,
               double tstart, double tstop, double *output, unsigned int verbose, 
               unsigned int llsq_model,
               unsigned int isweight, double *weights)         //Pure c does not require extern  "C" 
{

  int        ai, help=0, version=0;
  // unsigned int verbose;
  DFT        data, input, temp;
  int        save_stat=1, always_mid=0;
  double     LC=-1.0, Ca=-1.0, density=-1.0;
  double     fixed_Ic=-9.E99;
  int        ri=0, fi, pi, ret, inputtype=0;
  // unsigned int llsq_model;
  char       dfile, ifile, rfile,          
             sfile, tmp[1024], *cptr;          // pfile,  [FILENAME_MAX]
  // double     tstart, tstop;
  double     Ki, KiSD, Ic, IcSD, SWSS;
  double     istart=0.0;
  double     f, xm, ym, xs, ys;
//double    *w, *wx, *wy, *cx, *cy;
  RES        res;

  double    *t, *theta, *dv, *ci, *ici, *ct;
  int        dataNr=0, first, last;
  // double    *t0, *t1, *tac, *ctt, *output, *weights;
  int       voiNr = 1;
  // unsigned int    frameNr, isweight = 0;
  const char *debugfile1 = "debug1.txt";
  // FILE *pfile = fopen(debugfile1, "a+");
  const char *debugfile = "debug.txt";

  dftInit(&data); dftInit(&input); dftInit(&temp); resInit(&res);


  // /* 
  //  make sure interpolation and integration has been done before hand!
  // */
  // frameNr= *(unsigned int*) argv[0];
  // t0     =  (double*) argv[1];
  // t1     =  (double*) argv[2];
  // tac    =  (double*) argv[3];
  // ctt    =  (double*) argv[4];
  // tstart = *(double*) argv[5];
  // tstop  = *(double*) argv[6]; 
  // output =  (double*) argv[7];
  // verbose= *(unsigned int*) argv[8];
  // llsq_model = *(unsigned int*) argv[9]; 
  // isweight =  *(unsigned int*) argv[10];
  // weights = (double*) argv[11];

  // t0 = [0.,0.25,0.5,0.75,1.0,1.5,2.0,3.0,4.0,5.0,10,.0,15.0,20.0,25.0,30.,40.,45.,50.,60,70,80,90,100,110];
  // t1 = [0.25,0.5,0.75,1.0,1.5,2.0,3.0,4.0,5.0,10.0,15.0,20.0,25.0,30.,40.,45.,50.,60,70,80,90,100,110,120];   // from idl
  // tac = [2.99e-4,5.463e-3,4.679e-1,1.523e0,2.443563e+00,3.236325e+00,4.015189e+00,4.897500e+00,5.656930e+00, $
  // 7.241067e+00,9.052315e+00,1.013719e+01,1.076269e+01,1.107761e+01,1.086400e+01,1.065737e+01,1.043218e+01, $
  // 1.008137e+01,9.610046e+00,9.152024e+00,8.711683e+00,8.287338e+00,7.875306e+00];
  // ctt = [0.06415,1.879,84.8,66.15,41.86,26.86,21.36,19.87,17.15,13.54,10.89,9.305,8.147,7.294,6.655,6.167,5.786,5.479,$
  // 5.225,5.008,4.729,4.406,4.114,3.836,3.564,3.295];
  // data.frameNr = 26;

  /* Allocate memory for data */
  if(verbose>1) {printf("allocating memory\n"); }
  if(dftSetmem(&data, frameNr, voiNr)) {
    printf( "out of memory\n");}
  if(dftSetmem(&temp, frameNr, voiNr)) {
    printf("out of memory\n");}
  if(dftSetmem(&input, frameNr, voiNr)) {
    printf("out of memory\n");}

  if(dft_nr_of_NA(&data)>0) {  // check if file contains NAs (missing values)
    printf( "Error: missing values in %s", dfile);
    dftEmpty(&data); return(2);
  }
  if(dft_nr_of_NA(&temp)>0) {  // check if file contains NAs (missing values)
    printf( "Error: missing values in %s", dfile);
    dftEmpty(&temp); return(2);
  }
  if(dft_nr_of_NA(&input)>0) {  // check if file contains NAs (missing values)
    printf("Error: missing values in %s", dfile);
    dftEmpty(&input); return(2);
  }


  /* Set voiNr and frameNr, and type */
  data.voiNr=voiNr; data.frameNr=frameNr;
  data.isweight=isweight;
  data._type=DFT_FORMAT_PLAIN;
  char cnr[1] = {'1'};
  strcpy(cnr, data.studynr);
  char  cunit[6]={'k'};
  strcpy(cunit, data.unit);
  data.timeunit=2;
  data.timetype=3;
  temp.voiNr=voiNr; temp.frameNr=frameNr;
  temp._type=DFT_FORMAT_PLAIN;
  strcpy(cnr, data.studynr);
  strcpy(cunit, temp.unit); 
  temp.timeunit=2;
  temp.timetype=3;
  // if(strlen(data->radiopharmaceutical))
  //   printf("Radiopharmaceutical: %s\n", data->radiopharmaceutical);
  // if(strlen(data->isotope)) printf("Isotope: %s\n", data->isotope);
  // if(strlen(data->scanStartTime))
  //   printf("Scan start time: %s\n", data->scanStartTime);
  // if(strlen(data->injectionTime))
  //   printf("Injection time: %s\n", data->injectionTime);
  // if(data->decayCorrected==DFT_DECAY_CORRECTED)
  //   printf("Corrected for physical decay: yes\n");
  // else if(data->decayCorrected==DFT_DECAY_NOTCORRECTED)
  //   printf("Corrected for physical decay: no\n");
  // printf("_datasize = %d\n", data->_dataSize);
  // if(type==DFT_FORMAT_IFT) {
  //   data._type=DFT_FORMAT_PLAIN;
  // } else if(type==DFT_FORMAT_PMOD) {
  //   data._type=type; // keeping PMOD format
  //   // unless TAC names could not be read
  //   for(i=0; i<data.voiNr; i++) if(strlen(data.voi[i].name)<1) {
  //     data._type=DFT_FORMAT_PLAIN; break;} 
  // } else {
  //   data._type=type;
  // }


  for (int i=0; i<frameNr; i++) { 
    data.x1[i] = *(t0+i);
    data.x2[i] = t1[i];
    data.x[i]=0.5*(data.x1[i]+data.x2[i]);
    data.voi[ri].y[i]= tac[i];

    temp.x1[i] = t0[i];
    temp.x2[i] = t1[i];
    temp.x[i]=0.5*(temp.x1[i]+temp.x2[i]);
    temp.voi[ri].y[i]= ctt[i];
    if(data.isweight) data.w[fi]=weights[i];
    }


 
  double *ti1 = 0;
  double *ti2 = NULL;
  int verifypeak = 0;
  char *status;

  // /* Convert input time units to the same as in tissue data */
  // (void)dftTimeunitConversion(&temp, data.timeunit);
  // /* Check the tissue and plasma TAC concentration units */
  // ret=dftUnitConversion(&temp, petCunitId(data.unit));
  // if(ret!=0) {
  //   sprintf(status, "check the units of input and tissue data");
  // }

  /* Tell user what was the original input time range */
  if(temp.timetype==DFT_TIME_STARTEND) {
    if(ti1!=NULL) *ti1=temp.x1[0];
    if(ti2!=NULL) *ti2=temp.x2[temp.frameNr-1];
  } else {
    if(ti1!=NULL) *ti1=temp.x[0];
    if(ti2!=NULL) *ti2=temp.x[temp.frameNr-1];
  }

  /* Verify the peak if requested */
  if(verifypeak!=0) {
    ret=dftVerifyPeak(&temp, 0, verbose-2, status);
    //if(ret!=0) sprintf(status, "input TAC should start at time zero");
    if(ret>0) {dftEmpty(&temp); return 101;}
  }

  /* Interpolate and integrate data to pet times */
  ret=dftInterpolate(&temp, &data, &input, status, verbose);
  dftEmpty(&temp);
  if(ret!=0) return 4; if(ret==2) {printf("nothing to be done for interpolation!");}

  if(verbose>9) {
    printf("\nIDL input data:\n");
    dftPrint(&temp);
    printf("\nInput data:\n");
    dftPrint(&input);
    printf("\nTissue data:\n");
    dftPrint(&data);
  }


  /* Finishing reading data */
    // theta = (double *)malloc(frameNr*sizeof(double));
    // dv    = (double *)malloc(frameNr*sizeof(double));
    // ct    = (double *)malloc(frameNr*sizeof(double));
    // ci    = (double *)malloc(frameNr*sizeof(double));
    // ici   = (double *)malloc(frameNr*sizeof(double));
    // t     = (double *)malloc(frameNr*sizeof(double));
    // theta = (double *)malloc(frameNr*sizeof(double));

    // for (int i=0; i<data.frameNr; i++) { 
    //   ct[i]=data.voi[ri].y[i];
    //   ci[i] = input.voi[ri].y[i];
    //   ici[i]=input.voi[ri].y2[i];
    //   t[i]=input.x[i];
    //   theta[i] = data.voi[ri].y2[i];  
    //   dv[i] = data.voi[ri].y3[i];
    // }
    ct  = data.voi[ri].y;
    ci  = input.voi[ri].y;
    ici = input.voi[ri].y2;
    t   = input.x;
    theta = data.voi[ri].y2;  
    dv  = data.voi[ri].y3;

    if (verbose>9) {
        // FILE *ppfile = fopen(debugfile, "a+");
        printf( "CT %f %f supplied\n", ct[0],ct[frameNr-1]);
        printf( "CI %f %f supplied\n", ci[0],ci[frameNr-1]);
        printf( "ici %f %f supplied\n", ici[0],ici[frameNr-1]);
        printf( "theta %f %f supplied\n", theta[0],theta[frameNr-1]);
        printf( "dv %f %f supplied\n", dv[0],dv[frameNr-1]);
        printf( "t %f %f supplied\n", t[0],t[frameNr-1]);        
        // fclose(ppfile);
    }


  if(data.frameNr==1 && fixed_Ic<=-1.E99) {
    /* If static study, then automatically set IC to zero */
    fixed_Ic=0.0;
    if(verbose>=0)
      printf("Suggestion: for FUR calculation use regfur.");
  }

  if(dft_nr_of_NA(&data)>0) {  // check if file contains NAs (missing values)
    printf( "Error: missing values in %s", dfile);
    dftEmpty(&data); return(2);
  }

  if(data.frameNr==1 && data.timetype==DFT_TIME_MIDDLE) 
    data.x2[0]=data.x1[0]=data.x[0];
  // like if we had only frame mid times
  if(always_mid!=0) data.timetype=DFT_TIME_MIDDLE;


//   /* Sort the samples by time in case data is catenated from several curves */
//   (void)dftSortByFrame(&data);

//   /* Make sure that there is no overlap in frame times */
//   if(data.timetype==DFT_TIME_STARTEND) {
//     if(verbose>1) fprintf(pfile, "checking frame overlap in %s\n", dfile);
//     ret=dftDeleteFrameOverlap(&data);
//     if(ret) {
//       fprintf(pfile, "Error: %s has overlapping frame times.\n", dfile);
//       dftEmpty(&data);  return(2);
//     }
//   }


  /* Set time unit to min */
  ret=dftTimeunitConversion(&data, TUNIT_MIN);
  if(ret) printf( "Warning: check that regional data times are in minutes.");
  /* Get and check fit time range */
  dataNr=fittime_from_dft(&data, &tstart, &tstop, &first, &last, verbose-8);
  if(verbose>2) {
    printf("dataNr_in_range := %d\n", dataNr);
    printf("first_in_range := %d\n", first);
    printf("last_in_range := %d\n", last);
  }
  if(dataNr<1) {
    printf( "Error: data does not contain the specified time range.");
    dftEmpty(&data); return(2);
  } else if(dataNr<2 && fixed_Ic<=-1.E99) {
    printf( "Error: cannot make plot from less than 2 points.");
    dftEmpty(&data); return(2);
  } else if(dataNr==2 && fixed_Ic<=-1.E99) {
    printf( "Warning: only two samples in the time range.");
  }
  if(verbose>2) {
    printf("dataNr := %d\n", dataNr);
    printf("tstart := %g\ntstop := %g\n", tstart, tstop);
    printf("first := %d\nlast := %d\n", first, last);
  }


//   if(verbose>9) {
//     printf("\nInput data:\n");
//     dftPrint(&input);
//     printf("\nTissue data:\n");
//     dftPrint(&data);
//   }
//   if(inputtype==5) { // Reference region name was given
//     if(verbose>0) fprintf(pfile, "selected reference region := %s\n",
//                           input.voi[0].name);
//     for(ri=1; ri<input.voiNr; ri++)
//       fprintf(pfile, "Warning: reference region %s unused.\n",
//               input.voi[ri].name);
//   } else {
//     if(input.voiNr>1)
//       fprintf(pfile, "Warning: only the first of input curves is used.\n");
//   }

//   /* Check that original input data started early enough, otherwise AUC(0-T)
//      might be wrong */ 
//   if(istart>0.3) {
//     fprintf(pfile, "Warning: input TAC should start at time zero.\n");
//   }

  // /*
  //  *  Prepare the room for results
  //  */
  // if(verbose>1) printf("initializing result data\n");
  // ret=res_allocate_with_dft(&res, &data); if(ret!=0) {
  //   fprintf(pfile, "Error: cannot setup memory for results.\n");
  //   dftEmpty(&input); dftEmpty(&data);  return(4);
  // }
  // /* Copy titles & filenames */
  // tpcProgramName(argv[0], 1, 1, res.program, 128);
  // if(verbose>10) printf("res.program := '%s'\n", res.program);
  // strcpy(res.datafile, dfile);
  // if(inputtype==5) strcpy(res.refroi, input.voi[0].name);
  // else strcpy(res.plasmafile, ifile);
  // if(strlen(res.studynr)==0 || strcmp(res.studynr, ".")==0)
  //   studynr_from_fname2(dfile, res.studynr, 1);
  // /* Constants */
  // if(Ca>0.0) {res.density=density; res.lc=LC; res.concentration=Ca;}
  // /* Set data range */
  // sprintf(res.datarange, "%g - %g %s", tstart, tstop, petTunit(data.timeunit));
  // res.datanr=dataNr;
  // if(llsq_model==0) strcpy(res.fitmethod, "Traditional regression model");
  // else if(llsq_model==1) strcpy(res.fitmethod, "Iterative method");
  // else if(llsq_model==2) strcpy(res.fitmethod, "Perpendicular regression model");
  // else if(llsq_model==3) strcpy(res.fitmethod, "Median of two-point slopes");
  // /* Set current time to results */
  // res.time=time(NULL);
  // res.isweight=0; /*data.isweight;*/
  // /* Set parameter number, including also the extra "parameters" */
  // /* Set the parameter names */
  // res.parNr=3; if(Ca>0.0) res.parNr++;
  // pi=0;
  // if(Ca>0.0) {
  //   strcpy(res.parname[pi], "MR");
  //   strcpy(res.parunit[pi], petCunit(CUNIT_UMOL_PER_MIN_PER_100G));
  //   pi++;
  // }
  // strcpy(res.parname[pi], "Ki");
  // strcpy(res.parunit[pi], petCunit(CUNIT_ML_PER_ML_PER_MIN));
  // pi++;
  // strcpy(res.parname[pi], "Ic");
  // strcpy(res.parunit[pi], petCunit(CUNIT_ML_PER_ML));
  // pi++;
  // if(llsq_model==1 || llsq_model==3) {
  //   strcpy(res.parname[pi], "SqrtWSS");
  //   strcpy(res.parunit[pi], data.unit);
  // } else if(llsq_model==0) {
  //   strcpy(res.parname[pi], "r");
  //   strcpy(res.parunit[pi], petCunit(CUNIT_UNITLESS));
  // } else if(llsq_model==2) {
  //   strcpy(res.parname[pi], "SSD");
  //   strcpy(res.parunit[pi], petCunit(CUNIT_UNITLESS));
  // }




  /*
   *  Allocate memory for llsqwt() and weights
   */
  double w[data.frameNr];
  double wx[data.frameNr];
  double wy[data.frameNr];
  double cx[data.frameNr];
  double cy[data.frameNr];


  /*
   *  Calculate Ki for each region
   */

  /* One region at a time */
  //   for(ri=0; ri<data.voiNr; ri++) {
    // ri = 0;
    if(verbose>0) printf("calculating %s\n", data.voi[ri].name);


    /* Set data pointers */
    theta=data.voi[ri].y2;
    dv=data.voi[ri].y3;
    ci=input.voi[0].y; ici=input.voi[0].y2;
    ct=data.voi[ri].y; t=input.x;

    /* Calculate Patlak plot data */
    for(fi=data.frameNr-1; fi>=0; fi--) if(ci[fi]!=0.0) {
      if(verbose>8) {
        printf("%03d %8.3f : ici=%g ci=%g ct=%g\n",
          fi+1, t[fi], ici[fi], ci[fi], ct[fi] );
      }
      /* theta (x axis) */
      theta[fi]=ici[fi]/ci[fi]; wx[fi]=1.0;
      //if(theta[fi]>=0.0) wx[fi]=1.0; else wx[fi]=0.0;
      /* dv (y axis) */
      dv[fi]=ct[fi]/ci[fi]; wy[fi]=1.0;
      // check the close-to-zeroes in the first frames
      if(data.x[fi]<0.1*data.x[data.frameNr-1]) { 
        if(theta[fi]>theta[data.frameNr-1] || dv[fi]>dv[data.frameNr-1]) {
          if(verbose>2)
            printf("Possible close-to-zero plot point at %g -> set to zero.\n",
                   data.x[fi]);
          theta[fi]=dv[fi]=wx[fi]=wy[fi]=0.0;
        }
      }
    } else if(fixed_Ic>-1.E99) { // this only if input concentration is zero
      theta[fi]=ici[fi]; dv[fi]=ct[fi]; wx[fi]=wy[fi]=1.0;
    } else theta[fi]=dv[fi]=wx[fi]=wy[fi]=0.0;

    /* Set x weight to 0, if integral is still <=0, whether weighted or not */
    for(fi=input.frameNr-1; fi>=0; fi--) if(input.voi[0].y2[fi]<=0.0) break;
    for(; fi>=0; fi--) wx[fi]=0.0;

    if(verbose>6) {
      for(fi=first; fi<=last; fi++)
        printf("%03d %8.3f : %g %g  (%g %g)\n",
          fi+1, data.x[fi], theta[fi], dv[fi], wx[fi], wy[fi] );
    }



 /* Fit */    

    KiSD=Ki=Ic=IcSD=SWSS=0.0;

    if(fixed_Ic>-1.E99) { /* y axis is constrained to fixed_Ic */
 
      /* Calculate the means and SDs of plot data */
      ret=mean(&theta[first], &dv[first], dataNr, &xm, &xs, &ym, &ys);
      /* Calculate the slope through constrained intercept and plot mean */
      if(fixed_Ic>-1.E99) Ic=fixed_Ic; else Ic=0.0;
      Ki=(ym-Ic)/xm; if(xm!=0.0) KiSD=ys/xm; SWSS=1.0;

    } else if(llsq_model==1) {

      /* Calculation of LLSQ fit with errors in both coordinates */
      ret=llsqwt(
        /* input */
        &theta[first], &dv[first], dataNr, &wx[first], &wy[first],
        1.0e-10, /* allowed tolerance in slope estimation */
        /* input/output */
        &w[first], /* work vector; effective weights w[i] are returned in it */
        /* output ; SWSS is already normalized */
        &Ic, &Ki, &SWSS, &IcSD, &KiSD, cx, cy
      );
      if(verbose>5) {
        printf("%s:\n", data.voi[ri].name);
        for(fi=first; fi<=last; fi++)
          printf("%03d %8.3f : %g %g  (%g %g -> %g)\n",
            fi+1, data.x[fi], theta[fi], dv[fi], wx[fi], wy[fi], w[fi] );
      }

    } else if(llsq_model==2) {

      /* Remove plot points with zero weight */
      for(fi=first; fi<=last; fi++)
       if(wx[fi]<=0.0 || wy[fi]<=0.0) theta[fi]=dv[fi]=nan("");

      /* Calculation of perpendicular line fit */
      ret=llsqperp3(
        /* input */
        &theta[first], &dv[first], dataNr,
        /* output ; SSD is already normalized */
        &Ki, &Ic, &SWSS
      );

    } else if(llsq_model==0) {

      /* Remove plot points with zero weight */
      for(fi=first; fi<=last; fi++)
       if(wx[fi]<=0.0 || wy[fi]<=0.0) theta[fi]=dv[fi]=nan("");

      if(verbose>9) for(fi=first; fi<=last; fi++)
        printf(" %d  %g  %g\n", fi, theta[fi], dv[fi]);

      /* Calculation of linear regression using pearson() */
      ret=pearson3(
        /* input */
        &theta[first], &dv[first], dataNr,
        /* output */
        &Ki, &KiSD, &Ic, &IcSD, &SWSS, &f
      );
      if(verbose>9) printf("Ki=%g Ic=%g\n", Ki, Ic);

    } else if(llsq_model==3) {

      /* Remove plot points with zero weight */
      for(fi=first; fi<=last; fi++)
       if(wx[fi]<=0.0 || wy[fi]<=0.0) theta[fi]=dv[fi]=nan("");

      /* Calculation of median slope and ic */
      ret=medianline(
        &theta[first], &dv[first], dataNr, &Ki, &Ic
      );
    }

    // if(ret==0) {
    //   res.voi[ri].parameter[0]=Ki; if(save_stat) res.voi[ri].sd[0]=KiSD;
    //   res.voi[ri].parameter[1]=Ic; if(save_stat) res.voi[ri].sd[1]=IcSD;
    //   res.voi[ri].parameter[2]=SWSS;
    //   if(verbose>1) printf("Ki := %g (%g)\n", Ki, KiSD);

    // } else
    //   fprintf(pfile, "Error (%d) in linear fit of %s\n", ret, data.voi[ri].name);
 //   } /* next region */


  /* return results */
  output[0] = Ki;
  output[1] = Ic;
  output[2] = KiSD;
  output[3] = IcSD;
  output[4] = SWSS;           //correlation coefficient

  // clean memory! dftEmpty(&temp);
  resEmpty(&res); dftEmpty(&input); dftEmpty(&data); 


//  fclose(pfile);
  //  free(theta);
  //  free(dv);
  //  free(ci)  ;
  //  free(ici)  ;
  //  free(ct)  ;
  //  free(t)  ;  
}
