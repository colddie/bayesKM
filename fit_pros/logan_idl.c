/** @file logan.c
 *  @brief Regional Logan plot.
 *  @details Estimation of the tracer distribution volume or distribution
            volume ratio from regional PET data using Logan multiple-time
            graphical analysis. 
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


/******************************************************************************/
/* Local functions */
int best_logan_reed(
  double *x, double *y, double *wx, double *wy, int nr, int min_nr,
  double *slope, double *ic, double *nwss, double *sslope, double *sic,
  double *cx, double *cy, int *bnr, int verbose
);
int best_logan_regr(
  double *x, double *y, double *wx, double *wy, int nr, int min_nr,
  double *slope, double *ic, double *r, double *sslope, double *sic, int *bnr,
  int verbose
);
/******************************************************************************/


/*****************************************************************************/
/**
 *  Main
 */
int logan_idl(int argc,  float * argv[])         //Pure c does not require extern  "C" 
{
  int        ai, help=0, version=0;
  unsigned int verbose;
  DFT        data, input, temp;
  int        save_stat=1, always_mid=0;
  double     LC=-1.0, Ca=-1.0, density=-1.0;
  double     fixed_Ic=-9.E99;
  int        ri=0, fi, pi, ret, inputtype=0;
  unsigned int llsq_model;
  char       dfile, ifile, rfile,          
             sfile, tmp[1024], *cptr;          // pfile,  [FILENAME_MAX]
  double     tstart, tstop, DV, DVSD, Ic, IcSD, SWSS;
  double     istart=0.0;
  double     f, xm, ym, xs, ys, k2=-1.0;
  RES        res;
  int        bp_type; // 0=no, 1=DVR, 2=BPnd, 3=BPp
  int        dvr_minus_one; // 0=DVR reported, 1=BPnd reported with ref input

  double    *t, *theta, *dv, *ci, *ici, *ct, *ict;
  int        dataNr=0, first, last;
  double    *t0, *t1, *tac, *ctt, *output, *weights;
  int       voiNr = 1;
  unsigned int    frameNr, isweight = 0, logan_mode = 0;
  const char *debugfile1 = "debug1.txt";
  // FILE *pfile = fopen(debugfile1, "a+");
  const char *debugfile = "debug.txt";

  dftInit(&data); dftInit(&input); dftInit(&temp); resInit(&res);


  // dfile = 'dynamic_1tac.tac';
  // pfile = '';
  // rfile = '';
  // ifile = 'plasma.bld';
  // sfile = '';
  // tstart = 30;
  // tstop = 90;

  /* 
   make sure interpolation and integration has been done before hand!
  */
  frameNr= *(unsigned int*) argv[0];
  t0     =  (double*) argv[1];
  t1     =  (double*) argv[2];
  tac    =  (double*) argv[3];
  ctt    =  (double*) argv[4];
  tstart = *(double*) argv[5];
  tstop  = *(double*) argv[6]; 
  output =  (double*) argv[7];
  verbose= *(unsigned int*) argv[8];
  llsq_model = *(unsigned int*) argv[9]; 
  k2     = *(double*) argv[10]; 
  isweight =  *(unsigned int*) argv[11];
  weights = (double*) argv[12];
  logan_mode = *(unsigned int*) argv[13];
  // bp_type= *(unsigned int*) argv[11]; 
  // dvr_minus_one = *(unsigned int*) argv[12]; 


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
    printf( "out of memory\n"); return(1);}
  if(dftSetmem(&temp, frameNr, voiNr)) {
    printf("out of memory\n"); return(1);}
  if(dftSetmem(&input, frameNr, voiNr)) {
    printf("out of memory\n"); return(1);}

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
  char cnr[1] = "1";
  strcpy(cnr, data.studynr);
  char  cunit[6]="kBq/mL";
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
  // if(ret!=0) return 4; if(ret==2) printf('nothing to be done for interpolation!')

  if(verbose>9) {
    printf("\nIDL input data:\n");
    dftPrint(&temp);
    printf("\nInput data:\n");
    dftPrint(&input);
    printf("\nTissue data:\n");
    dftPrint(&data);
  }

  if(inputtype==5) { // Reference region name was given
    if(verbose>0) fprintf(stdout, "selected reference region := %s\n",
                          input.voi[0].name);
    for(ri=1; ri<input.voiNr; ri++)
      fprintf(stderr, "Warning: reference region %s unused.\n",
        input.voi[ri].name);
  } else {
    if(input.voiNr>1)
      fprintf(stderr, "Warning: only the first of input curves is used.\n");
  }
  // if(inputtype!=5 && k2>0.0) {
  //   fprintf(stderr, 
  //     "Error: reference region k2 must be used only with reference tissue.\n");
  //   dftEmpty(&input); dftEmpty(&data); return(3);
  // }
  
  /* Check that original input data started early enough, otherwise AUC(0-T)
     might be wrong */ 
  if(istart>0.3) {
    printf("Warning: input TAC should start at time zero.\n");
  }


 /* Integrate tissue data */
  if(verbose>1) printf("integrating tissue data\n");
  // for(ri=0; ri<data.voiNr; ri++) {
    if(data.timetype==DFT_TIME_STARTEND && always_mid==0)
      ret=petintegral(data.x1, data.x2, data.voi[ri].y, data.frameNr,
                      data.voi[ri].y2, NULL);
    else
      ret=integrate(data.x, data.voi[ri].y, data.frameNr, data.voi[ri].y2);

    
    if(ret) {
      printf( "Error in integration of tissue data. %d \n", ret);
      dftEmpty(&data); dftEmpty(&input); return(2);
    }
  // }

if (verbose>1) {
    printf("\nTissue data:\n");
    dftPrint(&data);
    printf("data.timetype %d DFT_TIME_STARTEND %d", data.timetype, DFT_TIME_STARTEND);
}
  

//  /*
//    *  Find the reference ROI for DVR calculation with plasma input
//    */
//   if(dvrname[0]) {
//     if(verbose>1) printf("selecting reference curve\n");
//     n=dftSelectRegions(&data, dvrname, 1);
//     if(n<=0) {        /* no vois are matching */
//       fprintf(stderr, "Error: Cannot find ref voi '%s'.\n", dvrname);
//       dftEmpty(&input); dftEmpty(&data); return(1);
//     }
//     if(n>1) {
//       n=dftSelectBestReference(&data);
//       if(n<0) {        /* no vois are matching */
//         fprintf(stderr, "Error: Cannot find ref voi '%s'.\n", dvrname);
//         dftEmpty(&input); dftEmpty(&data); return(1);
//       }
//       dvr_roi=n;
//       fprintf(stderr,"Warning: several ref regions match; %s is selected.\n",
//         data.voi[dvr_roi].name );
//     } else {
//       for(ri=0; ri<data.voiNr; ri++) if(data.voi[ri].sw) {dvr_roi=ri; break;}
//     }
//     if(verbose>0)
//       printf("Reference region: %s (voi=%d)\n", data.voi[dvr_roi].name, dvr_roi);
//   }



  /* Finishing reading data */
    // theta = (double *)malloc(frameNr*sizeof(double));
    // dv    = (double *)malloc(frameNr*sizeof(double));
    // ct    = (double *)malloc(frameNr*sizeof(double));
    // ict   = (double *)malloc(frameNr*sizeof(double));
    // ci    = (double *)malloc(frameNr*sizeof(double));
    // ici   = (double *)malloc(frameNr*sizeof(double));
    // t     = (double *)malloc(frameNr*sizeof(double));
    // // theta = (double *)malloc(frameNr*sizeof(double));

    // for (int i=0; i<data.frameNr; i++) { 
    //   ct[i]=data.voi[0].y[i];
    //   ict[i]=data.voi[ri].y2[i];
    //   ci[i] = input.voi[ri].y[i];
    //   ici[i]=input.voi[ri].y2[i];
    //   t[i]=input.x[i];
    //   theta[i] = data.voi[ri].y2[i];  
    //   dv[i] = data.voi[ri].y3[i];
    // }

    ct   = data.voi[0].y;
    ict  = data.voi[ri].y2;
    ci   = input.voi[ri].y;
    ici  = input.voi[ri].y2;
    t    = input.x;
    theta = data.voi[ri].y2;  
    dv   = data.voi[ri].y3;


    if (verbose>9) {
        // FILE *ppfile = fopen(debugfile, "a+");
        printf( "CT %f %f supplied\n", ct[0],ct[frameNr-1]);
        printf( "CT %f %f supplied\n", ict[0],ict[frameNr-1]);
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
      if(verbose>2) printf("calculating %s\n", data.voi[ri].name);

    /*
     *  Set axis weights;
     *  do this again for each region, because these can be set to zero below
     */
    for(fi=0; fi<data.frameNr; fi++) {
      if(data.isweight) wx[fi]=wy[fi]=data.w[fi]; else wx[fi]=wy[fi]=1.0; 
      if(data.timetype==DFT_TIME_STARTEND) 
        wx[fi]*=data.x2[fi]; else wx[fi]*=data.x[fi]; 
    }

    /* Set data pointers */
    ci=input.voi[0].y; ici=input.voi[0].y2;
    ct=data.voi[ri].y; ict=data.voi[ri].y2;
    theta=data.voi[ri].y2; /* Note that this is the same where integral is! */
    dv=data.voi[ri].y3;

    /* Calculate Logan plot data */
    for(fi=data.frameNr-1; fi>=0; fi--) if(ct[fi]!=0.0) {
      if(verbose>8) {
        printf("%03d %8.3f : ici=%g ci=%g ict=%g ct=%g\n",
          fi+1, data.x[fi], ici[fi], ci[fi], ict[fi], ct[fi] );
      }
      if (logan_mode==0) {
      /* dv (y axis) */
      dv[fi]=ict[fi]/ct[fi]; /*if(dv[fi]<0.0) wy[fi]=0;*/
      /* theta (x axis) */
      if(k2>0) theta[fi]=(ici[fi]+ci[fi]/k2)/ct[fi];
      else theta[fi]=ici[fi]/ct[fi];
      } else {
      /* dv (y axis) */
      dv[fi]=ict[fi]/ci[fi]; /*if(dv[fi]<0.0) wy[fi]=0;*/
      /* theta (x axis) */
      if(k2>0) theta[fi]=(ici[fi]+ci[fi]/k2)/ci[fi];
      else theta[fi]=ici[fi]/ci[fi];
      }

      /* check the close-to-zeroes in the first frames */
      if(data.x[fi]<0.1*data.x[data.frameNr-1]) { 
        if(theta[fi]>theta[data.frameNr-1] || dv[fi]>dv[data.frameNr-1]) {
          if(verbose>2)
            printf("Possible close-to-zero plot point at %g -> set to zero.\n",
                   data.x[fi]);
          theta[fi]=dv[fi]=wx[fi]=wy[fi]=0.0;
        }
      }
    } else theta[fi]=dv[fi]=wx[fi]=wy[fi]=0.0;

    if(verbose>6) {
      for(fi=first; fi<=last; fi++)
        printf("%03d %8.3f : %g %g  (%g %g)\n",
          fi+1, data.x[fi], theta[fi], dv[fi], wx[fi], wy[fi] );
    }

    /* Linear fit */    
    DVSD=DV=Ic=IcSD=SWSS=0.0; ret=0;
    if(llsq_model==1) {
      if(first==0) {
        /* Calculation of best LLSQ fit with errors in both coordinates */
        ret=best_logan_reed(
          &theta[first], &dv[first], &wx[first], &wy[first], dataNr, 5,
          &DV, &Ic, &SWSS, &DVSD, &IcSD, cx, cy, &fi, verbose-4
        );
        if(verbose>7) printf("Min NWSS with %d data points.\n", fi);
        res.voi[ri].sw=fi;
      } else {
        /* Calculation of LLSQ fit with errors in both coordinates */
        ret=llsqwt(
          /* input */
          &theta[first], &dv[first], dataNr, &wx[first], &wy[first],
          1.0e-10, /* allowed tolerance in slope estimation */
          /* input/output */
          &w[first], /* work vector; effective weights w[i] are returned in it */
          /* output ; SWSS is already normalized */
          &Ic, &DV, &SWSS, &IcSD, &DVSD, cx, cy
        );
      }
      if(verbose>6) {
        printf("%s:\n", data.voi[ri].name);
        for(fi=first; fi<=last; fi++)
          printf("%03d %8.3f : %g %g  (%g %g -> %g)\n",
            fi+1, data.x[fi], theta[fi], dv[fi], wx[fi], wy[fi], w[fi] );
      }
    } else if(llsq_model==2) {
      /* Check that theta is not negative */
      for(fi=first; fi<=last; fi++) if(theta[fi]<0.0) theta[fi]=nan("");
      /* Calculation of perpendicular line fit */
      ret=llsqperp3(
        /* input */
        &theta[first], &dv[first], dataNr,
        /* output ; SSD is already normalized */
        &DV, &Ic, &SWSS
      );
    } else if(llsq_model==0) {
      if(first==0) {
        /* Calculation of best regression line */
        ret=best_logan_regr(
          &theta[first], &dv[first], &wx[first], &wy[first], dataNr, 5,
          &DV, &Ic, &SWSS, &DVSD, &IcSD, &fi, verbose-4
        );
        res.voi[ri].sw=fi;
        if(verbose>9) printf("Dv=%g Ic=%g\n", DV, Ic);
      } else {
        /* Check that theta is not negative */
        for(fi=first; fi<=last; fi++) if(theta[fi]<0.0) theta[fi]=nan("");
        /* Calculation of linear regression using pearson() */
        ret=pearson3(
          /* input */
          &theta[first], &dv[first], dataNr,
          /* output */
          &DV, &DVSD, &Ic, &IcSD, &SWSS, &f
        );
      }
    } else if(llsq_model==3) {
      /* Check that theta is not negative */
      for(fi=first; fi<=last; fi++) if(theta[fi]<0.0) theta[fi]=nan("");
      /* Calculation of median slope and ic */
      ret=medianline(
        &theta[first], &dv[first], dataNr, &DV, &Ic
      );
    }
  //   if(ret==0) {
  //     res.voi[ri].parameter[0]=DV; if(save_stat) res.voi[ri].sd[0]=DVSD;
  //     if(inputtype==5 && dvr_minus_one!=0) {
  //       res.voi[ri].parameter[0]-=1.0;
  //     }
  //     res.voi[ri].parameter[1]=Ic; if(save_stat) res.voi[ri].sd[1]=IcSD;
  //     res.voi[ri].parameter[2]=SWSS;
  //     if(verbose>2) printf("DV := %g (%g)\n", DV, DVSD);
  //   } else
  //      fprintf(stderr, "Error (%d) in linear fit of %s\n",
  //                    ret, data.voi[ri].name);
  // // } /* next region */

  // /* Compute the DVR, BPnd or BPp if plasma input */
  // if(dvr_roi>=0) {
  //   for(ri=0; ri<data.voiNr; ri++) {
  //     if(bp_type==1)
  //       res.voi[ri].parameter[3]=
  //         res.voi[ri].parameter[0]/res.voi[dvr_roi].parameter[0];
  //     else if(bp_type==2)
  //       res.voi[ri].parameter[3]=
  //         res.voi[ri].parameter[0]/res.voi[dvr_roi].parameter[0] - 1.0;
  //     else
  //       res.voi[ri].parameter[3]=
  //         res.voi[ri].parameter[0]-res.voi[dvr_roi].parameter[0];
  //   }
  // }


  /* return results */
  output[0] = DV;
  output[1] = Ic;
  output[2] = DVSD;
  output[3] = IcSD;
  output[4] = SWSS;           //correlation coefficient


  // clean memory! dftEmpty(&temp);
  resEmpty(&res); dftEmpty(&input); dftEmpty(&data); 
  
//  fclose(pfile);
//  free(theta);
//  free(dv);
//   free(ci)  ;
//   free(ici)  ;
//   free(ct)  ;
//   free(t)  ;  
}


/******************************************************************************/
/// @endcond
/******************************************************************************/
/** Finds the best least-squares line to (x,y)-data, leaving points out
    from the beginning.
    @return Returns 0, if ok.
 */
int best_logan_reed(
  /** Plot x axis values */
  double *x,
  /** Plot y axis values */
  double *y,
  /** Weighting factors for x */
  double *wx,
  /** Weighting factors for y */
  double *wy,
  /** Nr of plot data points */
  int nr,
  /** Min nr of points to use in the fit; must be >=4 */
  int min_nr,
  /** Slope is returned in here */
  double *slope,
  /** Y axis intercept is returned in here */
  double *ic, 
  /** sqrt(WSS)/wsum is returned here */
  double *nwss,
  /** Expected sd of slope at calculated points */
  double *sslope,
  /** Expected sd of intercept at calculated points */
  double *sic,
  /** Calculated x data points */
  double *cx,
  /** Calculated y data points */
  double *cy,
  /** Number of points in the best fit (incl data with w=0) */
  int *bnr,
  /** Verbose level; if zero, then nothing is printed to stderr or stdout */
  int verbose
) {
  int from, to, ret, from_min, to_min;
  double *w, lic, lslope, lnwss, nwss_min;

  if(verbose>0) printf("best_logan_reed()\n");
  /* Check the data */
  if(x==NULL || y==NULL || nr<min_nr || nr<2) return(1);
  /* Check parameters */
  if(min_nr<4) return(2);

  /* Make room for work vector */
  w=(double*)malloc(nr*sizeof(double)); if(w==NULL) return(3);

  /* Search the plot range that gives the best nwss */
  nwss_min=9.99E+99; from_min=to_min=-1;
  for(from=0, to=nr-1; from<nr-min_nr; from++) {
    ret=llsqwt(x+from, y+from, (to-from)+1, wx+from, wy+from, 1.0E-10, w,
               &lic, &lslope, &lnwss, NULL, NULL, cx+from, cy+from);
    if(verbose>1) {
      printf("  range: %d-%d ; nwss=%g ; min=%g ; ret=%d\n",
        from, to, lnwss, nwss_min, ret);
    }
    if(ret==0 && lslope>0.0 && lnwss<nwss_min) {
      nwss_min=lnwss; from_min=from; to_min=to;}
  }
  if(from_min<0) {free(w); return(5);}

  /* Run llsqwt() again with that range, now with better resolution      */
  /* and this time compute also SD's                                     */
  from=from_min; to=to_min;
  ret=llsqwt(x+from, y+from, (to-from)+1, wx+from, wy+from, 1.0E-12, w,
             ic, slope, nwss, sic, sslope, cx+from, cy+from);
  free(w); if(ret) return(6);
  *bnr=(to-from)+1;

  return(0);
}
/******************************************************************************/

/******************************************************************************/
/** Find the best regression line to (x,y)-data, leaving points out
 *  from the beginning.
 *  @return Returns 0, if ok.
 */
int best_logan_regr(
  /** Plot x axis values */
  double *x,
  /** Plot y axis values */
  double *y,
  /** Weighting factors for x */
  double *wx,
  /** Weighting factors for y */
  double *wy,
  /** Nr of plot data points */
  int nr,
  /** Min nr of points to use in the fit; must be >=4 */
  int min_nr,
  /** Slope is returned in here */
  double *slope,
  /** Y axis intercept is returned in here */
  double *ic,
  /** Pearson's correlation coefficient is returned here */
  double *r,
  /** Expected sd of slope at calculated points */
  double *sslope,
  /** Expected sd of intercept at calculated points */
  double *sic,
  /** Number of points in the best fit (incl data with w=0) */
  int *bnr,
  /** Verbose level; if zero, then nothing is printed to stderr or stdout */
  int verbose
) {
  int fi, from, to, ret, from_min, to_min, n;
  double lic, lslope, lic_sd, lslope_sd, cv, cv_min, f;
  double *cx, *cy;

  if(verbose>0) printf("best_logan_regr()\n");
  /* Check the data */
  if(x==NULL || y==NULL || nr<min_nr || nr<2) return(1);
  /* Check parameters */
  if(min_nr<4) return(2);

  /* Construct a checked data with no negative values and weights>0 */
  cx=(double*)malloc(nr*sizeof(double));
  cy=(double*)malloc(nr*sizeof(double));
  if(cx==NULL || cy==NULL) return(3);
  for(fi=n=0; fi<nr; fi++)
    if(wx[fi]>0 && wy[fi]>0 && !isnan(x[fi]) && !isnan(y[fi])) {
      cx[n]=x[fi]; cy[n]=y[fi]; n++;}
  if(n<min_nr) {free(cx); free(cy); return(4);}

  /* Search the plot range that gives the lowest CV for slope */
  cv_min=9.99E+99; from_min=to_min=-1;
  for(from=0, to=n-1; from<n-min_nr; from++) {
    /* Calculation of linear regression using pearson() */
    ret=pearson(
      cx+from, cy+from, (to-from)+1,
      &lslope, &lslope_sd, &lic, &lic_sd, r, &f
    );
    if(ret==0 && lslope>0) {
      cv=lslope_sd/lslope;
    } else cv=9.99E+99;
    if(cv<cv_min) {
      cv_min=cv; from_min=from; to_min=to;}
  }
  if(from_min<0) {free(cx); free(cy); return(5);}

  /* Run pearson() again with that range */
  from=from_min; to=to_min;
  ret=pearson(
    cx+from, cy+from, (to-from)+1,
    slope, sslope, ic, sic, r, &f
  );
  free(cx); free(cy); if(ret) return(6);
  *bnr=(to-from)+1;

  return(0);
}
/******************************************************************************/
