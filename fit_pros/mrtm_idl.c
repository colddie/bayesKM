//mrtm_idl.c


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
#define NNLS_N 3
/*****************************************************************************/


/*****************************************************************************/
/**
 *  Main
 */
int mrtm_idl(int argc,  float * argv[])         //Pure c does not require extern  "C" 
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
  double     f, xm, ym, xs, ys, k2=-1.0, lambda=0;
  RES        res;
  int        bp_type; // 0=no, 1=DVR, 2=BPnd, 3=BPp
  int        dvr_minus_one; // 0=DVR reported, 1=BPnd reported with ref input

  double    *t, *theta, *dv, *ci, *ici, *ct, *ict;
  int        dataNr=0, first, last;
  double    *t0, *t1, *tac, *ctt, *output, *weights;
  int       voiNr = 1;
  unsigned int    frameNr, isweight = 0, logan_mode = 0, directbp=0;

  /* nnls */
  int      nnls_n, nnls_m, n, m, nnls_index[NNLS_N];
  double  *nnls_a[NNLS_N], *nnls_b, *nnls_zz, nnls_x[NNLS_N], *nnls_mat,
           nnls_wp[NNLS_N], *dptr, nnls_rnorm;

  const char *debugfile1 = "debug1.txt";
  // FILE *pfile = fopen(debugfile1, "a+");
  const char *debugfile = "debug.txt";

  dftInit(&data); dftInit(&input); resInit(&res);



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
  isweight =  *(unsigned int*) argv[9];
  weights = (double*) argv[10];
  directbp = *(unsigned int*) argv[11];


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
    if(data.isweight) data.w[i]=weights[i];
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
    // printf("\nIDL input data:\n");
    // dftPrint(&temp);
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
    // theta = (double *)malloc(frameNr*sizeof(double));

    // for (int i=0; i<data.frameNr; i++) { 
    //   ct[i]=data.voi[0].y[i];
    //   ict[i]=data.voi[ri].y2[i];
    //   ci[i] = input.voi[ri].y[i];
    //   ici[i]=input.voi[ri].y2[i];
    //   t[i]=input.x[i];
    //   theta[i] = data.voi[ri].y2[i];  
    //   dv[i] = data.voi[ri].y3[i];
    // }
    ct=data.voi[0].y;
    ict=data.voi[ri].y2;
    ci = input.voi[ri].y;
    ici=input.voi[ri].y2;
    t=input.x;
    theta = data.voi[ri].y2;  
    dv = data.voi[ri].y3;

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
   *  Calculate R1, k2',k2A for each region
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

  // /* Copy weights if available */
  // /* or set them to frame lengths */

  /*
   *  Allocate memory required by NNLS
   */
  if(verbose>1) printf("allocating memory for NNLS\n");
  nnls_n=NNLS_N; nnls_m=dataNr;
  nnls_mat=(double*)malloc(((nnls_n+2)*nnls_m)*sizeof(double));
  if(nnls_mat==NULL) {
    fprintf(stderr, "Error: cannot allocate memory for NNLS.\n");
    return(5);
  }
  for(n=0, dptr=nnls_mat; n<nnls_n; n++) {nnls_a[n]=dptr; dptr+=nnls_m;}
  nnls_b=dptr; dptr+=nnls_m; nnls_zz=dptr;

    /* Calculate the matrix for NNLS */
    /* Fill  A matrix: */
    /* function #1:  */
    for(m=0; m<nnls_m; m++) nnls_a[0][m]=ci[m];
    /* function #2:  */
    for(m=0; m<nnls_m; m++) nnls_a[1][m]=ici[m];
    /* function #3:  */
    for(m=0; m<nnls_m; m++) nnls_a[2][m]=-ict[m];
    /* Fill  B array:  */
    for(m=0; m<nnls_m; m++) nnls_b[m]=ct[m];


    /* Apply data weights */
    if(data.isweight) nnlsWght(nnls_n, nnls_m, nnls_a, nnls_b, data.w);
    if(verbose>6) {
      printf("Matrix A                     Array B\n");
      for(m=0; m<nnls_m; m++) {
        printf("%12.3f %12.3f %12.3f     %12.3f\n",
          nnls_a[0][m], nnls_a[1][m], nnls_a[2][m], nnls_b[m]);
      }
    }

    /* NNLS */
    ret=nnls(nnls_a, nnls_m, nnls_n, nnls_b, nnls_x, &nnls_rnorm,
              nnls_wp, nnls_zz, nnls_index);
    if(ret>1) { /* no solution is possible */
      printf('no solution available'); return(ret); // nosolution_nr++; continue;
    }
 
    for(n=0; n<nnls_n; n++) output[n]=nnls_x[n];
    // for(n=0; n<nnls_n; n++) tout.m[pi][yi][xi][n]=nnls_x[n];
    // /* Get R1 */
    // if(r1file[0]) r1out.m[pi][yi][xi][0]=nnls_x[0];
    // /* Get k2 */
    // if(k2file[0]) k2out.m[pi][yi][xi][0]=nnls_x[1];
    // /* Get theta3=p3+lambda */ //k2A, not used lambda here theta3 = k2/(1+BPnd)+lambda
    // if(t3file[0]) {
    //   if(nnls_wp[2]==0.0) t3out.m[pi][yi][xi][0]=nnls_x[2]+lambda;
    //   else t3out.m[pi][yi][xi][0]=0.0;
    // }
    // /* Get k2' (k2 of reference region), when BP>0 */
    // if(k2sfile[0] || model==MODEL_SRTM2) {
    //   if(nnls_x[2]>0.0) f=nnls_x[1]/nnls_x[2]; else f=0.0; // f=BP+1
    //   if(f>1.0 && nnls_x[0]>0.0 && nnls_x[1]>0.0)
    //     k2sout.m[pi][yi][xi][0]=nnls_x[1]/nnls_x[0];
    // }
  // }


 /* Estimate BP without division, unless using SRTM2 model */
 if (directbp) {
 /* Fill  A matrix: */
    /* function #1:  */
    for(m=0; m<nnls_m; m++) nnls_a[0][m]=ci[m];
    /* function #2:  */
    for(m=0; m<nnls_m; m++) nnls_a[1][m]=ici[m];
    /* function #3:  */
    for(m=0; m<nnls_m; m++) nnls_a[2][m]=-ct[m];
    /* Fill  B array:  */
    for(m=0; m<nnls_m; m++) nnls_b[m]=ict[m];


    /* Apply data weights */
    if(data.isweight) nnlsWght(nnls_n, nnls_m, nnls_a, nnls_b, data.w);
    if(verbose>6) {
      printf("Matrix A                     Array B\n");
      for(m=0; m<nnls_m; m++) {
        printf("%12.3f %12.3f %12.3f     %12.3f\n",
          nnls_a[0][m], nnls_a[1][m], nnls_a[2][m], nnls_b[m]);
      }
    }

    /* NNLS */
    ret=nnls(nnls_a, nnls_m, nnls_n, nnls_b, nnls_x, &nnls_rnorm,
              nnls_wp, nnls_zz, nnls_index);
    if(ret>1) { /* no solution is possible */
      printf('no solution available'); return(ret); // nosolution_nr++; continue;
    }
 
    for(n=0; n<nnls_n; n++) output[n]=nnls_x[n]-1;   // BP=DVR-1
    /* Get BP */
    // bpout.m[pi][yi][xi][0]=nnls_x[1];
 }
















//  fclose(pfile);
//  free(theta);
//  free(dv);
//   free(ci)  ;
//   free(ici)  ;
//   free(ct)  ;
//   free(t)  ;  
}
