/******************************************************************************
 * This file is not compiled into the library, but it contains main()
 * which is compiled to an executable, used to test the library functions. 
 *****************************************************************************/

/*****************************************************************************/
#include "libtpcmodel.h"
#include "libtpcmisc.h"
/*****************************************************************************/

/*****************************************************************************/
/* Global variables and array pointers for certain objective functions */
int fitframeNr;
double *simdata, *measdata, *pmin, *pmax, *p, *w;
/*****************************************************************************/
#define MAXVAL 100000;
/*****************************************************************************/
// Test function declarations:
int test_re(int VERBOSE);
int test_runst(int VERBOSE);
int test_normaldistr(int VERBOSE);
int test_integr(int VERBOSE);
int test_tgoRandomParametersST(int VERBOSE);
int test_powell(int VERBOSE);
int test_tgo(int VERBOSE);
int test_bobyqa1(int VERBOSE);
int test_constraints1(int VERBOSE);
int test_scales1(int VERBOSE);
int test_scales2(int VERBOSE);
int test_onedim1(int VERBOSE);
int test_onedim2(int VERBOSE);
int test_banana1(int VERBOSE);
int test_rastrigin(int VERBOSE);
int test_nptrange(int VERBOSE);
int test_bootstrap1(int VERBOSE);
double bobyqa_problem1(int n, double *x, void *func_data);
double bobyqa_problem2(int n, double *x, void *func_data);
double optfunc_dejong2(int n, double *x, void *func_data);
double optfunc_rastrigin(int n, double *x, void *func_data);
double func_deviation(int parNr, double *p, void *fdata);
/*****************************************************************************/

/*****************************************************************************/
static char *info[] = {
  "Usage: @P [options]",
  " ",
  "Options:",
  " -stdoptions", // List standard options like --help, -v, etc
  " -t, --test",
  "     Run all tests for library functions.",
  0};
/*****************************************************************************/

/*****************************************************************************/
/** Run unit tests to the library functions
 *  @author Vesa Oikonen
 *  @return 0 if all tests pass, otherwise >0.
 * */
int main(
  /** Nr of arguments */
  int argc,
  /** Pointer to arrays of argument string */
  char *argv[ ]
) {
  int i, help=0, version=0, verbose=1, error=0, test=0;
  int ret;
  char *cptr;

  if(argc==1) {tpcPrintUsage(argv[0], info, stdout); return(0);}
  for(i=1; i<argc; i++) {
    if(tpcProcessStdOptions(argv[i], &help, &version, &verbose)==0) continue;
    cptr=argv[i]; if(*cptr=='-') cptr++; if(*cptr=='-') cptr++;
    if(strncasecmp(cptr, "TEST", 1)==0) {
      test=1; continue;
    } else {
      error++; break;
    }
  }
  if(error>0) {
    fprintf(stderr, "Error: specify --help for usage.\n");
    return(1);
  }
  /* Print help or version? */
  if(help) {tpcPrintUsage(argv[0], info, stdout); return(0);}
  if(version) {tpcPrintBuild(argv[0], stdout); return(0);}

  if(test==0) return(0);

  if(verbose>0) printf("running tests for library functions...\n");
  drandSeed(1); //srand(time(0));
  i=10;
  i++; if((ret=test_re(verbose-1))!=0) {
    fprintf(stderr, "failed (%d).\n", ret); return(i);}
  i++; if((ret=test_runst(verbose-1))!=0) {
    fprintf(stderr, "failed (%d).\n", ret); return(i);}
  i++; if((ret=test_normaldistr(verbose-1))!=0) {
    fprintf(stderr, "failed (%d).\n", ret); return(i);}
  i++; if((ret=test_integr(verbose-1))!=0) {
    fprintf(stderr, "failed (%d).\n", ret); return(i);}
  i++; if((ret=test_tgoRandomParametersST(verbose-1))!=0) {
    fprintf(stderr, "failed (%d).\n", ret); return(i);}
  i++; if((ret=test_powell(verbose-1))!=0) {
    fprintf(stderr, "failed (%d).\n", ret); return(i);}
  i++; if((ret=test_bootstrap1(verbose-1))!=0) {
    fprintf(stderr, "failed (%d).\n", ret); return(i);}
  i++; if((ret=test_tgo(verbose-1))!=0) {
    fprintf(stderr, "failed (%d).\n", ret); return(i);}

  /* Bobyqa */
  i++; if((ret=test_bobyqa1(verbose-1))!=0) {
    fprintf(stderr, "failed (%d).\n", ret); return(i);}
  i++; if((ret=test_constraints1(verbose-1))!=0) {
    fprintf(stderr, "failed (%d).\n", ret); return(i);}
  i++; if((ret=test_scales1(verbose-1))!=0) {
    fprintf(stderr, "failed (%d).\n", ret); return(i);}
  i++; if((ret=test_scales2(verbose-1))!=0) {
    fprintf(stderr, "failed (%d).\n", ret); return(i);}
  i++; if((ret=test_onedim1(verbose-1))!=0) {
    fprintf(stderr, "failed (%d).\n", ret); return(i);}
  i++; if((ret=test_onedim2(verbose-1))!=0) {
    fprintf(stderr, "failed (%d).\n", ret); return(i);}
  i++; if((ret=test_banana1(verbose-1))!=0) {
    fprintf(stderr, "failed (%d).\n", ret); return(i);}
  i++; if((ret=test_rastrigin(verbose-1))!=0) {
    fprintf(stderr, "failed (%d).\n", ret); return(i);}
  i++; if((ret=test_nptrange(verbose-1))!=0) {
    fprintf(stderr, "failed (%d).\n", ret); return(i);}


  if(verbose>0) printf("\nAll tests passed.\n\n");
  return(0);
}
/*****************************************************************************/

/*****************************************************************************/
int test_bootstrap1(int VERBOSE)
{
  printf("test_bootstrap1()\n");
  if(VERBOSE) printf("\nOne parameter, data with Gaussian noise, no weights.\n");

  int i, j, ret;
  const int dataNr=50, parNr=1;
  const double SD=100.0;
  double observedData[dataNr], fittedData[dataNr], mean, meansd;
  double local_p[parNr], local_pmin[parNr], local_pmax[parNr];
  double sd[parNr], cl1[parNr], cl2[parNr];
  double local_w[dataNr];
  char temp[128];
  double tmpdata1[dataNr], tmpdata2[dataNr];
  simdata=tmpdata1; // used by objective function to simulate data in
  measdata=tmpdata2; // used by objective function to get measured data
                     // (in this case, bootstrapped data)
  const int repeats=200;
  double repmeans[repeats], bs_sdlist[repeats];
  
  drandSeed(1);
  //GAUSSDEV_SEED=(long int)time(NULL); srand(GAUSSDEV_SEED);

  /* Set global pointers and variables to be used by objective function */
  p=local_p; pmin=local_pmin; pmax=local_pmax; w=local_w;
  fitframeNr=dataNr;

  for(j=0; j<repeats; j++) {

    /* Create dataset */
    for(i=0; i<dataNr; i++) {
      observedData[i]=1000.0;
      local_w[i]=1.0;
    }
    /* Add noise */
    for(i=0; i<dataNr; i++) observedData[i] += SD*gaussdev2();
    /* Compute mean and sd from noisy data */
    mean=dmean(observedData, dataNr, &meansd);
    if(VERBOSE>1) printf("  simulated mean=%g and sd=%g\n", mean, meansd);
    /* 'fitted' data is the mean of noisy data */
    for(i=0; i<dataNr; i++) fittedData[i]=mean;
    /* 'fitted' parameter value is the mean */
    for(i=0; i<parNr; i++) local_p[i]=mean;
    repmeans[j]=mean;

    /* Set parameter limits */
    for(i=0; i<parNr; i++) {local_pmin[i]=0.0; local_pmax[i]=2000.0;}

    /* Try bootstrapping */
    // double matrix[300*parNr]; 
    ret=bootstrap(
           300, cl1, cl2, sd, p, pmin, pmax,
           dataNr, observedData, fittedData,
           measdata, parNr, w, func_deviation,
           temp, VERBOSE-2);
    if(ret!=0) {
      printf("Error %d in bootstrap() function: %s\n", ret, temp);
      return 11;
    }
    if(VERBOSE>1) {
      printf("  sd := %g\n", sd[0]);
      printf("  CL95%% := %g -%g\n", cl1[0], cl2[0]);
    }
    bs_sdlist[j]=sd[0];
  }

  /* Calculate the SD of means from repeats */
  mean=dmean(repmeans, repeats, &meansd);
  if(VERBOSE) printf("  simulated sd=%g\n", meansd);
  /* Calculate the mean of bootstrapped SDs from repeats */
  mean=dmean(bs_sdlist, repeats, NULL);
  if(VERBOSE) printf("  bootstrapped sd=%g\n", mean);
  /* Compare that they match */
  if(fabs(2.0*(meansd-mean)/(meansd+mean))>0.15) {
    printf("Error: SD from bootstrap() is too far from true SD.\n");
    if(VERBOSE==0) {
      printf("  simulated sd=%g\n", meansd);
      printf("  bootstrapped sd=%g\n", mean);
    }
    return 21;
  }

  if(VERBOSE) printf("   SUCCEEDED\n");
  return 0;
}
/******************************************************************************/

/******************************************************************************/
/* Check robust functions */
int test_re(int VERBOSE) {

  double data[100], est, var;
  int dataNr, n;
  int error_code = 0;

  printf("test_re()\n");
  /* For a single data point 0, all functions should give 0 as result */

  data[0]=data[1]=0.;
  dataNr=2;

  est=dmedian(data, dataNr);
  if(est<0.0001 && est>-0.0001){
    if(VERBOSE)
      printf("   SUCCESFULL: dmedian() operation succeeded for one data point.\n");
  } else{
    if(VERBOSE){printf("   FAILED: dmedian() failed for one data point.\n");}
    return (1);
  }

  est=mEstim(data, dataNr, 10, 0.5);
  if(est<0.0001 && est>-0.0001){
    if(VERBOSE)
      printf("   SUCCESFULL: mEstim() operation succeeded for one data point.\n");
  } else{
    if(VERBOSE){printf("   FAILED: mEstim() failed for one data point.\n");}
    return (1);
  }

  est=least_median_of_squares(data, dataNr);
  if(est<0.0001 && est>-0.0001){
    if(VERBOSE) {
      printf("   SUCCESFULL: least_median_of_squares() operation succeeded");
      printf(" for one data point.\n");
    }
  } else {
    if(VERBOSE)
      printf("   FAILED: least_median_of_squares() failed for one data point.\n");
    return (1);
  }


  error_code=least_trimmed_square(data, dataNr, &est, &var);
  if(error_code){
    if(est<0.0001 && est>-0.0001) {
      if(VERBOSE){
        printf("   SUCCESFULL: least_trimmed_square() succeeded with");
        printf(" error code: %i\n",error_code);
      }
    } else if(VERBOSE)
      printf("   FAILED: least_trimmed_square() failed with error code: %i\n",
             error_code);
    return(error_code);
  }
  else{
    if(est<0.0001 && est>-0.0001){
      if(VERBOSE){
        printf("   SUCCESFULL: least_trimmed_square() operation succeeded");
        printf(" for one data point.\n");
      }
    } else{
      if(VERBOSE)
        printf("   FAILED: least_trimmed_square() failed for one data point.\n");
      return (1);
    }
  }

  /* For a uniform distribution 1,2,3,4,5,6,7,8,9,10 all functions should
     give 5.5 as result */
  for(n=0; n<10; n++) data[n]=n+1;
  dataNr=10;

  est=dmedian(data, dataNr);
  if(est<5.50001 && est>5.49999){
    if(VERBOSE){
      printf("   SUCCESFULL: dmedian() operation succeeded for");
      printf(" uniform distribution.\n");
    }
  } else {
    if(VERBOSE){
       printf("   FAILED: dmedian() failed for uniform distribution.\n");}
    return (2);
  }

  est=mEstim(data, dataNr, 10, 0.5);
  if(est<5.50001 && est>5.49999) {
    if(VERBOSE){
      printf("   SUCCESFULL: mEstim() operation succeeded for");
      printf(" uniform distribution.\n");
    }
  } else {
    if(VERBOSE){
      printf("   FAILED: mEstim() failed for uniform distribution.\n"); }
    return (2);
  }


  /* LMS and LTS methods aren't able to handle uniform distributions */

  /* For skewed distribution 2.1,3.1,3.3,3.5,3.6,3.7,4.5,5.2,6.0,7.4
     median is 3.6, lms 3.8 */
  data[0]=2.1; data[1]=3.1; data[2]=3.3; data[3]=3.5; data[4]=3.6;
  data[5]=3.7; data[6]=4.5; data[7]=5.2; data[8]=6.0; data[9]=7.4;

  /* median is the average of data points 3.6 and 3.7*/

  est=dmedian(data, dataNr);
  if(est<3.6501 && est>3.6499){
    if(VERBOSE) {
      printf("   SUCCESFULL: dmedian() operation succeeded for");
      printf(" skewed distribution.\n");
    }
  } else{
    if(VERBOSE){
      printf("   FAILED: dmedian() failed for skewed distribution.\n");}
    return (3);
  }

  /* LMS finds the "shortest half" 3.1-4.5 and the result is midpoint of
     these values (3.8).*/
  est=least_median_of_squares(data, dataNr);
  if(est<3.8001 && est>3.7999) {
    if(VERBOSE) {
      printf("   SUCCESFULL: least_median_of_squares() operation succeeded for");
      printf(" skewed distribution.\n");
    }
  } else {
    if(VERBOSE){
      printf("   FAILED: least_median_of_squares() failed for");
      printf(" skewed distribution.\n");
    }
    return (3);
  }


  /* For skewed distribution 1, 1, 1, 1, 10
     hubers M-estimator is about 1.243, lts estimate is 1*/
  data[0]=1.; data[1]=1.; data[2]=1.; data[3]=1.; data[4]=10.;
  dataNr=5;

  /*Huber's M-estimator with one iteration and cutoffpoint 0.5 should
    give approximately 1.243 for this data according to algorithm described in
    Nonlinear_signal_prosessing_lecture_notes in
    www.cs.tut.fi/~eeroh/nonlin.html*/

  est=mEstim(data, dataNr, 1, 0.5);
  if(est<1.244 && est>1.242) {
    if(VERBOSE)
      printf("   SUCCESFULL: mEstim() succeeded for skewed distribution.\n");
  } else {
    if(VERBOSE){
      printf("   FAILED: mEstim() failed for skewed distribution.\n");}
    return (3);
  }

  /* LTS should maybe have more sophisticated test data.
    The results could also be compared to reqular least square*/

  error_code=least_trimmed_square(data, dataNr, &est, &var);
  if(error_code) {
    if(est<1.0001 && est>0.9999) {
      if(VERBOSE) {
        printf("   SUCCESFULL: least_trimmed_square() succeeded with");
        printf(" error code: %i\n", error_code);
      }
    } else if(VERBOSE) {
      printf("   FAILED: least_trimmed_square() failed with error code: %i\n",
             error_code);
    }
    return(error_code);
  } else {
    if(est<1.0001 && est>0.9999) {
      if(VERBOSE){
        printf("   SUCCESFULL: least_trimmed_square() operation succeeded");
        printf(" for skewed distribution.\n");
      }
    } else{
      if(VERBOSE){
        printf("   FAILED: least_trimmed_square() failed");
        printf(" for skewed distribution.\n");
      }
      return (3);
    }
  }
  return (error_code);
}
/******************************************************************************/

/******************************************************************************/
/** Function for testing runs_test
  \return
*/
int test_runst(int VERBOSE)
{
  double r1[12], r2[12];
  int Nr=12, level=5, runs=1, neg=0, pos=0, ret;

  printf("test_runst()\n");
  /*Test: independent residuals */

  r1[0]=8.0;
  r1[1]=50.3;
  r1[2]=162.4;
  r1[3]=379.4;
  r1[4]=225.9;
  r1[5]=100.2;
  r1[6]=87.4;
  r1[7]=89.2;
  r1[8]=85.6;
  r1[9]=73.1;
  r1[10]=61.2;
  r1[11]=61.3;

  r2[0]=6.03;
  r2[1]=37.7;
  r2[2]=140.1;
  r2[3]=311.0;
  r2[4]=192.6;
  r2[5]=98.7;
  r2[6]=86.2;
  r2[7]=91.4;
  r2[8]=85.3;
  r2[9]=76.6;
  r2[10]=58.9;
  r2[11]=57.1;

  /*test residuals function */
  
  residuals(r1, r2, Nr, &runs, &neg, &pos);


  if(runs!=5){
    if(VERBOSE){
      printf("   FAILED: residuals() failed to calculate nr of runs.\n");}
    return (2);
  }

  if(neg!=2){
    if(VERBOSE)
      printf("   FAILED: residuals() failed to calculate negative residuals.\n");
    return (2);
  }
  if(pos!=10){
    if(VERBOSE)
      printf("   FAILED: residuals() failed to calculate positive residuals.\n");
    return (3);
  }

  /*test runs_test function */

  ret=runs_test(r1, r2, Nr, level, NULL);
  if(ret>0){
    if(VERBOSE){
      printf("   FAILED: runs_test() failed with error code %d.\n", ret);}
    return (1);
  }
  if(ret==-1){
    if(VERBOSE){
      printf("   FAILED: runs_test() failed for independent residuals.\n");}
    return (1);
  }
  if(ret==0){
    if(VERBOSE)
      printf("   SUCCEEDED: runs_test() succeeded for independent residuals.\n");
  }


  /*Test: dependent residuals */

  r1[0]=8.0;
  r1[1]=50.3;
  r1[2]=162.4;
  r1[3]=379.4;
  r1[4]=225.9;
  r1[5]=100.2;
  r1[6]=85.1;
  r1[7]=83.2;
  r1[8]=85.0;
  r1[9]=73.1;
  r1[10]=58.3;
  r1[11]=57.0;

  r2[0]=6.03;
  r2[1]=37.7;
  r2[2]=140.1;
  r2[3]=311.0;
  r2[4]=192.6;
  r2[5]=98.7;
  r2[6]=86.2;
  r2[7]=84.4;
  r2[8]=85.3;
  r2[9]=76.6;
  r2[10]=58.9;
  r2[11]=57.1;

  /*test residuals function */
  
  residuals(r1, r2, Nr, &runs, &neg, &pos);

  if(runs!=2){
    if(VERBOSE){
      printf("   FAILED: residuals() failed to calculate nr of runs.\n");}
    return (2);
  }

  if(neg!=6){
    if(VERBOSE)
      printf("   FAILED: residuals() failed to calculate negative residuals.\n");
    return (2);
  }
  if(pos!=6){
    if(VERBOSE)
      printf("   FAILED: residuals() failed to calculate postive residuals.\n");
    return (3);
  }

  /*test runs_test function */

  ret=runs_test(r1, r2, Nr, level, NULL);
  if(ret>0){
    if(VERBOSE){
      printf("   FAILED: runs_test() failed with error code %d.\n", ret);}
    return (1);
  }
  if(ret==0){
    if(VERBOSE){
      printf("   FAILED: runs_test() failed for dependent residuals.\n");}
    return (1);
 
  }
  if(ret==-1){
    if(VERBOSE){
      printf("   SUCCEEDED: runs_test() succeeded for dependent residuals.\n");}
  }
  return (0);
}
/******************************************************************************/

/******************************************************************************/
int test_normaldistr(int VERBOSE)
{
  double x, ret;
  
  printf("test_normaldistr()\n");
  /*Test x=1.55 */
  x=1.55;
  ret=normal_pvalue_2(x);
  if(ret<0.121 || ret>0.122){
    if(VERBOSE){printf("   FAILED: normal_pvalue_2() failed for x=1.55.\n");}
    return (1);
  }
  ret=normal_pvalue_1(x);
  if(ret<0.060 || ret>0.061){
    if(VERBOSE){printf("   FAILED: normal_pvalue_1() failed for x=1.55.\n");}
    return (1);
  }

  /*Test x=0.5 */
  x=0.5;
  ret=normal_pvalue_2(x);
  if(ret<0.616 || ret>0.618){
    if(VERBOSE){printf("   FAILED: normal_pvalue_2() failed for x=0.5.\n");}
    return (1);
  }
  ret=normal_pvalue_1(x);
  if(ret<0.3084 || ret>0.3086){
    if(VERBOSE){printf("   FAILED: normal_pvalue_1() failed for x=0.5.\n");}
    return (1);
  }

  if(VERBOSE){
    printf("   SUCCEEDED: normal_pvalue_1() and normal_pvalue2() passed.\n");}

  return (0);
}
/******************************************************************************/

/******************************************************************************/
#if(0)
int test_polevl(int VERBOSE)
{
  double c2[2], c5[5], x, ret;

  printf("test_polevl()\n");
  /*test polynomial 1+1*x */

  c2[0]=c2[1]=1;
  x=1;

  ret=polevl(x, c2, 1);
  if(ret<1.99999 || ret>2.000001){
    if(VERBOSE)
      printf("   FAILED: polevl() failed for polynomial 1+1*x.\n");
    return (1);
  }

  /*test polynomial 1 + 1*0 + 1*0^2 + 1*0^3 + 1*0^4 */

  c5[0]=c5[1]=c5[2]=c5[3]=c5[4]=1.0;
  x=0;

  ret=polevl(x, c5, 4);

  if(ret<0.99999 || ret>1.00001){
    if(VERBOSE) 
      printf("   FAILED: polevl() failed for polynomial 1+1*0+1*0^2+...\n");
    return (2);
  }

  /*test polynomial 5+ 4*x + 3*x^2 + 2*x^3 + x^4 */
  c5[0]=1.0; c5[1]=2.0; c5[2]=3.0; c5[3]=4.0; c5[4]=5.0;
  x=2;

  ret=polevl(x, c5, 4);
  if(ret<56.99999 || ret>57.000001) {
    if(VERBOSE) {
      printf("   FAILED: polevl() failed for polynomial ");
      printf("5+ 4*x + 3*x^2 + 2*x^3 + x^4 .\n");
    }
    return (3);
  }

  if(VERBOSE){printf("   SUCCEEDED: polevl() passed.\n");}

  return (0);

}
#endif
/******************************************************************************/

/******************************************************************************/
#if(0)
int test_expx2(int VERBOSE)
{
  double x, ret;

  printf("test_expx2()\n");
  /*test exp(0*0)*/
  x=0;
  ret=expx2(x, 1);
  if(ret!=1.0){
    if(VERBOSE){printf("   FAILED: expx2() failed for exp(0*0).\n");}
    return (1);
  }

  /*test exp(-1*1)*/
  x=1;
  ret=expx2(x, -1);
  if(ret<0.3678 || ret>0.3679){
    if(VERBOSE){printf("   FAILED: expx2() failed for exp(-1*1).\n");}
    return (2);
  }

  /*test exp(0.25*0.25)*/
  x=0.25;
  ret=expx2(x, 1);
  if(ret<1.0644 || ret>1.0645){
    if(VERBOSE){printf("   FAILED: expx2() failed for exp(0.25*0.25).\n");}
    return (2);
  }
  if(VERBOSE){printf("   SUCCEEDED: expx2() passed.\n");}

  return (0);
}
#endif
/******************************************************************************/

/******************************************************************************/
int test_integr(int VERBOSE)
{
  double xin[10], x2in[10], yin[10], xnew1[10], yout[10];
  float f_xin[10], f_x2in[10], f_yin[10], f_xnew1[10], f_yout[10];
  int ret, i, nrin, nrout;

  printf("test_integr()\n");
  /* DATA 1:Five timepoints (1, 5), values 1.0, 2.0,...5.0
           new frames (1.5, 2.5,.., 5.5)*/

  for(i=0; i<5; i++){
     xin[i]=f_xin[i]=i+1;
     /*     x2in[i]=f_x2in[i]=i+2;*/
     yin[i]=f_yin[i]=i+1.0;
     xnew1[i]=f_xnew1[i]=i+1.5;
  }
  nrin=5;
  nrout=5;
  
  /* Function interpolate() and finterpolate*/

  ret=interpolate(xin, yin, nrin, xnew1, yout, NULL, NULL, nrout);
  if(ret!=0){
    if(VERBOSE){printf("\n   Test FAILED: interpolate() failed to execute.\n");}
    return (2);
  }
  /*  printf("%f %f %f %f %f\n", yout[0], yout[1], yout[2], yout[3], yout[4]);*/
  if(yout[0]!=1.5 || yout[1]!=2.5 || yout[2]!=3.5 || yout[3]!=4.5
     || yout[4]!=5.0)
  {
    if(VERBOSE){printf("\n   Test FAILED: interpolate() failed.\n");}
    return (2);
  }

  ret=finterpolate(f_xin, f_yin, nrin, f_xnew1, f_yout, NULL, NULL, nrout);
  if(ret!=0){
    if(VERBOSE){printf("   1 FAILED: f_interpolate() failed to execute.\n");}
    return (2);
  }
  if(f_yout[0]!=1.5 || f_yout[1]!=2.5 || f_yout[2]!=3.5 || f_yout[3]!=4.5 ||
     f_yout[4]!=5.0){
      if(VERBOSE){printf("   2 FAILED: f_interpolate() failed.\n");}
      return (2);
  }
  /* DATA 1.2:Three timepoints (3,5), values 3.0, 4.0, 5.0
             New frames (1.5,3.5)*/

  for(i=0; i<3; i++){
     xin[i]=f_xin[i]=i+3;
     /*     x2in[i]=f_x2in[i]=i+2;*/
     yin[i]=f_yin[i]=i+3;
 }
  xnew1[0]=f_xnew1[0]=1.5;
  xnew1[1]=f_xnew1[1]=3.5;
  nrin=3;
  nrout=2;

  /* Functions interpolate() and finterpolate*/

  ret=interpolate(xin, yin, nrin, xnew1, yout, NULL, NULL, nrout);
  if(ret!=0){
    if(VERBOSE){printf("   FAILED: interpolate() failed to execute.\n");}
    return (2);
  }
  if(yout[0]!=0.0 || yout[1]!=3.5 ){
      if(VERBOSE){printf("\n   Test FAILED: interpolate() failed.\n");}
      return (2);
  }
  /*  printf("%f %f %f %f %f %f\n",
      f_xin[0], f_xin[1], f_xin[2], f_yin[0], f_yin[1], f_yin[2]);*/

  ret=finterpolate(f_xin, f_yin, nrin, f_xnew1, f_yout, NULL, NULL, nrout);
  if(ret!=0){
    if(VERBOSE)
      printf("\n  3  Test FAILED: f_interpolate() failed to execute.\n");
    return (2);
  }
  if(f_yout[0]!=0.0 || f_yout[1]!=3.5 ){
      if(VERBOSE){printf("\n 4   Test FAILED: f_interpolate() failed.\n");}
      return (2);
  }

  /* DATA 2.1: Ten timepoints (1,10), all of value 1.0 */

  for(i=0; i<10; i++){
     xin[i]=f_xin[i]=i+1;
     x2in[i]=f_x2in[i]=i+2;
     yin[i]=f_yin[i]=1.0;
  }
  nrin=10;

  /* Functions integrate() and fintegrate */
  ret=integrate(xin, yin, nrin, yout);
  if(ret!=0){     
    if(VERBOSE){printf("   FAILED: integrate() failed to execute.\n");}
    return (2);
  }
  if(yout[9]!=9.5){
    if(VERBOSE){printf("   FAILED: integrate() failed.\n");}
    return (2);
  }

  ret=fintegrate(f_xin, f_yin, nrin, f_yout);
  if(ret!=0){     
    if(VERBOSE){printf("   FAILED: fintegrate() failed to execute.\n");}
    return (2);
  }
  if(yout[9]!=9.5){
    if(VERBOSE){printf("   FAILED: fintegrate() failed.\n");}
    return (2);
  }

  /* Functions petintegrate() and fpetintegrate */
  
  ret=petintegrate(xin, x2in, yin, nrin, yout, NULL);
  /*for(i=0;i<10;i++)printf("petintegrate %f\n", yout[i]);*/
  if(ret!=0){     
    if(VERBOSE){printf("   FAILED: petintegrate() failed to execute.\n");}
    return (2);
  }
  if(yout[9]>10.3334 || yout[9]<10.3332){
    if(VERBOSE){printf("   FAILED: petintegrate() failed.\n");}
    return (2);
  }

  ret=fpetintegrate(f_xin, f_x2in, f_yin, nrin, f_yout, NULL);
  if(ret!=0){     
    if(VERBOSE){printf("   FAILED: fpetintegrate() failed to execute.\n");}
    return (2);
  }
  if(f_yout[9]>10.3334 || f_yout[9]<10.3332){
    if(VERBOSE){printf("   FAILED: fpetintegrate() failed.\n");}
    return (2);
  }


  /* DATA 2.2: Ten timepoints (3,13), all of value 1.0 */

  for(i=0; i<10; i++){
     xin[i]=f_xin[i]=i+3;
     x2in[i]=f_x2in[i]=i+4;
     yin[i]=f_yin[i]=1.0;
  }
  nrin=10;

 /* Functions integrate() and fintegrate */
  ret=integrate(xin, yin, nrin, yout);
  if(ret!=0){     
    if(VERBOSE){printf("   FAILED: integrate() failed to execute.\n");}
    return (2);
  }
  if(yout[9]!=9.0){
    if(VERBOSE){printf("   FAILED: integrate() failed.\n");}
    return (2);
  }

  ret=fintegrate(f_xin, f_yin, nrin, f_yout);
  if(ret!=0){     
    if(VERBOSE){printf("   FAILED: fintegrate() failed to execute.\n");}
    return (2);
  }
  if(yout[9]!=9.0){
    if(VERBOSE){printf("   FAILED: fintegrate() failed.\n");}
    return (2);
  }

  /* Functions petintegrate() and fpetintegrate */
  ret=petintegrate(xin, x2in, yin, nrin, yout, NULL);
  if(ret!=0){     
    if(VERBOSE){printf("   FAILED: petintegrate() failed to execute.\n");}
    return (2);
  }
  if(yout[9]!=10.0){
    if(VERBOSE){printf("   FAILED: petintegrate() failed.\n");}
    return (2);
  }

  ret=fpetintegrate(f_xin, f_x2in, f_yin, nrin, f_yout, NULL);
  if(ret!=0){     
    if(VERBOSE){printf("   FAILED: fpetintegrate() failed execute.\n");}
    return (2);
  }
  if(yout[9]!=10.0){
    if(VERBOSE){printf("   FAILED: fpetintegrate() failed.\n");}
    return (2);
    }

  if(VERBOSE){printf("   SUCCEEDED: functions in integr.c passed.\n");}

  return(0);

}
/******************************************************************************/

/******************************************************************************/
int test_tgoRandomParametersST(int VERBOSE)
{
  int i, j;
  const int parNr=5, sampleNr=10000;
  double low[parNr], up[parNr], avg[parNr], medn[parNr],
         parmin[parNr], parmax[parNr];
  double parlist[sampleNr];
  TGO_POINT *points;

  printf("test_tgoRandomParametersST()\n");
  if(VERBOSE) printf("  sampleNr: %d\n", sampleNr);
  /* Allocate memory */
  points=(TGO_POINT*)malloc(sampleNr*sizeof(TGO_POINT));

  /* Set parameter limits */
  i=0; low[i]=0.0; up[i]=1000.0;
  i=1; low[i]=-10.0; up[i]=1000.0;
  i=2; low[i]=-1000.0; up[i]=10.0;
  i=3; low[i]=800.0; up[i]=1000.0;
  i=4; low[i]=1000.0; up[i]=1000.0;
  for(j=0; j<sampleNr; j++) points[j].topomin=0;
  tgoRandomParametersST(points, parNr, sampleNr, low, up);
  for(i=0; i<parNr; i++) {
    if(VERBOSE) printf("  Parameter %d:\n", i+1);
    for(j=0; j<sampleNr; j++) parlist[j]=points[j].par[i];
    avg[i]=dmean(parlist, sampleNr, NULL);
    medn[i]=dmedian(parlist, sampleNr);
    parmin[i]=parmax[i]=parlist[0];
    for(j=1; j<sampleNr; j++)
      if(parlist[j]<parmin[i]) parmin[i]=parlist[j];
      else if(parlist[j]>parmax[i]) parmax[i]=parlist[j];
    if(VERBOSE) {
      printf("    limits: [%g,%g]\n", low[i], up[i]);
      printf("    mean: %g\n", avg[i]);
      printf("    median: %g\n", medn[i]);
      printf("    range: %g - %g\n", parmin[i], parmax[i]);
    }
    if(parmin[i]<low[i] || parmax[i]>up[i]) {
      if(VERBOSE) printf("   FAILED with limits.\n\n");
      free(points); return 2;
    }
    if(low[i]<up[i] && fabs(avg[i])<fabs(medn[i])) {
      if( fabs(avg[i]-medn[i]) > 0.005*(up[i]-low[i]) ) {
        if(VERBOSE) printf("   FAILED with bias.\n\n");
        free(points); return 3;
      } else {
        if(VERBOSE) printf("   Warning of bias.\n\n");
      }
    }
    if(low[i]<0.0 && parmin[i]>=0.0) {
      if(VERBOSE) printf("   FAILED with negatives.\n\n");
      free(points); return 4;
    }
    if(up[i]>0.0 && parmax[i]<=0.0) {
      if(VERBOSE) printf("   FAILED with positives.\n\n");
      free(points); return 5;
    }
  }

  if(VERBOSE) printf("   SUCCEEDED\n");
  return 0;
}
/******************************************************************************/

/******************************************************************************/
/* test_powell() */

/** De Jong's second function or Rosenbrock's valley or Banana function.
 *  One minimum f(xi)=0. With n=2 minimum is at (1,1) and with n=3 at(1,1,1).
 *  With 4<=n<=7 there is one global minimum at (1,1,...,1) and local
 *  minimum somewhere near (-1,1,...1). For larger n the method breaks down.
 */
double optfunc_dejong2(int n, double *x, void *func_data)
{
  int i, last_n=0, dif;
  double d1, d2, f;
  double last_x[50];

  if(func_data==NULL) {} // just to prevent compiler warning
  dif=0; if(n!=last_n) dif=1;
  if(dif==0) for(i=0; i<n; i++) if(x[i]!=last_x[i]) {dif=1; break;}
  if(dif==0) {printf("objf called with the same parameters again!\n");}
  if(dif==1) {last_n=n; for(i=0; i<n; i++) last_x[i]=x[i];}

  f=0.0;
  for(i=0; i<n-1; i++) {

    //d1=fma(-x[i],x[i],x[i+1]); // would increase nr of evaluations
    d1=x[i+1]-x[i]*x[i];
    d2=1.0-x[i];
    f+=100.0*d1*d1 + d2*d2; 
  }

  return f;
}


int test_powell(int VERBOSE)
{
  double f, par[50], delta[50];
  int i, parNr, ret, iterNr;


  printf("test_powell(): Rosenbrock's valley with N=2\n");
  parNr=2; iterNr=100;
  for(i=0; i<parNr; i++) {
    par[i]=0.3*(double)(i+1); delta[i]=0.2;
  }
  ret=powell(par, delta, parNr, 0.002, &iterNr, &f, optfunc_dejong2, NULL, 0);
  if(ret>1) {
    if(VERBOSE) {printf("   FAILED: powell() returned error %d\n", ret);}
    return (1);
  }
  if(ret==1) {
    if(VERBOSE) {printf("   FAILED: powell() did not reach required tolerance\n");}
    return (1);
  }
  if(VERBOSE) {
    printf("powell() iterations: %d\n", iterNr);
    i=0; printf("estimated parameter values: %g", par[i]);
    for(i=1; i<parNr; i++) printf(", %g", par[i]);
    printf("\nestimated function minimum: %g\n", f);
  }
  for(i=0; i<parNr; i++) {
    if(fabs(par[i]-1.0)>1.0E-06) {
      if(VERBOSE) {printf("   FAILED: powell() did not reach required minimum\n");}
      return (1);
    }
  }

  
  printf("test_powell(): Rosenbrock's valley with N=3\n");
  parNr=3; iterNr=100;
  for(i=0; i<parNr; i++) {
    par[i]=0.3*(double)(i+1); delta[i]=0.2;
  }
  ret=powell(par, delta, parNr, 0.002, &iterNr, &f, optfunc_dejong2, NULL, 0);
  if(ret>1) {
    if(VERBOSE) {printf("   FAILED: powell() returned error %d\n", ret);}
    return (1);
  }
  if(ret==1) {
    if(VERBOSE) {printf("   FAILED: powell() did not reach required tolerance\n");}
    return (1);
  }
  if(VERBOSE) {
    printf("powell() iterations: %d\n", iterNr);
    i=0; printf("estimated parameter values: %g", par[i]);
    for(i=1; i<parNr; i++) printf(", %g", par[i]);
    printf("\nestimated function minimum: %g\n", f);
  }
  for(i=0; i<parNr; i++) {
    if(fabs(par[i]-1.0)>1.0E-06) {
      if(VERBOSE) {printf("   FAILED: powell() did not reach required minimum\n");}
      return (1);
    }
  }

  
  printf("test_powell(): Rosenbrock's valley with N=6\n");
  parNr=6; iterNr=100;
  for(i=0; i<parNr; i++) {
    par[i]=2.0; delta[i]=0.1;
  }
  ret=powell(par, delta, parNr, 0.002, &iterNr, &f, optfunc_dejong2, NULL, 0);
  if(ret>1) {
    if(VERBOSE) {printf("   FAILED: powell() returned error %d\n", ret);}
    return (1);
  }
  if(ret==1) {
    if(VERBOSE) {printf("   FAILED: powell() did not reach required tolerance\n");}
    return (1);
  }
  if(VERBOSE) {
    printf("powell() iterations: %d\n", iterNr);
    i=0; printf("estimated parameter values: %g", par[i]);
    for(i=1; i<parNr; i++) printf(", %g", par[i]);
    printf("\nestimated function minimum: %g\n", f);
  }
  for(i=0; i<parNr; i++) {
    if(fabs(par[i]-1.0)>1.0E-06) {
      if(VERBOSE) {printf("   FAILED: powell() did not reach required minimum\n");}
      return (1);
    }
  }
  
  printf("test_powell(): Rosenbrock's valley with N=6, starting from global min\n");
  parNr=6; iterNr=100;
  for(i=0; i<parNr; i++) {
    par[i]=1.0; delta[i]=0.1;
  }
  ret=powell(par, delta, parNr, 0.002, &iterNr, &f, optfunc_dejong2, NULL, 0);
  if(ret>1) {
    if(VERBOSE) {printf("   FAILED: powell() returned error %d\n", ret);}
    return (1);
  }
  if(ret==1) {
    if(VERBOSE) {printf("   FAILED: powell() did not reach required tolerance\n");}
    return (1);
  }
  if(VERBOSE) {
    printf("powell() iterations: %d\n", iterNr);
    i=0; printf("estimated parameter values: %g", par[i]);
    for(i=1; i<parNr; i++) printf(", %g", par[i]);
    printf("\nestimated function minimum: %g\n", f);
  }
  for(i=0; i<parNr; i++) {
    if(fabs(par[i]-1.0)>1.0E-06) {
      if(VERBOSE) {printf("   FAILED: powell() did not reach required minimum\n");}
      return (1);
    }
  }

  
  printf("test_powell(): Rosenbrock's valley with N=6, one parameters fixed\n");
  parNr=6; iterNr=100;
  for(i=0; i<parNr; i++) {
    par[i]=2.0; delta[i]=0.1;
  }
  par[3]=1.0; delta[3]=0.0;
  ret=powell(par, delta, parNr, 0.002, &iterNr, &f, optfunc_dejong2, NULL, 0);
  if(ret>1) {
    if(VERBOSE) {printf("   FAILED: powell() returned error %d\n", ret);}
    return (1);
  }
  if(ret==1) {
    if(VERBOSE) {printf("   FAILED: powell() did not reach required tolerance\n");}
    return (1);
  }
  if(VERBOSE) {
    printf("powell() iterations: %d\n", iterNr);
    i=0; printf("estimated parameter values: %g", par[i]);
    for(i=1; i<parNr; i++) printf(", %g", par[i]);
    printf("\nestimated function minimum: %g\n", f);
  }
  for(i=0; i<parNr; i++) {
    if(fabs(par[i]-1.0)>1.0E-06) {
      if(VERBOSE) {printf("   FAILED: powell() did not reach required minimum\n");}
      return (1);
    }
  }

  
  printf("test_powell(): ");
  printf("Rosenbrock's valley with N=6, starting close to local min\n");
  parNr=6; iterNr=100;
  for(i=0; i<parNr; i++) {
    if(i&1) par[i]=0.9; else par[i]=1.1; delta[i]=0.02;
  }
  par[0]=-0.95; delta[0]=0.03;
  if(VERBOSE) {
    i=0; printf("initial parameter values: %g", par[i]);
    for(i=1; i<parNr; i++) printf(", %g", par[i]);
    printf("\n");
  }
  ret=powell(par, delta, parNr, 0.00000001, &iterNr, &f, optfunc_dejong2, NULL, 0);
  if(ret>1) {
    if(VERBOSE) {printf("   FAILED: powell() returned error %d\n", ret);}
    return (1);
  }
  if(ret==1) {
    if(VERBOSE) {printf("   FAILED: powell() did not reach required tolerance\n");}
    return (1);
  }
  if(VERBOSE) {
    printf("powell() iterations: %d\n", iterNr);
    i=0; printf("estimated parameter values: %g", par[i]);
    for(i=1; i<parNr; i++) printf(", %g", par[i]);
    printf("\nestimated function minimum: %g\n", f);
  }
  if(fabs(f-3.97394)>0.001) {
    if(VERBOSE) {printf("   FAILED: powell() did not find the local minimum\n");}
    return (1);
  }



  if(VERBOSE){
    printf("   SUCCEEDED: powell() passed.\n");}

  return (0);
}
/******************************************************************************/

/******************************************************************************/
////// TGO ///////

/******************************************************************************/
int test_tgo(int VERBOSE)
{
  int    i, n, ret;
  double x[50];
  double xl[50], xu[50];
  double xtrue[50], ftrue;
  double d, f, flimit, xlimit;


  /*
   *  Testing tgo-powell with Banana function
   */
  printf("\ntest_tgo() with Banana function and Powell-Brent\n");
  TGO_LOCAL_OPT=0;
  TGO_SQUARED_TRANSF=1;
  TGO_LOCAL_INSIDE=0;
  n=7; ftrue=0.0; for(i=0; i<n; i++) xtrue[i]=1.0;
  for(i=0; i<n; i++) x[i]=0.0;
  i=0;    xl[i]=-10.0;     xu[i]=5.0;
  i++;    xl[i]=-5.0;      xu[i]=20.0;
  i++;    xl[i]=-1.0;      xu[i]=200.0;
  i++;    xl[i]=-1.0;      xu[i]=5.0;
  i++;    xl[i]=-1.0;      xu[i]=2.0;
  i++;    xl[i]=-1.0;      xu[i]=3.0;
  i++;    xl[i]= 0.0;      xu[i]=5.0;
  xlimit=5.0E-03; flimit=1.0E-04;
  if(0) {
    printf("n := %d\n", n);
    printf("Initial parameter values and limits:\n");
    for(i=0; i<n; i++) {
      printf("  x[%d]=%g   xl[%d]=%g    xu[%d]=%g\n",
        i, x[i], i, xl[i], i, xu[i]);
    }
  }
  ret=tgo(xl, xu, optfunc_dejong2, NULL, n, 10, &f, x, 1000, 1, 0);
  if(ret!=0) {
    if(VERBOSE) {printf("   FAILED: tgo() returned error %d\n", ret);}
    return (1);
  }
  if(VERBOSE) {
    printf("Optimized parameter values and true values:\n");
    for(i=0; i<n; i++) {
      printf("  x[%d]=%g   xtrue[%d]=%g\n",
        i, x[i], i, xtrue[i]);
    }
    printf("  min=%g   truemin=%g\n", f, ftrue);
  }
  ret=0;
  d=fabs(f-ftrue);
  if(d>flimit) {
    if(VERBOSE)
      fprintf(stderr, "Error: tgo() did not reach required minimum.\n");
    ret++;
  }
  for(i=0; i<n; i++) {
    d=fabs(x[i]-xtrue[i]);
    if(d>xlimit) {
      if(VERBOSE)
        fprintf(stderr, "Error: tgo() did not reach required x[%d].\n",i);
      ret++;
    }
  }
  if(ret!=0) return ret;
  if(VERBOSE)
    printf("tgo() optimization with Powell-Brent (n=%d) successful.\n", n);



  /*
   *  Testing tgo-bobyqa with Banana function
   */
  printf("\ntest_tgo() with Banana function and Bobyqa\n");
  TGO_LOCAL_OPT=1;
  TGO_SQUARED_TRANSF=1;
  TGO_LOCAL_INSIDE=0;
  n=7; ftrue=0.0; for(i=0; i<n; i++) xtrue[i]=1.0;
  for(i=0; i<n; i++) x[i]=0.0;
  i=0;    xl[i]=-10.0;     xu[i]=5.0;
  i++;    xl[i]=-5.0;      xu[i]=20.0;
  i++;    xl[i]=-1.0;      xu[i]=200.0;
  i++;    xl[i]=-1.0;      xu[i]=5.0;
  i++;    xl[i]=-1.0;      xu[i]=2.0;
  i++;    xl[i]=-1.0;      xu[i]=3.0;
  i++;    xl[i]= 0.0;      xu[i]=5.0;
  xlimit=5.0E-03; flimit=1.0E-04;
  if(0) {
    printf("n := %d\n", n);
    printf("Initial parameter values and limits:\n");
    for(i=0; i<n; i++) {
      printf("  x[%d]=%g   xl[%d]=%g    xu[%d]=%g\n",
        i, x[i], i, xl[i], i, xu[i]);
    }
  }
  ret=tgo(xl, xu, optfunc_dejong2, NULL, n, 10, &f, x, 1000, 1, 0);
  if(ret!=0) {
    if(VERBOSE) {printf("   FAILED: tgo() returned error %d\n", ret);}
    return (1);
  }
  if(VERBOSE) {
    printf("Optimized parameter values and true values:\n");
    for(i=0; i<n; i++) {
      printf("  x[%d]=%g   xtrue[%d]=%g\n",
        i, x[i], i, xtrue[i]);
    }
    printf("  min=%g   truemin=%g\n", f, ftrue);
  }
  ret=0;
  d=fabs(f-ftrue);
  if(d>flimit) {
    if(VERBOSE)
      fprintf(stderr, "Error: tgo() did not reach required minimum.\n");
    ret++;
  }
  for(i=0; i<n; i++) {
    d=fabs(x[i]-xtrue[i]);
    if(d>xlimit) {
      if(VERBOSE)
        fprintf(stderr, "Error: tgo() did not reach required x[%d].\n",i);
      ret++;
    }
  }
  if(ret!=0) return ret;
  if(VERBOSE)
    printf("tgo() optimization with Bobyqa (n=%d) successful.\n", n);



  /*
   *  Testing tgo-bobyqa with Banana function
   */
  printf("\ntest_tgo() with Generalized Rastrigin function and Powell-Brent\n");
  TGO_LOCAL_OPT=0;
  TGO_SQUARED_TRANSF=1;
  TGO_LOCAL_INSIDE=0;
  n=5; ftrue=0.0; for(i=0; i<n; i++) xtrue[i]=0.0;
  for(i=0; i<n; i++) x[i]=1.0;
  i=0;    xl[i]=-3.12;     xu[i]=2.12;
  i++;    xl[i]=-1.12;     xu[i]=3.12;
  i++;    xl[i]=-2.12;     xu[i]=3.12;
  i++;    xl[i]=-2.12;     xu[i]=2.12;
  i++;    xl[i]=-3.12;     xu[i]=1.12;
  xlimit=5.0E-03; flimit=1.0E-04;
  if(0) {
    printf("n := %d\n", n);
    printf("Initial parameter values and limits:\n");
    for(i=0; i<n; i++) {
      printf("  x[%d]=%g   xl[%d]=%g    xu[%d]=%g\n",
        i, x[i], i, xl[i], i, xu[i]);
    }
  }
  ret=tgo(xl, xu, optfunc_rastrigin, NULL, n, 10, &f, x, 1000, 1, 0);
  if(ret!=0) {
    if(VERBOSE) {printf("   FAILED: tgo() returned error %d\n", ret);}
    return (1);
  }
  if(VERBOSE) {
    printf("Optimized parameter values and true values:\n");
    for(i=0; i<n; i++) {
      printf("  x[%d]=%g   xtrue[%d]=%g\n",
        i, x[i], i, xtrue[i]);
    }
    printf("  min=%g   truemin=%g\n", f, ftrue);
  }
  ret=0;
  d=fabs(f-ftrue);
  if(d>flimit) {
    if(VERBOSE)
      fprintf(stderr, "Error: tgo() did not reach required minimum.\n");
    ret++;
  }
  for(i=0; i<n; i++) {
    d=fabs(x[i]-xtrue[i]);
    if(d>xlimit) {
      if(VERBOSE)
        fprintf(stderr, "Error: tgo() did not reach required x[%d].\n",i);
      ret++;
    }
  }
  if(ret!=0) return ret;
  if(VERBOSE)
    printf("tgo() optimization with Powell-Brent (n=%d) successful.\n", n);



  /*
   *  Testing tgo-bobyqa with Banana function
   */
  printf("\ntest_tgo() with Generalized Rastrigin function and Bobyqa\n");
  TGO_LOCAL_OPT=1;
  TGO_SQUARED_TRANSF=1;
  TGO_LOCAL_INSIDE=0;
  n=5; ftrue=0.0; for(i=0; i<n; i++) xtrue[i]=0.0;
  for(i=0; i<n; i++) x[i]=1.0;
  i=0;    xl[i]=-3.12;     xu[i]=2.12;
  i++;    xl[i]=-1.12;     xu[i]=3.12;
  i++;    xl[i]=-2.12;     xu[i]=3.12;
  i++;    xl[i]=-2.12;     xu[i]=2.12;
  i++;    xl[i]=-3.12;     xu[i]=1.12;
  xlimit=5.0E-03; flimit=1.0E-04;
  if(0) {
    printf("n := %d\n", n);
    printf("Initial parameter values and limits:\n");
    for(i=0; i<n; i++) {
      printf("  x[%d]=%g   xl[%d]=%g    xu[%d]=%g\n",
        i, x[i], i, xl[i], i, xu[i]);
    }
  }
  ret=tgo(xl, xu, optfunc_rastrigin, NULL, n, 10, &f, x, 1000, 1, 0);
  if(ret!=0) {
    if(VERBOSE) {printf("   FAILED: tgo() returned error %d\n", ret);}
    return (1);
  }
  if(VERBOSE) {
    printf("Optimized parameter values and true values:\n");
    for(i=0; i<n; i++) {
      printf("  x[%d]=%g   xtrue[%d]=%g\n",
        i, x[i], i, xtrue[i]);
    }
    printf("  min=%g   truemin=%g\n", f, ftrue);
  }
  ret=0;
  d=fabs(f-ftrue);
  if(d>flimit) {
    if(VERBOSE)
      fprintf(stderr, "Error: tgo() did not reach required minimum.\n");
    ret++;
  }
  for(i=0; i<n; i++) {
    d=fabs(x[i]-xtrue[i]);
    if(d>xlimit) {
      if(VERBOSE)
        fprintf(stderr, "Error: tgo() did not reach required x[%d].\n",i);
      ret++;
    }
  }
  if(ret!=0) return ret;
  if(VERBOSE)
    printf("tgo() optimization with Bobyqa (n=%d) successful.\n", n);


  return 0;
}
/******************************************************************************/
////// BOBYQA TESTS //////////
/******************************************************************************/

/******************************************************************************/
int test_bobyqa1(int VERBOSE)
{
  int error_code = 0;
  int i, j, m, n;
  bobyqa_result ret;

  double x[100], xl[100], xu[100], bdl, bdu, dx[100];
  double truex[100], *dp, truef=0.0, d;
  int npt, nevals, jcase;
  double twopi, minf, temp;
  double (*func)(int n, double *x, void *func_data);

  printf("\n=====================================\n");
  printf("\nTesting bobyqa with test problem 1...\n");
  printf("\n=====================================\n");

  twopi = atan(1.) * 8.;
  bdl = -1.;
  bdu = 1.;
  m = 5;
  func=bobyqa_problem1;

  if(VERBOSE>1)
    printf("Powell's Fortran code gave these results with n=20:\n");
  n=20; 
  x[0]=1.0; x[1]=1.0; x[2]=3.616077E-1; x[3]=1.0; x[4]=-3.616078E-1;
  x[5]=1.0; x[6]=-1.0; x[7]=1.0; x[8]=-1.0; x[9]=1.910563E-08; 
  x[10]=-1.0; x[11]=-1.0; x[12]=-3.616078E-1; x[13]=-1.0; x[14]=3.616080E-01; 
  x[15]=-1.0; x[16]=1.0; x[17]=-1.0; x[18]=1.0; x[19]=-1.376918E-07;
  if(VERBOSE>1) {
    printf("X is:\n");
    for(j=0; j<n; j++) {printf("%-16.5e", x[j]); if(j%5==4) printf("\n");}
  }
  minf=func(n, x, NULL);
  if(VERBOSE>1)
    printf("with these estimates the Least value of F = %.15e\n\n", minf);


  do {
    n = 2*m;
    for (i = 0; i < n; ++i) {
      xl[i] = bdl;
      xu[i] = bdu;
      dx[i] = 0.01*(bdu-bdl);
    }
    for (jcase = 1; jcase <= 2; ++jcase) {
      npt = n + 6;
      if (jcase == 2) npt = 2*n + 1;
      if(VERBOSE>1)
        printf("\n2D output with M = %d, N = %d and NPT = %d\n", m, n, npt);
      fflush(stdout);
      /* Set correct results, based on Powells Fortran program */
      if(m==5 && n==10 && npt==16) {
        dp=truex;
        *dp=0.2612470; *(++dp)=1.0; *(++dp)=-1.0; *(++dp)=1.0; *(++dp)=-1.0;
        *(++dp)=-1.0; *(++dp)=1.0; *(++dp)=-1.0; *(++dp)=1.0; *(++dp)=0.2612471;
        truef=5.680353888084283;
      } else if(m==5 && n==10 && npt==21) {
        dp=truex;
        *dp=1.0; *(++dp)=1.0; *(++dp)=-1.0; *(++dp)=1.0; *(++dp)=-1.0;
        *(++dp)=-1.0; *(++dp)=1.0; *(++dp)=-1.0; *(++dp)=1.0; *(++dp)=3.40121E-8;
        truef=5.601533972186465;
      } else if(m==10 && n==20 && npt==26) {
        dp=truex;
        *dp=1.0; *(++dp)=1.0; *(++dp)=0.3616078; *(++dp)=1.0; *(++dp)=-0.3616079;
        *(++dp)=1.0; *(++dp)=-1.0; *(++dp)=1.0; *(++dp)=-1.0; *(++dp)=1.100284E-7;
        *(++dp)=-1.0; *(++dp)=-1.0; *(++dp)=-0.3616080; *(++dp)=-1.0;
        *(++dp)=0.3616079;
        *(++dp)=-1.0; *(++dp)=1.0; *(++dp)=-1.0; *(++dp)=1.0;
        *(++dp)=-2.025736E-07;
        truef=3.220305336883057E+01;
      } else if(m==10 && n==20 && npt==41) {
        dp=truex;
        *dp=1.0; *(++dp)=1.0; *(++dp)=0.3616077; *(++dp)=1.0; *(++dp)=-0.3616079;
        *(++dp)=1.0; *(++dp)=-1.0; *(++dp)=1.0; *(++dp)=-1.0; 
        *(++dp)=1.910563E-08;
        *(++dp)=-1.0; *(++dp)=-1.0; *(++dp)=-0.3616078; *(++dp)=-1.0;
        *(++dp)=0.3616080;
        *(++dp)=-1.0; *(++dp)=1.0; *(++dp)=-1.0; *(++dp)=1.0;
        *(++dp)=-1.376918E-07;
        truef=3.220305336883060E+01;
      }
      /* Calculate initial guesses */
      for (j = 1; j <= m; ++j) {
        temp = (double)j * twopi / (double)m;
        x[2*j - 2] = cos(temp);
        x[2*j - 1] = sin(temp);
      }
      for(i=0; i<n; i++) { // check that initial guesses are inside limits
        if(x[i]<xl[i]) x[i]=xl[i];
        if(x[i]>xu[i]) x[i]=xu[i];
      }
      if(VERBOSE>1) {
        printf("Alkuarvaukset X is:\n");
        for(j=0; j<n; j++) {
          printf("%-16.5e", x[j]); if(j%5==4) printf("\n");
        }
      }
      minf=func(n, x, NULL);
      if(VERBOSE>1)
        printf("... joiden F = %.15e\n", minf);

      if(VERBOSE>1) {printf("\nseuraavaksi bobyqa()\n"); fflush(stdout);}
      ret=bobyqa(n, npt, x, xl, xu, dx, 0.0, 1.0E-08, 0.01, 1.0E-12, 1.0E-12,
                 1000, &nevals, &minf, func, NULL, NULL, 0);
      if(ret<1) {printf("error in bobyqa!\n"); return 1;}
      if(VERBOSE>1) {
        printf("At the return from BOBYQA   Number of function calls = %d\n",
               nevals);
        printf("Least value of F = %.15e    The corresponding X is:\n", minf);
        for(j=0; j<n; j++) {printf("%-16.5e", x[j]); if(j%5==4) printf("\n");}
      }
      /* Check that results are at least as good as with Fortran program */
      if(VERBOSE>1)
        printf("Fortran program gave F = %.15e\n", truef);
      if(minf>truef*1.00000000001) {
        printf("F is worse than with Powell's original Fortran SW\n");
        error_code=10; break;
      }
      for(j=0; j<n; j++) {
        d=x[j]-truex[j];
        if(fabs(d)>5.0E-06) {
          printf("fitted parameter differs too much from original Fortran SW\n");
          if(error_code!=0) {error_code=11; break;}
          printf("but F is better than with Fortran so it is ok\n");
        }
      }
      if(error_code!=0) break;
    }
    m += m;
  } while(m<=10 && error_code==0);

  
  if(error_code)
    printf("\n    Test FAILED: test_bobyqa1 failed with error code: %i\n",
           error_code);
  else
    printf("\n    Test SUCCESFULL: test_bobyqa1 exited with: %i\n", error_code); 

  return(error_code);
}
/******************************************************************************/

/******************************************************************************/
int test_constraints1(int VERBOSE)
{
  int error_code = 0;
  int i, j, m, n;
  bobyqa_result ret;

  double x[100], xl[100], xu[100], bdl, bdu, dx[100];
  double truex[100], *dp, truef=0.0, d, dmax;
  int npt, nevals;
  double twopi, minf, temp;
  double (*func)(int n, double *x, void *func_data);

  printf("\n=====================================\n");
  printf("\nTesting bobyqa with constraints 1...\n");
  printf("\n=====================================\n");

  twopi = atan(1.) * 8.;
  bdl = -1.;
  bdu = 1.;
  m=10;
  n=20; 
  npt = 2*n + 1;
  func=bobyqa_problem1;

  for (i = 0; i < n; ++i) {
    xl[i] = bdl;
    xu[i] = bdu;
    dx[i] = 0.1*(bdu-bdl);
  }
  //printf("2D output with M = %d, N = %d and NPT = %d\n", m, n, npt);
  /* Set correct results, based on Powells Fortran program */
  dp=truex;
  *dp=1.0; *(++dp)=1.0; *(++dp)=0.3616077; *(++dp)=1.0; *(++dp)=-0.3616078;
  *(++dp)=1.0; *(++dp)=-1.0; *(++dp)=1.0; *(++dp)=-1.0; *(++dp)=1.910563E-08;
  *(++dp)=-1.0; *(++dp)=-1.0; *(++dp)=-0.3616078; *(++dp)=-1.0; *(++dp)=0.3616080;
  *(++dp)=-1.0; *(++dp)=1.0; *(++dp)=-1.0; *(++dp)=1.0; *(++dp)=-1.376918E-07;
  if(VERBOSE>1) {
    printf("Powell's Fortran code gave these results with n=20:\n");
    printf("X is:\n");
    for(j=0; j<n; j++) {printf("%-16.5e", truex[j]); if(j%5==4) printf("\n");}
  }
  truef=func(n, truex, NULL);
  if(VERBOSE>1)
    printf("with these estimates the Least value of F = %.15e\n\n", truef);
  /* Calculate initial guesses */
  for (j = 1; j <= m; ++j) {
    temp = (double)j * twopi / (double)m;
    x[2*j - 2] = cos(temp);
    x[2*j - 1] = sin(temp);
  }
  for(i=0; i<n; i++) { // check that initial guesses ae inside limits
    if(x[i]<xl[i]) x[i]=xl[i];
    if(x[i]>xu[i]) x[i]=xu[i];
  }
  /* Fix certain parameters */
  i=0; x[i]=xl[i]=xu[i]=truex[i];
  i=6; x[i]=xl[i]=xu[i]=truex[i];
  i=17; x[i]=xl[i]=xu[i]=truex[i];
  i=18; x[i]=xl[i]=xu[i]=truex[i];
  npt = 2*(n-4) + 1; // recalc NPT based on fitted N
  if(VERBOSE>1) {
    printf("Alkuarvaukset X is:\n");
    for(j=0; j<n; j++) {
      printf("%-16.5e", x[j]); if(j%5==4) printf("\n");
    }
  }
  minf=func(n, x, NULL);
  if(VERBOSE>1)
    printf("... joiden F = %.15e\n", minf);

  if(VERBOSE>1) {printf("\nseuraavaksi bobyqa()\n"); fflush(stdout);}
  ret=bobyqa(n, npt, x, xl, xu, dx, 1.0E-08, 1.0E-06, 0.01, 1.0E-13, 1.0E-12,
             1000, &nevals, &minf, func, NULL, NULL, 0);
  if(ret<1) {printf("error in bobyqa!\n"); return 1;}
  if(VERBOSE>1) {
    printf("At the return from BOBYQA   Number of function calls = %d\n", nevals);
    printf("Least value of F = %.15e    The corresponding X is:\n", minf);
    for(j=0; j<n; j++) {printf("%-16.5e", x[j]); if(j%5==4) printf("\n");}
  }
  /* Check that results are at least as good as with Fortran program */
  if(minf>truef*1.00000000001) {
    printf("F is worse than with Powell's original Fortran SW\n");
    error_code=+10;
  }
  for(j=0, dmax=0.0; j<n; j++) {
    d=fabs(x[j]-truex[j]); if(d>dmax) dmax=d;
    if(d>5.0E-06) {
      if(VERBOSE>0) {
        printf("   FAILED: bobyqa() did not reach required par%d estimate\n", i+1);
        printf("           |estimate-true|=%.20E\n", d);
      }
      error_code+=11; //break;
    }
  }
  if(VERBOSE>0) printf("Max abs parameter difference: %.15E\n", dmax);

  if(error_code)
    printf("\n    Test FAILED: test_constraints1 failed with error code: %i\n",
           error_code);
  else
    printf("\n    Test SUCCESFULL: test_constraints1 exited with: %i\n",
           error_code); 

  return(error_code);
}
/******************************************************************************/

/******************************************************************************/
int test_scales1(int VERBOSE)
{
  int error_code = 0;
  int i, j, m, n;
  bobyqa_result ret;

  double x[100], xl[100], xu[100], bdl, bdu, dx[100];
  double truex[100], *dp, truef=0.0, d, dmax;
  int npt, jcase, nevals;
  double twopi, minf, temp;
  double (*func)(int n, double *x, void *func_data);

  printf("\n=======================================================\n");
  printf("\nTesting bobyqa with parameters of different scales 1...\n");
  printf("\n=======================================================\n");

  twopi = atan(1.) * 8.;
  bdl = -1.;
  bdu = 1.;
  m = 5;
  func=bobyqa_problem2;

  do {
    n = 2*m;
    for (i = 0; i < n; ++i) {
      xl[i] = bdl;
      xu[i] = bdu;
      if(i%5==4) {xl[i]*=10.; xu[i]*=10.;}
      else if(i%5==3) {xl[i]*=0.001; xu[i]*=0.001;}
      dx[i] = 0.08*(xu[i]-xl[i]);
    }
    for (jcase = 1; jcase <= 2; ++jcase) {
      npt = n + 6;
      if (jcase == 2) npt = 2*n + 1;
      if(VERBOSE>1)
        printf("\n2D output with M = %d, N = %d and NPT = %d\n", m, n, npt);
      fflush(stdout);
      /* Set correct results, based on Powells Fortran program, but
         with changed scales */
      if(m==5 && n==10 && npt==16) {
        dp=truex;
        *dp=0.2612470; *(++dp)=1.0; *(++dp)=-1.0; *(++dp)=1.0E-03; 
        *(++dp)=-1.0E+01;
        *(++dp)=-1.0; *(++dp)=1.0; *(++dp)=-1.0; *(++dp)=1.0E-03; 
        *(++dp)=2.612471E+00;
        truef=5.680353888084283;
      } else if(m==5 && n==10 && npt==21) {
        dp=truex;
        *dp=1.0; *(++dp)=1.0; *(++dp)=-1.0; *(++dp)=1.0E-03; *(++dp)=-1.0E+01;
        *(++dp)=-1.0; *(++dp)=1.0; *(++dp)=-1.0; *(++dp)=1.0E-03; 
        *(++dp)=3.40121E-03;
        truef=5.601533972186465;
      } else if(m==10 && n==20 && npt==26) {
        dp=truex;
        *dp=1.0; *(++dp)=1.0; *(++dp)=0.3616078; *(++dp)=1.0E-03; 
        *(++dp)=-3.616079E+00;
        *(++dp)=1.0; *(++dp)=-1.0; *(++dp)=1.0; *(++dp)=-1.0E-03; 
        *(++dp)=1.100284E-6;
        *(++dp)=-1.0; *(++dp)=-1.0; *(++dp)=-0.3616080; *(++dp)=-1.0E-03; 
        *(++dp)=3.616079E+00;
        *(++dp)=-1.0; *(++dp)=1.0; *(++dp)=-1.0; *(++dp)=1.0E-03;
        *(++dp)=-2.025736E-06;
        truef=3.220305336883057E+01;
      } else if(m==10 && n==20 && npt==41) {
        dp=truex;
        *dp=1.0; *(++dp)=1.0; *(++dp)=0.3616077; *(++dp)=1.0E-03;
        *(++dp)=-3.616079E+00;
        *(++dp)=1.0; *(++dp)=-1.0; *(++dp)=1.0; *(++dp)=-1.0E-03;
        *(++dp)=1.910563E-07;
        *(++dp)=-1.0; *(++dp)=-1.0; *(++dp)=-0.3616078; *(++dp)=-1.0E-03; 
        *(++dp)=3.616080E+00;
        *(++dp)=-1.0; *(++dp)=1.0; *(++dp)=-1.0; *(++dp)=1.0E-03; 
        *(++dp)=-1.376918E-06;
        truef=3.220305336883060E+01;
      }
      if(VERBOSE>1) {
        printf("Powell's Fortran code gave these results:\n");
        printf("X is:\n");
        for(j=0; j<n; j++) {printf("%-16.5e", truex[j]); if(j%5==4) printf("\n");}
      }
      truef=func(n, truex, NULL);
      if(VERBOSE>1)
        printf("with these estimates the Least value of F = %.15e\n\n", truef);
      /* Calculate initial guesses */
      for (j = 1; j <= m; ++j) {
        temp = (double) j * twopi / (double) m;
        x[2*j - 2] = cos(temp);
        x[2*j - 1] = sin(temp);
      }
      for(i=0; i<n; i++) {
        if(i%5==4) {
          x[i]*=10.;
        } else if(i%5==3) {
          x[i]*=0.001;
        }
      }
      for(i=0; i<n; i++) { // check that initial guesses are inside limits
        if(x[i]<xl[i]) x[i]=xl[i];
        if(x[i]>xu[i]) x[i]=xu[i];
      }
      if(VERBOSE>1) {
        printf("Alkuarvaukset X is:\n");
        for(j=0; j<n; j++) {
          printf("%-16.5e", x[j]); if(j%5==4) printf("\n");
        }
        fflush(stdout);
      }
      minf=func(n, x, NULL);
      if(VERBOSE>1)
        printf("... joiden F = %.15e\n", minf);

      if(VERBOSE>1) {printf("\nseuraavaksi bobyqa()\n"); fflush(stdout);}
      ret=bobyqa(n, npt, x, xl, xu, dx, 1.0E-08, 1.0E-08, 0.01, 1.0E-12, 1.0E-12,
                 1000, &nevals, &minf, func, NULL, NULL, 0);
      if(ret<1) {printf("error in bobyqa!\n"); return 1;}
      if(VERBOSE>1) {
        printf("At the return from BOBYQA   Number of function calls = %d\n",
               nevals);
        printf("Least value of F = %.15e    The corresponding X is:\n", minf);
        for(j=0; j<n; j++) {printf("%-16.5e", x[j]); if(j%5==4) printf("\n");}
      }
      /* Check that results are at least as good as with Fortran program */
      if(minf>truef*1.0000000001) {
        printf("F is worse than with Powell's original Fortran SW\n");
        printf("   minf= %.15E\n", minf);
        printf("   truef=%.15E\n", truef);
        error_code+=10; break;
      }
      for(j=0, dmax=0.0; j<n; j++) {
        d=fabs(x[j]-truex[j]); if(d>dmax) dmax=d;
        if(d>5.0E-03) {
          if(minf<0.99999999*truef) {
            printf("fitted parameter differs from original Fortran SW, ");
            printf("but min is better\n");
          } else {
            printf("fitted parameter differs too much from original Fortran SW\n");
            if(VERBOSE) {
              printf("   FAILED: bobyqa() did not reach required par%d estimate\n",
                     i+1);
              printf("           |estimate-true|=%.20E\n", d);
            }
            error_code+=11; //break;
          }
        }
      }
      if(VERBOSE) printf("Max abs parameter difference: %.15E\n", dmax);
      if(error_code!=0) break;
    }
    if(error_code!=0) break;
    m += m;
  } while(m <= 10);


  
  if(error_code)
    printf("\n    Test FAILED: test_scales1 failed with error code: %i\n",
           error_code);
  else
    printf("\n    Test SUCCESFULL: test_scales1 exited with: %i\n", error_code); 

  return(error_code);
}
/******************************************************************************/

/******************************************************************************/
int test_scales2(int VERBOSE)
{
  int error_code = 0;
  int i, j, m, n;
  bobyqa_result ret;

  double x[100], xl[100], xu[100], bdl, bdu, dx[100];
  double truex[100], *dp, truef=0.0, d;
  int npt, nevals;
  double twopi, temp, minf;
  double (*func)(int n, double *x, void *func_data);

  printf("\n=======================================================\n");
  printf("\nTesting bobyqa with constraints and fixed parameters...\n");
  printf("\n=======================================================\n");

  twopi = atan(1.) * 8.;
  bdl = -1.;
  bdu = 1.;
  m=10;
  n=20; 
  npt = 2*n + 1;
  func=bobyqa_problem2;

  for (i = 0; i < n; ++i) {
    xl[i] = bdl;
    xu[i] = bdu;
    if(i%5==4) {xl[i]*=10.; xu[i]*=10.;}
    else if(i%5==3) {xl[i]*=0.001; xu[i]*=0.001;}
    dx[i] = 0.07*(xu[i]-xl[i]);
  }
  /* Set correct results, based on Powells Fortran program */
  dp=truex;
  *dp=1.0; *(++dp)=1.0; *(++dp)=0.3616077; *(++dp)=1.0; *(++dp)=-0.3616078;
  *(++dp)=1.0; *(++dp)=-1.0; *(++dp)=1.0; *(++dp)=-1.0; *(++dp)=1.910563E-08;
  *(++dp)=-1.0; *(++dp)=-1.0; *(++dp)=-0.3616078; *(++dp)=-1.0; *(++dp)=0.3616080;
  *(++dp)=-1.0; *(++dp)=1.0; *(++dp)=-1.0; *(++dp)=1.0; *(++dp)=-1.376918E-07;
  for (i = 0; i < n; ++i) {
    if(i%5==4) truex[i]*=10.; else if(i%5==3) truex[i]*=0.001;}
  if(VERBOSE>1) {
    printf("Powell's Fortran code gave these results with n=20:\n");
    printf("X is:\n");
    for(j=0; j<n; j++) {printf("%-16.5e", truex[j]); if(j%5==4) printf("\n");}
  }
  truef=func(n, truex, NULL);
  if(VERBOSE>1)
    printf("with these estimates the Least value of F = %.15e\n\n", truef);
  /* Calculate initial guesses */
  for (j = 1; j <= m; ++j) {
    temp = (double) j * twopi / (double) m;
    x[2*j - 2] = cos(temp);
    x[2*j - 1] = sin(temp);
  }
  for(i=0; i<n; i++) {
    if(i%5==4) {
      x[i]*=10.;
    } else if(i%5==3) {
      x[i]*=0.001;
    }
  }
  for(i=0; i<n; i++) { // check that initial guesses ae inside limits
    if(x[i]<xl[i]) x[i]=xl[i];
    if(x[i]>xu[i]) x[i]=xu[i];
  }
  /* Fix certain parameters */
  i=0; x[i]=xl[i]=xu[i]=truex[i]; dx[i]=0.0;
  i=6; x[i]=xl[i]=xu[i]=truex[i]; dx[i]=0.0;
  i=17; x[i]=xl[i]=xu[i]=truex[i]; dx[i]=0.0;
  i=18; x[i]=xl[i]=xu[i]=truex[i]; dx[i]=0.0;
  npt = 2*(n-4) + 1; // recalc NPT based on fitted N
  if(VERBOSE>1) {
    printf("Alkuarvaukset X is:\n");
    for(j=0; j<n; j++) {
      printf("%-16.5e", x[j]); if(j%5==4) printf("\n");
    }
  }
  minf=func(n, x, NULL);
  if(VERBOSE>1)
    printf("... joiden F = %.15e\n", minf);

  if(VERBOSE>1) {printf("\nseuraavaksi bobyqa()\n"); fflush(stdout);}
  ret=bobyqa(n, npt, x, xl, xu, dx, 1.0E-06, 1.0E-12, 30.0, 1.0E-12, 1.0E-12,
             1000, &nevals, &minf, func, NULL, NULL, 0);
  if(ret<1) {printf("error in bobyqa!\n"); return 1;}
  if(VERBOSE>1) {
    printf("At the return from BOBYQA   Number of function calls = %d\n", nevals);
    printf("Least value of F = %.15e    The corresponding X is:\n", minf);
    for(j=0; j<n; j++) {printf("%-16.5e", x[j]); if(j%5==4) printf("\n");}
  }
  /* Check that results are at least as good as with Fortran program */
  if(minf>truef*1.00000000001) {
    printf("F is worse than with Powell's original Fortran SW\n");
    error_code=+10;
  }
  for(j=0; j<n; j++) {
    d=x[j]-truex[j];
    if(fabs(d)>5.0E-03) {
      if(minf<0.99999999*truef) {
        printf("fitted parameter differs from original Fortran SW, \n");
        printf("but min is better\n");
      } else {
        printf("fitted parameter differs too much from original Fortran SW\n");
        error_code+=11; break;
      }
    }
  }
  
  if(error_code)
    printf("\n    Test FAILED: test_scales2 failed with error code: %i\n",
           error_code);
  else
    printf("\n    Test SUCCESFULL: test_scales2 exited with: %i\n", error_code); 

  return(error_code);
}
/******************************************************************************/

/******************************************************************************/
int test_onedim1(int VERBOSE)
{
  int error_code = 0;
  int i, j, m, n;
  bobyqa_result ret;

  double x[100], xl[100], xu[100], bdl, bdu, dx[100];
  double truex[100], *dp, truef=0.0, d;
  int npt, nevals;
  double twopi, temp, minf;
  double (*func)(int n, double *x, void *func_data);

  printf("\n=======================================================\n");
  printf("\nTesting bobyqa with all parameters fixed except one...\n");
  printf("\n=======================================================\n");

  twopi = atan(1.) * 8.;
  bdl = -1.;
  bdu = 1.;
  m=10;
  n=20; 
  npt = 2*n + 1;
  func=bobyqa_problem2;

  for (i = 0; i < n; ++i) {
    xl[i] = bdl;
    xu[i] = bdu;
    if(i%5==4) {xl[i]*=10.; xu[i]*=10.;}
    else if(i%5==3) {xl[i]*=0.001; xu[i]*=0.001;}
    dx[i] = 0.07*(xu[i]-xl[i]);
  }
  /* Set correct results, based on Powells Fortran program */
  dp=truex;
  *dp=1.0; *(++dp)=1.0; *(++dp)=0.3616077; *(++dp)=1.0; *(++dp)=-0.3616078;
  *(++dp)=1.0; *(++dp)=-1.0; *(++dp)=1.0; *(++dp)=-1.0; *(++dp)=1.910563E-08;
  *(++dp)=-1.0; *(++dp)=-1.0; *(++dp)=-0.3616078; *(++dp)=-1.0; *(++dp)=0.3616080;
  *(++dp)=-1.0; *(++dp)=1.0; *(++dp)=-1.0; *(++dp)=1.0; *(++dp)=-1.376918E-07;
  for (i = 0; i < n; ++i) {
    if(i%5==4) truex[i]*=10.; else if(i%5==3) truex[i]*=0.001;}
  if(VERBOSE>1) {
    printf("Powell's Fortran code gave these results with n=20:\n");
    printf("X is:\n");
    for(j=0; j<n; j++) {printf("%-16.5e", truex[j]); if(j%5==4) printf("\n");}
  }
  truef=func(n, truex, NULL);
  if(VERBOSE>1)
    printf("with these estimates the Least value of F = %.15e\n\n", truef);
  /* Calculate initial guesses */
  for (j = 1; j <= m; ++j) {
    temp = (double)j * twopi / (double)m;
    x[2*j - 2] = cos(temp);
    x[2*j - 1] = sin(temp);
  }
  for(i=0; i<n; i++) {
    if(i%5==4) {
      x[i]*=10.;
    } else if(i%5==3) {
      x[i]*=0.001;
    }
  }
  for(i=0; i<n; i++) { // check that initial guesses are inside limits
    if(x[i]<xl[i]) x[i]=xl[i];
    if(x[i]>xu[i]) x[i]=xu[i];
  }
  /* Fix certain parameters */
  for(i=0; i<n; i++) if(i!=2) {
    x[i]=xl[i]=xu[i]=truex[i]; dx[i]=0.0;
  }
  npt = 2*(1) + 1; // recalc NPT based on fitted N
  printf("Alkuarvaukset X is:\n");
  for(j=0; j<n; j++) {
     printf("%-16.5e", x[j]); if(j%5==4) printf("\n");
  }
  minf=func(n, x, NULL);
  printf("... joiden F = %.15e\n", minf);

  printf("\nseuraavaksi bobyqa()\n"); fflush(stdout);
  ret=bobyqa(n, npt, x, xl, xu, dx, 1.0E-06, 1.0E-10, 30.0, 1.0E-12, 1.0E-12,
             1000, &nevals, &minf, func, NULL, NULL, 0);
  if(ret<1) {printf("error in bobyqa!\n"); return 1;}
  printf("At the return from BOBYQA   Number of function calls = %d\n", nevals);
  printf("Least value of F = %.15e    The corresponding X is:\n", minf);
  for(j=0; j<n; j++) {printf("%-16.5e", x[j]); if(j%5==4) printf("\n");}
  /* Check that results are at least as good as with Fortran program */
  if(minf>truef*1.00000000001) {
    printf("F is worse than with Powell's original Fortran SW\n");
    error_code=+10;
  }
  for(j=0; j<n; j++) {
    d=x[j]-truex[j];
    if(fabs(d)>5.0E-03) {
      if(minf<0.99999999*truef) {
        printf("fitted parameter differs from original Fortran SW, ");
        printf("but min is better\n");
      } else {
        printf("fitted parameter differs too much from original Fortran SW\n");
        error_code+=11; break;
      }
    }
  }
  
  if(error_code)
    printf("\n    Test FAILED: test_scales2 failed with error code: %i\n",
           error_code);
  else
    printf("\n    Test SUCCESFULL: test_scales2 exited with: %i\n", error_code); 

  return(error_code);
}
/******************************************************************************/

/******************************************************************************/
int test_onedim2(int VERBOSE)
{
  int error_code = 0;
  int i, j, m, n;
  bobyqa_result ret;

  double x[100], xl[100], xu[100], bdl, bdu, dx[100];
  double truex[100], *dp, truef=0.0, d;
  int npt, nevals;
  double twopi, temp, minf;
  double (*func)(int n, double *x, void *func_data);

  printf("\n=======================================================\n");
  printf("\nTesting bobyqa with all parameters fixed except one,\n");
  printf("and the fitted one has its minimum at the limit...\n");
  printf("\n=======================================================\n");

  twopi = atan(1.) * 8.;
  bdl = -1.;
  bdu = 1.;
  m=10;
  n=20; 
  npt = 2*n + 1;
  func=bobyqa_problem2;

  for (i = 0; i < n; ++i) {
    xl[i] = bdl;
    xu[i] = bdu;
    if(i%5==4) {xl[i]*=10.; xu[i]*=10.;}
    else if(i%5==3) {xl[i]*=0.001; xu[i]*=0.001;}
    dx[i] = 0.07*(xu[i]-xl[i]);
  }
  /* Set correct results, based on Powells Fortran program */
  dp=truex;
  *dp=1.0; *(++dp)=1.0; *(++dp)=0.3616077; *(++dp)=1.0; *(++dp)=-0.3616078;
  *(++dp)=1.0; *(++dp)=-1.0; *(++dp)=1.0; *(++dp)=-1.0; *(++dp)=1.910563E-08;
  *(++dp)=-1.0; *(++dp)=-1.0; *(++dp)=-0.3616078; *(++dp)=-1.0; *(++dp)=0.3616080;
  *(++dp)=-1.0; *(++dp)=1.0; *(++dp)=-1.0; *(++dp)=1.0; *(++dp)=-1.376918E-07;
  for (i = 0; i < n; ++i) {
    if(i%5==4) truex[i]*=10.; else if(i%5==3) truex[i]*=0.001;}
  if(VERBOSE>1) {
    printf("Powell's Fortran code gave these results with n=20:\n");
    printf("X is:\n");
    for(j=0; j<n; j++) {printf("%-16.5e", truex[j]); if(j%5==4) printf("\n");}
  }
  truef=func(n, truex, NULL);
  printf("with these estimates the Least value of F = %.15e\n\n", truef);
  /* Calculate initial guesses */
  for (j = 1; j <= m; ++j) {
    temp = (double)j * twopi / (double)m;
    x[2*j - 2] = cos(temp);
    x[2*j - 1] = sin(temp);
  }
  for(i=0; i<n; i++) {
    if(i%5==4) {
      x[i]*=10.;
    } else if(i%5==3) {
      x[i]*=0.001;
    }
  }
  for(i=0; i<n; i++) { // check that initial guesses are inside limits
    if(x[i]<xl[i]) x[i]=xl[i];
    if(x[i]>xu[i]) x[i]=xu[i];
  }
  /* Fix certain parameters */
  for(i=0; i<n; i++) if(i!=0) {
    x[i]=xl[i]=xu[i]=truex[i]; dx[i]=0.0;
  }
  npt = 2*(1) + 1; // recalc NPT based on fitted N
  printf("Alkuarvaukset X is:\n");
  for(j=0; j<n; j++) {
     printf("%-16.5e", x[j]); if(j%5==4) printf("\n");
  }
  minf=func(n, x, NULL);
  printf("... joiden F = %.15e\n", minf);

  printf("\nseuraavaksi bobyqa()\n"); fflush(stdout);
  ret=bobyqa(n, npt, x, xl, xu, dx, 1.0E-08, 1.0E-12, 30.0, 1.0E-15, 1.0E-15,
             1000, &nevals, &minf, func, NULL, NULL, 0);
  if(ret<1) {printf("error in bobyqa!\n"); return 1;}
  printf("At the return from BOBYQA   Number of function calls = %d\n", nevals);
  printf("Least value of F = %.15e    The corresponding X is:\n", minf);
  for(j=0; j<n; j++) {printf("%-16.5e", x[j]); if(j%5==4) printf("\n");}
  /* Check that results are at least as good as with Fortran program */
  if(minf>truef*1.00000000001) {
    printf("F is worse than with Powell's original Fortran SW\n");
    error_code=+10;
  }
  for(j=0; j<n; j++) {
    d=x[j]-truex[j];
    if(fabs(d)>5.0E-03) {
      if(minf<0.99999999*truef) {
        printf("fitted parameter differs from original Fortran SW, ");
        printf("but min is better\n");
      } else {
        printf("fitted parameter differs too much from original Fortran SW\n");
        error_code+=11; break;
      }
    }
  }
  
  if(error_code)
    printf("\n    Test FAILED: test_scales2 failed with error code: %i\n", error_code);
  else
    printf("\n    Test SUCCESFULL: test_scales2 exited with: %i\n", error_code); 

  return(error_code);
}
/******************************************************************************/

/******************************************************************************/
int test_banana1(int VERBOSE)
{
  double f, par[50], delta[50], parl[50], paru[50], dif;
  int i, parNr, ret, fevalNr, maxFevalNr, error_code=0;

  printf("\n=======================================================\n");
  printf("\nTesting bobyqa with De Jong's second function = Rosenbrock's valley\n");
  printf("= Banana function, at its basic form...\n");
  printf("\n=======================================================\n");

  drandSeed(1); //srand(time(0));

  printf("\nRosenbrock's valley with N=2\n");
  parNr=2; maxFevalNr=1000;
  for(i=0; i<parNr; i++) {
    par[i]=0.3*(double)(i+1); delta[i]=0.01; parl[i]=-1.0E+03; paru[i]=+1.0E+03;
  }
  if(VERBOSE) {
    i=0; printf("initial parameter values: %g", par[i]);
    for(i=1; i<parNr; i++) printf(", %g", par[i]);
    printf("\n");
    f=optfunc_dejong2(parNr, par, NULL);
    printf("function value with initial estimates: %g\n", f);
  }
  ret=bobyqa(parNr, 0, par, parl, paru, delta, 1.0E-06, 1.0E-06, 1.0E-10,
             1.0E-03, 1.0E-10, maxFevalNr,
             &fevalNr, &f, optfunc_dejong2, NULL, NULL, 0);
  if(ret<1) {printf("error in bobyqa!\n"); return 1;}
  if(VERBOSE) {
    printf("bobyqa() return code: %d\n", ret);
    printf("bobyqa() function call nr: %d\n", fevalNr);
    i=0; printf("estimated parameter values: %.10E", par[i]);
    for(i=1; i<parNr; i++) printf(", %.10E", par[i]);
    printf("\nestimated function minimum: %.10E\n", f);
  }
  if(f>1.0E-06) {
    if(VERBOSE) {printf("   FAILED: bobyqa() did not reach required minimum\n");}
    error_code=10;
  }
  for(i=0; i<parNr; i++) {
    dif=fabs(par[i]-1.0);
    if(dif>1.0E-05) {
      if(VERBOSE) {
        printf("   FAILED: bobyqa() did not reach required par%d estimate\n", i+1);
        printf("           |estimate-true|=%.20E\n", dif);
      }
      error_code=11;
    }
  }
  if(error_code) {
    printf("\n    Test FAILED: test_banana1 failed with error code: %i\n",
           error_code);
    return(error_code);
  }

  printf("\nRosenbrock's valley with N=3\n");
  parNr=3; maxFevalNr=1000;
  for(i=0; i<parNr; i++) {
    par[i]=0.3*(double)(i+1); delta[i]=0.1; parl[i]=-1.0E+03; paru[i]=+1.0E+03;
  }
  if(VERBOSE) {
    i=0; printf("initial parameter values: %g", par[i]);
    for(i=1; i<parNr; i++) printf(", %g", par[i]);
    printf("\n");
    f=optfunc_dejong2(parNr, par, NULL);
    printf("function value with initial estimates: %g\n", f);
  }
  ret=bobyqa(parNr, 0, par, parl, paru, delta, 1.0E-06, 1.0E-06, 1.0E-10,
             1.0E-06, 1.0E-10, maxFevalNr,
             &fevalNr, &f, optfunc_dejong2, NULL, NULL, 0);
  if(ret<1) {printf("error in bobyqa!\n"); return 1;}
  if(VERBOSE) {
    printf("bobyqa() return code: %d\n", ret);
    printf("bobyqa() function call nr: %d\n", fevalNr);
    i=0; printf("estimated parameter values: %.10E", par[i]);
    for(i=1; i<parNr; i++) printf(", %.10E", par[i]);
    printf("\nestimated function minimum: %.10E\n", f);
  }
  if(f>1.0E-06) {
    if(VERBOSE) {printf("   FAILED: bobyqa() did not reach required minimum\n");}
    error_code=20;
  }
  for(i=0; i<parNr; i++) {
    dif=fabs(par[i]-1.0);
    if(dif>1.0E-05) {
      if(VERBOSE) {
        printf("   FAILED: bobyqa() did not reach required par%d estimate\n", i+1);
        printf("           |estimate-true|=%.20E\n", dif);
      }
      error_code=21;
    }
  }
  if(error_code) {
    printf("\n    Test FAILED: test_banana1 failed with error code: %i\n",
           error_code);
    return(error_code);
  }


  printf("\nRosenbrock's valley with N=6\n");
  parNr=6; maxFevalNr=1000;
  for(i=0; i<parNr; i++) {
    par[i]=2.0; delta[i]=0.2; parl[i]=-1.0E+03; paru[i]=+1.0E+03;
  }
  if(VERBOSE) {
    i=0; printf("initial parameter values: %g", par[i]);
    for(i=1; i<parNr; i++) printf(", %g", par[i]);
    printf("\n");
    f=optfunc_dejong2(parNr, par, NULL);
    printf("function value with initial estimates: %g\n", f);
  }
  ret=bobyqa(parNr, 0, par, parl, paru, delta, 1.0E-06, 1.0E-06, 1.0E-12,
             1.0E-06, 1.0E-14, maxFevalNr,
             &fevalNr, &f, optfunc_dejong2, NULL, NULL, 0);
  if(ret<1) {printf("error in bobyqa!\n"); return 1;}
  if(VERBOSE) {
    printf("bobyqa() return code: %d\n", ret);
    printf("bobyqa() function call nr: %d\n", fevalNr);
    i=0; printf("estimated parameter values: %.10E", par[i]);
    for(i=1; i<parNr; i++) printf(", %.10E", par[i]);
    printf("\nestimated function minimum: %.10E\n", f);
  }
  if(f>1.0E-06) {
    if(VERBOSE) {printf("   FAILED: bobyqa() did not reach required minimum\n");}
    error_code=30;
  }
  for(i=0; i<parNr; i++) {
    dif=fabs(par[i]-1.0);
    if(dif>5.0E-05) {
      if(VERBOSE) {
        printf("   FAILED: bobyqa() did not reach required par%d estimate\n", i+1);
        printf("           |estimate-true|=%.20E\n", dif);
      }
      error_code=31;
    }
  }
  if(error_code) {
    printf("\n    Test FAILED: test_banana1 failed with error code: %i\n",
           error_code);
    return(error_code);
  }


  printf("\nRosenbrock's valley with N=6, starting from global min\n");
  parNr=6; maxFevalNr=1000;
  for(i=0; i<parNr; i++) {
    par[i]=1.0; delta[i]=0.1; parl[i]=-1.0E+03; paru[i]=+1.0E+03;
  }
  if(VERBOSE) {
    i=0; printf("initial parameter values: %g", par[i]);
    for(i=1; i<parNr; i++) printf(", %g", par[i]);
    printf("\n");
    f=optfunc_dejong2(parNr, par, NULL);
    printf("function value with initial estimates: %g\n", f);
  }
  ret=bobyqa(parNr, 0, par, parl, paru, delta, 1.0E-06, 1.0E-06, 1.0E-10,
             1.0E-06, 1.0E-10, maxFevalNr,
             &fevalNr, &f, optfunc_dejong2, NULL, NULL, 0);
  if(ret<1) {printf("error in bobyqa!\n"); return 1;}
  if(VERBOSE) {
    printf("bobyqa() return code: %d\n", ret);
    printf("bobyqa() function call nr: %d\n", fevalNr);
    i=0; printf("estimated parameter values: %.10E", par[i]);
    for(i=1; i<parNr; i++) printf(", %.10E", par[i]);
    printf("\nestimated function minimum: %.10E\n", f);
  }
  if(f>1.0E-10) {
    if(VERBOSE) {printf("   FAILED: bobyqa() did not reach required minimum\n");}
    error_code=40;
  }
  for(i=0; i<parNr; i++) {
    dif=fabs(par[i]-1.0);
    if(dif>1.0E-10) {
      if(VERBOSE) {
        printf("   FAILED: bobyqa() did not reach required par%d estimate\n", i+1);
        printf("           |estimate-true|=%.20E\n", dif);
      }
      error_code=41;
    }
  }
  if(error_code) {
    printf("\n    Test FAILED: test_banana1 failed with error code: %i\n",
           error_code);
    return(error_code);
  }


  printf("\nRosenbrock's valley with N=6, one parameter fixed\n");
  parNr=6; maxFevalNr=1000;
  for(i=0; i<parNr; i++) {
    par[i]=2.0; delta[i]=0.5; parl[i]=-1.0E+03; paru[i]=+1.0E+03;
  }
  i=3; par[i]=1.0; delta[i]=0.0;
  if(VERBOSE) {
    i=0; printf("initial parameter values: %g", par[i]);
    for(i=1; i<parNr; i++) printf(", %g", par[i]);
    printf("\n");
    f=optfunc_dejong2(parNr, par, NULL);
    printf("function value with initial estimates: %g\n", f);
  }
  ret=bobyqa(parNr, 0, par, parl, paru, delta, 1.0E-06, 1.0E-06, 1.0E-10,
             1.0E-06, 1.0E-10, maxFevalNr,
             &fevalNr, &f, optfunc_dejong2, NULL, NULL, 0);
  if(ret<1) {printf("error in bobyqa!\n"); return 1;}
  if(VERBOSE) {
    printf("bobyqa() return code: %d\n", ret);
    printf("bobyqa() function call nr: %d\n", fevalNr);
    i=0; printf("estimated parameter values: %.10E", par[i]);
    for(i=1; i<parNr; i++) printf(", %.10E", par[i]);
    printf("\nestimated function minimum: %.10E\n", f);
  }
  if(f>1.0E-06) {
    if(VERBOSE) {printf("   FAILED: bobyqa() did not reach required minimum\n");}
    error_code=50;
  }
  for(i=0; i<parNr; i++) {
    dif=fabs(par[i]-1.0);
    if(dif>1.0E-05) {
      if(VERBOSE) {
        printf("   FAILED: bobyqa() did not reach required par%d estimate\n", i+1);
        printf("           |estimate-true|=%.20E\n", dif);
      }
      error_code=51;
    }
  }
  if(error_code) {
    printf("\n    Test FAILED: test_banana1 failed with error code: %i\n",
           error_code);
    return(error_code);
  }


  printf("\nRosenbrock's valley with N=6, starting close to local min\n");
  parNr=6; maxFevalNr=1000;
  for(i=0; i<parNr; i++) {
    if(i&1) par[i]=0.9; else par[i]=1.1;
    delta[i]=0.5; parl[i]=-1.0E+03; paru[i]=+1.0E+03;
  }
  i=0; par[i]=-0.95; delta[i]=0.03;
  if(VERBOSE) {
    i=0; printf("initial parameter values: %g", par[i]);
    for(i=1; i<parNr; i++) printf(", %g", par[i]);
    printf("\n");
    f=optfunc_dejong2(parNr, par, NULL);
    printf("function value with initial estimates: %g\n", f);
  }
  ret=bobyqa(parNr, 0, par, parl, paru, delta, 1.0E-06, 1.0E-06, 1.0E-10,
             1.0E-08, 1.0E-10, maxFevalNr,
             &fevalNr, &f, optfunc_dejong2, NULL, NULL, 0);
  if(ret<1) {printf("error in bobyqa!\n"); return 1;}
  if(VERBOSE) {
    printf("bobyqa() return code: %d\n", ret);
    printf("bobyqa() function call nr: %d\n", fevalNr);
    i=0; printf("estimated parameter values: %.10E", par[i]);
    for(i=1; i<parNr; i++) printf(", %.10E", par[i]);
    printf("\nestimated function minimum: %.10E\n", f);
  }
  if(fabs(f-3.97394)>0.001) {
    if(VERBOSE) {printf("   FAILED: bobyqa() did not find the local minimum\n");}
    error_code=70;
  }
  for(i=0; i<parNr; i++) {
    if(i==0) dif=fabs(par[i]-(-1.0));
    else dif=fabs(par[i]-1.0);
    if(dif>0.3) {
      if(VERBOSE) {
        printf("   FAILED: bobyqa() did not reach required par%d estimate\n", i+1);
        printf("           for local min, |estimate-true|=%.20E\n", dif);
      }
      error_code=71;
    }
  }
  if(error_code) {
    printf("\n    Test FAILED: test_banana1 failed with error code: %i\n",
      error_code);
    return(error_code);
  }


  printf("\nRosenbrock's valley with N=6, trying to get rescue() called\n");
  printf("and to get ROUNOFF_LIMITED error because of too many iterations.\n");
  parNr=6; maxFevalNr=10000;
  for(i=0; i<parNr; i++) {
    if(i==0) par[i]=2.0; else if(i&1) par[i]=-4.1; else par[i]=3.0;
    delta[i]=0.5;
    if(i>0) parl[i]=-1.0E+01; else parl[i]=0.0;
    paru[i]=+5.0E+01;
  }
  if(VERBOSE) {
    i=0; printf("initial parameter values: %g", par[i]);
    for(i=1; i<parNr; i++) printf(", %g", par[i]);
    f=optfunc_dejong2(parNr, par, NULL);
    printf("function value with initial estimates: %g\n", f);
  }
  ret=bobyqa(parNr, 0, par, parl, paru, delta, 1.0E-100, 1.0E-100, 1.0E-100,
             1.0E-100, 1.0E-100, maxFevalNr,
             &fevalNr, &f, optfunc_dejong2, NULL, NULL, 0);
  if(ret<1) {printf("error in bobyqa!\n"); /*return 1;*/}
  if(ret!=BOBYQA_ROUNDOFF_LIMITED) {
    if(VERBOSE) {printf("   FAILED: bobyqa() did not stop in ROUNOFF error\n");}
    error_code=89;
  }
  if(VERBOSE) {
    printf("bobyqa() return code: %d\n", ret);
    printf("bobyqa() function call nr: %d\n", fevalNr);
    i=0; printf("estimated parameter values: %.10E", par[i]);
    for(i=1; i<parNr; i++) printf(", %.10E", par[i]);
    printf("\nestimated function minimum: %.10E\n", f);
  }
  if(fabs(f-0.0)>1.0E-10) {
    if(VERBOSE) {printf("   FAILED: bobyqa() did not find the global minimum\n");}
    error_code=80;
  }
  for(i=0; i<parNr; i++) {
    dif=fabs(par[i]-1.0);
    if(dif>1.0E-10) {
      if(VERBOSE) {
        printf("   FAILED: bobyqa() did not reach required par%d estimate\n", i+1);
        printf("           |estimate-true|=%.20E\n", dif);
      }
      error_code=81;
    }
  }
  if(error_code) {
    printf("\n    Test FAILED: test_banana1 failed with error code: %i\n",
           error_code);
    return(error_code);
  }



  printf("\n    Test SUCCESFULL: test_banana1 exited with: %i\n", error_code); 
  return(error_code);
}

/******************************************************************************/

/******************************************************************************/
int test_rastrigin(int VERBOSE)
{
  double f, fo, par[50], delta[50], parl[50], paru[50], dif;
  int i, parNr, ret, fevalNr, maxFevalNr, error_code=0;

  printf("\n===============================================================\n");
  printf("\nTesting bobyqa with Rastrigin function; this function should be\n");
  printf("used to test global optimization, not local optimization like\n");
  printf("bobyqa; therefore this now only tests that bobyqa() does not crash\n");
  printf("and finds at least some estimates.\n");
  printf("\n===============================================================\n");

  printf("\nTesting that global min is found with suitable start values\n");
  parNr=2; maxFevalNr=1000;
  for(i=0; i<parNr; i++) {
    par[i]=0.3; delta[i]=0.1; parl[i]=-5.12; paru[i]=+5.12;
  }
  if(VERBOSE) {
    i=0; printf("initial parameter values: %g", par[i]);
    for(i=1; i<parNr; i++) printf(", %g", par[i]);
    printf("\n");
    f=optfunc_rastrigin(parNr, par, NULL);
    printf("function value with initial estimates: %g\n", f);
  }
  ret=bobyqa(parNr, 0, par, parl, paru, delta, 1.0E-06, 1.0E-06, 1.0E-10,
             1.0E-03, 1.0E-10, maxFevalNr,
             &fevalNr, &f, optfunc_rastrigin, NULL, NULL, 0);
  if(ret<1) {printf("error in bobyqa!\n"); return 1;}
  if(VERBOSE) {
    printf("bobyqa() return code: %d\n", ret);
    printf("bobyqa() function call nr: %d\n", fevalNr);
    i=0; printf("estimated parameter values: %.10E", par[i]);
    for(i=1; i<parNr; i++) printf(", %.10E", par[i]);
    printf("\nestimated function minimum: %.10E\n", f);
  }
  if(f>1.0E-06) {
    if(VERBOSE) {printf("   FAILED: bobyqa() did not reach required minimum\n");}
    error_code=10;
  }
  for(i=0; i<parNr; i++) {
    dif=fabs(par[i]-0.0);
    if(dif>1.0E-08) {
      if(VERBOSE) {
        printf("   FAILED: bobyqa() did not reach required par%d estimate\n", i+1);
        printf("           |estimate-true|=%.20E\n", dif);
      }
      error_code=11;
    }
  }
  if(error_code) {
    printf("\n    Test FAILED: test_rastrigin failed with error code: %i\n",
           error_code);
    return(error_code);
  }


  printf("\nTesting that local min is found, depending on start values\n");
  parNr=2; maxFevalNr=1000;
  for(i=0; i<parNr; i++) {
    par[i]=4.3; delta[i]=0.1; parl[i]=-5.12; paru[i]=+5.12;
  }
  if(VERBOSE) {
    i=0; printf("initial parameter values: %g", par[i]);
    for(i=1; i<parNr; i++) printf(", %g", par[i]);
    printf("\n");
    f=optfunc_rastrigin(parNr, par, NULL);
    printf("function value with initial estimates: %g\n", f);
  }
  ret=bobyqa(parNr, 0, par, parl, paru, delta, 1.0E-06, 1.0E-06, 1.0E-10,
             1.0E-08, 1.0E-06, maxFevalNr,
             &fevalNr, &f, optfunc_rastrigin, NULL, NULL, 0);
  if(ret<1) {printf("error in bobyqa!\n"); return 1;}
  if(VERBOSE) {
    printf("bobyqa() return code: %d\n", ret);
    printf("bobyqa() function call nr: %d\n", fevalNr);
    i=0; printf("estimated parameter values: %.10E", par[i]);
    for(i=1; i<parNr; i++) printf(", %.10E", par[i]);
    printf("\nestimated function minimum: %.10E\n", f);
  }
  if(f<1.0) {
    if(VERBOSE) {printf("   FAILED: bobyqa() did not reach required minimum\n");}
    error_code=20;
  }
  fo=f;
  if(error_code==0) {
    par[0]-=1.0E-04; par[1]-=1.0E-04;
    f=optfunc_rastrigin(parNr, par, NULL);
    printf("close-by function minimum:  %.10E\n", f);
    if(f<fo) {
      if(VERBOSE) {printf("   FAILED: bobyqa() did not reach local minimum\n");}
      error_code=21;
    }
  }
  if(error_code==0) {
    par[0]+=2.0E-04; par[1]+=2.0E-04;
    f=optfunc_rastrigin(parNr, par, NULL);
    printf("close-by function minimum:  %.10E\n", f);
    if(f<fo) {
      if(VERBOSE) {printf("   FAILED: bobyqa() did not reach local minimum\n");}
      error_code=22;
    }
  }
  if(error_code==0) {
    par[0]-=2.0E-04; par[1]+=0.0E-04;
    f=optfunc_rastrigin(parNr, par, NULL);
    printf("close-by function minimum:  %.10E\n", f);
    if(f<fo) {
      if(VERBOSE) {printf("   FAILED: bobyqa() did not reach local minimum\n");}
      error_code=23;
    }
  }
  if(error_code==0) {
    par[0]+=0.0E-04; par[1]-=2.0E-04;
    f=optfunc_rastrigin(parNr, par, NULL);
    printf("close-by function minimum:  %.10E\n", f);
    if(f<fo) {
      if(VERBOSE) {printf("   FAILED: bobyqa() did not reach local minimum\n");}
      error_code=24;
    }
  }
  if(error_code) {
    printf("\n    Test FAILED: test_rastrigin failed with error code: %i\n",
           error_code);
    return(error_code);
  }

  printf("\n    Test SUCCESFULL: test_rastrigin exited with: %i\n", error_code); 
  return(error_code);
}

/******************************************************************************/

/******************************************************************************/
int test_nptrange(int VERBOSE)
{
  double f, par[50], delta[50], parl[50], paru[50], dif;
  int i, parNr, ret, fevalNr, maxFevalNr, error_code=0, npt;

  printf("\n=======================================================\n");
  printf("\nTesting bobyqa with a range of NPT parameters, because\n");
  printf("other tests use default npt=2*n+1.\n");
  printf("Banana function with N=6 is used as test function...\n");
  printf("\n=======================================================\n");

  printf("\nSet NPT to too high value (npt>(N+1)(N+2)/2):\n");
  parNr=6; maxFevalNr=1000;
  for(i=0; i<parNr; i++) {
    par[i]=2.0; delta[i]=0.5; parl[i]=-1.0E+03; paru[i]=+1.0E+03;
  }
  npt=(parNr+1)*(parNr+2)/2 + 1;
  if(VERBOSE) {
    printf("npt := %d\n", npt);
    i=0; printf("initial parameter values: %g", par[i]);
    for(i=1; i<parNr; i++) printf(", %g", par[i]);
    printf("\n");
    f=optfunc_dejong2(parNr, par, NULL);
    printf("function value with initial estimates: %g\n", f);
  }
  ret=bobyqa(parNr, npt, par, parl, paru, delta, 1.0E-08, 1.0E-08, 1.0E-12,
             1.0E-10, 1.0E-12, maxFevalNr,
             &fevalNr, &f, optfunc_dejong2, NULL, NULL, 0);
  if(ret!=BOBYQA_INVALID_ARGS) {
    if(VERBOSE)
      printf("   FAILED: bobyqa() did not give error about too high NPT\n");
    error_code=10;
  }
  if(error_code) {
    printf("\n    Test FAILED: test_nptrange failed with error code: %i\n",
           error_code);
    return(error_code);
  }


  printf("\nSet NPT to too low value (npt<(N+2):\n");
  parNr=6; maxFevalNr=1000;
  for(i=0; i<parNr; i++) {
    par[i]=2.0; delta[i]=0.5; parl[i]=-1.0E+03; paru[i]=+1.0E+03;
  }
  npt=(parNr+2) - 1;
  if(VERBOSE) {
    printf("npt := %d\n", npt);
    i=0; printf("initial parameter values: %g", par[i]);
    for(i=1; i<parNr; i++) printf(", %g", par[i]);
    printf("\n");
    f=optfunc_dejong2(parNr, par, NULL);
    printf("function value with initial estimates: %g\n", f);
  }
  ret=bobyqa(parNr, npt, par, parl, paru, delta, 1.0E-08, 1.0E-08, 1.0E-12,
             1.0E-10, 1.0E-12, maxFevalNr,
             &fevalNr, &f, optfunc_dejong2, NULL, NULL, 0);
  if(ret!=BOBYQA_INVALID_ARGS) {
    if(VERBOSE)
      printf("   FAILED: bobyqa() did not give error about too low NPT\n");
    error_code=20;
  }
  if(error_code) {
    printf("\n    Test FAILED: test_nptrange failed with error code: %i\n",
           error_code);
    return(error_code);
  }


  printf("\nSet NPT to its maximum (not recommended) (npt=(N+1)(N+2)/2):\n");
  parNr=6; maxFevalNr=1000;
  for(i=0; i<parNr; i++) {
    par[i]=2.0; delta[i]=0.5; parl[i]=-1.0E+03; paru[i]=+1.0E+03;
  }
  npt=(parNr+1)*(parNr+2)/2;
  if(VERBOSE) {
    printf("npt := %d\n", npt);
    i=0; printf("initial parameter values: %g", par[i]);
    for(i=1; i<parNr; i++) printf(", %g", par[i]);
    printf("\n");
    f=optfunc_dejong2(parNr, par, NULL);
    printf("function value with initial estimates: %g\n", f);
  }
  ret=bobyqa(parNr, npt, par, parl, paru, delta, 1.0E-08, 1.0E-10, 1.0E-12,
             1.0E-12, 1.0E-14, maxFevalNr,
             &fevalNr, &f, optfunc_dejong2, NULL, NULL, 0);
  if(ret==BOBYQA_INVALID_ARGS) {
    if(VERBOSE) printf("   FAILED: bobyqa() gave error about wrong arguments\n");
    error_code=30;
  }
  if(VERBOSE) {
    printf("bobyqa() return code: %d\n", ret);
  }
  if(ret<1) {printf("error in bobyqa!\n"); return 1;}
  if(VERBOSE) {
    printf("bobyqa() function call nr: %d\n", fevalNr);
    i=0; printf("estimated parameter values: %.10E", par[i]);
    for(i=1; i<parNr; i++) printf(", %.10E", par[i]);
    printf("\nestimated function minimum: %.10E\n", f);
  }
  if(f>1.0E-06) {
    if(VERBOSE) {printf("   FAILED: bobyqa() did not reach required minimum\n");}
    error_code=33;
  }
  for(i=0; i<parNr; i++) {
    dif=fabs(par[i]-1.0);
    if(dif>1.0E-05) {
      if(VERBOSE) {
        printf("   FAILED: bobyqa() did not reach required par%d estimate\n", i+1);
        printf("           |estimate-true|=%.20E\n", dif);
      }
      error_code=35;
    }
  }
  if(error_code) {
    printf("\n    Test FAILED: test_nptrange failed with error code: %i\n",
           error_code);
    return(error_code);
  }



  printf("\nSet NPT to its maximum (not recommended) (npt=(N+1)(N+2)/2)\n");
  printf("AND try to get rescue() called, too:\n\n");
  parNr=6; maxFevalNr=10000;
  for(i=0; i<parNr; i++) {
    if(i&1) par[i]=-3.0; else par[i]=+4.1;
    delta[i]=0.5; parl[i]=-1.0E+01; paru[i]=+5.0E+01;
  }
  npt=(parNr+1)*(parNr+2)/2;
  if(VERBOSE) {
    printf("npt := %d\n", npt);
    i=0; printf("initial parameter values: %g", par[i]);
    for(i=1; i<parNr; i++) printf(", %g", par[i]);
    printf("\n");
    f=optfunc_dejong2(parNr, par, NULL);
    printf("function value with initial estimates: %g\n", f);
  }
  ret=bobyqa(parNr, npt, par, parl, paru, delta, 1.0E-100, 1.0E-100, 1.0E-100,
             1.0E-100, 1.0E-100, maxFevalNr,
             &fevalNr, &f, optfunc_dejong2, NULL, NULL, 0);
  if(ret<1) {printf("error in bobyqa!\n"); /*return 1;*/}
  if(ret!=BOBYQA_ROUNDOFF_LIMITED && ret!=BOBYQA_MINF_MAX_REACHED) {
    /* Sometimes returns just BOBYQA_MINF_MAX_REACHED and not ROUNDOFF, 
       that is okayish */
    if(VERBOSE) printf("   FAILED: bobyqa() did not stop in ROUNDOFF error\n");
    if(VERBOSE) printf("   bobyqa() := %d\n", ret);
    error_code=50;
  }
  if(VERBOSE) {
    printf("bobyqa() return code: %d\n", ret);
    printf("bobyqa() function call nr: %d\n", fevalNr);
    i=0; printf("estimated parameter values: %.10E", par[i]);
    for(i=1; i<parNr; i++) printf(", %.10E", par[i]);
    printf("\n");
    printf("estimated function minimum: %.10E\n", f);
  }
  if(fabs(f-0.0)>1.0E-10) {
    if(VERBOSE) printf("   FAILED: bobyqa() did not find the global minimum\n");
    error_code=51;
  }
  for(i=0; i<parNr; i++) {
    dif=fabs(par[i]-1.0);
    if(dif>1.0E-10) {
      if(VERBOSE) {
        printf("   FAILED: bobyqa() did not reach required par%d estimate\n", i+1);
        printf("           |estimate-true|=%.20E\n", dif);
      }
      error_code=52;
    }
  }
  if(error_code) {
    printf("\n    Test FAILED: test_nptrange failed with error code: %i\n",
           error_code);
    return(error_code);
  }


  printf("\n    Test SUCCESFULL: test_nptrange exited with: %i\n", error_code); 
  return(error_code);
}

/******************************************************************************/

/******************************************************************************/
/* BOBYQA test problems: */

/** Test problem for BOBYQA, the objective function being the sum of 
    the reciprocals of all pairwise distances between the points P_I, 
    I=1,2,...,M in two dimensions, where M=N/2 and where the components 
    of P_I are X(2*I-1) and X(2*I). Thus each vector X of N variables 
    defines the M points P_I. The initial X gives equally spaced points 
    on a circle. Four different choices of the pairs (N,NPT) are tried, 
    namely (10,16), (10,21), (20,26) and (20,41). Convergence to a local 
    minimum that is not global occurs in both the N=10 cases. The details
    of the results are highly sensitive to computer rounding errors. The
    choice IPRINT=2 provides the current X and optimal F so far whenever
    RHO is reduced. The bound constraints of the problem require every 
    component of X to be in the interval [-1,1]. */
double bobyqa_problem1(int n, double *x, void *func_data)
{
  int i, j, k;
  double d1, d2, f, temp;


  if(func_data==NULL) {} // just to prevent compiler warning
  f=0.0;
  for (i = 3; i < n; i += 2) {
    k = i - 1;
    for (j = 1; j < k; j += 2) {
      d1 = x[i-1] - x[j-1];
      d2 = x[i] - x[j];
      temp = d1*d1 + d2*d2;
      temp = fmax(temp,1e-6);
      f += 1. / sqrt(temp);
    }
  }
  return f;
}

/** Same function as previous, expect that best fit parameters are of very
 *  different scales, some 1000-fold higher and some 1000-fold lower than
 *  others. */ 
double bobyqa_problem2(int n, double *x, void *func_data)
{
  int i, j, k;
  double d1, d2, f, temp;
  double localx[100];

#if(0)
  printf("func()\n"); fflush(stdout);
  printf("  x={%g", x[0]);
  for(j=1; j<n; j++) printf(", %g", x[j]); printf("}\n");
#endif
  if(func_data==NULL) {} // just to prevent compiler warning

  for(i=0; i<n; i++) {
    localx[i]=x[i];
    if(i%5==4) localx[i]*=0.1;
    else if(i%5==3) localx[i]*=1000.;
  }
  f=0.0;
  for (i = 3; i < n; i += 2) {
    k = i - 1;
    for (j = 1; j < k; j += 2) {
      d1 = localx[i-1] - localx[j-1];
      d2 = localx[i] - localx[j];
      temp = d1*d1 + d2*d2;
      temp = fmax(temp,1e-6);
      f += 1. / sqrt(temp);
    }
  }
  return f;
}

/** Generalized Rastrigin function,
 *  f(x) = A*n + SUM(i=1..n; (xi^2 - A*cos(2*PI*xi))),
 *  where A=10 and xi=[-5.12,5.12].
 *  It has global minimum value f(x)=0.0 at xi=0.0,
 *  and large nr of local minima around it.
 */
double optfunc_rastrigin(int n, double *x, void *func_data)
{
  int i, last_n=0, dif;
  double f;
  double last_x[50];

  if(func_data==NULL) {} // just to prevent compiler warning
  if(n<1 || n>50) return(nan(""));
  for(i=0; i<n; i++) if(x[i]<-5.12 || x[i]>5.12) return(nan(""));
  dif=0; if(n!=last_n) dif=1;
  if(dif==0) for(i=0; i<n; i++) if(x[i]!=last_x[i]) {dif=1; break;}
  if(dif==0) {printf("objf called with the same parameters again!\n");}
  if(dif==1) {last_n=n; for(i=0; i<n; i++) last_x[i]=x[i];}

  f=0.0;
  for(i=0; i<n; i++)
    f+=fma(x[i], x[i], -10.0*cos(2.0*M_PI*x[i]));
    //f+=x[i]*x[i]-10.0*cos(2.0*M_PI*x[i]);
  f+=10.0*n;

  return f;
}
/******************************************************************************/

/******************************************************************************/
/* Simple Objective functions working like in model fitting programs,
   requiring global arrays simdata, measdata, pmin, pmax, p, and w, and
   variable fitframeNr */

/* Returns WSS between parameter value and 1000 at each sample */
double func_deviation(int parNr, double *p, void *fdata)
{
  int fi;
  double wss=0.0, d;
  
  if(fdata==NULL) {} // just to prevent compiler warning
  if(parNr==0) {}    // just to prevent compiler warning

  /* Simulate the curve */
  for(fi=0; fi<fitframeNr; fi++)
    simdata[fi]=p[0];
  
  /* Calculate wss */
  for(fi=0, wss=0.; fi<fitframeNr; fi++) if(w[fi]>0.0) {
    /* residual */
    d=measdata[fi]-simdata[fi];
    /* error */
    wss+=w[fi]*d*d;
  }

  return(wss);  
}
  
/******************************************************************************/

/******************************************************************************/
