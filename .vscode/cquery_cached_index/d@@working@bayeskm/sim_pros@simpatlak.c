// simpatlak.c
// simulate the activity using patlak graphical model, work for FDG and other irreversible tracers
//           CROI(t)=Ki*integ(CP(t))+Int*CP(t)
// or when a reference region used
//           CROI(t)=Ki'*integ(Cref(t))+Int'*CP(t)



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



int simPatlak(unsigned int frameNr,
              double Ki,
              double Vb,
              double *t0,
              double *t1, 
              double *ctt,
              double tstart, 
              double tstop,
              double *output,
              unsigned int verbose)         //Pure c does not require extern  "C" 
{
 double  *ci, *ici, *t, *ct, *theta, *dv;
 DFT input, temp, data;  //data dummy
 int ri=0, voiNr=1;
 int ret;
 char *status;

//   frameNr= *(unsigned int*) argv[0];
//   Ki     =  *(double*) argv[1];
//   Vb     =  *(double*) argv[2];
//   t0     =  (double*) argv[3];
//   t1     =  (double*) argv[4];
//   ctt    =  (double*) argv[5];          // reference region?
//   tstart = *(double*) argv[6];
//   tstop  = *(double*) argv[7]; 
//   output =  (double*) argv[8];
//   verbose  = *(unsigned int*) argv[9];



dftInit(&input); dftInit(&temp); dftInit(&data);
if(dftSetmem(&input, frameNr, voiNr)) {
    printf("out of memory\n");}
if(dft_nr_of_NA(&input)>0) {  // check if file contains NAs (missing values)
    printf("Error: missing values in %s");
    dftEmpty(&input); return(2);
} 
if(dftSetmem(&temp, frameNr, voiNr)) {
    printf("out of memory\n");}
if(dft_nr_of_NA(&temp)>0) {  // check if file contains NAs (missing values)
    printf("Error: missing values in %s");
    dftEmpty(&temp); return(2);
} 
if(dftSetmem(&data, frameNr, voiNr)) {
    printf("out of memory\n");}
if(dft_nr_of_NA(&data)>0) {  // check if file contains NAs (missing values)
    printf("Error: missing values in %s");
    dftEmpty(&data); return(2);
} 

data.voiNr=voiNr; data.frameNr=frameNr;
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
strcpy(cunit, data.unit);
temp.timeunit=2;
temp.timetype=3;


for (int i=0; i<frameNr; i++) { 
    data.x1[i] = *(t0+i);
    data.x2[i] = t1[i];
    data.x[i]=0.5*(data.x1[i]+data.x2[i]);
    // data.voi[ri].y[i]= tac[i];
    temp.x1[i] = t0[i];
    temp.x2[i] = t1[i];
    temp.x[i]=0.5*(temp.x1[i]+temp.x2[i]);
    temp.voi[ri].y[i]= ctt[i];
}


/* Interpolate and integrate data to pet times */
ret=dftInterpolate(&temp, &data, &input, status, verbose);
dftEmpty(&temp);


ci    = (double *)malloc(frameNr*sizeof(double));
ici   = (double *)malloc(frameNr*sizeof(double));
t     = (double *)malloc(frameNr*sizeof(double));

for (int i=0; i<data.frameNr; i++) { 
    // ct[i]=data.voi[ri].y[i];
    ci[i]=input.voi[ri].y[i];
    ici[i]=input.voi[ri].y2[i];
    t[i]=input.x[i];
    // theta[i] = data.voi[ri].y2[i];  
    // dv[i] = data.voi[ri].y3[i];

    output[i] = Ki*ici[i] + Vb*ci[i] ;
}

    
if (verbose>9) {
    // FILE *ppfile = fopen(debugfile, "a+");
    printf( "CT %f %f supplied\n", ct[0],ct[frameNr-1]);
    printf( "CI %f %f supplied\n", ci[0],ci[frameNr-1]);
    printf( "ici %f %f supplied\n", ici[0],ici[frameNr-1]);
    printf( "t %f %f supplied\n", t[0],t[frameNr-1]);        
    // fclose(ppfile);
}



return(0);

}