/** @file simdispersion.c
 *  @brief Simulation of dispersion.
 */
/*****************************************************************************/
#include "tpcclibConfig.h"
/*****************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
/*****************************************************************************/
#include "tpccm.h"
/*****************************************************************************/

/*****************************************************************************/
/** Simulate the effect of dispersion on a time-activity curve.

    The units of rate constants must be related to the TAC time units;
    1/min and min, or 1/sec and sec.
   
    @return Function returns 0 when successful, else a value >= 1.
    @author Vesa Oikonen
    @sa simC1
 */
int simDispersion(
  /** Array of sample times. */
  double *x,
  /** Array of sample values, which will be replaced here by dispersion added values. */
  double *y,
  /** Nr of samples. */
  const int n,
  /** First dispersion time constant (zero if no dispersion); in same time unit as sample times. */ 
  const double tau1,
  /** 2nd dispersion time constant (zero if no dispersion); in same time unit as sample times. */ 
  const double tau2,
  /** Array for temporary data, for at least n samples; enter NULL to let function to 
      allocate and free the temporary space. */
  double *tmp
) {
  /* Check input */
  if(x==NULL || y==NULL || n<2) return 1;
  if(tau1<0.0 || tau2<0.0) return 2;

  /* Allocate memory if not allocated by user */
  double *buf;
  if(tmp!=NULL) buf=tmp; else buf=(double*)malloc(n*sizeof(double));
  if(buf==NULL) return 3;

  /* First dispersion */
  if(tau1>0.0) {
    double k=1.0/tau1;
    int ret=simC1(x, y, n, k, k, buf); 
    if(ret!=0) {
      if(tmp==NULL) free(buf);
      return 100+ret;
    }
    for(int i=0; i<n; i++) y[i]=buf[i];
  }

  /* Second dispersion */
  if(tau2>0.0) {
    double k=1.0/tau2;
    int ret=simC1(x, y, n, k, k, buf);
    if(ret!=0) {
      if(tmp==NULL) free(buf);
      return 200+ret;
    }
    for(int i=0; i<n; i++) y[i]=buf[i];
  }

  if(tmp==NULL) free(buf);
  return 0;
}
/*****************************************************************************/

/*****************************************************************************/
