/** @file sim1cm.c
 *  @brief Simulation of 1-tissue compartmental models.
 */
/*****************************************************************************/
/* #include "tpcclibConfig.h" */
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
/** Simulate myocardial tissue TAC using Iida's compartment model.
    
    @details
    Memory for ct must be allocated in the calling program.
    The units of rate constants must be related to the time unit;
    1/min and min, or 1/sec and sec.
   
    @return Function returns 0 when succesful, else a value >= 1.
    @author Vesa Oikonen
    @sa simC1, simOxygen
 */
int simMBF(
  /** Array of time values */
  double *t,
  /** Input activities */
  double *ci,
  /** Number of values in TACs */
  const int nr,
  /** Apparent k1 */
  const double k1,
  /** Apparent k2 */
  const double k2,
  /** Vfit */
  const double Vfit,
  /** Pointer for TAC array to be simulated; must be allocated */
  double *ct
) {
  int i;
  double dt2;
  double cii, ci_last, t_last;
  double ct_last, cti, cti_last;


  /* Check for data */
  if(nr<2) return(1);
  if(t==NULL || ci==NULL || ct==NULL) return(2);

  /* Calculate curves */
  t_last=0.0; cii=ci_last=0.0;
  cti=ct_last=cti_last=0.0;
  for(i=0; i<nr; i++) {
    /* delta time / 2 */
    dt2=0.5*(t[i]-t_last);
    /* calculate values */
    if(dt2<0.0) {
      return(5);
    } else if(dt2>0.0) {
      /* input integral */
      cii+=(ci[i]+ci_last)*dt2;
      /* Tissue compartment and its integral */
      ct[i] = (Vfit*ci[i] + k1*cii - k2*(cti_last+dt2*ct_last)) / (1.0 + dt2*k2);
      cti = cti_last + dt2*(ct_last+ct[i]);
    }
    /* set very small values to zero */
    if(fabs(ct[i])<1.0e-12) ct[i]=0.0;
    /* prepare to the next loop */
    t_last=t[i]; ci_last=ci[i];
    ct_last=ct[i]; cti_last=cti;
  }
  return(0);
}
/*****************************************************************************/

/*****************************************************************************/
/** Simulate tissue TAC using 1 tissue compartmental model and plasma TAC,
    at plasma TAC times.
     
    @details
    Memory for ct must be allocated in the calling program.
  
    The units of rate constants must be related to the time unit; 1/min and min,
    or 1/sec and sec.
   
    @sa simC1_i, simMBF, simC2, simOxygen
    @return Function returns 0 when succesful, else a value >= 1.
    @author Vesa Oikonen
 */
int simC1(
  /** Array of time values */
  double *t,
  /** Array of arterial activities */
  double *ca,
  /** Number of values in TACs */
  const int nr,
  /** Rate constant of the model */
  const double k1,
  /** Rate constant of the model */
  const double k2,
  /** Pointer for TAC array to be simulated; must be allocated */
  double *ct
) {
  int i;
  double dt2;
  double cai, ca_last, t_last;
  double ct1, ct1_last;
  double ct1i, ct1i_last;


  /* Check for data */
  if(nr<2) return 1;
  if(t==NULL || ca==NULL || ct==NULL) return 2;

  /* Check actual parameter number */
  if(!(k1>=0.0)) return 3;

  /* Calculate curves */
  t_last=0.0; if(t[0]<t_last) t_last=t[0]; 
  cai=ca_last=0.0;
  ct1_last=ct1i_last=0.0;
  ct1=ct1i=0.0;
  for(i=0; i<nr; i++) {
    /* delta time / 2 */
    dt2=0.5*(t[i]-t_last);
    /* calculate values */
    if(dt2<0.0) {
      return 5;
    } else if(dt2>0.0) {
      /* arterial integral */
      cai+=(ca[i]+ca_last)*dt2;
      /* tissue compartment and its integral */
      ct1 = (k1*cai - k2*(ct1i_last+dt2*ct1_last)) / (1.0 + dt2*k2);
      ct1i = ct1i_last + dt2*(ct1_last+ct1);
    }
    /* copy values to argument arrays; set very small values to zero */
    ct[i]=ct1; if(fabs(ct[i])<1.0e-12) ct[i]=0.0;
    /* prepare to the next loop */
    t_last=t[i]; ca_last=ca[i];
    ct1_last=ct1; ct1i_last=ct1i;
  }

  return 0;
}
/*****************************************************************************/

/*****************************************************************************/
/** Simulate tissue TAC using 1 tissue compartmental model and plasma TAC,
    at plasma TAC times.
     
    @details
    This version uses integral of arterial TAC as input function.
    Only advantage over simC1() is that the calculation of integral can be
    fully controlled and possibly more precice in some situations.
  
    The units of rate constants must be related to the time unit; 1/min and min,
    or 1/sec and sec.
   
    @sa simC1
    @return Function returns 0 when succesful, else a value >= 1.
    @author Vesa Oikonen
 */
int simC1_i(
  /** Array of time values */
  double *t,
  /** Array of AUC 0-t of arterial activities */
  double *cai,
  /** Number of values in TACs */
  const int nr,
  /** Rate constant of the model */
  const double k1,
  /** Rate constant of the model */
  const double k2,
  /** Pointer for TAC array to be simulated; must be allocated */
  double *ct
) {
  int i;
  double dt2;
  double t_last;
  double ct1, ct1_last;
  double ct1i, ct1i_last;


  /* Check for data */
  if(nr<2) return 1;
  if(t==NULL || cai==NULL || ct==NULL) return 2;

  /* Check actual parameter number */
  if(!(k1>=0.0)) return 3;

  /* Calculate curves */
  t_last=0.0; if(t[0]<t_last) t_last=t[0]; 
  ct1_last=ct1i_last=0.0;
  ct1=ct1i=0.0;
  for(i=0; i<nr; i++) {
    /* delta time / 2 */
    dt2=0.5*(t[i]-t_last);
    /* calculate values */
    if(dt2<0.0) {
      return 5;
    } else if(dt2>0.0) {
      /* tissue compartment and its integral */
      ct1 = (k1*cai[i] - k2*(ct1i_last+dt2*ct1_last)) / (1.0 + dt2*k2);
      ct1i = ct1i_last + dt2*(ct1_last+ct1);
    }
    /* copy values to argument arrays; set very small values to zero */
    ct[i]=ct1; if(fabs(ct[i])<1.0e-12) ct[i]=0.0;
    /* prepare to the next loop */
    t_last=t[i]; ct1_last=ct1; ct1i_last=ct1i;
  }

  return 0;
}
/*****************************************************************************/

/*****************************************************************************/
