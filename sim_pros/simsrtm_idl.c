/** @file sim2cm.c
 *  @brief Simulation of 2-tissue compartmental models.
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
/** Simulate tissue TAC using simplified reference tissue input compartment model.
 *   
 *  @details
 *  Memory for ct must be allocated in the calling program.
 *  The units of rate constants must be related to the time unit; 1/min and min,
 *  or 1/sec and sec.
 * 
 *  @return Function returns 0 when successful, else a value >= 1.
 *  @author Vesa Oikonen
 */
int simSRTM_idl(int argc, float * argv[])
{
  int i;
  double dt2;
  double cri, cr_last, t_last;
  double ct_last, cti, cti_last;
  
  /** Array of time values */
  double *t;
  /** Reference region activities */
  double *cr;
  /** Number of values in TACs */
  int nr;
  /** Ratio K1/K1' */
  double R1;
  /** Rate constant of the model */
  double k2;
  /** Binding potential */
  double BP;
  /** Pointer for TAC array to be simulated; must be allocated */
  double *ct;
    
  const char *debugfile  = "debug.txt";  


    /* read in parameters */	
    t               =  (double *)  argv[0];
    cr              =  (double *)  argv[1];
    nr              =  *(int*)    argv[2];
    R1              =  *(double*) argv[3];
    k2              =  *(double*) argv[4];
    BP              =  *(double*) argv[5];
    ct              =  (double *) argv[6];
  
    /* debug */
    //FILE *pfile = fopen(debugfile, "a+");
    //fprintf(pfile, "NIproj3d: 11 arguments required, %d supplied %f %f %f %f\n", nr,k1,k2,k3,k4);
    //fclose(pfile);


 /* Check for data */
  if(nr<2) return 1;
  if(ct==NULL) return 2;

  /* Calculate curves */
  t_last=0.0; if(t[0]<t_last) t_last=t[0];
  cri=cr_last=0.0; cti=ct_last=cti_last=0.0;
  for(i=0; i<nr; i++) {
    /* delta time / 2 */
    dt2=0.5*(t[i]-t_last);
    /* calculate values */
    if(dt2<0.0) {
      return 5;
    } else if(dt2>0.0) {
      /* reference integral */
      cri+=(cr[i]+cr_last)*dt2;
      /* Tissue compartment and its integral */
      ct[i] = ( R1*cr[i] + k2*cri - (k2/(1.0+BP))*(cti_last+dt2*ct_last) ) /
              ( 1.0 + dt2*(k2/(1.0+BP)) );
      cti = cti_last + dt2*(ct_last+ct[i]);
    }
    /* set very small values to zero */
    if(fabs(ct[i])<1.0e-12) ct[i]=0.0;
    /* prepare to the next loop */
    t_last=t[i]; cr_last=cr[i];
    ct_last=ct[i]; cti_last=cti;
  }

  return 0;
}