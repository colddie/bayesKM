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
/** Simulate tissue TAC using two-tissue compartment model and plasma TAC, 
    at plasma TAC times.
     
    @details
    Memory for ct must be allocated in the calling program.
    To retrieve the separate tissue compartment TACs, pointer to allocated
    memory for cta and/or ctb can be given; if compartmental TACs are not
    required, NULL pointer can be given instead.
  
    The units of rate constants must be related to the time unit; 1/min and min,
    or 1/sec and sec.
   
    @return Function returns 0 when succesful, else a value >= 1.
    @author Vesa Oikonen
    @sa simC2_i, simC1, simC3s, simC3p
 */
 
   int simMBF_idl(int argc, float * argv[])

{
  int i;
  double dt2;
  double cii, ci_last, t_last;
  double ct_last, cti, cti_last;

  double *t;
  /** Array of arterial activities */
  double *ci;
  /** Number of values in TACs */
   int nr;
  /** Rate constant of the model */
   double k1;
  /** Rate constant of the model */
   double k2;
  /** Blood volume to simulate spill-over */
   double Vfit;
  /** Pointer for TAC array to be simulated; must be allocated */
  double *ct;
    
  const char *debugfile  = "debug.txt";  


    /* read in parameters */	
    t               =  (double *)  argv[0];
    ci              =  (double *)  argv[1];
    nr              =  *(int*)    argv[2];
    k1              =  *(double*) argv[3];
    k2              =  *(double*) argv[4];
    Vfit            =  *(double*) argv[5]; 
    ct              =  (double *) argv[6];

    /* debug */
    //FILE *pfile = fopen(debugfile, "a+");
    //fprintf(pfile, "NIproj3d: 11 arguments required, %d supplied %f %f %f %f\n", nr,k1,k2,k3,k4);
    //fclose(pfile);


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
/** Simulate tissue TAC using two-tissue compartment model and plasma TAC, 
    at plasma TAC times.
     
    @details
    Memory for ct must be allocated in the calling program.
    To retrieve the separate tissue compartment TACs, pointer to allocated
    memory for cta and/or ctb can be given; if compartmental TACs are not
    required, NULL pointer can be given instead.
  
    The units of rate constants must be related to the time unit; 1/min and min,
    or 1/sec and sec.
     
    This version uses integral of arterial TAC as input function.
    Only advantage over simC2() is that the calculation of integral can be
    fully controlled and possibly more precice in some situations.
   
    @return Function returns 0 when succesful, else a value >= 1.
    @author Vesa Oikonen
    @sa simC2, simC1, simC3s, simC3p
 */

