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
/** Simulate tissue TAC using full reference tissue compartment model (original)
 *  and reference region TAC, at reference region TAC times.
 *   
 *  @details
 *  Memory for ct must be allocated in the calling program.
 *  To retrieve the separate tissue compartment TACs, pointer to allocated
 *  memory for cf and/or cb can be given; if compartmental TACs are not
 *  required, NULL pointer can be given instead.
 *
 *  The units of rate constants must be related to the time unit; 1/min and min,
 *  or 1/sec and sec.
 * 
 *  @return Function returns 0 when successful, else a value >= 1.
 *  @author Vesa Oikonen
 */
int simRTCM_idl(int argc, float * argv[])
{
  int i;
  double f, b, w, dt2;
  double cri, cr_last, t_last;
  double cf, cf_last, cb, cb_last;
  double cfi, cfi_last, cbi, cbi_last;
  
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
  /** Rate constant of the model */
  double k3;
  /** Rate constant of the model */
  double k4;
  /** Pointer for TAC array to be simulated; must be allocated */
  double *ct;
  /** Pointer for 1st compartment TAC to be simulated, or NULL */
  double *cta;
  /** Pointer for 2nd compartment TAC to be simulated, or NULL */
  double *ctb;
    
  const char *debugfile  = "debug.txt";  


    /* read in parameters */	
    t               =  (double *)  argv[0];
    cr              =  (double *)  argv[1];
    nr              =  *(int*)    argv[2];
    R1              =  *(double*) argv[3];
    k2              =  *(double*) argv[4];
    k3              =  *(double*) argv[5];
    k4              =  *(double*) argv[6];
    ct              =  (double *) argv[7];
    cta             =  NULL;
    ctb             =  NULL;
  
    /* debug */
    //FILE *pfile = fopen(debugfile, "a+");
    //fprintf(pfile, "NIproj3d: 11 arguments required, %d supplied %f %f %f %f\n", nr,k1,k2,k3,k4);
    //fclose(pfile);


 /* Check for data */
  if(nr<2) return 1;
  if(ct==NULL) return 2;

  /* Calculate curves */
  t_last=0.0; if(t[0]<t_last) t_last=t[0];
  cri=cr_last=0.0; cf_last=cb_last=cfi_last=cbi_last=cf=cb=cfi=cbi=0.0;
  for(i=0; i<nr; i++) {
    /* delta time / 2 */
    dt2=0.5*(t[i]-t_last);
    /* calculate values */
    if(dt2<0.0) {
      return 5;
    } else if(dt2>0.0) {
      /* reference integral */
      cri+=(cr[i]+cr_last)*dt2;
      /* partial results */
      f=cfi_last+dt2*cf_last;
      b=cbi_last+dt2*cb_last;
      w=k2 + k3 + k2*k4*dt2;
      /* 1st tissue compartment and its integral */
      cf = ( (1.0 + k4*dt2)*(R1*cr[i] + k2*cri) + k4*b - w*f ) /
           ( 1.0 + dt2*(w+k4) );
      cfi = cfi_last + dt2*(cf_last+cf);
      /* 2nd tissue compartment and its integral */
      cb = (k3*cfi - k4*b) / (1.0 + k4*dt2);
      cbi = cbi_last + dt2*(cb_last+cb);
    }
    /* copy values to argument arrays; set very small values to zero */
    ct[i]=cf+cb; if(fabs(ct[i])<1.0e-12) ct[i]=0.0;
    if(cta!=NULL) {cta[i]=cf; if(fabs(cta[i])<1.0e-12) cta[i]=0.0;}
    if(ctb!=NULL) {ctb[i]=cb; if(fabs(ctb[i])<1.0e-12) ctb[i]=0.0;}
    /* prepare to the next loop */
    t_last=t[i]; cr_last=cr[i];
    cf_last=cf; cfi_last=cfi;
    cb_last=cb; cbi_last=cbi;
  }

  return 0;
}

