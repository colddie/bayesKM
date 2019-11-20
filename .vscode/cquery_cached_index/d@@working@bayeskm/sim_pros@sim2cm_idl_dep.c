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
 
extern "C"   int simC2_idl(
  double *t;
  /** Array of arterial activities */
  double *ca;
  /** Number of values in TACs */
   int nr;
  /** Rate constant of the model */
   double k1;
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
)
{
  int i;
  double dt2, r, u, v;
  double cai, ca_last, t_last;
  double ct1, ct1_last, ct2, ct2_last;
  double ct1i, ct1i_last, ct2i, ct2i_last;
    
  const char *debugfile  = "debug.txt";  

  /* debug */
  //FILE *pfile = fopen(debugfile, "a+");
  //fprintf(pfile, "NIproj3d: 11 arguments required, %d supplied %f %f %f %f\n", nr,k1,k2,k3,k4);
  //fclose(pfile);
	
  
  /* Check for data */
  if(nr<2) return 1;
  if(t==NULL || ca==NULL || ct==NULL) return 2;

  /* Check parameters */
  if(k1<0.0) return 3;

  /* Calculate curves */
  t_last=0.0; if(t[0]<t_last) t_last=t[0];
  cai=ca_last=0.0;
  ct1_last=ct2_last=ct1i_last=ct2i_last=0.0;
  ct1=ct2=ct1i=ct2i=0.0;
  for(i=0; i<nr; i++) {
    /* delta time / 2 */
    dt2=0.5*(t[i]-t_last);
    /* calculate values */
    if(dt2<0.0) {
      return 5;
    } else if(dt2>0.0) {
      /* arterial integral */
      cai+=(ca[i]+ca_last)*dt2;
      /* Calculate partial results */
      r=1.0+k4*dt2;
      u=ct1i_last+dt2*ct1_last;
      v=ct2i_last+dt2*ct2_last;
      /* 1st tissue compartment and its integral */
      ct1 = ( k1*cai - (k2 + (k3/r))*u + (k4/r)*v )
            / ( 1.0 + dt2*(k2 + (k3/r)) );
      ct1i = ct1i_last + dt2*(ct1_last+ct1);
      /* 2nd tissue compartment and its integral */
      ct2 = (k3*ct1i - k4*v) / r;
      ct2i = ct2i_last + dt2*(ct2_last+ct2);
    }
    /* copy values to argument arrays; set very small values to zero */
    ct[i]=ct1+ct2; if(fabs(ct[i])<1.0e-12) ct[i]=0.0;
    if(cta!=NULL) {cta[i]=ct1; if(fabs(cta[i])<1.0e-12) cta[i]=0.0;}
    if(ctb!=NULL) {ctb[i]=ct2; if(fabs(ctb[i])<1.0e-12) ctb[i]=0.0;}
    /* prepare to the next loop */
    t_last=t[i]; ca_last=ca[i];
    ct1_last=ct1; ct1i_last=ct1i;
    ct2_last=ct2; ct2i_last=ct2i;
  }

  return 0;
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
extern "C" int simRTCM_idl(
  /** Array of time values */
  double *t,
  /** Reference region activities */
  double *cr,
  /** Number of values in TACs */
  int nr,
  /** Ratio K1/K1' */
  double R1,
  /** Rate constant of the model */
  double k2,
  /** Rate constant of the model */
  double k3,
  /** Rate constant of the model */
  double k4,
  /** Pointer for TAC array to be simulated; must be allocated */
  double *ct,
  /** Pointer for 1st compartment TAC to be simulated, or NULL */
  double *cta,
  /** Pointer for 2nd compartment TAC to be simulated, or NULL */
  double *ctb
)
{
  int i;
  double f, b, w, dt2;
  double cri, cr_last, t_last;
  double cf, cf_last, cb, cb_last;
  double cfi, cfi_last, cbi, cbi_last;
    
  const char *debugfile  = "debug.txt";  
  
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
extern "C" int simSRTM_idl(  /** Array of time values */
  double *t,
  /** Reference region activities */
  double *cr,
  /** Number of values in TACs */
  int nr,
  /** Ratio K1/K1' */
  double R1,
  /** Rate constant of the model */
  double k2,
  /** Binding potential */
  double BP,
  /** Pointer for TAC array to be simulated; must be allocated */
  double *ct
)
{
  int i;
  double dt2;
  double cri, cr_last, t_last;
  double ct_last, cti, cti_last;
  
  const char *debugfile  = "debug.txt";  

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