/// @file simulate.c
/// @brief Procedures for simulating PET time-activity curves.
/// @author Vesa Oikonen
/// @todo Remove setting of small values to zero, and add instead a separate
///       function to do that for whole DFT with user-defined limit.
///
/*****************************************************************************/
#include "libtpcmodel.h"
/*****************************************************************************/

/*****************************************************************************/
/** Simulates tissue TAC using 1-3 tissue compartment model (in series) and
    plasma TAC, at plasma TAC times.

    Memory for ct must be allocated in the calling program.
    To retrieve the separate tissue compartment TACs, pointer to allocated
    memory for cta, ctb and/or ctc can be given; if compartmental TACs are not
    required, NULL pointer can be given instead.

    The units of rate constants must be related to the time unit; 1/min and min,
    or 1/sec and sec.

\return Function returns 0 when successful, else a value >= 1.
*/

int SIMULATE_TEST;
int simC3s(
  /** Array of time values */
  double *t,
  /** Array of arterial activities */
  double *ca,
  /** Number of values in TACs */
  int nr,
  /** Rate constant of the model */
  double k1,
  /** Rate constant of the model */
  double k2,
  /** Rate constant of the model */
  double k3,
  /** Rate constant of the model */
  double k4,
  /** Rate constant of the model */
  double k5,
  /** Rate constant of the model */
  double k6,
  /** Pointer for TAC array to be simulated; must be allocated */
  double *ct,
  /** Pointer for 1st compartment TAC to be simulated, or NULL */
  double *cta,
  /** Pointer for 2nd compartment TAC to be simulated, or NULL */
  double *ctb,
  /** Pointer for 3rd compartment TAC to be simulated, or NULL */
  double *ctc
) {
  int i, model, parNr;
  double b, c, d, w, z, dt2;
  double cai, ca_last, t_last;
  double ct1, ct1_last, ct2, ct2_last, ct3, ct3_last;
  double ct1i, ct1i_last, ct2i, ct2i_last, ct3i, ct3i_last;


  /* Check for data */
  if(nr<2) return 1;
  if(ct==NULL) return 2;

  /* Check actual parameter number */
  if(k1<0.0) return 3;
  if(k3<=0.0) {k3=0.0; model=1; if(k2<=0.0) {k2=0.0; parNr=1;} else parNr=2;}
  else if(k5<=0.0) {k5=0.0; model=2; if(k4<=0.0) {k4=0.0;parNr=3;} else parNr=4;}
  else {model=3; if(k6<=0.0) {k6=0.0; parNr=5;} else parNr=6;}
  if(SIMULATE_TEST) printf("simulate(): model=%d parNr=%d\n", model, parNr);

  /* Calculate curves */
  t_last=0.0; if(t[0]<t_last) t_last=t[0]; 
  cai=ca_last=0.0;
  ct1_last=ct2_last=ct3_last=ct1i_last=ct2i_last=ct3i_last=0.0;
  ct1=ct2=ct3=ct1i=ct2i=ct3i=0.0;
  for(i=0; i<nr; i++) {
    /* delta time / 2 */
    dt2=0.5*(t[i]-t_last);
    /* calculate values */
    if(dt2<0.0) {
      return 5;
    } else if(dt2>0.0) {
      /* arterial integral */
      cai+=(ca[i]+ca_last)*dt2;
      /* partial results */
      b=ct1i_last+dt2*ct1_last;
      c=ct2i_last+dt2*ct2_last;
      d=ct3i_last+dt2*ct3_last;
      w=k4 + k5 - (k5*k6*dt2)/(1.0+k6*dt2);
      z=1.0+w*dt2;
      /* 1st tissue compartment and its integral */
      ct1 = (
          + k1*z*cai + (k3*k4*dt2 - (k2+k3)*z)*b
          + k4*c + k4*k6*dt2*d/(1.0+k6*dt2)
        ) / ( z*(1.0 + dt2*(k2+k3)) - k3*k4*dt2*dt2 );
      ct1i = ct1i_last + dt2*(ct1_last+ct1);
      /* 2nd tissue compartment and its integral */
      ct2 = (k3*ct1i - w*c + k6*d/(1.0+k6*dt2)) / z;
      ct2i = ct2i_last + dt2*(ct2_last+ct2);
      /* 3rd tissue compartment and its integral */
      ct3 = (k5*ct2i - k6*d) / (1.0 + k6*dt2);
      ct3i = ct3i_last + dt2*(ct3_last+ct3);
    }
    /* copy values to argument arrays; set very small values to zero */
    ct[i]=ct1+ct2+ct3; if(fabs(ct[i])<1.0e-12) ct[i]=0.0;
    if(cta!=NULL) {cta[i]=ct1; if(fabs(cta[i])<1.0e-12) cta[i]=0.0;}
    if(ctb!=NULL) {ctb[i]=ct2; if(fabs(ctb[i])<1.0e-12) ctb[i]=0.0;}
    if(ctc!=NULL) {ctc[i]=ct3; if(fabs(ctc[i])<1.0e-12) ctc[i]=0.0;}
    /* prepare to the next loop */
    t_last=t[i]; ca_last=ca[i];
    ct1_last=ct1; ct1i_last=ct1i;
    ct2_last=ct2; ct2i_last=ct2i;
    ct3_last=ct3; ct3i_last=ct3i;
  }

  return 0;
}
/*****************************************************************************/

/*****************************************************************************/
/** Simulates tissue TAC using 1-3 tissue compartment model (2nd and 3rd
    compartments in parallel) and plasma TAC, at plasma TAC times.

    Memory for ct must be allocated in the calling program.
    To retrieve the separate tissue compartment TACs, pointer to allocated
    memory for cta, ctb and/or ctc can be given; if compartmental TACs are not
    required, NULL pointer can be given instead.

    The units of rate constants must be related to the time unit; 1/min and min,
    or 1/sec and sec.

\return Function returns 0 when successful, else a value >= 1.
*/
int simC3p(
  /** Array of time values */
  double *t,
  /** Array of arterial activities */
  double *ca,
  /** Number of values in TACs */
  int nr,
  /** Rate constant of the model */
  double k1,
  /** Rate constant of the model */
  double k2,
  /** Rate constant of the model */
  double k3,
  /** Rate constant of the model */
  double k4,
  /** Rate constant of the model */
  double k5,
  /** Rate constant of the model */
  double k6,
  /** Pointer for TAC array to be simulated; must be allocated */
  double *ct,
  /** Pointer for 1st compartment TAC to be simulated, or NULL */
  double *cta,
  /** Pointer for 2nd compartment TAC to be simulated, or NULL */
  double *ctb,
  /** Pointer for 3rd compartment TAC to be simulated, or NULL */
  double *ctc
) {
  int i;
  double dt2, r, s, u, v, w;
  double cai, ca_last, t_last;
  double ct1, ct1_last, ct2, ct2_last, ct3, ct3_last;
  double ct1i, ct1i_last, ct2i, ct2i_last, ct3i, ct3i_last;


  /* Check for data */
  if(nr<2) return 1;
  if(ct==NULL) return 2;

  /* Check parameters */
  if(k1<0.0) return 3;

  /* Calculate curves */
  t_last=0.0; if(t[0]<t_last) t_last=t[0];
  cai=ca_last=0.0;
  ct1_last=ct2_last=ct3_last=ct1i_last=ct2i_last=ct3i_last=0.0;
  ct1=ct2=ct3=ct1i=ct2i=ct3i=0.0;
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
      s=1.0+k6*dt2;
      u=ct1i_last+dt2*ct1_last;
      v=ct2i_last+dt2*ct2_last;
      w=ct3i_last+dt2*ct3_last;
      /* 1st tissue compartment and its integral */
      ct1 = ( k1*cai - (k2 + (k3/r) + (k5/s))*u + (k4/r)*v + (k6/s)*w )
            / ( 1.0 + dt2*(k2 + (k3/r) + (k5/s)) );
      ct1i = ct1i_last + dt2*(ct1_last+ct1);
      /* 2nd tissue compartment and its integral */
      ct2 = (k3*ct1i - k4*v) / r;
      ct2i = ct2i_last + dt2*(ct2_last+ct2);
      /* 3rd tissue compartment and its integral */
      ct3 = (k5*ct1i - k6*w) / s;
      ct3i = ct3i_last + dt2*(ct3_last+ct3);
    }
    /* copy values to argument arrays; set very small values to zero */
    ct[i]=ct1+ct2+ct3; if(fabs(ct[i])<1.0e-12) ct[i]=0.0;
    if(cta!=NULL) {cta[i]=ct1; if(fabs(cta[i])<1.0e-12) cta[i]=0.0;}
    if(ctb!=NULL) {ctb[i]=ct2; if(fabs(ctb[i])<1.0e-12) ctb[i]=0.0;}
    if(ctc!=NULL) {ctc[i]=ct3; if(fabs(ctc[i])<1.0e-12) ctc[i]=0.0;}
    /* prepare to the next loop */
    t_last=t[i]; ca_last=ca[i];
    ct1_last=ct1; ct1i_last=ct1i;
    ct2_last=ct2; ct2i_last=ct2i;
    ct3_last=ct3; ct3i_last=ct3i;
  }

  return 0;
}
/*****************************************************************************/

/*****************************************************************************/
/** Simulates tissue TAC using 1-3 tissue compartment model (in series) and
    plasma TAC, at plasma TAC times, considering also arterial and venous
    vasculature.

    Memory for cpet must be allocated in the calling program.
    To retrieve the separate tissue compartment TACs, pointer to allocated
    memory for cta, ctb, ctc, ctab and/or ctvb can be given; if compartmental
    TACs are not required, NULL pointer can be given instead.

    The units of rate constants must be related to the time unit; 1/min and min,
    or 1/sec and sec.

    If blood flow is set to 0, function assumes that f>>k1, and Cvb=Cab.,

\return Function returns 0 when successful, else a value >= 1.
*/
int simC3vs(
  /** Array of time values */
  double *t,
  /** Array of arterial plasma activities */
  double *ca,
  /** Array of arterial blood activities */
  double *cb,
  /** Number of values in TACs */
  int nr,
  /** Rate constant of the model */
  double k1,
  /** Rate constant of the model */
  double k2,
  /** Rate constant of the model */
  double k3,
  /** Rate constant of the model */
  double k4,
  /** Rate constant of the model */
  double k5,
  /** Rate constant of the model */
  double k6,
  /** Blood flow; if 0, function assumes that f>>k1, and Cvb=Cab. */
  double f,
  /** Vascular volume fraction */
  double vb,
  /** Arterial fraction of vascular volume */
  double fa,
  /** Pointer for TAC array to be simulated; must be allocated */
  double *cpet,
  /** Pointer for 1st compartment TAC to be simulated, or NULL */
  double *cta,
  /** Pointer for 2nd compartment TAC to be simulated, or NULL */
  double *ctb,
  /** Pointer for 3rd compartment TAC to be simulated, or NULL */
  double *ctc,
  /** Pointer for arterial TAC in tissue, or NULL */
  double *ctab,
  /** Pointer for venous TAC in tissue, or NULL */
  double *ctvb
) {
  int i;
  double b, c, d, w, z, dt2, va, vv;
  double cai, ca_last, t_last, dct, cvb;
  double ct1, ct1_last, ct2, ct2_last, ct3, ct3_last;
  double ct1i, ct1i_last, ct2i, ct2i_last, ct3i, ct3i_last;


  /* Check for data */
  if(nr<2) return 1;
  if(cpet==NULL) return 2;

  /* Check parameters */
  if(k1<0.0) return 3;
  if(vb<0.0 || vb>=1.0) return 4;
  if(fa<=0.0 || fa>1.0) return 5;
  va=fa*vb; vv=(1.0-fa)*vb;

  /* Calculate curves */
  t_last=0.0; if(t[0]<t_last) t_last=t[0];
  cai=ca_last=0.0;
  ct1_last=ct2_last=ct3_last=ct1i_last=ct2i_last=ct3i_last=0.0;
  ct1=ct2=ct3=ct1i=ct2i=ct3i=0.0;
  for(i=0; i<nr; i++) {
    /* delta time / 2 */
    dt2=0.5*(t[i]-t_last);
    /* calculate values */
    if(dt2<0.0) {
      return 5;
    } else if(dt2>0.0) {
      /* arterial integral */
      cai+=(ca[i]+ca_last)*dt2;
      /* partial results */
      b=ct1i_last+dt2*ct1_last;
      c=ct2i_last+dt2*ct2_last;
      d=ct3i_last+dt2*ct3_last;
      w=k4 + k5 - (k5*k6*dt2)/(1.0+k6*dt2);
      z=1.0+w*dt2;
      /* 1st tissue compartment and its integral */
      ct1 = (
          + k1*z*cai + (k3*k4*dt2 - (k2+k3)*z)*b
          + k4*c + k4*k6*dt2*d/(1.0+k6*dt2)
        ) / ( z*(1.0 + dt2*(k2+k3)) - k3*k4*dt2*dt2 );
      ct1i = ct1i_last + dt2*(ct1_last+ct1);
      /* 2nd tissue compartment and its integral */
      ct2 = (k3*ct1i - w*c + k6*d/(1.0+k6*dt2)) / z;
      ct2i = ct2i_last + dt2*(ct2_last+ct2);
      /* 3rd tissue compartment and its integral */
      ct3 = (k5*ct2i - k6*d) / (1.0 + k6*dt2);
      ct3i = ct3i_last + dt2*(ct3_last+ct3);
    }
    /* Venous curve */
    if(f>0.) {dct = k1*ca[i] - k2*ct1; cvb = cb[i] - dct/f;} else cvb=cb[i];
    /* copy values to argument arrays; set very small values to zero */
    cpet[i]= va*cb[i] + vv*cvb + (1.0-vb)*(ct1+ct2+ct3);
    if(fabs(cpet[i])<1.0e-12) cpet[i]=0.0;
    if(cta!=NULL) {cta[i]=(1.0-vb)*ct1; if(fabs(cta[i])<1.0e-12) cta[i]=0.0;}
    if(ctb!=NULL) {ctb[i]=(1.0-vb)*ct2; if(fabs(ctb[i])<1.0e-12) ctb[i]=0.0;}
    if(ctc!=NULL) {ctc[i]=(1.0-vb)*ct3; if(fabs(ctc[i])<1.0e-12) ctc[i]=0.0;}
    if(ctab!=NULL) {ctab[i]=va*cb[i]; if(fabs(ctab[i])<1.0e-12) ctab[i]=0.0;}
    if(ctvb!=NULL) {ctvb[i]=vv*cvb; if(fabs(ctvb[i])<1.0e-12) ctvb[i]=0.0;}
    /* prepare to the next loop */
    t_last=t[i]; ca_last=ca[i];
    ct1_last=ct1; ct1i_last=ct1i;
    ct2_last=ct2; ct2i_last=ct2i;
    ct3_last=ct3; ct3i_last=ct3i;
  }

  return 0;
}
/*****************************************************************************/

/*****************************************************************************/
/** Simulates tissue TAC using 1-3 tissue compartment model (2nd and 3rd
    compartments in parallel) and plasma TAC, at plasma TAC times,
    considering also arterial and venous vasculature.

    Memory for cpet must be allocated in the calling program.
    To retrieve the separate tissue compartment TACs, pointer to allocated
    memory for cta, ctb, ctc, ctab and/or ctvb can be given; if compartmental
    TACs are not required, NULL pointer can be given instead.

    The units of rate constants must be related to the time unit; 1/min and min,
    or 1/sec and sec.

    If blood flow is set to 0, function assumes that f>>k1, and Cvb=Cab.,

\return Function returns 0 when successful, else a value >= 1.
*/
int simC3vp(
  /** Array of time values */
  double *t,
  /** Array of arterial plasma activities */
  double *ca,
  /** Array of arterial blood activities */
  double *cb,
  /** Number of values in TACs */
  int nr,
  /** Rate constant of the model */
  double k1,
  /** Rate constant of the model */
  double k2,
  /** Rate constant of the model */
  double k3,
  /** Rate constant of the model */
  double k4,
  /** Rate constant of the model */
  double k5,
  /** Rate constant of the model */
  double k6,
  /** Blood flow; if 0, function assumes that f>>k1, and Cvb=Cab. */
  double f,
  /** Vascular volume fraction */
  double vb,
  /** Arterial fraction of vascular volume */
  double fa,
  /** Pointer for TAC array to be simulated; must be allocated */
  double *cpet,
  /** Pointer for 1st compartment TAC to be simulated, or NULL */
  double *cta,
  /** Pointer for 2nd compartment TAC to be simulated, or NULL */
  double *ctb,
  /** Pointer for 3rd compartment TAC to be simulated, or NULL */
  double *ctc,
  /** Pointer for arterial TAC in tissue, or NULL */
  double *ctab,
  /** Pointer for venous TAC in tissue, or NULL */
  double *ctvb
) {
  int i;
  double dt2, r, s, u, v, w, va, vv;
  double cai, ca_last, t_last, dct, cvb;
  double ct1, ct1_last, ct2, ct2_last, ct3, ct3_last;
  double ct1i, ct1i_last, ct2i, ct2i_last, ct3i, ct3i_last;


  /* Check for data */
  if(nr<2) return 1;
  if(cpet==NULL) return 2;

  /* Check parameters */
  if(k1<0.0) return 3;
  if(vb<0.0 || vb>=1.0) return 4;
  if(fa<=0.0 || fa>1.0) return 5;
  va=fa*vb; vv=(1.0-fa)*vb;

  /* Calculate curves */
  t_last=0.0; if(t[0]<t_last) t_last=t[0];
  cai=ca_last=0.0;
  ct1_last=ct2_last=ct3_last=ct1i_last=ct2i_last=ct3i_last=0.0;
  ct1=ct2=ct3=ct1i=ct2i=ct3i=0.0;
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
      s=1.0+k6*dt2;
      u=ct1i_last+dt2*ct1_last;
      v=ct2i_last+dt2*ct2_last;
      w=ct3i_last+dt2*ct3_last;
      /* 1st tissue compartment and its integral */
      ct1 = ( k1*cai - (k2 + (k3/r) + (k5/s))*u + (k4/r)*v + (k6/s)*w )
            / ( 1.0 + dt2*(k2 + (k3/r) + (k5/s)) );
      ct1i = ct1i_last + dt2*(ct1_last+ct1);
      /* 2nd tissue compartment and its integral */
      ct2 = (k3*ct1i - k4*v) / r;
      ct2i = ct2i_last + dt2*(ct2_last+ct2);
      /* 3rd tissue compartment and its integral */
      ct3 = (k5*ct1i - k6*w) / s;
      ct3i = ct3i_last + dt2*(ct3_last+ct3);
    }
    /* Venous curve */
    if(f>0.) {dct = k1*ca[i] - k2*ct1; cvb = cb[i] - dct/f;} else cvb=cb[i];
    /* copy values to argument arrays; set very small values to zero */
    cpet[i]= va*cb[i] + vv*cvb + (1.0-vb)*(ct1+ct2+ct3);
    if(fabs(cpet[i])<1.0e-12) cpet[i]=0.0;
    if(cta!=NULL) {cta[i]=(1.0-vb)*ct1; if(fabs(cta[i])<1.0e-12) cta[i]=0.0;}
    if(ctb!=NULL) {ctb[i]=(1.0-vb)*ct2; if(fabs(ctb[i])<1.0e-12) ctb[i]=0.0;}
    if(ctc!=NULL) {ctc[i]=(1.0-vb)*ct3; if(fabs(ctc[i])<1.0e-12) ctc[i]=0.0;}
    if(ctab!=NULL) {ctab[i]=va*cb[i]; if(fabs(ctab[i])<1.0e-12) ctab[i]=0.0;}
    if(ctvb!=NULL) {ctvb[i]=vv*cvb; if(fabs(ctvb[i])<1.0e-12) ctvb[i]=0.0;}
    /* prepare to the next loop */
    t_last=t[i]; ca_last=ca[i];
    ct1_last=ct1; ct1i_last=ct1i;
    ct2_last=ct2; ct2i_last=ct2i;
    ct3_last=ct3; ct3i_last=ct3i;
  }

  return 0;
}
/*****************************************************************************/

/*****************************************************************************/
/** Simulates tissue TAC using 2 tissue compartment model (in series) and
    plasma TAC, at plasma TAC times.
    In contrary to the common model, kLoss represents a direct loss rate from
    the 2nd tissue compartment to venous plasma.

    Memory for ct must be allocated in the calling program.
    To retrieve the separate tissue compartment TACs, pointer to allocated
    memory for cta and ctb can be given; if compartmental TACs are not
    required, NULL pointer can be given instead.

    The units of rate constants must be related to the time unit; 1/min and min,
    or 1/sec and sec.

\return Function returns 0 when successful, else a value >= 1.
*/
int simC2l(
  /** Array of time values */
  double *t,
  /** Array of arterial activities */
  double *ca,
  /** Number of values in TACs */
  int nr,
  /** Rate constant of the model */
  double k1,
  /** Rate constant of the model */
  double k2,
  /** Rate constant of the model */
  double k3,
  /** Rate constant of the model */
  double kLoss,
  /** Pointer for TAC array to be simulated; must be allocated */
  double *ct,
  /** Pointer for 1st compartment TAC to be simulated, or NULL */
  double *cta,
  /** Pointer for 2nd compartment TAC to be simulated, or NULL */
  double *ctb
) {
  int i;
  double b, c, dt2;
  double cai, ca_last, t_last;
  double ct1, ct1_last, ct2, ct2_last;
  double ct1i, ct1i_last, ct2i, ct2i_last;


  /* Check for data */
  if(nr<2) return 1;
  if(ct==NULL) return 2;

  /* Check actual parameter number */
  if(k1<0.0) return 3;

  /* Calculate curves */
  t_last=0.0; if(t[0]<t_last) t_last=t[0];
  cai=ca_last=0.0;
  ct1_last=ct2_last=ct1i_last=ct2i_last=ct1=ct2=ct1i=ct2i=0.0;
  for(i=0; i<nr; i++) {
    /* delta time / 2 */
    dt2=0.5*(t[i]-t_last);
    /* calculate values */
    if(dt2<0.0) {
      return 5;
    } else if(dt2>0.0) {
      /* arterial integral */
      cai+=(ca[i]+ca_last)*dt2;
      /* partial results */
      b=ct1i_last+dt2*ct1_last;
      c=ct2i_last+dt2*ct2_last;
      /* 1st tissue compartment and its integral */
      ct1 = (k1*cai - (k2+k3)*b) / (1.0 + (k2+k3)*dt2 );
      ct1i = ct1i_last + dt2*(ct1_last+ct1);
      /* 2nd tissue compartment and its integral */
      ct2 = (k3*ct1i - kLoss*c) / (1.0 + kLoss*dt2);
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
/** Simulates tissue TAC using 2 tissue compartment model and plasma TAC,
    at plasma TAC times, considering also arterial and venous vasculature.
    The efflux from 2nd tissue compartment (at rate kL) goes directly to blood.

    Memory for cpet must be allocated in the calling program.
    To retrieve the separate tissue compartment TACs, pointer to allocated
    memory for cta, ctb, ctab and/or ctvb can be given; if compartmental
    TACs are not required, NULL pointer can be given instead.

    The units of rate constants must be related to the time unit; 1/min and min,
    or 1/sec and sec.

    If blood flow is set to 0, function assumes that f>>k1, and Cvb=Cab.,

\return Function returns 0 when successful, else a value >= 1.
*/
int simC2vl(
  /** Array of time values */
  double *t,
  /** Array of arterial plasma activities */
  double *ca,
  /** Array of arterial blood activities */
  double *cb,
  /** Number of values in TACs */
  int nr,
  /** Rate constant of the model */
  double k1,
  /** Rate constant of the model */
  double k2,
  /** Rate constant of the model */
  double k3,
  /** Rate constant of the model */
  double kL,
  /** Blood flow; if 0, function assumes that f>>k1, and Cvb=Cab. */
  double f,
  /** Vascular volume fraction */
  double vb,
  /** Arterial fraction of vascular volume */
  double fa,
  /** Pointer for TAC array to be simulated; must be allocated */
  double *cpet,
  /** Pointer for 1st compartment TAC to be simulated, or NULL */
  double *cta,
  /** Pointer for 2nd compartment TAC to be simulated, or NULL */
  double *ctb,
  /** Pointer for arterial TAC in tissue, or NULL */
  double *ctab,
  /** Pointer for venous TAC in tissue, or NULL */
  double *ctvb
) {
  int i;
  double dt2, b, c, va, vv;
  double cai, ca_last, t_last, dct, cvb;
  double ct1, ct1_last, ct2, ct2_last;
  double ct1i, ct1i_last, ct2i, ct2i_last;


  /* Check for data */
  if(nr<2) return 1;
  if(cpet==NULL) return 2;

  /* Check parameters */
  if(k1<0.0) return 3;
  if(vb<0.0 || vb>=1.0) return 4;
  if(fa<=0.0 || fa>1.0) return 5;
  va=fa*vb; vv=(1.0-fa)*vb;

  /* Calculate curves */
  t_last=0.0; if(t[0]<t_last) t_last=t[0];
  cai=ca_last=0.0;
  ct1_last=ct2_last=ct1i_last=ct2i_last=ct1=ct2=ct1i=ct2i=0.0;
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
      b=ct1i_last+dt2*ct1_last;
      c=ct2i_last+dt2*ct2_last;
      /* 1st tissue compartment and its integral */
      ct1 = (k1*cai - (k2+k3)*b) / (1.0 + (k2+k3)*dt2 );
      ct1i = ct1i_last + dt2*(ct1_last+ct1);
      /* 2nd tissue compartment and its integral */
      ct2 = (k3*ct1i - kL*c) / (1.0 + kL*dt2);
      ct2i = ct2i_last + dt2*(ct2_last+ct2);
    }
    /* Venous curve */
    if(f>0.) {dct = k1*ca[i] - k2*ct1 - kL*ct2; cvb = cb[i] - dct/f;}
    else cvb=cb[i];
    /* copy values to argument arrays; set very small values to zero */
    cpet[i]= va*cb[i] + vv*cvb + (1.0-vb)*(ct1+ct2);
    if(fabs(cpet[i])<1.0e-12) cpet[i]=0.0;
    if(cta!=NULL) {cta[i]=(1.0-vb)*ct1; if(fabs(cta[i])<1.0e-12) cta[i]=0.0;}
    if(ctb!=NULL) {ctb[i]=(1.0-vb)*ct2; if(fabs(ctb[i])<1.0e-12) ctb[i]=0.0;}
    if(ctab!=NULL) {ctab[i]=va*cb[i]; if(fabs(ctab[i])<1.0e-12) ctab[i]=0.0;}
    if(ctvb!=NULL) {ctvb[i]=vv*cvb; if(fabs(ctvb[i])<1.0e-12) ctvb[i]=0.0;}
    /* prepare to the next loop */
    t_last=t[i]; ca_last=ca[i];
    ct1_last=ct1; ct1i_last=ct1i;
    ct2_last=ct2; ct2i_last=ct2i;
  }

  return 0;
}
/*****************************************************************************/

/*****************************************************************************/
/** Simulates tissue TAC using 3 tissue compartmental model with two parallel
    compartments, and plasma TAC, at plasma TAC sample times, considering also
    arterial and venous vasculature.
    The efflux from 3rd tissue compartment (C) goes directly to blood
    at rate kLoss.

    Memory for cpet must be allocated in the calling program.
    To retrieve the separate tissue compartment TACs, pointer to allocated
    memory for cta, ctb, ctab and/or ctvb can be given; if compartmental
    TACs are not required, NULL pointer can be given instead.

    The units of rate constants must be related to the time unit; 1/min and min,
    or 1/sec and sec.

    If blood flow is set to 0, function assumes that f>>k1, and Cvb=Cab.,

\return Function returns 0 when successful, else a value >= 1.
*/
int simC3vpKLoss(
  /** Array of sample times */
  double *t,
  /** Array of arterial plasma activities */
  double *ca,
  /** Array of arterial blood activities */
  double *cb,
  /** Number of sample values in TACs */
  int nr,
  /** Rate constant of the model */
  double k1,
  /** Rate constant of the model */
  double k2,
  /** Rate constant of the model */
  double k3,
  /** Rate constant of the model */
  double k4,
  /** Rate constant of the model */
  double k5,
  /** Rate constant of the model */
  double k6,
  /** Rate constant of the model */
  double kLoss,
  /** Blood flow; if 0, function assumes that f>>k1, and Cvb=Cab. */
  double f,
  /** Vascular volume fraction */
  double vb,
  /** Arterial fraction of vascular volume */
  double fa,
  /** Pointer for TAC array to be simulated; must be allocated */
  double *cpet,
  /** Pointer for 1st tissue compartment TAC to be simulated, or NULL */
  double *cta,
  /** Pointer for 2nd compartment TAC to be simulated, or NULL */
  double *ctb,
  /** Pointer for 3rd compartment TAC to be simulated, or NULL */
  double *ctc,
  /** Pointer for arterial TAC in tissue, or NULL */
  double *ctab,
  /** Pointer for venous TAC in tissue, or NULL */
  double *ctvb
) {
  int i;
  double dt2, b, c, d, w, z, u, va, vv;
  double cai, ca_last, t_last, dct, cvb;
  double ct1, ct1_last, ct2, ct2_last, ct3, ct3_last;
  double ct1i, ct1i_last, ct2i, ct2i_last, ct3i, ct3i_last;
 
 
  /* Check for data */
  if(nr<2) return 1;
  if(cpet==NULL) return 2;
 
  /* Check parameters */
  if(k1<0.0) return 3;
  if(vb<0.0 || vb>=1.0) return 4;
  if(fa<=0.0 || fa>1.0) return 5;
  va=fa*vb; vv=(1.0-fa)*vb;
 
  /* Calculate curves */
  t_last=0.0; if(t[0]<t_last) t_last=t[0];
  cai=ca_last=0.0;
  ct1_last=ct2_last=ct3_last=ct1i_last=ct2i_last=ct3i_last=0.0;
  ct1=ct2=ct3=ct1i=ct2i=ct3i=0.0;
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
      w = 1.0 + k4*dt2;
      z = 1.0 + dt2*(k6 + kLoss); 
      u = k2 + k3 + k5 - k3*k4*dt2/w - k5*k6*dt2/z;
      b = ct1i_last+dt2*ct1_last;
      c = ct2i_last+dt2*ct2_last;
      d = ct3i_last+dt2*ct3_last;
      /* 1st tissue compartment and its integral */
      ct1 = ( k1*cai - u*b + k4*c/w + k6*d/z ) / ( 1.0 + dt2*u);
      ct1i = ct1i_last + dt2*(ct1_last+ct1);
      /* 2nd tissue compartment and its integral */
      ct2 = (k3*ct1i - k4*c) / w;
      ct2i = ct2i_last + dt2*(ct2_last+ct2);
      /* 3rd tissue compartment and its integral */
      ct3 = (k5*ct1i - (k6 + kLoss)*d) / z;
      ct3i = ct3i_last + dt2*(ct3_last+ct3);
    }
    /* Venous curve */
    if(f>0.) {dct = k1*ca[i] - k2*ct1 - kLoss*ct3; cvb = cb[i] - dct/f;}
    else cvb=cb[i];
    /* copy values to argument arrays; set very small values to zero */
    cpet[i]= va*cb[i] + vv*cvb + (1.0-vb)*(ct1+ct2+ct3);
    if(fabs(cpet[i])<1.0e-12) cpet[i]=0.0;
    if(cta!=NULL) {cta[i]=(1.0-vb)*ct1; if(fabs(cta[i])<1.0e-12) cta[i]=0.0;}
    if(ctb!=NULL) {ctb[i]=(1.0-vb)*ct2; if(fabs(ctb[i])<1.0e-12) ctb[i]=0.0;}
    if(ctc!=NULL) {ctc[i]=(1.0-vb)*ct3; if(fabs(ctc[i])<1.0e-12) ctc[i]=0.0;}
    if(ctab!=NULL) {ctab[i]=va*cb[i]; if(fabs(ctab[i])<1.0e-12) ctab[i]=0.0;}
    if(ctvb!=NULL) {ctvb[i]=vv*cvb; if(fabs(ctvb[i])<1.0e-12) ctvb[i]=0.0;}
    /* prepare to the next loop */
    t_last=t[i]; ca_last=ca[i];
    ct1_last=ct1; ct1i_last=ct1i;
    ct2_last=ct2; ct2i_last=ct2i;
    ct3_last=ct3; ct3i_last=ct3i;
  }
  
  return 0;
}
/*****************************************************************************/

/*****************************************************************************/
/** Simulates tissue TAC using reference tissue compartment model (original) and
    reference region TAC, at reference region TAC times.

    Memory for ct must be allocated in the calling program.
    To retrieve the separate tissue compartment TACs, pointer to allocated
    memory for cf and/or cb can be given; if compartmental TACs are not
    required, NULL pointer can be given instead.

    The units of rate constants must be related to the time unit; 1/min and min,
    or 1/sec and sec.

\return Function returns 0 when successful, else a value >= 1.
*/
int simRTCM(
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
) {
  int i;
  double f, b, w, dt2;
  double cri, cr_last, t_last;
  double cf, cf_last, cb, cb_last;
  double cfi, cfi_last, cbi, cbi_last;


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

/*****************************************************************************/
/** Simulates tissue TAC using reference tissue compartment model (simplified)
    and reference region TAC, at reference region TAC times.

    Memory for ct must be allocated in the calling program.

    The units of rate constants must be related to the time unit; 1/min and min,
    or 1/sec and sec.

\return Function returns 0 when successful, else a value >= 1.
*/
int simSRTM(
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
  /** Binding potential */
  double BP,
  /** Pointer for TAC array to be simulated; must be allocated */
  double *ct
) {
  int i;
  double dt2;
  double cri, cr_last, t_last;
  double ct_last, cti, cti_last;


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
    /* set invalid or very small values to zero */
    if(!(fabs(ct[i])>=1.0e-12)) ct[i]=0.0;
    /* prepare to the next loop */
    t_last=t[i]; cr_last=cr[i];
    ct_last=ct[i]; cti_last=cti;
  }

  return 0;
}
/*****************************************************************************/

/*****************************************************************************/
/** Simulates tissue TAC using reference tissue compartment model (transport
    limited in ref region) and reference region TAC, at reference region TAC
    times.

    Memory for ct must be allocated in the calling program.

    The units of rate constants must be related to the time unit; 1/min and min,
    or 1/sec and sec.

\return Function returns 0 when successful, else a value >= 1.
*/
int simTRTM(
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
  /** Pointer for TAC array to be simulated; must be allocated */
  double *ct
) {
  int i;
  double dt2;
  double cri, cr_last, t_last;
  double ct_last, cti, cti_last;


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
      ct[i] = ( R1*cr[i] + R1*k3*cri - (k2+k3)*(cti_last+dt2*ct_last) ) /
              ( 1.0 + dt2*(k2+k3) );
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
/*****************************************************************************/

/*****************************************************************************/
/** Simulation of TACs of parent tracer, and 1-2 of its metabolites in plasma
    using Huang's compartmental model.

    The units of model parameters must be related to the sample time unit;
    1/min and min, or 1/sec and sec.

    Pointers to memory for output TACs must be specified, or NULL if TAC is not
    needed.

\return Returns 0, if ok.
*/
int simHuangmet(
  /** Input: Sample times (preferably with short intervals)  */
  double *t,
  /** Input: Measured total plasma TAC */
  double *ctot,
  /** Input: Nr of samples */
  int nr,
  /** Input: Model parameters */
  double k01,
  /** Input: Model parameters */
  double k12,
  /** Input: Model parameters */
  double k21,
  /** Input: Model parameters */
  double k03,
  /** Input: Model parameters */
  double k34,
  /** Input: Model parameters */
  double k43,
  /** Output: unchanged (parent) tracer TAC */
  double *c0,
  /** Output: TAC of the 1st metabolite */
  double *c1,
  /** Output: TAC of the 2nd metabolite */
  double *c3
) {
  int i;
  double dt2, t_last, ictot, ctot_last;
  double c0_, c1_, c2_, c3_, c4_;
  double ic0_, ic1_, ic2_, ic3_, ic4_;
  double c0_last, c1_last, c2_last, c3_last, c4_last;
  double ic0_last, ic1_last, ic2_last, ic3_last, ic4_last;
  double a, b, am1, am2, am3, am4;


  /* Check input */
  if(t==NULL || ctot==NULL || nr<2) return(1);
  if(k01<0 || k12<0 || k21<0 || k03<0 || k34<0 || k43<0) return(2);
  if(t[0]<0.0) return(3);

  /* Compute the TACs */
  t_last=0.0; if(t[0]<t_last) t_last=t[0];
  ictot=ctot_last=0.0;
  c0_=c0_last=ic0_=ic0_last=0.0;
  c1_=c1_last=ic1_=ic1_last=0.0;
  c2_=c2_last=ic2_=ic2_last=0.0;
  c3_=c3_last=ic3_=ic3_last=0.0;
  c4_=c4_last=ic4_=ic4_last=0.0;
  if(SIMULATE_TEST)
    printf("%6.6s %4.4s %10.10s %10.10s %10.10s %10.10s %10.10s %10.10s\n",
      "t", "dt/2", "ictot", "C0", "C1", "C2", "C3", "C4");
  for(i=0; i<nr; i++) {
    /* delta time / 2 */
    dt2=0.5*(t[i]-t_last); if(dt2<0.0) return(5);
    if(dt2>0.0) {
      /* Compute temp constants */
      a=k01+k12-(k12*k21*dt2/(1.0+dt2*k21));
      b=k03+k34-(k34*k43*dt2/(1.0+dt2*k43));
      am1=ic1_last+dt2*c1_last;
      am2=ic2_last+dt2*c2_last;
      am3=ic3_last+dt2*c3_last;
      am4=ic4_last+dt2*c4_last;
    
      /* Compute the "input" i.e. ctot integral */
      ictot+=(ctot[i]+ctot_last)*dt2;
      /* Compute C1(t) and its integral */
      c1_= ( k01*(1.0-k03*dt2/(1.0+dt2*b))*ictot
            -(a-k01*k03*dt2/(1.0+dt2*b))*am1
            +(k21/(1.0+dt2*k21))*am2
            -(k01/(1.0+dt2*b))*am3
            -(k01*k43*dt2/((1.0+dt2*b)*(1.0+dt2*k43)))*am4
           ) / ( 1.0+dt2*(a-k01*k03*dt2/(1.0+dt2*b)) );
      ic1_= ic1_last + dt2*(c1_+c1_last);
      /* Compute C2(t) and its integral */
      c2_= (k12*ic1_-k21*am2)/(1.0+dt2*k21);
      ic2_= ic2_last + dt2*(c2_+c2_last);

      /* Compute C3(t) and its integral */
      c3_= ( k03*(1.0-k01*dt2/(1.0+dt2*a))*ictot
            -(b-k01*k03*dt2/(1.0+dt2*a))*am3
            +(k43/(1.0+dt2*k43))*am4
            -(k03/(1.0+dt2*a))*am1
            -(k03*k21*dt2/((1.0+dt2*a)*(1.0+dt2*k21)))*am2
           ) / ( 1.0+dt2*(b-k01*k03*dt2/(1.0+dt2*a)) );
      ic3_= ic3_last + dt2*(c3_+c3_last);
      /* Compute C4(t) and its integral */
      c4_= (k34*ic3_-k43*am4)/(1.0+dt2*k43);
      ic4_= ic4_last + dt2*(c4_+c4_last);

      /* Compute the C0(t) */
      c0_=ctot[i]-c1_-c3_; /*if(c0_<0.0) return(-1);*/
    }
    if(SIMULATE_TEST)
      printf("%6.2f %4.2f %10.2e %10.2e %10.2e %10.2e %10.2e %10.2e\n",
        t[i], dt2, ictot, c0_, c1_, c2_, c3_, c4_ );
    /* Set output data */
    if(c0) c0[i]=c0_; 
    if(c1) c1[i]=c1_; 
    if(c3) c3[i]=c3_;
    /* Prepare for the next sample */
    c0_last=c0_; c1_last=c1_; c2_last=c2_; c3_last=c3_; c4_last=c4_;
    ic0_last=ic0_; ic1_last=ic1_; ic2_last=ic2_; ic3_last=ic3_; ic4_last=ic4_;
    ctot_last=ctot[i]; t_last=t[i];
  } /* next sample */

  return(0);
}
/*****************************************************************************/

/*****************************************************************************/
/** Simulate parent tracer TAC using plasma metabolite model TPCMOD0009C,
 *  http://www.turkupetcentre.net/reports/tpcmod0009_app_c.pdf
\return Returns 0 if successful.
*/
int simTPCMOD0009c(
  /** Sample times */
  double *t,
  /** Total plasma TAC */
  double *ctot,
  /** Sample number */
  int nr,
  /** km */
  double km,
  /** k1m */
  double k1m,
  /** k2m */
  double k2m,
  /** k3m */
  double k3m,
  /** k4m */
  double k4m,
  /** Pointer to array where parent tracer TAC will be written;
   *  NULL, if not needed. */
  double *ca,
  /** Pointer to array where metabolized tracer TAC will be written;
   *  NULL, if not needed. */
  double *cm
) {
  int i;
  double t_last, dt2;
  double ctoti, ctot_last;
  double a, b;
  double ct1m, ct1mi, ct1m_last, ct1mi_last;
  double ct2m, ct2mi, ct2m_last, ct2mi_last;
  double cpm, cpmi, cpm_last, cpmi_last;
  

  /* Check for data */
  if(nr<2) return 1;
  if(t==NULL || ctot==NULL || ca==NULL) return 2;

  /* Calculate curves */
  t_last=0.0; if(t[0]<t_last) t_last=t[0]; 
  ctoti=ctot_last=0.0;
  ct1m_last=ct1mi_last=ct2m_last=ct2mi_last=cpm_last=cpmi_last=0.0;
  //ct1=ct2=ct3=ct1i=ct2i=ct3i=0.0;
  for(i=0; i<nr; i++) {
    /* delta time / 2 */
    dt2=0.5*(t[i]-t_last);
    /* calculate values */
    if(dt2<0.0) {
      return 5;
    } else {
      /* arterial integral */
      ctoti+=(ctot[i]+ctot_last)*dt2;
      /* partial results */
      a=k4m/(1.0+k4m*dt2);
      b=(k1m-km)/(1.0+dt2*(k1m+k3m-k3m*a*dt2));
      /* 1st tissue compartment and its integral */
      ct1m = (km*ctoti - k2m*(1.0-dt2*b)*(ct1mi_last+dt2*ct1m_last)
             + b*(cpmi_last+dt2*cpm_last) + a*b*dt2*(ct2mi_last+dt2*ct2m_last)
	     ) / (1.0 + dt2*k2m*(1.0-dt2*b));
      ct1mi = ct1mi_last + dt2*(ct1m_last+ct1m);
      /* Metabolite plasma compartment and its integral */
      cpm = (k2m*ct1mi + a*(ct2mi_last+dt2*ct2m_last)
             - (k1m+k3m-k3m*dt2*a)*(cpmi_last+dt2*cpm_last)
            ) / (1.0 + dt2*(k1m+k3m-k3m*dt2*a));
      cpmi = cpmi_last + dt2*(cpm_last+cpm);
      /* 2nd tissue compartment and its integral */
      ct2m = (k3m*cpmi - k4m*(ct2mi_last+dt2*ct2m_last)) / (1.0 + dt2*k4m);
      ct2mi = ct2mi_last + dt2*(ct2m_last+ct2m);
    }
    /* copy values to argument arrays; set very small values to zero */
    if(ca!=NULL) {ca[i]=ctot[i]-cpm; if(fabs(ca[i])<1.0e-12) ca[i]=0.0;}
    if(cm!=NULL) {cm[i]=cpm; if(fabs(cm[i])<1.0e-12) cm[i]=0.0;}
    /* prepare to the next loop */
    t_last=t[i]; ctot_last=ctot[i];
    ct1m_last=ct1m; ct1mi_last=ct1mi;
    ct2m_last=ct2m; ct2mi_last=ct2mi;
    cpm_last=cpm; cpmi_last=cpmi;
  } // next sample
  return 0;
}
/*****************************************************************************/

/*****************************************************************************/
/** Simulate myocardial tissue TAC using Iida's compartment model.
    Memory for ct must be allocated in the calling program.
    The units of rate constants must be related to the time unit;
    1/min and min, or 1/sec and sec.
    Function returns 0 when successful, else a value >= 1.
 */
int simMBF_v1(
  /** Array of time values */
  double *t,
  /** Input activities */
  double *ci,
  /** Number of values in TACs */
  int nr,
  /** Apparent k1 */
  double k1,
  /** Apparent k2 */
  double k2,
  /** Vfit */
  double Vfit,
  /** Pointer for TAC array to be simulated; must be allocated */
  double *ct
) {
  int i;
  double dt2;
  double cii, ci_last, t_last;
  double ct_last, cti, cti_last;


  /* Check for data */
  if(nr<2) return(1);
  if(ct==NULL) return(2);

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
/** Simulates tissue TAC using 1 tissue compartment model
    plasma TAC, at plasma TAC times.

    Memory for ct must be allocated in the calling program.

    The units of rate constants must be related to the time unit; 1/min and min,
    or 1/sec and sec.

\return Function returns 0 when successful, else a value >= 1.
*/
int simC1_v1(
  /** Array of time values */
  double *t,
  /** Array of arterial activities */
  double *ca,
  /** Number of values in TACs */
  int nr,
  /** Rate constant of the model */
  double k1,
  /** Rate constant of the model */
  double k2,
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
  if(ct==NULL) return 2;

  /* Check actual parameter number */
  if(k1<0.0) return 3;

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
/** Simulates tissue TAC using dual-input tissue compartment model
    (1-3 compartments in series for tracer1, and 1 compartment for tracer2)
    at plasma TAC times, considering also contribution of arterial and venous
    vasculature, but no exchange between compartments for tracer1 and tracer2.
    The units of rate constants must be related to the time unit; 1/min and min,
    or 1/sec and sec.
    If blood flow is set to 0, function assumes that f>>k1, and Cvb=Cab.

\return Function returns 0 when successful, else a value >= 1.
*/
int simC3DIvs(
  /** Array of time values */
  double *t,
  /** Array of arterial plasma activities of tracer1 */
  double *ca1,
  /** Array of arterial plasma activities of tracer2 */
  double *ca2,
  /** Array of arterial blood activities */
  double *cb,
  /** Number of values in TACs */
  int nr,
  /** Rate constant of the model for tracer1 (from plasma to C1) */
  double k1,
  /** Rate constant of the model for tracer1 (from C1 to plasma) */
  double k2,
  /** Rate constant of the model for tracer1 (from C1 to C2) */
  double k3,
  /** Rate constant of the model for tracer1 (from C2 to C1) */
  double k4,
  /** Rate constant of the model for tracer1 (from C2 to C3) */
  double k5,
  /** Rate constant of the model for tracer1 (from C3 to C2) */
  double k6,
  /** Rate constant of the model for tracer2 (from plasma to C4) */
  double k1b,
  /** Rate constant of the model for tracer2 (from C4 to plasma) */
  double k2b,
  /** Blood flow; if 0, function assumes that f>>k1, and Cvb=Cab. */
  double f,
  /** Vascular volume fraction */
  double vb,
  /** Arterial fraction of vascular volume */
  double fa,
  /** Pointer for TAC array to be simulated; must be allocated */
  double *scpet,
  /** Pointer for 1st tracer1 compartment TAC, or NULL if not needed */
  double *sct1,
  /** Pointer for 2nd tracer1 compartment TAC, or NULL if not needed */
  double *sct2,
  /** Pointer for 3rd tracer1 compartment TAC, or NULL if not needed */
  double *sct3,
  /** Pointer for 1st tracer2 compartment TAC, or NULL if not needed */
  double *sct1b,
  /** Pointer for arterial TAC in tissue, or NULL if not needed */
  double *sctab,
  /** Pointer for venous TAC in tissue, or NULL if not needed */
  double *sctvb
) {
  int i;
  double b, c, d, e, w, z, dt2, va, vv;
  double ca1i, ca1_last, ca2i, ca2_last, t_last, dct, cvb;
  double ct1, ct1_last, ct2, ct2_last, ct3, ct3_last;
  double ct1i, ct1i_last, ct2i, ct2i_last, ct3i, ct3i_last;
  double ct1b, ct1b_last, ct1bi, ct1bi_last;


  /* Check for data */
  if(nr<2) return 1;
  if(scpet==NULL) return 2;

  /* Check parameters */
  if(k1<0.0) return 3;
  if(vb<0.0 || vb>=1.0) return 4;
  if(fa<=0.0 || fa>1.0) return 5;
  va=fa*vb; vv=(1.0-fa)*vb;

  /* Calculate curves */
  t_last=0.0; if(t[0]<t_last) t_last=t[0];
  ca1i=ca1_last=ca2i=ca2_last=0.0;
  ct1_last=ct2_last=ct3_last=ct1i_last=ct2i_last=ct3i_last=0.0;
  ct1b_last=ct1bi_last=0.0;
  ct1=ct2=ct3=ct1b=ct1i=ct2i=ct3i=ct1bi=0.0;
  for(i=0; i<nr; i++) {
    /* delta time / 2 */
    dt2=0.5*(t[i]-t_last);
    /* calculate values */
    if(dt2<0.0) {
      return 5;
    } else if(dt2>0.0) {
      /* arterial integrals */
      ca1i+=(ca1[i]+ca1_last)*dt2;
      ca2i+=(ca2[i]+ca2_last)*dt2;
      /* partial results */
      b=ct1i_last+dt2*ct1_last;
      c=ct2i_last+dt2*ct2_last;
      d=ct3i_last+dt2*ct3_last;
      e=ct1bi_last+dt2*ct1b_last;
      w=k4 + k5 - (k5*k6*dt2)/(1.0+k6*dt2);
      z=1.0+w*dt2;
      /* 1st tissue compartment and its integral */
      ct1 = (
          + k1*z*ca1i + (k3*k4*dt2 - (k2+k3)*z)*b
          + k4*c + k4*k6*dt2*d/(1.0+k6*dt2)
        ) / ( z*(1.0 + dt2*(k2+k3)) - k3*k4*dt2*dt2 );
      ct1i = ct1i_last + dt2*(ct1_last+ct1);
      ct1b = (k1b*ca2i - k2b*e) / (1.0 + dt2*k2b);
      ct1bi = ct1bi_last + dt2*(ct1b_last+ct1b);
      /* 2nd tissue compartment and its integral */
      ct2 = (k3*ct1i - w*c + k6*d/(1.0+k6*dt2)) / z;
      ct2i = ct2i_last + dt2*(ct2_last+ct2);
      /* 3rd tissue compartment and its integral */
      ct3 = (k5*ct2i - k6*d) / (1.0 + k6*dt2);
      ct3i = ct3i_last + dt2*(ct3_last+ct3);
    }
    /* Venous curve */
    if(f>0.) {
      dct = k1*ca1[i] - k2*ct1 + k1b*ca2[i] - k2b*ct1b;
      cvb = cb[i] - dct/f;
    } else cvb=cb[i];
    /* copy values to argument arrays; set very small values to zero */
    scpet[i]= va*cb[i] + vv*cvb + (1.0-vb)*(ct1+ct2+ct3+ct1b);
    if(fabs(scpet[i])<1.0e-12) scpet[i]=0.0;
    if(sct1!=NULL) {sct1[i]=(1.0-vb)*ct1; if(fabs(sct1[i])<1.0e-12) sct1[i]=0.0;}
    if(sct2!=NULL) {sct2[i]=(1.0-vb)*ct2; if(fabs(sct2[i])<1.0e-12) sct2[i]=0.0;}
    if(sct3!=NULL) {sct3[i]=(1.0-vb)*ct3; if(fabs(sct3[i])<1.0e-12) sct3[i]=0.0;}
    if(sct1b!=NULL) {
      sct1b[i]=(1.0-vb)*ct1b; if(fabs(sct1b[i])<1.0e-12) sct1b[i]=0.0;}
    if(sctab!=NULL) {sctab[i]=va*cb[i]; if(fabs(sctab[i])<1.0e-12) sctab[i]=0.0;}
    if(sctvb!=NULL) {sctvb[i]=vv*cvb; if(fabs(sctvb[i])<1.0e-12) sctvb[i]=0.0;}
    /* prepare to the next loop */
    t_last=t[i]; ca1_last=ca1[i]; ca2_last=ca2[i];
    ct1_last=ct1; ct1i_last=ct1i;
    ct2_last=ct2; ct2i_last=ct2i;
    ct3_last=ct3; ct3i_last=ct3i;
    ct1b_last=ct1b; ct1bi_last=ct1bi;
  }

  return 0;
}
/*****************************************************************************/

/*****************************************************************************/
/** Simulates tissue TAC using dual-input tissue compartment model
    (compartments 2 and 3 in parallel for tracer1, and 1 compartment for
    tracer2) at plasma TAC sample times, considering also contribution of
    arterial and venous vasculature, and transfer of tracer1 to tracer2.
    The units of rate constants must be related to the time unit; 1/min and min,
    or 1/sec and sec.
    If blood flow is set to 0, function assumes that f>>k1, and Cvb=Cab.
    Reference: TPCMOD0001 Appendix C.
    Tested with program p2t_di -parallel.

\return Function returns 0 when successful, else a value >= 1.
*/
int simC4DIvp(
  /** Array of time values */
  double *t,
  /** Array of arterial plasma activities of tracer1 (parent tracer) */
  double *ca1,
  /** Array of arterial plasma activities of tracer2 (metabolite) */
  double *ca2,
  /** Array of (total) arterial blood activities */
  double *cb,
  /** Number of values in TACs */
  int nr,
  /** Rate constant of the model for tracer1 (from plasma to C1) */
  double k1,
  /** Rate constant of the model for tracer1 (from C1 to plasma) */
  double k2,
  /** Rate constant of the model for tracer1 (from C1 to C2) */
  double k3,
  /** Rate constant of the model for tracer1 (from C2 to C1) */
  double k4,
  /** Rate constant of the model for tracer1 (from C1 to C3) */
  double k5,
  /** Rate constant of the model for tracer1 (from C3 to C1) */
  double k6,
  /** Rate constant of the model for tracer1 (from C3 to plasma) */
  double k7,
  /** Rate constant of the model (from tracer1 in C1 to tracer2 in C4) */
  double km,
  /** Rate constant of the model for tracer2 (from plasma to C4) */
  double k1b,
  /** Rate constant of the model for tracer2 (from C4 to plasma) */
  double k2b,
  /** Blood flow; if 0, function assumes that f>>k1, and Cvb=Cab */
  double f,
  /** Vascular volume fraction (0<=Vb<1) */
  double vb,
  /** Arterial fraction of vascular volume (0<=fa<=1) */
  double fa,
  /** Pointer for TAC array to be simulated; must be allocated */
  double *scpet,
  /** Pointer for 1st tracer1 compartment TAC, or NULL if not needed */
  double *sct1,
  /** Pointer for 2nd tracer1 compartment TAC, or NULL if not needed */
  double *sct2,
  /** Pointer for 3rd tracer1 compartment TAC, or NULL if not needed */
  double *sct3,
  /** Pointer for 1st tracer2 compartment TAC, or NULL if not needed */
  double *sct1b,
  /** Pointer for arterial TAC in tissue, or NULL if not needed */
  double *sctab,
  /** Pointer for venous TAC in tissue, or NULL if not needed */
  double *sctvb,
  /** Verbose level; if zero, then nothing is printed into stdout or stderr */
  int verbose
) {
  int i;
  double b, c, d, e, pt, qt, dt2, va, vv;
  double ca1i, ca1_last, ca2i, ca2_last, t_last, dct, cvb;
  double ct1, ct1_last, ct2, ct2_last, ct3, ct3_last;
  double ct1i, ct1i_last, ct2i, ct2i_last, ct3i, ct3i_last;
  double ct1b, ct1b_last, ct1bi, ct1bi_last;


  if(verbose>0) {
    printf("simC4DIvp()\n");
    if(verbose>1) {
      printf("  k1 := %g\n", k1);
      printf("  k2 := %g\n", k2);
      printf("  k3 := %g\n", k3);
      printf("  k4 := %g\n", k4);
      printf("  k5 := %g\n", k5);
      printf("  k6 := %g\n", k6);
      printf("  k7 := %g\n", k7);
      printf("  km := %g\n", km);
      printf("  k1b := %g\n", k1b);
      printf("  k2b := %g\n", k2b);
      printf("  vb := %g\n", vb);
      printf("  fa := %g\n", fa);
      printf("  f := %g\n", f);
    }
  }

  /* Check for data */
  if(nr<2) return 1;
  if(scpet==NULL) return 2;

  /* Check parameters */
  if(k1<0.0 || k1b<0.0) return 3;
  if(vb<0.0 || vb>=1.0) return 4;
  if(fa<0.0 || fa>1.0) return 5;
  va=fa*vb; vv=(1.0-fa)*vb;

  /* Calculate curves */
  t_last=0.0; if(t[0]<t_last) t_last=t[0];
  ca1i=ca1_last=ca2i=ca2_last=0.0;
  ct1_last=ct2_last=ct3_last=ct1i_last=ct2i_last=ct3i_last=ct1b_last=
  ct1bi_last=0.0;
  ct1=ct2=ct3=ct1b=ct1i=ct2i=ct3i=ct1bi=0.0;
  for(i=0; i<nr; i++) {
    /* delta time / 2 */
    dt2=0.5*(t[i]-t_last);
    /* calculate values */
    if(dt2<0.0) {
      return 5;
    } else if(dt2>0.0) {
      /* arterial integrals */
      ca1i+=(ca1[i]+ca1_last)*dt2;
      ca2i+=(ca2[i]+ca2_last)*dt2;
      /* partial results */
      b=ct1i_last+dt2*ct1_last;
      c=ct2i_last+dt2*ct2_last;
      d=ct3i_last+dt2*ct3_last;
      e=ct1bi_last+dt2*ct1b_last;
      pt=k6+k7;
      qt=k2+k3+k5+km-(k3*k4*dt2)/(1.0+k4*dt2)-(k5*k6*dt2)/(1.0+pt*dt2);
      /* 1st tissue compartment and its integral */
      ct1 = (k1/(1.0+qt*dt2))*ca1i
	  - (qt/(1.0+qt*dt2))*b
          + (k4/((1.0+qt*dt2)*(1.0+k4*dt2)))*c
          + (k6/((1.0+qt*dt2)*(1.0+pt*dt2)))*d;
      ct1i = ct1i_last + dt2*(ct1_last+ct1);
      /* 2nd tissue compartment and its integral */
      ct2 = (k3/(1.0+k4*dt2))*ct1i
          - (k4/(1.0+k4*dt2))*c;
      ct2i = ct2i_last + dt2*(ct2_last+ct2);
      /* 3rd tissue compartment and its integral */
      ct3 = (k5/(1.0+pt*dt2))*ct1i
          - (pt/(1.0+pt*dt2))*d;
      ct3i = ct3i_last + dt2*(ct3_last+ct3);
      /* 4th tissue compartment (the 1st for tracer 2) and its integral */
      ct1b = (k1b/(1.0+k2b*dt2))*ca2i
           - (k2b/(1.0+k2b*dt2))*e
           + (km/(1.0+k2b*dt2))*ct1i;
      ct1bi = ct1bi_last + dt2*(ct1b_last+ct1b);
    }
    /* Venous curve */
    if(f>0.) {
      dct = k1*ca1[i] - k2*ct1 - k7*ct3 + k1b*ca2[i] - k2b*ct1b;
      cvb = cb[i] - dct/f;
    } else cvb=cb[i];
    /* copy values to argument arrays; set very small values to zero */
    scpet[i]= va*cb[i] + vv*cvb + (1.0-vb)*(ct1+ct2+ct3+ct1b);
    if(fabs(scpet[i])<1.0e-12) scpet[i]=0.0;
    if(sct1!=NULL) {
      sct1[i]=(1.0-vb)*ct1; if(fabs(sct1[i])<1.0e-12) sct1[i]=0.0;}
    if(sct2!=NULL) {
      sct2[i]=(1.0-vb)*ct2; if(fabs(sct2[i])<1.0e-12) sct2[i]=0.0;}
    if(sct3!=NULL) {
      sct3[i]=(1.0-vb)*ct3; if(fabs(sct3[i])<1.0e-12) sct3[i]=0.0;}
    if(sct1b!=NULL) {
      sct1b[i]=(1.0-vb)*ct1b; if(fabs(sct1b[i])<1.0e-12) sct1b[i]=0.0;}
    if(sctab!=NULL) {
      sctab[i]=va*cb[i]; if(fabs(sctab[i])<1.0e-12) sctab[i]=0.0;}
    if(sctvb!=NULL) {
      sctvb[i]=vv*cvb; if(fabs(sctvb[i])<1.0e-12) sctvb[i]=0.0;}
    /* prepare to the next loop */
    t_last=t[i]; ca1_last=ca1[i]; ca2_last=ca2[i];
    ct1_last=ct1; ct1i_last=ct1i;
    ct2_last=ct2; ct2i_last=ct2i;
    ct3_last=ct3; ct3i_last=ct3i;
    ct1b_last=ct1b; ct1bi_last=ct1bi;
  }
  
  if(verbose>2) {
    printf("AUC 0-%g:\n", t_last);
    printf(" ca1i := %g\n", ca1i);
    printf(" ca2i := %g\n", ca2i);
    printf(" ct1i := %g\n", ct1i_last);
    printf(" ct2i := %g\n", ct2i_last);
    printf(" ct3i := %g\n", ct3i_last);
    printf(" ct1bi := %g\n", ct1bi_last);
  }

  return 0;
}
/*****************************************************************************/

/*****************************************************************************/
/** Simulates tissue TAC using dual-input tissue compartment model
    (1-3 compartments in series for tracer1, and 1 compartment for tracer2)
    at plasma TAC times, considering also contribution of arterial and venous
    vasculature, and transfer of tracer1 to tracer2.
    The units of rate constants must be related to the time unit; 1/min and min,
    or 1/sec and sec.
    If blood flow is set to 0, function assumes that f>>k1, and Cvb=Cab.
    Reference: TPCMOD0001 Appendix B.
    Tested with program p2t_di -series.

\return Function returns 0 when successful, else a value >= 1.
*/
int simC4DIvs(
  /** Array of time values */
  double *t,
  /** Array of arterial plasma activities of tracer1 (parent tracer) */
  double *ca1,
  /** Array of arterial plasma activities of tracer2 (metabolite) */
  double *ca2,
  /** Array of (total) arterial blood activities */
  double *cb,
  /** Number of values in TACs */
  int nr,
  /** Rate constant of the model for tracer1 (from plasma to C1) */
  double k1,
  /** Rate constant of the model for tracer1 (from C1 to plasma) */
  double k2,
  /** Rate constant of the model for tracer1 (from C1 to C2) */
  double k3,
  /** Rate constant of the model for tracer1 (from C2 to C1) */
  double k4,
  /** Rate constant of the model for tracer1 (from C2 to C3) */
  double k5,
  /** Rate constant of the model for tracer1 (from C3 to C2) */
  double k6,
  /** Rate constant of the model for tracer1 (from C3 to plasma) */
  double k7,
  /** Rate constant of the model (from tracer1 in C1 to tracer2 in C4) */
  double km,
  /** Rate constant of the model for tracer2 (from plasma to C4) */
  double k1b,
  /** Rate constant of the model for tracer2 (from C4 to plasma) */
  double k2b,
  /** Blood flow; if 0, function assumes that f>>k1, and Cvb=Cab. */
  double f,
  /** Vascular volume fraction */
  double vb,
  /** Arterial fraction of vascular volume */
  double fa,
  /** Pointer for TAC array to be simulated; must be allocated */
  double *scpet,
  /** Pointer for 1st tracer1 compartment TAC, or NULL if not needed */
  double *sct1,
  /** Pointer for 2nd tracer1 compartment TAC, or NULL if not needed */
  double *sct2,
  /** Pointer for 3rd tracer1 compartment TAC, or NULL if not needed */
  double *sct3,
  /** Pointer for 1st tracer2 compartment TAC, or NULL if not needed */
  double *sct1b,
  /** Pointer for arterial TAC in tissue, or NULL if not needed */
  double *sctab,
  /** Pointer for venous TAC in tissue, or NULL if not needed */
  double *sctvb,
  /** Verbose level; if zero, then nothing is printed into stdout or stderr */
  int verbose
) {
  int i;
  double b, c, d, e, pt, qt, rt, dt2, va, vv;
  double ca1i, ca1_last, ca2i, ca2_last, t_last, dct, cvb;
  double ct1, ct1_last, ct2, ct2_last, ct3, ct3_last;
  double ct1i, ct1i_last, ct2i, ct2i_last, ct3i, ct3i_last;
  double ct1b, ct1b_last, ct1bi, ct1bi_last;


  if(verbose>0) {
    printf("simC4DIvs()\n");
    if(verbose>1) {
      printf("  k1 := %g\n", k1);
      printf("  k2 := %g\n", k2);
      printf("  k3 := %g\n", k3);
      printf("  k4 := %g\n", k4);
      printf("  k5 := %g\n", k5);
      printf("  k6 := %g\n", k6);
      printf("  k7 := %g\n", k7);
      printf("  km := %g\n", km);
      printf("  k1b := %g\n", k1b);
      printf("  k2b := %g\n", k2b);
      printf("  vb := %g\n", vb);
      printf("  fa := %g\n", fa);
      printf("  f := %g\n", f);
    }
  }

  /* Check for data */
  if(nr<2) return 1;
  if(scpet==NULL) return 2;

  /* Check parameters */
  if(k1<0.0 || k1b<0.0) return 3;
  if(vb<0.0 || vb>=1.0) return 4;
  if(fa<0.0 || fa>1.0) return 5;
  va=fa*vb; vv=(1.0-fa)*vb;

  /* Calculate curves */
  t_last=0.0; if(t[0]<t_last) t_last=t[0];
  ca1i=ca1_last=ca2i=ca2_last=0.0;
  ct1_last=ct2_last=ct3_last=ct1i_last=ct2i_last=ct3i_last=ct1b_last=
  ct1bi_last=0.0;
  ct1=ct2=ct3=ct1b=ct1i=ct2i=ct3i=ct1bi=0.0;
  for(i=0; i<nr; i++) {
    /* delta time / 2 */
    dt2=0.5*(t[i]-t_last);
    /* calculate values */
    if(dt2<0.0) {
      return 5;
    } else if(dt2>0.0) {
      /* arterial integrals */
      ca1i+=(ca1[i]+ca1_last)*dt2;
      ca2i+=(ca2[i]+ca2_last)*dt2;
      //printf("ca1[%d]=%g int=%g\n", i, ca1[i], ca1i);
      //printf("ca2[%d]=%g int=%g\n", i, ca2[i], ca2i);
      /* partial results */
      b=ct1i_last+dt2*ct1_last;
      c=ct2i_last+dt2*ct2_last;
      d=ct3i_last+dt2*ct3_last;
      e=ct1bi_last+dt2*ct1b_last;
      pt=k6+k7;
      qt=k4+k5-(k5*k6*dt2)/(1.0+pt*dt2);
      rt=k2+k3+km-(k3*k4*dt2)/(1.0+qt*dt2);
      /* 1st tissue compartment and its integral */
      ct1 = (k1/(1.0+rt*dt2))*ca1i
	  - (rt/(1.0+rt*dt2))*b
          + (k4/((1.0+qt*dt2)*(1.0+rt*dt2)))*c
          + ((k4*k6*dt2)/((1.0+pt*dt2)*(1.0+qt*dt2)*(1.0+rt*dt2)))*d;
      ct1i = ct1i_last + dt2*(ct1_last+ct1);
      /* 2nd tissue compartment and its integral */
      ct2 = (k3/(1.0+qt*dt2))*ct1i
          - (qt/(1.0+qt*dt2))*c
	  + (k6/((1.0+pt*dt2)*(1.0+qt*dt2)))*d;
      ct2i = ct2i_last + dt2*(ct2_last+ct2);
      /* 3rd tissue compartment and its integral */
      ct3 = (k5/(1.0+pt*dt2))*ct2i
          - (pt/(1.0+pt*dt2))*d;
      ct3i = ct3i_last + dt2*(ct3_last+ct3);
      /* 4th tissue compartment (the 1st for tracer 2) and its integral */
      ct1b = (k1b/(1.0+k2b*dt2))*ca2i
           - (k2b/(1.0+k2b*dt2))*e
           + (km/(1.0+k2b*dt2))*ct1i;
      ct1bi = ct1bi_last + dt2*(ct1b_last+ct1b);
    }
    /* Venous curve */
    if(f>0.) {
      dct = k1*ca1[i] - k2*ct1 - k7*ct3 + k1b*ca2[i] - k2b*ct1b;
      cvb = cb[i] - dct/f;
    } else cvb=cb[i];
    /* copy values to argument arrays; set very small values to zero */
    scpet[i]= va*cb[i] + vv*cvb + (1.0-vb)*(ct1+ct2+ct3+ct1b);
    if(fabs(scpet[i])<1.0e-12) scpet[i]=0.0;
    if(sct1!=NULL) {
      sct1[i]=(1.0-vb)*ct1; if(fabs(sct1[i])<1.0e-12) sct1[i]=0.0;}
    if(sct2!=NULL) {
      sct2[i]=(1.0-vb)*ct2; if(fabs(sct2[i])<1.0e-12) sct2[i]=0.0;}
    if(sct3!=NULL) {
      sct3[i]=(1.0-vb)*ct3; if(fabs(sct3[i])<1.0e-12) sct3[i]=0.0;}
    if(sct1b!=NULL) {
      sct1b[i]=(1.0-vb)*ct1b; if(fabs(sct1b[i])<1.0e-12) sct1b[i]=0.0;}
    if(sctab!=NULL) {
      sctab[i]=va*cb[i]; if(fabs(sctab[i])<1.0e-12) sctab[i]=0.0;}
    if(sctvb!=NULL) {
      sctvb[i]=vv*cvb; if(fabs(sctvb[i])<1.0e-12) sctvb[i]=0.0;}
    /* prepare to the next loop */
    t_last=t[i]; ca1_last=ca1[i]; ca2_last=ca2[i];
    ct1_last=ct1; ct1i_last=ct1i;
    ct2_last=ct2; ct2i_last=ct2i;
    ct3_last=ct3; ct3i_last=ct3i;
    ct1b_last=ct1b; ct1bi_last=ct1bi;
  }
  
  if(verbose>2) {
    printf("AUC 0-%g:\n", t_last);
    printf(" ca1i := %g\n", ca1i);
    printf(" ca2i := %g\n", ca2i);
    printf(" ct1i := %g\n", ct1i_last);
    printf(" ct2i := %g\n", ct2i_last);
    printf(" ct3i := %g\n", ct3i_last);
    printf(" ct1bi := %g\n", ct1bi_last);
  }

  return 0;
}
/*****************************************************************************/

/*****************************************************************************/
/** Simulate the effect of dispersion on a time-activity curve.
\return Returns 0 when successful, otherwise <>0.
 */
int simDispersion(
  /** Array of sample times */
  double *x,
  /** Array of sample values, which will be replaced here by dispersion added
   *  values */
  double *y,
  /** Nr of samples */
  int n,
  /** First dispersion time constant (zero if no dispersion); in same
   *  time unit as sample times */ 
  double tau1,
  /** 2nd dispersion time constant (zero if no dispersion); in same
   *  time unit as sample times */ 
  double tau2,
  /** Array for temporary data, for at least n samples */
  double *tmp
) {
  int i, ret;
  double k;

  /* Check input */
  if(x==NULL || y==NULL || tmp==NULL || n<2) return 1;
  if(tau1<0.0 || tau2<0.0) return 2;

  /* First dispersion */
  if(tau1>0.0) {
    k=1.0/tau1;
    ret=simC1_v1(x, y, n, k, k, tmp); if(ret!=0) return 100+ret;
    for(i=0; i<n; i++) y[i]=tmp[i];
  }

  /* Second dispersion */
  if(tau2>0.0) {
    k=1.0/tau2;
    ret=simC1_v1(x, y, n, k, k, tmp); if(ret!=0) return 200+ret;
    for(i=0; i<n; i++) y[i]=tmp[i];
  }

  return 0;
}
/*****************************************************************************/

/*****************************************************************************/
/** Simulate tissue and venous blood TACs using dual-input compartment model
 *  for [O-15]O2 (one tissue compartment for [O-15]O2, and another tissue
 *  compartment for its metabolite [O-15]H2O).
 *   
 *  @details
 *  The units of rate constants must be related to the time unit of the data;
 *  1/min and min, or 1/sec and sec.
 * 
 *  @return Function returns 0 when successful, else a value >= 1.
 *  @author Vesa Oikonen
 */
int simOxygen(
  /** Array of sample times */
  double *t,
  /** Array of arterial blood activities of tracer1 ([O-15]O2) */
  double *ca1,
  /** Array of arterial blood activities of tracer2 ([O-15]H2O) */
  double *ca2,
  /** Array of AUC 0-t of arterial tracer1 activities; NULL if not available */
  double *ca1i,
  /** Array of AUC 0-t of arterial tracer2 activities; NULL if not available */
  double *ca2i,
  /** Nr of samples (array lengths) */
  const int n,
  /** Rate constant of the model for tracer1 (from blood to C1) */
  const double k1a,
  /** Rate constant of the model for tracer1 (from C1 to blood) */
  const double k2a,
  /** Rate constant of the model (from tracer1 in C1 to tracer2 in C2) */
  const double km,
  /** Rate constant of the model for tracer2 (from blood to C2) */
  const double k1b,
  /** Rate constant of the model for tracer2 (from C2 to blood) */
  const double k2b,
  /** Vascular volume fraction [0-1) */
  const double vb,
  /** Arterial fraction of vascular volume [0-1] */
  const double fa,
  /** Pointer for TTAC array to be simulated; allocate in the calling program
   *  or set to NULL if not needed */
  double *scpet,
  /** Simulated TAC of tracer1 in tissue; allocate in the calling program or
   *  set to NULL if not needed */
  double *sct1,
  /** Simulated TAC of tracer2 in tissue; allocate in the calling program or
   *  set to NULL if not needed */
  double *sct2,
  /** Total arterial contribution to PET TTAC; allocate in the calling program
   *  or set to NULL if not needed */
  double *sctab,
  /** Venous tracer1 contribution to PET TAC; allocate in the calling program or
   *  set to NULL if not needed */
  double *sctvb1,
  /** Venous tarcer1 contribution to PET TAC; allocate in the calling program or
   *  set to NULL if not needed */
  double *sctvb2,
  /** Venous BTAC of tracer1; allocate in the calling program or
   *  set to NULL if not needed */
  double *scvb1,
  /** Venous BTAC of tracer2; allocate in the calling program or
   *  set to NULL if not needed */
  double *scvb2,
  /** Verbose level; if zero, then nothing is printed into stdout or stderr */
  const int verbose
) {
  if(verbose>0) {
    printf("simOxygen()\n");
    if(verbose>1) {
      printf("  k1a := %g\n", k1a);
      printf("  k2a := %g\n", k2a);
      printf("  km := %g\n", km);
      printf("  k1b := %g\n", k1b);
      printf("  k2b := %g\n", k2b);
      printf("  vb := %g\n", vb);
      printf("  fa := %g\n", fa);
      printf("  n := %d\n", n);
    }
  }

  /* Check for data */
  if(n<2) return 1;
  if(t==NULL) return 2;
  if(ca1==NULL || ca2==NULL) return 3;

  /* Check parameters */
  if(k1a<0.0 || k1b<0.0 || k2a<0.0 || k2b<0.0) return(4);
  if(vb<0.0 || vb>=1.0) return(5);
  if(fa<0.0 || fa>1.0) return(6);
  double va=fa*vb;       // arterial volume fraction in tissue
  double vv=(1.0-fa)*vb; // venous volume fraction in tissue

  /* Set initial condition */
  double t_last=0.0; // previous sample time
  if(t[0]<t_last) t_last=t[0];
  /* Concentrations, integrals, and previous values are zero */
  double cba1i=0.0, cba1_last=0.0;
  double cba2i=0.0, cba2_last=0.0;
  double ct1=0.0, ct1_last=0.0, ct1i=0.0, ct1i_last=0.0;
  double ct2=0.0, ct2_last=0.0, ct2i=0.0, ct2i_last=0.0;
  double cvb1=0.0, cvb2=0.0;

  /* Calculate curves */
  double p, q;
  double dt2; // delta t / 2
  for(int i=0; i<n; i++) {
    dt2=0.5*(t[i]-t_last);
    if(dt2>0.0) {
      /* arterial integrals */
      if(ca1i!=NULL) cba1i=ca1i[i]; else cba1i+=(ca1[i]+cba1_last)*dt2;
      if(ca2i!=NULL) cba2i=ca2i[i]; else cba2i+=(ca2[i]+cba2_last)*dt2;
      /* partial results */
      p=ct1i_last+dt2*ct1_last;
      q=ct2i_last+dt2*ct2_last;
      /* 1st tissue compartment and its integral */
      ct1 = (k1a*cba1i - (k2a+km)*p) / (1.0 + dt2*(k2a+km));
      ct1i = ct1i_last + dt2*(ct1_last+ct1);
      /* 2nd tissue compartment and its integral */
      ct2 = (km*ct1i + k1b*cba2i - k2b*q) / (1.0 + dt2*k2b);
      ct2i = ct2i_last + dt2*(ct2_last+ct2);
    }
    /* Venous BTACs */
    if(k1a>0.0 && k2a>0.0) cvb1=ct1/(k1a/k2a);
    else if(k2a>0.0) cvb1=0.0; else cvb1=ca1[i];
    if(k1b>0.0 && k2b>0.0) cvb2=ct2/(k1b/k2b);
    else if(k2b>0.0) cvb2=0.0; else cvb2=ca2[i];
    /* copy values to argument arrays; set very small values to zero */
    if(scpet!=NULL) {
      scpet[i]= va*(ca1[i]+ca2[i]) + vv*(cvb1+cvb2) + (1.0-vb)*(ct1+ct2);
      if(fabs(scpet[i])<1.0e-12) scpet[i]=0.0;
    }
    if(sct1!=NULL) {
      sct1[i]=(1.0-vb)*ct1; 
      if(fabs(sct1[i])<1.0e-12) sct1[i]=0.0;
    }
    if(sct2!=NULL) {
      sct2[i]=(1.0-vb)*ct2; 
      if(fabs(sct2[i])<1.0e-12) sct2[i]=0.0;
    }
    if(sctab!=NULL) {
      sctab[i]=va*(ca1[i]+ca2[i]); 
      if(fabs(sctab[i])<1.0e-12) sctab[i]=0.0;
    }
    if(sctvb1!=NULL) {
      sctvb1[i]=vv*cvb1; 
      if(fabs(sctvb1[i])<1.0e-12) sctvb1[i]=0.0;
    }
    if(sctvb2!=NULL) {
      sctvb2[i]=vv*cvb2; 
      if(fabs(sctvb2[i])<1.0e-12) sctvb2[i]=0.0;
    }
    if(scvb1!=NULL) scvb1[i]=cvb1;
    if(scvb2!=NULL) scvb2[i]=cvb2;
    /* prepare for the next loop */
    t_last=t[i]; 
    cba1_last=ca1[i]; cba2_last=ca2[i];
    ct1_last=ct1; ct1i_last=ct1i;
    ct2_last=ct2; ct2i_last=ct2i;
  } // next sample

  if(verbose>2) {
    printf("AUC 0-%g:\n", t_last);
    printf(" cba1i := %g\n", cba1i);
    printf(" cba2i := %g\n", cba2i);
    printf(" ct1i := %g\n", ct1i_last);
    printf(" ct2i := %g\n", ct2i_last);
  }

  return 0;
}
/*****************************************************************************/

/*****************************************************************************/
