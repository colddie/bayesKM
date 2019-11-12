/// @file integr.c
/// @author Vesa Oikonen, Kaisa Liukko
/// @brief linear interpolation and integration of PET and blood/plasma TACs.
///
#include "libtpcmodel.h"
/*****************************************************************************/

/*****************************************************************************/
/** @brief Linear interpolation and integration.


  It is assumed that both original and new interpolated data represent
  the actual values at specified time points (not from framed data).
  The integration is calculated dot-by-dot.

  If NULL is specified for newy[], newyi[] and/or newyii[], their values are
  not calculated.

  Original data and new x values must be sorted by ascending x.
  Subsequent x values can have equal values, enabling the use of step functions.
  Negative x (time) values can be processed.
  If necessary, the data is extrapolated assuming that:
    1) y[inf]=y[nr-1] and
    2) if x[0]>0, y[0]>0 and x[0]<=x[1]-x[0], then an imaginary line is drawn
       from (0,0) to (x[0],y[0]).
  @return Non-zero in case of an error.
*/

int INTEGR_TEST;

int interpolate(
  /** Input data x (time) values */
  double *x,
  /** Input data y values */
  double *y,
  /** Number of values in input data */
  int nr,
  /** Output data x values */
  double *newx,
  /** Interpolated (extrapolated) y values; NULL can be given if not needed */
  double *newy,
  /** Integrals; NULL can be given if not needed */
  double *newyi,
  /** 2nd integrals; NULL can be given if not needed */
  double *newyii,
  /** Nr of values in output data */
  int newnr
) {
  int i, j;
  double ty, tyi, tyii;
  double ox1, ox2, oy1, oy2, oi1, oi2, oii1, oii2, dt, ndt;


  if(INTEGR_TEST) printf("in interpolate()\n");
  /* Check for data */
  if(nr<1 || newnr<1) return 1;
  /* Check that programmer understood that outp data must've been allocated */
  if(newy==NULL && newyi==NULL && newyii==NULL) return 2;

  /* Initiate first two input samples */
  ox1=ox2=x[0];
  oy1=oy2=y[0];
  /* Extrapolate the initial phase with triangle... */
  if(ox1>0.0) {
    oy1=0.0;
    /* unless:
       -first input sample y<=0
       -initial gap is longer than input TAC sampling frequency
    */
    if(y[0]>0.0 && (nr>1 && x[0]<=x[1]-x[0])) ox1=0.0;
  }
  /* Integrals too */
  oi1=oii1=0.0;
  oi2=oi1 + (ox2-ox1)*(oy1+oy2)/2.0;
  oii2=oii1 + (ox2-ox1)*(oi1+oi2)/2.0;
  if(INTEGR_TEST>2) {
    printf("ox1=%g oy1=%g oi1=%g oii1=%g\n", ox1, oy1, oi1, oii1);
    printf("ox2=%g oy2=%g oi2=%g oii2=%g\n", ox2, oy2, oi2, oii2);
  }

  /* Set interpolated data before input data (even imaginary) to zero */
  j=0; while(j<newnr && newx[j]<ox1) {
    ty=tyi=tyii=0.0; if(INTEGR_TEST>4) printf("  ndt=%g\n", ox1-newx[j]);
    if(INTEGR_TEST>4) printf("  j=%d newx=%g ty=%g tyi=%g tyii=%g\n",
                             j, newx[j], ty, tyi, tyii);
    if(newy!=NULL) newy[j]=ty;
    if(newyi!=NULL) newyi[j]=tyi;
    if(newyii!=NULL) newyii[j]=tyii;
    j++;
  } 

  /* Set interpolated data between ox1 and ox2 */
  dt=ox2-ox1; if(dt>0.0) while(j<newnr && newx[j]<=ox2) {
    ndt=newx[j]-ox1; if(INTEGR_TEST>4) printf("  ndt=%g\n", ndt);
    ty=((oy2-oy1)/dt)*ndt + oy1; if(newy!=NULL) newy[j]=ty;
    tyi=oi1 + 0.5*(ty+oy1)*ndt; if(newyi!=NULL) newyi[j]=tyi;
    if(newyii!=NULL) newyii[j]=tyii= oii1 + 0.5*(tyi+oi1)*ndt;
    if(INTEGR_TEST>4)
      printf("  j=%d newx=%g ty=%g tyi=%g\n", j, newx[j], ty, tyi);
    j++;
  }

  /* Go through input data, sample-by-sample */
  for(i=0; i<nr && j<newnr; i++) {

    ox1=ox2; oy1=oy2; oi1=oi2; oii1=oii2;
    ox2=x[i]; oy2=y[i];
    oi2=oi1 + (ox2-ox1)*(oy1+oy2)/2.0;
    oii2=oii1 + (ox2-ox1)*(oi1+oi2)/2.0;
    if(INTEGR_TEST>3) {
      printf("ox1=%g oy1=%g oi1=%g oii1=%g\n", ox1, oy1, oi1, oii1);
      printf("ox2=%g oy2=%g oi2=%g oii2=%g\n", ox2, oy2, oi2, oii2);
    }
    
    /* Calculate input sample distance */
    dt=ox2-ox1; if(dt<0.0) return 3; else if(dt==0.0) continue;

    /* Any need for interpolation between ox1 and ox2? */
    while(j<newnr && newx[j]<=ox2) {
      ndt=newx[j]-ox1;
      ty=((oy2-oy1)/dt)*ndt + oy1; if(newy!=NULL) newy[j]=ty;
      tyi=oi1 + 0.5*(ty+oy1)*ndt; if(newyi!=NULL) newyi[j]=tyi;
      if(newyii!=NULL) newyii[j]=tyii= oii1 + 0.5*(tyi+oi1)*ndt;
      if(INTEGR_TEST>5)
        printf("  j=%d newx=%g ty=%g tyi=%g\n", j, newx[j], ty, tyi);
      j++;
    }

  } // next input sample

  /* Set interpolated data after input data, assuming steady input */
  while(j<newnr) {
    ndt=newx[j]-ox2;
    ty=oy2;
    tyi=oi2 + oy2*ndt;
    tyii=oii2 + 0.5*(tyi+oi2)*ndt;
    if(newy!=NULL) newy[j]=ty;
    if(newyi!=NULL) newyi[j]=tyi;
    if(newyii!=NULL) newyii[j]=tyii;
    if(INTEGR_TEST>5)
      printf("  j=%d newx=%g ty=%g tyi=%g\n", j, newx[j], ty, tyi);
    j++;
  } 

  if(INTEGR_TEST) printf("out interpolate()\n");
  return 0;
}

/** @brief float version of interpolate().

  It is assumed that both original and new interpolated data represent
  the actual values at specified time points (not from framed data).
  The integration is calculated dot-by-dot.

  If NULL is specified for newy[], newyi[] and/or newyii[], their values are
  not calculated.

  Original data and new x values must be sorted by ascending x.
  Subsequent x values can have equal values, enabling the use of step functions.
  Negative x (time) values can be processed.
  If necessary, the data is extrapolated assuming that:
    1) y[inf]=y[nr-1] and
    2) if x[0]>0, y[0]>0 and x[0]<=x[1]-x[0], then an imaginary line is drawn
       from (0,0) to (x[0],y[0]).
  @return Non-zero in case of an error.
*/
int finterpolate(
  /** Input data x (time) values */
  float *x,
  /** Input data y values */
  float *y,
  /** Number of values in input data */
  int nr,
  /** Output data x values */
  float *newx,
  /** Interpolated (extrapolated) y values; NULL can be given if not needed */
  float *newy,
  /** Integrals; NULL can be given if not needed */
  float *newyi,
  /** 2nd integrals; NULL can be given if not needed */
  float *newyii,
  /** Nr of values in output data */
  int newnr
) {
  int i, j;
  float ty, tyi, tyii;
  float ox1, ox2, oy1, oy2, oi1, oi2, oii1, oii2, dt, ndt;


  if(INTEGR_TEST) printf("in finterpolate()\n");
  /* Check for data */
  if(nr<1 || newnr<1) return 1;
  /* Check that programmer understood that outp data must've been allocated */
  if(newy==NULL && newyi==NULL && newyii==NULL) return 2;

  /* Initiate first two input samples */
  ox1=ox2=x[0];
  oy1=oy2=y[0];
  /* Extrapolate the initial phase with triangle... */
  if(ox1>0.0) {
    oy1=0.0;
    /* unless:
       -first input sample y<=0
       -initial gap is longer than input TAC sampling frequency
    */
    if(y[0]>0.0 && (nr>1 && x[0]<=x[1]-x[0])) ox1=0.0;
  }
  /* Integrals too */
  oi1=oii1=0.0;
  oi2=oi1 + (ox2-ox1)*(oy1+oy2)/2.0;
  oii2=oii1 + (ox2-ox1)*(oi1+oi2)/2.0;

  /* Set interpolated data before input data (even imaginary) to zero */
  j=0; while(j<newnr && newx[j]<ox1) {
    ty=tyi=tyii=0.0;
    if(newy!=NULL) newy[j]=ty;
    if(newyi!=NULL) newyi[j]=tyi;
    if(newyii!=NULL) newyii[j]=tyii;
    j++;
  } 

  /* Set interpolated data between ox1 and ox2 */
  dt=ox2-ox1; if(dt>0.0) while(j<newnr && newx[j]<=ox2) {
    ndt=newx[j]-ox1;
    ty=((oy2-oy1)/dt)*ndt + oy1; if(newy!=NULL) newy[j]=ty;
    tyi=oi1 + 0.5*(ty+oy1)*ndt; if(newyi!=NULL) newyi[j]=tyi;
    if(newyii!=NULL) newyii[j]=tyii= oii1 + 0.5*(tyi+oi1)*ndt;
    j++;
  }

  /* Go through input data, sample-by-sample */
  for(i=0; i<nr && j<newnr; i++) {

    ox1=ox2; oy1=oy2; oi1=oi2; oii1=oii2;
    ox2=x[i]; oy2=y[i];
    oi2=oi1 + (ox2-ox1)*(oy1+oy2)/2.0;
    oii2=oii1 + (ox2-ox1)*(oi1+oi2)/2.0;
    
    /* Calculate input sample distance */
    dt=ox2-ox1; if(dt<0.0) return 3; else if(dt==0.0) continue;

    /* Any need for interpolation between ox1 and ox2? */
    while(j<newnr && newx[j]<=ox2) {
      ndt=newx[j]-ox1;
      ty=((oy2-oy1)/dt)*ndt + oy1; if(newy!=NULL) newy[j]=ty;
      tyi=oi1 + 0.5*(ty+oy1)*ndt; if(newyi!=NULL) newyi[j]=tyi;
      if(newyii!=NULL) newyii[j]=tyii= oii1 + 0.5*(tyi+oi1)*ndt;
      j++;
    }

  } // next input sample

  /* Set interpolated data after input data, assuming steady input */
  while(j<newnr) {
    ndt=newx[j]-ox2;
    ty=oy2;
    tyi=oi2 + oy2*ndt;
    tyii=oii2 + 0.5*(tyi+oi2)*ndt;
    if(newy!=NULL) newy[j]=ty;
    if(newyi!=NULL) newyi[j]=tyi;
    if(newyii!=NULL) newyii[j]=tyii;
    j++;
  } 

  if(INTEGR_TEST) printf("out finterpolate()\n");
  return 0;
}
/*****************************************************************************/

/*****************************************************************************/
/** Linear integration from time 0 to x[0..nr-1].
    If x[0] is >0 and x[0]<=(x[1]-x[0]),
    then the beginning is interpolated from it to (0,0).
\return Returns 0 if OK, or 1, if error.
*/
int integrate(
  /** Original x values; duplicates are not allowed, all must be >=0.
      Data must be sorted by ascending x */
  double *x,
  /** Original y values */
  double *y,
  /** Nr of values */
  int nr,
  /** Array for integrals */
  double *yi
) {
  int j;

  if(nr==1 || x[0]<=(x[1]-x[0])) yi[0]=0.5*y[0]*x[0];
  else yi[0]=0; /*If the gap in the beginning is longer than time 
                 between first and second time point */
  for(j=1; j<nr; j++) yi[j]=yi[j-1]+0.5*(y[j]+y[j-1])*(x[j]-x[j-1]);

  return 0;
}

/** @brief float version of integrate().
 * 
 *  Linear integration from time 0 to x[0..nr-1].
 *  If x[0] is >0 and x[0]<=(x[1]-x[0]),
 *  then the beginning is interpolated from it to (0,0).
\return Returns 0 if OK, or 1, if error.
*/
int fintegrate(
  /** Original x values; duplicates are not allowed, all must be >=0.
      Data must be sorted by ascending x */
  float *x,
  /** Original y values */
  float *y,
  /** Nr of values */
  int nr,
  /** Array for integrals */
  float *yi
) {
  int j;

  if(nr==1 || x[0]<=(x[1]-x[0])) yi[0]=0.5*y[0]*x[0];
  else yi[0]=0; /*If the gap in the beginning is longer than time 
                 between first and second time point */
  for(j=1; j<nr; j++) yi[j]=yi[j-1]+0.5*(y[j]+y[j-1])*(x[j]-x[j-1]);

  return 0;
}
/*****************************************************************************/

/*****************************************************************************/
/** Calculates integrals of PET data at frame end times.

    Data does not have to be continuous, but it must be increasing in time.
    For faster performance in repetitive calls, allocate memory for integral[]
    even if it is not needed.
    If x1[0] is >0 AND x1[0] is <=(x2[0]-x1[0]),
    then the beginning is interpolated from (0, 0) to (x1[0], y1[0]*x1[0]/x)
    where x is midtime of the first frame.
\return Returns 0 if ok.
*/
int petintegrate(
  /** Array of frame start times */
  double *x1,
  /** Array of frame end times */
  double *x2,
  /** Array of y values (avg during each frame) */
  double *y,
  /** Nr of frames */
  int nr,
  /** Output: integral values at frame end times, or NULL */
  double *newyi,
  /** Output: 2nd integral values at frame end times, or NULL */
  double *newyii
) {
  int i, allocated=0;
  double *ti, x, a;


  if(INTEGR_TEST) printf("in petintegrate()\n");
  /* Check that data is increasing in time and is not (much) overlapping */
  if(nr<1 || x1[0]<0) return 1;
  for(i=0; i<nr; i++) if(x2[i]<x1[i]) return 2;
  for(i=1; i<nr; i++) if(x1[i]<=x1[i-1]) return 3;

  /* Allocate memory for temp data, if necessary */
  if(newyi==NULL) {
    allocated=1;
    ti=(double*)malloc(nr*sizeof(double)); if(ti==NULL) return 4;
  } else
    ti=newyi;

  /* Calculate the integral at the end of first frame */
  ti[0]=(x2[0]-x1[0])*y[0];
  /* If frame does not start at 0 time, add the area in the beginning */
  /* But only if the gap in the beginning is less than the lenght of the 
     first frame */
  if(x1[0]>0) {
    if(x1[0]<=x2[0]-x1[0]) {
       x=(x1[0]+x2[0])/2.0; 
       a=(x1[0]*(y[0]/x)*x1[0])/2.0; 
       ti[0]+=a;
    }
  }

  /* Calculate integrals at the ends of following frames */
  for(i=1; i<nr; i++) {
    /* Add the area of this frame to the previous integral */
    a=(x2[i]-x1[i])*y[i]; ti[i]=ti[i-1]+a;
    /* Check whether frames are contiguous */
    if(x1[i]==x2[i-1]) continue;
    /* When not, calculate the area of an imaginary frame */
    x=(x1[i]+x2[i-1])/2.0;
    a=(x1[i]-x2[i-1])*
      (y[i]-(y[i]-y[i-1])*(x2[i]+x1[i]-2.0*x)/(x2[i]+x1[i]-x2[i-1]-x1[i-1]));
    ti[i]+=a;
  }

  /* Copy integrals to output if required */
  if(allocated) for(i=0; i<nr; i++) newyi[i]=ti[i];

  /* Calculate 2nd integrals if required */
  if(newyii!=NULL) {
    newyii[0]=x2[0]*ti[0]/2.0;
    for(i=1; i<nr; i++)
      newyii[i]=newyii[i-1]+(x2[i]-x2[i-1])*(ti[i-1]+ti[i])/2.0;
  }
  
  /* Free memory */
  if(allocated) free((char*)ti);

  if(INTEGR_TEST) printf("out petintegrate()\n");
  return 0;
}

/** @brief Calculates integrals of PET data at frame end times.
    Float version of petintegrate().

    Data does not have to be continuous, but it must be increasing in time.
    For faster performance in repetitive calls, allocate memory for integral[]
    even if it is not needed.
    If x1[0] is >0 AND x1[0] is <=(x2[0]-x1[0]),
    then the beginning is interpolated from (0, 0) to (x1[0], y1[0]*x1[0]/x)
    where x is midtime of the first frame.
\return Returns 0 if ok.
*/
int fpetintegrate(
  /** Array of frame start times */
  float *x1,
  /** Array of frame end times */
  float *x2,
  /** Array of y values (avg during each frame) */
  float *y,
  /** Nr of frames */
  int nr,
  /** Output: integral values at frame end times, or NULL */
  float *newyi,
  /** Output: 2nd integral values at frame end times, or NULL */
  float *newyii
) {
  int i, allocated=0;
  float *ti, x, a;


  if(INTEGR_TEST) printf("in fpetintegrate()\n");
  /* Check that data is increasing in time and is not (much) overlapping */
  if(nr<1 || x1[0]<0) return 1;
  for(i=0; i<nr; i++) if(x2[i]<x1[i]) return 2;
  for(i=1; i<nr; i++) if(x1[i]<=x1[i-1]) return 3;

  /* Allocate memory for temp data, if necessary */
  if(newyi==NULL) {
    allocated=1;
    ti=(float*)malloc(nr*sizeof(float)); if(ti==NULL) return 4;
  } else
    ti=newyi;

  /* Calculate the integral at the end of first frame */
  ti[0]=(x2[0]-x1[0])*y[0];
  /* If frame does not start at 0 time, add the area in the beginning */
  /* But only if the gap in the beginning is less than the lenght of the 
     first frame */
  if(x1[0]>0) {
    if(x1[0]<=x2[0]-x1[0]) {
       x=(x1[0]+x2[0])/2.0; 
       a=(x1[0]*(y[0]/x)*x1[0])/2.0; 
       ti[0]+=a;
    }
  }

  /* Calculate integrals at the ends of following frames */
  for(i=1; i<nr; i++) {
    /* Add the area of this frame to the previous integral */
    a=(x2[i]-x1[i])*y[i]; ti[i]=ti[i-1]+a;
    /* Check whether frames are contiguous */
    if(x1[i]==x2[i-1]) continue;
    /* When not, calculate the area of an imaginary frame */
    x=(x1[i]+x2[i-1])/2.0;
    a=(x1[i]-x2[i-1])*
      (y[i]-(y[i]-y[i-1])*(x2[i]+x1[i]-2.0*x)/(x2[i]+x1[i]-x2[i-1]-x1[i-1]));
    ti[i]+=a;
  }

  /* Copy integrals to output if required */
  if(allocated) for(i=0; i<nr; i++) newyi[i]=ti[i];

  /* Calculate 2nd integrals if required */
  if(newyii!=NULL) {
    newyii[0]=x2[0]*ti[0]/2.0;
    for(i=1; i<nr; i++)
      newyii[i]=newyii[i-1]+(x2[i]-x2[i-1])*(ti[i-1]+ti[i])/2.0;
  }
  
  /* Free memory */
  if(allocated) free((char*)ti);

  if(INTEGR_TEST) printf("out fpetintegrate()\n");
  return 0;
}
/*****************************************************************************/

/*****************************************************************************/
/** @brief Interpolate and integrate TAC to PET frames.

  It is assumed, that original data is not from framed data, but that
  the values represent the actual value at specified time point, which
  allows the integration to be calculated dot-by-dot.

  If NULL is specified for *newy, *newyi and/or *newyii, their values are
  not calculated.

  Original data and new x values must be sorted by ascending x.
  Subsequent x values can have equal values, enabling the use of step functions.
  PET frames can overlap, but interpolation may then be slower.
  Negative x (time) values can be processed.
  If necessary, the data is extrapolated assuming that 1) y[inf]=y[nr-1] and
  2) if x[0]>0 and y[0]>0 and x[0]>x[1]-x[0], then an imaginary line is drawn
  from (0,0) to (x[0],y[0]).

  @return Non-zero in case of an error.
*/
int interpolate4pet(
  /** Times of original data */
  double *x,
  /** Values of original data */
  double *y,
  /** Number of original data values */
  int nr,
  /** PET frame start times; frames may overlap */
  double *newx1,
  /** PET frame end times; frames may overlap */
  double *newx2, 
  /** Mean value during PET frame, or NULL if not needed;
   *  calculation may be faster if newyi is calculated too */
  double *newy,
  /** Integral at frame mid time, or NULL if not needed */
  double *newyi,
  /** 2nd integral at frame mid time, or NULL if not needed */
  double *newyii,
  /** Number of PET frames */
  int newnr
) {
  int ret, fi, overlap=0, zeroframe=0;
  double petx[3], pety[3], petyi[3], petyii[3], fdur;

  if(INTEGR_TEST) printf("in interpolate4pet()\n");

  /* Check for data */
  if(nr<1 || newnr<1) return 1;
  /* Check that programmer understood that outp data must've been allocated */
  if(newy==NULL && newyi==NULL && newyii==NULL) return 2;
  /* Check that input data is not totally outside the input data */
  if(newx2[newnr-1]<=x[0] || newx1[0]>=x[nr-1]) return 3;

  /* Check frame lengths, also for overlap and zero length frames */
  for(fi=0; fi<newnr; fi++) {
    /* Calculate frame length */
    fdur=newx2[fi]-newx1[fi];
    /* If frame length is <0, that is an error */
    if(fdur<0.0) return 4;
    if(fdur==0.0) zeroframe++;
    /* Overlap? */
    if(fi>0 && newx2[fi-1]>newx1[fi]) overlap++;
  }
  if(INTEGR_TEST>1) {
    printf("overlap := %d\n", overlap);
    printf("zeroframe := %d\n", zeroframe);
  }

  /* Interpolate and integrate one frame at a time, if there is:
     -overlap in frames
     -frames of zero length
     -only few interpolated frames
     -if only integrals are needed
     -if neither of integrals is needed, then there is no temp memory to
      do otherwise
  */
  if(overlap>0 || zeroframe>0 || newnr<=3 ||
     newy==NULL || (newyi==NULL && newyii==NULL) )
  {

    if(INTEGR_TEST>1) printf("frame-by-frame interpolation/integration\n");
    for(fi=0; fi<newnr; fi++) {
      /* Set frame start, middle and end times */
      petx[0]=newx1[fi]; petx[2]=newx2[fi]; petx[1]=0.5*(petx[0]+petx[2]);
      /* Calculate frame length */
      fdur=petx[2]-petx[0];
      /* If frame length is <0, that is an error */
      if(fdur<0.0) return 4;
      /* If frame length is 0, then use direct interpolation */
      if(fdur==0.0) {
        ret=interpolate(x, y, nr, petx, pety, petyi, petyii, 1);
        if(ret) return 10+ret;
        if(newy!=NULL) newy[fi]=petyi[0];
        if(newyi!=NULL) newyi[fi]=petyi[0];
        if(newyii!=NULL) newyii[fi]=petyii[0];
        continue;
      }
      /* Calculate integrals at frame start, middle, and end */
      ret=interpolate(x, y, nr, petx, NULL, petyi, petyii, 3);
      if(ret) return 20+ret;
      /* Set output integrals, if required */
      if(newyi!=NULL) newyi[fi]=petyi[1];
      if(newyii!=NULL) newyii[fi]=petyii[1];
      /* Calculate frame mean, if required */
      if(newy!=NULL) newy[fi]=(petyi[2]-petyi[0])/fdur;
    } // next frame

  } else {

    if(INTEGR_TEST>1) printf("all-frames-at-once interpolation/integration\n");
    /* Set temp array */
    double *tp; if(newyii!=NULL) tp=newyii; else tp=newyi;
    /* Integrate at frame start times */
    ret=interpolate(x, y, nr, newx1, NULL, tp, NULL, newnr);
    if(ret) return 10+ret;
    /* Integrate at frame end times */
    ret=interpolate(x, y, nr, newx2, NULL, newy, NULL, newnr);
    if(ret) return 10+ret;
    /* Calculate average frame value */
    for(fi=0; fi<newnr; fi++)
      newy[fi]=(newy[fi]-tp[fi])/(newx2[fi]-newx1[fi]);
    /* Calculate integrals */
    if(newyi!=NULL || newyii!=NULL) {
      for(fi=0; fi<newnr; fi++) newx1[fi]+=0.5*(newx2[fi]-newx1[fi]);
      ret=interpolate(x, y, nr, newx1, NULL, newyi, newyii, newnr);
      if(ret) return 10+ret;
      for(fi=0; fi<newnr; fi++) newx1[fi]-=(newx2[fi]-newx1[fi]);
    }
  }

  if(INTEGR_TEST) printf("out interpolate4pet()\n");
  return 0;
}

/** @brief Interpolate and integrate TAC to PET frames.
  Float version of finterpolate4pet().

  It is assumed, that original data is not from framed data, but that
  the values represent the actual value at specified time point, which
  allows the integration to be calculated dot-by-dot.

  If NULL is specified for *newy, *newyi and/or *newyii, their values are
  not calculated.

  Original data and new x values must be sorted by ascending x.
  Subsequent x values can have equal values, enabling the use of step functions.
  PET frames can overlap, but interpolation may then be slower.
  Negative x (time) values can be processed.
  If necessary, the data is extrapolated assuming that 1) y[inf]=y[nr-1] and
  2) if x[0]>0 and y[0]>0 and x[0]>x[1]-x[0], then an imaginary line is drawn
  from (0,0) to (x[0],y[0]).

  @return Non-zero in case of an error.
*/
int finterpolate4pet(
  /** Times of original data */
  float *x,
  /** Values of original data */
  float *y,
  /** Number of original data values */
  int nr,
  /** PET frame start times; frames may overlap */
  float *newx1,
  /** PET frame end times; frames may overlap */
  float *newx2, 
  /** Mean value during PET frame, or NULL if not needed;
   *  calculation may be faster if newyi is calculated too */
  float *newy,
  /** Integral at frame mid time, or NULL if not needed */
  float *newyi,
  /** 2nd integral at frame mid time, or NULL if not needed */
  float *newyii,
  /** Number of PET frames */
  int newnr
) {
  int ret, fi, overlap=0, zeroframe=0;
  float petx[3], pety[3], petyi[3], petyii[3], fdur;

  if(INTEGR_TEST) printf("in finterpolate4pet()\n");

  /* Check for data */
  if(nr<1 || newnr<1) return 1;
  /* Check that programmer understood that outp data must've been allocated */
  if(newy==NULL && newyi==NULL && newyii==NULL) return 2;
  /* Check that input data is not totally outside the input data */
  if(newx2[newnr-1]<=x[0] || newx1[0]>=x[nr-1]) return 3;

  /* Check frame lengths, also for overlap and zero length frames */
  for(fi=0; fi<newnr; fi++) {
    /* Calculate frame length */
    fdur=newx2[fi]-newx1[fi];
    /* If frame length is <0, that is an error */
    if(fdur<0.0) return 4;
    if(fdur==0.0) zeroframe++;
    /* Overlap? */
    if(fi>0 && newx2[fi-1]>newx1[fi]) overlap++;
  }
  if(INTEGR_TEST>1) {
    printf("overlap := %d\n", overlap);
    printf("zeroframe := %d\n", zeroframe);
  }

  /* Interpolate and integrate one frame at a time, if there is:
     -overlap in frames
     -frames of zero length
     -only few interpolated frames
     -if only integrals are needed
     -if neither of integrals is needed, then there is no temp memory to
      do otherwise
  */
  if(overlap>0 || zeroframe>0 || newnr<=3 ||
     newy==NULL || (newyi==NULL && newyii==NULL) )
  {

    if(INTEGR_TEST>1) printf("frame-by-frame interpolation/integration\n");
    for(fi=0; fi<newnr; fi++) {
      /* Set frame start, middle and end times */
      petx[0]=newx1[fi]; petx[2]=newx2[fi]; petx[1]=0.5*(petx[0]+petx[2]);
      /* Calculate frame length */
      fdur=petx[2]-petx[0];
      /* If frame length is <0, that is an error */
      if(fdur<0.0) return 4;
      /* If frame length is 0, then use direct interpolation */
      if(fdur==0.0) {
        ret=finterpolate(x, y, nr, petx, pety, petyi, petyii, 1);
        if(ret) return 10+ret;
        if(newy!=NULL) newy[fi]=petyi[0];
        if(newyi!=NULL) newyi[fi]=petyi[0];
        if(newyii!=NULL) newyii[fi]=petyii[0];
        continue;
      }
      /* Calculate integrals at frame start, middle, and end */
      ret=finterpolate(x, y, nr, petx, NULL, petyi, petyii, 3);
      if(ret) return 20+ret;
      /* Set output integrals, if required */
      if(newyi!=NULL) newyi[fi]=petyi[1];
      if(newyii!=NULL) newyii[fi]=petyii[1];
      /* Calculate frame mean, if required */
      if(newy!=NULL) newy[fi]=(petyi[2]-petyi[0])/fdur;
    } // next frame

  } else {

    if(INTEGR_TEST>1) printf("all-frames-at-once interpolation/integration\n");
    /* Set temp array */
    float *tp; if(newyii!=NULL) tp=newyii; else tp=newyi;
    /* Integrate at frame start times */
    ret=finterpolate(x, y, nr, newx1, NULL, tp, NULL, newnr);
    if(ret) return 10+ret;
    /* Integrate at frame end times */
    ret=finterpolate(x, y, nr, newx2, NULL, newy, NULL, newnr);
    if(ret) return 10+ret;
    /* Calculate average frame value */
    for(fi=0; fi<newnr; fi++)
      newy[fi]=(newy[fi]-tp[fi])/(newx2[fi]-newx1[fi]);
    /* Calculate integrals */
    if(newyi!=NULL || newyii!=NULL) {
      for(fi=0; fi<newnr; fi++) newx1[fi]+=0.5*(newx2[fi]-newx1[fi]);
      ret=finterpolate(x, y, nr, newx1, NULL, newyi, newyii, newnr);
      if(ret) return 10+ret;
      for(fi=0; fi<newnr; fi++) newx1[fi]-=(newx2[fi]-newx1[fi]);
    }

  }

  if(INTEGR_TEST) printf("out finterpolate4pet()\n");
  return 0;
}
/*****************************************************************************/

/*****************************************************************************/
/** @brief Integrate PET TAC data to frame mid times.

    Any of output arrays may be set to NULL if that is not needed.
    Frames must be in ascending time order. Gaps and small overlap are allowed.
    If x1[0]>0 and x1[0]<=x2[0]-x1[0], then an imaginary line is 
    drawn from (0,0) to (x[0],y[0]).
\return Returns 0 if ok.
*/
int petintegral(
  /** frame start times */
  double *x1,
  /** frame end times */
  double *x2,
  /** avg value during frame */
  double *y,
  /** number of frames */
  int nr,
  /** integrals at frame mid time */
  double *ie,
  /** 2nd integrals at frame mid time */
  double *iie
) {
  int i;
  double x, last_x, last_x2, last_y, last_integral, box_integral, half_integral;
  double gap_integral, integral, integral2, frame_len, xdist, s;

  if(INTEGR_TEST) printf("in petintegral()\n");
  /* Check for data */
  if(nr<1 || nr<1) return(1);
  /* Check that programmer understood that output must've been allocated */
  if(ie==NULL && iie==NULL) return(2);

  /* Initiate values to zero */
  last_x=last_x2=last_y=last_integral=0.0;
  box_integral=integral=integral2=frame_len=0.0;

  for(i=0; i<nr; i++) {
    frame_len=x2[i]-x1[i]; if(frame_len<0.0) return(5);
    x=0.5*(x1[i]+x2[i]); xdist=x-last_x; 
    if(last_x>0.0 && xdist<=0.0) return(6);
    if(x<0) {
      if(ie!=NULL) ie[i]=integral;
      if(iie!=NULL) iie[i]=integral2;
      continue;
    }
    s=(y[i]-last_y)/xdist; /* slope */
    /*If there is a big gap in the beginning it is eliminated */
    if(i==0)if(x1[0]>x2[0]-x1[0]){last_x2=last_x=x1[0];}
    /*Integral of a possible gap between frames*/
    gap_integral=(x1[i]-last_x2)*(last_y+s*((last_x2+x1[i])/2.0-last_x)); 
    /*Integral from the beginning of the frame to the middle of the frame */
    half_integral=(x-x1[i])*(last_y+s*((x1[i]+x)/2.0-last_x));
    integral=box_integral+gap_integral+half_integral; 
    /* half_integral is not added to box because it is more accurate to 
       increment integrals of full frames */
    box_integral+=gap_integral+frame_len*y[i];
    integral2+=xdist*(integral+last_integral)*0.5;
    if(ie!=NULL) ie[i]=integral;
    if(iie!=NULL) iie[i]=integral2;
    last_x=x; last_x2=x2[i];
    last_y=y[i]; last_integral=integral;
  }

  if(INTEGR_TEST) printf("out petintegral()\n");
  return(0);
}

/** @brief Integrate PET TAC data to frame mid times.
    Float version of petintegral().

    Any of output arrays may be set to NULL if that is not needed.
    Frames must be in ascending time order. Gaps and small overlap are allowed.
    If x1[0]>0 and x1[0]<=x2[0]-x1[0], then an imaginary line is 
    drawn from (0,0) to (x[0],y[0]).
\return Returns 0 if ok.
*/
int fpetintegral(
  /** frame start times */
  float *x1,
  /** frame end times */
  float *x2,
  /** avg value during frame */
  float *y,
  /** number of frames */
  int nr,
  /** integrals at frame mid time */
  float *ie,
  /** 2nd integrals at frame mid time */
  float *iie
) {
  int i;
  float x, last_x, last_x2, last_y, last_integral, box_integral, half_integral;
  float gap_integral, integral, integral2, frame_len, xdist, s;

  if(INTEGR_TEST) printf("in fpetintegral()\n");
  /* Check for data */
  if(nr<1 || nr<1) return(1);
  /* Check that programmer understood that output must've been allocated */
  if(ie==NULL && iie==NULL) return(2);

  /* Initiate values to zero */
  last_x=last_x2=last_y=last_integral=0.0;
  box_integral=gap_integral=integral=integral2=frame_len=0.0;

  for(i=0; i<nr; i++) {
    frame_len=x2[i]-x1[i]; if(frame_len<0.0) return(5);
    x=0.5*(x1[i]+x2[i]); xdist=x-last_x; if(last_x>0.0 && xdist<=0.0) return(6);
    if(x<0) {
      if(ie!=NULL) ie[i]=integral;
      if(iie!=NULL) iie[i]=integral2;
      continue;
    }
    s=(y[i]-last_y)/xdist; /* slope */
    /*If there is a big gap in the beginning it is eliminated */
    if(i==0)if(x1[0]>x2[0]-x1[0]){last_x2=last_x=x1[0];}
    /*Integral of a possible gap between frames*/
    gap_integral=(x1[i]-last_x2)*(last_y+s*((last_x2+x1[i])/2.0-last_x));
    /*Integral from the beginning of the frame i to the middle of the frame */
    half_integral=(x-x1[i])*(last_y+s*((x1[i]+x)/2.0-last_x));
    integral=box_integral+gap_integral+half_integral; 
    /* half_integral is not added to box because it is more accurate to 
       sum integrals of full frames */
    box_integral+=gap_integral+frame_len*y[i];
    integral2+=xdist*(integral+last_integral)*0.5;
    if(ie!=NULL) ie[i]=integral;
    if(iie!=NULL) iie[i]=integral2;
    last_x=x; last_x2=x2[i];
    last_y=y[i]; last_integral=integral;
  }

  if(INTEGR_TEST) printf("out fpetintegral()\n");
  return(0);
}
/*****************************************************************************/

/*****************************************************************************/
/** @brief Integrate PET TAC data to frame end times.
 *
    Any of output arrays may be set to NULL if that is not needed.
    Frames must be in ascending time order. Gaps and small overlap are allowed.
    If x1[0]>0 and x1[0]<=x2[0]-x1[0], then an imaginary line is 
    drawn from (0,0) to (x[0],y[0]).
\return Returns 0 if ok.
*/
int petintegrate2fe(
  /** frame start times */
  double *x1,
  /** frame end times */
  double *x2,
  /** avg value during frame */
  double *y,
  /** number of frames */
  int nr,
  /** values at frame end time */
  double *e,
  /** integrals at frame end time */
  double *ie,
  /** 2nd integrals at frame end time */
  double *iie
) {
  int i;
  double x, last_x, last_x2, last_y, last_integral;
  double value, integral, integral2, frame_len, xdist, s;

  if(INTEGR_TEST) printf("in petintegrate2fe()\n");
  /* Check for data */
  if(nr<1 || nr<1) return(1);
  /* Check that programmer understood that output must've been allocated */
  if(e==NULL && ie==NULL && iie==NULL) return(2);

  /* Initiate values to zero */
  last_x=last_x2=last_y=last_integral=value=integral=integral2=frame_len=s=0.0;

  for(i=0; i<nr; i++) {
    frame_len=x2[i]-x1[i]; if(frame_len<0.0) return(5);
    x=0.5*(x1[i]+x2[i]); xdist=x-last_x; if(last_x>0.0 && xdist<=0.0) return(6);
    if(x<0) {
      if(e!=NULL) e[i]=value;
      if(ie!=NULL) ie[i]=integral;
      if(iie!=NULL) iie[i]=integral2;
      continue;
    }
    s=(y[i]-last_y)/xdist; /* slope between x[i-1] and x[i] */
    /* If there is a big gap in the beginning, it is eliminated */
    if(i==0 && x1[0]>x2[0]-x1[0]) {last_x2=last_x=x1[0];}
    integral+=(x1[i]-last_x2)*(last_y+s*((last_x2+x1[i])/2.0-last_x)); /*gap*/
    integral+=frame_len*y[i];
    integral2+=(x2[i]-last_x2)*(integral+last_integral)*0.5;
    if(e!=NULL && i>0) {
      value=last_y+s*(last_x2-last_x); /* value at previous frame end */
      e[i-1]=value;
    }
    if(ie!=NULL) ie[i]=integral;
    if(iie!=NULL) iie[i]=integral2;
    last_x=x; last_x2=x2[i]; last_y=y[i]; last_integral=integral;
  }
  if(e!=NULL) {
    value=last_y+s*(last_x2-last_x); /* Value for the last frame */
    e[i-1]=value;
  }

  if(INTEGR_TEST) printf("out petintegrate2fe()\n");
  return(0);
}

/** @brief Integrate PET TAC data to frame end times.
    Float version of petintegrate2fe().

    Any of output arrays may be set to NULL if that is not needed.
    Frames must be in ascending time order. Gaps and small overlap are allowed.
    If x1[0]>0 and x1[0]<=x2[0]-x1[0], then an imaginary line is 
    drawn from (0,0) to (x[0],y[0]).
\return Returns 0 if ok.
*/
int fpetintegrate2fe(
  /** frame start times */
  float *x1,
  /** frame end times */
  float *x2,
  /** avg value during frame */
  float *y,
  /** number of frames */
  int nr,
  /** values at frame end time */
  float *e,
  /** integrals at frame end time */
  float *ie,
  /** 2nd integrals at frame end time */
  float *iie
) {
  int i;
  float x, last_x, last_x2, last_y, last_integral;
  float value, integral, integral2, frame_len, xdist, s;

  if(INTEGR_TEST) printf("in fpetintegrate2fe()\n");
  /* Check for data */
  if(nr<1 || nr<1) return(1);
  /* Check that programmer understood that output must've been allocated */
  if(e==NULL && ie==NULL && iie==NULL) return(2);

  /* Initiate values to zero */
  last_x=last_x2=last_y=last_integral=value=integral=integral2=frame_len=s=0.0;

  for(i=0; i<nr; i++) {
    frame_len=x2[i]-x1[i]; if(frame_len<0.0) return(5);
    x=0.5*(x1[i]+x2[i]); xdist=x-last_x; if(last_x>0.0 && xdist<=0.0) return(6);
    if(x<0) {
      if(e!=NULL) e[i]=value;
      if(ie!=NULL) ie[i]=integral;
      if(iie!=NULL) iie[i]=integral2;
      continue;
    }
    s=(y[i]-last_y)/xdist; /* slope between x[i-1] and x[i] */
    /* If there is a big gap in the beginning, it is eliminated */
    if(i==0 && x1[0]>x2[0]-x1[0]) {last_x2=last_x=x1[0];}
    integral+=(x1[i]-last_x2)*(last_y+s*((last_x2+x1[i])/2.0-last_x)); /*gap*/
    integral+=frame_len*y[i];
    integral2+=(x2[i]-last_x2)*(integral+last_integral)*0.5;
    if(e!=NULL && i>0) {
      value=last_y+s*(last_x2-last_x); /* value at previous frame end */
      e[i-1]=value;
    }
    if(ie!=NULL) ie[i]=integral;
    if(iie!=NULL) iie[i]=integral2;
    last_x=x; last_x2=x2[i]; last_y=y[i]; last_integral=integral;
  }
  if(e!=NULL) {
    value=last_y+s*(last_x2-last_x); /* Value for the last frame */
    e[i-1]=value;
  }

  if(INTEGR_TEST) printf("out fpetintegrate2fe()\n");
  return(0);
}
/*****************************************************************************/

/*****************************************************************************/
