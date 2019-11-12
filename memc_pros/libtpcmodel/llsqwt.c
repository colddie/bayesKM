/// @file llsqwt.c
/// @author Vesa Oikonen
/// @brief Linear least-squares fit with errors in both coordinates.
///
/******************************************************************************/
#include "libtpcmodel.h"
/******************************************************************************/
/* Local function definitions */
/// @cond
int _medianline_cmp(const void *e1, const void *e2);
/// @endcond
/******************************************************************************/

/******************************************************************************/
/** Iterative method for linear least-squares fit with errors in both
    coordinates. This function is fully based on article [3].

  For n data-point pairs (x[i], y[i]) each point has its own weighting factors
  in (wx[i], wy[i]). This routine finds the values of the parameters m (slope)
  and c (intercept, ic) that yield the "best-fit" of the model equation
  Y = mX + c to the data, where X and Y are the predicted or calculated values
  of the data points.

  Weighting factors wx and wy must be assigned as the inverses of the variances
  or squares of the measurement uncertainties (SDs), i.e. w[i]=1/(sd[i])^2

  If true weights are unknown but yet the relative weights are correct,
  the slope, intercept and residuals (WSS) will be correct.
  The applied term S/(N-2) makes also the estimate of sigma (sd) of slope
  less dependent on the scaling of weights. The sigmas are not exact, since
  only the lowest-order terms in Taylor-series expansion are incorporated;
  anyhow sigmas are more accurate than the ones based on York algorithm.

  One or more data points can be excluded from the fit by setting either x or y
  weight to 0.

  References:
  1. York, D. Linear-squares fitting of a straight line. Can J Phys. 1966;44:1079-1086.
  2. Lybanon, M. A better least squares method when both variables have uncertainties. 
     Am J Phys. 1984;52:22-26 and 276-278.
  3. Reed BC. Linear least-squares fits with errors in both coordinates. II:
     Comments on parameter variances. Am J Phys. 1992;60:59-62.

  @sa llsqperp, best_llsqwt, medianline
  @return If successful, function returns value 0.
*/

int LLSQWT_TEST;
int llsqwt(
  /** coordinates of data points (of dimension n). */
  double *x,
  /** coordinates of data points (of dimension n). */
  double *y,
  /** number of data points. */
  int n,
  /** weighting factors in x. */
  double *wx,
  /** weighting factors in y. */
  double *wy,
  /** allowed tolerance in slope estimation. */
  double tol,
  /* Input/output */
  /** work vector (of dimension n); effective weights w[i] are returned in it. */
  double *w,
  /* Output */
  /** Estimated intercept. */
  double *ic,
  /** Estimated slope. */
  double *slope,
  /** sqrt(WSS)/wsum of the residuals. */
  double *nwss,
  /** expected sd of intercept at calculated points;
      If NULL, then not calculated and variable is left unchanged. */
  double *sic,
  /** Expected sd of slope at calculated points;
      If NULL, then not calculated and variable is left unchanged. */
  double *sslope,
  /** Estimated data points (X,y);
      If NULL, then not calculated and variable is left unchanged. */
  double *cx,
  /** Estimated data points (x,Y);
      If NULL, then not calculated and variable is left unchanged. */
  double *cy
) {
  int i, j, niter=0, nn, calcVar=1;
  double c, m, f, xb=0.0, yb=0.0, qa=0.0, qb=0.0, qc, ss=0.0, s1, s2, wsum=0.0;
  double xsum=0.0, x2sum=0.0, ysum=0.0, xysum=0.0, delta, discr, sqdis;
  double m_1, m_2, bcont, m2=0.0, w2;
  double u, v, AA, BB, CC, DD, EE, FF, GG, HH, JJ;
  double varc, varm, dmx, dmy, dcx, dcy;


  /*
   *  Lets look what we got for arguments
   */
  if(LLSQWT_TEST) {fprintf(stdout, "llsqwt()\n"); fflush(stdout);}
  /* Check that there is some data */
  if(n<2) return(1);
  if(LLSQWT_TEST>3) {
    for(i=0; i<n; i++) printf("%e +- %e    %e +- %e\n",
      x[i], wx[i], y[i], wy[i]);
  }
  if(tol<1.0e-100) return(1);
  if(w==NULL) return(1);
  /* Check if variances and fitted data will be calculated */
  if(sic==NULL || sslope==NULL || cx==NULL || cy==NULL) calcVar=0;
  else calcVar=1;

  /* If only 2 datapoints */
  if(n==2) {
    f=x[1]-x[0];
    if(f==0.0) {*slope=*ic=*sic=*sslope=*nwss=w[0]=w[1]=0.0; return(0);}
    *slope=(y[1]-y[0])/f; *ic=y[0]-(*slope)*x[0];
    *sic=*sslope=0.; w[0]=1.0; w[1]=1.0; *nwss=0.0;
    return(0);
  }

  /*
   *  Fit the LLSQ line
   */

  /* First estimation of the slope and intercept by unweighted regression */
  for(i=0; i<n; i++) {
    xsum+=x[i]; ysum+=y[i]; x2sum+=x[i]*x[i]; xysum+=x[i]*y[i];
  }
  delta=(double)n*x2sum - xsum*xsum;
  /* Initial guesses of the slope and intercept */
  if(delta==0.0) {
    if(LLSQWT_TEST) printf("x axis values contain only zeroes.\n");
    *slope=*ic=*sic=*sslope=*nwss=w[0]=w[1]=0.0;
    return(0);
  }
  m=((double)n*xysum - xsum*ysum)/delta;
  c=(x2sum*ysum - xsum*xysum)/delta;
  if(LLSQWT_TEST) printf("initial guesses: a=%e b=%e\n", c, m);

  /* Begin the iterations */
  bcont=m+2.0*tol;
  while(fabs(m-bcont)>tol && niter<20) {
    if(LLSQWT_TEST>2) printf(" %d. iteration, improvement=%g, tol=%g\n", niter, fabs(m-bcont), tol);
    bcont=m; niter++; 
    /* Compute the weighting factors */
    m2=m*m;
    for(i=0, nn=0; i<n; i++) {
      if(wx[i]<=0.0 || wy[i]<=0.0) w[i]=0.0;
      else {w[i]=wx[i]*wy[i]/(m2*wy[i]+wx[i]); nn++;}
      //printf("wx[i]=%g wy[i]=%g -> w[%d]=%g\n", wx[i], wy[i], i, w[i]);
    }
    if(nn<2) {
      if(LLSQWT_TEST) printf("less than two points with weight > 0.\n");
      *slope=*ic=*sic=*sslope=*nwss=0.0;
      return(0);   
    }
    /* Compute the barycentre coordinates */
    for(i=0, xb=yb=wsum=0.0; i<n; i++) {xb+=w[i]*x[i]; yb+=w[i]*y[i]; wsum+=w[i];}
    if(wsum<=0.0) return(2);
    xb/=wsum; yb/=wsum;
    if(LLSQWT_TEST>2) printf("barycentre: xb=%g yb=%g\n", xb, yb);
    /* Estimate the slope as either of the roots of the quadratic */
    /* compute the coefficients of the quadratic */
    for(i=0, qa=qb=qc=0.0; i<n; i++) if(w[i]>0.0) {
      u=x[i]-xb; v=y[i]-yb; w2=w[i]*w[i];
      qa += w2*u*v/wx[i];
      qb += w2*(u*u/wy[i] - v*v/wx[i]);
      qc += -w2*u*v/wy[i];
    }
    if(LLSQWT_TEST>2) printf("quadratic coefs: qa=%g qb=%g qc=%g\n", qa, qb, qc);
    if(qa==0.0) {
      m=0.0;
      /* Calculate WSS */
      for(i=0, ss=0.0; i<n; i++) {f=v=y[i]-yb; ss+=w[i]*f*f;}
    } else if(qa==1.0) { /* check if quadratic reduces to a linear form */
      m=-qc/qb;
      /* Calculate WSS */
      for(i=0, ss=0.0; i<n; i++) {u=x[i]-xb; v=y[i]-yb; f=v-m*u; ss+=w[i]*f*f;}
    } else {
      /* discriminant of quadratic */
      discr=qb*qb-4.0*qa*qc; if(discr<=0.0) sqdis=0.0; else sqdis=sqrt(discr);
      /* compute the two solutions of quadratic */
      m_1=(-qb+sqdis)/(2.0*qa); m_2=(-qb-sqdis)/(2.0*qa);
      /* Calculate WSS for both solutions */
      for(i=0, s1=s2=0.0; i<n; i++) {
        u=x[i]-xb; v=y[i]-yb;
        f=v-m_1*u; s1+=w[i]*f*f;
        f=v-m_2*u; s2+=w[i]*f*f;
      }
      /* choose the solution with lower WSS */
      if(s1<=s2) {m=m_1; ss=s1;} else {m=m_2; ss=s2;}
    }
    /* Calculate the intercept */
    c = yb - m*xb;
    if(LLSQWT_TEST>2) printf("iterated parameters: c=%e m=%e\n", c, m);

  } /* end of iteration loop */
  *ic=c; *slope=m; *nwss=sqrt(ss)/wsum; /* Set function return values */
  if(LLSQWT_TEST) {
    fprintf(stdout, "c=%14.7e m=%14.7e ss=%14.7e niter=%d\n", c, m, ss, niter);
  }

  /* Return here, if variances are not to be calculated */
  if(!calcVar) return(0);

  /*
   *  Calculate the fitted line
   */
  for(i=0; i<n; i++) {
    if(w[i]>0.0) {
      f=w[i]*(c + m*x[i] - y[i]); /* Lagrangian multiplier */ 
      cx[i]=x[i]-f*m/wx[i]; cy[i]=y[i]+f/wy[i];
    } else {
      cx[i]=x[i]; cy[i]=y[i];
    }
  }

  /*
   *  Estimate the variances of the parameters (Reed 1992)
   */
  /* varm = var of slope ; varc = var of intercept  */

  /* Use the true nr of data points, i.e. exclude the ones with zero weight */
  for(i=nn=0; i<n; i++) if(w[i]>0.0) nn++;
  if(nn<3) {*sslope=*sic=0.0; return(0);}

  /* Compute the barycentre coordinates again from computed data points */
  for(i=0, xb=yb=wsum=0.0; i<n; i++) {xb+=w[i]*cx[i]; yb+=w[i]*cy[i]; wsum+=w[i];}
  if(wsum<=0.0) return(2);
  xb/=wsum; yb/=wsum;
  if(LLSQWT_TEST) printf("barycentre: xb=%g yb=%g\n", xb, yb);

  /* common factors */
  HH=JJ=qa=qb=qc=0.0;
  for(i=0; i<n; i++) if(w[i]>0.0) {
    u=cx[i]-xb; v=cy[i]-yb; w2=w[i]*w[i];
    qa += w2*u*v/wx[i];
    qb += w2*(u*u/wy[i] - v*v/wx[i]);
    qc += -w2*u*v/wy[i];
    HH += w2*v/wx[i];
    JJ += w2*u/wx[i];
  }
  HH*=-2.0*m/wsum; JJ*=-2.0*m/wsum;
  if(LLSQWT_TEST>3)
    printf("quadratic coefs: qa=%g qb=%g qc=%g ; HH=%g JJ=%g\n", qa, qb, qc, HH, JJ);
  AA=BB=CC=0.0;
  for(i=0; i<n; i++) if(w[i]>0.0) {
    u=cx[i]-xb; v=cy[i]-yb; w2=w[i]*w[i];
    AA += w[i]*w[i]*w[i]*u*v / (wx[i]*wx[i]);
    BB -= w2 * (4.0*m*(w[i]/wx[i])*(u*u/wy[i]-v*v/wx[i]) - 2.0*v*HH/wx[i] + 2.0*u*JJ/wy[i]);
    CC -= (w2/wy[i]) * (4.0*m*w[i]*u*v/wx[i] + v*JJ + u*HH);
  }
  if(m!=0) AA = 4.0*m*AA - wsum*HH*JJ/m; else AA=0.0;
  if(LLSQWT_TEST>3) printf("AA=%g BB=%g CC=%g\n", AA, BB, CC);
  /* sigmas */
  varc=varm=0.0;
  for(j=0; j<n; j++) if(w[j]>0.0) {
    /* factors for this j */
    DD=EE=FF=GG=0.0;
    for(i=0; i<n; i++) if(w[i]>0.0) {
      u=cx[i]-xb; v=cy[i]-yb; w2=w[i]*w[i];
      if(i==j) delta=1.0; else delta=0.0; /* Kronecker delta */
      f = delta - w[j]/wsum;
      DD += (w2*v/wx[i])*f;
      EE += (w2*u/wy[i])*f;
      FF += (w2*v/wy[i])*f;
      GG += (w2*u/wx[i])*f;
    }
    EE*=2.0;
    /* derivatives at jth data point */
    f = (2.0*m*qa + qb - AA*m2 + BB*m - CC); m2=m*m;
    dmx = -(m2*DD + m*EE - FF) / f;
    dmy = -(m2*GG - 2.0*m*DD - EE/2.0) / f;
    dcx = (HH - m*JJ - xb)*dmx - m*w[j]/wsum;
    dcy = (HH - m*JJ - xb)*dmy + w[j]/wsum;
    /* sum terms for sigmas */
    varm += dmy*dmy/wy[j] + dmx*dmx/wx[j];
    varc += dcy*dcy/wy[j] + dcx*dcx/wx[j];
    if(LLSQWT_TEST>3)
      printf("DD=%g EE=%g FF=%g GG=%g dmx=%g dmy=%g dcx=%g dcy=%g\n",
             DD, EE, FF, GG, dmx, dmy, dcx, dcy);
  }
  varm *= ss/(double)(nn-2);
  varc *= ss/(double)(nn-2);
  if(LLSQWT_TEST>3) printf("varm=%g varc=%g\n", varm, varc);
  *sslope = sqrt(varm);
  *sic = sqrt(varc);
  if(LLSQWT_TEST) {fprintf(stdout, "sslope=%14.7e sic=%14.7e\n", *sslope, *sic);}
  return(0);
}
/******************************************************************************/

/******************************************************************************/
/** Finds the best least-squares line to (x,y)-data, leaving points out either
    from the beginning (mode=0) or from the end (mode=1).
    @sa llsqwt, llsqperp, medianline
    @return Returns 0, if ok.
 */
int best_llsqwt(
  /** Plot x axis values. */
  double *x,
  /** Plot y axis values. */
  double *y,
  /** Weighting factors for x. */
  double *wx,
  /** Weighting factors for y. */
  double *wy,
  /** Nr of plot data points. */
  int nr,
  /** Min nr of points to use in the fit; must be >=4. */
  int min_nr,
  /** Leave out points from beginning (0) or from end (1). */
  int mode,
  /** Slope is returned in here. */
  double *slope,
  /** Y axis intercept is returned in here. */
  double *ic,
  /** sqrt(WSS)/wsum is returned here. */
  double *nwss,
  /** Expected sd of slope at calculated points. */
  double *sslope,
  /** Expected sd of intercept at calculated points. */
  double *sic,
  /** Calculated x data points. */
  double *cx,
  /** Calculated y data points. */
  double *cy,
  /** Number of points in the best fit (incl data with w=0). */
  int *bnr
) {
  int from, to, ret, from_min, to_min;
  double *w, lic, lslope, lnwss=1.0E+100, nwss_min=1.0E+100;

  /* Check the data */
  if(x==NULL || y==NULL || nr<min_nr || nr<2) return(1);
  /* Check parameters */
  if(min_nr<4 || (mode!=0 && mode!=1)) return(2);

  /* Make room for work vector */
  w=(double*)malloc(nr*sizeof(double)); if(w==NULL) return(3);

  /* Search the plot range that gives the best nwss */
  nwss_min=9.99E+99; from_min=to_min=-1;
  if(mode==0) { /* Leave out plot data from the beginning */
    for(from=0, to=nr-1; from<nr-min_nr; from++) {
      ret=llsqwt(x+from, y+from, (to-from)+1, wx+from, wy+from, 1.0E-10, w,
                 &lic, &lslope, &lnwss, NULL, NULL, NULL, NULL);
      /*lnwss/=(double)((to-from)+1);*/
      if(LLSQWT_TEST) {
        printf("  range: %d-%d ; nwss=%g ; min=%g ; ret=%d\n", from, to, lnwss, nwss_min, ret);
      }
      if(ret==0 && lnwss<nwss_min) {nwss_min=lnwss; from_min=from; to_min=to;}
    }
  } else { /* Leave out plot data from the end */
    for(from=0, to=min_nr-1; to<nr; to++) {
      ret=llsqwt(x+from, y+from, (to-from)+1, wx+from, wy+from, 1.0E-10, w,
                 &lic, &lslope, &lnwss, NULL, NULL, NULL, NULL);
      /*lnwss/=(double)((to-from)+1);*/
      if(LLSQWT_TEST) {
        printf("  range: %d-%d ; nwss=%g ; min=%g ; ret=%d\n", from, to, lnwss, nwss_min, ret);
      }
      if(ret==0 && lnwss<nwss_min) {nwss_min=lnwss; from_min=from; to_min=to;}
    }
  }
  if(from_min<0) {free(w); return(4);}

  /* Run llsqwt() again with that range, now with better resolution      */
  /* and this time compute also SD's                                     */
  from=from_min; to=to_min;
  ret=llsqwt(x+from, y+from, (to-from)+1, wx+from, wy+from, 1.0E-15, w,
             ic, slope, nwss, sic, sslope, cx+from, cy+from);
  free(w);
  if(ret) return(5);
  *bnr=(to-from)+1;

  return(0);
}
/******************************************************************************/

/******************************************************************************/
/** Simple non-iterative perpendicular line fitting.
    This function is fully based on the article [1].
    
    References:
    1. Varga J & Szabo Z. Modified regression model for the Logan plot.
       J Cereb Blood Flow Metab 2002; 22:240-244.

    @sa llsqperp3, medianline, llsqwt, mtga_best_perp
    @return If successful, function returns value 0.
 */
int llsqperp(
  /** Coordinates of data points (dimension nr). */
  double *x,
  /** Coordinates of data points (dimension nr). */
  double *y,
  /** Number of data points. */
  int nr,
  /** Estimated slope. */
  double *slope,
  /** Estimated intercept. */
  double *ic,
  /** Sum of squared distances / nr. */
  double *ssd
) {
  int i, rnr;
  double qxx, qyy, qxy, mx, my, a, b, c, d, m1, m2, ssd1, ssd2;

  if(LLSQWT_TEST) {fprintf(stdout, "llsqperp()\n"); fflush(stdout);}
  /* Check the data */
  if(nr<2 || x==NULL || y==NULL) return(1);
  /* Calculate the means */
  for(i=0, mx=my=0; i<nr; i++) {mx+=x[i]; my+=y[i];}
  mx/=(double)nr; my/=(double)nr;
  /* Calculate the Q's */
  for(i=0, qxx=qyy=qxy=0; i<nr; i++) {
    a=x[i]-mx; b=y[i]-my; qxx+=a*a; qyy+=b*b; qxy+=a*b;
  }
  if(qxx<1.0E-100 || qyy<1.0E-100) return(2);
  /* Calculate the slope as real roots of quadratic equation */
  a=qxy; b=qxx-qyy; c=-qxy;
  rnr=quadratic(a, b, c, &m1, &m2);
  if(LLSQWT_TEST) {
    fprintf(stdout, "%d quadratic roots", rnr);
    if(rnr>0) fprintf(stdout, " %g", m1);
    if(rnr>1) fprintf(stdout, " %g", m2);
    fprintf(stdout, " ; traditional slope %g\n", qxy/qxx);
  }
  if(rnr==0) return(3);
  /* Calculate the sum of squared distances for the first root */
  a=m1; b=-1; c=my-m1*mx;
  for(i=0, ssd1=0; i<nr; i++) {
    /* calculate the distance from point (x[i],y[i]) to line ax+by+c=0 */
    //d=(a*x[i]+b*y[i]+c)/sqrt(a*a+b*b);
    d=(a*x[i]+b*y[i]+c)/hypot(a, b);
    ssd1+=d*d;
  }
  /* Also to the 2nd root, if there is one */
  if(rnr==2) {
    a=m2; b=-1; c=my-m2*mx;
    for(i=0, ssd2=0; i<nr; i++) {
      /* calculate the distance from point (x[i],y[i]) to line ax+by+c=0 */
      //d=(a*x[i]+b*y[i]+c)/sqrt(a*a+b*b);
      d=(a*x[i]+b*y[i]+c)/hypot(a, b);
      ssd2+=d*d;
    }
  } else ssd2=ssd1;
  /* If there were 2 roots, select the one with smaller ssd */
  if(rnr==2 && ssd2<ssd1) {ssd1=ssd2; m1=m2;}
  /* Set the slope and intercept */
  *slope=m1; *ic=my-m1*mx; *ssd=ssd1/(double)nr;

  return(0);
}
/******************************************************************************/

/******************************************************************************/
/** Simple non-iterative perpendicular line fitting. 
    This version of the function accepts data that contains NaN's.
    @sa llsqperp, medianline
    @return If successful, function returns value 0.
*/
int llsqperp3(
  /** Coordinates of data points (dimension nr). */
  double *x,
  /** Coordinates of data points (dimension nr). */
  double *y,
  /** Number of data points. */
  int nr,
  /** Estimated slope. */
  double *slope,
  /** Estimated intercept. */
  double *ic,
  /** Sum of squared distances / nr. */
  double *ssd
) {
  int i, j;
  double *nx, *ny;

  /* Allocate memory for new data pointers  */
  nx=(double*)calloc(nr, sizeof(double)); if(nx==NULL) return 1;
  ny=(double*)calloc(nr, sizeof(double));
  if(ny==NULL) {free((char*)nx); return 1;}

  /* Copy data to pointers  */
  for(i=0, j=0; i<nr; i++)
    if(!isnan(x[i]) && !isnan(y[i])) {nx[j]=x[i]; ny[j]=y[i]; j++;}

  /* Use llsqperp() */
  i=llsqperp(nx, ny, j, slope, ic, ssd);
  free((char*)nx); free((char*)ny);
  return(i);
}
/******************************************************************************/

/******************************************************************************/
/** Finds the real roots of a*x^2 + b*x + c = 0
    @return Returns the nr of roots, and the roots in m1 and m2.
 */
int quadratic(
  /** Input A */
  double a, 
  /** Input B */
  double b, 
  /** Input C */
  double c, 
  /** Output: Root 1 */
  double *m1, 
  /** Output: Root 2 */
  double *m2
) {
  double discriminant, r, r1, r2, sgnb, temp;

  if(a==0) {if(b==0) return(0); else {*m1=*m2=-c/b; return(1);}}
  discriminant=b*b - 4*a*c;
  if(discriminant>0) {
    if(b==0) {r=fabs(0.5*sqrt(discriminant)/a); *m1=-r; *m2=r;}
    else {
      sgnb=(b>0 ? 1:-1); temp=-0.5*(b + sgnb*sqrt(discriminant));
      r1=temp/a; r2=c/temp; if(r1<r2) {*m1=r1; *m2=r2;} else {*m1=r2; *m2=r1;}
    }
    return(2);
  } else if(discriminant==0) {
    *m1=-0.5*b/a; *m2=-0.5*b/a;
    return(2);
  } else {
    return(0);
  }
}
/******************************************************************************/

/******************************************************************************/
/** Median-based distribution-free estimation of slope and intercept.
    This method has no need for weighting and is insensitive to outliers.
    Note that this is not LMS !
    
    Reference (containing reference to the original idea):
    1. Siegel AF. Robust regression using repeated medians. Biometrika 1982;69(1):242-244.
  
    @sa llsqperp
    @return If successful, function returns value 0.
 */
int medianline(
  /** Coordinates of data points (dimension nr). */
  double *x,
  /** Coordinates of data points (dimension nr). */
  double *y,
  /** Number of data points. */
  int nr,
  /** Estimated slope. */
  double *slope,
  /** Estimated intercept. */
  double *ic
) {
  int i, j, snr;
  double *sp, *ip, d;

  if(LLSQWT_TEST) fprintf(stdout, "medianline()\n");
  /* Check the data */
  if(nr<2 || x==NULL || y==NULL) return(1);
  /* Allocate memory for slopes and intercepts */
  for(i=0, snr=0; i<nr-1; i++) for(j=i+1; j<nr; j++) snr++;
  sp=(double*)malloc(snr*sizeof(double));
  ip=(double*)malloc(snr*sizeof(double));
  if(sp==NULL || ip==NULL) return(2);
  /* Calculate the slopes and intercepts */
  for(i=0, snr=0; i<nr-1; i++) for(j=i+1; j<nr; j++) {
    if(isnan(x[i]) || isnan(x[j]) || isnan(y[i]) || isnan(y[j])) continue;
    d=x[j]-x[i]; if(d==0) continue;
    sp[snr]=(y[j]-y[i])/d; ip[snr]=y[i]-sp[snr]*x[i];
    snr++;
  }
  if(snr<2) {free(sp); free(ip); return(3);}
  /* Get the medians of slope and intercept */
  qsort(sp, snr, sizeof(double), _medianline_cmp);
  qsort(ip, snr, sizeof(double), _medianline_cmp);
  if(snr%2==1) {*slope=sp[snr/2]; *ic=ip[snr/2];}
  else {*slope=0.5*(sp[snr/2]+sp[(snr/2)-1]); *ic=0.5*(ip[snr/2]+ip[(snr/2)-1]);}
  free(sp); free(ip);
  return(0);
}
/******************************************************************************/

/******************************************************************************/
/// @cond
int _medianline_cmp(const void *e1, const void *e2)
{
  if(*((double*)e1) > *((double*)e2)) return(-1);
  else if(*((double*)e1) < *((double*)e2)) return(1);
  else return(0);
}
/// @endcond
/******************************************************************************/

/******************************************************************************/
