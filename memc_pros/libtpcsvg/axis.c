/// @file axis.c
/// @author Vesa Oikonen
/// @brief Axis for XY plots.
///
/*****************************************************************************/
#include "libtpcsvg.h"
/*****************************************************************************/

/*****************************************************************************/
/** Define suitable tick positions for XY plot.
   @return Returns 0 if successful, <>0 in case of error.
   @sa axis_check_range, svg_calculate_axes
*/
int axis_tick_positions(
  /** axis range minimum. */
  const double begin,
  /** axis range maximum. */
  const double end,
  /** Pointer to array where tick values will be written. Length must be at least tick_nr. */
  double *ticks,
  /** Input: Max allowed nr of ticks; Output: actual nr of ticks. */
  int *tick_nr,
  /** Output: Suggested scale factor (x10^sf). NULL, if not needed. */
  double *scale_factor,
  /** Output: tick value precision (nr of decimals). NULL, if not needed. */
  int *tick_decimals,
  /** Verbose level; if zero, then nothing is printed to stderr or stdout. */
  int verbose
) {
  int ti;
  double step, scale;

  if(verbose>0)
    printf("axis_tick_positions(%.10E, %.10E, ticks, %d, sf, td, %d)\n", 
      begin, end, *tick_nr, verbose);

  /* Check input */
  if(ticks==NULL || *tick_nr<1) return(1);
  //if(verbose>1) printf("  begin:=%.10E\n  end:=%.10E\n", begin, end);
  if(end<=begin) {
    *tick_nr=0;
    if(scale_factor!=NULL) *scale_factor=1.0;
    if(tick_decimals!=NULL) *tick_decimals=0;
    return(0);
  }
  /* Calculate the initial tick step size */
  step=(end-begin)/(double)(*tick_nr);
  if(verbose>1) printf("  tick_nr:=%d\n  step:=%20.10E\n", *tick_nr, step);

  /* Calculate a feasible step size and scale */
  scale=1.0;
  while(step<=0.5) {step*=10.0; scale/=10.0;}
  while(step>5.0) {step/=10.0; scale*=10.0;}
  if(verbose>1) printf("  scaled_step:=%.10E\n  scale:=%E\n", step, scale);
  if     (step<1.0) {step=1.0;}
  else if(step<2.0) {step=2.0;}
  else              {step=5.0;}
  if(verbose>1) {
    printf("  feasible step:=%g\n", step);
  }
  ti=0; ticks[ti]=step*scale*ceil(begin/(step*scale));
  while(ticks[ti]<=(end+(end-begin)*0.00001) && ti<*tick_nr) {
    ti++; ticks[ti]=ticks[0]+step*scale*(double)ti;
  }
  *tick_nr=ti;

  if(verbose>1) {
    printf("   final tick_nr := %d\n", *tick_nr);
    printf("    ticks: %.10E", ticks[0]);
    for(ti=1; ti<*tick_nr; ti++) printf(", %.10E", ticks[ti]);
    printf("\n");
  }

  /* Quit, if user did not like scales etc */
  if(scale_factor==NULL && tick_decimals==NULL) return(0);

  /* Calculate step scale */
  double step_scale=log10(scale);
  if(verbose>1) printf("  step_scale=%g\n", step_scale);

  /* Find the highest tick label */
  double tick_high;
  int tick_scale, scale_dif;
  if(fabs(ticks[0])>fabs(ticks[*tick_nr-1])) tick_high=ticks[0];
  else tick_high=ticks[*tick_nr-1];
  if(verbose>1) printf("  tick_high := %.10E\n", tick_high);
  tick_scale=0;
  while(fabs(tick_high)<1.0) {tick_high*=10.0; tick_scale--;}
  while(fabs(tick_high)>=10.0) {tick_high/=10.0; tick_scale++;}
  if(verbose>1) printf("  scaled_tick_high := %g\n  tick_scale := %d\n", tick_high, tick_scale);

  /* Calculate the difference between tick and step scales */
  scale_dif=tick_scale-step_scale;
  if(verbose>1) printf("  scale_dif := %d\n", scale_dif);
  /* That determines the tick number precision */
  if(tick_decimals!=NULL) {
    *tick_decimals=1+scale_dif;
    if(verbose>1) printf("  tick_decimals := %d\n", *tick_decimals);
  }

  /* Calculate the preferred tick scale factor, or quit, if it is not needed */
  if(scale_factor==NULL) return(0);
  *scale_factor=tick_scale;
  if(verbose>1) printf("  scale_factor := %g\n", *scale_factor);

  return(0);
}
/*****************************************************************************/

/*****************************************************************************/
/** Check and if necessary correct axis range: min must be smaller than max.
    Also, if the range is close to zero, then make it larger.
    Also, if the range is relatively large, then change it to start from zero.
   @sa axis_tick_positions, svg_calculate_axes
 */
void axis_check_range(
  /** axis range minimum */
  double *begin,
  /** axis range maximum */
  double *end,
  /** Verbose level; if zero, then nothing is printed to stderr or stdout */
  int verbose
) {
  double temp;

  if(verbose>0) printf("axis_check_range(%g, %g, %d)\n", *begin, *end, verbose);

  if(*begin<*end) {
    // fine now
  } else if(*begin>*end) {
    temp=*end; *end=*begin; *begin=temp;
  } else { // this means that begin==end
    if(*begin==0.0) {
      *begin=-1.0; *end=1.0;
    } else if(*begin<0.0) {
      *begin*=2.0; *end=0.0;
    } else {
      *begin=0.0; *end*=2.0;
    }
    if(verbose>2) {
      printf("  new begin := %g\n", *begin);
      printf("  new end := %g\n", *end);
    }
    return;
  }
  /* If range is very small */
  if((*end-*begin)<1.0e-98) {
    if(*begin>=1.0e-98 || *begin<0.0) *begin-=1.0e-98;
    else if(*begin>=0.0) *begin=0.0;
    if(*end<=-1.0e-98 || *end>0.0) *end+=1.0e-98;
    else if(*end<0.0) *end=0.0;
    if(verbose>2) {
      printf("   new begin := %g\n", *begin);
      printf("   new end := %g\n", *end);
    }
    return;
  }
  /* If range is relatively large */
  if(*begin>0.0 && *end>0.0) {
    if((*end-*begin)>3.3*(*begin)) *begin=0.0;
  } else if(*begin<0.0 && *end<0.0) {
    if((*end-*begin)>3.3*(-*end)) *end=0.0;
  }
  /* If data range is relatively small (compared to level) */
  temp=(*end-*begin)*2.0/(fabs(*end)+fabs(*begin));
  if(temp<0.01) {
    temp=0.5*(*end-*begin)*0.01/temp;
    if(*begin<0.0 || *begin>temp) *begin-=temp; else *begin=0.0;
    if(*end>0.0 || *end<-temp) *end+=temp; else *end=0.0;
  }

  return;
}
/*****************************************************************************/

/*****************************************************************************/
/** Remove initial exponential zeroes from the exponential part of string
    representation of values, for example, '1.01E-010' -> '1.01E-10'. */
void strRmExpZeroes(char *str)
{
  char *p, *r;
  int i, n;
  
  if(str==NULL || strlen(str)<3) return;
  // search the exp E
  p=strrchr(str, 'E'); if(p==NULL) p=strrchr(str, 'e'); if(p==NULL) return;
  // if only +,- and 0 after it, then remove also E
  r=p+1; if(strspn(r, "0+-")==strlen(r)) {*p=(char)0; return;}
  // remove next '+' or jump over '-'
  p++;
  if(*p=='+') {r=p; while(1) {*r=*(r+1); if(*r==(char)0) break; r++;}}
  else if(*p=='-') p++;
  // count initial zeroes
  n=strspn(p, "0"); if(n==0) return;
  // copy the rest of the string over zeroes
  r=p+n; n=strlen(r);
  for(i=0; i<=n; i++) {*p=*r; p++; r++;}
  return;
}
/*****************************************************************************/

/*****************************************************************************/
