/// @file mtga.c
/// @author Vesa Oikonen
/// @brief Multiple-time graphical analysis (including Patlak and Logan plots).
///
/*****************************************************************************/
#include "libtpcmodel.h"
/*****************************************************************************/

/*****************************************************************************/
/** Calculates Gjedde-Patlak plot x,y values from the measured input and ROI concentration TACs. 

    Plot will not include data where: 
    1) any of the values is not available (NaN), 
    2) integral is negative at this or later point, 
    3) divider is too close to zero, 4) plot x value is negative.

   @sa logan_data, img_patlak, llsqperp, mtga_best_perp
   @return Returns the number of acceptable Gjedde-Patlak plot data pairs,
    or <0 in case of an error. Note that zero or one return values prevent the 
    line fitting, but are here not considered as errors.
 */
int patlak_data(
  /** Nr of samples. */
  int data_nr,
  /** Array of input concentrations. */
  double *i,
  /** Array of integrals (from zero to sample time) of input concentrations;
      if reference region input, then remember to consider the frame length. */
  double *ii,
  /** Array of ROI concentrations. */
  double *c,
  /** Pointer to preallocated memory (at least size dnr) where MTGA plot 
      x values will be written. */
  double *x,
  /** Pointer to preallocated memory (at least size dnr) where MTGA plot 
      y values will be written. */
  double *y
) {
  int fi, plot_nr=0;
  double divider_limit=1.0E-12;

  if(data_nr<0 || i==NULL || ii==NULL || c==NULL || x==NULL || y==NULL)
    return -1;
  for(fi=0; fi<data_nr; fi++) {
    // check that all measured data is available
    if(isnan(i[fi]) || isnan(ii[fi]) || isnan(c[fi])) continue;
    if(!(i[fi]>-1.0E+20 && i[fi]<+1.0E+20)) continue; 
    if(!(ii[fi]>-1.0E+20 && ii[fi]<+1.0E+20)) continue; 
    if(!(c[fi]>-1.0E+20 && c[fi]<+1.0E+20)) continue; 
    // check that integral has been >=0 all the time
    if(ii[fi]<0.0) {plot_nr=0; continue;}
    // check that dividers are not too close to zero
    if(fabs(i[fi])<divider_limit) continue;
    // calculate plot x axis value
    x[plot_nr]=ii[fi]/i[fi];
    if(!(x[plot_nr]>-1.0E+20 && x[plot_nr]<+1.0E+20)) continue;
    if(x[plot_nr]<0.0) continue;
    // calculate plot y axis value
    y[plot_nr]=c[fi]/i[fi];
    if(!(y[plot_nr]>-1.0E+20 && y[plot_nr]<+1.0E+20)) continue;
    // so this plot data point is fine
    plot_nr++;
  }
  return(plot_nr);
}
/*****************************************************************************/

/*****************************************************************************/
/** Calculates Logan plot x,y values from the measured input and ROI concentration TACs. 

    Plot will not include data where: 
    1) any of the values is not available (NaN), 
    2) integral is negative at this or later point, 
    3) divider is too close to zero.

    @sa patlak_data, img_logan, llsqperp, mtga_best_perp
    @return Returns the number of acceptable Logan plot data pairs, or <0 in case of an error. 
    Note that zero or one return values prevent the line fitting, but are here not considered 
    as errors.
 */
int logan_data(
  /** Nr of samples. */
  int data_nr,
  /** Array of input concentrations. */
  double *i,
  /** Array of integrals (from zero to sample time) of input concentrations;
      if reference region input, then remember to consider the frame length. */
  double *ii,
  /** Array of ROI concentrations. */
  double *c,
  /** Array of ROI integrals (from zero to frame middle time). */
  double *ci,
  /** Reference region k2; set to <=0 if not needed. */
  double k2,
  /** Pointer to preallocated memory (at least size dnr) where MTGA plot x values will be written. */
  double *x,
  /** Pointer to preallocated memory (at least size dnr) where MTGA plot y values will be written. */
  double *y
) {
  int fi, plot_nr=0;
  double divider_limit=1.0E-18;

  if(data_nr<0 || i==NULL || ii==NULL || c==NULL || ci==NULL) return -1;
  if(x==NULL || y==NULL) return -1;

  for(fi=0; fi<data_nr; fi++) {
    // check that all measured data is available
    if(isnan(i[fi]) || isnan(ii[fi]) || isnan(c[fi]) || isnan(ci[fi])) continue;
    if(!(i[fi]>-1.0E+30 && i[fi]<+1.0E+30)) continue; 
    if(!(ii[fi]>-1.0E+30 && ii[fi]<+1.0E+30)) continue; 
    if(!(c[fi]>-1.0E+30 && c[fi]<+1.0E+30)) continue; 
    if(!(ci[fi]>-1.0E+30 && ci[fi]<+1.0E+30)) continue; 
    // check that integrals have been >=0 all the time
    if(ii[fi]<0.0) {plot_nr=0; continue;}
    if(ci[fi]<0.0) {plot_nr=0; continue;}
    // check that dividers are not too close to zero
    if(fabs(c[fi])<divider_limit) continue;
    // calculate plot x axis value
    if(k2>0.0) x[plot_nr]=(ii[fi]+i[fi]/k2)/c[fi];
    else x[plot_nr]=ii[fi]/c[fi];
    if(!(x[plot_nr]>-1.0E+30 && x[plot_nr]<+1.0E+30)) continue;
    // calculate plot y axis value
    y[plot_nr]=ci[fi]/c[fi];
    if(!(y[plot_nr]>-1.0E+30 && y[plot_nr]<+1.0E+30)) continue;
    // so this plot data point is fine
    plot_nr++;
  }
  return(plot_nr);
}
/*****************************************************************************/

/*****************************************************************************/
/** Finds the best regression line to (x,y)-data, leaving points out from the beginning, 
    because Gjedde-Patlak and Logan plots reach linearity at some later phase. 

    This function applies llsqperp() which is a non-iterative perpendicular line fitting routine.

    @sa llsqperp(), logan_data, patlak_data
    @return Returns 0, if successful, and <>0 in case of an error
 */
int mtga_best_perp(
  /** Plot x axis values. */
  double *x,
  /** Plot y axis values. */
  double *y,
  /** Nr of plot data points. */
  int nr,
  /** Slope is returned in here. */
  double *slope,
  /** Y axis intercept is returned in here. */
  double *ic,
  /** Sum of squared distances / fnr, or NULL if not needed. */
  double *ssd,
  /** Number of points in the best fit, or NULL if not needed. */
  int *fnr
) {
  int from, to, ret, from_min, to_min;
  double lic, lslope, lssd, ssd_min;

  /* Search the plot range that gives the lowest ssd */
  ssd_min=9.99E+99; from_min=to_min=-1;
  for(from=0, to=nr-1; ((to-from)+1)>=MTGA_BEST_MIN_NR; from++) {
    ret=llsqperp(x+from, y+from, (to-from)+1, &lslope, &lic, &lssd);
    if(ret==0 && lssd<ssd_min) {
      ssd_min=lssd; from_min=from; to_min=to;
      *slope=lslope; *ic=lic; if(ssd!=NULL) *ssd=lssd;
    }
  }
  if(from_min<0) {
    if(fnr!=NULL) *fnr=0;
    if(ssd!=NULL) *ssd=0.0;
    return(5);
  }
  if(fnr!=NULL) *fnr=(to_min-from_min)+1;
  return(0);
}
/*****************************************************************************/

/*****************************************************************************/
