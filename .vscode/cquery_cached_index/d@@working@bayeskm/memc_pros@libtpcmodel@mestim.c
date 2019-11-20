/// @file mestim.c
/// @author Kaisa Sederholm, Vesa Oikonen
/// @brief Calculating Hubers M-estimator for single data.
///
/*****************************************************************************/
#include "libtpcmodel.h"
/*****************************************************************************/

/*****************************************************************************/
/** Fit a constant (horisontal straight line) to the data with M-estimator.

    The algorithm is described in the lecture notes of Nonlinear signal 
    processing course of Tampere university of technology
    http://www.cs.tut.fi/~eeroh/nonlin.html

\return Returns Hubers M-estimator for single dataset.
 */
double mEstim(
  /** Data array */
  double *data,
  /** Number of data values */
  int nr,
  /** Number of iterations */
  int iterNr,
  /** cutoff point */
  double cutoff
) {

  int ii, in;
  double theta, sum1, sum2, help;

  theta=dmedian(data, nr);
  for(ii=0; ii<iterNr; ii++){
    sum1=sum2=0;
    for(in=0; in<nr; in++){
      if(data[in]<0.9999*theta || data[in]>1.0001*theta){
        help=huber(data[in]-theta, cutoff);
        sum1=sum1+data[in]*help/(data[in]-theta);
        sum2=sum2+help/(data[in]-theta);
      } else {
	sum1=sum1+cutoff*data[in];
        sum2=sum2+cutoff; 
      }
    }
    theta=sum1/sum2;
  }
  return theta;
}
/*****************************************************************************/

/*****************************************************************************/
/** Hubers function.
\return Returns x if |x|<b, and b otherwise.
*/
double huber(
  /** parameter x*/
  double x, 
  /** cutoff point*/
  double b
) {
  double help;

  if(x<-b) help=-b; else help=x;
  if(help<b) {return help;} else {return b;}
}
/*****************************************************************************/

/*****************************************************************************/
