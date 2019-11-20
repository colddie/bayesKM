/// @file aic.c
/// @brief Routines for model selection and weighting using Akaike's information criteria.
///
/// For reference, see Sederholm K. TPCMOD0016 in
/// http://pet.utu.fi/staff/vesoik/reports/tpcmod0000.html
///
/// @author Kaisa Sederholm, Vesa Oikonen
///
/****************************************************************************/
#include "libtpcmodel.h"
/****************************************************************************/

/****************************************************************************/
/** Computation of AICc in the special case of sum-of-squares optimization
    from the SS, nr of fitted points and nr of fitted parameters.

    If variance is different between the data points, weighted SS must be given.
    @sa parFreeNr, dftWSampleNr, aicWeights, aicWeightedAvg
    @return Returns the AIC value.
*/
double aicSS(
  /** Sum-of-Squares of the fit */
  double ss,
  /** Sample size, i.e. nr of fitted data points */
  const int n,
  /** Number of fitted model parameters; AICc calculation is valid only when (n-k)>1. */
  const int k
) {
  if(!(ss>=0.0) || n<1 || k<0 || (n-k)<2) return(nan(""));

  double aic=0.0, bias_adj, css;
  int dr, dv;
  
  dr=n-k-1; dv=2*k*(k+1);
  if(dr>0) bias_adj=((double)dv)/((double)dr); else bias_adj=0.0;
  if(ss<1.0e-50) css=1.0e-50; else css=ss; /* Because log(0) is an error */
  if(n>0) aic= n*log(css/(double)n) + 2.0*(double)k + bias_adj;
  return(aic);
}
/****************************************************************************/

/****************************************************************************/
/** @brief Calculate the number of free parameters.

    Model parameters can be fixed by setting lower and upper limit to equal values.
    This function simply checks the limits for each parameter.
    
    @return Returns the number of free parameters.
    @sa aicSS, dftWSampleNr
 */
int parFreeNr(
  /** Nr of parameters */
  const int n,
  /** Lower limits (array of length n) */
  double *pLower,
  /** Upper limits (array of length n) */
  double *pUpper
) {
  if(n<1 || pLower==NULL || pUpper==NULL) return(0);
  int nf=0;
  for(int i=0; i<n; i++) {
    double range=pUpper[i]-pLower[i];
    if(range>1.0E-10) nf++;
  }
  return(nf);
}
/*****************************************************************************/

/*****************************************************************************/
/** Computation of the Akaike weights for model averaging.
    Requires an array of AIC values, and an output array for weights.
   @sa aicSS, aicWeightedAvg, aicModel
   @return Returns 0, if OK.
*/
int aicWeights(
  /** Array of AICs */
  double *aic,
  /** Array of weights (output) */
  double *w,
  /** Lengths of arrays */
  int n
) {
  int i, mini;
  double minaic, sume;

  /* Check data */
  if(n<1 || aic==NULL || w==NULL) return(1);
  if(n==1) {w[0]=1.0; return(0);}
  /* Find out which model gave the smallest AIC */
  mini=0; for(i=1; i<n; i++) if(aic[i]<aic[mini]) mini=i;
  minaic=aic[mini];
  /* Compute relative weights for each model */
  for(i=0, sume=0.0; i<n; i++) {
    w[i]=exp(-0.5*(aic[i]-minaic));
    sume+=w[i];
  }
  if(sume==0.0) return(2);
  for(i=0; i<n; i++) w[i]/=sume;
  return(0);
}
/****************************************************************************/

/****************************************************************************/
/** Computation of the Akaike weighted model parameter average.
    Requires arrays of AIC weight values, and corresponding parameter values.
   @sa aicSS, aicModel, aicWeights
   @return Returns the weighted average.
*/
double aicWeightedAvg(
  /** Array of weights */
  double *w,
  /** Array of parameters */
  double *p,
  /** Lengths of arrays */
  int n
) {
  int i;
  double avg;

  /* Check data */
  if(n<1 || w==NULL || p==NULL) return(1);
  for(i=0, avg=0.0; i<n; i++) avg+=w[i]*p[i];
  return(avg);
}
/****************************************************************************/

/****************************************************************************/
/** Calculates a value describing the relative goodness of models, based on
    an array of model weights.
   @sa aicSS, aicWeightedAvg, aicWeights
   @return Returns the weighted average of model number.
*/
double aicModel(
  /** Array of weights */
  double *w,
  /** Length of array */
  int n
) {
  int i;
  double avg;

  if(n<1 || w==NULL) avg=0.0;
  else for(i=0, avg=0.0; i<n; i++) avg+=w[i]*(double)(i+1);
  return(avg);
}
/****************************************************************************/

/****************************************************************************/

