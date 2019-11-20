/// @file hholder.c
/// @author Kaisa Sederholm, Vesa Oikonen
/// @brief Implementation and use of Householder transform.
///
/// These routines are based on the code 
/// provided in the GSL library (http://sources.redhat.com/gsl/).
///
/*****************************************************************************/
#include "libtpcmodel.h"
/*****************************************************************************/

/*****************************************************************************/
/** This function prepares a Householder transformation
  P = I - tau h h^T 
  which can be used to zero all the elements of the input vector 
  except the first one that will get value beta. 
  On output the elements 1 - size-1 of the vector h are stored 
  in locations vector[1] - vector[size-1] of the input vector
  and value of beta is stored in location vector[0].
  
  @return The scalar tau is returned.
*/
double householder_transform(
  /** The N-vector to be transformed. */
  double *v,
  /** size of the vector. */
  int N
) {
  double vnorm, alpha, beta, tau;

  if(N<1) return 0.0; // tau = 0
  /* Euclidean norm of the vector starting from the second value */
  vnorm=0.0; for(int n=1; n<N; n++) vnorm+=v[n]*v[n];
  vnorm=sqrt(vnorm); if(isnan(vnorm) || vnorm==0.0) return 0.0; // tau = 0

  /* Computing the coefficient tau */
  alpha=v[0];
  beta= - (alpha >= 0.0 ? +1.0 : -1.0) * hypot(alpha, vnorm);
  tau=(beta-alpha)/beta ;

  /* Scale the Householder vector so that the first element will be 1.
   * (Scaling is also affecting the coefficient tau).
   * Without scaling, the first element would have value (alpha - beta). */
  {
    double s=alpha-beta;
    if(fabs(s)>DBL_MIN) {
      v[0]=beta;
      for(int n=1; n<N; n++) v[n]*=(1.0/s);
    } else {
      v[0]=beta;
      for(int n=1; n<N; n++) v[n]*=(doubleMachEps()/s);
      for(int n=1; n<N; n++) v[n]*=(1.0/doubleMachEps());
    }
  }

  return tau;
}
/*****************************************************************************/

/*****************************************************************************/
/** Applies a householder transformation defined by vector "vector" and
 *  scalar tau to the left-hand side of the matrix. 
 *  (I - tau vector vector^T)*matrix
 *  The result of the transform is stored in matrix.
 *
 *  @return Returns 0 if ok.
 */
int householder_hm(
  /** Coefficient defining householder transform. */
  double tau,
  /** Vector defining householder transform (of size rowNr). */
  double *vector,
  /** the matrix that is to be transformed. */
  double **matrix,
  /** Nr of rows in matrix. */
  int rowNr,
  /** Nr of columns in matrix. */
  int columnNr
) {
  int i, j;
  double wj;

  if(tau==0.0) return(0); // success
  if(rowNr<1 || columnNr<1) return(1);
  for(j=0; j<columnNr; j++) {
    /* Compute wj = vk Akj */
    wj=matrix[0][j];
    for(i=1; i<rowNr; i++)  /* note, computed for v(0) = 1 above */
      wj += vector[i]*matrix[i][j];
    /* Aij = Aij - tau vi wj */
    /* i = 0 */
    matrix[0][j]-=tau*wj;
    /* i = 1 .. M-1 */
    for(i=1; i<rowNr; i++) matrix[i][j]-=tau*vector[i]*wj;
  }
  return(0);
}
/*****************************************************************************/

/*****************************************************************************/
/** Applies a householder transformation defined by vector v and
 *  coefficient tau to vector w 
 *  w = (I - tau v v^T) w.
 * 
 *  @return Returns 0 if ok.
 */
int householder_hv(
  /** Coefficient defining householder transform. */
  double tau,
  /** Size of vectors v and w. */
  int size,
  /** Vector v. */
  double *v,
  /** Vector w. */
  double *w
) {
  int i;
  double d;

  if(tau==0) return(0); // success
  if(size<1) return(1);
  /* d = v'w */
  d=w[0]; for(i=1; i<size; i++) d+=v[i]*w[i];
  /* w = w - tau (v) (v'w) */
  w[0]-=tau*d;
  for(i=1; i<size; i++) w[i]-=tau*v[i]*d;
  return(0);
}
/*****************************************************************************/

/*****************************************************************************/
/** Calculates the euclidean norm of vector v[].
 *
\return Returns the euclidean norm of a vector.
 */
double householder_norm(
  /** Vector v */
  double *v,
  /** Size of vector v[] */
  int size
) {
  double help;
  int i;

  for(i=0, help=0; i<size; i++) help+=v[i]*v[i];
  return sqrt(help);
}
/*****************************************************************************/

/*****************************************************************************/
