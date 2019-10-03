/// @file qr.c
/// @author Kaisa Sederholm, Vesa Oikonen
/// @brief Routines needed in the use of QR 
///        decomposition when solving least squares problems.
///
/// These routines are based on the code of Gerard Jungman and 
/// Brian Gough provided in the GSL library (http://sources.redhat.com/gsl/).
///
/// @todo Try the other implementation of qr (in v2), and add svd as an alternative. 
///
/*****************************************************************************/
#include "libtpcmodel.h"
/*****************************************************************************/

/*****************************************************************************/
/** Algorithm QR.
 * 
 * Solves a matrix form least square problem 
 * min||A x - b|| => A x =~ b     (A is m*n matrix, m>=n)
 * using the QR decomposition for overdetermined systems.
 * Based on GNU Scientific Library, edited by Kaisa Liukko.
 *
 * Instead of pointers for working space, NULL can be given to let this
 * function to allocate and free the required memory.
 *
 * @return Function returns 0 if successful and 1 in case of invalid problem
 *         dimensions or memory allocation error.
 */
int qr(
  /** On entry, a[m][n] contains the m by n matrix A.
      On exit, a[][] contains the QR factorization. */
  double **A,
  /** Dimensions of matrix A are a[m][n]. */
  int m,
  /** Dimensions of matrix A are a[m][n]. */
  int n,
  /** B[] is an m-size vector containing the right-hand side vector b. */
  double *B,
  /** On exit, x[] will contain the solution vector x (size of n). */
  double *X,
  /** On exit, rnorm (pointer to double) contains the squared Euclidean norm of
   *  the residual vector (R^2); enter NULL if not needed. */
  double *rnorm,
  /** On exit, tau[] will contain the householder coefficients (size of n);
   *  enter NULL, if not needed. */
  double *tau,
  /** An m-size array of working space, res[]. On output contains
      residual b - Ax. Enter NULL to let qr() to handle it. */
  double *res,
  /** m*n array of working space. Enter NULL to let qr() to handle it. */ 
  double **wws,
  /** 2m-array of working space. Enter NULL to let qr() to handle it. */
  double *ws
) {
  int i;
  double *qrRes, *qrTau, **qrWws, *qrWs, *chain;

  /* Check the parameters and data */
  if(m<=0 || n<=0 || A==NULL || B==NULL || X==NULL) return(1);
  if(m<n) return(1);

  /* Allocate memory for working space, if required */
  if(tau!=NULL) qrTau=tau; else qrTau=(double*)calloc(n, sizeof(double));
  if(res!=NULL) qrRes=res; else qrRes=(double*)calloc(m, sizeof(double));
  if(wws!=NULL) {
    qrWws=wws; chain=(double*)NULL;
  } else {
    qrWws=(double**)malloc(m * sizeof(double*));
    chain=(double*)malloc(m*n * sizeof(double));
    for(i=0; i<m; i++) qrWws[i]=chain + i*n;
  }
  if(ws!=NULL) qrWs=ws; else qrWs=(double*)calloc(2*m, sizeof(double));
  if(qrTau==NULL || qrRes==NULL || qrWws==NULL || qrWs==NULL) return(1);

  /* Form the householder decomposition and solve the least square problem */
  if(qr_decomp(A, m, n, qrTau, qrWws, qrWs)) return(2);
  if(qr_solve(A, m, n, qrTau, B, X, qrRes, rnorm, qrWws, qrWs)) return(3);

  /* Free working space, if it was allocated here */
  if(tau==NULL) free(qrTau);
  if(res==NULL) free(qrRes);
  if(wws==NULL) {free(qrWws); free(chain);}
  if(ws==NULL) free(qrWs);
  for(i=0; i<n; i++) if(isnan(X[i])) return(4);
  return(0);
} /* qr */
/*****************************************************************************/

/*****************************************************************************/
/** Factorise a general M x N matrix A into A = Q R ,
    where Q is orthogonal (M x M) and R is upper triangular (M x N).

    Q is stored as a packed set of Householder vectors in the strict lower triangular 
    part of the input matrix A and a set of coefficients in vector tau.

    R is stored in the diagonal and upper triangle of the input matrix.

    The full matrix for Q can be obtained as the product Q = Q_1 Q_2 .. Q_k 
    and it's transform as the product Q^T = Q_k .. Q_2 Q_1 , where k = min(M,N) and
    Q_i = (I - tau_i * h_i * h_i^T)
    and where h_i is a Householder vector h_i = [1, A(i+1,i), A(i+2,i), ... , A(M,i)].
    This storage scheme is the same as in LAPACK.  
  
    NOTICE! The calling program must take care that pointer tau is 
    of size N or M, whichever is smaller. 

   @return Function returns 0 if ok.
 */
int qr_decomp(
  /** contains coefficient matrix A (m*n) as input and factorisation QR as output. */
  double **a,
  /** nr of rows in matrix A. */
  int M,
  /** nr of columns in matrix A. */
  int N,
  /** Vector for householder coefficients, of length N or M, whichever is smaller. */
  double *tau,
  /** m*n matrix of working space. */
  double **cchain,
  /** m size array of working space. */
  double *chain
) {
  //printf("qr_decomp()\n");
  int i, m, n, MNmin;
  double *subvector, **submatrix;

  /* Local variables */
  if(M<N) MNmin=M; else MNmin=N;
  if(MNmin<1 || a==NULL || tau==NULL || cchain==NULL || chain==NULL) return(1);

  subvector=chain;
  submatrix=cchain;

  for(i=0; i<MNmin; i++) {
    //printf("i=%d (MNmin=%d)\n", i, MNmin);
    /* Compute the Householder transformation to reduce the j-th column of the matrix A to 
       a multiple of the j-th unit vector. Householder vector h_i is saved in the lower triangular 
       part of the column and Householder coefficient tau_i in the vector tau. */
    for(m=i; m<M; m++) subvector[m-i]=a[m][i];
    tau[i] = householder_transform(subvector, M-i);
    for(m=i; m<M; m++) a[m][i]=subvector[m-i];

    /* Apply the transformation to the remaining columns to get upper triangular part of matrix R  */
    if(i+1 < N) {
      for(m=i; m<M; m++)
        for(n=i+1; n<N; n++)
          submatrix[m-i][n-i-1]=a[m][n];
      if(householder_hm(tau[i], subvector, submatrix, M-i, N-i)) return(2);
      for(m=i; m<M; m++)
        for(n=i+1; n<N; n++)
          a[m][n]=submatrix[m-i][n-i-1];
    }
  }

  return(0);
}
/*****************************************************************************/

/*****************************************************************************/
/** Find the least squares solution to the overdetermined system
 
     A x = b 
    
   for m >= n using the QR factorisation A = Q R. 
   qr_decomp() must be used prior to this function in order to form
   the QR factorisation of A.
   Solution is formed in the following order: 
   QR x = b  =>  R x = Q^T b  =>  x = R^-1 (Q^T b)
   
   NOTICE! The calling program must take care that pointers b, x
   and residual are of the right size.
   @return Function returns 0 if ok.
 */
int qr_solve(
  /** m*n matrix containing householder vectors of A. */
  double **QR,
  /** nr of rows in matrix A. */
  int M,
  /** nr of columns in matrix A. */
  int N,
  /** vector containing householder coefficients tau; length is M or N, whichever is smaller. */
  double *tau,
  /** Contains m-size vector b of A x = b. */
  double *b,
  /** solution vector x of length n. */
  double *x,
  /** residual vector of length m. */
  double *residual,
  /** norm^2 of the residual vector; enter NULL if not needed. */
  double *resNorm,
  /** m*n matrix of the working space. */
  double **cchain,
  /** 2m length array of the working space. */ 
  double *chain
) {
  //printf("qr_solve()\n");
  if(QR==NULL || tau==NULL || b==NULL || x==NULL || residual==NULL) return(1);
  if(cchain==NULL || chain==NULL) return(1);
  if(M<1 || N<1) return(2);

  int MNmin; if(M<N) MNmin=M; else MNmin=N;

  int m, n;
  double **Rmatrix=cchain;
  double *h=chain;
  double *w=chain+M;

  /* Get matrix R from the upper triangular part of QR
     First the rows N - M-1 are eliminated from R matrix*/
  for(m=0; m<N; m++) for(n=0; n<N; n++) Rmatrix[m][n]=QR[m][n];
  for(m=0; m<M; m++) residual[m]=b[m];

  /* Compute b = Q^T b */
  /* Form the product Q^T residual from householder vectors saved in the lower triangle of 
     QR matrix and householder coefficients saved in vector tau. */
  for(int i=0; i<MNmin; i++) {
    for(m=i; m<M; m++) h[m-i]=QR[m][i]; 
    for(m=i; m<M; m++) w[m-i]=residual[m];
    if(householder_hv(tau[i], M-i, h, w)) return(2);
    for(m=i; m<M; m++) residual[m]=w[m-i];
  }

  /* Solve R x = b by computing x = R^-1 b */
  for(n=0; n<N; n++) x[n]=residual[n];
  /* back-substitution */
  x[N-1]=x[N-1]/Rmatrix[N-1][N-1];
  for(int i=N-2; i>=0; i--) {
    for(int j=i+1; j<N; j++) x[i]-=Rmatrix[i][j]*x[j];
    x[i]/=Rmatrix[i][i];
  }

  /* Compute residual = b - A x = Q (Q^T b - R x) */
  for(n=0; n<N; n++) residual[n]=0.0;
  /* Compute residual= Q*residual */
  /* Form the product Q*residual from householder vectors saved in the lower triangle
     of QR matrix and householder coefficients saved in vector tau. */
  for(int i=MNmin-1; i>=0; i--) {
    for(int m=i; m<M; m++) h[m-i]=QR[m][i]; 
    for(int m=i; m<M; m++) w[m-i]=residual[m];
    if(householder_hv(tau[i], M-i, h, w)) return(2);
    for(int m=i; m<M; m++) residual[m]=w[m-i];
  }

  /* Compute norm^2 of the residual vector, if needed */
  if(resNorm!=NULL)
    for(m=0, *resNorm=0.0; m<M; m++) *resNorm +=residual[m]*residual[m];

  return(0);
}
/*****************************************************************************/

/*****************************************************************************/
/** Algorithm for weighting the problem that is given to QR algorithm.

    Square roots of weights are used because in QR the difference w*A-w*b is squared.
    @return Algorithm returns zero if successful, otherwise <>0.
*/
int qr_weight(
  /** Dimensions of matrix A are a[m][n]. */
  int N,
  /** Dimensions of matrix A are a[m][n]; size of vector B is m. */
  int M,
  /** Matrix a[m][n] for QR, contents will be weighted here. */
  double **A,
  /** B[] is an m-size vector for QR, contents will be weighted here. */
  double *b,
  /** Pointer to array of size m, which contains sample weights, used here
   *  to weight matrix A and vector b. */
  double *weight,
  /** m-sized vector for working space; enter NULL to allocate locally. */
  double *ws
) {
  int n, m;
  double *w;

  /* Check the arguments */
  if(N<1 || M<1 || A==NULL || b==NULL || weight==NULL) return(1);

  /* Allocate memory, if necessary */
  if(ws==NULL) {
    w=(double*)malloc(M*sizeof(double)); if(w==NULL) return(2);
  } else {
    w=ws;
  }

  /* Check that weights are not zero and get the square roots of them to w[]. */
  for(m=0; m<M; m++) {
    if(weight[m]<=1.0e-100) w[m]=1.0e-50;
    else w[m]=sqrt(weight[m]);
  }

  /* Multiply rows of matrix A and elements of vector b with weights. */
  for(m=0; m<M; m++) {
    for(n=0; n<N; n++) A[m][n]*=w[m];
    b[m]*=w[m];
  }

  if(ws==NULL) free(w);
  return(0);
}
/*****************************************************************************/

/*****************************************************************************/
/** @brief Solve over-determined least-squares problem A x ~ b using
    successive Householder rotations.

    This routine is based on the text and Fortran code in
    C.L. Lawson and R.J. Hanson, Solving Least Squares Problems,
    Prentice-Hall, Englewood Cliffs, New Jersey, 1974,
    and Fortran code by R.L. Parker and P.B. Stark.

   @return Returns 0 when successful, 1 if system is singular, and 2 in case of
    other errors, including that system is under-determined.
*/
int qrLH(
  /** Number of samples in matrix A and the length of vector b. */
  const unsigned int m, 
  /** Number of parameters in matrix A and the length of vector x. 
      The n must be smaller or equal to m. */
  const unsigned int n,
  /** Pointer to matrix A; matrix must be given as an n*m array,
      containing n consecutive m-length vectors. 
      Contents of A are modified in this routine. */
  double *a,
  /** Pointer to vector b of length n.
      Contents of b are modified in this routine. */
  double *b,
  /** Pointer to the result vector x of length n. */
  double *x,
  /** Pointer to a double value, in where the sum of squared residuals is written. */
  double *r2
) {
  /* Check the input */
  if(a==NULL || b==NULL || x==NULL || r2==NULL) return(2);
  if(n<1 || m<n) {*r2=nan(""); return(2);}

  /* Initiate output to zeroes, in case of exit because of singularity */
  for(unsigned int ni=0; ni<n; ni++) x[ni]=0.0;
  *r2=0.0;

  /* Rotates matrix A into upper triangular form */
  for(unsigned int ni=0; ni<n; ni++) {
    /* Find constants for rotation and diagonal entry */
    double sq=0.0;
    for(unsigned int mi=ni; mi<m; mi++) sq+=a[mi + ni*m]*a[mi + ni*m];
    if(sq==0.0) return(1);
    double qv1=-copysign(sqrt(sq), a[ni + ni*m]);
    double u1=a[ni + ni*m] - qv1;
    a[ni + ni*m]=qv1;
    unsigned int ni1=ni+1;
    /*  Rotate the remaining columns of sub-matrix. */
    for(unsigned int nj=ni1; nj<n; nj++) {
      double dot=u1*a[ni + nj*m];
      for(unsigned int mi=ni1; mi<m; mi++)
        dot+=a[mi + nj*m] * a[mi + ni*m];
      double c=dot/fabs(qv1*u1);
      for(unsigned int mi=ni1; mi<m; mi++)
        a[mi + nj*m]-=c*a[mi + ni*m];
      a[ni + nj*m]-=c*u1;
    }
    /* Rotate vector B */
    double dot=u1*b[ni];
    for(unsigned int mi=ni1; mi<m; mi++)
      dot+=b[mi]*a[mi + ni*m];
    double c=dot/fabs(qv1*u1);
    b[ni]-=c*u1;
    for(unsigned int mi=ni1; mi<m; mi++)
      b[mi]-=c*a[mi + ni*m];
  } // end of rotation loop

  /* Solve triangular system by back-substitution. */
  for(unsigned int ni=0; ni<n; ni++) {
    int k=n-ni-1;
    double s=b[k];
    for(unsigned int nj=k+1; nj<n; nj++) 
      s-=a[k + nj*m] * x[nj];
    if(a[k + k*m]==0.0) return(1);
    x[k]=s/a[k + k*m];
  }

  /* Calculate the sum of squared residuals. */
  *r2=0.0;
  for(unsigned int mi=n; mi<m; mi++)
    *r2 += b[mi]*b[mi];

  return(0);
}
/*****************************************************************************/

/*****************************************************************************/
