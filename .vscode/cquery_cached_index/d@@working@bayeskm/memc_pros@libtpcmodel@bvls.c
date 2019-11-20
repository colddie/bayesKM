/// @file bvls.c
/// @author Vesa Oikonen
/// @brief BVLS (Bounded-value least-squares).
///
///  Function bvls() is dependent on function qrLH() in qrlsq.c.
///
/*****************************************************************************/
#include "libtpcmodel.h"
/*****************************************************************************/

/*****************************************************************************/
/** @brief Bounded-value least-squares method to solve the linear problem
     A x ~ b , subject to limit1 <= x <= limit2.

    This routine is based on the text and Fortran code in
    C.L. Lawson and R.J. Hanson, Solving Least Squares Problems,
    Prentice-Hall, Englewood Cliffs, New Jersey, 1974,
    and Fortran codes by R.L. Parker and P.B. Stark, and by J. Burkardt.

    @sa nnls, qrLH, llsqWghtSquared, llsqWght
    @return Returns 0 when successful, -1 if iteration count exceeded the limit,
    1 in case of invalid data or settings, and 2 if problem cannot be solved.
*/
int bvls(
  /** Enter 0 to solve the problem from scratch, or <>0 to initialize the routine
      using caller's guess about which parameters are active within their bounds,
      which are at their lower bounds, and which at their upper bounds, as supplied 
      in the array istate ('warm start'). When key <> 0, the routine initially sets 
      the active components to the averages of their upper and lower bounds. */
  int key, 
  /** Number of samples in matrix A and the length of vector b. */
  const /*unsigned*/ int m, 
  /** Number of parameters in matrix A and the length of vector x. 
      The n must not be > m unless at least n-m variables are non-active
      (set to their bounds).. */
  const /*unsigned*/ int n,
  /** Pointer to matrix A; matrix must be given as an n*m array,
      containing n consecutive m-length vectors. 
      Contents of A are modified in this routine. */
  double *a,
  /** Pointer to vector b of length m.
      Contents of b are modified in this routine. */
  double *b,
  /** Array BL[0..n-1] of lower bounds for parameters. */
  double *bl,
  /** Array BU[0..n-1] of upper bounds for parameters. */
  double *bu,
  /** Pointer to the result vector x of length n. */
  double *x,
  /** Pointer to an array of length n, which will be used as working memory,
      and the minimum 2-norm || a.x-b ||, (R^2), will be written in w[0]. */
  double *w,
  /** Pointer to an array of length m*(n+2), or m*(m+2) if m<n, which will be used as 
      working memory. */
  double *act,
  /** Pointer to an array of length m, to be used as working memory. */
  double *zz,
  /** Pointer to an integer array of length n+1. If parameter key <>0, then
      the last position istate[n] must contain the total number of components at 
      their bounds (nbound, the `bound variables'). The absolute values of the first 
      nbound entries of istate[] are the indices of these `bound' components of x[].
      The sign of istate[0.. nbound-1] entries indicates whether x[|istate[i]|] is at its 
      upper (positive) or lower (negative) bound.
      Entries istate[nbound..n-1] contain the indices of the active components. */
  int *istate, 
  /** Number of performed iterations is returned in this variable. 
      Maximum nr of iterations is set with this variable, or set it to 0 to use
      the default maximum, 3*n. */
  int *iter,
  /** Verbose level; if zero, then nothing is printed to stderr or stdout */
  int verbose
) {
  if(verbose>0) {printf("bvls(%d, %d, %d, ...)\n", key, m, n); fflush(stdout);}
  /* Check the input */
  if(a==NULL || b==NULL || bl==NULL || bu==NULL || x==NULL || w==NULL ||
     act==NULL || zz==NULL || istate==NULL || iter==NULL || n<1 || m<1) 
  {
    if(verbose>0) fprintf(stderr, "Error: invalid input to BVLS.\n");
    return(1);
  }

  int maxIter=*iter; if(maxIter<3) maxIter=3*n;
  *iter=0; 
  const double eps=1.0E-13; // stopping rule

  /* Step 1. Initialize everything -- active and bound sets, initial values, etc. */
  if(verbose>1) {printf("step 1\n"); fflush(stdout);}

  /* Set mm to the smaller of matrix dimensions n and m. */
  int mm; if(m<n) mm=m; else mm=n;
  int aindx=0; // one-based index of X[] that determined the alpha; zero means not set.
  int aindxsign=0; // +1 if zz[aindx]> bu[aindx], -1 if zz[aindx]<bl[aindx]
  /* istateFromStep5 is the one-based index of the parameter that most wants to be active;
     zero value means its is currently not set. */
  int istateFromStep5 = 0;

  /* Check the consistency of given bounds bl[] and bu[]. */
  {
    double maxrange=0.0;
    for(int ni=0; ni<n; ni++) {
      double d=bu[ni]-bl[ni];
      if(verbose>3) printf("  bounds[%d]: %g %g\n", 1+ni, bl[ni], bu[ni]);
      if(d<0.0) {
        if(verbose>0) fprintf(stderr, "Error: inconsistent bounds in BVLS.\n"); 
        return(1);
      }
      maxrange=fmax(maxrange, d);
    }
    if(verbose>2) printf("  maxrange := %g\n", maxrange);
    if(maxrange<1.0E-10) {
      if(verbose>0) fprintf(stderr, "Error: no free variables in BVLS.\n");
      return(1);
    }
  }

  /* In a fresh initialization (key=0), bind all variables at their lower bounds. 
     If key<>0, use the supplied istate[] array to initialize the variables. */
  int nbound, nact; // number of bound and active parameters, respectively.
  if(key==0) {
    nbound=n;
    /* Write the indices; negative sign indicates lower limit */
    /* Note that these indices must start from 1, because 0 can not have sign */
    for(int ni=0; ni<nbound; ni++) istate[ni]=-(1+ni);
  } else {
    nbound=istate[n];
  }
  nact=n-nbound;
  if(nact>mm) {
    if(verbose>0) fprintf(stderr, "Error: too many active variables in BVLS starting solution.\n");
    return(2);
  }
  for(int ni=0; ni<nbound; ni++) {
    int i=abs(istate[ni])-1;
    if(istate[ni]<0) x[i]=bl[i]; else x[i]=bu[i];
  }

  /* In a warm start (key<>0, and nbound<n) initialize the active variables to the mean of
     their bounds. This is needed in case the initial QR results in active variables 
     out-of-bounds and Steps 8-11 get executed the first time through. */
  for(int ni=nbound; ni<n; ni++) {
    int i=abs(istate[ni])-1; // indices of active variables should not have signs, but let's be sure
    x[i]=0.5*(bl[i]+bu[i]);
  }

  /* Compute bnorm, the norm of the data vector b, for reference. */
  double bnorm=0.0;
  for(int mi=0; mi<m; mi++) bnorm+=b[mi]*b[mi];
  bnorm=sqrt(bnorm); if(verbose>2) printf("  initial_bnorm := %g\n", bnorm);


  /*
   *  Main loop
   */
  int skipStep2=0;
  double obj=0.0;
  int iact=0; // The component x[iact] is the one that most wants to become active

  for(*iter=1; *iter<=maxIter; ++(*iter)) {

    if(verbose>1) {printf("iteration %d\n", *iter); fflush(stdout);}

    if(!skipStep2) {
      /* Step 2. */ if(verbose>1) {printf("  step 2\n"); fflush(stdout);}

      /*  Initialize the negative gradient vector w(*). */
      for(int ni=0; ni<n; ni++) w[ni]=0.0;

      /* Compute the residual vector b-a.x , the negative gradient vector w[*], and 
         the current objective value obj = || a.x - b ||.
         The residual vector is stored in the mm+1'st column of act[*,*]. */
      obj=0.0;
      for(int mi=0; mi<m; mi++) {
        double ri=b[mi];
        for(int ni=0; ni<n; ni++) ri-=a[mi+ni*m]*x[ni];
        obj+=ri*ri;
        for(int ni=0; ni<n; ni++) w[ni]+=a[mi+ni*m]*ri;
        act[mi+mm*m]=ri;
      }
      if(verbose>3) {printf("    obj := %g\n", obj); fflush(stdout);}

      /* Converged?  Stop if the misfit << || b ||, or if all components are active 
        (unless this is the first iteration from a 'warm start'). */
      if((sqrt(obj)<=bnorm*eps) || ((*iter)>1 && nbound==0)) {
        if(verbose>1) {printf("bvls converged.\n"); fflush(stdout);}
        istate[n]=nbound;
        w[0]=sqrt(obj);
        return(0);
      }

      /* Add the contribution of the active components back into the residual. */
      for(int ni=nbound; ni<n; ni++) {
        int i=abs(istate[ni])-1;
        for(int mi=0; mi<m; mi++) act[mi+mm*m] += a[mi+i*m]*x[i];
      }
      if(verbose>9) {
        printf("Residual vector:\n");
        for(int mi=0; mi<m; mi++) printf("\t%g", act[mi+mm*m]);
        printf("\n");
      }

    }



    /* The first iteration in a 'warm start' requires immediate QR in Step 6
       but mostly we want go through Steps 3-5 first . */
    if(key!=0 && (*iter)==1) {
      if(verbose>1) printf("  'warm start' requires immediate QR in Step 6\n");
    } else {

      int it; // variable indicating the element in istate that most wants to be active

      do {
        /* Steps 3, 4. */ if(verbose>1) {printf("  steps 3 and 4\n"); fflush(stdout);}
        /* Find the bound element that most wants to be active. */
        double worst=0.0;
        it=1;
        for(int ni=0; ni<nbound; ni++) {
          int i=abs(istate[ni])-1;
          double bad; if(istate[ni] < 0) bad=-w[i]; else bad=+w[i];
          if(bad < worst) { it=ni+1; worst=bad; iact=i; }
        }

        /* Test whether the Kuhn-Tucker condition is met. */
        if(worst>=0.0) {
          if(verbose>1) {printf("Kuhn-Tucker condition is met.\n"); fflush(stdout);}
          istate[n]=nbound;
          w[0]=sqrt(obj);
          return(0);
        }

        /* The component x[iact] is the one that most wants to become active.
           If the last successful change in the active set was to move x[iact] to a bound, 
           don't let x[iact] in now: set the derivative of the misfit with respect to x[iact] 
           to zero and return to the Kuhn-Tucker test. */
        if(iact==(aindx-1)) w[aindx-1]=0.0;
      } while(iact==(aindx-1)); // Step 3 again

      /* Step 5. */ if(verbose>1) {printf("  step 5\n"); fflush(stdout);}

      /* Undo the effect of the new (potentially) active variable on the residual vector. */
      if(istate[it-1]==0) { // remove this if never happening
        if(verbose>0) fprintf(stderr, "Error: BVLS istate is zero!\n"); 
        return(1);
      }
      {
        double bnd;
        if(istate[it-1]>0) bnd=bu[iact]; else bnd=bl[iact];
        for(int mi=0; mi<m; mi++) act[mi+mm*m]+=bnd*a[mi+iact*m];
      }

      /* Set flag istateFromStep5, indicating that Step 6 was entered from Step 5.
         This forms the basis of a test for instability: the gradient calculation shows that x[iact] 
         wants to join the active set; if QR puts x[iact] beyond the bound from which it came, 
         the gradient calculation was in error and the variable should not have been introduced. */
      istateFromStep5=istate[it-1]; 

      /* Swap the indices (in istate) of the new active variable and the rightmost bound variable; 
         `unbind' that location by decrementing nbound. */
      istate[it-1]=istate[nbound-1];
      nbound--;  nact++;
      istate[nbound]=1+iact;
      if(mm<nact) {
        if(verbose>0) fprintf(stderr, "Error: too many free variables in BVLS.\n");
        return(2);
      }

    } // finalized steps 3-5

    do {

      skipStep2=0;

      /* Step 6. */ if(verbose>1) {printf("  step 6\n"); fflush(stdout);}

      /* Load array act with the appropriate columns of A for QR. For added stability, reverse 
         the column ordering so that the most recent addition to the active set is in the last 
         column. Also copy the residual vector from act[., mm] into act[., mm+1]. */
      for(int mi=0; mi<m; mi++) {
        act[mi+(mm+1)*m]=act[mi+mm*m]; // vector b for QR
        for(int ni=nbound; ni<n; ni++) {
          int i=abs(istate[ni])-1;
          act[mi+(nact+nbound-ni-1)*m]=a[mi+i*m];
        }
      }
      if(verbose>9) {
        printf("Matrix A for QR:\n");
        for(int ni=0; ni<nact; ni++) {
          for(int mi=0; mi<m; mi++) printf("\t%g", act[mi+ni*nact]);
          printf("\n");
        }
        printf("Vector B for QR:\n");
        for(int mi=0; mi<m; mi++) printf("\t%g", act[(mm+1)*m + mi]);
        printf("\n");
      }

      /* Test for linear dependence in QR, and for an instability that moves the variable 
         just introduced away from the feasible region (rather than into the region or 
         all the way through it). In either case, remove the latest vector introduced from 
         the active set and adjust the residual vector accordingly.
         Set the gradient component (w[iact]) to zero and return to the Kuhn-Tucker test. */
      double r2;
      if(qrLH(m, nact, act, &act[(mm+1)*m], zz, &r2) !=0 || 
         (istateFromStep5>0 && zz[nact-1]>bu[iact]) || 
         (istateFromStep5<0 && zz[nact-1]<bl[iact]) )
      {
        nbound++;
        if(bu[iact]>x[iact]) istate[nbound-1]=-istate[nbound-1];
        nact--;
        for(int mi=0; mi<m; mi++) act[mi+mm*m]-=x[iact]*a[mi+iact*m];
        istateFromStep5 = 0; // not from step 5
        w[iact]=0.0;
        skipStep2=1; // we want to skip Step 2 and go directly to Step 3
        if(verbose>3) {printf("    going from step 6 to step 3\n"); fflush(stdout);}
        break; // go to step 3
      }

      /* If Step 6 was entered from Step 5 and we are here, a new variable has been successfully 
         introduced into the active set; the last variable that was fixed at a bound is again 
         permitted to become active. */
      if(istateFromStep5!=0) aindx=0;
      istateFromStep5=0;

      /* Step 7. */ if(verbose>1) {printf("  step 7\n"); fflush(stdout);}
      /* Check for strict feasibility of the new QR solution. */
      int qr_solution_feasible=1;
      int indexHolder=0;
      if(verbose>8) printf("    nact=%d  nbound=%d\n", nact, nbound);
      for(int ni=0; ni<nact; ni++) {
        indexHolder=ni; // Loop in step 8 will start from this
        int i=abs(istate[ni+nbound])-1;
        if(verbose>8) {
          printf("      istate[%d]=%d\n", ni+nbound, 1+i);
          printf("      zz[%d]=%g  bl[%d]=%g  bu[%d]=%g\n", nact-ni-1, zz[nact-ni-1], i, bl[i], i, bu[i]);
        }
        if(zz[nact-ni-1]<bl[i] || zz[nact-ni-1]>bu[i]) {
          if(verbose>3) {printf("    new iterate is not feasible\n"); fflush(stdout);}
          qr_solution_feasible=0; break; // go to Step 8
        }
      }
      if(verbose>8) printf("    indexHolder=%d\n", indexHolder);
      if(qr_solution_feasible) {
        if(verbose>3) {printf("    new iterate is feasible\n"); fflush(stdout);}
        for(int ni=0; ni<nact; ni++) {
          int i=abs(istate[ni+nbound])-1;
          x[i]=zz[nact-ni-1];
        }
        /* New iterate is feasible; back to the top. */
        break; // Back to the start of the main loop
      }

      { // keep local variables alpha and alf local
        double alpha=2.0;
        /* Steps 8 and 9 */ if(verbose>1) {printf("  steps 8 and 9\n"); fflush(stdout);}
        double alf=alpha;
        for(int ni=indexHolder; ni<nact; ni++) {
          int i=abs(istate[ni+nbound])-1;
          if(zz[nact-ni-1] > bu[i]) alf=(bu[i]-x[i])/(zz[nact-ni-1]-x[i]);
          if(zz[nact-ni-1] < bl[i]) alf=(bl[i]-x[i])/(zz[nact-ni-1]-x[i]);
          if(alf<alpha) {
            alpha=alf;
            aindx=1+i;
            if((zz[nact-ni-1]-bl[i])<0.0) aindxsign=-1; else aindxsign=+1;
          }
        }
        /* Step 10 */ if(verbose>1) {printf("  step 10\n"); fflush(stdout);}
        for(int ni=0; ni<nact; ni++) {
          int i=abs(istate[ni+nbound])-1;
          x[i]+=alpha*(zz[nact-ni-1]-x[i]);
        }
      }

      /* Step 11 */ if(verbose>1) {printf("  step 11\n"); fflush(stdout);}
      /* Move the variable that determined alpha to the appropriate bound
         (aindx is its index; sj is + if zz[aindx]> bu[aindx], - if zz[aindx]<bl[aindx] ).
         If any other component of  x  is infeasible at this stage, it must be due to round-off.
         Bind every infeasible component and every component at a bound to the appropriate bound.
         Correct the residual vector for any variables moved to bounds. Since at least one 
         variable is removed from the active set in this step, Loop B
         (Steps 6-11) terminates after at most nact steps. */
      {
        int noldb=nbound;
        for(int ni=0; ni<nact; ni++) {
          int i=abs(istate[ni+noldb])-1;
          if((bu[i]-x[i]<=0.0) || (i==(aindx-1) && aindxsign>0)) {
            /* Move x[i] to its upper bound. */
            x[i]=bu[i];
            istate[ni+noldb]=istate[nbound]; istate[nbound]=+(1+i); nbound++;
            for(int mi=0; mi<m; mi++) act[mi+mm*m]-=bu[i]*a[mi+i*m];
          } else if( ((x[i]-bl[i])<=0.0) || (i==(aindx-1) && aindxsign<0)) {
            /* Move x(j) to its lower bound. */
            x[i]=bl[i];
            istate[ni+noldb]=istate[nbound]; istate[nbound]=-(1+i); nbound++;
            for(int mi=0; mi<m; mi++) act[mi+mm*m]-=bl[i]*a[mi+i*m];
          }
        }
        nact=n-nbound;
      }
      /* If there are still active variables left, repeat the QR; if not, go back to step 6. */
    } while(nact>0);

  } // main loop

  /* iterMax reached */
  if(verbose>0) fprintf(stderr, "Error: BVLS fails to converge.\n");
  return(-1);
}
/*****************************************************************************/

/*****************************************************************************/
/** Algorithm for weighting the problem that is given to a LLSQ algorithm.

    Square roots of weights are used because in LLSQ algorithms the difference
    w*A-w*b is squared.

    @sa llsqWghtSquared
    @return Algorithm returns zero if successful, 1 if arguments are inappropriate.
*/
int llsqWght(
  /** Matrix A dimension N (nr of parameters). */
  int N,
  /** Matrix A dimension M (nr of samples). */
  int M,
  /** Pointer to matrix A[N][M]; enter NULL to use following matrix format instead. */
  double **A,
  /** Pointer to matrix A[N*M]; enter NULL to use previous matrix format instead. */
  double *a,
  /** Vector B of length M. */
  double *b,
  /** Weights for each sample (array of length M). */
  double *weight
) {
  int n, m;
  double *w;

  /* Check the arguments */
  if(N<1 || M<1 || (A==NULL && a==NULL) || b==NULL || weight==NULL) return(1);

  /* Allocate memory */
  w=(double*)malloc(M*sizeof(double)); if(w==NULL) return(2);

  /* Check that weights are not zero and get the square roots of them to w[] */
  for(m=0; m<M; m++) {
    if(weight[m]<=1.0e-20) w[m]=0.0;
    else w[m]=sqrt(weight[m]);
  }
 
  /* Multiply rows of matrix A and elements of vector b with weights*/
  for(m=0; m<M; m++) {
    for(n=0; n<N; n++) {
      if(A!=NULL) A[n][m]*=w[m];
      if(a!=NULL) a[m+n*M]*=w[m];
    }
    b[m]*=w[m];
  }

  free(w);
  return(0);
}
/*****************************************************************************/

/*****************************************************************************/
/** Algorithm for weighting the problem that is given to a LLSQ algorithm.

    Square roots of weights are used because in LLSQ algorithms the difference
    w*A-w*b is squared.

    Here user must give squared weights; this makes calculation faster, when
    this function needs to be called many times. 

    @sa llsqWght
    @return Algorithm returns zero if successful, 1 if arguments are inappropriate.
*/
int llsqWghtSquared(
  /** Matrix A dimension N (nr of parameters). */
  int N,
  /** Matrix A dimension M (nr of samples). */
  int M,
  /** Pointer to matrix A[N][M]; enter NULL to use following matrix format instead. */
  double **A,
  /** Pointer to matrix A[N*M]; enter NULL to use previous matrix format instead. */
  double *a,
  /** Vector B of length M. */
  double *b,
  /** Squared weights for each sample (array of length M). */
  double *sweight
) {
  int n, m;

  /* Check the arguments */
  if(N<1 || M<1 || (A==NULL && a==NULL) || b==NULL || sweight==NULL) return(1);

  /* Multiply rows of matrix A and elements of vector b with weights*/
  for(m=0; m<M; m++) {
    for(n=0; n<N; n++) {
      if(A!=NULL) A[n][m]*=sweight[m];
      if(a!=NULL) a[m+n*M]*=sweight[m];
    }
    b[m]*=sweight[m];
  }

  return(0);
}
/*****************************************************************************/

/*****************************************************************************/
