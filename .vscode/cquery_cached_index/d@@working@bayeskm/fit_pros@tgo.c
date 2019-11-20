/// @file tgo.c
/// @author Kaisa Sederholm, Vesa Oikonen
/// @copyright (c) Turku PET Centre
/// @brief Topographical minimization algorithm.
///
/// TGO searches the global minimum of a function using clusterization. 
/// Calls a local minimization algorithm. 
/// Based on an algorithm by Aimo Torn and Sami Viitanen. See the article 
/// Topographical Global optimization in: C.A. Floudas and P.M. Pardalos (eds.) 
/// Recent advances in Global Optimization, Princeton University Press, 1992
/// or webpage www.abo.fi/~atorn/ProbAlg/Page53.html
///
/******************************************************************************/
#include "tgo.h"
/******************************************************************************/
#ifndef RAND_MAX
/** Max nr for random nr generator */
#define RAND_MAX 32767
#endif
/******************************************************************************/
#ifndef TGO_SAMPLNR  
/** Sample nr for TGO; must be even number. */
#define TGO_SAMPLNR 1000
#endif
/******************************************************************************/
/** Biased (1) or even (0) parameter distribution */
int TGO_SQUARED_TRANSF = 1;
/** Local optimization outside (0) or inside (1) iTGO */
int TGO_LOCAL_INSIDE = 1;
/** Local optimization is done using Powell-Brent (0) or Bobyqa (1). */
int TGO_LOCAL_OPT = 0;
/******************************************************************************/
/** Topographical minimization algorithm, that searches the global minimum of 
    a function using clusterization. Calls a local minimization 
    algorithm. Based on an algorithm by Aimo Torn and Sami Viitanen.
    @return Returns 0, if ok.
*/
extern "C" int tgo(
  /** Lower limits for the parameters */
  double *lowlim,
  /** Upper limits for the parameters */
  double *uplim,
  /** The object function */
  double (*objf)(int, double*, void*),
  /** Data to objective function; NULL if not needed */
  void *objfData,
  /** Dimension = nr of parameters */
  int dim,
  /** Nr of neighbours to investigate; enter large nr relative to samNr (below)
   *  if only few local minima are expected */
  int neighNr,
  /** Function value at global minimum */
  double *fmin,
  /** Global minimum = parameter estimates */
  double *gmin,
  /** Nr of points to sample in one iteration; enter larger samNr if nr of
   *  iterations (below) is small */
  int samNr,
  /** Nr of TGO iterations; enter 0 to use the default; enter 1 to run TGO
   *  just once ie TGO instead of iTGO. Large iteration nr is needed if
   *  samNr (above) would otherwise require too much memory. */
  int tgoNr,
  /** Verbose level; if zero, then nothing is printed into stdout or stderr */
  int verbose
) {
  int i, j, k, l, IDmin, itNr, samplNr, topoNr, badNr, nevals=0, ret;
  double *delta, temp, min, *tempp, deltaf, tol;
  TGO_POINT *sampled_points;
  int fixed_n, fitted_n;


  if(verbose>0) {printf("in tgo()\n"); fflush(stdout);}
  if(TGO_LOCAL_OPT==1) {
    if(verbose>0) printf("local optimization routine: bobyqa\n");
  } else {
    if(verbose>0) printf("local optimization routine: powell\n");
  }

  /* Check input */
  if(lowlim==NULL || uplim==NULL || objf==NULL || dim<=0) return(1);
  if(fmin==NULL || gmin==NULL) return(1);

  /* Check if any of parameters is fixed */
  for(i=0, fixed_n=0; i<dim; i++) if(uplim[i]<=lowlim[i]) fixed_n++;
  fitted_n=dim-fixed_n;
  if(verbose>1) printf("%d parameter(s) are fixed.\n", fixed_n);
  if(fitted_n<1) return(1);

  /* Continue input checking */
  if(samNr<=0) samplNr=TGO_SAMPLNR; else samplNr=samNr;
  if(samplNr&1) samplNr++; // If number is odd, then add 1
  if(neighNr>samplNr-1) neighNr=samplNr-1; // Make sure "neighNr" isn't too big
  if(tgoNr<1) tgoNr=fitted_n; // Set TGO iteration number, if not set by user
  if(verbose>1) {
    printf("samplNr := %d\n", samplNr);
    printf("neighNr := %d\n", neighNr);
    printf("tgoNr := %d\n", tgoNr);
    if(verbose>2) {
      printf("iTGO limits: [%g,%g]", lowlim[0], uplim[0]);
      for(i=1; i<dim; i++) printf(" [%g,%g]", lowlim[i], uplim[i]);
      printf("\n");
    }
#ifdef OMP_NUM_THREADS
    printf("OMP_NUM_THREADS := %d\n", OMP_NUM_THREADS);
#endif
    fflush(stdout);
  }

  /* Allocate memory */
  sampled_points=(TGO_POINT*)calloc(samplNr, sizeof(TGO_POINT));
  if(sampled_points==NULL) return(2);
  for(i=0; i<samplNr; i++) {
    sampled_points[i].topomin=0;
    sampled_points[i].fvalue=0.0;
  }
  delta=(double*)malloc(dim*sizeof(double));
  tempp=(double*)malloc(samplNr*sizeof(double));
  if(delta==NULL || tempp==NULL) {free(sampled_points); return(2);}

  /* Set seed for random number generator */
  //srand(15345); srand(time(NULL));
  drandSeed(1);

  /*
   *  Iterative TGO, or non-iterative if tgoNr==1
   */
  for(l=0; l<tgoNr; l++) {

    if(verbose>2) {printf("TGO Loop # %d: \n", l+1); fflush(stdout);}

    /*
     *  Sample N points in the feasible region and compute the object function
     *  values for those points which do not already have it.
     */
    if(TGO_SQUARED_TRANSF==1)
      tgoRandomParametersST(sampled_points, dim, samplNr, lowlim, uplim);
    else
      tgoRandomParameters(sampled_points, dim, samplNr, lowlim, uplim);
    badNr=0;
    for(i=0; i<samplNr; i++) if(sampled_points[i].topomin==0) {
      sampled_points[i].fvalue=objf(dim, sampled_points[i].par, objfData);
      /* If function return value was not normal then we'll try 
         later (twice) with new guesses */
      if(!isfinite(sampled_points[i].fvalue)) {
        badNr++;
        if(verbose>5) {
          printf("this point did not give normal return value:\n");
          for(k=0; k<dim; k++) printf("  %10.2e", sampled_points[i].par[k]);
          printf("\n");
        }
      }
    }
    if(verbose>4 && badNr>0) printf("Nr of bad points: %d\n", badNr);
    /* New guesses for bad points */
    k=0; while(k<2 && badNr>0) {
      badNr=0; k++;
      for(i=0; i<samplNr; i++) 
        if(sampled_points[i].topomin==0 && !isfinite(sampled_points[i].fvalue)) {
          /* sample a new random point */
          if(TGO_SQUARED_TRANSF==1)
            tgoRandomParametersST(sampled_points+i, dim, 1, lowlim, uplim);
          else
            tgoRandomParameters(sampled_points+i, dim, 1, lowlim, uplim);
          /* compute the object function value for that */
          sampled_points[i].fvalue=objf(dim, sampled_points[i].par, objfData);
          if(!isfinite(sampled_points[i].fvalue)) badNr++;
        }
      if(verbose>4 && badNr>0) printf("Nr of bad points: %d\n", badNr);
    }
    /* Print sampled points */
    if(verbose>6) {
      printf("Sampled points:\n");
      for(j=0; j<samplNr; j++) {
        printf("%d", j+1);
        for(i=0; i<dim; i++) printf(" %e ", sampled_points[j].par[i]);
        printf("=>%e\n", sampled_points[j].fvalue);
      }
       fflush(stdout);
    }
    /* Object functions values must be good for at least NeigNr points */
    if(l==0 && (samplNr-badNr)<=neighNr) {
      if(verbose>0) {
        printf("Error in TGO: invalid function return value from all points.\n");
        fflush(stdout);
      }
      free(sampled_points); free(delta); free(tempp);
      return(3);
    }

    /*
     *  For each point i find out if it is a "topografic minimum" 
     *  = better than k neighbour points
     */
    /* Save the distances to point i in the vector tempp */
    /* Find the closest neighbour from {x1,..,xn}/{xi} k times, */
    /* extracting it after comparing */
    for(i=0, topoNr=0; i<samplNr; i++) {
      sampled_points[i].topomin=0;
      /* If function value is not ok, then point cannot be minimum */
      if(!isfinite(sampled_points[i].fvalue)) continue;

      /* Compute the (scaled) distances */
      for(j=0; j<samplNr; j++) {
        tempp[j]=1.0E+99;
        if(i!=j) {
          for(k=0, tempp[j]=0.0; k<dim; k++) {
            deltaf=uplim[k]-lowlim[k]; if(deltaf<=0.0) continue;
            temp=sampled_points[i].par[k]-sampled_points[j].par[k];
            if(deltaf>1.0E-20) temp/=deltaf;
            if(isfinite(temp)) tempp[j] += temp*temp;
          }
          /* Distance is computed as square root, but it does not affect
             the order of distances, and sqrt() is relatively slow */
          /*tempp=sqrt(tempp);*/
        }
      }

      /* Find the closest neighbours */
      /* At the same time, collect info for max fvalue of the neighbours
         and for the mean distance to every direction */
      for(k=0; k<dim; k++) sampled_points[i].delta[k]=0.0; // Init delta array 
      sampled_points[i].fvalrange=sampled_points[i].fvalue;
      for(j=0; j<neighNr; j++) {
        min=tempp[0]; IDmin=0;
        for(k=1; k<samplNr; k++) {if(tempp[k]<min) {min=tempp[k]; IDmin=k;}}
        tempp[IDmin]=1e+99; // so that this will not be used again
        /* If point i is worse than any of the closest neighbours, then
           point i is not a topographic minimum; then stop this loop and go to
           the next point i+1 */	  
        if(isfinite(sampled_points[IDmin].fvalue) &&
          sampled_points[IDmin].fvalue<sampled_points[i].fvalue) break;

        /* Sum the distances to every direction for delta calculation */
        for(k=0; k<dim; k++)
          sampled_points[i].delta[k]+=
            fabs(sampled_points[i].par[k]-sampled_points[IDmin].par[k]); 
        if(isfinite(sampled_points[IDmin].fvalue) &&
           sampled_points[IDmin].fvalue>sampled_points[i].fvalrange)
          sampled_points[i].fvalrange=sampled_points[IDmin].fvalue;
      } // next neighbour point

      /* If this was NOT a topographic minimum, then continue with the next i */
      if(j!=neighNr) continue;
      /* otherwise mark this as topografic minimum (TM) */
      sampled_points[i].topomin=1; topoNr++;

      /* Compute the mean distance of neighbours from the TM in each
         dimension; local optimization delta will be based on that */
      for(k=0; k<dim; k++) sampled_points[i].delta[k]/=(double)neighNr;
      /* Compute the max range in fvalues */
      sampled_points[i].fvalrange-=sampled_points[i].fvalue;

    } /* next sample */
    if(verbose>2) {printf("  %d topographical minima\n", topoNr); fflush(stdout);}

    /* Check that if no topographical minimum was found, then set the smallest */
    /* minimum as 'topographical' minimum */
    if(topoNr==0) {
      min=sampled_points[0].fvalue; IDmin=0;
      for(k=1; k<samplNr; k++)
        if(!isfinite(min) || sampled_points[k].fvalue<min) {
          min=sampled_points[k].fvalue; IDmin=k;}
      sampled_points[IDmin].topomin=1;
      for(k=0; k<dim; k++)
        sampled_points[IDmin].delta[k]=0.1*(uplim[k]-lowlim[k]);
      sampled_points[i].fvalrange+=100.0*fabs(sampled_points[i].fvalue);
      if(verbose>2) {
        printf("  ; therefore minimum was set to point %d at %e\n", 
          IDmin, sampled_points[IDmin].fvalue);
         fflush(stdout);
      }
      topoNr=1;
    }
    if(verbose>3) { // Print the best TM 
      for(i=0, min=1e+99, IDmin=0; i<samplNr; i++)
        if(sampled_points[i].topomin==1) {
          if(isfinite(sampled_points[i].fvalue) && sampled_points[i].fvalue<min)
          {
            min=sampled_points[i].fvalue; IDmin=i;
          }
        }
      printf("  best topographical min:"); fflush(stdout);
      for(k=0; k<dim; k++) printf(" %e", sampled_points[IDmin].par[k]);
      printf(" => %e\n", sampled_points[IDmin].fvalue); fflush(stdout);
    }

    if(TGO_LOCAL_INSIDE==1) {
      /* Local optimization for each TM */
      if(verbose>2) printf("local optimization for each TM\n");
      for(i=0; i<samplNr; i++) if(sampled_points[i].topomin==1) {
        //tol=sampled_points[IDmin].fvalue*1.0E-02;
        //for(k=0; k<dim; k++) delta[k]=sampled_points[i].delta[k];
        for(k=0; k<dim; k++) delta[k]=0.1*sampled_points[i].delta[k];
        if(verbose>3) printf("point %d: original fvalue=%.2e\n",
                             i+1, sampled_points[i].fvalue);
        if(TGO_LOCAL_OPT==1) {
          tol=1.0E-08; //tol=1.0E-09;
          ret=bobyqa(dim, 0, sampled_points[i].par, lowlim, uplim, delta, 0.0, 1.0E-03, 
                     1.0E-10, tol, tol, 2000, &nevals, 
                     &sampled_points[i].fvalue, objf, objfData, NULL, verbose-3);
          if(ret<0 && verbose>0) {
            printf("bobyqa error %d\n", ret); fflush(stdout);}
          if(ret<0 && ret!=BOBYQA_ROUNDOFF_LIMITED) {
            free(sampled_points); free(delta); free(tempp);
            return(5);
          }
          if(verbose>3) {
            printf("  local opt => %.2e (nr of evals=%d)\n",
                   sampled_points[i].fvalue, nevals);
            fflush(stdout);
          }
        } else {
          itNr=40; //itNr=100;
          tol=1.0E-03; //tol=2.0E-03; //tol=1.0E-09;
          POWELL_LINMIN_MAXIT=30; // 100;
          ret=powell(sampled_points[i].par, delta, dim, tol, &itNr,
                     &sampled_points[i].fvalue, objf, objfData, verbose-3);
          if(ret>1 && verbose>0) {printf("powell error %d\n", ret); fflush(stdout);}
          if(ret>3) {
            free(sampled_points); free(delta); free(tempp);
            return(5);
          }
          if(verbose>3) {
            printf("  local opt => %.2e (itNr=%d)\n",
                   sampled_points[i].fvalue, itNr);
            fflush(stdout);
          }
        }
      }
    } // end of local optimizations inside this iTGO loop

  } /* end of tgo iterations */

  if(verbose>1) { // Print the final topographical minima
    if(verbose>2) printf("Final topographical minima and deltas\n");
    else printf("Final topographical minima\n");
    for(i=0; i<samplNr; i++) if(sampled_points[i].topomin==1) {
      k=0; printf("  %3d: %.2e", i, sampled_points[i].par[k]);
      for(k=1; k<dim; k++) printf(" %.2e", sampled_points[i].par[k]);
      printf(" => %.2e\n", sampled_points[i].fvalue); fflush(stdout);
      if(verbose>2) {
        k=0; printf("       %.2e", sampled_points[i].delta[k]);
        for(k=1; k<dim; k++) printf(" %.2e", sampled_points[i].delta[k]);
        printf(" => %.2e\n", sampled_points[i].fvalrange); fflush(stdout);
      }
    }
  }

  if(TGO_LOCAL_INSIDE==0) { 
    /*
     *  Use the points in TM as starting points for local optimization;
     *  this first local opt is done only if not done already inside iTGO 
     */
    if(verbose>2) {printf("Topographic minima:\n"); fflush(stdout);}
    for(i=0; i<samplNr; i++) if(sampled_points[i].topomin==1) {
      //tol=sampled_points[IDmin].fvalue*1.0E-02;
      for(k=0; k<dim; k++) {
        if(verbose>2) {printf("%e ", sampled_points[i].par[k]); fflush(stdout);}
        delta[k]=0.1*sampled_points[i].delta[k];
      }
      if(verbose>3) {printf("=> %e ", sampled_points[i].fvalue); fflush(stdout);}
      if(TGO_LOCAL_OPT==1) {
        tol=1.0E-09; //tol=1.0E-09;
        ret=bobyqa(dim, 0, sampled_points[i].par, lowlim, uplim, delta, 0.0, 1.0E-03, 
                   1.0E-10, tol, tol, 2000, &nevals, 
                   &sampled_points[i].fvalue, objf, objfData, NULL, verbose-3);
        if(ret<0 && verbose>0) {
          printf("bobyqa error %d\n", ret); fflush(stdout);}
        if(ret<0 && ret!=BOBYQA_ROUNDOFF_LIMITED) {
          free(sampled_points); free(delta); free(tempp);
          return(5);
        }
        if(verbose>2) {
          printf("local opt 1st round point %d => %e (nr of evals=%d)\n",
            i+1, sampled_points[i].fvalue, nevals);
          fflush(stdout);
        }
      } else {
        itNr=50; //itNr=100; 
        tol=1.0E-03; //tol=1.0E-04;
        POWELL_LINMIN_MAXIT=60; // POWELL_LINMIN_MAXIT=50; // 100;
        ret=powell(sampled_points[i].par, delta, dim, tol, &itNr,
                 &sampled_points[i].fvalue, objf, objfData, verbose-3);
        if(ret>1 && verbose>0) {printf("powell error %d\n", ret); fflush(stdout);}
        if(ret>3) {
          if(verbose>0) {printf("powell error %d\n", ret); fflush(stdout);}
          free(sampled_points); free(delta); free(tempp);
          return(5);
        }
        if(verbose>2) {
          printf("=> %e (itNr=%d) ", sampled_points[i].fvalue, itNr);
          fflush(stdout);
        }
      }
    }
    if(verbose>0) { // Print the topographical minima after local optimization
      printf("Final topographical minima after local optimization\n");
      for(i=0; i<samplNr; i++) if(sampled_points[i].topomin==1) {
        k=0; printf("  %3d: %.2e", i, sampled_points[i].par[k]);
        for(k=1; k<dim; k++) printf(" %.2e", sampled_points[i].par[k]);
        printf(" => %.2e\n", sampled_points[i].fvalue); fflush(stdout);
      }
    }
  }

  /* Rerun of local optimization with smaller tolerance and delta */
  for(i=0; i<samplNr; i++) if(sampled_points[i].topomin==1) {
    //tol=sampled_points[IDmin].fvalue*1.0E-04;
    for(k=0; k<dim; k++) delta[k]=0.1*sampled_points[i].delta[k];
    if(TGO_LOCAL_OPT==1) {
      tol=1.0E-10;
      ret=bobyqa(dim, 0, sampled_points[i].par, lowlim, uplim, delta, 0.0, 1.0E-05, 
                 1.0E-10, tol, tol, 1000, &nevals, 
                 &sampled_points[i].fvalue, objf, objfData, NULL, verbose-3);
      if(ret<0 && verbose>0) {
        printf("bobyqa error %d\n", ret); fflush(stdout);}
      if(ret<0 && ret!=BOBYQA_ROUNDOFF_LIMITED) {
        free(sampled_points); free(delta); free(tempp);
        return(5);
      }
      if(verbose>2) {
        printf("local opt 2nd round point %d => %e (nr of evals=%d)\n",
          i+1, sampled_points[i].fvalue, nevals);
        fflush(stdout);
      }
    } else {
      itNr=40; /*itNr=100;*/
      tol=1.0E-04; //tol=1.0E-03; //tol=1.0E-08;
      POWELL_LINMIN_MAXIT=60; // 100;
      ret=powell(sampled_points[i].par, delta, dim, tol, &itNr,
               &sampled_points[i].fvalue, objf, objfData, verbose-3);
      if(ret>1 && verbose>0) {printf("powell error %d\n", ret); fflush(stdout);}
      if(ret>3) {
        free(sampled_points); free(delta); free(tempp);
        return(6);
      }
      if(verbose>2) {
        printf("=> %e (itNr=%d)\n", sampled_points[i].fvalue, itNr); 
        fflush(stdout);
      }
    }
  }

  if(verbose>0) { // Print the topographical minima after 2nd local optimization
    printf("Final topographical minima after 2nd local optimization\n");
    for(i=0; i<samplNr; i++) if(sampled_points[i].topomin==1) {
      k=0; printf("  %3d: %.2e", i, sampled_points[i].par[k]);
      for(k=1; k<dim; k++) printf(" %.2e", sampled_points[i].par[k]);
      printf(" => %.2e\n", sampled_points[i].fvalue); fflush(stdout);
    }
  }


  /*
   *  Find the best locally optimized TM and run local opt again from
   *  this point with better accuracy
   */

  for(i=0, min=1e+99, IDmin=0; i<samplNr; i++) if(sampled_points[i].topomin==1) {
    if(isfinite(sampled_points[i].fvalue) && sampled_points[i].fvalue<min) {
      min=sampled_points[i].fvalue; IDmin=i;}
  }
  if(verbose>1) {
    printf("Best topographical minimum:");
    for(k=0; k<dim; k++) printf("%e ", sampled_points[IDmin].par[k]);
    printf("-> %e \n", sampled_points[IDmin].fvalue); fflush(stdout); 
  }

  /* Rerun of local optimization to the best point */
  deltaf=0.01; tol=5.0E-04; //tol=1.0E-15;
  do {
    for(k=0; k<dim; k++) delta[k]=deltaf*sampled_points[IDmin].delta[k];
    //for(k=0; k<dim; k++) delta[k]=deltaf*(uplim[k]-lowlim[k]);
    if(TGO_LOCAL_OPT==1) {
      ret=bobyqa(dim, 0, sampled_points[IDmin].par, lowlim, uplim, delta, 0.0, tol, // !!! 
                 1.0E-10, tol, tol, 5000, &nevals, 
                 &sampled_points[IDmin].fvalue, objf, objfData, NULL, verbose-3);
      if(ret<0 && verbose>0) {
        printf("bobyqa error %d\n", ret); fflush(stdout);}
      if(ret<0 && ret!=BOBYQA_ROUNDOFF_LIMITED) {
        free(sampled_points); free(delta); free(tempp);
        return(5);
      }
      if(verbose>2) {
        printf("local opt of the best point with tol=%g => %e (nr of evals=%d)\n",
          tol, sampled_points[IDmin].fvalue, nevals);
        fflush(stdout);
      }
      itNr=nevals;
      deltaf*=0.5; tol*=0.1;
    } else {
      itNr=100;
      POWELL_LINMIN_MAXIT=100; // 100;
      ret=powell(sampled_points[IDmin].par, delta, dim, tol, &itNr,
                 &sampled_points[IDmin].fvalue, objf, objfData, verbose-3);
      if(ret>1 && verbose>0) {printf("powell error %d\n", ret); fflush(stdout);}
      if(ret>3) {
        free(sampled_points); free(delta); free(tempp);
        return(7);
      }
      if(verbose>1) {
        printf("  powell once more with %d iteration(s) -> WSS=%e\n",
          itNr, sampled_points[IDmin].fvalue);
          fflush(stdout);
      }
      //deltaf*=0.5; tol*=0.5;
      deltaf*=0.5; tol*=0.25;
    }
  } while(itNr>1 && deltaf>1.0E-05);
//} while(itNr>1 && deltaf>1.0E-15);


  /* Save best point to gmin */
  for(k=0; k<dim; k++) gmin[k]=sampled_points[IDmin].par[k];
  *fmin=sampled_points[IDmin].fvalue;
  if(!isfinite(sampled_points[IDmin].fvalue)) {
    if(verbose>0) {printf("TGO error: valid minimum value was not reached.\n");
    fflush(stdout);}
    free(sampled_points); free(delta); free(tempp); return(9);
  }

  /* Exit TGO */
  free(sampled_points); free(delta); free(tempp);
  if(verbose>0) {printf("out of tgo\n"); fflush(stdout);}
  return(0);
} /* end tgo */
/******************************************************************************/

/******************************************************************************/
/** Create randomized parameters for TGO */
extern "C" void tgoRandomParameters(
  /** Pointer to list of TGO points */
  TGO_POINT *p,
  /** Nr of parameters in the point */
  int parNr,
  /** Nr of TGO points */
  int sNr,
  /** List of lower limits for each parameter. */ 
  double *low,
  /** List of upper limits for each parameter. */ 
  double *up
) {
  int i, j;
  double dif;

  for(j=0; j<parNr; j++) {
    dif=up[j]-low[j];
    if(dif<=0.0) {
      for(i=0; i<sNr; i++) if(p[i].topomin==0) p[i].par[j]=low[j];
    } else {
      for(i=0; i<sNr; i++) if(p[i].topomin==0) {
        //p[i].par[j]=((double)rand()/(double)RAND_MAX) * dif + low[j];
        p[i].par[j]= drand()*dif + low[j];
      }
    }
  }
}
/******************************************************************************/
/** Create randomized parameters for TGO with square-root transformation,
 *  that is, parameter distribution is biased towards low absolute value. */
extern "C" void tgoRandomParametersST(
  /** Pointer to list of TGO points */
  TGO_POINT *p,
  /** Nr of parameters in the point */
  int parNr,
  /** Nr of TGO points */
  int sNr,
  /** List of lower limits for each parameter. */ 
  double *low,
  /** List of upper limits for each parameter. */ 
  double *up
) {
  int i, j;
  double v, stl, stu, dif;

  for(j=0; j<parNr; j++) {
    dif=up[j]-low[j];
    if(dif<=0.0) {
      for(i=0; i<sNr; i++) if(p[i].topomin==0) p[i].par[j]=low[j];
    } else {
      stl=copysign(sqrt(fabs(low[j])),low[j]); if(!isnormal(stl)) stl=0.0;
      stu=copysign(sqrt(fabs(up[j])), up[j]); if(!isnormal(stu)) stu=0.0;
      dif=stu-stl;
      for(i=0; i<sNr; i++) if(p[i].topomin==0) {
        v=drand()*dif + stl;
        p[i].par[j]=copysign(v*v, v);
      }
    }
  }
}
/******************************************************************************/

/******************************************************************************/
