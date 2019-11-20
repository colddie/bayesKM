/// @file bootstrap.c
/// @author Kaisa Sederholm, Vesa Oikonen
/// @brief Procedure for counting the confidence intervals and standard
///        deviations for estimates of parameters of compartmental PET models.
///
/// This method is based on the article "Estimation of component and
/// parametric distributions in spectral analysis" by Turkheimer,
/// Sokoloff, Bertoldo, Lucignani, Reivich, Jaggi and Schmidt in
/// J Cereb Blood Flow Metabol. 1998; 18:1211-1222.
/// Standard deviation calculation is based on the book
/// "An introduction to bootstrap" by Efron & Tibshirani, 1993,
/// Chapman & Hall, New York.
///
/// @todo Enable optionally to use bobyqa instead of powell.
/// @test Test with weights.
///
/*****************************************************************************/
#include "libtpcmodel.h"
/*****************************************************************************/
/// @cond
#ifndef RAND_MAX
#define RAND_MAX 32767
#endif
/*****************************************************************************/
/** local function definitions */
static int bootstrapQSort(const void *par1, const void *par2);
/*****************************************************************************/
int bs_parNr, bs_frameNr;
double *bs_parameter;
double *bs_uplim;
double *bs_lowlim;
double *bs_weight;
double (*bs_func)(int, double*, void*);
/*****************************************************************************/
/// @endcond

/*****************************************************************************/
/** Bootstrap method.

  Original weights are assumed to be inversely proportional to variance.
  Square root is used, because bootstrap assumes them to be proportional
  to standard deviation.
  If only standard deviation is needed then cLim1 and cLim2 can be set to
  be NULL, and if only the confidence limits are wanted then SD can be set
  to be NULL.
  This function will not set seed for random number generator, therefore,
  make sure that it is set in your program, for example with srand(time(NULL));. 

\return Return values:
  - 0, if ok.
  - 1 - 3 if some of the given parameters is not qualified.
  - 4, if Powell fails, and 5, if out of memory.

*/
int bootstrapr(
  /** Bootstrap iteration number (>=100), set to zero to use the default (200). */
  int iterNr,
  /** Vector to put the lower confidence limits to; NULL if not needed. */
  double *cLim1,
  /** Vector to put the upper confidence limits to; NULL if not needed. */
  double *cLim2,
  /** Vector to put the standard deviations to; NULL if not needed. */
  double *SD,
  /** Best parameter estimates (preserved by this function). */
  double *parameter,
  /** Lower limits for the parameters. */
  double *lowlim,
  /** Upper limits for the parameters. */
  double *uplim,
  /** Nr of samples in tissue TAC data. */
  int frameNr,
  /** Measured tissue TAC values (not modified; local copy is used). */
  double *origTac,
  /** Best fitted tissue TAC values (not modified; local copy is used). */
  double *fitTac,
  /** Pointer to (empty) tissue TAC vector, where bootstrapped TACs will be written in each
   *  bootstrap iteration, and which will be used by objective function to calculate WSS. */
  double *bsTac,
  /** Nr of parameters. */
  int parNr,
  /** sample weights. */
  double *weight,
  /** The object function. */
  double (*objf)(int, double*, void*),
  /** Pointer to a string (allocated for at least 64 chars) where error message
      or other execution status will be written; enter NULL, if not needed. */   
  char *status,
  /** Verbose level; if zero, then nothing is printed to stderr or stdout. */
  int verbose,

  /** Full sampling matrix add later */
  double *bmatrix
) {

  int ret, i, j, powellItNr, lowindex, upindex;
  double fret=0.0, *parMean, *chainptr, *chain, *delta, *bsFitTac, **matrix;
  double help, *error, *wError, *biasEst, *unbiasPar, *estInOrder;
  char bserrmsg[64];

  if(verbose>0)
    printf("%s(%d, ..., %d, ..., %d, ..., %d)\n", __func__, iterNr, frameNr, parNr, verbose);

  /* Checking the given parameters */
  if(status!=0) strcpy(status, "");
  if(iterNr<100) iterNr=200;
  if(frameNr<1 || parNr<1 ) {
    strcpy(bserrmsg, "either framenumber or parameternumber is negative");
    if(verbose>0) fprintf(stderr, "Error: %s.\n", bserrmsg);
    if(status!=0) strcpy(status, bserrmsg);
    return(1);
  }
  if(bsTac==NULL || parameter==NULL || objf==NULL || origTac==NULL ||
     fitTac==NULL)
  {
    strcpy(bserrmsg, "some of the given parameters are not pointing anywhere");
    if(verbose>0) fprintf(stderr, "Error: %s.\n", bserrmsg);
    if(status!=0) strcpy(status, bserrmsg);
    return(2);
  }
  for(i=0; i<parNr; i++) {
    if(lowlim[i]>uplim[i]) {
      strcpy(bserrmsg, "given limit values are not qualified");
      if(verbose>0) fprintf(stderr, "Error: %s.\n", bserrmsg);
      if(status!=0) strcpy(status, bserrmsg);
      return(3);
    }
  }

  /* Allocating memory */
  if(verbose>1) printf("  allocating memory\n");
  bsFitTac=(double*)malloc(frameNr*sizeof(double));
  bs_parameter=(double*)malloc(parNr*sizeof(double));
  bs_weight=(double*)malloc(frameNr*sizeof(double));
  delta=(double*)malloc(parNr*sizeof(double));
  parMean=(double*)malloc(parNr*sizeof(double));
  chain=(double*)malloc((iterNr*parNr)*sizeof(double));
  matrix=(double**)malloc(parNr*sizeof(double*));
  // bmatrix = (double*)malloc(iterNr*parNr*sizeof(double));   // add later
  error=(double*)malloc(frameNr*sizeof(double));
  wError=(double*)malloc(frameNr*sizeof(double));
  biasEst=(double*)malloc(parNr*sizeof(double));
  unbiasPar=(double*)malloc(parNr*sizeof(double));
  if(bsFitTac==NULL || bs_parameter==NULL || bs_weight==NULL || delta==NULL ||
     parMean==NULL || chain==NULL || matrix==NULL || error==NULL ||
     wError==NULL || biasEst==NULL || unbiasPar==NULL) {
    strcpy(bserrmsg, "out of memory");
    if(verbose>0) fprintf(stderr, "Error: %s.\n", bserrmsg);
    if(status!=0) strcpy(status, bserrmsg);
    return(5);
  }
  for(i=0, chainptr=chain; i<parNr; i++) {
    matrix[i]=(chainptr + i*iterNr);
  }

  bs_parNr=parNr;
  bs_frameNr=frameNr;
  bs_func=objf;
  bs_lowlim=lowlim;
  bs_uplim=uplim;
  for(i=0; i<bs_parNr; i++) {
    bs_parameter[i]=parameter[i];
  }

  /* copy data and check the weights */
  /* calculate the errors and weighted errors */
  if(verbose>2) printf("  calculating errors and weighted errors");
  for(i=0; i<bs_frameNr; i++) {
    bsFitTac[i]=fitTac[i];
    if(weight[i]<=0.0) bs_weight[i]=1; else bs_weight[i]=weight[i];
    error[i]=origTac[i]-bsFitTac[i];
    wError[i]=error[i]/sqrt(bs_weight[i]);
  }
  if(verbose>3) {
    printf("  weighted errors:\n  ");
    for(i=0; i<bs_frameNr; i++) printf("%g ", wError[i]);
    printf("\n");
  }

  /*
   *  bootstrap iterations
   */
  if(verbose>1) printf("  bootstrap iterations\n");
  if(verbose>4) printf("Bootstrap matrix:\n");
  int powellfailNr=0;

  for(i=0; i<iterNr; i++) {

    /* sample a new error distribution (to the pointer error) */
    for(j=0; j<bs_frameNr; j++) {
      error[j]=wError[(int)((bs_frameNr)*(double)rand()/((double)(RAND_MAX)+1))];
      bsTac[j]=bsFitTac[j]+bs_weight[j]*error[j];
    }

    /* Powell local search */
    for(j=0; j<bs_parNr; j++) {
      delta[j]=0.01*(bs_uplim[j]-bs_lowlim[j]);
      bs_parameter[j]=parameter[j];
    }
    powellItNr=400;
    ret=powell(bs_parameter, delta, bs_parNr, 0.00001, &powellItNr, &fret, bs_func, NULL, 0);
    if(ret>1 && ret!=3)	{
      sprintf(bserrmsg, "error %d in powell()", ret);
      if(verbose>0) fprintf(stderr, "Error: %s.\n", bserrmsg);
      if(status!=0) strcpy(status, bserrmsg);
      free(bsFitTac);
      free(bs_parameter);
      free(bs_weight);
      free(delta);
      free(parMean);
      free(matrix);
      free(chain);
      free(error);
      free(wError);
      free(biasEst);
      free(unbiasPar);
      return(4);
    }
    if(ret==3) { // powell sometimes fails, do not worry if not too often
      powellfailNr++;
    }

    for(j=0; j<bs_parNr; j++) {
      matrix[j][i]=bs_parameter[j];
    }
    if(verbose>4) {
      for(j=0; j<bs_parNr; j++) printf("%g ", bs_parameter[j]);
      printf("\n");
    }

  } /* end of bootstrap iterations */
  if(powellfailNr>(iterNr/3)) {
    sprintf(bserrmsg, "error %d in powell()", ret);
    if(verbose>0) fprintf(stderr, "Error: too often %s.\n", bserrmsg);
    if(status!=0) strcpy(status, bserrmsg);
    free(bsFitTac);
    free(bs_parameter);
    free(bs_weight);
    free(delta);
    free(parMean);
    free(matrix);
    free(chain);
    free(error);
    free(wError);
    free(biasEst);
    free(unbiasPar);
    return(4);
  }

  /* Computing the mean of each parameter and estimates for bias */
  if(verbose>1) printf("  computing parameter bias\n");
  for(i=0; i<bs_parNr; i++) {
    for(j=0, help=0.0; j<iterNr; j++) help+=matrix[i][j];
    parMean[i]=help/(double)iterNr;
    biasEst[i]=parMean[i]-parameter[i];
    /*unbiasPar[i]=parameter[i]-biasEst[i];*/
    if(verbose>1) {
      printf("parMean[%d] := %g\n", i, parMean[i]);
      printf("parameter[%d] := %g\n", i, parameter[i]);
      printf("biasEst[%d] := %g\n", i, biasEst[i]);
    }
  }

  /* Computing the standard deviation for each parameter estimate */
  if(SD!=NULL) {
    if(verbose>1) printf("Standard deviations:\n");
    for(i=0; i<bs_parNr; i++) {
      for(j=0, help=0; j<iterNr; j++) {
        help+=(matrix[i][j]-parMean[i])*(matrix[i][j]-parMean[i])/(iterNr-1);
      }
      SD[i]=sqrt(help); if(verbose>1) printf("  %g\n", SD[i]);
    }
  }

  if(cLim1!=NULL && cLim2!=NULL) {
    if(verbose>1) printf("Confidence intervals:\n");
    lowindex=(int)temp_roundf(0.025*iterNr);
    upindex=(int)temp_roundf(0.975*iterNr)-1;
    /* Computing the 95% confidence intervals for each parameter estimate */
    for(i=0; i<bs_parNr; i++) {

      /* Arrange estimate samples in ascending order */
      estInOrder=matrix[i];
      qsort(estInOrder, iterNr, sizeof(double), bootstrapQSort);

      /* Get the indices within which 95% of samples lie in */
      if(fabs(estInOrder[lowindex]-biasEst[i])<1e-99) cLim1[i]=0.0;
      else cLim1[i]=estInOrder[lowindex]-biasEst[i];
      if(fabs(estInOrder[upindex]-biasEst[i])<1e-99) cLim2[i]=0.0;
      else cLim2[i]=estInOrder[upindex]-biasEst[i];
      if(verbose>1) {
        printf("  %g - %g\n", cLim1[i], cLim2[i]);
      }
    }
    if(verbose>6) {
      printf("Sorted matrix\n");
      for(j=0; j<iterNr; j++) {
        for(i=0; i<bs_parNr; i++) printf("  %12.3e", matrix[i][j]);
        printf("\n");
      }
      printf("lowindex := %d\nupindex := %d\n", lowindex, upindex);
    }
  }

  //
  for(j=0; j<iterNr; j++){
     for(i=0; i<bs_parNr; i++) {
       bmatrix[bs_parNr*j+i]=matrix[i][j]; } }
  //    printf("return full sampling matrix... %f \n", bmatrix[bs_parNr*j+i]); } }
  //
  
  
  free(bsFitTac);
  free(bs_parameter);
  free(bs_weight);
  free(delta);
  free(parMean);
  free(matrix);
  free(chain);
  free(error);
  free(wError);
  free(biasEst);
  free(unbiasPar);
  // free(bmatrix);

  if(verbose>0) {printf("  end of bootstrap()\n");}
  return(0);

} /* end bootstrap() */
/*****************************************************************************/

/*****************************************************************************/
/// @cond
int bootstrapQSort(const void *par1, const void *par2)
{
  if( *((double*)par1) < *((double*)par2)) return(-1);
  else if( *((double*)par1) > *((double*)par2)) return(1);
  else return(0);
}
/// @endcond
/*****************************************************************************/

/*****************************************************************************/
