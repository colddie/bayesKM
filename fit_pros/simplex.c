/// @file simplex.c
/// @author Vesa Oikonen
/// @brief Nelder-Mead algorithm (Downhill simplex) for function minimization.
///
/*****************************************************************************/
#include "libtpcmodel.h"
/*****************************************************************************/
/// @cond
/** Local functions */
void    _simplexGenNew(int M, double F);
/** Local variables for this routine                                          */
int     _simplexParNr, Worst, NewPnt;
double  _simplexP[MAX_PARAMETERS+3][MAX_PARAMETERS],
        _simplexC[MAX_PARAMETERS], _simplexR[MAX_PARAMETERS+3];
double  (*_simplexFunc)(double*);
/// @endcond
/*****************************************************************************/

/********************* Simplex main function *********************************/
/** Downhill simplex function minimization routine.
    Note that if any constraints are required for the parameter
    values they must be set in the function.
    @return Function returns the least calculated value of func.
    @sa powell, tgo, bobyqa, nlopt1D
*/
double simplex(
  /** Pointer to the function to be minimized. It must be defined in main program as: 
      double func(double *p); where p is the parameter array. */
  double (*_fun)(double*),
  /** The number of unknown parameters */
  int parNr,
  /** This double array contains the minimized parameters.
      Initial values must be set. */
  double *par,
  /** This double array contains the initial changes to parameters.
      To fix a parameter, set the corresponding delta to 0. */
  double *delta,
  /** Maximal error allowed (stopping rule #1) */
  double maxerr,
  /** Maximal nr of iterations allowed (stopping rule #2) */
  int maxiter
) {
  int         i, j, Meas, it;
  double      Max, Min, Max2, Min2, LastChi;
  int         NextBest, New2, Best;


  if(SIMPLEX_TEST>0) printf("in simplex()\n");
  /* SetUp */
  _simplexFunc=_fun;
  _simplexParNr=parNr; it=0; NewPnt=_simplexParNr+1;
  for(i=0; i<_simplexParNr; i++)
    for(Meas=0; Meas<_simplexParNr+3; Meas++) _simplexP[Meas][i]=par[i];
  if(SIMPLEX_TEST) {
    for(i=0; i<_simplexParNr; i++)
      printf("%12g   %12g\n", _simplexP[0][i], delta[i]);
    printf("ChiSqr of guesses: %f\n", (*_simplexFunc)(_simplexP[0]));
  }
  New2=NewPnt+1;
  for(Meas=0; Meas<=_simplexParNr; Meas++) {
    it++;
    _simplexR[Meas] = (*_simplexFunc)(_simplexP[Meas]);
    for (i=0; i<_simplexParNr; i++) {
      if(i==Meas) delta[i]= -delta[i];
      _simplexP[Meas+1][i] = _simplexP[Meas][i] + delta[i];
    }
  }

  /* Simplex minimization */
  LastChi = 1.0E30; NextBest=Best=0;
  do {
    for(j=0; j<100; j++) {
      /* Find the max and min response measured */
      Max=0.; Min=1.0E30;
      for (i=0; i<=_simplexParNr; i++) {
        if(_simplexR[i] > Max) {Max=_simplexR[i]; Worst=i;}
        if(_simplexR[i] < Min) {Min=_simplexR[i]; Best=i; }
      }
      /* Find 2nd best and 2nd worst, too */
      Max2=0.; Min2=1.0E30;
      for (i=0; i<=_simplexParNr; i++) {
        if((_simplexR[i] > Max2) && (_simplexR[i] < Max)) Max2=_simplexR[i];
        if((_simplexR[i] < Min2) && (_simplexR[i] > Min)) {
          Min2=_simplexR[i]; NextBest=i;}
      }
      /* Calculate centroid of all measurements */
      for(i=0; i<_simplexParNr; i++) {
        _simplexC[i]=0.;
        for(Meas=0; Meas<=_simplexParNr; Meas++)
          if(Meas!=Worst) _simplexC[i]+=_simplexP[Meas][i];
        _simplexC[i]/=(double)_simplexParNr;
      }
      /* Measure the response at the point reflected away from worst */
      for(i=0; i<_simplexParNr; i++)
        _simplexP[NewPnt][i] = 2.*_simplexC[i] - _simplexP[Worst][i];
      _simplexR[NewPnt]= (*_simplexFunc)(_simplexP[NewPnt]);
      it++;
      /* If this one is better than previous best, then expand in this
          direction */
      if(_simplexR[NewPnt] < _simplexR[Best]) {
        _simplexGenNew(New2,2.0); it++;
      } else {
        /* If this one is worse than previous worst, measure point halfway
           between worst and centroid */
        if(_simplexR[NewPnt] > _simplexR[Worst]) {
          _simplexGenNew(New2,-0.5); it++;
        } else {
          /* If newest response is worse than next best point
             but better than worst, measure response halfway
             between centroid and newest point */
          if((_simplexR[NextBest] < _simplexR[NewPnt]) &&
             (_simplexR[NewPnt] < _simplexR[Worst])) {
            _simplexGenNew(New2,0.5); it++;
          } else {
            /* If none of the above, keep the new point as best */
            for(i=0; i<_simplexParNr; i++)
              _simplexP[Worst][i] = _simplexP[NewPnt][i];
            _simplexR[Worst] = _simplexR[NewPnt];
          }
        }
      }
    }
    if(SIMPLEX_TEST>0) printf(" it=%i; ChiSqr=%f\n", it, _simplexR[Best]);
    if(SIMPLEX_TEST>1)
      for(i=0; i<_simplexParNr; i++) printf("     %12g\n", _simplexP[Best][i]);
    /* Check if fitting is not proceeding */
    if(_simplexR[Best] == LastChi) {
      for(i=0; i<_simplexParNr; i++) par[i]=_simplexP[Best][i];
      return _simplexR[Best];
    }
    LastChi = _simplexR[Best];
  } while ((_simplexR[Best]>maxerr) && (it<=maxiter));

  for(i=0; i<_simplexParNr; i++) par[i]=_simplexP[Best][i];
  if(SIMPLEX_TEST>0) printf("out simplex()\n");
  return _simplexR[Best];
}
/*****************************************************************************/
/// \cond
/** _simplexGenNew() */
void _simplexGenNew(
  /** M */
  int M, 
  /** F!=1.0 */
  double F
) {
  int i;

  for(i=0; i<_simplexParNr; i++)
    _simplexP[M][i] = _simplexC[i] + F*(_simplexC[i]-_simplexP[Worst][i]);
  _simplexR[M] = (*_simplexFunc)(_simplexP[M]);
  if (_simplexR[M] < _simplexR[NewPnt]) {
    /*_simplexP[M][M]*/
    for(i=0; i<_simplexParNr; i++) _simplexP[Worst][i] = _simplexP[M][i];    
    _simplexR[Worst] = _simplexR[M];
  } else {
    for (i=0; i<_simplexParNr; i++) _simplexP[Worst][i] = _simplexP[NewPnt][i];
    _simplexR[Worst] = _simplexR[NewPnt];
  }
}
/*****************************************************************************/
/// \endcond
