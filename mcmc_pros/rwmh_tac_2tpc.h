#ifndef _RWMH_TAC_2TPC_H_
#define _RWMH_TAC_2TPC_H_

#include "mcmc.hpp"
// #include "tpccm.h"
// #include "sim1cm.cpp"
// #include "sim2cm.cpp"
// #include "simrtcm.cpp"
// #include "simpct.cpp"
// #include "simPatlak.c"
// #include "simLogan.c"


extern "C" int rwmh_tac_2tpc(int argc, float * argv[]);


int simC1(
  /** Array of time values */
  double *t,
  /** Array of arterial activities */
  double *ca,
  /** Number of values in TACs */
  const int nr,
  /** Rate constant of the model */
  const double k1,
  /** Rate constant of the model */
  const double k2,
  /** Pointer for TAC array to be simulated; must be allocated */
  double *ct
);

int simC1_i(
  /** Array of time values */
  double *t,
  /** Array of AUC 0-t of arterial activities */
  double *cai,
  /** Number of values in TACs */
  const int nr,
  /** Rate constant of the model */
  const double k1,
  /** Rate constant of the model */
  const double k2,
  /** Pointer for TAC array to be simulated; must be allocated */
  double *ct
);


int simMBF(
  /** Array of time values */
  double *t,
  /** Input activities */
  double *ci,
  /** Number of values in TACs */
  const int nr,
  /** Apparent k1 */
  const double k1,
  /** Apparent k2 */
  const double k2,
  /** Vfit */
  const double Vfit,
  /** Pointer for TAC array to be simulated; must be allocated */
  double *ct
); 


int simC2(
  /** Array of time values */
  double *t,
  /** Array of arterial activities */
  double *ca,
  /** Number of values in TACs */
  const int nr,
  /** Rate constant of the model */
  const double k1,
  /** Rate constant of the model */
  const double k2,
  /** Rate constant of the model */
  const double k3,
  /** Rate constant of the model */
  const double k4,
  /** Pointer for TAC array to be simulated; must be allocated */
  double *ct,
  /** Pointer for 1st compartment TAC to be simulated, or NULL */
  double *cta,
  /** Pointer for 2nd compartment TAC to be simulated, or NULL */
  double *ctb
); 

int simC2_i(
  /** Array of time values */
  double *t,
  /** Array of AUC 0-t of arterial activities */
  double *cai,
  /** Number of values in TACs */
  const int nr,
  /** Rate constant of the model */
  const double k1,
  /** Rate constant of the model */
  const double k2,
  /** Rate constant of the model */
  const double k3,
  /** Rate constant of the model */
  const double k4,
  /** Pointer for TAC array to be simulated; must be allocated */
  double *ct,
  /** Pointer for 1st compartment TAC to be simulated, or NULL */
  double *cta,
  /** Pointer for 2nd compartment TAC to be simulated, or NULL */
  double *ctb
); 


int simRTCM(
  /** Array of time values */
  double *t,
  /** Reference region activities */
  double *cr,
  /** Number of values in TACs */
  const int nr,
  /** Ratio K1/K1' */
  const double R1,
  /** Rate constant of the model */
  const double k2,
  /** Rate constant of the model */
  const double k3,
  /** Rate constant of the model */
  const double k4,
  /** Pointer for TAC array to be simulated; must be allocated */
  double *ct,
  /** Pointer for 1st compartment TAC to be simulated, or NULL */
  double *cta,
  /** Pointer for 2nd compartment TAC to be simulated, or NULL */
  double *ctb
);

int simSRTM(
  /** Array of time values */
  double *t,
  /** Reference region activities */
  double *cr,
  /** Number of values in TACs */
  const int nr,
  /** Ratio K1/K1' */
  const double R1,
  /** Rate constant of the model */
  const double k2,
  /** Binding potential */
  const double BP,
  /** Pointer for TAC array to be simulated; must be allocated */
  double *ct
);

int simTRTM(
  /** Array of time values */
  double *t,
  /** Reference region activities */
  double *cr,
  /** Number of values in TACs */
  const int nr,
  /** Ratio K1/K1' */
  const double R1,
  /** Rate constant of the model */
  const double k2,
  /** Rate constant of the model */
  const double k3,
  /** Pointer for TAC array to be simulated; must be allocated */
  double *ct
);


int simpct
(
    double *ts,
    double *ctt,
    int    frameNr,
    double cbf,
    double mtt,
    double *tac
);

#endif  /** _RWMH_TAC_2TPC_H_ */