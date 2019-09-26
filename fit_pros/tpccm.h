/** @file tpccm.h
 *  @brief Header file for libtpccm.
 *  @details Header file for compartmental model library.
 *  @author Vesa Oikonen
 *  @copyright (c) Turku PET Centre
 */
#ifndef _TPCCM_H_
#define _TPCCM_H_
/*****************************************************************************/

/*****************************************************************************/
#include "tpcclibConfig.h"
/*****************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
/*****************************************************************************/

/*****************************************************************************/
/* sim1cm */
/*****************************************************************************/
int simMBF(
  double *t, double *ci, const int nr, 
  const double k1, const double k2, const double Vfit, double *ct
);
int simC1(
  double *t, double *ca, const int nr, 
  const double k1, const double k2, double *ct
);
int simC1_i(
  double *t, double *cai, const int nr, 
  const double k1, const double k2, double *ct
);
/*****************************************************************************/
/* sim2cm */
/*****************************************************************************/
int simC2(
  double *t, double *ca, const int nr, 
  const double k1, const double k2, const double k3, const double k4,
  double *ct, double *cta, double *ctb
);
int simC2_i(
  double *t, double *cai, const int nr, 
  const double k1, const double k2, const double k3, const double k4,
  double *ct, double *cta, double *ctb
);
/*****************************************************************************/
/* sim3cms */
/*****************************************************************************/
int simC3s(
  double *t, double *ca, const int nr, double k1, double k2,
  double k3, double k4, double k5, double k6,
  double *ct, double *cta, double *ctb, double *ctc
);
int simC3vs(
  double *t, double *ca, double *cb, const int nr,
  const double k1, const double k2, const double k3, const double k4, 
  const double k5, const double k6,
  const double f, const double vb, const double fa,
  double *cpet, double *cta, double *ctb, double *ctc,
  double *ctab, double *ctvb
);
/*****************************************************************************/
/* sim3cmp */
/*****************************************************************************/
int simC3p(
  double *t, double *ca, const int nr, const double k1, const double k2,
  const double k3, const double k4, const double k5, const double k6,
  double *ct, double *cta, double *ctb, double *ctc
);
int simC3vp(
  double *t, double *ca, double *cb, const int nr,
  const double k1, const double k2, const double k3, const double k4, 
  const double k5, const double k6,
  const double f, const double vb, const double fa,
  double *cpet, double *cta, double *ctb, double *ctc,
  double *ctab, double *ctvb
);
/*****************************************************************************/
/* simkloss */
/*****************************************************************************/
int simC2l(
  double *t, double *ca, const int nr, const double k1, const double k2,
  const double k3, const double kLoss, double *ct, double *cta, double *ctb
);
int simC2vl(
  double *t, double *ca, double *cb, const int nr,
  const double k1, const double k2, const double k3, const double kL,
  const double f, const double vb, const double fa,
  double *cpet, double *cta, double *ctb,
  double *ctab, double *ctvb
);
int simC3vpKLoss(
  double *t, double *ca, double *cb, const int nr,
  const double k1, const double k2, const double k3, const double k4, 
  const double k5, const double k6, const double kLoss, 
  const double f, const double vb, const double fa,
  double *cpet, double *cta, double *ctb, double *ctc,
  double *ctab, double *ctvb
);
/*****************************************************************************/
/* simrtcm */
/*****************************************************************************/
int simRTCM(
  double *t, double *cr, const int nr, const double R1, const double k2,
  const double k3, const double k4, double *ct, double *cta, double *ctb
);
int simSRTM(
  double *t, double *cr, const int nr, const double R1, const double k2,
  const double BP, double *ct
);
int simTRTM(
  double *t, double *cr, const int nr, const double R1, const double k2,
  double k3, double *ct
);
/*****************************************************************************/
/* simdicm */
/*****************************************************************************/
int simC3DIvs(
  double *t, double *ca1, double *ca2, double *cb, const int nr,
  const double k1, const double k2, const double k3, const double k4, 
  const double k5, const double k6, const double k1b, const double k2b, 
  const double f, const double vb, const double fa,
  double *scpet, double *sct1, double *sct2, double *sct3, double *sct1b,
  double *sctab, double *sctvb
);
int simC4DIvp(
  double *t, double *ca1, double *ca2, double *cb, const int nr,
  const double k1, const double k2, const double k3, const double k4, 
  const double k5, const double k6, const double k7, const double km, 
  const double k1b, const double k2b, 
  const double f, const double vb, const double fa,
  double *scpet, double *sct1, double *sct2, double *sct3,
  double *sct1b, double *sctab, double *sctvb,
  const int verbose
);
int simC4DIvs(
  double *t, double *ca1, double *ca2, double *cb, const int nr,
  const double k1, const double k2, const double k3, const double k4, 
  const double k5, const double k6, const double k7, const double km, 
  const double k1b, const double k2b, 
  const double f, const double vb, const double fa,
  double *scpet, double *sct1, double *sct2, double *sct3,
  double *sct1b, double *sctab, double *sctvb,
  const int verbose
);
/*****************************************************************************/
/* simdispersion */
/*****************************************************************************/
int simDispersion(
  double *x, double *y, const int n,
  const double tau1, const double tau2, double *tmp
);
/*****************************************************************************/
/* simoxygen */
/*****************************************************************************/
int simOxygen(
  double *t, double *ca1, double *ca2, double *ca1i, double *ca2i, const int n,
  const double k1a, const double k2a, const double km, 
  const double k1b, const double k2b, 
  const double vb, const double fa,
  double *scpet, double *sct1, double *sct2, double *sctab, 
  double *sctvb1, double *sctvb2, double *scvb1, double *scvb2,
  const int verbose
);
/*****************************************************************************/
/* convolut */
/*****************************************************************************/
int convolve1D(double *data, const int n, double *kernel, const int m, double *out);
int simIsSteadyInterval(double *x, const int n, double *f);
/*****************************************************************************/
/* simblood */
/*****************************************************************************/
/** Parameters of input CM for a single compound (parent or metabolite).
    @sa icmparcInit, icmparcAddMetabolites, icmparcAllocateTACs, icmparcFree
 */
typedef struct ICMPARC {
  /** Compound name. */
  char name[256];
  /** Compound infusion start time, from outside of system to BV. */
  double Ti;
  /** Compound infusion duration, from outside of system to BV. */
  double Tdur;
  /** Compound infusion rate (step function height). */
  double Irate;
  /** Rate constant from BV to BA. */
  double k_BV_BA;
  /** Rate constant for extraction from BA to U (out of system). */
  double k_BA_U;
  /** Rate constant from BA to TF. */
  double k_BA_TF;
  /** Rate constant from BA to TS. */
  double k_BA_TS;
  /** Rate constant from TF to BV. */
  double k_TF_BV;
  /** Rate constant from TS to BV. */
  double k_TS_BV;
  /** Number of metabolites. */
  unsigned int mNr;
  /** Pointer to list of metabolites, with list length mNr. */
  struct ICMPARC *metabolite;
  /** Pointer to parent compound. */
  struct ICMPARC *parent;
  /** Rate constant of formation from parent in BV. */
  double kp_BV;
  /** Rate constant of formation from parent in TF. */
  double kp_TF;
  /** Rate constant of formation from parent in TS. */
  double kp_TS;
  /** Optional storage for BV TAC integral. */
  double *ic_BV;
  /** Optional storage for TS TAC integral. */
  double *ic_TS;
  /** Optional storage for TF TAC integral. */
  double *ic_TF;
  /** Optional storage for BA TAC. */
  double *c_BA;
  /** Optional storage for BV TAC. */
  double *c_BV;
  /** Optional storage for TS TAC. */
  double *c_TS;
  /** Optional storage for TF TAC. */
  double *c_TF;
} ICMPARC;

void icmparcInit(ICMPARC *d);
int icmparcAddMetabolites(ICMPARC *d, unsigned const int mNr);
int icmparcAllocateTACs(ICMPARC *d, unsigned const int sNr, const int sub);
void icmparcFree(ICMPARC *d);
int simBTAC(double *t, const unsigned int nr, ICMPARC *p, double *cb);
/*****************************************************************************/

/*****************************************************************************/
#endif /* _TPCCM_H_ */
