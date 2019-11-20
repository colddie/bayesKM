/// @file o2.c
/// @author Vesa Oikonen
/// @brief Default parameters and helper functions for oxygen metabolism.
///
/*****************************************************************************/
#include "libtpcmodel.h"
/*****************************************************************************/

/*****************************************************************************/
/* Default variables for oxygen */
/** Arterial oxygen saturation fraction */
double SaO2=0.97;
/** Half-saturation pressure p50 (kPa) for hemoglobin */
double p50Hb=3.6;
/** Half-saturation pressure p50 (kPa) p50 for myoglobin */
double p50Mb=0.319;
/** Hill coefficient n for hemoglobin */
double nHb=2.7;
/** Hemoglobin concentration in blood (mg/g) */
double cHb=150.0;
/** Myoglobin concentration in muscle (mg/g) */
double cMb=4.7;
/*****************************************************************************/

/*****************************************************************************/
/** @brief Calculates K1/k2 ratio for [O-15]O2 in muscle, based on OER.
    @return Returns K1/k2.
 */
double mo2k1k2(
  /** Oxygen extraction ratio; assumptions hold only in the range OER=0.2-0.9 */
  const double OER,
  /** Arterial oxygen saturation fraction; you can enter DEFAULT_SAO2 */
  const double SaO2,
  /** Half-saturation pressure p50 (kPa) for hemoglobin; 
      you can enter DEFAULT_P50HB */
  const double p50Hb,
  /** Half-saturation pressure p50 (kPa) p50 for myoglobin; you can enter 
      DEFAULT_P50MB */
  const double p50Mb,
  /** Hill coefficient n for hemoglobin; you can enter DEFAULT_NHB */
  const double nHb,
  /** Hemoglobin concentration in blood (mg/g); you can enter DEFAULT_CHB */
  const double cHb,
  /** Myoglobin concentration in muscle (mg/g); you can enter DEFAULT_CMB */
  const double cMb,
  /** Verbose level; if zero, then nothing is printed into stdout or stderr */
  const int verbose
) {
  if(verbose>0) printf("mo2k1k2()\n");
  if(verbose>1) {
    printf("input:\n");
    printf("  OER := %g\n", OER);
    printf("  SaO2 := %g\n", SaO2);
    printf("  p50Hb := %g\n", p50Hb);
    printf("  p50Mb := %g\n", p50Mb);
    printf("  nHb := %g\n", nHb);
    printf("  cHb := %g\n", cHb);
    printf("  cMb := %g\n", cMb);
  }
  double SHb=(1.0-OER)*SaO2;
  if(verbose>2) printf("SHb := %g\n", SHb);

  double pO2=p50Hb*pow(SHb/(1.0-SHb), 1.0/nHb);
  if(verbose>2) printf("pO2 := %g\n", pO2);

  double SMb=pO2/(pO2+p50Mb);
  if(verbose>2) printf("SMb := %g\n", SMb);

  double rO2=cMb/cHb;
  if(verbose>2) printf("rO2 := %g\n", rO2);

  double k1k2=rO2*(SMb/SHb);
  if(verbose>1) printf("K1/k2 := %g\n", k1k2);

  return(k1k2);
}
/*****************************************************************************/

/*****************************************************************************/
/** Calculates the partial pressure of oxygen in muscle, based on OER and K1/k2.
    @return Returns pO2 (kPa).
 */
double mo2pO2(
  /** Oxygen extraction ratio; assumptions hold only in the range OER=0.2-0.9 */
  const double OER,
  /** K1/k2 ratio for [O-15]O2 */
  const double K1k2,
  /** Arterial oxygen saturation fraction; you can enter DEFAULT_SAO2 */
  const double SaO2,
  /** Half-saturation pressure p50 (kPa) p50 for myoglobin; you can enter 
      DEFAULT_P50MB */
  const double p50Mb,
  /** Hemoglobin concentration in blood (mg/g); you can enter DEFAULT_CHB */
  const double cHb,
  /** Myoglobin concentration in muscle (mg/g); you can enter DEFAULT_CMB */
  const double cMb,
  /** Verbose level; if zero, then nothing is printed into stdout or stderr */
  const int verbose
) {
  if(verbose>0) printf("mo2pO2()\n");
  if(verbose>1) {
    printf("input:\n");
    printf("  OER := %g\n", OER);
    printf("  K1/k2 := %g\n", K1k2);
    printf("  SaO2 := %g\n", SaO2);
    printf("  cHb := %g\n", cHb);
    printf("  cMb := %g\n", cMb);
    printf("  p50Mb := %g\n", p50Mb);
  }
  double SHb=(1.0-OER)*SaO2;
  if(verbose>2) printf("SHb := %g\n", SHb);

  double rO2=cMb/cHb;
  if(verbose>2) printf("rO2 := %g\n", rO2);

  double SMb=K1k2*(1.0-OER)*SaO2/rO2;
  if(verbose>2) printf("SMb := %g\n", SMb);

  double pO2=SMb*p50Mb/(1.0-SMb);
  if(verbose>1) printf("pO2 := %g\n", pO2);

  return(pO2);
}
/*****************************************************************************/

/*****************************************************************************/
