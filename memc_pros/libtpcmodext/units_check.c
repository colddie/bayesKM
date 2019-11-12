/// @file units_check.c
/// @brief Check and set data units for PET modelling.
/// @author Vesa Oikonen
///
/*****************************************************************************/

/*****************************************************************************/
#include "libtpcmodext.h"
/*****************************************************************************/

/*****************************************************************************/
/** Check that calibration units in IMG (PET image) and DFT (input TAC) are
 *  the same, and if not, then try to convert DFT calibration unit to IMG unit.
 *  If input unit is unknown, then assume it is the same as the PET unit.
\return Returns 0 if successful, >0 in case of error, and <0 in case of
    a warning or error message to user is suggested.
 */
int cunit_check_dft_vs_img(
  /** Pointer to DFT */
  DFT *dft,
  /** Pointer to IMG */
  IMG *img,
  /** Char pointer to string (at least of length 128) where possible
      error description or warning is copied; set to NULL if not necessary */
  char *errmsg,
  /** Verbose level; if zero, then nothing is printed to stderr or stdout */
  int verbose
) {
  int iunit, punit;

  if(verbose>0) printf("calibration_unit_check_dft_vs_img()\n");
  if(errmsg!=NULL) sprintf(errmsg, "program error");
  if(dft==NULL || img==NULL) return 1;

  iunit=petCunitId(dft->unit); // identify input file unit
  punit=img->unit;

  if(iunit==CUNIT_UNKNOWN) { // Input file unit is unknown
    // If PET unit is not known either, give a warning
    if(punit==CUNIT_UNKNOWN) {
      if(errmsg!=NULL) sprintf(errmsg, "unknown concentration units");
      return -1;
    } else { // Set to PET unit, and give a warning
      if(errmsg!=NULL)
        sprintf(errmsg, "unknown input concentration unit, now set to PET unit");
      strcpy(dft->unit, imgUnit(img->unit));
      return -2;
    }
  }

  // Input unit is known; if PET unit is not, then give a warning
  if(punit==CUNIT_UNKNOWN) {
    if(errmsg!=NULL) sprintf(errmsg, "unknown concentration units in PET data");
    return -3;
  }

  // Both units are known, so convert input data if necessary/possible
  if(iunit==CUNIT_KBQ_PER_ML) { // input is in units kBq/ml
    // If PET unit is the same, then everything is fine
    if(punit==CUNIT_KBQ_PER_ML) {
      if(errmsg!=NULL)
        sprintf(errmsg, "input and PET data have the same concentration units.\n");
      return 0;
    } else if(punit==CUNIT_BQ_PER_ML) { // image is in Bq/ml, convert input
      dftUnitConversion(dft, CUNIT_BQ_PER_ML);
      if(errmsg!=NULL)
        sprintf(errmsg, "input units converted to %s\n", dft->unit);
      return 0;
    } else { // image is in some other units, just give a warning for now
      if(errmsg!=NULL)
        sprintf(errmsg, "different concentration units in input and PET data");
      return -4;
    }
  } else if(iunit==CUNIT_BQ_PER_ML) { // input is in units Bq/ml
    // If PET unit is the same, then everything is fine
    if(punit==CUNIT_BQ_PER_ML) {
      if(errmsg!=NULL)
        sprintf(errmsg, "input and PET data have the same concentration units.\n");
      return 0;
    } else if(punit==CUNIT_KBQ_PER_ML) { // image is in kBq/ml, convert input
      dftUnitConversion(dft, CUNIT_KBQ_PER_ML);
      if(errmsg!=NULL)
        sprintf(errmsg, "input units converted to %s\n", dft->unit);
      return 0;
    } else { // image is in some other units, just give a warning for now
      if(errmsg!=NULL)
        sprintf(errmsg, "different concentration units in input and PET data");
      return -4;
    }
  } else { // input unit is known, but not kBq/ml or Bq/ml
    if(errmsg!=NULL)
      sprintf(errmsg, "check the concentration units in input and PET data");
    return -5;
  }

  return 0;
}
/*****************************************************************************/

/*****************************************************************************/
