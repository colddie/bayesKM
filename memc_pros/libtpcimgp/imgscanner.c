/// @file imgscanner.c
/// @author Vesa Oikonen
/// @brief Scanner specific parameters for IMG data.
///
/*****************************************************************************/
#include "libtpcimgp.h"
/*****************************************************************************/

/*****************************************************************************/
/** Sets scanner specific parameters in IMG data.
    If possible, set image zoom before calling this.
\return Returns 0, if ok.
 */
int imgSetScanner(
  /** IMG data which is filled with scanner specific information */
  IMG *img,
  /** SCANNER_ECAT931, SCANNER_ADVANCE, SCANNER_HRPLUS, SCANNER_HRRT, SCANNER_STEVCT_PET,
      as defined in libtpcimgio.h */
  int scanner_type
) {
  int rayNr;
  
  if(img->status<IMG_STATUS_OCCUPIED) return(1);
  img->scanner=scanner_type;
  /* Set zoom to 1.0, if it is not set */
  if(img->zoom<=0.0) img->zoom=1.0;
  /* Then the others */
  if(scanner_type==SCANNER_ECAT931) {
    // dimz 15
    rayNr=192;
    img->axialFOV=108.;
    img->transaxialFOV=600.826;
    img->sampleDistance=3.12932;
    img->sizez=6.75;
  } else if(scanner_type==SCANNER_ADVANCE) {
    // dimz 35
    rayNr=281;
    img->axialFOV=153.;
    img->transaxialFOV=550.;
    img->sampleDistance=1.970177;
    img->sizez=4.25;
  } else if(scanner_type==SCANNER_HRPLUS) {
    rayNr=288;
    img->axialFOV=155.2;
    img->transaxialFOV=583.;
    img->sampleDistance=2.25; /* bin size */
    img->sizez=2.425;
  } else if(scanner_type==SCANNER_HRRT) {
    rayNr=256;
    img->axialFOV=252.28;
    img->transaxialFOV=312.;
    img->sampleDistance=1.08; /* bin size */
    img->sizez=img->sizex=img->sizey=1.218750;
  } else if(scanner_type==SCANNER_STEVCT_PET) {
    // dimz 47
    rayNr=0;
    img->sizez=3.27;
    img->sizex=img->sizey=5.46875;
/*} else if(scanner_type==SCANNER_DMI_PET) { // Aino
    rayNr=544;
    img->axialFOV=200.;
    img->transaxialFOV=700.;
    img->sizez=nan("");
    img->sizex=img->sizey=nan("");
*/
  } else
    return(2);
  /* If this is image, then set also pixel sizes */
  if(img->type==IMG_TYPE_IMAGE) {
    if(scanner_type!=SCANNER_HRRT && scanner_type!=SCANNER_STEVCT_PET) {
      img->sizex=img->sizey=
        img->sampleDistance*(float)rayNr / ( (float)img->dimx*img->zoom );
    } else {
      img->sizex=img->sizey=
        img->transaxialFOV / ( (float)img->dimx*img->zoom );
    }
  }
  
  return(0);
}
/*****************************************************************************/

/*****************************************************************************/
