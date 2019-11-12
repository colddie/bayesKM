/// @file sif.c
/// @author Vesa Oikonen
/// @brief Routines for Scan Information Files (SIF).
///
/*****************************************************************************/
#include "libtpcimgio.h"
/*****************************************************************************/

/*****************************************************************************/
/*!
   Initiate SIF structure. This should be called once before first use.
   @sa sifEmpty, sifRead, sifSetmem
 */
void sifInit(
  /** Pointer to sif data struct. */
  SIF *data
) {
  if(SIF_TEST) printf("sifInit()\n");
  if(data==NULL) return;
  memset(data, 0, sizeof(SIF));
  data->frameNr=data->colNr=0;
}
/*****************************************************************************/

/*****************************************************************************/
/*!
   Free memory allocated for SIF. All contents are destroyed.
   @sa sifInit, sifSetmem
 */
void sifEmpty(
  /** Pointer to sif data struct. */
  SIF *data
) {
  if(SIF_TEST) printf("sifEmpty()\n");
  if(data==NULL) return;
  if(data->frameNr>0) {
    free((char*)(data->x1)); free((char*)(data->x2));
    free((char*)(data->prompts)); free((char*)(data->randoms));
    free((char*)(data->trues)); free((char*)(data->weights));
    data->frameNr=data->colNr=0;
  }
  data->scantime=(time_t)0; data->version=0;
  strcpy(data->studynr, ""); strcpy(data->isotope_name, "");
}
/*****************************************************************************/

/*****************************************************************************/
/*!
   Allocates memory for SIF data.
   @sa sifInit, sifEmpty
   @return 0 if ok, 1 failed memory allocation
 */
int sifSetmem(
  /** Pointer to initiated sif data struct. Any existing data is destroyed. */
  SIF *data,
  /** Number of PET time frames. */
  int frameNr
) {
  if(SIF_TEST) printf("sifSetmem()\n");
  if(data==NULL) return(1);
  /* Clear previous data, if necessary */
  if(data->frameNr>0) sifEmpty(data);
  if(frameNr<1) return(0);
  
  /* Allocate memory */
  data->x1=(double*)calloc(frameNr, sizeof(double));
  data->x2=(double*)calloc(frameNr, sizeof(double));
  data->prompts=(double*)calloc(frameNr, sizeof(double));
  data->randoms=(double*)calloc(frameNr, sizeof(double));
  data->trues=(double*)calloc(frameNr, sizeof(double));
  data->weights=(double*)calloc(frameNr, sizeof(double));
  if(data->x1==NULL || data->x2==NULL || data->prompts==NULL ||
     data->randoms==NULL || data->trues==NULL || data->weights==NULL) {
    strcpy(siferrmsg, "out of memory"); return(1);}
  data->frameNr=frameNr;

  return(0);
}
/*****************************************************************************/

/*****************************************************************************/
