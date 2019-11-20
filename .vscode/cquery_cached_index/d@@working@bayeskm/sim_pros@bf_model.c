/// @file bf_model.c
/// @brief Functions for calculation of basis functions for PET modelling.
/// @author Vesa Oikonen
///
/*****************************************************************************/

/*****************************************************************************/
#include "libtpcmodext.h"
/*****************************************************************************/

/*****************************************************************************/
/** Calculates set of basis functions for SRTM.
    @return Returns 0 if successful, otherwise non-zero.
 */
int bf_srtm(
  /** PET frame mid times */
  double *t,
  /** Non-decay corrected Cr(t) */
  double *cr,
  /** Nr of PET frames */
  int n,
  /** Nr of basis functions to calculate */
  int bfNr,
  /** theta3 min */
  double t3min,
  /** theta3 max */
  double t3max,
  /** data for basis functions is allocated and filled here */
  DFT *bf
) {
  int bi, fi, ret;
  double a, b, c;

  /* Check the parameters */
  if(t==NULL || cr==NULL || n<2 || bfNr<1 || t3min<1.0E-10 || t3min>=t3max) return(1);
  if(bf==NULL || bf->voiNr>0) return(1);
  
  /* Allocate meory for basis functions */
  ret=dftSetmem(bf, n, bfNr); if(ret) return(2);

  /* Copy and set information fields */
  bf->voiNr=bfNr; bf->frameNr=n;
  bf->_type=DFT_FORMAT_STANDARD;
  for(bi=0; bi<bf->voiNr; bi++) {
    snprintf(bf->voi[bi].voiname, 6, "B%5.5d", bi+1);
    strcpy(bf->voi[bi].hemisphere, ".");
    strcpy(bf->voi[bi].place, ".");
    strcpy(bf->voi[bi].name, bf->voi[bi].voiname);
  }
  for(fi=0; fi<bf->frameNr; fi++) bf->x[fi]=t[fi];

  /* Compute theta3 values to size fields */
  a=log10(t3min); b=log10(t3max); c=(b-a)/(double)(bfNr-1);
  for(bi=0; bi<bf->voiNr; bi++) {
    bf->voi[bi].size=pow(10.0, (double)bi*c+a);
  }
  
  /* Calculate the functions */
  for(bi=0; bi<bf->voiNr; bi++) {
    a=bf->voi[bi].size;
    ret=simC1_v1(t, cr, n, 1.0, a, bf->voi[bi].y);
    if(ret) return(4);
  }

  return(0);
}
/*****************************************************************************/

/*****************************************************************************/
/** Calculates set of basis functions for generic radiowater model.
    @return Returns 0 if successful, otherwise non-zero.
 */
int bfRadiowater(
  /** Arterial blood input TAC (not modified). */
  DFT *input,
  /** PET TACs (not modified, just to get frame times). */
  DFT *tissue,
  /** Place for basis functions (initiated DFT struct, allocated and filled here). */
  DFT *bf,
  /** Nr of basis functions to calculate. */
  int bfNr,
  /** Minimum of k2 (sec-1 or min-1, corresponding to TAC time units). */
  double k2min,
  /** Maximum of k2 (sec-1 or min-1, corresponding to TAC time units). */
  double k2max,
  /** Pointer to a string (allocated for at least 64 chars) where error message
      or other execution status will be written; enter NULL, if not needed. */   
  char *status,
  /** Verbose level; if zero, then nothing is printed to stderr or stdout. */
  int verbose
) {
  int bi, fi, ret;
  double a, b, c;

  if(verbose>0)
    printf("\nbfRadiowater(*inp, *tis, *bf, %d, %g, %g, status, %d)\n",
           bfNr, k2min, k2max, verbose);

  /* Check the parameters */
  if(input==NULL || tissue==NULL || bf==NULL) {
    if(status!=NULL) strcpy(status, "program error");
    return 1;
  }
  if(input->frameNr<3 || input->voiNr<1) {
    if(status!=NULL) strcpy(status, "no input data");
    return 2;
  }
  if(tissue->frameNr<1) {
    if(status!=NULL) strcpy(status, "no pet data");
    return 3;
  }
  if(input->timeunit!=tissue->timeunit) {
    if(status!=NULL) strcpy(status, "invalid time units");
    return 4;
  }
  if(bfNr<2) {
    if(status!=NULL) strcpy(status, "invalid nr of basis functions");
    return 5;
  }
  if(k2min<1.0E-10) k2min=1.0E-10; // range calculation does not work otherwise
  if(k2min>=k2max || k2min<0.0) {
    if(status!=NULL) strcpy(status, "invalid k2 range");
    return 6;
  }
  if(verbose>1) {
    printf("input timerange: %g - %g\n", input->x[0], input->x[input->frameNr-1]);
    printf("tissue timerange: %g - %g\n", tissue->x[0], tissue->x[tissue->frameNr-1]);
  }
  
  /* Allocate memory for basis functions */
  if(verbose>1) printf("allocating memory for basis functions\n");
  ret=dftSetmem(bf, tissue->frameNr, bfNr);
  if(ret) {
    if(status!=NULL) strcpy(status, "out of memory");
    return 10;
  }

  /* Copy and set information fields */
  bf->voiNr=bfNr; bf->frameNr=tissue->frameNr;
  bf->_type=tissue->_type;
  dftCopymainhdr2(tissue, bf, 1);
  for(bi=0; bi<bf->voiNr; bi++) {
    snprintf(bf->voi[bi].voiname, 6, "B%5.5d", bi+1);
    strcpy(bf->voi[bi].hemisphere, ".");
    strcpy(bf->voi[bi].place, ".");
    strcpy(bf->voi[bi].name, bf->voi[bi].voiname);
  }
  for(fi=0; fi<bf->frameNr; fi++) {
    bf->x[fi]=tissue->x[fi];
    bf->x1[fi]=tissue->x1[fi];
    bf->x2[fi]=tissue->x2[fi];
  }
  
  /* Compute the range of k2 values to size fields */
  if(verbose>1) printf("computing k2 values\n");
  a=log10(k2min); b=log10(k2max); c=(b-a)/(double)(bfNr-1);
  if(verbose>20) printf("a=%g b=%g, c=%g\n", a, b, c);
  for(bi=0; bi<bf->voiNr; bi++) {
    bf->voi[bi].size=pow(10.0, (double)bi*c+a);
  }
  if(verbose>2) {
    printf("final BF k2 range: %g - %g\n", 
           bf->voi[0].size, bf->voi[bf->voiNr-1].size);
  }
  
  /* Allocate memory for simulated TAC */
  double *sim;
  sim=(double*)malloc(input->frameNr*sizeof(double));
  if(sim==NULL) {
    if(status!=NULL) strcpy(status, "out of memory");
    dftEmpty(bf); return 11;
  }
  
  /* Calculate the basis functions at input time points */
  if(verbose>1) printf("computing basis functions at input sample times\n");
  for(bi=0; bi<bf->voiNr; bi++) {
    a=bf->voi[bi].size;
    ret=simC1_v1(input->x, input->voi[0].y, input->frameNr,
               1.0, a, sim);
    if(ret) {
      if(status!=NULL) strcpy(status, "simulation problem");
      free(sim); dftEmpty(bf);
      return(20);
    }
    if(verbose>100) {
      printf("\nk2 := %g\n", a);
      printf("simulated TAC:\n");
      for(fi=0; fi<input->frameNr; fi++)
        printf("  %12.6f  %12.3f\n", input->x[fi], sim[fi]);
    }
    /* interpolate to PET time frames */
    if(tissue->timetype==DFT_TIME_STARTEND)
      ret=interpolate4pet(input->x, sim, input->frameNr, tissue->x1, tissue->x2,
                          bf->voi[bi].y, NULL, NULL, bf->frameNr);
    else
      ret=interpolate(input->x, sim, input->frameNr, tissue->x,
                      bf->voi[bi].y, NULL, NULL, bf->frameNr);
    if(ret) {
      if(status!=NULL) strcpy(status, "simulation problem");
      free(sim); dftEmpty(bf);
      return(20);
    }

  } // next basis function

  free(sim);
  if(verbose>1) printf("bfRadiowater() done.\n\n");
  if(status!=NULL) strcpy(status, "ok");
  return(0);
}
/*****************************************************************************/

/*****************************************************************************/
/** Calculates set of basis functions for irreversible 2TCM.
   @return Returns 0 if successful, otherwise non-zero.
 */
int bfIrr2TCM(
  /** Arterial PTAC (not modified). */
  DFT *input,
  /** PET TTAC (not modified, just to get frame times). */
  DFT *tissue,
  /** Place for basis functions (initiated DFT struct, allocated and filled here). */
  DFT *bf,
  /** Nr of basis functions to calculate. */
  int bfNr,
  /** Minimum of theta=k2+k3 (sec-1 or min-1, corresponding to TAC time units). */
  double thetamin,
  /** Maximum of theta=k2+k3 (sec-1 or min-1, corresponding to TAC time units). */
  double thetamax,
  /** Pointer to a string (allocated for at least 64 chars) where error message
      or other execution status will be written; enter NULL, if not needed. */
  char *status,
  /** Verbose level; if zero, then nothing is printed to stderr or stdout. */
  int verbose
) {
  int bi, fi, ret;
  double a, b, c;

  if(verbose>0)
    printf("\nbfIrr2TCM(*inp, *tis, *bf, %d, %g, %g, status, %d)\n",
           bfNr, thetamin, thetamax, verbose);

  /* Check the parameters */
  if(input==NULL || tissue==NULL || bf==NULL) {
    if(status!=NULL) strcpy(status, "program error");
    else if(verbose>0) fprintf(stderr, "invalid function parameters\n");
    return(1);
  }
  if(input->frameNr<3 || input->voiNr<1) {
    if(status!=NULL) strcpy(status, "no input data");
    else if(verbose>0) fprintf(stderr, "invalid input data\n");
    return(2);
  }
  if(tissue->frameNr<1) {
    if(status!=NULL) strcpy(status, "no pet data");
    else if(verbose>0) fprintf(stderr, "invalid PET data\n");
    return(3);
  }
  if(input->timeunit!=tissue->timeunit) {
    if(status!=NULL) strcpy(status, "invalid time units");
    else if(verbose>0) fprintf(stderr, "invalid time units\n");
    return(4);
  }
  if(bfNr<2) {
    if(status!=NULL) strcpy(status, "invalid nr of basis functions");
    else if(verbose>0) fprintf(stderr, "invalid number of basis functions\n");
    return(5);
  }
  if(thetamin<0.0) thetamin=0.0;
  if(thetamin>=thetamax) {
    if(status!=NULL) strcpy(status, "invalid theta range");
    else if(verbose>0) fprintf(stderr, "invalid theta range\n");
    return(6);
  }
  if(verbose>1) {
    printf("input timerange: %g - %g\n", input->x[0], input->x[input->frameNr-1]);
    printf("tissue timerange: %g - %g\n", tissue->x[0], tissue->x[tissue->frameNr-1]);
  }
  
  /* Allocate memory for basis functions */
  if(verbose>1) printf("allocating memory for basis functions\n");
  ret=dftSetmem(bf, tissue->frameNr, bfNr);
  if(ret) {
    if(status!=NULL) strcpy(status, "out of memory");
    else if(verbose>0) fprintf(stderr, "out of memory\n");
    return(10);
  }

  /* Copy and set information fields */
  bf->voiNr=bfNr; bf->frameNr=tissue->frameNr;
  bf->_type=tissue->_type;
  dftCopymainhdr2(tissue, bf, 1);
  for(bi=0; bi<bf->voiNr; bi++) {
    snprintf(bf->voi[bi].voiname, 6, "B%5.5d", bi+1);
    strcpy(bf->voi[bi].hemisphere, ".");
    strcpy(bf->voi[bi].place, ".");
    strcpy(bf->voi[bi].name, bf->voi[bi].voiname);
  }
  for(fi=0; fi<bf->frameNr; fi++) {
    bf->x[fi]=tissue->x[fi];
    bf->x1[fi]=tissue->x1[fi];
    bf->x2[fi]=tissue->x2[fi];
  }
  
  /* Compute the range of theta values to size fields */
  if(verbose>1) printf("computing theta values\n");
  a=thetamin; b=thetamax; c=(b-a)/(double)(bfNr-1);
  if(verbose>20) printf("a=%g b=%g, c=%g\n", a, b, c);
  for(bi=0; bi<bf->voiNr; bi++) bf->voi[bi].size=(double)bi*c+a;
  if(verbose>2) {
    printf("final BF theta range: %g - %g\n", bf->voi[0].size, bf->voi[bf->voiNr-1].size);
    printf("theta step size: %g\n", c);
  }
  
  /* Allocate memory for simulated TAC */
  double *sim;
  sim=(double*)malloc(input->frameNr*sizeof(double));
  if(sim==NULL) {
    if(status!=NULL) strcpy(status, "out of memory");
    else if(verbose>0) fprintf(stderr, "out of memory\n");
    dftEmpty(bf); return(11);
  }
  
  /* Calculate the basis functions at input time points */
  if(verbose>1) printf("computing basis functions at input sample times\n");
  for(bi=0; bi<bf->voiNr; bi++) {
    a=bf->voi[bi].size;
    ret=simC1_v1(input->x, input->voi[0].y, input->frameNr, 1.0, a, sim);
    if(ret) {
      if(status!=NULL) strcpy(status, "simulation problem");
      else if(verbose>0) fprintf(stderr, "simulation problem\n");
      free(sim); dftEmpty(bf); return(20);
    }
    if(verbose>100) {
      printf("\ntheta := %g\n", a);
      printf("simulated TAC:\n");
      for(fi=0; fi<input->frameNr; fi++)
        printf("  %12.6f  %12.3f\n", input->x[fi], sim[fi]);
    }
    /* interpolate to PET time frames */
    if(tissue->timetype==DFT_TIME_STARTEND)
      ret=interpolate4pet(input->x, sim, input->frameNr, tissue->x1, tissue->x2,
                          bf->voi[bi].y, NULL, NULL, bf->frameNr);
    else
      ret=interpolate(input->x, sim, input->frameNr, tissue->x,
                      bf->voi[bi].y, NULL, NULL, bf->frameNr);
    if(ret) {
      if(status!=NULL) strcpy(status, "simulation problem");
      else if(verbose>0) fprintf(stderr, "simulation problem\n");
      free(sim); dftEmpty(bf); return(20);
    }

  } // next basis function

  free(sim);
  if(verbose>1) printf("bfIrr2TCM() done.\n\n");
  if(status!=NULL) strcpy(status, "ok");
  return(0);
}
/*****************************************************************************/

/*****************************************************************************/
