/// @file img_ana.c
/// @author Vesa Oikonen
/// @brief I/O routines for IMG data from/to Analyze 7.5 format.
///
/*****************************************************************************/
#include "libtpcimgio.h"
/*****************************************************************************/

/*****************************************************************************/
/*!
 *   Read Analyze 7.5 image.
 *
 *   Analyze database name must be given with path. Image and header files
 *   with .img and .hdr extensions must exist. Also SIF file with .sif
 *   extension is used, if it exists.
 *   anaFlipping() determines whether image is flipped in z-direction;
 *   image is always flipped in x,y-directions.
 *
 * @param dbname Analyze database name with path, with or without extension
 * @param img Pointer to initialized IMG structure
 * @return 0 if ok, and otherwise IMG status code; sets IMG->statmsg
 *  in case of error.
 */
int imgReadAnalyze(const char *dbname, IMG *img) {
  FILE *fp;
  int ret, fi, pi, xi, yi;
  float *fdata=NULL, *fptr;
  ANALYZE_DSR dsr;
  char datfile[FILENAME_MAX], hdrfile[FILENAME_MAX], siffile[FILENAME_MAX];
  int dimNr, dimx, dimy, dimz=1, dimt=1, pxlNr=0;
  SIF sif;


  if(IMG_TEST) printf("imgReadAnalyze(%s, *img)\n", dbname);

  /* Check the arguments */
  imgSetStatus(img, STATUS_OK);
  if(img==NULL || img->status!=IMG_STATUS_INITIALIZED) {
    imgSetStatus(img, STATUS_FAULT); return(2);}
  if(dbname==NULL || !dbname[0]) {imgSetStatus(img, STATUS_FAULT); return(1);}
  
  /* Make the image and header filenames */
  ret=anaExistsNew(dbname, hdrfile, datfile, siffile);
  if(ret==0) {imgSetStatus(img, STATUS_NOHEADERFILE); return(3);}
  if(ret==1 && IMG_TEST>0) printf("no SIF found for %s\n", dbname); 

  /* Read Analyze header file */
  ret=anaReadHeader(hdrfile, &dsr);
  if(ret) {
    if(ret==1) imgSetStatus(img, STATUS_FAULT);
    else if(ret==2) imgSetStatus(img, STATUS_NOHEADERFILE);
    else imgSetStatus(img, STATUS_UNSUPPORTED);
    return(3);
  }
  if(IMG_TEST) anaPrintHeader(&dsr, stdout);

  /* Open image datafile */
  if(IMG_TEST) fprintf(stdout, "reading image data %s\n", datfile);
  if((fp=fopen(datfile, "rb")) == NULL) {
    imgSetStatus(img, STATUS_NOIMGDATA); return(5);}

  /* Prepare IMG for Analyze image */
  /* Get the image dimensions from header */
  dimNr=dsr.dime.dim[0];
  if(dimNr<2) {fclose(fp); imgSetStatus(img, STATUS_INVALIDHEADER); return(4);}
  dimx=dsr.dime.dim[1]; dimy=dsr.dime.dim[2];
  if(dimNr>2) {dimz=dsr.dime.dim[3]; if(dimNr>3) dimt=dsr.dime.dim[4];}
  pxlNr=dimx*dimy*dimz;
  if(pxlNr<1) {fclose(fp); imgSetStatus(img, STATUS_INVALIDHEADER); return(4);}
  /* Allocate memory for IMG */
  ret=imgAllocate(img, dimz, dimy, dimx, dimt);
  if(ret) {fclose(fp); imgSetStatus(img, STATUS_NOMEMORY); return(11);}
  /* Copy information from Analyze header */
  img->type=IMG_TYPE_IMAGE;
  strlcpy(img->studyNr, dsr.hist.patient_id, MAX_STUDYNR_LEN+1);
  if(strcmp(img->studyNr, ".")==0) strcpy(img->studyNr, "");
  strcpy(img->patientName, dsr.hist.patient_id);
  img->sizex=dsr.dime.pixdim[1];
  img->sizey=dsr.dime.pixdim[2];
  img->sizez=dsr.dime.pixdim[3];
  img->xform[0]=NIFTI_XFORM_UNKNOWN; // qform
  img->xform[1]=NIFTI_XFORM_SCANNER_ANAT; // sform
  /*if(dsr.dime.funused2>1.E-5) img->zoom=dsr.dime.funused2;*/
  if(dsr.dime.funused3>1.E-5) img->isotopeHalflife=dsr.dime.funused3;
  for(pi=0; pi<dimz; pi++) img->planeNumber[pi]=pi+1;
  if(dsr.little) img->_fileFormat=IMG_ANA_L; else img->_fileFormat=IMG_ANA;
  /* Decay correction */
  if(strstr(dsr.hist.descrip, "Decay corrected.")!=NULL)
    img->decayCorrection=IMG_DC_CORRECTED;
  else if(strstr(dsr.hist.descrip, "No decay correction.")!=NULL)
    img->decayCorrection=IMG_DC_NONCORRECTED;
  else
    img->decayCorrection=IMG_DC_CORRECTED; // just assumed so

  /* Allocate memory for one image frame */
  fdata=malloc(pxlNr*sizeof(float));
  if(fdata==NULL) {fclose(fp); imgSetStatus(img, STATUS_NOMEMORY); return(12);}

  /* Read one image frame at a time */
  for(fi=0; fi<dimt; fi++) {
    fptr=fdata;
    ret=anaReadImagedata(fp, &dsr, fi+1, fptr);
    if(ret) {
      free(fdata); fclose(fp); imgSetStatus(img, STATUS_NOIMGDATA); return(7);}
    /* Copy pixel values to IMG */
    fptr=fdata;
    if(anaFlipping()==0) { /* no flipping in z-direction */
      for(pi=0; pi<img->dimz; pi++)
        for(yi=dimy-1; yi>=0; yi--)
          for(xi=dimx-1; xi>=0; xi--)
            img->m[pi][yi][xi][fi]=*fptr++;
    } else {
      for(pi=dimz-1; pi>=0; pi--)
        for(yi=dimy-1; yi>=0; yi--)
          for(xi=dimx-1; xi>=0; xi--)
            img->m[pi][yi][xi][fi]=*fptr++;
    }
  } /* next frame */
  free(fdata);
  fclose(fp);
  
  /* Try to read frame time information from SIF file */
  /* Check if SIF file is found */
  if(siffile[0]) {
    if(IMG_TEST) printf("reading SIF file %s\n", siffile);
    if(access(siffile, 0) == -1) {
      if(IMG_TEST) printf(" No SIF file; therefore unknown frame times.\n");
      return(0);
    }
  }
  /* If found, then read it */
  sifInit(&sif); ret=sifRead(siffile, &sif);
  if(ret) {imgSetStatus(img, STATUS_NOSIFDATA); return(21);}

  /* Copy SIF contents */
  ret=sif2img(&sif, img, 1, 1, 1, IMG_TEST-2);
  sifEmpty(&sif);
  if(ret!=0) {imgSetStatus(img, STATUS_WRONGSIFDATA); return(22);}
  
  return(0);
}
/*****************************************************************************/

/*****************************************************************************/
/**  Write Analyze 7.5 image.
 *
 *   Analyze database name must be given with path. Path must exist.
 *   Image and header files with .img and .hdr extensions are created.
 *   Existing files are overwritten.
 *   anaFlipping() determines whether image is flipped in z-direction;
 *   image is always flipped in x,y-directions.
 *   Byte order is determined based on _fileFormat field.
 * 
 * @param dbname analyze database name with path, without extension
 * @param img pointer to IMG data
 * @return 0 if ok, 1 invalid input, 2 invalid image status (image not occupied), 
 * 3 failed to resolve extreme values (min and max),
 * 12 failed to allocate temp memory, 14 failed to open file for writing,
 * 15 failed to write data, 21 failed to write header, and
 * sets IMG->statmsg in case of error
*/
int imgWriteAnalyze(const char *dbname, IMG *img) {
  FILE *fp;
  int ret, fi, pi, xi, yi, little;
  float g;
  ANALYZE_DSR dsr;
  char datfile[FILENAME_MAX], hdrfile[FILENAME_MAX], siffile[FILENAME_MAX];
  const char *cptr;
  int pxlNr=0;
  struct tm tm; //*st;
  short int *sdata, *sptr, smin, smax;
  SIF sif;


  if(IMG_TEST) printf("imgWriteAnalyze(%s, *img)\n", dbname);

  /* Check the arguments */
  imgSetStatus(img, STATUS_OK);
  if(img==NULL || img->status!=IMG_STATUS_OCCUPIED) {
    imgSetStatus(img, STATUS_FAULT); return(2);}
  if(dbname==NULL || !dbname[0]) {imgSetStatus(img, STATUS_FAULT); return(1);}
  
  /* Make the image and header filenames */  
  strcpy(datfile, dbname); strcat(datfile, ".img");
  strcpy(hdrfile, dbname); strcat(hdrfile, ".hdr");
  strcpy(siffile, dbname); strcat(siffile, ".sif");


  /*
   *  Fill Analyze header
   */
  if(img->_fileFormat==IMG_ANA_L) dsr.little=1; else dsr.little=0;
  /* Header key */
  memset(&dsr.hk, 0, sizeof(ANALYZE_HEADER_KEY));
  memset(&dsr.dime, 0, sizeof(ANALYZE_HEADER_IMGDIM));
  memset(&dsr.hist, 0, sizeof(ANALYZE_HEADER_HISTORY));
  dsr.hk.sizeof_hdr=348;
  strcpy(dsr.hk.data_type, "");
  cptr=strrchr(dbname, '/'); if(cptr==NULL) cptr=strrchr(dbname, '\\');
  if(cptr!=NULL) cptr++;
  if(cptr==NULL) cptr=dbname;
  strncpy(dsr.hk.db_name, cptr, 17);
  dsr.hk.extents=16384;
  dsr.hk.regular='r';
  /* Image dimension */
  dsr.dime.dim[0]=4;
  dsr.dime.dim[1]=img->dimx;
  dsr.dime.dim[2]=img->dimy;
  dsr.dime.dim[3]=img->dimz;
  dsr.dime.dim[4]=img->dimt;
  dsr.dime.datatype=ANALYZE_DT_SIGNED_SHORT;
  dsr.dime.bitpix=16;
  dsr.dime.pixdim[0]=0.0;
  dsr.dime.pixdim[1]=img->sizex;
  dsr.dime.pixdim[2]=img->sizey;
  dsr.dime.pixdim[3]=img->sizez;
  dsr.dime.pixdim[4]=0.0;
  dsr.dime.funused1=0.0; /* Scale factor is set later */
  /* dsr.dime.funused2=img->zoom; */ /* Reconstruction zoom */
  dsr.dime.funused3=img->isotopeHalflife;
  /* Data history */
  if(img->decayCorrection==IMG_DC_CORRECTED)
    strcpy(dsr.hist.descrip, "Decay corrected.");
  else if(img->decayCorrection==IMG_DC_NONCORRECTED) 
    strcpy(dsr.hist.descrip, "No decay correction.");
  else
    strcpy(dsr.hist.descrip, "");
  if(strlen(img->studyNr)>0 && strcmp(img->studyNr, ".")!=0)
    memcpy(dsr.hist.scannum, img->studyNr, 10);
  else
    strcpy(dsr.hist.scannum, "");
  gmtime_r((time_t*)&img->scanStart, &tm);
  if(!strftime(dsr.hist.exp_date, 10, "%Y-%m-%d", &tm))
    memcpy(dsr.hist.exp_date, "1900-01-01", 10);
  if(!strftime(dsr.hist.exp_time, 10, "%H:%M:%S", &tm))
    strlcpy(dsr.hist.exp_time, "00:00:00", 10);

  /*
   *  Scale data to short int range
   *  Determine and set scale factor and cal_min & cal_max
   */
  if(IMG_TEST) printf("scaling data to short ints\n");
  ret=imgMinMax(img, &dsr.dime.cal_min, &dsr.dime.cal_max);
  if(ret) {imgSetStatus(img, STATUS_FAULT); return(3);}
  if(IMG_TEST) printf("min=%g max=%g\n", dsr.dime.cal_min, dsr.dime.cal_max);
  if(fabs(dsr.dime.cal_min)>fabs(dsr.dime.cal_max)) g=fabs(dsr.dime.cal_min);
  else g=fabs(dsr.dime.cal_max);
  if(isnan(g)) {imgSetStatus(img, STATUS_FAULT); return 3;}
  if(g<1E-20) g=1.0; else g=32767./g; dsr.dime.funused1=1.0/g;
  if(IMG_TEST) printf("scale_factor=%g\n", dsr.dime.funused1);

  /* Allocate memory for short int array */
  pxlNr=(img->dimx)*(img->dimy)*(img->dimz);
  sdata=malloc(pxlNr*sizeof(short int));
  if(sdata==NULL) {imgSetStatus(img, STATUS_NOMEMORY); return 12;}

  /* Open image data file for write */
  if((fp=fopen(datfile, "wb")) == NULL) {
    imgSetStatus(img, STATUS_CANTWRITEIMGFILE);
    free(sdata);
    return 14;
  }

  /* Copy and write image matrix data to short int array */
  /* Data is written one frame at a time */
  smin=smax=temp_roundf(g*img->m[0][0][0][0]);
  for(fi=0; fi<img->dimt; fi++) {
    sptr=sdata;
    if(anaFlipping()==0) {
      for(pi=0; pi<img->dimz; pi++)
        for(yi=img->dimy-1; yi>=0; yi--)
          for(xi=img->dimx-1; xi>=0; xi--) {
            *sptr=temp_roundf(g*img->m[pi][yi][xi][fi]);
            if(*sptr>smax) smax=*sptr; else if(*sptr<smin) smin=*sptr;
            sptr++;
          }
    } else {
      for(pi=img->dimz-1; pi>=0; pi--)
        for(yi=img->dimy-1; yi>=0; yi--)
          for(xi=img->dimx-1; xi>=0; xi--) {
            *sptr=temp_roundf(g*img->m[pi][yi][xi][fi]);
            if(*sptr>smax) smax=*sptr; else if(*sptr<smin) smin=*sptr;
            sptr++;
          }
    }
    /* Change byte order if necessary */
    little=little_endian();
    if(little!=dsr.little)
      swabip(sdata, pxlNr*sizeof(short int));
    /* Write image data */
    if(fwrite(sdata, 2, pxlNr, fp) != (unsigned int)pxlNr) {
      imgSetStatus(img, STATUS_CANTWRITEIMGFILE);
      free(sdata); fclose(fp);
      return 15;
    }
  }
  /* Done writing */
  fclose(fp);
  free(sdata);

  if(IMG_TEST) printf("smin=%d smax=%d\n", smin, smax);

  /* Set header glmin & glmax */
  dsr.dime.glmin=smin; dsr.dime.glmax=smax;
  
  /* Write Analyze header */
  ret=anaWriteHeader(hdrfile, &dsr);
  if(ret) {
    imgSetStatus(img, STATUS_CANTWRITEHEADERFILE);
    return 21;
  }
  imgSetStatus(img, STATUS_OK);

  /* Otherwise ready, but check if SIF should/can be written */
  sifInit(&sif);
  /* Try to read existing SIF */
  ret=sifRead(siffile, &sif);
  if(ret==0) { // SIF could be read
    if(sif.frameNr==img->dimt) {
      /* If size matches, then update the contents, but keep counts, in case
         previous SIF comes with actual count info from scanner */
      ret=img2sif(img, &sif, 1, 1, 0, IMG_TEST-2);
    } else {
      /* otherwise create SIF contents */
      ret=img2sif(img, &sif, 1, 1, 2, IMG_TEST-2);
    }
  } else {
    /* otherwise create SIF contents */
    ret=img2sif(img, &sif, 1, 1, 2, IMG_TEST-2);
  }
  if(ret!=0) {
    /* SIF data could not be made: do not give error, just do not write it */
    if(IMG_TEST>0) printf("SIF contents could not be filled.\n");
    return 0;
  }
  /* Write SIF */
  ret=sifWrite(&sif, siffile);  
  if(ret!=0) {
    /* SIF could not be written: do not give error, just do not write it */
    if(IMG_TEST>0)
      fprintf(stderr, "Error: SIF could not be written (%d).\n", ret);
  }

  imgSetStatus(img, STATUS_OK);
  return 0;
}
/*****************************************************************************/

/*****************************************************************************/
/*!
 * Fill IMG struct header information from Analyze 7.5 database files.
 *
 *  SIF file is read if available. Information concerning separate frames
 *  or planes is not filled though.
 *
 * @param dbname name of Analyze database, may contain filename extension
 * @param img pointer to the initiated IMG data
 * @return errstatus, which is STATUS_OK (0) when call was successful,
 * and >0 in case of an error.
 */
int imgReadAnalyzeHeader(const char *dbname, IMG *img) {
  char hdrfile[FILENAME_MAX], siffile[FILENAME_MAX];
  ANALYZE_DSR ana_header;
  SIF sif;
  double f;
  int ret;

  if(IMG_TEST) printf("\nimgReadAnalyzeHeader(%s, *img)\n", dbname);
  
  /* Check the input */
  if(img==NULL) return STATUS_FAULT;
  if(img->status!=IMG_STATUS_INITIALIZED) return STATUS_FAULT;
  imgSetStatus(img, STATUS_FAULT);
  if(dbname==NULL) return STATUS_FAULT;

  /* Determine the names of hdr and sif files */
  ret=anaDatabaseExists(dbname, hdrfile, NULL, siffile);
  if(ret==0) return STATUS_NOFILE;
  
  /* Read Analyze header file */
  ret=anaReadHeader(hdrfile, &ana_header);
  if(ret!=0) {
    if(IMG_TEST>1) printf("anaReadHeader() return value := %d\n", ret);
    if(ret==1) return STATUS_FAULT;
    else if(ret==2) return STATUS_NOHEADERFILE;
    else return STATUS_UNSUPPORTED;
    return(STATUS_FAULT);
  }
  /* and set IMG contents */
  ret=imgGetAnalyzeHeader(img, &ana_header);
  if(ret!=0) {
    imgSetStatus(img, ret);
    return(ret);
  }

  /* If SIF does not exist, then that's it */
  if(!siffile[0]) {
    imgSetStatus(img, STATUS_OK);
    return STATUS_OK;
  }

  /* SIF is available, so read that too */
  sifInit(&sif); ret=0;
  if(sifRead(siffile, &sif)!=0) return STATUS_OK;
  /* Copy scan time */
  img->scanStart=sif.scantime;
  /* Study number, if not yet defined */
  if(!img->studyNr[0] && strlen(sif.studynr)>1 )
    strlcpy(img->studyNr, sif.studynr, MAX_STUDYNR_LEN+1);
  /* Isotope half-life, if not yet defined */
  f=hlFromIsotope(sif.isotope_name);
  if(img->isotopeHalflife<=0.0 && f>0.0) img->isotopeHalflife=60.0*f;
  sifEmpty(&sif);

  return STATUS_OK;
}
/*****************************************************************************/

/*****************************************************************************/
/*!
 * Copy Analyze 7.5 header information into IMG.
 * 
 * @param img image structure
 * @param h Analyze header structure
 * @return errstatus, which is STATUS_OK (0) when call was successful,
 * and >0 in case of an error.
 */
int imgGetAnalyzeHeader(IMG *img, ANALYZE_DSR *h) {
  int dimNr, dimx, dimy, dimz=1, dimt=1, pxlNr=0;

  if(IMG_TEST) printf("\nimgGetAnalyzeHeader(*img, *dsr)\n");
  
  /* Check the input */
  
  if(img==NULL) return STATUS_FAULT;
  if(img->status!=IMG_STATUS_INITIALIZED && img->status!=IMG_STATUS_OCCUPIED)
    return STATUS_FAULT;
  imgSetStatus(img, STATUS_FAULT);
  if(h==NULL) return STATUS_FAULT;
    
  imgSetStatus(img, STATUS_INVALIDHEADER);

  /* Get the image dimensions from header */
  dimNr=h->dime.dim[0];
  if(dimNr<2) return STATUS_INVALIDHEADER;
  dimx=h->dime.dim[1]; dimy=h->dime.dim[2];
  if(dimNr>2) {dimz=h->dime.dim[3]; if(dimNr>3) dimt=h->dime.dim[4];}
  pxlNr=dimx*dimy*dimz;
  if(pxlNr<1) return STATUS_INVALIDHEADER;
  img->dimx=dimx; img->dimy=dimy; img->dimz=dimz; img->dimt=dimt;

  /* Copy information from Analyze header */
  img->type=IMG_TYPE_IMAGE;
  strlcpy(img->studyNr, h->hist.patient_id, MAX_STUDYNR_LEN+1);
  if(strcmp(img->studyNr, ".")==0) strcpy(img->studyNr, "");
  strcpy(img->patientName, h->hist.patient_id);
  img->sizex=h->dime.pixdim[1];
  img->sizey=h->dime.pixdim[2];
  img->sizez=h->dime.pixdim[3];
  /*if(h->dime.funused2>1.E-5) img->zoom=h->dime.funused2;*/
  if(h->dime.funused3>1.E-5) img->isotopeHalflife=h->dime.funused3;
  if(h->little) img->_fileFormat=IMG_ANA_L; else img->_fileFormat=IMG_ANA;
  
  /* Decay correction */
  if(strstr(h->hist.descrip, "Decay corrected.")!=NULL)
    img->decayCorrection=IMG_DC_CORRECTED;
  else if(strstr(h->hist.descrip, "No decay correction.")!=NULL)
    img->decayCorrection=IMG_DC_NONCORRECTED;
  else
    img->decayCorrection=IMG_DC_CORRECTED; // just assumed so

  imgSetStatus(img, STATUS_OK);
  return STATUS_OK;
}
/*****************************************************************************/

/*****************************************************************************/
/**
 * Copy header information in IMG struct into Analyze 7.5 header struct.
 *
 *  Min, max, and scale factor are set here and they apply to all frames.
 *
 * @return errstatus, which is STATUS_OK (0) when call was successful,
 * and >0 in case of an error.
 */
int imgSetAnalyzeHeader(
  /** pointer to IMG struct from which header information is read */
  IMG *img,
  /** Analyze 7.5 database name */
  const char *dbname,
  /** pointer to Analyze header struct to be filled */
  ANALYZE_DSR *dsr,
  /** minimum pixel value in all frames that will be written */
  float fmin,
  /** maximum pixel value in all frames that will be written */
  float fmax
) {
  struct tm tm; //*st;
  char *cptr;
  float g;

  if(IMG_TEST) printf("\nimgSetAnalyzeHeader(*img, %s, *dsr, %g, %g)\n",
                      dbname, fmin, fmax);
  
  /* Check the input */
  if(img==NULL) return STATUS_FAULT;
  if(img->status!=IMG_STATUS_INITIALIZED && img->status!=IMG_STATUS_OCCUPIED)
    return STATUS_FAULT;
  imgSetStatus(img, STATUS_FAULT);
  if(dsr==NULL) return STATUS_FAULT;

  /* Byte order */
  if(img->_fileFormat==IMG_ANA_L) dsr->little=1; else dsr->little=0;
  /* Header key */
  memset(&dsr->hk, 0, sizeof(ANALYZE_HEADER_KEY));
  memset(&dsr->dime, 0, sizeof(ANALYZE_HEADER_IMGDIM));
  memset(&dsr->hist, 0, sizeof(ANALYZE_HEADER_HISTORY));
  dsr->hk.sizeof_hdr=348;
  strcpy(dsr->hk.data_type, "");
  cptr=strrchr(dbname, '/'); if(cptr==NULL) cptr=strrchr(dbname, '\\');
  if(cptr!=NULL) cptr++;
  if(cptr==NULL) cptr=(char*)dbname;
  strncpy(dsr->hk.db_name, cptr, 17);
  dsr->hk.extents=16384;
  dsr->hk.regular='r';
  /* Image dimension */
  dsr->dime.dim[0]=4;
  dsr->dime.dim[1]=img->dimx;
  dsr->dime.dim[2]=img->dimy;
  dsr->dime.dim[3]=img->dimz;
  dsr->dime.dim[4]=img->dimt;
  dsr->dime.datatype=ANALYZE_DT_SIGNED_SHORT;
  dsr->dime.bitpix=16;
  dsr->dime.pixdim[0]=0.0;
  dsr->dime.pixdim[1]=img->sizex;
  dsr->dime.pixdim[2]=img->sizey;
  dsr->dime.pixdim[3]=img->sizez;
  dsr->dime.pixdim[4]=0.0;
  dsr->dime.funused1=0.0; /* Scale factor is set later */
  /* dsr.dime.funused2=img->zoom; */ /* Reconstruction zoom */
  dsr->dime.funused3=img->isotopeHalflife;
  /* Data history */
  if(img->decayCorrection==IMG_DC_CORRECTED)
    strcpy(dsr->hist.descrip, "Decay corrected.");
  else if(img->decayCorrection==IMG_DC_NONCORRECTED)
    strcpy(dsr->hist.descrip, "No decay correction.");
  else
    strcpy(dsr->hist.descrip, "");
  if(strlen(img->studyNr)>0 && strcmp(img->studyNr, ".")!=0)
    memcpy(dsr->hist.scannum, img->studyNr, 10);
  else
    strcpy(img->studyNr, "");
  gmtime_r((time_t*)&img->scanStart, &tm);
  if(!strftime(dsr->hist.exp_date, 10, "%Y%m%d", &tm))
    strcpy(dsr->hist.exp_date, "19000101");
  if(!strftime(dsr->hist.exp_time, 10, "%H:%M:%S", &tm))
    strncpy(dsr->hist.exp_time, "00:00:00", 10);

  /* Determine and set scale factor and cal_min & cal_max */
  if(fmin<fmax) {
    dsr->dime.cal_min=fmin; dsr->dime.cal_max=fmax;
  } else { /* not given in function call, try to find those here */
    if(img->status==IMG_STATUS_OCCUPIED &&
       imgMinMax(img, &dsr->dime.cal_min, &dsr->dime.cal_max)==0) {}
    else return STATUS_FAULT;
  }
  if(fabs(dsr->dime.cal_min) > fabs(dsr->dime.cal_max)) g=fabs(dsr->dime.cal_min);
  else g = fabs(dsr->dime.cal_max);
  /* if(fabs(dsr->dime.cal_min)>fabs(dsr->dime.cal_max))
     g=fabs(dsr->dime.cal_min); */
  /* else g=fabs(dsr->dime.cal_max); */
  if(g<1E-20) g=1.0; else g=32767./g; dsr->dime.funused1=1.0/g;
  /* Set header glmin & glmax */
  dsr->dime.glmin=temp_roundf(fmin*g); dsr->dime.glmax=temp_roundf(fmax*g);
  /* printf("glmin=%d\n", dsr->dime.glmin); */
  /* printf("glmax=%d\n", dsr->dime.glmax); */

  imgSetStatus(img, STATUS_OK);
  return STATUS_OK;
}
/*****************************************************************************/

/*****************************************************************************/
/*!
 * Read the first frame from an Analyze 7.5 database into IMG data structure.
 *
 * @param fname Name of Analyze database from which IMG contents will be read
 * @param img pointer to the initiated but not preallocated IMG data
 * @return errstatus, which is STATUS_OK (0) when call was successful,
 * and >0 in case of an error.
 */
int imgReadAnalyzeFirstFrame(const char *fname, IMG *img) {
  int ret=0;

  if(IMG_TEST) printf("\nimgReadAnalyzeFirstFrame(%s, *img)\n", fname);
  /* Check the input */
  if(img==NULL) return STATUS_FAULT;
  if(img->status!=IMG_STATUS_INITIALIZED) return STATUS_FAULT;
  imgSetStatus(img, STATUS_FAULT);
  if(fname==NULL) return STATUS_FAULT;

  /* Read header information from file */
  ret=imgReadAnalyzeHeader(fname, img); if(ret) return(ret);
  if(IMG_TEST>3) imgInfo(img);

  /* Allocate memory for one frame */
  img->dimt=1;
  ret=imgAllocate(img, img->dimz, img->dimy, img->dimx, img->dimt);
  if(ret) return STATUS_NOMEMORY;

  /* Read the first frame */
  ret=imgReadAnalyzeFrame(fname, 1, img, 0); if(ret) return(ret); 

  /* All went well */
  imgSetStatus(img, STATUS_OK);
  return STATUS_OK;
}
/*****************************************************************************/

/*****************************************************************************/
/*!
 *  Read a specified frame from an Analyze 7.5 database into preallocated IMG
 *  data structure. 
 *
 *  Analyze database consists of two or three files in the same 
 *  directory: fname.hdr, fname.img, and optionally fname.sif.
 *  IMG header is assumed to be filled correctly before calling this function,
 *  except for information concerning separate planes and this frame,
 *  which is filled here.
 *  If frame does not exist, then and only then STATUS_NOMATRIX is returned.
 *
 * @param fname name of Analyze database from which IMG contents will be read
 * @param frame_to_read frame which will be read [1..frameNr]
 * @param img pointer to the IMG data. Place for the frame must be preallocated
 * @param frame_index IMG frame index [0..dimt-1] where data will be placed
 * @return errstatus, which is STATUS_OK (0) when call was successful,
 * and >0 in case of an error.
 */ 
int imgReadAnalyzeFrame(
  const char *fname, int frame_to_read, IMG *img, int frame_index
) {
  FILE *fp;
  int ret, pi, xi, yi;
  float *fdata=NULL, *fptr;
  ANALYZE_DSR dsr;
  char datfile[FILENAME_MAX], hdrfile[FILENAME_MAX], siffile[FILENAME_MAX];
  SIF sif;


  if(IMG_TEST) printf("\nimgReadAnalyzeFrame(%s, %d, *img, %d)\n",
    fname, frame_to_read, frame_index);
    
  /* Check the input */
  if(img==NULL) return STATUS_FAULT;
  if(img->status!=IMG_STATUS_OCCUPIED) return STATUS_FAULT;
  if(fname==NULL) return STATUS_FAULT;
  if(frame_index<0 || frame_index>img->dimt-1) return STATUS_FAULT;
  if(frame_to_read<1) return STATUS_FAULT;
  imgSetStatus(img, STATUS_FAULT);
  
  /* Determine the names of hdr, data and sif files */
  ret=anaDatabaseExists(fname, hdrfile, datfile, siffile);
  if(ret==0) return STATUS_NOFILE;

  /* Read Analyze header file */
  ret=anaReadHeader(hdrfile, &dsr);
  if(ret!=0) {
    if(ret==1) return STATUS_FAULT;
    else if(ret==2) return STATUS_NOHEADERFILE;
    else return STATUS_UNSUPPORTED;
    return(STATUS_FAULT);
  }

  /* Open image datafile */
  if(IMG_TEST>2) fprintf(stdout, "reading image data %s\n", datfile);
  if((fp=fopen(datfile, "rb")) == NULL) return STATUS_NOIMGDATA;

  /* Allocate memory for one image frame */
  fdata=malloc(img->dimx*img->dimy*img->dimz*sizeof(float));
  if(fdata==NULL) {fclose(fp); return STATUS_NOMEMORY;}

  /* Read the required image frame */
  fptr=fdata;
  ret=anaReadImagedata(fp, &dsr, frame_to_read, fptr);
  fclose(fp);
  if(ret==3) {free(fdata); return STATUS_NOMATRIX;} /* no more frames */
  if(ret!=0) {free(fdata); return STATUS_UNSUPPORTED;}

  /* Copy pixel values to IMG */
  fptr=fdata;
  if(anaFlipping()==0) { /* no flipping in z-direction */
    for(pi=0; pi<img->dimz; pi++)
      for(yi=img->dimy-1; yi>=0; yi--)
        for(xi=img->dimx-1; xi>=0; xi--)
          img->m[pi][yi][xi][frame_index]=*fptr++;
  } else {
    for(pi=img->dimz-1; pi>=0; pi--)
      for(yi=img->dimy-1; yi>=0; yi--)
        for(xi=img->dimx-1; xi>=0; xi--)
          img->m[pi][yi][xi][frame_index]=*fptr++;
  }
  free(fdata);

  /* Set decay correction factor to zero */
  img->decayCorrFactor[frame_index]=0.0;

  imgSetStatus(img, STATUS_OK); /* If the rest is failed, no problem */

  /* 
   *  Try to read frame time information from SIF file
   */
  sifInit(&sif);
  if(sifRead(siffile, &sif)!=0) return STATUS_OK;
  /* Frame information */
  if(sif.frameNr>=frame_to_read) {
    img->start[frame_index]=sif.x1[frame_to_read-1];
    img->end[frame_index]=sif.x2[frame_to_read-1];
    img->mid[frame_index]=0.5*(img->start[frame_index]+img->end[frame_index]);
    img->prompts[frame_index]=sif.prompts[frame_to_read-1];
    img->randoms[frame_index]=sif.randoms[frame_to_read-1];
  }
  sifEmpty(&sif);

  return STATUS_OK;
}
/*****************************************************************************/

/*****************************************************************************/
/*!
 * Write one PET frame from IMG data struct into Analyze 7.5 database file.
 *  This function can be called repeatedly to write all frames one at a time
 *  to conserve memory. This function does not write SIF.
 *
 * @param dbname name of file where IMG contents will be written.
 *  If file does not exist, it is created.
 *  Make sure to delete existing file, unless you want to add data
 * @param frame_to_write PET frame number (1..frameNr) which will be written:
 *  If set to 0, frame data will be written to an existing or new PET file as
 *  a new frame, never overwriting existing data.
 *  If >0, then frame data is written as specified frame number, overwriting
 *  any data existing with the same frame number.
 * @param img pointer to the IMG data struct
 * @param frame_index IMG frame index (0..dimt-1) which will be written
 * @param fmin minimum pixel value in all frames that will be written;
 *  used only when writing the first frame
 * @param fmax maximum pixel value in all frames that will be written;
 *  used only when writing the first frame
 * @return errstatus, which is STATUS_OK (0) when call was successful,
 * and >0 in case of an error.
 */
int imgWriteAnalyzeFrame(
  const char *dbname, int frame_to_write, IMG *img, int frame_index,
  float fmin, float fmax
) {
  IMG test_img;
  int ret=0, pxlNr, zi, xi, yi, little;
  FILE *fp;
  short int *sdata=NULL, *sptr;
  char datfile[FILENAME_MAX], hdrfile[FILENAME_MAX], siffile[FILENAME_MAX];
  ANALYZE_DSR dsr;
  float scale_factor=1.0;


  if(IMG_TEST) printf("\nimgWriteAnalyzeFrame(%s, %d, *img, %d, %g, %g)\n",
    dbname, frame_to_write, frame_index, fmin, fmax);

  /*
   *  Check the input 
   */
  if(dbname==NULL) return STATUS_FAULT;
  if(img==NULL) return STATUS_FAULT;
  if(img->status!=IMG_STATUS_OCCUPIED) return STATUS_FAULT;
  if(frame_to_write<0) return STATUS_FAULT;
  if(frame_index<0 || frame_index>=img->dimt) return STATUS_FAULT;
  if(img->_fileFormat!=IMG_ANA_L && img->_fileFormat!=IMG_ANA)
    return STATUS_FAULT;

  /*
   *  If database does not exist, then create it with new header,
   *  and if it does exist, then read and check header information.
   *  Create or edit header to contain correct frame nr.
   *  Determine the global scaling factor.   
   */
  imgInit(&test_img);
  if(anaDatabaseExists(dbname, hdrfile, datfile, siffile)==0) { // not existing

    /* Create database filenames */
    sprintf(hdrfile, "%s.hdr", dbname);
    sprintf(datfile, "%s.img", dbname);
    sprintf(siffile, "%s.sif", dbname);

    /* Set main header */
    imgSetAnalyzeHeader(img, dbname, &dsr, fmin, fmax);
    if(frame_to_write==0) frame_to_write=1;
    dsr.dime.dim[4]=frame_to_write;
    scale_factor=dsr.dime.funused1;
    if(fabs(scale_factor)>1.0E-20) scale_factor=1.0/scale_factor;

    /* Write Analyze header */
    ret=anaWriteHeader(hdrfile, &dsr);
    if(ret && IMG_TEST) printf("anaWriteHeader() := %d\n", ret);
    if(ret) return STATUS_CANTWRITEHEADERFILE;

    /* Remove datafile if necessary */
    if(access(datfile, 0) != -1) remove(datfile);

  } else { /* database does exist */
  
    /* Read header information for checking */
    ret=imgReadAnalyzeHeader(dbname, &test_img);
    if(ret!=0) {imgEmpty(&test_img); return ret;}
    /* Check that file format is the same */
    if(img->_fileFormat!=test_img._fileFormat || img->type!=test_img.type) {
      imgEmpty(&test_img); return STATUS_WRONGFILETYPE;}
    /* Check that matrix sizes are the same */
    if(img->dimz!=test_img.dimz || img->dimx!=test_img.dimx ||
       img->dimy!=test_img.dimy) {
      imgEmpty(&test_img); return STATUS_VARMATSIZE;}
    imgEmpty(&test_img);

    /* Read the header, set new frame number, and write it back */
    /* Get also the scale factor */
    if((ret=anaReadHeader(hdrfile, &dsr))!=0) return STATUS_NOMAINHEADER;
    scale_factor=1.0/dsr.dime.funused1;
    if(frame_to_write==0) frame_to_write=dsr.dime.dim[4]+1;
    if(dsr.dime.dim[4]<frame_to_write) {
      if(dsr.dime.dim[4]+1<frame_to_write) return STATUS_MISSINGMATRIX;
      dsr.dime.dim[4]=frame_to_write;
    }
    if((ret=anaWriteHeader(hdrfile, &dsr))!=0) return STATUS_NOWRITEPERM;
  }
  if(IMG_TEST>2) {
    printf("frame_to_write := %d\n", frame_to_write);
    printf("hdrfile := %s\n", hdrfile);
    printf("datfile := %s\n", datfile);
    printf("siffile := %s\n", siffile);
  }

  /* Allocate memory for matrix short int data (one plane) */
  pxlNr=img->dimx*img->dimy;
  sdata=(short int*)malloc(pxlNr*sizeof(short int));
  if(sdata==NULL) return STATUS_NOMEMORY;

  /* Open datafile, not removing possible old contents */
  if(frame_to_write==1) fp=fopen(datfile, "wb"); else fp=fopen(datfile, "r+b");
  if(fp==NULL) {free(sdata); return STATUS_CANTWRITEIMGFILE;}
  /* Move file pointer to the place of current frame */
  if(fseek(fp, (frame_to_write-1)*pxlNr*img->dimz*sizeof(short int),
           SEEK_SET)!=0) {
    free(sdata); fclose(fp); return STATUS_MISSINGMATRIX;}
  little=little_endian();
  /* Copy, scale and write data plane-by-plane */
  if(anaFlipping()==0) {
    for(zi=0; zi<img->dimz; zi++) {
      sptr=sdata;
     /*printf("plane := %d\n  scale_factor := %g\n", zi+1, scale_factor);*/
      for(yi=img->dimy-1; yi>=0; yi--) for(xi=img->dimx-1; xi>=0; xi--) {
        *sptr=temp_roundf(scale_factor*img->m[zi][yi][xi][frame_index]); sptr++;
      }
      /* Change byte order if necessary */
      sptr=sdata; if(little!=dsr.little) swabip(sptr, pxlNr*sizeof(short int));
      /* Write image data */
      sptr=sdata;
      if(fwrite(sptr, sizeof(short int), pxlNr, fp) != (unsigned int)pxlNr) {
        free(sdata); fclose(fp); return STATUS_CANTWRITEIMGFILE;
      }
    }
  } else {
    for(zi=img->dimz-1; zi>=0; zi--) {
      sptr=sdata;
      for(yi=img->dimy-1; yi>=0; yi--) for(xi=img->dimx-1; xi>=0; xi--) {
        *sptr=temp_roundf(scale_factor*img->m[zi][yi][xi][frame_index]); sptr++;
      }
      /* Change byte order if necessary */
      sptr=sdata; if(little!=dsr.little) swabip(sptr, pxlNr*sizeof(short int));
      /* Write image data */
      sptr=sdata;
      if(fwrite(sptr, sizeof(short int), pxlNr, fp) != (unsigned int)pxlNr) {
        free(sdata); fclose(fp); return STATUS_CANTWRITEIMGFILE;
      }
    }
  }
  free(sdata);
  fclose(fp);

  return STATUS_OK;
}
/*****************************************************************************/

/*****************************************************************************/
