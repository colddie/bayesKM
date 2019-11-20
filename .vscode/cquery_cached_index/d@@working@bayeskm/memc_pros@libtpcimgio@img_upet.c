/// @file img_upet.c
/// @author Vesa Oikonen
/// @brief I/O routines for IMG data from/to Siemens Inveon format.
///
/*****************************************************************************/
#include "libtpcimgio.h"
/*****************************************************************************/

/*****************************************************************************/
/** Read MicroPET image and write ECAT 7 image volume frame-by-frame.
 
    @return Returns errstatus, which is STATUS_OK (0) when call was successful,
    and >0 in case of an error.  
*/
int imgMicropetToEcat7(
  /** MicroPET image filename. */
  char *upetname,
  /** ECAT image filename. */
  char *ecatfile,
  /** Verbose level */
  int verbose
) {
  char upetheader[FILENAME_MAX], upetimage[FILENAME_MAX];
  int n, ret;
  int acquisition_mode, data_type;
  char tmp[MAX_MICROPET_LINE_LEN];


  if(verbose>1) printf("\nimgMicropetToEcat7(%s, %s, %d)\n",
    upetname, ecatfile, verbose);
  /* Check the arguments */
  if(upetname==NULL || ecatfile==NULL) return STATUS_FAULT;
  ret=upetExists(upetname, upetheader, upetimage, verbose-1);
  if(ret!=2) return STATUS_NOFILE;

  /*
   *  Open Micropet Header and binary data files
   */
  FILE *fph, *fpi;
  if((fph=fopen(upetheader, "r"))==NULL) return(STATUS_NOHEADERFILE);
  if((fpi=fopen(upetimage, "rb"))==NULL) {fclose(fph); return(STATUS_NOIMGDATA);}


  /*
   *  Check that image format is (currently) supported
   */
  rewind(fph);
  if(verbose>1) printf("checking that image format is supported\n");
  n=-1; if(upetHeaderReadParameter(fph, "file_type", tmp)==0)
    (void)sscanf(tmp, "%d", &n);
  if(verbose>2) printf("file_type := %d\n", n);
  if(n!=5) {fclose(fph); fclose(fpi); return(STATUS_UNSUPPORTED);}
  acquisition_mode=-1;
  if(upetHeaderReadParameter(fph, "acquisition_mode", tmp)==0)
    (void)sscanf(tmp, "%d", &acquisition_mode);
  if(verbose>2) printf("acquisition_mode := %d\n", acquisition_mode);
  if(acquisition_mode!=2 && acquisition_mode!=3 && acquisition_mode!=9) {
    fclose(fph); fclose(fpi); return(STATUS_UNSUPPORTED);}
  data_type=-1;
  if(upetHeaderReadParameter(fph, "data_type", tmp)==0)
    (void)sscanf(tmp, "%d", &data_type);
  if(verbose>2) printf("data_type := %d\n", data_type);
  if(data_type!=4 && data_type!=2) {
    fclose(fph); fclose(fpi); return(STATUS_UNSUPPORTED);}

  /*
   *  Convert PET or CT image
   */
  if(acquisition_mode==2 || acquisition_mode==3) {
    ret=imgMicropetPETToEcat7(fph, fpi, ecatfile, verbose);
  } else if(acquisition_mode==9) {
    ret=imgMicropetCTToEcat7(fph, fpi, ecatfile, verbose);
  } else {
    ret=STATUS_UNSUPPORTED;
  }
  fclose(fph); fclose(fpi);
  return ret;
}
/*****************************************************************************/

/*****************************************************************************/
/** Read MicroPET static or dynamic PET image and write ECAT 7 image volume 
    frame-by-frame.
    @return Returns errstatus, which is STATUS_OK (0) when call was successful,
    and >0 in case of an error.  
*/
int imgMicropetPETToEcat7(
  /** MicroPET header file pointer. */
  FILE *fph,
  /** MicroPET image datafile pointer. */
  FILE *fpi,
  /** ECAT image filename. */
  char *ecatfile,
  /** Verbose level. */
  int verbose
) {
  int n, zi, xi, yi, ti, ret;


  if(verbose>1) printf("imgMicropetPETToEcat7(*fph, *fpi, %s, %d)\n",
                             ecatfile, verbose);
  /* Check input */
  if(fph==NULL || fpi==NULL || ecatfile==NULL) return STATUS_FAULT;

  /* Remove existing ECAT file */
  if(access(ecatfile, 0)!=-1 && remove(ecatfile)!=0) {
    return(STATUS_CANNOTERASE);
  }

  /*
   *  Read image dimensions from header
   */
  int zdim, xdim, ydim, tdim;
  ret=upetGetImageDimensions(fph, &zdim, &xdim, &ydim, &tdim);
  if(ret) {return(STATUS_INVALIDHEADER);}
  if(verbose>1) {
    printf("z_dim := %d\n", zdim);
    printf("x_dim := %d\n", xdim);
    printf("y_dim := %d\n", ydim);
    printf("t_dim := %d\n", tdim);
  }

  /*
   *  Read and write image frame-by-frame
   */
  IMG img; imgInit(&img);
  /* Allocate memory for one frame */
  ret=imgAllocate(&img, zdim, ydim, xdim, 1);
  if(ret) {return(STATUS_NOMEMORY);}
  /* Fill header with what we now can */
  float calibration_factor;
  ret=imgGetMicropetMainHeader(fph, &img, &calibration_factor, verbose-2);
  if(ret) {
    if(verbose>2) printf("ret := %d\n", ret);
    imgEmpty(&img); return(STATUS_INVALIDHEADER);
  }
  if(verbose>1) printf("calibration_factor := %g\n", calibration_factor);
  img._fileFormat=IMG_E7;
  img.type=IMG_TYPE_IMAGE;
  studynr_from_fname(ecatfile, img.studyNr);
  upetScanStart(fph, &img.scanStart);
  /* Allocate memory for the binary data */
  int pxlnr=xdim*ydim*zdim;
  char *mdata, *mptr;
  mdata=(char*)malloc(pxlnr*sizeof(float)); if(mdata==NULL) {
    imgEmpty(&img); return(STATUS_NOMEMORY);
  }
  /* Frame-by-frame */
  float *fptr;
  for(ti=0; ti<tdim; ti++) {
    if(verbose>3) {printf("ti=%d\n", ti); fflush(stdout);}
    /* Read frame information from MicroPET header into IMG */
    ret=imgGetMicropetFrameHeader(fph, &img, ti, verbose-2);
    if(ret) {
      if(verbose==0) {fprintf(stdout, "\n"); fflush(stdout);}
      free(mdata); imgEmpty(&img);
      return(STATUS_INVALIDHEADER);
    }
    /* Read floats */
    mptr=mdata;
    if((n=fread(mptr, 4, pxlnr, fpi)) < pxlnr) {
      if(verbose==0) {fprintf(stdout, "\n"); fflush(stdout);}
      free(mdata); imgEmpty(&img);
      return(STATUS_NOMATRIX);
    }
    /* Copy floats to IMG */
    mptr=mdata;
    for(zi=0; zi<zdim; zi++)
      for(yi=0; yi<ydim; yi++)
        for(xi=0; xi<xdim; xi++) {
          fptr=(float*)mptr;
          img.m[zi][yi][xi][0]=(*fptr)*img.weight[0]*calibration_factor;
          mptr+=4;
        }
    /* Write frame */
    ret=imgWriteFrame(ecatfile, ti+1, &img, 0); //printf("ret := %d\n", ret);
    if(ret!=STATUS_OK) break;
    if(verbose>1) {printf("    frame written.\n"); fflush(stdout);}
    else if(verbose==0) {fprintf(stdout, "."); fflush(stdout);}
  };
  free(mdata); imgEmpty(&img);
  if(verbose==0) {fprintf(stdout, "\n"); fflush(stdout);}
  if(verbose==0 && ret==STATUS_NOMATRIX) {
    fprintf(stdout, "  %d frame(s) processed.\n", ti);
  }
  if(ret!=STATUS_OK && ret!=STATUS_NOMATRIX) {
    remove(ecatfile); return ret;
  }

  return STATUS_OK;
}
/*****************************************************************************/

/*****************************************************************************/
/** Read MicroPET CT image and write ECAT 7 image volume.
    @return Returns errstatus, which is STATUS_OK (0) when call was successful,
    and >0 in case of an error.  
*/
int imgMicropetCTToEcat7(
  /** MicroPET header file pointer */
  FILE *fph,
  /** MicroPET image datafile pointer */
  FILE *fpi,
  /** ECAT image filename */
  char *ecatfile,
  /** Verbose level */
  int verbose
) {
  IMG img;
  int n, pxlnr, zi, xi, yi, zdim, xdim, ydim, ret;
  float f, scale_factor;
  char *mdata, *mptr;
  char tmp[MAX_MICROPET_LINE_LEN];
  short int *si;


  if(MICROPET_TEST>1) printf("imgMicropetCTToEcat7(*fph, *fpi, %s, %d)\n",
                             ecatfile, verbose);
  /* Check input */
  if(fph==NULL || fpi==NULL || ecatfile==NULL) return STATUS_FAULT;

  /*
   *  Read image dimensions from header
   */
  ret=upetGetImageDimensions(fph, &zdim, &xdim, &ydim, NULL);
  if(ret) {return(STATUS_INVALIDHEADER);}
  if(MICROPET_TEST>1) {
    printf("z_dim := %d\n", zdim);
    printf("x_dim := %d\n", xdim);
    printf("y_dim := %d\n", ydim);
  }

  /* Read scale factor */
  rewind(fph);
  if(upetHeaderReadParameter(fph, "scale_factor", tmp)!=0) {
    return(STATUS_INVALIDHEADER);}
  scale_factor=-1; (void)sscanf(tmp, "%f", &scale_factor);
  if(scale_factor<=0) return(STATUS_INVALIDHEADER);
  if(MICROPET_TEST>1) {
    printf("scale_factor := %g\n", scale_factor);
  }

  /* Remove existing ECAT file */
  if(access(ecatfile, 0)!=-1 && remove(ecatfile)!=0) {
    return(STATUS_CANNOTERASE);
  }

  /*
   *  Read and write image
   */
  imgInit(&img);
  /* Allocate memory for one frame */
  ret=imgAllocate(&img, zdim, ydim, xdim, 1);
  if(ret) {return(STATUS_NOMEMORY);}
  /* Fill header with what we now can */
  ret=imgGetMicropetMainHeader(fph, &img, NULL, verbose-2);
  if(ret) {
    if(MICROPET_TEST) printf("ret := %d\n", ret);
    imgEmpty(&img); return(STATUS_INVALIDHEADER);
  }
  img._fileFormat=IMG_E7;
  img.type=IMG_TYPE_IMAGE;
  studynr_from_fname(ecatfile, img.studyNr);
  upetScanStart(fph, &img.scanStart);
  /* Allocate memory for the binary data */
  pxlnr=xdim*ydim;
  mdata=(char*)malloc(pxlnr*sizeof(short int)); if(mdata==NULL) {
    imgEmpty(&img); return(STATUS_NOMEMORY);
  }
  /* Read image data, plane-by-plane */
  for(zi=0; zi<zdim; zi++) {
    mptr=mdata;
    if((n=fread(mptr, 2, pxlnr, fpi)) < pxlnr) {
      if(verbose==0) {fprintf(stdout, "\n"); fflush(stdout);}
      free(mdata); imgEmpty(&img);
      return(STATUS_NOMATRIX);
    }
    /* Copy short ints to IMG */
    mptr=mdata;
    for(yi=0; yi<ydim; yi++)
      for(xi=0; xi<xdim; xi++) {
        si=(short int*)mptr;
        f=(float)*si*scale_factor;
        if(f>=0.0) img.m[zi][yi][xi][0]=f; else img.m[zi][yi][xi][0]=0.0;
        mptr+=2;
      }
    if(MICROPET_TEST>1) printf("   plane %d\n", zi+1);
    else if(verbose==0) {fprintf(stdout, "."); fflush(stdout);}
  }
  free(mdata);
  if(verbose==0) {fprintf(stdout, "\n"); fflush(stdout);}
  /* Save ECAT 7 image volume */
  ret=imgWrite(ecatfile, &img);
  if(ret!=0) {imgEmpty(&img); return(STATUS_CANNOTWRITE);}

  imgEmpty(&img);
  return STATUS_OK;
}
/*****************************************************************************/

/*****************************************************************************/
/** Read main header information from MicroPET header into one-frame-IMG.
    @return Returns 0 when successful and >0 in case of an error.  
 */
int imgGetMicropetMainHeader(
  /** MicroPET header file pointer */
  FILE *fp,
  /** Pointer to allocated IMG struct */
  IMG *img,
  /** Calibration factor / Branching fraction */
  float *calibration_factor,
  /** Verbose level */
  int verbose
) {
  char tmp[MAX_MICROPET_LINE_LEN];
  int n;
  float f;


  if(verbose>0) printf("imgGetMicropetMainHeader(*fp, *img, *f)\n");
  if(fp==NULL) return 1;
  if(img==NULL) return 2;

  /* scanner model */
  rewind(fp);
  if(verbose>1) printf("  reading 'model'\n");
  if(upetHeaderReadParameter(fp, "model", tmp)!=0) return 11;
  n=-1; (void)sscanf(tmp, "%d", &n); if(n<0) return 11;
  img->scanner=n;

  /* zoom */
  rewind(fp);
  if(verbose>1) printf("  reading 'zoom'\n");
  if(upetHeaderReadParameter(fp, "zoom", tmp)!=0) return 11;
  f=-1; (void)sscanf(tmp, "%f", &f); if(f<0) return 11;
  img->zoom=f;


  /* pixel size x */
  rewind(fp);
  if(verbose>1) printf("  reading 'pixel_size_x'\n");
  if(upetHeaderReadParameter(fp, "pixel_size_x", tmp)==0) {
    f=-1; (void)sscanf(tmp, "%f", &f); if(f<0) return 12;
  } else {
    rewind(fp);
    if(verbose>1) printf("  reading 'pixel_size'\n");
    if(upetHeaderReadParameter(fp, "pixel_size", tmp)!=0) return 12;
    f=-1; (void)sscanf(tmp, "%f", &f); if(f<0) return 12;
    f*=10.;
  }
  img->sizex=f;

  /* pixel size y */
  rewind(fp);
  if(verbose>1) printf("  reading 'pixel_size_y'\n");
  if(upetHeaderReadParameter(fp, "pixel_size_y", tmp)==0) {
    f=-1; (void)sscanf(tmp, "%f", &f); if(f<0) return 13;
  } else {
    rewind(fp);
    if(verbose>1) printf("  reading 'pixel_size'\n");
    if(upetHeaderReadParameter(fp, "pixel_size", tmp)!=0) return 13;
    f=-1; (void)sscanf(tmp, "%f", &f); if(f<0) return 12;
    f*=10.;
  }
  img->sizey=f;

  /* pixel size z; note that this is replaced below by transaxial_bin_size */
  rewind(fp);
  if(verbose>1) printf("  reading 'pixel_size_z'\n");
  if(upetHeaderReadParameter(fp, "pixel_size_z", tmp)==0) {
    f=-1; (void)sscanf(tmp, "%f", &f); if(f<0) return 14;
  } else {
    rewind(fp);
    if(verbose>1) printf("  reading 'axial_plane_size'\n");
    if(upetHeaderReadParameter(fp, "axial_plane_size", tmp)!=0) {
      if(verbose>0) printf("  cannot find z pixel size\n");
      f=0.0;
    } else {
      f=-1; (void)sscanf(tmp, "%f", &f); if(f<0) return 14;
      f*=10.;
    }
  }
  img->sizez=f;

  /* transaxial_bin_size */
  rewind(fp);
  if(verbose>1) printf("  reading 'transaxial_bin_size'\n");
  if(upetHeaderReadParameter(fp, "transaxial_bin_size", tmp)==0) {
    f=-1; (void)sscanf(tmp, "%f", &f); if(f>0) img->sizez=10.0*f;
  }

  /* isotope halflife */
  rewind(fp);
  if(verbose>1) printf("  reading 'isotope_half_life'\n");
  if(upetHeaderReadParameter(fp, "isotope_half_life", tmp)==0) {
    f=-1; (void)sscanf(tmp, "%f", &f); if(f<0) return 15;
    img->isotopeHalflife=f;
  }

  /* branching_fraction */
  rewind(fp);
  if(verbose>1) printf("  reading 'isotope_branching_fraction'\n");
  if(upetHeaderReadParameter(fp, "isotope_branching_fraction", tmp)==0) {
    f=-1; (void)sscanf(tmp, "%f", &f); if(f<0) return 16;
    img->branchingFraction=f;
  }

  /* decay correction applied */
  rewind(fp);
  if(verbose>1) printf("  reading 'decay_correction_applied'\n");
  if(upetHeaderReadParameter(fp, "decay_correction_applied", tmp)==0) {
    n=-1; (void)sscanf(tmp, "%d", &n); if(n<0) return 17;
    if(n==0) img->decayCorrection=IMG_DC_NONCORRECTED;
    else img->decayCorrection=IMG_DC_CORRECTED;
  }

  /* calibration units */
  rewind(fp);
  if(verbose>1) printf("  reading 'calibration_units'\n");
  if(upetHeaderReadParameter(fp, "calibration_units", tmp)==0) {
    n=-1; (void)sscanf(tmp, "%d", &n); if(n<0) return 18;
    switch(n) {
      case 1: img->unit=CUNIT_NCI_PER_ML; break;
      case 2: img->unit=CUNIT_BQ_PER_ML; break;
      case 0:
      default: img->unit=CUNIT_UNKNOWN; break;
    }
  }

  /* calibration factor */
  rewind(fp);
  if(verbose>1) printf("  reading 'calibration_factor'\n");
  if(calibration_factor!=NULL &&
     upetHeaderReadParameter(fp, "calibration_factor", tmp)==0)
  {
    f=-1; (void)sscanf(tmp, "%f", &f); if(f<=0.0) return 19;
    *calibration_factor=f;
    if(img->branchingFraction>0.0) *calibration_factor/=img->branchingFraction;
  }

  /* FOV */
  rewind(fp);
  if(verbose>1) printf("  reading 'radial_fov'\n");
  if(upetHeaderReadParameter(fp, "radial_fov", tmp)==0) {
    f=-1; (void)sscanf(tmp, "%f", &f); if(f<0) return 20;
    img->transaxialFOV=10.0*f;
  }

  return 0;
}
/*****************************************************************************/

/*****************************************************************************/
/** Read frame information from MicroPET header into one-frame-IMG. 
    @return Returns 0 when successful, otherwise >0.
 */
int imgGetMicropetFrameHeader(
  /** File pointer to Concorde/MicroPET header. */
  FILE *fp,
  /** Pointer to IMG struct, allocated for one frame; frame information is
     written in frame 0. */
  IMG *img,
  /** Frame index [0..tdim-1]. */
  int frame_index,
  /** Verbose level. */
  int verbose
) {
  char tmp[MAX_MICROPET_LINE_LEN];
  float f;


  if(verbose>0)
    printf("imgGetMicropetFrameHeader(*fp, *img, %d)\n", frame_index);
  if(fp==NULL) return 1;
  if(img==NULL) return 2;
  if(frame_index<0) return 3;

  /* Search required frame from the beginning of header file */
  rewind(fp);

  /* Find correct frame index */
  sprintf(tmp, "frame %d", frame_index);
  if(verbose>1) printf("  reading '%s'\n", tmp);
  if(upetHeaderReadParameter(fp, tmp, tmp)!=0) return 5;

  /* frame start time */
  if(verbose>1) printf("  reading 'frame_start'\n");
  if(upetHeaderReadParameter(fp, "frame_start", tmp)!=0) return 11;
  f=-1.0; (void)sscanf(tmp, "%f", &f); if(f<0.0) return 11;
  img->start[0]=f;

  /* frame duration; can be 0 at least in single-frame image */
  if(verbose>1) printf("  reading 'frame_duration'\n");
  if(upetHeaderReadParameter(fp, "frame_duration", tmp)!=0) return 12;
  f=-1.0; (void)sscanf(tmp, "%f", &f); if(f<0.0) return 12;
  img->end[0]=img->start[0]+f;
  img->mid[0]=0.5*(img->end[0]+img->start[0]);

  /* scale factor (written in 'weight' since there is no better place) */
  if(verbose>1) printf("  reading 'scale_factor'\n");
  if(upetHeaderReadParameter(fp, "scale_factor", tmp)!=0) return 13;
  f=-1; (void)sscanf(tmp, "%f", &f); if(f<0) return 13;
  img->weight[0]=f;

  /* decay correction */
  if(verbose>1) printf("  reading 'decay_correction'\n");
  if(upetHeaderReadParameter(fp, "decay_correction", tmp)!=0) return 14;
  f=-1; (void)sscanf(tmp, "%f", &f); if(f<0) return 14;
  img->decayCorrFactor[0]=f;

  return 0;
}
/*****************************************************************************/

/*****************************************************************************/
/** Read Scan Information from Concorde/MicroPET header file.
    @return Returns 0 if successful. 
 */
int imgGetMicropetSIF(
  /** File pointer to Concorde/MicroPET header file */
  FILE *fp,
  /** Pointer to initiated but non-allocated SIF struct;
   *  Studynr should be filled afterwards. */
  SIF *sif
) {
  char tmp[MAX_MICROPET_LINE_LEN], tmp2[64], tmp3[64];
  int n, i, ret;


  if(fp==NULL) return 1;
  if(sif==NULL) return 2;


  /* Get frame number */
  rewind(fp);
  if(upetHeaderReadParameter(fp, "total_frames", tmp)!=0) return 11;
  n=-1; (void)sscanf(tmp, "%d", &n); if(n<1) return 11;

  /* Allocate memory for SIF */
  ret=sifSetmem(sif, n); if(ret!=0) return 4;
  sif->frameNr=n;
  sif->colNr=4;
  sif->version=1;

  /* Scan time */
  upetScanStart(fp, &sif->scantime);

  /* Isotope */
  rewind(fp);
  if(upetHeaderReadParameter(fp, "isotope", tmp)!=0) return 13;
  strlcpy(sif->isotope_name, tmp, 8);

  /* Frames */
  for(i=0; i<sif->frameNr; i++) {
    /* Find correct frame index */
    sprintf(tmp, "frame %d", i);
    if(upetHeaderReadParameter(fp, tmp, tmp)!=0) return 21;
    /* frame start time */
    if(upetHeaderReadParameter(fp, "frame_start", tmp)!=0) return 22;
    n=-1; (void)sscanf(tmp, "%d", &n); if(n<0) return 22;
    sif->x1[i]=n;
    /* frame duration */
    if(upetHeaderReadParameter(fp, "frame_duration", tmp)!=0) return 23;
    n=-1; (void)sscanf(tmp, "%d", &n); if(n<0) return 23;
    sif->x2[i]=sif->x1[i]+n;
    /* prompts */
    if(upetHeaderReadParameter(fp, "prompts", tmp)!=0) return 24;
    n=-1; (void)sscanf(tmp, "%s %s %d", tmp2, tmp3, &n); if(n<0) return 24;
    sif->prompts[i]=n;
    /* delays */
    if(upetHeaderReadParameter(fp, "delays", tmp)!=0) return 25;
    n=-1; (void)sscanf(tmp, "%s %s %d", tmp2, tmp3, &n); if(n<0) return 25;
    sif->randoms[i]=n;
    /* trues */
    sif->trues[i]=sif->prompts[i]-sif->randoms[i];
  }
  return 0;
}
/*****************************************************************************/

/*****************************************************************************/
/** Copy microPET header information from IFT struct to IMG struct.
   @return Returns errstatus, which is STATUS_OK (0) when call was successful,
    and >0 in case of an error.
*/
int imgGetMicropetHeader(
  /** Pointer to initiated IMG struct */
  IMG *img
) {
  int i, n;
  float f;
  char key[MAX_MICROPET_LINE_LEN], tmp2[64], tmp3[64];
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wmissing-field-initializers"
  struct tm scanstart={0};
#pragma clang diagnostic pop
  
  if(MICROPET_TEST) printf("\nimgGetMicropetHeader(*img)\n");
  
  /* Check the input */
  if(img==NULL) return STATUS_FAULT;
  if(img->status!=IMG_STATUS_INITIALIZED && img->status!=IMG_STATUS_OCCUPIED)
    return STATUS_FAULT;
  imgSetStatus(img, STATUS_FAULT);
  if(img->ift.keyNr<10) return STATUS_FAULT;
    
  imgSetStatus(img, STATUS_INVALIDHEADER);

  /* Check image format */
  i=iftGetIntValue(&img->ift, 0, "file_type", &n);
  if(i<0) return STATUS_INVALIDHEADER;
  if(n!=5) {imgSetStatus(img, STATUS_UNSUPPORTED); return STATUS_UNSUPPORTED;}

  i=iftGetIntValue(&img->ift, 0, "acquisition_mode", &n);
  if(i<0) return STATUS_INVALIDHEADER;
  if(n!=2 && n!=3) { // CT (9) not yet supported
    imgSetStatus(img, STATUS_UNSUPPORTED); return STATUS_UNSUPPORTED;}

  i=iftGetIntValue(&img->ift, 0, "data_type", &n);
  if(i<0) return STATUS_INVALIDHEADER;
  if(n!=4 && n!=2) {
    imgSetStatus(img, STATUS_UNSUPPORTED); return STATUS_UNSUPPORTED;}

  /* scanner model */
  i=iftGetIntValue(&img->ift, 0, "model", &img->scanner);
  if(i<0) return STATUS_INVALIDHEADER;
  if(img->scanner<0) return STATUS_INVALIDHEADER;

  /* image dimensions */
  i=iftGetIntValue(&img->ift, 0, "total_frames", &n);
  if(i<0) return STATUS_INVALIDHEADER;
  img->dimt=n; if(img->dimt<1) return STATUS_INVALIDHEADER;
  i=iftGetIntValue(&img->ift, 0, "x_dimension", &n);
  if(i<0) return STATUS_INVALIDHEADER;
  img->dimx=n; if(img->dimx<1) return STATUS_INVALIDHEADER;
  i=iftGetIntValue(&img->ift, 0, "y_dimension", &n);
  if(i<0) return STATUS_INVALIDHEADER;
  img->dimy=n; if(img->dimy<1) return STATUS_INVALIDHEADER;
  i=iftGetIntValue(&img->ift, 0, "z_dimension", &n);
  if(i<0) return STATUS_INVALIDHEADER;
  img->dimz=n; if(img->dimz<1) return STATUS_INVALIDHEADER;

  /* zoom */
  i=iftGetFloatValue(&img->ift, 0, "zoom", &img->zoom);
  if(i<0) return STATUS_INVALIDHEADER;
  if(img->zoom<0.0) return STATUS_INVALIDHEADER;

  /* pixel size x */
  i=iftGetFloatValue(&img->ift, 0, "pixel_size_x", &img->sizex);
  if(i>=0) {
    if(img->sizex<0.0) return STATUS_INVALIDHEADER;
  } else {
    i=iftGetFloatValue(&img->ift, 0, "pixel_size", &img->sizex);
    if(i<0 || img->sizex<0.0) return STATUS_INVALIDHEADER;
    img->sizex*=10.0;
  }

  /* pixel size y */
  i=iftGetFloatValue(&img->ift, 0, "pixel_size_y", &img->sizey);
  if(i>=0) {
    if(img->sizey<0.0) return STATUS_INVALIDHEADER;
  } else {
    i=iftGetFloatValue(&img->ift, 0, "pixel_size", &img->sizey);
    if(i<0 || img->sizey<0.0) return STATUS_INVALIDHEADER;
    img->sizey*=10.0;
  }

  /* pixel size z, replaced by transaxial_bin_size, if available */
  /* transaxial_bin_size */
  i=iftGetFloatValue(&img->ift, 0, "pixel_size_z", &img->sizez);
  if(i>=0) {
    if(img->sizez<0.0) return STATUS_INVALIDHEADER;
  } else {
    i=iftGetFloatValue(&img->ift, 0, "axial_plane_size", &img->sizez);
    if(i<0 || img->sizez<0.0) return STATUS_INVALIDHEADER;
    img->sizez*=10.0;
  }
  i=iftGetFloatValue(&img->ift, 0, "transaxial_bin_size", &f);
  if(i>=0 && f>0.0) img->sizez=10.0*f;

  /* isotope halflife */
  i=iftGetFloatValue(&img->ift, 0, "isotope_half_life", &img->isotopeHalflife);
  if(i<0 || img->isotopeHalflife<0.0) return STATUS_INVALIDHEADER;

  /* branching_fraction */
  i=iftGetFloatValue(&img->ift, 0, "isotope_branching_fraction", &img->branchingFraction);
  if(i<0 || img->branchingFraction<0.0) return STATUS_INVALIDHEADER;

  /* decay correction applied */
  i=iftGetIntValue(&img->ift, 0, "decay_correction_applied", &n);
  if(i<0 || n<0.0) return STATUS_INVALIDHEADER;
  if(n==0) img->decayCorrection=IMG_DC_NONCORRECTED;
  else img->decayCorrection=IMG_DC_CORRECTED;

  /* calibration units */
  i=iftGetIntValue(&img->ift, 0, "calibration_units", &n);
  if(i<0 || n<0.0) return STATUS_INVALIDHEADER;
  switch(n) {
    case 1: img->unit=CUNIT_NCI_PER_ML; break;
    case 2: img->unit=CUNIT_BQ_PER_ML; break;
    case 0:
    default: img->unit=CUNIT_UNKNOWN; break;
  }

  /* calibration factor */
  i=iftGetFloatValue(&img->ift, 0, "calibration_factor", &img->calibrationFactor);
  if(i<0 || img->calibrationFactor<0.0) return STATUS_INVALIDHEADER;
  if(img->branchingFraction>0.0) img->calibrationFactor/=img->branchingFraction;

  /* FOV */
  i=iftGetFloatValue(&img->ift, 0, "radial_fov", &img->transaxialFOV);
  if(i<0) return STATUS_INVALIDHEADER;
  img->transaxialFOV*=10.0;

  /* General */
  img->_fileFormat=IMG_MICROPET;
  img->type=IMG_TYPE_IMAGE;

  /* Studynumber, if possible */
  strcpy(key, "study"); n=1;
  if((i=iftGet(&img->ift, key))>=0) 
    n=studynr_from_fname2(img->ift.item[i].value, img->studyNr, 0);
  if(n!=0) {
    strcpy(key, "file_name");
    if((i=iftGet(&img->ift, key))>=0)
      n=studynr_from_fname2(img->ift.item[i].value, img->studyNr, 0);
  }
  if(n!=0 && MICROPET_TEST>1) printf("Valid studyNr could not be read.\n");

  /* Scan start */
  strcpy(key, "scan_time");
  if((i=iftGet(&img->ift, key))<0) return STATUS_INVALIDHEADER;
  n=sscanf(img->ift.item[i].value, "%s %s %d %d:%d:%d %d",
           tmp2, tmp3, &scanstart.tm_mday,
           &scanstart.tm_hour, &scanstart.tm_min, &scanstart.tm_sec, &i);
  if(n==7) {
    scanstart.tm_year=i-1900;
    if(strcasecmp(tmp3, "Jan")==0)      scanstart.tm_mon=0;
    else if(strcasecmp(tmp3, "Feb")==0) scanstart.tm_mon=1;
    else if(strcasecmp(tmp3, "Mar")==0) scanstart.tm_mon=2;
    else if(strcasecmp(tmp3, "Apr")==0) scanstart.tm_mon=3;
    else if(strcasecmp(tmp3, "May")==0) scanstart.tm_mon=4;
    else if(strcasecmp(tmp3, "Jun")==0) scanstart.tm_mon=5;
    else if(strcasecmp(tmp3, "Jul")==0) scanstart.tm_mon=6;
    else if(strcasecmp(tmp3, "Aug")==0) scanstart.tm_mon=7;
    else if(strcasecmp(tmp3, "Sep")==0) scanstart.tm_mon=8;
    else if(strcasecmp(tmp3, "Oct")==0) scanstart.tm_mon=9;
    else if(strcasecmp(tmp3, "Nov")==0) scanstart.tm_mon=10;
    else if(strcasecmp(tmp3, "Dec")==0) scanstart.tm_mon=11;
    scanstart.tm_isdst=-1;
    img->scanStart=timegm(&scanstart); //mktime(&scanstart);
    if(img->scanStart<0) img->scanStart=0;
  } else return STATUS_INVALIDHEADER;


  imgSetStatus(img, STATUS_OK);
  return STATUS_OK;
}
/******************************************************************************/

/******************************************************************************/
/** Fill IMG struct header information from microPET database files.

    Information concerning separate frames or planes is not filled though.

   @return Returns errstatus, which is STATUS_OK (0) when call was successful,
    and >0 in case of an error.
 */
int imgReadMicropetHeader(
  /** Name of microPET database, may contain filename extension */
  const char *dbname,
  /** Pointer to the initiated IMG data */
  IMG *img
) {
  char hdrfile[FILENAME_MAX];
  int ret;

  if(IMG_TEST) printf("\nimgReadMicropetHeader(%s, *img)\n", dbname);
  
  /* Check the input */
  if(img==NULL) return STATUS_FAULT;
  if(img->status!=IMG_STATUS_INITIALIZED) return STATUS_FAULT;
  imgSetStatus(img, STATUS_FAULT);
  if(dbname==NULL) return STATUS_FAULT;

  /* Determine the names of hdr and sif files */
  ret=upetExists(dbname, hdrfile, NULL, IMG_TEST-1);
  if(ret==0) return STATUS_NOFILE;
  
  /* Read microPET header file into IFT */
  iftEmpty(&img->ift);
  ret=defRead(&img->ift, hdrfile);
  if(ret!=0) {
    if(IMG_TEST>1) printf("defRead() return value := %d\n", ret);
    return(STATUS_FAULT);
  }
  /* and set IMG contents */
  ret=imgGetMicropetHeader(img);
  if(ret!=0) {
    imgSetStatus(img, ret);
    return(ret);
  }

  return STATUS_OK;
}
/******************************************************************************/

/******************************************************************************/
/** Read a specified frame from microPET image into preallocated IMG
    data structure. 
    
    MicroPET image consists of two files in the same directory: 
    fname.hdr and fname.img.
    IMG header is assumed to be filled correctly before calling this function,
    except for information concerning separate planes and this frame,
    which is filled here.
    If frame does not exist, then and only then STATUS_NOMATRIX is returned.
  
  @return Returns errstatus, which is STATUS_OK (0) when call was successful,
   and >0 in case of an error.
 */ 
int imgReadMicropetFrame(
  /** Name of microPET image (hdr or img file, or without extension)
      from which IMG contents will be read */
  const char *fname,
  /** Frame which will be read [1..frameNr] */
  int frame_to_read,
  /** Pointer to the IMG data. Place for the frame must be preallocated */
  IMG *img,
  /** IMG frame index [0..dimt-1] where data will be placed */
  int frame_index
) {
  FILE *fp;
  int i, fi, ret, zi, xi, yi;
  float *fdata=NULL, *fptr, f;
  char datfile[FILENAME_MAX], hdrfile[FILENAME_MAX];
  char key[MAX_MICROPET_LINE_LEN], value[MAX_MICROPET_LINE_LEN];


  if(IMG_TEST) printf("\nimgReadMicropetFrame(%s, %d, *img, %d)\n",
    fname, frame_to_read, frame_index);
    
  /* Check the input */
  if(img==NULL) return STATUS_FAULT;
  if(img->status!=IMG_STATUS_OCCUPIED) return STATUS_FAULT;
  if(fname==NULL) return STATUS_FAULT;
  if(frame_index<0 || frame_index>img->dimt-1) return STATUS_FAULT;
  if(frame_to_read<1) return STATUS_FAULT;
  imgSetStatus(img, STATUS_FAULT);
  
  /* Determine the names of hdr and data files */
  ret=upetExists(fname, hdrfile, datfile, IMG_TEST-1);
  if(ret<2) {imgSetStatus(img, STATUS_NOFILE); return STATUS_NOFILE;}

  /* Read microPET header file into IFT, if not available already */
  imgSetStatus(img, STATUS_INVALIDHEADER);
  if(img->ift.keyNr<10) {
    iftEmpty(&img->ift);
    ret=defRead(&img->ift, hdrfile);
    if(ret!=0) {
      if(IMG_TEST>1) printf("defRead() return value := %d\n", ret);
      return(STATUS_INVALIDHEADER);
    }
    if(IMG_TEST>3) printf("ift.keyNr := %d\n", img->ift.keyNr);
  }

  /* Read frame header information */
  strcpy(key, "frame"); sprintf(value, "%d", frame_to_read-1);
  fi=iftGetFullmatchFrom(&img->ift, 0, key, value);
  if(fi<0) {imgSetStatus(img, STATUS_NOMATRIX); return STATUS_NOMATRIX;}
  i=iftGetFloatValue(&img->ift, fi+1, "frame_start", &f);
  if(i<0 || isnan(f)) {return STATUS_INVALIDHEADER;}
  img->start[frame_index]=f;
  i=iftGetFloatValue(&img->ift, fi+1, "frame_duration", &f);
  if(i<0 || isnan(f) || f<0.0) {return STATUS_INVALIDHEADER;}
  img->end[frame_index]=img->start[frame_index]+f;
  img->mid[frame_index]=0.5*(img->end[frame_index]+img->start[frame_index]);
  /* decay correction */
  i=iftGetFloatValue(&img->ift, fi+1, "decay_correction", &f);
  if(i<0 || isnan(f) || f<0.0) {return STATUS_INVALIDHEADER;}
  img->decayCorrFactor[frame_index]=f;
  //printf("   decayCorrFactor=%g\n", img->decayCorrFactor[0]);
  /* set plane numbers */
  for(zi=0; zi<img->dimz; zi++) img->planeNumber[zi]=zi+1;
  /* prompts and randoms (delays), without checking the result */
  i=iftGetFloatValue(&img->ift, fi+1, "prompts_rate", &img->prompts[frame_index]);
  i=iftGetFloatValue(&img->ift, fi+1, "delays_rate", &img->randoms[frame_index]);

  /* Open image datafile */
  if(IMG_TEST>2) fprintf(stdout, "reading image data %s\n", datfile);
  if((fp=fopen(datfile, "rb"))==NULL) {
    imgSetStatus(img, STATUS_NOIMGDATA); return STATUS_NOIMGDATA;}

  /* Allocate memory for one image frame */
  fdata=malloc(img->dimx*img->dimy*img->dimz*sizeof(float));
  if(fdata==NULL) {
    fclose(fp); imgSetStatus(img, STATUS_NOMEMORY);
    return STATUS_NOMEMORY;
  }

  /* Read the required image frame */
  fptr=fdata;
  ret=upetReadImagedata(fp, &img->ift, frame_to_read, fptr);
  if(ret!=0 && IMG_TEST) fprintf(stdout, "upetReadImagedata() := %d\n", ret);
  fclose(fp);
  if(ret==3) { /* no more frames */
    free(fdata); imgSetStatus(img, STATUS_NOMATRIX); return STATUS_NOMATRIX;}
  if(ret!=0) {
    free(fdata); imgSetStatus(img, STATUS_UNSUPPORTED); return STATUS_UNSUPPORTED;}

  /* Copy pixel values to IMG */
  fptr=fdata;
  for(zi=0; zi<img->dimz; zi++)
    for(yi=0; yi<img->dimy; yi++)
      for(xi=0; xi<img->dimx; xi++)
        img->m[zi][yi][xi][frame_index]=*fptr++;
  free(fdata);

  imgSetStatus(img, STATUS_OK); /* If the rest is failed, no problem */
  return STATUS_OK;
}
/******************************************************************************/

/******************************************************************************/
/** Read the first frame from a microPET image into IMG data structure.
   @return Returns errstatus, which is STATUS_OK (0) when call was successful,
    and >0 in case of an error.
 */
int imgReadMicropetFirstFrame(
  /** Name of microPET image (hdr or img file, or without extension)
      from which IMG contents will be read */
  const char *fname,
  /** Pointer to the initiated but not preallocated IMG data */
  IMG *img
) {
  int ret=0;

  if(IMG_TEST) printf("\nimgReadMicropetFirstFrame(%s, *img)\n", fname);
  /* Check the input */
  if(img==NULL) return STATUS_FAULT;
  if(img->status!=IMG_STATUS_INITIALIZED) return STATUS_FAULT;
  imgSetStatus(img, STATUS_FAULT);
  if(fname==NULL) return STATUS_FAULT;

  /* Read header information from file */
  ret=imgReadMicropetHeader(fname, img);
  if(IMG_TEST>1) printf("imgReadMicropetHeader() := %s\n", img->statmsg);
  if(ret) return(ret);
  if(IMG_TEST>3) imgInfo(img);

  /* Allocate memory for one frame */
  img->dimt=1;
  ret=imgAllocate(img, img->dimz, img->dimy, img->dimx, img->dimt);
  if(ret) return STATUS_NOMEMORY;

  /* Read the first frame */
  ret=imgReadMicropetFrame(fname, 1, img, 0);
  if(IMG_TEST>1) printf("imgReadMicropetFrame() := %s\n", img->statmsg);
  if(ret) return(ret); 

  /* All went well */
  imgSetStatus(img, STATUS_OK);
  return STATUS_OK;
}
/******************************************************************************/

/******************************************************************************/
/** Read the whole dynamic microPET image into IMG data structure. 

    Note that microPET images are often too large for 32-bit systems.

    @return Returns errstatus, which is STATUS_OK (0) when call was successful,
     and >0 in case of an error.
 */
int imgReadMicropet(
  /** Name of microPET image (hdr or img file, or without extension)
      from which IMG contents will be read */
  const char *fname,
  /** Pointer to the initiated but not preallocated IMG data */
  IMG *img
) {
  int fi=0, ret=0;

  if(IMG_TEST) printf("\nimgReadMicropet(%s, *img)\n", fname);
  /* Check the input */
  if(img==NULL) return STATUS_FAULT;
  if(img->status!=IMG_STATUS_INITIALIZED) return STATUS_FAULT;
  imgSetStatus(img, STATUS_FAULT);
  if(fname==NULL) return STATUS_FAULT;

  /* Read the first frame into IMG struct */
  fi=0;
  if(IMG_TEST>2) printf("reading frame %d\n", fi+1);
  ret=imgReadMicropetFirstFrame(fname, img);
  if(ret!=STATUS_OK) {
    if(IMG_TEST>0) printf("imgReadMicropetFirstFrame() := %s\n", img->statmsg);
    imgEmpty(img); return ret;
  }
  /* Read rest of the frames */
  do {
    fi++;
    if(IMG_TEST>2) printf("reading frame %d\n", fi+1);
    ret=imgReadMicropetFrame(fname, fi+1, img, fi);
  } while(ret==0);
  if(ret!=STATUS_OK && ret!=STATUS_NOMATRIX) {
    if(IMG_TEST>0) printf("imgReadMicropetFrame() := %s\n", img->statmsg);
    imgEmpty(img); return ret;
  }
  if(IMG_TEST>1) printf("%d frame(s) were read.\n", fi);

  /* All went well */
  imgSetStatus(img, STATUS_OK);
  return STATUS_OK;
}
/******************************************************************************/

/******************************************************************************/
