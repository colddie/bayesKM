/// @file img_nii.c
/// @author Vesa Oikonen
/// @brief NIfTI-1 PET image I/O routines for IMG data.
///
///  Function are not intended to support all NIfTI files or file
///  properties, but only those that have been found necessary in
///  Turku PET Centre.
///
/*****************************************************************************/
#include "libtpcimgio.h"
/*****************************************************************************/

/*****************************************************************************/
/** Read Nifti-1 image.

    Nifti database name must be given with path. Either one file with
    extension .nii, or image and header files with extensions .img and .hdr must exist.
    Also SIF file with .sif extension is used, if it exists.
  
    @return 0 if ok, and otherwise IMG status code (<>0); sets IMG->statmsg in case of an error.
    @sa imgInit, imgReadNiftiFirstFrame, imgReadNiftiHeader, imgWriteNifti
 */
int imgReadNifti(
  /** Nifti database name with path, with or without extension */
  const char *filename,
  /** Pointer to initialized IMG structure */
  IMG *img,
  /** Verbose level; if zero, then nothing is printed to stderr or stdout */
  int verbose
) {
  int fi, ret;

  if(verbose>0) {printf("imgReadNifti(%s, ...)\n", filename); fflush(stdout);}

  /* Check the arguments */
  if(img==NULL || img->status!=IMG_STATUS_INITIALIZED) {
    if(img!=NULL) imgSetStatus(img, STATUS_FAULT);
    if(verbose>0) fprintf(stderr, "Error: invalid IMG argument\n");
    return(STATUS_FAULT);
  }
  if(filename==NULL || !filename[0]) {
    imgSetStatus(img, STATUS_FAULT);
    if(verbose>0) fprintf(stderr, "Error: invalid filename\n");
    return(STATUS_FAULT);
  }
  imgSetStatus(img, STATUS_OK);

  /* Read header information from file */
  ret=imgReadNiftiHeader(filename, img, verbose-1); if(ret) {return(ret);}
  if(verbose>10) imgInfo(img);

  /* Allocate memory for all frames */
  ret=imgAllocate(img, img->dimz, img->dimy, img->dimx, img->dimt);
  if(ret) {
    imgSetStatus(img, STATUS_OK);
    return STATUS_NOMEMORY;
  }

  /* Read one frame at a time */
  for(fi=0; fi<img->dimt; fi++) {
    ret=imgReadNiftiFrame(filename, 1+fi, img, fi, verbose-1);
    if(ret) return(ret);
  }

  /* All went well */
  imgSetStatus(img, STATUS_OK);
  return(STATUS_OK);
}
/*****************************************************************************/

/*****************************************************************************/
/** Read the first frame of Nifti-1 image into IMG data structure.

    Nifti database name must be given with path. Either one file with
    extension .nii, or image and header files with extensions .img and .hdr must exist.
    Also SIF file with .sif extension is used, if it exists.
 
    @return 0 if ok, and otherwise IMG status code (<>0); sets IMG->statmsg in case of an error.
 */
int imgReadNiftiFirstFrame(
  /** Nifti database name with path, with or without extension */
  const char *filename,
  /** Pointer to initialized IMG structure */
  IMG *img,
  /** Verbose level; if zero, then nothing is printed to stderr or stdout */
  int verbose
) {
  int ret;

  if(verbose>0) {printf("imgReadNiftiFirstFrame(%s, ...)\n", filename); fflush(stdout);}

  /* Check the arguments */
  if(img==NULL || img->status!=IMG_STATUS_INITIALIZED) {
    if(img!=NULL) imgSetStatus(img, STATUS_FAULT);
    if(verbose>0) fprintf(stderr, "Error: invalid IMG argument\n");
    return(STATUS_FAULT);
  }
  if(filename==NULL || !filename[0]) {
    imgSetStatus(img, STATUS_FAULT);
    if(verbose>0) fprintf(stderr, "Error: invalid filename\n");
    return(STATUS_FAULT);
  }
  imgSetStatus(img, STATUS_OK);

  /* Read header information from file */
  ret=imgReadNiftiHeader(filename, img, verbose-1); if(ret) {return(ret);}
  if(verbose>10) imgInfo(img);

  /* Allocate memory for one frame */
  img->dimt=1;
  ret=imgAllocate(img, img->dimz, img->dimy, img->dimx, img->dimt);
  if(ret) {
    imgSetStatus(img, STATUS_OK);
    return STATUS_NOMEMORY;
  }

  /* Read the first frame */
  ret=imgReadNiftiFrame(filename, 1, img, 0, verbose-1); if(ret) return(ret); 

  /* All went well */
  imgSetStatus(img, STATUS_OK);
  return(STATUS_OK);
}
/*****************************************************************************/

/*****************************************************************************/
/** Fill IMG struct header information from Nifti database files.

    Nifti database name must be given with path. Either one file with
    extension .nii, or image and header files with extensions .img and .hdr
    must exist.
    Also SIF file with .sif extension is used, if it exists; however,
    information concerning separate frames or planes is not filled.
  
    @return 0 if ok, and otherwise IMG status code (<>0); sets IMG->statmsg in case of an error.
 */
int imgReadNiftiHeader(
  /** Nifti database name with path, with or without extension */
  const char *filename,
  /** Pointer to initialized IMG structure */
  IMG *img,
  /** Verbose level; if zero, then nothing is printed to stderr or stdout */
  int verbose
) {
  char basefile[FILENAME_MAX], tmp[256];
  char datfile[FILENAME_MAX], hdrfile[FILENAME_MAX], siffile[FILENAME_MAX];
  NIFTI_DSR dsr;
  int ret;
  SIF sif;
  double f;

  if(verbose>0) {
    printf("imgReadNiftiHeader(%s, ...)\n", filename); fflush(stdout);}

  /* Check the arguments */
  if(img==NULL || img->status!=IMG_STATUS_INITIALIZED) {
    if(img!=NULL) imgSetStatus(img, STATUS_FAULT);
    if(verbose>0) fprintf(stderr, "Error: invalid IMG argument\n");
    return(STATUS_FAULT);
  }
  if(filename==NULL || !filename[0]) {
    imgSetStatus(img, STATUS_FAULT);
    if(verbose>0) fprintf(stderr, "Error: invalid filename\n");
    return(STATUS_FAULT);
  }
  imgSetStatus(img, STATUS_OK);

  /* Extract the base file name without extensions */
  strcpy(basefile, filename); niftiRemoveFNameExtension(basefile);
  if(strlen(basefile)<1) {
    imgSetStatus(img, STATUS_FAULT);
    if(verbose>0) fprintf(stderr, "Error: invalid filename\n");
    return(STATUS_FAULT);
  }

  /* Make the image and header filenames, and read Nifti header */
  ret=niftiExists(basefile, hdrfile, datfile, siffile, &dsr, verbose-2, tmp);
  if(ret==0) {
    imgSetStatus(img, STATUS_NOFILE);
    if(verbose>0) fprintf(stderr, "Error: %s\n", tmp);
    return(STATUS_NOFILE);
  }
  if(ret==1 && verbose>1) printf("no SIF found for %s\n", basefile); 

  /* and set IMG contents */
  ret=imgGetNiftiHeader(img, &dsr, verbose-2);
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
  if(verbose>1) {printf("reading SIF %s\n", siffile); fflush(stdout);}
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

  return(STATUS_OK);
}
/*****************************************************************************/

/*****************************************************************************/
/** Copy Nifti header information into IMG.
  @return Returns IMG status, which is STATUS_OK (0) when call was successful,
   and >0 in case of an error.
 */
int imgGetNiftiHeader(
  /** Pointer to IMG struct */
  IMG *img,
  /** Pointer to Nifti header contents */
  NIFTI_DSR *dsr,
  /** Verbose level; if zero, then nothing is printed to stderr or stdout */
  int verbose
) {
  int i, dimNr, dimx, dimy, dimz=1, dimt=1, pxlNr=0;
  float f;

  if(verbose>0) {printf("imgGetNiftiHeader()\n"); fflush(stdout);}
  
  /* Check the input */
  if(img==NULL) return STATUS_FAULT;
  if(img->status!=IMG_STATUS_INITIALIZED && img->status!=IMG_STATUS_OCCUPIED)
    return STATUS_FAULT;
  imgSetStatus(img, STATUS_FAULT);
  if(dsr==NULL) return STATUS_FAULT;
    
  imgSetStatus(img, STATUS_INVALIDHEADER);

  /* Get the image dimensions from header */
  dimNr=dsr->h.dim[0];
  if(dimNr<2 || dimNr>4) {
    if(verbose>0) fprintf(stderr, "Error: Nifti image dimension %d is not supported\n", dimNr);
    return STATUS_UNSUPPORTED;
  }
  dimx=dsr->h.dim[1]; dimy=dsr->h.dim[2];
  if(dimNr>2) {dimz=dsr->h.dim[3]; if(dimNr>3) dimt=dsr->h.dim[4];}
  pxlNr=dimx*dimy*dimz;
  if(pxlNr<1) {
    if(verbose>0) fprintf(stderr, "Error: invalid Nifti image dimensions.\n");
    return STATUS_INVALIDHEADER;
  }
  img->dimx=dimx; img->dimy=dimy; img->dimz=dimz; img->dimt=dimt;

  /* Copy information from header */
  img->type=IMG_TYPE_IMAGE;
  strcpy(img->studyNr, "");
  strcpy(img->patientName, "");
  /* NIfTI file format */
  if(strcmp(dsr->h.magic, "ni1")==0) img->_fileFormat=IMG_NIFTI_1D;
  else if(strcmp(dsr->h.magic, "n+1")==0) img->_fileFormat=IMG_NIFTI_1S;
  else {
    if(verbose>0) fprintf(stderr, "Error: invalid Nifti magic number.\n");
    return STATUS_INVALIDHEADER;
  }
  /* NIfTI datatype is not needed, because currently saved as floats anyway */
  //img->_dataType=dsr->h.datatype;
  /* Pixel x,y,z sizes, converting units if necessary */
  f=1.0;
  if(dsr->h.xyzt_units & NIFTI_UNITS_METER) f=1000.;
  else if(dsr->h.xyzt_units & NIFTI_UNITS_MICRON) f=0.001;
  else if(dsr->h.xyzt_units & NIFTI_UNITS_MM) f=1.0;
  if(verbose>2) printf("pixel size conversion factor := %g\n", f);
  for(i=1; i<=3; i++) {
    if(i==1) img->sizex=f*dsr->h.pixdim[i];
    else if(i==2) img->sizey=f*dsr->h.pixdim[i];
    else if(i==3) img->sizez=f*dsr->h.pixdim[i];
  }
  /* Orientation, quaternion, and transformation parameters */
  img->xform[0]=dsr->h.qform_code;
  img->xform[1]=dsr->h.sform_code;
  img->quatern[0]=dsr->h.quatern_b;
  img->quatern[1]=dsr->h.quatern_c;
  img->quatern[2]=dsr->h.quatern_d;
  img->quatern[3]=dsr->h.qoffset_x;
  img->quatern[4]=dsr->h.qoffset_y;
  img->quatern[5]=dsr->h.qoffset_z;
  for(i=0; i<4; i++) img->quatern[6+i]=dsr->h.srow_x[i];
  for(i=0; i<4; i++) img->quatern[10+i]=dsr->h.srow_y[i];
  for(i=0; i<4; i++) img->quatern[14+i]=dsr->h.srow_z[i];

  /* Assumptions */
  img->decayCorrection=IMG_DC_CORRECTED;

  imgSetStatus(img, STATUS_OK);
  return STATUS_OK;
}
/*****************************************************************************/

/*****************************************************************************/
/** Read a specified frame from a Nifti database files into preallocated IMG
    data structure.

    IMG header is assumed to be filled correctly before calling this function,
    except for information concerning separate planes and this frame,
    which is filled here.

    @return 0 if ok, and otherwise IMG status code (<>0); sets IMG->statmsg
     in case of an error.
     If frame does not exist, then and only then STATUS_NOMATRIX is returned.
 */ 
int imgReadNiftiFrame(
  /** Pointer to string that contains the name of NIfTI database with path,
      with or without filename extension. Either one file with
      extension .nii, or image and header files with extensions .img and .hdr
      must exist.
      Also SIF file with .sif extension is used, if it exists; however,
      information concerning separate frames or planes is not filled. */
  const char *filename,
  /** Frame which will be read from database [1..frameNr] */
  int frame_to_read,
  /** Pointer to the IMG data. Place for the frame must be preallocated */
  IMG *img,
  /** IMG frame index [0..dimt-1] where data will be placed */
  int frame_index,
  /** Verbose level; if zero, then nothing is printed to stderr or stdout */
  int verbose
) {
  char basefile[FILENAME_MAX], tmp[256];
  char datfile[FILENAME_MAX], hdrfile[FILENAME_MAX], siffile[FILENAME_MAX];
  NIFTI_DSR dsr;
  int ret, zi, yi, xi, fi;
  SIF sif;
  FILE *fp;
  float *fdata=NULL, *fptr;


  if(verbose>0) {
    printf("\nimgReadNiftiFrame(%s, %d, *img, %d, %d)\n",
           filename, frame_to_read, frame_index, verbose);
    fflush(stdout);
  }

  /* Check the input */
  if(img==NULL || img->status!=IMG_STATUS_OCCUPIED) {
    if(img!=NULL) imgSetStatus(img, STATUS_FAULT);
    if(verbose>0) fprintf(stderr, "Error: invalid IMG argument\n");
    return(STATUS_FAULT);
  }
  if(filename==NULL || !filename[0]) {
    imgSetStatus(img, STATUS_FAULT);
    if(verbose>0) fprintf(stderr, "Error: invalid filename\n");
    return(STATUS_FAULT);
  }
  if(frame_index<0 || frame_index>img->dimt-1 || frame_to_read<1) {
    imgSetStatus(img, STATUS_FAULT);
    if(verbose>0) fprintf(stderr, "Error: invalid frame settings\n");
    return STATUS_FAULT;
  }

  /* Extract the base file name without extensions */
  strcpy(basefile, filename); niftiRemoveFNameExtension(basefile);
  if(strlen(basefile)<1) {
    imgSetStatus(img, STATUS_FAULT);
    if(verbose>0) fprintf(stderr, "Error: invalid filename\n");
    return(STATUS_FAULT);
  }

  /* Make the image and header filenames, and read Nifti header */
  ret=niftiExists(basefile, hdrfile, datfile, siffile, &dsr, verbose-2, tmp);
  if(ret==0) {
    imgSetStatus(img, STATUS_NOFILE);
    if(verbose>0) fprintf(stderr, "Error: %s\n", tmp);
    return(STATUS_NOFILE);
  }
  if(ret==1 && verbose>1) {printf("no SIF found for %s\n", basefile); fflush(stdout);}

  /* Open image datafile */
  if(verbose>2) {
    fprintf(stdout, "reading image data %s\n", datfile); fflush(stdout);}
  imgSetStatus(img, STATUS_NOIMGDATA);
  if((fp=fopen(datfile, "rb")) == NULL) return STATUS_NOIMGDATA;

  /* Allocate memory for one image frame */
  imgSetStatus(img, STATUS_NOMEMORY);
  fdata=calloc(img->dimx*img->dimy*img->dimz, sizeof(float));
  if(fdata==NULL) {fclose(fp); return STATUS_NOMEMORY;}

  /* Read the required image frame */
  fptr=fdata;
  ret=niftiReadImagedata(fp, &dsr, frame_to_read, fptr, verbose-1, tmp);
  if(verbose>1) printf("niftiReadImagedata() -> %s\n", tmp);
  fclose(fp);
  if(ret==-1) { /* no more frames */
    free(fdata); imgSetStatus(img, STATUS_NOMATRIX);
    return STATUS_NOMATRIX;
  }
  if(ret!=0) {
    free(fdata); imgSetStatus(img, STATUS_UNSUPPORTED);
    return STATUS_UNSUPPORTED;
  }

  /* Copy pixel values to IMG */
  fptr=fdata;
  for(zi=0; zi<img->dimz; zi++)
    for(yi=0; yi<img->dimy; yi++)
      for(xi=0; xi<img->dimx; xi++)
        img->m[zi][yi][xi][frame_index]=*fptr++;
  free(fdata);

  /* Set decay correction factor to zero */
  img->decayCorrFactor[frame_index]=0.0;

  /* Set plane numbers */
  for(fi=0; fi<img->dimz; fi++) img->planeNumber[fi]=fi+1;

  /* 
   *  Try to read frame time information from SIF file
   */
  imgSetStatus(img, STATUS_OK); /* If the rest is failed, no problem */
  sifInit(&sif);
  if(sifRead(siffile, &sif)!=0) {
    if(verbose>1) {fprintf(stdout, "  cannot read SIF (%s)\n", siffile); fflush(stdout);}
    return STATUS_OK;
  }
  /* Frame information */
  if(verbose>3) {fprintf(stdout, "  setting frame times\n"); fflush(stdout);}
  if(sif.frameNr>=frame_to_read) {
    img->start[frame_index]=sif.x1[frame_to_read-1];
    img->end[frame_index]=sif.x2[frame_to_read-1];
    img->mid[frame_index]=0.5*(img->start[frame_index]+img->end[frame_index]);
    img->prompts[frame_index]=sif.prompts[frame_to_read-1];
    img->randoms[frame_index]=sif.randoms[frame_to_read-1];
  }
  sifEmpty(&sif);

  imgSetStatus(img, STATUS_OK);
  return STATUS_OK;
}
/*****************************************************************************/

/*****************************************************************************/
/** Copy header information in IMG struct into NIfTI header struct.

    Min, max, and scale factor are set here and they apply to all frames.

   @return errstatus, which is STATUS_OK (0) when call was successful,
    and >0 in case of an error.
 */
int imgSetNiftiHeader(
  /** pointer to IMG struct from which header information is read */
  IMG *img,
  /** NIfTI database name */
  const char *dbname,
  /** pointer to NIfTI header struct to be filled */
  NIFTI_DSR *dsr,
  /** minimum pixel value in all frames that will be written */
  float fmin,
  /** maximum pixel value in all frames that will be written */
  float fmax,
  /** Verbose level; if zero, then nothing is printed to stderr or stdout */
  int verbose
) {
  int i;
  char *cptr;

  if(verbose>0) {
    printf("\nimgSetNiftiHeader(*img, %s, *dsr, %g, %g, ...)\n", dbname, fmin, fmax);
    fflush(stdout);
  }
  /* Check the input */
  if(dbname==NULL || !dbname[0]) {
    if(verbose>0) fprintf(stderr, "Error: invalid filename\n");
    return(STATUS_FAULT);
  }
  if(img==NULL || (img->status!=IMG_STATUS_INITIALIZED && img->status!=IMG_STATUS_OCCUPIED)) {
    if(verbose>0) fprintf(stderr, "Error: invalid IMG argument\n");
    return(STATUS_FAULT);
  }
  if(dsr==NULL) {
    if(verbose>0) fprintf(stderr, "Error: invalid header struct\n");
    return STATUS_FAULT;
  }

  /* Set NIfTI byte order to current machines byte order */
  dsr->byte_order=little_endian();

  /* Initiate header struct with zeroes */
  memset(&dsr->h, 0, sizeof(NIFTI_1_HEADER));
  memset(&dsr->e, 0, sizeof(NIFTI_EXTENDER));

  /* Set header */
  dsr->h.sizeof_hdr=NIFTI_HEADER_SIZE;
  strcpy(dsr->h.data_type, "");
  cptr=strrchr(dbname, '/'); if(cptr==NULL) cptr=strrchr(dbname, '\\');
  if(cptr!=NULL) cptr++;
  if(cptr==NULL) cptr=(char*)dbname;
  strncpy(dsr->h.db_name, cptr, 17);
  dsr->h.extents=16384; // not used in NIfTI, but required for Analyze compatibility
  dsr->h.regular='r'; // not used in NIfTI, but required for Analyze compatibility
  dsr->h.dim_info='\0'; // MRI slice ordering

  /* Image dimension */
  for(i=0; i<8; i++) dsr->h.dim[i]=1;
  dsr->h.dim[0]=4;
  dsr->h.dim[1]=img->dimx;
  dsr->h.dim[2]=img->dimy;
  dsr->h.dim[3]=img->dimz;
  dsr->h.dim[4]=img->dimt;
  dsr->h.intent_p1=0.0;
  dsr->h.intent_p2=0.0;
  dsr->h.intent_p3=0.0;
  dsr->h.intent_code=NIFTI_INTENT_NONE;
  dsr->h.datatype=NIFTI_DT_FLOAT; // data as floats, so no need to scale
  dsr->h.bitpix=32;
  dsr->h.slice_start=0;
  for(i=0; i<8; i++) dsr->h.pixdim[i]=0.0;
  // https://nifti.nimh.nih.gov/nifti-1/documentation/nifti1fields/nifti1fields_pages/qsform.html
  dsr->h.pixdim[0]=1.0; // Set to either 1.0 or -1.0
  dsr->h.pixdim[1]=img->sizex;
  dsr->h.pixdim[2]=img->sizey;
  dsr->h.pixdim[3]=img->sizez;
  if(img->_fileFormat==IMG_NIFTI_1D)
    dsr->h.vox_offset=0;
  else
    dsr->h.vox_offset=352;
  dsr->h.scl_slope=1.0; // data as floats, so no need to scale
  dsr->h.scl_inter=0.0; // data as floats, so no need to scale
  dsr->h.slice_end=0;
  dsr->h.slice_code=0;
  dsr->h.xyzt_units=NIFTI_UNITS_MM+NIFTI_UNITS_SEC;
  dsr->h.cal_max=fmax;
  dsr->h.cal_min=0.0;
  dsr->h.slice_duration=0.0;
  dsr->h.toffset=0.0;
  dsr->h.glmax=fmax; // unused in NIfTI
  dsr->h.glmin=fmin; // unused in NIfTI

  strlcpy(dsr->h.descrip, img->studyNr, 80);
  strcpy(dsr->h.aux_file, "");

  dsr->h.qform_code=img->xform[0];
  dsr->h.sform_code=img->xform[1];
  dsr->h.quatern_b=img->quatern[0];
  dsr->h.quatern_c=img->quatern[1];
  dsr->h.quatern_d=img->quatern[2];
  dsr->h.qoffset_x=img->quatern[3];
  dsr->h.qoffset_y=img->quatern[4];
  dsr->h.qoffset_z=img->quatern[5];
  for(i=0; i<4; i++) dsr->h.srow_x[i]=img->quatern[6+i];
  for(i=0; i<4; i++) dsr->h.srow_y[i]=img->quatern[10+i];
  for(i=0; i<4; i++) dsr->h.srow_z[i]=img->quatern[14+i];
  strcpy(dsr->h.intent_name, "");

  if(img->_fileFormat==IMG_NIFTI_1D) strcpy(dsr->h.magic, "ni1");
  else strcpy(dsr->h.magic, "n+1");

  /* Extension is left as 0 0 0 0 */

  return STATUS_OK;
}
/*****************************************************************************/

/*****************************************************************************/
/** Write one PET frame from IMG data struct into NIfTI file.

    This function can be called repeatedly to write all frames one at a time
    to conserve memory. This function does not write SIF.
    Single or dual file format is determined based on _fileFormat field in
    IMG struct. Byte order is the not changed.

    @return Returns errstatus, which is STATUS_OK (0) when call was successful,
     and >0 in case of an error.
 */
int imgWriteNiftiFrame(
  /** Name of NIfTI file where IMG contents will be written.
      If file does not exist, it is created.
      Make sure to delete existing file, unless you want to add data. */
  const char *dbname,
  /** PET frame number (1..frameNr) which will be written:
      If set to 0, frame data will be written to an existing or new PET file as
      a new additional frame, never overwriting existing data.
      If >0, then frame data is written as specified frame number, overwriting
      any data existing with the same frame number */
  int frame_to_write,
  /** pointer to the IMG data struct. */
  IMG *img,
  /** IMG frame index (0..dimt-1) which will be written. */
  int frame_index,
  /** minimum pixel value in all frames that will be written;
      used only when writing the first frame. */
  float fmin,
  /** maximum pixel value in all frames that will be written;
      used only when writing the first frame. */
  float fmax,
  /** Verbose level; if zero, then nothing is printed to stderr or stdout */
  int verbose
) {
  IMG test_img;
  NIFTI_DSR dsr;
  int ret=0, fileis=0, voxNr;
  char imgfile[FILENAME_MAX], hdrfile[FILENAME_MAX], siffile[FILENAME_MAX];
  char tmp[256];
  FILE *fp;
  float *fdata=NULL, *fptr;
  int xi, yi, zi;


  if(verbose>0) {
    printf("\nimgWriteNiftiFrame(%s, %d, *img, %d, %g, %g, ...)\n",
           dbname, frame_to_write, frame_index, fmin, fmax);
    fflush(stdout);
  }
  /*
   *  Check the input 
   */
  if(dbname==NULL || !dbname[0]) {
    imgSetStatus(img, STATUS_FAULT);
    if(verbose>0) fprintf(stderr, "Error: invalid filename\n");
    return(STATUS_FAULT);
  }
  if(img==NULL || img->status!=IMG_STATUS_OCCUPIED) {
    if(img!=NULL) imgSetStatus(img, STATUS_FAULT);
    if(verbose>0) fprintf(stderr, "Error: invalid IMG argument\n");
    return(STATUS_FAULT);
  }
  if(frame_index<0 || frame_index>img->dimt-1 || frame_to_write<0) {
    imgSetStatus(img, STATUS_FAULT);
    if(verbose>0) fprintf(stderr, "Error: invalid frame settings\n");
    return STATUS_FAULT;
  }
  if(img->_fileFormat!=IMG_NIFTI_1D && img->_fileFormat!=IMG_NIFTI_1S) {
    imgSetStatus(img, STATUS_FAULT);
    if(verbose>0) fprintf(stderr, "Error: invalid file format setting\n");
    return STATUS_FAULT;
  }

  /*
   *  If NIfTI does not exist, then create it with new header,
   *  and if it does exist, then read and check header information.
   *  Create or edit header to contain correct frame nr.
   *  Determine the global scaling factors.   
   */
  imgInit(&test_img);
  fileis=niftiExists(dbname, hdrfile, imgfile, siffile, &dsr, verbose-2, NULL);
  if(fileis==0) { // not existing
    if(verbose>1) {
      printf("  writing 1st frame to a new file\n"); fflush(stdout);}

    /* Create database filenames */
    niftiCreateFNames(dbname, hdrfile, imgfile, siffile, img->_fileFormat);

    /* Set main header */
    ret=imgSetNiftiHeader(img, dbname, &dsr, fmin, fmax, verbose-1);
    if(ret!=0) {
      if(verbose>0) fprintf(stderr, "Error: cannot read NIfTI header\n");
      return STATUS_INVALIDHEADER;
    }
    if(frame_to_write==0) frame_to_write=1;
    dsr.h.dim[4]=1;

    /* Write NIfTI header */
    ret=niftiWriteHeader(hdrfile, &dsr, verbose-1, tmp);
    if(ret!=0) {
      if(verbose>0) fprintf(stderr, "Error in niftiWriteHeader(): %s\n", tmp);
      return STATUS_CANTWRITEHEADERFILE;
    }

    /* Remove dual format datafile if necessary */
    if(img->_fileFormat==IMG_NIFTI_1D)
      if(access(imgfile, 0) != -1) {
        if(verbose>0) printf("  removing %s\n", imgfile);
        remove(imgfile);
      }

  } else { /* NIfTI exists */
    if(verbose>1) printf("  adding frame to an existing file\n");

    /* Read header information for checking */
    ret=imgGetNiftiHeader(&test_img, &dsr, verbose-1);
    if(ret!=0) {
      if(verbose>0) fprintf(stderr, "Error: cannot read NIfTI header\n");
      imgEmpty(&test_img); return ret;
    }
    /* Check that file format is the same */
    if(img->_fileFormat!=test_img._fileFormat) {
      if(verbose>0) {
        fprintf(stderr, "Error: different file format\n");
        printf("  new._fileFormat:=%d\n", img->_fileFormat);
        printf("  prev._fileFormat:=%d\n", test_img._fileFormat);
      }
      imgEmpty(&test_img); return STATUS_WRONGFILETYPE;
    }
    /* Check also data type, if available */
    if(img->_dataType>0 && test_img._dataType>0 &&
       img->_dataType!=test_img._dataType)
    {
      if(verbose>0) {
        fprintf(stderr, "Error: different datatype\n");
        printf("  new._dataType:=%d\n", img->_dataType);
        printf("  prev._dataType:=%d\n", test_img._dataType);
      }
      imgEmpty(&test_img); return STATUS_WRONGFILETYPE;
    }
    /* Check that matrix sizes are the same */
    if(img->dimz!=test_img.dimz || img->dimx!=test_img.dimx ||
       img->dimy!=test_img.dimy) {
      if(verbose>0) fprintf(stderr, "Error: different matrix size\n");
      imgEmpty(&test_img); return STATUS_VARMATSIZE;}
    imgEmpty(&test_img);

    /* Set new frame number in NIfTI header */
    if(frame_to_write==0) frame_to_write=dsr.h.dim[4]+1;
    if(dsr.h.dim[4]<frame_to_write) {
      if(dsr.h.dim[4]+1<frame_to_write) {
        if(verbose>0) fprintf(stderr, "Error: missing matrix\n");
        return STATUS_MISSINGMATRIX;
      }
      dsr.h.dim[4]=frame_to_write;
    }
    /* and save the updated header */
    if((ret=niftiWriteHeader(hdrfile, &dsr, verbose-1, tmp))!=0) {
      if(verbose>0) fprintf(stderr, "Error: %s.\n", tmp);
      return STATUS_NOWRITEPERM;
    }
  }
  if(verbose>2) {
    printf("frame_to_write := %d\n", frame_to_write);
    printf("vox_offset := %d\n", (int)dsr.h.vox_offset);
    printf("hdrfile := %s\n", hdrfile);
    printf("imgfile := %s\n", imgfile);
    printf("siffile := %s\n", siffile);
    printf("magic := %s\n", dsr.h.magic);
  }

  /* Open voxel data file, not removing possible old contents like header */
  if(img->_fileFormat==IMG_NIFTI_1D && frame_to_write==1)
    fp=fopen(imgfile, "wb");
  else
    fp=fopen(imgfile, "r+b");
  if(fp==NULL) {
    if(verbose>0) fprintf(stderr, "Error: cannot open %s for write.\n",imgfile);
    return STATUS_CANTWRITEIMGFILE;
  }

  /* Move file pointer to the place of current frame */
  voxNr=img->dimz*img->dimy*img->dimx;
  if(fseek(fp, (int)dsr.h.vox_offset+(frame_to_write-1)*voxNr*sizeof(float),
           SEEK_SET)!=0)
  {
    if(verbose>0) fprintf(stderr, "Error: invalid file write position.\n");
    fclose(fp); return STATUS_MISSINGMATRIX;
  }

  /* Allocate memory for matrix float data (one plane) */
  fdata=(float*)calloc(voxNr, sizeof(float));
  if(fdata==NULL) {
    if(verbose>0) fprintf(stderr, "Error: out of memory.\n");
    fclose(fp); return STATUS_NOMEMORY;
  }

  /* Write voxel values as floats */
  for(zi=0, fptr=fdata; zi<img->dimz; zi++)
    for(yi=0; yi<img->dimy; yi++)
      for(xi=0; xi<img->dimx; xi++, fptr++)
        *fptr=img->m[zi][yi][xi][frame_index];
  fptr=fdata;
  if(fwrite(fptr, sizeof(float), voxNr, fp) != (unsigned int)voxNr) {
    if(verbose>0) fprintf(stderr, "Error: disk full or no write permission.\n");
    free(fdata); fclose(fp); return STATUS_CANTWRITEIMGFILE;
  }
  free(fdata);
  fclose(fp);

  return STATUS_OK;
}
/*****************************************************************************/

/*****************************************************************************/
/** Write NIfTI-1 image.

   Nifti database name can be given with path, however, path is not
   created in here.
   IMG field _fileFormat determines whether NIfTI is written in single
   file format (*.nii) or dual file format (*.hdr and *.img).
   Optionally SIF file with .sif extension is saved to store frame times.
 
   @return 0 if ok, and otherwise IMG status code (<>0); sets IMG->statmsg
    in case of an error.
   @sa imgInit, imgReadNiftiFirstFrame, imgReadNiftiHeader, imgReadNifti
 */
int imgWriteNifti(
  /** Nifti database name with path, with or without extension */
  const char *dbname,
  /** Pointer to IMG structure containing the image data to be written */
  IMG *img,
  /** Save (1) or do not save (0) SIF */
  int save_sif,
  /** Verbose level; if zero, then nothing is printed to stderr or stdout */
  int verbose
) {
  int ret, fi;
  char imgfile[FILENAME_MAX], hdrfile[FILENAME_MAX], siffile[FILENAME_MAX];
  float fmin, fmax;
  SIF sif;

  if(verbose>0) {
    printf("imgWriteNifti(%s, *img, %d, ...)\n", dbname, save_sif);
    fflush(stdout);
  }

  /*
   *  Check the input 
   */
  if(dbname==NULL || !dbname[0]) {
    imgSetStatus(img, STATUS_FAULT);
    if(verbose>0) fprintf(stderr, "Error: invalid filename\n");
    return(STATUS_FAULT);
  }
  if(img==NULL || img->status!=IMG_STATUS_OCCUPIED) {
    if(img!=NULL) imgSetStatus(img, STATUS_FAULT);
    if(verbose>0) fprintf(stderr, "Error: invalid IMG argument\n");
    return(STATUS_FAULT);
  }
  if(img->_fileFormat!=IMG_NIFTI_1D && img->_fileFormat!=IMG_NIFTI_1S) {
    imgSetStatus(img, STATUS_FAULT);
    if(verbose>0) fprintf(stderr, "Error: invalid file format setting\n");
    return STATUS_FAULT;
  }

  /* Create the NIfTI filename(s) */
  ret=niftiCreateFNames(dbname, hdrfile, imgfile, siffile, img->_fileFormat);
  if(ret!=0) {
    if(verbose>0) fprintf(stderr, "  Error: invalid NIfTI name %s\n", dbname);
    imgSetStatus(img, STATUS_FAULT);
    return STATUS_FAULT;
  }
  /* Delete previous NIfTI */
  ret=niftiRemove(dbname, 0, verbose-1);
  if(ret!=0) {
    if(verbose>0) fprintf(stderr, "  Error: cannot delete previous NIfTI.\n");
    imgSetStatus(img, STATUS_CANNOTERASE);
    return STATUS_CANNOTERASE;
  }

  /*
   *  Get the global min and max pixel values;
   *  those are currently only saved in header, but may be later
   *  needed for scaling pixels to short ints
   */
  if(verbose>1) {fprintf(stdout, "  searching min and max\n"); fflush(stdout);}
  ret=imgMinMax(img, &fmin, &fmax);
  if(ret) {
    if(verbose>0) fprintf(stderr, "  Error: %s\n", imgStatus(ret));
    imgSetStatus(img, STATUS_NOIMGDATA);
    return STATUS_NOIMGDATA;
  }
  if(verbose>1) {
    printf("    global_min := %g\n    global_max := %g\n", fmin, fmax);
    fflush(stdout);
  }
  /*
   *  Write the image frames
   */
  for(fi=0, ret=0; fi<img->dimt; fi++) {
    ret=imgWriteNiftiFrame(dbname, fi+1, img, fi, fmin, fmax, verbose-2);
    if(ret!=STATUS_OK) break;
    if(verbose>4) {printf("    frame written.\n"); fflush(stdout);}
  } // next frame
  //printf("ret := %d\n", ret);
  if(ret!=STATUS_OK) {
    niftiRemove(dbname, img->_fileFormat, verbose-3);
    if(verbose>0) fprintf(stderr, "Error: %s.\n", imgStatus(ret));
    return ret;
  }

  /* If SIF is not needed, then thats it */
  if(save_sif==0) {
    imgSetStatus(img, STATUS_OK);
    return 0;
  }

  /* Copy contents from IMG to SIF and save SIF */
  sifInit(&sif);
  /* Try to read existing SIF */
  ret=sifRead(siffile, &sif);
  if(ret==0) { // SIF could be read
    if(sif.frameNr==img->dimt) {
      /* If size matches, then update the contents, but keep counts, in case
         previous SIF comes with actual count info from scanner */
      ret=img2sif(img, &sif, 1, 1, 0, verbose-3);
    } else {
      /* otherwise create SIF contents */
      ret=img2sif(img, &sif, 1, 1, 2, verbose-3);
    }
  } else {
    /* otherwise create SIF contents */
    ret=img2sif(img, &sif, 1, 1, 2, verbose-3);
  }
  if(ret!=0) {
    if(verbose>0) fprintf(stderr, "  Error: cannot create SIF contents.\n");
    imgSetStatus(img, STATUS_CANNOTWRITE);
    sifEmpty(&sif);
    return STATUS_CANNOTWRITE;
  }
  /* Write SIF */
  ret=sifWrite(&sif, siffile);  
  if(ret!=0) {
    if(verbose>0) fprintf(stderr, "  Error: cannot write %s\n", siffile);
    imgSetStatus(img, STATUS_CANNOTWRITE);
    sifEmpty(&sif);
    return STATUS_CANNOTWRITE;
  }
  sifEmpty(&sif);

  imgSetStatus(img, STATUS_OK);
  return 0;
}
/*****************************************************************************/

/*****************************************************************************/
