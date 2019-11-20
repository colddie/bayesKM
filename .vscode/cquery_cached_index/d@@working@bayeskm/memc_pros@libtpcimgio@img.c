/// @file img.c
/// @author Vesa Oikonen
/// @brief Basic tools for working with IMG data struct.
///
/******************************************************************************/
#include "libtpcimgio.h"
/******************************************************************************/
/** Status (error) messages from image processing; see enum img_status_msg */
static const char *imgmsg[] = {
  /*  0, STATUS_OK                   */ "ok",
  /*  1, STATUS_FAULT                */ "fault in calling routine",
  /*  2, STATUS_NOMEMORY             */ "out of memory",
  /*  3, STATUS_NOFILE               */ "cannot open file",
  /*  4, STATUS_UNKNOWNFORMAT        */ "unknown file format",
  /*  5, STATUS_UNSUPPORTED          */ "unsupported file type",
  /*  6, STATUS_MISSINGMATRIX        */ "missing matrix/matrices",
  /*  7, STATUS_NOWRITEPERM          */ "no permission to write",
  /*  8, STATUS_DISKFULL             */ "disk full?",
  /*  9, STATUS_NOMATLIST            */ "cannot read matrix list",
  /* 10, STATUS_INVALIDMATLIST       */ "invalid matrix list",
  /* 11, STATUS_VARMATSIZE           */ "variable matrix size",
  /* 12, STATUS_NOMAINHEADER         */ "cannot read mainheader",
  /* 13, STATUS_NOSUBHEADER          */ "cannot read subheader",
  /* 14, STATUS_NOMATRIX             */ "cannot read matrix",
  /* 15, STATUS_UNSUPPORTEDAXIALCOMP */ "axial compression is not supported",
  /* 16, STATUS_NOIMGDATAFILE        */ "image datafile does not exist",
  /* 17, STATUS_NOHEADERFILE         */ "header file does not exist",
  /* 18, STATUS_INVALIDHEADER        */ "invalid header contents",
  /* 19, STATUS_NOIMGDATA            */ "cannot read image data",
  /* 20, STATUS_NOSIFDATA            */ "cannot read sif data",
  /* 21, STATUS_WRONGSIFDATA         */ "wrong sif data",
  /* 22, STATUS_CANTWRITEIMGFILE     */ "cannot write image datafile",
  /* 23, STATUS_CANTWRITEHEADERFILE  */ "cannot write header file",
  /* 24, STATUS_WRONGFILETYPE        */ "wrong file type",
  /* 25, STATUS_CANNOTERASE          */ "cannot erase file",
  /* 26, STATUS_CANNOTREAD           */ "cannot read data",
  /* 27, STATUS_CANNOTWRITE          */ "cannot write data",
  /* 28, STATUS_UNSUPPORTEDPOLARMAP  */ "polar map is not supported",
  /* 29, STATUS_INVALIDPOLARMAP      */ "invalid polar map",
  0
};
#if 0
/** Status (error) messages from image processing; see enum img_status_msg */
static const char *_imgStatusMessage[] = 
{
  /*  0, IMG_ERR_OK      */ "ok",
  /*  1, IMG_ERR_CALLING */ "fault in calling routine",
  /*  2, IMG_ERR_OOM     */ "out of memory"
};
#endif
/******************************************************************************/

/******************************************************************************/
/** Call this once before any use of IMG data.

   @sa imgAllocate, imgEmpty 
 */
void imgInit(
  /** Pointer to IMG struct. */
  IMG *image
) {
  int i=0;

  if(IMG_TEST) printf("imgInit()\n");
  if(image==NULL) return;
  memset(image, 0, sizeof(IMG));
  /*if(image->status!=IMG_STATUS_UNINITIALIZED) return;*/
  image->status=IMG_STATUS_INITIALIZED;
  imgSetStatus(image, STATUS_OK);
  image->type=0;
  image->unit=0;
  image->calibrationFactor=0.0;
  image->zoom=0.0;
  image->radiopharmaceutical[0]=(char)0;
  image->isotopeHalflife=0.0;
  image->decayCorrection=(char)IMG_DC_UNKNOWN;
  image->branchingFraction=0.0;
  image->unit=0;
  image->scanStart=0;
  image->orientation=0;
  image->axialFOV=image->transaxialFOV=image->sampleDistance=0.0;
  image->studyNr[0]=image->patientName[0]=(char)0;
  image->sizex=image->sizey=image->sizez=0;
  image->_dataType=0;
  image->_fileFormat=0;
  image->scanner=0;
  image->modality=0;
  for(i=0; i<2; i++) image->xform[i]=NIFTI_XFORM_UNKNOWN;
  for(i=0; i<18; i++) image->quatern[i]=0.0;
  for(i=0; i<12; i++) image->mt[i]=0.0;
  iftInit(&image->ift);
  image->polarmap_num_rings=0;
  for(i=0; i<MAX_POLARMAP_NUM_RINGS; i++) {
    image->polarmap_sectors_per_ring[i]=0;
    image->polarmap_ring_position[i]=0.0;
    image->polarmap_ring_angle[i]=0;
  }
  image->polarmap_start_angle=0;
  image->dimt=image->dimx=image->dimy=image->dimz=0;
  image->gapx=image->gapy=image->gapz=0.0;
  image->resolutionx=image->resolutiony=image->resolutionz=0.0;
  image->m=(float****)NULL;
  image->_header=(float*)NULL;
  image->pixel=(float*)NULL;
  image->column=(float**)NULL;
  image->row=(float***)NULL;
  image->plane=(float****)NULL;
  image->planeNumber=(int*)NULL;
  image->start=image->end=image->mid=(float*)NULL;
  image->isWeight=0;
  image->weight=image->sd=image->prompts=image->randoms=(float*)NULL;
  image->decayCorrFactor=(float*)NULL;
  image->errstatus=STATUS_OK;
}
/******************************************************************************/

/******************************************************************************/
/** Free memory that is allocated for IMG.
    @sa imgInit, imgAllocate, imgRead
 */
void imgEmpty(
  /** Pointer to IMG struct */
  IMG *image
) {
  int i=0;

  if(IMG_TEST) printf("imgEmpty()\n");
  if(image==NULL || image->status<IMG_STATUS_OCCUPIED) return;
  /* Free up memory */
  if(image->_pxl!=NULL) free(image->_pxl);
  if(image->_col!=NULL) free(image->_col);
  if(image->_row!=NULL) free(image->_row);
  if(image->_pln!=NULL) free(image->_pln);
  if(image->dimz>0) {free(image->planeNumber);}
  if(image->dimt>0) free(image->_header);
  /* Set variables */
  imgSetStatus(image, STATUS_OK);
  image->type=0;
  image->unit=0;
  image->calibrationFactor=0;
  image->zoom=0.0;
  image->radiopharmaceutical[0]=(char)0;
  image->isotopeHalflife=0.0;
  image->decayCorrection=(char)IMG_DC_UNKNOWN;
  image->branchingFraction=0.0;
  image->unit=0;
  image->scanStart=0;
  image->orientation=0;
  image->axialFOV=image->transaxialFOV=image->sampleDistance=0.0;
  image->studyNr[0]=image->patientName[0]=image->patientID[0]=(char)0;
  image->userProcessCode[0]=image->studyDescription[0]=(char)0;
  image->sizex=image->sizey=image->sizez=0;
  image->gapx=image->gapy=image->gapz=0.0;
  image->resolutionx=image->resolutiony=image->resolutionz=0.0;
  image->_dataType=0;
  image->_fileFormat=0;
  image->scanner=0;
  image->modality=0;
  for(i=0; i<2; i++) image->xform[i]=NIFTI_XFORM_UNKNOWN;
  for(i=0; i<18; i++) image->quatern[i]=0.0;
  for(i=0; i<12; i++) image->mt[i]=0.0;
  iftEmpty(&image->ift);
  image->polarmap_num_rings=0;
  for(i=0; i<MAX_POLARMAP_NUM_RINGS; i++) {
    image->polarmap_sectors_per_ring[i]=0;
    image->polarmap_ring_position[i]=0.0;
    image->polarmap_ring_angle[i]=0;
  }
  image->polarmap_start_angle=0;
  image->dimt=image->dimx=image->dimy=image->dimz=0;
  image->m=(float****)NULL;
  image->_header=(float*)NULL;
  image->pixel=(float*)NULL;
  image->column=(float**)NULL;
  image->row=(float***)NULL;
  image->plane=(float****)NULL;
  image->planeNumber=(int*)NULL;
  image->start=image->end=image->mid=(float*)NULL;
  image->isWeight=0;
  image->weight=image->sd=(float*)NULL;
  image->decayCorrFactor=(float*)NULL;
  /* Set status */
  image->status=IMG_STATUS_INITIALIZED;
  image->errstatus=STATUS_OK;
}
/******************************************************************************/

/******************************************************************************/
/** Allocates memory for img data. Old contents are not saved.
 
   @sa imgAllocateWithHeader, imgDup, imgEmpty, imgInit
   @return 0 if ok, 1 image is not initialized,
   2 invalid input dimension(s), 3 failed to allocate header, 
   4 - 8 failed to allocate image data
 */
int imgAllocate(
  /** Pointer to initialized image struct; old contents are deleted. */
  IMG *image,
  /** Nr of image planes (dim z) to allocate. */
  int planes,
  /** Nr of image rows (dim y) to allocate. */
  int rows,
  /** Nr of image columns (dim x) to allocate. */
  int columns,
  /** Nr of image time frames (dim t) to allocate. */
  int frames
) {
  unsigned short int zi, ri, ci;
  float ***rptr, **cptr, *pptr;

  if(IMG_TEST) printf("imgAllocate(*image, %d, %d, %d, %d)\n",
    planes, rows, columns, frames);
  /* Check arguments */
  if(image==NULL) return(1);
  imgSetStatus(image, STATUS_FAULT);
  if(image->status==IMG_STATUS_UNINITIALIZED) return(1);
  if(planes<1 || rows<1 || columns<1 || frames<1) return(2);
  if(image->status>=IMG_STATUS_OCCUPIED) imgEmpty(image);
  /* Allocate memory for header data */
  imgSetStatus(image, STATUS_NOMEMORY);
  image->_header=(float*)calloc(8*frames, sizeof(float));
  if(image->_header==NULL) return(3);
  image->planeNumber=(int*)calloc(planes, sizeof(int));
  if(image->planeNumber==NULL) {free(image->_header); return(4);}
  /* Allocate memory for image data */
  image->_pln=(float****)calloc(planes, sizeof(float***));
  if(image->_pln==NULL) {
    free(image->_header); free(image->planeNumber);
    return(5);
  }
  image->_row=(float***)calloc(planes*rows, sizeof(float**));
  if(image->_row==NULL) {
    free(image->_header); free(image->planeNumber);
    free(image->_pln); return(6);
  }
  image->_col=(float**)calloc(planes*rows*columns, sizeof(float*));
  if(image->_col==NULL) {
    free(image->_header); free(image->planeNumber);
    free(image->_pln); free(image->_row); return(7);
  }
  image->_pxl=(float*)calloc(planes*rows*columns*frames, sizeof(float));
  if(image->_pxl==NULL) {
    free(image->_header); free(image->planeNumber);
    free(image->_pln); free(image->_row); free(image->_col); return(8);
  }
  /* Set data pointers */
  rptr=image->_row; cptr=image->_col; pptr=image->_pxl;
  for(zi=0; zi<planes; zi++) {
    image->_pln[zi]=rptr;
    for(ri=0; ri<rows; ri++) {
      *rptr++=cptr;
      for(ci=0; ci<columns; ci++) {
        *cptr++=pptr; pptr+=frames;
      }
    }
  }
  image->m=image->_pln;
  image->plane=image->_pln;
  image->column=image->_col;
  image->row=image->_row;
  image->pixel=image->_pxl;
  /* Set header pointers */
  image->start=          image->_header+0*frames;
  image->end=            image->_header+1*frames;
  image->mid=            image->_header+2*frames;
  image->weight=         image->_header+3*frames;
  image->sd=             image->_header+4*frames;
  image->prompts=        image->_header+5*frames;
  image->randoms=        image->_header+6*frames;
  image->decayCorrFactor=image->_header+7*frames;
  /* Ok */
  image->dimz=planes; image->dimy=rows; image->dimx=columns; image->dimt=frames;
  imgSetStatus(image, STATUS_OK);
  image->status=IMG_STATUS_OCCUPIED;
  return(0);
}
/******************************************************************************/

/******************************************************************************/
/** This functions just combines imgAllocate() and imgCopyhdr().
    @sa imgAllocate, imgCopyhdr, imgDup
    @return Returns 0 if successful, otherwise returns <>0.
 */
int imgAllocateWithHeader(
  /** Pointer to IMG struct which will be allocated here. */
  IMG *image,
  /** Image matrix dimensions; z. */
  int planes,
  /** Image matrix dimensions; y. */
  int rows,
  /** Image matrix dimensions; x. */
  int columns,
  /** Image matrix dimensions; t. */
  int frames,
  /** Pointer to IMG struct where header contents will be copied from. */
  IMG *image_from
) {
  int ret;
  ret=imgAllocate(image, planes, rows, columns, frames); if(ret) return ret;
  ret=imgCopyhdr(image_from, image); return ret;
}
/******************************************************************************/

/******************************************************************************/
/** Duplicate IMG struct. Any existing contents in img2 will be deleted.
    @sa imgAllocate, imgAllocateWithHeader, imgExtractRange
    @return Returns 0 if successful.
 */
int imgDup(
  /** Pointer to IMG struct which will be duplicated */
  IMG *img1,
  /** Pointer to initiated IMG struct into which the duplicate will be made */
  IMG *img2
) {
  int i, n, ret;

  if(img1==NULL || img2==NULL) return(1);
  /* Delete any old contents */
  imgEmpty(img2);
  /* Allocate memory and copy headers */
  ret=imgAllocateWithHeader(img2, img1->dimz, img1->dimy, img1->dimx, img1->dimt, img1);
  if(ret!=0) return(10+ret);
  /* Copy voxel contents */
  n=img1->dimx*img1->dimy*img1->dimz*img1->dimt;
  for(i=0; i<n; i++) img2->pixel[i]=img1->pixel[i];
  return 0;
}
/******************************************************************************/

/******************************************************************************/
/*!
 * Return pointer to string describing the image error status message
 *
 * @param status_index index of img_status_string
 * @return pointer to string
 */
char *imgStatus(int status_index) {
  int n=0;
  while(imgmsg[n]!=0) n++;
  if(status_index<0 || status_index>n-1) return((char*)imgmsg[STATUS_FAULT]);
  else return((char*)imgmsg[status_index]);
}
/******************************************************************************/

/******************************************************************************/
/*!
 * Set the IMG image error status message pointer (statmsg) and index
 *
 * @param img Pointer to IMG struct where statmsg will be set
 * @param status_index Index of img_status_string
 */
void imgSetStatus(IMG *img, int status_index)
{
  int n=0;
  if(img==NULL) return;
  while(imgmsg[n]!=0) n++;
  if(status_index<0 || status_index>n-1) img->errstatus=STATUS_FAULT;
  else img->errstatus=status_index;
  img->statmsg=imgmsg[img->errstatus];
}
/******************************************************************************/

/******************************************************************************/
/** Prints img information to stdout; mainly for testing purposes.
 */
void imgInfo(
  /** Pointer to IMG data. */
  IMG *image
) {
  int i;
  char buf[64];

  if(IMG_TEST) printf("imgInfo()\n");
  if(image==NULL) {
    fprintf(stdout, "image := NULL\n"); return;
  } else if(image->status<=IMG_STATUS_UNINITIALIZED) {
    fprintf(stdout, "image_status := not initialized\n"); return;
  } else if(image->status==IMG_STATUS_INITIALIZED) {
    fprintf(stdout, "image_status := initialized but empty\n"); /* return; */
  } else if(image->status==IMG_STATUS_ERROR) {
    fprintf(stdout, "image_status := error\n");
  }
  fprintf(stdout, "image_error_status := %s\n", image->statmsg);
  fprintf(stdout, "image_type := %d\n", image->type);
  fprintf(stdout, "saved_data_type := %d\n", image->_dataType);
  fprintf(stdout, "file_format := %d\n", image->_fileFormat);
  fprintf(stdout, "scanner := %d\n", image->scanner);
  fprintf(stdout, "modality := %d\n", image->modality);

  fprintf(stdout, "qform := %d\n", image->xform[0]);
  fprintf(stdout, "sform := %d\n", image->xform[1]);
  fprintf(stdout, "quatern_b := %g\n", image->quatern[0]);
  fprintf(stdout, "quatern_c := %g\n", image->quatern[1]);
  fprintf(stdout, "quatern_d := %g\n", image->quatern[2]);
  fprintf(stdout, "quatern_x_shift := %g\n", image->quatern[3]);
  fprintf(stdout, "quatern_y_shift := %g\n", image->quatern[4]);
  fprintf(stdout, "quatern_z_shift := %g\n", image->quatern[5]);
  for(i=0; i<4; i++)
    fprintf(stdout, "srow_x[%d] := %g\n", 1+i, image->quatern[6+i]);
  for(i=0; i<4; i++)
    fprintf(stdout, "srow_y[%d] := %g\n", 1+i, image->quatern[10+i]);
  for(i=0; i<4; i++)
    fprintf(stdout, "srow_z[%d] := %g\n", 1+i, image->quatern[14+i]);
  for(i=0; i<12; i++)
    fprintf(stdout, "matrix_transformation[%d] := %g\n", 1+i, image->mt[i]);

  fprintf(stdout, "ift.keyNr := %d\n", image->ift.keyNr);
  fprintf(stdout, "identification_code := %.*s\n",
                   MAX_STUDYNR_LEN, image->studyNr);
  fprintf(stdout, "data_unit := %s\n", imgUnit((int)image->unit));
  fprintf(stdout, "image_zoom := %g\n", image->zoom);
  fprintf(stdout, "radiopharmaceutical := %.32s\n", image->radiopharmaceutical);
  fprintf(stdout, "isotope_halflife := %e [sec]\n", image->isotopeHalflife);
  fprintf(stdout, "branching_fraction := %f\n", image->branchingFraction);
  fprintf(stdout, "calibration_factor := %e\n", image->calibrationFactor);
  if(!ctime_r_int(&image->scanStart, buf)) strcpy(buf, "1900-01-01 00:00:00");
  fprintf(stdout, "scan_start_time := %s\n", buf);
  fprintf(stdout, "patient_name := %s\n", image->patientName);
  fprintf(stdout, "patient_id := %s\n", image->patientID);
  fprintf(stdout, "patient_orientation := %d\n", image->orientation);
  fprintf(stdout, "FOV_axial := %g [mm]\n", image->axialFOV);
  fprintf(stdout, "FOV_transaxial := %g [mm]\n", image->transaxialFOV);
  fprintf(stdout, "sample_distance := %g [mm]\n", image->sampleDistance);
  fprintf(stdout, "pixel_size_x := %g [mm]\n", image->sizex);
  fprintf(stdout, "pixel_size_y := %g [mm]\n", image->sizey);
  fprintf(stdout, "pixel_size_z := %g [mm]\n", image->sizez);
  fprintf(stdout, "dimension_x := %d\n", image->dimx);
  fprintf(stdout, "dimension_y := %d\n", image->dimy);
  fprintf(stdout, "dimension_z := %d\n", image->dimz);
  fprintf(stdout, "dimension_t := %d\n", image->dimt);
  /* Polar map */
  fprintf(stdout, "polarmap_num_rings := %d\n", image->polarmap_num_rings);
  if(image->polarmap_num_rings>0) {
    fprintf(stdout, "polarmap_sectors_per_ring :=");
    for(i=0; i<image->polarmap_num_rings; i++)
      fprintf(stdout, " %d", image->polarmap_sectors_per_ring[i]);
    fprintf(stdout, "\n");
    fprintf(stdout, "polarmap_ring_position :=");
    for(i=0; i<image->polarmap_num_rings; i++)
      fprintf(stdout, " %g", image->polarmap_ring_position[i]);
    fprintf(stdout, "\n");
    fprintf(stdout, "polarmap_ring_angle :=");
    for(i=0; i<image->polarmap_num_rings; i++)
      fprintf(stdout, " %d", image->polarmap_ring_angle[i]);
    fprintf(stdout, "\n");
    fprintf(stdout, "polarmap_start_angle := %d\n", image->polarmap_start_angle);
  }
  /* Check if the rest is available */
  if(image->status==IMG_STATUS_INITIALIZED) return;

  fprintf(stdout, "actual_plane_numbers := %d", image->planeNumber[0]);
  for(i=1; i<image->dimz; i++) fprintf(stdout, " %d", image->planeNumber[i]);
  fprintf(stdout, "\n");
  fprintf(stdout, "Frame times (sec):\n");
  for(i=0; i<image->dimt; i++) fprintf(stdout, "  %e %e %e\n",
    image->start[i], image->end[i], image->mid[i]);
  if(image->isWeight) fprintf(stdout, "Frames are weighted.\n");
  else fprintf(stdout, "Frames are not weighted.\n");
  if(image->decayCorrection==IMG_DC_CORRECTED) {
    fprintf(stdout, "Decay correction factors for each frame:\n");
    for(i=0; i<image->dimt; i++)
      fprintf(stdout, "%03i  %e\n", i+1, image->decayCorrFactor[i]);
  } else
    fprintf(stdout, "Image is not decay corrected.\n");
  return;
}
/******************************************************************************/

/******************************************************************************/
/** Copy the header fields from one image struct to another.

    Does not copy memory addresses or IMG sizes.
    Frame times, decay correction factors etc are copied, when possible.
    Plane numbers are copied, when possible.

    @sa imgAllocateWithHeader, imgDup
    @return 0 if successful, 1 invalid input, 2 pointers are to the same image
 */
int imgCopyhdr(
  /** Pointer to input IMG data. */
  IMG *image1,
  /** Pointer to output IMG data. */
  IMG *image2
) {
  int i;

  if(IMG_TEST) printf("imgCopyhdr()\n");
  /* check */
  if(image1==NULL || image2==NULL) return(1);
  if(image1==image2) return(2);
  /* copy */
  image2->type=image1->type;
  image2->unit=image1->unit;
  image2->calibrationFactor=image1->calibrationFactor;
  strcpy(image2->studyNr, image1->studyNr);
  strcpy(image2->patientName, image1->patientName);
  strcpy(image2->patientID, image1->patientID);
  strcpy(image2->userProcessCode, image1->userProcessCode);
  strcpy(image2->studyDescription, image1->studyDescription);
  image2->zoom=image1->zoom;
  strcpy(image2->radiopharmaceutical, image1->radiopharmaceutical);
  image2->isotopeHalflife=image1->isotopeHalflife;
  image2->decayCorrection=image1->decayCorrection;
  image2->branchingFraction=image1->branchingFraction;
  image2->scanStart=image1->scanStart;
  image2->axialFOV=image1->axialFOV;
  image2->transaxialFOV=image1->transaxialFOV;
  image2->sampleDistance=image1->sampleDistance;
  image2->sizex=image1->sizex;
  image2->sizey=image1->sizey;
  image2->sizez=image1->sizez;
  image2->gapx=image1->gapx;
  image2->gapy=image1->gapy;
  image2->gapz=image1->gapz;
  image2->resolutionx=image1->resolutionx;
  image2->resolutiony=image1->resolutiony;
  image2->resolutionz=image1->resolutionz;
  image2->_dataType=image1->_dataType;
  image2->_fileFormat=image1->_fileFormat;
  image2->orientation=image1->orientation;
  image2->scanner=image1->scanner;
  image2->modality=image1->modality;
  for(i=0; i<2; i++) image2->xform[i]=image1->xform[i];
  for(i=0; i<18; i++) image2->quatern[i]=image1->quatern[i];
  for(i=0; i<12; i++) image2->mt[i]=image1->mt[i];
  i=iftdup(&image1->ift, &image2->ift); if(i!=0) return(8);
  image2->polarmap_num_rings=image1->polarmap_num_rings;
  for(i=0; i<MAX_POLARMAP_NUM_RINGS; i++) {
    image2->polarmap_sectors_per_ring[i]=image1->polarmap_sectors_per_ring[i];
    image2->polarmap_ring_position[i]=image1->polarmap_ring_position[i];
    image2->polarmap_ring_angle[i]=image1->polarmap_ring_angle[i];
  }
  image2->polarmap_start_angle=image1->polarmap_start_angle;
  if(image1->dimz==image2->dimz) for(i=0; i<image1->dimz; i++) {
    image2->planeNumber[i]=image1->planeNumber[i];
  }
  if(image1->dimt==image2->dimt) for(i=0; i<image1->dimt; i++) {
    image2->start[i]=image1->start[i]; image2->end[i]=image1->end[i];
    image2->mid[i]=image1->mid[i];
    image2->weight[i]=image1->weight[i]; image2->sd[i]=image1->sd[i];
    image2->prompts[i]=image1->prompts[i];
    image2->randoms[i]=image1->randoms[i];
    image2->decayCorrFactor[i]=image1->decayCorrFactor[i];
  }
  image2->isWeight=image1->isWeight;
  return(0);
}
/******************************************************************************/

/******************************************************************************/
/** Extract a smaller 4D image from inside an IMG.

   Any existing data is overwritten.
  
   @sa imgInit, imgEmpty, imgDup
   @return 0 if ok, 1 invalid input, 2 failed to allocate memory for target image
 */
int imgExtractRange(
  /** Source image structure, 'occupied' (has allocated data). */
  IMG *img1,
  /** Image range structure. */
  IMG_RANGE r,
  /** Target image structure 'initialized' (has not allocated data). */
  IMG *img2
) {
  int zi, yi, xi, fi, zj, yj, xj, fj;

  if(IMG_TEST) {
    printf("imgExtractRange(*img1, r, *img2)\n");
    printf("  z=[%d,%d] y=[%d,%d] x=[%d,%d] f=[%d,%d]\n",
      r.z1, r.z2, r.y1, r.y2, r.x1, r.x2, r.f1, r.f2);
  }
  /* Check arguments */
  if(img2==NULL) return(1); else imgSetStatus(img2, STATUS_FAULT);
  if(img1->status!=IMG_STATUS_OCCUPIED) return(1);
  if(img2->status==IMG_STATUS_UNINITIALIZED) return(1);
  if(r.z1<1 || r.z2>img1->dimz || r.z1>r.z2) return(1);
  if(r.y1<1 || r.y2>img1->dimy || r.y1>r.y2) return(1);
  if(r.x1<1 || r.x2>img1->dimx || r.x1>r.x2) return(1);
  if(r.f1<1 || r.f2>img1->dimt || r.f1>r.f2) return(1);

  /* Allocate memory unless the same size was previously allocated */
  imgSetStatus(img2, STATUS_NOMEMORY);
  zi=r.z2-r.z1+1; yi=r.y2-r.y1+1; xi=r.x2-r.x1+1; fi=r.f2-r.f1+1;
  if(img2->status>=IMG_STATUS_OCCUPIED)
    if(img2->dimz!=zi || img2->dimy!=yi || img2->dimx!=xi || img2->dimt!=fi)
      imgEmpty(img2);
  if(img2->status!=IMG_STATUS_OCCUPIED) {
    if(imgAllocate(img2, zi, yi, xi, fi)!=0) return(2);
  }

  /* Copy data */
  imgCopyhdr(img1, img2);
  for(fi=r.f1-1, fj=0; fi<r.f2; fi++, fj++) {
    img2->start[fj]=img1->start[fi];
    img2->end[fj]=img1->end[fi];
    img2->mid[fj]=img1->mid[fi];
    img2->weight[fj]=img1->weight[fi];
    img2->sd[fj]=img1->sd[fi];
    img2->prompts[fj]=img1->prompts[fi];
    img2->randoms[fj]=img1->randoms[fi];
    img2->decayCorrFactor[fj]=img1->decayCorrFactor[fi];
  }
  for(zi=r.z1-1, zj=0; zi<r.z2; zi++, zj++) {
    img2->planeNumber[zj]=img1->planeNumber[zi];
  }
  for(zi=r.z1-1, zj=0; zi<r.z2; zi++, zj++)
    for(yi=r.y1-1, yj=0; yi<r.y2; yi++, yj++)
      for(xi=r.x1-1, xj=0; xi<r.x2; xi++, xj++)
        for(fi=r.f1-1, fj=0; fi<r.f2; fi++, fj++)
          img2->m[zj][yj][xj][fj]=img1->m[zi][yi][xi][fi];

  imgSetStatus(img2, STATUS_OK);
  return(0);
}
/******************************************************************************/

/******************************************************************************/
/** Verify that IMG contains frame times.
    @sa imgExistentCounts, imgFramesCheck
    @return Returns nonzero if frame times are there, and 0 if not.
 */
int imgExistentTimes(
  /** Pointer to IMG struct. */
  IMG *img
) {
  int fi;
  if(img==NULL || img->status!=IMG_STATUS_OCCUPIED || img->dimt<1) return 0;
  for(fi=0; fi<img->dimt; fi++) if(img->end[fi]>0.00000001) return 1;
  return 0;  
}
/*****************************************************************************/

/*****************************************************************************/
/** Verify that IMG contains prompts and randoms.

    @sa imgExistentTimes
    @return Returns 0 if neither prompts or randoms can be found, 
     1 or 2 if prompts or randoms can be found, and 3 if both prompts and randoms are there.
 */
int imgExistentCounts(
  /** Pointer to IMG struct */
  IMG *img
) {
  int fi, p=0, r=0;
  float v1, v2;
  if(img==NULL || img->status!=IMG_STATUS_OCCUPIED || img->dimt<1) return 0;
  /* If just one frame, then value > 0 is fine */
  if(img->dimt==1) {
    if(img->prompts[0]>0.00000001) p=1;
    if(img->randoms[0]>0.00000001) r=2;
    return(p+r);
  }
  /* Otherwise, check also that frames have different count level */
  for(fi=1; fi<img->dimt; fi++) {
    v1=img->prompts[fi]-img->prompts[fi-1]; if(fabs(v1)>0.001) p=1;
    v2=img->randoms[fi]-img->randoms[fi-1]; if(fabs(v2)>0.001) r=2;
    if((p+r)>2) break;
  }
  return(p+r);
}
/*****************************************************************************/

/*****************************************************************************/
