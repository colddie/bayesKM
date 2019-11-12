/// @file img_e7.c
/// @author Vesa Oikonen
/// @brief ECAT 7 I/O routines for IMG data.
///
/*****************************************************************************/
#include "libtpcimgio.h"
/*****************************************************************************/

/*****************************************************************************/
/** Read ECAT 7 image, volume or 2D sinogram.
 
 @return 0 if ok, 1 invalid input, 2 image status is not 'initialized', 
  3 failed to open file fiel for reading, 4 recognize file,
  5 file type not supported, 6 invalid matrix list,
  7 invalid number of matrixes/frames, 8 variable matrix size,
  9 failed to read header, 11 failed to allocate memory for data, 
  13 failed to read data.
 */
int imgReadEcat7(
  /** Name of the input ECAT 7 file. */
  const char *fname,
  /** Data structure in which the file is read; must be initialized
      before using this function. */
  IMG *img
) {
  FILE *fp;
  int ret, m, i, fi, pi, xi, yi, frame, plane, prev_frame, prev_plane;
  int dimx, dimy, dimz, planeNr, frameNr, blkNr=0, pxlNr;
  ECAT7_mainheader main_header;
  ECAT7_imageheader image_header;
  ECAT7_2Dscanheader scan2d_header;
  ECAT7_scanheader scan_header;
  ECAT7_polmapheader polmap_header;
  ECAT7_MATRIXLIST mlist;
  ECAT7_Matval matval;
  float *fdata=NULL, *fptr;


  if(IMG_TEST) printf("imgReadEcat7(%s, *img)\n", fname);
  /* Check the arguments */
  if(fname==NULL) return(1);
  if(img==NULL || img->status!=IMG_STATUS_INITIALIZED) {
    imgSetStatus(img, STATUS_FAULT); return(2);}

  /* Open file for read */
  if((fp=fopen(fname, "rb"))==NULL) {imgSetStatus(img,STATUS_NOFILE); return(3);}
  
  /* Read main header */
  ret=ecat7ReadMainheader(fp, &main_header);
  if(ret) {fclose(fp); imgSetStatus(img, STATUS_UNKNOWNFORMAT); return(4);}
  /* Check for magic number */
  if(strncmp(main_header.magic_number, ECAT7V_MAGICNR, 7)!=0) {
    fclose(fp); imgSetStatus(img, STATUS_UNKNOWNFORMAT); return(4);
  }
  
  /* Check if file_type is supported */
  if(imgEcat7Supported(&main_header)==0) {
    fclose(fp); imgSetStatus(img, STATUS_UNSUPPORTED); return(5);
  }

  /* Read matrix list */
  ecat7InitMatlist(&mlist);
  ret=ecat7ReadMatlist(fp, &mlist, IMG_TEST-1);
  if(ret || mlist.matrixNr<1 || ecat7CheckMatlist(&mlist)) {
    fclose(fp); imgSetStatus(img, STATUS_INVALIDMATLIST); return(6);}
  ecat7SortMatlistByPlane(&mlist);
  if(IMG_TEST>2) ecat7PrintMatlist(&mlist);

  /* Calculate the number of planes, frames and gates */
  /* Check that all planes have equal nr of frames (gates) */
  /* and check that frames (gates) are consequentally numbered */
  prev_plane=plane=-1; prev_frame=frame=-1; frameNr=planeNr=0; ret=0;
  for(m=0; m<mlist.matrixNr; m++) {
    ecat7_id_to_val(mlist.matdir[m].id, &matval); plane=matval.plane;
    if(main_header.num_frames>=main_header.num_gates) frame=matval.frame;
    else frame=matval.gate;
    if(plane!=prev_plane) {frameNr=1; planeNr++;}
    else {frameNr++; if(prev_frame>0 && frame!=prev_frame+1) {ret=1; break;}}
    prev_plane=plane; prev_frame=frame;
    /* Calculate and check the size of one data matrix */
    if(m==0) {
      blkNr=mlist.matdir[m].endblk-mlist.matdir[m].strtblk;
    } else if(blkNr!=mlist.matdir[m].endblk-mlist.matdir[m].strtblk) {
      ret=2; break;
    }
  } /* next matrix */
  if(IMG_TEST>2) printf("frameNr=%d planeNr=%d\n", frameNr, planeNr);
  if(ret==1 || (frameNr*planeNr != mlist.matrixNr)) {
    fclose(fp); imgSetStatus(img, STATUS_MISSINGMATRIX);
    ecat7EmptyMatlist(&mlist); return(7);
  }
  if(ret==2) {
    fclose(fp); imgSetStatus(img, STATUS_VARMATSIZE);
    ecat7EmptyMatlist(&mlist); return(8);
  }

  /* Read the first subheader to get planeNr from volumes and to get x&y dim */
  m=0; dimz=1; imgSetStatus(img, STATUS_NOSUBHEADER);
  switch(main_header.file_type) {
    case ECAT7_IMAGE8:
    case ECAT7_IMAGE16:
    case ECAT7_VOLUME8:
    case ECAT7_VOLUME16:
      img->type=IMG_TYPE_IMAGE;
      ret=ecat7ReadImageheader(fp, mlist.matdir[m].strtblk, &image_header);
      dimx=image_header.x_dimension; dimy=image_header.y_dimension;
      if(image_header.num_dimensions>2 && image_header.z_dimension>1)
        planeNr=dimz=image_header.z_dimension;
      break;
    case ECAT7_2DSCAN:
      img->type=IMG_TYPE_RAW;
      ret=ecat7Read2DScanheader(fp, mlist.matdir[m].strtblk, &scan2d_header);
      dimx=scan2d_header.num_r_elements; dimy=scan2d_header.num_angles;
      if(scan2d_header.num_dimensions>2 && scan2d_header.num_z_elements>1)
        planeNr=dimz=scan2d_header.num_z_elements;
      break;
    case ECAT7_3DSCAN:
    case ECAT7_3DSCAN8:
    case ECAT7_3DSCANFIT:
      img->type=IMG_TYPE_RAW;
      ret=ecat7ReadScanheader(fp, mlist.matdir[m].strtblk, &scan_header);
      dimx=scan_header.num_r_elements; dimy=scan_header.num_angles;
      for(i=dimz=0; i<64; i++) dimz+=scan_header.num_z_elements[i];
      planeNr=dimz;
      /*if(scan_header.axial_compression!=0) {
        img->statmsg=imgmsg[STATUS_UNSUPPORTEDAXIALCOMP]; ret=-1;}*/
      break;
    case ECAT7_POLARMAP:
      img->type=IMG_TYPE_POLARMAP;
      ret=ecat7ReadPolmapheader(fp, mlist.matdir[m].strtblk, &polmap_header);
      planeNr=dimz=dimy=1; dimx=0;
      for(i=0; i<polmap_header.num_rings; i++)
        dimx+=polmap_header.sectors_per_ring[i];
      break;
    default: dimx=dimy=dimz=planeNr=0; ret=-1;
  }
  pxlNr=dimx*dimy;
  if(ret || pxlNr<1 || planeNr<1) {
    fclose(fp);  ecat7EmptyMatlist(&mlist); return(9);}
  imgSetStatus(img, STATUS_OK);

  /* Allocate memory for IMG data */
  ret=imgAllocate(img, planeNr, dimy, dimx, frameNr);
  if(ret) {
    fclose(fp); imgSetStatus(img, STATUS_NOMEMORY);
    ecat7EmptyMatlist(&mlist); return(11);
  }
  /* Copy information from mainheader */
  imgGetEcat7MHeader(img, &main_header);
  /* Set fileFormat */
  switch(main_header.file_type) {
    case ECAT7_IMAGE8:
    case ECAT7_IMAGE16:
      img->_fileFormat=IMG_E7_2D; break;
    case ECAT7_VOLUME8:
    case ECAT7_VOLUME16:
      img->_fileFormat=IMG_E7; break;
    case ECAT7_2DSCAN:
      img->_fileFormat=IMG_E7_2D; break;
    case ECAT7_3DSCAN:
    case ECAT7_3DSCAN8:
    case ECAT7_3DSCANFIT:
      img->_fileFormat=IMG_E7; break;
    case ECAT7_POLARMAP:
      img->_fileFormat=IMG_POLARMAP; break;
    default:
      img->_fileFormat=IMG_UNKNOWN; break;
  }

  if(dimz>1) {
    /* Read ECAT volume matrices */
    fi=0;
    for(m=0; m<mlist.matrixNr; m++) {
      /* Get matrix values */
      ecat7_id_to_val(mlist.matdir[m].id, &matval);
      /* Read subheader and data */
      if(img->type==IMG_TYPE_IMAGE)
        ret=ecat7ReadImageMatrix(fp, mlist.matdir[m].strtblk,
              mlist.matdir[m].endblk, &image_header, &fdata);
      else
        ret=ecat7ReadScanMatrix(fp, mlist.matdir[m].strtblk,
              mlist.matdir[m].endblk, &scan_header, &fdata);
      if(ret || fdata==NULL) {
        if(IMG_TEST) printf("ecat7ReadXMatrix()=%d\n%s\n", ret, ecat7errmsg);
        fclose(fp); imgSetStatus(img, STATUS_NOMATRIX);
        ecat7EmptyMatlist(&mlist); return(13);
      }
      /* Copy subheader information */
      if(img->type==IMG_TYPE_POLARMAP) {
        img->_dataType=polmap_header.data_type;
        img->start[fi]=polmap_header.frame_start_time/1000.;
        img->end[fi]=img->start[fi]+polmap_header.frame_duration/1000.;
        img->mid[fi]=0.5*(img->start[fi]+img->end[fi]);
        img->sizex=0.001*polmap_header.pixel_size;
      } else if(img->type==IMG_TYPE_IMAGE) {
        img->_dataType=image_header.data_type;
        img->start[fi]=image_header.frame_start_time/1000.;
        img->end[fi]=img->start[fi]+image_header.frame_duration/1000.;
        img->mid[fi]=0.5*(img->start[fi]+img->end[fi]);
        if(image_header.decay_corr_fctr>1.0) {
          img->decayCorrFactor[fi]=image_header.decay_corr_fctr;
          img->decayCorrection=IMG_DC_CORRECTED;
        } else {
          img->decayCorrFactor[fi]=0.0;
          img->decayCorrection=IMG_DC_UNKNOWN;
        }
        img->zoom=image_header.recon_zoom;
        img->sizex=10.*image_header.x_pixel_size;
        img->sizey=10.*image_header.y_pixel_size;
        img->sizez=10.*image_header.z_pixel_size;
        img->resolutionx=10.*image_header.x_resolution;
        img->resolutiony=10.*image_header.y_resolution;
        img->resolutionz=10.*image_header.z_resolution;
        img->xform[0]=NIFTI_XFORM_UNKNOWN; // qform
        img->xform[1]=NIFTI_XFORM_SCANNER_ANAT; // sform
        img->quatern[6]=img->sizex; img->quatern[9]=img->sizex;
        img->quatern[11]=img->sizey; img->quatern[13]=img->sizey;
        img->quatern[16]=img->sizez; img->quatern[17]=img->sizez;
        img->mt[0]=image_header.mt_1_1;
        img->mt[1]=image_header.mt_1_2;
        img->mt[2]=image_header.mt_1_3;
        img->mt[3]=image_header.mt_1_4;
        img->mt[4]=image_header.mt_2_1;
        img->mt[5]=image_header.mt_2_2;
        img->mt[6]=image_header.mt_2_3;
        img->mt[7]=image_header.mt_2_4;
        img->mt[8]=image_header.mt_3_1;
        img->mt[9]=image_header.mt_3_2;
        img->mt[10]=image_header.mt_3_3;
        img->mt[11]=image_header.mt_3_4;
      } else {
        img->_dataType=scan_header.data_type;
        img->start[fi]=scan_header.frame_start_time/1000.;
        img->end[fi]=img->start[fi]+scan_header.frame_duration/1000.;
        img->mid[fi]=0.5*(img->start[fi]+img->end[fi]);
        if(scan_header.x_resolution>0.0)
          img->sampleDistance=10.0*scan_header.x_resolution;
        else
          img->sampleDistance=10.0*main_header.bin_size;
        /* also, correct for dead-time */
        if(scan_header.deadtime_correction_factor>0.0)
          for(i=0, fptr=fdata; i<dimz*pxlNr; i++, fptr++)
            *fptr*=scan_header.deadtime_correction_factor;
        img->prompts[fi]=(float)scan_header.prompts;
        img->randoms[fi]=scan_header.delayed;
      }
      /* Copy matrix data through volume planes */
      for(pi=0; pi<dimz; pi++) {
        for(yi=0, fptr=fdata+pi*pxlNr; yi<dimy; yi++) for(xi=0; xi<dimx; xi++)
          img->m[pi][yi][xi][fi]= *fptr++;
      }
      free(fdata); fi++;
    } /* next matrix */
    /* Set plane numbers */
    for(pi=0; pi<dimz; pi++) img->planeNumber[pi]=pi+1;
  } else {
    /* Read separate matrices */
    prev_plane=plane=-1; prev_frame=frame=-1; pi=fi=-1;
    for(m=0; m<mlist.matrixNr; m++) {
      ecat7_id_to_val(mlist.matdir[m].id, &matval); plane=matval.plane;
      if(main_header.num_frames>=main_header.num_gates) frame=matval.frame;
      else frame=matval.gate;
      if(plane!=prev_plane) {fi=0; pi++;} else fi++;
      /* Read subheader and data */
      if(img->type==IMG_TYPE_POLARMAP)
        ret=ecat7ReadPolarmapMatrix(fp, mlist.matdir[m].strtblk,
              mlist.matdir[m].endblk, &polmap_header, &fdata);
      else if(img->type==IMG_TYPE_IMAGE)
        ret=ecat7ReadImageMatrix(fp, mlist.matdir[m].strtblk,
              mlist.matdir[m].endblk, &image_header, &fdata);
      else
        ret=ecat7Read2DScanMatrix(fp, mlist.matdir[m].strtblk,
              mlist.matdir[m].endblk, &scan2d_header, &fdata);
      if(ret || fdata==NULL) {
        fclose(fp); imgSetStatus(img, STATUS_NOMATRIX);
        ecat7EmptyMatlist(&mlist); return(13);
      }
      /* Copy subheader information */
      if(fi==0) img->planeNumber[pi]=plane;
      if(img->type==IMG_TYPE_POLARMAP) {
        img->_dataType=polmap_header.data_type;
        img->start[fi]=polmap_header.frame_start_time/1000.;
        img->end[fi]=img->start[fi]+polmap_header.frame_duration/1000.;
        img->mid[fi]=0.5*(img->start[fi]+img->end[fi]);
        img->sizex=0.001*polmap_header.pixel_size;
      } else if(img->type==IMG_TYPE_IMAGE) {
        img->_dataType=image_header.data_type;
        img->start[fi]=image_header.frame_start_time/1000.;
        img->end[fi]=img->start[fi]+image_header.frame_duration/1000.;
        img->mid[fi]=0.5*(img->start[fi]+img->end[fi]);
        if(image_header.decay_corr_fctr>1.0) {
          img->decayCorrFactor[fi]=image_header.decay_corr_fctr;
          img->decayCorrection=IMG_DC_CORRECTED;
        } else {
          img->decayCorrFactor[fi]=0.0;
          img->decayCorrection=IMG_DC_UNKNOWN;
        }
        img->zoom=image_header.recon_zoom;
        img->sizex=10.*image_header.x_pixel_size;
        img->sizey=10.*image_header.y_pixel_size;
        img->sizez=10.*image_header.z_pixel_size;
        img->resolutionx=10.*image_header.x_resolution;
        img->resolutiony=10.*image_header.y_resolution;
        img->resolutionz=10.*image_header.z_resolution;
        img->xform[0]=NIFTI_XFORM_UNKNOWN; // qform
        img->xform[1]=NIFTI_XFORM_SCANNER_ANAT; // sform
        img->quatern[6]=img->sizex; img->quatern[9]=img->sizex;
        img->quatern[11]=img->sizey; img->quatern[13]=img->sizey;
        img->quatern[16]=img->sizez; img->quatern[17]=img->sizez;
        img->mt[0]=image_header.mt_1_1;
        img->mt[1]=image_header.mt_1_2;
        img->mt[2]=image_header.mt_1_3;
        img->mt[3]=image_header.mt_1_4;
        img->mt[4]=image_header.mt_2_1;
        img->mt[5]=image_header.mt_2_2;
        img->mt[6]=image_header.mt_2_3;
        img->mt[7]=image_header.mt_2_4;
        img->mt[8]=image_header.mt_3_1;
        img->mt[9]=image_header.mt_3_2;
        img->mt[10]=image_header.mt_3_3;
        img->mt[11]=image_header.mt_3_4;
      } else {
        img->_dataType=scan2d_header.data_type;
        img->start[fi]=scan2d_header.frame_start_time/1000.;
        img->end[fi]=img->start[fi]+scan2d_header.frame_duration/1000.;
        img->mid[fi]=0.5*(img->start[fi]+img->end[fi]);
        if(scan_header.x_resolution>0.0)
          img->sampleDistance=10.0*scan_header.x_resolution;
        else
          img->sampleDistance=10.0*main_header.bin_size;
        /* also, correct for dead-time */
        if(scan2d_header.deadtime_correction_factor>0.0)
          for(i=0, fptr=fdata; i<pxlNr; i++, fptr++)
            *fptr*=scan2d_header.deadtime_correction_factor;
        img->prompts[fi]=(float)scan_header.prompts;
        img->randoms[fi]=scan_header.delayed;
      }
      /* Copy matrix data */
      for(yi=0, fptr=fdata; yi<dimy; yi++) for(xi=0; xi<dimx; xi++)
        img->m[pi][yi][xi][fi]= *fptr++;
      free(fdata);
      /* prepare for the next matrix */
      prev_plane=plane; prev_frame=frame;
    } /* next matrix */
  }
  fclose(fp); ecat7EmptyMatlist(&mlist);

  /* Calibrate */
  if(main_header.ecat_calibration_factor>0.0)
    for(pi=0; pi<img->dimz; pi++)
      for(yi=0; yi<img->dimy; yi++) for(xi=0; xi<img->dimx; xi++)
        for(fi=0; fi<img->dimt; fi++)
          img->m[pi][yi][xi][fi]*=main_header.ecat_calibration_factor;

  imgSetStatus(img, STATUS_OK);
  return(0);
}
/*****************************************************************************/

/*****************************************************************************/
/*!
 * Write ECAT 7 3D image volume or 3D sinogram.
 *
 * @param fname output filename
 * @param img pointer to IMG data
 * @return 0 if ok, 1 invalid input, 2 image status is not 'occupied', 
 * 3 failed to allocate memory for data, 6 failed to create file, 
 * 7 failed to write data, 8 unsupported image type
 * sets IMG->statmsg in case of error
 */
int imgWriteEcat7(const char *fname, IMG *img) {
  ECAT7_mainheader main_header;
  ECAT7_imageheader image_header;
  ECAT7_scanheader scan_header;
  FILE *fp;
  int fi, pi, xi, yi, pxlNr, matrixId, ret;
  float *fdata, *fptr;


  if(IMG_TEST) printf("imgWriteEcat7(%s, *img)\n", fname);
  if(IMG_TEST>1 && ECAT7_TEST==0) ECAT7_TEST=1;
  /* Check the arguments */
  if(fname==NULL) return(1);
  if(img==NULL || img->status!=IMG_STATUS_OCCUPIED) {
    imgSetStatus(img, STATUS_FAULT); return(2);}
  if(img->type!=IMG_TYPE_RAW && img->type!=IMG_TYPE_IMAGE) {
    imgSetStatus(img, STATUS_FAULT); return(2);}

  /* Initiate headers */
  memset(&main_header, 0, sizeof(ECAT7_mainheader));
  memset(&image_header,0, sizeof(ECAT7_imageheader));
  memset(&scan_header, 0, sizeof(ECAT7_scanheader));

  /* Set main header */
  imgSetEcat7MHeader(img, &main_header);
  main_header.bin_size=img->sampleDistance/10.0;

  /* Allocate memory for matrix float data */
  pxlNr=img->dimx*img->dimy*img->dimz;
  fdata=(float*)malloc(pxlNr*sizeof(float));
  if(fdata==NULL) {imgSetStatus(img, STATUS_NOMEMORY); return(3);}

  /* Open file, write main header and initiate matrix list */
  fp=ecat7Create(fname, &main_header);
  if(fp==NULL) {free(fdata); imgSetStatus(img, STATUS_NOWRITEPERM); return(6);}

  /* Set (most of) subheader contents */
  if(img->type==IMG_TYPE_RAW) {
    scan_header.x_resolution=img->sampleDistance/10.0;
    scan_header.num_dimensions=4;
    if(img->dimz==239) {
      scan_header.num_z_elements[0]=63;
      scan_header.num_z_elements[1]=106;
      scan_header.num_z_elements[2]=70;
    } else {
      scan_header.num_z_elements[0]=img->dimz;
    }
    scan_header.storage_order=1;
    scan_header.data_type=ECAT7_SUNI2;
    scan_header.num_r_elements=img->dimx;
    scan_header.num_angles=img->dimy;
  } else if(img->type==IMG_TYPE_IMAGE) {
    image_header.num_dimensions=3;
    image_header.z_dimension=img->dimz;
    image_header.data_type=ECAT7_SUNI2;
    image_header.x_dimension=img->dimx;
    image_header.y_dimension=img->dimy;
    image_header.recon_zoom=img->zoom;
    image_header.x_pixel_size=0.1*img->sizex;
    image_header.y_pixel_size=0.1*img->sizey;
    image_header.z_pixel_size=0.1*img->sizez;
    image_header.x_resolution=0.1*img->resolutionx;
    image_header.y_resolution=0.1*img->resolutiony;
    image_header.z_resolution=0.1*img->resolutionz;
    img->quatern[6]=img->sizex; img->quatern[9]=img->sizex;
    img->quatern[11]=img->sizey; img->quatern[13]=img->sizey;
    img->quatern[16]=img->sizez; img->quatern[17]=img->sizez;
    image_header.mt_1_1=img->mt[0];
    image_header.mt_1_2=img->mt[1];
    image_header.mt_1_3=img->mt[2];
    image_header.mt_1_4=img->mt[3];
    image_header.mt_2_1=img->mt[4];
    image_header.mt_2_2=img->mt[5];
    image_header.mt_2_3=img->mt[6];
    image_header.mt_2_4=img->mt[7];
    image_header.mt_3_1=img->mt[8];
    image_header.mt_3_2=img->mt[9];
    image_header.mt_3_3=img->mt[10];
    image_header.mt_3_4=img->mt[11];
  }

  /* Write each matrix */
  for(fi=0; fi<img->dimt; fi++) {

    /* Create new matrix id (i.e. matnum) */
    matrixId=ecat7_val_to_id(fi+1, 1, 1, 0, 0);

    /* Copy matrix pixel values to fdata */
    fptr=fdata;
    for(pi=0; pi<img->dimz; pi++)
      for(yi=0; yi<img->dimy; yi++) for(xi=0; xi<img->dimx; xi++)
        *fptr++=img->m[pi][yi][xi][fi];

    /* Write subheader and data */
    fptr=fdata;
    if(img->type==IMG_TYPE_RAW) {
      scan_header.frame_start_time=(int)temp_roundf(1000.*img->start[fi]);
      scan_header.frame_duration=
        (int)temp_roundf(1000.*(img->end[fi]-img->start[fi]));
      scan_header.prompts=temp_roundf(img->prompts[fi]);
      scan_header.delayed=temp_roundf(img->randoms[fi]);
      /*ecat7PrintScanheader(&scan_header, stdout);*/
      ret=ecat7WriteScanMatrix(fp, matrixId, &scan_header, fptr);
    } else if(img->type==IMG_TYPE_IMAGE) {
      image_header.frame_start_time=(int)temp_roundf(1000.*img->start[fi]);
      image_header.frame_duration=
        (int)temp_roundf(1000.*(img->end[fi]-img->start[fi]));
      if(img->decayCorrection==IMG_DC_CORRECTED)
        image_header.decay_corr_fctr=img->decayCorrFactor[fi];
      else
        image_header.decay_corr_fctr=0.0;
      /*ecat7PrintImageheader(&image_header, stdout);*/
      ret=ecat7WriteImageMatrix(fp, matrixId, &image_header, fptr);
    } else {
      free(fdata); fclose(fp); imgSetStatus(img, STATUS_UNSUPPORTED); return(8);
    }
    if(ret) {
      if(IMG_TEST) {printf("matrixId=%d ret=%d\n", matrixId, ret);}
      free(fdata); fclose(fp); imgSetStatus(img, STATUS_DISKFULL); return(7);
    }

  } /* next matrix */
  free(fdata); fclose(fp);

  imgSetStatus(img, STATUS_OK);
  return(0);
}
/*****************************************************************************/

/*****************************************************************************/
/*!
 * Write ECAT 7 2D image or 2D sinogram.
 *
 * @param fname output filename
 * @param img pointer to image structure
 * @return 0 if ok, 1 invalid input, 2 image status is not 'occupied', 
 * 3 failed to allocate memory for data, 6 faield to create file, 
 * 7 failed to write data, 8 image type not supported,
 * sets IMG->statmsg in case of error
 */
int imgWrite2DEcat7(const char *fname, IMG *img) {
  ECAT7_mainheader main_header;
  ECAT7_imageheader image_header;
  ECAT7_2Dscanheader scan2d_header;
  FILE *fp;
  int fi, pi, xi, yi, pxlNr, matrixId, ret;
  float *fdata, *fptr;


  if(IMG_TEST) printf("imgWrite2DEcat7(%s, *img)\n", fname);
  if(IMG_TEST>1 && ECAT7_TEST==0) ECAT7_TEST=1;
  /* Check the arguments */
  if(fname==NULL) {return(1);}
  if(img==NULL || img->status!=IMG_STATUS_OCCUPIED) {
    imgSetStatus(img, STATUS_FAULT); return(2);}

  /* Initiate headers */
  memset(&main_header, 0, sizeof(ECAT7_mainheader));
  memset(&image_header,0, sizeof(ECAT7_imageheader));
  memset(&scan2d_header, 0, sizeof(ECAT7_2Dscanheader));

  /* Set main header */
  imgSetEcat7MHeader(img, &main_header);
  main_header.bin_size=img->sampleDistance/10.0;
  if(img->type==IMG_TYPE_RAW) main_header.file_type=ECAT7_2DSCAN;
  else main_header.file_type=ECAT7_IMAGE16;
  main_header.num_planes=img->dimz;

  /* Allocate memory for matrix float data */
  pxlNr=img->dimx*img->dimy;
  fdata=(float*)malloc(pxlNr*sizeof(float));
  if(fdata==NULL) {imgSetStatus(img, STATUS_NOMEMORY); return(3);}

  /* Open file, write main header and initiate matrix list */
  fp=ecat7Create(fname, &main_header);
  if(fp==NULL) {free(fdata); imgSetStatus(img, STATUS_NOWRITEPERM); return(6);}

  /* Set (most of) subheader contents */
  if(img->type==IMG_TYPE_RAW) {
    scan2d_header.num_dimensions=2;
    scan2d_header.num_z_elements=1;
    scan2d_header.data_type=ECAT7_SUNI2;
    scan2d_header.num_r_elements=img->dimx;
    scan2d_header.num_angles=img->dimy;
  } else if(img->type==IMG_TYPE_IMAGE) {
    image_header.num_dimensions=2;
    image_header.z_dimension=1;
    image_header.data_type=ECAT7_SUNI2;
    image_header.x_dimension=img->dimx;
    image_header.y_dimension=img->dimy;
    image_header.recon_zoom=img->zoom;
    image_header.x_pixel_size=0.1*img->sizex;
    image_header.y_pixel_size=0.1*img->sizey;
    image_header.z_pixel_size=0.1*img->sizez;
    image_header.x_resolution=0.1*img->resolutionx;
    image_header.y_resolution=0.1*img->resolutiony;
    image_header.z_resolution=0.1*img->resolutionz;
    img->quatern[6]=img->sizex; img->quatern[9]=img->sizex;
    img->quatern[11]=img->sizey; img->quatern[13]=img->sizey;
    img->quatern[16]=img->sizez; img->quatern[17]=img->sizez;
    image_header.mt_1_1=img->mt[0];
    image_header.mt_1_2=img->mt[1];
    image_header.mt_1_3=img->mt[2];
    image_header.mt_1_4=img->mt[3];
    image_header.mt_2_1=img->mt[4];
    image_header.mt_2_2=img->mt[5];
    image_header.mt_2_3=img->mt[6];
    image_header.mt_2_4=img->mt[7];
    image_header.mt_3_1=img->mt[8];
    image_header.mt_3_2=img->mt[9];
    image_header.mt_3_3=img->mt[10];
    image_header.mt_3_4=img->mt[11];
  }

  /* Write each matrix */
  for(fi=0; fi<img->dimt; fi++) for(pi=0; pi<img->dimz; pi++) {

    /* Create new matrix id (i.e. matnum) */
    matrixId=ecat7_val_to_id(fi+1, img->planeNumber[pi], 1, 0, 0);

    /* Copy matrix pixel values to fdata */
    fptr=fdata;
    for(yi=0; yi<img->dimy; yi++) for(xi=0; xi<img->dimx; xi++)
      *fptr++=img->m[pi][yi][xi][fi];

    /* Write subheader and data */
    fptr=fdata;
    if(img->type==IMG_TYPE_RAW) {
      scan2d_header.frame_start_time=(int)temp_roundf(1000.*img->start[fi]);
      scan2d_header.frame_duration=
       (int)temp_roundf(1000.*(img->end[fi]-img->start[fi]));
      scan2d_header.prompts=temp_roundf(img->prompts[fi]);
      scan2d_header.delayed=temp_roundf(img->randoms[fi]);
      ret=ecat7Write2DScanMatrix(fp, matrixId, &scan2d_header, fptr);
    } else if(img->type==IMG_TYPE_IMAGE) {
      image_header.frame_start_time=(int)temp_roundf(1000.*img->start[fi]);
      image_header.frame_duration=
        (int)temp_roundf(1000.*(img->end[fi]-img->start[fi]));
      if(img->decayCorrection==IMG_DC_CORRECTED)
        image_header.decay_corr_fctr=img->decayCorrFactor[fi];
      else
        image_header.decay_corr_fctr=0.0;
      ret=ecat7WriteImageMatrix(fp, matrixId, &image_header, fptr);
    } else {
      free(fdata); fclose(fp); imgSetStatus(img, STATUS_UNSUPPORTED); return(8);
    }
    if(ret) {
      if(IMG_TEST) {printf("matrixId=%d ret=%d\n", matrixId, ret);}
      free(fdata); fclose(fp); imgSetStatus(img, STATUS_DISKFULL); return(7);
    }

  } /* next matrix */
  free(fdata); fclose(fp);

  imgSetStatus(img, STATUS_OK);
  return(0);
}
/*****************************************************************************/

/*****************************************************************************/
/*!
 * Write ECAT 7 polar map.
 *
 * @param fname output filename
 * @param img pointer to image structure
 * @return 0 if ok, 1 invalid input, 2 image status is not 'occupied',
 * 3 failed to allocate memory for data, 6 faield to create file,
 * 7 failed to write data, 8 image type not supported,
 * sets IMG->statmsg in case of error
 */
int imgWritePolarmap(const char *fname, IMG *img) {
  ECAT7_mainheader main_header;
  ECAT7_polmapheader polmap_header;
  FILE *fp;
  int fi, pi, xi, yi, pxlNr, matrixId, ret;
  float *fdata, *fptr;


  if(IMG_TEST) printf("imgWritePolarmap(%s, *img)\n", fname);
  if(IMG_TEST>1 && ECAT7_TEST==0) ECAT7_TEST=1;
  /* Check the arguments */
  if(fname==NULL) return(1);
  if(img==NULL || img->status!=IMG_STATUS_OCCUPIED) {
    imgSetStatus(img, STATUS_FAULT); return(2);}
  if(img->type!=IMG_TYPE_POLARMAP) {
    imgSetStatus(img, STATUS_FAULT); return(2);}

  /* Initiate headers */
  memset(&main_header, 0, sizeof(ECAT7_mainheader));
  memset(&polmap_header,0, sizeof(ECAT7_polmapheader));

  /* Set main header */
  imgSetEcat7MHeader(img, &main_header);
  main_header.bin_size=img->sampleDistance/10.0;

  /* Allocate memory for matrix float data */
  pxlNr=img->dimx*img->dimy*img->dimz;
  fdata=(float*)malloc(pxlNr*sizeof(float));
  if(fdata==NULL) {imgSetStatus(img, STATUS_NOMEMORY); return(3);}

  /* Open file, write main header and initiate matrix list */
  fp=ecat7Create(fname, &main_header);
  if(fp==NULL) {free(fdata); imgSetStatus(img, STATUS_NOWRITEPERM); return(6);}

  /* Set (most of) subheader contents */
  imgSetEcat7SHeader(img, &polmap_header);

  /* Write each matrix */
  for(fi=0; fi<img->dimt; fi++) {

    /* Create new matrix id (i.e. matnum) */
    matrixId=ecat7_val_to_id(fi+1, 1, 1, 0, 0);

    /* Copy matrix pixel values to fdata */
    fptr=fdata;
    for(pi=0; pi<img->dimz; pi++)
      for(yi=0; yi<img->dimy; yi++) for(xi=0; xi<img->dimx; xi++)
        *fptr++=img->m[pi][yi][xi][fi];

    /* Write subheader and data */
    fptr=fdata;
    polmap_header.frame_start_time=(int)temp_roundf(1000.*img->start[fi]);
    polmap_header.frame_duration=
      (int)temp_roundf(1000.*(img->end[fi]-img->start[fi]));
    /*ecat7PrintImageheader(&image_header, stdout);*/
    ret=ecat7WritePolarmapMatrix(fp, matrixId, &polmap_header, fptr);
    if(ret) {
      if(IMG_TEST) {printf("matrixId=%d ret=%d\n", matrixId, ret);}
      free(fdata); fclose(fp); imgSetStatus(img, STATUS_DISKFULL); return(7);
    }

  } /* next matrix */
  free(fdata); fclose(fp);

  imgSetStatus(img, STATUS_OK);
  return(0);
}
/*****************************************************************************/

/*****************************************************************************/
/*!
 * Copy ECAT 7 main header information into IMG
 *
 * @param img target structure
 * @param h source structure
 */
void imgGetEcat7MHeader(IMG *img, ECAT7_mainheader *h) 
{
  img->scanner=h->system_type;
  imgUnitFromEcat7(img, h);
  strlcpy(img->radiopharmaceutical, h->radiopharmaceutical, 32);
  img->isotopeHalflife=h->isotope_halflife;
  img->scanStart=h->scan_start_time;
  img->axialFOV=10.0*h->distance_scanned;
  img->transaxialFOV=10.0*h->transaxial_fov;
  strlcpy(img->studyNr, h->study_type, MAX_STUDYNR_LEN+1);
  strlcpy(img->patientName, h->patient_name, 32);
  strncpy(img->patientID, h->patient_id, 16);
  img->sizez=10.0*h->plane_separation;
  switch(h->file_type) {
    case ECAT7_IMAGE8:
    case ECAT7_IMAGE16:
    case ECAT7_VOLUME8:
    case ECAT7_VOLUME16: img->type=IMG_TYPE_IMAGE; break;
    case ECAT7_POLARMAP: img->type=IMG_TYPE_POLARMAP; break;
    default: img->type=IMG_TYPE_RAW;
  }
  img->orientation=h->patient_orientation;
  strlcpy(img->studyDescription, h->study_description, 32);
  strncpy(img->userProcessCode, h->user_process_code, 10);
  img->userProcessCode[10]=(char)0;
  /* If valid study number is found in user_process_code, then take it */  
  if(!img->studyNr[0] && studynr_validity_check(img->userProcessCode))
    strlcpy(img->studyNr, img->userProcessCode, MAX_STUDYNR_LEN+1);
  img->branchingFraction=h->branching_fraction;
}
/*****************************************************************************/

/*****************************************************************************/
/*!
 * Copy information from IMG to ECAT 7 main header
 *
 * @param img source structure
 * @param h target structure
 */
void imgSetEcat7MHeader(IMG *img, ECAT7_mainheader *h) 
{
  h->sw_version=72;
  if(img->type==IMG_TYPE_POLARMAP) {
    strcpy(h->magic_number, ECAT7V_MAGICNR);
    h->file_type=ECAT7_POLARMAP;
  } else if(img->type==IMG_TYPE_RAW) {
    strcpy(h->magic_number, ECAT7S_MAGICNR);
    if(img->_fileFormat==IMG_E7_2D) h->file_type=ECAT7_2DSCAN;
    else h->file_type=ECAT7_3DSCAN;
  } else {
    strcpy(h->magic_number, ECAT7V_MAGICNR);
    if(img->_fileFormat==IMG_E7_2D) h->file_type=ECAT7_IMAGE16;
    else h->file_type=ECAT7_VOLUME16;
  }
  h->system_type=img->scanner;
  h->scan_start_time=img->scanStart;
  h->isotope_halflife=img->isotopeHalflife;
  imgUnitToEcat7(img, h);
  h->ecat_calibration_factor=1.0;
  h->transaxial_fov=img->transaxialFOV/10.0;
  h->num_planes=img->dimz; /* h->num_planes=1; */
  h->num_frames=img->dimt;  
  h->num_gates=1;
  h->num_bed_pos=0;
  h->distance_scanned=img->axialFOV/10.0;
  h->plane_separation=img->sizez/10.0;
  strncpy(h->radiopharmaceutical, img->radiopharmaceutical, 32);
  strcpy(h->isotope_name, imgIsotope(img));
  strlcpy(h->study_type, img->studyNr, 12);
  strcpy(h->patient_name, img->patientName);
  strcpy(h->patient_id, img->patientID);
  h->patient_orientation=img->orientation;
  strcpy(h->study_description, img->studyDescription);
  strlcpy(h->user_process_code, img->userProcessCode, 10);
  h->branching_fraction=img->branchingFraction;
}
/*****************************************************************************/

/*****************************************************************************/
/*!
 * Return the IMG fileformat based on ECAT7 file_type.
 *
 * @param h Ecat7 main header
 * @return IMG._fileFormat code value
 */
int imgGetEcat7Fileformat(ECAT7_mainheader *h) {
  int fileFormat=IMG_UNKNOWN;
  switch(h->file_type) {
    case ECAT7_IMAGE8:
    case ECAT7_IMAGE16:
      fileFormat=IMG_E7_2D; break;
    case ECAT7_VOLUME8:
    case ECAT7_VOLUME16:
      fileFormat=IMG_E7; break;
    case ECAT7_2DSCAN:
      fileFormat=IMG_E7_2D; break;
    case ECAT7_3DSCAN:
    case ECAT7_3DSCAN8:
    case ECAT7_3DSCANFIT:
      fileFormat=IMG_E7; break;
    case ECAT7_POLARMAP:
      fileFormat=IMG_POLARMAP; break;
    default:
      fileFormat=IMG_UNKNOWN; break;
  }
  return fileFormat;
}
/*****************************************************************************/
/*!
  Fill IMG struct header information from an image or sinogram file
   in ECAT 7 format. Information concerning separate frames or planes is not
   filled.
  
  @param fname image or sinogram filename
  @param img pointer to initialized IMG structure
  @return errstatus, which is STATUS_OK (0) when call was successful,
  and >0 in case of an error.
 */
int imgReadEcat7Header(const char *fname, IMG *img) {
  FILE *fp;
  int ret, m, i;
  int planeNr, frameNr, blkNr=0;
  ECAT7_mainheader main_header;
  ECAT7_imageheader image_header;
  ECAT7_2Dscanheader scan2d_header;
  ECAT7_scanheader scan_header;
  ECAT7_polmapheader polmap_header;
  ECAT7_MATRIXLIST mlist;


  if(IMG_TEST) printf("\nimgReadEcat7Header(%s, *img)\n", fname);
  
  /* Check the arguments */
  if(img==NULL) return STATUS_FAULT;
  if(img->status!=IMG_STATUS_INITIALIZED) return STATUS_FAULT;
  imgSetStatus(img, STATUS_FAULT);
  if(fname==NULL) return STATUS_FAULT;
    
  /* Open the file */
  if((fp=fopen(fname, "rb")) == NULL) return STATUS_NOFILE;

  /* Read main header */
  ret=ecat7ReadMainheader(fp, &main_header);
  if(ret) {fclose(fp); return STATUS_NOMAINHEADER;}
  /* Check for magic number */
  if(strncmp(main_header.magic_number, ECAT7V_MAGICNR, 7)!=0) {
    fclose(fp); return STATUS_UNKNOWNFORMAT;}
  /* Check if file_type is supported */
  if(imgEcat7Supported(&main_header)==0) {fclose(fp); return STATUS_UNSUPPORTED;}
  /* Copy main header information into IMG; sets also img.type */
  imgGetEcat7MHeader(img, &main_header);
  if(IMG_TEST>7) printf("img.type := %d\n", img->type);
  img->_fileFormat=imgGetEcat7Fileformat(&main_header);
  if(IMG_TEST>7) printf("img._fileFormat := %d\n", img->_fileFormat);
  if(img->_fileFormat==IMG_UNKNOWN) {fclose(fp); return STATUS_UNSUPPORTED;}

  /* Read matrix list */
  ecat7InitMatlist(&mlist);
  ret=ecat7ReadMatlist(fp, &mlist, IMG_TEST-1);
  if(ret || mlist.matrixNr<1 || ecat7CheckMatlist(&mlist)) {
    fclose(fp); return STATUS_INVALIDMATLIST;}
  /* Make sure that plane and frame numbers are continuous */
  ecat7GatherMatlist(&mlist, 1, 1, 1, 1);
  /* Get plane and frame numbers and check that volume is full */
  ret=ecat7GetPlaneAndFrameNr(&mlist, &main_header, &planeNr, &frameNr);
  if(ret) {ecat7EmptyMatlist(&mlist); fclose(fp); return ret;}
  img->dimz=planeNr;
  img->dimt=frameNr;
  /* Get and check the size of data matrices */
  ret=ecat7GetMatrixBlockSize(&mlist, &blkNr);
  if(ret) {ecat7EmptyMatlist(&mlist); fclose(fp); return ret;}

  /* Read one subheader */
  if(IMG_TEST>5) printf("main_header.file_type := %d\n", main_header.file_type);
  m=0;
  switch(main_header.file_type) {
    case ECAT7_IMAGE8:
    case ECAT7_IMAGE16:
    case ECAT7_VOLUME8:
    case ECAT7_VOLUME16:
      ret=ecat7ReadImageheader(fp, mlist.matdir[m].strtblk, &image_header);
      break;
    case ECAT7_2DSCAN:
      ret=ecat7Read2DScanheader(fp, mlist.matdir[m].strtblk, &scan2d_header);
      break;
    case ECAT7_3DSCAN:
    case ECAT7_3DSCAN8:
    case ECAT7_3DSCANFIT:
      ret=ecat7ReadScanheader(fp, mlist.matdir[m].strtblk, &scan_header);
      break;
    case ECAT7_POLARMAP:
      ret=ecat7ReadPolmapheader(fp, mlist.matdir[m].strtblk, &polmap_header);
      break;
    default: ret=-1;
  }
  /* Free locally allocated memory and close the file */
  ecat7EmptyMatlist(&mlist); fclose(fp);
  /* Check whether subheader was read */
  if(ret) return STATUS_NOSUBHEADER;

  /* Get the following information in the subheader:
     dimensions x, y and z; datatype;
     image decay correction on/off, zoom, pixel size and resolution;
     sinogram sample distance; matrix transformation parameters.
   */
  switch(main_header.file_type) {
    case ECAT7_IMAGE8:
    case ECAT7_IMAGE16:
    case ECAT7_VOLUME8:
    case ECAT7_VOLUME16:
      img->dimx=image_header.x_dimension; img->dimy=image_header.y_dimension;
      if(image_header.num_dimensions>2 && image_header.z_dimension>1)
        /*planeNr=*/ img->dimz=image_header.z_dimension;
      img->_dataType=image_header.data_type;
      if(image_header.decay_corr_fctr>1.0) img->decayCorrection=IMG_DC_CORRECTED;
      img->zoom=image_header.recon_zoom;
      img->sizex=10.*image_header.x_pixel_size;
      img->sizey=10.*image_header.y_pixel_size;
      img->sizez=10.*image_header.z_pixel_size;
      img->resolutionx=10.*image_header.x_resolution;
      img->resolutiony=10.*image_header.y_resolution;
      img->resolutionz=10.*image_header.z_resolution;
      img->xform[0]=NIFTI_XFORM_UNKNOWN;
      img->xform[1]=NIFTI_XFORM_SCANNER_ANAT;
      img->quatern[6]=img->sizex; img->quatern[9]=img->sizex;
      img->quatern[11]=img->sizey; img->quatern[13]=img->sizey;
      img->quatern[16]=img->sizez; img->quatern[17]=img->sizez;
      img->mt[0]=image_header.mt_1_1;
      img->mt[1]=image_header.mt_1_2;
      img->mt[2]=image_header.mt_1_3;
      img->mt[3]=image_header.mt_1_4;
      img->mt[4]=image_header.mt_2_1;
      img->mt[5]=image_header.mt_2_2;
      img->mt[6]=image_header.mt_2_3;
      img->mt[7]=image_header.mt_2_4;
      img->mt[8]=image_header.mt_3_1;
      img->mt[9]=image_header.mt_3_2;
      img->mt[10]=image_header.mt_3_3;
      img->mt[11]=image_header.mt_3_4;
      break;
    case ECAT7_2DSCAN:
      img->dimx=scan2d_header.num_r_elements;
      img->dimy=scan2d_header.num_angles;
      if(scan2d_header.num_dimensions>2 && scan2d_header.num_z_elements>1)
        planeNr=img->dimz=scan2d_header.num_z_elements;
      img->_dataType=scan2d_header.data_type;
      if(scan2d_header.x_resolution>0.0)
        img->sampleDistance=10.0*scan2d_header.x_resolution;
      else
        img->sampleDistance=10.0*main_header.bin_size;
      break;
    case ECAT7_3DSCAN:
    case ECAT7_3DSCAN8:
    case ECAT7_3DSCANFIT:
      img->dimx=scan_header.num_r_elements; img->dimy=scan_header.num_angles;
      for(i=img->dimz=0; i<64; i++) img->dimz+=scan_header.num_z_elements[i];
      /* planeNr=img->dimz; */
      img->_dataType=scan_header.data_type;
      if(scan_header.x_resolution>0.0)
        img->sampleDistance=10.0*scan_header.x_resolution;
      else
        img->sampleDistance=10.0*main_header.bin_size;
      break;
    case ECAT7_POLARMAP:
      img->dimy=img->dimz=1;
      img->polarmap_num_rings=polmap_header.num_rings;
      if(img->polarmap_num_rings>MAX_POLARMAP_NUM_RINGS)
        return STATUS_INVALIDPOLARMAP;
      for(i=0; i<img->polarmap_num_rings; i++) {
        img->polarmap_sectors_per_ring[i]=polmap_header.sectors_per_ring[i];
        img->polarmap_ring_position[i]=polmap_header.ring_position[i];
        img->polarmap_ring_angle[i]=polmap_header.ring_angle[i];
      }
      img->polarmap_start_angle=polmap_header.start_angle;
      for(i=0, img->dimx=0; i<img->polarmap_num_rings; i++)
        img->dimx+=img->polarmap_sectors_per_ring[i];
      /* pixel_size actually contains volume, in cubic cm */
      img->sizex=img->sizey=img->sizez=0.001*polmap_header.pixel_size;
      break;
  }

  imgSetStatus(img, STATUS_OK);
  return STATUS_OK;
}
/*****************************************************************************/

/*****************************************************************************/
/*!
 * Check whether read functions in IMG library support this ECAT 7.x file_type.
 *
 * @param h Ecat7 main header
 * @return 1 if supported, 0 if not.
 */
int imgEcat7Supported(ECAT7_mainheader *h) {
  if(h==NULL) return(0);
  if(h->file_type==ECAT7_VOLUME8) return(1);
  if(h->file_type==ECAT7_VOLUME16) return(1);
  if(h->file_type==ECAT7_IMAGE8) return(1);
  if(h->file_type==ECAT7_IMAGE16) return(1);
  if(h->file_type==ECAT7_2DSCAN) return(1);
  if(h->file_type==ECAT7_3DSCAN) return(1);
  if(h->file_type==ECAT7_3DSCAN8) return(1);
  if(h->file_type==ECAT7_3DSCANFIT) return(1);
  if(h->file_type==ECAT7_POLARMAP) return(1);
  return(0);
}
/*****************************************************************************/

/*****************************************************************************/
/*!
 * Read the first frame from an ECAT 7 file into IMG data structure.
 *
 * @param fname name of file from which IMG contents will be read
 * @param img pointer to the initiated but not preallocated IMG data
 * @return errstatus, which is STATUS_OK (0) when call was successful,
 * and >0 in case of an error.
 */
int imgReadEcat7FirstFrame(const char *fname, IMG *img) {
  int ret=0;

  if(IMG_TEST) printf("\nimgReadEcat7FirstFrame(%s, *img)\n", fname);
  /* Check the input */
  if(img==NULL) return STATUS_FAULT;
  if(img->status!=IMG_STATUS_INITIALIZED) return STATUS_FAULT;
  imgSetStatus(img, STATUS_FAULT);
  if(fname==NULL) return STATUS_FAULT;

  /* Read header information from file */
  ret=imgReadEcat7Header(fname, img); if(ret) return(ret);
  if(IMG_TEST>3) imgInfo(img);

  /* Allocate memory for one frame */
  img->dimt=1;
  ret=imgAllocate(img, img->dimz, img->dimy, img->dimx, img->dimt);
  if(ret) return STATUS_NOMEMORY;

  /* Read the first frame */
  ret=imgReadEcat7Frame(fname, 1, img, 0); if(ret) return(ret); 

  /* All went well */
  imgSetStatus(img, STATUS_OK);
  return STATUS_OK;
}
/*****************************************************************************/

/*****************************************************************************/
/*!
 * Read a specified frame from an ECAT 7 file into preallocated IMG
 *  data structure. IMG header is assumed to be filled correctly before
 *  calling this function, except for information concerning separate planes
 *  and this frame, which is filled here.
 *  If frame does not exist, then and only then STATUS_NOMATRIX is returned.
 *
 * @param fname name of file from which IMG contents will be read
 * @param frame_to_read frame which will be read [1..frameNr]
 * @param img pointer to the IMG data. Place for the frame must be preallocated
 * @param frame_index IMG frame index [0..dimt-1] where data will be placed
 * @return errstatus, which is STATUS_OK (0) when call was successful,
 * and >0 in case of an error.
 */
int imgReadEcat7Frame(
  const char *fname, int frame_to_read, IMG *img, int frame_index
) {
  FILE *fp;
  int ret, m, i, pi, xi, yi, frame, plane, seqplane, pxlNr;
  /* int blkNr=0; */
  ECAT7_mainheader main_header;
  ECAT7_imageheader image_header;
  ECAT7_2Dscanheader scan2d_header;
  ECAT7_scanheader scan_header;
  ECAT7_polmapheader polmap_header;
  ECAT7_MATRIXLIST mlist;
  ECAT7_Matval matval;
  float *fdata=NULL, *fptr;


  if(IMG_TEST) printf("\nimgReadEcat7Frame(%s, %d, *img, %d)\n",
    fname, frame_to_read, frame_index);
    
  /* Check the input */
  if(img==NULL) return STATUS_FAULT;
  if(img->status!=IMG_STATUS_OCCUPIED) return STATUS_FAULT;
  if(fname==NULL) return STATUS_FAULT;
  if(frame_index<0 || frame_index>img->dimt-1) return STATUS_FAULT;
  if(frame_to_read<1) return STATUS_FAULT;

  /* Open file for read */
  if((fp=fopen(fname, "rb")) == NULL) return STATUS_NOFILE;

  /* Read main header */
  ret=ecat7ReadMainheader(fp, &main_header);
  if(ret) {fclose(fp); return STATUS_NOMAINHEADER;}
  /* Check for magic number */
  if(strncmp(main_header.magic_number, ECAT7V_MAGICNR, 7)!=0) {
    fclose(fp); return STATUS_UNKNOWNFORMAT;}
  /* Check if file_type is supported */
  if(imgEcat7Supported(&main_header)==0) {
    fclose(fp); return STATUS_UNSUPPORTED;}
  /* Read matrix list and nr and sort it */
  ecat7InitMatlist(&mlist);
  ret=ecat7ReadMatlist(fp, &mlist, IMG_TEST-1);
  if(ret) {fclose(fp); return STATUS_NOMATLIST;}
  if(mlist.matrixNr<=0 || ecat7CheckMatlist(&mlist)) {
    fclose(fp); ecat7EmptyMatlist(&mlist); return STATUS_INVALIDMATLIST;}
  /* Make sure that plane and frame numbers are continuous */
  ecat7GatherMatlist(&mlist, 1, 1, 1, 1);
  ecat7SortMatlistByFrame(&mlist); /* printf("matlist sorted\n"); */
  /* Calculate and check the size of one data matrix */
  /*ret=ecat7GetMatrixBlockSize(&mlist, &blkNr);
    if(ret) {fclose(fp); return ret;}*/

  /* Read all matrices that belong to the required frame */
  /*blkNr=-1;*/
  ret=0; seqplane=-1; pxlNr=img->dimx*img->dimy;
  for(m=0; m<mlist.matrixNr; m++) {
    /* get frame and plane */
    ecat7_id_to_val(mlist.matdir[m].id, &matval); plane=matval.plane;
    if(main_header.num_frames>=main_header.num_gates) frame=matval.frame;
    else frame=matval.gate; /* printf("frame=%d plane=%d\n", frame, plane); */
    if(frame!=frame_to_read) continue; /* do not process other frames */
    /*if(img->dimz>1) seqplane++; else seqplane=plane-1; */
    if(img->_fileFormat==IMG_E7_2D) seqplane=plane-1; else seqplane++;
    
    /* Read subheader and data */
    if(IMG_TEST>4) printf("reading matrix %d,%d\n", frame, plane);
    if(img->type==IMG_TYPE_IMAGE) { /* 2D or 3D image */
      ret=ecat7ReadImageMatrix(fp, mlist.matdir[m].strtblk,
              mlist.matdir[m].endblk, &image_header, &fdata);
    } else if(img->type==IMG_TYPE_POLARMAP) {  /* polarmap */
      ret=ecat7ReadPolarmapMatrix(fp, mlist.matdir[m].strtblk,
              mlist.matdir[m].endblk, &polmap_header, &fdata);
    } else if(img->dimz>1) { /* 3D sinogram */
      ret=ecat7ReadScanMatrix(fp, mlist.matdir[m].strtblk,
            mlist.matdir[m].endblk, &scan_header, &fdata);
    } else { /* 2D sinogram */
      ret=ecat7Read2DScanMatrix(fp, mlist.matdir[m].strtblk,
            mlist.matdir[m].endblk, &scan2d_header, &fdata);
    }
    if(ret || fdata==NULL) {
      fclose(fp); ecat7EmptyMatlist(&mlist); return STATUS_NOMATRIX;}

    /* Copy information concerning this frame and make correction to data */
    if(img->type==IMG_TYPE_IMAGE) {
      img->start[frame_index]=image_header.frame_start_time/1000.;
      img->end[frame_index]=
        img->start[frame_index]+image_header.frame_duration/1000.;
      img->mid[frame_index]=0.5*(img->start[frame_index]+img->end[frame_index]);
      if(image_header.decay_corr_fctr>1.0) {
        img->decayCorrFactor[frame_index]=image_header.decay_corr_fctr;
        img->decayCorrection=IMG_DC_CORRECTED;
      } else {
        img->decayCorrFactor[frame_index]=0.0;
        img->decayCorrection=IMG_DC_UNKNOWN;
      }
    } else if(img->type==IMG_TYPE_POLARMAP) {  /* polarmap */
      img->start[frame_index]=polmap_header.frame_start_time/1000.;
      img->end[frame_index]=
        img->start[frame_index]+polmap_header.frame_duration/1000.;
      img->mid[frame_index]=0.5*(img->start[frame_index]+img->end[frame_index]);
    } else if(img->dimz>1) {
      img->start[frame_index]=scan_header.frame_start_time/1000.;
      img->end[frame_index]=
        img->start[frame_index]+scan_header.frame_duration/1000.;
      img->mid[frame_index]=0.5*(img->start[frame_index]+img->end[frame_index]);
      /* also, correct for dead-time */
      if(scan_header.deadtime_correction_factor>0.0)
        for(i=0, fptr=fdata; i<img->dimz*pxlNr; i++, fptr++)
          *fptr*=scan_header.deadtime_correction_factor;
      img->prompts[frame_index]=(float)scan_header.prompts;
      img->randoms[frame_index]=scan_header.delayed;
    } else {
      img->start[frame_index]=scan2d_header.frame_start_time/1000.;
      img->end[frame_index]=
        img->start[frame_index]+scan2d_header.frame_duration/1000.;
      img->mid[frame_index]=0.5*(img->start[frame_index]+img->end[frame_index]);
      /* also, correct for dead-time */
      if(scan2d_header.deadtime_correction_factor>0.0)
        for(i=0, fptr=fdata; i<pxlNr; i++, fptr++)
          *fptr*=scan2d_header.deadtime_correction_factor;
      img->prompts[frame_index]=(float)scan2d_header.prompts;
      img->randoms[frame_index]=scan2d_header.delayed;
    }
    /* Copy matrix data through volume planes */
    if(img->_fileFormat!=IMG_E7_2D) {
    /* if(img->dimz>1) { */
      for(pi=0; pi<img->dimz; pi++) {
        if(IMG_TEST>5)
          printf("  putting data into m[%d][][][%d]\n", pi, frame_index);
        for(yi=0, fptr=fdata+pi*pxlNr; yi<img->dimy; yi++)
          for(xi=0; xi<img->dimx; xi++)
            img->m[pi][yi][xi][frame_index]= *fptr++;
      }
    } else {
        if(IMG_TEST>5)
          printf("  putting data into m[%d][][][%d]\n", seqplane, frame_index);
        for(yi=0, fptr=fdata; yi<img->dimy; yi++)
          for(xi=0; xi<img->dimx; xi++)
            img->m[seqplane][yi][xi][frame_index]= *fptr++;
        img->planeNumber[seqplane]=plane;
    }
    free(fdata);
    /* Calibrate */
    if(main_header.ecat_calibration_factor>0.0)
      for(pi=0; pi<img->dimz; pi++)
        for(yi=0; yi<img->dimy; yi++) for(xi=0; xi<img->dimx; xi++)
          img->m[pi][yi][xi][frame_index]*=main_header.ecat_calibration_factor;
  } /* next matrix */
  if(IMG_TEST>3) printf("end of matrices.\n");
  ecat7EmptyMatlist(&mlist);
  fclose(fp);

  /* seqplane is <0 only if this frame did not exist at all; return -1 */
  if(IMG_TEST>4) printf("last_seqplane := %d.\n", seqplane);
  if(seqplane<0) return STATUS_NOMATRIX;

  /* check that correct number of planes was read */
  if(seqplane>0 && (seqplane+1 != img->dimz)) return STATUS_MISSINGMATRIX;

  /* All went well */
  imgSetStatus(img, STATUS_OK);
  return STATUS_OK;
}
/*****************************************************************************/

/*****************************************************************************/
/*!
 * Write one PET frame from IMG data struct into ECAT 7 image or sinogram
 *  file; format is specified in IMG struct. This function can be called
 *  repeatedly to write all frames one at a time to conserve memory.
 *  However, file with just mainheader and matrix list without any previous
 *  frame is not accepted.
 *
 * @param fname name of file where IMG contents will be written.
 *  If file does not exist, it is created.
 *  Make sure to delete existing file, unless you want to add data
 * @param frame_to_write PET frame number (1..frameNr) which will be written:
 *  If set to 0, frame data will be written to an existing or new PET file as
 *  a new frame, never overwriting existing data.
 *  If >0, then frame data is written as specified frame number, overwriting
 *  any data existing with the same frame number
 * @param img pointer to the IMG data struct
 * @param frame_index IMG frame index (0..dimt-1) which will be written
 * @return errstatus, which is STATUS_OK (0) when call was successful,
 * and >0 in case of an error.
 */
int imgWriteEcat7Frame(
  const char *fname, int frame_to_write, IMG *img, int frame_index
) {
  IMG test_img;
  int ret=0, pxlNr, zi, xi, yi, matrixId;
  ECAT7_mainheader main_header;
  ECAT7_imageheader image_header;
  ECAT7_scanheader scan_header;
  ECAT7_2Dscanheader scan2d_header;
  ECAT7_polmapheader polmap_header;
  void *sub_header=NULL;
  FILE *fp;
  float *fdata=NULL, *fptr;


  if(IMG_TEST) printf("\nimgWriteEcat7Frame(%s, %d, *img, %d)\n",
    fname, frame_to_write, frame_index);
    
  /*
   *  Check the input 
   */
  if(fname==NULL) return STATUS_FAULT;
  if(img==NULL) return STATUS_FAULT;
  if(img->status!=IMG_STATUS_OCCUPIED) return STATUS_FAULT;
  if(frame_to_write<0) return STATUS_FAULT;
  if(frame_index<0 || frame_index>=img->dimt) return STATUS_FAULT;
  if(img->_fileFormat!=IMG_E7 &&
     img->_fileFormat!=IMG_POLARMAP &&
     img->_fileFormat!=IMG_E7_2D) return STATUS_FAULT;

  /* Initiate headers */
  memset(&main_header, 0, sizeof(ECAT7_mainheader));
  memset(&image_header,0, sizeof(ECAT7_imageheader));
  memset(&scan_header, 0, sizeof(ECAT7_scanheader));
  memset(&scan2d_header, 0, sizeof(ECAT7_2Dscanheader));
  memset(&polmap_header,0, sizeof(ECAT7_polmapheader));
  imgInit(&test_img);


  /*
   *  If file does not exist, then create it with new mainheader,
   *  and if it does exist, then read and check header information.
   *  Create or edit mainheader to contain correct frame nr.   
   *  In any case, leave file pointer open for write.   
   */
  if(access(fname, 0) == -1) { /* file does not exist */

    /* Set main header */
    imgSetEcat7MHeader(img, &main_header);
    main_header.bin_size=img->sampleDistance/10.0;
    if(frame_to_write==0) frame_to_write=1;
    main_header.num_frames=frame_to_write;

    /* Open file, write main header and initiate matrix list */
    fp=ecat7Create(fname, &main_header); if(fp==NULL) return STATUS_NOWRITEPERM;

  } else { /* file exists */
  
    /* Read header information for checking */
    ret=imgReadEcat7Header(fname, &test_img); if(ret!=0) return ret;
    /* Check that file format is the same */
    if(img->_fileFormat!=test_img._fileFormat || img->type!=test_img.type)
      return STATUS_WRONGFILETYPE;
    /* Check that matrix sizes are the same */
    if(img->dimz!=test_img.dimz || img->dimx!=test_img.dimx ||
       img->dimy!=test_img.dimy)
      return STATUS_VARMATSIZE;
    imgEmpty(&test_img);

    /* Open the file for read and write */
    if((fp=fopen(fname, "r+b"))==NULL) return STATUS_NOWRITEPERM;

    /* Read the mainheader, set new frame number, and write it back */
    if((ret=ecat7ReadMainheader(fp, &main_header))!=0)
      return STATUS_NOMAINHEADER;
    if(frame_to_write==0) frame_to_write=main_header.num_frames+1;
    if(main_header.num_frames<frame_to_write)
      main_header.num_frames=frame_to_write;
    if((ret=ecat7WriteMainheader(fp, &main_header))!=0)
      return STATUS_NOWRITEPERM;
    
  }
  if(IMG_TEST>2) {
    printf("frame_to_write := %d\n", frame_to_write);
  }

  /* Allocate memory for matrix float data */
  pxlNr=img->dimx*img->dimy*img->dimz; /* for 2D too! */
  fdata=(float*)malloc(pxlNr*sizeof(float));
  if(fdata==NULL) {fclose(fp); return STATUS_NOMEMORY;}
  
  /* Set matrix subheader */
  if(img->_fileFormat==IMG_E7) {
    if(img->type==IMG_TYPE_RAW) sub_header=(void*)&scan_header;
    else sub_header=&image_header;
  } else if(img->_fileFormat==IMG_E7_2D) {
    if(img->type==IMG_TYPE_RAW) sub_header=(void*)&scan_header;
    else sub_header=(void*)&image_header;
  } else if(img->_fileFormat==IMG_POLARMAP) {
    sub_header=(void*)&polmap_header;
  } else {fclose(fp); return STATUS_FAULT;}
  imgSetEcat7SHeader(img, sub_header);

  /* Copy matrix pixel values to fdata */
  fptr=fdata;
  for(zi=0; zi<img->dimz; zi++)
    for(yi=0; yi<img->dimy; yi++) for(xi=0; xi<img->dimx; xi++)
      *fptr++=img->m[zi][yi][xi][frame_index];

  /* Write subheader and data, and set the rest of subheader contents */
  fptr=fdata; ret=-1;
  if(img->_fileFormat==IMG_E7) {
    /* Create new matrix id (i.e. matnum) */
    matrixId=ecat7_val_to_id(frame_to_write, 1, 1, 0, 0);
    if(img->type==IMG_TYPE_RAW) {
      scan_header.frame_start_time=
       (int)temp_roundf(1000.*img->start[frame_index]);
      scan_header.frame_duration=
        (int)temp_roundf(1000.*(img->end[frame_index]-img->start[frame_index]));
      scan_header.prompts=temp_roundf(img->prompts[frame_index]);
      scan_header.delayed=temp_roundf(img->randoms[frame_index]);
      ret=ecat7WriteScanMatrix(fp, matrixId, &scan_header, fptr);
    } else {
      image_header.frame_start_time=
       (int)temp_roundf(1000.*img->start[frame_index]);
      image_header.frame_duration=
        (int)temp_roundf(1000.*(img->end[frame_index]-img->start[frame_index]));
      if(img->decayCorrection==IMG_DC_CORRECTED)
        image_header.decay_corr_fctr=img->decayCorrFactor[frame_index];
      else
        image_header.decay_corr_fctr=0.0;
      /*ecat7PrintImageheader(&image_header, stdout);*/
      ret=ecat7WriteImageMatrix(fp, matrixId, &image_header, fptr);
    }
  } else if(img->_fileFormat==IMG_E7_2D) {
    for(zi=0; zi<img->dimz; zi++, fptr+=img->dimx*img->dimy) {
      /* Create new matrix id (i.e. matnum) */
      matrixId=ecat7_val_to_id(frame_to_write, img->planeNumber[zi], 1, 0, 0);
      if(img->type==IMG_TYPE_RAW) {
        scan2d_header.frame_start_time=
	  (int)temp_roundf(1000.*img->start[frame_index]);
        scan2d_header.frame_duration=
          (int)temp_roundf(1000.*(img->end[frame_index]-img->start[frame_index]));
        scan2d_header.prompts=temp_roundf(img->prompts[frame_index]);
        scan2d_header.delayed=temp_roundf(img->randoms[frame_index]);
        ret=ecat7Write2DScanMatrix(fp, matrixId, &scan2d_header, fptr);
      } else {
        image_header.frame_start_time=
	  (int)temp_roundf(1000.*img->start[frame_index]);
        image_header.frame_duration=
	  (int)temp_roundf(1000.*(img->end[frame_index]-img->start[frame_index]));
        if(img->decayCorrection==IMG_DC_CORRECTED)
          image_header.decay_corr_fctr=img->decayCorrFactor[frame_index];
        else
          image_header.decay_corr_fctr=0.0;
        ret=ecat7WriteImageMatrix(fp, matrixId, &image_header, fptr);
      }
      if(ret!=0) break;
    } /* next plane */
  } else if(img->_fileFormat==IMG_POLARMAP) {
    /* Create new matrix id (i.e. matnum) */
    matrixId=ecat7_val_to_id(frame_to_write, 1, 1, 0, 0);
    polmap_header.frame_start_time=
     (int)temp_roundf(1000.*img->start[frame_index]);
    polmap_header.frame_duration=
        (int)temp_roundf(1000.*(img->end[frame_index]-img->start[frame_index]));
    ret=ecat7WritePolarmapMatrix(fp, matrixId, &polmap_header, fptr);
  }
  free(fdata); fclose(fp);
  if(ret) return STATUS_DISKFULL;

  return STATUS_OK;
}
/*****************************************************************************/

/*****************************************************************************/
/*!
 * Set ECAT7 subheader based on IMG contents
 *
 * @param img image structure
 * @param h Ecat7 image, scan, 2D scan or polar map header
 */
void imgSetEcat7SHeader(IMG *img, void *h) {
  ECAT7_imageheader *image_header;
  ECAT7_scanheader *scan_header;
  ECAT7_2Dscanheader *scan2d_header;
  ECAT7_polmapheader *polmap_header;
  int i;

  if(img->type==IMG_TYPE_POLARMAP) {
    polmap_header=(ECAT7_polmapheader*)h;
    polmap_header->data_type=ECAT7_SUNI2;
    polmap_header->num_rings=img->polarmap_num_rings;
    if(polmap_header->num_rings>32) polmap_header->num_rings=32;
    for(i=0; i<polmap_header->num_rings; i++) {
      polmap_header->sectors_per_ring[i]=img->polarmap_sectors_per_ring[i];
      polmap_header->ring_position[i]=img->polarmap_ring_position[i];
      polmap_header->ring_angle[i]=img->polarmap_ring_angle[i];
    }
    polmap_header->start_angle=258;
    polmap_header->pixel_size=1000.0*img->sizex;  
    polmap_header->quant_units=0; /* default, see main header */
  } else if(img->type==IMG_TYPE_RAW) {
    if(img->_fileFormat==IMG_E7_2D) {
      scan2d_header=(ECAT7_2Dscanheader*)h;
      scan2d_header->num_dimensions=2;
      scan2d_header->num_z_elements=1;
      scan2d_header->data_type=ECAT7_SUNI2;
      scan2d_header->num_r_elements=img->dimx;
      scan2d_header->num_angles=img->dimy;
    } else {
      scan_header=(ECAT7_scanheader*)h;
      scan_header->x_resolution=img->sampleDistance/10.0;
      scan_header->num_dimensions=4;
      if(img->dimz==239) {
        scan_header->num_z_elements[0]=63;
        scan_header->num_z_elements[1]=106;
        scan_header->num_z_elements[2]=70;
      } else {
        scan_header->num_z_elements[0]=img->dimz;
      }
      scan_header->storage_order=1;
      scan_header->data_type=ECAT7_SUNI2;
      scan_header->num_r_elements=img->dimx;
      scan_header->num_angles=img->dimy;
    }
  } else {
    image_header=(ECAT7_imageheader*)h;
    image_header->data_type=ECAT7_SUNI2;
    image_header->x_dimension=img->dimx;
    image_header->y_dimension=img->dimy;
    image_header->recon_zoom=img->zoom;
    image_header->x_pixel_size=0.1*img->sizex;
    image_header->y_pixel_size=0.1*img->sizey;
    image_header->z_pixel_size=0.1*img->sizez;
    image_header->x_resolution=0.1*img->resolutionx;
    image_header->y_resolution=0.1*img->resolutiony;
    image_header->z_resolution=0.1*img->resolutionz;
    image_header->mt_1_1=img->mt[0];
    image_header->mt_1_2=img->mt[1];
    image_header->mt_1_3=img->mt[2];
    image_header->mt_1_4=img->mt[3];
    image_header->mt_2_1=img->mt[4];
    image_header->mt_2_2=img->mt[5];
    image_header->mt_2_3=img->mt[6];
    image_header->mt_2_4=img->mt[7];
    image_header->mt_3_1=img->mt[8];
    image_header->mt_3_2=img->mt[9];
    image_header->mt_3_3=img->mt[10];
    image_header->mt_3_4=img->mt[11];
    if(img->_fileFormat==IMG_E7_2D) {
      image_header->num_dimensions=2;
      image_header->z_dimension=1;
    } else {
      image_header->num_dimensions=3;
      image_header->z_dimension=img->dimz;
    }
    for(i=0; i<49; i++) image_header->fill_user[i]=0;
  }
}
/*****************************************************************************/

/*****************************************************************************/
