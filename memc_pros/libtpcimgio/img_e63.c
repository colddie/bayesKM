/// @file img_e63.c
/// @author Vesa Oikonen
/// @brief ECAT 6.3 I/O routines for IMG data.
///
/*****************************************************************************/
#include "libtpcimgio.h"
/*****************************************************************************/

/*****************************************************************************/
/** Read all matrices in ECAT file to memory.

    Sinograms are dead-time corrected.

 @sa imgInit

 @return 0 if ok, 1 invalid input, 3 failed to open file for reading,
  4 failed to read main header, 5 failed to read matrix list, 
  6 matrix not found, 7 variable matrix sizes, 8 failed to read matrix
  sub header, 9 failed to allocate memory for data,
  10 failed to read sub header, 11 failed to read matrix data.
 */
int ecat63ReadAllToImg(
  /** Name of the input ECAT 6.3 file. */
  const char *fname,
  /** Data structure in which the file is read; must be initialized
      before using this function. */
  IMG *img
) {
  int                 i, j, m, ret, blkNr=0, dim_x=0, dim_y=0;
  int                 frameNr, planeNr, del_nr=0;
  int                 frame, plane, prev_frame, prev_plane, seqplane;
  FILE               *fp;
  ECAT63_mainheader   main_header;
  MATRIXLIST          mlist;
  MatDir              tmpMatdir;
  Matval              matval;
  ECAT63_imageheader  image_header;
  ECAT63_scanheader   scan_header;
  ECAT63_attnheader   attn_header;
  ECAT63_normheader   norm_header;
  float               scale=1.0;
  short int          *sptr;
  char               *mdata=NULL, *mptr;
  /*int                *iptr;*/


  if(IMG_TEST) printf("ecat63ReadAllToImg(%s, *img)\n", fname);
  /* Check the arguments */
  if(fname==NULL || img==NULL || img->status!=IMG_STATUS_INITIALIZED) {
    if(img) imgSetStatus(img, STATUS_FAULT);
    return 1;
  }

  /* Open input CTI file for read */
  if((fp=fopen(fname, "rb")) == NULL) {
    imgSetStatus(img, STATUS_NOFILE);
    return 3;
  }

  /* Read main header */
  if((ret=ecat63ReadMainheader(fp, &main_header))) {
    imgSetStatus(img, STATUS_NOMAINHEADER);
    return 4;
  }
  if(IMG_TEST>5) ecat63PrintMainheader(&main_header, stdout);

  /* Read matrix list and nr */
  ecat63InitMatlist(&mlist);
  ret=ecat63ReadMatlist(fp, &mlist, ECAT63_TEST-1);
  if(ret) {
    imgSetStatus(img, STATUS_NOMATLIST);
    return 5;
  }
  if(mlist.matrixNr<=0) {
    strcpy(ecat63errmsg, "matrix list is empty");
    return 5;
  }
  /* Sort matrix list */
  for(i=0; i<mlist.matrixNr-1; i++) for(j=i+1; j<mlist.matrixNr; j++) {
    if(mlist.matdir[i].matnum>mlist.matdir[j].matnum) {
      tmpMatdir=mlist.matdir[i];
      mlist.matdir[i]=mlist.matdir[j]; mlist.matdir[j]=tmpMatdir;
    }
  }
  if(IMG_TEST>6) {
    printf("matrix list after sorting:\n");
    ecat63PrintMatlist(&mlist);
  }

  /* Trim extra frames */
  if(main_header.num_frames>0) {
    del_nr=ecat63DeleteLateFrames(&mlist, main_header.num_frames);
    if(IMG_TEST && del_nr>0)
      printf("  %d entries in matrix list will not be used.\n", del_nr);
  }
  /* Calculate the number of planes, frames and gates */
  /* Check that all planes have equal nr of frames (gates) */
  /* and check that frames (gates) are consequentially numbered */
  prev_plane=plane=-1; prev_frame=frame=-1; frameNr=planeNr=0; ret=0;
  for(m=0; m<mlist.matrixNr; m++) if(mlist.matdir[m].matstat==1) {
    /* get frame and plane */
    mat_numdoc(mlist.matdir[m].matnum, &matval);
    plane=matval.plane;
    if(main_header.num_frames>=main_header.num_gates) frame=matval.frame;
    else frame=matval.gate;
    if(IMG_TEST>2) printf("frame=%d plane=%d\n", frame, plane);
    /* did plane change? */
    if(plane!=prev_plane) {
      frameNr=1; planeNr++;
    } else {
      frameNr++;
      /* In each plane, frame (gate) numbers must be continuous */
      if(prev_frame>0 && frame!=prev_frame+1) {ret=1; break;}
    }
    prev_plane=plane; prev_frame=frame;
    /* Calculate and check the size of one data matrix */
    if(m==0) {
      blkNr=mlist.matdir[m].endblk-mlist.matdir[m].strtblk;
    } else if(blkNr!=mlist.matdir[m].endblk-mlist.matdir[m].strtblk) {
      ret=2; break;
    }
  } /* next matrix */
  if(IMG_TEST) printf("frameNr=%d planeNr=%d\n", frameNr, planeNr);
  if(ret==1 || (frameNr*planeNr != mlist.matrixNr-del_nr)) {
    strcpy(ecat63errmsg, "missing matrix");
    ecat63EmptyMatlist(&mlist); fclose(fp);
    return(6); /* this number is used in imgRead() */
  }
  if(ret==2) {
    strcpy(ecat63errmsg, "matrix sizes are different");
    ecat63EmptyMatlist(&mlist); fclose(fp); return(7);
  }
  /* Read x,y-dimensions from 1st matrix subheader */
  m=0; if(main_header.file_type==IMAGE_DATA) {
    ret=ecat63ReadImageheader(fp, mlist.matdir[m].strtblk, &image_header, IMG_TEST-4, NULL);
    dim_x=image_header.dimension_1; dim_y=image_header.dimension_2;
  } else if(main_header.file_type==RAW_DATA) {
    ret=ecat63ReadScanheader(fp, mlist.matdir[m].strtblk, &scan_header, IMG_TEST-4, NULL);
    dim_x=scan_header.dimension_1; dim_y=scan_header.dimension_2;
  } else if(main_header.file_type==ATTN_DATA) {
    ret=ecat63ReadAttnheader(fp, mlist.matdir[m].strtblk, &attn_header, IMG_TEST-10, NULL);
    dim_x=attn_header.dimension_1; dim_y=attn_header.dimension_2;
  } else if(main_header.file_type==NORM_DATA) {
    ret=ecat63ReadNormheader(fp, mlist.matdir[m].strtblk, &norm_header, IMG_TEST-10, NULL);
    dim_x=norm_header.dimension_1; dim_y=norm_header.dimension_2;
  }
  /*pxlNr=dim_x*dim_y;*/
  if(ret) {
    sprintf(ecat63errmsg, "cannot read matrix %u subheader",
            mlist.matdir[m].matnum);
    ecat63EmptyMatlist(&mlist); fclose(fp); return(8);
  }

  /* Allocate memory for ECAT data matrix */
  if(IMG_TEST>1) printf("allocating memory for %d blocks\n", blkNr);
  mdata=(char*)malloc(blkNr*MatBLKSIZE); if(mdata==NULL) {
    strcpy(ecat63errmsg, "out of memory");
    fclose(fp); ecat63EmptyMatlist(&mlist); return(8);
  }
  /* Allocate memory for whole img data */
  ret=imgAllocate(img, planeNr, dim_y, dim_x, frameNr);
  if(ret) {
    sprintf(ecat63errmsg, "out of memory (%d)", ret);
    fclose(fp); ecat63EmptyMatlist(&mlist); free(mdata); return(9);
  }

  /* Fill img info with ECAT main and sub header information */
  img->scanner=main_header.system_type;
  img->unit=main_header.calibration_units; /* this continues below */
  strlcpy(img->radiopharmaceutical, main_header.radiopharmaceutical, 32);
  img->isotopeHalflife=main_header.isotope_halflife;
  img->scanStart=ecat63Scanstarttime(&main_header);
  if(img->scanStart==-1) {
    img->scanStart=0;
    if(IMG_TEST>0) printf("invalid scan_start_time in mainheader\n");
  }
  img->axialFOV=10.*main_header.axial_fov;
  img->transaxialFOV=10.*main_header.transaxial_fov;
  strlcpy(img->studyNr, main_header.study_name, MAX_STUDYNR_LEN+1);
  strlcpy(img->patientName, main_header.patient_name, 32);
  strlcpy(img->patientID, main_header.patient_id, 16);
  img->sizez=10.*main_header.plane_separation;
  if(main_header.file_type==IMAGE_DATA) {
    img->type=IMG_TYPE_IMAGE;
    if(img->unit<1) img->unit=image_header.quant_units;
    img->zoom=image_header.recon_scale;
    if(image_header.decay_corr_fctr>1.0) img->decayCorrection=IMG_DC_CORRECTED;
    else img->decayCorrection=IMG_DC_UNKNOWN;
    img->sizex=img->sizey=10.*image_header.pixel_size;
    if(img->sizez<10.*image_header.slice_width)
      img->sizez=10.*image_header.slice_width;
    img->xform[0]=NIFTI_XFORM_UNKNOWN; // qform
    img->xform[1]=NIFTI_XFORM_SCANNER_ANAT; // sform
    img->quatern[6]=img->sizex; img->quatern[9]=img->sizex;
    img->quatern[11]=img->sizey; img->quatern[13]=img->sizey;
    img->quatern[16]=img->sizez; img->quatern[17]=img->sizez;
  } else if(main_header.file_type==RAW_DATA) {
    img->type=IMG_TYPE_RAW;
    img->decayCorrection=IMG_DC_NONCORRECTED;
  } else if(main_header.file_type==ATTN_DATA) {
    img->type=IMG_TYPE_ATTN;
    img->decayCorrection=IMG_DC_NONCORRECTED;
  } else if(main_header.file_type==NORM_DATA) {
    img->type=IMG_TYPE_RAW;
    img->decayCorrection=IMG_DC_NONCORRECTED;
  }
  strlcpy(img->studyDescription, main_header.study_description, 32);
  strncpy(img->userProcessCode, main_header.user_process_code, 10);
  img->userProcessCode[10]=(char)0;
  /* If valid study number is found in user_process_code, then take it */  
  if(!img->studyNr[0] && studynr_validity_check(img->userProcessCode))
    strlcpy(img->studyNr, img->userProcessCode, MAX_STUDYNR_LEN+1);
    
  /* Set _fileFormat */
  img->_fileFormat=IMG_E63;

  /* Read one ECAT matrix at a time and put them to img */
  prev_plane=plane=-1; seqplane=-1;
  for(m=0; m<mlist.matrixNr; m++) if(mlist.matdir[m].matstat==1) {
    /* get plane */
    mat_numdoc(mlist.matdir[m].matnum, &matval);
    plane=matval.plane;
    /* did plane change? */
    if(plane!=prev_plane) {seqplane++; frame=1;} else frame++;
    prev_plane=plane;
    img->planeNumber[seqplane]=plane;
    /* Read subheader */
    if(main_header.file_type==IMAGE_DATA) {
      ret=ecat63ReadImageheader(fp, mlist.matdir[m].strtblk, &image_header, IMG_TEST-10, NULL);
      if(ret==0 && (dim_x!=image_header.dimension_1 || dim_y!=image_header.dimension_2)) ret=1;
      scale=image_header.quant_scale;
      if(image_header.ecat_calibration_fctr>0.0 &&
         fabs(main_header.calibration_factor/image_header.ecat_calibration_fctr - 1.0)>0.0001)
        scale*=image_header.ecat_calibration_fctr;
      if(img->unit==0 && image_header.quant_units>0)
        img->unit=image_header.quant_units;
      img->_dataType=image_header.data_type;
      img->start[frame-1]=image_header.frame_start_time/1000.;
      img->end[frame-1]=img->start[frame-1]+image_header.frame_duration/1000.;
      img->mid[frame-1]=0.5*(img->start[frame-1]+img->end[frame-1]);
      if(image_header.decay_corr_fctr>1.0)
        img->decayCorrFactor[frame-1]=image_header.decay_corr_fctr;
      else
        img->decayCorrFactor[frame-1]=0.0;
      /* Dead-time correction is assumed to be included in scale factor */
    } else if(main_header.file_type==RAW_DATA) {
      ret=ecat63ReadScanheader(fp, mlist.matdir[m].strtblk, &scan_header, IMG_TEST-10, NULL);
      if(ret==0 && (dim_x!=scan_header.dimension_1 || dim_y!=scan_header.dimension_2)) ret=1;
      scale=scan_header.scale_factor;
      /* Dead-time correction is done here */
      if(scan_header.loss_correction_fctr>0.0)
        scale*=scan_header.loss_correction_fctr;
      img->_dataType=scan_header.data_type;
      img->start[frame-1]=scan_header.frame_start_time/1000.;
      img->end[frame-1]=img->start[frame-1]+scan_header.frame_duration/1000.;
      img->mid[frame-1]=0.5*(img->start[frame-1]+img->end[frame-1]);
      img->sampleDistance=10.0*scan_header.sample_distance;
      img->prompts[frame-1]=(float)scan_header.prompts;
      img->randoms[frame-1]=(float)scan_header.delayed;
    } else if(main_header.file_type==ATTN_DATA) {
      ret=ecat63ReadAttnheader(fp, mlist.matdir[m].strtblk, &attn_header, IMG_TEST-10, NULL);
      if(ret==0 && (dim_x!=attn_header.dimension_1 || dim_y!=attn_header.dimension_2)) ret=1;
      img->_dataType=attn_header.data_type;
      scale=attn_header.scale_factor;
      img->sampleDistance=10.0*attn_header.sample_distance;
    } else if(main_header.file_type==NORM_DATA) {
      ret=ecat63ReadNormheader(fp, mlist.matdir[m].strtblk, &norm_header, IMG_TEST-10, NULL);
      if(ret==0 && (dim_x!=norm_header.dimension_1 || dim_y!=norm_header.dimension_2)) ret=1;
      scale=norm_header.scale_factor;
      img->_dataType=norm_header.data_type;
    } else img->_dataType=-1;
    if(ret) {
      sprintf(ecat63errmsg, "cannot read matrix %u subheader", mlist.matdir[m].matnum);
      ecat63EmptyMatlist(&mlist); free(mdata); fclose(fp); return(10);
    }
    if(main_header.calibration_factor>0.0) scale*=main_header.calibration_factor;
    if(IMG_TEST>2) printf("scale=%e datatype=%d\n", scale, img->_dataType);
    /* Read the pixel data */
    ret=ecat63ReadMatdata(fp, mlist.matdir[m].strtblk+1,
         mlist.matdir[m].endblk-mlist.matdir[m].strtblk,
         mdata, img->_dataType);
    if(ret) {
      strcpy(ecat63errmsg, "cannot read matrix data");
      ecat63EmptyMatlist(&mlist); free(mdata); fclose(fp); return(11);
    }
    if(img->_dataType==BYTE_TYPE) {
      for(i=0, mptr=mdata; i<dim_y; i++) for(j=0; j<dim_x; j++, mptr++) {
        img->m[seqplane][i][j][frame-1]=scale*(float)(*mptr);
      }
    } else if(img->_dataType==VAX_I2 || img->_dataType==SUN_I2) {
      for(i=0, mptr=mdata; i<dim_y; i++) for(j=0; j<dim_x; j++, mptr+=2) {
        sptr=(short int*)mptr;
        img->m[seqplane][i][j][frame-1]=scale*(float)(*sptr);
      }
    } else if(img->_dataType==VAX_I4 || img->_dataType==SUN_I4) {
      for(i=0, mptr=mdata; i<dim_y; i++) for(j=0; j<dim_x; j++, mptr+=4) {
        /*iptr=(int*)mptr;*/
        img->m[seqplane][i][j][frame-1]=1.0; /*scale*(float)(*iptr);*/
      }
    } else if(img->_dataType==VAX_R4 || img->_dataType==IEEE_R4) {
      for(i=0, mptr=mdata; i<dim_y; i++) for(j=0; j<dim_x; j++, mptr+=4) {
        memcpy(&img->m[seqplane][i][j][frame-1], mptr, 4);
        img->m[seqplane][i][j][frame-1]*=scale;
      }
    }
  } /* next matrix */

  /* Set the unit */
  imgUnitFromEcat(img, img->unit);

  /* Free memory and close file */
  free(mdata);
  ecat63EmptyMatlist(&mlist);
  fclose(fp);

  /* For saving, only 2-byte VAX data types are allowed */
  if(img->_dataType==VAX_I4 || img->_dataType==VAX_R4) img->_dataType=VAX_I2;

  return(0);
}
/*****************************************************************************/

/*****************************************************************************/
/*!
 * Write all matrices in memory to the ECAT file.
 *
 * @param fname name of the output ECAT 6.3 file, If ECAT file exists, 
 * it is renamed as <code>filename%</code>
 * @param img data structure from which the data is written
 * @return 0 if ok, 1 invalid data, 2 image status is not 'occupied', 
 * 3 failed to create file, 4 failed to allocate memory for data, 
 * 9 failed to write data
 */
int ecat63WriteAllImg(const char *fname, IMG *img) {
  int                 frame, plane, n, i, j, m, ret=0;
  float               f, fmax, fmin, g, scale;
  short int          *sdata, *sptr, smin, smax;
  FILE               *fp;
  ECAT63_mainheader   main_header;
  ECAT63_imageheader  image_header;
  ECAT63_scanheader   scan_header;
  struct tm           tm;

  if(IMG_TEST) printf("ecat63WriteAllImg(%s, *img)\n", fname);
  /* Check the arguments */
  if(fname==NULL) {strcpy(ecat63errmsg, "invalid ECAT filename"); return(1);}
  if(img==NULL || img->status!=IMG_STATUS_OCCUPIED) {
    strcpy(ecat63errmsg, "image data is empty"); return(2);}
  if(img->_dataType<1) img->_dataType=VAX_I2;

  /* Initiate headers */
  memset(&main_header, 0, sizeof(ECAT63_mainheader));
  memset(&image_header,0, sizeof(ECAT63_imageheader));
  memset(&scan_header, 0, sizeof(ECAT63_scanheader));

  /* Make a main header */
  main_header.sw_version=2;
  main_header.num_planes=img->dimz;  
  main_header.num_frames=img->dimt;  
  main_header.num_gates=1;
  main_header.num_bed_pos=1;
  if(img->type==IMG_TYPE_IMAGE) main_header.file_type=IMAGE_DATA;
  else if(img->type==IMG_TYPE_RAW) main_header.file_type=RAW_DATA;
  else {strcpy(ecat63errmsg, "invalid filetype"); return(1);}
  main_header.data_type=img->_dataType;
  if(img->scanner>0) main_header.system_type=img->scanner;
  else main_header.system_type=ECAT63_SYSTEM_TYPE_DEFAULT;
  main_header.calibration_factor=1.0;
  main_header.axial_fov=img->axialFOV/10.0;
  main_header.transaxial_fov=img->transaxialFOV/10.0;
  main_header.plane_separation=img->sizez/10.0;
  main_header.calibration_units=imgUnitToEcat6(img);
  strncpy(main_header.radiopharmaceutical, img->radiopharmaceutical, 32);
  if(gmtime_r((time_t*)&img->scanStart, &tm)!=NULL) {
    main_header.scan_start_year=tm.tm_year+1900;
    main_header.scan_start_month=tm.tm_mon+1;
    main_header.scan_start_day=tm.tm_mday;
    main_header.scan_start_hour=tm.tm_hour;
    main_header.scan_start_minute=tm.tm_min;
    main_header.scan_start_second=tm.tm_sec;
    if(IMG_TEST>2) {
      printf("  img->scanStart := %ld\n", (long int)img->scanStart);
      printf("  -> tm_year := %d\n", tm.tm_year);
      printf("  -> tm_hour := %d\n", tm.tm_hour);
    } 
  } else {
    main_header.scan_start_year=1900;
    main_header.scan_start_month=1;
    main_header.scan_start_day=1;
    main_header.scan_start_hour=0;
    main_header.scan_start_minute=0;
    main_header.scan_start_second=0;
    if(IMG_TEST>0) printf("invalid scan_start_time in IMG\n");
  }
  main_header.isotope_halflife=img->isotopeHalflife;
  strcpy(main_header.isotope_code, imgIsotope(img));
  strlcpy(main_header.study_name, img->studyNr, 12);
  strcpy(main_header.patient_name, img->patientName);
  strcpy(main_header.patient_id, img->patientID);
  strlcpy(main_header.user_process_code, img->userProcessCode, 10);
  strncpy(main_header.study_description, img->studyDescription, 32);
  if(IMG_TEST) ecat63PrintMainheader(&main_header, stdout);

  /* Allocate memory for matrix data */
  sdata=(short int*)malloc(img->dimx*img->dimy*sizeof(short int));
  if(sdata==NULL) {strcpy(ecat63errmsg, "out of memory"); return(4);}

  /* Open output ECAT file */
  fp=ecat63Create(fname, &main_header);
  if(fp==NULL) {strcpy(ecat63errmsg, "cannot write ECAT file"); return(3);}

  /* Set the subheader contents, as far as possible */
  switch(main_header.file_type) {
    case RAW_DATA:
      scan_header.data_type=main_header.data_type;
      scan_header.dimension_1=img->dimx;
      scan_header.dimension_2=img->dimy;
      scan_header.frame_duration_sec=1;
      scan_header.scale_factor=1.0;
      scan_header.frame_start_time=0;
      scan_header.frame_duration=1000;
      /* Deadtime correction was done when reading, set to 1.0 to prevent
         second correction when reading again. */
      scan_header.loss_correction_fctr=1.0;
      /*if(IMG_TEST) ecat63PrintScanheader(&scan_header);*/
      break;
    case IMAGE_DATA:
      image_header.data_type=main_header.data_type;
      image_header.num_dimensions=2;
      image_header.dimension_1=img->dimx;
      image_header.dimension_2=img->dimy;
      image_header.recon_scale=img->zoom;
      image_header.quant_scale=1.0;
      image_header.slice_width=img->sizez/10.;
      image_header.pixel_size=img->sizex/10.;
      image_header.frame_start_time=0;
      image_header.frame_duration=1000;
      image_header.plane_eff_corr_fctr=1.0;
      image_header.decay_corr_fctr=1.0;
      image_header.loss_corr_fctr=1.0;
      image_header.ecat_calibration_fctr=1.0;
      image_header.well_counter_cal_fctr=1.0;
      image_header.quant_units=main_header.calibration_units;
      /*if(IMG_TEST) ecat63PrintImageheader(&image_header);*/
      break;
  }

  /* Write one matrix at a time */
  n=0;
  for(plane=1; plane<=img->dimz;plane++) for(frame=1;frame<=img->dimt;frame++) {
    /* Scale data to short ints */
    /* Search min and max */
    fmin=fmax=f=img->m[plane-1][0][0][frame-1];
    for(i=0; i<img->dimy; i++) for(j=0; j<img->dimx; j++) {
      f=img->m[plane-1][i][j][frame-1];
      if(f>fmax) fmax=f; else if(f<fmin) fmin=f;
    }
    /* Calculate scaling factor */
    if(fabs(fmin)>fabs(fmax)) g=fabs(fmin); else g=fabs(fmax);
    if(g!=0) scale=32766./g; else scale=1.0;
    /*printf("fmin=%e fmax=%e g=%e scale=%e\n", fmin, fmax, g, scale);*/
    /* Scale matrix data to shorts */
    sptr=sdata;
    for(i=0; i<img->dimy; i++) for(j=0; j<img->dimx; j++) {
      *sptr=(short int)temp_roundf(scale*img->m[plane-1][i][j][frame-1]);
      sptr++;
    }
    /* Calculate and set subheader min&max and scale */
    smin=(short int)temp_roundf(scale*fmin);
    smax=(short int)temp_roundf(scale*fmax);
    if(main_header.file_type==RAW_DATA) {
      scan_header.scan_min=smin; scan_header.scan_max=smax;
      scan_header.scale_factor=1.0/scale;
      scan_header.frame_start_time=(int)temp_roundf(1000.*img->start[frame-1]);
      scan_header.frame_duration=
        (int)temp_roundf(1000.*(img->end[frame-1]-img->start[frame-1]));
      scan_header.sample_distance=(img->sampleDistance)/10.0;
      scan_header.prompts=temp_roundf(img->prompts[frame-1]);
      scan_header.delayed=temp_roundf(img->randoms[frame-1]);
    } else if(main_header.file_type==IMAGE_DATA) {
      image_header.image_min=smin; image_header.image_max=smax;
      image_header.quant_scale=1.0/scale;
      image_header.frame_start_time=(int)temp_roundf(1000.*img->start[frame-1]);
      image_header.frame_duration=
        (int)temp_roundf(1000.*(img->end[frame-1]-img->start[frame-1]));
      /* Set decay correction factor */
      if(img->decayCorrection==IMG_DC_CORRECTED)
        image_header.decay_corr_fctr=img->decayCorrFactor[frame-1];
      else
        image_header.decay_corr_fctr=0.0;
    } 
    /* Write matrix data */
    m=mat_numcod(frame, img->planeNumber[plane-1], 1, 0, 0);
    /*m=mat_numcod(frame, plane, 1, 0, 0);*/
    sptr=sdata;
    if(IMG_TEST) printf("  writing matnum=%d\n", m);
    if(main_header.file_type==RAW_DATA) {
      if(IMG_TEST) ecat63PrintScanheader(&scan_header, stdout);
      ret=ecat63WriteScan(fp, m, &scan_header, sptr);
    } else if(main_header.file_type==IMAGE_DATA) {
      if(IMG_TEST) ecat63PrintImageheader(&image_header, stdout);
      ret=ecat63WriteImage(fp, m, &image_header, sptr);
    }
    if(ret) {
      sprintf(ecat63errmsg, "cannot write data on pl%02d fr%02d (%d).",
        plane, frame, ret);
      fclose(fp); remove(fname); free(sdata);
      return(9);
    }
    n++;
  } /* next matrix */
  if(IMG_TEST) printf("  %d matrices written in %s\n", n, fname);

  /* Close file and free memory */
  fclose(fp); free(sdata);

  return(0);
}
/*****************************************************************************/

/*****************************************************************************/
/*!
 * Reads one CTI ECAT 6.3 plane (all frames or gates) at a time to memory.
 * Img data must be initialized before this procedure.
 * Existing img->_dataType is not changed.
 * If img data structure is empty, reads the first plane.
 * If img data structure contains data, reads the next plane.
 * Any existing data in img is cleared and replaced by the new plane.
 *
 * @param fname name of the input ECAT 6.3 file
 * @param img data structure in which the file is read
 * @return 0 if ok, 1 next plane was requested but not found any more,
 * 2 invalid input data, 3 failed to open file,
 * 4 failed to read main header, 5 failed to read matrix list,
 * 6 invalid matrix data, 7 failed to read matrix sub header,
 * 8 failed to allocate memory, 9 failed to read matrix data
 */
int ecat63ReadPlaneToImg(const char *fname, IMG *img) {
  int                 i, j, m, ret, blkNr=0, dim_x=0, dim_y=0, del_nr=0;
  int                 frameNr, planeNr, datatype=0, firstm=0, isFirst=0;
  int                 frame, plane, prev_frame, prev_plane, prev_frameNr=0;
  FILE               *fp;
  ECAT63_mainheader   main_header;
  MATRIXLIST          mlist;
  MatDir              tmpMatdir;
  Matval              matval;
  ECAT63_imageheader  image_header;
  ECAT63_scanheader   scan_header;
  ECAT63_attnheader   attn_header;
  ECAT63_normheader   norm_header;
  float               scale=1.0;
  short int          *sptr;
  char               *mdata=NULL, *mptr;
  int                *iptr;


  if(IMG_TEST) printf("ecat63ReadPlaneToImg(%s, *img)\n", fname);
  /* Check the arguments */
  if(fname==NULL) {strcpy(ecat63errmsg, "invalid ECAT filename"); return(2);}
  if(img==NULL || img->status==IMG_STATUS_UNINITIALIZED) {
    strcpy(ecat63errmsg, "image data not initialized"); return(2);}

  /* Open input CTI file for read */
  if((fp=fopen(fname, "rb")) == NULL) {
    strcpy(ecat63errmsg, "cannot open ECAT file"); return(3);}

  /* Read main header */
  if((ret=ecat63ReadMainheader(fp, &main_header))) {
    sprintf(ecat63errmsg, "cannot read main header (%d)", ret);
    fclose(fp); return(4);
  }
  if(IMG_TEST) ecat63PrintMainheader(&main_header, stdout);

  /* Read matrix list and nr */
  ecat63InitMatlist(&mlist);
  ret=ecat63ReadMatlist(fp, &mlist, IMG_TEST);
  if(ret) {
    sprintf(ecat63errmsg, "cannot read matrix list (%d)", ret);
    fclose(fp); return(5);
  }
  if(mlist.matrixNr<=0) {
    strcpy(ecat63errmsg, "matrix list is empty"); fclose(fp); return(5);}
  /* Sort matrix list */
  for(i=0; i<mlist.matrixNr-1; i++) for(j=i+1; j<mlist.matrixNr; j++) {
    if(mlist.matdir[i].matnum>mlist.matdir[j].matnum) {
      tmpMatdir=mlist.matdir[i];
      mlist.matdir[i]=mlist.matdir[j]; mlist.matdir[j]=tmpMatdir;
    }
  }
  /* Trim extra frames */
  if(main_header.num_frames>0) {
    del_nr=ecat63DeleteLateFrames(&mlist, main_header.num_frames);
    if(IMG_TEST && del_nr>0)
      printf("  %d entries in matrix list will not be used.\n", del_nr);
  }

  /* Decide the plane to read */
  if(img->status==IMG_STATUS_OCCUPIED) {
    /* img contains data */
    isFirst=0; prev_frameNr=img->dimt;
    /* get the current plane in there */
    prev_plane=img->planeNumber[img->dimz-1];
    /* Clear all data in img: not here but only after finding new data */
  } else {
    /* img does not contain any data */
    isFirst=1;
    /* set "current plane" to -1 */
    prev_plane=-1;
  }
  /* from sorted matrix list, find first plane larger than the current one */
  for(m=0, plane=-1; m<mlist.matrixNr; m++) if(mlist.matdir[m].matstat==1) {
    mat_numdoc(mlist.matdir[m].matnum, &matval);
    if(matval.plane>prev_plane) {plane=matval.plane; break;}
  }
  /* If not found, return an error or 1 if this was not the first one */
  if(plane<0) {
    fclose(fp); ecat63EmptyMatlist(&mlist);
    if(isFirst) {
      sprintf(ecat63errmsg, "ECAT file contains no matrices");
      return(7);
    } else {
      sprintf(ecat63errmsg, "ECAT file contains no more planes");
      if(IMG_TEST) printf("%s\n", ecat63errmsg);
      return(1);
    }
  }
  if(IMG_TEST) printf("Next plane to read: %d\n", plane);
  /* Clear all data in img */
  imgEmpty(img);

  /* Calculate the number of frames and gates */
  prev_frame=frame=-1; frameNr=0; ret=0; planeNr=1;
  for(m=0; m<mlist.matrixNr; m++) if(mlist.matdir[m].matstat==1) {
    /* get frame and plane; work only with current plane */
    mat_numdoc(mlist.matdir[m].matnum, &matval);
    if(matval.plane<plane) continue; else if(matval.plane>plane) break;
    if(main_header.num_frames>=main_header.num_gates) frame=matval.frame;
    else frame=matval.gate;
    frameNr++;
    /* frame (gate) numbers must be continuous */
    if(prev_frame>0 && frame!=prev_frame+1) {ret=1; break;}
    prev_frame=frame;
    /* Calculate and check the size of one data matrix */
    if(frameNr==1) {
      firstm=m;
      blkNr=mlist.matdir[m].endblk-mlist.matdir[m].strtblk;
    } else if(blkNr!=mlist.matdir[m].endblk-mlist.matdir[m].strtblk) {
      ret=2; break;
    }
    prev_frame=frame;
  } /* next matrix */
  if(ret==1) {
    strcpy(ecat63errmsg, "missing matrix");
    ecat63EmptyMatlist(&mlist); fclose(fp); return(6);
  }
  if(ret==2) {
    strcpy(ecat63errmsg, "matrix sizes are different");
    ecat63EmptyMatlist(&mlist); fclose(fp); return(6);
  }
  /* Check that frameNr is equal to the dimt in IMG */
  if(!isFirst && frameNr!=prev_frameNr) {
    strcpy(ecat63errmsg, "different frame nr in different planes");
    ecat63EmptyMatlist(&mlist); fclose(fp); return(6);
  }

  /* Read x,y-dimensions from 1st matrix subheader on current plane */
  m=firstm; if(main_header.file_type==IMAGE_DATA) {
    ret=ecat63ReadImageheader(fp, mlist.matdir[m].strtblk, &image_header, IMG_TEST-10, NULL);
    dim_x=image_header.dimension_1; dim_y=image_header.dimension_2;
  } else if(main_header.file_type==RAW_DATA) {
    ret=ecat63ReadScanheader(fp, mlist.matdir[m].strtblk, &scan_header, IMG_TEST-10, NULL);
    dim_x=scan_header.dimension_1; dim_y=scan_header.dimension_2;
  } else if(main_header.file_type==ATTN_DATA) {
    ret=ecat63ReadAttnheader(fp, mlist.matdir[m].strtblk, &attn_header, IMG_TEST-10, NULL);
    dim_x=attn_header.dimension_1; dim_y=attn_header.dimension_2;
  } else if(main_header.file_type==NORM_DATA) {
    ret=ecat63ReadNormheader(fp, mlist.matdir[m].strtblk, &norm_header, IMG_TEST-10, NULL);
    dim_x=norm_header.dimension_1; dim_y=norm_header.dimension_2;
  }
  /*pxlNr=dim_x*dim_y;*/
  if(ret) {
    sprintf(ecat63errmsg, "cannot read matrix %u subheader",
            mlist.matdir[m].matnum);
    ecat63EmptyMatlist(&mlist); fclose(fp); return(7);
  }

  /* Allocate memory for ECAT data matrix */
  if(IMG_TEST) printf("allocating memory for %d blocks\n", blkNr);
  mdata=(char*)malloc(blkNr*MatBLKSIZE); if(mdata==NULL) {
    strcpy(ecat63errmsg, "out of memory");
    fclose(fp); ecat63EmptyMatlist(&mlist); return(8);
  }
  /* Allocate memory for whole img data */
  ret=imgAllocate(img, planeNr, dim_y, dim_x, frameNr);
  if(ret) {
    sprintf(ecat63errmsg, "out of memory (%d)", ret);
    fclose(fp); ecat63EmptyMatlist(&mlist); free(mdata); return(8);
  }

  /* Fill img info with ECAT main and sub header information */
  img->scanner=main_header.system_type;
  img->unit=main_header.calibration_units; /* this continues below */
  strlcpy(img->radiopharmaceutical, main_header.radiopharmaceutical, 32);
  img->isotopeHalflife=main_header.isotope_halflife;
  img->scanStart=ecat63Scanstarttime(&main_header);
  if(img->scanStart==-1) {
    img->scanStart=0;
    if(IMG_TEST>0) printf("invalid scan_start_time in mainheader\n");
  }
  img->axialFOV=10.*main_header.axial_fov;
  img->transaxialFOV=10.*main_header.transaxial_fov;
  strncpy(img->studyNr, main_header.study_name, MAX_STUDYNR_LEN);
  strlcpy(img->patientName, main_header.patient_name, 32);
  strlcpy(img->patientID, main_header.patient_id, 16);
  img->sizez=10.*main_header.plane_separation;
  if(main_header.file_type==IMAGE_DATA) {
    img->type=IMG_TYPE_IMAGE;
    if(img->unit<1) img->unit=image_header.quant_units;
    img->zoom=image_header.recon_scale;
    if(image_header.decay_corr_fctr>1.0) img->decayCorrection=IMG_DC_CORRECTED;
    else img->decayCorrection=IMG_DC_UNKNOWN;
    img->sizex=img->sizey=10.*image_header.pixel_size;
    if(img->sizez<10.*image_header.slice_width)
      img->sizez=10.*image_header.slice_width;
    img->xform[0]=NIFTI_XFORM_UNKNOWN; // qform
    img->xform[1]=NIFTI_XFORM_SCANNER_ANAT; // sform
    img->quatern[6]=img->sizex; img->quatern[9]=img->sizex;
    img->quatern[11]=img->sizey; img->quatern[13]=img->sizey;
    img->quatern[16]=img->sizez; img->quatern[17]=img->sizez;
  } else if(main_header.file_type==RAW_DATA) {
    img->type=IMG_TYPE_RAW;
    img->decayCorrection=IMG_DC_NONCORRECTED;
  } else if(main_header.file_type==ATTN_DATA) {
    img->type=IMG_TYPE_ATTN;
    img->decayCorrection=IMG_DC_NONCORRECTED;
  } else if(main_header.file_type==NORM_DATA) {
    img->type=IMG_TYPE_RAW;
    img->decayCorrection=IMG_DC_NONCORRECTED;
  }
  img->planeNumber[0]=plane;
  strlcpy(img->studyDescription, main_header.study_description, 32);
  strlcpy(img->userProcessCode, main_header.user_process_code, 10);
  //img->userProcessCode[10]=(char)0;
  /* If valid study number is found in user_process_code, then take it */  
  if(!img->studyNr[0] && studynr_validity_check(img->userProcessCode))
    strlcpy(img->studyNr, img->userProcessCode, MAX_STUDYNR_LEN+1);
  /* Set _fileFormat */
  img->_fileFormat=IMG_E63;

  /* Read one ECAT matrix at a time and put them to img */
  for(m=firstm, frame=1; m<mlist.matrixNr; m++) if(mlist.matdir[m].matstat==1) {
    /* get plane */
    mat_numdoc(mlist.matdir[m].matnum, &matval);
    if(matval.plane>plane) break; /* Quit when current plane is read */
    /* Read subheader */
    if(main_header.file_type==IMAGE_DATA) {
      ret=ecat63ReadImageheader(fp, mlist.matdir[m].strtblk, &image_header, IMG_TEST-10, NULL);
      if(ret==0 && (dim_x!=image_header.dimension_1 || dim_y!=image_header.dimension_2)) ret=1;
      scale=image_header.quant_scale;
      if(image_header.ecat_calibration_fctr>0.0 &&
         fabs(main_header.calibration_factor/image_header.ecat_calibration_fctr - 1.0)>0.0001)
        scale*=image_header.ecat_calibration_fctr;
      if(img->unit==0 && image_header.quant_units>0)
        img->unit=image_header.quant_units;
      datatype=image_header.data_type;
      img->start[frame-1]=image_header.frame_start_time/1000.;
      img->end[frame-1]=img->start[frame-1]+image_header.frame_duration/1000.;
      img->mid[frame-1]=0.5*(img->start[frame-1]+img->end[frame-1]);
      if(image_header.decay_corr_fctr>1.0)
        img->decayCorrFactor[frame-1]=image_header.decay_corr_fctr;
      else
        img->decayCorrFactor[frame-1]=0.0;
    } else if(main_header.file_type==RAW_DATA) {
      ret=ecat63ReadScanheader(fp, mlist.matdir[m].strtblk, &scan_header, IMG_TEST-10, NULL);
      if(ret==0 && (dim_x!=scan_header.dimension_1 ||dim_y!=scan_header.dimension_2)) ret=1;
      scale=scan_header.scale_factor;
      if(scan_header.loss_correction_fctr>0.0) // dead-time correction
        scale*=scan_header.loss_correction_fctr;
      datatype=scan_header.data_type;
      img->start[frame-1]=scan_header.frame_start_time/1000.;
      img->end[frame-1]=img->start[frame-1]+scan_header.frame_duration/1000.;
      img->mid[frame-1]=0.5*(img->start[frame-1]+img->end[frame-1]);
      img->sampleDistance=10.0*scan_header.sample_distance;
      img->prompts[frame-1]=(float)scan_header.prompts;
      img->randoms[frame-1]=(float)scan_header.delayed;
    } else if(main_header.file_type==ATTN_DATA) {
      ret=ecat63ReadAttnheader(fp, mlist.matdir[m].strtblk, &attn_header, IMG_TEST-10, NULL);
      if(ret==0 && (dim_x!=attn_header.dimension_1 ||dim_y!=attn_header.dimension_2)) ret=1;
      img->sampleDistance=10.0*attn_header.sample_distance;
      scale=attn_header.scale_factor;
      datatype=attn_header.data_type;
    } else if(main_header.file_type==NORM_DATA) {
      ret=ecat63ReadNormheader(fp, mlist.matdir[m].strtblk, &norm_header, IMG_TEST-10, NULL);
      if(ret==0 && (dim_x!=norm_header.dimension_1 ||dim_y!=norm_header.dimension_2)) ret=1;
      scale=norm_header.scale_factor;
      datatype=norm_header.data_type;
    } else datatype=-1;
    if(ret) {
      sprintf(ecat63errmsg, "cannot read matrix %u subheader",
        mlist.matdir[m].matnum);
      ecat63EmptyMatlist(&mlist); free(mdata); fclose(fp); return(7);
    }
    if(main_header.calibration_factor>0.0)
      scale*=main_header.calibration_factor;
    if(IMG_TEST>2) printf("scale=%e datatype=%d\n", scale, datatype);
    /* Read the pixel data */
    ret=ecat63ReadMatdata(fp, mlist.matdir[m].strtblk+1,
         mlist.matdir[m].endblk-mlist.matdir[m].strtblk,
         mdata, datatype);
    if(ret) {
      strcpy(ecat63errmsg, "cannot read matrix data");
      ecat63EmptyMatlist(&mlist); free(mdata); fclose(fp); return(9);
    }
    if(datatype==BYTE_TYPE) {
      for(i=0, mptr=mdata; i<dim_y; i++) for(j=0; j<dim_x; j++, mptr++) {
        img->m[0][i][j][frame-1]=scale*(float)(*mptr);
      }
    } else if(datatype==VAX_I2 || datatype==SUN_I2) {
      for(i=0, mptr=mdata; i<dim_y; i++) for(j=0; j<dim_x; j++, mptr+=2) {
        sptr=(short int*)mptr;
        img->m[0][i][j][frame-1]=scale*(float)(*sptr);
      }
    } else if(datatype==VAX_I4 || datatype==SUN_I4) {
      for(i=0, mptr=mdata; i<dim_y; i++) for(j=0; j<dim_x; j++, mptr+=4) {
        iptr=(int*)mptr;
        img->m[0][i][j][frame-1]=scale*(float)(*iptr);
      }
    } else if(datatype==VAX_R4 || datatype==IEEE_R4) {
      for(i=0, mptr=mdata; i<dim_y; i++) for(j=0; j<dim_x; j++, mptr+=4) {
        memcpy(&img->m[0][i][j][frame-1], mptr, 4);
        img->m[0][i][j][frame-1]*=scale;
      }
    }
    frame++;
  } /* next matrix (frame) */
  /* Set the unit */
  imgUnitFromEcat(img, img->unit);

  /* Free memory and close file */
  free(mdata);
  ecat63EmptyMatlist(&mlist);
  fclose(fp);

  /* Set _dataType if it has not been set before */
  if(img->_dataType<1) img->_dataType=datatype;
  /* For saving, only 2-byte VAX data types are allowed */
  if(img->_dataType==VAX_I4 || img->_dataType==VAX_R4) img->_dataType=VAX_I2;

  return(0);
}
/*****************************************************************************/

/*****************************************************************************/
/**
   Adds all matrices in memory to the ECAT file.
   If ECAT file does not exist, it is created.

   Please note that existing ECAT file is NOT saved as bak file.
  
   @return 0 if ok, 1 invalid input, 2 image status is not 'occupied', 
   3 failed to open file for reading, 4 failed to allocate memory for data, 
   9 failed to write data, 21 invalid matrix list,
   22 failed to write main header
 */
int ecat63AddImg(
  /** Name of the output ECAT 6.3 file. */
  const char *fname,
  /** Pointer to data structure from which the data is written. */
  IMG *img
) {
  int                 n, i, j, m, ret=0, add=0;
  int                 frameNr, planeNr;
  int                 frame, plane, prev_plane;
  float               f, fmax, fmin, g, scale;
  short int          *sdata, *sptr, smin, smax;
  FILE               *fp;
  ECAT63_mainheader   main_header;
  ECAT63_imageheader  image_header;
  ECAT63_scanheader   scan_header;
  MATRIXLIST          mlist;
  MatDir              tmpMatdir;
  Matval              matval;
  struct tm           tm;

  if(IMG_TEST) printf("ecat63AddImg(%s, *img)\n", fname);
  if(IMG_TEST>4) imgInfo(img);
  /* Check the arguments */
  if(fname==NULL) {strcpy(ecat63errmsg, "invalid ECAT filename"); return(1);}
  if(img==NULL || img->status!=IMG_STATUS_OCCUPIED) {
    strcpy(ecat63errmsg, "image data is empty"); return(2);}
  if(img->_dataType<1) img->_dataType=VAX_I2;

  /* Initiate headers */
  memset(&main_header, 0, sizeof(ECAT63_mainheader));
  memset(&image_header,0, sizeof(ECAT63_imageheader));
  memset(&scan_header, 0, sizeof(ECAT63_scanheader));

  /* Make a main header */
  main_header.sw_version=2;
  main_header.num_planes=img->dimz;  
  main_header.num_frames=img->dimt;  
  main_header.num_gates=1;
  main_header.num_bed_pos=1;
  if(img->type==IMG_TYPE_IMAGE) main_header.file_type=IMAGE_DATA;
  else if(img->type==IMG_TYPE_RAW) main_header.file_type=RAW_DATA;
  else {strcpy(ecat63errmsg, "invalid filetype"); return(1);}
  main_header.data_type=img->_dataType;
  if(img->scanner>0) main_header.system_type=img->scanner;
  else main_header.system_type=ECAT63_SYSTEM_TYPE_DEFAULT;
  main_header.calibration_factor=1.0;
  main_header.axial_fov=img->axialFOV/10.0;
  main_header.transaxial_fov=img->transaxialFOV/10.0;
  main_header.plane_separation=img->sizez/10.0;
  main_header.calibration_units=imgUnitToEcat6(img);
  strncpy(main_header.radiopharmaceutical, img->radiopharmaceutical, 32);
  if(gmtime_r((time_t*)&img->scanStart, &tm)!=NULL) {
    main_header.scan_start_year=tm.tm_year+1900;
    main_header.scan_start_month=tm.tm_mon+1;
    main_header.scan_start_day=tm.tm_mday;
    main_header.scan_start_hour=tm.tm_hour;
    main_header.scan_start_minute=tm.tm_min;
    main_header.scan_start_second=tm.tm_sec;
    if(IMG_TEST>2) {
      printf("  img->scanStart := %ld\n", (long int)img->scanStart);
      printf("  -> tm_year := %d\n", tm.tm_year);
      printf("  -> tm_hour := %d\n", tm.tm_hour);
    } 
  } else {
    main_header.scan_start_year=1900;
    main_header.scan_start_month=1;
    main_header.scan_start_day=1;
    main_header.scan_start_hour=0;
    main_header.scan_start_minute=0;
    main_header.scan_start_second=0;
    if(IMG_TEST>0) printf("invalid scan_start_time in IMG\n");
  }
  main_header.isotope_halflife=img->isotopeHalflife;
  strcpy(main_header.isotope_code, imgIsotope(img));
  strlcpy(main_header.study_name, img->studyNr, 12);
  strcpy(main_header.patient_name, img->patientName);
  strcpy(main_header.patient_id, img->patientID);
  strlcpy(main_header.study_description, img->studyDescription, 32);
  strlcpy(main_header.user_process_code, img->userProcessCode, 10);

  /* Check if the ECAT file exists */
  if(access(fname, 0) != -1) {
    if(IMG_TEST) printf("Opening existing ECAT file.\n");
    add=1;
    /* Open the ECAT file */
    if((fp=fopen(fname, "rb+")) == NULL) {
      strcpy(ecat63errmsg, "cannot create ECAT file"); return(3);}
    /* Read main header */
    if((ret=ecat63ReadMainheader(fp, &main_header))) {
      sprintf(ecat63errmsg, "cannot read main header (%d)", ret);
      fclose(fp); return(3);
    }
    fflush(fp);
    /* Check that filetypes are matching */
    if((main_header.file_type==IMAGE_DATA && img->type==IMG_TYPE_IMAGE) ||
       (main_header.file_type==RAW_DATA && img->type==IMG_TYPE_RAW)) {
      /* ok */
    } else {
      sprintf(ecat63errmsg, "cannot add different filetype");
      fclose(fp); return(3);
    }
  } else {
    if(IMG_TEST) printf("ECAT file does not exist.\n");
    add=0;
    /* Create output ECAT file */
    fp=ecat63Create(fname, &main_header);
    if(fp==NULL) {strcpy(ecat63errmsg, "cannot create ECAT file"); return(3);}
  }
  if(IMG_TEST) ecat63PrintMainheader(&main_header, stdout);

  /* Allocate memory for matrix data */
  sdata=(short int*)malloc(img->dimx*img->dimy*sizeof(short int));
  if(sdata==NULL) {strcpy(ecat63errmsg, "out of memory"); return(4);}

  /* Set the subheader contents, as far as possible */
  switch(main_header.file_type) {
    case RAW_DATA:
      scan_header.data_type=main_header.data_type;
      scan_header.dimension_1=img->dimx;
      scan_header.dimension_2=img->dimy;
      scan_header.frame_duration_sec=1;
      scan_header.scale_factor=1.0;
      scan_header.frame_start_time=0;
      scan_header.frame_duration=1000;
      scan_header.loss_correction_fctr=1.0;
      scan_header.sample_distance=(img->sampleDistance)/10.0;
      /*if(IMG_TEST) ecat63PrintScanheader(&scan_header);*/
      break;
    case IMAGE_DATA:
      image_header.data_type=main_header.data_type;
      image_header.num_dimensions=2;
      image_header.dimension_1=img->dimx;
      image_header.dimension_2=img->dimy;
      image_header.recon_scale=img->zoom;
      image_header.quant_scale=1.0;
      image_header.slice_width=img->sizez/10.;
      image_header.pixel_size=img->sizex/10.;
      image_header.frame_start_time=0;
      image_header.frame_duration=1000;
      image_header.plane_eff_corr_fctr=1.0;
      image_header.decay_corr_fctr=0.0;
      image_header.loss_corr_fctr=1.0;
      image_header.ecat_calibration_fctr=1.0;
      image_header.well_counter_cal_fctr=1.0;
      image_header.quant_units=main_header.calibration_units;
      /*if(IMG_TEST) ecat63PrintImageheader(&image_header);*/
      break;
  }

  /* Write one matrix at a time */
  n=0;
  for(plane=1; plane<=img->dimz;plane++) for(frame=1;frame<=img->dimt;frame++) {
    /* Scale data to short ints */
    /* Search min and max */
    fmin=fmax=f=img->m[plane-1][0][0][frame-1];
    for(i=0; i<img->dimy; i++) for(j=0; j<img->dimx; j++) {
      f=img->m[plane-1][i][j][frame-1];
      if(f>fmax) fmax=f; else if(f<fmin) fmin=f;
    }
    /* Calculate scaling factor */
    if(fabs(fmin)>fabs(fmax)) g=fabs(fmin); else g=fabs(fmax);
    if(g!=0) scale=32766./g; else scale=1.0;
    /*printf("fmin=%e fmax=%e g=%e scale=%e\n", fmin, fmax, g, scale);*/
    /* Scale matrix data to shorts */
    sptr=sdata;
    for(i=0; i<img->dimy; i++) for(j=0; j<img->dimx; j++) {
      *sptr=(short int)temp_roundf(scale*img->m[plane-1][i][j][frame-1]);
      sptr++;
    }
    /* Calculate and set subheader min&max and scale */
    smin=(short int)temp_roundf(scale*fmin);
    smax=(short int)temp_roundf(scale*fmax);
    if(main_header.file_type==RAW_DATA) {
      scan_header.scan_min=smin; scan_header.scan_max=smax;
      scan_header.scale_factor=1.0/scale;
      scan_header.frame_start_time=(int)temp_roundf(1000.*img->start[frame-1]);
      scan_header.frame_duration=
        (int)temp_roundf(1000.*(img->end[frame-1]-img->start[frame-1]));
      scan_header.prompts=temp_roundf(img->prompts[frame-1]);
      scan_header.delayed=temp_roundf(img->randoms[frame-1]);
    } else if(main_header.file_type==IMAGE_DATA) {
      image_header.image_min=smin; image_header.image_max=smax;
      image_header.quant_scale=1.0/scale;
      image_header.frame_start_time=(int)temp_roundf(1000.*img->start[frame-1]);
      image_header.frame_duration=
        (int)temp_roundf(1000.*(img->end[frame-1]-img->start[frame-1]));
      /* Set decay correction factor */
      if(img->decayCorrection==IMG_DC_CORRECTED)
        image_header.decay_corr_fctr=img->decayCorrFactor[frame-1];
      else
        image_header.decay_corr_fctr=0.0;
    } 
    /* Write matrix data */
    m=mat_numcod(frame, img->planeNumber[plane-1], 1, 0, 0);
    sptr=sdata;
    if(IMG_TEST) printf("  writing matnum=%d\n", m);
    if(main_header.file_type==RAW_DATA) {
      /*if(IMG_TEST) ecat63PrintScanheader(&scan_header);*/
      ret=ecat63WriteScan(fp, m, &scan_header, sptr);
    } else if(main_header.file_type==IMAGE_DATA) {
      /*if(IMG_TEST) ecat63PrintImageheader(&image_header);*/
      ret=ecat63WriteImage(fp, m, &image_header, sptr);
    }
    if(ret) {
      sprintf(ecat63errmsg, "cannot write data on pl%02d fr%02d (%d).",
        plane, frame, ret);
      fclose(fp); remove(fname); free(sdata);
      return(9);
    }
    n++;
  } /* next matrix */
  if(IMG_TEST) printf("  %d matrices written in %s\n", n, fname);

  /* Free matrix memory */
  free(sdata);

  /* If matrices were added, set main header */
  if(add==1) {
    fflush(fp);
    /* Read matrix list */
    ecat63InitMatlist(&mlist);
    ret=ecat63ReadMatlist(fp, &mlist, IMG_TEST);
    if(ret) {
      sprintf(ecat63errmsg, "cannot read matrix list (%d)", ret);
      fclose(fp); return(21);
    }
    if(mlist.matrixNr<=0) {
      strcpy(ecat63errmsg, "matrix list is empty"); fclose(fp); return(21);}
    /* Sort matrix list */
    for(i=0; i<mlist.matrixNr-1; i++) for(j=i+1; j<mlist.matrixNr; j++) {
      if(mlist.matdir[i].matnum>mlist.matdir[j].matnum) {
        tmpMatdir=mlist.matdir[i];
        mlist.matdir[i]=mlist.matdir[j]; mlist.matdir[j]=tmpMatdir;
      }
    }
    /* Calculate the number of planes and frames */
    prev_plane=plane=-1; /*prev_frame=frame=-1;*/ frameNr=planeNr=0;
    for(m=0; m<mlist.matrixNr; m++) {
      mat_numdoc(mlist.matdir[m].matnum, &matval);
      plane=matval.plane; frame=matval.frame;
      if(plane!=prev_plane) {frameNr=1; planeNr++;} else {frameNr++;}
      prev_plane=plane; /*prev_frame=frame;*/
    } /* next matrix */
    ecat63EmptyMatlist(&mlist);
    main_header.num_planes=planeNr; main_header.num_frames=frameNr;
    /* Write main header */
    ret=ecat63WriteMainheader(fp, &main_header);
    if(ret) {
      sprintf(ecat63errmsg, "cannot write mainheader (%d)", ret);
      fclose(fp); return(22);
    }
  }

  /* Close file and free memory */
  fclose(fp);

  return(0);
}
/*****************************************************************************/

/*****************************************************************************/
/*!
 * Check whether read functions in IMG library support this ECAT 6.3 file_type.
 *
 * @param h Ecat 6.3 main header
 * @return 1 if supported, 0 if not.
 */
int imgEcat63Supported(ECAT63_mainheader *h) {
  if(h==NULL) return(0);
  if(h->file_type==IMAGE_DATA) return(1);
  if(h->file_type==RAW_DATA)   return(1);
  if(h->file_type==ATTN_DATA)  return(1);
  if(h->file_type==NORM_DATA)  return(1);
  return(0);
}
/*****************************************************************************/

/*****************************************************************************/
/*!
 * Copy ECAT 6.3 main header information into IMG
 *
 * @param img target image structure
 * @param h source Ecat 6.3 main header
 */
void imgGetEcat63MHeader(IMG *img, ECAT63_mainheader *h) 
{
  if(IMG_TEST>0) printf("imgGetEcat63MHeader()\n");
  img->_dataType=h->data_type; /* again in subheaders*/
  /* For saving IMG data, only 2-byte VAX data types are allowed,
     so change it now */
  if(img->_dataType==VAX_I4 || img->_dataType==VAX_R4) img->_dataType=VAX_I2;
  img->scanner=h->system_type;
  strlcpy(img->radiopharmaceutical, h->radiopharmaceutical, 32);
  img->isotopeHalflife=h->isotope_halflife;
  img->scanStart=ecat63Scanstarttime(h); 
  if(img->scanStart==-1) {
    img->scanStart=0;
    if(IMG_TEST>0) printf("invalid scan_start_time in mainheader\n");
  }
  if(IMG_TEST>1) {
    char buf1[32], buf2[32];
    printf(" %s -> %s\n", ecat63ScanstarttimeInt(h, buf1), 
                          ctime_r_int(&img->scanStart, buf2));
  }
  img->axialFOV=10.0*h->axial_fov;
  img->transaxialFOV=10.0*h->transaxial_fov;
  strlcpy(img->studyNr, h->study_name, MAX_STUDYNR_LEN+1);
  strlcpy(img->patientName, h->patient_name, 32);
  strlcpy(img->patientID, h->patient_id, 16);
  img->sizez=10.0*h->plane_separation;
  switch(h->file_type) {
    case IMAGE_DATA: img->type=IMG_TYPE_IMAGE; break;
    case RAW_DATA:
    case ATTN_DATA:
    case NORM_DATA: 
    default: img->type=IMG_TYPE_RAW;
  }
  strlcpy(img->studyDescription, h->study_description, 32);
  strncpy(img->userProcessCode, h->user_process_code, 10);
  img->userProcessCode[10]=(char)0;
  /* If valid study number is found in user_process_code, then take it */  
  if(!img->studyNr[0] && studynr_validity_check(img->userProcessCode))
    strlcpy(img->studyNr, img->userProcessCode, MAX_STUDYNR_LEN+1);
  /* Set calibration unit (this may have to be read from subheader later) */
  imgUnitFromEcat(img, h->calibration_units);
}
/*****************************************************************************/

/*****************************************************************************/
/*!
 * Copy information from IMG struct into ECAT 6.3 main header
 *
 * @param img source image structure
 * @param h target Ecat 6.3 main header
 */
void imgSetEcat63MHeader(IMG *img, ECAT63_mainheader *h) 
{
  struct tm tm;
  memset(&tm, 0, sizeof(struct tm));

  if(IMG_TEST>0) printf("imgSetEcat63MHeader()\n");
  if(IMG_TEST>2) {
    char buf[32];
    if(!ctime_r_int(&img->scanStart, buf)) strcpy(buf, "1900-01-01 00:00:00");
    fprintf(stdout, "  scan_start_time := %s\n", buf);
  }
  h->sw_version=2;
  h->num_planes=img->dimz;  
  h->num_frames=img->dimt;  
  h->num_gates=1;
  h->num_bed_pos=1;
  if(img->type==IMG_TYPE_IMAGE) h->file_type=IMAGE_DATA;
  else h->file_type=RAW_DATA;
  h->data_type=VAX_I2;
  if(img->scanner>0) h->system_type=img->scanner;
  else h->system_type=ECAT63_SYSTEM_TYPE_DEFAULT;
  h->calibration_factor=1.0;
  h->axial_fov=img->axialFOV/10.0;
  h->transaxial_fov=img->transaxialFOV/10.0;
  h->plane_separation=img->sizez/10.0;
  h->calibration_units=imgUnitToEcat6(img);
  strncpy(h->radiopharmaceutical, img->radiopharmaceutical, 32);
  if(gmtime_r(&img->scanStart, &tm)!=NULL) {
    h->scan_start_year=tm.tm_year+1900;
    h->scan_start_month=tm.tm_mon+1;
    h->scan_start_day=tm.tm_mday;
    h->scan_start_hour=tm.tm_hour;
    h->scan_start_minute=tm.tm_min;
    h->scan_start_second=tm.tm_sec;
    if(IMG_TEST>2) {
      printf("  img->scanStart := %ld\n", (long int)img->scanStart);
      printf("  -> tm_year := %d\n", tm.tm_year);
      printf("  -> tm_hour := %d\n", tm.tm_hour);
    } 
  } else {
    h->scan_start_year=1900;
    h->scan_start_month=1;
    h->scan_start_day=1;
    h->scan_start_hour=0;
    h->scan_start_minute=0;
    h->scan_start_second=0;
    if(IMG_TEST>0) printf("invalid scan_start_time in IMG\n");
  }
  h->isotope_halflife=img->isotopeHalflife;
  strcpy(h->isotope_code, imgIsotope(img));
  strlcpy(h->study_name, img->studyNr, 12);
  strcpy(h->patient_name, img->patientName);
  strcpy(h->patient_id, img->patientID);
  strlcpy(h->user_process_code, img->userProcessCode, 10);
  strlcpy(h->study_description, img->studyDescription, 32);
}
/*****************************************************************************/

/*****************************************************************************/
/*!
 * Return the IMG fileformat based on ECAT 6.3 file_type.
 *
 * @param h Ecat 6.3 main header
 * @return IMG._fileFormat value.
 */
int imgGetEcat63Fileformat(ECAT63_mainheader *h) {
  int fileFormat=IMG_UNKNOWN;
  switch(h->file_type) {
    case IMAGE_DATA:
    case RAW_DATA:
    case ATTN_DATA:
    case NORM_DATA:
      fileFormat=IMG_E63; break;
    default:
      fileFormat=IMG_UNKNOWN; break;
  }
  return fileFormat;
}
/*****************************************************************************/

/*****************************************************************************/
/**
   Fill IMG struct header information from an image or sinogram file
   in ECAT 6.3 format. Information concerning separate frames or planes is not
   filled.
 
   @note ECAT 6.3 files do not have a magic number, therefore, do
    not use this function to determine if your file is in this format, at least
    test all other possible formats before calling this.
   
   @param fname image or sinogram filename
   @param img pointer to initialized IMG structure
   @return errstatus, which is STATUS_OK (0) when call was successful,
   and >0 in case of an error.
 */
int imgReadEcat63Header(
  const char *fname,
  IMG *img
) {
  int                 m, blkNr=0, ret=0;
  int                 frameNr, planeNr /*, del_nr=0*/;
  FILE               *fp;
  ECAT63_mainheader   main_header;
  MATRIXLIST          mlist;
  ECAT63_imageheader  image_header;
  ECAT63_scanheader   scan_header;
  ECAT63_attnheader   attn_header;
  ECAT63_normheader   norm_header;


  if(IMG_TEST>0) printf("\nimgReadEcat63Header(%s, *img)\n", fname);
  
  /* Check the arguments */
  if(img==NULL) return STATUS_FAULT;
  if(img->status!=IMG_STATUS_INITIALIZED) return STATUS_FAULT;
  imgSetStatus(img, STATUS_FAULT);
  if(fname==NULL) return STATUS_FAULT;
    
  /* Open the file */
  if((fp=fopen(fname, "rb")) == NULL) return STATUS_NOFILE;
  
  /* Read main header */
  if((ret=ecat63ReadMainheader(fp, &main_header))) {
    fclose(fp); return STATUS_NOMAINHEADER;}
  /* Check if file_type is supported */
  if(imgEcat63Supported(&main_header)==0) {
    fclose(fp); return STATUS_UNSUPPORTED;}
  /* Copy main header information into IMG; sets also img.type */
  imgGetEcat63MHeader(img, &main_header);
  if(IMG_TEST>7) printf("img.type := %d\n", img->type);
  img->_fileFormat=imgGetEcat63Fileformat(&main_header);
  if(IMG_TEST>7) printf("img._fileFormat := %d\n", img->_fileFormat);
  if(img->_fileFormat==IMG_UNKNOWN) {fclose(fp); return STATUS_UNSUPPORTED;}

  /* Read matrix list and nr and sort it */
  ecat63InitMatlist(&mlist);
  ret=ecat63ReadMatlist(fp, &mlist, IMG_TEST-1);
  if(ret) {fclose(fp); return STATUS_NOMATLIST;}
  if(mlist.matrixNr<=0) {
    ecat63EmptyMatlist(&mlist); fclose(fp); return STATUS_INVALIDMATLIST;}
  /* Make sure that plane and frame numbers are continuous */
  ecat63GatherMatlist(&mlist, 1, 1, 1, 1);
  ecat63SortMatlistByPlane(&mlist); // otherwise frameNr can't be computed below
  /* Trim extra frames */
  if(main_header.num_frames>0)
    /*del_nr=*/ecat63DeleteLateFrames(&mlist, main_header.num_frames);
  /* Get plane and frame numbers and check that volume is full */
  ret=ecat63GetPlaneAndFrameNr(&mlist, &main_header, &planeNr, &frameNr);
  if(ret) {ecat63EmptyMatlist(&mlist); fclose(fp); return ret;}
  img->dimz=planeNr;
  img->dimt=frameNr;
  /* Calculate and check the size of one data matrix */
  ret=ecat63GetMatrixBlockSize(&mlist, &blkNr);
  if(ret) {ecat63EmptyMatlist(&mlist); fclose(fp); return ret;}

  /* Read one subheader */
  if(IMG_TEST>5) printf("main_header.file_type := %d\n", main_header.file_type);
  m=0;
  switch(main_header.file_type) {
    case IMAGE_DATA:
      ret=ecat63ReadImageheader(fp, mlist.matdir[m].strtblk, &image_header, IMG_TEST-10, NULL);
      break;
    case RAW_DATA:
      ret=ecat63ReadScanheader(fp, mlist.matdir[m].strtblk, &scan_header, IMG_TEST-10, NULL);
      break;
    case ATTN_DATA:
      ret=ecat63ReadAttnheader(fp, mlist.matdir[m].strtblk, &attn_header, IMG_TEST-10, NULL);
      break;
    case NORM_DATA:
      ret=ecat63ReadNormheader(fp, mlist.matdir[m].strtblk, &norm_header, IMG_TEST-10, NULL);
      break;
    default: ret=-1;
  }
  /* Free locally allocated memory and close the file */
  ecat63EmptyMatlist(&mlist); fclose(fp);
  /* Check whether subheader was read */
  if(ret) return STATUS_NOSUBHEADER;

  /* Get the following information in the subheader:
     dimensions x, y and z; datatype;
     image decay correction on/off, zoom, and pixel size;
     sinogram sample distance.
   */
  switch(main_header.file_type) {
    case IMAGE_DATA:
      img->dimx=image_header.dimension_1; img->dimy=image_header.dimension_2;
      if(img->unit<1) imgUnitFromEcat(img, image_header.quant_units);
      img->_dataType=image_header.data_type;
      img->zoom=image_header.recon_scale;
      if(image_header.decay_corr_fctr>1.0) img->decayCorrection=IMG_DC_CORRECTED;
      else img->decayCorrection=IMG_DC_UNKNOWN;
      img->sizex=img->sizey=10.*image_header.pixel_size;
      if(img->sizez<10.*image_header.slice_width)
        img->sizez=10.*image_header.slice_width;
      break;
    case RAW_DATA:
      img->dimx=scan_header.dimension_1; img->dimy=scan_header.dimension_2;
      img->type=IMG_TYPE_RAW;
      img->_dataType=scan_header.data_type;
      img->decayCorrection=IMG_DC_NONCORRECTED;
      img->sampleDistance=10.0*scan_header.sample_distance;
      break;
    case ATTN_DATA:
      img->dimx=attn_header.dimension_1; img->dimy=attn_header.dimension_2;
      img->type=IMG_TYPE_ATTN;
      img->decayCorrection=IMG_DC_NONCORRECTED;
      img->_dataType=attn_header.data_type;
      break;
    case NORM_DATA:
      img->dimx=norm_header.dimension_1; img->dimy=norm_header.dimension_2;
      img->type=IMG_TYPE_RAW;
      img->decayCorrection=IMG_DC_NONCORRECTED;
      img->_dataType=norm_header.data_type;
      break;
  }

  /* For saving IMG data, only 2-byte VAX data types are allowed,
     so change it now */
  if(img->_dataType==VAX_I4 || img->_dataType==VAX_R4) img->_dataType=VAX_I2;

  imgSetStatus(img, STATUS_OK);
  return STATUS_OK;
}
/*****************************************************************************/

/*****************************************************************************/
/*!
 * Read the first frame from an ECAT 6.3 file into IMG data structure.
 *
 * @param fname name of file from which IMG contents will be read
 * @param img pointer to the initiated but not preallocated IMG data
 * @return errstatus, which is STATUS_OK (0) when call was successful,
 * and >0 in case of an error.
 */
int imgReadEcat63FirstFrame(const char *fname, IMG *img) {
  int ret=0;

  if(IMG_TEST) printf("\nimgReadEcat63FirstFrame(%s, *img)\n", fname);
  /* Check the input */
  if(img==NULL) return STATUS_FAULT;
  if(img->status!=IMG_STATUS_INITIALIZED) return STATUS_FAULT;
  imgSetStatus(img, STATUS_FAULT);
  if(fname==NULL) return STATUS_FAULT;

  /* Read header information from file */
  ret=imgReadEcat63Header(fname, img); if(ret) return(ret);
  if(IMG_TEST>4) imgInfo(img);

  /* Allocate memory for one frame */
  img->dimt=1;
  ret=imgAllocate(img, img->dimz, img->dimy, img->dimx, img->dimt);
  if(ret) return STATUS_NOMEMORY;

  /* Read the first frame */
  ret=imgReadEcat63Frame(fname, 1, img, 0); if(ret) return(ret); 

  /* All went well */
  imgSetStatus(img, STATUS_OK);
  return STATUS_OK;
}
/*****************************************************************************/

/*****************************************************************************/
/**
   Read a specified frame from an ECAT 6.3 file into preallocated IMG
   data structure.

   IMG header is assumed to be filled correctly before
   calling this function, except for information concerning separate planes
   and this frame, which is filled here.
 
  @param fname name of file from which IMG contents will be read
  @param frame_to_read frame which will be read (1..frameNr)
  @param img pointer to the IMG data. Place for the frame must be preallocated
  @param frame_index IMG frame index (0..dimt-1) where data will be placed

  @return errstatus, which is STATUS_OK (0) when call was successful,
  and >0 in case of an error. If frame does not exist, then and only then
  STATUS_NOMATRIX is returned.
 */
int imgReadEcat63Frame(
  const char *fname, int frame_to_read, IMG *img, int frame_index
) {
  int                 m, ret=0, blkNr=0, frame, plane, seqplane=-1, xi, yi;
  int                 local_data_type=0;
  FILE               *fp;
  ECAT63_mainheader   main_header;
  MATRIXLIST          mlist;
  Matval              matval;
  ECAT63_imageheader  image_header;
  ECAT63_scanheader   scan_header;
  ECAT63_attnheader   attn_header;
  ECAT63_normheader   norm_header;
  float               scale=1.0;
  short int          *sptr;
  char               *mdata=NULL, *mptr;
  /*int                *iptr;*/


  if(IMG_TEST) printf("\nimgReadEcat63Frame(%s, %d, *img, %d)\n",
    fname, frame_to_read, frame_index);
    
  /* Check the input */
  if(img==NULL) return STATUS_FAULT;
  if(img->status!=IMG_STATUS_OCCUPIED) return STATUS_FAULT;
  imgSetStatus(img, STATUS_FAULT);
  if(fname==NULL) return STATUS_FAULT;
  if(frame_index<0 || frame_index>img->dimt-1) return STATUS_FAULT;
  if(frame_to_read<1) return STATUS_FAULT;

  /* Open the file */
  if((fp=fopen(fname, "rb")) == NULL) return STATUS_NOFILE;
  
  /* Read main header */
  if((ret=ecat63ReadMainheader(fp, &main_header))) {
    fclose(fp); return STATUS_NOMAINHEADER;}

  /* Read matrix list and nr and sort it */
  ecat63InitMatlist(&mlist);
  ret=ecat63ReadMatlist(fp, &mlist, IMG_TEST-1);
  if(ret) {fclose(fp); return STATUS_NOMATLIST;}
  if(mlist.matrixNr<=0) {
    fclose(fp); ecat63EmptyMatlist(&mlist); return STATUS_INVALIDMATLIST;}
  /* Make sure that plane and frame numbers are continuous */
  ecat63GatherMatlist(&mlist, 1, 1, 1, 1);
  ecat63SortMatlistByFrame(&mlist);
  
  /* Calculate and check the size of one data matrix */
  ret=ecat63GetMatrixBlockSize(&mlist, &blkNr);
  if(ret) {ecat63EmptyMatlist(&mlist); fclose(fp); return ret;}
  /* And allocate memory for the raw data blocks */
  if(IMG_TEST>6) printf("allocating memory for %d blocks\n", blkNr);
  mdata=(char*)malloc(blkNr*MatBLKSIZE);
  if(mdata==NULL) {
    fclose(fp); ecat63EmptyMatlist(&mlist); return STATUS_NOMEMORY;}

  /* Read all matrices that belong to the required frame */
  ret=0; seqplane=-1;
  for(m=0; m<mlist.matrixNr; m++) if(mlist.matdir[m].matstat==1) {
    /* get frame and plane */
    mat_numdoc(mlist.matdir[m].matnum, &matval);
    plane=matval.plane;
    if(main_header.num_frames>=main_header.num_gates) frame=matval.frame;
    else frame=matval.gate; /* printf("frame=%d plane=%d\n", frame, plane); */
    if(frame!=frame_to_read) continue; /* do not process other frames */
    seqplane++; 
    if(img->planeNumber[seqplane]<1) {
      img->planeNumber[seqplane]=plane;
    } else if(img->planeNumber[seqplane]!=plane) {
      fclose(fp); ecat63EmptyMatlist(&mlist); free(mdata);
      return STATUS_MISSINGMATRIX;    
    }
    
    /* Read the subheader */
    if(IMG_TEST>4) printf("reading subheader for matrix %d,%d\n", frame, plane);
    if(main_header.file_type==IMAGE_DATA) {
      ret=ecat63ReadImageheader(fp, mlist.matdir[m].strtblk, &image_header, IMG_TEST-10, NULL);
      scale=image_header.quant_scale;
      if(image_header.ecat_calibration_fctr>0.0 && 
         fabs(main_header.calibration_factor/image_header.ecat_calibration_fctr - 1.0)>0.0001)
        scale*=image_header.ecat_calibration_fctr;
      local_data_type=image_header.data_type;
      img->start[frame_index]=image_header.frame_start_time/1000.;
      img->end[frame_index]=img->start[frame_index]+image_header.frame_duration/1000.;
      img->mid[frame_index]=0.5*(img->start[frame_index]+img->end[frame_index]);
      if(image_header.decay_corr_fctr>1.0) {
        img->decayCorrFactor[frame_index]=image_header.decay_corr_fctr;
        img->decayCorrection=IMG_DC_CORRECTED;
      } else {
        img->decayCorrFactor[frame_index]=0.0;
        img->decayCorrection=IMG_DC_UNKNOWN;
      }
      img->xform[0]=NIFTI_XFORM_UNKNOWN; // qform
      img->xform[1]=NIFTI_XFORM_SCANNER_ANAT; // sform
      img->quatern[6]=img->sizex; img->quatern[9]=img->sizex;
      img->quatern[11]=img->sizey; img->quatern[13]=img->sizey;
      img->quatern[16]=img->sizez; img->quatern[17]=img->sizez;
    } else if(main_header.file_type==RAW_DATA) {
      ret=ecat63ReadScanheader(fp, mlist.matdir[m].strtblk, &scan_header, IMG_TEST-10, NULL);
      scale=scan_header.scale_factor;
      if(scan_header.loss_correction_fctr>0.0) scale*=scan_header.loss_correction_fctr;
      local_data_type=scan_header.data_type;
      img->start[frame_index]=scan_header.frame_start_time/1000.;
      img->end[frame_index]=img->start[frame_index]+scan_header.frame_duration/1000.;
      img->mid[frame_index]=0.5*(img->start[frame_index]+img->end[frame_index]);
      img->sampleDistance=10.0*scan_header.sample_distance;
      img->prompts[frame_index]=(float)scan_header.prompts;
      img->randoms[frame_index]=(float)scan_header.delayed;
    } else if(main_header.file_type==ATTN_DATA) {
      ret=ecat63ReadAttnheader(fp, mlist.matdir[m].strtblk, &attn_header, IMG_TEST-10, NULL);
      img->sampleDistance=10.0*attn_header.sample_distance;
      scale=attn_header.scale_factor;
      local_data_type=attn_header.data_type;
    } else if(main_header.file_type==NORM_DATA) {
      ret=ecat63ReadNormheader(fp, mlist.matdir[m].strtblk, &norm_header, IMG_TEST-10, NULL);
      scale=norm_header.scale_factor;
      local_data_type=norm_header.data_type;
    } else
      img->_dataType=-1;
    if(ret) {
      ecat63EmptyMatlist(&mlist); free(mdata); fclose(fp);
      return STATUS_NOSUBHEADER;
    }
    if(main_header.calibration_factor>0.0) scale*=main_header.calibration_factor;
    if(IMG_TEST>6) printf("scale=%e datatype=%d\n", scale, img->_dataType);
    /* Read the pixel data */
    if(IMG_TEST>4) printf("reading matrix data\n");
    ret=ecat63ReadMatdata(fp, mlist.matdir[m].strtblk+1,
         mlist.matdir[m].endblk-mlist.matdir[m].strtblk,
         mdata, local_data_type);
    if(ret) {
      ecat63EmptyMatlist(&mlist); free(mdata); fclose(fp);
      return STATUS_MISSINGMATRIX;
    }
    if(IMG_TEST>4) printf("converting matrix data to floats\n");
    if(local_data_type==BYTE_TYPE) {
      for(yi=0, mptr=mdata; yi<img->dimy; yi++)
        for(xi=0; xi<img->dimx; xi++, mptr++) {
          img->m[seqplane][yi][xi][frame_index]=scale*(float)(*mptr);
        }
    } else if(local_data_type==VAX_I2 || local_data_type==SUN_I2) {
      for(yi=0, mptr=mdata; yi<img->dimy; yi++)
        for(xi=0; xi<img->dimx; xi++, mptr+=2) {
          sptr=(short int*)mptr;
          img->m[seqplane][yi][xi][frame_index]=scale*(float)(*sptr);
        }
    } else if(local_data_type==VAX_I4 || local_data_type==SUN_I4) {
      for(yi=0, mptr=mdata; yi<img->dimy; yi++)
        for(xi=0; xi<img->dimx; xi++, mptr+=4) {
          /*iptr=(int*)mptr;*/
          img->m[seqplane][yi][xi][frame_index]=1.0; /* scale*(float)(*iptr); */
        }
    } else if(local_data_type==VAX_R4 || local_data_type==IEEE_R4) {
      for(yi=0, mptr=mdata; yi<img->dimy; yi++)
        for(xi=0; xi<img->dimx; xi++, mptr+=4) {
          memcpy(&img->m[seqplane][yi][xi][frame_index], mptr, 4);
          img->m[seqplane][yi][xi][frame_index]*=scale;
        }
    }
  } /* next matrix */
  if(IMG_TEST>3) printf("end of matrices.\n");
  
  free(mdata);
  ecat63EmptyMatlist(&mlist); 
  fclose(fp);

  /* seqplane is <0 only if this frame did not exist at all; return -1 */
  if(IMG_TEST>4) printf("last_seqplane := %d.\n", seqplane);
  if(seqplane<0) {imgSetStatus(img, STATUS_NOMATRIX); return STATUS_NOMATRIX;}

  /* check that correct number of planes was read */
  if(seqplane+1 != img->dimz) {
    imgSetStatus(img, STATUS_MISSINGMATRIX);
    return STATUS_MISSINGMATRIX;
  }

  /* For saving IMG data, only 2-byte VAX data types are allowed,
     so change it now */
  if(img->_dataType==VAX_I4 || img->_dataType==VAX_R4) img->_dataType=VAX_I2;

  /* All went well */
  imgSetStatus(img, STATUS_OK);
  return STATUS_OK;
}
/*****************************************************************************/

/*****************************************************************************/
/*!
 * Write one PET frame from IMG data struct into ECAT 6.3 image or sinogram
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
int imgWriteEcat63Frame(
  const char *fname, int frame_to_write, IMG *img, int frame_index
) {
  IMG test_img;
  int ret=0, pxlNr, zi, xi, yi, matrixId;
  ECAT63_mainheader   main_header;
  ECAT63_imageheader  image_header;
  ECAT63_scanheader   scan_header;
  void *sub_header=NULL;
  FILE *fp;
  float *fdata=NULL, *fptr;


  if(IMG_TEST>0) printf("\nimgWriteEcat63Frame(%s, %d, *img, %d)\n",
    fname, frame_to_write, frame_index);
  if(IMG_TEST>1) ECAT63_TEST=IMG_TEST-1;
  if(IMG_TEST>4) {
    char buf[32];
    if(!ctime_r_int(&img->scanStart, buf)) strcpy(buf, "1900-01-01 00:00:00");
    fprintf(stdout, "  scan_start_time := %s\n", buf);
  }

  /*
   *  Check the input 
   */
  if(fname==NULL) return STATUS_FAULT;
  if(img==NULL) return STATUS_FAULT;
  if(img->status!=IMG_STATUS_OCCUPIED) return STATUS_FAULT;
  if(frame_to_write<0) return STATUS_FAULT;
  if(frame_index<0 || frame_index>=img->dimt) return STATUS_FAULT;
  if(img->_fileFormat!=IMG_E63) return STATUS_FAULT;

  /* Initiate headers */
  memset(&main_header, 0, sizeof(ECAT63_mainheader));
  memset(&image_header,0, sizeof(ECAT63_imageheader));
  memset(&scan_header, 0, sizeof(ECAT63_scanheader));
  imgInit(&test_img);


  /*
   *  If file does not exist, then create it with new mainheader,
   *  and if it does exist, then read and check header information.
   *  Create or edit mainheader to contain correct frame nr.   
   *  In any case, leave file pointer open for write.   
   */
  if(access(fname, 0) == -1) { /* file does not exist*/

    if(IMG_TEST>1) printf("  new file\n");
    /* Set main header */
    imgSetEcat63MHeader(img, &main_header);
    if(IMG_TEST>6) {
      ecat63PrintMainheader(&main_header, stdout);
    }
    if(frame_to_write==0) frame_to_write=1;
    main_header.num_frames=frame_to_write;

    /* Open file, write main header and initiate matrix list */
    fp=ecat63Create(fname, &main_header); if(fp==NULL) return STATUS_NOWRITEPERM;

  } else { /* file exists*/
  
    if(IMG_TEST>1) printf("  existing file\n");
    /* Read header information for checking */
    ret=imgReadEcat63Header(fname, &test_img); if(ret!=0) return ret;
    if(IMG_TEST>1) {
      char buf[32];
      if(!ctime_r_int(&test_img.scanStart, buf)) 
        strcpy(buf, "1900-01-01 00:00:00");
      fprintf(stdout, "scan_start_time := %s\n", buf);
    }
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
    if((ret=ecat63ReadMainheader(fp, &main_header))!=0)
      return STATUS_NOMAINHEADER;
    if(frame_to_write==0) frame_to_write=main_header.num_frames+1;
    if(main_header.num_frames<frame_to_write)
      main_header.num_frames=frame_to_write;
    if((ret=ecat63WriteMainheader(fp, &main_header))!=0)
      return STATUS_NOWRITEPERM;
    if(IMG_TEST>0) {
      char tmp[32];
      printf("  scan_start_time := %s\n", 
             ecat63ScanstarttimeInt(&main_header, tmp));
    }
    
  }
  if(IMG_TEST>2) {
    printf("frame_to_write := %d\n", frame_to_write);
  }

  /* Allocate memory for matrix float data */
  pxlNr=img->dimx*img->dimy*img->dimz;
  fdata=(float*)malloc(pxlNr*sizeof(float));
  if(fdata==NULL) {fclose(fp); return STATUS_NOMEMORY;}
  
  /* Set matrix subheader */
  if(img->type==IMG_TYPE_RAW) sub_header=(void*)&scan_header;
  else if(img->type==IMG_TYPE_IMAGE) sub_header=&image_header;
  else {fclose(fp); return STATUS_FAULT;}
  imgSetEcat63SHeader(img, sub_header);

  /* Copy matrix pixel values to fdata */
  fptr=fdata;
  for(zi=0; zi<img->dimz; zi++)
    for(yi=0; yi<img->dimy; yi++) for(xi=0; xi<img->dimx; xi++)
      *fptr++=img->m[zi][yi][xi][frame_index];

  /* Write subheader and data, and set the rest of subheader contents */
  fptr=fdata; ret=-1;
  for(zi=0; zi<img->dimz; zi++, fptr+=img->dimx*img->dimy) {
    /* Create new matrix id (i.e. matnum) */
    matrixId=mat_numcod(frame_to_write, img->planeNumber[zi], 1, 0, 0);
    if(img->type==IMG_TYPE_RAW) {
      scan_header.frame_start_time=
        (int)temp_roundf(1000.*img->start[frame_index]);
      scan_header.frame_duration=
        (int)temp_roundf(1000.*(img->end[frame_index]-img->start[frame_index]));
      scan_header.prompts=temp_roundf(img->prompts[frame_index]);
      scan_header.delayed=temp_roundf(img->randoms[frame_index]);
      ret=ecat63WriteScanMatrix(fp, matrixId, &scan_header, fptr);
    } else {
      image_header.frame_start_time=
        (int)temp_roundf(1000.*img->start[frame_index]);
      image_header.frame_duration=
        (int)temp_roundf(1000.*(img->end[frame_index]-img->start[frame_index]));
      if(img->decayCorrection==IMG_DC_CORRECTED)
        image_header.decay_corr_fctr=img->decayCorrFactor[frame_index];
      else
        image_header.decay_corr_fctr=0.0;
      ret=ecat63WriteImageMatrix(fp, matrixId, &image_header, fptr);
    }
    if(ret!=0) break;
  } /* next plane*/
  free(fdata); fclose(fp);
  if(ret) return STATUS_DISKFULL;

  return STATUS_OK;
}
/*****************************************************************************/

/*****************************************************************************/
/** Copies Ecat6.3 image or scan sub header information from IMG struct.

 @param img source image structure
 @param h target sub header structure
 */
void imgSetEcat63SHeader(IMG *img, void *h) {
  ECAT63_imageheader *image_header;
  ECAT63_scanheader *scan_header;

  if(img->type==IMG_TYPE_RAW) {
    scan_header=(ECAT63_scanheader*)h;
    scan_header->data_type=VAX_I2;
    scan_header->dimension_1=img->dimx;
    scan_header->dimension_2=img->dimy;
    scan_header->frame_duration_sec=1;
    scan_header->scale_factor=1.0;
    scan_header->frame_start_time=0;
    scan_header->frame_duration=1000;
    scan_header->loss_correction_fctr=1.0;
    scan_header->sample_distance=img->sampleDistance/10.0;
  } else {
    image_header=(ECAT63_imageheader*)h;
    image_header->data_type=VAX_I2;
    image_header->num_dimensions=2;
    image_header->dimension_1=img->dimx;
    image_header->dimension_2=img->dimy;
    image_header->recon_scale=img->zoom;
    image_header->quant_scale=1.0;
    image_header->slice_width=img->sizez/10.;
    image_header->pixel_size=img->sizex/10.;
    image_header->frame_start_time=0;
    image_header->frame_duration=1000;
    image_header->plane_eff_corr_fctr=1.0; // now we don't know which plane
    image_header->decay_corr_fctr=1.0;
    image_header->loss_corr_fctr=1.0;
    image_header->ecat_calibration_fctr=1.0;
    image_header->well_counter_cal_fctr=1.0;
    image_header->quant_units=imgUnitToEcat6(img);
  }
}
/*****************************************************************************/

/*****************************************************************************/
