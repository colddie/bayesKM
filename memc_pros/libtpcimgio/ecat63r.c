/// @file ecat63r.c
/// @author Vesa Oikonen
/// @brief Procedures for reading ECAT 6.3 format.
///
///  Assumptions:
///  1. Data_type in headers specifies, whether ints, long ints and floats in
///     header and matrix data are in VAX format (1, 2, 3 and 4) or in IEEE
///     (i386, SUN) format.
///  2. Data is automatically converted to little or big endian when read,
///     according to the current platform.
///  3. Data is automatically converted out from the VAX format when read,
///     but header data_type remains as original.
///
/*****************************************************************************/
#include "libtpcimgio.h"
/*****************************************************************************/

/*****************************************************************************/
/** Read ECAT 6.3 main header.

 @return 0 if ok, 1 invalid input, 2 failed to find subheader block, 
  3 invalid magic number (should be "ECAT63") at start of file, 
  5 invalid data type, 6 invalid calibration factor, 7 invalid file type.
 */
int ecat63ReadMainheader(
  /** file pointer */
  FILE *fp, 
  /** target Ecat 6.3 main header structure */
  ECAT63_mainheader *h
) {
  unsigned char buf[MatBLKSIZE];
  int i;
  int little; /* 1 if current platform is little endian (i386), else 0 */
  int vaxdata=1; /* 1 if data is in VAX format, else 0 */

  if(ECAT63_TEST>0) {printf("ecat63ReadMainheader()\n"); fflush(stdout);}
  if(fp==NULL || h==NULL) return(1);
  little=little_endian(); if(ECAT63_TEST>1) {printf("little := %d\n", little); fflush(stdout);}
  /* Seek the first block */
  fseek(fp, 0, SEEK_SET); if(ftell(fp)!=0) return(1);
  /* Read the main header block */
  if(fread(buf, MatBLKSIZE, 1, fp)<1) return(2);

  /* Copy char data to header structure */
  memcpy(h->ecat_format, buf+0, 14); 
  memcpy(h->fill1, buf+14, 14);
  memcpy(h->original_file_name, buf+28, 20);
  if(ECAT63_TEST>10) {
    printf("original_file_name := '%.20s'\n", h->original_file_name); fflush(stdout);}
  memcpy(h->node_id, buf+56, 10);
  memcpy(h->isotope_code, buf+78, 8);
  memcpy(h->radiopharmaceutical, buf+90, 32);
  memcpy(h->study_name, buf+162, 12);
  memcpy(h->patient_id, buf+174, 16);
  memcpy(h->patient_name, buf+190, 32);
  h->patient_sex=buf[222];
  memcpy(h->patient_age, buf+223, 10);
  memcpy(h->patient_height, buf+233, 10);
  memcpy(h->patient_weight, buf+243, 10);
  h->patient_dexterity=buf[253];
  memcpy(h->physician_name, buf+254, 32);
  memcpy(h->operator_name, buf+286, 32);
  memcpy(h->study_description, buf+318, 32);
  memcpy(h->facility_name, buf+356, 20);
  memcpy(h->user_process_code, buf+462, 10);

  /* Copy short ints; in big endian platform, change byte order */
  if(!little) swabip(buf+50, 2);
  memcpy(&h->data_type, buf+50, 2); 
  if(ECAT63_TEST>10) {printf("main_header.data_type=%d\n", h->data_type); fflush(stdout);}
  if(h->data_type<1) {
    if(ECAT63_TEST>1) {printf("invalid data_type; assuming VAX_I2\n"); fflush(stdout);}
    h->data_type=VAX_I2;
  }
  if(h->data_type>4) vaxdata=0;

  if(!little) swabip(buf+48, 2);
  memcpy(&h->sw_version, buf+48, 2);
  if(!little) swabip(buf+52, 2);
  memcpy(&h->system_type, buf+52, 2);
  if(!little) swabip(buf+54, 2);
  memcpy(&h->file_type, buf+54, 2);
  if(ECAT63_TEST>10) {printf("main_header.file_type=%d\n", h->file_type); fflush(stdout);}
  if(!little) swabip(buf+66, 2);
  memcpy(&h->scan_start_day, buf+66, 2);
  if(!little) swabip(buf+68, 2);
  memcpy(&h->scan_start_month, buf+68, 2);
  if(!little) swabip(buf+70, 2);
  memcpy(&h->scan_start_year, buf+70, 2);
  if(!little) swabip(buf+72, 2);
  memcpy(&h->scan_start_hour, buf+72, 2);
  if(!little) swabip(buf+74, 2);
  memcpy(&h->scan_start_minute, buf+74, 2);
  if(!little) swabip(buf+76, 2);
  memcpy(&h->scan_start_second, buf+76, 2);
  if(!little) swabip(buf+134, 2);
  memcpy(&h->rot_source_speed, buf+134, 2);
  if(!little) swabip(buf+136, 2);
  memcpy(&h->wobble_speed, buf+136, 2);
  if(!little) swabip(buf+138, 2);
  memcpy(&h->transm_source_type, buf+138, 2);
  if(!little) swabip(buf+148, 2);
  memcpy(&h->transaxial_samp_mode, buf+148, 2);
  if(!little) swabip(buf+150, 2);
  memcpy(&h->coin_samp_mode, buf+150, 2);
  if(!little) swabip(buf+152, 2);
  memcpy(&h->axial_samp_mode, buf+152, 2);
  if(!little) swabip(buf+158, 2);
  memcpy(&h->calibration_units, buf+158, 2);
  if(!little) swabip(buf+160, 2);
  memcpy(&h->compression_code, buf+160, 2);
  if(!little) swabip(buf+350, 2);
  memcpy(&h->acquisition_type, buf+350, 2);
  if(!little) swabip(buf+352, 2);
  memcpy(&h->bed_type, buf+352, 2);
  if(!little) swabip(buf+354, 2);
  memcpy(&h->septa_type, buf+354, 2);
  if(!little) swabip(buf+376, 2);
  memcpy(&h->num_planes, buf+376, 2);
  if(!little) swabip(buf+378, 2);
  memcpy(&h->num_frames, buf+378, 2);
  if(!little) swabip(buf+380, 2);
  memcpy(&h->num_gates, buf+380, 2);
  if(!little) swabip(buf+382, 2);
  memcpy(&h->num_bed_pos, buf+382, 2);
  if(!little) swabip(buf+452, 2);
  memcpy(&h->lwr_sctr_thres, buf+452, 2);
  if(!little) swabip(buf+454, 2);
  memcpy(&h->lwr_true_thres, buf+454, 2);
  if(!little) swabip(buf+456, 2);
  memcpy(&h->upr_true_thres, buf+456, 2);
  if(!little) swabip(buf+472, 40);
  memcpy(h->fill2, buf+472, 40);

  /* Copy floats */
  h->isotope_halflife=ecat63rFloat(buf+86, vaxdata, little);
  h->gantry_tilt=ecat63rFloat(buf+122, vaxdata, little);
  h->gantry_rotation=ecat63rFloat(buf+126, vaxdata, little);
  h->bed_elevation=ecat63rFloat(buf+130, vaxdata, little);
  h->axial_fov=ecat63rFloat(buf+140, vaxdata, little);
  h->transaxial_fov=ecat63rFloat(buf+144, vaxdata, little);
  h->calibration_factor=ecat63rFloat(buf+154, vaxdata, little);
  h->init_bed_position=ecat63rFloat(buf+384, vaxdata, little);
  for(i=0; i<15; i++) h->bed_offset[i]=ecat63rFloat(buf+388+i*4, vaxdata, little);
  h->plane_separation=ecat63rFloat(buf+448, vaxdata, little);
  h->collimator=ecat63rFloat(buf+458, vaxdata, little);

  /* Check file format and platform */
  if(ECAT63_TEST>1) {printf("ecat_format='%.14s'\n", h->ecat_format); fflush(stdout);}
  /* if format is not specified, ECAT63 is assumed */
  if(h->ecat_format[0]==(char)0) {
    strcpy(h->ecat_format, "ECAT63");
  }
  /* only ECAT63 format is approved here */
  if(ECAT63_TEST>1) {printf("ecat_format='%.14s'\n", h->ecat_format); fflush(stdout);}
  if(strncmp(h->ecat_format, "ECAT63", 6)!=0) return(3);

  /* Check that most important contents are ok */
  if(ECAT63_TEST>3) {
    printf("  mhdr.data_type := %s\n", ecat63Datatype(h->data_type));
    fflush(stdout);
  }
  if(h->data_type<BYTE_TYPE || h->data_type>SUN_I4) {
    if(ECAT63_TEST>1) {printf("Invalid data types; probable conversion error.\n"); fflush(stdout);}
    return(5);
  }
  if(h->calibration_factor<0.0 || h->calibration_factor>1.0e12) {
    if(ECAT63_TEST>1) {
      printf("Invalid calibration factor; possible conversion error.\n"); fflush(stdout);}
    return(6);
  }
  if(h->file_type!=RAW_DATA && h->file_type!=IMAGE_DATA &&
     h->file_type!=ATTN_DATA && h->file_type!=NORM_DATA) {
    if(ECAT63_TEST>1) {printf("Invalid file types; probable conversion error.\n"); fflush(stdout);}
    return(7);
  }

  return(0);
}
/*****************************************************************************/

/*****************************************************************************/
/** Read ECAT 6.3 image header.
   @return 0 if ok, 1 invalid input, 2 failed to find block, 
   3 failed to read block, 4 invalid data type, 5 invalid calibration factor, 
   6 invalid frame duration.
 */
int ecat63ReadImageheader(
  /** File pointer to ECAT file; position is not important. */
  FILE *fp, 
  /** Block number from matrix list [2..number of blocks]. */
  int blk, 
  /** Pointer to the header struct where contents are put. */
  ECAT63_imageheader *h,
  /** Verbose level; if zero, then only warnings are printed into stderr. */
  int verbose,
  /** Pointer to error message, at least 128 characters; enter NULL, if not needed. */
  char *errmsg 
) {
  unsigned char buf[MatBLKSIZE];
  int i;
  int little; /* 1 if current platform is little endian (i386), else 0 */
  int vaxdata=1; /* 1 if data is in VAX format, else 0 */

  if(verbose>0) {printf("ecat63ReadImageheader(fp, %d ih)\n", blk); fflush(stdout);}
  if(fp==NULL || blk<2 || h==NULL) {
    if(errmsg!=NULL) strcpy(errmsg, "invalid input");
    if(verbose>0) {fprintf(stderr, "Invalid input.\n"); fflush(stderr);}
    return(1);
  }
  little=little_endian();
  /* Seek the subheader block */
  fseek(fp, (blk-1)*MatBLKSIZE, SEEK_SET); 
  if(ftell(fp)!=(blk-1)*MatBLKSIZE) {
    if(errmsg!=NULL) sprintf(errmsg, "failed to find block %d", blk);
    if(verbose>0) {fprintf(stderr, "Failed to find block %d.\n", blk); fflush(stderr);}
    return(2);
  }
  /* Read the subheader block */
  if(fread(buf, MatBLKSIZE, 1, fp)<1) {
    if(errmsg!=NULL) sprintf(errmsg, "failed to read block %d", blk);
    if(verbose>0) {fprintf(stderr, "Failed to read block %d.\n", blk); fflush(stderr);}
    return(3);
  }

  /* Copy char data to header structure */
  memcpy(h->fill1, buf+0, 126); memcpy(h->annotation, buf+420, 40);

  /* Copy short ints */
  /* in big endian platform, change byte order temporarily */
  if(!little) swabip(buf, MatBLKSIZE);
  memcpy(&h->data_type, buf+126, 2); if(h->data_type>4) vaxdata=0;
  if(verbose>10) printf("data_type=%d\n", h->data_type);
  memcpy(&h->num_dimensions, buf+128, 2);
  memcpy(&h->dimension_1, buf+132, 2); memcpy(&h->dimension_2, buf+134, 2);
  memcpy(&h->image_min, buf+176, 2); memcpy(&h->image_max, buf+178, 2);
  memcpy(&h->slice_location, buf+200, 2); memcpy(&h->recon_start_hour, buf+202, 2);
  memcpy(&h->recon_start_min, buf+204, 2); memcpy(&h->recon_start_sec, buf+206, 2);
  memcpy(&h->filter_code, buf+236, 2); memcpy(&h->processing_code, buf+376, 2);
  memcpy(&h->quant_units, buf+380, 2); memcpy(&h->recon_start_day, buf+382, 2);
  memcpy(&h->recon_start_month, buf+384, 2); memcpy(&h->recon_start_year, buf+386, 2);
  memcpy(h->fill2, buf+460, 52);
  /* Change back the byte order */
  if(!little) swabip(buf, MatBLKSIZE);

  /* Copy ints */
  h->frame_duration=ecat63rInt(buf+192, vaxdata, little);
  h->frame_start_time=ecat63rInt(buf+196, vaxdata, little);
  h->recon_duration=ecat63rInt(buf+208, vaxdata, little);
  h->scan_matrix_num=ecat63rInt(buf+238, vaxdata, little);
  h->norm_matrix_num=ecat63rInt(buf+242, vaxdata, little);
  h->atten_cor_mat_num=ecat63rInt(buf+246, vaxdata, little);

  /* Copy floats */
  h->x_origin=ecat63rFloat(buf+160, vaxdata, little);
  h->y_origin=ecat63rFloat(buf+164, vaxdata, little);
  h->recon_scale=ecat63rFloat(buf+168, vaxdata, little);
  h->quant_scale=ecat63rFloat(buf+172, vaxdata, little);
  h->pixel_size=ecat63rFloat(buf+184, vaxdata, little);
  h->slice_width=ecat63rFloat(buf+188, vaxdata, little);
  h->image_rotation=ecat63rFloat(buf+296, vaxdata, little);
  h->plane_eff_corr_fctr=ecat63rFloat(buf+300, vaxdata, little);
  h->decay_corr_fctr=ecat63rFloat(buf+304, vaxdata, little);
  h->loss_corr_fctr=ecat63rFloat(buf+308, vaxdata, little);
  h->intrinsic_tilt=ecat63rFloat(buf+312, vaxdata, little);
  h->ecat_calibration_fctr=ecat63rFloat(buf+388, vaxdata, little);
  h->well_counter_cal_fctr=ecat63rFloat(buf+392, vaxdata, little);
  for(i=0; i<6; i++) h->filter_params[i]=ecat63rFloat(buf+396+i*4, vaxdata, little);

  /* Check that most important contents are ok */
  if(h->data_type<BYTE_TYPE || h->data_type>SUN_I4) {
    if(errmsg!=NULL) strcpy(errmsg, "invalid data types; probable conversion error");
    if(verbose>0) {fprintf(stderr, "Invalid data types; probable conversion error.\n"); fflush(stderr);}
    if(verbose>1) {printf("data_type := %d\n", h->data_type); fflush(stdout);}
    return(4);
  }
  if(h->ecat_calibration_fctr<0.0 || h->ecat_calibration_fctr>1.0e10) {
    if(errmsg!=NULL) strcpy(errmsg, "invalid calibration factor; probable conversion error");
    if(verbose>0) {fprintf(stderr, "Invalid calibration factor; probable conversion error.\n"); fflush(stderr);}
    return(5);
  }
  if(h->frame_duration<0.0 || h->frame_duration>1.0e12) {
    if(errmsg!=NULL) strcpy(errmsg, "invalid frame duration; probable conversion error");
    if(verbose>0) {fprintf(stderr, "Invalid frame duration; probable conversion error.\n"); fflush(stderr);}
    return(6);
  }
  if(errmsg!=NULL) strcpy(errmsg, "ok");

  return(0);
}
/*****************************************************************************/

/*****************************************************************************/
/** Read ECAT 6.3 attenuation header
   @return 0 if ok, 1 failed to find block, 3 failed to read block, 
    4 invalid data type, 5 invalid scale factor.
 */
int ecat63ReadAttnheader(
  /** File pointer to ECAT file; position is not important. */
  FILE *fp, 
  /** Block number from matrix list [2..number of blocks]. */
  int blk, 
  /** Pointer to the header struct where contents are put. */
  ECAT63_attnheader *h,
  /** Verbose level; if zero, then only warnings are printed into stderr. */
  int verbose,
  /** Pointer to error message, at least 128 characters; enter NULL, if not needed. */
  char *errmsg 
) {
  unsigned char buf[MatBLKSIZE];
  int little; /* 1 if current platform is little endian (i386), else 0 */
  int vaxdata=1; /* 1 if data is in VAX format, else 0 */

  if(ECAT63_TEST) printf("ecat63ReadAttnheader(fp, %d, ah)\n", blk);
  if(fp==NULL || blk<2 || h==NULL) {
    if(errmsg!=NULL) strcpy(errmsg, "invalid input");
    if(verbose>0) fprintf(stderr, "Invalid input.\n");
    return(1);
  }
  little=little_endian();
  /* Seek the subheader block */
  fseek(fp, (blk-1)*MatBLKSIZE, SEEK_SET);
  if(ftell(fp)!=(blk-1)*MatBLKSIZE) {
    if(errmsg!=NULL) sprintf(errmsg, "failed to find block %d", blk);
    if(verbose>0) fprintf(stderr, "Failed to find block %d.\n", blk);
    return(2);
  }
  /* Read the subheader block */
  if(fread(buf, MatBLKSIZE, 1, fp)<1) {
    if(errmsg!=NULL) sprintf(errmsg, "failed to read block %d", blk);
    if(verbose>0) fprintf(stderr, "Failed to read block %d.\n", blk);
    return(3);
  }

  /* Copy short ints */
  /* in big endian platform, change byte order temporarily */
  if(!little) swabip(buf, MatBLKSIZE);
  memcpy(&h->data_type, buf+126, 2); if(h->data_type>4) vaxdata=0;
  /*printf("data_type=%d\n", h->data_type);*/
  memcpy(&h->attenuation_type, buf+128, 2);
  memcpy(&h->dimension_1, buf+132, 2); memcpy(&h->dimension_2, buf+134, 2);
  /* Change back the byte order */
  if(!little) swabip(buf, MatBLKSIZE);

  /* Copy floats */
  h->scale_factor=ecat63rFloat(buf+182, vaxdata, little);
  h->x_origin=ecat63rFloat(buf+186, vaxdata, little);
  h->y_origin=ecat63rFloat(buf+190, vaxdata, little);
  h->x_radius=ecat63rFloat(buf+194, vaxdata, little);
  h->y_radius=ecat63rFloat(buf+198, vaxdata, little);
  h->tilt_angle=ecat63rFloat(buf+202, vaxdata, little);
  h->attenuation_coeff=ecat63rFloat(buf+206, vaxdata, little);
  h->sample_distance=ecat63rFloat(buf+210, vaxdata, little);

  /* Check that most important contents are ok */
  if(h->data_type<BYTE_TYPE || h->data_type>SUN_I4) {
    if(errmsg!=NULL) strcpy(errmsg, "invalid data types; probable conversion error");
    if(verbose>0) fprintf(stderr, "Invalid data types; probable conversion error.\n");
    if(verbose>1) printf("data_type := %d\n", h->data_type);
    return(4);
  }
  if(h->scale_factor<=0.0 || h->scale_factor>1.0e8) {
    if(errmsg!=NULL) strcpy(errmsg, "invalid scale factor; probable conversion error");
    if(verbose>0) fprintf(stderr, "Invalid scale factor; probable conversion error.\n");
    return(5);
  }
  if(errmsg!=NULL) strcpy(errmsg, "ok");

  return(0);
}
/*****************************************************************************/

/*****************************************************************************/
/** Read ECAT 6.3 scan header.
   @return 0 if ok, 1 invalid input, 2 failed to find block, 
   3 failed to read block, 4 invalid data type, 5 invalid scale factor, 
   6 invalid frame duration.
 */
int ecat63ReadScanheader(
  /** File pointer to ECAT file; position is not important. */
  FILE *fp, 
  /** Block number from matrix list [2..number of blocks]. */
  int blk, 
  /** Pointer to the header struct where contents are put. */
  ECAT63_scanheader *h,
  /** Verbose level; if zero, then only warnings are printed into stderr. */
  int verbose,
  /** Pointer to error message, at least 128 characters; enter NULL, if not needed. */
  char *errmsg 
) {
  unsigned char buf[MatBLKSIZE];
  int i;
  int little; /* 1 if current platform is little endian (i386), else 0 */
  int vaxdata=1; /* 1 if data is in VAX format, else 0 */

  if(ECAT63_TEST) printf("ecat63ReadScanheader(fp, %d, sh)\n", blk);
  if(fp==NULL || blk<2 || h==NULL) {
    if(errmsg!=NULL) strcpy(errmsg, "invalid input");
    if(verbose>0) fprintf(stderr, "Invalid input.\n");
    return(1);
  }
  little=little_endian();
  /* Seek the subheader block */
  fseek(fp, (blk-1)*MatBLKSIZE, SEEK_SET); 
  if(ftell(fp)!=(blk-1)*MatBLKSIZE) {
    if(errmsg!=NULL) sprintf(errmsg, "failed to find block %d", blk);
    if(verbose>0) fprintf(stderr, "Failed to find block %d.\n", blk);
    return(2);
  }
  /* Read the subheader block */
  if(fread(buf, MatBLKSIZE, 1, fp)<1) {
    if(errmsg!=NULL) sprintf(errmsg, "failed to read block %d", blk);
    if(verbose>0) fprintf(stderr, "Failed to read block %d.\n", blk);
    return(3);
  }

  /* Copy char data to header structure */
  memcpy(h->fill1, buf+0, 126);

  /* Copy short ints */
  /* in big endian platform, change byte order temporarily */
  if(!little) swabip(buf, MatBLKSIZE);
  memcpy(&h->data_type, buf+126, 2); if(h->data_type>4) vaxdata=0;
  /*printf("data_type=%d\n", h->data_type);*/
  memcpy(&h->dimension_1, buf+132, 2); memcpy(&h->dimension_2, buf+134, 2);
  memcpy(&h->smoothing, buf+136, 2); memcpy(&h->processing_code, buf+138, 2);
  memcpy(&h->frame_duration_sec, buf+170, 2);
  memcpy(&h->scan_min, buf+192, 2); memcpy(&h->scan_max, buf+194, 2);
  memcpy(h->fill2, buf+468, 44);
  /* Change back the byte order */
  if(!little) swabip(buf, MatBLKSIZE);

  /* Copy ints */
  h->gate_duration=ecat63rInt(buf+172, vaxdata, little);
  h->r_wave_offset=ecat63rInt(buf+176, vaxdata, little);
  h->prompts=ecat63rInt(buf+196, vaxdata, little);
  h->delayed=ecat63rInt(buf+200, vaxdata, little);
  h->multiples=ecat63rInt(buf+204, vaxdata, little);
  h->net_trues=ecat63rInt(buf+208, vaxdata, little);
  h->total_coin_rate=ecat63rInt(buf+452, vaxdata, little);
  h->frame_start_time=ecat63rInt(buf+456, vaxdata, little);
  h->frame_duration=ecat63rInt(buf+460, vaxdata, little);

  /* Copy floats */
  h->sample_distance=ecat63rFloat(buf+146, vaxdata, little);
  h->isotope_halflife=ecat63rFloat(buf+166, vaxdata, little);
  h->scale_factor=ecat63rFloat(buf+182, vaxdata, little);
  for(i=0; i<16; i++) h->cor_singles[i]=ecat63rFloat(buf+316+i*4, vaxdata, little);
  for(i=0; i<16; i++) h->uncor_singles[i]=ecat63rFloat(buf+380+i*4, vaxdata, little);
  h->tot_avg_cor=ecat63rFloat(buf+444, vaxdata, little);
  h->tot_avg_uncor=ecat63rFloat(buf+448, vaxdata, little);
  h->loss_correction_fctr=ecat63rFloat(buf+464, vaxdata, little);

  /* Check that most important contents are ok */
  if(h->data_type<BYTE_TYPE || h->data_type>SUN_I4) {
    if(errmsg!=NULL) strcpy(errmsg, "invalid data types; probable conversion error");
    if(verbose>0) fprintf(stderr, "Invalid data types; probable conversion error.\n");
    if(verbose>1) printf("data_type := %d\n", h->data_type);
    return(4);
  }
  if(h->scale_factor<=0.0 || h->scale_factor>1.0e8) {
    if(errmsg!=NULL) strcpy(errmsg, "invalid calibration factor; probable conversion error");
    if(verbose>0) fprintf(stderr, "Invalid calibration factor; probable conversion error.\n");
    return(5);
  }
  if(h->frame_duration<0.0 || h->frame_duration>1.0e12) {
    if(errmsg!=NULL) strcpy(errmsg, "invalid frame duration; probable conversion error");
    if(verbose>0) fprintf(stderr, "Invalid frame duration; probable conversion error.\n");
    return(6);
  }
  if(errmsg!=NULL) strcpy(errmsg, "ok");

  return(0);
}
/*****************************************************************************/

/*****************************************************************************/
/** Read ECAT 6.3 normalization header.

    Note that ECAT 6.3 normalization data is usually stored in scan file format,
    not in normalization format.

   @return 0 if ok, 1 invalid input, 2 failed to find block, 
    3 failed to read block, 4 invalid data type, 5 invalid scale factor.
 */
int ecat63ReadNormheader(
  /** File pointer to ECAT file; position is not important. */
  FILE *fp, 
  /** Block number from matrix list [2..number of blocks]. */
  int blk, 
  /** Pointer to the header struct where contents are put. */
  ECAT63_normheader *h,
  /** Verbose level; if zero, then only warnings are printed into stderr. */
  int verbose,
  /** Pointer to error message, at least 128 characters; enter NULL, if not needed. */
  char *errmsg 
) {
  unsigned char buf[MatBLKSIZE];
  int little; /* 1 if current platform is little endian (i386), else 0 */
  int vaxdata=1; /* 1 if data is in VAX format, else 0 */

  if(ECAT63_TEST) printf("ecat63ReadNormheader(fp, %d, nh)\n", blk);
  if(fp==NULL || blk<2 || h==NULL) {
    if(errmsg!=NULL) strcpy(errmsg, "invalid input");
    if(verbose>0) fprintf(stderr, "Invalid input.\n");
    return(1);
  }
  little=little_endian();
  /* Seek the subheader block */
  fseek(fp, (blk-1)*MatBLKSIZE, SEEK_SET); 
  if(ftell(fp)!=(blk-1)*MatBLKSIZE) {
    if(errmsg!=NULL) sprintf(errmsg, "failed to find block %d", blk);
    if(verbose>0) fprintf(stderr, "Failed to find block %d.\n", blk);
    return(2);
  }
  /* Read the subheader block */
  if(fread(buf, MatBLKSIZE, 1, fp)<1) {
    if(errmsg!=NULL) sprintf(errmsg, "failed to read block %d", blk);
    if(verbose>0) fprintf(stderr, "Failed to read block %d.\n", blk);
    return(3);
  }

  /* Copy short ints */
  /* in big endian platform, change byte order temporarily */
  if(!little) swabip(buf, MatBLKSIZE);
  memcpy(&h->data_type, buf+126, 2); if(h->data_type>4) vaxdata=0;
  if(verbose>10) printf("data_type=%d\n", h->data_type);
  memcpy(&h->dimension_1, buf+132, 2); memcpy(&h->dimension_2, buf+134, 2);
  memcpy(&h->norm_hour, buf+186, 2); memcpy(&h->norm_minute, buf+188, 2);
  memcpy(&h->norm_second, buf+190, 2); memcpy(&h->norm_day, buf+192, 2);
  memcpy(&h->norm_month, buf+194, 2); memcpy(&h->norm_year, buf+196, 2);
  /* Change back the byte order */
  if(!little) swabip(buf, MatBLKSIZE);

  /* Copy floats */
  h->scale_factor=ecat63rFloat(buf+182, vaxdata, little);
  h->fov_source_width=ecat63rFloat(buf+198, vaxdata, little);

  /* Check that most important contents are ok */
  if(h->data_type<BYTE_TYPE || h->data_type>SUN_I4) {
    if(errmsg!=NULL) strcpy(errmsg, "invalid data types; probable conversion error");
    if(verbose>0) fprintf(stderr, "Invalid data types; probable conversion error.\n");
    if(verbose>1) printf("data_type := %d\n", h->data_type);
    return(4);
  }
  if(h->scale_factor<=0.0 || h->scale_factor>1.0e8) {
    if(errmsg!=NULL) strcpy(errmsg, "invalid scale factor; probable conversion error");
    if(verbose>0) fprintf(stderr, "Invalid scale factor; probable conversion error.\n");
    return(5);
  }
  if(errmsg!=NULL) strcpy(errmsg, "ok");

  return(0);
}
/*****************************************************************************/

/*****************************************************************************/
/*!
 * Read ECAT 6.3 matrix data and convert byte order if necessary
 *  Remember to allocate memory for full blocks!
 *  There are differences here when compared to ecat7.c
 *
 * @param fp file pointer from where data is read
 * @param strtblk starting block [>= 1]
 * @param blkNr number of block to be read [>= 0]
 * @param data pointer to block where data is read
 * @param dtype data type code
 * @return 0 if ok, 1 invalid input, 2 failed to read data,
 * 9 failed to find starting block from file,
 */
int ecat63ReadMatdata(FILE *fp, int strtblk, int blkNr, char *data, int dtype) {
  int i, n, little, err=0;
  char *cptr;
  float f;


  if(ECAT63_TEST) printf("ecat63ReadMatdata(fp, %d, %d, data, %d)\n", strtblk, blkNr, dtype);
  /* Check the arguments */
  if(blkNr<=0 || strtblk<1 || data==NULL) return(1);
  /* Seek the first data block */
  fseek(fp, (strtblk-1)*MatBLKSIZE, SEEK_SET);
  if(ftell(fp)!=(strtblk-1)*MatBLKSIZE) return(9);
  /* Read the data blocks */
  if(fread(data, MatBLKSIZE, blkNr, fp) < (unsigned int)blkNr) return(2);
  /* Translate data if necessary */
  little=little_endian();
  switch(dtype) {
    case BYTE_TYPE: /* byte format...no translation necessary */
      break;
    case VAX_I2:    /* byte conversion necessary on big endian platform */
      if(!little) {cptr=data; swabip(cptr, blkNr*MatBLKSIZE);}
      break;
    case VAX_I4:
      for(i=0, cptr=data; i<blkNr*MatBLKSIZE; i+=4, cptr+=4) {
        n=ecat63rInt(cptr, 1, little); memcpy(cptr, &n, 4);
      }
      break;
    case VAX_R4:
      for(i=0, cptr=data; i<blkNr*MatBLKSIZE; i+=4, cptr+=4) {
        f=ecat63rFloat(cptr, 1, little); memcpy(cptr, &f, 4);
      }
      break;
    case IEEE_R4:   /* IEEE float ; byte conversion necessary on big end platforms */
    case SUN_I4:    /* SUN int ; byte conversion necessary on big end platforms */
      if(!little) swawbip(data, blkNr*MatBLKSIZE);
      break;
    case SUN_I2:    /* SUN short ; byte conversion necessary on big end platforms */
      if(!little) swabip(data, blkNr*MatBLKSIZE);
      break;
    default:  /* if something else, for now think it as an error */
      err=2;
      break;
  }
  return(err);
}
/*****************************************************************************/

/*****************************************************************************/
/** Read ECAT63 image matrix header and data.

    If only header is to be read, set last_block=first_block.
    Note: data is not calibrated with factor in main header.

   @return 0 if ok, 1 invalid input, 5 failed to read sub header, 
    6 invalid (x,y,z) dimensions, 8 failed to allocate memory for meta-data,
    9 failed to read matrix data, 11 failed to allocate memory for image data
 */
int ecat63ReadImageMatrix(
  /** ECAT file pointer */
  FILE *fp, 
  /** Subheader record number */
  int first_block, 
  /** Last data block number */
  int last_block, 
  /** Ptr to subheader data which is filled */
  ECAT63_imageheader *h, 
  /** Ptr to the address of the matrix data; any old contents are not freed,
      therefore you must free the pointer after data is copied elsewhere. */
  float **fdata
) {
  int i, ret, blockNr, pxlNr;
  char *mdata, *mptr, errmsg[128];
  float *_fdata, *fptr;
  short int *sptr;
  int *iptr;
  
  
  if(ECAT63_TEST) printf("ecat63ReadImageMatrix(fp, %d, %d, hdr, fdata)\n",
    first_block, last_block);
  if(fp==NULL || first_block<=MatFirstDirBlk || h==NULL) {
    sprintf(ecat63errmsg, "invalid function parameter.\n");
    return(1);
  }
  *fdata=(float*)NULL;
  
  /* Read subheader */
  ret=ecat63ReadImageheader(fp, first_block, h, ECAT63_TEST-2, errmsg);
  if(ret) {
    strcpy(ecat63errmsg, errmsg);
    return(5);
  }
  if(ECAT63_TEST>4) ecat63PrintImageheader(h, stdout);
  pxlNr=h->dimension_1*h->dimension_2;
  if(pxlNr<=0) {
    sprintf(ecat63errmsg, "invalid matrix dimension.\n");
    return(6);  
  }
  
  /* Read matrix data */
  blockNr=last_block-first_block; if(blockNr<1) return(0);
  mdata=(char*)malloc(blockNr*MatBLKSIZE);
  if(mdata==NULL) {
    sprintf(ecat63errmsg, "cannot allocate memory.\n");
    return(8);  
  }
  mptr=mdata;
  ret=ecat63ReadMatdata(fp, first_block+1, blockNr, mptr, h->data_type);
  if(ret || mdata==NULL) {
    sprintf(ecat63errmsg, "cannot read matrix data (%d).\n", ret);
    free(mdata); return(9);
  }
  
  /* Allocate memory for float data */
  _fdata=(float*)malloc(pxlNr*sizeof(float));
  if(_fdata==NULL) {
    sprintf(ecat63errmsg, "cannot allocate memory.\n");
    free(mdata); return(11);  
  }

  /* Convert matrix data to floats */
  if(h->ecat_calibration_fctr>0.0) h->quant_scale*=h->ecat_calibration_fctr;
  fptr=_fdata; mptr=mdata;
  if(h->data_type==BYTE_TYPE) {
    for(i=0; i<pxlNr; i++, mptr++, fptr++)
      *fptr=h->quant_scale*(float)(*mptr);
  } else if(h->data_type==VAX_I2 || h->data_type==SUN_I2) {
    for(i=0; i<pxlNr; i++, mptr+=2, fptr++) {
      sptr=(short int*)mptr;
      *fptr=h->quant_scale*(float)(*sptr);
      if(!(*fptr>-1.0E+22 && *fptr<1.0E+22)) *fptr=0.0;
    }
  } else if(h->data_type==VAX_I4 || h->data_type==SUN_I4) {
    for(i=0; i<pxlNr; i++, mptr+=4, fptr++) {
      iptr=(int*)mptr;
      *fptr=h->quant_scale*(float)(*iptr);
      if(!(*fptr>-1.0E+22 && *fptr<1.0E+22)) *fptr=0.0;
    }
  } else if(h->data_type==VAX_R4 || h->data_type==IEEE_R4) {
    memcpy(fptr, mptr, pxlNr*4);
    for(i=0; i<pxlNr; i++, fptr++) {
      *fptr *= h->quant_scale;
      if(!(*fptr>-1.0E+22 && *fptr<1.0E+22)) *fptr=0.0;
    }
  }
  free(mdata);
  *fdata=_fdata;

  return(0);
}
/*****************************************************************************/

/*****************************************************************************/
/** Read ECAT63 scan matrix header and data.

    If only header is to be read, set last_block=first_block.
    Note: data is not calibrated with factor in main header.
 
   @return 0 if ok, 1 invalid input, 5 failed to read sub header,
    6 invalid (x,y,z) dimensions, 8 failed to allocate memory for meta-data,
    9 failed to read matrix data, 11 failed to allocate memory for the data
 */
int ecat63ReadScanMatrix(
  /** ECAT file pointer. */
  FILE *fp, 
  /** Subheader record number. */
  int first_block, 
  /** Last data block number. */
  int last_block, 
  /** Pointer to subheader data which is filled. */
  ECAT63_scanheader *h, 
  /** Double pointer to the address of the matrix data; any old contents are 
      not freed, therefore you must free the pointer after data is copied 
      elsewhere. */
  float **fdata
) {
  int i, ret, blockNr, pxlNr;
  char *mdata, *mptr, errmsg[128];
  float *_fdata, *fptr;
  short int *sptr;
  int *iptr;
  
  
  if(ECAT63_TEST) printf("ecat63ReadScanMatrix(fp, %d, %d, hdr, fdata)\n",
    first_block, last_block);
  if(fp==NULL || first_block<=MatFirstDirBlk || h==NULL) {
    sprintf(ecat63errmsg, "invalid function parameter.\n");
    return(1);
  }
  *fdata=(float*)NULL;
  
  /* Read subheader */
  ret=ecat63ReadScanheader(fp, first_block, h, ECAT63_TEST-2, errmsg);
  if(ret) {
    strcpy(ecat63errmsg, errmsg);
    return(5);
  }
  if(ECAT63_TEST>4) ecat63PrintScanheader(h, stdout);
  pxlNr=h->dimension_1*h->dimension_2;
  if(pxlNr<=0) {
    sprintf(ecat63errmsg, "invalid matrix dimension.\n");
    return(6);  
  }
  
  /* Read matrix data */
  blockNr=last_block-first_block; if(blockNr<1) return(0);
  mdata=(char*)malloc(blockNr*MatBLKSIZE);
  if(mdata==NULL) {
    sprintf(ecat63errmsg, "cannot allocate memory.\n");
    return(8);  
  }
  mptr=mdata;
  ret=ecat63ReadMatdata(fp, first_block+1, blockNr, mptr, h->data_type);
  if(ret || mdata==NULL) {
    sprintf(ecat63errmsg, "cannot read matrix data (%d).\n", ret);
    free(mdata); return(9);
  }
  
  /* Allocate memory for float data */
  _fdata=(float*)malloc(pxlNr*sizeof(float));
  if(_fdata==NULL) {
    sprintf(ecat63errmsg, "cannot allocate memory.\n");
    free(mdata); return(11);  
  }

  /* Convert matrix data to floats */
  fptr=_fdata; mptr=mdata;
  if(h->data_type==BYTE_TYPE) {
    for(i=0; i<pxlNr; i++, mptr++, fptr++)
      *fptr=h->scale_factor*(float)(*mptr);
  } else if(h->data_type==VAX_I2 || h->data_type==SUN_I2) {
    for(i=0; i<pxlNr; i++, mptr+=2, fptr++) {
      sptr=(short int*)mptr;
      *fptr=h->scale_factor*(float)(*sptr);
    }
  } else if(h->data_type==VAX_I4 || h->data_type==SUN_I4) {
    for(i=0; i<pxlNr; i++, mptr+=4, fptr++) {
      iptr=(int*)mptr;
      *fptr=h->scale_factor*(float)(*iptr);
    }
  } else if(h->data_type==VAX_R4 || h->data_type==IEEE_R4) {
    memcpy(fptr, mptr, pxlNr*4);
    for(i=0; i<pxlNr; i++, fptr++) *fptr *= h->scale_factor;
  }
  free(mdata);
  *fdata=_fdata;

  return(0);
}
/*****************************************************************************/

/*****************************************************************************/
/** Read ECAT63 attenuation matrix header and data.

    If only header is to be read, set last_block=first_block.
    Note: data is not calibrated with factor in main header.
 
   @return 0 if ok, 1 invalid input, 5 failed to read sub header,
    6 invalid (x,y,z) dimensions, 8 failed to allocate memory for meta-data,
    9 failed to read matrix data, 11 failed to allocate memory for the data
 */
int ecat63ReadAttnMatrix(
  /** ECAT file pointer. */
  FILE *fp, 
  /** Subheader record number. */
  int first_block, 
  /** Last data block number. */
  int last_block, 
  /** Pointer to subheader data which is filled. */
  ECAT63_attnheader *h, 
  /** Float pointer to the address of the matrix data; any old contents are 
      not freed, therefore you must free the pointer after data is copied 
      elsewhere. */
  float **fdata
) {
  int i, ret, blockNr, pxlNr;
  char *mdata, *mptr, errmsg[128];
  float *_fdata, *fptr;
  short int *sptr;
  int *iptr;
  
  
  if(ECAT63_TEST) printf("ecat63ReadAttnMatrix(fp, %d, %d, hdr, fdata)\n",
    first_block, last_block);
  if(fp==NULL || first_block<=MatFirstDirBlk || h==NULL) {
    sprintf(ecat63errmsg, "invalid function parameter.\n");
    return(1);
  }
  *fdata=(float*)NULL;
  
  /* Read subheader */
  ret=ecat63ReadAttnheader(fp, first_block, h, ECAT63_TEST-2, errmsg);
  if(ret) {
    strcpy(ecat63errmsg, errmsg);
    return(5);
  }
  if(ECAT63_TEST>4) ecat63PrintAttnheader(h, stdout);
  pxlNr=h->dimension_1*h->dimension_2;
  if(pxlNr<=0) {
    sprintf(ecat63errmsg, "invalid matrix dimension.\n");
    return(6);  
  }
  
  /* Read matrix data */
  blockNr=last_block-first_block; if(blockNr<1) return(0);
  mdata=(char*)malloc(blockNr*MatBLKSIZE);
  if(mdata==NULL) {
    sprintf(ecat63errmsg, "cannot allocate memory.\n");
    return(8);  
  }
  mptr=mdata;
  ret=ecat63ReadMatdata(fp, first_block+1, blockNr, mptr, h->data_type);
  if(ret || mdata==NULL) {
    sprintf(ecat63errmsg, "cannot read matrix data (%d).\n", ret);
    free(mdata); return(9);
  }
  
  /* Allocate memory for float data */
  _fdata=(float*)malloc(pxlNr*sizeof(float));
  if(_fdata==NULL) {
    sprintf(ecat63errmsg, "cannot allocate memory.\n");
    free(mdata); return(11);  
  }

  /* Convert matrix data to floats */
  fptr=_fdata; mptr=mdata;
  if(h->data_type==BYTE_TYPE) {
    for(i=0; i<pxlNr; i++, mptr++, fptr++)
      *fptr=h->scale_factor*(float)(*mptr);
  } else if(h->data_type==VAX_I2 || h->data_type==SUN_I2) {
    for(i=0; i<pxlNr; i++, mptr+=2, fptr++) {
      sptr=(short int*)mptr;
      *fptr=h->scale_factor*(float)(*sptr);
    }
  } else if(h->data_type==VAX_I4 || h->data_type==SUN_I4) {
    for(i=0; i<pxlNr; i++, mptr+=4, fptr++) {
      iptr=(int*)mptr;
      *fptr=h->scale_factor*(float)(*iptr);
    }
  } else if(h->data_type==VAX_R4 || h->data_type==IEEE_R4) {
    memcpy(fptr, mptr, pxlNr*4);
    for(i=0; i<pxlNr; i++, fptr++) *fptr *= h->scale_factor;
  }
  free(mdata);
  *fdata=_fdata;

  return(0);
}
/*****************************************************************************/

/*****************************************************************************/
/*!
 * Reading ECAT 6.3 floats
 *
 * @param bufi pointer to 32-bit long data block
 * @param isvax 1 for VAX format
 * @param islittle 1 for little endian
 * @return read float value
 */
float ecat63rFloat(void *bufi, int isvax, int islittle) {
  union {unsigned int ul; float f;} t;

  memcpy(&t.ul, bufi, 4); if(t.ul==0) {return(0.0);}
  if(isvax) { /* if input is in VAX format */
    /* Swap words on i386 and bytes on SUN */
    if(islittle) swawip(&t.ul, 4); else swabip(&t.ul, 4);
    t.ul-=(2L<<23); /* subtract 2 from exp */
  } else { /* input is in i386 format */
    if(!islittle) swawbip(&t.ul, 4); /* Switch words and bytes on SUN */
  }
  return(t.f);
}

/*!
 * Reading and writing ECAT 6.3 32-bit ints.
 *   32-bit int format is same in VAX and i386
 *
 * @param bufi pointer to 32-bit long data block
 * @param isvax 1 for VAX format
 * @param islittle 1 for littel endian
 * @return read data as interger number
 */
int ecat63rInt(void *bufi, int isvax, int islittle) {
  int i;

  if(isvax==0) {}  // just to prevent compiler warning
  /* Swap both words and bytes on SUN */
  memcpy(&i, bufi, 4); if(!islittle) swawbip(&i, 4);
  return(i);
}
/*****************************************************************************/

/*****************************************************************************/
/*!
 * Returns the nr of bytes required for storage of one pixel of specified
 * data_type
 *
 * @param data_type data type code
 * @return number of bytes
 */
int ecat63pxlbytes(short int data_type) {
  int byteNr=0;
  switch(data_type) {
    case BYTE_TYPE: byteNr=1; break;
    case VAX_I2:
    case SUN_I2: byteNr=2; break;
    case VAX_I4:
    case VAX_R4:
    case IEEE_R4:
    case SUN_I4: byteNr=4; break;
  }
  return(byteNr);
}
/*****************************************************************************/

/*****************************************************************************/

