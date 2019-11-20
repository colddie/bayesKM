/// @file ecat63w.c
/// @author Vesa Oikonen
/// @brief Procedures for writing ECAT 6.3 matrix data.
///
///  Assumptions:
///  1. All data is always saved in little endian byte order (i386 and VAX).
///  2. Data is automatically saved in one of the little endian formats
///     as specified in header data_type.
///  3. VAX data can be saved correctly only in 2-byte formats.
///
/*****************************************************************************/
#include "libtpcimgio.h"
/*****************************************************************************/

/*****************************************************************************/
/*!
 * Write ECAT 6.3 main header.
 *
 * @param fp target file pointer
 * @param h Ecat 6.3 main header
 * @return 0, if ok, 1 invalid input, 2 failed to find block,
 * 3 failed to write block
 */
int ecat63WriteMainheader(FILE *fp, ECAT63_mainheader *h) {
  char buf[MatBLKSIZE];
  int i, little, tovax;


  if(ECAT63_TEST) printf("ecat63WriteMainheader()\n");
  little=little_endian();
  /* Clear buf */
  memset(buf, 0, MatBLKSIZE);
  /* Check arguments */
  if(fp==NULL || h->data_type<1 || h->data_type>7) return(1);
  if(h->data_type==VAX_I2 || h->data_type==VAX_I4 || h->data_type==VAX_R4)
    tovax=1; else tovax=0;

  /* Copy short ints to buf */
  memcpy(buf+50, &h->data_type, 2); memcpy(buf+48, &h->sw_version, 2);
  memcpy(buf+52, &h->system_type, 2); memcpy(buf+54, &h->file_type, 2);
  memcpy(buf+66, &h->scan_start_day, 2); memcpy(buf+68, &h->scan_start_month, 2);
  memcpy(buf+70, &h->scan_start_year, 2); memcpy(buf+72, &h->scan_start_hour, 2);
  memcpy(buf+74, &h->scan_start_minute, 2); 
  memcpy(buf+76, &h->scan_start_second, 2);
  memcpy(buf+134, &h->rot_source_speed, 2); memcpy(buf+136, &h->wobble_speed, 2);
  memcpy(buf+138, &h->transm_source_type, 2); 
  memcpy(buf+148, &h->transaxial_samp_mode, 2);
  memcpy(buf+150, &h->coin_samp_mode, 2); memcpy(buf+152, &h->axial_samp_mode, 2);
  memcpy(buf+158, &h->calibration_units, 2); 
  memcpy(buf+160, &h->compression_code, 2);
  memcpy(buf+350, &h->acquisition_type, 2); memcpy(buf+352, &h->bed_type, 2);
  memcpy(buf+354, &h->septa_type, 2); memcpy(buf+376, &h->num_planes, 2);
  memcpy(buf+378, &h->num_frames, 2); memcpy(buf+380, &h->num_gates, 2);
  memcpy(buf+382, &h->num_bed_pos, 2); memcpy(buf+452, &h->lwr_sctr_thres, 2);
  memcpy(buf+454, &h->lwr_true_thres, 2); memcpy(buf+456, &h->upr_true_thres, 2);
  memcpy(buf+472, h->fill2, 40);
  /* big to little endian if necessary */
  if(!little) swabip(buf, MatBLKSIZE);

  /* Copy floats to buf */
  ecat63wFloat(&h->isotope_halflife, buf+86, tovax, little);
  ecat63wFloat(&h->gantry_tilt, buf+122, tovax, little);
  ecat63wFloat(&h->gantry_rotation, buf+126, tovax, little);
  ecat63wFloat(&h->bed_elevation, buf+130, tovax, little);
  ecat63wFloat(&h->axial_fov, buf+140, tovax, little);
  ecat63wFloat(&h->transaxial_fov, buf+144, tovax, little);
  ecat63wFloat(&h->calibration_factor, buf+154, tovax, little);
  ecat63wFloat(&h->init_bed_position, buf+384, tovax, little);
  for(i=0; i<15; i++) ecat63wFloat(&h->bed_offset[i], buf+388+4*i, tovax, little);
  ecat63wFloat(&h->plane_separation, buf+448, tovax, little);
  ecat63wFloat(&h->init_bed_position, buf+458, tovax, little);

  /* Copy chars */
  /*memcpy(buf+0, h->ecat_format, 14);*/
  memcpy(buf+14, h->fill1, 14);
  memcpy(buf+28, h->original_file_name, 20); memcpy(buf+56, h->node_id, 10);
  memcpy(buf+78, h->isotope_code, 8); memcpy(buf+90, h->radiopharmaceutical, 32);
  memcpy(buf+162, h->study_name, 12); memcpy(buf+174, h->patient_id, 16);
  memcpy(buf+190, h->patient_name, 32); buf[222]=h->patient_sex;
  memcpy(buf+223, h->patient_age, 10); memcpy(buf+233, h->patient_height, 10);
  memcpy(buf+243, h->patient_weight, 10); buf[253]=h->patient_dexterity;
  memcpy(buf+254, h->physician_name, 32); memcpy(buf+286, h->operator_name, 32);
  memcpy(buf+318, h->study_description, 32); memcpy(buf+356, h->facility_name, 20);
  memcpy(buf+462, h->user_process_code, 10);

  /* Write main header */
  fseek(fp, 0*MatBLKSIZE, SEEK_SET); if(ftell(fp)!=0*MatBLKSIZE) return(2);
  if(fwrite(buf, 1, 1*MatBLKSIZE, fp) != 1*MatBLKSIZE) return(3);

  return(0);
}
/*****************************************************************************/

/*****************************************************************************/
/*!
 * Write ECAT 6.3 image header
 *
 * @param fp target file pointer
 * @param block block number [>= 3]
 * @param h Ecat 6.3 image header
 * @return 0, if ok, 1 invalid input, 2 failed to find block, 
 * 3 failed to write block
 */
int ecat63WriteImageheader(FILE *fp, int block, ECAT63_imageheader *h) {
  char buf[MatBLKSIZE];
  int i, little, tovax;


  if(ECAT63_TEST) printf("ecat63WriteImageheader(fp, %d, ih)\n", block);
  little=little_endian();
  /* Clear buf */
  memset(buf, 0, MatBLKSIZE);
  /* Check arguments */
  if(fp==NULL || block<3 || h->data_type<1 || h->data_type>7) return(1);
  if(h->data_type==VAX_I2 || h->data_type==VAX_I4 || h->data_type==VAX_R4)
    tovax=1; else tovax=0;

  /* Copy short ints to buf */
  memcpy(buf+126, &h->data_type, 2); 
  memcpy(buf+128, &h->num_dimensions, 2);
  memcpy(buf+132, &h->dimension_1, 2); 
  memcpy(buf+134, &h->dimension_2, 2);
  memcpy(buf+176, &h->image_min, 2); 
  memcpy(buf+178, &h->image_max, 2);
  memcpy(buf+200, &h->slice_location, 2); 
  memcpy(buf+202, &h->recon_start_hour, 2);
  memcpy(buf+204, &h->recon_start_min, 2); 
  memcpy(buf+206, &h->recon_start_sec, 2);
  memcpy(buf+236, &h->filter_code, 2); 
  memcpy(buf+376, &h->processing_code, 2);
  memcpy(buf+380, &h->quant_units, 2); 
  memcpy(buf+382, &h->recon_start_day, 2);
  memcpy(buf+384, &h->recon_start_month, 2); 
  memcpy(buf+386, &h->recon_start_year, 2); 
  memcpy(buf+460, h->fill2, 52);
  /* big to little endian if necessary */
  if(!little) swabip(buf, MatBLKSIZE);

  /* Copy floats to buf */
  ecat63wFloat(&h->x_origin, buf+160, tovax, little);
  ecat63wFloat(&h->y_origin, buf+164, tovax, little);
  ecat63wFloat(&h->recon_scale, buf+168, tovax, little);
  ecat63wFloat(&h->quant_scale, buf+172, tovax, little);
  ecat63wFloat(&h->pixel_size, buf+184, tovax, little);
  ecat63wFloat(&h->slice_width, buf+188, tovax, little);
  ecat63wFloat(&h->image_rotation, buf+296, tovax, little);
  ecat63wFloat(&h->plane_eff_corr_fctr, buf+300, tovax, little);
  ecat63wFloat(&h->decay_corr_fctr, buf+304, tovax, little);
  ecat63wFloat(&h->loss_corr_fctr, buf+308, tovax, little);
  ecat63wFloat(&h->intrinsic_tilt, buf+312, tovax, little);
  ecat63wFloat(&h->ecat_calibration_fctr, buf+388, tovax, little);
  ecat63wFloat(&h->well_counter_cal_fctr, buf+392, tovax, little);
  for(i=0; i<6; i++) 
    ecat63wFloat(&h->filter_params[i], buf+396+4*i, tovax, little);

  /* Copy ints to buf */
  ecat63wInt(&h->frame_duration, buf+192, tovax, little);
  ecat63wInt(&h->frame_start_time, buf+196, tovax, little);
  ecat63wInt(&h->scan_matrix_num, buf+238, tovax, little);
  ecat63wInt(&h->norm_matrix_num, buf+242, tovax, little);
  ecat63wInt(&h->atten_cor_mat_num, buf+246, tovax, little);

  /* Copy chars */
  memcpy(buf+0, h->fill1, 126); 
  memcpy(buf+420, h->annotation, 40);

  /* Write subheader */
  fseek(fp, (block-1)*MatBLKSIZE, SEEK_SET); 
  if(ftell(fp)!=(block-1)*MatBLKSIZE) return(2);
  if(fwrite(buf, 1, 1*MatBLKSIZE, fp) != 1*MatBLKSIZE) return(3);

  return(0);
}
/*****************************************************************************/

/*****************************************************************************/
/*!
 * Write ECAT 6.3 attenuation header
 *
 * @param fp target file pointer
 * @param block block number [>=3]
 * @param h Ecat 6.3 attenuation header
 * @return 0 if ok, 1 invalid input, 2 failed to find block, 
 * 3 failed to write block
 */
int ecat63WriteAttnheader(FILE *fp, int block, ECAT63_attnheader *h) {
  unsigned char buf[MatBLKSIZE];
  int little, tovax;

  if(ECAT63_TEST) printf("ecat63WriteAttnheader(fp, %d, ah)\n", block);
  little=little_endian();
  /* Clear buf */
  memset(buf, 0, MatBLKSIZE);
  /* Check arguments */
  if(fp==NULL || block<3 || h->data_type<1 || h->data_type>7) return(1);
  if(h->data_type==VAX_I2 || h->data_type==VAX_I4 || h->data_type==VAX_R4)
    tovax=1; else tovax=0;

  /* Copy short ints to buf */
  memcpy(buf+126, &h->data_type, 2); 
  memcpy(buf+128, &h->attenuation_type, 2);
  memcpy(buf+132, &h->dimension_1, 2); 
  memcpy(buf+134, &h->dimension_2, 2);
  /* big to little endian if necessary */
  if(!little) swabip(buf, MatBLKSIZE);

  /* Copy floats to buf */
  ecat63wFloat(&h->scale_factor, buf+182, tovax, little);
  ecat63wFloat(&h->x_origin, buf+186, tovax, little);
  ecat63wFloat(&h->y_origin, buf+190, tovax, little);
  ecat63wFloat(&h->x_radius, buf+194, tovax, little);
  ecat63wFloat(&h->y_radius, buf+198, tovax, little);
  ecat63wFloat(&h->tilt_angle, buf+202, tovax, little);
  ecat63wFloat(&h->attenuation_coeff, buf+206, tovax, little);
  ecat63wFloat(&h->sample_distance, buf+210, tovax, little);

  /* Write subheader */
  fseek(fp, (block-1)*MatBLKSIZE, SEEK_SET);
  if(ftell(fp)!=(block-1)*MatBLKSIZE) return(2);
  if(fwrite(buf, 1, 1*MatBLKSIZE, fp) != 1*MatBLKSIZE) return(3);

  return(0);
}
/*****************************************************************************/

/*****************************************************************************/
/*!
 * Write ECAT 6.3 scan header
 * 
 * @param fp target file pointer
 * @param block block number [>=3]
 * @param h Ecat 6.3 scan header
 * @return 0 if ok, 1 invalid input, 2 failed to find block,
 * 3 failed to write block
 */
int ecat63WriteScanheader(FILE *fp, int block, ECAT63_scanheader *h) {
  unsigned char buf[MatBLKSIZE];
  int i, little, tovax;


  if(ECAT63_TEST) printf("ecat63WriteScanheader(fp, %d, ih)\n", block);
  little=little_endian();
  /* Clear buf */
  memset(buf, 0, MatBLKSIZE);
  /* Check arguments */
  if(fp==NULL || block<3 || h->data_type<1 || h->data_type>7) return(1);
  if(h->data_type==VAX_I2 || h->data_type==VAX_I4 || h->data_type==VAX_R4)
    tovax=1; else tovax=0;

  /* Copy short ints to buf */
  memcpy(buf+126, &h->data_type, 2);
  memcpy(buf+132, &h->dimension_1, 2); memcpy(buf+134, &h->dimension_2, 2);
  memcpy(buf+136, &h->smoothing, 2); memcpy(buf+138, &h->processing_code, 2);
  memcpy(buf+170, &h->frame_duration_sec, 2);
  memcpy(buf+192, &h->scan_min, 2); memcpy(buf+194, &h->scan_max, 2);
  memcpy(buf+468, h->fill2, 44);
  /* big to little endian if necessary */
  if(!little) swabip(buf, MatBLKSIZE);

  /* Copy floats to buf */
  ecat63wFloat(&h->sample_distance, buf+146, tovax, little);
  ecat63wFloat(&h->isotope_halflife, buf+166, tovax, little);
  ecat63wFloat(&h->scale_factor, buf+182, tovax, little);
  for(i=0; i<16; i++) {
    ecat63wFloat(&h->cor_singles[i], buf+316+4*i, tovax, little);
    ecat63wFloat(&h->uncor_singles[i], buf+380+4*i, tovax, little);
  }
  ecat63wFloat(&h->tot_avg_cor, buf+444, tovax, little);
  ecat63wFloat(&h->tot_avg_uncor, buf+448, tovax, little);
  ecat63wFloat(&h->loss_correction_fctr, buf+464, tovax, little);

  /* Copy ints to buf */
  ecat63wInt(&h->gate_duration, buf+172, tovax, little);
  ecat63wInt(&h->r_wave_offset, buf+176, tovax, little);
  ecat63wInt(&h->prompts, buf+196, tovax, little);
  ecat63wInt(&h->delayed, buf+200, tovax, little);
  ecat63wInt(&h->multiples, buf+204, tovax, little);
  ecat63wInt(&h->net_trues, buf+208, tovax, little);
  ecat63wInt(&h->total_coin_rate, buf+452, tovax, little);
  ecat63wInt(&h->frame_start_time, buf+456, tovax, little);
  ecat63wInt(&h->frame_duration, buf+460, tovax, little);

  /* Copy chars */
  memcpy(buf+0, h->fill1, 126);

  /* Write subheader */
  fseek(fp, (block-1)*MatBLKSIZE, SEEK_SET); 
  if(ftell(fp)!=(block-1)*MatBLKSIZE) return(2);
  if(fwrite(buf, 1, 1*MatBLKSIZE, fp) != 1*MatBLKSIZE) return(3);

  return(0);
}
/*****************************************************************************/

/*****************************************************************************/
/*!
 * Write ECAT 6.3 normalization header
 *
 * @param fp target file pointer
 * @param block block number [>=3]
 * @param h Ecat 6.3 normalization header
 * @return 0 if ok, 1 invalid input, 2 failed to find block,
 * 3 failed to write block
 */
int ecat63WriteNormheader(FILE *fp, int block, ECAT63_normheader *h)
{
  unsigned char buf[MatBLKSIZE];
  int little, tovax;

  if(ECAT63_TEST) printf("ecat63WriteNormheader(fp, %d, nh)\n", block);
  little=little_endian();
  /* Clear buf */
  memset(buf, 0, MatBLKSIZE);
  /* Check arguments */
  if(fp==NULL || block<3 || h->data_type<1 || h->data_type>7) return(1);
  if(h->data_type==VAX_I2 || h->data_type==VAX_I4 || h->data_type==VAX_R4)
    tovax=1; else tovax=0;

  /* Copy short ints to buf */
  memcpy(buf+126, &h->data_type, 2);
  memcpy(buf+132, &h->dimension_1, 2); 
  memcpy(buf+134, &h->dimension_2, 2);
  memcpy(buf+372, &h->norm_hour, 2);
  memcpy(buf+376, &h->norm_minute, 2);
  memcpy(buf+380, &h->norm_second, 2);
  memcpy(buf+384, &h->norm_day, 2);
  memcpy(buf+388, &h->norm_month, 2);
  memcpy(buf+392, &h->norm_year, 2);
  /* big to little endian if necessary */
  if(!little) swabip(buf, MatBLKSIZE);

  /* Copy floats to buf */
  ecat63wFloat(&h->scale_factor, buf+182, tovax, little);
  ecat63wFloat(&h->fov_source_width, buf+198, tovax, little);

  /* Write subheader */
  fseek(fp, (block-1)*MatBLKSIZE, SEEK_SET);
  if(ftell(fp)!=(block-1)*MatBLKSIZE) return(2);
  if(fwrite(buf, 1, 1*MatBLKSIZE, fp) != 1*MatBLKSIZE) return(3);

  return(0);
}
/*****************************************************************************/

/*****************************************************************************/
/*!
 * Create a new ECAT 6.3 file and return file pointer
 *  or NULL in case of an error.
 *  If file exists, it is renamed as fname% if possible.
 *  Directory list is written in big endian byte order.
 *
 * @param fname file name
 * @param h Ecat 6.3 main header
 * @return opened file pointer, or NULL in case of failure
 */
FILE *ecat63Create(const char *fname, ECAT63_mainheader *h) {
  FILE *fp;
  char tmp[FILENAME_MAX];
  int buf[MatBLKSIZE/4];

  if(ECAT63_TEST) printf("ecat63Create()\n");
  /* Check the arguments */
  if(fname==NULL || h==NULL) return(NULL);
  /* Check if file exists; backup, if necessary */
  if(access(fname, 0) != -1) {
    strcpy(tmp, fname); strcat(tmp, BACKUP_EXTENSION);
    if(access(tmp, 0) != -1) remove(tmp);
    if(ECAT63_TEST) printf("Renaming %s -> %s\n", fname, tmp);
    rename(fname, tmp);
  }
  /* Open file */
  fp=fopen(fname, "wb+"); if(fp==NULL) return(fp);
  /* Write main header */
  if(ecat63WriteMainheader(fp, h)) return(NULL);
  /* Construct an empty matrix list ; convert to little endian if necessary */
  memset(buf, 0, MatBLKSIZE);  
  buf[0]=31; buf[1]=2; if(!little_endian()) swawbip(buf, MatBLKSIZE);
  /* Write data buffer */
  fseek(fp, (MatFirstDirBlk-1)*MatBLKSIZE, SEEK_SET);
  if(ftell(fp)!=(MatFirstDirBlk-1)*MatBLKSIZE) return(NULL);
  if(fwrite(buf, 4, MatBLKSIZE/4, fp) != MatBLKSIZE/4) return(NULL);
  /* OK, then return file pointer */
  return(fp);
}
/*****************************************************************************/

/*****************************************************************************/
/*!
 * Write ECAT 6.3 image matrix header and data
 *
 * @param fp target file pointer
 * @param matnum matrix number [1..number of matrixes]
 * @param h Ecat 6.3 image header
 * @param data pointer to data that is written
 * @return 0 if ok, 1 invalid input or invalid image dimensions,
 * 2 failed to resolve data type
 * 3 too little data size, 4 failed to resolve next block size in file
 */
int ecat63WriteImage(FILE *fp, int matnum, ECAT63_imageheader *h, void *data) {
  int nxtblk, blkNr, data_size, pxlNr, pxlSize, ret;

  if(ECAT63_TEST) printf("ecat63WriteImage(fp, %d, ih, data)\n", matnum);
  if(fp==NULL || matnum<1 || h==NULL || data==NULL) return(1);
  /* nr of pixels */
  pxlNr=h->dimension_1*h->dimension_2; if(pxlNr<1) return(2);
  /* mem taken by one pixel */
  switch(h->data_type) {
    case BYTE_TYPE: pxlSize=1;
      break;
    case VAX_I2:
    case SUN_I2: pxlSize=2;
      break;
    case VAX_I4: return(3);
    case VAX_R4: return(3);
    case IEEE_R4:
    case SUN_I4: pxlSize=4;
      break;
    default: return(2);
  }
  /* mem taken by all pixels */
  data_size=pxlNr*pxlSize;
  /* block nr taken by all pixels */
  blkNr=(data_size+MatBLKSIZE-1)/MatBLKSIZE; if(blkNr<1) return(3);
  /* Get block number for matrix header and data */
  nxtblk=ecat63Matenter(fp, matnum, blkNr); if(nxtblk<1) return(4);
  if(ECAT63_TEST) printf("  block=%d\n", nxtblk);
  /* Write header */
  ret=ecat63WriteImageheader(fp, nxtblk, h); if(ret) return(40+ret);
  /* Write matrix data */
  ret=ecat63WriteMatdata(fp, nxtblk+1, data, pxlNr, pxlSize);
  if(ret) return(50+ret);
  return 0;
}
/*****************************************************************************/

/*****************************************************************************/
/*!
 * Write ECAT 6.3 sinogram matrix header and data
 *
 * @param fp target file pointer
 * @param matnum matrix number [1..number of matrixes]
 * @param h Ecat 6.3 scan header
 * @param data pointer to data that is written
 * @return 0 if ok, 1 invalid input or invalid image dimensions,
 * 2 failed to resolve data type
 * 3 too little data size, 4 failed to resolve next block size in file
 */
int ecat63WriteScan(FILE *fp, int matnum, ECAT63_scanheader *h, void *data) {
  int nxtblk, blkNr, data_size, pxlNr, pxlSize, ret;

  if(ECAT63_TEST) printf("ecat63WriteScan(fp, %d, sh, data)\n", matnum);
  if(fp==NULL || matnum<1 || h==NULL || data==NULL) return(1);
  /* nr of pixels */
  pxlNr=h->dimension_1*h->dimension_2; if(pxlNr<1) return(1);
  /* mem taken by one pixel */
  switch(h->data_type) {
    case BYTE_TYPE: pxlSize=1;
      break;
    case VAX_I2:
    case SUN_I2: pxlSize=2;
      break;
    case VAX_I4: return(3);
    case VAX_R4: return(3);
    case IEEE_R4:
    case SUN_I4: pxlSize=4;
      break;
    default: return(2);
  }
  /* mem taken by all pixels */
  data_size=pxlNr*pxlSize;
  /* block nr taken by all pixels */
  blkNr=(data_size+MatBLKSIZE-1)/MatBLKSIZE; if(blkNr<1) return(3);
  /* Get block number for matrix header and data */
  nxtblk=ecat63Matenter(fp, matnum, blkNr); if(nxtblk<1) return(4);
  if(ECAT63_TEST) printf("  block=%d\n", nxtblk);
  /* Write header */
  ret=ecat63WriteScanheader(fp, nxtblk, h); if(ret) return(40+ret);
  /* Write matrix data */
  ret=ecat63WriteMatdata(fp, nxtblk+1, data, pxlNr, pxlSize);
  if(ret) return(50+ret);
  return 0;
}
/*****************************************************************************/

/*****************************************************************************/
/*!
 * Write ECAT 6.3 normalization matrix header and data
 *
 * @param fp target file pointer
 * @param matnum matrix number [1..number of matrixes]
 * @param h Ecat 6.3 normalization header
 * @param data pointer to data that is written
 * @return 0 if ok, 1 invalid input or invalid image dimensions,
 * 2 failed to resolve data type
 * 3 too little data size, 4 failed to resolve next block size in file
 */
int ecat63WriteNorm(FILE *fp, int matnum, ECAT63_normheader *h, void *data) {
  int nxtblk, blkNr, data_size, pxlNr, pxlSize, ret;

  if(ECAT63_TEST) printf("ecat63WriteNorm(fp, %d, nh, data)\n", matnum);
  if(fp==NULL || matnum<1 || h==NULL || data==NULL) return(1);
  /* nr of pixels */
  pxlNr=h->dimension_1*h->dimension_2; if(pxlNr<1) return(1);
  /* mem taken by one pixel */
  switch(h->data_type) {
    case BYTE_TYPE: pxlSize=1;
      break;
    case VAX_I2:
    case SUN_I2: pxlSize=2;
      break;
    case VAX_I4: return(3);
    case VAX_R4: return(3);
    case IEEE_R4:
    case SUN_I4: pxlSize=4;
      break;
    default: return(2);
  }
  /* mem taken by all pixels */
  data_size=pxlNr*pxlSize;
  /* block nr taken by all pixels */
  blkNr=(data_size+MatBLKSIZE-1)/MatBLKSIZE; if(blkNr<1) return(3);
  /* Get block number for matrix header and data */
  nxtblk=ecat63Matenter(fp, matnum, blkNr); if(nxtblk<1) return(4);
  if(ECAT63_TEST) printf("  block=%d\n", nxtblk);
  /* Write header */
  ret=ecat63WriteNormheader(fp, nxtblk, h); if(ret) return(40+ret);
  /* Write matrix data */
  ret=ecat63WriteMatdata(fp, nxtblk+1, data, pxlNr, pxlSize);
  if(ret) return(50+ret);
  return 0;
}
/*****************************************************************************/

/*****************************************************************************/
/*!
 * Write ECAT 6.3 attenuation matrix header and data
 *
 * @param fp target file pointer
 * @param matnum matrix number [1..number of matrixes]
 * @param h Ecat 6.3 attenuation header
 * @param data pointer to data that is written
 * @return 0 if ok, 1 invalid input or invalid image dimensions,
 * 2 failed to resolve data type
 * 3 too little data size, 4 failed to resolve next block size in file
 */
int ecat63WriteAttn(FILE *fp, int matnum, ECAT63_attnheader *h, void *data) {
  int nxtblk, blkNr, data_size, pxlNr, pxlSize, ret;

  if(ECAT63_TEST) printf("ecat63WriteAttn(fp, %d, ah, data)\n", matnum);
  if(fp==NULL || matnum<1 || h==NULL || data==NULL) return(1);
  /* nr of pixels */
  pxlNr=h->dimension_1*h->dimension_2; if(pxlNr<1) return(1);
  /* mem taken by one pixel */
  switch(h->data_type) {
    case BYTE_TYPE: pxlSize=1;
      break;
    case VAX_I2:
    case SUN_I2: pxlSize=2;
      break;
    case VAX_I4: return(3);
    case VAX_R4: return(3);
    case IEEE_R4:
    case SUN_I4: pxlSize=4;
      break;
    default: return(2);
  }
  /* mem taken by all pixels */
  data_size=pxlNr*pxlSize;
  /* block nr taken by all pixels */
  blkNr=(data_size+MatBLKSIZE-1)/MatBLKSIZE; if(blkNr<1) return(3);
  /* Get block number for matrix header and data */
  nxtblk=ecat63Matenter(fp, matnum, blkNr); if(nxtblk<1) return(4);
  if(ECAT63_TEST) printf("  block=%d\n", nxtblk);
  /* Write header */
  ret=ecat63WriteAttnheader(fp, nxtblk, h); if(ret) return(40+ret);
  /* Write matrix data */
  ret=ecat63WriteMatdata(fp, nxtblk+1, data, pxlNr, pxlSize);
  if(ret) return(50+ret);
  return 0;
}
/*****************************************************************************/

/*****************************************************************************/
/*!
 * Write ECAT 6.3 matrix data to a specified file position.
 *   Data does not need to be allocated for full blocks.
 *   Data must be represented in current machines byte order, and it is
 *   always saved in big endian byte order.
 *   Give also nr of pixels and byte size of one pixel.
 *
 * @param fp target file pointer
 * @param strtblk starting image block [>=1]
 * @param data pointer to data that is written
 * @param pxlNr number of items to be written [>=1]
 * @param pxlSize size of one data item in bytes [>=1]
 * @return 0 if ok, 1 invalid input, 2 failed to find starting block, 
 * 3 failed to write data
 */
int ecat63WriteMatdata(
  FILE *fp, int strtblk, char *data, int pxlNr, int pxlSize) 
{
  unsigned char buf[MatBLKSIZE];
  char *dptr;
  int i, blkNr, dataSize, byteNr;

  if(ECAT63_TEST) 
    printf("ecat63WriteMatdata(fp, %d, data, %d, %d)\n", strtblk, pxlNr, pxlSize);
  if(fp==NULL || strtblk<1 || data==NULL || pxlNr<1 || pxlSize<1) return(1);
  memset(buf, 0, MatBLKSIZE);
  dataSize=pxlNr*pxlSize; if(dataSize<1) return(1);
  /* block nr taken by all pixels */
  blkNr=(dataSize+MatBLKSIZE-1)/MatBLKSIZE; if(blkNr<1) return(1);
  if(ECAT63_TEST>1) printf("    blkNr=%d\n", blkNr);
  /* Search the place for writing */
  fseek(fp, (strtblk-1)*MatBLKSIZE, SEEK_SET); 
  if(ftell(fp)!=(strtblk-1)*MatBLKSIZE) return(2);
  /* Save blocks one at a time */
  for(i=0, dptr=data; i<blkNr && dataSize>0; i++) {
    byteNr=(dataSize<MatBLKSIZE)?dataSize:MatBLKSIZE;
    memcpy(buf, dptr, byteNr);
    /* Change matrix byte order in big endian platforms */
    if(!little_endian()) {
      if(pxlSize==2) swabip(buf, byteNr);
      else if(pxlSize==4) swawbip(buf, byteNr);
    }
    /* Write block */
    if(fwrite(buf, 1, MatBLKSIZE, fp)!=MatBLKSIZE) return(3);
    /* Prepare for the next block */
    dptr+=byteNr;
    dataSize-=byteNr;
  } /* next block */
  return(0);
}
/*****************************************************************************/

/*****************************************************************************/
/*!
 * Check if pixel float values need to be scaled to be saved as short ints,
 *   or if they are already all very close to integers.
 *
 * @param amax absolute maximum value
 * @param data float array
 * @param nr number of float values in flaot array
 * @return 1 if scaling is necessary, and 0 if not.
*/
int ecat63_is_scaling_needed(float amax, float *data, int nr) {
  int i;
  double d;

  if(nr<1 || data==NULL) return(0);
  /* scaling is necessary if all values are between -1 - 1 */
  if(amax<0.9999) return(1);
  /* Lets check first if at least the max value is close to integers or not */
  if(modf(amax, &d)>0.0001) return(1);
  /* if it is, then check all pixels */
  for(i=0; i<nr; i++) if(modf(*data++, &d)>0.0001) return(1);
  return(0);
}
/*****************************************************************************/

/*****************************************************************************/
/*!
 * Write ECAT 6.3 image matrix header and data
 *
 * @param fp target file pointer
 * @param matnum
 * @param h Ecat 6.3 image header
 * @param fdata
 * @return 0 if ok, 1 invalid input, 3 invalid matrix dimensions, 
 * 4 invalid block number, 5 failed to allocate memory, 
 * 8 failed to resolve new matrix block number, 
 * 10 failed to write image sub header, 
 * 13 failed to write matrix data
 */
int ecat63WriteImageMatrix(
  FILE *fp, int matnum, ECAT63_imageheader *h, float *fdata) 
{
  int i, nxtblk, blkNr, data_size, pxlNr, ret;
  float *fptr, fmin, fmax, g, f;
  char *mdata, *mptr;
  short int *sptr;



  if(ECAT63_TEST) printf("ecat63WriteImageMatrix(fp, %d, h, data)\n", matnum);
  if(fp==NULL || matnum<1 || h==NULL || fdata==NULL) {
    sprintf(ecat63errmsg, "invalid function parameter.\n");
    return(1);
  }
  /* nr of pixels */
  pxlNr=h->dimension_1*h->dimension_2;
  if(pxlNr<1) {
    sprintf(ecat63errmsg, "invalid matrix dimension.\n");
    return(3);
  }
  /* How much memory is needed for ALL pixels */
  data_size=pxlNr*ecat63pxlbytes(h->data_type);
  /* block nr taken by all pixels */
  blkNr=(data_size+MatBLKSIZE-1)/MatBLKSIZE; if(blkNr<1) {
    sprintf(ecat63errmsg, "invalid block number.\n");
    return(4);
  }
  /* Allocate memory for matrix data */
  mdata=(char*)calloc(blkNr, MatBLKSIZE); if(mdata==NULL) {
    sprintf(ecat63errmsg, "out of memory.\n");
    return(5);
  }
  /* Search for min and max for calculation of scale factor */
  fMinMaxFin(fdata, pxlNr, &fmin, &fmax);
/*
  fptr=fdata; fmin=fmax=*fptr;
  for(i=0; i<pxlNr; i++, fptr++) {
    if(*fptr>fmax) fmax=*fptr; else if(*fptr<fmin) fmin=*fptr;
  }
*/
  if(fabs(fmin)>fabs(fmax)) g=fabs(fmin); else g=fabs(fmax);
  if(g>0) f=32766./g; else f=1.0;
  /* Check if pixels values can be left as such with scale_factor = 1 */
  fptr=fdata;
  if(f>=1.0 && ecat63_is_scaling_needed(g, fptr, pxlNr)==0) f=1.0;
  /* Scale matrix data to shorts */
  h->quant_scale=1.0/f;
  sptr=(short int*)mdata; fptr=fdata;
  for(i=0; i<pxlNr; i++, sptr++, fptr++) *sptr=(short int)temp_roundf(f*(*fptr));
  /* Set header short min & max */
  h->image_min=(short int)temp_roundf(f*fmin);
  h->image_max=(short int)temp_roundf(f*fmax);
  /* Get block number for matrix header and data */
  nxtblk=ecat63Matenter(fp, matnum, blkNr); if(nxtblk<1) {
    sprintf(ecat63errmsg, "cannot determine matrix block (%d).\n", -nxtblk);
    free(mdata); return(8);
  }
  if(ECAT63_TEST>2) printf("  block=%d fmin=%g fmax=%g\n", nxtblk, fmin, fmax);
  /* Write header */
  ret=ecat63WriteImageheader(fp, nxtblk, h); if(ret) {
    sprintf(ecat63errmsg, "cannot write subheader (%d).\n", ret);
    free(mdata); return(10);
  }
  /* Write matrix data */
  mptr=mdata;
  ret=ecat63WriteMatdata(fp, nxtblk+1, mptr, pxlNr, ecat63pxlbytes(h->data_type));
  free(mdata);
  if(ret) {
    sprintf(ecat63errmsg, "cannot write matrix data (%d).\n", ret);
    return(13);
  }
  return(0);
}
/*****************************************************************************/

/*****************************************************************************/
/*!
 * Write ECAT 6.3 sinogram matrix header and data
 *
 * @param fp target file pointer
 * @param matnum matrix number [1..number of matrixes]
 * @param h Ecat 6.3 scan header
 * @param fdata matrix data
 * @return 0 if ok, 1 invalid input, 3 invalid matrix dimension, 
 * 4 invalid block number, 5 failed to allocate memory for data, 
 * 8 failed to resolve next block number, 10 cannot write sub header, 
 * 13 failed to write data
 */
int ecat63WriteScanMatrix(
  FILE *fp, int matnum, ECAT63_scanheader *h, float *fdata) 
{
  int i, nxtblk, blkNr, data_size, pxlNr, ret;
  float *fptr, fmin, fmax, g, f;
  char *mdata, *mptr;
  short int *sptr;


  if(ECAT63_TEST) printf("ecat63WriteScanMatrix(fp, %d, h, data)\n", matnum);
  if(fp==NULL || matnum<1 || h==NULL || fdata==NULL) {
    sprintf(ecat63errmsg, "invalid function parameter.\n");
    return(1);
  }
  /* nr of pixels */
  pxlNr=h->dimension_1*h->dimension_2;
  if(pxlNr<1) {
    sprintf(ecat63errmsg, "invalid matrix dimension.\n");
    return(3);
  }
  /* How much memory is needed for ALL pixels */
  data_size=pxlNr*ecat63pxlbytes(h->data_type);
  /* block nr taken by all pixels */
  blkNr=(data_size+MatBLKSIZE-1)/MatBLKSIZE; if(blkNr<1) {
    sprintf(ecat63errmsg, "invalid block number.\n");
    return(4);
  }
  /* Allocate memory for matrix data */
  mdata=(char*)calloc(blkNr, MatBLKSIZE); if(mdata==NULL) {
    sprintf(ecat63errmsg, "out of memory.\n");
    return(5);
  }
  /* Search for min and max for calculation of scale factor */
  fMinMaxFin(fdata, pxlNr, &fmin, &fmax);
/*
  fptr=fdata; fmin=fmax=*fptr;
  for(i=0; i<pxlNr; i++, fptr++) {
    if(*fptr>fmax) fmax=*fptr; else if(*fptr<fmin) fmin=*fptr;
  }
*/
  if(fabs(fmin)>fabs(fmax)) g=fabs(fmin); else g=fabs(fmax);
  if(g>0) f=32766./g; else f=1.0;
  /* Check if pixels values can be left as such with scale_factor = 1 */
  fptr=fdata;
  if(f>=1.0 && ecat63_is_scaling_needed(g, fptr, pxlNr)==0) f=1.0;
  /* Scale matrix data to shorts */
  h->scale_factor=1.0/f;
  sptr=(short int*)mdata; fptr=fdata;
  for(i=0; i<pxlNr; i++, sptr++, fptr++) *sptr=(short int)temp_roundf(f*(*fptr));
  /* Set header short min & max */
  h->scan_min=(short int)temp_roundf(f*fmin);
  h->scan_max=(short int)temp_roundf(f*fmax);
  /* Get block number for matrix header and data */
  nxtblk=ecat63Matenter(fp, matnum, blkNr); if(nxtblk<1) {
    sprintf(ecat63errmsg, "cannot determine matrix block (%d).\n", -nxtblk);
    free(mdata); return(8);
  }
  if(ECAT63_TEST>2) printf("  block=%d fmin=%g fmax=%g\n", nxtblk, fmin, fmax);
  /* Write header */
  ret=ecat63WriteScanheader(fp, nxtblk, h); if(ret) {
    sprintf(ecat63errmsg, "cannot write subheader (%d).\n", ret);
    free(mdata); return(10);
  }
  /* Write matrix data */
  mptr=mdata;
  ret=ecat63WriteMatdata(fp, nxtblk+1, mptr, pxlNr, ecat63pxlbytes(h->data_type));
  free(mdata);
  if(ret) {
    sprintf(ecat63errmsg, "cannot write matrix data (%d).\n", ret);
    return(13);
  }
  return(0);
}
/*****************************************************************************/

/*****************************************************************************/
/*!
 * Writing ECAT 6.3 floats
 *
 * @param bufi pointer to 4-byte long input (float data)
 * @param bufo pointer to 4-byte long output
 * @param tovax 1 for VAX format
 * @param islittle 1 for little endian
 */
void ecat63wFloat(float *bufi, void *bufo, int tovax, int islittle) {
  unsigned int ul;

  memcpy(&ul, bufi, 4); if(ul==0) {memcpy(bufo, bufi, 4); return;}
  if(tovax) { /* If VAX format is needed */
    ul+=(2L<<23); /* increase exp by 2 */
    /* Swap words on i386 and bytes on SUN */
    if(islittle) swawip(&ul, 4); else swabip(&ul, 4);
  } else {
    if(!islittle) swawbip(&ul, 4); /* Switch words and bytes on SUN */
  }
  memcpy(bufo, &ul, 4);
}
/*!
 * Writing ECAT 6.3 32-bit ints.
 *  32-bit int format is same in VAX and i386
 *
 * @param bufi pointer to 4-byte long input (integer data)
 * @param bufo pointer to 4-byte long output
 * @param tovax 1 for VAX format
 * @param islittle 1 for little endian
 */
void ecat63wInt(int *bufi, void *bufo, int tovax, int islittle) {
  int i;

  if(tovax==0) {} // to prevent compiler warning
  /* Swap both words and bytes on SUN */
  memcpy(&i, bufi, 4); if(!islittle) swawbip(&i, 4);
  memcpy(bufo, &i, 4);
}
/*****************************************************************************/

/*****************************************************************************/
/** @brief Convert scan_start_time in ECAT 6.3 main header into a 
 *  struct tm.
 *  @author Vesa Oikonen
 *  @return Returns pointer to the struct_tm, or null in case of an error.
 */
struct tm* ecat63ScanstarttimeToTm(
  /** Pointer to ECAT 6.3 main header */
  const ECAT63_mainheader *h,
  /** Pointer to pre-allocated struct tm. */
  struct tm *tm
) {
  if(tm==NULL || h==NULL) return(NULL);
  memset(tm, 0, sizeof(struct tm));
  tm->tm_mday=h->scan_start_day;
  tm->tm_mon=h->scan_start_month-1;
  tm->tm_year=h->scan_start_year-1900;
  tm->tm_hour=h->scan_start_hour;
  tm->tm_min=h->scan_start_minute;
  tm->tm_sec=h->scan_start_second;
  tm->tm_isdst=-1;
  if(timegm(tm)==-1) return(NULL);
  return(tm);
}
/*****************************************************************************/

/*****************************************************************************/
/** @brief Get calendar time from ECAT 6.3 main header.
 *  @author Vesa Oikonen
 *  @return Returns time_t, or -1 in case of an error.
 */
time_t ecat63Scanstarttime(
  /** Pointer to ECAT 6.3 main header */
  const ECAT63_mainheader *h
) {
  if(h==NULL) return((time_t)-1);
  struct tm tm;
  tm.tm_mday=h->scan_start_day;
  tm.tm_mon=h->scan_start_month-1;
  tm.tm_year=h->scan_start_year-1900;
  tm.tm_hour=h->scan_start_hour;
  tm.tm_min=h->scan_start_minute;
  tm.tm_sec=h->scan_start_second;
  tm.tm_isdst=-1;
  return(timegm(&tm));
}
/*****************************************************************************/

/*****************************************************************************/
