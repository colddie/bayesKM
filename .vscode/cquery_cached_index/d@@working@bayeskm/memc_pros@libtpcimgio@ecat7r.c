/// @file ecat7r.c
/// @author Vesa Oikonen
/// @brief Functions for reading ECAT 7.x format.
///
/*****************************************************************************/
#include "libtpcimgio.h"
/*****************************************************************************/

/*****************************************************************************/
/**
 * Read ECAT 7.x main header
 * @return 0 if ok, 1 == invalid parameters, 2 == first header block not found,
 * 3 == header block not read properly
 */
int ecat7ReadMainheader(
  /** input file pointer */
  FILE *fp, 
  /** Ecat7 main header */
  ECAT7_mainheader *h
) {
  unsigned char buf[MatBLKSIZE];
  int little; /* 1 if current platform is little endian (i386), else 0 */
  struct tm st;

  if(ECAT7_TEST) printf("ecat7ReadMainheader()\n");
  if(fp==NULL || h==NULL) return(1);
  little=little_endian(); if(ECAT7_TEST) printf("little=%d\n", little);

  /* Seek the first block */
  fseek(fp, 0, SEEK_SET); if(ftell(fp)!=0) return(2); 
  /* Read the header block */
  if(fread(buf, MatBLKSIZE, 1, fp)<1) return(3);
  
  /* Copy the header fields and swap if necessary */
  memcpy(&h->magic_number, buf+0, 14);
  memcpy(&h->original_file_name, buf+14, 32);
  if(little) swabip(buf+46, 2);
  memcpy(&h->sw_version, buf+46, 2);
  if(little) swabip(buf+48, 2);
  memcpy(&h->system_type, buf+48, 2);
  if(little) swabip(buf+50, 2);
  memcpy(&h->file_type, buf+50, 2);
  memcpy(&h->serial_number, buf+52, 10);
  if(little) swawbip(buf+62, 4);
  memcpy(&h->scan_start_time, buf+62, 4);
  //printf("ecat7ReadMainheader(): scan_start_time := %d\n", h->scan_start_time);
  memcpy(&h->isotope_name, buf+66, 8);
  if(little) swawbip(buf+74, 4);
  memcpy(&h->isotope_halflife, buf+74, 4);
  memcpy(&h->radiopharmaceutical, buf+78, 32);
  if(little) swawbip(buf+110, 4);
  memcpy(&h->gantry_tilt, buf+110, 4);
  if(little) swawbip(buf+114, 4);
  memcpy(&h->gantry_rotation, buf+114, 4);
  if(little) swawbip(buf+118, 4);
  memcpy(&h->bed_elevation, buf+118, 4);
  if(little) swawbip(buf+122, 4);
  memcpy(&h->intrinsic_tilt, buf+122, 4);
  if(little) swabip(buf+126, 2);
  memcpy(&h->wobble_speed, buf+126, 2);
  if(little) swabip(buf+128, 2);
  memcpy(&h->transm_source_type, buf+128, 2);
  if(little) swawbip(buf+130, 4);
  memcpy(&h->distance_scanned, buf+130, 4);
  if(little) swawbip(buf+134, 4);
  memcpy(&h->transaxial_fov, buf+134, 4);
  if(little) swabip(buf+138, 2);
  memcpy(&h->angular_compression, buf+138, 2);
  if(little) swabip(buf+140, 2);
  memcpy(&h->coin_samp_mode, buf+140, 2);
  if(little) swabip(buf+142, 2);
  memcpy(&h->axial_samp_mode, buf+142, 2);
  if(little) swawbip(buf+144, 4);
  memcpy(&h->ecat_calibration_factor, buf+144, 4);
  if(little) swabip(buf+148, 2);
  memcpy(&h->calibration_units, buf+148, 2);
  if(little) swabip(buf+150, 2);
  memcpy(&h->calibration_units_label, buf+150, 2);
  if(little) swabip(buf+152, 2);
  memcpy(&h->compression_code, buf+152, 2);
  memcpy(&h->study_type, buf+154, 12);
  memcpy(&h->patient_id, buf+166, 16);
  memcpy(&h->patient_name, buf+182, 32);
  memcpy(&h->patient_sex, buf+214, 1);
  memcpy(&h->patient_dexterity, buf+215, 1);
  if(little) swawbip(buf+216, 4);
  memcpy(&h->patient_age, buf+216, 4);
  if(little) swawbip(buf+220, 4);
  memcpy(&h->patient_height, buf+220, 4);
  if(little) swawbip(buf+224, 4);
  memcpy(&h->patient_weight, buf+224, 4);
  if(little) swawbip(buf+228, 4);
  memcpy(&h->patient_birth_date, buf+228, 4);
  memcpy(&h->physician_name, buf+232, 32);
  memcpy(&h->operator_name, buf+264, 32);
  memcpy(&h->study_description, buf+296, 32);
  if(little) swabip(buf+328, 2);
  memcpy(&h->acquisition_type, buf+328, 2);
  if(little) swabip(buf+330, 2);
  memcpy(&h->patient_orientation, buf+330, 2);
  memcpy(&h->facility_name, buf+332, 20);
  if(little) swabip(buf+352, 2);
  memcpy(&h->num_planes, buf+352, 2);
  if(little) swabip(buf+354, 2);
  memcpy(&h->num_frames, buf+354, 2);
  if(little) swabip(buf+356, 2);
  memcpy(&h->num_gates, buf+356, 2);
  if(little) swabip(buf+358, 2);
  memcpy(&h->num_bed_pos, buf+358, 2);
  if(little) swawbip(buf+360, 4);
  memcpy(&h->init_bed_position, buf+360, 4);
  if(little) swawbip(buf+364, 15*4);
  memcpy(h->bed_position, buf+364, 15*4);
  if(little) swawbip(buf+424, 4);
  memcpy(&h->plane_separation, buf+424, 4);
  if(little) swabip(buf+428, 2);
  memcpy(&h->lwr_sctr_thres, buf+428, 2);
  if(little) swabip(buf+430, 2);
  memcpy(&h->lwr_true_thres, buf+430, 2);
  memcpy(&h->upr_true_thres, buf+432, 2); if(little) swabip(&h->upr_true_thres,2);
  memcpy(&h->user_process_code, buf+434, 10);
  if(little) swabip(buf+444, 2);
  memcpy(&h->acquisition_mode, buf+444, 2);
  if(little) swawbip(buf+446, 4);
  memcpy(&h->bin_size, buf+446, 4);
  if(little) swawbip(buf+450, 4);
  memcpy(&h->branching_fraction, buf+450, 4);
  if(little) swawbip(buf+454, 4);
  memcpy(&h->dose_start_time, buf+454, 4);
  if(little) swawbip(buf+458, 4);
  memcpy(&h->dosage, buf+458, 4);
  if(little) swawbip(buf+462, 4);
  memcpy(&h->well_counter_corr_factor, buf+462,4);
  memcpy(&h->data_units, buf+466, 32);
  if(little) swabip(buf+498, 2);
  memcpy(&h->septa_state, buf+498, 2);
  memcpy(&h->fill_cti, buf+500, 12);

  /* Patient birth date can have been saved in two different int formats,
     either YYYYMMDD or as seconds from start of year 1970. In latter case
     the number can be negative. */
  /* Seconds from start of year 1970 are converted to YYYYMMDD format */
  if(isdate4(h->patient_birth_date, NULL, NULL, NULL)!=0) {
    time_t t=h->patient_birth_date; gmtime_r(&t, &st);
    h->patient_birth_date=10000*(st.tm_year+1900)+100*(st.tm_mon+1)+st.tm_mday;
  }

  return(0);
}
/*****************************************************************************/

/*****************************************************************************/
/*!
 * Read ECAT 7.x image header
 *
 * @param fp 	input file pointer
 * @param blk block number [1..number of blocks]
 * @param h	Ecat7 image header
 * @return 0 if ok, 1 == invalid parameters, 2 == first header block not found,
 * 3 == header block not read properly
 */
int ecat7ReadImageheader(FILE *fp, int blk, ECAT7_imageheader *h) {
  unsigned char buf[MatBLKSIZE];
  int little; /* 1 if current platform is little endian (i386), else 0 */

  if(ECAT7_TEST) printf("ecat7ReadImageheader()\n");
  if(fp==NULL || h==NULL) return(1);
  little=little_endian(); if(ECAT7_TEST) printf("little=%d\n", little);

  /* Seek the subheader block */
  fseek(fp, (blk-1)*MatBLKSIZE, SEEK_SET);
  if(ftell(fp)!=(blk-1)*MatBLKSIZE) return(2);
  /* Read the header block */
  if(fread(buf, MatBLKSIZE, 1, fp)<1) return(3);

  /* Copy the header fields and swap if necessary */
  if(little) swabip(buf+0, 2);
  memcpy(&h->data_type, buf+0, 2);
  if(little) swabip(buf+2, 2);
  memcpy(&h->num_dimensions, buf+2, 2);
  if(little) swabip(buf+4, 2);
  memcpy(&h->x_dimension, buf+4, 2);
  if(little) swabip(buf+6, 2);
  memcpy(&h->y_dimension, buf+6, 2);
  if(little) swabip(buf+8, 2);
  memcpy(&h->z_dimension, buf+8, 2);
  if(little) swawbip(buf+10, 4);
  memcpy(&h->x_offset, buf+10, 4);
  if(little) swawbip(buf+14, 4);
  memcpy(&h->y_offset, buf+14, 4);
  if(little) swawbip(buf+18, 4);
  memcpy(&h->z_offset, buf+18, 4);
  if(little) swawbip(buf+22, 4);
  memcpy(&h->recon_zoom, buf+22, 4);
  if(little) swawbip(buf+26, 4);
  memcpy(&h->scale_factor, buf+26, 4);
  if(little) swabip(buf+30, 2);
  memcpy(&h->image_min, buf+30, 2);
  if(little) swabip(buf+32, 2);
  memcpy(&h->image_max, buf+32, 2);
  if(little) swawbip(buf+34, 4);
  memcpy(&h->x_pixel_size, buf+34, 4);
  if(little) swawbip(buf+38, 4);
  memcpy(&h->y_pixel_size, buf+38, 4);
  if(little) swawbip(buf+42, 4);
  memcpy(&h->z_pixel_size, buf+42, 4);
  if(little) swawbip(buf+46, 4);
  memcpy(&h->frame_duration, buf+46, 4);
  if(little) swawbip(buf+50, 4);
  memcpy(&h->frame_start_time, buf+50, 4);
  if(little) swabip(buf+54, 2);
  memcpy(&h->filter_code, buf+54, 2);
  if(little) swawbip(buf+56, 4);
  memcpy(&h->x_resolution, buf+56, 4);
  if(little) swawbip(buf+60, 4);
  memcpy(&h->y_resolution, buf+60, 4);
  if(little) swawbip(buf+64, 4);
  memcpy(&h->z_resolution, buf+64, 4);
  if(little) swawbip(buf+68, 4);
  memcpy(&h->num_r_elements, buf+68, 4);
  if(little) swawbip(buf+72, 4);
  memcpy(&h->num_angles, buf+72, 4);
  if(little) swawbip(buf+76, 4);
  memcpy(&h->z_rotation_angle, buf+76, 4);
  if(little) swawbip(buf+80, 4);
  memcpy(&h->decay_corr_fctr, buf+80, 4);
  if(little) swawbip(buf+84, 4);
  memcpy(&h->processing_code, buf+84, 4);
  if(little) swawbip(buf+88, 4);
  memcpy(&h->gate_duration, buf+88, 4);
  if(little) swawbip(buf+92, 4);
  memcpy(&h->r_wave_offset, buf+92, 4);
  if(little) swawbip(buf+96, 4);
  memcpy(&h->num_accepted_beats, buf+96, 4);
  if(little) swawbip(buf+100, 4);
  memcpy(&h->filter_cutoff_frequency, buf+100, 4);
  if(little) swawbip(buf+104, 4);
  memcpy(&h->filter_resolution, buf+104, 4);
  if(little) swawbip(buf+108, 4);
  memcpy(&h->filter_ramp_slope, buf+108, 4);
  if(little) swabip(buf+112, 2);
  memcpy(&h->filter_order, buf+112, 2);
  if(little) swawbip(buf+114, 4);
  memcpy(&h->filter_scatter_fraction, buf+114, 4);
  if(little) swawbip(buf+118, 4);
  memcpy(&h->filter_scatter_slope, buf+118, 4);
  memcpy(&h->annotation, buf+122, 40);
  if(little) swawbip(buf+162, 4);
  memcpy(&h->mt_1_1, buf+162, 4);
  if(little) swawbip(buf+166, 4);
  memcpy(&h->mt_1_2, buf+166, 4);
  if(little) swawbip(buf+170, 4);
  memcpy(&h->mt_1_3, buf+170, 4);
  if(little) swawbip(buf+174, 4);
  memcpy(&h->mt_2_1, buf+174, 4);
  if(little) swawbip(buf+178, 4);
  memcpy(&h->mt_2_2, buf+178, 4);
  if(little) swawbip(buf+182, 4);
  memcpy(&h->mt_2_3, buf+182, 4);
  if(little) swawbip(buf+186, 4);
  memcpy(&h->mt_3_1, buf+186, 4);
  if(little) swawbip(buf+190, 4);
  memcpy(&h->mt_3_2, buf+190, 4);
  if(little) swawbip(buf+194, 4);
  memcpy(&h->mt_3_3, buf+194, 4);
  if(little) swawbip(buf+198, 4);
  memcpy(&h->rfilter_cutoff, buf+198, 4);
  if(little) swawbip(buf+202, 4);
  memcpy(&h->rfilter_resolution, buf+202, 4);
  if(little) swabip(buf+206, 2);
  memcpy(&h->rfilter_code, buf+206, 2);
  if(little) swabip(buf+208, 2);
  memcpy(&h->rfilter_order, buf+208, 2);
  if(little) swawbip(buf+210, 4);
  memcpy(&h->zfilter_cutoff, buf+210, 4);
  if(little) swawbip(buf+214, 4);
  memcpy(&h->zfilter_resolution, buf+214, 4);
  if(little) swabip(buf+218, 2);
  memcpy(&h->zfilter_code, buf+218, 2);
  if(little) swabip(buf+220, 2);
  memcpy(&h->zfilter_order, buf+220, 2);
  if(little) swawbip(buf+222, 4);
  memcpy(&h->mt_1_4, buf+222, 4);
  if(little) swawbip(buf+226, 4);
  memcpy(&h->mt_2_4, buf+226, 4);
  if(little) swawbip(buf+230, 4);
  memcpy(&h->mt_3_4, buf+230, 4);
  if(little) swabip(buf+234, 2);
  memcpy(&h->scatter_type, buf+234, 2);
  if(little) swabip(buf+236, 2);
  memcpy(&h->recon_type, buf+236, 2);
  if(little) swabip(buf+238, 2);
  memcpy(&h->recon_views, buf+238, 2);
  memcpy(&h->fill_cti, buf+240, 174);
  memcpy(&h->fill_user, buf+414, 96);

  return(0);
}
/*****************************************************************************/

/*****************************************************************************/
/*!
 * Read ECAT 7.x attenuation header
 *
 * @param fp 	input file pointer
 * @param blk 	block number [1..number of blocks]
 * @param h	Ecat7 attenuation header
 * @return 0 if ok, 1 == invalid parameters, 2 == first header block not found,
 * 3 == header block not read properly
 */
int ecat7ReadAttenheader(FILE *fp, int blk, ECAT7_attenheader *h) {
  unsigned char buf[MatBLKSIZE];
  int little; /* 1 if current platform is little endian (i386), else 0 */

  if(ECAT7_TEST) printf("ecat7ReadAttenheader()\n");
  if(fp==NULL || h==NULL) return(1);
  little=little_endian();

  /* Seek the subheader block */
  fseek(fp, (blk-1)*MatBLKSIZE, SEEK_SET);
  if(ftell(fp)!=(blk-1)*MatBLKSIZE) return(2);
  /* Read the header block */
  if(fread(buf, MatBLKSIZE, 1, fp)<1) return(3);
  /* Copy the header fields and swap if necessary */
  if(little) swabip(buf+0, 2);
  memcpy(&h->data_type, buf+0, 2);
  if(little) swabip(buf+2, 2);
  memcpy(&h->num_dimensions, buf+2, 2);
  if(little) swabip(buf+4, 2);
  memcpy(&h->attenuation_type, buf+4, 2);
  if(little) swabip(buf+6, 2);
  memcpy(&h->num_r_elements, buf+6, 2);
  if(little) swabip(buf+8, 2);
  memcpy(&h->num_angles, buf+8, 2);
  if(little) swabip(buf+10, 2);
  memcpy(&h->num_z_elements, buf+10, 2);
  if(little) swabip(buf+12, 2);
  memcpy(&h->ring_difference, buf+12, 2);
  if(little) swawbip(buf+14, 4);
  memcpy(&h->x_resolution, buf+14, 4);
  if(little) swawbip(buf+18, 4);
  memcpy(&h->y_resolution, buf+18, 4);
  if(little) swawbip(buf+22, 4);
  memcpy(&h->z_resolution, buf+22, 4);
  if(little) swawbip(buf+26, 4);
  memcpy(&h->w_resolution, buf+26, 4);
  if(little) swawbip(buf+30, 4);
  memcpy(&h->scale_factor, buf+30, 4);
  if(little) swawbip(buf+34, 4);
  memcpy(&h->x_offset, buf+34, 4);
  if(little) swawbip(buf+38, 4);
  memcpy(&h->y_offset, buf+38, 4);
  if(little) swawbip(buf+42, 4);
  memcpy(&h->x_radius, buf+42, 4);
  if(little) swawbip(buf+46, 4);
  memcpy(&h->y_radius, buf+46, 4);
  if(little) swawbip(buf+50, 4);
  memcpy(&h->tilt_angle, buf+50, 4);
  if(little) swawbip(buf+54, 4);
  memcpy(&h->attenuation_coeff, buf+54, 4);
  if(little) swawbip(buf+58, 4);
  memcpy(&h->attenuation_min, buf+58, 4);
  if(little) swawbip(buf+62, 4);
  memcpy(&h->attenuation_max, buf+62, 4);
  if(little) swawbip(buf+66, 4);
  memcpy(&h->skull_thickness, buf+66, 4);
  if(little) swabip(buf+70, 2);
  memcpy(&h->num_additional_atten_coeff, buf+70, 2);
  if(little) swawbip(buf+72, 8*4);
  memcpy(h->additional_atten_coeff, buf+72, 8*4);
  if(little) swawbip(buf+104, 4);
  memcpy(&h->edge_finding_threshold, buf+104, 4);
  if(little) swabip(buf+108, 2);
  memcpy(&h->storage_order, buf+108, 2);
  if(little) swabip(buf+110, 2);
  memcpy(&h->span, buf+110, 2);
  if(little) swabip(buf+112, 64*2);
  memcpy(h->z_elements, buf+112, 64*2);
  if(little) swabip(buf+240, 86*2);
  memcpy(h->fill_cti, buf+240, 86*2);
  if(little) swabip(buf+412, 50*2);
  memcpy(h->fill_user, buf+412, 50*2);
  return(0);
}
/*****************************************************************************/

/*****************************************************************************/
/*!
 * Read ECAT 7.x polar map header
 *
 * @param fp 	input file pointer
 * @param blk 	block number [1..number of blocks]
 * @param h	Ecat7 polar map header
 * @return 0 if ok, 1 == invalid parameters, 2 == first header block not found,
 * 3 == header block not read properly
 */
int ecat7ReadPolmapheader(FILE *fp, int blk, ECAT7_polmapheader *h) {
  unsigned char buf[MatBLKSIZE];
  int little; /* 1 if current platform is little endian (i386), else 0 */

  if(ECAT7_TEST) printf("ecat7ReadPolarmapheader()\n");
  if(fp==NULL || h==NULL) return(1);
  little=little_endian();

  /* Seek the subheader block */
  fseek(fp, (blk-1)*MatBLKSIZE, SEEK_SET);
  if(ftell(fp)!=(blk-1)*MatBLKSIZE) return(2);
  /* Read the header block */
  if(fread(buf, MatBLKSIZE, 1, fp)<1) return(3);
  /* Copy the header fields and swap if necessary */
  if(little) swabip(buf+0, 2);
  memcpy(&h->data_type, buf+0, 2);
  if(little) swabip(buf+2, 2);
  memcpy(&h->polar_map_type, buf+2, 2);
  if(little) swabip(buf+4, 2);
  memcpy(&h->num_rings, buf+4, 2);
  if(little) swabip(buf+6, 32*2);
  memcpy(h->sectors_per_ring, buf+6, 32*2);
  if(little) swawbip(buf+70, 32*4);
  memcpy(h->ring_position, buf+70, 32*4);
  if(little) swabip(buf+198, 32*2);
  memcpy(h->ring_angle, buf+198, 32*2);
  if(little) swabip(buf+262, 2);
  memcpy(&h->start_angle, buf+262, 2);
  if(little) swabip(buf+264, 3*2);
  memcpy(h->long_axis_left, buf+264, 3*2);
  if(little) swabip(buf+270, 3*2);
  memcpy(h->long_axis_right, buf+270, 3*2);
  if(little) swabip(buf+276, 2);
  memcpy(&h->position_data, buf+276, 2);
  if(little) swabip(buf+278, 2);
  memcpy(&h->image_min, buf+278, 2);
  if(little) swabip(buf+280, 2);
  memcpy(&h->image_max, buf+280, 2);
  if(little) swawbip(buf+282, 4);
  memcpy(&h->scale_factor, buf+282, 4);
  if(little) swawbip(buf+286, 4);
  memcpy(&h->pixel_size, buf+286, 4);
  if(little) swawbip(buf+290, 4);
  memcpy(&h->frame_duration, buf+290, 4);
  if(little) swawbip(buf+294, 4);
  memcpy(&h->frame_start_time, buf+294, 4);
  if(little) swabip(buf+298, 2);
  memcpy(&h->processing_code, buf+298, 2);
  if(little) swabip(buf+300, 2);
  memcpy(&h->quant_units, buf+300, 2);
  memcpy(h->annotation, buf+302, 40);
  if(little) swawbip(buf+342, 4);
  memcpy(&h->gate_duration, buf+342, 4);
  if(little) swawbip(buf+346, 4);
  memcpy(&h->r_wave_offset, buf+346, 4);
  if(little) swawbip(buf+350, 4);
  memcpy(&h->num_accepted_beats, buf+350, 4);
  memcpy(h->polar_map_protocol, buf+354, 20);
  memcpy(h->database_name, buf+374, 30);
  if(little) swabip(buf+404, 27*2);
  memcpy(h->fill_cti, buf+404, 27*2);
  return(0);
}
/*****************************************************************************/

/*****************************************************************************/
/*!
 * Read ECAT 7.x 3D normalization header
 *
 * @param fp 	input file pointer
 * @param blk 	block number [1..number of blocks] 
 * @param h	Ecat7 normalization header
 * @return 0 if ok, 1 == invalid parameters, 2 == first header block not found,
 * 3 == header block not read properly
 */
int ecat7ReadNormheader(FILE *fp, int blk, ECAT7_normheader *h) {
  unsigned char buf[MatBLKSIZE];
  int little; /* 1 if current platform is little endian (i386), else 0 */

  if(ECAT7_TEST) printf("ecat7ReadNormheader()\n");
  if(fp==NULL || h==NULL) return(1);
  little=little_endian();

  /* Seek the subheader block */
  fseek(fp, (blk-1)*MatBLKSIZE, SEEK_SET);
  if(ftell(fp)!=(blk-1)*MatBLKSIZE) return(2); 
  /* Read the header block */
  if(fread(buf, MatBLKSIZE, 1, fp)<1) return(3);
  /* Copy the header fields and swap if necessary */
  if(little) swabip(buf+0, 2);
  memcpy(&h->data_type, buf+0, 2);
  if(little) swabip(buf+2, 2);
  memcpy(&h->num_r_elements, buf+2, 2);
  if(little) swabip(buf+4, 2);
  memcpy(&h->num_transaxial_crystals, buf+4, 2);
  if(little) swabip(buf+6, 2);
  memcpy(&h->num_crystal_rings, buf+6, 2);
  if(little) swabip(buf+8, 2);
  memcpy(&h->crystals_per_ring, buf+8, 2);
  if(little) swabip(buf+10, 2);
  memcpy(&h->num_geo_corr_planes, buf+10, 2);
  if(little) swabip(buf+12, 2);
  memcpy(&h->uld, buf+12, 2);
  if(little) swabip(buf+14, 2);
  memcpy(&h->lld, buf+14, 2);
  if(little) swabip(buf+16, 2);
  memcpy(&h->scatter_energy, buf+16, 2);
  if(little) swawbip(buf+18, 4);
  memcpy(&h->norm_quality_factor, buf+18, 4);
  if(little) swabip(buf+22, 2);
  memcpy(&h->norm_quality_factor_code, buf+22, 2);
  if(little) swawbip(buf+24, 32*4);
  memcpy(h->ring_dtcor1, buf+24, 32*4);
  if(little) swawbip(buf+152, 32*4);
  memcpy(h->ring_dtcor2, buf+152, 32*4);
  if(little) swawbip(buf+280, 8*4);
  memcpy(h->crystal_dtcor, buf+280, 8*4);
  if(little) swabip(buf+312, 2);
  memcpy(&h->span, buf+312, 2);
  if(little) swabip(buf+314, 2);
  memcpy(&h->max_ring_diff, buf+314, 2);
  if(little) swabip(buf+316, 48*2);
  memcpy(h->fill_cti, buf+316, 48*2);
  if(little) swabip(buf+412, 50*2);
  memcpy(h->fill_user, buf+412, 50*2);
  return(0);
}
/*****************************************************************************/

/*****************************************************************************/
/*!
 * Read ECAT 7.x 3D scan header (512 bytes)
 *
 * @param fp 	input file pointer
 * @param blk 	block number [1..number of blocks]
 * @param h	Ecat7 scan header
 * @return 0 if ok, 1 == invalid parameters, 2 == first header block not found,
 * 3 == header block not read properly
 */
int ecat7ReadScanheader(FILE *fp, int blk, ECAT7_scanheader *h) {
  unsigned char buf[2*MatBLKSIZE];
  int little; /* 1 if current platform is little endian (i386), else 0 */

  if(ECAT7_TEST) printf("ecat7ReadScanheader()\n");
  if(fp==NULL || h==NULL) return(1);
  little=little_endian(); if(ECAT7_TEST) printf("little=%d\n", little);

  /* Seek the subheader block */
  fseek(fp, (blk-1)*MatBLKSIZE, SEEK_SET);
  if(ftell(fp)!=(blk-1)*MatBLKSIZE) return(2);
  /* Read the header block */
  if(fread(buf, MatBLKSIZE, 2, fp)<1) return(3);

  /* Copy the header fields and swap if necessary */
  if(little) swabip(buf+0, 2);
  memcpy(&h->data_type, buf+0, 2);
  if(little) swabip(buf+2, 2);
  memcpy(&h->num_dimensions, buf+2, 2);
  if(little) swabip(buf+4, 2);
  memcpy(&h->num_r_elements, buf+4, 2);
  if(little) swabip(buf+6, 2);
  memcpy(&h->num_angles, buf+6, 2);
  if(little) swabip(buf+8, 2);
  memcpy(&h->corrections_applied, buf+8, 2);
  if(little) swabip(buf+10, 64*2);
  memcpy(h->num_z_elements, buf+10, 64*2);
  if(little) swabip(buf+138, 2);
  memcpy(&h->ring_difference, buf+138, 2);
  if(little) swabip(buf+140, 2);
  memcpy(&h->storage_order, buf+140, 2);
  if(little) swabip(buf+142, 2);
  memcpy(&h->axial_compression, buf+142, 2);
  if(little) swawbip(buf+144, 4);
  memcpy(&h->x_resolution, buf+144, 4);
  if(little) swawbip(buf+148, 4);
  memcpy(&h->v_resolution, buf+148, 4);
  if(little) swawbip(buf+152, 4);
  memcpy(&h->z_resolution, buf+152, 4);
  if(little) swawbip(buf+156, 4);
  memcpy(&h->w_resolution, buf+156, 4);
  if(little) swabip(buf+160, 6*2);
  memcpy(h->fill_gate, buf+160, 6*2);
  if(little) swawbip(buf+172, 4);
  memcpy(&h->gate_duration, buf+172, 4);
  if(little) swawbip(buf+176, 4);
  memcpy(&h->r_wave_offset, buf+176, 4);
  if(little) swawbip(buf+180, 4);
  memcpy(&h->num_accepted_beats, buf+180, 4);
  if(little) swawbip(buf+184, 4);
  memcpy(&h->scale_factor, buf+184, 4);
  if(little) swabip(buf+188, 2);
  memcpy(&h->scan_min, buf+188, 2);
  if(little) swabip(buf+190, 2);
  memcpy(&h->scan_max, buf+190, 2);
  if(little) swawbip(buf+192, 4);
  memcpy(&h->prompts, buf+192, 4);
  if(little) swawbip(buf+196, 4);
  memcpy(&h->delayed, buf+196, 4);
  if(little) swawbip(buf+200, 4);
  memcpy(&h->multiples, buf+200, 4);
  if(little) swawbip(buf+204, 4);
  memcpy(&h->net_trues, buf+204, 4);
  if(little) swawbip(buf+208, 4);
  memcpy(&h->tot_avg_cor, buf+208, 4);
  if(little) swawbip(buf+212, 4);
  memcpy(&h->tot_avg_uncor, buf+212, 4);
  if(little) swawbip(buf+216, 4);
  memcpy(&h->total_coin_rate, buf+216, 4);
  if(little) swawbip(buf+220, 4);
  memcpy(&h->frame_start_time, buf+220, 4);
  if(little) swawbip(buf+224, 4);
  memcpy(&h->frame_duration, buf+224, 4);
  if(little) swawbip(buf+228, 4);
  memcpy(&h->deadtime_correction_factor, buf+228,4);
  if(little) swabip(buf+232, 90*2);
  memcpy(h->fill_cti, buf+232, 90*2);
  if(little) swabip(buf+412, 50*2);
  memcpy(h->fill_user, buf+412, 50*2);
  if(little) swawbip(buf+512, 128*4);
  memcpy(h->uncor_singles, buf+512, 128*4);
  return(0);
}
/*****************************************************************************/

/*****************************************************************************/
/*!
 * Read ECAT 7.x 2D scan header
 *
 * @param fp 	input file pointer
 * @param blk 	block number [1..number of blocks]
 * @param h	Ecat7 2D scan header
 * @return 0 if ok, 1 == invalid parameters, 2 == first header block not found,
 * 3 == header block not read properly
 */
int ecat7Read2DScanheader(FILE *fp, int blk, ECAT7_2Dscanheader *h) {
  unsigned char buf[MatBLKSIZE];
  int little; /* 1 if current platform is little endian (i386), else 0 */

  if(ECAT7_TEST) printf("ecat7Read2DScanheader()\n");
  if(fp==NULL || h==NULL) return(1);
  little=little_endian();

  /* Seek the subheader block */
  fseek(fp, (blk-1)*MatBLKSIZE, SEEK_SET);
  if(ftell(fp)!=(blk-1)*MatBLKSIZE) return(2);
  /* Read the header block */
  if(fread(buf, MatBLKSIZE, 1, fp)<1) return(3);
  /* Copy the header fields and swap if necessary */
  if(little) swabip(buf+0, 2);
  memcpy(&h->data_type, buf+0, 2);
  if(little) swabip(buf+2, 2);
  memcpy(&h->num_dimensions, buf+2, 2);
  if(little) swabip(buf+4, 2);
  memcpy(&h->num_r_elements, buf+4, 2);
  if(little) swabip(buf+6, 2);
  memcpy(&h->num_angles, buf+6, 2);
  if(little) swabip(buf+8, 2);
  memcpy(&h->corrections_applied, buf+8, 2);
  if(little) swabip(buf+10, 2);
  memcpy(&h->num_z_elements, buf+10, 2);
  if(little) swabip(buf+12, 2);
  memcpy(&h->ring_difference, buf+12, 2);
  if(little) swawbip(buf+14, 4);
  memcpy(&h->x_resolution, buf+14, 4);
  if(little) swawbip(buf+18, 4);
  memcpy(&h->y_resolution, buf+18, 4);
  if(little) swawbip(buf+22, 4);
  memcpy(&h->z_resolution, buf+22, 4);
  if(little) swawbip(buf+26, 4);
  memcpy(&h->w_resolution, buf+26, 4);
  if(little) swabip(buf+30, 6*2);
  memcpy(h->fill_gate, buf+30, 6*2);
  if(little) swawbip(buf+42, 4);
  memcpy(&h->gate_duration, buf+42, 4);
  if(little) swawbip(buf+46, 4);
  memcpy(&h->r_wave_offset, buf+46, 4);
  if(little) swawbip(buf+50, 4);
  memcpy(&h->num_accepted_beats, buf+50, 4);
  if(little) swawbip(buf+54, 4);
  memcpy(&h->scale_factor, buf+54, 4);
  if(little) swabip(buf+58, 2);
  memcpy(&h->scan_min, buf+58, 2);
  if(little) swabip(buf+60, 2);
  memcpy(&h->scan_max, buf+60, 2);
  if(little) swawbip(buf+62, 4);
  memcpy(&h->prompts, buf+62, 4);
  if(little) swawbip(buf+66, 4);
  memcpy(&h->delayed, buf+66, 4);
  if(little) swawbip(buf+70, 4);
  memcpy(&h->multiples, buf+70, 4);
  if(little) swawbip(buf+74, 4);
  memcpy(&h->net_trues, buf+74, 4);
  if(little) swawbip(buf+78, 16*4);
  memcpy(h->cor_singles, buf+78, 16*4);
  if(little) swawbip(buf+142, 16*4);
  memcpy(h->uncor_singles, buf+142, 16*4);
  if(little) swawbip(buf+206, 4);
  memcpy(&h->tot_avg_cor, buf+206, 4);
  if(little) swawbip(buf+210, 4);
  memcpy(&h->tot_avg_uncor, buf+210, 4);
  if(little) swawbip(buf+214, 4);
  memcpy(&h->total_coin_rate, buf+214, 4);
  if(little) swawbip(buf+218, 4);
  memcpy(&h->frame_start_time, buf+218, 4);
  if(little) swawbip(buf+222, 4);
  memcpy(&h->frame_duration, buf+222, 4);
  if(little) swawbip(buf+226, 4);
  memcpy(&h->deadtime_correction_factor, buf+226,4);
  if(little) swabip(buf+230, 8*2);
  memcpy(h->physical_planes, buf+230, 8*2);
  if(little) swabip(buf+246, 83*2);
  memcpy(h->fill_cti, buf+246, 83*2);
  if(little) swabip(buf+412, 50*2);
  memcpy(h->fill_user, buf+412, 50*2);
  return(0);
}
/*****************************************************************************/

/*****************************************************************************/
/*!
 * Read ECAT 7.x 2D normalization header
 *
 * @param fp 	input file pointer
 * @param blk 	block number [1..number of blocks]
 * @param h	Ecat7 normalization header
 * @return 0 if ok, 1 == invalid parameters, 2 == first header block not found,
 * 3 == header block not read properly
 */
int ecat7Read2DNormheader(FILE *fp, int blk, ECAT7_2Dnormheader *h) {
  unsigned char buf[MatBLKSIZE];
  int little; /* 1 if current platform is little endian (i386), else 0 */

  if(ECAT7_TEST) printf("ecat7Read2Dnormheader()\n");
  if(fp==NULL || h==NULL) return(1);
  little=little_endian();

  /* Seek the subheader block */
  fseek(fp, (blk-1)*MatBLKSIZE, SEEK_SET);
  if(ftell(fp)!=(blk-1)*MatBLKSIZE) return(2); 
  /* Read the header block */
  if(fread(buf, MatBLKSIZE, 1, fp)<1) return(3);
  /* Copy the header fields and swap if necessary */
  if(little) swabip(buf+0, 2);
  memcpy(&h->data_type, buf+0, 2);
  if(little) swabip(buf+2, 2);
  memcpy(&h->num_dimensions, buf+2, 2);
  if(little) swabip(buf+4, 2);
  memcpy(&h->num_r_elements, buf+4, 2);
  if(little) swabip(buf+6, 2);
  memcpy(&h->num_angles, buf+6, 2);
  if(little) swabip(buf+8, 2);
  memcpy(&h->num_z_elements, buf+8, 2);
  if(little) swabip(buf+10, 2);
  memcpy(&h->ring_difference, buf+10, 2);
  if(little) swawbip(buf+12, 4);
  memcpy(&h->scale_factor, buf+12, 4);
  if(little) swawbip(buf+16, 4);
  memcpy(&h->norm_min, buf+16, 4);
  if(little) swawbip(buf+20, 4);
  memcpy(&h->norm_max, buf+20, 4);
  if(little) swawbip(buf+24, 4);
  memcpy(&h->fov_source_width, buf+24, 4);
  if(little) swawbip(buf+28, 4);
  memcpy(&h->norm_quality_factor, buf+28, 4);
  if(little) swabip(buf+32, 2);
  memcpy(&h->norm_quality_factor_code, buf+32, 2);
  if(little) swabip(buf+34, 2);
  memcpy(&h->storage_order, buf+34, 2);
  if(little) swabip(buf+36, 2);
  memcpy(&h->span, buf+36, 2);
  if(little) swabip(buf+38, 64*2);
  memcpy(h->fill_cti, buf+38, 64*2);
  if(little) swabip(buf+166, 123*2);
  memcpy(h->fill_cti, buf+166, 123*2);
  if(little) swabip(buf+412, 50*2);
  memcpy(h->fill_user, buf+412, 50*2);
  return(0);
}
/*****************************************************************************/

/*****************************************************************************/
/*!
 * Read ECAT7 matrix data and convert byte order if necessary
 * Remember to allocate memory for full blocks!
 * There are differences here when compared to ecat63.c
 *
 * @param fp 		input file pointer
 * @param start_block 	starting block index
 * @param block_nr 	number of blocks to be read
 * @param data		target buffer 
 * @param dtype		data type of target buffer
 * @return 0 if ok, 1 == invalid parameters, 9 == start block not found,
 * 2 == data blocks read properly
 */
int ecat7ReadMatrixdata(
  FILE *fp, int start_block, int block_nr, char *data, int dtype
) {
  int i, n, little, err=0;
  char *cptr;
  float f;

  if(ECAT7_TEST) printf("ecat7ReadMatrixdata(fp, %d, %d, data, %d)\n",
    start_block, block_nr, dtype);
  /* Check the arguments */
  if(block_nr<=0 || start_block<1 || data==NULL) return(1);
  /* Seek the first data block */
  fseek(fp, (start_block-1)*MatBLKSIZE, SEEK_SET);
  if(ftell(fp)!=(start_block-1)*MatBLKSIZE) return(9);
  /* Read the data blocks */
  if(fread(data, MatBLKSIZE, block_nr, fp) < (unsigned int)block_nr) return(2);
  /* Translate data if necessary */
  little=little_endian();
  switch(dtype) {
    case ECAT7_BYTE: /* byte format...no translation necessary */
      break;
    case ECAT7_VAXI2:  /* byte conversion necessary on big endian platform */
      if(!little) {cptr=data; swabip(cptr, block_nr*MatBLKSIZE);}
      break;
    case ECAT7_VAXI4:
      for(i=0, cptr=data; i<block_nr*MatBLKSIZE; i+=4, cptr+=4) {
        n=ecat7rInt(cptr, 1, little); memcpy(cptr, &n, 4);
      }
      break;
    case ECAT7_VAXR4:
      for(i=0, cptr=data; i<block_nr*MatBLKSIZE; i+=4, cptr+=4) {
        f=ecat7rFloat(cptr, 1, little); memcpy(cptr, &f, 4);
      }
      break;
    case ECAT7_IEEER4: /* IEEE float ; byte conversion necessary on
                          little endian platforms */
    case ECAT7_SUNI4:  /* SUN int ; byte conversion necessary on
                          little endian platforms */
      if(little) swawbip(data, block_nr*MatBLKSIZE);
      break;
    case ECAT7_SUNI2:  /* SUN short ; byte conversion necessary on
                          little endian platforms */
      if(little) swabip(data, block_nr*MatBLKSIZE);
      break;
    default:  /* if something else, for now think it as an error */ 
      err=2;
      break;
  }
  return(err);
}
/*****************************************************************************/

/*****************************************************************************/
/*!
 * Read ECAT7 image matrix header and data. If only header is to be read,
 * set last_block=first_block.
 * Note: data is not calibrated with factor in main header.
 *
 * @param fp ECAT file pointer
 * @param first_block Subheader record number
 * @param last_block Last data block number
 * @param h Ptr to subheader data which is filled
 * @param fdata Ptr to the address of the matrix data
 * @return 0 if ok, 1 invalid input, 5 failed to read subheader,
 * 6 invalid image (x,y,z) dimensions, 8 failed to allocate memory for meta-data,
 * 9 failed to read matrix data, 11 failed to allocate memory for voxel data
 */
int ecat7ReadImageMatrix(
  FILE *fp, int first_block, int last_block, ECAT7_imageheader *h, float **fdata
) {
  int i, ret, blockNr, pxlNr;
  char *mdata, *mptr;
  float *_fdata, *fptr;
  short int *sptr;
  int *iptr;
  
  
  if(ECAT7_TEST) printf("ecat7ReadImageMatrix(fp, %d, %d, hdr, fdata)\n",
    first_block, last_block);
  if(fp==NULL || first_block<=MatFirstDirBlk || h==NULL) {
    sprintf(ecat7errmsg, "invalid function parameter.\n");
    return(1);
  }
  *fdata=(float*)NULL;
  
  /* Read subheader */
  ret=ecat7ReadImageheader(fp, first_block, h);
  if(ret) {
    sprintf(ecat7errmsg, "cannot read subheader (%d).\n", ret);
    return(5);
  }
  if(ECAT7_TEST>4) ecat7PrintImageheader(h, stdout);
  pxlNr=h->x_dimension*h->y_dimension;
  if(h->num_dimensions>2) pxlNr*=h->z_dimension;
  if(pxlNr<=0) {
    sprintf(ecat7errmsg, "invalid matrix dimension.\n");
    return(6);  
  }
  
  /* Read matrix data */
  blockNr=last_block-first_block; if(blockNr<1) return(0);
  mdata=(char*)malloc(blockNr*MatBLKSIZE);
  if(mdata==NULL) {
    sprintf(ecat7errmsg, "cannot allocate memory.\n");
    return(8);  
  }
  mptr=mdata;
  ret=ecat7ReadMatrixdata(fp, first_block+1, blockNr, mptr, h->data_type);
  if(ret || mdata==NULL) {
    sprintf(ecat7errmsg, "cannot read matrix data (%d).\n", ret);
    free(mdata); return(9);
  }
  
  /* Allocate memory for float data */
  _fdata=(float*)malloc(pxlNr*sizeof(float));
  if(_fdata==NULL) {
    sprintf(ecat7errmsg, "cannot allocate memory.\n");
    free(mdata); return(11);  
  }

  /* Convert matrix data to floats */
  fptr=_fdata; mptr=mdata;
  if(h->data_type==ECAT7_BYTE) {
    for(i=0; i<pxlNr; i++, mptr++, fptr++)
      *fptr=h->scale_factor*(float)(*mptr);
  } else if(h->data_type==ECAT7_VAXI2 || h->data_type==ECAT7_SUNI2) {
    for(i=0; i<pxlNr; i++, mptr+=2, fptr++) {
      sptr=(short int*)mptr;
      *fptr=h->scale_factor*(float)(*sptr);
      if(!(*fptr>-1.0E+22 && *fptr<1.0E+22)) *fptr=0.0;
    }
  } else if(h->data_type==ECAT7_VAXI4 || h->data_type==ECAT7_SUNI4) {
    for(i=0; i<pxlNr; i++, mptr+=4, fptr++) {
      iptr=(int*)mptr;
      *fptr=h->scale_factor*(float)(*iptr);
      if(!(*fptr>-1.0E+22 && *fptr<1.0E+22)) *fptr=0.0;
    }
  } else if(h->data_type==ECAT7_VAXR4 || h->data_type==ECAT7_IEEER4) {
    memcpy(fptr, mptr, pxlNr*4);
    for(i=0; i<pxlNr; i++, fptr++) {
      *fptr *= h->scale_factor;
      if(!(*fptr>-1.0E+22 && *fptr<1.0E+22)) *fptr=0.0;
    }
  }
  free(mdata);
  *fdata=_fdata;

  return(0);
}
/*****************************************************************************/

/*****************************************************************************/
/*!
 * Read ECAT7 2D sinogram matrix header and data
 *   Memory for fdata[] is allocated here, remember to free memory after usage.
 *   Note: data is not calibrated with factor in main header.
 *   Note: data is not multiplied with deadtime_correction_factor.
 *
 * @param fp ECAT file pointer
 * @param first_block Subheader record number
 * @param last_block Last data block number
 * @param h Ptr to subheader data which is filled
 * @param fdata Ptr to the address of the matrix data
 * @return 0 if ok, 1 invalid input, 5 failed to read scan header,
 * 6 invalid image (x,y,z) dimensions, 8 failed to allocate memory for meta-data,
 * 9 failed to read matrix data, 11 failed to allocate memory for voxel data
 */
int ecat7Read2DScanMatrix(FILE *fp, int first_block, int last_block,
			  ECAT7_2Dscanheader *h, float **fdata) {
  int i, ret, blockNr, pxlNr;
  char *mdata, *mptr;
  float *_fdata, *fptr;
  short int *sptr;
  int *iptr;


  if(ECAT7_TEST) printf("ecat7Read2DScanMatrix(fp, %d, %d, hdr, fdata)\n",
    first_block, last_block);
  if(fp==NULL || first_block<=MatFirstDirBlk || h==NULL) {
    sprintf(ecat7errmsg, "invalid function parameter.\n");
    return(1);
  }
  *fdata=(float*)NULL;
  
  /* Read subheader */
  ret=ecat7Read2DScanheader(fp, first_block, h);
  if(ret) {
    sprintf(ecat7errmsg, "cannot read subheader (%d).\n", ret);
    return(5);
  }
  if(ECAT7_TEST>4) ecat7Print2DScanheader(h, stdout);
  pxlNr=h->num_r_elements*h->num_angles;
  if(h->num_dimensions>2) pxlNr*=h->num_z_elements;
  if(pxlNr<=0) {
    sprintf(ecat7errmsg, "invalid matrix dimension.\n");
    return(6);  
  }
  
  /* Read matrix data */
  blockNr=last_block-first_block; if(blockNr<1) return(0);
  mdata=(char*)malloc(blockNr*MatBLKSIZE);
  if(mdata==NULL) {
    sprintf(ecat7errmsg, "cannot allocate memory.\n");
    return(8);  
  }
  mptr=mdata;
  ret=ecat7ReadMatrixdata(fp, first_block+1, blockNr, mptr, h->data_type);
  if(ret || mdata==NULL) {
    sprintf(ecat7errmsg, "cannot read matrix data (%d).\n", ret);
    free(mdata); return(9);
  }
  
  /* Allocate memory for float data */
  _fdata=(float*)malloc(pxlNr*sizeof(float));
  if(_fdata==NULL) {
    sprintf(ecat7errmsg, "cannot allocate memory.\n");
    free(mdata); return(11);  
  }

  /* Convert matrix data to floats */
  fptr=_fdata; mptr=mdata;
  if(h->data_type==ECAT7_BYTE) {
    for(i=0; i<pxlNr; i++, mptr++, fptr++)
      *fptr=h->scale_factor*(float)(*mptr);
  } else if(h->data_type==ECAT7_VAXI2 || h->data_type==ECAT7_SUNI2) {
    for(i=0; i<pxlNr; i++, mptr+=2, fptr++) {
      sptr=(short int*)mptr;
      *fptr=h->scale_factor*(float)(*sptr);
    }
  } else if(h->data_type==ECAT7_VAXI4 || h->data_type==ECAT7_SUNI4) {
    for(i=0; i<pxlNr; i++, mptr+=4, fptr++) {
      iptr=(int*)mptr;
      *fptr=h->scale_factor*(float)(*iptr);
    }
  } else if(h->data_type==ECAT7_VAXR4 || h->data_type==ECAT7_IEEER4) {
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
 * Read ECAT7 3D sinogram matrix header and data.
 *   Memory for fdata[] is allocated here, remember to free memory after usage.
 *   Note: data is converted to floats with scale_factor in the scan matrix header.
 *   Note: data is not calibrated with ecat_calibration_factor in main header.
 *   Note: data is not multiplied with deadtime_correction_factor.
 *
 * @param fp ECAT file pointer
 * @param first_block Subheader record number
 * @param last_block Last data block number
 * @param h Ptr to subheader data which is filled
 * @param fdata Ptr to the address of the matrix data
 * @return 0 if ok, 1 invalid input, 5 failed to read scan header,
 * 6 invalid image (x,y,z) dimensions, 8 failed to allocate memory for meta-data,
 * 9 failed to read matrix data, 11 failed to allocate memory for voxel data
 */
int ecat7ReadScanMatrix(
  FILE *fp, int first_block, int last_block, ECAT7_scanheader *h, float **fdata
) {
  int i, ret, blockNr, trueblockNr, pxlNr, dimz;
  char *mdata, *mptr;
  float *_fdata, *fptr;
  short int *sptr;
  int *iptr;
  
  
  if(ECAT7_TEST) printf("ecat7ReadScanMatrix(fp, %d, %d, hdr, fdata)\n",
    first_block, last_block);
  if(fp==NULL || first_block<=MatFirstDirBlk || h==NULL) {
    sprintf(ecat7errmsg, "invalid function parameter.\n");
    return(1);
  }
  *fdata=(float*)NULL;
  
  /* Read subheader */
  ret=ecat7ReadScanheader(fp, first_block, h);
  if(ret) {
    sprintf(ecat7errmsg, "cannot read subheader (%d).\n", ret);
    return(5);
  }
  if(ECAT7_TEST>4) ecat7PrintScanheader(h, stdout);
  pxlNr=h->num_r_elements*h->num_angles;
  for(i=dimz=0; i<64; i++) dimz+=h->num_z_elements[i];
  pxlNr*=dimz;
  if(pxlNr<=0) {
    sprintf(ecat7errmsg, "invalid matrix dimension.\n");
    return(6);  
  }
  trueblockNr=pxlNr*ecat7pxlbytes(h->data_type);
  trueblockNr=(trueblockNr+MatBLKSIZE-1)/MatBLKSIZE;

  /* Read matrix data; note that header takes 2 blocks */
  blockNr=last_block-first_block-1; if(blockNr<1) return(0);
  if(blockNr<trueblockNr) trueblockNr=blockNr;
  mdata=(char*)malloc(blockNr*MatBLKSIZE);
  if(mdata==NULL) {
    sprintf(ecat7errmsg, "cannot allocate memory.\n");
    return(8);  
  }
  mptr=mdata; /* note that only true block nr is read! */
  ret=ecat7ReadMatrixdata(fp, first_block+2, trueblockNr, mptr, h->data_type);
  if(ret || mdata==NULL) {
    sprintf(ecat7errmsg, "cannot read matrix data (%d).\n", ret);
    free(mdata); return(9);
  }
  
  /* Allocate memory for float data */
  _fdata=(float*)malloc(pxlNr*sizeof(float));
  if(_fdata==NULL) {
    sprintf(ecat7errmsg, "cannot allocate memory.\n");
    free(mdata); return(11);  
  }

  /* Convert matrix data to floats */
  fptr=_fdata; mptr=mdata;
  if(h->data_type==ECAT7_BYTE) {
    for(i=0; i<pxlNr; i++, mptr++, fptr++)
      *fptr=h->scale_factor*(float)(*mptr);
  } else if(h->data_type==ECAT7_VAXI2 || h->data_type==ECAT7_SUNI2) {
    for(i=0; i<pxlNr; i++, mptr+=2, fptr++) {
      sptr=(short int*)mptr;
      *fptr=h->scale_factor*(float)(*sptr);
    }
  } else if(h->data_type==ECAT7_VAXI4 || h->data_type==ECAT7_SUNI4) {
    for(i=0; i<pxlNr; i++, mptr+=4, fptr++) {
      iptr=(int*)mptr;
      *fptr=h->scale_factor*(float)(*iptr);
    }
  } else if(h->data_type==ECAT7_VAXR4 || h->data_type==ECAT7_IEEER4) {
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
 * Read ECAT7 polar map matrix header and data.
 * If only header is to be read, set last_block=first_block.
 * Note: data is not calibrated with factor in main header.
 *
 * @param fp ECAT file pointer
 * @param first_block Subheader record number
 * @param last_block Last data block number
 * @param h Ptr to subheader data which is filled
 * @param fdata Ptr to the address of the matrix data
 * @return 0 if ok, 1 invalid input, 5 failed to read scan header,
 * 6 invalid image (x,y,z) dimensions, 8 failed to allocate memory for meta-data,
 * 9 failed to read matrix data, 11 failed to allocate memory for voxel data
 */
int ecat7ReadPolarmapMatrix(
  FILE *fp, int first_block, int last_block, ECAT7_polmapheader *h, float **fdata
) {
  int i, ret, blockNr, pxlNr;
  char *mdata, *mptr;
  float *_fdata, *fptr;
  short int *sptr;
  int *iptr;


  if(ECAT7_TEST) printf("ecat7ReadPolarmapMatrix(fp, %d, %d, hdr, fdata)\n",
    first_block, last_block);
  if(fp==NULL || first_block<=MatFirstDirBlk || h==NULL) return 1;
  *fdata=(float*)NULL;

  /* Read subheader */
  ret=ecat7ReadPolmapheader(fp, first_block, h);
  if(ret) {
    sprintf(ecat7errmsg, "cannot read subheader (%d).\n", ret);
    return 2;
  }
  if(ECAT7_TEST>4) ecat7PrintPolmapheader(h, stdout);
  for(i=pxlNr=0; i<h->num_rings; i++) pxlNr+=h->sectors_per_ring[i];
  if(pxlNr<=0) return 3;

  /* Read matrix data */
  blockNr=last_block-first_block; if(blockNr<1) return 0;
  mdata=(char*)malloc(blockNr*MatBLKSIZE);
  if(mdata==NULL) return 4;
  mptr=mdata;
  ret=ecat7ReadMatrixdata(fp, first_block+1, blockNr, mptr, h->data_type);
  if(ret || mdata==NULL) {
    if(mdata!=NULL) free(mdata);
    return 5;
  }

  /* Allocate memory for float data */
  _fdata=(float*)malloc(pxlNr*sizeof(float));
  if(_fdata==NULL) {
    sprintf(ecat7errmsg, "cannot allocate memory.\n");
    free(mdata); return 4;
  }

  /* Convert matrix data to floats */
  fptr=_fdata; mptr=mdata;
  if(h->data_type==ECAT7_BYTE) {
    for(i=0; i<pxlNr; i++, mptr++, fptr++)
      *fptr=h->scale_factor*(float)(*mptr);
  } else if(h->data_type==ECAT7_VAXI2 || h->data_type==ECAT7_SUNI2) {
    for(i=0; i<pxlNr; i++, mptr+=2, fptr++) {
      sptr=(short int*)mptr;
      *fptr=h->scale_factor*(float)(*sptr);
    }
  } else if(h->data_type==ECAT7_VAXI4 || h->data_type==ECAT7_SUNI4) {
    for(i=0; i<pxlNr; i++, mptr+=4, fptr++) {
      iptr=(int*)mptr;
      *fptr=h->scale_factor*(float)(*iptr);
    }
  } else if(h->data_type==ECAT7_VAXR4 || h->data_type==ECAT7_IEEER4) {
    memcpy(fptr, mptr, pxlNr*4);
    for(i=0; i<pxlNr; i++, fptr++) *fptr *= h->scale_factor;
  }
  free(mdata);
  *fdata=_fdata;

  return 0;
}
/*****************************************************************************/

/*****************************************************************************/
/*!
 * Read ECAT7 floats
 *
 * @param bufi		pointer to 32-bit data block
 * @param isvax 	!= 0 for VAX format
 * @param islittle	!= 0 for little endian conversion
 * @return data in bufi as float value
 */
float ecat7rFloat(void *bufi, int isvax, int islittle) {
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
 * Reading and writing ECAT7 32-bit ints
 * 32-bit int format is same in VAX and i386
 *
 * @param bufi		pointer to one 32-bit data block
 * @param isvax         ignored
 * @param islittle	!= 0 for little endian conversion
 * @return converted 32-bit integer
 */
int ecat7rInt(void *bufi, int isvax, int islittle) {
  int i;

  if(isvax==0) {} // prevent compiler warning
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
 * @param data_type	defined value for data type
 * @return number of bytes (1,2 or 4) or 0 if type not recognized
 */
int ecat7pxlbytes(short int data_type) {
  int byteNr=0;
  switch(data_type) {
    case ECAT7_BYTE: byteNr=1; break;
    case ECAT7_VAXI2:
    case ECAT7_SUNI2: byteNr=2; break;
    case ECAT7_VAXI4:
    case ECAT7_VAXR4:
    case ECAT7_IEEER4:
    case ECAT7_SUNI4: byteNr=4; break;
  }
  return(byteNr);
}
/*****************************************************************************/

/*****************************************************************************/
