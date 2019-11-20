/// @file ecat7w.c
/// @author Vesa Oikonen
/// @brief Functions for writing ECAT 7.x format.
///
/*****************************************************************************/
#include "libtpcimgio.h"
/*****************************************************************************/

/*****************************************************************************/
/**
 * Write ECAT 7.x main header.
 * @details Writes header always in big endian byte order.
 * @return 0 in case of success, 1 == invalid parameters, 4 == file pointer is
 * at wrong position, 5 == writing of MatBLKSIZE bytes was not success
 */
int ecat7WriteMainheader(
  /** output file pointer */
  FILE *fp, 
  /** Ecat7 main header */
  ECAT7_mainheader *h
) {
  unsigned char buf[MatBLKSIZE];
  int little;

  if(ECAT7_TEST) printf("ecat7WriteMainheader()\n");
  /* Check arguments */
  if(fp==NULL || h==NULL) return(1);
  little=little_endian();
  /* Clear buf */
  memset(buf, 0, MatBLKSIZE);

  /* Copy header contents into buffer and change byte order if necessary */
  memcpy(buf+0, &h->magic_number, 14);
  memcpy(buf+14, &h->original_file_name, 32);
  memcpy(buf+46, &h->sw_version, 2); if(little) swabip(buf+46, 2);
  memcpy(buf+48, &h->system_type, 2); if(little) swabip(buf+48, 2);
  memcpy(buf+50, &h->file_type, 2); if(little) swabip(buf+50, 2);
  memcpy(buf+52, &h->serial_number, 10);
  //printf("ecat7WriteMainheader(): scan_start_time := %d\n", h->scan_start_time);
  memcpy(buf+62, &h->scan_start_time, 4); if(little) swawbip(buf+62, 4);
  memcpy(buf+66, &h->isotope_name, 8);
  memcpy(buf+74, &h->isotope_halflife, 4); if(little) swawbip(buf+74, 4);
  memcpy(buf+78, &h->radiopharmaceutical, 32);
  memcpy(buf+110, &h->gantry_tilt, 4); if(little) swawbip(buf+110, 4);
  memcpy(buf+114, &h->gantry_rotation, 4); if(little) swawbip(buf+114, 4);
  memcpy(buf+118, &h->bed_elevation, 4); if(little) swawbip(buf+118, 4);
  memcpy(buf+122, &h->intrinsic_tilt, 4); if(little) swawbip(buf+122, 4);
  memcpy(buf+126, &h->wobble_speed, 2); if(little) swabip(buf+126, 2);
  memcpy(buf+128, &h->transm_source_type, 2); if(little) swabip(buf+128, 2);
  memcpy(buf+130, &h->distance_scanned, 4); if(little) swawbip(buf+130, 4);
  memcpy(buf+134, &h->transaxial_fov, 4); if(little) swawbip(buf+134, 4);
  memcpy(buf+138, &h->angular_compression, 2); if(little) swabip(buf+138, 2);
  memcpy(buf+140, &h->coin_samp_mode, 2); if(little) swabip(buf+140, 2);
  memcpy(buf+142, &h->axial_samp_mode, 2); if(little) swabip(buf+142, 2);
  memcpy(buf+144, &h->ecat_calibration_factor, 4); if(little) swawbip(buf+144, 4);
  memcpy(buf+148, &h->calibration_units, 2); if(little) swabip(buf+148, 2);
  memcpy(buf+150, &h->calibration_units_label, 2); if(little) swabip(buf+150, 2);
  memcpy(buf+152, &h->compression_code, 2); if(little) swabip(buf+152, 2);
  memcpy(buf+154, &h->study_type, 12);
  memcpy(buf+166, &h->patient_id, 16);
  memcpy(buf+182, &h->patient_name, 32);
  memcpy(buf+214, &h->patient_sex, 1);
  memcpy(buf+215, &h->patient_dexterity, 1);
  memcpy(buf+216, &h->patient_age, 4); if(little) swawbip(buf+216, 4);
  memcpy(buf+220, &h->patient_height, 4); if(little) swawbip(buf+220, 4);
  memcpy(buf+224, &h->patient_weight, 4); if(little) swawbip(buf+224, 4);
  memcpy(buf+228, &h->patient_birth_date, 4); if(little) swawbip(buf+228, 4);
  memcpy(buf+232, &h->physician_name, 32);
  memcpy(buf+264, &h->operator_name, 32);
  memcpy(buf+296, &h->study_description, 32);
  memcpy(buf+328, &h->acquisition_type, 2); if(little) swabip(buf+328, 2);
  memcpy(buf+330, &h->patient_orientation, 2); if(little) swabip(buf+330, 2);
  memcpy(buf+332, &h->facility_name, 20);
  memcpy(buf+352, &h->num_planes, 2); if(little) swabip(buf+352, 2);
  memcpy(buf+354, &h->num_frames, 2); if(little) swabip(buf+354, 2);
  memcpy(buf+356, &h->num_gates, 2); if(little) swabip(buf+356, 2);
  memcpy(buf+358, &h->num_bed_pos, 2); if(little) swabip(buf+358, 2);
  memcpy(buf+360, &h->init_bed_position, 4); if(little) swawbip(buf+360, 4);
  memcpy(buf+364, h->bed_position, 15*4); if(little) swawbip(buf+364, 15*4);
  memcpy(buf+424, &h->plane_separation, 4); if(little) swawbip(buf+424, 4);
  memcpy(buf+428, &h->lwr_sctr_thres, 2); if(little) swabip(buf+428, 2);
  memcpy(buf+430, &h->lwr_true_thres, 2); if(little) swabip(buf+430, 2);
  memcpy(buf+432, &h->upr_true_thres, 2); if(little) swabip(buf+432, 2);
  memcpy(buf+434, &h->user_process_code, 10);
  memcpy(buf+444, &h->acquisition_mode, 2); if(little) swabip(buf+444, 2);
  memcpy(buf+446, &h->bin_size, 4); if(little) swawbip(buf+446, 4);
  memcpy(buf+450, &h->branching_fraction, 4); if(little) swawbip(buf+450, 4);
  memcpy(buf+454, &h->dose_start_time, 4); if(little) swawbip(buf+454, 4);
  memcpy(buf+458, &h->dosage, 4); if(little) swawbip(buf+458, 4);
  memcpy(buf+462, &h->well_counter_corr_factor, 4); if(little) swawbip(buf+462, 4);
  memcpy(buf+466, &h->data_units, 32);
  memcpy(buf+498, &h->septa_state, 2); if(little) swabip(buf+498, 2);
  memcpy(buf+500, &h->fill_cti, 12);

  /* Write main header */
  fseek(fp, 0*MatBLKSIZE, SEEK_SET); if(ftell(fp)!=0*MatBLKSIZE) return(4);
  if(fwrite(buf, 1, 1*MatBLKSIZE, fp) != 1*MatBLKSIZE) return(5);

  return(0);
}
/*****************************************************************************/

/*****************************************************************************/
/*!
 * Write ECAT 7.x image header. Changes data type to big endian.
 *
 * @param fp	output file pointer
 * @param blk	header block number, blk >= 2
 * @param h	Ecat7 image header
 * @return 0 in case of success, 1 == invalid parameters, 4 == file pointer is
 * at wrong position, 5 == writing of MatBLKSIZE bytes was not success
 */
int ecat7WriteImageheader(FILE *fp, int blk, ECAT7_imageheader *h) {
  unsigned char buf[MatBLKSIZE];
  int little; /* 1 if current platform is little endian (i386), else 0 */

  if(ECAT7_TEST) printf("ecat7WriteImageheader()\n");
  if(fp==NULL || blk<2 || h==NULL) return(1);
  little=little_endian(); if(ECAT7_TEST) printf("little=%d\n", little);
  /* Clear buf */
  memset(buf, 0, MatBLKSIZE);
  if(h->data_type==ECAT7_VAXI2) h->data_type=ECAT7_SUNI2;
  else if(h->data_type==ECAT7_VAXI4) h->data_type=ECAT7_SUNI4;
  else if(h->data_type==ECAT7_VAXR4) h->data_type=ECAT7_IEEER4;
  
  /* Copy the header fields and swap if necessary */
  memcpy(buf+0, &h->data_type, 2); if(little) swabip(buf+0, 2);
  memcpy(buf+2, &h->num_dimensions, 2); if(little) swabip(buf+2, 2);
  memcpy(buf+4, &h->x_dimension, 2); if(little) swabip(buf+4, 2);
  memcpy(buf+6, &h->y_dimension, 2); if(little) swabip(buf+6, 2);
  memcpy(buf+8, &h->z_dimension, 2); if(little) swabip(buf+8, 2);
  memcpy(buf+10, &h->x_offset, 4); if(little) swawbip(buf+10, 4);
  memcpy(buf+14, &h->y_offset, 4); if(little) swawbip(buf+14, 4);
  memcpy(buf+18, &h->z_offset, 4); if(little) swawbip(buf+18, 4);
  memcpy(buf+22, &h->recon_zoom, 4); if(little) swawbip(buf+22, 4);
  memcpy(buf+26, &h->scale_factor, 4); if(little) swawbip(buf+26, 4);
  memcpy(buf+30, &h->image_min, 2); if(little) swabip(buf+30, 2);
  memcpy(buf+32, &h->image_max, 2); if(little) swabip(buf+32, 2);
  memcpy(buf+34, &h->x_pixel_size, 4); if(little) swawbip(buf+34, 4);
  memcpy(buf+38, &h->y_pixel_size, 4); if(little) swawbip(buf+38, 4);
  memcpy(buf+42, &h->z_pixel_size, 4); if(little) swawbip(buf+42, 4);
  memcpy(buf+46, &h->frame_duration, 4); if(little) swawbip(buf+46, 4);
  memcpy(buf+50, &h->frame_start_time, 4); if(little) swawbip(buf+50, 4);
  memcpy(buf+54, &h->filter_code, 2); if(little) swabip(buf+54, 2);
  memcpy(buf+56, &h->x_resolution, 4); if(little) swawbip(buf+56, 4);
  memcpy(buf+60, &h->y_resolution, 4); if(little) swawbip(buf+60, 4);
  memcpy(buf+64, &h->z_resolution, 4); if(little) swawbip(buf+64, 4);
  memcpy(buf+68, &h->num_r_elements, 4); if(little) swawbip(buf+68, 4);
  memcpy(buf+72, &h->num_angles, 4); if(little) swawbip(buf+72, 4);
  memcpy(buf+76, &h->z_rotation_angle, 4); if(little) swawbip(buf+76, 4);
  memcpy(buf+80, &h->decay_corr_fctr, 4); if(little) swawbip(buf+80, 4);
  memcpy(buf+84, &h->processing_code, 4); if(little) swawbip(buf+84, 4);
  memcpy(buf+88, &h->gate_duration, 4); if(little) swawbip(buf+88, 4);
  memcpy(buf+92, &h->r_wave_offset, 4); if(little) swawbip(buf+92, 4);
  memcpy(buf+96, &h->num_accepted_beats, 4); if(little) swawbip(buf+96, 4);
  memcpy(buf+100, &h->filter_cutoff_frequency, 4); if(little) swawbip(buf+100, 4);
  memcpy(buf+104, &h->filter_resolution, 4); if(little) swawbip(buf+104, 4);
  memcpy(buf+108, &h->filter_ramp_slope, 4); if(little) swawbip(buf+108, 4);
  memcpy(buf+112, &h->filter_order, 2); if(little) swabip(buf+112, 2);
  memcpy(buf+114, &h->filter_scatter_fraction, 4); if(little) swawbip(buf+114, 4);
  memcpy(buf+118, &h->filter_scatter_slope, 4); if(little) swawbip(buf+118, 4);
  memcpy(buf+122, &h->annotation, 40);
  memcpy(buf+162, &h->mt_1_1, 4); if(little) swawbip(buf+162, 4);
  memcpy(buf+166, &h->mt_1_2, 4); if(little) swawbip(buf+166, 4);
  memcpy(buf+170, &h->mt_1_3, 4); if(little) swawbip(buf+170, 4);
  memcpy(buf+174, &h->mt_2_1, 4); if(little) swawbip(buf+174, 4);
  memcpy(buf+178, &h->mt_2_2, 4); if(little) swawbip(buf+178, 4);
  memcpy(buf+182, &h->mt_2_3, 4); if(little) swawbip(buf+182, 4);
  memcpy(buf+186, &h->mt_3_1, 4); if(little) swawbip(buf+186, 4);
  memcpy(buf+190, &h->mt_3_2, 4); if(little) swawbip(buf+190, 4);
  memcpy(buf+194, &h->mt_3_3, 4); if(little) swawbip(buf+194, 4);
  memcpy(buf+198, &h->rfilter_cutoff, 4); if(little) swawbip(buf+198, 4);
  memcpy(buf+202, &h->rfilter_resolution, 4); if(little) swawbip(buf+202, 4);
  memcpy(buf+206, &h->rfilter_code, 2); if(little) swabip(buf+206, 2);
  memcpy(buf+208, &h->rfilter_order, 2); if(little) swabip(buf+208, 2);
  memcpy(buf+210, &h->zfilter_cutoff, 4); if(little) swawbip(buf+210, 4);
  memcpy(buf+214, &h->zfilter_resolution, 4); if(little) swawbip(buf+214, 4);
  memcpy(buf+218, &h->zfilter_code, 2); if(little) swabip(buf+218, 2);
  memcpy(buf+220, &h->zfilter_order, 2); if(little) swabip(buf+220, 2);
  memcpy(buf+222, &h->mt_1_4, 4); if(little) swawbip(buf+222, 4);
  memcpy(buf+226, &h->mt_2_4, 4); if(little) swawbip(buf+226, 4);
  memcpy(buf+230, &h->mt_3_4, 4); if(little) swawbip(buf+230, 4);
  memcpy(buf+234, &h->scatter_type, 2); if(little) swabip(buf+234, 2);
  memcpy(buf+236, &h->recon_type, 2); if(little) swabip(buf+236, 2);
  memcpy(buf+238, &h->recon_views, 2); if(little) swabip(buf+238, 2);
  memcpy(buf+240, &h->fill_cti, 87);
  memcpy(buf+414, &h->fill_user, 48);

  /* Write header */
  fseek(fp, (blk-1)*MatBLKSIZE, SEEK_SET); if(ftell(fp)!=(blk-1)*MatBLKSIZE) return(4);
  if(fwrite(buf, 1, 1*MatBLKSIZE, fp) != 1*MatBLKSIZE) return(5);

  return(0);
}
/*****************************************************************************/

/*****************************************************************************/
/*!
 * Write ECAT 7.x attenuation header
 *
 * @param fp	output file pointer
 * @param blk	header block number, blk >= 2
 * @param h	Ecat7 attenuation header
 * @return 0 in case of success, 1 == invalid parameters, 4 == file pointer is
 * at wrong position, 5 == writing of MatBLKSIZE bytes was not success
*/
int ecat7WriteAttenheader(FILE *fp, int blk, ECAT7_attenheader *h) {
  unsigned char buf[MatBLKSIZE];
  int little; /* 1 if current platform is little endian (i386), else 0 */

  if(ECAT7_TEST) printf("ecat7WriteAttenheader()\n");
  if(fp==NULL || blk<2 || h==NULL) return(1);
  little=little_endian(); if(ECAT7_TEST) printf("little=%d\n", little);
  /* Clear buf */
  memset(buf, 0, MatBLKSIZE);
  if(h->data_type==ECAT7_VAXI2) h->data_type=ECAT7_SUNI2;
  else if(h->data_type==ECAT7_VAXI4) h->data_type=ECAT7_SUNI4;
  else if(h->data_type==ECAT7_VAXR4) h->data_type=ECAT7_IEEER4;
  
  /* Copy the header fields and swap if necessary */
  memcpy(buf+0, &h->data_type, 2); if(little) swabip(buf+0, 2);
  memcpy(buf+2, &h->num_dimensions, 2); if(little) swabip(buf+2, 2);
  memcpy(buf+4, &h->attenuation_type, 2); if(little) swabip(buf+4, 2);
  memcpy(buf+6, &h->num_r_elements, 2); if(little) swabip(buf+6, 2);
  memcpy(buf+8, &h->num_angles, 2); if(little) swabip(buf+8, 2);
  memcpy(buf+10, &h->num_z_elements, 2); if(little) swabip(buf+10, 2);
  memcpy(buf+12, &h->ring_difference, 2); if(little) swabip(buf+12, 2);
  memcpy(buf+14, &h->x_resolution, 4); if(little) swawbip(buf+14, 4);
  memcpy(buf+18, &h->y_resolution, 4); if(little) swawbip(buf+18, 4);
  memcpy(buf+22, &h->z_resolution, 4); if(little) swawbip(buf+22, 4);
  memcpy(buf+26, &h->w_resolution, 4); if(little) swawbip(buf+26, 4);
  memcpy(buf+30, &h->scale_factor, 4); if(little) swawbip(buf+30, 4);
  memcpy(buf+34, &h->x_offset, 4); if(little) swawbip(buf+34, 4);
  memcpy(buf+38, &h->y_offset, 4); if(little) swawbip(buf+38, 4);
  memcpy(buf+42, &h->x_radius, 4); if(little) swawbip(buf+42, 4);
  memcpy(buf+46, &h->y_radius, 4); if(little) swawbip(buf+46, 4);
  memcpy(buf+50, &h->tilt_angle, 4); if(little) swawbip(buf+50, 4);
  memcpy(buf+54, &h->attenuation_coeff, 4); if(little) swawbip(buf+54, 4);
  memcpy(buf+58, &h->attenuation_min, 4); if(little) swawbip(buf+58, 4);
  memcpy(buf+62, &h->attenuation_max, 4); if(little) swawbip(buf+62, 4);
  memcpy(buf+66, &h->skull_thickness, 4); if(little) swawbip(buf+66, 4);
  memcpy(buf+70, &h->num_additional_atten_coeff, 2); if(little) swabip(buf+70, 2);
  memcpy(buf+72, h->additional_atten_coeff, 8*4); if(little) swawbip(buf+72, 8*4);
  memcpy(buf+104, &h->edge_finding_threshold, 4); if(little) swawbip(buf+104, 4);
  memcpy(buf+108, &h->storage_order, 2); if(little) swabip(buf+108, 2);
  memcpy(buf+110, &h->span, 2); if(little) swabip(buf+110, 2);
  memcpy(buf+112, h->z_elements, 64*2); if(little) swabip(buf+112, 64*2);
  memcpy(buf+240, h->fill_cti, 86*2); if(little) swabip(buf+240, 86*2);
  memcpy(buf+412, h->fill_user, 50*2); if(little) swabip(buf+412, 50*2);

  /* Write header */
  fseek(fp, (blk-1)*MatBLKSIZE, SEEK_SET);
  if(ftell(fp)!=(blk-1)*MatBLKSIZE) return(4);
  if(fwrite(buf, 1, 1*MatBLKSIZE, fp) != 1*MatBLKSIZE) return(5);

  return(0);
}
/*****************************************************************************/

/*****************************************************************************/
/*!
 * Write ECAT 7.x polar map header
 *
 * @param fp	output file pointer
 * @param blk	header block number, blk >= 2
 * @param h	Ecat7 polar map header
 * @return 0 in case of success, 1 == invalid parameters, 4 == file pointer is
 * at wrong position, 5 == writing of MatBLKSIZE bytes was not success
*/
int ecat7WritePolmapheader(FILE *fp, int blk, ECAT7_polmapheader *h) {
  unsigned char buf[MatBLKSIZE];
  int little; /* 1 if current platform is little endian (i386), else 0 */

  if(ECAT7_TEST) printf("ecat7WritePolmapheader()\n");
  if(fp==NULL || blk<2 || h==NULL) return(1);
  little=little_endian(); if(ECAT7_TEST) printf("little=%d\n", little);
  /* Clear buf */
  memset(buf, 0, MatBLKSIZE);
  if(h->data_type==ECAT7_VAXI2) h->data_type=ECAT7_SUNI2;
  else if(h->data_type==ECAT7_VAXI4) h->data_type=ECAT7_SUNI4;
  else if(h->data_type==ECAT7_VAXR4) h->data_type=ECAT7_IEEER4;
  
  /* Copy the header fields and swap if necessary */
  memcpy(buf+0, &h->data_type, 2); if(little) swabip(buf+0, 2);
  memcpy(buf+2, &h->polar_map_type, 2); if(little) swabip(buf+2, 2);
  memcpy(buf+4, &h->num_rings, 2); if(little) swabip(buf+4, 2);
  memcpy(buf+6, h->sectors_per_ring, 32*2); if(little) swabip(buf+6, 32*2);
  memcpy(buf+70, h->ring_position, 32*4); if(little) swawbip(buf+70, 32*4);
  memcpy(buf+198, h->ring_angle, 32*2); if(little) swabip(buf+198, 32*2);
  memcpy(buf+262, &h->start_angle, 2); if(little) swabip(buf+262, 2);
  memcpy(buf+264, h->long_axis_left, 3*2); if(little) swabip(buf+264, 3*2);
  memcpy(buf+270, h->long_axis_right, 3*2); if(little) swabip(buf+270, 3*2);
  memcpy(buf+276, &h->position_data, 2); if(little) swabip(buf+276, 2);
  memcpy(buf+278, &h->image_min, 2); if(little) swabip(buf+278, 2);
  memcpy(buf+280, &h->image_max, 2); if(little) swabip(buf+280, 2);
  memcpy(buf+282, &h->scale_factor, 4); if(little) swawbip(buf+282, 4);
  memcpy(buf+286, &h->pixel_size, 4); if(little) swawbip(buf+286, 4);
  memcpy(buf+290, &h->frame_duration, 4); if(little) swawbip(buf+290, 4);
  memcpy(buf+294, &h->frame_start_time, 4); if(little) swawbip(buf+294, 4);
  memcpy(buf+298, &h->processing_code, 2); if(little) swabip(buf+298, 2);
  memcpy(buf+300, &h->quant_units, 2); if(little) swabip(buf+300, 2);
  memcpy(buf+302, h->annotation, 40);
  memcpy(buf+342, &h->gate_duration, 4); if(little) swawbip(buf+342, 4);
  memcpy(buf+346, &h->r_wave_offset, 4); if(little) swawbip(buf+346, 4);
  memcpy(buf+350, &h->num_accepted_beats, 4); if(little) swawbip(buf+350, 4);
  memcpy(buf+354, h->polar_map_protocol, 20);
  memcpy(buf+374, h->database_name, 30);
  memcpy(buf+404, h->fill_cti, 27*2); if(little) swabip(buf+404, 27*2);

  /* Write header */
  fseek(fp, (blk-1)*MatBLKSIZE, SEEK_SET);
  if(ftell(fp)!=(blk-1)*MatBLKSIZE) return(4);
  if(fwrite(buf, 1, 1*MatBLKSIZE, fp) != 1*MatBLKSIZE) return(5);

  return(0);
}
/*****************************************************************************/

/*****************************************************************************/
/*!
 * Write ECAT 7.x 3D normalization header
 *
 * @param fp	output file pointer
 * @param blk	header block number, blk >= 2
 * @param h	Ecat7 normalization header
 * @return 0 in case of success, 1 == invalid parameters, 4 == file pointer is
 * at wrong position, 5 == writing of MatBLKSIZE bytes was not success
 */
int ecat7WriteNormheader(FILE *fp, int blk, ECAT7_normheader *h) {
  unsigned char buf[MatBLKSIZE];
  int little; /* 1 if current platform is little endian (i386), else 0 */

  if(ECAT7_TEST) printf("ecat7WriteNormheader()\n");
  if(fp==NULL || blk<2 || h==NULL) return(1);
  little=little_endian(); if(ECAT7_TEST) printf("little=%d\n", little);
  /* Clear buf */
  memset(buf, 0, MatBLKSIZE);
  if(h->data_type==ECAT7_VAXI2) h->data_type=ECAT7_SUNI2;
  else if(h->data_type==ECAT7_VAXI4) h->data_type=ECAT7_SUNI4;
  else if(h->data_type==ECAT7_VAXR4) h->data_type=ECAT7_IEEER4;
  
  /* Copy the header fields and swap if necessary */
  memcpy(buf+0, &h->data_type, 2); if(little) swabip(buf+0, 2);
  memcpy(buf+2, &h->num_r_elements, 2); if(little) swabip(buf+2, 2);
  memcpy(buf+4, &h->num_transaxial_crystals, 2); if(little) swabip(buf+4, 2);
  memcpy(buf+6, &h->num_crystal_rings, 2); if(little) swabip(buf+6, 2);
  memcpy(buf+8, &h->crystals_per_ring, 2); if(little) swabip(buf+8, 2);
  memcpy(buf+10, &h->num_geo_corr_planes, 2); if(little) swabip(buf+10, 2);
  memcpy(buf+12, &h->uld, 2); if(little) swabip(buf+12, 2);
  memcpy(buf+14, &h->lld, 2); if(little) swabip(buf+14, 2);
  memcpy(buf+16, &h->scatter_energy, 2); if(little) swabip(buf+16, 2);
  memcpy(buf+18, &h->norm_quality_factor, 4); if(little) swawbip(buf+18, 4);
  memcpy(buf+22, &h->norm_quality_factor_code, 2); if(little) swabip(buf+22, 2);
  memcpy(buf+24, h->ring_dtcor1, 32*4); if(little) swawbip(buf+24, 32*4);
  memcpy(buf+152, h->ring_dtcor2, 32*4); if(little) swawbip(buf+152, 32*4);
  memcpy(buf+280, h->crystal_dtcor, 8*4); if(little) swawbip(buf+280, 8*4);
  memcpy(buf+312, &h->span, 2); if(little) swabip(buf+312, 2);
  memcpy(buf+314, &h->max_ring_diff, 2); if(little) swabip(buf+314, 2);
  memcpy(buf+316, h->fill_cti, 48*2); if(little) swabip(buf+316, 48*2);
  memcpy(buf+412, h->fill_user, 50*2); if(little) swabip(buf+412, 50*2);

  /* Write header */
  fseek(fp, (blk-1)*MatBLKSIZE, SEEK_SET);
  if(ftell(fp)!=(blk-1)*MatBLKSIZE) return(4);
  if(fwrite(buf, 1, 1*MatBLKSIZE, fp) != 1*MatBLKSIZE) return(5);

  return(0);
}
/*****************************************************************************/

/*****************************************************************************/
/*!
 * Write ECAT 7.x 3D scan header (512 bytes)
 * Changes data type to big endian.
 *
 * @param fp		pointer to output file
 * @param blk		block number, blk >= 2
 * @param h		Ecat7 scan header
 * @return 0 in case of success, 1 == invalid parameters, 4 == file pointer is
 * at wrong position, 5 == writing of MatBLKSIZE bytes was not success
 */
int ecat7WriteScanheader(FILE *fp, int blk, ECAT7_scanheader *h) {
  unsigned char buf[2*MatBLKSIZE];
  int little; /* 1 if current platform is little endian (i386), else 0 */

  if(ECAT7_TEST) printf("ecat7WriteScanheader()\n");
  if(fp==NULL || blk<2 || h==NULL) return(1);
  little=little_endian(); if(ECAT7_TEST) printf("little=%d\n", little);
  /* Clear buf */
  memset(buf, 0, 2*MatBLKSIZE);
  if(h->data_type==ECAT7_VAXI2) h->data_type=ECAT7_SUNI2;
  else if(h->data_type==ECAT7_VAXI4) h->data_type=ECAT7_SUNI4;
  else if(h->data_type==ECAT7_VAXR4) h->data_type=ECAT7_IEEER4;

  /* Copy the header fields and swap if necessary */
  memcpy(buf+0, &h->data_type, 2); if(little) swabip(buf+0, 2);
  memcpy(buf+2, &h->num_dimensions, 2); if(little) swabip(buf+2, 2);
  memcpy(buf+4, &h->num_r_elements, 2); if(little) swabip(buf+4, 2);
  memcpy(buf+6, &h->num_angles, 2); if(little) swabip(buf+6, 2);
  memcpy(buf+8, &h->corrections_applied, 2); if(little) swabip(buf+8, 2);
  memcpy(buf+10, h->num_z_elements, 64*2); if(little) swabip(buf+10, 64*2);
  memcpy(buf+138, &h->ring_difference, 2); if(little) swabip(buf+138, 2);
  memcpy(buf+140, &h->storage_order, 2); if(little) swabip(buf+140, 2);
  memcpy(buf+142, &h->axial_compression, 2); if(little) swabip(buf+142, 2);
  memcpy(buf+144, &h->x_resolution, 4); if(little) swawbip(buf+144, 4);
  memcpy(buf+148, &h->v_resolution, 4); if(little) swawbip(buf+148, 4);
  memcpy(buf+152, &h->z_resolution, 4); if(little) swawbip(buf+152, 4);
  memcpy(buf+156, &h->w_resolution, 4); if(little) swawbip(buf+156, 4);
  memcpy(buf+160, h->fill_gate, 6*2); if(little) swabip(buf+160, 6*2);
  memcpy(buf+172, &h->gate_duration, 4); if(little) swawbip(buf+172, 4);
  memcpy(buf+176, &h->r_wave_offset, 4); if(little) swawbip(buf+176, 4);
  memcpy(buf+180, &h->num_accepted_beats, 4); if(little) swawbip(buf+180, 4);
  memcpy(buf+184, &h->scale_factor, 4); if(little) swawbip(buf+184, 4);
  memcpy(buf+188, &h->scan_min, 2); if(little) swabip(buf+188, 2);
  memcpy(buf+190, &h->scan_max, 2); if(little) swabip(buf+190, 2);
  memcpy(buf+192, &h->prompts, 4); if(little) swawbip(buf+192, 4);
  memcpy(buf+196, &h->delayed, 4); if(little) swawbip(buf+196, 4);
  memcpy(buf+200, &h->multiples, 4); if(little) swawbip(buf+200, 4);
  memcpy(buf+204, &h->net_trues, 4); if(little) swawbip(buf+204, 4);
  memcpy(buf+208, &h->tot_avg_cor, 4); if(little) swawbip(buf+208, 4);
  memcpy(buf+212, &h->tot_avg_uncor, 4); if(little) swawbip(buf+212, 4);
  memcpy(buf+216, &h->total_coin_rate, 4); if(little) swawbip(buf+216, 4);
  memcpy(buf+220, &h->frame_start_time, 4); if(little) swawbip(buf+220, 4);
  memcpy(buf+224, &h->frame_duration, 4); if(little) swawbip(buf+224, 4);
  memcpy(buf+228, &h->deadtime_correction_factor, 4); if(little) swawbip(buf+228, 4);
  memcpy(buf+232, h->fill_cti, 90*2); if(little) swabip(buf+232, 90*2);
  memcpy(buf+412, h->fill_user, 50*2); if(little) swabip(buf+412, 50*2);
  memcpy(buf+512, h->uncor_singles, 128*4); if(little) swawbip(buf+512, 128*4);

  /* Write 3D scan header */
  fseek(fp, (blk-1)*MatBLKSIZE, SEEK_SET); if(ftell(fp)!=(blk-1)*MatBLKSIZE) return(4);
  if(fwrite(buf, 1, 2*MatBLKSIZE, fp) != 2*MatBLKSIZE) return(5);

  return(0);
}
/*****************************************************************************/

/*****************************************************************************/
/*!
 * Write ECAT 7.x 2D scan header
 *
 * @param fp		output file pointer
 * @param blk		header block number, blk >= 2
 * @param h		Ecat7 2D scan header
 * @return 0 in case of success, 1 == invalid parameters, 4 == file pointer is
 * at wrong position, 5 == writing of MatBLKSIZE bytes was not success
 */
int ecat7Write2DScanheader(FILE *fp, int blk, ECAT7_2Dscanheader *h) {
  unsigned char buf[MatBLKSIZE];
  int little; /* 1 if current platform is little endian (i386), else 0 */

  if(ECAT7_TEST) printf("ecat7Write2DScanheader()\n");
  if(fp==NULL || blk<2 || h==NULL) return(1);
  little=little_endian(); if(ECAT7_TEST) printf("little=%d\n", little);
  /* Clear buf */
  memset(buf, 0, MatBLKSIZE);
  if(h->data_type==ECAT7_VAXI2) h->data_type=ECAT7_SUNI2;
  else if(h->data_type==ECAT7_VAXI4) h->data_type=ECAT7_SUNI4;
  else if(h->data_type==ECAT7_VAXR4) h->data_type=ECAT7_IEEER4;

  /* Copy the header fields and swap if necessary */
  memcpy(buf+0, &h->data_type, 2); if(little) swabip(buf+0, 2);
  memcpy(buf+2, &h->num_dimensions, 2); if(little) swabip(buf+2, 2);
  memcpy(buf+4, &h->num_r_elements, 2); if(little) swabip(buf+4, 2);
  memcpy(buf+6, &h->num_angles, 2); if(little) swabip(buf+6, 2);
  memcpy(buf+8, &h->corrections_applied, 2); if(little) swabip(buf+8, 2);
  memcpy(buf+10, &h->num_z_elements, 2); if(little) swabip(buf+10, 2);
  memcpy(buf+12, &h->ring_difference, 2); if(little) swabip(buf+12, 2);
  memcpy(buf+14, &h->x_resolution, 4); if(little) swawbip(buf+14, 4);
  memcpy(buf+18, &h->y_resolution, 4); if(little) swawbip(buf+18, 4);
  memcpy(buf+22, &h->z_resolution, 4); if(little) swawbip(buf+22, 4);
  memcpy(buf+26, &h->w_resolution, 4); if(little) swawbip(buf+26, 4);
  memcpy(buf+30, h->fill_gate, 6*2); if(little) swabip(buf+30, 6*2);
  memcpy(buf+42, &h->gate_duration, 4); if(little) swawbip(buf+42, 4);
  memcpy(buf+46, &h->r_wave_offset, 4); if(little) swawbip(buf+46, 4);
  memcpy(buf+50, &h->num_accepted_beats, 4); if(little) swawbip(buf+50, 4);
  memcpy(buf+54, &h->scale_factor, 4); if(little) swawbip(buf+54, 4);
  memcpy(buf+58, &h->scan_min, 2); if(little) swabip(buf+58, 2);
  memcpy(buf+60, &h->scan_max, 2); if(little) swabip(buf+60, 2);
  memcpy(buf+62, &h->prompts, 4); if(little) swawbip(buf+62, 4);
  memcpy(buf+66, &h->delayed, 4); if(little) swawbip(buf+66, 4);
  memcpy(buf+70, &h->multiples, 4); if(little) swawbip(buf+70, 4);
  memcpy(buf+74, &h->net_trues, 4); if(little) swawbip(buf+74, 4);
  memcpy(buf+78, h->cor_singles, 16*4); if(little) swawbip(buf+78, 16*4);
  memcpy(buf+142, h->uncor_singles, 16*4); if(little) swawbip(buf+142, 16*4);
  memcpy(buf+206, &h->tot_avg_cor, 4); if(little) swawbip(buf+206, 4);
  memcpy(buf+210, &h->tot_avg_uncor, 4); if(little) swawbip(buf+210, 4);
  memcpy(buf+214, &h->total_coin_rate, 4); if(little) swawbip(buf+214, 4);
  memcpy(buf+218, &h->frame_start_time, 4); if(little) swawbip(buf+218, 4);
  memcpy(buf+222, &h->frame_duration, 4); if(little) swawbip(buf+222, 4);
  memcpy(buf+226, &h->deadtime_correction_factor, 4); if(little) swawbip(buf+226, 4);
  memcpy(buf+230, h->physical_planes, 8*2); if(little) swabip(buf+230, 8*2);
  memcpy(buf+246, h->fill_cti, 83*2); if(little) swabip(buf+246, 83*2);
  memcpy(buf+412, h->fill_user, 50*2); if(little) swabip(buf+412, 50*2);

  /* Write header */
  fseek(fp, (blk-1)*MatBLKSIZE, SEEK_SET);
  if(ftell(fp)!=(blk-1)*MatBLKSIZE) return(4);
  if(fwrite(buf, 1, 1*MatBLKSIZE, fp) != 1*MatBLKSIZE) return(5);

  return(0);
}
/*****************************************************************************/

/*****************************************************************************/
/*!
 * Write ECAT 7.x 2D normalization header.
 *
 * @param fp	file pointer
 * @param blk   header block number, blk >= 2
 * @param h	Ecat7 2D normalization header
 * @return 0 in case of success, 1 == invalid parameters, 4 == file pointer is
 * at wrong position, 5 == writing of MatBLKSIZE bytes was not success
 */
int ecat7Write2DNormheader(FILE *fp, int blk, ECAT7_2Dnormheader *h) {
  unsigned char buf[MatBLKSIZE];
  int little; /* 1 if current platform is little endian (i386), else 0 */

  if(ECAT7_TEST) printf("ecat7Write2DNormheader()\n");
  if(fp==NULL || blk<2 || h==NULL) return(1);
  little=little_endian(); if(ECAT7_TEST) printf("little=%d\n", little);
  /* Clear buf */
  memset(buf, 0, MatBLKSIZE);
  if(h->data_type==ECAT7_VAXI2) h->data_type=ECAT7_SUNI2;
  else if(h->data_type==ECAT7_VAXI4) h->data_type=ECAT7_SUNI4;
  else if(h->data_type==ECAT7_VAXR4) h->data_type=ECAT7_IEEER4;

  /* Copy the header fields and swap if necessary */
  memcpy(buf+0, &h->data_type, 2); if(little) swabip(buf+0, 2);
  memcpy(buf+2, &h->num_dimensions, 2); if(little) swabip(buf+2, 2);
  memcpy(buf+4, &h->num_r_elements, 2); if(little) swabip(buf+4, 2);
  memcpy(buf+6, &h->num_angles, 2); if(little) swabip(buf+6, 2);
  memcpy(buf+8, &h->num_z_elements, 2); if(little) swabip(buf+8, 2);
  memcpy(buf+10, &h->ring_difference, 2); if(little) swabip(buf+10, 2);
  memcpy(buf+12, &h->scale_factor, 4); if(little) swawbip(buf+12, 4);
  memcpy(buf+16, &h->norm_min, 4); if(little) swawbip(buf+16, 4);
  memcpy(buf+20, &h->norm_max, 4); if(little) swawbip(buf+20, 4);
  memcpy(buf+24, &h->fov_source_width, 4); if(little) swawbip(buf+24, 4);
  memcpy(buf+28, &h->norm_quality_factor, 4); if(little) swawbip(buf+28, 4);
  memcpy(buf+32, &h->norm_quality_factor_code, 2); if(little) swabip(buf+32, 2);
  memcpy(buf+34, &h->storage_order, 2); if(little) swabip(buf+34, 2);
  memcpy(buf+36, &h->span, 2); if(little) swabip(buf+36, 2);
  memcpy(buf+38, h->fill_cti, 64*2); if(little) swabip(buf+38, 64*2);
  memcpy(buf+166, h->fill_cti, 123*2); if(little) swabip(buf+166, 123*2);
  memcpy(buf+412, h->fill_user, 50*2); if(little) swabip(buf+412, 50*2);

  /* Write header */
  fseek(fp, (blk-1)*MatBLKSIZE, SEEK_SET);
  if(ftell(fp)!=(blk-1)*MatBLKSIZE) return(4);
  if(fwrite(buf, 1, 1*MatBLKSIZE, fp) != 1*MatBLKSIZE) return(5);

  return(0);
}
/*****************************************************************************/

/*****************************************************************************/
/*!
 * Create a new ECAT 7.x file. If file exists, it is renamed as fname% if possible.
 *  Directory list is written in big endian byte order.
 *
 * @param fname		filename
 * @param h		Ecat7 main header
 * @return file pointer or NULL in case of an error.
 */
FILE *ecat7Create(const char *fname, ECAT7_mainheader *h) {
  FILE *fp;
  char tmp[FILENAME_MAX];
  int buf[MatBLKSIZE/4];

  if(ECAT7_TEST) printf("ecat7Create(%s, h)\n", fname);
  /* Check the arguments */
  if(fname==NULL || h==NULL) return(NULL);
  /* Check if file exists; backup, if necessary */
  if(access(fname, 0) != -1) {
    strcpy(tmp, fname); strcat(tmp, BACKUP_EXTENSION);
    if(access(tmp, 0) != -1) remove(tmp);
    if(ECAT7_TEST) printf("Renaming %s -> %s\n", fname, tmp);
    rename(fname, tmp);
  }
  /* Open file */
  fp=fopen(fname, "wb+"); if(fp==NULL) return(fp);
  /* Write main header */
  if(ecat7WriteMainheader(fp, h)) return(NULL);
  /* Construct an empty matrix list ; convert to little endian if necessary */
  memset(buf, 0, MatBLKSIZE);
  buf[0]=31; buf[1]=MatFirstDirBlk; if(little_endian()) swawbip(buf, MatBLKSIZE);
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
 * Check if pixel float values need to be scaled to be saved as short ints,
 * or if they are already all very close to integers.
 *
 * @param amax		absolute maximum value
 * @param data		float array
 * @param nr		float array size
 * @return 1, if scaling is necessary, and 0 if not.
 */
int ecat7_is_scaling_needed(float amax, float *data, int nr) {
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
 * Write ECAT 7.x image or volume matrix header and data
 *
 * @param fp		output file pointer
 * @param matrix_id	coded matrix id
 * @param h		Ecat7 image header
 * @param fdata 	float data to be written
 * @return 0 if ok.
 */
int ecat7WriteImageMatrix(FILE *fp, int matrix_id, ECAT7_imageheader *h, float *fdata) {
  int i, nxtblk, blkNr, data_size, pxlNr, ret;
  float *fptr, fmin, fmax, g, f;
  char *mdata, *mptr;
  short int *sptr;
  
  
  if(ECAT7_TEST) printf("ecat7WriteImageMatrix(fp, %d, h, data)\n", matrix_id);
  if(fp==NULL || matrix_id<1 || h==NULL || fdata==NULL) {
    sprintf(ecat7errmsg, "invalid function parameter.\n");
    return(1);
  }
  if(h->data_type!=ECAT7_SUNI2) {
    sprintf(ecat7errmsg, "invalid data_type.\n");
    return(2);
  }
  /* nr of pixels */
  pxlNr=h->x_dimension*h->y_dimension;
  if(h->num_dimensions>2) pxlNr*=h->z_dimension;
  if(pxlNr<1) {
    sprintf(ecat7errmsg, "invalid matrix dimension.\n");
    return(3);
  }
  /* How much memory is needed for ALL pixels */
  data_size=pxlNr*ecat7pxlbytes(h->data_type);
  /* block nr taken by all pixels */
  blkNr=(data_size+MatBLKSIZE-1)/MatBLKSIZE; if(blkNr<1) {
    sprintf(ecat7errmsg, "invalid block number.\n");
    return(4);
  }
  /* Allocate memory for matrix data */
  mdata=(char*)calloc(blkNr, MatBLKSIZE); if(mdata==NULL) {
    sprintf(ecat7errmsg, "out of memory.\n");
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
  if(f>=1.0 && ecat7_is_scaling_needed(g, fptr, pxlNr)==0) f=1.0;
  /* Scale matrix data to shorts */
  h->scale_factor=1.0/f;
  sptr=(short int*)mdata; fptr=fdata;
  for(i=0; i<pxlNr; i++, sptr++, fptr++) *sptr=(short int)temp_roundf(f*(*fptr));
  /* Set header short min & max */
  h->image_min=(short int)temp_roundf(f*fmin);
  h->image_max=(short int)temp_roundf(f*fmax);
  /* Get block number for matrix header and data */
  nxtblk=ecat7EnterMatrix(fp, matrix_id, blkNr); if(nxtblk<1) {
    sprintf(ecat7errmsg, "cannot determine matrix block (%d).\n", -nxtblk);
    free(mdata); return(8);
  }
  if(ECAT7_TEST>2) printf("  block=%d fmin=%g fmax=%g\n", nxtblk, fmin, fmax);
  /* Write header */
  ret=ecat7WriteImageheader(fp, nxtblk, h); if(ret) {
    sprintf(ecat7errmsg, "cannot write subheader (%d).\n", ret);
    free(mdata); return(10);
  }
  /* Write matrix data */
  mptr=mdata;
  ret=ecat7WriteMatrixdata(fp, nxtblk+1, mptr, pxlNr, ecat7pxlbytes(h->data_type));
  free(mdata);
  if(ret) {
    sprintf(ecat7errmsg, "cannot write matrix data (%d).\n", ret);
    return(13);
  }
  return(0);
}
/*****************************************************************************/

/*****************************************************************************/
/*!
 * Write ECAT 7.x 2D sinogram matrix header and data
 *
 * @param fp		output file pointer
 * @param matrix_id     coded matrix id
 * @param h		Ecat7 2D image scan header
 * @param fdata		float data to be written
 * @returns 0 if ok.
 */
int ecat7Write2DScanMatrix(FILE *fp, int matrix_id, ECAT7_2Dscanheader *h, float *fdata) {
  int i, nxtblk, blkNr, data_size, pxlNr, ret;
  float *fptr, fmin, fmax, g, f;
  char *mdata, *mptr;
  short int *sptr;


  if(ECAT7_TEST) printf("ecat7Write2DScanMatrix(fp, %d, h, data)\n", matrix_id);
  if(fp==NULL || matrix_id<1 || h==NULL || fdata==NULL) {
    sprintf(ecat7errmsg, "invalid function parameter.\n");
    return(1);
  }
  if(h->data_type!=ECAT7_SUNI2) {
    sprintf(ecat7errmsg, "invalid data_type.\n");
    return(2);
  }
  /* nr of pixels */
  pxlNr=h->num_r_elements*h->num_angles;
  if(h->num_dimensions>2) pxlNr*=h->num_z_elements;
  if(pxlNr<1) {
    sprintf(ecat7errmsg, "invalid matrix dimension.\n");
    return(3);
  }
  /* How much memory is needed for ALL pixels */
  data_size=pxlNr*ecat7pxlbytes(h->data_type);
  /* block nr taken by all pixels */
  blkNr=(data_size+MatBLKSIZE-1)/MatBLKSIZE; if(blkNr<1) {
    sprintf(ecat7errmsg, "invalid block number.\n");
    return(4);
  }
  /* Allocate memory for matrix data */
  mdata=(char*)calloc(blkNr, MatBLKSIZE); if(mdata==NULL) {
    sprintf(ecat7errmsg, "out of memory.\n");
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
  if(f>=1.0 && ecat7_is_scaling_needed(g, fptr, pxlNr)==0) f=1.0;
  /* Scale matrix data to shorts */
  h->scale_factor=1.0/f;
  sptr=(short int*)mdata; fptr=fdata;
  for(i=0; i<pxlNr; i++, sptr++, fptr++) *sptr=(short int)temp_roundf(f*(*fptr));
  /* Set header short min & max */
  h->scan_min=(short int)temp_roundf(f*fmin);
  h->scan_max=(short int)temp_roundf(f*fmax);
  /* Get block number for matrix header and data */
  nxtblk=ecat7EnterMatrix(fp, matrix_id, blkNr); if(nxtblk<1) {
    sprintf(ecat7errmsg, "cannot determine matrix block (%d).\n", -nxtblk);
    free(mdata); return(8);
  }
  if(ECAT7_TEST>2) printf("  block=%d fmin=%g fmax=%g\n", nxtblk, fmin, fmax);
  /* Write header */
  ret=ecat7Write2DScanheader(fp, nxtblk, h); if(ret) {
    sprintf(ecat7errmsg, "cannot write subheader (%d).\n", ret);
    free(mdata); return(10);
  }
  /* Write matrix data */
  mptr=mdata;
  ret=ecat7WriteMatrixdata(fp, nxtblk+1, mptr, pxlNr, ecat7pxlbytes(h->data_type));
  free(mdata);
  if(ret) {
    sprintf(ecat7errmsg, "cannot write matrix data (%d).\n", ret);
    return(13);
  }
  return(0);
}
/*****************************************************************************/

/*****************************************************************************/
/*!
 * Write ECAT 7.x 3D sinogram matrix header and data
 *
 * @param fp		output file pointer
 * @param matrix_id	coded matrix id
 * @param h		Ecat7 scan header
 * @param fdata		float data
 * @return 0 if ok.
 */
int ecat7WriteScanMatrix(FILE *fp, int matrix_id, ECAT7_scanheader *h, float *fdata) {
  int i, nxtblk, blkNr, data_size, pxlNr, dimz, ret;
  float *fptr, fmin, fmax, g, f;
  char *mdata, *mptr;
  short int *sptr;


  if(ECAT7_TEST) printf("ecat7WriteScanMatrix(fp, %d, h, data)\n", matrix_id);
  if(fp==NULL || matrix_id<1 || h==NULL || fdata==NULL) {
    sprintf(ecat7errmsg, "invalid function parameter.\n");
    return(1);
  }
  if(h->data_type!=ECAT7_SUNI2) {
    sprintf(ecat7errmsg, "invalid data_type.\n");
    return(2);
  }
  /* nr of pixels */
  pxlNr=h->num_r_elements*h->num_angles;
  for(i=dimz=0; i<64; i++) dimz+=h->num_z_elements[i];
  pxlNr*=dimz;
  if(pxlNr<1) {
    sprintf(ecat7errmsg, "invalid matrix dimension.\n");
    return(3);
  }
  /* How much memory is needed for ALL pixels */
  data_size=pxlNr*ecat7pxlbytes(h->data_type);
  /* block nr taken by all pixels */
  blkNr=(data_size+MatBLKSIZE-1)/MatBLKSIZE; if(blkNr<1) {
    sprintf(ecat7errmsg, "invalid block number.\n");
    return(4);
  }
  /* Allocate memory for matrix data */
  mdata=(char*)calloc(blkNr, MatBLKSIZE); if(mdata==NULL) {
    sprintf(ecat7errmsg, "out of memory.\n");
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
  if(f>=1.0 && ecat7_is_scaling_needed(g, fptr, pxlNr)==0) f=1.0;
  /* Scale matrix data to shorts */
  h->scale_factor=1.0/f;
  sptr=(short int*)mdata; fptr=fdata;
  for(i=0; i<pxlNr; i++, sptr++, fptr++) *sptr=(short int)temp_roundf(f*(*fptr));
  /* Set header short min & max */
  h->scan_min=(short int)temp_roundf(f*fmin);
  h->scan_max=(short int)temp_roundf(f*fmax);
  /* Get block number for matrix header and data */
  /* Note that one extra block (blkNr+1) is needed for 3D scan header */
  nxtblk=ecat7EnterMatrix(fp, matrix_id, blkNr+1); if(nxtblk<1) {
    sprintf(ecat7errmsg, "cannot determine matrix block (%d).\n", -nxtblk);
    free(mdata); return(8);
  }
  if(ECAT7_TEST>2) printf("  block=%d fmin=%g fmax=%g\n", nxtblk, fmin, fmax);
  /* Write header */
  ret=ecat7WriteScanheader(fp, nxtblk, h); if(ret) {
    sprintf(ecat7errmsg, "cannot write subheader (%d).\n", ret);
    free(mdata); return(10);
  }
  /* Write matrix data */
  /* Note that 3D scan header takes TWO blocks */
  mptr=mdata;
  ret=ecat7WriteMatrixdata(fp, nxtblk+2, mptr, pxlNr, ecat7pxlbytes(h->data_type));
  free(mdata);
  if(ret) {
    sprintf(ecat7errmsg, "cannot write matrix data (%d).\n", ret);
    return(13);
  }
  return(0);
}
/*****************************************************************************/

/*****************************************************************************/
/*!
 * Write ECAT 7.x polarmap matrix header and data
 *
 * @param fp 		output file pointer
 * @param matrix_id	coded matrix information
 * @param h		Ecat7 polar map header
 * @param fdata		float data
 * @return 0 if ok.
 */
int ecat7WritePolarmapMatrix(FILE *fp, int matrix_id, ECAT7_polmapheader *h, float *fdata) {
  int i, nxtblk, blkNr, data_size, pxlNr, ret;
  float *fptr, fmin, fmax, g, f;
  char *mdata, *mptr;
  short int *sptr;
  
  
  if(ECAT7_TEST) printf("ecat7WritePolarmapMatrix(fp, %d, h, data)\n", matrix_id);
  if(fp==NULL || matrix_id<1 || h==NULL || fdata==NULL) {
    sprintf(ecat7errmsg, "invalid function parameter.\n");
    return(1);
  }
  if(h->data_type!=ECAT7_SUNI2) {
    sprintf(ecat7errmsg, "invalid data_type.\n");
    return(2);
  }
  /* nr of pixels */
  for(i=pxlNr=0; i<h->num_rings; i++) pxlNr+=h->sectors_per_ring[i];
  if(pxlNr<1) {
    sprintf(ecat7errmsg, "invalid matrix dimension.\n");
    return(3);
  }
  /* How much memory is needed for ALL pixels */
  data_size=pxlNr*ecat7pxlbytes(h->data_type);
  /* block nr taken by all pixels */
  blkNr=(data_size+MatBLKSIZE-1)/MatBLKSIZE; if(blkNr<1) {
    sprintf(ecat7errmsg, "invalid block number.\n");
    return(4);
  }
  /* Allocate memory for matrix data */
  mdata=(char*)calloc(blkNr, MatBLKSIZE); if(mdata==NULL) {
    sprintf(ecat7errmsg, "out of memory.\n");
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
  if(f>=1.0 && ecat7_is_scaling_needed(g, fptr, pxlNr)==0) f=1.0;
  /* Scale matrix data to shorts */
  h->scale_factor=1.0/f;
  sptr=(short int*)mdata; fptr=fdata;
  for(i=0; i<pxlNr; i++, sptr++, fptr++) *sptr=(short int)temp_roundf(f*(*fptr));
  /* Set header short min & max */
  h->image_min=(short int)temp_roundf(f*fmin);
  h->image_max=(short int)temp_roundf(f*fmax);
  /* Get block number for matrix header and data */
  nxtblk=ecat7EnterMatrix(fp, matrix_id, blkNr); if(nxtblk<1) {
    sprintf(ecat7errmsg, "cannot determine matrix block (%d).\n", -nxtblk);
    free(mdata); return(8);
  }
  if(ECAT7_TEST>2) printf("  block=%d fmin=%g fmax=%g\n", nxtblk, fmin, fmax);
  /* Write header */
  ret=ecat7WritePolmapheader(fp, nxtblk, h); if(ret) {
    sprintf(ecat7errmsg, "cannot write subheader (%d).\n", ret);
    free(mdata); return(10);
  }
  /* Write matrix data */
  mptr=mdata;
  ret=ecat7WriteMatrixdata(fp, nxtblk+1, mptr, pxlNr, ecat7pxlbytes(h->data_type));
  free(mdata);
  if(ret) {
    sprintf(ecat7errmsg, "cannot write matrix data (%d).\n", ret);
    return(13);
  }
  return(0);
}
/*****************************************************************************/

/*****************************************************************************/
/*!
 * Write ECAT 7.x matrix data to a specified file position.
 * Data does not need to be allocated for full blocks.
 * Data must be represented in current machines byte order, and it is
 * always saved in big endian byte order.
 *
 * @param fp           	Pointer to an opened ECAT file
 * @param start_block   Block number where matrix data is written
 * @param data          Pointer to matrix data
 * @param pxl_nr        Number of pixels
 * @param pxl_size      Size of data for one pixel in bytes
 * @return >0 in case of an error.
 */
int ecat7WriteMatrixdata(FILE *fp, int start_block, char *data, int pxl_nr, int pxl_size) {
  unsigned char buf[MatBLKSIZE];
  char *dptr;
  int i, blkNr, dataSize, byteNr, little;

  if(ECAT7_TEST) printf("ecat7WriteMatrixdata(fp, %d, data, %d, %d)\n",
                        start_block, pxl_nr, pxl_size);
  if(fp==NULL || start_block<1 || data==NULL || pxl_nr<1 || pxl_size<1) return(1);
  little=little_endian(); memset(buf, 0, MatBLKSIZE);
  dataSize=pxl_nr*pxl_size;
  /* block nr taken by all pixels */
  blkNr=(dataSize+MatBLKSIZE-1)/MatBLKSIZE; if(blkNr<1) return(1);
  if(ECAT7_TEST>2) printf("    blkNr=%d\n", blkNr);
  /* Search the place for writing */
  fseek(fp, (start_block-1)*MatBLKSIZE, SEEK_SET);
  if(ftell(fp)!=(start_block-1)*MatBLKSIZE) return(2);
  /* Save blocks one at a time */
  for(i=0, dptr=data; i<blkNr && dataSize>0; i++) {
    byteNr=(dataSize<MatBLKSIZE)?dataSize:MatBLKSIZE;
    memcpy(buf, dptr, byteNr);
    /* Change matrix byte order in little endian platforms */
    if(little) {
      if(pxl_size==2) swabip(buf, byteNr);
      else if(pxl_size==4) swawbip(buf, byteNr);
    }
    /* Write block */
    if(fwrite(buf, 1, MatBLKSIZE, fp)!=MatBLKSIZE) return(3);
    /* Prepare for the next block */
    dptr+=byteNr; dataSize-=byteNr;
  } /* next block */
  return(0);
}
/*****************************************************************************/

/*****************************************************************************/

