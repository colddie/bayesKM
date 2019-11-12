/// @file ecat7p.c
/// @author Vesa Oikonen
/// @brief Printing ECAT 7.x (header) contents.
///
/*****************************************************************************/
#include "libtpcimgio.h"
/*****************************************************************************/

/*****************************************************************************/
/*!
 * Print ECAT 7.x main header contents to specified file pointer
 *
 * @param h Ecat7 main header
 * @param fp target file pointer
 */
void ecat7PrintMainheader(ECAT7_mainheader *h, FILE *fp) {
  char tmp[64];
  time_t t;

  if(ECAT7_TEST) fprintf(stdout, "ecat7PrintMainheader()\n");
  fprintf(fp, "magic_number := %.14s\n", h->magic_number);
  fprintf(fp, "original_file_name := %.32s\n", h->original_file_name);
  fprintf(fp, "sw_version := %d\n", h->sw_version);
  fprintf(fp, "system_type := %d\n", h->system_type);
  fprintf(fp, "file_type := %d (%s)\n", h->file_type, ecat7filetype(h->file_type) );
  fprintf(fp, "serial_number := %.10s\n", h->serial_number);
  t=(time_t)h->scan_start_time;
  if(!ctime_r_int(&t, tmp)) strcpy(tmp, "1900-01-01 00:00:00");
  fprintf(fp, "scan_start_time := %s\n", tmp);
  fprintf(fp, "isotope_name := %.8s\n", h->isotope_name);
  fprintf(fp, "isotope_halflife := %E sec\n", h->isotope_halflife);
  fprintf(fp, "radiopharmaceutical := %.32s\n", h->radiopharmaceutical);
  fprintf(fp, "gantry_tilt := %g\n", h->gantry_tilt);
  fprintf(fp, "gantry_rotation := %g\n", h->gantry_rotation);
  fprintf(fp, "bed_elevation := %g\n", h->bed_elevation);
  fprintf(fp, "intrinsic_tilt := %g\n", h->intrinsic_tilt);
  fprintf(fp, "wobble_speed := %d\n", h->wobble_speed);
  fprintf(fp, "transm_source_type := %d\n", h->transm_source_type);
  fprintf(fp, "distance_scanned := %g\n", h->distance_scanned);
  fprintf(fp, "transaxial_fov := %g\n", h->transaxial_fov);
  fprintf(fp, "angular_compression := %d\n", h->angular_compression);
  fprintf(fp, "coin_samp_mode := %d\n", h->coin_samp_mode);
  fprintf(fp, "axial_samp_mode := %d\n", h->axial_samp_mode);
  fprintf(fp, "ecat_calibration_factor := %E\n", h->ecat_calibration_factor);
  fprintf(fp, "calibration_units := %d\n", h->calibration_units);
  fprintf(fp, "calibration_units_label := %d\n", h->calibration_units_label);
  fprintf(fp, "compression_code := %d\n", h->compression_code);
  fprintf(fp, "study_type := %.12s\n", h->study_type);
  fprintf(fp, "patient_id := %.16s\n", h->patient_id);
  fprintf(fp, "patient_name := %.32s\n", h->patient_name);
  fprintf(fp, "patient_sex := %c\n", (h->patient_sex!=0)?h->patient_sex:(char)32);
  fprintf(fp, "patient_dexterity := %c\n", (h->patient_dexterity!=0)?h->patient_dexterity:(char)32 );
  fprintf(fp, "patient_age := %g\n", h->patient_age);
  fprintf(fp, "patient_height := %g\n", h->patient_height);
  fprintf(fp, "patient_weight := %g\n", h->patient_weight);
  fprintf(fp, "patient_birth_date := %d\n", h->patient_birth_date);
  fprintf(fp, "physician_name := %.32s\n", h->physician_name);
  fprintf(fp, "operator_name := %.32s\n", h->operator_name);
  fprintf(fp, "study_description := %.32s\n", h->study_description);
  fprintf(fp, "acquisition_type := %d (%s)\n", h->acquisition_type,
          ecat7acquisitiontype(h->acquisition_type));
  fprintf(fp, "patient_orientation := %d\n", h->patient_orientation);
  fprintf(fp, "facility_name := %.20s\n", h->facility_name);
  fprintf(fp, "num_planes := %d\n", h->num_planes);
  fprintf(fp, "num_frames := %d\n", h->num_frames);
  fprintf(fp, "num_gates := %d\n", h->num_gates);
  fprintf(fp, "num_bed_pos := %d\n", h->num_bed_pos);
  fprintf(fp, "init_bed_position := %g\n", h->init_bed_position);
  fprintf(fp, "bed_position :=");
  for(int i=0; i<15; i++) fprintf(fp, " %g", h->bed_position[i]);
  fprintf(fp, "\n");
  fprintf(fp, "plane_separation := %g cm\n", h->plane_separation);
  fprintf(fp, "lwr_sctr_thres := %d\n", h->lwr_sctr_thres);
  fprintf(fp, "lwr_true_thres := %d\n", h->lwr_true_thres);
  fprintf(fp, "upr_true_thres := %d\n", h->upr_true_thres);
  fprintf(fp, "user_process_code := %.10s\n", h->user_process_code);
  fprintf(fp, "acquisition_mode := %d\n", h->acquisition_mode);
  fprintf(fp, "bin_size := %g cm\n", h->bin_size);
  fprintf(fp, "branching_fraction := %g\n", h->branching_fraction);
  t=(time_t)h->dose_start_time;
  if(!ctime_r_int(&t, tmp)) strcpy(tmp, "1900-01-01 00:00:00");
  fprintf(fp, "dose_start_time := %s\n", tmp);
  fprintf(fp, "dosage := %g\n", h->dosage);
  fprintf(fp, "well_counter_corr_factor := %E\n", h->well_counter_corr_factor);
  fprintf(fp, "data_units := %.32s\n", h->data_units);
  fprintf(fp, "septa_state := %d\n", h->septa_state);
  fprintf(fp, "fill_cti :=");
  for(int i=0; i<6; i++) fprintf(fp, " %d", h->fill_cti[i]);
  fprintf(fp, "\n");
}
/*****************************************************************************/

/*****************************************************************************/
/*!
 * Print ECAT 7.x image header contents to specified file pointer.
 *
 * @param h Ecat7 image header
 * @param fp target file pointer
 */
void ecat7PrintImageheader(ECAT7_imageheader *h, FILE *fp) 
{
  if(ECAT7_TEST) fprintf(stdout, "ecat7PrintImageheader()\n");
  fprintf(fp, "data_type := %d (%s)\n", h->data_type, ecat7datatype(h->data_type) );
  fprintf(fp, "num_dimensions := %d\n", h->num_dimensions);
  fprintf(fp, "x_dimension := %d\n", h->x_dimension);
  fprintf(fp, "y_dimension := %d\n", h->y_dimension);
  fprintf(fp, "z_dimension := %d\n", h->z_dimension);
  fprintf(fp, "x_offset := %g\n", h->x_offset);
  fprintf(fp, "y_offset := %g\n", h->y_offset);
  fprintf(fp, "z_offset := %g\n", h->z_offset);
  fprintf(fp, "recon_zoom := %g\n", h->recon_zoom);
  fprintf(fp, "scale_factor := %E\n", h->scale_factor);
  fprintf(fp, "image_min := %d\n", h->image_min);
  fprintf(fp, "image_max := %d\n", h->image_max);
  fprintf(fp, "x_pixel_size := %g\n", h->x_pixel_size);
  fprintf(fp, "y_pixel_size := %g\n", h->y_pixel_size);
  fprintf(fp, "z_pixel_size := %g\n", h->z_pixel_size);
  fprintf(fp, "frame_duration := %d\n", h->frame_duration);
  fprintf(fp, "frame_start_time := %d\n", h->frame_start_time);
  fprintf(fp, "filter_code := %d\n", h->filter_code);
  fprintf(fp, "x_resolution := %g\n", h->x_resolution);
  fprintf(fp, "y_resolution := %g\n", h->y_resolution);
  fprintf(fp, "z_resolution := %g\n", h->z_resolution);
  fprintf(fp, "num_r_elements := %g\n", h->num_r_elements);
  fprintf(fp, "num_angles := %g\n", h->num_angles);
  fprintf(fp, "z_rotation_angle := %g\n", h->z_rotation_angle);
  fprintf(fp, "decay_corr_fctr := %g\n", h->decay_corr_fctr);
  fprintf(fp, "processing_code := %d\n", h->processing_code);
  fprintf(fp, "gate_duration := %d\n", h->gate_duration);
  fprintf(fp, "r_wave_offset := %d\n", h->r_wave_offset);
  fprintf(fp, "num_accepted_beats := %d\n", h->num_accepted_beats);
  fprintf(fp, "filter_cutoff_frequency := %E\n", h->filter_cutoff_frequency);
  fprintf(fp, "filter_resolution := %E\n", h->filter_resolution);
  fprintf(fp, "filter_ramp_slope := %E\n", h->filter_ramp_slope);
  fprintf(fp, "filter_order := %d\n", h->filter_order);
  fprintf(fp, "filter_scatter_fraction := %E\n", h->filter_scatter_fraction);
  fprintf(fp, "filter_scatter_slope := %E\n", h->filter_scatter_slope);
  fprintf(fp, "annotation := %.40s\n", h->annotation);
  fprintf(fp, "mt_1_1 := %g\n", h->mt_1_1);
  fprintf(fp, "mt_1_2 := %g\n", h->mt_1_2);
  fprintf(fp, "mt_1_3 := %g\n", h->mt_1_3);
  fprintf(fp, "mt_2_1 := %g\n", h->mt_2_1);
  fprintf(fp, "mt_2_2 := %g\n", h->mt_2_2);
  fprintf(fp, "mt_2_3 := %g\n", h->mt_2_3);
  fprintf(fp, "mt_3_1 := %g\n", h->mt_3_1);
  fprintf(fp, "mt_3_2 := %g\n", h->mt_3_2);
  fprintf(fp, "mt_3_3 := %g\n", h->mt_3_3);
  fprintf(fp, "rfilter_cutoff := %g\n", h->rfilter_cutoff);
  fprintf(fp, "rfilter_resolution := %g\n", h->rfilter_resolution);
  fprintf(fp, "rfilter_code := %d\n", h->rfilter_code);
  fprintf(fp, "rfilter_order := %d\n", h->rfilter_order);
  fprintf(fp, "zfilter_cutoff := %g\n", h->zfilter_cutoff);
  fprintf(fp, "zfilter_resolution := %g\n", h->zfilter_resolution);
  fprintf(fp, "zfilter_code := %d\n", h->zfilter_code);
  fprintf(fp, "zfilter_order := %d\n", h->zfilter_order);
  fprintf(fp, "mt_1_4 := %g\n", h->mt_1_4);
  fprintf(fp, "mt_2_4 := %g\n", h->mt_2_4);
  fprintf(fp, "mt_3_4 := %g\n", h->mt_3_4);
  fprintf(fp, "scatter_type := %d\n", h->scatter_type);
  fprintf(fp, "recon_type := %d\n", h->recon_type);
  fprintf(fp, "recon_views := %d\n", h->recon_views);
  fprintf(fp, "fill_cti :=");
  for(int i=0; i<87; i++) fprintf(fp, " %d", h->fill_cti[i]);
  fprintf(fp, "\n");
  fprintf(fp, "fill_user :=");
  for(int i=0; i<49; i++) fprintf(fp, " %d", h->fill_user[i]);
  fprintf(fp, "\n");
}
/*****************************************************************************/

/*****************************************************************************/
/*!
 * Print ECAT 7.x 3D sinogram header contents to specified file pointer
 *
 * @param h Ecat7 scan header
 * @param fp target file pointer
 */
void ecat7PrintScanheader(ECAT7_scanheader *h, FILE *fp) 
{
  if(ECAT7_TEST) fprintf(stdout, "ecat7PrintScanheader()\n");
  fprintf(fp, "data_type := %d (%s)\n", h->data_type,
    ecat7datatype(h->data_type) );
  fprintf(fp, "num_dimensions := %d\n", h->num_dimensions);
  fprintf(fp, "num_r_elements := %d\n", h->num_r_elements);
  fprintf(fp, "num_angles := %d\n", h->num_angles);
  fprintf(fp, "corrections_applied := %d\n", h->corrections_applied);
  fprintf(fp, "num_z_elements :=");
  for(int i=0; i<64; i++) fprintf(fp, " %d", h->num_z_elements[i]);
  fprintf(fp, "\n");
  fprintf(fp, "ring_difference := %d\n", h->ring_difference);
  fprintf(fp, "storage_order := %d\n", h->storage_order);
  fprintf(fp, "axial_compression := %d (span)\n", h->axial_compression);
  fprintf(fp, "x_resolution := %g cm\n", h->x_resolution);
  fprintf(fp, "v_resolution := %g rad\n", h->v_resolution);
  fprintf(fp, "z_resolution := %g cm\n", h->z_resolution);
  fprintf(fp, "w_resolution := %g\n", h->w_resolution);
  fprintf(fp, "gate_duration := %d\n", h->gate_duration);
  fprintf(fp, "r_wave_offset := %d\n", h->r_wave_offset);
  fprintf(fp, "num_accepted_beats := %d\n", h->num_accepted_beats);
  fprintf(fp, "scale_factor := %E\n", h->scale_factor);
  fprintf(fp, "scan_min := %d\n", h->scan_min);
  fprintf(fp, "scan_max := %d\n", h->scan_max);
  fprintf(fp, "prompts := %d\n", h->prompts);
  fprintf(fp, "delayed := %d\n", h->delayed);
  fprintf(fp, "multiples := %d\n", h->multiples);
  fprintf(fp, "net_trues := %d\n", h->net_trues);
  fprintf(fp, "tot_avg_cor := %g\n", h->tot_avg_cor);
  fprintf(fp, "tot_avg_uncor := %g\n", h->tot_avg_uncor);
  fprintf(fp, "total_coin_rate := %d\n", h->total_coin_rate);
  fprintf(fp, "frame_start_time := %d\n", h->frame_start_time);
  fprintf(fp, "frame_duration := %d\n", h->frame_duration);
  fprintf(fp, "deadtime_correction_factor := %g\n", h->deadtime_correction_factor);
  fprintf(fp, "uncor_singles :=");
  for(int i=0; i<128; i++) fprintf(fp, " %g", h->uncor_singles[i]);
  fprintf(fp, "\n");
}
/*****************************************************************************/

/*****************************************************************************/
/*!
 * Print ECAT 7.x attenuation header contents to specified file pointer
 *
 * @param h Ecat7 attenuation header
 * @param fp target file pointer
 */
void ecat7PrintAttenheader(ECAT7_attenheader *h, FILE *fp) 
{
  if(ECAT7_TEST) fprintf(stdout, "ecat7PrintAttenheader()\n");
  fprintf(fp, "data_type := %d (%s)\n", h->data_type,
    ecat7datatype(h->data_type) );
  fprintf(fp, "num_dimensions := %d\n", h->num_dimensions);
  fprintf(fp, "attenuation_type := %d\n", h->attenuation_type);
  fprintf(fp, "num_r_elements := %d\n", h->num_r_elements);
  fprintf(fp, "num_angles := %d\n", h->num_angles);
  fprintf(fp, "num_z_elements := %d\n", h->num_z_elements);
  fprintf(fp, "ring_difference := %d\n", h->ring_difference);
  fprintf(fp, "x_resolution := %g\n", h->x_resolution);
  fprintf(fp, "y_resolution := %g\n", h->y_resolution);
  fprintf(fp, "z_resolution := %g\n", h->z_resolution);
  fprintf(fp, "w_resolution := %g\n", h->w_resolution);
  fprintf(fp, "scale_factor := %E\n", h->scale_factor);
  fprintf(fp, "x_offset := %g\n", h->x_offset);
  fprintf(fp, "y_offset := %g\n", h->y_offset);
  fprintf(fp, "x_radius := %g\n", h->x_radius);
  fprintf(fp, "y_radius := %g\n", h->y_radius);
  fprintf(fp, "tilt_angle := %g\n", h->tilt_angle);
  fprintf(fp, "attenuation_coeff := %E\n", h->attenuation_coeff);
  fprintf(fp, "attenuation_min := %E\n", h->attenuation_min);
  fprintf(fp, "attenuation_max := %E\n", h->attenuation_max);
  fprintf(fp, "skull_thickness := %g\n", h->skull_thickness);
  fprintf(fp, "num_additional_atten_coeff := %d\n", h->num_additional_atten_coeff);
  fprintf(fp, "additional_atten_coeff :=");
  for(int i=0; i<8; i++) fprintf(fp, " %E", h->additional_atten_coeff[i]);
  fprintf(fp, "\n");
  fprintf(fp, "edge_finding_threshold := %g\n", h->edge_finding_threshold);
  fprintf(fp, "storage_order := %d\n", h->storage_order);
  fprintf(fp, "span := %d\n", h->span);
  fprintf(fp, "z_elements :=");
  for(int i=0; i<64; i++) fprintf(fp, " %d", h->z_elements[i]);
  fprintf(fp, "\n");
}
/*****************************************************************************/

/*****************************************************************************/
/*!
 * Print ECAT 7.x polar map header contents to specified file pointer
 *
 * @param h Ecat7 polar map header
 * @param fp target file pointer
 */
void ecat7PrintPolmapheader(ECAT7_polmapheader *h, FILE *fp) {
  int i;

  if(ECAT7_TEST) fprintf(stdout, "ecat7PrintPolmapheader()\n");
  fprintf(fp, "data_type := %d (%s)\n", h->data_type,
    ecat7datatype(h->data_type) );
  fprintf(fp, "polar_map_type := %d\n", h->polar_map_type);
  fprintf(fp, "num_rings := %d\n", h->num_rings);
  fprintf(fp, "sectors_per_ring :=");
  for(i=0; i<32; i++) fprintf(fp, " %d", h->sectors_per_ring[i]);
  fprintf(fp, "\n");
  fprintf(fp, "ring_position :=");
  for(i=0; i<32; i++) fprintf(fp, " %g", h->ring_position[i]);
  fprintf(fp, "\n");
  fprintf(fp, "ring_angle :=");
  for(i=0; i<32; i++) fprintf(fp, " %d", h->ring_angle[i]);
  fprintf(fp, "\n");
  fprintf(fp, "start_angle := %d\n", h->start_angle);
  fprintf(fp, "long_axis_left :=");
  for(i=0; i<3; i++) fprintf(fp, " %d", h->long_axis_left[i]);
  fprintf(fp, "\n");
  fprintf(fp, "long_axis_right :=");
  for(i=0; i<3; i++) fprintf(fp, " %d", h->long_axis_right[i]);
  fprintf(fp, "\n");
  fprintf(fp, "position_data := %d\n", h->position_data);
  fprintf(fp, "image_min := %d\n", h->image_min);
  fprintf(fp, "image_max := %d\n", h->image_max);
  fprintf(fp, "scale_factor := %E\n", h->scale_factor);
  fprintf(fp, "pixel_size := %g\n", h->pixel_size);
  fprintf(fp, "frame_duration := %d\n", h->frame_duration);
  fprintf(fp, "frame_start_time := %d\n", h->frame_start_time);
  fprintf(fp, "processing_code := %d\n", h->processing_code);
  fprintf(fp, "quant_units := %d\n", h->quant_units);
  fprintf(fp, "annotation := %.40s\n", h->annotation);
  fprintf(fp, "gate_duration := %d\n", h->gate_duration);
  fprintf(fp, "r_wave_offset := %d\n", h->r_wave_offset);
  fprintf(fp, "num_accepted_beats := %d\n", h->num_accepted_beats);
  fprintf(fp, "polar_map_protocol := %.20s\n", h->polar_map_protocol);
  fprintf(fp, "database_name := %.30s\n", h->database_name);
}
/*****************************************************************************/

/*****************************************************************************/
/*!
 * Prints ECAT 7.x normalization header contents to specified file pointer
 *
 * @param h Ecat7 normalization header
 * @param fp tager file pointer
 */
void ecat7PrintNormheader(ECAT7_normheader *h, FILE *fp) 
{
  int i;

  if(ECAT7_TEST) fprintf(stdout, "ecat7PrintNormheader()\n");
  fprintf(fp, "data_type := %d (%s)\n", h->data_type,
    ecat7datatype(h->data_type) );
  fprintf(fp, "num_r_elements := %d\n", h->num_r_elements);
  fprintf(fp, "num_transaxial_crystals := %d\n", h->num_transaxial_crystals);
  fprintf(fp, "num_crystal_rings := %d\n", h->num_crystal_rings);
  fprintf(fp, "crystals_per_ring := %d\n", h->crystals_per_ring);
  fprintf(fp, "num_geo_corr_planes := %d\n", h->num_geo_corr_planes);
  fprintf(fp, "uld := %d\n", h->uld);
  fprintf(fp, "lld := %d\n", h->lld);
  fprintf(fp, "scatter_energy := %d\n", h->scatter_energy);
  fprintf(fp, "norm_quality_factor := %g\n", h->norm_quality_factor);
  fprintf(fp, "norm_quality_factor_code := %d\n", h->norm_quality_factor_code);
  fprintf(fp, "ring_dtcor1 :=");
  for(i=0; i<32; i++) fprintf(fp, " %E", h->ring_dtcor1[i]);
  fprintf(fp, "\n");
  fprintf(fp, "ring_dtcor2 :=");
  for(i=0; i<32; i++) fprintf(fp, " %E", h->ring_dtcor2[i]);
  fprintf(fp, "\n");
  fprintf(fp, "crystal_dtcor :=");
  for(i=0; i<8; i++) fprintf(fp, " %E", h->crystal_dtcor[i]);
  fprintf(fp, "\n");
  fprintf(fp, "span := %d\n", h->span);
  fprintf(fp, "max_ring_diff := %d\n", h->max_ring_diff);
}
/*****************************************************************************/

/*****************************************************************************/
/*!
 * Print ECAT 7.x 2D sinogram header contents to specified file pointer
 *
 * @param h Ecat7 2D scan header
 * @param fp target file pointer
 */
void ecat7Print2DScanheader(ECAT7_2Dscanheader *h, FILE *fp) 
{
  int i;

  if(ECAT7_TEST) fprintf(stdout, "ecat7Print2DScanheader()\n");
  fprintf(fp, "data_type := %d (%s)\n", h->data_type,
    ecat7datatype(h->data_type) );
  fprintf(fp, "num_dimensions := %d\n", h->num_dimensions);
  fprintf(fp, "num_r_elements := %d\n", h->num_r_elements);
  fprintf(fp, "num_angles := %d\n", h->num_angles);
  fprintf(fp, "corrections_applied := %d\n", h->corrections_applied);
  fprintf(fp, "num_z_elements := %d\n", h->num_z_elements);
  fprintf(fp, "ring_difference := %d\n", h->ring_difference);
  fprintf(fp, "x_resolution := %g\n", h->x_resolution);
  fprintf(fp, "y_resolution := %g\n", h->y_resolution);
  fprintf(fp, "z_resolution := %g\n", h->z_resolution);
  fprintf(fp, "w_resolution := %g\n", h->w_resolution);
  fprintf(fp, "gate_duration := %d\n", h->gate_duration);
  fprintf(fp, "r_wave_offset := %d\n", h->r_wave_offset);
  fprintf(fp, "num_accepted_beats := %d\n", h->num_accepted_beats);
  fprintf(fp, "scale_factor := %E\n", h->scale_factor);
  fprintf(fp, "scan_min := %d\n", h->scan_min);
  fprintf(fp, "scan_max := %d\n", h->scan_max);
  fprintf(fp, "prompts := %d\n", h->prompts);
  fprintf(fp, "delayed := %d\n", h->delayed);
  fprintf(fp, "multiples := %d\n", h->multiples);
  fprintf(fp, "net_trues := %d\n", h->net_trues);
  fprintf(fp, "cor_singles :=");
  for(i=0; i<16; i++) fprintf(fp, " %g", h->cor_singles[i]);
  fprintf(fp, "\n");
  fprintf(fp, "uncor_singles :=");
  for(i=0; i<16; i++) fprintf(fp, " %g", h->uncor_singles[i]);
  fprintf(fp, "\n");
  fprintf(fp, "tot_avg_cor := %g\n", h->tot_avg_cor);
  fprintf(fp, "tot_avg_uncor := %g\n", h->tot_avg_uncor);
  fprintf(fp, "total_coin_rate := %d\n", h->total_coin_rate);
  fprintf(fp, "frame_start_time := %d\n", h->frame_start_time);
  fprintf(fp, "frame_duration := %d\n", h->frame_duration);
  fprintf(fp, "deadtime_correction_factor := %E\n", h->deadtime_correction_factor);
  fprintf(fp, "physical_planes :=");
  for(i=0; i<8; i++) fprintf(fp, " %d", h->physical_planes[i]);
  fprintf(fp, "\n");
}
/*****************************************************************************/

/*****************************************************************************/
/*!
 * Print ECAT 7.x 2D normalization header contents to specified file pointer
 *
 * @param h Ecat7 2D normalization header
 * @param fp target file pointer
 */
void ecat7Print2DNormheader(ECAT7_2Dnormheader *h, FILE *fp) 
{
  int i;

  if(ECAT7_TEST) fprintf(stdout, "ecat7Print2DNormheader()\n");
  fprintf(fp, "data_type := %d (%s)\n", h->data_type,
    ecat7datatype(h->data_type) );
  fprintf(fp, "num_dimensions := %d\n", h->num_dimensions);
  fprintf(fp, "num_r_elements := %d\n", h->num_r_elements);
  fprintf(fp, "num_angles := %d\n", h->num_angles);
  fprintf(fp, "num_z_elements := %d\n", h->num_z_elements);
  fprintf(fp, "ring_difference := %d\n", h->ring_difference);
  fprintf(fp, "scale_factor := %E\n", h->scale_factor);
  fprintf(fp, "norm_min := %g\n", h->norm_min);
  fprintf(fp, "norm_max := %g\n", h->norm_max);
  fprintf(fp, "fov_source_width := %g\n", h->fov_source_width);
  fprintf(fp, "norm_quality_factor := %g\n", h->norm_quality_factor);
  fprintf(fp, "norm_quality_factor_code := %d\n", h->norm_quality_factor_code);
  fprintf(fp, "storage_order := %d\n", h->storage_order);
  fprintf(fp, "span := %d\n", h->span);
  fprintf(fp, "z_elements :=");
  for(i=0; i<64; i++) fprintf(fp, " %d", h->z_elements[i]);
  fprintf(fp, "\n");
}
/*****************************************************************************/

/*****************************************************************************/
/*!
 * Returns pointer to a string describing the ECAT7 file_type
 *
 * @param file_type file type code
 * @return pointer to static string
 */
char* ecat7filetype(short int file_type) {
  static char *info[] = {
  "unknown", "2D sinogram", "image-16", "attenuation correction",
  "2D normalization", "polar map", "volume 8", "volume 16",
  "projection 8", "projection 16", "image 8", "3D sinogram 16",
  "3D sinogram 8", "3D normalization", "3D sinogram fit",
  0};
  if(file_type>=0 && file_type<=14) return((char*)info[file_type]);
  else return((char*)info[0]);
}
/*****************************************************************************/

/*****************************************************************************/
/*!
 * Returns pointer to a string describing the ECAT7 acquisition_type
 *
 * @param acquisition_type acquisition type code
 * @return pointer to static string
 */
char* ecat7acquisitiontype(short int acquisition_type) {
  static char *info[] = {
  "undefined", "blank", "transmission", "static emission",
  "dynamic emission", "gated emission", "transmission rectilinear",
  "emission rectilinear",
  0};
  if(acquisition_type>=0 && acquisition_type<=7)
    return((char*)info[acquisition_type]);
  else return((char*)info[0]);
}
/*****************************************************************************/

/*****************************************************************************/
/*!
 * Returns pointer to a string describing the ECAT7 data_type
 *
 * @param data_type data type code
 * @return pointer to static string
 */
char* ecat7datatype(short int data_type) {
  static char *info[] = {
  "unknown", "byte", "VAX 2 byte integer", "VAX 4 byte integer",
  "VAX 4 byte float", "IEEE 4 byte float", "SUN 2 byte integer",
  "SUN 4 byte integer",
  0};
  if(data_type>=0 && data_type<=7) return((char*)info[data_type]);
  else return((char*)info[0]);
}
/*****************************************************************************/

/*****************************************************************************/
/** Print ECAT7 subheader contents into specified file pointer.
\return Returns 0 when successful.
 */
int ecat7PrintSubheader(
  /** ECAT7 mainheader (not printed but needed here) */
  ECAT7_mainheader mh,
  /** File pointer to ECAT7 file */
  FILE *fp,
  /** ECAT7 plane */
  int plane,
  /** ECAT7 frame */
  int frame,
  /** Output is written to this file pointer; it can be stdout */
  FILE *ofp
) {
  int                     mi, ret, nr=0;
  ECAT7_imageheader       image_header;
  ECAT7_scanheader        scan_header;
  ECAT7_2Dscanheader      scan2D_header;
  ECAT7_2Dnormheader      norm2D_header;
  ECAT7_normheader        norm_header;
  ECAT7_attenheader       atten_header;
  ECAT7_polmapheader      polmap_header;
  static ECAT7_MATRIXLIST mlist;
  ECAT7_Matval            matval;


  /*
   *  Read matrix list
   */
  ecat7InitMatlist(&mlist);
  ret=ecat7ReadMatlist(fp, &mlist, ECAT7_TEST);
  if(ret) {
    fprintf(stderr, "Error (%d): cannot read matrix list.\n", ret);
    return(2);
  }
  if(mlist.matrixNr<=0) {
    fprintf(stderr, "Error: matrix list is empty.\n");
    return(2);
  }
  if(ECAT7_TEST>1) ecat7PrintMatlist(&mlist);

  /*
   *  Read and print subheaders one at a time
   */
  for(mi=nr=0; mi<mlist.matrixNr; mi++) {
    /* Get frame nr */
    ecat7_id_to_val(mlist.matdir[mi].id, &matval);
    /* Check if this is supposed to be listed or not */
    if(frame>=0 && frame!=matval.frame) continue;
    if(plane>=0 && plane!=matval.plane) continue;
    fprintf(fp, "Matrix: plane %d frame %d gate %d bed %d\n",
      matval.plane, matval.frame, matval.gate, matval.bed);
    /* Read and print subheader */
    ret=0;
    switch(mh.file_type) {
      case ECAT7_ATTEN:
        ret=ecat7ReadAttenheader(fp, mlist.matdir[mi].strtblk, &atten_header);
        if(ret==0) ecat7PrintAttenheader(&atten_header, ofp);
        break;
      case ECAT7_3DNORM:
        ret=ecat7ReadNormheader(fp, mlist.matdir[mi].strtblk, &norm_header);
        if(ret==0) ecat7PrintNormheader(&norm_header, ofp);
        break;
      case ECAT7_IMAGE8:
      case ECAT7_IMAGE16:
      case ECAT7_VOLUME8:
      case ECAT7_VOLUME16:
        ret=ecat7ReadImageheader(fp, mlist.matdir[mi].strtblk, &image_header);
        if(ret==0) ecat7PrintImageheader(&image_header, ofp);
        break;
      case ECAT7_3DSCAN:
      case ECAT7_3DSCAN8:
      case ECAT7_3DSCANFIT:
        ret=ecat7ReadScanheader(fp, mlist.matdir[mi].strtblk, &scan_header);
        if(ret==0) ecat7PrintScanheader(&scan_header, ofp);
        break;
      case ECAT7_POLARMAP:
        ret=ecat7ReadPolmapheader(fp, mlist.matdir[mi].strtblk, &polmap_header);
        if(ret==0) ecat7PrintPolmapheader(&polmap_header, ofp);
        break;
      case ECAT7_2DSCAN:
        ret=ecat7Read2DScanheader(fp, mlist.matdir[mi].strtblk, &scan2D_header);
        if(ret==0) ecat7Print2DScanheader(&scan2D_header, ofp);
        break;
      case ECAT7_2DNORM:
        ret=ecat7Read2DNormheader(fp, mlist.matdir[mi].strtblk, &norm2D_header);
        if(ret==0) ecat7Print2DNormheader(&norm2D_header, ofp);
        break;
      default:
        fprintf(stderr, "Error: matrix filetype %d is not yet supported.\n",
          mh.file_type);
        ecat7EmptyMatlist(&mlist);
        return(8);
    }
    if(ret) {
      fprintf(stderr, "Error %d in reading subheader.\n", ret);
      ecat7EmptyMatlist(&mlist); return(5);
    }
    nr++; // counter
  } /* next matrix */
  ecat7EmptyMatlist(&mlist);
  
  if(nr==0 && (plane>=0 || frame>=0)) {
    fprintf(stderr, "Error: specified matrices not found.\n");
    return(11);
  }

  return(0);
}
/*****************************************************************************/

/*****************************************************************************/

