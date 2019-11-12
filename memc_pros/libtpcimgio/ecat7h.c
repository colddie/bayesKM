/// @file ecat7h.c
/// @author Vesa Oikonen
/// @brief Procedures for editing ECAT 7.x header contents.
///
/*****************************************************************************/
#include "libtpcimgio.h"
/*****************************************************************************/

/*****************************************************************************/
/** Edit ECAT 7 main header.
\return Returns 0, if ok, and 1 or 2, if field name or or value is invalid.
 */
int ecat7EditMHeader(
  /** Pointer to ECAT 7 mainheader structure */
  ECAT7_mainheader *h,
  /** Field name to be changed */
  char *field,
  /** New value for the field */
  char *value,
  /** Verbose level; if <=0, then nothing is printed into stdout */
  int verbose
) {
  int yy, mm, dd;
  short int si;
  float f;

  if(verbose>0) printf("ecat7EditMHeader('%s', '%s')\n", field, value);
  si=atoi(value); f=atof(value);
  if(strcmp(field, "magic_number")==0) {
    strlcpy(h->magic_number, value, 14);
  } else if(strcmp(field, "original_file_name")==0) {
    strlcpy(h->original_file_name, value, 32);
  } else if(strcmp(field, "sw_version")==0) {
    if(si<=0) return(2); else h->sw_version=si;
  } else if(strcmp(field, "system_type")==0) {
    if(si<0) return(2); else h->system_type=si;
  } else if(strcmp(field, "file_type")==0) {
    if(si<0) return(2); else h->file_type=si;
  } else if(strcmp(field, "serial_number")==0) {
    strlcpy(h->serial_number, value, 10);
  } else if(strcmp(field, "scan_start_time")==0) {
    struct tm stm;
    time_t t;
    if(get_datetime(value, &stm, verbose-1)!=0) return(2);
    if(verbose>3) printf("  year=%d\n", stm.tm_year);
    if(verbose>1) printf("  hour=%d\n", stm.tm_hour);
    /* ECAT 7 main header saves int, not time_t (long int),
       therefore checking whether time can be saved correctly */
    t=timegm(&stm);
    h->scan_start_time=timegm(&stm);
    if(t!=(time_t)h->scan_start_time) {
      /* no it can't be saved as required; set negative time as time base */
      if(t<(time_t)0) h->scan_start_time=0;
      else return(2); /* if too far in future, return error */
    }
    if(verbose>1) printf("  scan_start_time := %d\n", h->scan_start_time);
    if(h->scan_start_time==-1) h->scan_start_time=0; //return(2);
  } else if(strcmp(field, "isotope_name")==0) {
    strlcpy(h->isotope_name, value, 8);
  } else if(strcmp(field, "isotope_halflife")==0) {
    if(f<=1.0E-3) return(2); else h->isotope_halflife=f;
  } else if(strcmp(field, "radiopharmaceutical")==0) {
    strlcpy(h->radiopharmaceutical, value, 32);
  } else if(strcmp(field, "gantry_tilt")==0) {
    h->gantry_tilt=f;
  } else if(strcmp(field, "gantry_rotation")==0) {
    h->gantry_rotation=f;
  } else if(strcmp(field, "bed_elevation")==0) {
    h->bed_elevation=f;
  } else if(strcmp(field, "intrinsic_tilt")==0) {
    h->intrinsic_tilt=f;
  } else if(strcmp(field, "wobble_speed")==0) {
    h->wobble_speed=si;
  } else if(strcmp(field, "transm_source_type")==0) {
    h->transm_source_type=si;
  } else if(strcmp(field, "distance_scanned")==0) {
    h->distance_scanned=f;
  } else if(strcmp(field, "transaxial_fov")==0) {
    h->transaxial_fov=f;
  } else if(strcmp(field, "angular_compression")==0) {
    h->angular_compression=si;
  } else if(strcmp(field, "coin_samp_mode")==0) {
    h->coin_samp_mode=si;
  } else if(strcmp(field, "axial_samp_mode")==0) {
    h->axial_samp_mode=si;
  } else if(strcmp(field, "ecat_calibration_factor")==0) {
    h->ecat_calibration_factor=f;
  } else if(strcmp(field, "calibration_units")==0) {
    h->calibration_units=si;
  } else if(strcmp(field, "calibration_units_label")==0) {
    h->calibration_units_label=si;
  } else if(strcmp(field, "compression_code")==0) {
    h->compression_code=si;
  } else if(strcmp(field, "study_type")==0) {
    strlcpy(h->study_type, value, 12);
  } else if(strcmp(field, "patient_id")==0) {
    strlcpy(h->patient_id, value, 16);
  } else if(strcmp(field, "patient_name")==0) {
    strlcpy(h->patient_name, value, 32);
  } else if(strcmp(field, "patient_sex")==0) {
    h->patient_sex=value[0];
  } else if(strcmp(field, "patient_dexterity")==0) {
    h->patient_dexterity=value[0];
  } else if(strcmp(field, "patient_age")==0) {
    h->patient_age=f;
  } else if(strcmp(field, "patient_height")==0) {
    h->patient_height=f;
  } else if(strcmp(field, "patient_weight")==0) {
    h->patient_weight=f;
  } else if(strcmp(field, "patient_birth_date")==0) {
    struct tm st;
    time_t timet;
    timet=time(NULL); gmtime_r(&timet, &st);
    if(sscanf(value, "%d-%d-%d", &yy, &mm, &dd)!=3) return(2);
    st.tm_mday=dd; st.tm_mon=mm-1; st.tm_year=yy-1900;
    st.tm_hour=12; st.tm_min=0; st.tm_sec=0; st.tm_isdst=-1;
    h->patient_birth_date=timegm(&st);
  } else if(strcmp(field, "physician_name")==0) {
    strlcpy(h->physician_name, value, 32);
  } else if(strcmp(field, "operator_name")==0) {
    strlcpy(h->operator_name, value, 32);
  } else if(strcmp(field, "study_description")==0) {
    strlcpy(h->study_description, value, 32);
  } else if(strcmp(field, "acquisition_type")==0) {
    h->acquisition_type=si;
  } else if(strcmp(field, "patient_orientation")==0) {
    h->patient_orientation=si;
  } else if(strcmp(field, "facility_name")==0) {
    strlcpy(h->facility_name, value, 20);
  } else if(strcmp(field, "num_planes")==0) {
    h->num_planes=si;
  } else if(strcmp(field, "num_frames")==0) {
    h->num_frames=si;
  } else if(strcmp(field, "num_gates")==0) {
    h->num_gates=si;
  } else if(strcmp(field, "num_bed_pos")==0) {
    h->num_bed_pos=si;
  } else if(strcmp(field, "init_bed_position")==0) {
    h->init_bed_position=f;
  } else if(strcmp(field, "bed_position")==0) {
    sscanf(value, "%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f",
      h->bed_position+0, h->bed_position+1, h->bed_position+2,
      h->bed_position+3, h->bed_position+4, h->bed_position+5,
      h->bed_position+6, h->bed_position+7, h->bed_position+8,
      h->bed_position+9, h->bed_position+10, h->bed_position+11,
      h->bed_position+12, h->bed_position+13, h->bed_position+14
    );
  } else if(strcmp(field, "plane_separation")==0) {
    h->plane_separation=f;
  } else if(strcmp(field, "lwr_sctr_thres")==0) {
    h->lwr_sctr_thres=si;
  } else if(strcmp(field, "lwr_true_thres")==0) {
    h->lwr_true_thres=si;
  } else if(strcmp(field, "upr_true_thres")==0) {
    h->upr_true_thres=si;
  } else if(strcmp(field, "user_process_code")==0) {
    strlcpy(h->user_process_code, value, 10);
  } else if(strcmp(field, "acquisition_mode")==0) {
    h->acquisition_mode=si;
  } else if(strcmp(field, "bin_size")==0) {
    h->bin_size=f;
  } else if(strcmp(field, "branching_fraction")==0) {
    h->branching_fraction=f;
  } else if(strcmp(field, "dose_start_time")==0) {
    struct tm stm;
    if(get_datetime(value, &stm, verbose-1)!=0) return(2);
    h->dose_start_time=timegm(&stm);
  } else if(strcmp(field, "dosage")==0) {
    h->dosage=f;
  } else if(strcmp(field, "well_counter_corr_factor")==0) {
    h->well_counter_corr_factor=f;
  } else if(strcmp(field, "data_units")==0) {
    strlcpy(h->data_units, value, 32);
  } else if(strcmp(field, "septa_state")==0) {
    h->septa_state=si;
  } else if(strncasecmp(field, "FILL_CTI", 8)==0) {
    char *cptr;
    cptr=strtok(value, " \t,;\n\r");
    for(int i=0; i<6; i++) {
      if(cptr==NULL) break;
      h->fill_cti[i]=(short int)atoi(cptr);
      cptr=strtok(NULL, " \t,;\n\r");
    }
  } else
    return(1);

  return(0);
}
/*****************************************************************************/

/*****************************************************************************/
/** Edit ECAT 7 3Dscan header.
\return Returns 0, if ok, and 1 or 2, if field name or
    or value is invalid.
 */
int ecat7EditSHeader(
  /** Pointer to ECAT 7 3D scan header structure */
  ECAT7_scanheader *h,
  /** Field name to be changed */
  char *field,
  /** New value for the field */
  char *value,
  /** Verbose level; if <=0, then nothing is printed into stdout */
  int verbose
) {
  int i, ii;
  short int si;
  float f;
  char *cptr;

  if(verbose>0) printf("ecat7EditSHeader('%s', '%s')\n", field, value);

  si=atoi(value); ii=atoi(value); f=atof(value);

  if(strcasecmp(field, "DATA_TYPE")==0) {
    h->data_type=si;
  } else if(strcasecmp(field, "NUM_DIMENSIONS")==0) {
    h->num_dimensions=si;
  } else if(strcasecmp(field, "NUM_R_ELEMENTS")==0) {
    h->num_r_elements=si;
  } else if(strcasecmp(field, "NUM_ANGLES")==0) {
    h->num_angles=si;
  } else if(strcasecmp(field, "CORRECTIONS_APPLIED")==0) {
    h->corrections_applied=si;
  } else if(strncasecmp(field, "NUM_Z_ELEMENTS", 14)==0) {
    cptr=strtok(value, " \t,;\n\r");
    for(i=0; i<64; i++) {
      if(cptr==NULL) break;
      h->num_z_elements[i]=(short int)atoi(cptr);
      cptr=strtok(NULL, " \t,;\n\r");
    }
  } else if(strcasecmp(field, "RING_DIFFERENCE")==0) {
    h->ring_difference=si;
  } else if(strcasecmp(field, "STORAGE_ORDER")==0) {
    h->storage_order=si;
  } else if(strcasecmp(field, "AXIAL_COMPRESSION")==0) {
    h->axial_compression=si;
  } else if(strcasecmp(field, "X_RESOLUTION")==0) {
    h->x_resolution=f;
  } else if(strcasecmp(field, "V_RESOLUTION")==0) {
    h->v_resolution=f;
  } else if(strcasecmp(field, "Z_RESOLUTION")==0) {
    h->z_resolution=f;
  } else if(strcasecmp(field, "W_RESOLUTION")==0) {
    h->w_resolution=f;
  } else if(strncasecmp(field, "FILL_GATE", 9)==0) {
    cptr=strtok(value, " \t,;\n\r");
    for(i=0; i<6; i++) {
      if(cptr==NULL) break;
      h->fill_gate[i]=(short int)atoi(cptr);
      cptr=strtok(NULL, " \t,;\n\r");
    }
  } else if(strcasecmp(field, "GATE_DURATION")==0) {
    h->gate_duration=ii;
  } else if(strcasecmp(field, "R_WAVE_OFFSET")==0) {
    h->r_wave_offset=ii;
  } else if(strcasecmp(field, "NUM_ACCEPTED_BEATS")==0) {
    h->num_accepted_beats=ii;
  } else if(strcasecmp(field, "SCALE_FACTOR")==0) {
    h->scale_factor=f;
  } else if(strcasecmp(field, "SCAN_MIN")==0) {
    h->scan_min=si;
  } else if(strcasecmp(field, "SCAN_MAX")==0) {
    h->scan_max=si;
  } else if(strcasecmp(field, "PROMPTS")==0) {
    h->prompts=ii;
  } else if(strcasecmp(field, "DELAYED")==0) {
    h->delayed=ii;
  } else if(strcasecmp(field, "MULTIPLES")==0) {
    h->multiples=ii;
  } else if(strcasecmp(field, "NET_TRUES")==0) {
    h->net_trues=ii;
  } else if(strcasecmp(field, "TOT_AVG_COR")==0) {
    h->tot_avg_cor=f;
  } else if(strcasecmp(field, "TOT_AVG_UNCOR")==0) {
    h->tot_avg_uncor=f;
  } else if(strcasecmp(field, "TOTAL_COIN_RATE")==0) {
    h->total_coin_rate=ii;
  } else if(strcasecmp(field, "FRAME_START_TIME")==0) {
    h->frame_start_time=ii;
  } else if(strcasecmp(field, "FRAME_DURATION")==0) {
    h->frame_duration=ii;
  } else if(strcasecmp(field, "DEADTIME_CORRECTION_FACTOR")==0) {
    h->deadtime_correction_factor=f;
  } else if(strncasecmp(field, "FILL_CTI", 8)==0) {
    cptr=strtok(value, " \t,;\n\r");
    for(i=0; i<90; i++) {
      if(cptr==NULL) break;
      h->fill_cti[i]=(short int)atoi(cptr);
      cptr=strtok(NULL, " \t,;\n\r");
    }
  } else if(strncasecmp(field, "FILL_USER", 9)==0) {
    cptr=strtok(value, " \t,;\n\r");
    for(i=0; i<50; i++) {
      if(cptr==NULL) break;
      h->fill_user[i]=(short int)atoi(cptr);
      cptr=strtok(NULL, " \t,;\n\r");
    }
  } else if(strncasecmp(field, "UNCOR_SINGLES", 13)==0) {
    cptr=strtok(value, " \t,;\n\r");
    for(i=0; i<128; i++) {
      if(cptr==NULL) break;
      h->fill_user[i]=(float)atof(cptr);
      cptr=strtok(NULL, " \t,;\n\r");
    }
  } else
    return(1);
  return(0);
}
/*****************************************************************************/

/*****************************************************************************/
/** Edit ECAT 7 image volume header.
\return Returns 0, if ok, and 1 or 2, if field name or
    or value is invalid.
 */
int ecat7EditVHeader(
  /** Pointer to ECAT 7 image volume header structure */
  ECAT7_imageheader *h,
  /** Field name to be changed */
  char *field,
  /** New value for the field */
  char *value,
  /** Verbose level; if <=0, then nothing is printed into stdout */
  int verbose
) {
  int ii;
  short int si;
  float f;

  if(verbose>0) printf("ecat7EditVHeader('%s', '%s')\n", field, value);

  si=atoi(value); ii=atoi(value); f=atof(value);

  if(strcasecmp(field, "DATA_TYPE")==0) {
    h->data_type=si;
  } else if(strcasecmp(field, "NUM_DIMENSIONS")==0) {
    h->num_dimensions=si;
  } else if(strcasecmp(field, "X_DIMENSION")==0) {
    h->x_dimension=si;
  } else if(strcasecmp(field, "Y_DIMENSION")==0) {
    h->y_dimension=si;
  } else if(strcasecmp(field, "Z_DIMENSION")==0) {
    h->z_dimension=si;
  } else if(strcasecmp(field, "X_OFFSET")==0) {
    h->x_offset=f;
  } else if(strcasecmp(field, "Y_OFFSET")==0) {
    h->y_offset=f;
  } else if(strcasecmp(field, "Z_OFFSET")==0) {
    h->z_offset=f;
  } else if(strcasecmp(field, "RECON_ZOOM")==0) {
    h->recon_zoom=f;
  } else if(strcasecmp(field, "SCALE_FACTOR")==0) {
    h->scale_factor=f;
  } else if(strcasecmp(field, "IMAGE_MIN")==0) {
    h->image_min=si;
  } else if(strcasecmp(field, "IMAGE_MAX")==0) {
    h->image_max=si;
  } else if(strcasecmp(field, "X_PIXEL_SIZE")==0) {
    h->x_pixel_size=f;
  } else if(strcasecmp(field, "Y_PIXEL_SIZE")==0) {
    h->y_pixel_size=f;
  } else if(strcasecmp(field, "Z_PIXEL_SIZE")==0) {
    h->z_pixel_size=f;
  } else if(strcasecmp(field, "FRAME_DURATION")==0) {
    h->frame_duration=ii;
  } else if(strcasecmp(field, "FRAME_START_TIME")==0) {
    h->frame_start_time=ii;
  } else if(strcasecmp(field, "FILTER_CODE")==0) {
    h->filter_code=si;
  } else if(strcasecmp(field, "X_RESOLUTION")==0) {
    h->x_resolution=f;
  } else if(strcasecmp(field, "Y_RESOLUTION")==0) {
    h->y_resolution=f;
  } else if(strcasecmp(field, "Z_RESOLUTION")==0) {
    h->z_resolution=f;
  } else if(strcasecmp(field, "NUM_R_ELEMENTS")==0) {
    h->num_r_elements=f;
  } else if(strcasecmp(field, "NUM_ANGLES")==0) {
    h->num_angles=f;
  } else if(strcasecmp(field, "Z_ROTATION_ANGLE")==0) {
    h->z_rotation_angle=f;
  } else if(strcasecmp(field, "DECAY_CORR_FCTR")==0) {
    h->decay_corr_fctr=f;
  } else if(strcasecmp(field, "PROCESSING_CODE")==0) {
    h->processing_code=ii;
  } else if(strcasecmp(field, "GATE_DURATION")==0) {
    h->gate_duration=ii;
  } else if(strcasecmp(field, "R_WAVE_OFFSET")==0) {
    h->r_wave_offset=ii;
  } else if(strcasecmp(field, "NUM_ACCEPTED_BEATS")==0) {
    h->num_accepted_beats=ii;
  } else if(strcasecmp(field, "FILTER_CUTOFF_FREQUENCY")==0) {
    h->filter_cutoff_frequency=f;
  } else if(strcasecmp(field, "FILTER_RESOLUTION")==0) {
    h->filter_resolution=f;
  } else if(strcasecmp(field, "FILTER_RAMP_SLOPE")==0) {
    h->filter_ramp_slope=f;
  } else if(strcasecmp(field, "FILTER_ORDER")==0) {
    h->filter_order=f;
  } else if(strcasecmp(field, "FILTER_SCATTER_FRACTION")==0) {
    h->filter_scatter_fraction=f;
  } else if(strcasecmp(field, "FILTER_SCATTER_SLOPE")==0) {
    h->filter_scatter_slope=si;
  } else if(strcasecmp(field, "ANNOTATION")==0) {
    strlcpy(h->annotation, value, 40);
  } else if(strcasecmp(field, "MT_1_1")==0) {
    h->mt_1_1=f;
  } else if(strcasecmp(field, "MT_1_2")==0) {
    h->mt_1_2=f;
  } else if(strcasecmp(field, "MT_1_3")==0) {
    h->mt_1_3=f;
  } else if(strcasecmp(field, "MT_2_1")==0) {
    h->mt_2_1=f;
  } else if(strcasecmp(field, "MT_2_2")==0) {
    h->mt_2_2=f;
  } else if(strcasecmp(field, "MT_2_3")==0) {
    h->mt_2_3=f;
  } else if(strcasecmp(field, "MT_3_1")==0) {
    h->mt_3_1=f;
  } else if(strcasecmp(field, "MT_3_2")==0) {
    h->mt_3_2=f;
  } else if(strcasecmp(field, "MT_3_3")==0) {
    h->mt_3_3=f;
  } else if(strcasecmp(field, "RFILTER_CUTOFF")==0) {
    h->rfilter_cutoff=f;
  } else if(strcasecmp(field, "RFILTER_RESOLUTION")==0) {
    h->rfilter_resolution=f;
  } else if(strcasecmp(field, "RFILTER_CODE")==0) {
    h->rfilter_code=si;
  } else if(strcasecmp(field, "RFILTER_ORDER")==0) {
    h->rfilter_order=si;
  } else if(strcasecmp(field, "ZFILTER_CUTOFF")==0) {
    h->zfilter_cutoff=f;
  } else if(strcasecmp(field, "ZFILTER_RESOLUTION")==0) {
    h->zfilter_resolution=f;
  } else if(strcasecmp(field, "ZFILTER_CODE")==0) {
    h->zfilter_code=si;
  } else if(strcasecmp(field, "ZFILTER_ORDER")==0) {
    h->zfilter_order=si;
  } else if(strcasecmp(field, "MT_1_4")==0) {
    h->mt_1_4=f;
  } else if(strcasecmp(field, "MT_2_4")==0) {
    h->mt_2_4=f;
  } else if(strcasecmp(field, "MT_3_4")==0) {
    h->mt_3_4=f;
  } else if(strcasecmp(field, "SCATTER_TYPE")==0) {
    h->scatter_type=si;
  } else if(strcasecmp(field, "RECON_TYPE")==0) {
    h->recon_type=si;
  } else if(strcasecmp(field, "RECON_VIEWS")==0) {
    h->recon_views=si;
  } else if(strncasecmp(field, "FILL_CTI", 8)==0) {
    char *cptr;
    cptr=strtok(value, " \t,;\n\r");
    for(int i=0; i<87; i++) {
      if(cptr==NULL) break;
      h->fill_cti[i]=(short int)atoi(cptr);
      cptr=strtok(NULL, " \t,;\n\r");
    }
  } else if(strncasecmp(field, "FILL_USER", 9)==0) {
    char *cptr;
    cptr=strtok(value, " \t,;\n\r");
    for(int i=0; i<50; i++) {
      if(cptr==NULL) break;
      h->fill_user[i]=(short int)atoi(cptr);
      cptr=strtok(NULL, " \t,;\n\r");
    }
  } else
    return(1);

  return(0);
}
/*****************************************************************************/

/*****************************************************************************/

