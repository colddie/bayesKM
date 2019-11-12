/// @file ecat7ift.c
/// @author Vesa Oikonen
/// @brief Procedures for reading and writing ECAT 7.x headers with IFT struct.
///
/*****************************************************************************/
#include "libtpcmisc.h"
#include "libtpcimgio.h"
/*****************************************************************************/

/*****************************************************************************/
/** Initiate image data inside ECAT_MATRIX struct */
void ematInitiate(
  /** Pointer to struct */
  ECAT_MATRIX *emat
) {
  emat->mnum=0;
  iftInit(&emat->sh);
  emat->f=(float*)NULL;
}
/** Initiate ECAT_HEADERS struct */
void ehdrInitiate(
  /** Pointer to struct */
  ECAT_HEADERS *ehdr)
{
  iftInit(&ehdr->mh);
  ehdr->nr=0;
  ehdr->m=(ECAT_MATRIX*)NULL;  
}
/*****************************************************************************/

/*****************************************************************************/
/** Free memory allocated in ECAT_MATRIX */
void ematEmpty(
  /** Pointer to ECAT matrix list */
  ECAT_MATRIX *emat
) {
  if(emat==NULL) return;
  emat->mnum=0;
  iftEmpty(&emat->sh);
  if(emat->f!=NULL) free(emat->f);
  emat->f=NULL;
}
/*****************************************************************************/

/*****************************************************************************/
/** Free memory allocated in ECAT_HEADERS */
void ehdrEmpty(
  /** Pointer to header list struct */
  ECAT_HEADERS *ehdr
) {
  int i;

  if(ehdr==NULL) return;
  iftEmpty(&ehdr->mh);
  for(i=0; i<ehdr->nr; i++) {ematEmpty(&ehdr->m[i]);}
  free(ehdr->m); ehdr->m=NULL; ehdr->nr=0;
}
/*****************************************************************************/

/*****************************************************************************/
/** Removes previous matrix contents but preserves the main header.
 *  @return Returns STATUS_OK when succesful. */
int ehdrAllocate(
  /** Pointer to list of ECAT headers */
  ECAT_HEADERS *ehdr, 
  /** Nr of headers to allocate */
  int nr
) {
  int i;

  if(ehdr==NULL || nr<1) return STATUS_FAULT;
  for(i=0; i<ehdr->nr; i++) {ematEmpty(&ehdr->m[i]);}
  free(ehdr->m); ehdr->m=NULL; ehdr->nr=0;
  ehdr->m=(ECAT_MATRIX*)malloc(nr*sizeof(ECAT_MATRIX));
  if(ehdr->m==NULL) return STATUS_NOMEMORY;
  for(i=0; i<nr; i++) ematInitiate(&ehdr->m[i]);
  ehdr->nr=nr;
  return STATUS_OK;
}
/*****************************************************************************/

/*****************************************************************************/
/** Copy ECAT 7 mainheader into IFT struct. 
 *  @return Returns STATUS_OK when succesful. */
int ecat7MHeaderToIFT(
  /** Pointer to source ECAT 7 main header */
  ECAT7_mainheader *h,
  /** Pointer to initiated IFT struct */
  IFT *ift,
  /** Verbose level */
  int verbose
) {
  int i;
  char tmp[1024], tmp2[32];
  time_t t;

  if(verbose>0) printf("ecat7MHeaderToIFT(mh, ift)\n");
  if(h==NULL || ift==NULL) return STATUS_FAULT;
  if(strncmp(h->magic_number, ECAT7V_MAGICNR, 7)!=0) return STATUS_UNKNOWNFORMAT;

  if(iftPut(ift, "magic_number", h->magic_number, NULL)!=0) return STATUS_UNSUPPORTED;
  if(iftPut(ift, "original_file_name", h->original_file_name, NULL)!=0) return STATUS_UNSUPPORTED;
  sprintf(tmp, "%d", h->sw_version);
  if(iftPut(ift, "sw_version", tmp, NULL)!=0) return STATUS_UNSUPPORTED;
  sprintf(tmp, "%d", h->system_type);
  if(iftPut(ift, "system_type", tmp, NULL)!=0) return STATUS_UNSUPPORTED;
  sprintf(tmp, "%d", h->file_type);
  if(iftPut(ift, "file_type", tmp, NULL)!=0) return STATUS_UNSUPPORTED;
  if(iftPut(ift, "serial_number", h->serial_number, NULL)!=0) return STATUS_UNSUPPORTED;
  t=h->scan_start_time;
  if(!ctime_r_int(&t, tmp)) strcpy(tmp, "1970-01-01 00:00:00");
  if(iftPut(ift, "scan_start_time", tmp, NULL)!=0) return STATUS_UNSUPPORTED;
  if(iftPut(ift, "isotope_name", h->isotope_name, NULL)!=0) return STATUS_UNSUPPORTED;
  sprintf(tmp, "%g", h->isotope_halflife);
  if(iftPut(ift, "isotope_halflife", tmp, NULL)!=0) return STATUS_UNSUPPORTED;
  if(iftPut(ift, "radiopharmaceutical", h->radiopharmaceutical, NULL)!=0) return STATUS_UNSUPPORTED;
  sprintf(tmp, "%g", h->gantry_tilt);
  if(iftPut(ift, "gantry_tilt", tmp, NULL)!=0) return STATUS_UNSUPPORTED;
  sprintf(tmp, "%g", h->gantry_rotation);
  if(iftPut(ift, "gantry_rotation", tmp, NULL)!=0) return STATUS_UNSUPPORTED;
  sprintf(tmp, "%g", h->bed_elevation);
  if(iftPut(ift, "bed_elevation", tmp, NULL)!=0) return STATUS_UNSUPPORTED;
  sprintf(tmp, "%g", h->intrinsic_tilt);
  if(iftPut(ift, "intrinsic_tilt", tmp, NULL)!=0) return STATUS_UNSUPPORTED;
  sprintf(tmp, "%d", h->wobble_speed);
  if(iftPut(ift, "wobble_speed", tmp, NULL)!=0) return STATUS_UNSUPPORTED;
  sprintf(tmp, "%d", h->transm_source_type);
  if(iftPut(ift, "transm_source_type", tmp, NULL)!=0) return STATUS_UNSUPPORTED;
  sprintf(tmp, "%g", h->distance_scanned);
  if(iftPut(ift, "distance_scanned", tmp, NULL)!=0) return STATUS_UNSUPPORTED;
  sprintf(tmp, "%g", h->transaxial_fov);
  if(iftPut(ift, "transaxial_fov", tmp, NULL)!=0) return STATUS_UNSUPPORTED;
  sprintf(tmp, "%d", h->angular_compression);
  if(iftPut(ift, "angular_compression", tmp, NULL)!=0) return STATUS_UNSUPPORTED;
  sprintf(tmp, "%d", h->coin_samp_mode);
  if(iftPut(ift, "coin_samp_mode", tmp, NULL)!=0) return STATUS_UNSUPPORTED;
  sprintf(tmp, "%d", h->axial_samp_mode);
  if(iftPut(ift, "axial_samp_mode", tmp, NULL)!=0) return STATUS_UNSUPPORTED;
  sprintf(tmp, "%E", h->ecat_calibration_factor);
  if(iftPut(ift, "ecat_calibration_factor", tmp, NULL)!=0) return STATUS_UNSUPPORTED;
  sprintf(tmp, "%d", h->calibration_units);
  if(iftPut(ift, "calibration_units", tmp, NULL)!=0) return STATUS_UNSUPPORTED;
  sprintf(tmp, "%d", h->calibration_units_label);
  if(iftPut(ift, "calibration_units_label", tmp, NULL)!=0) return STATUS_UNSUPPORTED;
  sprintf(tmp, "%d", h->compression_code);
  if(iftPut(ift, "compression_code", tmp, NULL)!=0) return STATUS_UNSUPPORTED;
  if(iftPut(ift, "study_type", h->study_type, NULL)!=0) return STATUS_UNSUPPORTED;
  if(iftPut(ift, "patient_id", h->patient_id, NULL)!=0) return STATUS_UNSUPPORTED;
  if(iftPut(ift, "patient_name", h->patient_name, NULL)!=0) return STATUS_UNSUPPORTED;
  sprintf(tmp, "%c", h->patient_sex);
  if(iftPut(ift, "patient_sex", tmp, NULL)!=0) return STATUS_UNSUPPORTED;
  sprintf(tmp, "%c", h->patient_dexterity);
  if(iftPut(ift, "patient_dexterity", tmp, NULL)!=0) return STATUS_UNSUPPORTED;
  sprintf(tmp, "%g", h->patient_age);
  if(iftPut(ift, "patient_age", tmp, NULL)!=0) return STATUS_UNSUPPORTED;
  sprintf(tmp, "%g", h->patient_height);
  if(iftPut(ift, "patient_height", tmp, NULL)!=0) return STATUS_UNSUPPORTED;
  sprintf(tmp, "%g", h->patient_weight);
  if(iftPut(ift, "patient_weight", tmp, NULL)!=0) return STATUS_UNSUPPORTED;
  t=h->patient_birth_date;
  if(!ctime_r_int(&t, tmp)) strcpy(tmp, "1970-01-01 00:00:00");
  if(iftPut(ift, "patient_birth_date", tmp, NULL)!=0) return STATUS_UNSUPPORTED;
  if(iftPut(ift, "physician_name", h->physician_name, NULL)!=0) return STATUS_UNSUPPORTED;
  if(iftPut(ift, "operator_name", h->operator_name, NULL)!=0) return STATUS_UNSUPPORTED;
  if(iftPut(ift, "study_description", h->study_description, NULL)!=0) return STATUS_UNSUPPORTED;
  sprintf(tmp, "%d", h->acquisition_type);
  if(iftPut(ift, "acquisition_type", tmp, NULL)!=0) return STATUS_UNSUPPORTED;
  sprintf(tmp, "%d", h->patient_orientation);
  if(iftPut(ift, "patient_orientation", tmp, NULL)!=0) return STATUS_UNSUPPORTED;
  if(iftPut(ift, "facility_name", h->facility_name, NULL)!=0) return STATUS_UNSUPPORTED;
  sprintf(tmp, "%d", h->num_planes);
  if(iftPut(ift, "num_planes", tmp, NULL)!=0) return STATUS_UNSUPPORTED;
  sprintf(tmp, "%d", h->num_frames);
  if(iftPut(ift, "num_frames", tmp, NULL)!=0) return STATUS_UNSUPPORTED;
  sprintf(tmp, "%d", h->num_gates);
  if(iftPut(ift, "num_gates", tmp, NULL)!=0) return STATUS_UNSUPPORTED;
  sprintf(tmp, "%d", h->num_bed_pos);
  if(iftPut(ift, "num_bed_pos", tmp, NULL)!=0) return STATUS_UNSUPPORTED;
  sprintf(tmp, "%g", h->init_bed_position);
  if(iftPut(ift, "init_bed_position", tmp, NULL)!=0) return STATUS_UNSUPPORTED;
  sprintf(tmp, "%g", h->bed_position[0]);
  for(i=1; i<15; i++) {sprintf(tmp2, " %g", h->bed_position[i]); strcat(tmp, tmp2);}
  if(iftPut(ift, "bed_position", tmp, NULL)!=0) return STATUS_UNSUPPORTED;
  sprintf(tmp, "%g", h->plane_separation);
  if(iftPut(ift, "plane_separation", tmp, NULL)!=0) return STATUS_UNSUPPORTED;
  sprintf(tmp, "%d", h->lwr_sctr_thres);
  if(iftPut(ift, "lwr_sctr_thres", tmp, NULL)!=0) return STATUS_UNSUPPORTED;
  sprintf(tmp, "%d", h->lwr_true_thres);
  if(iftPut(ift, "lwr_true_thres", tmp, NULL)!=0) return STATUS_UNSUPPORTED;
  sprintf(tmp, "%d", h->upr_true_thres);
  if(iftPut(ift, "upr_true_thres", tmp, NULL)!=0) return STATUS_UNSUPPORTED;
  if(iftPut(ift, "user_process_code", h->user_process_code, NULL)!=0) return STATUS_UNSUPPORTED;
  sprintf(tmp, "%d", h->acquisition_mode);
  if(iftPut(ift, "acquisition_mode", tmp, NULL)!=0) return STATUS_UNSUPPORTED;
  sprintf(tmp, "%g", h->bin_size);
  if(iftPut(ift, "bin_size", tmp, NULL)!=0) return STATUS_UNSUPPORTED;
  sprintf(tmp, "%g", h->branching_fraction);
  if(iftPut(ift, "branching_fraction", tmp, NULL)!=0) return STATUS_UNSUPPORTED;
  t=h->dose_start_time;
  if(!ctime_r_int(&t, tmp)) strcpy(tmp, "1970-01-01 00:00:00");
  if(iftPut(ift, "dose_start_time", tmp, NULL)!=0) return STATUS_UNSUPPORTED;
  sprintf(tmp, "%g", h->dosage);
  if(iftPut(ift, "dosage", tmp, NULL)!=0) return STATUS_UNSUPPORTED;
  sprintf(tmp, "%E", h->well_counter_corr_factor);
  if(iftPut(ift, "well_counter_corr_factor", tmp, NULL)!=0) return STATUS_UNSUPPORTED;
  if(iftPut(ift, "data_units", h->data_units, NULL)!=0) return STATUS_UNSUPPORTED;
  sprintf(tmp, "%d", h->septa_state);
  if(iftPut(ift, "septa_state", tmp, NULL)!=0) return STATUS_UNSUPPORTED;

  sprintf(tmp, "%d", h->fill_cti[0]);
  for(i=1; i<6; i++) {sprintf(tmp2, " %d", h->fill_cti[i]); strcat(tmp, tmp2);}
  if(iftPut(ift, "fill_cti", tmp, NULL)!=0) return STATUS_UNSUPPORTED;

  return STATUS_OK;
}
/*****************************************************************************/

/*****************************************************************************/
/** Copy ECAT 7 main header from IFT struct to header struct. 
*  @return Returns STATUS_OK when succesful. */
int ecat7MainheaderFromIFT(
  /** Pointer to target ECAT 7 main header */
  ECAT7_mainheader *h,
  /** Pointer to source IFT struct */
  IFT *ift,
  /** Verbose level */
  int verbose
) {
  int ii, ret;

  if(verbose>0) printf("ecat7MainheaderFromIFT(mh, ift)\n");
  if(h==NULL || ift==NULL) return STATUS_FAULT;
  if(verbose>5) iftWrite(ift, "stdout");

  for(ii=ret=0; ii<ift->keyNr; ii++) {
    if(verbose>2) 
      printf("  key := %s\n  value := %s\n", ift->item[ii].key, ift->item[ii].value);
    ret=ecat7EditMHeader(h, ift->item[ii].key, ift->item[ii].value, verbose-1);
    if(ret!=0) {
      if(verbose>0) fprintf(stderr, "Error with key '%s'\n", ift->item[ii].key);
      break;
    }
  } 
  if(ret!=0) return STATUS_FAULT;
  return STATUS_OK;
}
/*****************************************************************************/

/*****************************************************************************/
/** Copy ECAT 7 image header into IFT struct. 
*  @return Returns STATUS_OK when succesful. */
int ecat7ImageheaderToIFT(
  /** Pointer to source ECAT 7 image header */
  ECAT7_imageheader *h,
  /** Pointer to initiated IFT struct */
  IFT *ift,
  /** Verbose level */
  int verbose
) {
  int i;
  char tmp[1024], tmp2[32];

  if(verbose>0) printf("ecat7ImageheaderToIFT(h, ift)\n");
  if(h==NULL || ift==NULL) return STATUS_FAULT;
  if(h->data_type<=0) return STATUS_UNKNOWNFORMAT;

  sprintf(tmp, "%d", h->data_type);
  if(iftPut(ift, "data_type", tmp, NULL)!=0) return STATUS_UNSUPPORTED;
  sprintf(tmp, "%d", h->num_dimensions);
  if(iftPut(ift, "num_dimensions", tmp, NULL)!=0) return STATUS_UNSUPPORTED;
  sprintf(tmp, "%d", h->x_dimension);
  if(iftPut(ift, "x_dimension", tmp, NULL)!=0) return STATUS_UNSUPPORTED;
  sprintf(tmp, "%d", h->y_dimension);
  if(iftPut(ift, "y_dimension", tmp, NULL)!=0) return STATUS_UNSUPPORTED;
  sprintf(tmp, "%d", h->z_dimension);
  if(iftPut(ift, "z_dimension", tmp, NULL)!=0) return STATUS_UNSUPPORTED;
  sprintf(tmp, "%g", h->x_offset);
  if(iftPut(ift, "x_offset", tmp, NULL)!=0) return STATUS_UNSUPPORTED;
  sprintf(tmp, "%g", h->y_offset);
  if(iftPut(ift, "y_offset", tmp, NULL)!=0) return STATUS_UNSUPPORTED;
  sprintf(tmp, "%g", h->z_offset);
  if(iftPut(ift, "z_offset", tmp, NULL)!=0) return STATUS_UNSUPPORTED;
  sprintf(tmp, "%g", h->recon_zoom);
  if(iftPut(ift, "recon_zoom", tmp, NULL)!=0) return STATUS_UNSUPPORTED;
  sprintf(tmp, "%E", h->scale_factor);
  if(iftPut(ift, "scale_factor", tmp, NULL)!=0) return STATUS_UNSUPPORTED;
  sprintf(tmp, "%d", h->image_min);
  if(iftPut(ift, "image_min", tmp, NULL)!=0) return STATUS_UNSUPPORTED;
  sprintf(tmp, "%d", h->image_max);
  if(iftPut(ift, "image_max", tmp, NULL)!=0) return STATUS_UNSUPPORTED;
  sprintf(tmp, "%g", h->x_pixel_size);
  if(iftPut(ift, "x_pixel_size", tmp, NULL)!=0) return STATUS_UNSUPPORTED;
  sprintf(tmp, "%g", h->y_pixel_size);
  if(iftPut(ift, "y_pixel_size", tmp, NULL)!=0) return STATUS_UNSUPPORTED;
  sprintf(tmp, "%g", h->z_pixel_size);
  if(iftPut(ift, "z_pixel_size", tmp, NULL)!=0) return STATUS_UNSUPPORTED;
  sprintf(tmp, "%d", h->frame_duration);
  if(iftPut(ift, "frame_duration", tmp, NULL)!=0) return STATUS_UNSUPPORTED;
  sprintf(tmp, "%d", h->frame_start_time);
  if(iftPut(ift, "frame_start_time", tmp, NULL)!=0) return STATUS_UNSUPPORTED;
  sprintf(tmp, "%d", h->filter_code);
  if(iftPut(ift, "filter_code", tmp, NULL)!=0) return STATUS_UNSUPPORTED;
  sprintf(tmp, "%g", h->x_resolution);
  if(iftPut(ift, "x_resolution", tmp, NULL)!=0) return STATUS_UNSUPPORTED;
  sprintf(tmp, "%g", h->y_resolution);
  if(iftPut(ift, "y_resolution", tmp, NULL)!=0) return STATUS_UNSUPPORTED;
  sprintf(tmp, "%g", h->z_resolution);
  if(iftPut(ift, "z_resolution", tmp, NULL)!=0) return STATUS_UNSUPPORTED;
  sprintf(tmp, "%g", h->num_r_elements);
  if(iftPut(ift, "num_r_elements", tmp, NULL)!=0) return STATUS_UNSUPPORTED;
  sprintf(tmp, "%g", h->num_angles);
  if(iftPut(ift, "num_angles", tmp, NULL)!=0) return STATUS_UNSUPPORTED;
  sprintf(tmp, "%g", h->z_rotation_angle);
  if(iftPut(ift, "z_rotation_angle", tmp, NULL)!=0) return STATUS_UNSUPPORTED;
  sprintf(tmp, "%g", h->decay_corr_fctr);
  if(iftPut(ift, "decay_corr_fctr", tmp, NULL)!=0) return STATUS_UNSUPPORTED;
  sprintf(tmp, "%d", h->processing_code);
  if(iftPut(ift, "processing_code", tmp, NULL)!=0) return STATUS_UNSUPPORTED;
  sprintf(tmp, "%d", h->gate_duration);
  if(iftPut(ift, "gate_duration", tmp, NULL)!=0) return STATUS_UNSUPPORTED;
  sprintf(tmp, "%d", h->r_wave_offset);
  if(iftPut(ift, "r_wave_offset", tmp, NULL)!=0) return STATUS_UNSUPPORTED;
  sprintf(tmp, "%d", h->num_accepted_beats);
  if(iftPut(ift, "num_accepted_beats", tmp, NULL)!=0) return STATUS_UNSUPPORTED;
  sprintf(tmp, "%E", h->filter_cutoff_frequency);
  if(iftPut(ift, "filter_cutoff_frequency", tmp, NULL)!=0) return STATUS_UNSUPPORTED;
  sprintf(tmp, "%E", h->filter_resolution);
  if(iftPut(ift, "filter_resolution", tmp, NULL)!=0) return STATUS_UNSUPPORTED;
  sprintf(tmp, "%E", h->filter_ramp_slope);
  if(iftPut(ift, "filter_ramp_slope", tmp, NULL)!=0) return STATUS_UNSUPPORTED;
  sprintf(tmp, "%d", h->filter_order);
  if(iftPut(ift, "filter_order", tmp, NULL)!=0) return STATUS_UNSUPPORTED;
  sprintf(tmp, "%E", h->filter_scatter_fraction);
  if(iftPut(ift, "filter_scatter_fraction", tmp, NULL)!=0) return STATUS_UNSUPPORTED;
  sprintf(tmp, "%E", h->filter_scatter_slope);
  if(iftPut(ift, "filter_scatter_slope", tmp, NULL)!=0) return STATUS_UNSUPPORTED;
  if(iftPut(ift, "annotation", h->annotation, NULL)!=0) return STATUS_UNSUPPORTED;
  sprintf(tmp, "%g", h->mt_1_1);
  if(iftPut(ift, "mt_1_1", tmp, NULL)!=0) return STATUS_UNSUPPORTED;
  sprintf(tmp, "%g", h->mt_1_2);
  if(iftPut(ift, "mt_1_2", tmp, NULL)!=0) return STATUS_UNSUPPORTED;
  sprintf(tmp, "%g", h->mt_1_3);
  if(iftPut(ift, "mt_1_3", tmp, NULL)!=0) return STATUS_UNSUPPORTED;
  sprintf(tmp, "%g", h->mt_2_1);
  if(iftPut(ift, "mt_2_1", tmp, NULL)!=0) return STATUS_UNSUPPORTED;
  sprintf(tmp, "%g", h->mt_2_2);
  if(iftPut(ift, "mt_2_2", tmp, NULL)!=0) return STATUS_UNSUPPORTED;
  sprintf(tmp, "%g", h->mt_2_3);
  if(iftPut(ift, "mt_2_3", tmp, NULL)!=0) return STATUS_UNSUPPORTED;
  sprintf(tmp, "%g", h->mt_3_1);
  if(iftPut(ift, "mt_3_1", tmp, NULL)!=0) return STATUS_UNSUPPORTED;
  sprintf(tmp, "%g", h->mt_3_2);
  if(iftPut(ift, "mt_3_2", tmp, NULL)!=0) return STATUS_UNSUPPORTED;
  sprintf(tmp, "%g", h->mt_3_3);
  if(iftPut(ift, "mt_3_3", tmp, NULL)!=0) return STATUS_UNSUPPORTED;
  sprintf(tmp, "%g", h->rfilter_cutoff);
  if(iftPut(ift, "rfilter_cutoff", tmp, NULL)!=0) return STATUS_UNSUPPORTED;
  sprintf(tmp, "%g", h->rfilter_resolution);
  if(iftPut(ift, "rfilter_resolution", tmp, NULL)!=0) return STATUS_UNSUPPORTED;
  sprintf(tmp, "%d", h->rfilter_code);
  if(iftPut(ift, "rfilter_code", tmp, NULL)!=0) return STATUS_UNSUPPORTED;
  sprintf(tmp, "%d", h->rfilter_order);
  if(iftPut(ift, "rfilter_order", tmp, NULL)!=0) return STATUS_UNSUPPORTED;
  sprintf(tmp, "%g", h->zfilter_cutoff);
  if(iftPut(ift, "zfilter_cutoff", tmp, NULL)!=0) return STATUS_UNSUPPORTED;
  sprintf(tmp, "%g", h->zfilter_resolution);
  if(iftPut(ift, "zfilter_resolution", tmp, NULL)!=0) return STATUS_UNSUPPORTED;
  sprintf(tmp, "%d", h->zfilter_code);
  if(iftPut(ift, "zfilter_code", tmp, NULL)!=0) return STATUS_UNSUPPORTED;
  sprintf(tmp, "%d", h->zfilter_order);
  if(iftPut(ift, "zfilter_order", tmp, NULL)!=0) return STATUS_UNSUPPORTED;
  sprintf(tmp, "%g", h->mt_1_4);
  if(iftPut(ift, "mt_1_4", tmp, NULL)!=0) return STATUS_UNSUPPORTED;
  sprintf(tmp, "%g", h->mt_2_4);
  if(iftPut(ift, "mt_2_4", tmp, NULL)!=0) return STATUS_UNSUPPORTED;
  sprintf(tmp, "%g", h->mt_3_4);
  if(iftPut(ift, "mt_3_4", tmp, NULL)!=0) return STATUS_UNSUPPORTED;
  sprintf(tmp, "%d", h->scatter_type);
  if(iftPut(ift, "scatter_type", tmp, NULL)!=0) return STATUS_UNSUPPORTED;
  sprintf(tmp, "%d", h->recon_type);
  if(iftPut(ift, "recon_type", tmp, NULL)!=0) return STATUS_UNSUPPORTED;
  sprintf(tmp, "%d", h->recon_views);
  if(iftPut(ift, "recon_views", tmp, NULL)!=0) return STATUS_UNSUPPORTED;

  sprintf(tmp, "%d", h->fill_cti[0]);
  for(i=1; i<87; i++) {sprintf(tmp2, " %d", h->fill_cti[i]); strcat(tmp, tmp2);}
  if(iftPut(ift, "fill_cti", tmp, NULL)!=0) return STATUS_UNSUPPORTED;

  sprintf(tmp, "%d", h->fill_user[0]);
  for(i=1; i<49; i++) {sprintf(tmp2, " %d", h->fill_user[i]); strcat(tmp, tmp2);}
  if(iftPut(ift, "fill_user", tmp, NULL)!=0) return STATUS_UNSUPPORTED;

/*
  if(iftPut(ift, "", h->, NULL)!=0) return STATUS_UNSUPPORTED;

  sprintf(tmp, "%g", h->);
  if(iftPut(ift, "", tmp, NULL)!=0) return STATUS_UNSUPPORTED;

  sprintf(tmp, "%d", h->);
  if(iftPut(ift, "", tmp, NULL)!=0) return STATUS_UNSUPPORTED;
  */
  return STATUS_OK;
}
/*****************************************************************************/

/*****************************************************************************/
/** Copy ECAT 7 scan header into IFT struct. 
*  @return Returns STATUS_OK when succesful. */
int ecat7ScanheaderToIFT(
  /** Pointer to source ECAT 7 scan header */
  ECAT7_scanheader *h,
  /** Pointer to initiated IFT struct. Any previous contents are preserved. */
  IFT *ift,
  /** Verbose level */
  int verbose
) {
  int i;
  char tmp[1024], tmp2[32];

  if(verbose>0) printf("ecat7ScanheaderToIFT(h, ift)\n");
  if(h==NULL || ift==NULL) return STATUS_FAULT;
  if(h->data_type<=0) return STATUS_UNKNOWNFORMAT;

  sprintf(tmp, "%d", h->data_type);
  if(iftPut(ift, "data_type", tmp, NULL)!=0) return STATUS_UNSUPPORTED;
  sprintf(tmp, "%d", h->num_dimensions);
  if(iftPut(ift, "num_dimensions", tmp, NULL)!=0) return STATUS_UNSUPPORTED;
  sprintf(tmp, "%d", h->num_r_elements);
  if(iftPut(ift, "num_r_elements", tmp, NULL)!=0) return STATUS_UNSUPPORTED;
  sprintf(tmp, "%d", h->num_angles);
  if(iftPut(ift, "num_angles", tmp, NULL)!=0) return STATUS_UNSUPPORTED;
  sprintf(tmp, "%d", h->corrections_applied);
  if(iftPut(ift, "corrections_applied", tmp, NULL)!=0) return STATUS_UNSUPPORTED;
  sprintf(tmp, "%d", h->num_z_elements[0]);
  for(i=1; i<64; i++) {sprintf(tmp2, " %d", h->fill_cti[i]); strcat(tmp, tmp2);}
  if(iftPut(ift, "num_z_elements", tmp, NULL)!=0) return STATUS_UNSUPPORTED;
  sprintf(tmp, "%d", h->ring_difference);
  if(iftPut(ift, "ring_difference", tmp, NULL)!=0) return STATUS_UNSUPPORTED;
  sprintf(tmp, "%d", h->storage_order);
  if(iftPut(ift, "storage_order", tmp, NULL)!=0) return STATUS_UNSUPPORTED;
  sprintf(tmp, "%d", h->axial_compression);
  if(iftPut(ift, "axial_compression", tmp, NULL)!=0) return STATUS_UNSUPPORTED;
  sprintf(tmp, "%g", h->x_resolution);
  if(iftPut(ift, "x_resolution", tmp, NULL)!=0) return STATUS_UNSUPPORTED;
  sprintf(tmp, "%g", h->v_resolution);
  if(iftPut(ift, "v_resolution", tmp, NULL)!=0) return STATUS_UNSUPPORTED;
  sprintf(tmp, "%g", h->z_resolution);
  if(iftPut(ift, "z_resolution", tmp, NULL)!=0) return STATUS_UNSUPPORTED;
  sprintf(tmp, "%g", h->w_resolution);
  if(iftPut(ift, "w_resolution", tmp, NULL)!=0) return STATUS_UNSUPPORTED;
  sprintf(tmp, "%d", h->gate_duration);
  if(iftPut(ift, "gate_duration", tmp, NULL)!=0) return STATUS_UNSUPPORTED;
  sprintf(tmp, "%d", h->r_wave_offset);
  if(iftPut(ift, "r_wave_offset", tmp, NULL)!=0) return STATUS_UNSUPPORTED;
  sprintf(tmp, "%d", h->num_accepted_beats);
  if(iftPut(ift, "num_accepted_beats", tmp, NULL)!=0) return STATUS_UNSUPPORTED;
  sprintf(tmp, "%E", h->scale_factor);
  if(iftPut(ift, "scale_factor", tmp, NULL)!=0) return STATUS_UNSUPPORTED;
  sprintf(tmp, "%d", h->scan_min);
  if(iftPut(ift, "scan_min", tmp, NULL)!=0) return STATUS_UNSUPPORTED;
  sprintf(tmp, "%d", h->scan_max);
  if(iftPut(ift, "scan_max", tmp, NULL)!=0) return STATUS_UNSUPPORTED;
  sprintf(tmp, "%d", h->prompts);
  if(iftPut(ift, "prompts", tmp, NULL)!=0) return STATUS_UNSUPPORTED;
  sprintf(tmp, "%d", h->delayed);
  if(iftPut(ift, "delayed", tmp, NULL)!=0) return STATUS_UNSUPPORTED;
  sprintf(tmp, "%d", h->multiples);
  if(iftPut(ift, "multiples", tmp, NULL)!=0) return STATUS_UNSUPPORTED;
  sprintf(tmp, "%d", h->net_trues);
  if(iftPut(ift, "net_trues", tmp, NULL)!=0) return STATUS_UNSUPPORTED;
  sprintf(tmp, "%g", h->tot_avg_cor);
  if(iftPut(ift, "tot_avg_cor", tmp, NULL)!=0) return STATUS_UNSUPPORTED;
  sprintf(tmp, "%g", h->tot_avg_uncor);
  if(iftPut(ift, "tot_avg_uncor", tmp, NULL)!=0) return STATUS_UNSUPPORTED;
  sprintf(tmp, "%d", h->total_coin_rate);
  if(iftPut(ift, "total_coin_rate", tmp, NULL)!=0) return STATUS_UNSUPPORTED;
  sprintf(tmp, "%d", h->frame_start_time);
  if(iftPut(ift, "frame_start_time", tmp, NULL)!=0) return STATUS_UNSUPPORTED;
  sprintf(tmp, "%d", h->frame_duration);
  if(iftPut(ift, "frame_duration", tmp, NULL)!=0) return STATUS_UNSUPPORTED;
  sprintf(tmp, "%g", h->deadtime_correction_factor);
  if(iftPut(ift, "deadtime_correction_factor", tmp, NULL)!=0) return STATUS_UNSUPPORTED;
  sprintf(tmp, "%g", h->uncor_singles[0]);
  for(i=1; i<128; i++) {sprintf(tmp2, " %g", h->uncor_singles[i]); strcat(tmp, tmp2);}
  if(iftPut(ift, "uncor_singles", tmp, NULL)!=0) return STATUS_UNSUPPORTED;

/*
  if(iftPut(ift, "", h->, NULL)!=0) return STATUS_UNSUPPORTED;

  sprintf(tmp, "%g", h->);
  if(iftPut(ift, "", tmp, NULL)!=0) return STATUS_UNSUPPORTED;

  sprintf(tmp, "%d", h->);
  if(iftPut(ift, "", tmp, NULL)!=0) return STATUS_UNSUPPORTED;
  */
  return STATUS_OK;
}
/*****************************************************************************/

/*****************************************************************************/
/** Read ECAT 7 subheader from file and store in IFT struct. 
*  @return Returns STATUS_OK when succesful. */
int ecat7ReadSubheaderToIFT(
  /** File pointer to opened ECAT7 file */
  FILE *fp,
  /** ECAT7 mainheader */
  ECAT7_mainheader *h,
  /** Subheader location */
  int strtblk,
  /** Preallocated location for header data */
  IFT *ift,
  /** Verbose level */
  int verbose
) {
  ECAT7_imageheader       image_header;
  ECAT7_scanheader        scan_header;
  ECAT7_2Dscanheader      scan2D_header;
  ECAT7_2Dnormheader      norm2D_header;
  ECAT7_normheader        norm_header;
  ECAT7_attenheader       atten_header;
  ECAT7_polmapheader      polmap_header;
  int                     ret;

  if(verbose>0) printf("ecat7ReadSubheaderToIFT(fp, mh, %d, ift)\n", strtblk);
  /* Check input */
  if(fp==NULL || h==NULL || strtblk<3 || ift==NULL) return STATUS_FAULT;

  /* Read subheader and copy into IFT */
  //if(iftPut(ift, "test_key", "test_value", NULL)!=0) return STATUS_UNSUPPORTED;
  switch(h->file_type) {
    case ECAT7_ATTEN:
      ret=ecat7ReadAttenheader(fp, strtblk, &atten_header);
      if(ret!=0) return STATUS_NOSUBHEADER;
      return STATUS_UNSUPPORTED;
      break;
    case ECAT7_3DNORM:
      ret=ecat7ReadNormheader(fp, strtblk, &norm_header);
      if(ret!=0) return STATUS_NOSUBHEADER;
      return STATUS_UNSUPPORTED;
      break;
    case ECAT7_IMAGE8:
    case ECAT7_IMAGE16:
    case ECAT7_VOLUME8:
    case ECAT7_VOLUME16:
      ret=ecat7ReadImageheader(fp, strtblk, &image_header);
      if(ret!=0) return STATUS_NOSUBHEADER;
      ret=ecat7ImageheaderToIFT(&image_header, ift, verbose); 
      if(ret!=STATUS_OK) return ret;
      break;
    case ECAT7_3DSCAN:
    case ECAT7_3DSCAN8:
    case ECAT7_3DSCANFIT:
      ret=ecat7ReadScanheader(fp, strtblk, &scan_header);
      if(ret!=0) return STATUS_NOSUBHEADER;
      ret=ecat7ScanheaderToIFT(&scan_header, ift, verbose); 
      if(ret!=STATUS_OK) return ret;
      ret=ecat7WriteScanheader(fp, strtblk, &scan_header);
      break;
    case ECAT7_POLARMAP:
      ret=ecat7ReadPolmapheader(fp, strtblk, &polmap_header);
      if(ret!=0) return STATUS_NOSUBHEADER;
      return STATUS_UNSUPPORTED;
      break;
    case ECAT7_2DSCAN:
      ret=ecat7Read2DScanheader(fp, strtblk, &scan2D_header);
      if(ret!=0) return STATUS_NOSUBHEADER;
      return STATUS_UNSUPPORTED;
      break;
    case ECAT7_2DNORM:
      ret=ecat7Read2DNormheader(fp, strtblk, &norm2D_header);
      if(ret!=0) return STATUS_NOSUBHEADER;
      return STATUS_UNSUPPORTED;
      break;
    default:
      return STATUS_UNSUPPORTED;
  }
  return STATUS_OK;
}
/*****************************************************************************/

/*****************************************************************************/
/** Write ECAT 7 subheader from IFT struct into ECAT 7 file. 
*  @return Returns STATUS_OK when succesful. */
int ecat7WriteSubheaderFromIFT(
  /** File pointer to opened ECAT7 file */
  FILE *fp,
  /** ECAT7 mainheader */
  ECAT7_mainheader *h,
  /** Subheader location */
  int strtblk,
  /** Header data */
  IFT *ift,
  /** Verbose level */
  int verbose
) {
  ECAT7_imageheader       image_header;
  ECAT7_scanheader        scan_header;
  ECAT7_2Dscanheader      scan2D_header;
  ECAT7_2Dnormheader      norm2D_header;
  ECAT7_normheader        norm_header;
  ECAT7_attenheader       atten_header;
  ECAT7_polmapheader      polmap_header;
  int ii, ret;

  if(verbose>0) 
    printf("ecat7WriteSubheaderFromIFT(fp, mh, %d, ift)\n", strtblk);
  /* Check input */
  if(fp==NULL || h==NULL || strtblk<3 || ift==NULL) return STATUS_FAULT;

  /* Read subheader, set its contents from IFT, write subheader back */
  switch(h->file_type) {
    case ECAT7_ATTEN:
      ret=ecat7ReadAttenheader(fp, strtblk, &atten_header);
      if(ret!=0) return STATUS_NOSUBHEADER;
      return STATUS_UNSUPPORTED;
      break;
    case ECAT7_3DNORM:
      ret=ecat7ReadNormheader(fp, strtblk, &norm_header);
      if(ret!=0) return STATUS_NOSUBHEADER;
      return STATUS_UNSUPPORTED;
      break;
    case ECAT7_IMAGE8:
    case ECAT7_IMAGE16:
    case ECAT7_VOLUME8:
    case ECAT7_VOLUME16:
/*
      ret=ecat7ImageheaderToIFT(&image_header, ift); if(ret!=STATUS_OK) return ret;
*/
      ret=ecat7ReadImageheader(fp, strtblk, &image_header);
      if(ret!=0) return STATUS_NOSUBHEADER;
      for(ii=ret=0; ii<ift->keyNr; ii++) {
        if(verbose>7) 
          printf("  key := %s\n  value := %s\n", 
                 ift->item[ii].key, ift->item[ii].value);
        ret=ecat7EditVHeader(&image_header, ift->item[ii].key, 
                             ift->item[ii].value, verbose-2);
        if(ret!=0) {
	  if(verbose>0) 
            fprintf(stderr, "Error with key '%s'\n", ift->item[ii].key);
	  break;
	}
      } 
      if(ret!=0) return STATUS_FAULT;
      ret=ecat7WriteImageheader(fp, strtblk, &image_header);
      if(ret!=0) return STATUS_NOSUBHEADER;
      break;
    case ECAT7_3DSCAN:
    case ECAT7_3DSCAN8:
    case ECAT7_3DSCANFIT:
      ret=ecat7ReadScanheader(fp, strtblk, &scan_header);
      if(ret!=0) return STATUS_NOSUBHEADER;
      for(ii=ret=0; ii<ift->keyNr; ii++) {
        if(verbose>7) 
          printf("  key := %s\n  value := %s\n", 
            ift->item[ii].key, ift->item[ii].value);
        ret=ecat7EditSHeader(&scan_header, ift->item[ii].key, 
                             ift->item[ii].value, verbose-2);
        if(ret!=0) {
	  if(verbose>0) 
            fprintf(stderr, "Error with key '%s'\n", ift->item[ii].key);
	  break;
	}
      } 
      if(ret!=0) return STATUS_FAULT;
      ret=ecat7WriteScanheader(fp, strtblk, &scan_header);
      if(ret!=0) return STATUS_NOSUBHEADER;
      break;
    case ECAT7_POLARMAP:
      ret=ecat7ReadPolmapheader(fp, strtblk, &polmap_header);
      if(ret!=0) return STATUS_NOSUBHEADER;
      return STATUS_UNSUPPORTED;
      break;
    case ECAT7_2DSCAN:
      ret=ecat7Read2DScanheader(fp, strtblk, &scan2D_header);
      if(ret!=0) return STATUS_NOSUBHEADER;
      return STATUS_UNSUPPORTED;
      break;
    case ECAT7_2DNORM:
      ret=ecat7Read2DNormheader(fp, strtblk, &norm2D_header);
      if(ret!=0) return STATUS_NOSUBHEADER;
      return STATUS_UNSUPPORTED;
      break;
    default:
      return STATUS_UNSUPPORTED;
  }
  return STATUS_OK;
}
/*****************************************************************************/

/*****************************************************************************/
/** Read ECAT7 header contents (both main header and subheaders).
    @return Returns errstatus, which is STATUS_OK (0) when call was successful,
    and >0 in case of an error.  
 */
int ecat7ReadHeaders(
  /** Image or sinogram filename */
  const char *fname,
  /** Pointer to empty headers structure */
  ECAT_HEADERS *ehdr,
  /** Verbose level */
  int verbose
) {
  ECAT7_mainheader main_header;
  ECAT7_MATRIXLIST mlist;
  FILE *fp;
  int ret, mi;


  if(verbose>0) printf("ecat7ReadHeaders(%s, ehdr)\n", fname);
  /* Check the arguments */
  if(ehdr==NULL) return STATUS_FAULT;
  if(fname==NULL) return STATUS_FAULT;

  /* Open the file */
  if(verbose>1) printf("open %s\n", fname);
  if((fp=fopen(fname, "rb")) == NULL) return STATUS_NOFILE;

  /* Read main header */
  ret=ecat7ReadMainheader(fp, &main_header);
  if(ret) {fclose(fp); return STATUS_NOMAINHEADER;}
  /* Check for magic number */
  if(verbose>1) printf("check magic number in %s\n", fname);
  if(strncmp(main_header.magic_number, ECAT7V_MAGICNR, 7)!=0) {
    fclose(fp); return STATUS_UNKNOWNFORMAT;}
  /* Copy main header information into IFT */
  ret=ecat7MHeaderToIFT(&main_header, &ehdr->mh, verbose);
  if(ret!=STATUS_OK) {fclose(fp); return ret;}
  if(verbose>5) iftWrite(&ehdr->mh, "stdout");

  /* Read matrix list */
  ecat7InitMatlist(&mlist);
  ret=ecat7ReadMatlist(fp, &mlist, verbose-1);
  if(ret || mlist.matrixNr<1 || ecat7CheckMatlist(&mlist)) {
    fclose(fp); return STATUS_INVALIDMATLIST;}
  /* Allocate space for each matrix */
  ret=ehdrAllocate(ehdr, mlist.matrixNr);
  if(ret!=STATUS_OK) {fclose(fp); ecat7EmptyMatlist(&mlist); return ret;}
  /* Read one subheader at a time */ 
  for(mi=0; mi<mlist.matrixNr; mi++) {
    ehdr->m[mi].mnum=mlist.matdir[mi].id;
    ecat7_id_to_val(mlist.matdir[mi].id, &ehdr->m[mi].matval);
    if(verbose>2) {
        printf("frame := %d\n", ehdr->m[mi].matval.frame);
        printf("plane := %d\n", ehdr->m[mi].matval.plane);
        printf("gate := %d\n", ehdr->m[mi].matval.gate);
        printf("data := %d\n", ehdr->m[mi].matval.data);
        printf("bed := %d\n", ehdr->m[mi].matval.bed);
    }
    ret=ecat7ReadSubheaderToIFT(fp, &main_header, mlist.matdir[mi].strtblk, 
                                &ehdr->m[mi].sh, verbose);
    if(ret!=STATUS_OK) {fclose(fp); ecat7EmptyMatlist(&mlist); return ret;}
    //iftWrite(&ehdr->m[mi].sh, "stdout");
  } // next matrix
  fclose(fp); ecat7EmptyMatlist(&mlist);
  
  return STATUS_OK;
}
/*****************************************************************************/

/*****************************************************************************/
/** Write ECAT7 header contents (both main header and subheaders).
    @return Returns errstatus, which is STATUS_OK (0) when call was successful,
    and >0 in case of an error.
 */
int ecat7WriteHeaders(
  /** Image or sinogram filename */
  const char *fname,
  /** Pointer to filled headers structure */
  ECAT_HEADERS *ehdr,
  /** Verbose level */
  int verbose
) {
  ECAT7_mainheader main_header;
  ECAT7_MATRIXLIST mlist;
  FILE *fp;
  int ret, mi;


  if(verbose>0) printf("ecat7WriteHeaders(%s, ehdr)\n", fname);
  /* Check the arguments */
  if(ehdr==NULL) return STATUS_FAULT;
  if(fname==NULL) return STATUS_FAULT;

  /* Open the file */
  if((fp=fopen(fname, "r+b")) == NULL) return STATUS_NOFILE;

  /* Read main header */
  ret=ecat7ReadMainheader(fp, &main_header);
  if(ret) {fclose(fp); return STATUS_NOMAINHEADER;}
  /* Check for magic number */
  if(strncmp(main_header.magic_number, ECAT7V_MAGICNR, 7)!=0) {
    fclose(fp); return STATUS_UNKNOWNFORMAT;}
  /* Copy IFT contents into main header */
  //iftWrite(&ehdr->mh, "stdout");
  ret=ecat7MainheaderFromIFT(&main_header, &ehdr->mh, verbose);
  if(ret!=STATUS_OK) {fclose(fp); return ret;}
  //if(HEADER_TEST>5) iftWrite(&ehdr->mh, "stdout");
  /* Write mainheader */
  ret=ecat7WriteMainheader(fp, &main_header);
  if(ret!=0) {fclose(fp); return STATUS_CANNOTWRITE;}
  /* Read matrix list */
  ecat7InitMatlist(&mlist);
  ret=ecat7ReadMatlist(fp, &mlist, verbose-1);
  if(ret || mlist.matrixNr<1 || ecat7CheckMatlist(&mlist)) {
    fclose(fp); return STATUS_INVALIDMATLIST;}
#if(0)
  /* Allocate space for each matrix */
  ret=ehdrAllocate(ehdr, mlist.matrixNr);
  if(ret!=STATUS_OK) {fclose(fp); ecat7EmptyMatlist(&mlist); return ret;}
#endif
  /* Edit one subheader at a time */
  for(mi=0; mi<mlist.matrixNr; mi++) {
    ehdr->m[mi].mnum=mlist.matdir[mi].id;
    ecat7_id_to_val(mlist.matdir[mi].id, &ehdr->m[mi].matval);
    if(verbose>2) {
        printf("frame := %d\n", ehdr->m[mi].matval.frame);
        printf("plane := %d\n", ehdr->m[mi].matval.plane);
        printf("gate := %d\n", ehdr->m[mi].matval.gate);
        printf("data := %d\n", ehdr->m[mi].matval.data);
        printf("bed := %d\n", ehdr->m[mi].matval.bed);
    }
    ret=ecat7WriteSubheaderFromIFT(fp, &main_header, mlist.matdir[mi].strtblk, 
      &ehdr->m[mi].sh, verbose);
    if(ret!=STATUS_OK) {fclose(fp); ecat7EmptyMatlist(&mlist); return ret;}
    //iftWrite(&ehdr->m[mi].sh, "stdout");
  } // next matrix
  fclose(fp); ecat7EmptyMatlist(&mlist);

  return STATUS_OK;
}
/*****************************************************************************/

/*****************************************************************************/
