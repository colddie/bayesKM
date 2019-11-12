/// @file ecat63h.c
/// @author Vesa Oikonen
/// @brief Procedures for editing ECAT 6.3 header contents.
///
/*****************************************************************************/
#include "libtpcimgio.h"
/*****************************************************************************/

/*****************************************************************************/
/** Copy ECAT 6.3 main header from one struct into another.
    @return Returns 0 if ok.
 */
int ecat63CopyMainheader(
  /** Pointer to source header */
  ECAT63_mainheader *h1,
  /** Pointer to target header */
  ECAT63_mainheader *h2
) {
  int i;

  if(h1==NULL || h2==NULL) return(1);
  memset(h2, 0, sizeof(ECAT63_mainheader));

  for(i=0; i<14; i++) h2->ecat_format[i]=h1->ecat_format[i];
  for(i=0; i<14; i++) h2->fill1[i]=h1->fill1[i];
  for(i=0; i<20; i++) h2->original_file_name[i]=h1->original_file_name[i];
  h2->sw_version=h1->sw_version;
  h2->data_type=h1->data_type;
  h2->system_type=h1->system_type;
  h2->file_type=h1->file_type;
  for(i=0; i<10; i++) h2->node_id[i]=h1->node_id[i];
  h2->scan_start_day=h1->scan_start_day;
  h2->scan_start_month=h1->scan_start_month;
  h2->scan_start_year=h1->scan_start_year;
  h2->scan_start_hour=h1->scan_start_hour;
  h2->scan_start_minute=h1->scan_start_minute;
  h2->scan_start_second=h1->scan_start_second;
  for(i=0; i<8; i++) h2->isotope_code[i]=h1->isotope_code[i];
  h2->isotope_halflife=h1->isotope_halflife;
  for(i=0; i<32; i++) h2->radiopharmaceutical[i]=h1->radiopharmaceutical[i];
  h2->gantry_tilt=h1->gantry_tilt;
  h2->gantry_rotation=h1->gantry_rotation;
  h2->bed_elevation=h1->bed_elevation;
  h2->rot_source_speed=h1->rot_source_speed;
  h2->wobble_speed=h1->wobble_speed;
  h2->transm_source_type=h1->transm_source_type;
  h2->axial_fov=h1->axial_fov;
  h2->transaxial_fov=h1->transaxial_fov;
  h2->transaxial_samp_mode=h1->transaxial_samp_mode;
  h2->coin_samp_mode=h1->coin_samp_mode;
  h2->axial_samp_mode=h1->axial_samp_mode;
  h2->calibration_factor=h1->calibration_factor;
  h2->calibration_units=h1->calibration_units;
  h2->compression_code=h1->compression_code;
  for(i=0; i<12; i++) h2->study_name[i]=h1->study_name[i];
  for(i=0; i<16; i++) h2->patient_id[i]=h1->patient_id[i];
  for(i=0; i<32; i++) h2->patient_name[i]=h1->patient_name[i];
  h2->patient_sex=h1->patient_sex;
  for(i=0; i<10; i++) h2->patient_age[i]=h1->patient_age[i];
  for(i=0; i<10; i++) h2->patient_height[i]=h1->patient_height[i];
  for(i=0; i<10; i++) h2->patient_weight[i]=h1->patient_weight[i];
  h2->patient_dexterity=h1->patient_dexterity;
  for(i=0; i<32; i++) h2->physician_name[i]=h1->physician_name[i];
  for(i=0; i<32; i++) h2->operator_name[i]=h1->operator_name[i];
  for(i=0; i<32; i++) h2->study_description[i]=h1->study_description[i];
  h2->acquisition_type=h1->acquisition_type;
  h2->bed_type=h1->bed_type;
  h2->septa_type=h1->septa_type;
  for(i=0; i<20; i++) h2->facility_name[i]=h1->facility_name[i];
  h2->num_planes=h1->num_planes;
  h2->num_frames=h1->num_frames;
  h2->num_gates=h1->num_gates;
  h2->num_bed_pos=h1->num_bed_pos;
  h2->init_bed_position=h1->init_bed_position;
  for(i=0; i<15; i++) h2->bed_offset[i]=h1->bed_offset[i];
  h2->plane_separation=h1->plane_separation;
  h2->lwr_sctr_thres=h1->lwr_sctr_thres;
  h2->lwr_true_thres=h1->lwr_true_thres;
  h2->upr_true_thres=h1->upr_true_thres;
  h2->collimator=h1->collimator;
  for(i=0; i<10; i++) h2->user_process_code[i]=h1->user_process_code[i];
  for(i=0; i<20; i++) h2->fill2[i]=h1->fill2[i];
  return(0);
}
/*****************************************************************************/

/*****************************************************************************/
/** Copy ECAT 6.3 scan sub header from one struct into another.
    @return Returns 0 if ok.
 */
int ecat63CopyScanheader(
  /** Pointer to source header */
  ECAT63_scanheader *h1,
  /** Pointer to target header */
  ECAT63_scanheader *h2
) {
  int i;

  if(h1==NULL || h2==NULL) return(1);
  memset(h2, 0, sizeof(ECAT63_scanheader));

  for(i=0; i<126; i++) h2->fill1[i]=h1->fill1[i];
  h2->data_type=h1->data_type;
  h2->dimension_1=h1->dimension_1;
  h2->dimension_2=h1->dimension_2;
  h2->smoothing=h1->smoothing;
  h2->processing_code=h1->processing_code;
  h2->sample_distance=h1->sample_distance;
  h2->isotope_halflife=h1->isotope_halflife;
  h2->frame_duration_sec=h1->frame_duration_sec;
  h2->gate_duration=h1->gate_duration;
  h2->r_wave_offset=h1->r_wave_offset;
  h2->scale_factor=h1->scale_factor;
  h2->scan_min=h1->scan_min;
  h2->scan_max=h1->scan_max;
  h2->prompts=h1->prompts;
  h2->delayed=h1->delayed;
  h2->multiples=h1->multiples;
  h2->net_trues=h1->net_trues;
  for(i=0; i<16; i++) h2->cor_singles[i]=h1->cor_singles[i];
  for(i=0; i<16; i++) h2->uncor_singles[i]=h1->uncor_singles[i];
  h2->tot_avg_cor=h1->tot_avg_cor;
  h2->tot_avg_uncor=h1->tot_avg_uncor;
  h2->total_coin_rate=h1->total_coin_rate;
  h2->frame_start_time=h1->frame_start_time;
  h2->frame_duration=h1->frame_duration;
  h2->loss_correction_fctr=h1->loss_correction_fctr;
  for(i=0; i<22; i++) h2->fill2[i]=h1->fill2[i];
  return(0);
}
/*****************************************************************************/

/*****************************************************************************/
/** Edit ECAT 6.3 main header.
    @return Returns 0, if ok, and 1 or 2, if field name or or value is invalid.
 */
int ecat63EditMHeader(
  /** Pointer to ECAT 6.3 mainheader structure */
  ECAT63_mainheader *h,
  /** Field name to be changed */
  char *field,
  /** New value for the field */
  char *value,
  /** Verbose level; if <=0, then nothing is printed into stdout */
  int verbose
) {
//  int yy, mm, dd;
  short int si;
  float f;

  if(verbose>0) printf("ecat63EditMHeader('%s', '%s')\n", field, value);
  si=atoi(value); f=atof(value);

  if(strcmp(field, "ecat_format")==0 || strcmp(field, "magic_number")==0) {
    strlcpy(h->ecat_format, value, 14);
  } else if(strcmp(field, "fill1")==0) {
    strlcpy(h->fill1, value, 14);
  } else if(strcmp(field, "original_file_name")==0) {
    strlcpy(h->original_file_name, value, 20);
  } else if(strcmp(field, "sw_version")==0) {
    if(si<=0) return(2); else h->sw_version=si;
  } else if(strcmp(field, "data_type")==0) {
    if(si<0) return(2); else h->data_type=si;
  } else if(strcmp(field, "system_type")==0) {
    if(si<0) return(2); else h->system_type=si;
  } else if(strcmp(field, "file_type")==0) {
    if(si<0) return(2); else h->file_type=si;
  } else if(strcmp(field, "node_id")==0 || strcmp(field, "serial_number")==0) {
    strlcpy(h->node_id, value, 10);
  } else if(strcmp(field, "scan_start_day")==0) {
    if(si<0) return(2); else h->scan_start_day=si;
  } else if(strcmp(field, "scan_start_month")==0) {
    if(si<0) return(2); else h->scan_start_month=si;
  } else if(strcmp(field, "scan_start_year")==0) {
    if(si<0) return(2); else h->scan_start_year=si;
  } else if(strcmp(field, "scan_start_hour")==0) {
    if(si<0) return(2); else h->scan_start_hour=si;
  } else if(strcmp(field, "scan_start_minute")==0) {
    if(si<0) return(2); else h->scan_start_minute=si;
  } else if(strcmp(field, "scan_start_second")==0) {
    if(si<0) return(2); else h->scan_start_second=si;
  } else if(strcmp(field, "scan_start_time")==0) {
    int YYYY, MM, DD, hh, mm, ss, n;
    n=sscanf(value, "%d-%d-%d %d:%d:%d", &YYYY, &MM, &DD, &hh, &mm, &ss);
    if(n!=6 && n!=3) return(2);
    h->scan_start_year=YYYY;
    h->scan_start_month=MM;
    h->scan_start_day=DD;
    if(n==6) {
      h->scan_start_hour=hh;
      h->scan_start_minute=mm;
      h->scan_start_second=ss;
    }
  } else if(strcmp(field, "isotope_code")==0 || strcmp(field, "isotope_name")==0) {
    strlcpy(h->isotope_code, value, 8);
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
  } else if(strcmp(field, "rot_source_speed")==0) {
    h->rot_source_speed=si;
  } else if(strcmp(field, "wobble_speed")==0) {
    h->wobble_speed=si;
  } else if(strcmp(field, "transm_source_type")==0) {
    h->transm_source_type=si;
  } else if(strcmp(field, "axial_fov")==0) {
    h->axial_fov=f;
  } else if(strcmp(field, "transaxial_fov")==0) {
    h->transaxial_fov=f;
  } else if(strcmp(field, "transaxial_samp_mode")==0) {
    h->transaxial_samp_mode=si;
  } else if(strcmp(field, "coin_samp_mode")==0) {
    h->coin_samp_mode=si;
  } else if(strcmp(field, "axial_samp_mode")==0) {
    h->axial_samp_mode=si;
  } else if(strcmp(field, "calibration_factor")==0) {
    h->calibration_factor=f;
  } else if(strcmp(field, "calibration_units")==0) {
    h->calibration_units=si;
  } else if(strcmp(field, "compression_code")==0) {
    h->compression_code=si;
  } else if(strcmp(field, "study_name")==0) {
    strlcpy(h->study_name, value, 12);
  } else if(strcmp(field, "patient_id")==0) {
    strlcpy(h->patient_id, value, 16);
  } else if(strcmp(field, "patient_name")==0) {
    strlcpy(h->patient_name, value, 32);
  } else if(strcmp(field, "patient_sex")==0) {
    h->patient_sex=value[0];
  } else if(strcmp(field, "patient_age")==0) {
    strlcpy(h->patient_name, value, 10);
  } else if(strcmp(field, "patient_height")==0) {
    strlcpy(h->patient_height, value, 10);
  } else if(strcmp(field, "patient_weight")==0) {
    strlcpy(h->patient_weight, value, 10);
  } else if(strcmp(field, "patient_dexterity")==0) {
    h->patient_dexterity=value[0];
  } else if(strcmp(field, "physician_name")==0) {
    strlcpy(h->physician_name, value, 32);
  } else if(strcmp(field, "operator_name")==0) {
    strlcpy(h->operator_name, value, 32);
  } else if(strcmp(field, "study_description")==0) {
    strlcpy(h->study_description, value, 32);
  } else if(strcmp(field, "acquisition_type")==0) {
    h->acquisition_type=si;
  } else if(strcmp(field, "bed_type")==0) {
    h->bed_type=si;
  } else if(strcmp(field, "septa_type")==0) {
    h->septa_type=si;
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
  } else if(strcmp(field, "bed_offset")==0) {
    sscanf(value, "%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f",
      h->bed_offset+0, h->bed_offset+1, h->bed_offset+2,
      h->bed_offset+3, h->bed_offset+4, h->bed_offset+5,
      h->bed_offset+6, h->bed_offset+7, h->bed_offset+8,
      h->bed_offset+9, h->bed_offset+10, h->bed_offset+11,
      h->bed_offset+12, h->bed_offset+13, h->bed_offset+14
    );
  } else if(strcmp(field, "plane_separation")==0) {
    h->plane_separation=f;
  } else if(strcmp(field, "lwr_sctr_thres")==0) {
    h->lwr_sctr_thres=si;
  } else if(strcmp(field, "lwr_true_thres")==0) {
    h->lwr_true_thres=si;
  } else if(strcmp(field, "upr_true_thres")==0) {
    h->upr_true_thres=si;
  } else if(strcmp(field, "collimator")==0) {
    h->collimator=f;
  } else if(strcmp(field, "user_process_code")==0) {
    strlcpy(h->user_process_code, value, 10);
  } else if(strcasecmp(field, "FILL2")==0) {
    char *cptr;
    cptr=strtok(value, " \t,;\n\r");
    for(int i=0; i<20; i++) {
      if(cptr==NULL) break;
      h->fill2[i]=(short int)atoi(cptr);
      cptr=strtok(NULL, " \t,;\n\r");
    }
  } else
    return(1);

  return(0);
}
/*****************************************************************************/

/*****************************************************************************/
