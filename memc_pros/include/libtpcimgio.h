/// @file libtpcimgio.h
/// @brief Header file for libtpcimgio.
/// @author Vesa Oikonen
///

#ifdef __cplusplus
extern "C" {
#endif


#ifndef _LIBTPCIMGIO_H
#define _LIBTPCIMGIO_H
/*****************************************************************************/

/*****************************************************************************/
#include "tpcclibConfig.h"
/*****************************************************************************/
//#include <omp.h>
/*****************************************************************************/
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <strings.h>
#include <ctype.h>
#include <float.h>
#include <unistd.h>
#include <locale.h>
#include <time.h>
/*****************************************************************************/
#include "libtpcmisc.h"
/*****************************************************************************/

/*****************************************************************************/
#ifndef BACKUP_EXTENSION
/** Backup file extension */
#define BACKUP_EXTENSION ".bak"
#endif 
/*****************************************************************************/

/*****************************************************************************/
/* Analyze */

/// @cond
#define ANALYZE_HEADER_KEY_SIZE 40
#define ANALYZE_HEADER_IMGDIM_SIZE 108
#define ANALYZE_HEADER_HISTORY_SIZE 200
#define ANALYZE_FLIP_DEFAULT 1
/// @endcond

/** Analyze 7.5 datatype */
#define ANALYZE_DT_NONE 0
/** Analyze 7.5 datatype */
#define ANALYZE_DT_UNKNOWN 0
/** Analyze 7.5 datatype, 1 bit */
#define ANALYZE_DT_BINARY 1
/** Analyze 7.5 datatype, 8 bits */
#define ANALYZE_DT_UNSIGNED_CHAR 2
/** Analyze 7.5 datatype, 16 bits */
#define ANALYZE_DT_SIGNED_SHORT 4
/** Analyze 7.5 datatype, 32 bits */
#define ANALYZE_DT_SIGNED_INT 8
/** Analyze 7.5 datatype, 32 bits */
#define ANALYZE_DT_FLOAT 16
/** Analyze 7.5 datatype, 64 bits (two floats) */
#define ANALYZE_DT_COMPLEX 32
/** Analyze 7.5 datatype, 64 bits */
#define ANALYZE_DT_DOUBLE 64
/** Analyze 7.5 datatype */
#define ANALYZE_DT_RGB 128
/** Analyze 7.5 datatype */
#define ANALYZE_DT_ALL 255

/** Verbose prints from Analyze functions */
int ANALYZE_TEST;

/** Analyze header key */
typedef struct {
  /** Size of header file in bytes; required. */
  int sizeof_hdr;
  /** data_type */
  char data_type[10];
  /** db_name */
  char db_name[18];
  /** Required field, should be 16384 if image file is created as contiguous 
      a minimum extent size. */
  int extents;
  /** session_error */
  short int session_error;
  /** Required field; 'r' to indicate that all images and volumes are
      the same size. */
  char regular;
  /** hkey_un0 */
  char hkey_un0;
} ANALYZE_HEADER_KEY;

/** Analyze header imgdim */
typedef struct {
  /** Array of the image dimensions; nr of dimensions (usually 4), x dimension,
      y dimension, z dimension, and time dimension (nr of image volumes). */
  short int dim[8];
  /** unused */
  short int unused8;
  /** unused */
  short int unused9;
  /** unused */
  short int unused10;
  /** unused */
  short int unused11;
  /** unused */
  short int unused12;
  /** unused */
  short int unused13;
  /** unused */
  short int unused14;
  /** datatype: ANALYZE_DT_NONE, ANALYZE_DT_UNKNOWN, ANALYZE_DT_BINARY (1 bit),
      ANALYZE_DT_UNSIGNED_CHAR (8 bits), ANALYZE_DT_SIGNED_SHORT (16 bits),
      ANALYZE_DT_SIGNED_INT (32 bits), ANALYZE_DT_FLOAT (32 bits), 
      ANALYZE_DT_COMPLEX (64 bits), ANALYZE_DT_DOUBLE (64 bits), 
      ANALYZE_DT_RGB.
      Notice that not all SW supports all data types.
  */
  short int datatype;
  /** Nr of bits per pixel; 1, 8, 16, 32, or 64. */
  short int bitpix;
  /** Unused */
  short int dim_un0;
  /** Pixel dimensions in mm or ms, corresponding to dim[]. */
  float pixdim[8];
  /** Byte offset in the .img file at which voxels start; negative value
      specifies that the absolute value is applied for every image in the file.
   */
  float vox_offset;
  /** funused1 */
  float funused1;
  /** funused2 */
  float funused2;
  /** funused3 */
  float funused3;
  /** Maximum of calibrated values. */
  float cal_max;
  /** Minimum of calibrated values. */
  float cal_min;
  /** compressed */
  float compressed;
  /** verified */
  float verified;
  /** Maximum pixel value of the whole database. */
  int glmax;
  /** Minimum pixel value of the whole database. */
  int glmin;
} ANALYZE_HEADER_IMGDIM;

/** Analyze header history */
typedef struct {
  /** descrip */
  char descrip[80];
  /** aux_file */
  char aux_file[24];
  /** Slice orientation for the dataset; 0=transverse unflipped, 1=coronal 
      unflipped, 2=sagittal unflipped, 3=transverse flipped, 4=coronal flipped,
      5=sagittal flipped. */
  char orient;
  /** originator */
  char originator[10];
  /** generated */
  char generated[10];
  /** scannum */
  char scannum[10];
  /** patient_id */
  char patient_id[10];
  /** exp_date */
  char exp_date[10];
  /** exp_time */
  char exp_time[10];
  /** hist_un0 */
  char hist_un0[3];
  /** views */
  int views;
  /** vols_added */
  int vols_added;
  /** start_field */
  int start_field;
  /** field_skip */
  int field_skip;
  /** omax */
  int omax;
  /** omin */
  int omin;
  /** smax */
  int smax;
  /** smin */
  int smin;
} ANALYZE_HEADER_HISTORY;


/** Analyze header combo */
typedef struct {
  /** hk */
  ANALYZE_HEADER_KEY hk;
  /** dime */
  ANALYZE_HEADER_IMGDIM dime;
  /** hist */
  ANALYZE_HEADER_HISTORY hist;
  /** Originally 1=little endian, 0=big endian; not stored in file. */
  int little;
} ANALYZE_DSR;
/*****************************************************************************/

/*****************************************************************************/
/* ECAT 7.x */

#ifndef MatBLKSIZE
/** ECAT matrix block size */
#define MatBLKSIZE 512
#endif
#ifndef MatFirstDirBlk
/** ECAT matrix directory start block */
#define MatFirstDirBlk 2
#endif

/// @cond
#define ECAT7V_MAGICNR "MATRIX72v"
#define ECAT7S_MAGICNR "MATRIX7011"
#define ECAT7_SW_VERSION 72
/// @endcond

/** ECAT7 matrix data type */
#define	ECAT7_BYTE      1
/** ECAT7 matrix data type */
#define	ECAT7_VAXI2     2
/** ECAT7 matrix data type */
#define ECAT7_VAXI4     3
/** ECAT7 matrix data type */
#define ECAT7_VAXR4     4
/** ECAT7 matrix data type */
#define ECAT7_IEEER4    5
/** ECAT7 matrix data type */
#define	ECAT7_SUNI2     6
/** ECAT7 matrix data type */
#define	ECAT7_SUNI4     7
/** ECAT matrix data type */
#define M68K_I2 SUN_I2
/** ECAT matrix data type */
#define M68K_I4 SUN_I4

/** ECAT7 matrix filetype */
#define ECAT7_UNKNOWN   0
/** ECAT7 matrix filetype */
#define ECAT7_2DSCAN    1
/** ECAT7 matrix filetype */
#define ECAT7_IMAGE16   2
/** ECAT7 matrix filetype */
#define ECAT7_ATTEN     3
/** ECAT7 matrix filetype */
#define ECAT7_2DNORM    4
/** ECAT7 matrix filetype */
#define ECAT7_POLARMAP  5
/** ECAT7 matrix filetype */
#define ECAT7_VOLUME8   6
/** ECAT7 matrix filetype */
#define ECAT7_VOLUME16  7
/** ECAT7 matrix filetype */
#define ECAT7_PROJ      8
/** ECAT7 matrix filetype */
#define ECAT7_PROJ16    9
/** ECAT7 matrix filetype */
#define ECAT7_IMAGE8    10
/** ECAT7 matrix filetype */
#define ECAT7_3DSCAN    11
/** ECAT7 matrix filetype */
#define ECAT7_3DSCAN8   12
/** ECAT7 matrix filetype */
#define ECAT7_3DNORM    13
/** ECAT7 matrix filetype */
#define ECAT7_3DSCANFIT 14

/** ECAT7 patient orientation */
#define ECAT7_Feet_First_Prone           0
/** ECAT7 patient orientation */
#define ECAT7_Head_First_Prone           1
/** ECAT7 patient orientation */
#define ECAT7_Feet_First_Supine          2
/** ECAT7 patient orientation */
#define ECAT7_Head_First_Supine          3
/** ECAT7 patient orientation */
#define ECAT7_Feet_First_Decubitus_Right 4
/** ECAT7 patient orientation */
#define ECAT7_Head_First_Decubitus_Right 5
/** ECAT7 patient orientation */
#define ECAT7_Feet_First_Decubitus_Left  6
/** ECAT7 patient orientation */
#define ECAT7_Head_First_Decubitus_Left  7
/** ECAT7 patient orientation */
#define ECAT7_Unknown_Orientation        8

/** Error message from ECAT 7 functions */
char ecat7errmsg[128];
/** Verbose prints from ECAT 7 functions */
int ECAT7_TEST;

/** ECAT7 main header; 512 bytes, if packed and not aligned; tpcclib does
    not require packing, but packing may be necessary if used as library 
    in other SW. */
typedef struct ecat7_mainheader {  /* 512 bytes, if packed (not aligned) */
  /** Unix file type identification number */
  char      magic_number[14];
  /** Scan file's creation number */
  char      original_file_name[32];
  /** sw_version */
  short int sw_version;
  /** Scanner model */
  short int system_type;
  /** Matrix file type */
  short int file_type;
  /** Serial number of the gantry */
  char      serial_number[10];
  /** Date and time when acquisition was started (sec from base time) */
  int       scan_start_time;
  /** String representation of the isotope */
  char      isotope_name[8];
  /** Half-life of isotope (sec) */
  float     isotope_halflife;
  /** String representation of the tracer name */
  char      radiopharmaceutical[32];
  /** Angle (degrees) */
  float     gantry_tilt;
  /** Angle (degrees) */
  float     gantry_rotation;
  /** Bed height from lowest point (cm) */
  float     bed_elevation;
  /** intrinsic_tilt */
  float     intrinsic_tilt;
  /** wobble_speed */
  short int wobble_speed;
  /** transm_source_type */
  short int transm_source_type;
  /** Total distance scanned (cm) */
  float     distance_scanned;
  /** Diameter of transaxial view (cm) */
  float     transaxial_fov;
  /** 0=no mash, 1=mash of 2, 2=mash of 4 */
  short int angular_compression;
  /** 0=Net trues, 1=Prompts and Delayed, 3=Prompts, Delayed, and Multiples */
  short int coin_samp_mode;
  /** 0=Normal, 1=2X, 2=3X */
  short int axial_samp_mode;
  /** ECAT calibration factor */
  float     ecat_calibration_factor;
  /** 0=Uncalibrated; 1=Calibrated; 2=Processed */
  short int calibration_units;
  /** Whether data_units[] is filled or not? */
  short int calibration_units_label;
  /** compression_code */
  short int compression_code;
  /** study_type */
  char      study_type[12];
  /** patient_id */
  char      patient_id[16];
  /** patient_name */
  char      patient_name[32];
  /** patient_sex */
  char      patient_sex;
  /** patient_dexterity */
  char      patient_dexterity;
  /** Patient age (years) */
  float     patient_age;
  /** Patient height (cm) */
  float     patient_height;
  /** Patient weight (kg) */
  float     patient_weight;
  /** YYYYMMDD. In HR+ files this field may contain birth date as seconds from
      time zero, thus negative number when born before 1970, but those are
      converted to YYYYMMDD when file is read. */
  int       patient_birth_date;
  /** physician_name */
  char      physician_name[32];
  /** operator_name */
  char      operator_name[32];
  /** study_description */
  char      study_description[32];
  /** 0=Undefined; 1=Blank; 2=Transmission; 3=Static emission; 4=Dynamic emission;
      5=Gated emission; 6=Transmission rectilinear; 7=Emission rectilinear */
  short int acquisition_type;
  /** patient_orientation */
  short int patient_orientation;
  /** facility_name */
  char      facility_name[20];
  /** num_planes */
  short int num_planes;
  /** Highest frame number in partially reconstruction files */
  short int num_frames;
  /** num_gates */
  short int num_gates;
  /** num_bed_pos */
  short int num_bed_pos;
  /** init_bed_position */
  float     init_bed_position;
  /** bed_position */
  float     bed_position[15];
  /** Physical distance between adjacent planes (cm) */
  float     plane_separation;
  /** lwr_sctr_thres */
  short int lwr_sctr_thres;
  /** lwr_true_thres */
  short int lwr_true_thres;
  /** upr_true_thres */
  short int upr_true_thres;
  /** user_process_code */
  char      user_process_code[10];
  /** acquisition_mode */
  short int acquisition_mode;
  /** Width of view sample (cm) */
  float     bin_size;
  /** Fraction of decay by positron emission (included in calibration factor
      in TPC) */
  float     branching_fraction;
  /** Time of injection */
  int       dose_start_time;
  /** Radiopharmaceutical dosage at time of injection (Bq/cc) */
  float     dosage;
  /** well_counter_corr_factor */
  float     well_counter_corr_factor;
  /** Free text field; fixed strings: "ECAT counts/sec", "Bq/cc" */
  char      data_units[32];
  /** septa_state */
  short int septa_state;
  /** fill */
  short int fill_cti[6];   
} ECAT7_mainheader;

/** ECAT7 image header (512 bytes) */
typedef struct ecat7_imageheader {
  /** data_type */
  short int data_type;
  /** num_dimensions */
  short int num_dimensions;
  /** x_dimension */
  short int x_dimension;
  /** y_dimension */
  short int y_dimension;
  /** z_dimension */
  short int z_dimension;
  /** cm */
  float     x_offset;
  /** cm */
  float     y_offset;
  /** cm */
  float     z_offset;
  /** Reconstruction magnification factor */
  float     recon_zoom;
  /** scale_factor */
  float     scale_factor;
  /** image_min */
  short int image_min;
  /** image_max */
  short int image_max;
  /** X dimension pixel size (cm) */
  float     x_pixel_size;
  /** Y dimension pixel size (cm) */
  float     y_pixel_size;
  /** Z dimension pixel size (cm) */
  float     z_pixel_size;
  /** msec */
  int       frame_duration;
  /** Offset from first frame (msec) */
  int       frame_start_time;
  /** filter_code */
  short int filter_code;
  /** cm */
  float     x_resolution;
  /** cm */
  float     y_resolution;
  /** cm */
  float     z_resolution;
  /** Number R elements from sinogram */
  float     num_r_elements;
  /** Nr of angles from sinogram */
  float     num_angles;
  /** Rotation in the xy plane (degrees) */
  float     z_rotation_angle;
  /** decay_corr_fctr */
  float     decay_corr_fctr;
  /** processing_code */
  int       processing_code;
  /** gate_duration */
  int       gate_duration;
  /** r_wave_offset */
  int       r_wave_offset;
  /** num_accepted_beats */
  int       num_accepted_beats;
  /** filter_cutoff_frequency */
  float     filter_cutoff_frequency;
  /** filter_resolution */
  float     filter_resolution;
  /** filter_ramp_slope */
  float     filter_ramp_slope;
  /** filter_order */
  short int filter_order;
  /** filter_scatter_fraction */
  float     filter_scatter_fraction;
  /** filter_scatter_slope */
  float     filter_scatter_slope;
  /** annotation */
  char      annotation[40];
  /** X Rotation angle, matrix transformation element (1,1) */
  float     mt_1_1;
  /** Y Rotation angle, matrix transformation element (1,2) */
  float     mt_1_2;
  /** Z Rotation angle, matrix transformation element (1,3) */
  float     mt_1_3;
  /** X translation, matrix transformation element (2,1) */
  float     mt_2_1;
  /** Y translation, matrix transformation element (2,2) */
  float     mt_2_2;
  /** Z translation, matrix transformation element (2,3) */
  float     mt_2_3;
  /** X scale factor, matrix transformation element (3,1) */
  float     mt_3_1;
  /** Y scale factor, matrix transformation element (3,2) */
  float     mt_3_2;
  /** Z scale factor, matrix transformation element (3,3) */
  float     mt_3_3;
  /** rfilter_cutoff */
  float     rfilter_cutoff;
  /** rfilter_resolution */
  float     rfilter_resolution;
  /** rfilter_code */
  short int rfilter_code;
  /** rfilter_order */
  short int rfilter_order;
  /** zfilter_cutoff */
  float     zfilter_cutoff;
  /** zfilter_resolution */
  float     zfilter_resolution;
  /** zfilter_code */
  short int zfilter_code;
  /** zfilter_order */
  short int zfilter_order;
  /** Matrix transformation element (1,4) */
  float     mt_1_4;
  /** Matrix transformation element (2,4) */
  float     mt_2_4;
  /** Matrix transformation element (3,4) */
  float     mt_3_4;
  /** scatter_type */
  short int scatter_type;
  /** recon_type */
  short int recon_type;
  /** recon_views */
  short int recon_views;
  /** fill */
  short int fill_cti[87];
  /** fill */
  short int fill_user[49];
} ECAT7_imageheader;

/** ECAT7 3D sinogram header (1024 bytes) */
typedef struct ecat7_scanheader {
  /** data_type */
  short int data_type;
  /** num_dimensions */
  short int num_dimensions;
  /** Total elements collected (r dimension ) */
  short int num_r_elements;
  /** Total views collected (theta dimension) */
  short int num_angles;
  /** Bit 0 - Norm; Bit 1 - Atten; Bit 2 - Smooth */
  short int corrections_applied;
  /** Nr of elements in z dimension for each ring difference segment */
  short int num_z_elements[64];
  /** Max ring difference (d dimension) in this frame */
  short int ring_difference;
  /** RThetaZD or RZThetaD */
  short int storage_order;
  /** Span */
  short int axial_compression;
  /** Resolution in r dimension (cm) */
  float     x_resolution;
  /** Resolution in Theta dimension (rad) */
  float     v_resolution;
  /** Resolution in z dimension (cm) */
  float     z_resolution;
  /** w_resolution */
  float     w_resolution;
  /** fill_gate */
  short int fill_gate[6];
  /** gate_duration */
  int       gate_duration;
  /** Time from start of first gate (msec) */
  int       r_wave_offset;
  /** num_accepted_beats */
  int       num_accepted_beats;
  /** scale_factor */
  float     scale_factor;
  /** scan_min */
  short int scan_min;
  /** scan_max */
  short int scan_max;
  /** prompts */
  int       prompts;
  /** delayed */
  int       delayed;
  /** multiples */
  int       multiples;
  /** net_trues */
  int       net_trues;
  /** tot_avg_cor */
  float     tot_avg_cor;
  /** tot_avg_uncor */
  float     tot_avg_uncor;
  /** total_coin_rate */
  int       total_coin_rate;
  /** Time offset from first frame (msec) */
  int       frame_start_time;
  /** Total duration of current frame (msec) */
  int       frame_duration;
  /** deadtime_correction_factor */
  float     deadtime_correction_factor;
  /** fill_cti */
  short int fill_cti[90];
  /** fill_user */
  short int fill_user[50];
  /** uncor_singles */
  float     uncor_singles[128];
} ECAT7_scanheader;

/** ECAT7 2D sinogram header (512 bytes) */
typedef struct ecat7_2Dscanheader {
  /** data_type */
  short int data_type;
  /** num_dimensions */
  short int num_dimensions;
  /** num_r_elements */
  short int num_r_elements;
  /** num_angles */
  short int num_angles;
  /** corrections_applied */
  short int corrections_applied;
  /** num_z_elements */
  short int num_z_elements;
  /** ring_difference */
  short int ring_difference;
  /** x_resolution */
  float     x_resolution;
  /** y_resolution */
  float     y_resolution;
  /** z_resolution */
  float     z_resolution;
  /** w_resolution */
  float     w_resolution;
  /** fill_gate */
  short int fill_gate[6];
  /** gate_duration */
  int       gate_duration;
  /** r_wave_offset */
  int       r_wave_offset;
  /** num_accepted_beats */
  int       num_accepted_beats;
  /** scale_factor */
  float     scale_factor;
  /** scan_min */
  short int scan_min;
  /** scan_max */
  short int scan_max;
  /** prompts */
  int       prompts;
  /** delayed */
  int       delayed;
  /** multiples */
  int       multiples;
  /** net_trues */
  int       net_trues;
  /** cor_singles */
  float     cor_singles[16];
  /** uncor_singles */
  float     uncor_singles[16];
  /** tot_avg_cor */
  float     tot_avg_cor;
  /** tot_avg_uncor */
  float     tot_avg_uncor;
  /** total_coin_rate */
  int       total_coin_rate;
  /** frame_start_time */
  int       frame_start_time;
  /** frame_duration */
  int       frame_duration;
  /** deadtime_correction_factor */
  float     deadtime_correction_factor;
  /** physical_planes */
  short int physical_planes[8];
  /** fill */
  short int fill_cti[83];
  /** fill */
  short int fill_user[50];
} ECAT7_2Dscanheader;

/** ECAT7 2D normalization header */
typedef struct ecat7_2Dnormheader {
  /** data_type */
  short int data_type;
  /** num_dimensions */
  short int num_dimensions;
  /** num_r_elements */
  short int num_r_elements;
  /** num_angles */
  short int num_angles;
  /** num_z_elements */
  short int num_z_elements;
  /** ring_difference */
  short int ring_difference;
  /** scale_factor */
  float     scale_factor;
  /** norm_min */
  float     norm_min;
  /** norm_max */
  float     norm_max;
  /** fov_source_width */
  float     fov_source_width;
  /** norm_quality_factor */
  float     norm_quality_factor;
  /** norm_quality_factor_code */
  short int norm_quality_factor_code;
  /** storage_order */
  short int storage_order;
  /** span */
  short int span;
  /** z_elements */
  short int z_elements[64];
  /** fill */
  short int fill_cti[123];
  /** fill */
  short int fill_user[50];
} ECAT7_2Dnormheader;

/** ECAT7 attenuation header */
typedef struct ecat7_attenheader {
  /** data_type */
  short int data_type;
  /** num_dimensions */
  short int num_dimensions;
  /** attenuation_type */
  short int attenuation_type;
  /** Total elements collected (x dim) */
  short int num_r_elements;
  /** Total views collected (y dim) */
  short int num_angles;
  /** Total elements collected (z dim) */
  short int num_z_elements;
  /** Max acceptance angle */
  short int ring_difference;
  /** Resolution in the x dimension (cm) */
  float     x_resolution;
  /** Resolution in the y dimension (cm) */
  float     y_resolution;
  /** Resolution in the z dimension (cm) */
  float     z_resolution;
  /** TBD */
  float     w_resolution;
  /** Attenuation scale factor */
  float     scale_factor;
  /** Ellipse offset in x axis from center (cm) */
  float     x_offset;
  /** Ellipse offset in y axis from center (cm) */
  float     y_offset;
  /** Ellipse radius in x axis (cm) */
  float     x_radius;
  /** Ellipse radius in y axis (cm) */
  float     y_radius;
  /** Tilt angle of the ellipse (degrees) */
  float     tilt_angle;
  /** Mu-absorption coefficient (cm-1) */
  float     attenuation_coeff;
  /** attenuation_min */
  float     attenuation_min;
  /** attenuation_max */
  float     attenuation_max;
  /** skull_thickness */
  float     skull_thickness;
  /** num_additional_atten_coeff */
  short int num_additional_atten_coeff;
  /** additional_atten_coeff */
  float     additional_atten_coeff[8];
  /** edge_finding_threshold */
  float     edge_finding_threshold;
  /** Data storage order (RThetaZD, RZThetaD) */
  short int storage_order;
  /** Axial compression specifier (nr of ring differences spanned by a segment) */
  short int span;
  /** Nr of 'planes' in z direction for each ring difference segment */
  short int z_elements[64];
  /** fill */
  short int fill_cti[86];
  /** fill */
  short int fill_user[50];
} ECAT7_attenheader;

/** ECAT7 3D normalization header */
typedef struct ecat7_normheader {
  /** data_type */
  short int data_type;
  /** Total elements collected (y dimension) */
  short int num_r_elements;
  /** Transaxial crystals per block */
  short int num_transaxial_crystals;
  /** Nr of crystal rings */
  short int num_crystal_rings;
  /** crystals_per_ring */
  short int crystals_per_ring;
  /** Nr of rows in plane geometric correction array */
  short int num_geo_corr_planes;
  /** Upper energy limit */
  short int uld;
  /** Lower energy limit */
  short int lld;
  /** Scatter energy threshold */
  short int scatter_energy;
  /** norm_quality_factor */
  float     norm_quality_factor;
  /** norm_quality_factor_code */
  short int norm_quality_factor_code;
  /** ring_dtcor1 */
  float     ring_dtcor1[32];
  /** ring_dtcor2 */
  float     ring_dtcor2[32];
  /** crystal_dtcor */
  float     crystal_dtcor[8];
  /** Nr of ring differences included in each segment */
  short int span;
  /** Maximum ring difference acquired */
  short int max_ring_diff;
  /** fill */
  short int fill_cti[48];
  /** fill */
  short int fill_user[50];
} ECAT7_normheader;

/** ECAT7 polarmap header */
typedef struct ecat7_polmapheader {
  /** data_type */
  short int data_type;
  /** polar_map_type */
  short int polar_map_type;         
  /** num_rings */
  short int num_rings;
  /** sectors_per_ring */
  short int sectors_per_ring[32];   
  /** ring_position */
  float     ring_position[32];           
  /** ring_angle */
  short int ring_angle[32];      
  /** start_angle */
  short int start_angle;   
  /** long_axis_left */
  short int long_axis_left[3];         
  /** long_axis_right */
  short int long_axis_right[3];   
  /** position_data */
  short int position_data;     
  /** image_min */
  short int image_min; 
  /** image_max */
  short int image_max;
  /** scale_factor */
  float     scale_factor;
  /** pixel_size */
  float     pixel_size;   
  /** frame_duration */
  int       frame_duration;
  /** frame_start_time */
  int       frame_start_time;        
  /** processing_code */
  short int processing_code;          
  /** quant_units */
  short int quant_units;
  /** annotation */
  char      annotation[40];
  /** gate_duration */
  int       gate_duration;
  /** r_wave_offset */
  int       r_wave_offset;
  /** num_accepted_beats */
  int       num_accepted_beats;
  /** polar_map_protocol */
  char      polar_map_protocol[20];
  /** database_name */
  char      database_name[30];   
  /** fill */
  short int fill_cti[27];       
  /** fill */
  short int fill_user[27];
} ECAT7_polmapheader;

/** ECAT7 matrix directory */
typedef struct {
  /** id */
  int id;
  /** start block */
  int strtblk;
  /** end block */
  int endblk;
  /** status */
  int status;
} ECAT7_MatDir;

/** ECAT7 matrix list */
typedef struct {
  /** matrixNr */
  int matrixNr;
  /** matrixSpace */
  int matrixSpace;
  /** matdir */
  ECAT7_MatDir *matdir;
} ECAT7_MATRIXLIST;

/** ECAT7 matrix values */
typedef struct {
  /** frame */
  int frame;
  /** plane */
  int plane;
  /** gate */
  int gate;
  /** data */
  int data;
  /** bed */
  int bed;
} ECAT7_Matval;
/*****************************************************************************/

/*****************************************************************************/
/* ECAT 6.3 */

/** ECAT 6.3 matrix block size */
#define MatBLKSIZE 512
/** ECAT 6.3 first directory block */
#define MatFirstDirBlk 2

/** ECAT 6.3 data type */
#define	BYTE_TYPE   1
/** ECAT 6.3 data type */
#define	VAX_I2      2
/** ECAT 6.3 data type */
#define VAX_I4      3
/** ECAT 6.3 data type */
#define VAX_R4      4
/** ECAT 6.3 data type */
#define IEEE_R4     5
/** ECAT 6.3 data type */
#define	SUN_I2      6
/** ECAT 6.3 data type */
#define	SUN_I4      7

/** ECAT 6.3 file type */
#define	RAW_DATA    1
/** ECAT 6.3 file type */
#define	IMAGE_DATA  2
/** ECAT 6.3 file type */
#define	ATTN_DATA   3
/** ECAT 6.3 file type */
#define	NORM_DATA   4

/** Default system type in ECAT 6.3 file header */
#define ECAT63_SYSTEM_TYPE_DEFAULT 931
/** Error message from ECAT 6.3 functions */
char ecat63errmsg[128];
/** Verbose prints from ECAT 6.3 functions */
int ECAT63_TEST;

/** ECAT 6.3 matrix directory */
typedef struct {
  /** matnum */
  int matnum;
  /** start block */
  int strtblk;
  /** end block */
  int endblk;
  /** matstat */
  int matstat;
} MatDir;

/** ECAT 6.3 matrix list */
typedef struct {
  /** matrixNr */
  int matrixNr;
  /** matrixSpace */
  int matrixSpace;
  /** matdir */
  MatDir *matdir;
} MATRIXLIST;

/** ECAT 6.3 matrix value */
typedef struct {
  /** frame */
  int frame;
  /** plane */
  int plane;
  /** gate */
  int gate;
  /** data */
  int data;
  /** bed */
  int bed;
} Matval;

/** CTI ECAT 6.3 main header */
typedef struct ecat63_mainheader {
  /** User reserved space */
  char      ecat_format[14];
  /** User reserved space */
  char      fill1[14];
  /** Scan file's creation name */
  char      original_file_name[20];
  /** Enumerated type (VER_PRE5, VER_5, etc) */
  short int sw_version;
  /** Enumerated type (DTYPE_BYTES, DTYPE_12, etc) */
  short int data_type;
  /** Enumerated type (MODEL_911_01, MODEL_911_02, etc) */
  short int system_type;
  /** Enumerated type (FTYPE_SCAN, FTYPE_IMAGE, etc) */
  short int file_type;
  /** Unique ID of the ECAT system used */
  char      node_id[10];
  /** Day acquisition was started */
  short int scan_start_day;
  /** Month acquisition was started */
  short int scan_start_month;
  /** Year acquisition was started */
  short int scan_start_year;
  /** Hour acquisition was started */
  short int scan_start_hour;
  /** Minute acquisition was started */
  short int scan_start_minute;
  /** Second acquisition was started */
  short int scan_start_second;
  /** Isotope specifier */
  char      isotope_code[8];
  /** Half-life of isotope specified (in sec) */
  float     isotope_halflife;
  /** Radiopharmaceutical (free format ASCII) */
  char      radiopharmaceutical[32];
  /** Gantry tilt angle in degrees */
  float     gantry_tilt;
  /** Gantry rotation angle in degrees */
  float     gantry_rotation;
  /** Bed height from the lowest point (in cm) */
  float     bed_elevation;
  /** Rot source speed (revolutions per min, 0 if not rotating) */
  short int rot_source_speed;
  /** Wobble speed (revolutions per min, 0 if not wobbled) */
  short int wobble_speed;
  /** Transmission source type; enumerated type (SRC_NONE, SRC_RRS, etc) */
  short int transm_source_type;
  /** Distance from first to last plane (in cm) */
  float     axial_fov;
  /** Diameter of transaxial view (in cm) */
  float     transaxial_fov;
  /** Transaxial sampling mode; enumerated type 
      (XSAMP_STAT, XSAMP_STAT_3D, etc) */
  short int transaxial_samp_mode;
  /** Enumerated type (CSAMP_NET_TRUES, etc) */
  short int coin_samp_mode;
  /** Enumerated type (ASAMP_NORM, ASAMP_NORM_2X, etc) */
  short int axial_samp_mode;
  /** Quantification scale factor */
  float     calibration_factor;
  /** Calibration units; enumerated type (UNIT_UCIML, etc) */
  short int calibration_units;
  /** Enumerated type (COMP_NONE, etc) */
  short int compression_code;
  /** Study descriptor */
  char      study_name[12];
  /** Patient identification descriptor */
  char      patient_id[16];
  /** Patient name */
  char      patient_name[32];
  /** Patient sex (SEX_MALE, SEX_FEMALE) */
  char      patient_sex;
  /** Patient age (free format) */
  char      patient_age[10];
  /** Patient height (free format) */
  char      patient_height[10];
  /** Patient weight (free format) */
  char      patient_weight[10];
  /** Patient dexterity (DEXT_RT, DEXT_LF, DEXT_AMB) */
  char      patient_dexterity;
  /** Physician name (free format) */
  char      physician_name[32];
  /** Operator name (free format) */
  char      operator_name[32];
  /** Study description (free format) */
  char      study_description[32];
  /** Acquisition type (ACQ_RECTTR, ACQ_DYEM, etc) */
  short int acquisition_type;
  /** Bed type (BED_CTI, BED_SIEMENS) */
  short int bed_type;
  /** Septa type (SEPTA_NONE, SEPTA_3MM, etc) */
  short int septa_type;
  /** Name of facility (free format) */
  char      facility_name[20];
  /** Number of planes of data collected 
      (not necessarily the number saved in file) */
  short int num_planes;
  /** Number of frames of data collected 
      (not necessarily the number saved in file) */
  short int num_frames;
  /** Number of gates of data collected 
      (not necessarily the number saved in file) */
  short int num_gates;
  /** Number of bed positions of data collected 
      (not necessarily the number saved in file) */
  short int num_bed_pos;
  /** Absolute bed location of bed position 0 (in cm) */
  float     init_bed_position;
  /** Bed offset from init_bed_position (in cm) */
  float     bed_offset[15];
  /** Distance between adjacent planes (in cm) */
  float     plane_separation;
  /** Lowest threshold setting for scatter (in KeV) */
  short int lwr_sctr_thres;
  /** Lower threshold setting for trues (in KeV) */
  short int lwr_true_thres;
  /** Upper threshold setting for trues (in KeV) */
  short int upr_true_thres;
  /** Collimator position, if applicable */
  float     collimator;
  /** Data processing code (defined by user) */
  char      user_process_code[10];
  /** User reserved space */
  short int fill2[20];
} ECAT63_mainheader;

/** ECAT 6.3 image header */
typedef struct ecat63_imageheader {
  /** User reserved space (126 bytes) */
  char      fill1[126];
  /** Data type (DTYPE_BYTES, DTYPE_I2, etc) */
  short int data_type;
  /** Number of dimensions */
  short int num_dimensions;
  /** Unused (2 bytes) */
  short int unused1;
  /** Dimension along x axis */
  short int dimension_1;
  /** Dimension along y axis */
  short int dimension_2;
  /** Unused (24 bytes) */
  short int unused2[12];
  /** Offset in x axis for recon target (in cm) */
  float     x_origin;
  /** Offset in y axis for recon target (in cm) */
  float     y_origin;
  /** Reconstruction magnification factor (zoom) */
  float     recon_scale;
  /** Quantification scale factor (in quant_units) */
  float     quant_scale;
  /** Image minimum pixel value */
  short int image_min;
  /** Image maximum pixel value */
  short int image_max;
  /** Unused (4 bytes) */
  short int unused3[2];
  /** Pixel size (in cm) */
  float     pixel_size;
  /** Axial slice thickness (in cm) */
  float     slice_width;
  /** Total duration of current frame (in msec) */
  int       frame_duration;
  /** Frame start time as offset from the first frame */
  int       frame_start_time;
  /** Location offset from initial bed position (in cm) */
  short int slice_location;
  /** Hour when reconstruction began */
  short int recon_start_hour;
  /** Minute when reconstruction began */
  short int recon_start_min;
  /** Second when reconstruction began */
  short int recon_start_sec;
  /** Duration of reconstruction (in msec) */
  int       recon_duration;
  /** Unused (24 bytes) */
  short int unused4[12];
  /** Enumerated filter code (FILT_NONE, FILT_RAMP, etc) */
  short int filter_code;
  /** File index to corresponding scan data */
  int       scan_matrix_num;
  /** File index to corresponding normalization data */
  int       norm_matrix_num;
  /** File index to attenuation correction data */
  int       atten_cor_mat_num;
  /** Unused (46 bytes) */
  short int unused5[23];
  /** Angle image was rotated in reconstruction (in degrees) */
  float     image_rotation;
  /** Plane efficiency factor applied */
  float     plane_eff_corr_fctr;
  /** Isotope decay compensation applied to data */
  float     decay_corr_fctr;
  /** Loss correction factor (dead time) applied */
  float     loss_corr_fctr;
  /** intrinsic_tilt (previously unused space) */
  float     intrinsic_tilt;
  /** Unused (60 bytes) */
  short int unused6[30];
  /** Bit encoded processing code (PROC_DECAY_MASK, etc) */
  short int processing_code;
  /** Unused (2 bytes) */
  short int unused7;
  /** Enumerated quantification units (UNIT_MCIML, UNIT_NONE, etc) */
  short int quant_units;
  /** Day image was reconstructed */
  short int recon_start_day;
  /** Month image was reconstructed */
  short int recon_start_month;
  /** Year image was reconstructed */
  short int recon_start_year;
  /** ECAT calibration factor */
  float     ecat_calibration_fctr;
  /** Well counter calibration factor */
  float     well_counter_cal_fctr;
  /** Filter cut-off frequency, DC component, ramp slope */
  float     filter_params[6];
  /** Free format annotation */
  char      annotation[40];
  /** User reserved space (52 bytes) */
  short int fill2[26];
} ECAT63_imageheader;

/** ECAT 6.3 sinogram header */
typedef struct ecat63_scanheader {
  /** User reserved space */
  char      fill1[126];
  /** Enumerated file data type */
  short int data_type;
  /** Unused (4 bytes) */
  short int unused1[2];
  /** Total views collected (y dimension) */
  short int dimension_1;
  /** Total elements collected (x dimension) */
  short int dimension_2;
  /** Smoothing; 0=not smoothed, 1= 9x9 smoothing */
  short int smoothing;
  /** Processing applied to scan data */
  short int processing_code;
  /** Unused (6 bytes) */
  short int unused2[3];
  /** Actual distance of view sample (in cm) */
  float     sample_distance;
  /** Unused (16 bytes) */
  short int unused3[8];
  /** Half-life of isotope (in sec) */
  float     isotope_halflife;
  /** Frame duration (in sec) */
  short int frame_duration_sec;
  /** Gating segment length (in msec) */
  int       gate_duration;
  /** Time from start of the first gate (in msec) */
  int       r_wave_offset;
  /** Unused (2 bytes) */
  short int unused4;
  /** Scale factor; should be 1 if data is stored in floats. */
  float     scale_factor;
  /** Unused (6 bytes) */
  short int unused5[3];
  /** Minimum value in sinogram */
  short int scan_min;
  /** Maximum value in sinogram */
  short int scan_max;
  /** Total prompts collected in this frame/gate */
  int       prompts;
  /** Total delays collected in this frame/gate */
  int       delayed;
  /** Total multiplies collected in this frame/gate */
  int       multiples;
  /** Total net trues (prompts-randoms) collected in this frame/gate */
  int       net_trues;
  /** Unused (104 bytes) */
  short int unused6[52];
  /** Total singles with loss correction factoring */
  float     cor_singles[16];
  /** Total singles without loss correction factoring */
  float     uncor_singles[16];
  /** Mean value of loss-corrected singles */
  float     tot_avg_cor;
  /** Mean value of singles (not loss corrected) */
  float     tot_avg_uncor;
  /** Measured coincidence rate from IPCP */
  int       total_coin_rate;
  /** Time offset from first frame time (in msec) */
  int       frame_start_time;
  /** Total duration of current frame (in msec) */
  int       frame_duration;
  /** Loss correction factor applied to the sinogram */
  float     loss_correction_fctr;
  /** Unused (44 bytes) */
  short int fill2[22];
} ECAT63_scanheader;

/** ECAT 6.3 normalization header */
typedef struct ecat63_normheader {
  /** User reserved space */
  char      fill1[126];
  /** Enumerated data type */
  short int data_type;
  /** Unused (4 bytes) */
  short int unused1[2];
  /** dim */
  short int dimension_1;
  /** dim */
  short int dimension_2;
  /** Unused (46 bytes) */
  short int unused2[23];
  /** Normalization scale factor, 
      may contain plane efficiency correction factor */
  float     scale_factor;
  /** Unused (12 bytes) */
  short int unused3[6];
  /** Width of normalization source (in cm) */
  float     fov_source_width;
  /** Unused (170 bytes) */
  short int unused4[85];
  /** Hour of normalization scan */
  short int norm_hour;
  /** Unused (2 bytes) */
  short int unused5;
  /** Minute of normalization scan */
  short int norm_minute;
  /** Unused (2 bytes) */
  short int unused6;
  /** Second of normalization scan */
  short int norm_second;
  /** Unused (2 bytes) */
  short int unused7;
  /** Day of normalization scan */
  short int norm_day;
  /** Unused (2 bytes) */
  short int unused8;
  /** Month of normalization scan */
  short int norm_month;
  /** Unused (2 bytes) */
  short int unused9;
  /** Year of normalization scan */
  short int norm_year;
  /** Unused (2 bytes) */
  short int unused10;
  /** Unused (116 bytes) */
  short int unused11[58];
} ECAT63_normheader;

/** ECAT 6.3 attenuation header */
typedef struct ecat63_attnheader {
  /** User reserved space */
  char      fill1[126];
  /** Enumerated data type */
  short int data_type;
  /** Enumerated attenuation type */
  short int attenuation_type;
  /** Unused (2 bytes) */
  short int unused1;
  /** Total elements collected (x dimension) */
  short int dimension_1;
  /** Total views collected (y dimension) */
  short int dimension_2;
  /** Unused (46 bytes) */
  short int unused2[23];
  /** Attenuation scale factor */
  float     scale_factor;
  /** Ellipse offset in x axis from center (in cm) */
  float     x_origin;
  /** Ellipse offset in y axis from center (in cm) */
  float     y_origin;
  /** Ellipse radius in x axis (in cm) */
  float     x_radius;
  /** Ellipse radius in y axis (in cm) */
  float     y_radius;
  /** Tilt angle of the ellipse (in degrees) */
  float     tilt_angle;
  /** Mu-absorption coefficient (in 1/cm) */
  float     attenuation_coeff;
  /** sample_distance */
  float     sample_distance;
  /** Unused (298 bytes) */
  short int unused3[149];
} ECAT63_attnheader;

/** ECAT 6.3 matrix directory node */
typedef struct matdirnode {
  /** matnum */
  int    matnum;
  /** start block */
  int    strtblk;
  /** end block */
  int    endblk;
  /** matstat */
  int    matstat;
  /** next */
  struct matdirnode *next;
} MatDirNode ;

/** ECAT 6.3 matrix directory list */
typedef struct matdirlist {
  /** nmats */
  int         nmats;
  /** first */
  MatDirNode *first;
  /** last */
  MatDirNode *last;
} MatDirList;

/** ECAT 6.3 matrix data */
typedef struct matrixdata {
  /** mat_type */
  int    mat_type;
  /** shptr */
  char  *shptr;
  /** data_ptr */
  char  *data_ptr;
  /** nviews */
  int    nviews;
  /** nelements */
  int    nelements;
  /** nblks */
  int    nblks;
  /** data_type */
  int    data_type;
} MatrixData ;

/** ECAT 6.3 matrix file */
typedef struct matrix_file {
  /** main header */
  ECAT63_mainheader *mhptr;
  /** dir list */
  MatDirList *dirlist;
  /** File pointer */
  FILE *fptr ;
} Matrix_file;
/*****************************************************************************/

/*****************************************************************************/
/* SIF */

/** Error message from SIF functions */
char siferrmsg[128];
/** Verbose prints from SIF functions */
int SIF_TEST;

/** Scan Information File Data structure */
typedef struct {
  /** Scan time */
  time_t scantime;
  /** Number of frames */
  int frameNr;
  /** Number of columns (usually 4) */
  int colNr;
  /** SIF version */
  int version;
  /** Study number */
  char studynr[MAX_STUDYNR_LEN+1];
  /** String representation of the isotope */
  char isotope_name[8];
  /** Frame start time (sec) */
  double *x1;
  /** Frame end time (sec) */
  double *x2;
  /** Prompts */
  double *prompts;
  /** Randoms */
  double *randoms;
  /** Trues = Prompts-randoms, but at least 1 */
  double *trues;
  /** Weights = (Frame duration)^2 / trues */
  double *weights;
} SIF;
/*****************************************************************************/

/*****************************************************************************/
/* IMG */

/** Definition for img struct status */
#define IMG_STATUS_UNINITIALIZED 0
/** Definition for img struct status */
#define IMG_STATUS_INITIALIZED   1
/** Definition for img struct status */
#define IMG_STATUS_OCCUPIED      2
/** Definition for img struct status */
#define IMG_STATUS_ERROR         3

/** Definition for img error status message */
#define IMG_ERR_OK      0
/** Definition for img error status message */
#define IMG_ERR_CALLING 1
/** Definition for img error status message */
#define IMG_ERR_OOM     2

/** Definition for image type */
#define IMG_TYPE_UNKNOWN  0
/** Definition for image type */
#define IMG_TYPE_IMAGE    1
/** Definition for 'image type' sinogram */
#define IMG_TYPE_RAW      2
/** Definition for image type */
#define IMG_TYPE_POLARMAP 3
/** Definition for 'image type' attenuation data */
#define IMG_TYPE_ATTN     4

/** Definition for file format */
#define IMG_UNKNOWN   0
/** Definition for file format */
#define IMG_E63       1
/** Definition for file format */
#define IMG_E7        2
/** Definition for file format */
#define IMG_E7_2D     3
/** Definition for file format */
#define IMG_POLARMAP  9
/** Definition for file format: big endian */
#define IMG_ANA       11
/** Definition for file format: little endian variant */
#define IMG_ANA_L     12 /* little endian variant */
/** Definition for file format */
#define IMG_INTERFILE 21
/** Definition for file format: dual file format */
#define IMG_NIFTI_1D  31 /* dual file format */   
/** Definition for file format: single file format */
#define IMG_NIFTI_1S  32 /* single file format */
/** Definition for file format */
#define IMG_MICROPET  41
/** Definition for file format */
#define IMG_FLAT      61
/** Definition for file format */
#define IMG_DICOM     100

/** Definition for physical decay correction */
#define IMG_DC_UNKNOWN 0
/** Definition for physical decay correction */
#define IMG_DC_CORRECTED 1
/** Definition for physical decay correction */
#define IMG_DC_NONCORRECTED 2

/** Definition for modality */
#define IMG_MODALITY_UNKNOWN 0
/** Definition for modality */
#define IMG_MODALITY_PET 1
/** Definition for modality */
#define IMG_MODALITY_MRI 2
/** Definition for modality */
#define IMG_MODALITY_CT 3
/** Definition for modality */
#define IMG_MODALITY_SPECT 4

/** Definition for scanner model (system type) */
#define SCANNER_UNKNOWN 0
/** Definition for scanner model (system type) */
#define SCANNER_ECAT931 12
/** Definition for scanner model (system type) */
#define SCANNER_ADVANCE 12096
/** Definition for scanner model (system type) */
#define SCANNER_HRPLUS 3
/** Definition for scanner model (system type) */
#define SCANNER_HRRT 4
/* these may change later */
/** Definition for scanner model (system type): nonstandard */
#define SCANNER_MRI 5
/** Definition for scanner model (system type): nonstandard */
#define SCANNER_STEVCT_PET 6
/** Definition for scanner model (system type): nonstandard */
#define SCANNER_STEVCT_CT 7
/** Definition for scanner model (system type): nonstandard */
#define SCANNER_DMI_PET 8
/* Concorde/MicropET scanners */
/** Definition for scanner model (system type): nonstandard */
#define SCANNER_PRIMATE 2000
/** Definition for scanner model (system type): nonstandard */
#define SCANNER_RODENT 2001
/** Definition for scanner model (system type): nonstandard */
#define SCANNER_MICROPET2 2002
/** Definition for scanner model (system type): nonstandard */
#define SCANNER_FOCUS_220 2500
/** Definition for scanner model (system type): nonstandard */
#define SCANNER_FOCUS_120 2501
/** Definition for scanner model (system type): nonstandard */
#define SCANNER_INVEON_DEDICATED_PET 5000
/** Definition for scanner model (system type): nonstandard */
#define SCANNER_INVEON_MM_PET 5500
/** Definition for scanner model (system type): nonstandard */
#define SCANNER_MR_PET_HEAD_INSERT 6000
/** Definition for scanner model (system type): nonstandard */
#define SCANNER_TUEBINGEN_PET_MR 8000

/** Maximum nr of rings in polar map (based on ECAT 7 polar map header) */
#define MAX_POLARMAP_NUM_RINGS 32

/** Definitions for IMG struct status message */
enum {STATUS_OK,STATUS_FAULT,STATUS_NOMEMORY,STATUS_NOFILE,STATUS_UNKNOWNFORMAT,
    STATUS_UNSUPPORTED,STATUS_MISSINGMATRIX,STATUS_NOWRITEPERM,STATUS_DISKFULL,
    STATUS_NOMATLIST,STATUS_INVALIDMATLIST,STATUS_VARMATSIZE,STATUS_NOMAINHEADER,
    STATUS_NOSUBHEADER, STATUS_NOMATRIX, STATUS_UNSUPPORTEDAXIALCOMP,
    STATUS_NOIMGDATAFILE, STATUS_NOHEADERFILE, STATUS_INVALIDHEADER,
    STATUS_NOIMGDATA, STATUS_NOSIFDATA, STATUS_WRONGSIFDATA, 
    STATUS_CANTWRITEIMGFILE, STATUS_CANTWRITEHEADERFILE, STATUS_WRONGFILETYPE, 
    STATUS_CANNOTERASE, STATUS_CANNOTREAD, STATUS_CANNOTWRITE,
    STATUS_UNSUPPORTEDPOLARMAP, STATUS_INVALIDPOLARMAP};

/** Verbose printing from IMG functions */
int IMG_TEST;

/** Struct for IMG pixel position*/
typedef struct {
  /** [1..dimx] */
  int x;
  /** [1..dimy] */
  int y;
  /** [1..dimz] */
  int z;
  /** [1..dimt] */
  int f;
} IMG_PIXEL;

/** Struct for a list of IMG pixel positions */
typedef struct {
  /** Number of stored pixels in the list */
  int        pxlNr;
  /** Number of allocated pixels in the list */
  int        _pxlNr;
  /** List of pixels */
  IMG_PIXEL *p;
} IMG_PIXELS;

/** Struct for IMG 4D pixel volume */
typedef struct {
  /** [1..dimx] */
  int x1;
  /** [1..dimx] */
  int x2;
  /** [1..dimy] */
  int y1;
  /** [1..dimy] */
  int y2;
  /** [1..dimz] */
  int z1;
  /** [1..dimz] */
  int z2;
  /** [1..dimt] */
  int f1;
  /** [1..dimt] */
  int f2;
} IMG_RANGE;

/** Struct for IMG 4D voxel volume */
typedef struct {
  /** IMG column [1..dimx] (from left to right) */
  int x;
  /** IMG row [1..dimy] (from top to bottom) */
  int y;
  /** IMG plane [1..dimz] (from up to down) */
  int z;
  /** IMG time frame [1..dimt] */
  int t;
} VOXEL_4D;

/** 4D IMG data structure for dynamic image */
typedef struct {

  /*
   *  State of image
   */
  /** Image status (note that this is different from errstatus below):
   *  IMG_STATUS_UNINITIALIZED, IMG_STATUS_INITIALIZED,
   *  IMG_STATUS_OCCUPIED, IMG_STATUS_ERROR   */
  char status;
  /** Pointer to _imgStatusMessage, describing current status */
  const char *statmsg;

  /*
   *  Information on the study
   */
  /** for calibration units see imgUnit() in img.c */
  char unit;
  /** Calibration factor (included in pixel values) */
  float calibrationFactor;
  /** study identification code, i.e. (consequential) study number */
  char studyNr[MAX_STUDYNR_LEN+1];
  /** patient name */
  char patientName[32];
  /** patient id, e.g. 311206-123H */
  char patientID[16];
  /** Name of radiopharmaceutical */
  char radiopharmaceutical[32];
  /** Half-life of isotope (sec) */
  float isotopeHalflife;
  /** Decay correction: IMG_DC_UNKNOWN (0), IMG_DC_CORRECTED (1),
   *  IMG_DC_NONCORRECTED (2) */
  char decayCorrection;
  /** Branching fraction (included in pixel values and in calibrationFactor) */
  float branchingFraction;
  /** Scan start time and date */
  time_t scanStart;
  /** Patient orientation (see ECAT 7.2 format) */
  int orientation;
  /** User process code (which may contain valid study number) */
  char userProcessCode[11];
  /** Study description (currently free text field for user to fill) */
  char studyDescription[32];

  /*
   *  Information on the image
   */
  /** IMG_TYPE_IMAGE, IMG_TYPE_RAW */
  char type;
  /** Reconstruction zoom factor */
  float zoom;
  /** Scanner axial FOV (mm) */
  float axialFOV;
  /** Scanner transaxial FOV (mm) */
  float transaxialFOV;
  /** Scanner sample distance (mm) */
  float sampleDistance;
  /** Pixel size (mm) */
  float sizex;
  /** Pixel size (mm) */
  float sizey;
  /** Pixel size (mm) */
  float sizez;
  /** Gaps between pixels in x direction (mm); negative value means overlap */
  float gapx;
  /** Gaps between pixels in y direction (mm); negative value means overlap */
  float gapy;
  /** Gaps between pixels in z direction (mm); negative value means overlap */
  float gapz;
  /** Image resolution in x direction (mm) */
  float resolutionx;
  /** Image resolution in y direction (mm) */
  float resolutiony;
  /** Image resolution in z direction (mm) */
  float resolutionz;
  /** Saved file-format specific data type; default 0 is always ok */
  int _dataType;
  /** File format: IMG_UNKNOWN, IMG_E63, IMG_E7, IMG_E7_2D, ...
      default 0 is always ok */
  int _fileFormat;
  /** Scanner type */
  int scanner;
  /** Modality */
  int modality;
  /** XForm codes Q and S (as in NIfTI-1) */
  short int xform[2];
  /** Quaternion parameters b, c, d, and x, y, z shift, and affine transform
   *  parameters for the 1st, 2nd and 3rd row, x[4], y[4], and z[4]
   *  (as in NIfTI-1) */
  float quatern[18];
  /** Matrix transformation parameters (1,1), (1,2), (1,3), (1,4), ... (3,4)
      as in ECAT 7 image subheader */
  float mt[12];
  /** IFT struct to store any additional header information */
  IFT ift;  

  /*
   *  Definitions for polar map
   */
  /** If data is not a polar map, polarmap_num_rings=0.
      If data is a polar map, polarmap_num_rings is between 1 and 
      MAX_POLARMAP_NUM_RINGS.
   */
  int polarmap_num_rings;
  /** Number of sectors in each polar map ring; defined only in polar map data.
      In case of polar map, dimz=dimy=1, dimx= sum of sectors in each ring.
      Polar map can contain dynamic data (time frames), in that case dimz>1.
   */
  int polarmap_sectors_per_ring[MAX_POLARMAP_NUM_RINGS];
  /** Polar map: fractional distance aong the long axis from base to apex,
   *  as defined in ECAT 7 header */
  float polarmap_ring_position[MAX_POLARMAP_NUM_RINGS];
  /** Polar map ring angle relative to long axis (90 degrees along cylinder,
   *  decreasing to 0 at the apex), as defined in ECAT 7 header */
  short int polarmap_ring_angle[MAX_POLARMAP_NUM_RINGS];
  /** Polar map start angle for rings, as defined in ECAT 7 header */
  short int polarmap_start_angle;

  /*
   *  Image data
   */
  /* Dimensions */
  /** Dimension of Time (t) */
  unsigned short int dimt;
  /** Dimension of Column (c/x) */
  unsigned short int dimx;
  /** Dimension of Row (r/y) */
  unsigned short int dimy;
  /** Dimension of Plane (p/z) */
  unsigned short int dimz;
/// @cond
  /** 'Hidden' pointer for actual data */
  float *_pxl;
  /** 'Hidden' pointer for actual data */
  float **_col;
  /** 'Hidden' pointer for actual data */
  float ***_row;
  /** 'Hidden' pointer for actual data */
  float ****_pln;
  /** 'Hidden' pointer for actual data */
  float *_header;
/// @endcond
  /* Pointers for data to be used */
  /** Pointer to image data in matrix format m[plane][row][col][frame] */
  float ****m;
  /** Pointer to image data in matrix format plane[plane][row][col][frame] */
  float ****plane;
  /** Pointer to image data in matrix format row[row][col][frame] */
  float ***row;
  /** Pointer to image data in matrix format column[col][frame] */
  float **column;
  /** Pointer to image data in matrix format pixel[frame] */
  float *pixel;
  /** Plane numbers (numbers need not be contiguous with each other) */
  int *planeNumber;

  /*
   *  Frame times
   */
  /** Frame start time (sec) */
  float *start;
  /** Frame end time (sec) */
  float *end;
  /** Frame mid time (sec) */
  float *mid;

  /*
   *  Frame weights
   */
  /** Weights: 0=not weighted, 1=weighted, 2=also SD known */
  char isWeight;
  /** Frame weight factor */
  float *weight;
  /** Frame S.D. for weighting */
  float *sd;
  /** Prompts / frame */
  float *prompts;
  /** Randoms (delayed) / frame */
  float *randoms;
  
  /*
   *  Decay correction factors for each frame
   */
  /** Decay correction factor for each frame; included in pixel values */
  float *decayCorrFactor;

  /** Error status: STATUS_OK, STATUS_FAULT, STATUS_NOMEMORY, etc */
  int errstatus;

} IMG;
/*****************************************************************************/

/*****************************************************************************/
/* IMG units (deprecated) */
/// @cond
#define IMGUNIT_UNKNOWN               CUNIT_UNKNOWN
#define IMGUNIT_CPS                   CUNIT_CPS
#define IMGUNIT_COUNTS                CUNIT_COUNTS
#define IMGUNIT_KBQ_PER_ML            CUNIT_KBQ_PER_ML
#define IMGUNIT_SEC_KBQ_PER_ML        CUNIT_SEC_KBQ_PER_ML
#define IMGUNIT_PER_SEC               CUNIT_PER_SEC
#define IMGUNIT_PER_MIN               CUNIT_PER_MIN
#define IMGUNIT_ML_PER_ML             CUNIT_ML_PER_ML
#define IMGUNIT_ML_PER_DL             CUNIT_ML_PER_DL
#define IMGUNIT_ML_PER_ML_PER_MIN     CUNIT_ML_PER_ML_PER_MIN
#define IMGUNIT_ML_PER_DL_PER_MIN     CUNIT_ML_PER_DL_PER_MIN
#define IMGUNIT_UNITLESS              CUNIT_UNITLESS
#define IMGUNIT_NCI_PER_ML            CUNIT_NCI_PER_ML
#define IMGUNIT_MBQ_PER_ML            CUNIT_MBQ_PER_ML
#define IMGUNIT_BQ_PER_ML             CUNIT_BQ_PER_ML
#define IMGUNIT_UCI_PER_ML            CUNIT_UCI_PER_ML
#define IMGUNIT_UMOL_PER_MIN_PER_100G CUNIT_UMOL_PER_MIN_PER_100G
#define IMGUNIT_MG_PER_MIN_PER_100G   CUNIT_MG_PER_MIN_PER_100G
/// @endcond
/*****************************************************************************/

/*****************************************************************************/
/* Definitions for VOL, struct for 3D image volume */

/** Verbose prints from Vol functions */
int VOL_TEST;

/** 3D pixel position */
typedef struct {
  /** [1..dimx] */
  int x;
  /** [1..dimy] */
  int y;
  /** [1..dimz] */
  int z;
} VOL_PIXEL;

/** 3D pixel volume */
typedef struct {
  /** [1..dimx] */
  int x1;
  /** [1..dimx] */
  int x2;
  /** [1..dimy] */
  int y1;
  /** [1..dimy] */
  int y2;
  /** [1..dimz] */
  int z1;
  /** [1..dimz] */
  int z2;
} VOL_RANGE;  

/** 3D volume data structure - 4-byte float voxels */
typedef struct {
  /** Volume status */
  char status;
  /** Pointer to _imgStatusMessage, describing current status */
  char *statmsg;
  /** Orientation */
  int orientation;
  /** Pixel size in dimension x (mm) */
  float sizex;
  /** Pixel size in dimension y (mm) */
  float sizey;
  /** Pixel size in dimension z (mm) */
  float sizez;
  /** Dimension of Column (c/x) */
  unsigned short int dimx;
  /** Dimension of Row (r/y) */
  unsigned short int dimy;
  /** Dimension of Plane (p/z) */
  unsigned short int dimz;
/// @cond
  /** Hidden pointer for actual data */
  float *_vxl;
  /** Hidden pointer for actual data */
  float *_col;
  /** Hidden pointer for actual data */
  float **_row;
  /** Hidden pointer for actual data */
  float ***_pln;
/// @endcond
  /** Pointer for data to be used */
  float ***v;
  /** Pointer for data to be used */
  float ***plane;
  /** Pointer for data to be used */
  float **row;
  /** Pointer for data to be used */
  float *column;
  /** Pointer for data to be used */
  float *voxel;
} VOL;

/** 3D volume data structure - 2-byte short int voxels */
typedef struct {
  /** Volume status */
  char status;
  /** Pointer to _imgStatusMessage, describing current status */
  char *statmsg;
  /** Orientation */
  int orientation;
  /** Pixel size in dimension x (mm) */
  float sizex;
  /** Pixel size in dimension y (mm) */
  float sizey;
  /** Pixel size in dimension z (mm) */
  float sizez;
  /** Dimension of Column (c/x) */
  unsigned short int dimx;
  /** Dimension of Row (r/y) */
  unsigned short int dimy;
  /** Dimension of Plane (p/z) */
  unsigned short int dimz;
  /** Scaling factor */
  float scale_factor;
/// @cond
  /** Hidden pointer for actual data */
  short int *_vxl;
  /** Hidden pointer for actual data */
  short int *_col;
  /** Hidden pointer for actual data */
  short int **_row;
  /** Hidden pointer for actual data */
  short int ***_pln;
/// @endcond
  /** Pointers for data to be used */
  short int ***v;
  /** plane */
  short int ***plane;
  /** row */
  short int **row;
  /** column */
  short int *column;
  /** voxel */
  short int *voxel;
} SVOL;
/*****************************************************************************/

/*****************************************************************************/
/* MicroPET */
#ifndef MAX_MICROPET_LINE_LEN
/** Max line length in micropet header */ 
#define MAX_MICROPET_LINE_LEN 1024
#endif 
/** Verbose prints from micropet functions */
int MICROPET_TEST;
/*****************************************************************************/

/*****************************************************************************/
/* NifTI */

/** NIFTI1 header size */
#define NIFTI_HEADER_SIZE 348
/** NIFTI1 header size */
#define NIFTI_HEADER_EXTENDER_SIZE 4

/** NIFTI1 units: unknown */
#define NIFTI_UNITS_UNKNOWN 0
/** NIFTI1 units: meter */
#define NIFTI_UNITS_METER   1
/** NIFTI1 units: millimetre */
#define NIFTI_UNITS_MM      2
/** NIFTI1 units: micrometer */
#define NIFTI_UNITS_MICRON  4
/** NIFTI1 units: seconds */
#define NIFTI_UNITS_SEC     8
/** NIFTI1 units: milliseconds */
#define NIFTI_UNITS_MSEC   16
/** NIFTI1 units: microseconds */
#define NIFTI_UNITS_USEC   24
/** NIFTI1 units: Hertz */
#define NIFTI_UNITS_HERTZ   32
/** NIFTI1 units: parts per million */
#define NIFTI_UNITS_PPM     40
/** NIFTI1 units: radians per second */
#define NIFTI_UNITS_RADS    48


/** NIFTI1 datatype (same as Analyze datatypes) */
#define NIFTI_DT_NONE 0
/** NIFTI1 datatype (same as Analyze datatypes) */
#define NIFTI_DT_UNKNOWN 0
/** NIFTI1 datatype 1 bit (same as Analyze datatypes) */
#define NIFTI_DT_BINARY 1
/** NIFTI1 datatype 8 bits (same as Analyze datatypes) */
#define NIFTI_DT_UNSIGNED_CHAR 2
/** NIFTI1 datatype 16 bits (same as Analyze datatypes) */
#define NIFTI_DT_SIGNED_SHORT 4
/** NIFTI1 datatype 32 bits (same as Analyze datatypes) */
#define NIFTI_DT_SIGNED_INT 8
/** NIFTI1 datatype 32 bits (same as Analyze datatypes) */
#define NIFTI_DT_FLOAT 16
/** NIFTI1 datatype 64 bits (same as Analyze datatypes) */
#define NIFTI_DT_COMPLEX 32
/** NIFTI1 datatype 64 bits (same as Analyze datatypes) */
#define NIFTI_DT_DOUBLE 64
/** NIFTI1 datatype 24 bits (same as Analyze datatypes) */
#define NIFTI_DT_RGB 128
/** NIFTI1 datatype (same as Analyze datatypes) */
#define NIFTI_DT_ALL 255
/** NIFTI1 datatype 8 bits*/
#define NIFTI_DT_SIGNED_CHAR 256
/** NIFTI1 datatype 16 bits*/
#define NIFTI_DT_UNSIGNED_SHORT 512
/** NIFTI1 datatype 32 bits*/
#define NIFTI_DT_UNSIGNED_INT 768
/** NIFTI1 datatype 64 bits*/
#define NIFTI_DT_LONG_LONG 1024
/** NIFTI1 datatype 64 bits*/
#define NIFTI_DT_UNSIGNED_LONG_LONG 1280
/** NIFTI1 datatype 128 bits*/
#define NIFTI_DT_LONG_DOUBLE 1536
/** NIFTI1 datatype 128 bits*/
#define NIFTI_DT_DOUBLE_PAIR 1792
/** NIFTI1 datatype 256 bits*/
#define NIFTI_DT_LONG_DOUBLE_PAIR 2048
/** NIFTI1 datatype 32 bits*/
#define NIFTI_DT_RGBA 2304


/** NIFTI1 intent dataset */
#define NIFTI_INTENT_NONE 0
/** NIFTI1 intent statistics */
#define NIFTI_INTENT_CORREL 2
/** NIFTI1 intent statistics */
#define NIFTI_INTENT_TTEST 3
/** NIFTI1 intent statistics */
#define NIFTI_INTENT_FTEST 4 
/** NIFTI1 intent statistics */
#define NIFTI_INTENT_ZSCORE 5
/** NIFTI1 intent statistics */
#define NIFTI_INTENT_CHISQ 6
/** NIFTI1 intent statistics */
#define NIFTI_INTENT_BETA 7
/** NIFTI1 intent statistics */
#define NIFTI_INTENT_BINOM 8
/** NIFTI1 intent statistics */
#define NIFTI_INTENT_GAMMA 9
/** NIFTI1 intent statistics */
#define NIFTI_INTENT_POISSON 10
/** NIFTI1 intent statistics */
#define NIFTI_INTENT_NORMAL 11
/** NIFTI1 intent statistics */
#define NIFTI_INTENT_FTEST_NONC 12
/** NIFTI1 intent statistics */
#define NIFTI_INTENT_CHISQ_NONC 13
/** NIFTI1 intent statistics */
#define NIFTI_INTENT_LOGISTIC 14
/** NIFTI1 intent statistics */
#define NIFTI_INTENT_LAPLACE 15
/** NIFTI1 intent statistics */
#define NIFTI_INTENT_UNIFORM 16
/** NIFTI1 intent statistics */
#define NIFTI_INTENT_TTEST_NONC 17
/** NIFTI1 intent statistics */
#define NIFTI_INTENT_WEIBULL 18
/** NIFTI1 intent statistics */
#define NIFTI_INTENT_CHI 19
/** NIFTI1 intent statistics */
#define NIFTI_INTENT_INVGAUSS 20
/** NIFTI1 intent statistics */
#define NIFTI_INTENT_EXTVAL 21
/** NIFTI1 intent statistics */
#define NIFTI_INTENT_PVAL 22
/** NIFTI1 intent statistics */
#define NIFTI_INTENT_LOGPVAL 23 
/** NIFTI1 intent statistics */
#define NIFTI_INTENT_LOG10PVAL 24
/** NIFTI1 intent other */
#define NIFTI_INTENT_ESTIMATE 1001
/** NIFTI1 intent other */
#define NIFTI_INTENT_LABEL 1002
/** NIFTI1 intent other */
#define NIFTI_INTENT_NEURONAME 1003
/** NIFTI1 intent other */
#define NIFTI_INTENT_GENMATRIX 1004
/** NIFTI1 intent other */
#define NIFTI_INTENT_SYMMATRIX 1005
/** NIFTI1 intent other */
#define NIFTI_INTENT_DISPVECT 1006
/** NIFTI1 intent other */
#define NIFTI_INTENT_VECTOR 1007
/** NIFTI1 intent other */
#define NIFTI_INTENT_POINTSET 1008
/** NIFTI1 intent other */
#define NIFTI_INTENT_TRIANGLE 1009
/** NIFTI1 intent other */
#define NIFTI_INTENT_QUATERNION 1010
/** NIFTI1 intent other */
#define NIFTI_INTENT_DIMLESS 1011

/** NIFTI1 Coordinate System: Arbitrary coordinates. */
#define NIFTI_XFORM_UNKNOWN 0
/** NIFTI1 Coordinate System: Scanner-based anatomical coordinates. */
#define NIFTI_XFORM_SCANNER_ANAT 1
/** NIFTI1 Coordinate System: Coordinates aligned to another file or "truth". */
#define NIFTI_XFORM_ALIGNED_ANAT 2
/** NIFTI1 Coordinate System: Coordinates aligned to the Talairach space. */
#define NIFTI_XFORM_TALAIRACH 3
/** NIFTI1 Coordinate System: Coordinates aligned to the MNI space. */
#define NIFTI_XFORM_MNI_152 4



/** Nifti-1 header, 348 bytes */
typedef struct {
  /** Size of the header. Must be 348 for NIFTI-1, and 540 for NIFTI-2 (byte offset 0) */
  int sizeof_hdr;
  /** Unused. Needed for compatibility with Analyze (byte offset 4) */
  char data_type[10];
  /** Unused. Needed for compatibility with Analyze (byte offset 14) */
  char db_name[18];
  /** Unused. Value 16384 needed for compatibility with Analyze (byte offset 32) */
  int extents;
  /** Unused. Needed for compatibility with Analyze (byte offset 36) */
  short int session_error;
  /** Unused. Value 'r' needed for compatibility with Analyze (byte offset 38) */
  char regular;
  /** MRI slice ordering, encoding directions(phase, frequency, slice). (byte offset 39) */
  char dim_info;

  /** Data array dimensions; dim[0] is for the nr of dimensions,
   *  1,2,3 are for space, 4 is for time, 5 is for storing multiple values
   *  at each spatiotemporal voxel. (byte offset 40) */
  short int dim[8];
  /** 1st intent parameter, dependent on intent_code (byte offset 56) */
  float intent_p1;
  /** 2nd intent parameter, dependent on intent_code (byte offset 60) */
  float intent_p2;
  /** 3rd intent parameter, dependent on intent_code (byte offset 64) */
  float intent_p3;
  /** NIFTI_INTENT_* (byte offset 68). */
  short int intent_code;
  /** Data type (byte offset 70) */
  short int datatype;
  /** Nr of bits per voxel (byte offset 72) */
  short int bitpix;
  /** First slice index (byte offset 74) */
  short int slice_start;
  /** Grid spacings starting from pixdim[1]; pixdim[0] contains orientation
   *  (byte offset 76) */
  float pixdim[8];
  /** Offset into .nii file (byte offset 108) */
  float vox_offset;
  /** Data scaling: slope (byte offset 112); pixel values should be scaled
   *  as scl_slope*x + scl_inter */
  float scl_slope;
  /** Data scaling: offset (byte offset 116); pixel values should be scaled
   *  as scl_slope*x + scl_inter */
  float scl_inter;
  /** Last slice index (byte offset 120) */
  short int slice_end;
  /** Slice timing order (byte offset 122) */
  char slice_code;
  /** Units of pixdim[1..4], combination of NIFTI_UNITS_* (byte offset 123). */
  char xyzt_units;
  /** Max display intensity (byte offset 124) */
  float cal_max;
  /** Min display intensity (byte offset 128) */
  float cal_min;
  /** Time for 1 slice (byte offset 132) */
  float slice_duration;
  /** Time axis shift (byte offset 136) */
  float toffset;
  /** Unused. Needed for compatibility with Analyze. (byte offset 140) */
  int glmax;
  /** Unused. Needed for compatibility with Analyze. (byte offset 144) */
  int glmin;

  /** Free text field for study description (byte offset 148) */
  char descrip[80];
  /** Auxiliary file name (byte offset 228) */
  char aux_file[24];
  /** Use the quaternion fields, NIFTI_XFORM_UNKNOWN, NIFTI_XFORM_SCANNER_ANAT, or
      NIFTI_XFORM_ALIGNED_ANAT. (byte offset 252) */
  short int qform_code;
  /** Use of the affine fields, NIFTI_XFORM_* code (byte offset 254) */
  short int sform_code;
  /** Quaternion b parameter (byte offset 256) */
  float quatern_b;
  /** Quaternion c parameter (byte offset 260) */
  float quatern_c;
  /** Quaternion d parameter (byte offset 264) */
  float quatern_d;
  /** Quaternion x shift (byte offset 268) */
  float qoffset_x;
  /** Quaternion y shift (byte offset 272) */
  float qoffset_y;
  /** Quaternion z shift (byte offset 276) */
  float qoffset_z;
  /** 1st row affine transformation (byte offset 280) */
  float srow_x[4];
  /** 2nd row affine transformation (byte offset 296) */
  float srow_y[4];
  /** 2rd row affine transformation (byte offset 312) */
  float srow_z[4];
  /** Name or Meaning of data (Offset 0) (byte offset 328) */
  char intent_name[16];
  /** Magic string, "ni1\0" (dual file) or "n+1\0" (single file) (byte offset 344). 
      If zero, file should be treated as Analyze.
  */
  char magic[4];
} NIFTI_1_HEADER;

/** This structure represents a 4-byte string that should follow the
    binary nifti_1_header data in a NIFTI-1 header file. 
 */
typedef struct {
 /** If the char values are {1,0,0,0}, the file is expected to contain
  *  extensions, values of {0,0,0,0} imply the file does not contain extensions.
  *  Other sequences of values are not currently defined. */
  char extension[4];
} NIFTI_EXTENDER;

/** Combination of all NIFTI headers */
typedef struct {
  /** NIfTI-1 header */
  NIFTI_1_HEADER h;
  /** Extension; obligatory in .nii file */ 
  NIFTI_EXTENDER e;
  /** Specifies whether data is stored on disk as little endian (1),
   *  or big endian (0). */
  int byte_order;
} NIFTI_DSR;
/*****************************************************************************/

/*****************************************************************************/
/* analyze */
int anaExists(const char *dbname);
int anaExistsNew(
  const char *dbname, char *hdrfile, char *imgile, char *siffile
);
int anaRemove(const char *dbname);
void anaRemoveFNameExtension(char *fname);
int anaDatabaseExists(
  const char *dbname, char *hdrfile, char *imgfile, char *siffile
);
int anaMakeSIFName(const char *dbname, char *siffile);
int anaFlipping();
int anaReadHeader(char *filename, ANALYZE_DSR *h);
int anaReadImagedata(FILE *fp, ANALYZE_DSR *h, int frame, float *data);
int anaWriteHeader(char *filename, ANALYZE_DSR *h);
int anaPrintHeader(ANALYZE_DSR *h, FILE *fp);
int anaEditHeader(ANALYZE_DSR *h, char *field, char *value);
/*****************************************************************************/

/*****************************************************************************/
/* ecat7 */

/* Read functions */
int ecat7ReadMainheader(FILE *fp, ECAT7_mainheader *h);
int ecat7ReadImageheader(FILE *fp, int blk, ECAT7_imageheader *h);
int ecat7ReadAttenheader(FILE *fp, int blk, ECAT7_attenheader *h);
int ecat7ReadPolmapheader(FILE *fp, int blk, ECAT7_polmapheader *h);
int ecat7ReadNormheader(FILE *fp, int blk, ECAT7_normheader *h);
int ecat7ReadScanheader(FILE *fp, int blk, ECAT7_scanheader *h);
int ecat7Read2DScanheader(FILE *fp, int blk, ECAT7_2Dscanheader *h);
int ecat7Read2DNormheader(FILE *fp, int blk, ECAT7_2Dnormheader *h);
int ecat7ReadMatrixdata(FILE *fp, int start_block, int block_nr,
  char *data, int dtype);
float ecat7rFloat(void *bufi, int isvax, int islittle);
int ecat7rInt(void *bufi, int isvax, int islittle);
int ecat7ReadImageMatrix(FILE *fp, int first_block, int last_block,
  ECAT7_imageheader *h, float **fdata);
int ecat7Read2DScanMatrix(FILE *fp, int first_block, int last_block,
  ECAT7_2Dscanheader *h, float **fdata);
int ecat7ReadScanMatrix(FILE *fp, int first_block, int last_block,
  ECAT7_scanheader *h, float **fdata);
int ecat7ReadPolarmapMatrix(FILE *fp, int first_block, int last_block,
  ECAT7_polmapheader *h, float **fdata);
int ecat7pxlbytes(short int data_type);

/* Matrix list functions */
void ecat7InitMatlist(ECAT7_MATRIXLIST *mlist);
void ecat7EmptyMatlist(ECAT7_MATRIXLIST *mlist);
int ecat7ReadMatlist(FILE *fp, ECAT7_MATRIXLIST *ml, int verbose);
void ecat7PrintMatlist(ECAT7_MATRIXLIST *ml);
int ecat7EnterMatrix(FILE *fp, int matrix_id, int block_nr);
int ecat7_val_to_id(int frame, int plane, int gate, int data, int bed);
void ecat7_id_to_val(int matrix_id, ECAT7_Matval *matval);
void ecat7SortMatlistByPlane(ECAT7_MATRIXLIST *ml);
void ecat7SortMatlistByFrame(ECAT7_MATRIXLIST *ml);
int ecat7CheckMatlist(ECAT7_MATRIXLIST *ml);
int ecat7DeleteLateFrames(ECAT7_MATRIXLIST *ml, int frame_nr);
int ecat7GetPlaneAndFrameNr(ECAT7_MATRIXLIST *mlist, ECAT7_mainheader *h,
  int *plane_nr, int *frame_nr);
int ecat7GetMatrixBlockSize(ECAT7_MATRIXLIST *mlist, int *blk_nr);
int ecat7GetNums(ECAT7_MATRIXLIST *ml, ECAT7_mainheader *mh, FILE *fp,
  short int *num_planes, short int *num_frames, short int *num_gates,
  short int *num_bed_pos);
int ecat7GatherMatlist(ECAT7_MATRIXLIST *ml, short int do_planes,
  short int do_frames, short int do_gates, short int do_beds);

/* Write functions */
int ecat7WriteMainheader(FILE *fp, ECAT7_mainheader *h);
int ecat7WriteImageheader(FILE *fp, int blk, ECAT7_imageheader *h);
int ecat7WriteAttenheader(FILE *fp, int blk, ECAT7_attenheader *h);
int ecat7WritePolmapheader(FILE *fp, int blk, ECAT7_polmapheader *h);
int ecat7WriteNormheader(FILE *fp, int blk, ECAT7_normheader *h);
int ecat7WriteScanheader(FILE *fp, int blk, ECAT7_scanheader *h);
int ecat7Write2DScanheader(FILE *fp, int blk, ECAT7_2Dscanheader *h);
int ecat7Write2DNormheader(FILE *fp, int blk, ECAT7_2Dnormheader *h);
int ecat7WritePolarmapMatrix(FILE *fp, int matrix_id,
  ECAT7_polmapheader *h, float *fdata);
int ecat7WriteMatrixdata(FILE *fp, int start_block, char *data,
  int pxl_nr, int pxl_size);
/*void ecat7wFloat(float *bufi, void *bufo, int tovax, int islittle);*/
/*void ecat7wInt(int *bufi, void *bufo, int tovax, int islittle);*/
FILE *ecat7Create(const char *fname, ECAT7_mainheader *h);
int ecat7WriteImageMatrix(FILE *fp, int matrix_id, ECAT7_imageheader *h,
  float *fdata);
int ecat7Write2DScanMatrix(FILE *fp, int matrix_id, ECAT7_2Dscanheader *h,
  float *fdata);
int ecat7WriteScanMatrix(FILE *fp, int matrix_id, ECAT7_scanheader *h,
  float *fdata);
int ecat7_is_scaling_needed(float amax, float *data, int nr);

/* Printing functions */
void ecat7PrintMainheader(ECAT7_mainheader *h, FILE *fp);
void ecat7PrintImageheader(ECAT7_imageheader *h, FILE *fp);
void ecat7PrintScanheader(ECAT7_scanheader *h, FILE *fp);
void ecat7PrintAttenheader(ECAT7_attenheader *h, FILE *fp);
void ecat7PrintPolmapheader(ECAT7_polmapheader *h, FILE *fp);
void ecat7PrintNormheader(ECAT7_normheader *h, FILE *fp);
void ecat7Print2DScanheader(ECAT7_2Dscanheader *h, FILE *fp);
void ecat7Print2DNormheader(ECAT7_2Dnormheader *h, FILE *fp);
int ecat7PrintSubheader(
  ECAT7_mainheader mh, FILE *fp, int plane, int frame, FILE *ofp);

/* Descriptive strings for printing */
char* ecat7filetype(short int file_type);
char* ecat7acquisitiontype(short int acquisition_type);
char* ecat7datatype(short int data_type);

/* Header edit functions */
int ecat7EditMHeader(ECAT7_mainheader *h, char *field, char *value, int verbose);
int ecat7EditSHeader(ECAT7_scanheader *h, char *field, char *value, int verbose);
int ecat7EditVHeader(ECAT7_imageheader *h, char *field, char *value, int verbose);
/*****************************************************************************/

/*****************************************************************************/
/* ECAT 6.3 */

/* ECAT 6.3 Read functions */
int ecat63ReadMainheader(FILE *fp, ECAT63_mainheader *h);
int ecat63ReadImageheader(FILE *fp, int blk, ECAT63_imageheader *h, int verbose, char *errmsg);
int ecat63ReadScanheader(FILE *fp, int blk, ECAT63_scanheader *h, int verbose, char *errmsg);
int ecat63ReadAttnheader(FILE *fp, int blk, ECAT63_attnheader *h, int verbose, char *errmsg);
int ecat63ReadNormheader(FILE *fp, int blk, ECAT63_normheader *h, int verbose, char *errmsg);
int ecat63ReadMatdata(FILE *fp, int strtblk, int blkNr, char *data, int dtype);
int ecat63ReadImageMatrix(
  FILE *fp, int strtblk, int lastblk, ECAT63_imageheader *h, float **f);
int ecat63ReadScanMatrix(
  FILE *fp, int strtblk, int lastblk, ECAT63_scanheader *h, float **f);
int ecat63ReadAttnMatrix(
  FILE *fp, int strtblk, int lastblk, ECAT63_attnheader *h, float **f);
float ecat63rFloat(void *bufi, int isvax, int islittle);
int ecat63rInt(void *bufi, int isvax, int islittle);
int ecat63pxlbytes(short int data_type);

/* ECAT 6.3 Matrix list functions */
void ecat63InitMatlist(MATRIXLIST *mlist);
void ecat63EmptyMatlist(MATRIXLIST *mlist);
int ecat63ReadMatlist(FILE *fp, MATRIXLIST *ml, int verbose);
void ecat63PrintMatlist(MATRIXLIST *ml);
int mat_numcod(int frame, int plane, int gate, int data, int bed);
void mat_numdoc(int matnum, Matval *matval);
int ecat63Matenter(FILE *fp, int matnum, int blkNr);
void ecat63SortMatlistByPlane(MATRIXLIST *ml);
void ecat63SortMatlistByFrame(MATRIXLIST *ml);
int ecat63CheckMatlist(MATRIXLIST *ml);
int ecat63DeleteLateFrames(MATRIXLIST *ml, int frame_nr);
int ecat63GetMatrixBlockSize(MATRIXLIST *mlist, int *blk_nr);
int ecat63GetPlaneAndFrameNr(
  MATRIXLIST *mlist, ECAT63_mainheader *h, int *plane_nr, int *frame_nr);
int ecat63GetNums(
  MATRIXLIST *ml, short int *num_planes, short int *num_frames, 
  short int *num_gates, short int *num_bed_pos);
int ecat63GatherMatlist(
  MATRIXLIST *ml, short int do_planes, short int do_frames, 
  short int do_gates, short int do_beds);

/* ECAT 6.3 Write functions */
int ecat63WriteMainheader(FILE *fp, ECAT63_mainheader *h);
int ecat63WriteImageheader(FILE *fp, int block, ECAT63_imageheader *h);
int ecat63WriteScanheader(FILE *fp, int block, ECAT63_scanheader *h);
int ecat63WriteAttnheader(FILE *fp, int block, ECAT63_attnheader *h);
int ecat63WriteNormheader(FILE *fp, int block, ECAT63_normheader *h);
FILE *ecat63Create(const char *fname, ECAT63_mainheader *h);
int ecat63WriteMatdata(
  FILE *fp, int strtblk, char *data, int pxlNr, int pxlSize);
int ecat63WriteScan(FILE *fp, int matnum, ECAT63_scanheader *h, void *data);
int ecat63WriteImage(FILE *fp, int matnum, ECAT63_imageheader *h, void *data);
int ecat63WriteNorm(FILE *fp, int matnum, ECAT63_normheader *h, void *data);
int ecat63WriteAttn(FILE *fp, int matnum, ECAT63_attnheader *h, void *data);
int ecat63WriteImageMatrix(
  FILE *fp, int matnum, ECAT63_imageheader *h, float *fdata);
int ecat63WriteScanMatrix(
  FILE *fp, int matnum, ECAT63_scanheader *h, float *fdata);
void ecat63wFloat(float *bufi, void *bufo, int tovax, int islittle);
void ecat63wInt(int *bufi, void *bufo, int tovax, int islittle);
int ecat63_is_scaling_needed(float amax, float *data, int nr);
struct tm* ecat63ScanstarttimeToTm(const ECAT63_mainheader *h, struct tm *tm);
time_t ecat63Scanstarttime(const ECAT63_mainheader *h);

/* ECAT 6.3 Printing functions */
void ecat63PrintMainheader(ECAT63_mainheader *h, FILE *fp);
void ecat63PrintImageheader(ECAT63_imageheader *h, FILE *fp);
void ecat63PrintScanheader(ECAT63_scanheader *h, FILE *fp);
void ecat63PrintAttnheader(ECAT63_attnheader *h, FILE *fp);
void ecat63PrintNormheader(ECAT63_normheader *h, FILE *fp);
int ecat6PrintSubheader(ECAT63_mainheader mh, FILE *fp,
  int plane, int frame, FILE *ofp);
char *ecat63Datatype(short int dtype);
char *ecat63Unit(short int dunit);
void float2parts(float *buf);
char* ecat63ScanstarttimeInt(const ECAT63_mainheader *h, char *buf);

/* Header edit functions */
int ecat63CopyMainheader(ECAT63_mainheader *h1, ECAT63_mainheader *h2);
int ecat63CopyScanheader(ECAT63_scanheader *h1, ECAT63_scanheader *h2);
int ecat63EditMHeader(
  ECAT63_mainheader *h, char *field, char *value, int verbose);
/*****************************************************************************/

/*****************************************************************************/
/* IMG data functions */

/* Initialization and memory handling of img data */
void imgInit(IMG *image);
void imgEmpty(IMG *image);
int imgAllocate(
  IMG *image, int planes, int rows, int columns, int frames);
int imgAllocateWithHeader(
  IMG *image, int planes, int rows, int columns, int frames, IMG *image_from);
int imgDup(IMG *img1, IMG *img2);

/* Retrieving image information */
char *imgStatus(int status_index);
void imgSetStatus(IMG *img, int status_index);
void imgInfo(IMG *image);
int imgCopyhdr(IMG *image1, IMG *image2);
int imgExtractRange(IMG *img1, IMG_RANGE r, IMG *img2);
int imgExistentTimes(IMG *img);
int imgExistentCounts(IMG *img);
/*****************************************************************************/

/*****************************************************************************/
/* IMG vs SIF */
int sif2img(
  SIF *sif, IMG *img, int copy_header, int copy_frames, int copy_counts,
  int verbose
);
int img2sif(
  IMG *img, SIF *sif, int copy_header, int copy_frames, int copy_counts,
  int verbose
);
/*****************************************************************************/

/*****************************************************************************/
/* IMG decay correction */
int imgDecayCorrection(IMG *img, int mode);
char *imgIsotope(IMG *img);
int imgSetDecayCorrFactors(IMG *image, int mode);
int imgBranchingCorrection(IMG *image, int mode, int verbose, char *status);
/*****************************************************************************/

/*****************************************************************************/
/* IMG file i/o */

/* General */
int imgRead(const char *fname, IMG *img);
//int imgReadMainHeader(const char *fname, IMG *img);
int imgWrite(const char *fname, IMG *img);
int imgReadHeader(const char *fname, IMG *img, int format);
//int imgReadNextFrame(char *fname, IMG *img);
int imgReadFrame(
  const char *fname, int frame_to_read, IMG *img, int frame_index
);
int imgWriteFrame(
  const char *fname, int frame_to_write, IMG *img, int frame_index
);
void imgFormatFromFName(IMG *img, const char *fname);
int imgFormatDetermine(
  const char *fname, char *basename, char *hdrfile,
  char *imgfile, char *siffile, int *file_format, int *scanner, int *type,
  int *modality, int verbose
);

/* ECAT 6.3 and IMG */
int ecat63ReadAllToImg(const char *fname, IMG *img);
int ecat63WriteAllImg(const char *fname, IMG *img);
int ecat63ReadPlaneToImg(const char *fname, IMG *img);
int ecat63AddImg(const char *fname, IMG *img);
void imgGetEcat63MHeader(IMG *img, ECAT63_mainheader *h);
void imgSetEcat63MHeader(IMG *img, ECAT63_mainheader *h);
int imgEcat63Supported(ECAT63_mainheader *h);
int imgGetEcat63Fileformat(ECAT63_mainheader *h);
int imgReadEcat63Header(const char *fname, IMG *img);
int imgReadEcat63FirstFrame(const char *fname, IMG *img);
int imgReadEcat63Frame(
  const char *fname, int frame_to_read, IMG *img, int frame_index);
int imgWriteEcat63Frame(
  const char *fname, int frame_to_write, IMG *img, int frame_index);
void imgSetEcat63SHeader(IMG *img, void *h);

/* ECAT 7.x and IMG */
int imgReadEcat7(const char *fname, IMG *img);
int imgWriteEcat7(const char *fname, IMG *img);
int imgWrite2DEcat7(const char *fname, IMG *img);
int imgWritePolarmap(const char *fname, IMG *img);

void imgGetEcat7MHeader(IMG *img, ECAT7_mainheader *h);
void imgSetEcat7MHeader(IMG *img, ECAT7_mainheader *h);
int imgReadEcat7Header(const char *fname, IMG *img);
int imgEcat7Supported(ECAT7_mainheader *h);
int imgReadEcat7Frame(
  const char *fname, int frame_to_read, IMG *img, int frame_index);
int imgReadEcat7FirstFrame(const char *fname, IMG *img);
int imgGetEcat7Fileformat(ECAT7_mainheader *h);
int imgWriteEcat7Frame(
  const char *fname, int frame_to_write, IMG *img, int frame_index);
void imgSetEcat7SHeader(IMG *img, void *h);

/* Analyze format and IMG */
int imgReadAnalyze(const char *dbname, IMG *img);
int imgWriteAnalyze(const char *dbname, IMG *img);
int imgReadAnalyzeHeader(const char *dbname, IMG *img);
int imgGetAnalyzeHeader(IMG *img, ANALYZE_DSR *h);
int imgSetAnalyzeHeader(
  IMG *img, const char *dbname, ANALYZE_DSR *h, float fmin, float fmax);
int imgReadAnalyzeFrame(
  const char *dbname, int frame_to_read, IMG *img, int frame_index);
int imgReadAnalyzeFirstFrame(const char *fname, IMG *img);
int imgWriteAnalyzeFrame(
  const char *fname, int frame_to_write, IMG *img, int frame_index,
  float fmin, float fmax
);

/* MicroPET format and IMG */
int imgMicropetToEcat7(char *upetname, char *ecatfile, int verbose);
int imgMicropetPETToEcat7(
  FILE *fph, FILE *fpi, char *ecatfile, int verbose);
int imgMicropetCTToEcat7(
  FILE *fph, FILE *fpi, char *ecatfile, int verbose);
int imgGetMicropetMainHeader(
  FILE *fp, IMG *img, float *calibration_factor, int verbose);
int imgGetMicropetFrameHeader(
  FILE *fp, IMG *img, int frame_index, int verbose);
int imgGetMicropetSIF(FILE *fp, SIF *sif);
int imgGetMicropetHeader(IMG *img);
int imgReadMicropetHeader(const char *dbname, IMG *img);
int imgReadMicropetFrame(const char *fname, int frame_to_read,
  IMG *img, int frame_index);
int imgReadMicropetFirstFrame(const char *fname, IMG *img);
int imgReadMicropet(const char *fname, IMG *img);

/* NifTI and IMG */
int imgReadNifti(const char *filename, IMG *img, int verbose);
int imgReadNiftiFirstFrame(const char *filename, IMG *img, int verbose);
int imgReadNiftiHeader(const char *filename, IMG *img, int verbose);
int imgGetNiftiHeader(IMG *img, NIFTI_DSR *h, int verbose);
int imgReadNiftiFrame(
  const char *filename, int frame_to_read, IMG *img, int frame_index,
  int verbose
);
int imgSetNiftiHeader(
  IMG *img, const char *dbname, NIFTI_DSR *dsr, float fmin, float fmax,
  int verbose
);
int imgWriteNiftiFrame(
  const char *dbname, int frame_to_write, IMG *img, int frame_index,
  float fmin, float fmax, int verbose
);
int imgWriteNifti(
  const char *dbname, IMG *img, int save_sif, int verbose
);
/*****************************************************************************/

/*****************************************************************************/
/* imgcomp */
int imgMatch(IMG *img1, IMG *img2, float accuracy);
int imgMatchMatrix(IMG *img1, IMG *img2, double accuracy);
int imgMatchHeader(IMG *img1, IMG *img2);
int imgMatchTransform(IMG *img1, IMG *img2);
int imgMatchFrames(IMG *img1, IMG *img2);
int imgMatchPlanes(IMG *img1, IMG *img2);
int imgMaxDifference(
  IMG *img1, IMG *img2,
  VOXEL_4D *absdiff, float *abs_max, VOXEL_4D *reldiff, float *rel_max
);
/*****************************************************************************/

/*****************************************************************************/
/* imgminmax */
int imgMax(IMG *img, float *maxvalue);
int imgAbsMax(IMG *img, float *maxvalue);
int imgRangeMinMax(
  IMG *img, IMG_RANGE *r, IMG_PIXEL *maxp, float *maxv,
  IMG_PIXEL *minp, float *minv
);
int imgMinMax(IMG *img, float *minvalue, float *maxvalue);
int imgFrameMinMax(IMG *img, int frame, float *minvalue, float *maxvalue);
int imgReadMinMax(const char *fname, float *fmin, float *fmax);
int imgSmoothMax(IMG *img, float *maxvalue, IMG_PIXEL *p);
int imgGetPeak(IMG *img, float beforeTime, IMG_PIXEL *p, int verbose);
int imgGetMaxFrame(IMG *img, IMG *mimg, int verbose);
int imgGetMaxTime(IMG *img, IMG *mimg, const int w, int verbose);
int imgAvg(IMG *img, IMG_RANGE *r, float *avg);
float f_kth_smallest(float *data, int n, int k);
float fmedian(float *data, int n);
float fmean(float *data, int n, float *sd);
void fMinMaxFin(float *data, int n, float *fmin, float *fmax);
/*****************************************************************************/

/*****************************************************************************/
/* imgunits */
int imgUnitId(char *unit);
void imgUnitFromEcat(IMG *img, int ecat_unit);
void imgUnitFromEcat7(IMG *img, ECAT7_mainheader *h);
int imgUnitToEcat6(IMG *img);
void imgUnitToEcat7(IMG *img, ECAT7_mainheader *h);
char *imgUnit(int dunit);
int imgSetUnit(IMG *img, char *unit);
/*****************************************************************************/

/*****************************************************************************/
/* interfile */
int interfile_read(char headerName[256], char searchWord[256],
  char returnValue[256], char errorMessage[300]);
int interfileIsHeader(const char *hdrfile, char *imgfile);
int interfileExists(
  const char *fname, char *hdrfile, char *imgfile, int verbose);
/*****************************************************************************/

/*****************************************************************************/
/* micropet */
int upetHeaderReadParameter(FILE *fp, char *parameter, char *value);
int upetIsHeader(char *hdrfile);
int upetExists(const char *upetname, char *hdrfile, char *imgfile, int verbose);
int upetGetImageDimensions(FILE *fp, int *z, int *x, int *y, int *f);
int upetScanStart(FILE *fp, time_t *scant);
int upetReadImagedata(FILE *fp, IFT *ift, int frame, float *data);
/*****************************************************************************/

/*****************************************************************************/
/* nifti */
int niftiExists(
  const char *dbname, char *hdrfile, char *imgile, char *siffile,
  NIFTI_DSR *header, int verbose, char *status
);
int niftiCreateFNames(
  const char *filename, char *hdrfile, char *imgfile, char *siffile,
  int fileformat
);
int niftiRemove(
  const char *dbname, int fileformat, int verbose
);
int niftiReadHeader(
  char *filename, NIFTI_DSR *h, int verbose, char *status
);
int niftiPrintHeader(NIFTI_DSR *h, FILE *fp);
void niftiRemoveFNameExtension(char *fname);
int niftiReadImagedata(
  FILE *fp, NIFTI_DSR *h, int frame, float *data, int verbose, char *status
);
int niftiWriteHeader(
  char *filename, NIFTI_DSR *dsr, int verbose, char *status
);
/*****************************************************************************/

/*****************************************************************************/
/* sif */
//void libsif_printdate(FILE *fp);
int sifRead(char *filename, SIF *data);
int sifWrite(SIF *data, char *filename);
void sifPrint(SIF *data);
void sifEmpty(SIF *data);
void sifInit(SIF *data);
int sifSetmem(SIF *data, int frameNr);
void sifWeight(SIF *data, double halflife);
void sifWeightByFrames(SIF *data, double halflife);
void sifWeightNorm(SIF *data);
void sifModerateTrues(SIF *sif, double limit);
void sifModerateWeights(SIF *sif, double limit);
int sifExistentCounts(SIF *sif);
/* Deprecated */
/// @cond
#define readSIF sifRead
#define writeSIF sifWrite
#define printSIF sifPrint
#define emptySIF sifEmpty
#define weightSIF sifWeight
#define initSIF sifInit
#define setmemSIF sifSetmem
/// @endcond
/*****************************************************************************/

/*****************************************************************************/
/* ecat7ift */
/** Definitions for ECAT matrix contents */
typedef struct {
  /** Matrix number */
  int mnum;
  /** Struct containing frame, plane, gate, data, bed */
  ECAT7_Matval matval;
  /** Subheader */
  IFT sh;
  /** Pointer to pixel data */
  float *f;
} ECAT_MATRIX;

/** Definitions for ECAT headers contents */
typedef struct {
  /** Main header */
  IFT mh;
  /** Number of ECAT matrices */
  int nr;
  /** Pointer to the list of ECAT matrices */
  ECAT_MATRIX *m;
} ECAT_HEADERS;
/* Functions */
void ematInitiate(ECAT_MATRIX *emat);
void ehdrInitiate(ECAT_HEADERS *ehdr);
void ematEmpty(ECAT_MATRIX *emat);
void ehdrEmpty(ECAT_HEADERS *ehdr);
int ehdrAllocate(ECAT_HEADERS *ehdr, int nr);
int ecat7MHeaderToIFT(ECAT7_mainheader *h, IFT *ift, int verbose);
int ecat7MainheaderFromIFT(ECAT7_mainheader *h, IFT *ift, int verbose);
int ecat7ImageheaderToIFT(ECAT7_imageheader *h, IFT *ift, int verbose);
int ecat7ScanheaderToIFT(ECAT7_scanheader *h, IFT *ift, int verbose);
int ecat7ReadSubheaderToIFT(FILE *fp, ECAT7_mainheader *h, int strtblk, 
  IFT *ift, int verbose);
int ecat7WriteSubheaderFromIFT(FILE *fp, ECAT7_mainheader *h, int strtblk, 
  IFT *ift, int verbose);
int ecat7ReadHeaders(const char *fname, ECAT_HEADERS *ehdr, int verbose);
int ecat7WriteHeaders(const char *fname, ECAT_HEADERS *ehdr, int verbose);
/*****************************************************************************/

/*****************************************************************************/
/* ird */
int string_to_xyzf(const char *str, IMG_PIXEL *v);
int irdReorder(IMG_RANGE *img_range);
int irdRead(char *irdfile, IMG_RANGE *img_range, char *status);
int irdCheck(IMG_RANGE *r, IMG *img);
/*****************************************************************************/

/*****************************************************************************/
/* vol */
/* 4-byte floats */
void volInit(VOL *vol);
void volEmpty(VOL *vol);
int volAllocate(VOL *vol, int planes, int rows, int columns);
int img2vol(IMG *img, VOL *vol, int frame);
int vol2img(VOL *vol, IMG *img, int frame);
void volInfo(VOL *vol, FILE *fp);
void volContents(VOL *vol, VOL_RANGE r, FILE *fp);
int volMax(
  VOL *vol, VOL_RANGE *r,
  VOL_PIXEL *maxp, float *maxv, VOL_PIXEL *minp, float *minv
);
int volAvg(VOL *vol, VOL_RANGE *r, float *avg);
/* 2-byte short ints */
void svolInit(SVOL *svol);
void svolEmpty(SVOL *svol);
int svolAllocate(SVOL *svol, int planes, int rows, int columns);
int img2svol(IMG *img, SVOL *svol, int frame);
int svol2img(SVOL *svol, IMG *img, int frame);
void svolInfo(SVOL *svol, FILE *fp);
/* General */
int vrdReorder(VOL_RANGE *vol_range);
int vrdVxlNr(VOL_RANGE *vol_range);
int vrd2vol(VOL_RANGE *r, VOL *vol, float in, float out, char *status);
int vrdRead(char *vdffile, VOL_RANGE *vol_range, char *status);
int string_to_xyz(char *str, int *x, int *y, int *z);
/*****************************************************************************/

/*****************************************************************************/
/* pixel */
void pxlInit(IMG_PIXELS *pxl);
void pxlFree(IMG_PIXELS *pxl);
int pxlAllocate(IMG_PIXELS *pxl, int pxlNr);
int pxlAllocateMore(IMG_PIXELS *pxl, int pxlNr);
int pxlMakeRoom(IMG_PIXELS *list, int i, int n);
int pxlAdd(IMG_PIXELS *list, IMG_PIXEL *pxl);
int pxlGet(IMG_PIXELS *list, int i, IMG_PIXEL *pxl);
int pxlAddFromMask(IMG_PIXELS *list, IMG *img);
void pxlMove(IMG_PIXELS *list, int from, int to);
int pxlRm(IMG_PIXELS *list, int index);
int pxlRmDuplicates(IMG_PIXELS *list);
int pxlWrite(IMG_PIXELS *pxl, FILE *fp, char *status);
int pxlRead(IMG_PIXELS *pxl, const char *fname, char *status);
/*****************************************************************************/

/*****************************************************************************/
/* DICOM */

/** DICOM tag. Tag group and element are shown in hex format. */
typedef struct DCMTAG {
  /** Group: even numbers defined in DICOM standard, odd numbers are
      vendor specific.
  */
  unsigned short int group;
  /** Element */
  unsigned short int element;
} DCMTAG;

/** @public DICOM Transfer Syntax UID. 

    In case of implicit VR, elements do NOT contain VR !
    Reference: DICOM PS3.5 2017a chapter 10.

    Items must be the same and in the same order as the dcm_truid list.
*/
typedef enum {
  DCM_TRUID_UNKNOWN, ///< Unknown Transfer Syntax UID  
  DCM_TRUID_LEI,     ///< Little Endian Implicit VR (DICOM default)
  DCM_TRUID_LEE,     ///< Little Endian Explicit VR
  DCM_TRUID_BEE,     ///< Big Endian Explicit VR
  DCM_TRUID_JPEG50,  ///< Lossy JPEG 8-bit compression
  DCM_TRUID_JPEG51,  ///< Lossy JPEG 12-bit compression
  DCM_TRUID_JPEG70,  ///< Lossless JPEG
  DCM_TRUID_JPEG80,  ///< Lossless JPEG-LS
  DCM_TRUID_JPEG81,  ///< Lossy JPEG-LS
  DCM_TRUID_JPEG90,  ///< Lossless JPEG 2000
  DCM_TRUID_JPEG91,  ///< JPEG 2000
  DCM_TRUID_JPEG92,  ///< Lossless multicomponent JPEG 2000
  DCM_TRUID_JPEG93,  ///< Multicomponent JPEG 2000
  DCM_TRUID_MPEG100, ///< MPEG-2
  DCM_TRUID_MPEG102, ///< MPEG-4
  DCM_TRUID_MPEG103, ///< MPEG-4 BD-compatible
  DCM_TRUID_RLE,     ///< Lossless RLE
  DCM_TRUID_RFC,     ///< RFC 2557
  DCM_TRUID_XML,     ///< XML encoding
  DCM_TRUID_INVALID  ///< Invalid Transfer Syntax UID
} dcmtruid;

/** @public DICOM value representation (VR). 

    Reference: DICOM PS3.5 2017a chapter 6.2.
 
    Items must be the same and in the same order as the dcm_vr list.
*/
typedef enum {
  DCM_VR_AE,        ///< DICOM application entity, max 16 bytes.
  DCM_VR_AS,        ///< DICOM age string, 4 bytes fixed.
  DCM_VR_AT,        ///< DICOM attribute tag, 4 bytes fixed.
  DCM_VR_CS,        ///< DICOM code (control) string, max 16 bytes.
  DCM_VR_DA,        ///< DICOM date in format YYYYMMDD, 8 bytes fixed. @note In old standard 10 bytes fixed, in format YYYY.MM.DD
  DCM_VR_DS,        ///< DICOM decimal string, max 16 bytes.
  DCM_VR_DT,        ///< DICOM date time, max 26 bytes.
  DCM_VR_FL,        ///< DICOM floating point single precision, 4 bytes fixed.
  DCM_VR_FD,        ///< DICOM floating point double precision, 8 bytes fixed.
  DCM_VR_IS,        ///< DICOM integer string, max 12 bytes.
  DCM_VR_LO,        ///< DICOM long string, max 64 chars.
  DCM_VR_LT,        ///< DICOM long text, max 10240 chars.
  DCM_VR_OB,        ///< DICOM other byte string, even bytes, endian insensitive.
  DCM_VR_OD,        ///< DICOM other double (64-bit) stream, endian sensitive.
  DCM_VR_OF,        ///< DICOM other float (32-bit) stream, endian sensitive.
  DCM_VR_OL,        ///< DICOM other long (32-bit) stream, endian sensitive.
  DCM_VR_OW,        ///< DICOM other word (16-bit) stream, even bytes, endian sensitive.
  DCM_VR_PN,        ///< DICOM person name, max 64 chars per component group.
  DCM_VR_SH,        ///< DICOM short string, max 16 chars.
  DCM_VR_SL,        ///< DICOM signed long (32-bit integer), 4 bytes fixed.
  DCM_VR_SQ,        ///< DICOM sequence of zero or more elements (used for nested data).
  DCM_VR_SS,        ///< DICOM signed short (16-bit integer), 2 bytes fixed.
  DCM_VR_ST,        ///< DICOM short text, max 1024 chars.
  DCM_VR_TM,        ///< DICOM time HHMMSS.FFFFFF, max 14 bytes. @note In old standard 16 bytes max.
  DCM_VR_UC,        ///< DICOM unlimited characters.
  DCM_VR_UI,        ///< DICOM unique identifier (UID), max 64 bytes.
  DCM_VR_UL,        ///< DICOM unsigned long (32-bit) integer, 4 bytes fixed.
  DCM_VR_UN,        ///< DICOM unknown, any valid length of another VR.
  DCM_VR_UR,        ///< DICOM URI or URL, string of characters.
  DCM_VR_US,        ///< DICOM unsigned short (16-bit) integer, 2 bytes fixed.
  DCM_VR_UT,        ///< DICOM unlimited text, character string.
  DCM_VR_INVALID    ///< Invalid DICOM value representation.
} dcmvr;

/** Data struct for one DICOM item; may be recursive. */
typedef struct DCMITEM {
  /** Nr of File pointer, NULL if not opened. */
  FILE *fp;
  /** File position of the start of this element; 0, if not set. */
  fpos_t pos;
  /** Enumerated Transfer Syntax UID. */
  dcmtruid truid;
  /** Item tag. */
  DCMTAG tag;
  /** Enumerated Value Representation. */
  dcmvr vr;
  /** Value Length. */
  unsigned int vl;
  /** Pointer to linked list of child elements; NULL if none. */
  struct DCMITEM *child_item;
  /** Pointer to linked list of parent elements; NULL if none. */
  struct DCMITEM *parent_item;
  /** Pointer to next item ; NULL if none. */
  struct DCMITEM *next_item;
  /** Pointer to previous item ; NULL if none. */
  struct DCMITEM *prev_item;
  /** Pointer to raw data value (no byte conversions etc); 
      NULL if not available. */
  char *rd;
} DCMITEM;

/** Main data struct for one DICOM file */
typedef struct DCMFILE {
  /** DICOM filename. */
  char filename[FILENAME_MAX];
  /** File pointer, NULL if not opened. */
  FILE *fp;
  /** Enumerated Transfer Syntax UID. */
  dcmtruid truid;
  /** Pointer to linked list of DICOM elements; recursive; NULL if none. */
  DCMITEM *item;
} DCMFILE;
 
/** Data struct for one DICOM image matrix.
   @sa DCMML
 */
typedef struct DCMMATRIX {
  /** File name; must be allocated and freed as necessary;
      note that one file may contain many matrices. */
  char *filename;
  /** Acquisition date. */ 
  char acqDate[16];
  /** Acquisition time. */ 
  char acqTime[16];
  /** Frame [1..frameNr]; note that one frame may contain several bed positions,
      each with their own frame times, but those will be on different planes. */
  unsigned int frame;
  /** Plane [1..planeNr]. */
  unsigned int plane;
  /** Frame start time (sec). */
  double frameStart;
  /** Frame duration (sec). */
  double frameDur;
} DCMMATRIX;

/** Main data struct for all DICOM matrices. */
typedef struct DCMML {
  /** Nr of matrices. */
  unsigned int nr;
  /** Allocate size of matrix list. */
  unsigned int anr;
  /** Pointer to matrix list. */
  DCMMATRIX *m;
} DCMML;
/*****************************************************************************/

/*****************************************************************************/
/* dcm */
int dcmVerifyMagic(const char *filename, FILE *fp);
unsigned char dcmVRReserved(dcmvr id);
dcmvr dcmVRId(const char *s);
char *dcmVRName(dcmvr id);
size_t dcmVRVLength(dcmvr id);
char *dcmVRDescr(dcmvr id);
char *dcmDA2intl(const char *orig, char *intl);
char *dcmTM2intl(const char *orig, char *intl);
char *dcmDT2intl(const char *orig, char *intl);
unsigned int dcmSOPIdentify(const char *s);
char *dcmSOPName(unsigned int i);
char *dcmSOPUID(unsigned int i);
char *dcmSOPUIDName(const char *s);
dcmtruid dcmTrUID(const char *s);
char *dcmTrUIDDescr(dcmtruid id);
char *dcmTrUIDString(dcmtruid id);
dcmtruid dcmReadTransferSyntaxUID(FILE *fp);
int dcmReadFileTag(FILE *fp, DCMTAG *tag);
int dcmWriteFileTag(FILE *fp, DCMTAG *tag);
int dcmWriteFileSQDelimItem(FILE *fp);
int dcmWriteFileSQItemDelimTag(FILE *fp);
dcmvr dcmReadFileVR(FILE *fp, char *vrstr);
unsigned int dcmReadFileVL(FILE *fp, unsigned int n);
int dcmReadFileVRVL(FILE *fp, dcmvr *vr, unsigned int *vl, unsigned int *n);
int dcmWriteFileVRVL(FILE *fp, dcmvr vr, unsigned int vl, unsigned int *n);

void dcmfileInit(DCMFILE *d);
void dcmitemFree(DCMITEM *d);
void dcmfileFree(DCMFILE *d);
unsigned short int dcmfileMaxDepth(DCMFILE *df);
unsigned short int dcmitemMaxDepth(DCMITEM *d);
unsigned short int dcmitemParentNr(DCMITEM *d);
char *dcmValueString(DCMITEM *d);
long int dcmitemGetInt(DCMITEM *d);
double dcmitemGetReal(DCMITEM *d);
DCMITEM *dcmFindTag(DCMITEM *d, const short int omit, DCMTAG *tag, const int verbose);
void dcmitemPrint(DCMITEM *d);
void dcmTagSet(DCMTAG *tag, unsigned short int group, unsigned short int element);
int dcmAddItem(
  DCMFILE *dcm, DCMITEM *d, short int aschild, DCMTAG tag, dcmvr vr, unsigned int vl,
  char *rd, const int verbose
);

int dcmFileReadNextElement(
  DCMFILE *dcm, DCMITEM *prev_item, DCMITEM *parent_item, const short int sub,
  const short int headerOnly, int verbose
);
int dcmFileRead(
  const char *filename, DCMFILE *dcm, const short int headerOnly,
  int verbose
);
int dcmFileWrite(const char *filename, DCMFILE *dcm, int verbose);

/*****************************************************************************/

/*****************************************************************************/
#endif // _LIBTPCIMGIO_H

#ifdef __cplusplus
}
#endif