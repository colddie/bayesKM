/// @file libtpcsvg.h
/// @brief Header file for libtpcsvg.
/// @author Vesa Oikonen
///

#ifdef __cplusplus
extern "C" {
#endif


#ifndef _LIBTPCSVG_H
#define _LIBTPCSVG_H
/*****************************************************************************/

/*****************************************************************************/
#include "tpcclibConfig.h"
/*****************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
/*****************************************************************************/

/*****************************************************************************/
/** Width of SVG viewport */
#define SVG_VIEWPORT_WIDTH 10000
/** Height of SVG viewport */
#define SVG_VIEWPORT_HEIGHT 10000
/** Max nr of axis ticks */
#define MAX_TICK_NR 12
/** Max length of axis tick labels */
#define MAX_TICK_LABEL_LEN 20
/** Max length of SVG legend texts */
#define MAX_SVG_LEGEND_LEN 30
/*****************************************************************************/

/*****************************************************************************/
/* axis */
int axis_tick_positions(
  const double begin, const double end, double *ticks, int *tick_nr, 
  double *scale_factor, int *tick_decimals, int verbose
);
void axis_check_range(double *begin, double *end, int verbose);
void strRmExpZeroes(char *str);
/*****************************************************************************/

/*****************************************************************************/
/** SVG plot symbol type */
typedef enum {RECTANGLE,CIRCLE,UPTRIANGLE,DOWNTRIANGLE,DIAMOND,
              LEFTTRIANGLE,RIGHTTRIANGLE} svgSymbolType;
/** SVG plot symbol open/filled */
typedef enum {SYMBOLOPEN,SYMBOLFILLED} svgSymbolFill;
/** SVG plot colour */
typedef enum {BLACK,RED,BLUE,GREEN,PURPLE,OLIVE,AQUA,FUCHSIA,GRAY,LIME,MAROON,
              NAVY,SILVER,TEAL,YELLOW} svgColor;
/*****************************************************************************/
/** Position of viewport
  @sa svg_init_viewport_pos
 */
struct svg_viewport_pos {
  /** Viewport exists? */
  int is;
  /** x pos */
  int x;
  /** y pos */
  int y;
  /** width */
  int w;
  /** height */
  int h;
  /** character size */
  int chr_size;
};

/** SVG plot coordinates
  @sa svg_init_coord
 */
struct svg_coord {
  /** min */
  double min;
  /** max */
  double max;
  /** scale */
  double scale;
  /** origo */
  double origo;
  /** Nr of ticks */
  int tick_nr;
  /** tick value */
  double tick[MAX_TICK_NR];
  /** scale of tick values */
  double tickscale;
  /** nr of decimals in tick value */
  int tick_decimals;
  /** tick labels */
  char tick_label[MAX_TICK_NR][MAX_TICK_LABEL_LEN+1];
  /** upper margin */
  int upper_margin;
  /** is min fixed? */
  int fixed_min;
  /** is max fixed? */
  int fixed_max;
};
/** Viewport for plotting data in SVG format.
  @sa svg_init_viewports, svg_define_viewports
 */
struct svg_viewports {
  /** Main viewport */
  struct svg_viewport_pos main_viewport;
  /** Main title; relative to the main viewport */
  struct svg_viewport_pos main_title_viewport;
  /** Y axis title; relative to the main viewport */
  struct svg_viewport_pos yaxis_title_viewport;
  /** X axis title; relative to the main viewport */
  struct svg_viewport_pos xaxis_title_viewport;
  /** Legend area; relative to the main viewport */
  struct svg_viewport_pos label_area_viewport;
  /** Plot area; relative to the main viewport */
  struct svg_viewport_pos plot_area_viewport;
  /** Coordinate area; relative to the plot_area */
  struct svg_viewport_pos coordinate_area_viewport;
  /** X coordinate area contents */
  struct svg_coord x;
  /** Y coordinate area contents */
  struct svg_coord y;
};
/*****************************************************************************/
/** Struct for storing data for one legend item */
typedef struct svg_legend {
  /** Plot type: 1=line, 2=symbols, 0=both line and symbols */
  int plot_type;
  /** Symbol type: RECTANGLE,CIRCLE,UPTRIANGLE,DOWNTRIANGLE,DIAMOND,
                   LEFTTRIANGLE, RIGHTTRIANGLE */
  svgSymbolType symbol_type;
  /** Symbol filling: SYMBOLOPEN, SYMBOLFILLED */
  svgSymbolFill symbol_fill;
  /** SVG color index */
  svgColor color;
  /** Legend text */
  char text[MAX_SVG_LEGEND_LEN+1];
} SVG_LEGEND;

/** Struct for storing all legends */
typedef struct svg_legends {
/// @cond
  /** Switch to show if struct contents are initiated (1) or not (0) */
  int _init; 
/// @endcond
  /** Nr of legends */
  int n;
  /** Pointer to one legend */
  SVG_LEGEND *l;
} SVG_LEGENDS;
/*****************************************************************************/

/*****************************************************************************/
/* svg_file.c */
FILE *svg_initiate(
  const char *filename, const double height, const double width, 
  struct svg_viewports *vp, char *errmsg, int verbose
);
int svg_close(
  FILE *fp, char *errmsg, int verbose
);
FILE *svg_xhtml_initiate(
  const char *filename, const char *XHTML_title, char *errmsg, int verbose
);
int svg_xhtml_close(
  FILE *fp, char *errmsg, int verbose
);
int svg_xhtml_svg_open(
  FILE *fp, const double height, const double width, struct svg_viewports *vp,
  char *errmsg, int verbose
);
int svg_xhtml_svg_close(
  FILE *fp, char *errmsg, int verbose
);
int svg_write(
  FILE *fp, const char *svg_string, char *errmsg, int verbose
);
/*****************************************************************************/
/* svg_vport.c */
void svg_init_viewport_pos(struct svg_viewport_pos *p);
void svg_init_coord(struct svg_coord *p);
void svg_init_viewports(struct svg_viewports *p);
int svg_define_viewports(
  const int main_viewport_width, const int main_viewport_height,
  const int is_main_title, const int is_yaxis_title, const int is_xaxis_title, 
  const int is_label_area, struct svg_viewports *vp, int verbose
);
/*****************************************************************************/
/* svg_plot.c */
int svg_start_plot_viewport(
  FILE *fp, struct svg_viewports *vp, char *errmsg, int verbose
);
int svg_end_plot_viewport(
  FILE *fp, /*struct svg_viewports *vp,*/ char *errmsg, int verbose
);
int svg_start_coordinate_viewport(
  FILE *fp, struct svg_viewports *vp, char *errmsg, int verbose
);
int svg_end_coordinate_viewport(
  FILE *fp, char *errmsg, int verbose
);
int svg_calculate_axes(
  struct svg_viewports *vp, int verbose
);
int svg_write_axes(
  FILE *fp, struct svg_viewports *vp, char *errmsg, int verbose
);
int svg_write_xticks(
  FILE *fp, struct svg_viewports *vp, char *errmsg, int verbose
);
int svg_write_yticks(
  FILE *fp, struct svg_viewports *vp, char *errmsg, int verbose
);
int svg_write_tac(
  FILE *fp, struct svg_viewports *vp, const int plot_type, const char *tac_id, 
  const char *tac_title, double *x, double *y, const int data_nr, 
  const char *color, const svgSymbolType symbol_type, const svgSymbolFill symbol_fill,
  char *errmsg, int verbose
);
int get_line_intersection(
  const double a1x, const double a1y, const double a2x, const double a2y,
  const double b1x, const double b1y, const double b2x, const double b2y,
  double *ix, double *iy, int verbose
);
int check_intersection_with_viewport(
  const double x1, const double y1, const double x2, const double y2,
  struct svg_viewport_pos *cavp,
  double *nx1, double *ny1, double *nx2, double *ny2, int verbose
);
/*****************************************************************************/
/* svg_title.c */
int svg_create_main_title(
  FILE *fp, const char *main_title_text, const char *sub_title_text,
  struct svg_viewports *vp, char *errmsg, int verbose
);
int svg_create_xaxis_title(
  FILE *fp, const char *title_text, struct svg_viewports *vp,
  char *errmsg, int verbose
);
int svg_create_yaxis_title(
  FILE *fp, const char *title_text, struct svg_viewports *vp,
  char *errmsg, int verbose
);
/*****************************************************************************/
/* svg_defs.c */
int svg_define_symbols(FILE *fp, char *errmsg, int verbose);
char *svgColorName(const svgColor index);
char *svgSymbolName(const svgSymbolType index);
/*****************************************************************************/
/* svg_legends.c */
void svg_init_legends(
  SVG_LEGENDS *legends
);
void svg_legend_empty(
  SVG_LEGENDS *legends
);
int svg_legend_add(
  SVG_LEGENDS *legends, const int plot_type, const int symbol_type, const svgSymbolFill symbol_fill, 
  const int color, const char *text
);
int svg_create_legends(
  FILE *fp, struct svg_viewports *vp, SVG_LEGENDS *legends,
  char *errmsg, int verbose
);
/*****************************************************************************/

/*****************************************************************************/
#endif

#ifdef __cplusplus
}
#endif
