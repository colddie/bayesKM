/// @file svg_vport.c
/// @author Vesa Oikonen
/// @brief Create SVG plot viewports for TPC SVG C library.
///
/*****************************************************************************/
#include "libtpcsvg.h"
#include <math.h>
/*****************************************************************************/
/** Write inline SVG (1) or separate SVG file (0) */
extern int SVG_INLINE;
/*****************************************************************************/

/*****************************************************************************/
/** Initiate struct svg_viewport_pos contents to all-zeroes before use. */
void svg_init_viewport_pos(
  struct svg_viewport_pos *p
) {
  p->is=0;
  p->x=0;
  p->y=0;
  p->w=0;
  p->h=0;
  p->chr_size=0;
}
/** Initiate struct svg_coord contents to all-zeroes before use. */
void svg_init_coord(
  struct svg_coord *p
) {
  p->min=0;
  p->max=0;
  p->scale=0.0;
  p->origo=0.0;
  p->tick_nr=0;
  for(int i=0; i<MAX_TICK_NR; i++) p->tick[i]=0.0;
  p->tickscale=0.0;
  p->tick_decimals=0.0;
  for(int i=0; i<MAX_TICK_NR; i++) p->tick_label[i][0]=(char)0;
  p->upper_margin=0;
  p->fixed_min=0;
  p->fixed_max=0;
}
/** Initiate struct svg_viewports contents to all-zeroes before use. */
void svg_init_viewports(
  struct svg_viewports *p
) {
  svg_init_viewport_pos(&p->main_viewport);
  svg_init_viewport_pos(&p->main_title_viewport);
  svg_init_viewport_pos(&p->yaxis_title_viewport);
  svg_init_viewport_pos(&p->xaxis_title_viewport);
  svg_init_viewport_pos(&p->label_area_viewport);
  svg_init_viewport_pos(&p->plot_area_viewport);
  svg_init_viewport_pos(&p->coordinate_area_viewport);
  svg_init_coord(&p->x);
  svg_init_coord(&p->y);
}
/*****************************************************************************/

/*****************************************************************************/
/** Define the viewport positions for further use. All measures are in pixels.

    Axis tick labels (y axis at least) should be set before calling this, so that
    enough room for the labels can be reserved.
   @return Returns pointer to the file if successful and NULL in case of an error.
 */
int svg_define_viewports(
  /** Main viewport width (zero if default is used) */
  const int main_viewport_width,
  /** Main viewport height (zero if default is used) */
  const int main_viewport_height,
  /** Is there main title? no=0, yes<>0 */
  const int is_main_title,
  /** Is there y axis title? no=0, yes<>0 */
  const int is_yaxis_title,
  /** Is there x axis title? no=0, yes<>0 */
  const int is_xaxis_title,
  /** Is there label area? no=0, yes<>0 */
  const int is_label_area,
  /** Pointer to structure which will be filled with viewport positions and sizes */
  struct svg_viewports *vp,
  /** Verbose level; if zero, then nothing is printed to stderr or stdout */
  int verbose
) {
  int ti, m, n;

  if(verbose>0)
    printf("svg_define_viewports(%d, %d, %d, %d, %d, %d, vp, %d)\n",
      main_viewport_width, main_viewport_height, is_main_title, is_yaxis_title,
      is_xaxis_title, is_label_area, verbose);

  if(vp==NULL) return(1);
  /* Set main viewport */
  vp->main_viewport.is=1; vp->main_viewport.x=vp->main_viewport.y=0;
  if(main_viewport_width<1) vp->main_viewport.w=SVG_VIEWPORT_WIDTH;
  if(main_viewport_height<1) vp->main_viewport.h=SVG_VIEWPORT_HEIGHT;
  /* Set viewport for main title(s) */
  if(is_main_title==0) {
    vp->main_title_viewport.is=0;
    vp->main_title_viewport.x=vp->main_title_viewport.y=0;
    vp->main_title_viewport.w=vp->main_title_viewport.h=0;
    vp->main_title_viewport.chr_size=0;
  } else {
    vp->main_title_viewport.is=1;
    vp->main_title_viewport.x=vp->main_title_viewport.y=0;
    vp->main_title_viewport.w=vp->main_viewport.w;
    vp->main_title_viewport.h=vp->main_viewport.w/12;
    vp->main_title_viewport.chr_size=5*vp->main_title_viewport.h/10;
  }
  /* Set viewport for x axis title */
  if(is_xaxis_title==0) {
    vp->xaxis_title_viewport.is=0;
    vp->xaxis_title_viewport.x=0;
    vp->xaxis_title_viewport.y=vp->main_title_viewport.h;
    vp->xaxis_title_viewport.w=vp->main_title_viewport.w;
    vp->xaxis_title_viewport.h=0;
    vp->xaxis_title_viewport.chr_size=0;
  } else {
    vp->xaxis_title_viewport.is=1;
    vp->xaxis_title_viewport.x=0;
    vp->xaxis_title_viewport.w=vp->main_viewport.w;
    vp->xaxis_title_viewport.h=vp->main_viewport.h/18;
    vp->xaxis_title_viewport.y=vp->main_viewport.h-vp->xaxis_title_viewport.h;
    vp->xaxis_title_viewport.chr_size=7*vp->xaxis_title_viewport.h/10;
  }
  /* Set viewport for y axis title */
  if(is_yaxis_title==0) {
    vp->yaxis_title_viewport.is=0;
    vp->yaxis_title_viewport.x=0;
    vp->yaxis_title_viewport.y=vp->main_title_viewport.h;
    vp->yaxis_title_viewport.w=0;
    vp->yaxis_title_viewport.h=vp->main_viewport.h-vp->main_title_viewport.h-vp->xaxis_title_viewport.h;
    vp->yaxis_title_viewport.chr_size=0;
  } else {
    vp->yaxis_title_viewport.is=1;
    vp->yaxis_title_viewport.x=0;
    vp->yaxis_title_viewport.y=vp->main_title_viewport.h;
    vp->yaxis_title_viewport.w=vp->main_viewport.w/18;
    vp->yaxis_title_viewport.h=vp->main_viewport.h-vp->main_title_viewport.h-vp->xaxis_title_viewport.h;
    if(vp->xaxis_title_viewport.is)
      vp->yaxis_title_viewport.chr_size=vp->xaxis_title_viewport.chr_size;
    else
      vp->yaxis_title_viewport.chr_size=7*vp->yaxis_title_viewport.w/10;
  }
  /* Set viewport for label area */
  if(is_label_area==0) {
    vp->label_area_viewport.is=0;
    vp->label_area_viewport.x=vp->main_viewport.w;
    vp->label_area_viewport.y=vp->main_title_viewport.h;
    vp->label_area_viewport.w=0;
    vp->label_area_viewport.h=vp->main_viewport.h-vp->main_title_viewport.h-vp->xaxis_title_viewport.h;
  } else {
    vp->label_area_viewport.is=1;
    vp->label_area_viewport.x=3*vp->main_viewport.w/4;
    vp->label_area_viewport.y=vp->main_title_viewport.h;
    vp->label_area_viewport.w=vp->main_viewport.w-vp->label_area_viewport.x;
    vp->label_area_viewport.h=vp->main_viewport.h-vp->main_title_viewport.h-vp->xaxis_title_viewport.h;
  }
  /* Set viewport for plot area */
  vp->plot_area_viewport.is=1;
  vp->plot_area_viewport.x=vp->yaxis_title_viewport.w;
  vp->plot_area_viewport.y=vp->main_title_viewport.h;
  vp->plot_area_viewport.w=vp->main_viewport.w-vp->yaxis_title_viewport.w-vp->label_area_viewport.w;
  vp->plot_area_viewport.h=vp->main_viewport.h-vp->main_title_viewport.h-vp->xaxis_title_viewport.h;
  /* Set plot area contents (inside plot area) */
  vp->coordinate_area_viewport.is=1;
  for(ti=m=0; ti<vp->y.tick_nr; ti++) {n=strlen(vp->y.tick_label[ti]); if(n>m) m=n;}
  if(verbose>2) printf("max_yaxis_label_len=%d\n", m);
  if(m<3) vp->coordinate_area_viewport.x=vp->plot_area_viewport.w/14;
  else if(m<5) vp->coordinate_area_viewport.x=vp->plot_area_viewport.w/10;
  else if(m<7) vp->coordinate_area_viewport.x=vp->plot_area_viewport.w/8;
  else vp->coordinate_area_viewport.x=vp->plot_area_viewport.w/5;
  vp->coordinate_area_viewport.y=0;
  vp->coordinate_area_viewport.w=vp->plot_area_viewport.w-vp->coordinate_area_viewport.x;
  vp->coordinate_area_viewport.h=19*vp->plot_area_viewport.h/20;
  /* Calculate the character size for tick labels etc */
  vp->plot_area_viewport.chr_size=vp->coordinate_area_viewport.chr_size=
    ceil(0.67*(double)(vp->plot_area_viewport.h-vp->coordinate_area_viewport.h));

  if(verbose>3) printf("coordinate_area_viewport.h := %d\n", vp->coordinate_area_viewport.h);

  return(0);
}
/*****************************************************************************/

/*****************************************************************************/
