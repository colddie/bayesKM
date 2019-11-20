/// @file svg_plot.c
/// @author Vesa Oikonen
/// @brief Create SVG plot contents for TPC SVG C library.
///
/*****************************************************************************/
#include "libtpcsvg.h"
/*****************************************************************************/
/** Write inline SVG (1) or separate SVG file (0) */
extern int SVG_INLINE;
/*****************************************************************************/

/*****************************************************************************/
/** Check whether two lines, each drawn between two points, intersect each other.

    If either end of lines intersects with the other line, that is NOT counted as intersection.
   @return Returns 1 in case of intersection, and 0 if they do not cross.
   @sa check_intersection_with_viewport
 */ 
int get_line_intersection(
  /** x,y coordinates of line a at points 1 and 2; x coordinate of point 1. */
  const double a1x,
  /** y coordinate of point 1. */
  const double a1y,
  /** x coordinate of point 2. */
  const double a2x,
  /** y coordinate of point 2. */
  const double a2y,
  /** x,y coordinates of line b at points 1 and 2; x coordinate of point 1. */
  const double b1x,
  /** y coordinate of point 1. */
  const double b1y,
  /** x coordinate of point 2. */
  const double b2x,
  /** y coordinate of point 2. */
  const double b2y,
  /** Pointers for intersection point coordinates; NULL if not needed. */
  double *ix,
  /** Intersection y coordinate */
  double *iy,
  /** Verbose level; if zero, then nothing is printed to stderr or stdout. */
  int verbose
) {
  double sax, say, sbx, sby, s, t, d;

  if(verbose>0) {
    printf("get_line_intersection(%g, %g, %g, %g, %g, %g, %g, %g, ix, iy, %d)\n",
      a1x, a1y, a2x, a2y, b1x, b1y, b2x, b2y, verbose);  
    fflush(stdout);
  }

  if(a1x==a2x && a1y==a2y) return(0);
  sax=a2x-a1x; say=a2y-a1y;
  sbx=b2x-b1x; sby=b2y-b1y;
  d=-sbx*say+sax*sby; if(d==0.0) return(0);
  s=(-say*(a1x-b1x) + sax*(a1y-b1y)) / d;
  t=(+sbx*(a1y-b1y) - sby*(a1x-b1x)) / d;
//if(s>=0.0 && s<=1.0 && t>=0.0 && t<=1.0) {
  if(s>0.0 && s<1.0 && t>0.0 && t<1.0) {
    if(verbose>3) printf("s=%g t=%g\n", s, t);
    if(ix!=NULL) *ix=a1x+(t*sax);
    if(iy!=NULL) *iy=a1y+(t*say);
    return(1);
  } 
  return(0);
}
/*****************************************************************************/

/*****************************************************************************/
/** Check the intersections between specified line and viewport borders.
   @return Returns the number (0-2) of crossings.
   @sa get_line_intersection
 */
int check_intersection_with_viewport(
  /** x,y coordinates of line at points 1 and 2; x coordinate of point 1. */
  const double x1, 
  /** y coordinate of point 1. */
  const double y1, 
  /** x coordinate of point 2. */
  const double x2,
  /** y coordinate of point 2. */
  const double y2,
  /** Pointer to coordinate area viewport. */
  struct svg_viewport_pos *cavp,
  /** Pointers for (possibly) modified line coordinates; NULL if not needed. */
  double *nx1,
  /** new y coordinate of point 1. */
  double *ny1,
  /** new x coordinate of point 2. */
  double *nx2,
  /** new y coordinate of point 2. */
  double *ny2,
  /** Verbose level; if zero, then nothing is printed to stderr or stdout. */
  int verbose
) {
  int i=0, n, crossings=0;
  double ix, iy;
  double mx[2], my[2];

  if(verbose>0) {
    printf("check_intersection_with_viewport(%g, %g, %g, %g, cavp, nx1, ny1, nx2, ny2, %d)\n",
      x1, y1, x2, y2, verbose); 
    fflush(stdout);
  }

  mx[0]=x1; my[0]=y1; mx[1]=x2; my[1]=y2;
  /* Upper border */
  n=get_line_intersection(x1, y1, x2, y2, 0.0, 0.0, (double)cavp->w, 0.0, &ix, &iy, verbose);
  if(n>0) {
    if(crossings==0) {
      if(y1<iy) i=0; else i=1; // which point was out? change that one
    } else { // both are out, change the other than first time
      if(i==0) i=1; else i=0;
    }
    mx[i]=ix; my[i]=iy; crossings++;
    if(verbose>3)
      printf("line between (%g,%g) and (%g,%g) would cross upper border at (%g,%g)\n",
        x1, y1, x2, y2, ix, iy);
  }
  /* Lower border */
  n=get_line_intersection(x1, y1, x2, y2,
          0.0, (double)cavp->h, (double)cavp->w, (double)cavp->h, &ix, &iy, verbose);
  if(n>0) {
    if(crossings==0) {
      if(y1>iy) i=0; else i=1; // which point was out? change that one
    } else { // both are out, change the other than first time
      if(i==0) i=1; else i=0;
    }
    mx[i]=ix; my[i]=iy; crossings++;
    if(verbose>3)
      printf("line between (%g,%g) and (%g,%g) would cross lower border at (%g,%g)\n",
        x1, y1, x2, y2, ix, iy);
  }
  /* Left border */
  n=get_line_intersection(x1, y1, x2, y2,
          0.0, 0.0, 0.0, cavp->h,
          &ix, &iy, verbose);
  if(n>0) {
    if(crossings==0) {
      if(x1<ix) i=0; else i=1; // which point was out? change that one
    } else { // both are out, change the other than first time
      if(i==0) i=1; else i=0;
    }
    mx[i]=ix; my[i]=iy; crossings++;
    if(verbose>3)
      printf("line between (%g,%g) and (%g,%g) would cross left border at (%g,%g)\n",
        x1, y1, x2, y2, ix, iy);
  }
  /* Right border */
  n=get_line_intersection(x1, y1, x2, y2,
          (double)cavp->w, 0.0, (double)cavp->w, (double)cavp->h,
          &ix, &iy, verbose);
  if(n>0) {
    if(crossings==0) {
      if(x1>ix) i=0; else i=1; // which point was out? change that one
    } else { // both are out, change the other than first time
      if(i==0) i=1; else i=0;
    }
    mx[i]=ix; my[i]=iy; crossings++;
    if(verbose>3)
      printf("line between (%g,%g) and (%g,%g) would cross right border at (%g,%g)\n",
        x1, y1, x2, y2, ix, iy);
  }
  if(verbose>3 && crossings>0) printf("crossings=%d\n", crossings);

  if(nx1!=NULL) *nx1=mx[0];
  if(ny1!=NULL) *ny1=my[0];
  if(nx2!=NULL) *nx2=mx[1];
  if(ny2!=NULL) *ny2=my[1];
  if(verbose>2 && crossings>0)
    printf("modified line (%g,%g) -> (%g,%g)\n", mx[0], my[0], mx[1], my[1]);

  return(crossings);
}
/*****************************************************************************/

/*****************************************************************************/
/** Start plot area viewport.
   @return Returns 0 if successful, <>0 in case of error.
   @sa svg_end_plot_viewport, svg_define_viewports
*/
int svg_start_plot_viewport(
  /** SVG graphics file pointer. */
  FILE *fp,
  /** Struct containing the viewport sizes. */
  struct svg_viewports *vp,
  /** Char pointer to string (at least of length 128) where possible
      error description is copied; set to NULL if not necessary. */
  char *errmsg,
  /** Verbose level; if zero, then nothing is printed to stderr or stdout. */
  int verbose
) {
  char tmp[1024], line[128], ilc[9];

  if(verbose>0) {
    printf("svg_start_plot_viewport(fp, vp, errmsg, %d)\n", verbose);
    fflush(stdout);
  }

  /* Check the input */
  if(fp==NULL) {
    if(errmsg!=NULL) sprintf(errmsg, "file was closed too early");
    return(1);
  }
  if(SVG_INLINE) strcpy(ilc, "svg:"); else strcpy(ilc, "");

  /* Create a new viewport for plot area */
  strcpy(tmp, "\n  <!-- Plot area viewport -->\n");
  strcat(tmp, "  <"); strcat(tmp, ilc); strcat(tmp, "svg");
  sprintf(line, " x=\"%dpx\" y=\"%dpx\" width=\"%dpx\" height=\"%d\"",
    vp->plot_area_viewport.x, vp->plot_area_viewport.y,
    vp->plot_area_viewport.w, vp->plot_area_viewport.h); strcat(tmp, line);
  sprintf(line, "\n      viewBox=\"0 0 %d %d\"",
    vp->plot_area_viewport.w, vp->plot_area_viewport.h); strcat(tmp, line);
  strcat(tmp, "\n      preserveAspectRatio=\"xMidYMid meet\"");
  strcat(tmp, ">\n");
  if(svg_write(fp, tmp, errmsg, verbose-5)!=0) {return(4);}

#if(0)
  sprintf(tmp, "    <%srect width=\"%dpx\" height=\"%dpx\" stroke=\"none\" fill=\"yellow\" fill-opacity=\"0.1\" />\n",
    ilc, vp->plot_area_viewport.w, vp->plot_area_viewport.h);
  if(svg_write(fp, tmp, errmsg, verbose-5)!=0) {return(99);}
#endif

#if(0) // moved elsewhere
  /* Create symbols for later use */
  if(svg_define_symbols(fp, errmsg)!=0) {return(6);}
#endif

  return(0);
}
/*****************************************************************************/

/*****************************************************************************/
/** End plot viewport.
   @return Returns 0 if successful, <>0 in case of error.
   @sa svg_start_plot_viewport
*/
int svg_end_plot_viewport(
  /** SVG graphics file pointer */
  FILE *fp,
  // /** Struct containing the viewport sizes */
  // struct svg_viewports *vp,
  /** Char pointer to string (at least of length 128) where possible
      error description is copied; set to NULL if not necessary */
  char *errmsg,
  /** Verbose level; if zero, then nothing is printed to stderr or stdout */
  int verbose
) {
  char tmp[1024];

  if(verbose>0) {
    printf("svg_end_plot_viewport(fp, vp, errmsg, %d)\n", verbose);
    fflush(stdout);
  }

  /* Check the input */
  if(fp==NULL) {
    if(errmsg!=NULL) sprintf(errmsg, "file was closed too early");
    return(1);
  }

  /* Write the end tag */
  if(SVG_INLINE) strcpy(tmp, "  </svg:svg>\n"); else strcpy(tmp, "  </svg>\n");

  if(svg_write(fp, tmp, errmsg, verbose-5)!=0) {return(2);}

  return(0);
}
/*****************************************************************************/

/*****************************************************************************/
/** Start coordinate area viewport.
   @return Returns 0 if successful, <>0 in case of error.
   @sa svg_end_coordinate_viewport, svg_define_viewports
*/
int svg_start_coordinate_viewport(
  /** SVG graphics file pointer. */
  FILE *fp,
  /** Struct containing the viewport sizes. */
  struct svg_viewports *vp,
  /** Char pointer to string (at least of length 128) where possible
      error description is copied; set to NULL if not necessary. */
  char *errmsg,
  /** Verbose level; if zero, then nothing is printed to stderr or stdout. */
  int verbose
) {
  char tmp[1024], line[128];

  if(verbose>0) {
    printf("svg_start_coordinate_viewport(fp, vp, errmsg, %d)\n", verbose);
    fflush(stdout);
  }

  /* Check the input */
  if(fp==NULL) {
    if(errmsg!=NULL) sprintf(errmsg, "file was closed too early");
    return(1);
  }

  /* Create a new viewport for plot area */
  strcpy(tmp, "\n  <!-- Coordinate area viewport -->\n");
  if(SVG_INLINE) strcat(tmp, "  <svg:"); else strcat(tmp, "  <");
  sprintf(line, "svg x=\"%dpx\" y=\"%dpx\" width=\"%dpx\" height=\"%d\"",
    vp->coordinate_area_viewport.x, vp->coordinate_area_viewport.y,
    vp->coordinate_area_viewport.w, vp->coordinate_area_viewport.h); strcat(tmp, line);
  sprintf(line, "\n      viewBox=\"0 0 %d %d\"",
    vp->coordinate_area_viewport.w, vp->coordinate_area_viewport.h); strcat(tmp, line);
  strcat(tmp, "\n      preserveAspectRatio=\"xMidYMid meet\"");
  strcat(tmp, ">\n");
  if(svg_write(fp, tmp, errmsg, verbose-5)!=0) {return(4);}

#if(0)
  sprintf(tmp, "    <rect width=\"%dpx\" height=\"%dpx\" stroke=\"none\" fill=\"silver\" fill-opacity=\"0.3\" />\n",
    vp->coordinate_area_viewport.w, vp->coordinate_area_viewport.h);
  if(svg_write(fp, tmp, errmsg, verbose-5)!=0) {return(99);}
#endif

  return(0);
}
/*****************************************************************************/

/*****************************************************************************/
/** End coordinate area viewport.
   @return Returns 0 if successful, <>0 in case of error.
   @sa svg_start_coordinate_viewport
*/
int svg_end_coordinate_viewport(
  /** SVG graphics file pointer */
  FILE *fp,
  /** Char pointer to string (at least of length 128) where possible
      error description is copied; set to NULL if not necessary. */
  char *errmsg,
  /** Verbose level; if zero, then nothing is printed to stderr or stdout. */
  int verbose
) {
  char tmp[1024];

  if(verbose>0) {
    printf("svg_end_coordinate_viewport(fp, errmsg, %d)\n", verbose);
    fflush(stdout);
  }

  /* Check the input */
  if(fp==NULL) {
    if(errmsg!=NULL) sprintf(errmsg, "file was closed too early");
    return(1);
  }

  if(SVG_INLINE) strcpy(tmp, "  </svg:svg>\n\n");
  else strcpy(tmp, "  </svg>\n\n");
  if(svg_write(fp, tmp, errmsg, verbose-5)!=0) {return(2);}

  return(0);
}
/*****************************************************************************/

/*****************************************************************************/
/** Calculate the axis tick positions.
    Before calling this, viewport must be filled with curve min and max values.
    This routine checks that max>min, changing the values if necessary.
   @return Returns 0 if successful, <>0 in case of error.
   @sa svg_define_viewports, svg_write_axes, svg_write_xticks, svg_write_yticks, axis_tick_positions
 */
int svg_calculate_axes(
  /** Struct containing the viewport sizes. */
  struct svg_viewports *vp,
  /** Verbose level; if zero, then nothing is printed to stderr or stdout. */
  int verbose
) {
  int ret, ti, prec;
  double orig_min, orig_max;

  if(verbose>0) {
    printf("svg_calculate_axes(vp, %d)\n", verbose);
    fflush(stdout);
  }

  /* x axis */
  orig_min=vp->x.min; orig_max=vp->x.max;
  axis_check_range(&vp->x.min, &vp->x.max, verbose);
  if(vp->x.fixed_min) vp->x.min=orig_min;
  if(vp->x.fixed_max) vp->x.max=orig_max;
  if(verbose>1) printf("x-range %g - %g -> %g - %g\n", orig_min, orig_max, vp->x.min, vp->x.max);
  vp->x.tick_nr=MAX_TICK_NR;
  if(vp->label_area_viewport.is) vp->x.tick_nr=1+vp->x.tick_nr/2; // Legends reduce the width of x axis
  ret=axis_tick_positions(
    vp->x.min, vp->x.max, vp->x.tick, &vp->x.tick_nr, &vp->x.tickscale,
    &vp->x.tick_decimals, verbose);
  if(ret!=0) return(ret+100);

  /* create tick labels to be written later */
  prec=vp->x.tick_decimals-1-(int)vp->x.tickscale; if(prec<0) prec=0;
  for(ti=0; ti<vp->x.tick_nr; ti++) {
    if(vp->x.tickscale<-2 || vp->x.tickscale>3) {
      sprintf(vp->x.tick_label[ti], "%.*E", vp->x.tick_decimals-1, vp->x.tick[ti]);
      strRmExpZeroes(vp->x.tick_label[ti]);
    } else if(vp->x.tickscale<=0)
      sprintf(vp->x.tick_label[ti], "%.*f", prec, vp->x.tick[ti]);
    else
      sprintf(vp->x.tick_label[ti], "%.*f", prec, vp->x.tick[ti]);
  }

  /* y axis */
  orig_min=vp->y.min; orig_max=vp->y.max;
  axis_check_range(&vp->y.min, &vp->y.max, verbose);
  if(vp->y.fixed_min) vp->y.min=orig_min;
  if(vp->y.fixed_max) vp->y.max=orig_max;
  if(verbose>1) printf("y-range %g - %g -> %g - %g\n", orig_min, orig_max, vp->y.min, vp->y.max);
  vp->y.tick_nr=MAX_TICK_NR; // 10
  ret=axis_tick_positions(
    vp->y.min, vp->y.max, vp->y.tick, &vp->y.tick_nr, &vp->y.tickscale,
    &vp->y.tick_decimals, verbose);
  if(ret!=0) return(ret+200);
  /* create tick labels to be written later */
  prec=vp->y.tick_decimals-1-(int)vp->y.tickscale; if(prec<0) prec=0;
  for(ti=0; ti<vp->y.tick_nr; ti++) {
    if(vp->y.tickscale<-2 || vp->y.tickscale>3) {
      sprintf(vp->y.tick_label[ti], "%.*E", vp->y.tick_decimals-1, vp->y.tick[ti]);
      strRmExpZeroes(vp->y.tick_label[ti]);
    } else if(vp->y.tickscale<=0)
      sprintf(vp->y.tick_label[ti], "%.*f", prec, vp->y.tick[ti]);
    else
      sprintf(vp->y.tick_label[ti], "%.*f", prec, vp->y.tick[ti]);
  }

  return(0);
}
/*****************************************************************************/

/*****************************************************************************/
/** Draw the axes into SVG plot coordinate area.
   @return Returns 0 if successful, <>0 in case of error.
   @sa svg_define_viewports, svg_calculate_axes, svg_write_xticks, svg_write_yticks, axis_tick_positions
*/
int svg_write_axes(
  /** SVG graphics file pointer. */
  FILE *fp,
  /** Struct containing the viewport sizes. */
  struct svg_viewports *vp,
  /** Pointer to string where error message is written; NULL if not needed. */
  char *errmsg,
  /** Verbose level; if zero, then nothing is printed to stderr or stdout. */
  int verbose
) {
  int linew=20, coordw=10;
  char tmp[1024], line[128], ilc[9];
  double f;

  if(verbose>0) {
    printf("svg_write_axes(fp, vp, errmsg, %d)\n", verbose);
    fflush(stdout);
  }

  /* Check the input */
  if(fp==NULL) {
    if(errmsg!=NULL) sprintf(errmsg, "file was closed too early");
    return(1);
  }
  if(vp->x.min>=vp->x.max || vp->y.min>=vp->y.max) {
    if(errmsg!=NULL) sprintf(errmsg, "invalid plot range");
    if(verbose>1) {
      printf("vp->x.min=%g vp->x.max=%g\n", vp->x.min, vp->x.max);
      printf("vp->y.min=%g vp->y.max=%g\n", vp->y.min, vp->y.max);
    }
    return(2);
  }

  if(SVG_INLINE) strcpy(ilc, "svg:"); else strcpy(ilc, "");

  /* Draw the lines around the plot */
  sprintf(tmp, "    <%s", ilc);
  sprintf(line, "polyline fill=\"none\" stroke=\"%s\" stroke-width=\"%d\"\n", "black", linew);
  strcat(tmp, line);
#if(0)
  sprintf(line, "      points=\"%d,%d %d,%d %d,%d %d,%d %d,%d\" />\n",
    linew/2, linew/2,
    linew/2, vp->coordinate_area_viewport.h-linew,
    vp->coordinate_area_viewport.w-linew/2, vp->coordinate_area_viewport.h-linew,
    vp->coordinate_area_viewport.w-linew/2, linew/2,
    linew/2, linew/2);
  strcat(tmp, line);
#else
  sprintf(line, "      points=\"%d,%d %d,%d %d,%d %d,%d %d,%d\" />\n",
    linew/2, linew/2,
    linew/2, vp->coordinate_area_viewport.h-linew/2,
    vp->coordinate_area_viewport.w-linew/2, vp->coordinate_area_viewport.h-linew/2,
    vp->coordinate_area_viewport.w-linew/2, linew/2,
    linew/2, linew/2);
  strcat(tmp, line);
#endif
  if(svg_write(fp, tmp, errmsg, verbose-5)!=0) {return(4);}

  /* Set the upper margins for both axes */
  vp->y.upper_margin=0.04*(double)vp->coordinate_area_viewport.h;
#if(1)
  if(vp->label_area_viewport.is)
    vp->x.upper_margin=0.02*(double)vp->coordinate_area_viewport.w;
  else // If no legends, then more room may be needed to fit x tick labels
    vp->x.upper_margin=0.08*(double)vp->coordinate_area_viewport.w;
#else
  vp->x.upper_margin=0.08*(double)vp->coordinate_area_viewport.w;
#endif

  /* Calculate the scale factors */
  f=vp->coordinate_area_viewport.w-vp->x.upper_margin;
  vp->x.scale=f/(vp->x.max-vp->x.min);
  if(verbose>0) printf("xscalef:=%g (%g vs %g-%g)\n", vp->x.scale, f, vp->x.min, vp->x.max);
  f=vp->coordinate_area_viewport.h-vp->y.upper_margin;
  vp->y.scale=f/(vp->y.max-vp->y.min);
  if(verbose>1) printf("yscalef:=%g (%g vs %g-%g)\n", vp->y.scale, f, vp->y.min, vp->y.max);

  /* Calculate the origo in plot coordinates */
  vp->x.origo=-vp->x.scale*vp->x.min;
  if(verbose>1) printf("x.origo := %g\n", vp->x.origo);
  vp->y.origo=-vp->y.scale*vp->y.min;
  if(verbose>1) printf("y.origo := %g\n", vp->y.origo);

  /* Draw the x=0 line, if necessary */
  if(vp->x.origo>0 && vp->x.origo<vp->coordinate_area_viewport.w/*-vp->x.upper_margin*/) {
    if(verbose>1) printf("drawing x=0 line\n");
    sprintf(tmp, "    <%s", ilc);
    sprintf(line, "line fill=\"none\" stroke=\"%s\" stroke-width=\"%d\"\n", "black", coordw);
    strcat(tmp, line);
    sprintf(line, "      x1=\"%g\" x2=\"%g\" y1=\"%d\" y2=\"%d\" />\n",
            vp->x.origo, vp->x.origo, 0, vp->coordinate_area_viewport.h);
    strcat(tmp, line);
    if(svg_write(fp, tmp, errmsg, verbose-5)!=0) {return(6);}
  }

  /* Draw the y=0 line, if necessary */
  if(vp->y.origo>0 && vp->y.origo<vp->coordinate_area_viewport.h
     /*-vp->y.upper_margin*/) {
    if(verbose>1) printf("drawing y=0 line\n");
    sprintf(tmp, "    <%s", ilc);
    sprintf(line, "line fill=\"none\" stroke=\"%s\" stroke-width=\"%d\"\n", "black", coordw);
    strcat(tmp, line);
    sprintf(line, "      x1=\"%d\" x2=\"%d\" y1=\"%g\" y2=\"%g\" />\n",
      0, vp->coordinate_area_viewport.w,
      vp->coordinate_area_viewport.h-vp->y.origo,
      vp->coordinate_area_viewport.h-vp->y.origo);
    strcat(tmp, line);
    if(svg_write(fp, tmp, errmsg, verbose-5)!=0) {return(7);}
  }

  return(0);
}
/*****************************************************************************/

/*****************************************************************************/
/** Create SVG plot x axis ticks.
   @return Returns 0 if successful, <>0 in case of error.
   @sa svg_define_viewports, svg_write_axes, svg_write_yticks, svg_calculate_axes
*/
int svg_write_xticks(
  /** SVG graphics file pointer. */
  FILE *fp,
  /** Struct containing the viewport sizes. */
  struct svg_viewports *vp,
  /** Char pointer to string (at least of length 128) where possible
      error description is copied; set to NULL if not necessary. */
  char *errmsg,
  /** Verbose level; if zero, then nothing is printed to stderr or stdout. */
  int verbose
) {
  char tmp[1024], line[128], ilc[9];
  int ti;
  double pos, ypos, yheight;

  if(verbose>0) {printf("svg_write_xticks(fp, vp, errmsg, %d)\n", verbose); fflush(stdout);}

  /* Check the input */
  if(vp->x.tick_nr<1 || vp->plot_area_viewport.h==vp->coordinate_area_viewport.h) return(0);
  if(fp==NULL) {if(errmsg!=NULL) sprintf(errmsg, "file was closed too early"); return(1);}

  if(SVG_INLINE) strcpy(ilc, "svg:"); else strcpy(ilc, "");

  yheight=vp->plot_area_viewport.h-vp->coordinate_area_viewport.h;
  ypos=vp->coordinate_area_viewport.h;

  strcpy(tmp, "\n    <!-- X axis ticks inside plot area -->\n");
#if(0)
  sprintf(tmp, "      <%srect x=\"0px\" y=\"%dpx\" width=\"%dpx\" height=\"%gpx\" stroke=\"none\" fill=\"aqua\" fill-opacity=\"0.3\" />\n",
    ilc, vp->coordinate_area_viewport.h, vp->plot_area_viewport.w, yheight);
#endif
  if(svg_write(fp, tmp, errmsg, verbose-5)!=0) {return(99);}

  /* Write ticks */
  sprintf(tmp, "    <%sg", ilc);
  sprintf(line, " stroke=\"%s\"", "black"); strcat(tmp, line);
  sprintf(line, " stroke-width=\"%g\"", 20.); strcat(tmp, line);
  sprintf(line, " fill=\"%s\"", "none"); strcat(tmp, line);
  strcat(tmp, ">\n"); 
  if(svg_write(fp, tmp, errmsg, verbose-5)!=0) {return(5);}

  if(verbose>9) {
    printf("vp->plot_area_viewport.w := %d\n", vp->plot_area_viewport.w);
    printf("vp->coordinate_area_viewport.w := %d\n", vp->coordinate_area_viewport.w);
  }
  for(ti=0; ti<vp->x.tick_nr; ti++) {
    pos=vp->x.origo+vp->x.scale*vp->x.tick[ti];
    pos += vp->plot_area_viewport.w-vp->coordinate_area_viewport.w;
    if(verbose>1) printf("ti=%d: x tick pos=%g\n", ti, pos);
    sprintf(tmp, "      <%s", ilc);
    sprintf(line, "line x1=\"%g\" x2=\"%g\" y1=\"%g\" y2=\"%g\" />\n",
            pos, pos, ypos, ypos+yheight/8);
    strcat(tmp, line);
    if(svg_write(fp, tmp, errmsg, verbose-5)!=0) {return(6);}
  }
  // try to write one extra tick, if there is space
  pos=vp->x.origo+vp->x.scale*(2.0*vp->x.tick[vp->x.tick_nr-1] - vp->x.tick[vp->x.tick_nr-2]);
  pos += vp->plot_area_viewport.w-vp->coordinate_area_viewport.w;
  if(pos<vp->plot_area_viewport.w) {
    if(verbose>1) printf("extra ti=%d: x tick pos=%g\n", ti, pos);
    sprintf(tmp, "      <%s", ilc);
    sprintf(line, "line x1=\"%g\" x2=\"%g\" y1=\"%g\" y2=\"%g\" />\n",
            pos, pos, ypos, ypos+yheight/8);
    strcat(tmp, line);
    if(svg_write(fp, tmp, errmsg, verbose-5)!=0) {return(6);}
  }
  sprintf(tmp, "    </%sg>\n", ilc);
  if(svg_write(fp, tmp, errmsg, verbose-5)!=0) {return(7);}

  /* Write ticks labels */
  sprintf(tmp, "    <%sg", ilc);
  sprintf(line, " font-family=\"Sans-serif\""); strcat(tmp, line);
  sprintf(line, " text-anchor=\"middle\""); strcat(tmp, line);
  sprintf(line, " font-size=\"%d\"", vp->coordinate_area_viewport.chr_size); strcat(tmp, line);
  sprintf(line, " fill=\"black\""); strcat(tmp, line);
  strcat(tmp, ">\n"); 
  if(svg_write(fp, tmp, errmsg, verbose-5)!=0) {return(8);}
  for(ti=0; ti<vp->x.tick_nr; ti++) {
    pos=vp->x.origo+vp->x.scale*vp->x.tick[ti];
    pos += vp->plot_area_viewport.w-vp->coordinate_area_viewport.w;
    if(verbose>1) printf("ti=%d: x tick pos=%g\n", ti, pos);
    sprintf(tmp, "      <%s", ilc);
    sprintf(line, "text x=\"%g\" y=\"%g\"", pos, ypos+0.92*(double)yheight); strcat(tmp, line);
    sprintf(line, ">"); strcat(tmp, line);
    strcat(tmp, vp->x.tick_label[ti]);
    strcat(tmp, "</"); strcat(tmp, ilc); strcat(tmp, "text>\n");
    if(svg_write(fp, tmp, errmsg, verbose-5)!=0) {return(9);}
  }
  sprintf(tmp, "    </%sg>\n", ilc);
  if(svg_write(fp, tmp, errmsg, verbose-5)!=0) {return(10);}

  return(0);
}
/*****************************************************************************/

/*****************************************************************************/
/** Create SVG plot y axis ticks.
   @return Returns 0 if successful, <>0 in case of error.
   @sa svg_define_viewports, svg_write_axes, svg_write_xticks, svg_calculate_axes, axis_tick_positions
*/
int svg_write_yticks(
  /** SVG graphics file pointer. */
  FILE *fp,
  /** Struct containing the viewport sizes. */
  struct svg_viewports *vp,
  /** Char pointer to string (at least of length 128) where possible
      error description is copied; set to NULL if not necessary. */
  char *errmsg,
  /** Verbose level; if zero, then nothing is printed to stderr or stdout. */
  int verbose
) {
  char tmp[1024], line[128], ilc[9];
  int ti;
  double pos, xwidth;

  if(verbose>0) {printf("svg_write_yticks(fp, vp, errmsg, %d)\n", verbose); fflush(stdout);}

  /* Check the input */
  if(vp->y.tick_nr<1 || vp->plot_area_viewport.w==vp->coordinate_area_viewport.w) return(0);
  if(fp==NULL) {if(errmsg!=NULL) sprintf(errmsg, "file was closed too early"); return(1);}

  if(SVG_INLINE) strcpy(ilc, "svg:"); else strcpy(ilc, "");

  xwidth=vp->plot_area_viewport.w-vp->coordinate_area_viewport.w;

  strcpy(tmp, "\n    <!-- Y axis ticks inside plot area -->\n");
#if(0)
  sprintf(tmp,
    "      <%srect width=\"%gpx\" height=\"%dpx\" stroke=\"none\" fill=\"aqua\" fill-opacity=\"0.3\" />\n",
    ilc, xwidth, vp->plot_area_viewport.h);
#endif
  if(svg_write(fp, tmp, errmsg, verbose-5)!=0) {return(99);}

  if(verbose>0) {
    printf("vp->y.tick_nr=%d\n", vp->y.tick_nr);
    printf("vp->y.tickscale=%g\n", vp->y.tickscale);
    printf("vp->y.tick_decimals=%d\n", vp->y.tick_decimals);
  }

  /* Write ticks */
  sprintf(tmp, "    <%sg", ilc);
  sprintf(line, " stroke=\"%s\"", "black"); strcat(tmp, line);
  sprintf(line, " stroke-width=\"%g\"", 20.); strcat(tmp, line);
  sprintf(line, " fill=\"%s\"", "none"); strcat(tmp, line);
  strcat(tmp, ">\n"); 
  if(svg_write(fp, tmp, errmsg, verbose-5)!=0) {return(5);}
  for(ti=0; ti<vp->y.tick_nr; ti++) {
    pos=vp->coordinate_area_viewport.h-(vp->y.origo+vp->y.scale*vp->y.tick[ti]);
    if(verbose>1) printf("ti=%d: y tick pos=%g\n", ti, pos);
    sprintf(tmp, "      <%s", ilc);
    sprintf(line, "line x1=\"%g\" x2=\"%g\" y1=\"%g\" y2=\"%g\" />\n",
      xwidth,
      xwidth-(double)vp->coordinate_area_viewport.chr_size/8.0,
      pos, pos); strcat(tmp, line);
    if(svg_write(fp, tmp, errmsg, verbose-5)!=0) {return(6);}
  }
  sprintf(tmp, "    </%sg>\n", ilc);
  if(svg_write(fp, tmp, errmsg, verbose-5)!=0) {return(7);}

  /* Write ticks labels */
  sprintf(tmp, "    <%sg", ilc);
  sprintf(line, " font-family=\"Sans-serif\""); strcat(tmp, line);
  sprintf(line, " text-anchor=\"end\""); strcat(tmp, line);
  sprintf(line, " font-size=\"%d\"", vp->coordinate_area_viewport.chr_size); strcat(tmp, line);
  sprintf(line, " fill=\"black\""); strcat(tmp, line);
  strcat(tmp, ">\n"); 
  if(svg_write(fp, tmp, errmsg, verbose-5)!=0) {return(8);}
  for(ti=0; ti<vp->y.tick_nr; ti++) {
    pos=vp->coordinate_area_viewport.h-(vp->y.origo+vp->y.scale*vp->y.tick[ti]);
    if(verbose>1) printf("ti=%d: y tick pos=%g\n", ti, pos);
    sprintf(tmp, "      <%s", ilc);
    sprintf(line, "text x=\"%g\" y=\"%g\"",
      0.92*xwidth, pos+0.4*(double)vp->coordinate_area_viewport.chr_size);
    strcat(tmp, line);
    sprintf(line, ">"); strcat(tmp, line);
    strcat(tmp, vp->y.tick_label[ti]);
    strcat(tmp, "</"); strcat(tmp, ilc); strcat(tmp, "text>\n");
    if(svg_write(fp, tmp, errmsg, verbose-5)!=0) {return(9);}
  }
  sprintf(tmp, "    </%sg>\n", ilc);
  if(svg_write(fp, tmp, errmsg, verbose-5)!=0) {return(10);}

  return(0);
}
/*****************************************************************************/

/*****************************************************************************/
/** Draw sample curve in an SVG file.
   @return Returns 0 if successful, <>0 in case of error.
   @sa svg_initiate, svg_write, svg_define_viewports
*/
int svg_write_tac(
  /** SVG graphics file pointer */
  FILE *fp,
  /** Struct containing the viewport sizes */
  struct svg_viewports *vp,
  /** Plot type: 1=line, 2=symbols, 0=both line and symbols */
  const int plot_type,
  /** Unique ID for the curve */
  const char *tac_id,
  /** Title of the curve, which may be shown in the graph */
  const char *tac_title,
  /** Pointer to the polyline data x array (original quantities) */
  double *x,
  /** Pointer to the polyline data y array (original quantities) */
  double *y,
  /** Nr of data points in the array (half of array length) */
  const int data_nr,
  /** SVG color name as a string, e.g. aqua,black,blue,fuchsia,gray,
      green,lime,maroon,navy,olive,purple,red,silver,teal,yellow.
      Note that this string is not tested */
  const char *color,
  /** Symbol type: RECTANGLE,CIRCLE,UPTRIANGLE,DOWNTRIANGLE,DIAMOND,
                   LEFTTRIANGLE, RIGHTTRIANGLE */
  const svgSymbolType symbol_type,
  /** Symbol filling: SYMBOLOPEN, SYMBOLFILLED */
  const svgSymbolFill symbol_fill,
  /** Char pointer to string (at least of length 128) where possible
      error description is copied; set to NULL if not necessary */
  char *errmsg,
  /** Verbose level; if zero, then nothing is printed to stderr or stdout */
  int verbose
) {
  char tmp[1024], line[128], ilc[9];
  int i, j;
  double px, py, size=100.0, trgsize=140.0, circsize=120.0;

  if(verbose>0) {
    printf("svg_write_tac(fp, vp, %d, %s, %s, x, y, %d, %s, %d, %d, errmsg, %d)\n",
      plot_type, tac_id, tac_title, data_nr, color, (int)symbol_type, (int)symbol_fill, verbose);
    fflush(stdout);
  }

  /* Check the input */
  if(data_nr<1) return(0);
  if(fp==NULL) {if(errmsg!=NULL) sprintf(errmsg, "file was closed too early"); return(1);}
  if(color==NULL || strlen(color)<2) {
    if(errmsg!=NULL) sprintf(errmsg, "invalid color");
    return(1);
  }

  if(SVG_INLINE) strcpy(ilc, "svg:"); else strcpy(ilc, "");

  /* Initiate the curve object group */
  sprintf(tmp, "\n    <!-- %s : %s -->\n    <%sg", tac_id, tac_title, ilc);
  sprintf(line, " stroke=\"%s\"", color); strcat(tmp, line);
  sprintf(line, " stroke-width=\"%g\"", 0.25*size); strcat(tmp, line);
  sprintf(line, " fill=\"%s\"", color); strcat(tmp, line);
  if(symbol_fill==SYMBOLOPEN) sprintf(line, " fill-opacity=\"0.02\"");
  else sprintf(line, " fill-opacity=\"0.67\"");
  strcat(tmp, line); strcat(tmp, ">\n");
  if(svg_write(fp, tmp, errmsg, verbose-5)!=0) {return(11);}
  /* Set the curve title */
  sprintf(tmp, "      <%s", ilc);
  sprintf(line, "title>"); strcat(tmp, line);
  strcat(tmp, tac_title);
  strcat(tmp, "</"); strcat(tmp, ilc); strcat(tmp, "title>\n");
  if(svg_write(fp, tmp, errmsg, verbose-5)!=0) {return(12);}

  /* Plot the line, if required */
  strcpy(tmp, ""); strcpy(line, "");
  if(plot_type==0 || plot_type==1) {
    int lineon=0, path_started=0, prev_exists=0, cross_nr;
    double prev_px=-1, prev_py=-1, nx1, ny1, nx2, ny2;
    /* Start a new print line */
    strcpy(line, "\n       ");
    /* Write line coordinates */
    for(i=j=0; i<data_nr; i++) {
      if(isnan(x[i]) || isnan(y[i])) {lineon=0; continue;}
      /* Print recent line coordinates in file */
      if(j>=5) { /* line end */
        if(path_started==0) {
          sprintf(tmp, "      <%spath fill=\"none\" d=\"", ilc);
          if(svg_write(fp, tmp, errmsg, verbose-5)!=0) {return(14);}
          path_started=1;
        }
        if(svg_write(fp, line, errmsg, verbose-5)!=0) {return(15);}
        /* Start a new line */
        j=0; strcpy(line, "\n       ");
      }
      /* Compute the point coordinates in viewport */
      if(verbose>3) printf("x[%d]=%g y[%d]=%g\n", i, x[i], i, y[i]);
      px=vp->x.origo+x[i]*vp->x.scale;
      py=vp->coordinate_area_viewport.h-(vp->y.origo+vp->y.scale*y[i]);
      /* Make sure that prev point exists */
      if(prev_exists==0) {prev_px=px; prev_py=py; prev_exists=1;}
      /* Check if line would cross viewport border(s) */
      cross_nr=check_intersection_with_viewport(prev_px, prev_py, px, py,
                &vp->coordinate_area_viewport, &nx1, &ny1, &nx2, &ny2, verbose);
      if(verbose>2 &&cross_nr>0)
        printf("new line coordinates (%g,%g) -> (%g,%g)\n", nx1, ny1, nx2, ny2);
      if(cross_nr==2) {
        // Move to nx1,ny1 and draw line from nx1,ny1 to nx2,ny2
        if(j>0) strcat(line, " ");
        sprintf(tmp, "M%.0f %.0f L%.0f %.0f", nx1, ny1, nx2, ny2);
        if(verbose>4) printf("  write %s\n", tmp);
        strcat(line, tmp); j+=2;
        lineon=0; 
        // Proceed to next sample 
        prev_px=px; prev_py=py;
        continue;
      } else if(cross_nr==1) {
        if(nx1!=prev_px || ny1!=prev_py) lineon=0;
        // Draw line from nx1,ny1 to nx2,ny2
        if(j>0) strcat(line, " ");
        if(lineon==0) {sprintf(tmp, "M%.0f %.0f L%.0f %.0f", nx1, ny1, nx2, ny2); j+=2;}
        else if(lineon==1) {sprintf(tmp, "L%.0f %.0f", nx2, ny2); j++;}
        else {sprintf(tmp, "%.0f %.0f", nx2, ny2); j++;}
        if(verbose>4) printf("  write %s\n", tmp);
        strcat(line, tmp);
        if(nx2!=px || ny2!=py) lineon=0; else lineon++;
        // Proceed to next sample 
        prev_px=px; prev_py=py;
        continue;
      } else {
        nx1=prev_px; ny1=prev_py; nx2=px; ny2=py;
      }
      /* Draw line if coordinates are within viewport */
      if(nx1>=0 && nx1<=vp->coordinate_area_viewport.w+1 &&
         nx2>=0 && nx2<=vp->coordinate_area_viewport.w+1 &&
         ny1>=0 && ny1<=vp->coordinate_area_viewport.h+1 &&
         ny2>=0 && ny2<=vp->coordinate_area_viewport.h+1)
      {
        if(j>0) strcat(line, " ");
        if(lineon==0) {
          if(nx1!=nx2 && ny1!=ny2) sprintf(tmp, "M%.0f %.0f L%.0f %.0f", nx1, ny1, nx2, ny2); 
          else sprintf(tmp, "M%.0f %.0f", nx2, ny2); 
          j+=2;
        }
        else if(lineon==1) {sprintf(tmp, "L%.0f %.0f", nx2, ny2); j++;}
        else {sprintf(tmp, "%.0f %.0f", nx2, ny2); j++;}
        if(verbose>4) printf("  write %s\n", tmp);
        strcat(line, tmp);
        lineon++;
      } else {
        lineon=0;
      }
      prev_px=px; prev_py=py;
    }
    /* Write into file the remaining (if any) points */
    if(j>0) {
      if(path_started==0) {
        sprintf(tmp, "      <%spath fill=\"none\" d=\"", ilc);
        if(svg_write(fp, tmp, errmsg, verbose-5)!=0) {return(14);}
      }
      if(svg_write(fp, line, errmsg, verbose-5)!=0) {return(16);}
    }
    /* Close line */
    strcpy(tmp, "\" />\n");
    if(svg_write(fp, tmp, errmsg, verbose-5)!=0) {return(17);}
  }

  /* Plot the symbols, if required */
  strcpy(tmp, "");
  if(plot_type==0 || plot_type==2) {
    int prev_px=-1, prev_py=-1;
    for(i=0; i<data_nr; i++) {
      if(isnan(x[i]) || isnan(y[i])) continue;
      px=vp->x.origo+x[i]*vp->x.scale;
      py=vp->coordinate_area_viewport.h-(vp->y.origo+vp->y.scale*y[i]);
      /* Do not plot points outside viewport */
      if(px<0 || py<0) continue;
      if(px>vp->coordinate_area_viewport.w+1 || py>vp->coordinate_area_viewport.h+1) continue;
      /* Do not plot 2nd time the same point */
      if(px==prev_px && py==prev_py) continue;
      prev_px=px; prev_py=py;
      /* Draw the symbol */
      sprintf(tmp, "      <%suse ", ilc);
      switch(symbol_type) {
        case RECTANGLE:
          sprintf(line, "xlink:href=\"#sym-rect\" x=\"%g\" y=\"%g\" width=\"%g\" height=\"%g\"",
            px-0.5*size, py-0.5*size, size, size);
          break;
        case UPTRIANGLE:
          sprintf(line, "xlink:href=\"#sym-uptr\" x=\"%g\" y=\"%g\" width=\"%g\" height=\"%g\"",
            px-0.5*trgsize, py-0.5*trgsize, trgsize, trgsize);
          break;
        case DOWNTRIANGLE:
          sprintf(line, "xlink:href=\"#sym-dotr\" x=\"%g\" y=\"%g\" width=\"%g\" height=\"%g\"",
            px-0.5*trgsize, py-0.5*trgsize, trgsize, trgsize);
          break;
        case DIAMOND:
          sprintf(line, "xlink:href=\"#sym-diam\" x=\"%g\" y=\"%g\" width=\"%g\" height=\"%g\"",
            px-0.5*trgsize, py-0.5*trgsize, trgsize, trgsize);
          break;
        case LEFTTRIANGLE:
          sprintf(line, "xlink:href=\"#sym-letr\" x=\"%g\" y=\"%g\" width=\"%g\" height=\"%g\"",
            px-0.5*trgsize, py-0.5*trgsize, trgsize, trgsize);
          break;
        case RIGHTTRIANGLE:
          sprintf(line, "xlink:href=\"#sym-ritr\" x=\"%g\" y=\"%g\" width=\"%g\" height=\"%g\"",
            px-0.5*trgsize, py-0.5*trgsize, trgsize, trgsize);
          break;
        case CIRCLE:
        default:
          sprintf(line, "xlink:href=\"#sym-circ\" x=\"%g\" y=\"%g\" width=\"%g\" height=\"%g\"",
            px-0.5*circsize, py-0.5*circsize, circsize, circsize);
          break;
      }
      strcat(tmp, line);
      strcat(tmp, " />\n");
      if(svg_write(fp, tmp, errmsg, verbose-5)!=0) {return(18);}
    }
  }

  /* Close the curve object group */
  sprintf(tmp, "    </%sg>\n", ilc);
  if(svg_write(fp, tmp, errmsg, verbose-5)!=0) {return(19);}

  return(0);
}
/*****************************************************************************/

/*****************************************************************************/
