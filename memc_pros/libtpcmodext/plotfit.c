/// @file plotfit.c
/// @brief Plot measured and fitted TACs in SVG format.
/// @author Vesa Oikonen
///
/*****************************************************************************/

/*****************************************************************************/
#include "libtpcmodext.h"
/*****************************************************************************/

/*****************************************************************************/
/** Writes specified range of plots of original and fitted TACs in SVG 1.1 format.
   @return Returns 0 if successful, otherwise nonzero.
   @sa plot_fit_svg, plot_svg
 */
int plot_fitrange_svg(
  /** Measured data points */
  DFT *dft1,
  /** Fitted data points. Times can be different but unit must be the same. */
  DFT *dft2,
  /** String for plot main title, or NULL */
  char *main_title,
  /** Start time; NaN if determined from data */
  double x1,
  /** End time; NaN if determined from data */
  double x2,
  /** Minimum y value; NaN if determined from data */
  double y1,
  /** Maximum y value; NaN if determined from data */
  double y2,
  /** SVG filename; existing file is backed up */ 
  char *fname,
  /** Verbose level; set to zero to not to print any comments */
  int verbose
) {
  int ret, n, ri, si, ei;
  char x_title[64], y_title[64], tac_id[32], tac_title[64];
  double minx, maxx, miny, maxy, tx1, tx2, ty1, ty2, f;
  struct svg_viewports viewports; svg_init_viewports(&viewports);
  int max_color_nr, color_nr;
  int max_symbol_nr, symbol_nr;
  SVG_LEGENDS legends; svg_init_legends(&legends);
  FILE *fp_svg=NULL;

  if(verbose>0) {
    printf("plot_fitrange_svg(dft1, dft2, mt, x1, x2, y1, y2, fn, %d)\n", verbose);
  }

  /* Check data */
  if(dft1==NULL || dft1->voiNr<1) return(1);
  if(dft2==NULL) return(1);
  if(dft2->voiNr!=dft1->voiNr) return(1);

  /* Check if file exists; backup, if necessary */
  ret=backupExistingFile(fname, NULL, NULL); if(ret) return(2);

  int is_label=0; if(dft1->voiNr>1) is_label=1;

  /* Determine the plot min and max x values */
  ret=dftMinMax(dft1, &tx1, &tx2, NULL, NULL); if(ret) return(3);
  minx=tx1; maxx=tx2;
  ret=dftMinMax(dft2, &tx1, &tx2, NULL, NULL); if(ret) return(3);
  if(minx>tx1) minx=tx1; 
  if(maxx<tx2) maxx=tx2;
  if(minx>0.0) {
    f=maxx-minx; minx-=0.05*f; 
    if(minx<0.0) minx=0.0;
  }
  if(!isnan(x1)) minx=x1; 
  if(!isnan(x2)) maxx=x2;

  /* Determine the plot min and max y values */
  ret=dftMaxY(dft1, minx, maxx, &ty1, &ty2); if(ret) return(3);
  miny=ty1; maxy=ty2;
  ret=dftMaxY(dft2, minx, maxx, &ty1, &ty2); if(ret) return(3);
  if(miny>ty1) miny=ty1; 
  if(maxy<ty2) maxy=ty2;
  if(miny>0.0) {
    f=maxy-miny; miny-=0.05*f; 
    if(miny<0.0) miny=0.0;
  }
  if(!isnan(y1)) miny=y1; 
  if(!isnan(y2)) maxy=y2;

  if(verbose>1) printf("minx:=%g\nmaxx:=%g\nminy:=%g\nmaxy:=%g\n", minx, maxx, miny, maxy);

  /* Calculate the axis ticks */
  viewports.label_area_viewport.is=is_label; // needed for x axis ticks
  viewports.x.min=minx; viewports.x.max=maxx;
  viewports.y.min=miny; viewports.y.max=maxy;
  if(isnan(x1) || isnan(x2)) viewports.x.fixed_min=0; 
  if(isnan(y1) || isnan(y2)) viewports.y.fixed_min=0; 
  else viewports.y.fixed_min=1; 
  ret=svg_calculate_axes(&viewports, verbose-3); if(ret) return(4);

  /* Set x and y axis titles based on activity and time units */
  strcpy(x_title, "");
  if(dft1->timeunit==DFTTIME_SEC || dft1->timeunit==DFTTIME_MIN) 
    sprintf(x_title, "Time (%s)", dftTimeunit(dft1->timeunit));
  else if(dft1->timeunit!=DFTTIME_UNKNOWN)
    strcpy(x_title, dftTimeunit(dft1->timeunit));
  strcpy(y_title, "");
  if(dftUnitId(dft1->unit)!=DFTUNIT_UNKNOWN)
    strcpy(y_title, dftUnit(dftUnitId(dft1->unit)));

  /* Set the plot window and window area sizes */
  ret=svg_define_viewports(0, 0, strlen(main_title), strlen(y_title),
    strlen(x_title), is_label, &viewports, verbose-3);
  if(ret) return(5);

  /* Initiate graphics file */
  fp_svg=svg_initiate(fname, 0, 0, &viewports, NULL, verbose-3);
  if(fp_svg==NULL) return(6);

  /* Put the graph titles into their own viewports */
  ret=svg_create_main_title(fp_svg, main_title, "", &viewports, NULL,verbose-3);
  if(ret) return(7);
  ret=svg_create_yaxis_title(fp_svg, y_title, &viewports, NULL, verbose-3);
  if(ret) return(8);
  ret=svg_create_xaxis_title(fp_svg, x_title, &viewports, NULL, verbose-3);
  if(ret) return(9);

  /* Put the plot into its own viewport */
  ret=svg_start_plot_viewport(fp_svg, &viewports, NULL, verbose-3);
  if(ret) return(10);

  /*  Start coordinate area viewport */
  ret=svg_start_coordinate_viewport(fp_svg, &viewports, NULL, verbose-3);
  if(ret) return(11);

  /* Write plot axes */
  ret=svg_write_axes(fp_svg, &viewports, NULL, verbose-3);
  if(ret) return(12);

  /*
   *  Draw the plots
   */
  max_color_nr=0; while(svgColorName(max_color_nr)!=NULL) max_color_nr++;
  if(max_color_nr==0) max_color_nr=1; // remainder works only if 2nd operator>0
  if(verbose>3) printf("max_color_nr := %d\n", max_color_nr);
  max_symbol_nr=0; while(svgSymbolName(max_symbol_nr)!=NULL) max_symbol_nr++;
  if(max_symbol_nr==0) max_symbol_nr=1; // remainder works only if 2nd operator>0
  if(verbose>3) printf("max_symbol_nr := %d\n", max_symbol_nr);
  if(dft1->voiNr==1) color_nr=0; else color_nr=1;
  symbol_nr=0;
  for(ri=0, n=0; ri<dft1->voiNr; ri++) {
    sprintf(tac_id, "plot_%d", n);
/*
    if(strlen(dft1->studynr)>0 && strcmp(dft1->studynr, ".")!=0)
      sprintf(tac_title, "%s: %s", dft1->studynr, dft1->voi[ri].name);
    else strcpy(tac_title, dft1->voi[ri].name);
*/
    rnameRmDots(dft1->voi[ri].name, tac_title);
    /* Draw the fitted line */
    for(si=0; si<dft2->frameNr; si++) if(!isnan(dft2->voi[ri].y[si])) break;
    for(ei=dft2->frameNr-1; ei>=si; ei--) if(!isnan(dft2->voi[ri].y[ei])) break;
    if((ei-si)>0) {
      ret=svg_write_tac(fp_svg, &viewports, 1, tac_id, tac_title,
            dft2->x+si, dft2->voi[ri].y+si, 1+ei-si,
            svgColorName(color_nr%max_color_nr), symbol_nr%max_symbol_nr, SYMBOLFILLED,
            NULL, verbose-3);
      if(ret) {svg_legend_empty(&legends); return(21);}
    }
    /* Draw the measured points */
    ret=svg_write_tac(fp_svg, &viewports, 2, tac_id, tac_title,
            dft1->x, dft1->voi[ri].y, dft1->frameNr,
            svgColorName(color_nr%max_color_nr), symbol_nr%max_symbol_nr, SYMBOLFILLED,
            NULL, verbose-3);
    if(ret) {svg_legend_empty(&legends); return(22);}
    /* Set legend too, if requested */
    if(is_label!=0) {
      svg_legend_add(&legends, 0, symbol_nr%max_symbol_nr, SYMBOLFILLED,
                     color_nr%max_color_nr, tac_title);
    }
    /* Prepare for the next plot */
    color_nr++; n++;
    if(color_nr==max_color_nr) {symbol_nr++; color_nr=0;}
    if(symbol_nr==max_symbol_nr) symbol_nr=0;
  }

  /* Close the coordinate viewport */
  ret=svg_end_coordinate_viewport(fp_svg, NULL, verbose-3);
  if(ret) {svg_legend_empty(&legends); return(91);}

  /* Write the axis ticks */
  if(svg_write_xticks(fp_svg, &viewports, NULL, verbose-3)!=0) {
    svg_legend_empty(&legends); return(92);}
  if(svg_write_yticks(fp_svg, &viewports, NULL, verbose-3)!=0) {
    svg_legend_empty(&legends); return(93);}

  /* Close the plot viewport */
  ret=svg_end_plot_viewport(fp_svg, NULL, verbose-3);
  if(ret) {svg_legend_empty(&legends); return(94);}

  /* Make the plot legends into their own viewport */
  if(viewports.label_area_viewport.is!=0) {
    if(verbose>2) printf("creating plot legends\n");
    ret=svg_create_legends(fp_svg, &viewports, &legends, NULL, verbose-3);
    if(ret) {svg_legend_empty(&legends); return(95);}
  }
  svg_legend_empty(&legends);

  /* Close the SVG file */
  ret=svg_close(fp_svg, NULL, verbose-3); if(ret) return(101);

  return(0);
}
/*****************************************************************************/

/*****************************************************************************/
/** Writes plots of original and fitted TACs in SVG 1.1 format.
    Data must not contain NaNs.
   @return Returns 0 if successful, otherwise nonzero.
   @sa plot_svg, plot_fitrange_svg
 */
int plot_fit_svg(
  /** Measured data points */
  DFT *dft1,
  /** Fitted data points. Times can be different but unit must be the same. */
  DFT *dft2,
  /** String for plot main title, or "" */
  char *main_title,
  /** SVG filename; existing file is backed up */ 
  char *fname,
  /** Verbose level; if zero, then nothing is printed to stderr or stdout */
  int verbose
) {
  int ret, n, ri, si, ei;
  char x_title[64], y_title[64], tac_id[32], tac_title[64];
  double minx, maxx, miny, maxy, tx1, tx2, ty1, ty2, f;
  struct svg_viewports viewports; svg_init_viewports(&viewports);
  int max_color_nr, color_nr;
  int max_symbol_nr, symbol_nr;
  SVG_LEGENDS legends; svg_init_legends(&legends);
  FILE *fp_svg=NULL;


  if(verbose>0) {
    printf("plot_fit_svg(dft1, dft2, mt, fn, %d)\n", verbose);
  }

  /* Check data */
  if(dft1==NULL || dft1->voiNr<1) return(1);
  if(dft2==NULL) return(1);
  if(dft2->voiNr!=dft1->voiNr) return(1);

  int is_label=0; if(dft1->voiNr>1) is_label=1;

  /* Check if file exists; backup, if necessary */
  ret=backupExistingFile(fname, NULL, NULL); if(ret) return(2);

  /* Determine the plot min and max values */
  minx=maxx=miny=maxy=0;
  ret=dftMinMax(dft1, &tx1, &tx2, &ty1, &ty2); if(ret) return(3);
  minx=tx1; maxx=tx2; miny=ty1; maxy=ty2;
  ret=dftMinMax(dft2, &tx1, &tx2, &ty1, &ty2); if(ret) return(3);
  if(minx>tx1) minx=tx1;
  if(maxx<tx2) maxx=tx2;
  if(miny>ty1) miny=ty1;
  if(maxy<ty2) maxy=ty2;
  if(verbose>1) 
    printf("minx:=%g\nmaxx:=%g\nminy:=%g\nmaxy:=%g\n", minx, maxx, miny, maxy);
  if(miny>0.0) {f=maxy-miny; miny-=0.01*f;}

  /* Calculate the axis ticks */
  viewports.label_area_viewport.is=is_label; // needed for x axis ticks
  viewports.x.fixed_min=0; viewports.y.fixed_min=0;
  viewports.x.min=minx; viewports.x.max=maxx;
  viewports.y.min=miny; viewports.y.max=maxy;
  ret=svg_calculate_axes(&viewports, verbose-3); if(ret) return(4);

  /* Set x and y axis titles based on activity and time units */
  if(verbose>2) printf("set x title\n");
  strcpy(x_title, "");
  if(dft1->timeunit==DFTTIME_SEC || dft1->timeunit==DFTTIME_MIN) 
    sprintf(x_title, "Time (%s)", dftTimeunit(dft1->timeunit));
  else if(dft1->timeunit!=DFTTIME_UNKNOWN)
    strcpy(x_title, dftTimeunit(dft1->timeunit));
  if(verbose>2) printf("set y title\n");
  strcpy(y_title, "");
  if(dftUnitId(dft1->unit)!=DFTUNIT_UNKNOWN)
    strcpy(y_title, dftUnit(dftUnitId(dft1->unit)));

  /* Set the plot window and window area sizes */
  if(verbose>2) printf("set window sizes\n");
  ret=svg_define_viewports(0, 0, strlen(main_title), strlen(y_title),
    strlen(x_title), is_label, &viewports, verbose-3);
  if(ret) return(5);

  /* Initiate graphics file */
  fp_svg=svg_initiate(fname, 0, 0, &viewports, NULL, verbose-3);
  if(fp_svg==NULL) return(6);

  /* Put the graph titles into their own viewports */
  ret=svg_create_main_title(fp_svg, main_title, "", &viewports, NULL,verbose-3);
  if(ret) return(7);
  ret=svg_create_yaxis_title(fp_svg, y_title, &viewports, NULL, verbose-3);
  if(ret) return(8);
  ret=svg_create_xaxis_title(fp_svg, x_title, &viewports, NULL, verbose-3);
  if(ret) return(9);

  /* Put the plot into its own viewport */
  ret=svg_start_plot_viewport(fp_svg, &viewports, NULL, verbose-3);
  if(ret) return(10);

  /*  Start coordinate area viewport */
  ret=svg_start_coordinate_viewport(fp_svg, &viewports, NULL, verbose-3);
  if(ret) return(11);

  /* Write plot axes */
  ret=svg_write_axes(fp_svg, &viewports, NULL, verbose-3);
  if(ret) {
    if(verbose>0) printf("svg_write_axes() := %d\n", ret);
    return(12);
  }

  /*
   *  Draw the plots
   */
  max_color_nr=0; while(svgColorName(max_color_nr)!=NULL) max_color_nr++;
  if(max_color_nr<1) max_color_nr=1; // remainder works only if 2nd operator>0
  if(verbose>3) printf("max_color_nr := %d\n", max_color_nr);
  max_symbol_nr=0; while(svgSymbolName(max_symbol_nr)!=NULL) max_symbol_nr++;
  if(max_symbol_nr<1) max_symbol_nr=1; // remainder works only if 2nd operator>0
  if(verbose>3) printf("max_symbol_nr := %d\n", max_symbol_nr);
  if(dft1->voiNr==1) color_nr=0; else color_nr=1;
  symbol_nr=0;
  for(ri=0, n=0; ri<dft1->voiNr; ri++) {
    sprintf(tac_id, "plot_%d", n);
/*
    if(strlen(dft1->studynr)>0 && strcmp(dft1->studynr, ".")!=0)
      sprintf(tac_title, "%s: %s", dft1->studynr, dft1->voi[ri].name);
    else strcpy(tac_title, dft1->voi[ri].name);
*/
    rnameRmDots(dft1->voi[ri].name, tac_title);
    /* Draw the fitted line */
    for(si=0; si<dft2->frameNr; si++) if(!isnan(dft2->voi[ri].y[si])) break;
    for(ei=dft2->frameNr-1; ei>=si; ei--) if(!isnan(dft2->voi[ri].y[ei])) break;
    if((ei-si)>0) {
      ret=svg_write_tac(fp_svg, &viewports, 1, tac_id, tac_title,
            dft2->x+si, dft2->voi[ri].y+si, 1+ei-si,
            svgColorName(color_nr%max_color_nr), symbol_nr%max_symbol_nr, SYMBOLFILLED,
            NULL, verbose-3);
      if(ret) {svg_legend_empty(&legends); return(21);}
    }
    /* Draw the measured points */
    ret=svg_write_tac(fp_svg, &viewports, 2, tac_id, tac_title,
            dft1->x, dft1->voi[ri].y, dft1->frameNr,
            svgColorName(color_nr%max_color_nr), symbol_nr%max_symbol_nr, SYMBOLFILLED,
            NULL, verbose-3);
    if(ret) {svg_legend_empty(&legends); return(22);}
    /* Set legend too, if requested */
    if(is_label!=0) {
      svg_legend_add(&legends, 0, symbol_nr%max_symbol_nr, SYMBOLFILLED,
                     color_nr%max_color_nr, tac_title);
    }
    /* Prepare for the next plot */
    color_nr++; n++;
    if(color_nr==max_color_nr) {symbol_nr++; color_nr=0;}
    if(symbol_nr==max_symbol_nr) symbol_nr=0;
  }

  /* Close the coordinate viewport */
  ret=svg_end_coordinate_viewport(fp_svg, NULL, verbose-3);
  if(ret) {svg_legend_empty(&legends); return(91);}

  /* Write the axis ticks */
  if(svg_write_xticks(fp_svg, &viewports, NULL, verbose-3)!=0) {
    svg_legend_empty(&legends); return(92);}
  if(svg_write_yticks(fp_svg, &viewports, NULL, verbose-3)!=0) {
    svg_legend_empty(&legends); return(93);}

  /* Close the plot viewport */
  ret=svg_end_plot_viewport(fp_svg, NULL, verbose-3);
  if(ret) {svg_legend_empty(&legends); return(94);}

  /* Make the plot legends into their own viewport */
  if(viewports.label_area_viewport.is!=0) {
    if(verbose>2) printf("creating plot legends\n");
    ret=svg_create_legends(fp_svg, &viewports, &legends, NULL, verbose-3);
    if(ret) {svg_legend_empty(&legends); return(95);}
  }
  svg_legend_empty(&legends);

  /* Close the SVG file */
  ret=svg_close(fp_svg, NULL, verbose-3); if(ret) return(101);

  return(0);
}
/*****************************************************************************/

/*****************************************************************************/
