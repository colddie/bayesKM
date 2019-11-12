/// @file plotdata.c
/// @brief Write linear plot data and fitted lines in HTML tables or SVG plots.
/// @author Vesa Oikonen
///
/*****************************************************************************/

/*****************************************************************************/
#include "libtpcmodext.h"
/*****************************************************************************/

/*****************************************************************************/
/** Writes graphical analysis plots in SVG 1.1 format.
    Assumes that line slope and ic are in res->parameter[0] and [1].
   @return Returns 0 if successful, otherwise nonzero.
   @sa plot_fit_svg, plot_fitrange_svg
 */
int plot_svg(
  /** Plot points: X in y2, Y in y3 */
  DFT *dft,
  /** Results containing parameters of line */
  RES *res, 
  /** First sample (starting from 0) used in linear fit */
  int first,
  /** last sample (starting from 0) used in linear fit */
  int last,
  /** String for plot main title, or NULL */
  char *main_title,
  /** String for X axis title, or NULL */
  char *x_title,
  /** String for Y axis title, or NULL */
  char *y_title,
  /** SVG filename; existing file is renamed as *.bak */ 
  char *fname,
  /** Verbose level; if zero, then nothing is printed to stderr or stdout */
  int verbose
) {
  int n, ri, fi, ret;
  char tac_id[32], tac_title[64];
  double maxPlotX=0, maxy;
  double px[2], py[2];
  struct svg_viewports viewports; svg_init_viewports(&viewports);
  int max_color_nr, color_nr;
  int max_symbol_nr, symbol_nr;
  SVG_LEGENDS legends; svg_init_legends(&legends);
  FILE *fp_svg=NULL;

  if(verbose>0) {
    printf("plot_svg(dft, res, %d, %d, mt, xt, yt, fn, %d)\n", first, last, verbose);
  }

  /* Check data */
  if(dft==NULL || dft->voiNr<1) return(1);
  if(res==NULL || res->voiNr!=dft->voiNr) return(1);
  if(first>last) return(1);

  int is_label=0; if(dft->voiNr>1) is_label=1;

  /* Check if file exists; backup, if necessary */
  backupExistingFile(fname, NULL, NULL);

  /* Search the largest plot x-value to be used as line end point */
  for(ri=0; ri<dft->voiNr; ri++) for(fi=0; fi<dft->frameNr; fi++)
    if(dft->voi[ri].y2[fi]>maxPlotX) maxPlotX=dft->voi[ri].y2[fi];
  /* Get maxy */
  maxy=0.0;
  for(ri=0; ri<dft->voiNr; ri++) for(fi=0; fi<dft->frameNr; fi++)
    if(dft->voi[ri].y3[fi]>maxy) maxy=dft->voi[ri].y3[fi];

  /* Calculate the axis ticks */
  viewports.label_area_viewport.is=is_label; // needed for x axis ticks
  viewports.x.fixed_min=0; viewports.y.fixed_min=0;
  viewports.x.min=0.0; viewports.x.max=maxPlotX;
  viewports.y.min=0.0; viewports.y.max=maxy;
  ret=svg_calculate_axes(&viewports, verbose-3);
  if(ret) return(2);

  /* Set the plot window and window area sizes */
  ret=svg_define_viewports(0, 0, strlen(main_title), strlen(y_title),
                           strlen(x_title), is_label, &viewports, verbose-3);
  if(ret) return(3);

  /* Initiate graphics file */
  fp_svg=svg_initiate(fname, 0, 0, &viewports, NULL, verbose-3);
  if(fp_svg==NULL) return(4);

  /* Put the graph titles into their own viewports */
  ret=svg_create_main_title(fp_svg, main_title, "", &viewports, NULL,verbose-3);
  if(ret) return(5);
  ret=svg_create_yaxis_title(fp_svg, y_title, &viewports, NULL, verbose-3);
  if(ret) return(6);
  ret=svg_create_xaxis_title(fp_svg, x_title, &viewports, NULL, verbose-3);
  if(ret) return(7);

  /* Put the plot into its own viewport */
  ret=svg_start_plot_viewport(fp_svg, &viewports, NULL, verbose-3);
  if(ret) return(8);

  /*  Start coordinate area viewport */
  ret=svg_start_coordinate_viewport(fp_svg, &viewports, NULL, verbose-3);
  if(ret) return(9);

  /* Write plot axes */
  ret=svg_write_axes(fp_svg, &viewports, NULL, verbose-3);
  if(ret) return(10);

  /*
   *  Draw the plots
   */
  max_color_nr=0; while(svgColorName(max_color_nr)!=NULL) max_color_nr++;
  if(max_color_nr==0) max_color_nr=1; // remainder works only if 2nd operator>0
  if(verbose>3) printf("max_color_nr := %d\n", max_color_nr);
  max_symbol_nr=0; while(svgSymbolName(max_symbol_nr)!=NULL) max_symbol_nr++;
  if(max_symbol_nr==0) max_symbol_nr=1; // remainder works only if 2nd operator>0
  if(verbose>3) printf("max_symbol_nr := %d\n", max_symbol_nr);
  if(dft->voiNr==1) color_nr=0; else color_nr=1;
  symbol_nr=0;
  for(ri=0, n=0; ri<dft->voiNr; ri++) {
    sprintf(tac_id, "plot_%d", n);
    //printf("ri=%d color_nr=%d symbol_nr=%d\n", ri, color_nr, symbol_nr);
/*
    if(strlen(dft->studynr)>0 && strcmp(dft->studynr, ".")!=0)
      sprintf(tac_title, "%s: %s", dft->studynr, dft->voi[ri].name);
    else strcpy(tac_title, dft->voi[ri].name);
*/
    rnameRmDots(dft->voi[ri].name, tac_title);
    /*printf("tac_title := %s ; color := %s ; symbol := %d\n",
      tac_title, svgcolor[color_nr%max_color_nr], symbol_nr%(max_symbol_nr+1));*/
    /* Draw symbols */
    if(dft->frameNr<150) {
      /* plot symbols only if less than 150 samples */
      ret=svg_write_tac(fp_svg, &viewports, 2, tac_id, tac_title,
            dft->voi[ri].y2, dft->voi[ri].y3, dft->frameNr,
            svgColorName(color_nr%max_color_nr), symbol_nr%max_symbol_nr, SYMBOLFILLED,
            NULL, verbose-3);
    } else {
      /* plot samples as line if 150 samples or more */
      ret=svg_write_tac(fp_svg, &viewports, 1, tac_id, tac_title,
            dft->voi[ri].y2, dft->voi[ri].y3, dft->frameNr,
            svgColorName(color_nr%max_color_nr), symbol_nr%max_symbol_nr, SYMBOLFILLED,
            NULL, verbose-3);
    }
    if(ret) {svg_legend_empty(&legends); return(21);}
    /* Draw the line */
    px[0]=0.0; py[0]=res->voi[ri].parameter[1];
    px[1]=maxPlotX;
    py[1]=maxPlotX*res->voi[ri].parameter[0]+res->voi[ri].parameter[1];
    sprintf(tac_id, "line_%d", n);
    ret=svg_write_tac(fp_svg, &viewports, 1, tac_id, tac_title,
            px, py, 2,
            svgColorName(color_nr%max_color_nr), symbol_nr%max_symbol_nr, SYMBOLFILLED,
            NULL, verbose-3);
    if(ret) {svg_legend_empty(&legends); return(22);}
    /* Set legend too, if requested */
    if(is_label!=0) {
      svg_legend_add(&legends, 0, symbol_nr%max_symbol_nr, SYMBOLFILLED, color_nr%max_color_nr, tac_title);
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
/******************************************************************************/

/******************************************************************************/
/** Write plot and line fit data in XHTML 1.1. Strict table format.
    Assumes that line slope and ic are in res->parameter[0] and [1].
   @return Returns 0 if successful, otherwise nonzero.
   @sa plotdata_as_dft
 */
int plotdata(
  /** Plot points: X in y2, Y in y3 */
  DFT *dft,
  /** Results containing parameters of line */
  RES *res,
  /** First sample (starting from 0) used in linear fit */
  int first,
  /** last sample (starting from 0) used in linear fit */
  int last,
  /** String for plot main title, or NULL */
  char *mtitle,
  /** String for X axis title, or NULL */
  char *xtitle,
  /** String for Y axis title, or NULL */
  char *ytitle,
  /** Filename for plot data; existing file is renamed as *%.
   *  If extension is .dft, plot data (excluding lines) is written in DFT
   *  format with x values as separate columns before corresponding y values. */ 
  char *fname
) {
  int n, ri, row, fi;
  char tmp[FILENAME_MAX], *cptr=NULL;
  FILE *fp;
  double maxPlotX=0, maxRegX=0, maxFitX=0, f;


  /* Check data */
  if(dft==NULL || dft->voiNr<1) return(1);
  if(res==NULL || res->voiNr!=dft->voiNr) return(1);
  if(first>last) return(1);

  /* Get filename extension to determine output type */
  cptr=strrchr(fname, '.');
  if(cptr!=NULL && strcasecmp(cptr, ".DFT")==0) {
    /* Write in DFT if required */
    return(plotdata_as_dft(dft, fname));
  }

  /* Check if file exists; backup, if necessary */
  backupExistingFile(fname, NULL, NULL);

  /* Search the largest plot x-value to be used as line end point */
  for(ri=0; ri<dft->voiNr; ri++) for(fi=0; fi<dft->frameNr; fi++)
    if(dft->voi[ri].y2[fi]>maxPlotX) maxPlotX=dft->voi[ri].y2[fi];

  /* Open output file */
  if((fp = fopen(fname, "w")) == NULL) return(3);

  /* Write XHTML doctype and head */
//  n=fprintf(fp, "<!DOCTYPE html PUBLIC \"-//W3C//DTD XHTML 1.1//EN\" \"ht
//tp://www.w3.org/TR/xhtml11/DTD/xhtml11.dtd\">\n");
  n=fprintf(fp, "<!DOCTYPE html PUBLIC \"-//W3C//DTD XHTML 1.1//EN\"");
  if(n<10) return(4);
  n=fprintf(fp, " \"http://www.w3.org/TR/xhtml11/DTD/xhtml11.dtd\">\n");
  if(n<10) return(4);
  n=fprintf(fp, "<html xmlns=\"http://www.w3.org/1999/xhtml\" xml:lang=\"en\">\n\n");
  if(n<20) return(4);
  /* Write XHTML header */
  n=fprintf(fp, "<head>\n"); if(n<6) return(4);
  fprintf(fp, "  <title>Graphical analysis plot</title>\n");
//  fprintf(fp, "  <meta http-equiv=\"content-type\" cont
//ent=\"text/html; charset=iso-8859-1\" />\n");
  fprintf(fp, "  <meta http-equiv=\"content-type\" content=\"text/html;");
  fprintf(fp, " charset=iso-8859-1\" />\n");

  fprintf(fp, "  <meta http-equiv=\"content-language\" content=\"en-gb\" />\n");
  fprintf(fp, "  <meta name=\"ProgId\" content=\"Excel.Sheet\" />\n");
/*
  fprintf(fp, "  <link rel=\"icon\" href=\"http://www.turkupetcentre.
net/favicon.ico\" type=\"image/x-icon\" />\n");
  fprintf(fp, "  <link rel=\"shortcut icon\" href=\"http://www.turkupetcentre.
net/favicon.ico\" type=\"image/x-icon\" />\n");
*/
  fprintf(fp, "  <link rel=\"icon\" href=\"http://www.turkupetcentre.net/");
  fprintf(fp, "favicon.ico\" type=\"image/x-icon\" />\n");
  fprintf(fp, "  <link rel=\"shortcut icon\" href=\"http://www.turkupet");
  fprintf(fp, "centre.net/favicon.ico\" type=\"image/x-icon\" />\n");

  fprintf(fp, "  <style type=\"text/css\">\n");
  fprintf(fp, "    thead {background-color:#999999; color:black;}\n");
/*  fprintf(fp, "    table {text-align:left; width:100%%; 
border-collapse:collapse; empty-cells:show;}\n");*/
  fprintf(fp, "    table {text-align:left; width:100%%;");
  fprintf(fp, " border-collapse:collapse; empty-cells:show;}\n");

  fprintf(fp, "    td {border:1px solid black;}\n");
  fprintf(fp, "    <!--table\n");
  fprintf(fp, "    	{mso-displayed-decimal-separator:\"\\.\";\n");
  fprintf(fp, "    	mso-displayed-thousand-separator:\" \";}\n");
  fprintf(fp, "    -->\n");
  fprintf(fp, "  </style>\n");
  fprintf(fp, "</head>\n");

  /* Start writing the body of the HTML file */
  fprintf(fp, "\n<body>\n");

  /* Start the div for tables */
  fprintf(fp, "\n<div id=\"tables\">\n");

  /* Write information on the graphical analysis */
  fprintf(fp, "<table>\n");
  fprintf(fp, "<tbody>\n");
  fprintf(fp, "<tr><th>Main title</th><th>%s</th></tr>\n", mtitle);
  fprintf(fp, "<tr><th>X title</th><th>%s</th></tr>\n", xtitle);
  fprintf(fp, "<tr><th>Y title</th><th>%s</th></tr>\n", ytitle);
  if(ctime_r_int(&res->time, tmp))
    fprintf(fp, "<tr><th>Date</th><th>%s</th></tr>\n", tmp);
  fprintf(fp, "</tbody>\n");
  fprintf(fp, "</table>\n");

  /* Write the plots, each to their own table */
  for(ri=0; ri<dft->voiNr; ri++) {
    /* Search the largest regional plot x-value to be used as line end points */
    for(fi=0, maxRegX=maxFitX=0; fi<dft->frameNr; fi++) {
      if(dft->voi[ri].y2[fi]>maxRegX) maxRegX=dft->voi[ri].y2[fi];
      if(fi>=first && fi<=last && dft->voi[ri].y2[fi]>maxFitX)
        maxFitX=dft->voi[ri].y2[fi];
    }
    /* Begin a new table */
    fprintf(fp, "<table>\n");
    /* Write the title row */
    fprintf(fp, "<thead>\n");
    fprintf(fp, "<tr><th>%s %s %s</th>", dft->voi[ri].voiname,
            dft->voi[ri].hemisphere, dft->voi[ri].place);
    fprintf(fp, "<th>symbol open</th><th>symbol filled</th><th>text</th>");
    fprintf(fp, "<th>X</th><th>line</th>");
    fprintf(fp, "</tr>\n");
    fprintf(fp, "</thead>\n");
    /* Write the plot rows */
    fprintf(fp, "<tbody>\n");
    row=0;
    for(fi=0; fi<(dft->frameNr>2?dft->frameNr:2); fi++)
      if(!isnan(dft->voi[ri].y2[fi]) && !isnan(dft->voi[ri].y3[fi]))
    {
        fprintf(fp, "<tr>");
        if(fi<dft->frameNr) {
          fprintf(fp, "<th>%g</th>", dft->voi[ri].y2[fi]); /* x-axis value */
          fprintf(fp, "<th>%g</th>", dft->voi[ri].y3[fi]); /* y-axis value */
        } else {
          fprintf(fp, "<th> </th>");
          fprintf(fp, "<th> </th>");
        }
      /* If included in the fit, y-axis value again */
      if(fi>=first && fi<=last)
        fprintf(fp, "<th>%g</th>", dft->voi[ri].y3[fi]);
      else fprintf(fp, "<th></th>");
      if(fi<dft->frameNr)
        fprintf(fp, "<th>%g</th>", dft->x[fi]); /* Frame time as text */
      else
        fprintf(fp, "<th> </th>");
      /* Line points */
      if(row==0) { /* line start */
        fprintf(fp, "<th>0</th>"); /* x-axis value */
        fprintf(fp, "<th>%g</th>", res->voi[ri].parameter[1]); /* y-axis value */
      } else if(row==1) { /* line point at the end of fitted range */
        fprintf(fp, "<th>%g</th>", maxFitX); /* x-axis value */
        f=maxFitX*res->voi[ri].parameter[0]+res->voi[ri].parameter[1];
        fprintf(fp, "<th>%g</th>", f); /* y-axis value */
      } else if(row==2) { /* line "mid" point at the end of regional data */
        fprintf(fp, "<th>%g</th>", maxRegX); /* x-axis value */
        f=maxRegX*res->voi[ri].parameter[0]+res->voi[ri].parameter[1];
        fprintf(fp, "<th>%g</th>", f); /* y-axis value */
      } else if(row==3) { /* line end point at the end of all plots */
        fprintf(fp, "<th>%g</th>", maxPlotX); /* x-axis value */
        f=maxPlotX*res->voi[ri].parameter[0]+res->voi[ri].parameter[1];
        fprintf(fp, "<th>%g</th>", f); /* y-axis value */
      }
      fprintf(fp, "</tr>\n");
      row++;
    }
    fprintf(fp, "</tbody>\n");

    /* End the data table */
    fprintf(fp, "</table>\n");

  } /* next region plot */

  /* End the div for tables */
  fprintf(fp, "</div>\n");

  /* Stop writing the body of the HTML file, and end the file */
  n=fprintf(fp, "</body></html>\n");
  if(n==0) return(4);

  /* Close file */
  fclose(fp);

  return(0);
}
/******************************************************************************/

/*****************************************************************************/
/**  Write plot data in DFT format with x values as separate columns before corresponding y values.
    @return Returns 0 if successful, otherwise nonzero.
    @sa plotdata
 */
int plotdata_as_dft(
  /** Plot points: X in y2, Y in y3 */
  DFT *dft,
  /** Filename for plot data */
  char *fname
) {
  int ri, rj, fi, ret;
  DFT plot;

  /* Check input data */
  if(dft==NULL || dft->voiNr<1) return(1);
  /* Create the plot data */
  dftInit(&plot);
  ret=dftSetmem(&plot, dft->frameNr, 2*dft->voiNr); if(ret) return(ret);
  ret=dftCopymainhdr(dft, &plot); if(ret) {dftEmpty(&plot); return(ret);}
  for(ri=rj=0; ri<dft->voiNr; ri++) {
    /* x */
    strcpy(plot.voi[rj].voiname, "X");
    strcpy(plot.voi[rj].name, plot.voi[rj].voiname);
    for(fi=0; fi<dft->frameNr; fi++) plot.voi[rj].y[fi]=dft->voi[ri].y2[fi];
    rj++;
    /* y */
    dftCopyvoihdr(dft, ri, &plot, rj);
    for(fi=0; fi<dft->frameNr; fi++) plot.voi[rj].y[fi]=dft->voi[ri].y3[fi];
    rj++;
  }
  for(fi=0; fi<dft->frameNr; fi++) {
    plot.x[fi]=dft->x[fi]; plot.x1[fi]=dft->x1[fi]; plot.x2[fi]=dft->x2[fi];
  }
  plot.voiNr=2*dft->voiNr; plot.frameNr=dft->frameNr;
  /* Save plot data */
  strcpy(plot.comments, "");
  ret=dftWrite(&plot, fname);
  dftEmpty(&plot);
  return(ret);
}
/******************************************************************************/

/******************************************************************************/
