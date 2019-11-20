/// @file svg_title.c
/// @author Vesa Oikonen
/// @brief Create SVG plot titles for TPC SVG C library.
///
/*****************************************************************************/
#include "libtpcsvg.h"
/*****************************************************************************/
/** Write inline SVG (1) or separate SVG file (0) */
extern int SVG_INLINE;
/*****************************************************************************/

/*****************************************************************************/
/** Create SVG plot main title.
\return Returns 0 if successful, <>0 in case of error.
*/
int svg_create_main_title(
  /** SVG graphics file pointer */
  FILE *fp,
  /** Text for main title */
  const char *main_title_text,
  /** Text for sub title */
  const char *sub_title_text,
  /** Struct containing the viewport sizes */
  struct svg_viewports *vp,
  /** Char pointer to string (at least of length 128) where possible
      error description is copied; set to NULL if not necessary */
  char *errmsg,
  /** Verbose level; if zero, then nothing is printed to stderr or stdout */
  int verbose
) {
  char tmp[1024], line[128], ilc[9];
  double main_pos, sub_pos;

  if(verbose>0)
    printf("svg_create_main_title(fp, mtt, stt, vp, errmsg, %d)\n", verbose);
  if(verbose>1) {
    printf("main_title_text := '%s'\n", main_title_text);
    printf("sub_title_text := '%s'\n", sub_title_text);
  }

  /* Check the input */
  if(vp->main_title_viewport.is==0) return(0);
  if(fp==NULL) {
    if(errmsg!=NULL) sprintf(errmsg, "file was closed too early");
    return(1);
  }

  if(SVG_INLINE) strcpy(ilc, "svg:"); else strcpy(ilc, "");

  /* Set title of the SVG document */
  sprintf(tmp, "\n  <%stitle>%s</%stitle>\n", ilc, main_title_text, ilc);
  if(svg_write(fp, tmp, errmsg, verbose-5)!=0) {return(4);}

  /* Create a new viewport for main title */
  strcpy(tmp, "\n  <!-- Main title viewport -->\n");
  sprintf(line, "  <%ssvg x=\"%dpx\" y=\"%dpx\" width=\"%dpx\" height=\"%d\"",
    ilc, vp->main_title_viewport.x, vp->main_title_viewport.y,
    vp->main_title_viewport.w, vp->main_title_viewport.h); strcat(tmp, line);
  sprintf(line, "\n      viewBox=\"0 0 %d %d\"",
    vp->main_title_viewport.w, vp->main_title_viewport.h); strcat(tmp, line);
  strcat(tmp, "\n      preserveAspectRatio=\"xMidYMid meet\"");
  strcat(tmp, ">\n");
  if(svg_write(fp, tmp, errmsg, verbose-5)!=0) {return(4);}


#if(0)
  sprintf(tmp, "    <rect width=\"%dpx\" height=\"%dpx\" stroke=\"none\" fill=\"lime\" fill-opacity=\"0.3\" />\n",
    vp->main_title_viewport.w, vp->main_title_viewport.h);
  if(svg_write(fp, tmp, errmsg, verbose-5)!=0) {return(99);}
#endif

  /* Determine the y positions for titles */
  if(strlen(main_title_text)>0) { /* main title exists */
    if(strlen(sub_title_text)==0) { /* no subtitle */
      main_pos=0.75*(double)vp->main_title_viewport.h;
      sub_pos=vp->main_title_viewport.h;
    } else { /* also subtitle */
      main_pos=0.52*(double)vp->main_title_viewport.h;
      sub_pos=0.9*(double)vp->main_title_viewport.h;
    }
  } else { /* no main title */
    if(strlen(sub_title_text)==0) { /* no subtitle either */
      main_pos=0.5*(double)vp->main_title_viewport.h;
      sub_pos=vp->main_title_viewport.h;
    } else { /* only subtitle */
      main_pos=0.0;
      sub_pos=0.4*(double)vp->main_title_viewport.h;
    }
  }

  /* Set main title text */
  sprintf(tmp, "    <%stext x=\"%d\" y=\"%g\"\n",
    ilc, vp->main_title_viewport.w/2, main_pos);
  sprintf(line, "        font-family=\"Sans-serif\""); strcat(tmp, line);
  sprintf(line, " text-anchor=\"middle\""); strcat(tmp, line);
  sprintf(line, " font-size=\"%d\"",
    vp->main_title_viewport.chr_size); strcat(tmp, line);
  sprintf(line, " fill=\"black\""); strcat(tmp, line);
  sprintf(line, ">\n"); strcat(tmp, line);
  sprintf(line, "      %s\n", main_title_text); strcat(tmp, line);
  sprintf(line, "    </%stext>\n", ilc); strcat(tmp, line);
  if(svg_write(fp, tmp, errmsg, verbose-5)!=0) {return(6);}

  /* Set sub title text */
  sprintf(tmp, "    <%stext x=\"%d\" y=\"%g\"\n",
    ilc, vp->main_title_viewport.w/2, sub_pos);
  sprintf(line, "        font-family=\"Sans-serif\""); strcat(tmp, line);
  sprintf(line, " text-anchor=\"middle\""); strcat(tmp, line);
  sprintf(line, " font-size=\"%d\"",
    2*vp->main_title_viewport.chr_size/3); strcat(tmp, line);
  sprintf(line, " fill=\"black\""); strcat(tmp, line);
  sprintf(line, ">\n"); strcat(tmp, line);
  sprintf(line, "      %s\n", sub_title_text); strcat(tmp, line);
  sprintf(line, "    </%stext>\n", ilc); strcat(tmp, line);
  if(svg_write(fp, tmp, errmsg, verbose-5)!=0) {return(8);}

  /* Close the view port */
  sprintf(tmp, "  </%ssvg>\n", ilc);
  if(svg_write(fp, tmp, errmsg, verbose-5)!=0) {return(9);}

  return(0);
}
/*****************************************************************************/

/*****************************************************************************/
/** Create SVG plot x axis title.
\return Returns 0 if successful, <>0 in case of error.
*/
int svg_create_xaxis_title(
  /** SVG graphics file pointer */
  FILE *fp,
  /** Text for x axis title */
  const char *title_text,
  /** Struct containing the viewport sizes */
  struct svg_viewports *vp,
  /** Char pointer to string (at least of length 128) where possible
      error description is copied; set to NULL if not necessary */
  char *errmsg,
  /** Verbose level; if zero, then nothing is printed to stderr or stdout */
  int verbose
) {
  char tmp[1024], line[128], ilc[9];

  if(verbose>0)
    printf("svg_create_xaxis_title(fp, tt, vp, errmsg, %d)\n", verbose);
  if(verbose>1)
    printf("title_text := '%s'\n", title_text);

  /* Check the input */
  if(vp->xaxis_title_viewport.is==0) return(0);
  if(fp==NULL) {
    if(errmsg!=NULL) sprintf(errmsg, "file was closed too early");
    return(1);
  }

  if(SVG_INLINE) strcpy(ilc, "svg:"); else strcpy(ilc, "");

  /* Create a new viewport for x axis title */
  strcpy(tmp, "\n  <!-- X axis title viewport -->\n");
  sprintf(line, "  <%ssvg x=\"%dpx\" y=\"%dpx\" width=\"%dpx\" height=\"%d\"",
    ilc, vp->xaxis_title_viewport.x, vp->xaxis_title_viewport.y,
    vp->xaxis_title_viewport.w, vp->xaxis_title_viewport.h); strcat(tmp, line);
  sprintf(line, "\n      viewBox=\"0 0 %d %d\"",
    vp->xaxis_title_viewport.w, vp->xaxis_title_viewport.h); strcat(tmp, line);
  strcat(tmp, "\n      preserveAspectRatio=\"xMidYMid meet\"");
  strcat(tmp, ">\n");
  if(svg_write(fp, tmp, errmsg, verbose-5)!=0) {return(4);}


#if(0)
  sprintf(tmp, "    <rect width=\"%dpx\" height=\"%dpx\" stroke=\"none\" fill=\"green\" fill-opacity=\"0.3\" />\n",
    vp->xaxis_title_viewport.w, vp->xaxis_title_viewport.h);
  if(svg_write(fp, tmp, errmsg, verbose-5)!=0) {return(99);}
#endif


  /* Set the text contents */
  sprintf(tmp, "    <%stext x=\"%d\" y=\"%g\"\n",
    ilc, vp->xaxis_title_viewport.w/2, 0.75*(double)vp->xaxis_title_viewport.h);
  sprintf(line, "        font-family=\"Sans-serif\""); strcat(tmp, line);
  sprintf(line, " text-anchor=\"middle\""); strcat(tmp, line);
  sprintf(line, " font-size=\"%d\"",
    vp->xaxis_title_viewport.chr_size); strcat(tmp, line);
  sprintf(line, " fill=\"black\""); strcat(tmp, line);
  sprintf(line, ">\n"); strcat(tmp, line);
  sprintf(line, "      %s\n", title_text); strcat(tmp, line);
  sprintf(line, "    </%stext>\n", ilc); strcat(tmp, line);
  if(svg_write(fp, tmp, errmsg, verbose-5)!=0) {return(6);}

  /* Close the view port */
  sprintf(tmp, "  </%ssvg>\n", ilc);
  if(svg_write(fp, tmp, errmsg, verbose-5)!=0) {return(9);}

  return(0);
}
/*****************************************************************************/

/*****************************************************************************/
/** Create SVG plot y axis title.
\return Returns 0 if successful, <>0 in case of error.
*/
int svg_create_yaxis_title(
  /** SVG graphics file pointer */
  FILE *fp,
  /** Text for y axis title */
  const char *title_text,
  /** Struct containing the viewport sizes */
  struct svg_viewports *vp,
  /** Char pointer to string (at least of length 128) where possible
      error description is copied; set to NULL if not necessary */
  char *errmsg,
  /** Verbose level; if zero, then nothing is printed to stderr or stdout */
  int verbose
) {
  char tmp[1024], line[128], ilc[9];
  double xpos, ypos;

  if(verbose>0)
    printf("svg_create_yaxis_title(fp, tt, vp, errmsg, %d)\n", verbose);
  if(verbose>1)
    printf("title_text := '%s'\n", title_text);

  /* Check the input */
  if(vp->yaxis_title_viewport.is==0) return(0);
  if(fp==NULL) {
    if(errmsg!=NULL) sprintf(errmsg, "file was closed too early");
    return(1);
  }

  if(SVG_INLINE) strcpy(ilc, "svg:"); else strcpy(ilc, "");

  /* Create a new viewport for y axis title */
  strcpy(tmp, "\n  <!-- Y axis title viewport -->\n");
  sprintf(line, "  <%ssvg x=\"%dpx\" y=\"%dpx\" width=\"%dpx\" height=\"%d\"",
    ilc, vp->yaxis_title_viewport.x, vp->yaxis_title_viewport.y,
    vp->yaxis_title_viewport.w, vp->yaxis_title_viewport.h); strcat(tmp, line);
  sprintf(line, "\n      viewBox=\"0 0 %d %d\"",
    vp->yaxis_title_viewport.w, vp->yaxis_title_viewport.h); strcat(tmp, line);
  strcat(tmp, "\n      preserveAspectRatio=\"xMidYMid meet\"");
  strcat(tmp, ">\n");
  if(svg_write(fp, tmp, errmsg, verbose-5)!=0) {return(4);}


#if(0)
  sprintf(tmp, "    <rect width=\"%dpx\" height=\"%dpx\" stroke=\"none\" fill=\"red\" fill-opacity=\"0.3\" />\n",
    vp->yaxis_title_viewport.w, vp->yaxis_title_viewport.h);
  if(svg_write(fp, tmp, errmsg, verbose-5)!=0) {return(99);}
#endif

  /* Set the text contents */
  xpos=0.75*(double)vp->yaxis_title_viewport.w;
  ypos=0.5*(double)vp->yaxis_title_viewport.h;
  sprintf(tmp, "    <%stext x=\"%g\" y=\"%g\"\n", ilc, xpos, ypos);
  sprintf(line, "        font-family=\"Sans-serif\""); strcat(tmp, line);
  sprintf(line, " text-anchor=\"middle\""); strcat(tmp, line);
  sprintf(line, " font-size=\"%d\"",
    vp->yaxis_title_viewport.chr_size); strcat(tmp, line);
  sprintf(line, " fill=\"black\""); strcat(tmp, line);
  sprintf(line, " transform=\"rotate(270,%g,%g)\"", xpos, ypos); strcat(tmp, line);
  sprintf(line, ">\n"); strcat(tmp, line);
  sprintf(line, "      %s\n", title_text); strcat(tmp, line);
  sprintf(line, "    </%stext>\n", ilc); strcat(tmp, line);
  if(svg_write(fp, tmp, errmsg, verbose-5)!=0) {return(6);}

  /* Close the view port */
  sprintf(tmp, "  </%ssvg>\n", ilc);
  if(svg_write(fp, tmp, errmsg, verbose-5)!=0) {return(9);}

  return(0);
}
/*****************************************************************************/

/*****************************************************************************/

