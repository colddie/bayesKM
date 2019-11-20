/// @file svg_legend.c
/// @author Vesa Oikonen
/// @brief Functions for drawing legends to SVG plots.
///
/*****************************************************************************/
#include "libtpcsvg.h"
/*****************************************************************************/
/** Write inline SVG (1) or separate SVG file (0) */
extern int SVG_INLINE;
/*****************************************************************************/

/*****************************************************************************/
/** Initiate SVG plot legends struct contents; call this once before usage */
void svg_init_legends(
  /** Pointer to legends struct */
  SVG_LEGENDS *legends
) {
  legends->_init=1;
  legends->n=0;
  legends->l=NULL;
}
/*****************************************************************************/

/*****************************************************************************/
/** Empty the legends struct contents and free the allocated memory. */
void svg_legend_empty(
  /** Pointer to legends struct */
  SVG_LEGENDS *legends
) {
  if(legends==NULL) return;
  if(legends->_init==0 || legends->n<1) return;
  free(legends->l);
  legends->n=0;
}
/*****************************************************************************/

/*****************************************************************************/
/** Add information of one legend item to legends struct.
 *  Memory will be allocated here.
\return Returns 0 when successful, otherwise <>0.
 */
int svg_legend_add(
  /** Pointer to legends structure. */
  SVG_LEGENDS *legends,
  /** Plot type: 1=line, 2=symbols, 0=both line and symbols. */
  const int plot_type,
  /** Symbol type: RECTANGLE,CIRCLE,UPTRIANGLE,DOWNTRIANGLE,DIAMOND,
                   LEFTTRIANGLE, RIGHTTRIANGLE */
  const int symbol_type,
  /** Symbol filling: SYMBOLOPEN, SYMBOLFILLED */
  const svgSymbolFill symbol_fill,
  /** SVG color index. */
  const int color,
  /** Pointer to Legend text. */
  const char *text
) {
  if(legends==NULL || legends->_init!=1) return 1;
  if(legends->n==0)
    legends->l=malloc((legends->n+1)*sizeof(SVG_LEGEND));
  else
    legends->l=realloc(legends->l, (legends->n+1)*sizeof(SVG_LEGEND));
  if(legends->l==NULL) {legends->n=0; return 2;}
  legends->n++;
  legends->l[legends->n-1].plot_type=plot_type;
  legends->l[legends->n-1].symbol_type=symbol_type;
  legends->l[legends->n-1].symbol_fill=symbol_fill;
  legends->l[legends->n-1].color=color;
  strncpy(legends->l[legends->n-1].text, text, MAX_SVG_LEGEND_LEN);
  legends->l[legends->n-1].text[MAX_SVG_LEGEND_LEN]=(char)0;
  return 0;
}
/*****************************************************************************/

/*****************************************************************************/
/** Create SVG plot legends.
\return Returns 0 if successful, <>0 in case of error.
*/
int svg_create_legends(
  /** SVG graphics file pointer. */
  FILE *fp,
  /** Struct containing the viewport sizes. */
  struct svg_viewports *vp,
  /** Pointer to struct containing legends. */
  SVG_LEGENDS *legends,
  /** Char pointer to string (at least of length 128) where possible
      error description is copied; set to NULL if not necessary. */
  char *errmsg,
  /** Verbose level; if zero, then nothing is printed to stderr or stdout. */
  int verbose
) {
  char tmp[1024], line[128], ilc[9];
  double xpos, ypos, ygap, size=100.0, trgsize=140.0, circsize=120.0;
  int ti, len, maxlen, text_space;
  const int defNr=24;
  const double hrratio=1.80; // average character height/width 

  if(verbose>0)
    printf("svg_create_legends(fp, vp, legends, errmsg, %d)\n", verbose);

  /* Check the input */
  if(vp->label_area_viewport.is==0) return(0);
  if(legends==NULL || legends->n<1) return(0);
  if(fp==NULL) {
    if(errmsg!=NULL) sprintf(errmsg, "file was closed too early");
    return(1);
  }

  if(SVG_INLINE) strcpy(ilc, "svg:"); else strcpy(ilc, "");

  /* Create a new viewport for plot legends */
  strcpy(tmp, "\n  <!-- Legends viewport -->\n");
  sprintf(line, "  <%ssvg x=\"%dpx\" y=\"%dpx\" width=\"%dpx\" height=\"%d\"",
    ilc, vp->label_area_viewport.x, vp->label_area_viewport.y,
    vp->label_area_viewport.w, vp->label_area_viewport.h); strcat(tmp, line);
  sprintf(line, "\n      viewBox=\"0 0 %d %d\"",
    vp->label_area_viewport.w, vp->label_area_viewport.h); strcat(tmp, line);
  strcat(tmp, "\n      preserveAspectRatio=\"xMidYMid meet\"");
  strcat(tmp, ">\n");
  if(svg_write(fp, tmp, errmsg, verbose-5)!=0) {return(4);}


#if(0) // during testing show the vieport with background color
  sprintf(tmp, "    <%srect width=\"%dpx\" height=\"%dpx\" stroke=\"none\" fill=\"red\" fill-opacity=\"0.3\" />\n",
    ilc, vp->label_area_viewport.w, vp->label_area_viewport.h);
  if(svg_write(fp, tmp, errmsg, verbose-5)!=0) {return(99);}
#endif

  /* Set the space for legend text */
  text_space=(8*vp->label_area_viewport.w)/10;

  /* Determine font size for legends */
  /* first, based on legend nr */
  ti=legends->n; if(ti<defNr) ti=defNr;
  vp->label_area_viewport.chr_size=
    (double)vp->label_area_viewport.h/(double)(ti+1);
  //printf("vp->label_area_viewport.chr_size := %d\n", vp->label_area_viewport.chr_size);
  /* then, make smaller based on max length, if necessary */
  ti=0; maxlen=len=strlen(legends->l[ti].text);
  for(ti=1; ti<legends->n; ti++) {
    len=strlen(legends->l[ti].text); if(len>maxlen) maxlen=len;}
  //printf("legend max length := %d\n", maxlen);
  if(vp->label_area_viewport.chr_size*maxlen>hrratio*(double)text_space)
    vp->label_area_viewport.chr_size=hrratio*(double)text_space/(double)maxlen;
  //printf("vp->label_area_viewport.chr_size := %d\n", vp->label_area_viewport.chr_size);
  /* Set line gap, if there is space for that */
  ygap=0; if(legends->n<=2*defNr/3) ygap=vp->label_area_viewport.chr_size/3;

  /* Write legend texts as a group */
  xpos=(double)(vp->label_area_viewport.w - text_space);
  ypos=1.5*(double)vp->label_area_viewport.chr_size;
  sprintf(tmp, "    <%sg", ilc);
  sprintf(line, " font-family=\"Sans-serif\""); strcat(tmp, line);
  sprintf(line, " text-anchor=\"Start\""); strcat(tmp, line);
  sprintf(line, " font-size=\"%d\"", vp->label_area_viewport.chr_size); strcat(tmp, line);
  sprintf(line, " fill=\"black\""); strcat(tmp, line);
  strcat(tmp, ">\n"); 
  if(svg_write(fp, tmp, errmsg, verbose-5)!=0) {return(8);}
  /* Write one legend text at a time */
  for(ti=0; ti<legends->n; ti++) {
    //printf("ti=%d text='%s'\n", ti, legends->l[ti].text);
    sprintf(tmp, "      <%s", ilc);
    sprintf(line, "text x=\"%g\" y=\"%g\"", xpos, ypos); strcat(tmp, line);
    sprintf(line, ">"); strcat(tmp, line);
    strcat(tmp, legends->l[ti].text);
    strcat(tmp, "</"); strcat(tmp, ilc); strcat(tmp, "text>\n");
    if(svg_write(fp, tmp, errmsg, verbose-5)!=0) {return(9);}
    ypos+=ygap+(double)vp->label_area_viewport.chr_size;
  }
  sprintf(tmp, "    </%sg>\n", ilc);
  if(svg_write(fp, tmp, errmsg, verbose-5)!=0) {return(10);}

  /* Write legend symbols */
  xpos=0.5*(double)(vp->label_area_viewport.w - text_space);
  ypos=1.25*(double)vp->label_area_viewport.chr_size;
  for(ti=0; ti<legends->n; ti++) {
    sprintf(tmp, "    <%sg", ilc);
    sprintf(line, " stroke=\"%s\"", svgColorName(legends->l[ti].color)); strcat(tmp, line);
    sprintf(line, " fill=\"%s\"", svgColorName(legends->l[ti].color)); strcat(tmp, line);
    if(legends->l[ti].symbol_fill==SYMBOLOPEN) sprintf(line, " fill-opacity=\"0.02\"");
    else sprintf(line, " fill-opacity=\"0.67\"");
    strcat(tmp, line);
    sprintf(line, " stroke-width=\"25\""); strcat(tmp, line);
    strcat(tmp, ">\n"); 
    if(svg_write(fp, tmp, errmsg, verbose-5)!=0) {return(8);}

    /* Line */
    if(legends->l[ti].plot_type==0 || legends->l[ti].plot_type==1) {
      sprintf(tmp, "      <%sline", ilc);
      sprintf(line, " x1=\"%g\"", 0.25*xpos); strcat(tmp, line);
      sprintf(line, " y1=\"%g\"", ypos); strcat(tmp, line);
      sprintf(line, " x2=\"%g\"", 1.75*xpos); strcat(tmp, line);
      sprintf(line, " y2=\"%g\"", ypos); strcat(tmp, line);
      strcat(tmp, " />\n");
      if(svg_write(fp, tmp, errmsg, verbose-5)!=0) {return(18);}
    }

    /* Symbol */
    if(legends->l[ti].plot_type==0 || legends->l[ti].plot_type==2) {
      sprintf(tmp, "      <%suse ", ilc);
      switch(legends->l[ti].symbol_type) {
        case RECTANGLE:
          sprintf(line, "xlink:href=\"#sym-rect\" x=\"%g\" y=\"%g\" width=\"%g\" height=\"%g\"",
            xpos-0.5*size, ypos-0.5*size, size, size);
          break;
        case UPTRIANGLE:
          sprintf(line, "xlink:href=\"#sym-uptr\" x=\"%g\" y=\"%g\" width=\"%g\" height=\"%g\"",
            xpos-0.5*trgsize, ypos-0.5*trgsize, trgsize, trgsize);
          break;
        case DOWNTRIANGLE:
          sprintf(line, "xlink:href=\"#sym-dotr\" x=\"%g\" y=\"%g\" width=\"%g\" height=\"%g\"",
            xpos-0.5*trgsize, ypos-0.5*trgsize, trgsize, trgsize);
          break;
        case DIAMOND:
          sprintf(line, "xlink:href=\"#sym-diam\" x=\"%g\" y=\"%g\" width=\"%g\" height=\"%g\"",
            xpos-0.5*trgsize, ypos-0.5*trgsize, trgsize, trgsize);
          break;
        case LEFTTRIANGLE:
          sprintf(line, "xlink:href=\"#sym-letr\" x=\"%g\" y=\"%g\" width=\"%g\" height=\"%g\"",
            xpos-0.5*trgsize, ypos-0.5*trgsize, trgsize, trgsize);
          break;
        case RIGHTTRIANGLE:
          sprintf(line, "xlink:href=\"#sym-ritr\" x=\"%g\" y=\"%g\" width=\"%g\" height=\"%g\"",
            xpos-0.5*trgsize, ypos-0.5*trgsize, trgsize, trgsize);
          break;
        case CIRCLE:
        default:
          sprintf(line, "xlink:href=\"#sym-circ\" x=\"%g\" y=\"%g\" width=\"%g\" height=\"%g\"",
            xpos-0.5*circsize, ypos-0.5*circsize, circsize, circsize);
          break;
      }
      strcat(tmp, line);
      strcat(tmp, " />\n");
      if(svg_write(fp, tmp, errmsg, verbose-5)!=0) {return(18);}
    }


    sprintf(tmp, "    </%sg>\n", ilc);
    if(svg_write(fp, tmp, errmsg, verbose-5)!=0) {return(10);}
    ypos+=ygap+(double)vp->label_area_viewport.chr_size;
  } // next legend

  /* Close the view port */
  sprintf(tmp, "  </%ssvg>\n", ilc);
  if(svg_write(fp, tmp, errmsg, verbose-5)!=0) {return(9);}

  return(0);
}
/*****************************************************************************/

/*****************************************************************************/
