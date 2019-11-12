/// @file svg_defs.c
/// @author Vesa Oikonen
/// @brief Defined objects for TPC SVG C library.
///
/*****************************************************************************/
#include "libtpcsvg.h"
/*****************************************************************************/
/** Write inline SVG (1) or separate SVG file (0) */
extern int SVG_INLINE;
/*****************************************************************************/

/*****************************************************************************/
/** SVG colors used by libtpcsvg; corresponding svgColor enums as comments */
static char *svgcolor[] = {
  /* BLACK   */ "black",
  /* RED     */ "red",
  /* BLUE    */ "blue",
  /* GREEN   */ "green",
  /* PURPLE  */ "purple",
  /* OLIVE   */ "olive",
  /* AQUA    */ "aqua", 
  /* FUCHSIA */ "fuchsia",
  /* GRAY    */ "gray",
  /* LIME    */ "lime",
  /* MAROON  */ "maroon",
  /* NAVY    */ "navy",
  /* SILVER  */ "silver",
  /* TEAL    */ "teal", 
  /* YELLOW  */ "yellow",
  0
};
/*****************************************************************************/

/*****************************************************************************/
/** Return pointer to string describing the color, or NULL if outside
 *  of limits */
char *svgColorName(
  /** index of color */
  const svgColor i
) {
  unsigned int n=0;
  while(svgcolor[n]!=0) n++;
  if(/*i<0 ||*/ i>n-1) return(NULL); //return(svgcolor[BLACK]);
  else return(svgcolor[i]);
}
/*****************************************************************************/

/*****************************************************************************/
/** Plot symbols used by libtpcsvg; corresponding svgSymbolType enums
 *  as comments */
static char *svgsymbol[] = {
  /* RECTANGLE     */ "rect",
  /* CIRCLE        */ "circ",
  /* UPTRIANGLE    */ "uptr",
  /* DOWNTRIANGLE  */ "dotr",
  /* DIAMOND       */ "diam",
  /* LEFTTRIANGLE  */ "letr",
  /* RIGHTTRIANGLE */ "ritr", 
  0
};
/*****************************************************************************/

/*****************************************************************************/
/** Return pointer to string describing the symbol, or NULL if outside
 *  of limits */
char *svgSymbolName(
  /** index of symbol */
  const svgSymbolType i
) {
  unsigned int n=0;
  while(svgsymbol[n]!=0) n++;
  if(/*i<0 ||*/ i>n-1) return(NULL);
  else return(svgsymbol[i]);
}
/*****************************************************************************/

/*****************************************************************************/
/** Define the curve symbols for SVG graphics file.
\return Returns 0 if successful, <>0 in case of error.
 */
int svg_define_symbols(
  /** SVG graphics file pointer */
  FILE *fp,
  /** Char pointer to string (at least of length 128) where possible
      error description is copied; set to NULL if not necessary */
  char *errmsg,
  /** Verbose level; if zero, then nothing is printed to stderr or stdout */
  int verbose
) {
  char tmp[2048], ilc[9], line[256], allsame[128];

  if(verbose>0) printf("svg_define_symbols(fp, errmsg, %d)\n", verbose);

  if(SVG_INLINE) strcpy(ilc, "svg:"); else strcpy(ilc, "");
  strcpy(allsame, "viewBox=\"0 0 120 120\" preserveAspectRatio=\"xMinYMin meet\"");

  strcpy(tmp, "");
  sprintf(line, "  <%sdefs>\n", ilc); strcat(tmp, line);

  sprintf(line, "    <%ssymbol id=\"sym-rect\" %s >\n", ilc, allsame);
  strcat(tmp, line);
  sprintf(line, "      <%srect x=\"10\" y=\"10\" width=\"100\" height=\"100\" />\n", ilc);
  strcat(tmp, line);
  sprintf(line, "    </%ssymbol>\n", ilc); strcat(tmp, line);

  sprintf(line, "    <%ssymbol id=\"sym-circ\" %s >\n", ilc, allsame);
  strcat(tmp, line);
  sprintf(line, "      <%scircle cx=\"60\" cy=\"60\" r=\"50\" />\n", ilc);
  strcat(tmp, line);
  sprintf(line, "    </%ssymbol>\n", ilc); strcat(tmp, line);

  sprintf(line, "    <%ssymbol id=\"sym-uptr\" %s >\n", ilc, allsame);
  strcat(tmp, line);
  sprintf(line, "      <%spolygon points=\" 10 17, 110 17, 60 103 \" />\n", ilc);
  strcat(tmp, line);
  sprintf(line, "    </%ssymbol>\n", ilc); strcat(tmp, line);

  sprintf(line, "    <%ssymbol id=\"sym-dotr\" %s >\n", ilc, allsame);
  strcat(tmp, line);
  sprintf(line, "      <%spolygon points=\" 10 103, 110 103, 60 17 \" />\n", ilc);
  strcat(tmp, line);
  sprintf(line, "    </%ssymbol>\n", ilc); strcat(tmp, line);

  sprintf(line, "    <%ssymbol id=\"sym-letr\" %s >\n", ilc, allsame);
  strcat(tmp, line);
  sprintf(line, "      <%spolygon points=\" 103 10, 103 110, -103 60 \" />\n", ilc);
  strcat(tmp, line);
  sprintf(line, "    </%ssymbol>\n", ilc); strcat(tmp, line);

  sprintf(line, "    <%ssymbol id=\"sym-ritr\" %s >\n", ilc, allsame);
  strcat(tmp, line);
  sprintf(line, "      <%spolygon points=\" 17 10, 17 110, 103 60 \" />\n", ilc);
  strcat(tmp, line);
  sprintf(line, "    </%ssymbol>\n", ilc); strcat(tmp, line);

  sprintf(line, "    <%ssymbol id=\"sym-diam\" %s >\n", ilc, allsame);
  strcat(tmp, line);
  sprintf(line, "      <%spolygon points=\" 60 10, 110 60, 60 110, 10 60 \" />\n", ilc);
  strcat(tmp, line);
  sprintf(line, "    </%ssymbol>\n", ilc); strcat(tmp, line);

  sprintf(line, "  </%sdefs>\n", ilc); strcat(tmp, line);

  return(svg_write(fp, tmp, errmsg, verbose-5));
}
/*****************************************************************************/

/*****************************************************************************/
