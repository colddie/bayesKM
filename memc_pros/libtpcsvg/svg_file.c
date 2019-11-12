/// @file svg_file.c
/// @brief File io for TPC SVG C library.
/// @author Vesa Oikonen
/// @todo If Inkscape, Gimp and Batik start to support plot data outside of
///       viewport, then decide whether to write that data in SVG again.
///       Embedded SVG in HTML5.
///
/*****************************************************************************/
#include "libtpcsvg.h"
/*****************************************************************************/
/** Write inline SVG (1) or separate SVG file (0) */
int SVG_INLINE = 0;
/*****************************************************************************/

/*****************************************************************************/
/** Initiate a new SVG graphics file.

    If file with same name exists, it is overwritten without backup.
   @return Returns pointer to the file if successful and NULL in case of an error.
   @sa svg_close, svg_xhtml_initiate
 */
FILE *svg_initiate(
  /** File name for SVG graphics */
  const char *filename,
  /** Plot height in cm; 0, if not predefined */
  const double height,
  /** Plot width in cm; 0, if not predefined */
  const double width,
  /** Struct containing the viewport sizes */
  struct svg_viewports *vp,
  /** Char pointer to string (at least of length 128) where possible
      error description is copied; set to NULL if not necessary */
  char *errmsg,
  /** Verbose level; if zero, then nothing is printed to stderr or stdout */
  int verbose
) {
  FILE *fp;
  char tmp[1024], *cptr;

  if(verbose>0)
    printf("svg_initiate(%s, %g, %g, vp, errmsg, %d)\n",
      filename, height, width, verbose);

  /* Check input */
  if(filename==NULL || vp==NULL) return(NULL);
  if(vp->main_viewport.w<3 || vp->main_viewport.h<3) return(NULL);

  /* Open file for write */
  fp=fopen(filename, "w");
  if(fp==NULL) {
    if(errmsg!=NULL) strcpy(errmsg, "cannot open file for write");
    return(fp);
  }

  strcpy(tmp, "<?xml version=\"1.0\" encoding=\"utf-8\"?>\n");
  //strcat(tmp, "<!DOCTYPE svg PUBLIC \"-//W3C//DTD SVG 1.1//EN\"\n");
  //strcat(tmp, " \"http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd\">\n");
  strcat(tmp, "<svg version=\"1.1\" baseProfile=\"full\"\n");
  strcat(tmp, "     xmlns=\"http://www.w3.org/2000/svg\"\n");
  strcat(tmp, "     xmlns:xlink=\"http://www.w3.org/1999/xlink\"\n");
  strcat(tmp, "     xmlns:ev=\"http://www.w3.org/2001/xml-events\"");
  if(svg_write(fp, tmp, errmsg, verbose-5)!=0) {
    fclose(fp); remove(filename); return(NULL);}

  sprintf(tmp, "\n    viewBox=\"0 0 %d %d\"",
    vp->main_viewport.w, vp->main_viewport.h);
  strcat(tmp,  "\n    preserveAspectRatio=\"xMinYMin meet\"");
  if(svg_write(fp, tmp, errmsg, verbose-5)!=0) {
    fclose(fp); remove(filename); return(NULL);}

  if(width>0.0) {
    sprintf(tmp, "\n     width=\"%gcm\"", width);
    if(svg_write(fp, tmp, errmsg, verbose-5)!=0) {
      fclose(fp); remove(filename); return(NULL);}
  }
  if(height>0.0) {
    sprintf(tmp, "\n     height=\"%gcm\"", height);
    if(svg_write(fp, tmp, errmsg, verbose-5)!=0) {
      fclose(fp); remove(filename); return(NULL);}
  }

  strcpy(tmp, ">\n");
  if(svg_write(fp, tmp, errmsg, verbose-5)!=0) {
    fclose(fp); remove(filename); return(NULL);}

  /* Write file name as title element */
  cptr=strrchr(filename, '/'); if(cptr==NULL) cptr=strrchr(filename, '\\');
  if(cptr!=NULL) cptr++; else cptr=(char*)filename;
  if(cptr!=NULL) {
    sprintf(tmp, "  <title>%s</title>\n", cptr);
    if(svg_write(fp, tmp, errmsg, verbose-5)!=0) {
      fclose(fp); remove(filename); return(NULL);}
  }

  /* Create plot symbols for possible later use */
  if(svg_define_symbols(fp, errmsg, verbose)!=0) {
    fclose(fp); remove(filename); return(NULL);}

  return(fp);
}
/*****************************************************************************/

/*****************************************************************************/
/** Close SVG graphics file.
   @return Returns 0 if successful, <>0 in case of error.
   @sa svg_initiate
*/
int svg_close(
  /** SVG graphics file pointer */
  FILE *fp,
  /** Char pointer to string (at least of length 128) where possible
      error description is copied; set to NULL if not necessary */
  char *errmsg,
  /** Verbose level; if zero, then nothing is printed to stderr or stdout */
  int verbose
) {
  if(verbose>0) printf("svg_close(fp, errmsg, %d)\n", verbose);

  if(svg_write(fp, "</svg>\n", errmsg, verbose-5)!=0) {fclose(fp); return(1);}
  fflush(fp); fclose(fp);
  return(0);
}
/*****************************************************************************/

/*****************************************************************************/
/** Initiate a new XHTML file for one or more inline SVG graphics files.

    If file with same name exists, it is overwritten without backup.
   @return Returns pointer to the file if successful and NULL in case of an error.
   @sa svg_initiate, svg_xhtml_close, svg_xhtml_svg_open
 */
FILE *svg_xhtml_initiate(
  /** File name for SVG graphics */
  const char *filename,
  /** XHTML title; if NULL, then filename is used */
  const char *XHTML_title,
  /** Char pointer to string (at least of length 128) where possible
      error description is copied; set to NULL if not necessary */
  char *errmsg,
  /** Verbose level; if zero, then nothing is printed to stderr or stdout */
  int verbose
) {
  FILE *fp;
  char tmp[2048], line[256];

  if(verbose>0)
    printf("svg_xhtml_initiate(%s, %s, errmsg, %d)\n",
      filename, XHTML_title, verbose);

  SVG_INLINE = 1;

  /* Check input */
  if(filename==NULL) return(NULL);
  if(XHTML_title==NULL) XHTML_title=filename;

  /* Open file for write */
  fp=fopen(filename, "w");
  if(fp==NULL) {
    if(errmsg!=NULL) strcpy(errmsg, "cannot open file for write");
    return(fp);
  }

  strcpy(tmp, "<?xml version=\"1.0\" encoding=\"iso-8859-1\"?>\n");
  if(svg_write(fp, tmp, errmsg, verbose-5)!=0) {
    fclose(fp); remove(filename); return(NULL);}

  strcpy(tmp, "<!DOCTYPE html PUBLIC\n");
  strcat(tmp, "     \"-//W3C//DTD XHTML 1.0 Strict//EN\"\n");
  strcat(tmp, "     \"http://www.w3.org/TR/xhtml1/DTD/xhtml1-strict.dtd\">\n");
  if(svg_write(fp, tmp, errmsg, verbose-5)!=0) {
    fclose(fp); remove(filename); return(NULL);}

  strcpy(tmp, "<html xmlns=\"http://www.w3.org/1999/xhtml\"\n");
  strcat(tmp, "     xmlns:svg=\"http://www.w3.org/2000/svg\"\n");
  strcat(tmp, "     xmlns:xlink=\"http://www.w3.org/1999/xlink\"\n");
  strcat(tmp, "     xmlns:ev=\"http://www.w3.org/2001/xml-events\"\n");
  strcat(tmp, "     xml:lang=\"en\" lang=\"en\">\n");
  if(svg_write(fp, tmp, errmsg, verbose-5)!=0) {
    fclose(fp); remove(filename); return(NULL);}

  strcpy(tmp,   "<head>\n");
  sprintf(line, "  <title>%s</title>\n", XHTML_title); strcat(tmp, line);
  sprintf(line, "  <meta http-equiv=\"content-type\" content=\"text/html; charset=iso-8859-1\" />\n"); strcat(tmp, line);
  sprintf(line, "  <meta http-equiv=\"content-language\" content=\"en-gb\" />\n"); strcat(tmp, line);
  strcat(tmp,   "  <object id=\"AdobeSVG\" classid=\"clsid:78156a80-c6a1-4bbf-8e6a-3cd390eeb4e2\"> </object>\n");
  strcat(tmp,   "  <?import namespace=\"svg\" urn=\"http://www.w3.org/2000/svg\" implementation=\"#AdobeSVG\"?>\n");
  strcat(tmp,   "</head>\n\n");
  if(svg_write(fp, tmp, errmsg, verbose-5)!=0) {
    fclose(fp); remove(filename); return(NULL);}

  strcpy(tmp, "<body>\n\n");
  if(svg_write(fp, tmp, errmsg, verbose-5)!=0) {
    fclose(fp); remove(filename); return(NULL);}

  /* Create plot symbols for possible later use */
  if(svg_define_symbols(fp, errmsg, verbose)!=0) {
    fclose(fp); remove(filename); return(NULL);}

  return(fp);
}
/*****************************************************************************/

/*****************************************************************************/
/** Close XHTML file containing inline SVG.
   @return Returns 0 if successful, <>0 in case of error.
   @sa svg_xhtml_initiate
*/
int svg_xhtml_close(
  /** SVG graphics file pointer */
  FILE *fp,
  /** Char pointer to string (at least of length 128) where possible
      error description is copied; set to NULL if not necessary */
  char *errmsg,
  /** Verbose level; if zero, then nothing is printed to stderr or stdout */
  int verbose
) {
  if(verbose>0) printf("svg_xhtml_close(fp, errmsg, %d)\n", verbose);

  if(svg_write(fp, "</body>\n</html>\n", errmsg, verbose-5)!=0) {
    fclose(fp); return(1);}
  fflush(fp); fclose(fp);

  SVG_INLINE = 0;

  return(0);
}
/*****************************************************************************/

/*****************************************************************************/
/** Open a new SVG inline XHTML file.
    @return Returns 0 if successful, <>0 in case of error.
   @sa svg_xhtml_svg_close, svg_xhtml_initiate
 */
int svg_xhtml_svg_open(
  /** SVG graphics file pointer */
  FILE *fp,
  /** Plot height in cm; 0, if not predefined */
  const double height,
  /** Plot width in cm; 0, if not predefined */
  const double width,
  /** Struct containing the viewport sizes */
  struct svg_viewports *vp,
  /** Char pointer to string (at least of length 128) where possible
      error description is copied; set to NULL if not necessary */
  char *errmsg,
  /** Verbose level; if zero, then nothing is printed to stderr or stdout */
  int verbose
) {
  char tmp[1024];

  if(verbose>0)
    printf("svg_xhtml_svg_open(fp, %g, %g, vp, errmsg, %d)\n",
      height, width, verbose);

  /* Check input */
  if(fp==NULL || vp==NULL) return(1);
  if(vp->main_viewport.w<3 || vp->main_viewport.h<3) return(1);

  strcpy(tmp, "<svg:svg version=\"1.1\" baseProfile=\"full\"");
  if(svg_write(fp, tmp, errmsg, verbose-5)!=0) return(2);

  sprintf(tmp, "\n    viewBox=\"0 0 %d %d\"", vp->main_viewport.w, vp->main_viewport.h);
  strcat(tmp,  "\n    preserveAspectRatio=\"xMinYMin meet\"");
  if(svg_write(fp, tmp, errmsg, verbose-5)!=0) return(3);

  if(width>0.0) {
    sprintf(tmp, "\n     width=\"%gcm\"", width);
    if(svg_write(fp, tmp, errmsg, verbose-5)!=0) return(4);
  }
  if(height>0.0) {
    sprintf(tmp, "\n     height=\"%gcm\"", height);
    if(svg_write(fp, tmp, errmsg, verbose-5)!=0) return(5);
  }

  strcpy(tmp, ">\n");
  if(svg_write(fp, tmp, errmsg, verbose-5)!=0) return(6);

  return(0);
}
/*****************************************************************************/

/*****************************************************************************/
/** Close SVG graphics inline in XHTML file. Leaves the file open.
   @return Returns 0 if successful, <>0 in case of error.
   @sa svg_xhtml_svg_open, svg_write
*/
int svg_xhtml_svg_close(
  /** SVG graphics file pointer */
  FILE *fp,
  /** Char pointer to string (at least of length 128) where possible
      error description is copied; set to NULL if not necessary */
  char *errmsg,
  /** Verbose level; if zero, then nothing is printed to stderr or stdout */
  int verbose
) {
  if(verbose>0) printf("svg_xhtml_svg_close(fp, errmsg, %d)\n", verbose);

  if(svg_write(fp, "</svg:svg>\n", errmsg, verbose-5)!=0) return(1);
  return(0);
}
/*****************************************************************************/

/*****************************************************************************/
/** Write given string into open SVG file
   @return Returns 0 if successful, <>0 in case of an error.
   @sa svg_initiate, svg_xhtml_initiate, svg_xhtml_svg_open
*/
int svg_write(
  /** SVG graphics file pointer */
  FILE *fp,
  /** Char pointer to NULL terminated string to be written into file */
  const char *svg_string,
  /** Char pointer to string (at least of length 128) where possible
      error description is copied; set to NULL if not necessary */
  char *errmsg,
  /** Verbose level; if zero, then nothing is printed to stderr or stdout */
  int verbose
) {
  int pnr, len;

  if(verbose>0) printf("svg_write(fp, svg_string, errmsg, %d)\n", verbose);
  if(verbose>1) printf("svg_string := %s\n", svg_string);

  len=strlen(svg_string); if(len<1) return(0);
  if(fp==NULL) {
    if(errmsg!=NULL) sprintf(errmsg, "file was closed too early");
    return(1);
  }
  pnr=fprintf(fp, "%s", svg_string);
  if(pnr<len) {
    if(errmsg!=NULL) sprintf(errmsg, "cannot write into file");
    return(2);
  }
  return(0);
}
/*****************************************************************************/

/*****************************************************************************/

