/** @file proginfo.c
 *  @author Vesa Oikonen
 *  @brief Functions for printing usage and build information from executables.
 */
/*****************************************************************************/
#include "libtpcmisc.h"
/*****************************************************************************/

/*****************************************************************************/
static char *tpclicense[] = {
  "This program comes with ABSOLUTELY NO WARRANTY.",
  "This is free software, and you are welcome to redistribute it",
  "under GNU General Public License.",
  0};
static char *tpcstdoptions[] = {
  " -h, --help",
  "     Display usage information on standard output and exit.",
  " -v, --version",
  "     Display version and compile information on standard output and exit.",
  " -d[n], --debug[=n], --verbose[=n]",
  "     Set the level (n) of debugging messages and listings.",
  " -q, --quiet",
  "     Suppress displaying normal results on standard output.",
  " -s, --silent",
  "     Suppress displaying anything except errors.",
  0};
/*****************************************************************************/

/*****************************************************************************/
/** Check if given command-line argument string is one of the standard
 *  command-line options of this project. 
 *  @return 0 if string was identified as standard option, otherwise 1.
 *  @author Vesa Oikonen
 */
int tpcProcessStdOptions(
  /** Pointer to command-line option string */
  const char *s,
  /** If option string is either -h or --help, then this variable is set to 1 */
  int *print_usage,  
  /** If option string is either -v, -V, or --version, then this variable is
   *  set to 1 */
  int *print_version,  
  /** The level of debugging messages and listings:
   *  - If option string is -d, --debug or --verbose, then +1 is added
   *    to this variable.
   *  - If option string is -d[n], --debug[=n], or --verbose[=n], then +n
   *    is added to this variable.
   *  - If options string is -q or --quiet, then this variable is set to 0.         
   *  - If options string is -s or --silent, then this variable is set to -1.         
   */
  int *verbose_level
) {
  int n;

  char *cptr;
  /* Check that string is option, starting with '-' or '--' */
  if(s==NULL || strlen(s)<2 || s[0]!='-') return 1;
  /* Set pointer to the character after the first '-' */
  cptr=(char*)s+1;
  /* If also the next character is '-', then try the long forms of options */
  if(*cptr=='-') {
    cptr++; if(strlen(cptr)<1) return 1;
    if(strcasecmp(cptr, "help")==0) {*print_usage=1; return 0;}
    if(strcasecmp(cptr, "helphtml")==0) {*print_usage=2; return 0;}
    if(strcasecmp(cptr, "version")==0) {*print_version=1; return 0;}
    if(strcasecmp(cptr, "debug")==0) {*verbose_level+=1; return 0;}
    if(strcasecmp(cptr, "verbose")==0) {*verbose_level+=1; return 0;}
    if(strncasecmp(cptr, "debug=", 6)==0) {
      if(!isdigit(cptr[6])) return 1;
      n=atoi(cptr+6); *verbose_level+=n; return 0;
    }
    if(strncasecmp(cptr, "verbose=", 8)==0) {
      if(!isdigit(cptr[8])) return 1;
      n=atoi(cptr+8); *verbose_level+=n; return 0;
    }
    if(strcasecmp(cptr, "quiet")==0) {*verbose_level=0; return 0;}
    if(strcasecmp(cptr, "silent")==0) {*verbose_level=-1; return 0;}
    return 1;
  }
  /* So it is the short form, if anything */
  if(strcmp(cptr, "h")==0) {*print_usage=1; return 0;}
  if(strcasecmp(cptr, "v")==0) {*print_version=1; return 0;}
  if(strcmp(cptr, "d")==0) {*verbose_level+=1; return 0;}
  if(strncmp(cptr, "d", 1)==0 && strlen(cptr)>1) {
    if(!isdigit(cptr[1])) return 1;
    n=atoi(cptr+1); *verbose_level+=n; return 0;
  }
  if(strcmp(cptr, "q")==0) {*verbose_level=0; return 0;}
  if(strcmp(cptr, "s")==0) {*verbose_level=-1; return 0;}
  return 1;
}
/*****************************************************************************/

/*****************************************************************************/
/** Process program name and optionally version into given string from argv[0].
 */
void tpcProgramName(
  /** Set to argv[0] */
  const char *program,
  /** Add version (1) or do not add (0) */
  int version,
  /** Add copyright (1) or do not add (0) */
  int copyright,
  /** Pointer to string where program name is written */
  char *prname,
  /** Length of prname string, including trailing zero */
  int n
) {
  char *tmp;

  /* Check the input */
  if(prname==NULL || n<1) return;
  prname[0]=(char)0;
  
  /* Remove path and extension */
  if(strlen(program)>0) tmp=strdup(program); else tmp=strdup("unknown"); 
  filenameRmPath(tmp); filenameRmExtension(tmp);
  /* Copy it if possible */
  n-=strlen(tmp); if(n>0) strcpy(prname, tmp); else {free(tmp); return;}
  free(tmp);

  /* Add version, if required */
  if(version!=0) {
    /* Create string with version number */
    char v[256];
    sprintf(v, "%d.%d.%d", tpcclib_VERSION_MAJOR,
                           tpcclib_VERSION_MINOR, tpcclib_VERSION_PATCH);
    /* Copy it if possible */
    n--; // space
    n-=strlen(v); if(n>0) {strcat(prname, " "); strcat(prname, v);} 
  }
  
  /* Add copyright, if required */
  if(copyright!=0) {
    /* Copy it if possible */
    n--; // space
    n-=strlen(tpcclib_COPYRIGHT); 
    if(n>0) {strcat(prname, " "); strcat(prname, tpcclib_COPYRIGHT);} 
  }
  
  return;
}
/*****************************************************************************/

/*****************************************************************************/
/** Print program usage given as argument, plus program name,
 *  tpcclib version, and default copyright text.
 *  - Any string \@P, separated by space characters, is replaced by program name
 *    in the output. It can only be used once per line.
 *  - When line contains text 'stdoptions', in place of that the description of
 *    standard command-line options '-h, -v, -q, -s etc' are displayed.
 *  @author Vesa Oikonen
 */
void tpcPrintUsage(
  /** Program name */
  const char *program,
  /** Program usage text */
  char *text[],
  /** File pointer where to print; usually stdout */
  FILE *fp
) {
  int i;
  char *cptr, *bprogram;

  /* Print program name, version, and copyright */
  if(strlen(program)>0) bprogram=strdup(program);
  else bprogram=strdup("unknown"); 
  filenameRmPath(bprogram); filenameRmExtension(bprogram);
  fprintf(fp, "\n  %s - tpcclib %d.%d.%d %s\n \n", bprogram, 
           tpcclib_VERSION_MAJOR, tpcclib_VERSION_MINOR, tpcclib_VERSION_PATCH,
           tpcclib_COPYRIGHT);
  /* Print usage */
  i=0; while(text[i]!=0) {
    /* If line contains string 'stdoptions' then print instead of it the
       description of standard command-line options */
    if(strstr(text[i], "stdoptions")) {
      int j=0;
      while(tpcstdoptions[j]!=0) fprintf(fp, "  %s\n", tpcstdoptions[j++]);        
      i++; continue;
    }
    /* If line does not contain program name, then just print it as it is */
    cptr=strstr(text[i], " @P ");
    if(cptr==NULL) {fprintf(fp, "  %s\n", text[i++]); continue;}
    /* Replace '@P' with program name */
    char *s; s=strdup(text[i]);
    s[strlen(text[i])-strlen(cptr)]=(char)0;
    fprintf(fp, "  %s %s %s\n", s, bprogram, cptr+4);
    free(s); i++;
  }
  fprintf(fp, " \n");
  /* Print licence info */
  i=0; while(tpclicense[i]!=0) fprintf(fp, "  %s\n", tpclicense[i++]);
  fprintf(fp, "\n");
  free(bprogram);
}
/*****************************************************************************/

/*****************************************************************************/
/** Write program usage given as argument, plus program name,
 *  tpcclib version, and default copyright text, into HTML file.
 *  - Any string \@P, separated by space characters, is replaced by program name
 *    in the output. It can only be used once per line.
 *  - When line contains text 'stdoptions', in place of that the description of
 *    standard command-line options '-h, -v, -q, -s etc' are displayed.
 *  @return Returns 0 when successful. 
 *  @author Vesa Oikonen
 */
int tpcHtmlUsage(
  /** Program name, may contain extension and path */
  const char *program,
  /** Program usage text */
  char *text[],
  /** Path name where to create file programname.html;
   *  path may contain trailing '/' or '\\'.   */
  const char *path
) {
  unsigned int len, i, j;
  char *bprogram, *fname, *cptr, *line;
  FILE *fp;

  if(program==NULL || text==NULL || strlen(program)<1) return 1;
  
  /* Clean program name */
  bprogram=strdup(program);
  filenameRmPath(bprogram); filenameRmExtension(bprogram);
  
  /* Make filename */
  fname=calloc(strlen(path)+1+strlen(bprogram)+5, sizeof(char));
  if(fname==NULL) {free(bprogram); return 1;}
  strcpy(fname, path); len=strlen(fname); 
  if(len>0 && (fname[len-1]=='/' || fname[len-1]=='\\')) fname[len-1]=(char)0;
  len=strlen(fname); if(len>0) strcat(fname, "/"); 
  strcat(fname, bprogram); strcat(fname, ".html");
  
  //printf("fname := '%s'\n", fname);

  /* Open file for write */
  fp=stdout;

  /* Write HTML header */
  len=fprintf(fp, "<!DOCTYPE html>\n");
  if(len<10) {free(bprogram); free(fname); return 2;}
  fprintf(fp, "<html>\n");
  fprintf(fp, "<head>\n");
  fprintf(fp, "  <meta charset=\"UTF-8\">\n");
  fprintf(fp, "  <title>%s</title>\n", bprogram);
  fprintf(fp, "  <style type=\"text/css\">\n");
  fprintf(fp, "    * {font-family: monospace;}\n");
  fprintf(fp, "    footer {\n");
  fprintf(fp, "      border:1px solid gray;\n");
  fprintf(fp, "      font-size: smaller;\n");
  fprintf(fp, "    }\n");
  fprintf(fp, "    img {border-width: 0px;}\n");
  fprintf(fp, "  </style>\n");
  fprintf(fp, "</head>\n\n");
  
  /* Write HTML body */
  fprintf(fp, "<body>\n");
  
  /* Write program name, version and copyright as title; */
  /* replace (c) with html code when necessary */
  fprintf(fp, "<h2>%s - tpcclib %d.%d.%d ", bprogram, tpcclib_VERSION_MAJOR,
          tpcclib_VERSION_MINOR, tpcclib_VERSION_PATCH);
  line=tpcclib_COPYRIGHT; len=strlen(line);
  for(j=0; j<len; j++) {
    if(strncasecmp(line+j, "(C)", 3)==0) {fputs("&copy;", fp); j+=2; continue;}
    fputc(line[j], fp);
  }
  fputs("</h2>\n\n", fp);

  /* Print usage */
  fprintf(fp, "<pre>\n");
  i=0; while(text[i]!=0) {
    line=text[i];
    /* If line contains string 'stdoptions' then print instead of it the
       description of standard command-line options */
    if(strstr(line, "stdoptions")) {
      int j=0;
      while(tpcstdoptions[j]!=0) fprintf(fp, "%s\n", tpcstdoptions[j++]);        
      i++; continue;
    }

    /* Process "See also" line: add links to other programs */
    if(strstr(line, "See also: ")!=NULL) {
      /* copy until the first ':' */
      j=0; while(line[j]!='\0') {
        fputc(line[j], fp);
        j++; if(line[j-1]==':') break;
      }
      /* the rest of line with token */
      char *tline; unsigned int n=0;
      tline=strdup(line+j); cptr=strtok(tline, ", :;\t\n\r");
      while(cptr!=NULL) {
        if(n>0) fputc(',', fp);
        fprintf(fp, " <a href=\"./%s.html\">%s</a>", cptr, cptr);
        cptr=strtok(NULL, ", :;\t\n\r"); 
	n++;
      }
      fputs("\n", fp);
      free(tline); 
      i++; continue;
    }

    
    /* Print the line one character at the time */
    len=strlen(line); j=0;
    while(j<len) {
      
      /* If WWW Address follows, then add a link */
      if(strncasecmp(line+j, "http://", 7)==0) {
        unsigned int li;
        cptr=line+j; len=strcspn(cptr, " ),;");
        fputs("<a href=\"", fp);
        for(li=0; li<len; li++) fputc(line[j+li], fp);
        fputs("\">", fp);
        for(li=0; li<len; li++) fputc(line[j+li], fp);
        fputs("</a>", fp);
        j+=len; continue;
      }
      
      /* If necessary, replace '@P' with program name */
      if(strncmp(line+j, " @P ", 4)==0) {
        fprintf(fp, " %s ", bprogram);
        j+=4; continue;
      }
      /* Replace (c) or (C) with html code */
      if(strncasecmp(line+j, "(C)", 3)==0) {
        fputs("&copy;", fp);
        j+=3; continue;
      }

      /* Replace <, >, and & characters with html codes */
      if(line[j]=='<') {fputs("&lt;", fp); j++; continue;}
      if(line[j]=='>') {fputs("&gt;", fp); j++; continue;}
      if(line[j]=='&') {fputs("&amp;", fp); j++; continue;}

      /* Just write the char normally */
      fputc(line[j], fp); j++;
    }
    fprintf(fp, "\n");
    i++; continue;    

  }
  fprintf(fp, "</pre>\n");

  /* Write footer */
  fprintf(fp, "\n<footer>\n");

  /* Licence info */
  fprintf(fp, "<div>\n");
  /* Icon with link to GNU-GPL */
  fprintf(fp, 
    "<a href=\"http://www.gnu.org/licenses/gpl-3.0-standalone.html\">\n");
  fprintf(fp, "<img alt=\"GNU GPL\" ");
  fprintf(fp, 
    "style=\"width:88px; height:31px; float:left; margin: 5px 20px 5px 5px;\"");
  fprintf(fp, 
    "\n src=\"http://www.turkupetcentre.net/petanalysis/pic/gplv3-88x31.png\"></a>\n");
  /* Text */
  fprintf(fp, "<p>");
  i=0; while(tpclicense[i]!=0) fprintf(fp, "%s<br>\n", tpclicense[i++]);
  fprintf(fp, "</p>\n");
  fprintf(fp, "</div>\n");
  fprintf(fp, "</footer>\n");
           
  /* Close HTML */
  fprintf(fp, "</body>\n");
  fprintf(fp, "</html>\n");


  free(bprogram); free(fname);
  return 0;
}

/*****************************************************************************/

/*****************************************************************************/
/** Print tpctools build information.
 *  @author Vesa Oikonen
 */
void tpcPrintBuild(
  /** Program name; enter NULL, if not to be printed */
  const char *program,
  /** File pointer where to print; usually stdout */
  FILE *fp
) {
  fprintf(fp, "\n");
  if(program!=NULL) {
    /* Print program name */
    char *s; s=strdup(program); 
    filenameRmPath(s); filenameRmExtension(s);
    fprintf(fp, " Program: %s\n", s);
    free(s);
  }
  /* Build time */
  fprintf(fp, " Build: %s %s\n",__DATE__,__TIME__);
  /* tpcclib (and program) version */
  fprintf(fp, " tpcclib version: %d.%d.%d\n", tpcclib_VERSION_MAJOR,
           tpcclib_VERSION_MINOR, tpcclib_VERSION_PATCH);
  /* Compiler information */
#if defined(__STDC_VERSION__)
  fprintf(fp, " Version of C: %ld\n", __STDC_VERSION__);
#endif
#if defined(__GNUC__) && defined(__VERSION__)
  fprintf(fp, " GNU C version: %s\n", __VERSION__);
#endif
#if defined(__clang__) && defined(__clang_version__)
  fprintf(fp, " Clang/LLVM version: %s\n", __clang_version__);
#endif
#ifdef _OPENMP
  fprintf(fp, " OpenMP version: %d\n", _OPENMP);
#endif
  /* Platform information */
#if defined(__x86_64__) || defined(__LP64__) || defined(__ppc64__) || \
    defined(__LLP64__) || defined(__ILP64__) 
  fprintf(fp, " Architecture: 64-bit\n");
#else
  fprintf(fp, " Architecture: 32-bit\n");
#endif
}
/*****************************************************************************/

/*****************************************************************************/
