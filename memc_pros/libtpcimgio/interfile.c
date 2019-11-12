/// @file interfile.c
/// @author Roman Krais, Vesa Oikonen
/// @brief Function(s) for interfile headers.
///
/******************************************************************************/
#include "libtpcimgio.h"
/******************************************************************************/

/******************************************************************************/
/*!
 * The function searches the keyword in the header and passes the value
 * belonging to that value back to the main program.
 * The name of the header (string 'headerName') and the requested keyword
 * (string 'searchWord') are passed to the function. It passes back the
 * value of the keyword (string 'returnValue') and possibly an error message
 * or warning (string 'errorMessage'). So the values are passed back as strings.
 * The interpretation (conversion to integer, float, time etc) is up
 * to the programmer.
 *
 * The interfile header has to comply to the following rules:
 * - first line in the file is '!INTERFILE'
 * - maximal length of a line is 512 characters
 * - A line has two fields sperated by ':=' (keyword := value)
 * - maximal length of keyword and value is 256 characters.
 * - no header entries after a line  '!END OF INTERFILE'
 * - a line starting with a semicolon ';' is a comment
 *
 * @param headerName header file name
 * @param searchWord keyword to look for
 * @param returnValue value for keyword in header
 * @param errorMessage error message/warnings. In case there is a error message 
 *        it will be returnd as string in the
 * variable 'errmsg'.
 * @return 0 if ok, 1  keyword appears more than once in the interfile header
 * (value of last occurence of keyword is returned), 
 * 2 keyword not found in interfile header (returned value is empty 
 *   (i.e. contains '/0's only)),
 * 3  interfile header cold not be opened for reading (returned value is empty 
 *    (i.e. contains '/0's only)),
 * 4  wrong file format?! (No '!INTERFILE' in the first line) (returned value 
 *    is empty (i.e. contains '/0's only))
 */
int interfile_read(char headerName[256], char searchWord[256], 
                   char returnValue[256], char errorMessage[300]
) {
  short int  i, pos;
  short int  count=0;    /* counter: How often appears keyword in the header? */
  int        n;
  char       *c[1];
  char       keyword[256], value[256];
  char       line[512];  /* max length of a line accepted in interfile header */
  FILE       *interfileHeader;

                                                        /* initialise strings */
  for (i=0;i<256;i++) returnValue[i] = '\0';
  for (i=0;i<300;i++) errorMessage[i] = '\0';

  /* open interfile header for reading */
  if ((interfileHeader = fopen(headerName,"r"))==NULL) {
    strcpy(errorMessage,headerName);
    strcat(errorMessage," could not be opened for reading");
    return 3;
  }

  /* check from first line if file is really interfile header */
  n=fread(&c,1,1,interfileHeader); if(n<1) {
    strcpy(errorMessage,"wrong file header format?! No '!INTERFILE' at start of ");
    strcat(errorMessage,headerName);
    fclose(interfileHeader);
    return 4;
  }
  i=0;
  memcpy(&line[i],c,1);
  while (memcmp(c,"\n",1) && memcmp(c,"\r",1)) {
    i++;
    n=fread(&c,1,1,interfileHeader); if(n<1) {
      strcpy(errorMessage,"wrong file header format?! No '!INTERFILE' at start of ");
      strcat(errorMessage,headerName);
      fclose(interfileHeader);
      return 4;
    }
    memcpy(&line[i],c,1);
  }
  if (memcmp(line,"!INTERFILE",10)) {
    strcpy(errorMessage,"wrong file header format?! No '!INTERFILE' at start of ");
    strcat(errorMessage,headerName);
    fclose(interfileHeader);
    return 4;
  }

 /* read file line by line */
 while (fread(&c,1,1,interfileHeader) == 1) {
    for (i=0;i<512;i++) line[i] = '\0';     /* initialise line */
    for (i=0;i<256;i++) keyword[i] = '\0';  /* initialise keyword */
    for (i=0;i<256;i++) value[i] = '\0';    /* initialise value */
    i=0;
    /* \n = end of line, \r = carriage return. Lines in  ASCII files */
    /* on Sun-Solaris end with \n, on Intel-Windows with \r\n        */
    while (memcmp(c,"\r",1) && memcmp(c,"\n",1) && i<512) {
      memcpy(&line[i],c,1);
      n=fread(&c,1,1,interfileHeader); if(n<1) {
        strcpy(errorMessage,"wrong file header format: ");
        strcat(errorMessage,headerName);
        fclose(interfileHeader);
        return 4;
      }
      i++;
    }
    /* comments are not processed */
    if (strncmp(&line[0],";",1)) {
      /* get keyword and value from line */
      /* find position of the field seperator ':=' */
      for (pos=1; pos<512; pos++)
        if (line[pos] == '=' && line[pos-1] == ':') break; 
      /* now get the first and the second field */
      for (i=0;i<pos-2 && i<256;i++) keyword[i] = line[i];
      for (i=pos+2;i<256+pos+2 && i<512;i++) {
        if (!memcmp(&line[i],"\0",1) || !memcmp(&line[i],"\r",1) || !memcmp(&line[i],"\n",1)) 
          break;      /* stop at the end of "line" */
        value[i-pos-2] = line[i];
      }
      if (!memcmp(keyword,"!END OF INTERFILE",17)) break;     /* are we done? */
      /* check if we found the keyword */
      else if (!strcmp(keyword,searchWord)) {
        strcpy(returnValue,value);
        count++;
      }
    }
  }
  fclose(interfileHeader);   /* done with reading */
  if (count == 0) {
    strcpy(errorMessage,"keyword '");
    strcat(errorMessage,searchWord);
    strcat(errorMessage,"' not found in header");
    return 2;
  }
  if (count > 1) {
    strcpy(errorMessage,"keyword '");
    strcat(errorMessage,searchWord);
    strcat(errorMessage,"' appears more than once in header");
    return 1;
  }
  return 0;
}
/******************************************************************************/

/******************************************************************************/
/** Verify that given file is a valid Interfile header file.
\return Returns 0 if not, and 1 if it is a valid header file. If image data
    filename is requested and found in header, then 2 is returned. 
 */
int interfileIsHeader(
  /** Interfile header filename, with correct extension */ 
  const char *hdrfile,
  /** Pointer to allocated string where Interfile image data filename is
   *  written if found in header; enter NULL if not needed. */
  char *imgfile
) {
  char /*key[256], value[256],*/ temp[FILENAME_MAX];
  IFT ift;
  int li, ret;

  //printf("\ninterfileIsHeader(%s, *imgfile)\n", hdrfile);
  if(hdrfile==NULL || strlen(hdrfile)<1) return 0;
  if(imgfile!=NULL) strcpy(imgfile, "");
  iftInit(&ift); strcpy(temp, hdrfile);
  ret=iftRead(&ift, temp, 0);
  if(ret!=0 || ift.keyNr<2) {iftEmpty(&ift); return(0);}
  /* Check that file starts with !INTERFILE */
  //iftWriteItem(&ift, 0, stdout);
  if(ift.item[0].type!='!') {iftEmpty(&ift); return(0);}
  if(strcasecmp(ift.item[0].value, "INTERFILE")!=0) {iftEmpty(&ift); return(0);}
  /* If imgfile was not requested, then we are done */
  if(imgfile==NULL) {iftEmpty(&ift); return(1);}
  /* Find filename for image data */
  //iftWriteItem(&ift, 1, stdout);
  li=iftGetFrom(&ift, 1, "name of data file");
  if(li<0 || strlen(ift.item[li].value)<1) {iftEmpty(&ift); return(1);}
  strcpy(imgfile, ift.item[li].value);
  iftEmpty(&ift);
  return(2);
}
/******************************************************************************/

/******************************************************************************/
/** Check if specified image filename is a Interfile data.
\return Returns 0 if it is not, 1 if it is, and both image and header is found.
 */
int interfileExists(
  /** Filename, either header file, image file, or base name without extensions.
    */
  const char *fname,
  /** If fname is Interfile, then header filename will be
     written in this char pointer (space needs to allocated by caller);
     NULL if not needed. */
  char *hdrfile,
  /** If fname is Interfile, then image filename will be
     written in this char pointer (space needs to allocated by caller);
     NULL if not needed. Note that filename is stored without path. */
  char *imgfile,
  /** Verbose level; if zero, then nothing is printed into stdout or stderr */
  int verbose
) {
  char *cptr, basefile[FILENAME_MAX], temp[FILENAME_MAX], temp2[FILENAME_MAX];

  if(fname==NULL || strlen(fname)==0) return(0);
  if(verbose>0)
    printf("\ninterfileExists(%s, *str, *str, %d)\n", fname, verbose);

  /* Construct the base file name wo extensions */
  strcpy(basefile, fname);
  cptr=strrchr(basefile, '.');
  if(cptr!=NULL) {
    if(strcasecmp(cptr, ".HDR")==0 || strcasecmp(cptr, ".I")==0 )
      *cptr=(char)0;
  } 
  cptr=strrchr(basefile, '.');
  if(cptr!=NULL) {
    if(strcasecmp(cptr, ".I")==0)
      *cptr=(char)0;
  }
  if(verbose>1) printf("\n  basefile := %s\n", basefile);

  /* Header file exists? */
  strcpy(temp, basefile); strcat(temp, ".i.hdr");
  if(access(temp, 0) == -1) {
    strcpy(temp, basefile); strcat(temp, ".hdr");
    if(access(temp, 0) == -1) {
      if(verbose>0) printf("\n  hdr file not found or accessible.\n");
      return(0);
    }
  }
  /* Is this Interfile header file? */
  if(interfileIsHeader(temp, temp2) < 2) {
    if(verbose>0)
      printf("\n  %s was not identified as Interfile header file.\n", temp);
    return(0);
  }
  /* Preserve header filename */
  if(hdrfile!=NULL) strcpy(hdrfile, temp);

  /* Image file exists? Add path from hdr file */
  cptr=strrchr(temp, '/'); if(cptr==NULL) cptr=strrchr(temp, '\\');
  if(cptr!=NULL) {
    cptr++; *cptr=(char)0; strcat(temp, temp2); strcpy(temp2, temp);} 
  if(strlen(temp2)<1 || access(temp2, 0) == -1) {
    if(verbose>0) printf("\n  %s not found or accessible.\n", temp2);
    return(0);
  }
  /* Preserve image filename */
  if(imgfile!=NULL) strcpy(imgfile, temp2);

  return 1;
}
/******************************************************************************/

/******************************************************************************/
