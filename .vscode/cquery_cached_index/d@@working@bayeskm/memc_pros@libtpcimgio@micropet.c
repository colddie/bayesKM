/// @file micropet.c
/// @author Vesa Oikonen
/// @brief Procedures for reading Siemens Inveon images.
///
/******************************************************************************/
#include "libtpcimgio.h"
/******************************************************************************/

/******************************************************************************/
/** Read specified parameter value from Concorde/MicroPET header.
   @return Returns 0 if parameter was found, even if value is empty, and 1,
    if parameter was not found, and <0 in case of other errors.
 */
int upetHeaderReadParameter(
  /** File pointer to Concorde/MicroPET header; parameter is read starting
      from file pointer forward, therefore rewind file pointer before calling
      this routine if you want to search parameter from beginning. */ 
  FILE *fp,
  /** Pointer to string which contains the header parameter name. */
  char *parameter,
  /** Pointer to allocated string where parameter value will be written;
      memory for at least MAX_MICROPET_LINE_LEN chars must be allocated;
      NULL if not needed. */
  char *value
) {
  char *cptr, tmp[MAX_MICROPET_LINE_LEN];

  if(fp==NULL) return -1;
  if(parameter==NULL || strlen(parameter)<1) return -2;
  do {
    if(fgets(tmp, MAX_MICROPET_LINE_LEN-1, fp)==NULL) return 1;
    if(tmp[0]=='#') continue;
    if(strncasecmp(tmp, parameter, strlen(parameter))!=0) continue;
    /* Get parameter value, if one exists */
    cptr=tmp+strlen(parameter)+1;
    if(strlen(cptr)>0) {
      if(value!=0) strnCopyClean(value, cptr, MAX_MICROPET_LINE_LEN-1);
    } else {
      if(value!=0) strcpy(value, "");
    }
    /* In any case, we found the parameter */
    if(MICROPET_TEST>9) printf("%s := %s\n", parameter, value);
    return 0;
  } while(1);
  return 0;
}
/******************************************************************************/

/******************************************************************************/
/** Verify that given file is a valid Concorde/microPET file header file.
    @return Returns 0 if not, and 1 if it is a valid header file. 
 */
int upetIsHeader(
  /** Concorde/microPET file header filename, with correct extension */ 
  char *hdrfile
) {
  char tmp[MAX_MICROPET_LINE_LEN];
  FILE *fp;
  int ret;

  if(hdrfile==NULL || strlen(hdrfile)<5) return 0;
  if(MICROPET_TEST>1) printf("\nupetIsHeader(%s)\n", hdrfile);
  if((fp=fopen(hdrfile, "r"))==NULL) return 0;
  /* Check that first line starts with '#' */
  if(fgets(tmp, MAX_MICROPET_LINE_LEN-1, fp)==NULL) {fclose(fp); return 0;}
  if(tmp[0]!='#') {fclose(fp); return 0;}
  /* Check that certain header parameters do exist */
  ret=upetHeaderReadParameter(fp, "version", tmp);
  if(ret!=0) {fclose(fp); return 0;}
  ret=upetHeaderReadParameter(fp, "model", tmp);
  if(ret!=0) {fclose(fp); return 0;}
  ret=upetHeaderReadParameter(fp, "institution", tmp);
  if(ret!=0) {fclose(fp); return 0;}
  fclose(fp);
  return 1;
}
/******************************************************************************/

/******************************************************************************/
/** Check if specified image filename is a Concorde/microPET file
    @return Returns 0 if it is not, and >0 if it is:
     1 if header but not image is found, or 2 if both image and header are found.
 */
int upetExists(
  /** Filename, either header file, image file, or base name without extensions. */
  const char *upetname,
  /** If upetname is a Concorde/microPET file, then header filename will be
     written in this char pointer (space needs to allocated by caller);
     NULL if not needed. */
  char *hdrfile,
  /** If upetname is a Concorde/microPET file, then image filename will be
     written in this char pointer (space needs to allocated by caller);
     NULL if not needed. */
  char *imgfile,
  /** Verbose level; if <=0, then nothing is printed into stdout */
  int verbose
) {
  char *cptr, basefile[FILENAME_MAX], temp[FILENAME_MAX];

  if(upetname==NULL || strlen(upetname)==0) return(0);
  if(verbose>0) printf("\nupetExists(%s, *str, *str)\n", upetname);

  /* Construct the base file name wo extensions */
  strlcpy(basefile, upetname, FILENAME_MAX);
  cptr=strrchr(basefile, '.');
  if(cptr!=NULL) {
    if(strncasecmp(cptr, ".HDR", 4)==0 || strncasecmp(cptr, ".IMG", 4)==0 )
      *cptr=(char)0;
  } 
  cptr=strrchr(basefile, '.');
  if(cptr!=NULL) {
    if(strncasecmp(cptr, ".IMG", 4)==0 )
      *cptr=(char)0;
  }
  if(verbose>1) printf("\n  basefile := %s\n", basefile);

  /* Header file exists? */
  strlcpy(temp, basefile, FILENAME_MAX-4); strcat(temp, ".hdr");
  if(access(temp, 0) == -1) {
    strlcpy(temp, basefile, FILENAME_MAX-8); strcat(temp, ".img.hdr");
    if(access(temp, 0) == -1) {
      if(verbose>0) printf("\n  hdr file not found or accessible.\n");
      return(0);
    }
  }
  /* Is this microPET header file? */
  if(upetIsHeader(temp)==0) {
    if(verbose>0)
      printf("\n  %s was not identified as microPET header file.\n", temp);
    return(0);
  }
  /* Preserve header filename */
  if(hdrfile!=NULL) strlcpy(hdrfile, temp, FILENAME_MAX);

  /* Image file exists? */
  strlcpy(temp, basefile, FILENAME_MAX); strcat(temp, ".img");
  if(access(temp, 0) == -1) {
    if(verbose>0) printf("\n  %s not found or accessible.\n", temp);
    return(1);
  }
  /* Preserve image filename */
  if(imgfile!=NULL) strlcpy(imgfile, temp, FILENAME_MAX);

  return(2);
}
/******************************************************************************/

/******************************************************************************/
/** Read image dimensions from header.
    @return Returns 0 when successful.
 */
int upetGetImageDimensions(
  /** File pointer to MicroPET image header */
  FILE *fp,
  /** Pointers to dimensions: planes */
  int *z,
  /** Pointers to dimensions: columns */
  int *x,
  /** Pointers to dimensions: rows */
  int *y,
  /** Pointers to dimensions: frames; if not existent (CT), enter NULL */
  int *f
) {
  char tmp[MAX_MICROPET_LINE_LEN];

  if(fp==NULL) return 1;
  *z=*x=*y=0; if(f!=NULL) *f=0;
  rewind(fp);

  if(f!=NULL) {
    if(upetHeaderReadParameter(fp, "time_frames", tmp)!=0) {
      if(upetHeaderReadParameter(fp, "total_frames", tmp)!=0) return 11;
    }  
    *f=-1; (void)sscanf(tmp, "%d", f);
  }
  if(upetHeaderReadParameter(fp, "x_dimension", tmp)!=0) return 12;
  *x=-1; (void)sscanf(tmp, "%d", x);
  if(upetHeaderReadParameter(fp, "y_dimension", tmp)!=0) return 13;
  *y=-1; (void)sscanf(tmp, "%d", y);
  if(upetHeaderReadParameter(fp, "z_dimension", tmp)!=0) return 14;
  *z=-1; (void)sscanf(tmp, "%d", z);
  if(*z<1 || *x<1 || *y<1) return 2;
  if(f!=NULL && *f<1) return 2;
  return 0;
}
/******************************************************************************/

/******************************************************************************/
/** Read scan start time from Concorde/MicroPET header.
    @return Returns 0 when successful.
 */
int upetScanStart(
  /** File pointer to Concorde/MicroPET header */
  FILE *fp,
  /** Pointer to time_t where time and date will be saved */
  time_t *scant
) {
  char tmp[MAX_MICROPET_LINE_LEN], tmp2[64], tmp3[64];
  int n, i;
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wmissing-field-initializers"
  struct tm scanstart={0};
#pragma clang diagnostic pop

  if(fp==NULL || scant==NULL) return 1;

  rewind(fp);
  if(upetHeaderReadParameter(fp, "scan_time", tmp)!=0) return 2;
  n=sscanf(tmp, "%s %s %d %d:%d:%d %d", tmp2, tmp3, &scanstart.tm_mday,
           &scanstart.tm_hour, &scanstart.tm_min, &scanstart.tm_sec, &i);
  if(n==7) {
    scanstart.tm_year=i-1900;
    if(strcasecmp(tmp3, "Jan")==0)      scanstart.tm_mon=0;
    else if(strcasecmp(tmp3, "Feb")==0) scanstart.tm_mon=1;
    else if(strcasecmp(tmp3, "Mar")==0) scanstart.tm_mon=2;
    else if(strcasecmp(tmp3, "Apr")==0) scanstart.tm_mon=3;
    else if(strcasecmp(tmp3, "May")==0) scanstart.tm_mon=4;
    else if(strcasecmp(tmp3, "Jun")==0) scanstart.tm_mon=5;
    else if(strcasecmp(tmp3, "Jul")==0) scanstart.tm_mon=6;
    else if(strcasecmp(tmp3, "Aug")==0) scanstart.tm_mon=7;
    else if(strcasecmp(tmp3, "Sep")==0) scanstart.tm_mon=8;
    else if(strcasecmp(tmp3, "Oct")==0) scanstart.tm_mon=9;
    else if(strcasecmp(tmp3, "Nov")==0) scanstart.tm_mon=10;
    else if(strcasecmp(tmp3, "Dec")==0) scanstart.tm_mon=11;
    scanstart.tm_isdst=-1;
    //*scant=mktime(&scanstart); if(*scant<0) return 4;
    *scant=timegm(&scanstart); if(*scant<0) return 4;
  } else return 5;

  return 0;
}
/******************************************************************************/

/******************************************************************************/
/** Reads microPET image data, scaling values to floats if necessary.

    Reads only one frame at a time!

    @return 0 when successful, otherwise <>0.
 */
int upetReadImagedata(
  /** file opened previously in binary mode */
  FILE *fp,
  /** microPET header in IFT struct */
  IFT *ift,
  /** frame number to read [1..number of frames] */
  int frame,
  /** pointer to image float data allocated previously */
  float *data
) {
  int dimx, dimy, dimz=1, dimt=1, pxlNr=0;
  int i, fi, n, data_type, data_bytes, little, start_pos, rawSize;
  char *mdata, *mptr;
  float f, cf, bf, *fptr;
  short int *sptr;
  int *iptr;
  char key[MAX_MICROPET_LINE_LEN], value[MAX_MICROPET_LINE_LEN];


  if(MICROPET_TEST) printf("upetReadImagedata(fp, ift, %d, data)\n", frame);

  /* Check the arguments */
  if(frame<=0 || fp==NULL || ift==NULL || data==NULL) return(1);

  /* Get the image dimensions from header */
  i=iftGetIntValue(ift, 0, "time_frames", &dimt);
  if(i<0 || dimt<1) {
    i=iftGetIntValue(ift, 0, "total_frames", &dimt);
    if(i<0 || dimt<1) return 4;
  }
  i=iftGetIntValue(ift, 0, "x_dimension", &dimx);
  if(i<0 || dimx<1) return 4;
  i=iftGetIntValue(ift, 0, "y_dimension", &dimy);
  if(i<0 || dimy<1) return 4;
  i=iftGetIntValue(ift, 0, "z_dimension", &dimz);
  if(i<0 || dimz<1) return 4;
  pxlNr=dimx*dimy*dimz; if(pxlNr<1) return(4);
  if(frame>dimt) return(3); // do not change this return value

  /* Get the data type (little=Intel, big=Sun) */
  i=iftGetIntValue(ift, 0, "data_type", &data_type);
  if(i<0 || data_type<1) return 5;
  switch(data_type) {
    case 1: data_bytes=1; little=little_endian(); break;
    case 2: data_bytes=2; little=1; break;
    case 3: data_bytes=4; little=1; break;
    case 4: data_bytes=4; little=1; break;
    case 5: data_bytes=4; little=0; break;
    case 6: data_bytes=2; little=0; break;
    case 7: data_bytes=4; little=0; break;
    default: return 5;
  }
  rawSize=data_bytes*pxlNr;
  if(MICROPET_TEST>0) {
    printf("  data_type=%d\n  data_bytes=%d\n", data_type, data_bytes);
    printf("  pxlNr=%d\n  rawSize=%d\n", pxlNr, rawSize);
  }

  /* Allocate memory for the binary data */
  mdata=(char*)malloc(rawSize); if(mdata==NULL) return(11);

  /* Seek the start of current frame data */
  start_pos=(frame-1)*rawSize;
  if(MICROPET_TEST>2) printf("  start_pos=%d\n", start_pos);
  fseek(fp, start_pos, SEEK_SET);
  if(ftell(fp)!=start_pos) {
    if(MICROPET_TEST>5) printf("could not move to start_pos\n");
    free(mdata); return(7);
  }

  /* Read the data */
  mptr=mdata;
  if((n=fread(mptr, rawSize, 1, fp)) < 1) {
    if(MICROPET_TEST>5)
      printf("could read only %d bytes when request was %d\n", n, rawSize);
    free(mdata); return(8);
  }

  /* Convert byte order if necessary */
  mptr=mdata;
  if(little!=little_endian()) {
    if(MICROPET_TEST>0) {printf("byte conversion\n"); fflush(stdout);}
    switch(data_bytes) {
      case 1: /* no conversion needed */ break;
      case 2: swabip(mptr, rawSize); break;
      case 4: swawbip(mptr, rawSize); break;
      default: 
        if(MICROPET_TEST>5)
          printf("unsupported data_type := %d\n", data_type);
        free(mdata); return(5);
    }
  }

  /* Get scale factor for this frame */
  strcpy(key, "frame"); sprintf(value, "%d", frame-1);
  fi=iftGetFullmatchFrom(ift, 0, key, value); if(fi<0) {free(mdata); return(6);}
  i=iftGetFloatValue(ift, fi+1, "scale_factor", &f);
  if(i<0 || f<=0.0) {free(mdata); return(6);}
  if(MICROPET_TEST>5) printf("  scale_factor := %g\n", f);

  /* calibration factor */
  i=iftGetFloatValue(ift, 0, "calibration_factor", &cf);
  if(i<0 || cf<0.0) {free(mdata); return 7;}
  if(MICROPET_TEST>5) printf("  calibration_factor := %g\n", cf);
  /* branching_fraction */
  i=iftGetFloatValue(ift, 0, "isotope_branching_fraction", &bf);
  if(i<0 || bf<0.0) {free(mdata); return 7;}
  if(MICROPET_TEST>5) printf("  branching_fraction := %g\n", bf);
  if(cf>0.0) {f=cf; if(bf>0.0) f/=bf;}
  if(MICROPET_TEST>5) printf("  f := %g\n", f);
 
  /* Copy data to float pixel values */
  mptr=mdata; fptr=data;
  switch(data_type) {
    case 1: // unsigned char
      for(i=0; i<pxlNr; i++, mptr++, fptr++) *fptr=f*(float)(*mptr);
      break;
    case 2: // short int
    case 6:
      for(i=0; i<pxlNr; i++, mptr+=2, fptr++) {
        sptr=(short int*)mptr; *fptr=f*(float)(*sptr);
      }
      break;
    case 3: // int
    case 7:
      for(i=0; i<pxlNr; i++, mptr+=4, fptr++) {
        iptr=(int*)mptr; *fptr=f*(float)(*iptr);
      }
      break;
    case 4: // float
    case 5:
      memcpy(fptr, mptr, pxlNr*data_bytes);
      for(i=0; i<pxlNr; i++, fptr++) *fptr*=f;
      break;
    default:
      if(MICROPET_TEST>5)
        printf("unsupported data_type := %d\n", data_type);
      free(mdata); return(5);
  }

  free(mdata);
  if(MICROPET_TEST>1) printf("anaReadImagedata() succeeded\n");
  return(0);
}
/******************************************************************************/

/******************************************************************************/
