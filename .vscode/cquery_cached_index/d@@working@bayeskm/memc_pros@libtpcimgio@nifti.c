/// @file nifti.c
/// @author Vesa Oikonen
/// @brief Procedures for reading and writing NIfTI-1 PET images.
///
/// Function are not intended to support all NIfTI files or file
/// properties, but only those that have been found necessary in
/// Turku PET Centre. For full NIfTI support, use other libraries
/// e.g. niftilib <http://niftilib.sourceforge.net/>
///
/// NIfTI-1 and NIfTI-2 documentation and source codes in 
/// <http://nifti.nimh.nih.gov/>
///
/// Procedures in this file are not dependent on IMG struct. 
///
/******************************************************************************/
#include "libtpcimgio.h"
/******************************************************************************/

/*****************************************************************************/
/** Remove any extensions from Nifti file name, leaving only base file name.
 */
void niftiRemoveFNameExtension(
  /** Full name of file. */
  char *fname
) {
  char *cptr;
  cptr=strrchr(fname, '.'); if(cptr==NULL) return;
  if(strcasecmp(cptr, ".")==0 || strcasecmp(cptr, ".img")==0 ||
     strcasecmp(cptr, ".hdr")==0 || strcasecmp(cptr, ".sif")==0 ||
     strcasecmp(cptr, ".nii")==0)
       *cptr=(char)0;
  /* Remove also double extensions, e.g. from data.img.hdr */
  cptr=strrchr(fname, '.'); if(cptr==NULL) return;
  if(strcasecmp(cptr, ".img")==0 || strcasecmp(cptr, ".nii")==0) *cptr=(char)0;
}
/*****************************************************************************/

/*****************************************************************************/
/** Construct the file names for NIfTI image.
    @return Returns 0 if successful, otherwise <>0.
 */
int niftiCreateFNames(
  /** Filename, either header file, image file, or base name without extensions,
      but possibly with path name. This string is never modified. */
  const char *filename,
  /** Header filename will be written in this char pointer
      (space needs to allocated by caller);
      If header and image are combined, then this will be the name of combined
      file; enter NULL if not needed. */
  char *hdrfile,
  /** Image filename will be written in this char pointer
      (space needs to allocated by caller);
      If header and image are combined, then this will be the name of combined
      file; enter NULL if not needed. */
  char *imgfile,
  /** SIF filename will be written in this char pointer
      (space needs to allocated by caller); enter NULL if not needed. */
  char *siffile,
  /** NIfTI file format, either IMG_NIFTI_1D (31) or IMG_NIFTI_1S (32). */
  int fileformat
) {
  int n;
  char basename[FILENAME_MAX];

  if(hdrfile!=NULL) strcpy(hdrfile, "");
  if(imgfile!=NULL) strcpy(imgfile, "");
  if(siffile!=NULL) strcpy(siffile, "");
  n=strlen(filename); if(n<1 || n>=FILENAME_MAX) {return(1);}
  strlcpy(basename, filename, FILENAME_MAX); 
  niftiRemoveFNameExtension(basename);
  
  /* Create database filenames */
  if(fileformat==IMG_NIFTI_1D) {
    if(hdrfile!=NULL) snprintf(hdrfile, FILENAME_MAX, "%s.hdr", basename);
    if(imgfile!=NULL) snprintf(imgfile, FILENAME_MAX, "%s.img", basename);
  } else if(fileformat==IMG_NIFTI_1S) {
    if(hdrfile!=NULL) snprintf(hdrfile, FILENAME_MAX, "%s.nii", basename);
    if(imgfile!=NULL) snprintf(imgfile, FILENAME_MAX, "%s.nii", basename);
  } else
    return(2);
  if(siffile!=NULL) snprintf(siffile, FILENAME_MAX, "%s.sif", basename);
  return(0);
}
/*****************************************************************************/

/*****************************************************************************/
/** Remove header and voxel data files or the single .nii file belonging to
    specified NIfTI database. 

    SIF is not deleted in any case.
    Validity of NIfTI is not verified, therefore this can be used to
    delete any files with similar name as NIfTI would have.

    @return Returns 0 when call was successful, otherwise <>0. Call is considered
    successful, if files do not exist initially.
 */
int niftiRemove(
  /** NIfTI database name with path, possibly with filename extension */
  const char *dbname,
  /** NIfTI file format, either IMG_NIFTI_1D (31) or IMG_NIFTI_1S (32),
   *  or IMG_UNKNOWN (0) in case both are to be deleted. */
  int fileformat,
  /** Verbose level; if zero, then nothing is printed to stderr or stdout */
  int verbose
) {
  int ret=0, errNr=0;
  char imgfile[FILENAME_MAX], hdrfile[FILENAME_MAX], siffile[FILENAME_MAX];

  if(verbose>0) {
    printf("niftiRemove(%s, %d, ...)\n", dbname, fileformat);
    fflush(stdout);
  }

  ret=niftiCreateFNames(dbname, hdrfile, imgfile, siffile, fileformat);
  if(ret==0 && fileformat==IMG_NIFTI_1D) { // dual format
    if(access(hdrfile, 0)!=-1) {
      if(verbose>1) {printf("  removing %s\n", hdrfile); fflush(stdout);}
      if(remove(hdrfile)!=0) errNr++;
    }
    if(access(imgfile, 0)!=-1) {
      if(verbose>1) {printf("  removing %s\n", imgfile); fflush(stdout);}
      if(remove(imgfile)!=0) errNr++;
    }
  } else if(ret==0 && fileformat==IMG_NIFTI_1S) { // single format
    if(access(imgfile, 0)!=-1) {
      if(verbose>1) {printf("  removing %s\n", imgfile); fflush(stdout);}
      if(remove(imgfile)!=0) errNr++;
    }
  } else { // dual and single formats
    ret=niftiCreateFNames(dbname, hdrfile, imgfile, siffile, IMG_NIFTI_1D);
    if(ret!=0) return 1;
    if(access(hdrfile, 0)!=-1) {
      if(verbose>1) {printf("  removing %s\n", hdrfile); fflush(stdout);}
      if(remove(hdrfile)!=0) errNr++;
    }
    if(access(imgfile, 0)!=-1) {
      if(verbose>1) {printf("  removing %s\n", imgfile); fflush(stdout);}
      if(remove(imgfile)!=0) errNr++;
    }

    ret=niftiCreateFNames(dbname, hdrfile, imgfile, siffile, IMG_NIFTI_1S);
    if(ret!=0) return 1;
    if(access(imgfile, 0)!=-1) {
      if(verbose>1) {printf("  removing %s\n", imgfile); fflush(stdout);}
      if(remove(imgfile)!=0) errNr++;
    }
  }
  return errNr;
}
/*****************************************************************************/

/*****************************************************************************/
/** Verify if specified filename is a NIfTI file.
   @return Returns 0 if it is not, 1 if header and image data are found
   (either as one combined file or as separate files,
   and 2, if sif file is found too.
 */
int niftiExists(
  /** Filename, either header file, image file, or base name without extensions.
      This string is never modified. */
  const char *filename,
  /** If filename refers to a Nifti file, then header filename will be
      written in this char pointer (space needs to allocated by caller);
      If header and image are combined, then this will be the name of combined
      file; enter NULL if not needed. */
  char *hdrfile,
  /** If filename refers to a Nifti file, then image filename will be
      written in this char pointer (space needs to allocated by caller);
      If header and image are combined, then this will be the name of combined
      file; enter NULL if not needed. */
  char *imgfile,
  /** If filename refers to a Nifti file, and if SIF exists, then SIF filename
      will be written in this char pointer (space needs to allocated by caller);
      NULL if not needed. */
  char *siffile,
  /** Pointer to Nifti header, which is filled in this function; enter NULL, if not needed. */
  NIFTI_DSR *header,
  /** Verbose level; if zero, then nothing is printed to stderr or stdout. */
  int verbose,
  /** Pointer to a string (allocated for at least 64 chars) where error message
      or other execution status will be written; enter NULL, if not needed. */
  char *status
) {
  char basefile[FILENAME_MAX], temp[FILENAME_MAX], localhdrfile[FILENAME_MAX];
  NIFTI_DSR *dsr, local_dsr;
  int ret, combined=0; 

  if(filename==NULL || strlen(filename)==0) return(0);
  if(verbose>0) {printf("\nniftiExists(%s, ...)\n", filename); fflush(stdout);}
  if(status!=NULL) strcpy(status, "OK");
  if(header==NULL) dsr=&local_dsr; else dsr=header;

  localhdrfile[0]=(char)0;

  /* Construct the base file name wo extensions */
  strlcpy(basefile, filename, FILENAME_MAX); 
  niftiRemoveFNameExtension(basefile);
  if(verbose>1) printf("\n  basefile := %s\n", basefile);

  /* Combined header and image file exists? */
  strcpy(temp, basefile); strcat(temp, ".nii");
  if(access(temp, 0) == -1) {
    if(verbose>0) printf("  %s not found or accessible.\n", temp);
  } else {
    /* Preserve header and image filenames */
    strcpy(localhdrfile, temp);
    if(hdrfile!=NULL) strcpy(hdrfile, temp);
    if(imgfile!=NULL) strcpy(imgfile, temp);
    combined=1;
    if(verbose>1) printf("  %s is accessible.\n", temp);
  }

  /* If not, then check if header file exists */
  if(combined==0) {
    strcpy(temp, basefile); strcat(temp, ".hdr");
    if(access(temp, 0) == -1) {
      strcpy(temp, basefile); strcat(temp, ".img.hdr");
      if(access(temp, 0) == -1) {
        if(verbose>0) printf("  hdr file not found or accessible.\n");
        if(status!=NULL) strcpy(status, "file not accessible");
        return(0);
      }
    }
    /* Preserve header filename */
    strcpy(localhdrfile, temp);
    if(hdrfile!=NULL) strcpy(hdrfile, temp);
    if(verbose>1) printf("  %s is accessible.\n", temp);
  }

  /* If not combined, then does image file exists? */
  if(combined==0) {
    strcpy(temp, basefile); strcat(temp, ".img");
    if(access(temp, 0) == -1) {
      if(verbose>0) printf("  %s not found or accessible.\n", temp);
      if(status!=NULL) strcpy(status, "file not accessible");
      return(0);
    }
    /* Preserve image filename */
    if(imgfile!=NULL) strcpy(imgfile, temp);
    if(verbose>1) printf("  %s is accessible.\n", temp);
  }

  /* Is this Nifti file? */
  if((ret=niftiReadHeader(localhdrfile, dsr, verbose, temp))!=0) {
    if(status!=NULL) strcpy(status, "file is not Nifti");
    if(verbose>0) {
      printf("  %s was not identified as Nifti header file (%d).\n", localhdrfile, ret);
      printf("  %s\n", temp);
    }
    return(0);
  }
  if(verbose>1) printf("  %s is identified as Nifti.\n", localhdrfile);
  if(verbose>10) niftiPrintHeader(dsr, stdout);

  /* SIF exists? */
  strcpy(temp, basefile); strcat(temp, ".sif");
  if(verbose>3) printf("  checking if %s exists\n", temp);
  if(access(temp, 0) == -1) {
    strcpy(temp, basefile); strcat(temp, ".img.sif");
    if(verbose>3) printf("  checking if %s exists\n", temp);
    if(access(temp, 0) == -1) {
      strcpy(temp, basefile); strcat(temp, ".nii.sif");
      if(verbose>3) printf("  checking if %s exists\n", temp);
      if(access(temp, 0) == -1) {
        if(verbose>0) printf("\n  SIF not found or accessible.\n");
        if(siffile!=NULL) strcpy(siffile, "");
        if(status!=NULL) {
          if(combined) strcpy(status, "combined Nifti file is accessible");
          else strcpy(status, "Nifti files are accessible");
        }
        return(1); // but otherwise ok
      }
    }
  }
  /* Preserve SIF filename */
  if(siffile!=NULL) strcpy(siffile, temp);
  if(verbose>1) printf("  %s is accessible.\n", temp);
  if(status!=NULL) {
    if(combined) strcpy(status, "combined Nifti file and SIF are accessible");
    else strcpy(status, "Nifti files and SIF are accessible");
  }
  return(2);
}
/*****************************************************************************/

/*****************************************************************************/
/** Read Nifti header contents. Currently, does not read Nifti-1 header extension.
    @return Returns 0, if successful, otherwise >0.
 */
int niftiReadHeader(
  /** Name of file to read (including path and extension) */
  char *filename,
  /** Pointer to previously allocated header structure */
  NIFTI_DSR *dsr,
  /** Verbose level; if zero, then nothing is printed to stderr or stdout */
  int verbose,
  /** Pointer to a string (allocated for at least 64 chars) where error message
      or other execution status will be written; enter NULL, if not needed */     
  char *status
) {
  FILE *fp;
  int little; // 1 if current platform is little endian (x86), else 0
  int same_order, extender=0;
  unsigned char buf[NIFTI_HEADER_SIZE];
  short int s;
  int n;

  /* Check arguments */
  if(filename==NULL || strlen(filename)==0 || dsr==NULL) return(1);
  if(verbose>0) {
    printf("\nniftiReadHeader(%s, ...)\n", filename); fflush(stdout);}
  if(status!=NULL) strcpy(status, "OK");
  little=little_endian(); if(verbose>3) printf("  little := %d\n", little);

  /* Open file */
  fp=fopen(filename, "rb"); if(fp==NULL) {
    if(status!=NULL) strcpy(status, "cannot open file");
    if(verbose>0) fprintf(stderr, "Error: cannot open file %s\n", filename);
    return(2);
  }

  /* Read Nifti header */
  if(fread(buf, NIFTI_HEADER_SIZE, 1, fp)<1) {
    if(status!=NULL) strcpy(status, "complete Nifti header not found");
    if(verbose>0)
      fprintf(stderr, "Error: invalid Nifti header file %s\n", filename);
    return(3);
  }
  /* Read nifti1 extender */
  for(n=0; n<4; n++) dsr->e.extension[n]=(char)0;
  if(fread(dsr->e.extension, 4, 1, fp)<1) {
    if(status!=NULL) strcpy(status, "complete Nifti header not found");
    if(verbose>1)
      fprintf(stdout, "Nifti header extender not found in %s\n", filename);
    extender=0;
  } else {
    extender=1;
  }
  /* Close file */
  fclose(fp);

  /* Read Nifti Magic number */
  memcpy(dsr->h.magic, buf+344, 4);
  if(strcasecmp(dsr->h.magic, "ni1")==0) {
    if(verbose>1) {printf("  separate hdr and img files.\n"); fflush(stdout);}
  } else if(strcasecmp(dsr->h.magic, "n+1")==0) {
    if(verbose>1) printf("  combined hdr and img data.\n");
  } else {
    if(status!=NULL) strcpy(status, "Nifti magic number not found");
    if(verbose>0) {
      fprintf(stderr, "Error: not a Nifti header file %s\n", filename);
      fflush(stderr);
    }
    if(verbose>2) {
      printf("magic := {%d, %d, %d, %d}\n", dsr->h.magic[0], dsr->h.magic[1],
              dsr->h.magic[2], dsr->h.magic[3]);
    } 
    return(4);
  }
  if(verbose>1) printf("  Nifti Magic number := %s\n", dsr->h.magic);
  /* Check that 4-byte header extender was found, if magic number is n+1 */
  if(strcasecmp(dsr->h.magic, "n+1")==0 && extender==0) {
    if(status!=NULL) strcpy(status, "Nifti header extender not found");
    if(verbose>0) {
      fprintf(stderr, "Error: not valid Nifti n+1 header file %s\n", filename);
      fflush(stderr);
    }
    return(5);
  }
  
  /* Determine from dim[0] if file is big or little endian */
  memcpy(&s, buf+40, 2); if(verbose>10) printf("  s := %d\n", s);
  if(s>0 && s<8) { // same order in file and current machine
    dsr->byte_order=little;
    same_order=1;
  } else {
    swabip(&s, 2); if(verbose>10) printf("  s := %d\n", s);
    if(s>0 && s<8) { // opposite order in file and in current machine
      if(little==1) dsr->byte_order=0; else dsr->byte_order=1;
      same_order=0;
    } else {
      if(status!=NULL) strcpy(status, "invalid Nifti byte order");
      if(verbose>0)
        fprintf(stderr, "Error: not a valid Nifti header file %s\n", filename);
      return(6);
    }
  }
  if(verbose>1) printf("  Nifti byte order := %d\n", dsr->byte_order);

  /* Size of header */
  memcpy(&n, buf+0, 4); if(!same_order) swawbip(&n, 4);
  if(n!=348) {
    if(status!=NULL) strcpy(status, "invalid Nifti sizeof_hdr");
    if(verbose>0)
      fprintf(stderr, "Error: not a valid Nifti header file %s\n", filename);
    return(7);
  }
  dsr->h.sizeof_hdr=n;

  /*  */
  memcpy(&dsr->h.data_type, buf+4, 10);
  memcpy(&dsr->h.db_name, buf+14, 18);
  if(!same_order) swawbip(buf+32, 4);
  memcpy(&dsr->h.extents, buf+32, 4);
  if(!same_order) swabip(buf+36, 2);
  memcpy(&dsr->h.session_error, buf+36, 2);
  memcpy(&dsr->h.regular, buf+38, 1);
  memcpy(&dsr->h.dim_info, buf+39, 1);

  /* dim */
  if(!same_order) swabip(buf+40, 16);
  memcpy(dsr->h.dim, buf+40, 16);
  /* intent parameters */
  if(!same_order) swawbip(buf+56, 4);
  memcpy(&dsr->h.intent_p1, buf+56, 4);
  if(!same_order) swawbip(buf+60, 4);
  memcpy(&dsr->h.intent_p2, buf+60, 4);
  if(!same_order) swawbip(buf+64, 4);
  memcpy(&dsr->h.intent_p3, buf+64, 4);
  if(!same_order) swabip(buf+68, 2);
  memcpy(&dsr->h.intent_code, buf+68, 2);

  /*  */
  if(!same_order) swabip(buf+70, 2);
  memcpy(&dsr->h.datatype, buf+70, 2);
  if(!same_order) swabip(buf+72, 2);
  memcpy(&dsr->h.bitpix, buf+72, 2);
  if(!same_order) swabip(buf+74, 2);
  memcpy(&dsr->h.slice_start, buf+74, 2);
  if(!same_order) swawbip(buf+76, 32);
  memcpy(dsr->h.pixdim, buf+76, 32);
  if(!same_order) swawbip(buf+108, 4);
  memcpy(&dsr->h.vox_offset, buf+108, 4);
  if(!same_order) swawbip(buf+112, 4);
  memcpy(&dsr->h.scl_slope, buf+112, 4);
  if(!same_order) swawbip(buf+116, 4);
  memcpy(&dsr->h.scl_inter, buf+116, 4);
  if(!same_order) swabip(buf+120, 2);
  memcpy(&dsr->h.slice_end, buf+120, 2);
  memcpy(&dsr->h.slice_code, buf+122, 1);
  memcpy(&dsr->h.xyzt_units, buf+123, 1);
  if(!same_order) swawbip(buf+124, 4);
  memcpy(&dsr->h.cal_max, buf+124, 4);
  if(!same_order) swawbip(buf+128, 4);
  memcpy(&dsr->h.cal_min, buf+128, 4);
  if(!same_order) swawbip(buf+132, 4);
  memcpy(&dsr->h.slice_duration,buf+132,4);
  if(!same_order) swawbip(buf+136, 4);
  memcpy(&dsr->h.toffset, buf+136, 4);
  if(!same_order) swawbip(buf+140, 4);
  memcpy(&dsr->h.glmax, buf+140, 4);
  if(!same_order) swawbip(buf+144, 4);
  memcpy(&dsr->h.glmin, buf+144, 4);

  /* study description */
  memcpy(&dsr->h.descrip, buf+148, 80);
  /* Auxiliary filename */
  memcpy(&dsr->h.aux_file, buf+228, 24);

  /* Transformation parameters */
  if(!same_order) swabip(buf+252, 2);
  memcpy(&dsr->h.qform_code, buf+252, 2);
  if(!same_order) swabip(buf+254, 2);
  memcpy(&dsr->h.sform_code, buf+254, 2);
  if(!same_order) swawbip(buf+256, 4);
  memcpy(&dsr->h.quatern_b, buf+256, 4);
  if(!same_order) swawbip(buf+260, 4);
  memcpy(&dsr->h.quatern_c, buf+260, 4);
  if(!same_order) swawbip(buf+264, 4);
  memcpy(&dsr->h.quatern_d, buf+264, 4);
  if(!same_order) swawbip(buf+268, 4);
  memcpy(&dsr->h.qoffset_x, buf+268, 4);
  if(!same_order) swawbip(buf+272, 4);
  memcpy(&dsr->h.qoffset_y, buf+272, 4);
  if(!same_order) swawbip(buf+276, 4);
  memcpy(&dsr->h.qoffset_z, buf+276, 4);
  if(!same_order) swawbip(buf+280, 16);
  memcpy(dsr->h.srow_x, buf+280, 16);
  if(!same_order) swawbip(buf+296, 16);
  memcpy(dsr->h.srow_y, buf+296, 16);
  if(!same_order) swawbip(buf+312, 16);
  memcpy(dsr->h.srow_z, buf+312, 16);

  memcpy(&dsr->h.intent_name, buf+328, 16);

  if(status!=NULL) strcpy(status, "complete Nifti header was read");
  if(verbose>0) fflush(stdout);
  return(0);
}
/*****************************************************************************/

/*****************************************************************************/
/** Print the contents of Nifti header to specified file pointer.
    @return Returns 0 if ok, 1 if invalid input.
 */
int niftiPrintHeader(
  /** Pointer to combined Nifti header. */
  NIFTI_DSR *dsr,
  /** File pointer where header information is printed. */
  FILE *fp
) {
  int i;
  char *cptr, tmp[256];

  /* Check input */
  if(fp==NULL || dsr==NULL) return(1);

  fprintf(fp, "Nifti header:\n");

  /* Byte order */
  if(dsr->byte_order==0) strcpy(tmp, "big"); else strcpy(tmp, "little");
  fprintf(fp, "byte_order := %s endian\n", tmp);
  /* sizeof_hdr */
  fprintf(fp, "sizeof_hdr := %d\n", dsr->h.sizeof_hdr);
  /* data_type */
  strncpy(tmp, dsr->h.data_type, 10); tmp[10]=(char)0;
  cptr=tmp; while(*cptr) {if(!isprint(cptr[0])) cptr[0]=' '; cptr++;}
  fprintf(fp, "data_type := %s\n", tmp);
  /* db_name */
  strncpy(tmp, dsr->h.db_name, 18); tmp[18]=(char)0;
  cptr=tmp; while(*cptr) {if(!isprint(*cptr)) *cptr=' '; cptr++;}
  fprintf(fp, "db_name := %s\n", tmp);
  /* extents */
  fprintf(fp, "extents := %d\n", dsr->h.extents);
  /* session_error */
  fprintf(fp, "session_error := %d\n", dsr->h.session_error);
  /* regular */
  fprintf(fp, "regular := %d\n", dsr->h.regular);
  /* dim_info */
  fprintf(fp, "dim_info := %d\n", dsr->h.dim_info);

  /* Data array dimensions */
  i=0; fprintf(fp, "dim := {%d", dsr->h.dim[i]);
  for(i=1; i<8; i++) fprintf(fp, ", %d", dsr->h.dim[i]);
  fprintf(fp, "}\n");

  /* Intent parameters */
  fprintf(fp, "intent_p1 := %g\n", dsr->h.intent_p1);
  fprintf(fp, "intent_p2 := %g\n", dsr->h.intent_p2);
  fprintf(fp, "intent_p3 := %g\n", dsr->h.intent_p3);
  fprintf(fp, "intent_code := %d\n", dsr->h.intent_code);

  /*  */
  fprintf(fp, "datatype := %d\n", dsr->h.datatype);
  fprintf(fp, "bitpix := %d\n", dsr->h.bitpix);
  fprintf(fp, "slice_start := %d\n", dsr->h.slice_start);
  i=0; fprintf(fp, "pixdim := {%g", dsr->h.pixdim[i]);
  for(i=1; i<8; i++) fprintf(fp, ", %g", dsr->h.pixdim[i]);
  fprintf(fp, "}\n");
  fprintf(fp, "vox_offset := %g\n", dsr->h.vox_offset);
  fprintf(fp, "scl_slope := %g\n", dsr->h.scl_slope);
  fprintf(fp, "scl_inter := %g\n", dsr->h.scl_inter);
  fprintf(fp, "slice_end := %d\n", dsr->h.slice_end);
  fprintf(fp, "slice_code := %d\n", dsr->h.slice_code);
  fprintf(fp, "xyzt_units := %d\n", dsr->h.xyzt_units);
  fprintf(fp, "cal_max := %g\n", dsr->h.cal_max);
  fprintf(fp, "cal_min := %g\n", dsr->h.cal_min);
  fprintf(fp, "slice_duration := %g\n", dsr->h.slice_duration);
  fprintf(fp, "toffset := %g\n", dsr->h.toffset);
  fprintf(fp, "glmax := %d\n", dsr->h.glmax);
  fprintf(fp, "glmin := %d\n", dsr->h.glmin);

  /* Study description */
  strncpy(tmp, dsr->h.descrip, 80); tmp[80]=(char)0;
  cptr=tmp; while(*cptr) {if(!isprint(*cptr)) *cptr=' '; cptr++;}
  fprintf(fp, "descrip := %s\n", tmp);
  strncpy(tmp, dsr->h.aux_file, 24); tmp[24]=(char)0;
  cptr=tmp; while(*cptr) {if(!isprint(*cptr)) *cptr=' '; cptr++;}
  fprintf(fp, "aux_file := %s\n", tmp);

  /* Transformation parameters */
  fprintf(fp, "qform_code := %d\n", dsr->h.qform_code);
  fprintf(fp, "sform_code := %d\n", dsr->h.sform_code);
  fprintf(fp, "quatern_b := %g\n", dsr->h.quatern_b);
  fprintf(fp, "quatern_c := %g\n", dsr->h.quatern_c);
  fprintf(fp, "quatern_d := %g\n", dsr->h.quatern_d);
  fprintf(fp, "qoffset_x := %g\n", dsr->h.qoffset_x);
  fprintf(fp, "qoffset_y := %g\n", dsr->h.qoffset_y);
  fprintf(fp, "qoffset_z := %g\n", dsr->h.qoffset_z);
  i=0; fprintf(fp, "srow_x := {%g", dsr->h.srow_x[i]);
  for(i=1; i<4; i++) fprintf(fp, ", %g", dsr->h.srow_x[i]);
  fprintf(fp, "}\n");
  i=0; fprintf(fp, "srow_y := {%g", dsr->h.srow_y[i]);
  for(i=1; i<4; i++) fprintf(fp, ", %g", dsr->h.srow_y[i]);
  fprintf(fp, "}\n");
  i=0; fprintf(fp, "srow_z := {%g", dsr->h.srow_z[i]);
  for(i=1; i<4; i++) fprintf(fp, ", %g", dsr->h.srow_z[i]);
  fprintf(fp, "}\n");

  strncpy(tmp, dsr->h.intent_name, 16); tmp[16]=(char)0;
  cptr=tmp; while(*cptr) {if(!isprint(*cptr)) *cptr=' '; cptr++;}
  fprintf(fp, "intent_name := %s\n", tmp);

  /* Nifti magic number */
  fprintf(fp, "magic := %s\n", dsr->h.magic);

  /* Nifti header extender */
  i=0; fprintf(fp, "extension := {%d", dsr->e.extension[i]);
  for(i=1; i<4; i++) fprintf(fp, ", %d", dsr->e.extension[i]);
  fprintf(fp, "}\n");

  fflush(fp);
  return(0);
}
/*****************************************************************************/

/*****************************************************************************/
/** Read Nifti image data, convert byte order if necessary,
    and scale values to floats. Reads only one frame at a time!
    @return Returns 0 if successful, >1 in case of an error, and specifically
     -1 in case that contents after the last image frame was requested.
 */
int niftiReadImagedata(
  /** File pointer to start of image data file, opened previously in binary mode. */
  FILE *fp,
  /** Pointer to previously filled Nifti header structure */
  NIFTI_DSR *dsr,
  /** Frame number to read [1..number of frames]. */
  int frame,
  /** Pointer to image float data allocated previously for dimz*dimy*dimx floats. */
  float *data,
  /** Verbose level; if zero, then nothing is printed to stderr or stdout. */
  int verbose,
  /** Pointer to a string (allocated for at least 64 chars) where error message
      or other execution status will be written; enter NULL, if not needed. */
  char *status
) {
  int dimNr, dimx, dimy, dimz=1, dimt=1, pxlNr=0;
  int i, n, little, start_pos, rawSize;
  char *mdata, *mptr;
  float *fptr, ss, si;
  short int *sptr;
  int *iptr;
  double d;


  if(verbose>0) {
    printf("niftiReadImagedata(fp, h, %d, data, %d)\n", frame, verbose);
    fflush(stdout);
  }
  /* Check the arguments */
  if(status!=NULL) sprintf(status, "invalid function input");
  if(frame<=0 || fp==NULL || dsr==NULL || data==NULL) return(1);

  /* Get the image data start location from header, in case of single file
     format */
  if(strcasecmp(dsr->h.magic, "n+1")==0) start_pos=(int)dsr->h.vox_offset;
  else start_pos=0;
  if(start_pos<0) start_pos=-start_pos;
  if(verbose>2) printf("  image_start_pos := %d\n", start_pos);
  /* edit it later to move to the correct frame */

  /* Get the image dimensions from header */
  if(status!=NULL) sprintf(status, "invalid image dimensions");
  dimNr=dsr->h.dim[0]; if(dimNr<2 || dimNr>4) return(2);
  dimx=dsr->h.dim[1];
  dimy=dsr->h.dim[2];
  if(dimNr>2) dimz=dsr->h.dim[3];
  if(dimNr>3) dimt=dsr->h.dim[4];
  if(frame>dimt) return(-1);
  pxlNr=dimx*dimy*dimz; if(pxlNr<1) return(4);

  // data_type is unused in Nifti
  /* Check that datatype is supported */
  if(verbose>1) printf("  verifying datatype\n");
  n=0;
  if(dsr->h.datatype & NIFTI_DT_RGB) n+=NIFTI_DT_RGB;
  if(dsr->h.datatype & NIFTI_DT_COMPLEX) n+=NIFTI_DT_COMPLEX;
  if(dsr->h.datatype & NIFTI_DT_BINARY) n+=NIFTI_DT_BINARY;
  if(dsr->h.datatype==NIFTI_DT_UNKNOWN) n+=512;
  if(n!=0) {
    if(verbose>0) printf("datatype error %d\n", n);
    if(status!=NULL) sprintf(status, "unsupported pixel datatype %d", dsr->h.datatype);
    return(6);
  }

  /* Allocate memory for the binary data */
  if(verbose>1) printf("  allocating memory for binary data\n");
  if(status!=NULL) sprintf(status, "invalid pixel data format");
  if(dsr->h.bitpix<8) return(5); // We don't support bit data
  rawSize=pxlNr*(dsr->h.bitpix/8); if(rawSize<1) return(6);
  if(verbose>1) printf("  pxlNr=%d  rawSize=%d\n", pxlNr, rawSize);
  if(status!=NULL) sprintf(status, "out of memory");
  mdata=(char*)malloc(rawSize); if(mdata==NULL) return(11);

  /* Seek the start of current frame data */
  if(verbose>1) printf("  seeking file position\n");
  start_pos+=(frame-1)*rawSize;
  if(verbose>2) printf("start_pos=%d\n", start_pos);
  fseek(fp, start_pos, SEEK_SET);
  if(ftell(fp)!=start_pos) {
    if(status!=NULL) sprintf(status, "could not move to start_pos %d", start_pos);
    free(mdata); return(7);
  }

  /* Read the data */
  if(verbose>1) printf("  reading binary data\n");
  mptr=mdata;
  if((n=fread(mptr, rawSize, 1, fp)) < 1) {
    if(status!=NULL) sprintf(status, "could read only %d bytes when request was %d", n, rawSize);
    free(mdata); return(8);
  }

  /* Convert byte order if necessary */
  little=little_endian(); mptr=mdata;
  if(little!=dsr->byte_order) {
    if(verbose>0) printf("byte conversion\n");
    switch(dsr->h.bitpix) {
      case 8: /* no conversion needed */ break;
      case 16: swabip(mptr, rawSize); break;
      case 32: swawbip(mptr, rawSize); break;
      case 64: swawbip(mptr, rawSize); break;
      default:
        if(verbose>5) printf("unsupported nifti bitpix := %d\n", dsr->h.bitpix);
        sprintf(status, "unsupported nifti bitpix := %d", dsr->h.bitpix);
        free(mdata); return(5);
    }
  }

  /* Get scaling factors */
  ss=dsr->h.scl_slope; if(ss==0) ss=1.0;
  si=dsr->h.scl_inter;

  /* Copy data to float pixel values */
  if(verbose>1) printf("  conversion to floating point voxel values\n");
  mptr=mdata;
  switch(dsr->h.datatype) {
    case NIFTI_DT_UNSIGNED_CHAR:
      if(dsr->h.bitpix!=8) {
        if(status!=NULL)
          sprintf(status, "invalid combination of datatype and bitpix (%d, %d)",
            dsr->h.datatype, dsr->h.bitpix);
        free(mdata); return(5);
      }
      fptr=data;
      for(i=0; i<pxlNr; i++, mptr++, fptr++) *fptr=si+ss*(float)(*mptr);
      break;
    case NIFTI_DT_UNSIGNED_SHORT:
      if(dsr->h.bitpix!=16) {
        if(status!=NULL)
          sprintf(status, "invalid combination of datatype and bitpix (%d, %d)",
            dsr->h.datatype, dsr->h.bitpix);
        free(mdata); return(5);
      }
      fptr=data;
      for(i=0; i<pxlNr; i++, mptr+=2, fptr++) {
        unsigned short int *uptr=(unsigned short int*)mptr; *fptr=si+ss*(float)(*uptr);
      }
      break;
    case NIFTI_DT_SIGNED_SHORT:
      if(dsr->h.bitpix!=16) {
        if(status!=NULL)
          sprintf(status, "invalid combination of datatype and bitpix (%d, %d)",
            dsr->h.datatype, dsr->h.bitpix);
        free(mdata); return(5);
      }
      fptr=data;
      for(i=0; i<pxlNr; i++, mptr+=2, fptr++) {
        sptr=(short int*)mptr; *fptr=si+ss*(float)(*sptr);
      }
      break;
    case NIFTI_DT_SIGNED_INT:
      if(dsr->h.bitpix!=16 && dsr->h.bitpix!=32) {
        if(status!=NULL)
          sprintf(status, "invalid combination of datatype and bitpix (%d, %d)",
            dsr->h.datatype, dsr->h.bitpix);
        free(mdata); return(5);
      }
      fptr=data;
      if(dsr->h.bitpix==16) {
        for(i=0; i<pxlNr; i++, mptr+=2, fptr++) {
          iptr=(int*)mptr; *fptr=si+ss*(float)(*iptr);
        }
      } else if(dsr->h.bitpix==32) {
        for(i=0; i<pxlNr; i++, mptr+=4, fptr++) {
          iptr=(int*)mptr; *fptr=si+ss*(float)(*iptr);
        }
      }
      break;
    case NIFTI_DT_FLOAT: // 16
      if(dsr->h.bitpix==32) {
        fptr=data; memcpy(fptr, mptr, pxlNr*4);
        for(i=0, fptr=data; i<pxlNr; i++, fptr++) *fptr*=ss;
        for(i=0, fptr=data; i<pxlNr; i++, fptr++) *fptr+=si;
      } else {
        if(status!=NULL)
          sprintf(status, "invalid combination of datatype and bitpix (%d, %d)",
            dsr->h.datatype, dsr->h.bitpix);
        free(mdata); return(5);
      }
      break;
    case NIFTI_DT_DOUBLE:
      if(dsr->h.bitpix!=64) {
        if(status!=NULL)
          sprintf(status, "invalid combination of datatype and bitpix (%d, %d)",
            dsr->h.datatype, dsr->h.bitpix);
        free(mdata); return(5);
      }
      for(i=0, fptr=data; i<pxlNr; i++, mptr+=8, fptr++) {
        memcpy(&d, mptr, 8); *fptr=si+ss*d;
      }
      break;
    default:
      if(status!=NULL)
        sprintf(status, "unsupported pixel datatype %d", dsr->h.datatype);
      free(mdata); return(5);
  }

  if(verbose>1) {printf("  data read successfully.\n"); fflush(stdout);}
  free(mdata);
  if(status!=NULL) sprintf(status, "ok");
  return 0;
}
/*****************************************************************************/

/*****************************************************************************/
/** Write NIfTI-1 header contents.

    Currently, does not write header extension.
    Header field 'byte_order' is used to determine the required byte order.
   @return Returns 0, if successful, otherwise >0.
 */
int niftiWriteHeader(
  /** Name of file to write (including path and extension). */
  char *filename,
  /** Pointer to previously allocated header structure. */
  NIFTI_DSR *dsr,
  /** Verbose level; if zero, then nothing is printed to stderr or stdout. */
  int verbose,
  /** Pointer to a string (allocated for at least 64 chars) where error message
      or other execution status will be written; enter NULL, if not needed. */
  char *status
) {
  FILE *fp;
  int little; // 1 if current platform is little endian (x86), else 0
  int same_order;
  unsigned char buf1[NIFTI_HEADER_SIZE];
  unsigned char buf2[NIFTI_HEADER_EXTENDER_SIZE];
  unsigned char *bptr;


  if(verbose>0) {
    printf("\nniftiWriteHeader(%s, ...)\n", filename); fflush(stdout);
  }

  /* Check arguments */
  if(status!=NULL) strcpy(status, "invalid function input");
  if(filename==NULL || strlen(filename)==0 || dsr==NULL) return(1);
  /* Check magic number */
  if(strcmp(dsr->h.magic, "ni1")!=0 && strcmp(dsr->h.magic, "n+1")!=0)
    return(1);

  /* Check if byte swapping is needed */
  little=little_endian(); if(verbose>3) printf("  little := %d\n", little);
  if(little==dsr->byte_order) same_order=1; else same_order=0;

  /* Make sure that buffers are all zeroes to begin with */
  memset(buf1, 0, sizeof(NIFTI_HEADER_SIZE));
  memset(buf2, 0, sizeof(NIFTI_HEADER_EXTENDER_SIZE));

  /* Copy header contents into buffer */
  if(verbose>2) printf("  setting write buffer\n");
  bptr=buf1+0;
  memcpy(bptr, &dsr->h.sizeof_hdr, 4); if(!same_order) swawbip(bptr, 4);
  bptr=buf1+4;
  memcpy(bptr, &dsr->h.data_type, 10);
  bptr=buf1+14;
  memcpy(bptr, &dsr->h.db_name, 18);
  bptr=buf1+32;
  memcpy(bptr, &dsr->h.extents, 4); if(!same_order) swawbip(bptr, 4);
  bptr=buf1+36;
  memcpy(bptr, &dsr->h.session_error, 2); if(!same_order) swabip(bptr, 2);
  bptr=buf1+38;
  memcpy(bptr, &dsr->h.regular, 1);
  bptr=buf1+39;
  memcpy(bptr, &dsr->h.dim_info, 1);

  bptr=buf1+40;
  memcpy(bptr, dsr->h.dim, 16); if(!same_order) swabip(bptr, 16);
  bptr=buf1+56;
  memcpy(bptr, &dsr->h.intent_p1, 4); if(!same_order) swawbip(bptr, 4);
  bptr=buf1+60;
  memcpy(bptr, &dsr->h.intent_p2, 4); if(!same_order) swawbip(bptr, 4);
  bptr=buf1+64;
  memcpy(bptr, &dsr->h.intent_p3, 4); if(!same_order) swawbip(bptr, 4);
  bptr=buf1+68;
  memcpy(bptr, &dsr->h.intent_code, 2); if(!same_order) swabip(bptr, 2);
  bptr=buf1+70;
  memcpy(bptr, &dsr->h.datatype, 2); if(!same_order) swabip(bptr, 2);
  bptr=buf1+72;
  memcpy(bptr, &dsr->h.bitpix, 2); if(!same_order) swabip(bptr, 2);
  bptr=buf1+74;
  memcpy(bptr, &dsr->h.slice_start, 2); if(!same_order) swabip(bptr, 2);
  bptr=buf1+76;
  memcpy(bptr, dsr->h.pixdim, 32); if(!same_order) swawbip(bptr, 32);
  bptr=buf1+108;
  memcpy(bptr, &dsr->h.vox_offset, 4); if(!same_order) swawbip(bptr, 4);
  bptr=buf1+112;
  memcpy(bptr, &dsr->h.scl_slope, 4); if(!same_order) swawbip(bptr, 4);
  bptr=buf1+116;
  memcpy(bptr, &dsr->h.scl_inter, 4); if(!same_order) swawbip(bptr, 4);
  bptr=buf1+120;
  memcpy(bptr, &dsr->h.slice_end, 2); if(!same_order) swabip(bptr, 2);
  bptr=buf1+122;
  memcpy(bptr, &dsr->h.slice_code, 1);
  bptr=buf1+123;
  memcpy(bptr, &dsr->h.xyzt_units, 1);
  bptr=buf1+124;
  memcpy(bptr, &dsr->h.cal_max, 4); if(!same_order) swawbip(bptr, 4);
  bptr=buf1+128;
  memcpy(bptr, &dsr->h.cal_min, 4); if(!same_order) swawbip(bptr, 4);
  bptr=buf1+132;
  memcpy(bptr, &dsr->h.slice_duration, 4); if(!same_order) swawbip(bptr, 4);
  bptr=buf1+136;
  memcpy(bptr, &dsr->h.toffset, 4); if(!same_order) swawbip(bptr, 4);
  bptr=buf1+140;
  memcpy(bptr, &dsr->h.glmax, 4); if(!same_order) swawbip(bptr, 4);
  bptr=buf1+144;
  memcpy(bptr, &dsr->h.glmin, 4); if(!same_order) swawbip(bptr, 4);

  bptr=buf1+148;
  memcpy(bptr, dsr->h.descrip, 80);
  bptr=buf1+228;
  memcpy(bptr, dsr->h.aux_file, 24);
  bptr=buf1+252;
  memcpy(bptr, &dsr->h.qform_code, 2); if(!same_order) swabip(bptr, 2);
  bptr=buf1+254;
  memcpy(bptr, &dsr->h.sform_code, 2); if(!same_order) swabip(bptr, 2);
  bptr=buf1+256;
  memcpy(bptr, &dsr->h.quatern_b, 4); if(!same_order) swawbip(bptr, 4);
  bptr=buf1+260;
  memcpy(bptr, &dsr->h.quatern_c, 4); if(!same_order) swawbip(bptr, 4);
  bptr=buf1+264;
  memcpy(bptr, &dsr->h.quatern_d, 4); if(!same_order) swawbip(bptr, 4);
  bptr=buf1+268;
  memcpy(bptr, &dsr->h.qoffset_x, 4); if(!same_order) swawbip(bptr, 4);
  bptr=buf1+272;
  memcpy(bptr, &dsr->h.qoffset_y, 4); if(!same_order) swawbip(bptr, 4);
  bptr=buf1+276;
  memcpy(bptr, &dsr->h.qoffset_z, 4); if(!same_order) swawbip(bptr, 4);
  bptr=buf1+280;
  memcpy(bptr, dsr->h.srow_x, 16); if(!same_order) swawbip(bptr, 16);
  bptr=buf1+296;
  memcpy(bptr, dsr->h.srow_y, 16); if(!same_order) swawbip(bptr, 16);
  bptr=buf1+312;
  memcpy(bptr, dsr->h.srow_z, 16); if(!same_order) swawbip(bptr, 16);
  bptr=buf1+328;
  memcpy(bptr, dsr->h.intent_name, 16);
  bptr=buf1+344;
  memcpy(bptr, dsr->h.magic, 4);

  /* Open header file for write; do not delete old contents, since this
     function may be called to update single format NIfTI */
  if(strcmp(dsr->h.magic, "ni1")==0) { // dual file format
    if(verbose>2) printf("  creating NIfTI header %s\n", filename);
    fp=fopen(filename, "wb");
  } else if(access(filename, 0)==-1) { // single file format, not exists
    if(verbose>2) printf("  creating NIfTI header %s\n", filename);
    fp=fopen(filename, "wb");
  } else { // single file format, exists already
    if(verbose>2) printf("  opening NIfTI header %s\n", filename);
    fp=fopen(filename, "r+b");
  }
  if(fp==NULL) {
    if(status!=NULL) strcpy(status, "cannot open Nifti header for write");
    return(2);
  }

  /* Write header */
  if(verbose>2) printf("  writing NIfTI header\n");
  if(fwrite(buf1, 1, NIFTI_HEADER_SIZE, fp) != NIFTI_HEADER_SIZE) {
    if(status!=NULL) strcpy(status, "cannot write Nifti header");
    fclose(fp); return(3);
  }

  /* Write extender, if necessary (leave the contents 0 0 0 0 for now) */
  if(verbose>2) printf("  writing NIfTI extender\n");
  if(fwrite(buf2, 1, NIFTI_HEADER_EXTENDER_SIZE, fp)!= NIFTI_HEADER_EXTENDER_SIZE) {
    if(status!=NULL) strcpy(status, "cannot write Nifti header extender");
    fclose(fp); return(3);
  }

  fclose(fp);

  if(verbose>2) {printf("  complete Nifti header was written\n"); fflush(stdout);}
  if(status!=NULL) strcpy(status, "complete Nifti header was written");
  return(0);
}
/*****************************************************************************/

/*****************************************************************************/
