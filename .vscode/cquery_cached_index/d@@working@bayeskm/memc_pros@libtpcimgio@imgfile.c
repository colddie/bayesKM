/// @file imgfile.c
/// @author Vesa Oikonen, Harri Merisaari, Calle Laakkonen
/// @brief I/O routines for IMG data.
///
/// Currently supported file formats:
///   - ECAT 6.3 images and sinograms
///   - ECAT 7.x 2D and 3D images (volumes) and sinograms
///   - Analyze 7.5 images (subset)
///   - NIfTI-1 images (subset)
///   - microPET images (only reading)
///
/*****************************************************************************/
#include "libtpcimgio.h"
/*****************************************************************************/

/*****************************************************************************/
/*!
   Read an image or sinogram file in ECAT 6.3 or ECAT 7.x format,
   or image in NIfTI-1, Analyze 7.5, or microPET format.
  
   @sa imgReadFrame, imgWrite
   @return 0 if ok, 1 invalid input, 2 image status is not 'initialized', 
   4 unrecognised format, 5 unsupported Ecat7 type,
   sets IMG->statmsg in case of error.
 */
int imgRead(
  /** Input filename. */
  const char *fname, 
  /** Pointer to initialized IMG structure.
      @sa imgInit, imgEmpty
   */
  IMG *img
) {
  FILE *fp;
  int ret;
  ECAT7_mainheader ecat7_main_header;
  ECAT63_mainheader ecat63_main_header;
  char temp[FILENAME_MAX];
  //char *cptr;

  if(IMG_TEST) {printf("imgRead(%s, *img)\n", fname); fflush(stdout);}
  /* Check the arguments */
  if(fname==NULL) {img->statmsg=imgStatus(STATUS_FAULT); return(1);}
  if(img==NULL || img->status!=IMG_STATUS_INITIALIZED) {
    img->statmsg=imgStatus(STATUS_FAULT); return(2);}

  /* Check if we have NIfTI file, which may be in single file format,
     or dual file format which has similar names as microPET and Analyze */
  if(niftiExists(fname, NULL, NULL, NULL, NULL, IMG_TEST-3, NULL)>0) {
    /* Read NIfTI image */
    ret=imgReadNifti(fname, img, IMG_TEST);
    if(IMG_TEST) {printf("imgReadNifti() := %d\n", ret); fflush(stdout);}
    if(ret==STATUS_OK) {
      if(IMG_TEST) printf("%s identified as supported NIfTI.\n", fname);
      img->statmsg=imgStatus(STATUS_OK);
      return(STATUS_OK);
    }
    img->statmsg=imgStatus(ret); return(4);
  }

  /* Check if we have microPET or Analyze file, which consist of separate header
     and data files, and have similar names */
  if(upetExists(fname, NULL, NULL, IMG_TEST-3)==2) {
    /* Read microPET image */
    ret=imgReadMicropet(fname, img); if(ret!=STATUS_OK) return(3);
    if(IMG_TEST) printf("%s identified as microPET format.\n", fname);
    return(0);
  }
  if(anaExistsNew(fname, temp, NULL, NULL)!=0) {
    anaRemoveFNameExtension(temp);
    /* Read Analyze image */
    ret=imgReadAnalyze(temp, img);
    if(IMG_TEST) {printf("imgReadAnalyze() := %d\n", ret); fflush(stdout);}
    if(ret==STATUS_OK) {
      if(IMG_TEST) printf("%s identified as supported Analyze 7.5 format.\n",
        fname);
      img->statmsg=imgStatus(STATUS_OK);
      return(0);
    }
    if(ret==STATUS_NOSIFDATA || ret==STATUS_WRONGSIFDATA) {
      img->statmsg=imgStatus(ret); return(0);}
    img->statmsg=imgStatus(ret); return(4);
  }

  /* Check if we have an ECAT file */
  /* Open file for read */
  if((fp=fopen(fname, "rb")) == NULL) {
    img->statmsg=imgStatus(STATUS_NOFILE); return(4);
  }
  /* Try to read ECAT 7.x main header */
  ret=ecat7ReadMainheader(fp, &ecat7_main_header);
  if(ret) {fclose(fp); img->statmsg=imgStatus(STATUS_UNKNOWNFORMAT); return(4);}
  /* If header could be read, check for magic number */
  if(strncmp(ecat7_main_header.magic_number, ECAT7V_MAGICNR, 7)==0) {
    /* This is ECAT 7.x file */
    /* Check if file type is supported */
    if(imgEcat7Supported(&ecat7_main_header)==0) {
      fclose(fp); img->statmsg=imgStatus(STATUS_UNSUPPORTED); return(5);
    }
    fclose(fp);
    /* Read file */
    if(IMG_TEST) printf("%s identified as supported ECAT 7.x %s format\n",
      fname, ecat7filetype(ecat7_main_header.file_type));
    ret=imgReadEcat7(fname, img);
    if(ret) {if(IMG_TEST) printf("imgReadEcat7()=%d\n", ret); return(6);}
  } else {
    /* Check if file is in ECAT 6.3 format */
    ret=ecat63ReadMainheader(fp, &ecat63_main_header);
    fclose(fp);
    if(ret==0) {
      /* It seems to be ECAT 6.3, so read it */
      if(IMG_TEST) printf("%s identified as supported ECAT 6.3 %s format\n",
        fname, ecat7filetype(ecat63_main_header.file_type));
      ret=ecat63ReadAllToImg(fname, img);
      if(ret) {
        if(IMG_TEST) fprintf(stderr, "ecat63ReaddAllToImg: %s\n", ecat63errmsg);
        if(ret==6) img->statmsg=imgStatus(STATUS_MISSINGMATRIX);
        else img->statmsg=imgStatus(STATUS_UNSUPPORTED);
        return(6);
      }
    } else {img->statmsg=imgStatus(STATUS_UNKNOWNFORMAT); return(4);}
  }
  img->statmsg=imgStatus(STATUS_OK);
  return(0);
}
/*****************************************************************************/

/*****************************************************************************/
/*!
   Write an image or sinogram file.

   Format depends on _fileFormat or file name extension.
  
   @sa imgInit, imgEmpty, imgRead
   @return 0 if ok, 1 invalid input, 2 invalid image type or status, 
   5 failed to write file, sets IMG->statmsg in case of error
 */
int imgWrite(
  /** File name for output image. */
  const char *fname,
  /** Pointer to IMG data. */
  IMG *img
) {
  int ret;

  if(IMG_TEST) printf("imgWrite(%s, *img)\n", fname);
  /* Check the arguments */
  if(fname==NULL) return(1);
  if(img==NULL || img->status!=IMG_STATUS_OCCUPIED) {
    imgSetStatus(img, STATUS_FAULT); return(2);}
  if(img->type!=IMG_TYPE_RAW &&
     img->type!=IMG_TYPE_IMAGE &&
     img->type!=IMG_TYPE_POLARMAP) {
    imgSetStatus(img, STATUS_FAULT); return(2);}

  /* If _fileFormat is not defined, then determine it from the file name */
  if(img->_fileFormat==IMG_UNKNOWN) {
    if(IMG_TEST>1) printf("  file format determined based on file name\n");
    imgFormatFromFName(img, fname);
    if(IMG_TEST>1) printf("  _fileFormat := %d\n", img->_fileFormat);
  }

  /* Write */
  if(img->_fileFormat==IMG_E63) {
    ret=ecat63WriteAllImg(fname, img);
    switch(ret) {
      case 0: break;
      case 4: imgSetStatus(img, STATUS_NOMEMORY); break;
      case 3: imgSetStatus(img, STATUS_NOWRITEPERM); break;
      case 9: imgSetStatus(img, STATUS_DISKFULL); break;
      default: imgSetStatus(img, STATUS_FAULT);
    }
    if(ret) return(7);
  } else if(img->_fileFormat==IMG_ANA || img->_fileFormat==IMG_ANA_L) {
    ret=imgWriteAnalyze(fname, img); if(ret) return(5);
  } else if(img->_fileFormat==IMG_NIFTI_1S || img->_fileFormat==IMG_NIFTI_1D) {
    ret=imgWriteNifti(fname, img, 1, IMG_TEST-1); if(ret) return(5);
  } else if(img->_fileFormat==IMG_E7_2D) {
    ret=imgWrite2DEcat7(fname, img); if(ret) return(5);
  } else if(img->_fileFormat==IMG_POLARMAP) {
    ret=imgWritePolarmap(fname, img); if(ret) return(5);
  } else {
    ret=imgWriteEcat7(fname, img); if(ret) return(5);
  }
  imgSetStatus(img, STATUS_OK);
  return(0);
}
/*****************************************************************************/

/*****************************************************************************/
/*!
   Fill IMG struct header information from an image or sinogram file
   in ECAT 6.3, ECAT 7.x or Analyze 7.5 format.

   Information concerning separate frames or planes is not filled.
 
   @sa imgInit, imgEmpty, imgReadFrame
   @return errstatus, which is STATUS_OK (0) when call was successful,
   and >0 in case of an error.
 */
int imgReadHeader(
  /** Image or sinogram file name */
  const char *fname,
  /** Pointer to initialized but not allocated IMG structure. */
  IMG *img,
  /** Image file format, if known (IMG_E63, IMG_E7, ...); if not known, set to
      IMG_UNKNOWN or 0. */
  int format
) {
  int ret;

  if(IMG_TEST) {
    printf("\nimgReadHeader(%s, *img, %d)\n", fname, format);
    fflush(stdout);
  }

  /* Check the arguments */
  if(fname==NULL) return STATUS_FAULT;
  if(img==NULL) return STATUS_FAULT;
  if(img->status!=IMG_STATUS_INITIALIZED) return STATUS_FAULT;

  /* If user did not know the file format, then try to decide it */
  if(format==IMG_UNKNOWN) {
    int scanner, imgtype, modality;
    ret=imgFormatDetermine(fname, NULL, NULL, NULL, NULL, 
           &format, &scanner, &imgtype, &modality, IMG_TEST-3);
    if(ret!=0) {
      imgSetStatus(img, ret);
      return(ret);
    }
    if(format==IMG_UNKNOWN) {
      imgSetStatus(img, STATUS_UNSUPPORTED);
      return(STATUS_UNSUPPORTED);
    }
  }

  /* Read the header for the file format */
  ret=STATUS_UNSUPPORTED;
  if(format==IMG_ANA || format==IMG_ANA_L) {
    /* Read Analyze header information */
    ret=imgReadAnalyzeHeader(fname, img);
  } else if(format==IMG_NIFTI_1S || format==IMG_NIFTI_1D) {
    /* Read NIfTI image header information */
    ret=imgReadNiftiHeader(fname, img, IMG_TEST-2);
  } else if(format==IMG_MICROPET) {
    /* Read microPET header information */
    ret=imgReadMicropetHeader(fname, img);
  } else if(format==IMG_E7 || format==IMG_E7_2D) {
    ret=imgReadEcat7Header(fname, img);
  } else if(format==IMG_E63 || format==IMG_POLARMAP) {
    ret=imgReadEcat63Header(fname, img);
  }  
  imgSetStatus(img, ret);
  return(ret);
}
/*****************************************************************************/

/*****************************************************************************/
/*!
   Read one time frame from a supported PET image or sinogram file into
   IMG data structure.

   This functions can be called repeatedly to read all
   the frames one at a time to conserve memory.
  
   @sa imgRead, imgReadHeader, imgWriteFrame
   @return errstatus, which is STATUS_OK (0) when call was successful,
   and >0 in case of an error.  Specifically, return value STATUS_NOMATRIX
   signals that frame does not exist, i.e. all frames have been read.
   IMG.statmsg can be set using ERROR_STATUS.
 */
int imgReadFrame(
  /** Name of file from which IMG contents will be read.
      Currently supported file formats are ECAT 6.3 images and sinograms,
      ECAT 7.x 2D and 3D images and sinograms, and NIfTI-1, Analyze 7.5, and
      microPET 3D and 4D images. */
  const char *fname,
  /** Frame which will be read [1..frameNr]. */
  int frame_to_read,
  /** Pointer to initiated or occupied IMG data.
      If occupied, then new frame is tested to match the previous file type,
      dimensions, and other fundamental information contained in the IMG.
      If not occupied, then memory is allocated here. */
  IMG *img,
  /** IMG frame index (0..dimt-1) where data will be placed.
      If index is >0, then the memory for that frame must be allocated before
      calling this function. */
  int frame_index
) {
  IMG test_img;
  int ret=0;

  if(IMG_TEST) {
    printf("\nimgReadFrame(%s, %d, *img, %d)\n",
           fname, frame_to_read, frame_index);
    fflush(stdout);
  }
  /*
   *  Check the input 
   */
  if(fname==NULL) return STATUS_FAULT;
  if(img==NULL) return STATUS_FAULT;
  if(img->status!=IMG_STATUS_INITIALIZED && img->status!=IMG_STATUS_OCCUPIED)
    return STATUS_FAULT;
  if(frame_to_read<1) return STATUS_FAULT;
  if(frame_index<0) return STATUS_FAULT;
  /* if frame_index>0, then there must be sufficient memory allocated for it */
  if(frame_index>0) {
    if(img->status!=IMG_STATUS_OCCUPIED) return STATUS_FAULT;
    if(frame_index>img->dimt-1) return STATUS_FAULT;
  }

  /*
   *  If IMG is preallocated, check that fundamental header information
   *  is compatible with old and new contents.
   *  If not allocated, then read the header contents, and allocate it
   */
  imgInit(&test_img);
  if(img->status==IMG_STATUS_OCCUPIED) {
    ret=imgReadHeader(fname, &test_img, img->_fileFormat); 
    imgSetStatus(&test_img, ret);
    if(IMG_TEST>1)
      printf("imgReadHeader() return message := %s\n", test_img.statmsg);
    if(ret) return(ret);
    if(IMG_TEST>3) imgInfo(&test_img);
    /* Test that file format and type are the same */
    ret=0;
    if(img->type!=test_img.type) ret++;
    if(img->_fileFormat!=test_img._fileFormat) ret++;
    /* Test that x, y, and z dimensions are the same */
    if(img->dimx!=test_img.dimx) ret++;
    if(img->dimy!=test_img.dimy) ret++;
    if(img->dimz!=test_img.dimz) ret++; 
    imgEmpty(&test_img); if(ret>0) return STATUS_INVALIDHEADER;
  } else {
    ret=imgReadHeader(fname, img, IMG_UNKNOWN);
    imgSetStatus(img, ret);
    if(IMG_TEST>1)
      printf("imgReadHeader() return message := %s\n", img->statmsg);
    if(ret) return(ret);
    if(IMG_TEST>3) imgInfo(img);
    /* Allocate memory for one frame */
    img->dimt=1;
    ret=imgAllocate(img, img->dimz, img->dimy, img->dimx, img->dimt);
    if(ret) return STATUS_NOMEMORY;
  }

  /*
   *  Read the frame data and corresponding information like frame time
   *  if available   
   */
  switch(img->_fileFormat) {
    case IMG_E7:
    case IMG_E7_2D:
    case IMG_POLARMAP:
      ret=imgReadEcat7Frame(fname, frame_to_read, img, frame_index);
      if(IMG_TEST>1) printf("imgReadEcat7Frame() return value := %d\n", ret);
      break;
    case IMG_E63:
      ret=imgReadEcat63Frame(fname, frame_to_read, img, frame_index);
      if(IMG_TEST>1) printf("imgReadEcat63Frame() return value := %d\n", ret);
      break;
    case IMG_ANA:
    case IMG_ANA_L:
      ret=imgReadAnalyzeFrame(fname, frame_to_read, img, frame_index);
      if(IMG_TEST>1) printf("imgReadAnalyzeFrame() return value := %d\n", ret);
      break;
    case IMG_NIFTI_1S:
    case IMG_NIFTI_1D:
      ret=imgReadNiftiFrame(fname, frame_to_read, img, frame_index, 0);
      if(IMG_TEST>1) printf("imgReadNiftiFrame() return value := %d\n", ret);
      break;
    case IMG_MICROPET:
      ret=imgReadMicropetFrame(fname, frame_to_read, img, frame_index);
      if(IMG_TEST>1) printf("imgReadAnalyzeFrame() return value := %d\n", ret);
      break;
    default:
      ret=STATUS_UNSUPPORTED;
  }
  imgSetStatus(img, ret);
  return ret;
}
/*****************************************************************************/

/*****************************************************************************/
/** Write one PET frame from IMG data struct into a supported PET image or
    sinogram file. 

    This function can be called repeatedly to write all
    frames one at a time to conserve memory.

    Currently supported file formats are ECAT 6.3 images and sinograms, and
    ECAT 7.x 2D and 3D images and sinograms, and NIfTI-1 images.
    Analyze 7.5 images are NOT supported (because global min and max would
    be needed). SIF for NIfTI-1 files is not written by this function.  

    @sa imgReadFrame, imgReadHeader, imgWrite
    @return Returns errstatus, which is STATUS_OK (0) when call was successful,
    and >0 in case of an error.
 */
int imgWriteFrame(
  /** Name of file where IMG contents will be written.
      If file exists, data is either overwritten or catenated as a new frame,
      depending on the following arguments.
      If file does not exist, it is created. */
  const char *fname,
  /** PET frame number (1..frameNr) which will be written:
      If set to 0, frame data will be written to an existing or new PET file as
      a new frame, never overwriting existing data.
      If >0, then frame data is written as specified frame number, overwriting
      any data existing with the same frame number */
  int frame_to_write,
  /** Pointer to the IMG data struct */
  IMG *img,
  /** IMG frame index (0..dimt-1) which will be written */
  int frame_index
) {
  int ret=0;

  if(IMG_TEST>0) {
    printf("\nimgWriteFrame(%s, %d, *img, %d)\n",
           fname, frame_to_write, frame_index);
    fflush(stdout);
  }
  if(IMG_TEST>3) {
    char buf[32];
    if(!ctime_r_int(&img->scanStart, buf)) strcpy(buf, "1900-01-01 00:00:00");
    fprintf(stdout, "  scan_start_time := %s\n", buf);
  }
    
  /*
   *  Check the input 
   */
  if(fname==NULL) return STATUS_FAULT;
  if(img==NULL) return STATUS_FAULT;
  if(img->status!=IMG_STATUS_OCCUPIED) return STATUS_FAULT;
  if(frame_to_write<0) return STATUS_FAULT;
  if(frame_index<0 || frame_index>=img->dimt) return STATUS_FAULT;


  /*
   *  Call separate function for each supported file format
   */
  imgFormatFromFName(img, fname);
  switch(img->_fileFormat) {
    case IMG_E7:
    case IMG_E7_2D:
    case IMG_POLARMAP:
      ret=imgWriteEcat7Frame(fname, frame_to_write, img, frame_index);
      break;
    case IMG_E63:
      ret=imgWriteEcat63Frame(fname, frame_to_write, img, frame_index);
      break;
    case IMG_ANA:
    case IMG_ANA_L:
      ret=STATUS_UNSUPPORTED;
      /* Not supported because would require global min&max values
       * if saved in short ints which is now the only possibility
       * ret=imgWriteAnaFrame(fname, frame_to_write, img, frame_index);
       */
      break;
    case IMG_NIFTI_1D:
    case IMG_NIFTI_1S:
#if(0)
      /* Not supported because would require global min&max values
       * if saved in short ints which is now the only possibility
       */
      ret=STATUS_UNSUPPORTED;
#endif
      /* Nifti is currently always written as floats, therefore
         global min and max pixel values are not yet needed */
      ret=imgWriteNiftiFrame(fname, frame_to_write, img, frame_index,
                             0, 0, IMG_TEST-2);
      break;
    default:
      ret=STATUS_UNSUPPORTED;
  }
  imgSetStatus(img, ret);
  return ret;
}
/*****************************************************************************/

/*****************************************************************************/
/**
   Determine IMG _fileFormat from file name extension, if not already defined.
   Default is ECAT 7 image volume, if nothing else can be guessed.

   Note that NIfTI dual file format, Analyze, and microPET files can not be
   separated by file naming, thus all of these formats will be identified as
   Analyze by this function (default may be changed to NIfTI in future).
  
   @sa imgRead, imgWrite
 */
void imgFormatFromFName(
  /** Target image structure where file format is saved;
      should have IMG_UNKNOWN as file type. */
  IMG *img,
  /** Name of file that is used to determine format. */
  const char *fname
) {
  char *cptr=NULL, *cptr2=NULL, temp[FILENAME_MAX];

  if(IMG_TEST>2) printf("imgFormatFromFName(img, %s)\n", fname);
  if(img->_fileFormat!=IMG_UNKNOWN && img->_fileFormat>0) {
    if(IMG_TEST>3)
      printf("  _fileFormat := %d, not changed\n", img->_fileFormat);
    return;
  }
  img->_fileFormat=IMG_E7; /* default */
  /* get extensions */
  strcpy(temp, fname); cptr=strrchr(temp, '.');
  if(cptr!=NULL) {
    *cptr=(char)0; cptr++;
    cptr2=strrchr(temp, '.'); if(cptr2!=NULL) {*cptr2=(char)0; cptr2++;}
  }
  if(cptr2!=NULL) {
    if(strcasecmp(cptr2, "i.hdr")==0) { img->_fileFormat=IMG_INTERFILE; return;}
    if(strcasecmp(cptr2, "i.img")==0) { img->_fileFormat=IMG_INTERFILE; return;}
  }
  if(cptr!=NULL) {
    if(strcasecmp(cptr, "hdr")==0) { img->_fileFormat=IMG_ANA; return;}
    if(strcasecmp(cptr, "polmap")==0) { img->_fileFormat=IMG_POLARMAP; return;}
    if(strcasecmp(cptr, "img")==0 ||
       strcasecmp(cptr, "scn")==0 ||
       strcasecmp(cptr, "nrm")==0 ||
       strcasecmp(cptr, "atn")==0) {
      img->_fileFormat=IMG_E63; return;
    }
    if(strcasecmp(cptr, "dcm")==0) { img->_fileFormat=IMG_DICOM; return;}
    if(strcasecmp(cptr, "i")==0) { img->_fileFormat=IMG_INTERFILE; return;}
    if(strcasecmp(cptr, "nii")==0) { img->_fileFormat=IMG_NIFTI_1S; return;}
  } else { /* no extension at all */
    img->_fileFormat=IMG_ANA;
  }
}
/*****************************************************************************/

/*****************************************************************************/
/** Determine the file format and type of an existing PET image data.

    Note that this function only identifies those formats and types that were
    needed in TPC at the time of writing this function.
    Note also that this function does not care whether the format is fully or
    not at all supported by other library functions.
    
    @sa imgRead, imgWrite, imgReadHeader
    @return Returns 0 when no errors occurred, which does not mean that format
    was identified.
 */
int imgFormatDetermine(
  /** Full name of image file. In case of Analyze or microPET image, hdr or img
   *  file name, or basename without extension is accepted. Pathname is not
   *  accepted. */
  const char *fname,
  /** Pointer to allocated string where image file name without extensions
   *  is written (only for Analyze and microPET); enter NULL, if not needed. */ 
  char *basename,
  /** Pointer to allocated string where header file name without extensions
   *  is written (only for Analyze and microPET); enter NULL, if not needed. */ 
  char *hdrfile,
  /** Pointer to allocated string where image data file name without extensions
   *  is written (only for Analyze and microPET); enter NULL, if not needed. */ 
  char *imgfile,
  /** Pointer to allocated string where SIF file name without extensions
   *  is written (only for Analyze), if available; enter NULL, if not needed. */ 
  char *siffile,
  /** Pointer to int where image format ID is written; 0=unknown, for other
   *  formats see definitions in img.h */
  int *file_format,
  /** Pointer to int where scanner ID is written; 0=unknown, for other
   *  formats see definitions in img.h */
  int *scanner,
  /** Pointer to int where image type is written; 0=unknown, 1=image,
   *  2=sinogram, 3=polarmap, please check definitions in img.h */
  int *type,
  /** Pointer to int where modality is written; 0=unknown, 1=PET,
   *  2=MRI, 3=CT, ..., please check definitions in img.h */
  int *modality,
  /** Verbose level; if zero, then nothing is printed into stdout or stderr */
  int verbose
) {
  char *cptr, temp[FILENAME_MAX];
  int ret, fformat=IMG_UNKNOWN;
  NIFTI_DSR nifti_dsr;
  IMG img;

  if(verbose>0) {
    printf("imgFormatDetermine(\"%s\", ...)\n", fname);
    fflush(stdout);
  }
  /* Initiate results */
  if(basename!=NULL) strcpy(basename, "");
  if(hdrfile!=NULL) strcpy(hdrfile, "");
  if(imgfile!=NULL) strcpy(imgfile, "");
  if(siffile!=NULL) strcpy(siffile, "");
  if(file_format!=NULL) *file_format=IMG_UNKNOWN;
  if(scanner!=NULL) *scanner=SCANNER_UNKNOWN;
  if(type!=NULL) *type=IMG_TYPE_UNKNOWN;
  if(modality!=NULL) *modality=IMG_MODALITY_UNKNOWN;
  if(fname==NULL) return STATUS_FAULT;
  if(strlen(fname)<1) return STATUS_NOFILE;

  /* Check the image data exists and is accessible */
  strcpy(temp, fname);
  if(access(temp, 0) == -1) {
    if(verbose>1) printf("  file is not directly accessible.\n");
    /* Try to add .nii to the name */
    sprintf(temp, "%s.nii", fname);
    if(access(temp, 0) == -1) {
      if(verbose>1)
        printf("  file is not accessible with .nii extension.\n");
      /* Try to add .img to the name */
      sprintf(temp, "%s.img", fname);
      if(access(temp, 0) == -1) {
        if(verbose>1)
          printf("  file is not accessible with .img extension.\n");
        /* Try to add .hdr to the name */
        sprintf(temp, "%s.hdr", fname);
        if(access(temp, 0) == -1) sprintf(temp, "%s.i.hdr", fname);
        if(access(temp, 0) == -1) sprintf(temp, "%s.img.hdr", fname);
        if(access(temp, 0) == -1) {
          if(verbose>1)
            printf("  file is not accessible with .hdr extension.\n");
          /* Try to add .dcm to the name */
          sprintf(temp, "%s.dcm", fname);
          if(access(temp, 0) == -1) {
            if(verbose>1)
              printf("  file is not accessible with .dcm extension.\n");
            return STATUS_NOFILE;
          }
        }
      }
    }
  }
  if(verbose>1) {printf("'%s' is accessible.\n", temp); fflush(stdout);}

  /* DICOM is identified from the file name extension, and not processed any further */
  cptr=strrchr(temp, '.');
  if(cptr!=NULL && strcasecmp(cptr, ".DCM")==0) {
    fformat=IMG_DICOM; if(file_format!=NULL) *file_format=fformat;
    if(verbose>1) printf("file was identified to be in DICOM format.\n");
    return STATUS_OK;
  }

  /* Try to read it as ECAT file first, because images which consist of 
     more than one file may reside in the same folder with
     other formats */
  imgInit(&img);
  /* Is this an ECAT7 file */
  ret=imgReadEcat7Header(fname, &img);
  if(ret==STATUS_OK) {
    fformat=img._fileFormat; if(file_format!=NULL) *file_format=fformat;
    if(verbose>1) printf("file was identified to be in ECAT7 format.\n");
  } else if(ret==STATUS_VARMATSIZE) {
    fformat=img._fileFormat; if(file_format!=NULL) *file_format=fformat;
    if(verbose>1) printf("file is ECAT7 but matrix sizes are different\n.");
  } else if(ret==STATUS_UNKNOWNFORMAT || ret==STATUS_NOFILE) {
    /* If main header was read but format was not identified as Ecat7, 
       it might be in Ecat6 format */
    ret=imgReadEcat63Header(fname, &img);
    /* if necessary, try also with .img added */
    if(ret!=STATUS_OK && ret!=STATUS_VARMATSIZE && ret!=STATUS_MISSINGMATRIX) {
      sprintf(temp, "%s.img", fname); ret=imgReadEcat63Header(temp, &img);}
    if(ret==STATUS_OK) {
      /* Is this an ECAT6 file; however this is rather uncertain step, because
         ECAT6 files don't contain any magic number */
      fformat=img._fileFormat; if(file_format!=NULL) *file_format=fformat;
      if(verbose>1) printf("file was identified to be in ECAT6 format.\n");
    } else if(ret==STATUS_VARMATSIZE || ret==STATUS_MISSINGMATRIX) {
      fformat=img._fileFormat; if(file_format!=NULL) *file_format=fformat;
      if(verbose>1) printf("file is ECAT63 but matrix sizes are different\n.");
    }
  }

  /* If format was not yet identified, then try to read it as NIfTI;
     this must be done before trying Analyze, because dual format NIfTI
     is compatible with Analyze format */
  if(fformat==IMG_UNKNOWN &&
     niftiExists(fname, hdrfile, imgfile, siffile, &nifti_dsr,
                 verbose-2, NULL)>0) {
    if(verbose>1) printf("file was identified to be in NIfTI format.\n");
    /* Read NIfTI header information into IMG struct */
    ret=imgGetNiftiHeader(&img, &nifti_dsr, verbose-2);
    if(ret==STATUS_OK) {
      fformat=img._fileFormat; if(file_format!=NULL) *file_format=fformat;
    }
  }

  /* If format was not yet identified, then try to read it as microPET */
  if(fformat==IMG_UNKNOWN && upetExists(fname, hdrfile, imgfile, verbose-1)==2) {
    fformat=IMG_MICROPET; if(file_format!=NULL) *file_format=fformat;
    if(verbose>1) printf("file was identified to be in microPET format\n.");
    /* Read header information from file */
    ret=imgReadMicropetHeader(fname, &img);
    if(ret==STATUS_OK) {
      fformat=img._fileFormat; if(file_format!=NULL) *file_format=fformat;
    }
  }

  /* If format was not yet identified, then try to read it as Analyze */
  if(fformat==IMG_UNKNOWN &&
     anaDatabaseExists(fname, hdrfile, imgfile, siffile)>0) {
    if(verbose>1) printf("file was identified to be in Analyze format.\n");
    fformat=IMG_ANA; if(file_format!=NULL) *file_format=fformat;
    /* Read Analyze header information */
    ret=imgReadAnalyzeHeader(hdrfile, &img);
    if(ret==STATUS_OK) {
      fformat=img._fileFormat; if(file_format!=NULL) *file_format=fformat;
    }
  }

  /* If format was not yet identified, check if it is DICOM (without .dcm extension
     that was checked previously) */
  if(dcmVerifyMagic(fname, NULL)) {
    fformat=IMG_DICOM; if(file_format!=NULL) *file_format=fformat;
    if(verbose>1) printf("file was identified to be in DICOM format.\n");
    return STATUS_OK;
  }

  /* If format was not yet identified, check if it is Interfile */
  if(fformat==IMG_UNKNOWN &&
     interfileExists(fname, hdrfile, imgfile, verbose-1)!=0) {
    fformat=IMG_INTERFILE; if(file_format!=NULL) *file_format=fformat;
    if(verbose>1) printf("file was identified to be in Interfile format\n.");
  }

  /* Fill other information */
  if(scanner!=NULL) *scanner=img.scanner;
  if(type!=NULL) *type=img.type;
  if(modality!=NULL) *modality=img.modality;

  /* Quit */
  imgEmpty(&img);
  if(verbose>1) printf("fformat := %d\n", fformat);
  return STATUS_OK;
}
/*****************************************************************************/

/*****************************************************************************/
