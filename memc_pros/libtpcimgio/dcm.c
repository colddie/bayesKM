/** @file dcm.c
    @brief IO functions for DICOM files.
 */
/******************************************************************************/
#include "libtpcimgio.h"
/******************************************************************************/

/******************************************************************************/
/*! @cond PRIVATE */
/** One item for table of DICOM value representations (VRs). */
typedef struct DCM_VR {
  /** VR id */
  dcmvr vr;
  /** Pointer to the NULL terminated VR string. */
  char name[8];
  /** Nr of reserved bytes following VR; 0 or 2. */
  unsigned char res;
  /** Endian sensitive or not. */
  short int es;
  /** Corresponding max value length in bytes; 0 if not defined. */
  size_t s;
  /** Description of the VR; NULL terminated string. */
  char descr[64];
} DCM_VR;

/** Table of DICOM value representations (VRs). 

    Items must be the same and in the same order as the dcmdr enums in tpcdcm.h.
    Table can not be accessed directly outside the c file.
*/
static DCM_VR dcm_vr[]={
  {DCM_VR_AE, "AE", 0, 0,    16, "application entity"},
  {DCM_VR_AS, "AS", 0, 0,     4, "age string"},
  {DCM_VR_AT, "AT", 0, 0,     4, "attribute tag"},
  {DCM_VR_CS, "CS", 0, 0,    16, "code string"},
  {DCM_VR_DA, "DA", 0, 0,     8, "date"},
  {DCM_VR_DS, "DS", 0, 0,    16, "decimal string"}, // actual files may contain longer strings too
  {DCM_VR_DT, "DT", 0, 0,    26, "date and time"},
  {DCM_VR_FL, "FL", 0, 1,     4, "floating point single precision"},
  {DCM_VR_FD, "FD", 0, 1,     8, "floating point double precision"},
  {DCM_VR_IS, "IS", 0, 0,    12, "integer string"},
  {DCM_VR_LO, "LO", 0, 0,    64, "long string"},
  {DCM_VR_LT, "LT", 0, 0, 10240, "long text"},
  {DCM_VR_OB, "OB", 2, 0,     0, "other byte (8-bit) stream"}, // even bytes
  {DCM_VR_OD, "OD", 2, 1,     0, "other double (64-bit) stream"},
  {DCM_VR_OF, "OF", 2, 1,     0, "other float (32-bit) stream"},
  {DCM_VR_OL, "OL", 2, 1,     0, "other long (32-bit) stream"},
  {DCM_VR_OW, "OW", 2, 1,     0, "other word (16-bit) stream"},
  {DCM_VR_PN, "PN", 0, 0,    64, "person name"},
  {DCM_VR_SH, "SH", 0, 0,    16, "short string"},
  {DCM_VR_SL, "SL", 0, 1,     4, "signed long (32-bit integer)"},
  {DCM_VR_SQ, "SQ", 2, 0,     0, "sequence of elements (used for nested data)"},
  {DCM_VR_SS, "SS", 0, 1,     2, "signed short (16-bit integer)"},
  {DCM_VR_ST, "ST", 0, 0,  1024, "short text"},
  {DCM_VR_TM, "TM", 0, 0,    14, "time"},
  {DCM_VR_UC, "UC", 2, 0,     0, "unlimited characters"},
  {DCM_VR_UI, "UI", 0, 0,    64, "UID"},
  {DCM_VR_UL, "UL", 0, 1,     4, "unsigned long (32-bit integer)"},
  {DCM_VR_UN, "UN", 2, 0,     0, "unknown, any valid length of another VR"},
  {DCM_VR_UR, "UR", 2, 0,    64, "URI or URL string"},
  {DCM_VR_US, "US", 0, 1,     2, "unsigned short (16-bit integer)"},
  {DCM_VR_UT, "UT", 2, 0,     0, "unlimited text"},
  // This MUST be kept as the last element
  {DCM_VR_INVALID, "INVALID", 0, 0, 0, "invalid value representation"}
};
/*! @endcond */
/******************************************************************************/

/******************************************************************************/
/*! @cond PRIVATE */
/** One item for table of DICOM Transfer Syntax UIDs. */
typedef struct DCM_TRUID_ITEM {
  /** Enumerated Id */
  dcmtruid id;
  /** Pointer to the NULL terminated UID string. */
  char uid[64];
  /** Description of the VR; NULL terminated string. */
  char descr[64];
} DCM_TRUID_ITEM;

/** Table of DICOM Transfer Syntax UIDs. 

    Items must be the same and in the same order as the dcmtruid enums 
    in tpcdcm.h. Table can not be accessed directly outside the c file.
*/
static DCM_TRUID_ITEM dcm_truid[]={
  {DCM_TRUID_UNKNOWN, "1.2.840.10008.1.2", "unknown"}, //< default
  {DCM_TRUID_LEI, "1.2.840.10008.1.2", "implicit VR little endian"},
  {DCM_TRUID_LEE, "1.2.840.10008.1.2.1", "explicit VR little endian"},
  {DCM_TRUID_BEE, "1.2.840.10008.1.2.2", "explicit VR big endian"},
  {DCM_TRUID_JPEG50, "1.2.840.10008.1.2.4.50", "lossy JPEG 8-bit compression"},
  {DCM_TRUID_JPEG51, "1.2.840.10008.1.2.4.51", "lossy JPEG 12-bit compression"},
  {DCM_TRUID_JPEG70, "1.2.840.10008.1.2.4.70", "lossless JPEG"},
  {DCM_TRUID_JPEG80, "1.2.840.10008.1.2.4.80", "lossless JPEG-LS"},
  {DCM_TRUID_JPEG81, "1.2.840.10008.1.2.4.81", "lossy JPEG-LS"},
  {DCM_TRUID_JPEG90, "1.2.840.10008.1.2.4.90", "lossless JPEG 2000"},
  {DCM_TRUID_JPEG91, "1.2.840.10008.1.2.4.91", "JPEG 2000"},
  {DCM_TRUID_JPEG92, "1.2.840.10008.1.2.4.92", "lossless multicomponent JPEG 2000"},
  {DCM_TRUID_JPEG93, "1.2.840.10008.1.2.4.93", "multicomponent JPEG 2000"},
  {DCM_TRUID_MPEG100, "1.2.840.10008.1.2.4.100", "MPEG-2"},
  {DCM_TRUID_MPEG102, "1.2.840.10008.1.2.4.102", "MPEG-4"},
  {DCM_TRUID_MPEG103, "1.2.840.10008.1.2.4.103", "MPEG-4 BD-compatible"},
  {DCM_TRUID_RLE, "1.2.840.10008.1.2.5", "lossless RLE"},
  {DCM_TRUID_RFC, "1.2.840.10008.1.2.6.1", "RFC 2557"},
  {DCM_TRUID_XML, "1.2.840.10008.1.2.6.2", "XML"},
  {DCM_TRUID_INVALID, "", "invalid"}
};
/*! @endcond */
/******************************************************************************/

/******************************************************************************/
/*! @cond PRIVATE */
/** One item for table of DICOM SOPs. */
typedef struct DCM_SOP_ITEM {
  /** Pointer to the NULL terminated SOP UID string. 
      See DICOM tag (0008,0016) for SOP Class UID. */
  char uid[64];
  /** SOP Name; NULL terminated string. */
  char name[64];
} DCM_SOP_ITEM;

/** Table of DICOM Storage SOPs. 
    @note Table can not be accessed directly outside the c file.
*/
static DCM_SOP_ITEM dcm_sop[]={
  {"invalid", "invalid SOP"}, // do not change this
  {"1.2.840.10008.5.1.4.1.1.1", "Computed Radiography Image Storage"},
  {"1.2.840.10008.5.1.4.1.1.12.1", "X-Ray Angiographic Image Storage"},
  {"1.2.840.10008.5.1.4.1.1.128", "Positron Emission Tomography Image Storage"},
  {"1.2.840.10008.5.1.4.1.1.130", "Enhanced PET Image Storage"},
  {"1.2.840.10008.5.1.4.1.1.128.1", "Legacy Converted Enhanced PET Image Storage"},
  {"1.2.840.10008.5.1.4.1.1.2", "CT Image Storage"},
  {"1.2.840.10008.5.1.4.1.1.20", "NM Image Storage"},
  {"1.2.840.10008.5.1.4.1.1.30", "Parametric Map Storage"},
  {"1.2.840.10008.5.1.4.1.1.3.1", "Ultrasound Multiframe Image Storage"},
  {"1.2.840.10008.5.1.4.1.1.4", "MR Image Storage"},
  {"1.2.840.10008.5.1.4.1.1.4.1", "Enhanced MR Image Storage"},
  {"1.2.840.10008.5.1.4.1.1.4.2", "MR Spectroscopy Storage"},
  {"1.2.840.10008.5.1.4.1.1.6.1", "Ultrasound Image Storage"},
  {"1.2.840.10008.5.1.4.1.1.66", "Raw Data Storage"},
  {"1.2.840.10008.5.1.4.1.1.66.1", "Spatial Registration Storage"},
  {"1.2.840.10008.5.1.4.1.1.66.2", "Spatial Fiducials Storage"},
  {"1.2.840.10008.5.1.4.1.1.66.3", "Deformable Spatial Registration Storage"}, 
  {"1.2.840.10008.5.1.4.1.1.66.4", "Segmentation Storage"},
  {"unknown", "unknown SOP"} // do not change this 
};
/*! @endcond */
/******************************************************************************/

/******************************************************************************/
/** Verify that given file (either file name or file pointer) appears to be
    DICOM file, based on the magic number.
    @sa dcmReadFile, dcmReadTransferSyntaxUID
    @return Returns 1 if DICOM magic number can be found, or 0 if not.
    @author Vesa Oikonen
 */
int dcmVerifyMagic(
  /** Name of file to open and to verify for the magic number; enter NULL to
      use the file pointer (next argument) instead. */
  const char *filename,
  /** File pointer of file to check, opened with fp=fopen(filename, "rb"); 
      enter NULL to open file (previous argument) locally. 
      Previously opened file pointer is first rewound to start; if DICOM magic
      number is found, then file pointer is left to the end of magic number,
      and if not, it is rewound to the file start.
  */ 
  FILE *fp
) {
  FILE *lfp;

  if(filename==NULL && fp==NULL) return(0);
  if(fp!=NULL) {
    lfp=fp; rewind(lfp);
  } else {
    lfp=fopen(filename, "rb");
  }
  if(lfp==NULL) return(0);

  /* Skip the first 128 bytes */
  if(fseek(lfp, 128, SEEK_SET)) {if(fp==NULL) fclose(lfp); return(0);}

  /* Read the next 4 bytes as characters */
  char buf[5];
  size_t n=fread(&buf, 1, 4, lfp);
  buf[4]=(char)0;
  if(n!=(size_t)4) {rewind(lfp); return(0);}

  /* Check the magic number */
  if(strncmp(buf, "DICM", 4)==0) {
    if(fp==NULL) fclose(lfp);
    return(1);
  } else {
    if(fp==NULL) fclose(lfp); else rewind(lfp);
    return(0);
  }
}
/******************************************************************************/

/******************************************************************************/
/** Is the explicit VR (2 bytes) followed by reserved 2 bytes?
    If yes, then the following Value Length is also given as 32-byte integer,
    if no, then as 16-bit integer.
    @return Returns 0, if not, and 2, if it is.
 */ 
unsigned char dcmVRReserved(
  /** VR id (DCM_VR_AE, ...). */ 
  dcmvr id
) {
  unsigned short int i=0;
  while(dcm_vr[i].vr!=DCM_VR_INVALID) {
    if(id==dcm_vr[i].vr) return(dcm_vr[i].res);
    i++;
  }
  return(2);
}
/*****************************************************************************/

/*****************************************************************************/
/** Identify the DICOM VR based on the two-character long string.
 *  @return Returns the VR id.
 *  @sa dcmVRName
 */
dcmvr dcmVRId(
  /** VR string. Two first characters are used. 
      String does not need to be null-terminated. */
  const char *s
) {
  if(s==NULL) return(DCM_VR_INVALID);
  char buf[3]; buf[0]=s[0]; buf[1]=s[1]; buf[2]=(char)0;

  /* Identify the VR */
  unsigned short int i=0;
  while(dcm_vr[i].vr!=DCM_VR_INVALID) {
    if(strncmp(dcm_vr[i].name, buf, 2)==0) return(dcm_vr[i].vr);
    i++;
  }
  return(DCM_VR_INVALID);
}
/*****************************************************************************/

/*****************************************************************************/
/** Get the DICOM VR name.
 *  @return Returns pointer to the name string.
 *  @sa dcmIdentifyVR
 */
char *dcmVRName(
  /** VR id (DCM_VR_AE, ...). */ 
  dcmvr id
) {
  unsigned short int i=0;
  while(dcm_vr[i].vr!=DCM_VR_INVALID) {
    if(id==dcm_vr[i].vr) return(dcm_vr[i].name);
    i++;
  }
  return(dcm_vr[DCM_VR_INVALID].name);
}
/*****************************************************************************/

/*****************************************************************************/
/** Get the DICOM VR max value length in bytes; 0 if not defined.
 *  @return Returns the length in bytes.
 *  @sa dcmIdentifyVR
 */
size_t dcmVRVLength(
  /** VR id (DCM_VR_AE, ...). */ 
  dcmvr id
) {
  unsigned short int i=0;
  while(dcm_vr[i].vr!=DCM_VR_INVALID) {
    if(id==dcm_vr[i].vr) return(dcm_vr[i].s);
    i++;
  }
  return(dcm_vr[DCM_VR_INVALID].s);
}
/*****************************************************************************/

/*****************************************************************************/
/** Get the DICOM VR description.
 *  @return Returns pointer to the description string.
 *  @sa dcmIdentifyVR
 */
char *dcmVRDescr(
  /** VR id (DCM_VR_AE, ...). */ 
  dcmvr id
) {
  unsigned short int i=0;
  while(dcm_vr[i].vr!=DCM_VR_INVALID) {
    if(id==dcm_vr[i].vr) return(dcm_vr[i].descr);
    i++;
  }
  return(dcm_vr[DCM_VR_INVALID].descr);
}
/*****************************************************************************/

/*****************************************************************************/
/** Convert DICOM date 'DA' to international format YYYY-MM-DD.
 *  @return Returns pointer to the date string, or NULL in case of an error.
 */
char *dcmDA2intl(
  /** Pointer to original DICOM string. */ 
  const char *orig,
  /** Pointer to string where date in international format will be written;
      must be allocated for at least 11 characters. */
  char *intl 
) {
  if(orig==NULL || intl==NULL) return(NULL);
  if(strnlen(orig, 10)<8) return(NULL);
  if(isdigit(orig[4])) { // modern format YYYYMMDD
    sprintf(intl, "%4.4s-%2.2s-%2.2s", orig, orig+4, orig+6);
  } else { // old format YYYY.MM.DD
    sprintf(intl, "%4.4s-%2.2s-%2.2s", orig, orig+5, orig+8);
  }
  if(isdate(intl)) {intl[0]=(char)0; return(NULL);}
  return(intl);
}
/*****************************************************************************/

/*****************************************************************************/
/** Convert DICOM time 'TM' to international format hh:mm:ss.
 *  @return Returns pointer to the time string, or NULL in case of an error.
 */
char *dcmTM2intl(
  /** Pointer to original DICOM string. */ 
  const char *orig,
  /** Pointer to string where time in international format will be written;
      must be allocated for at least 9 characters. */
  char *intl 
) {
  if(orig==NULL || intl==NULL) return(NULL);
  if(strnlen(orig, 14)<6) return(NULL);
  if(isdigit(orig[2])) { // modern format hhmmss.fract
    sprintf(intl, "%2.2s:%2.2s:%2.2s", orig, orig+2, orig+4);
  } else { // old format hh.mm.ss
    sprintf(intl, "%2.2s:%2.2s:%2.2s", orig, orig+3, orig+6);
  }
  if(istime(intl)) {intl[0]=(char)0; return(NULL);}
  return(intl);
}
/*****************************************************************************/

/*****************************************************************************/
/** Convert DICOM datetime 'DT' to international format YYYY-MM-DD hh:mm:ss.
 *  @return Returns pointer to the time string, or NULL in case of an error.
 */
char *dcmDT2intl(
  /** Pointer to original DICOM string. 
      Should be in format YYYYMMDDhhmmss.FFFFFF+hhmm */ 
  const char *orig,
  /** Pointer to string where date and time in international format will be
      written; must be allocated for at least 20 characters. */
  char *intl 
) {
  if(orig==NULL || intl==NULL) return(NULL);
  if(strnlen(orig, 26)<14) return(NULL);
  sprintf(intl, "%4.4s-%2.2s-%2.2s %2.2s:%2.2s:%2.2s",
          orig, orig+4, orig+6, orig+8, orig+10, orig+12);
  if(isdatetime(intl, NULL)) {intl[0]=(char)0; return(NULL);}
  return(intl);
}
/*****************************************************************************/

/*****************************************************************************/
/** Identify the DICOM SOP UID.
 *  @return Returns the SOP UID list index.
 *  @sa dcmSOPName, dcmSOPUID
 */
unsigned int dcmSOPIdentify(
  /** SOP UID string. */
  const char *s
) {
  if(s==NULL || strnlen(s, 3)<3) return(0);

  /* Identify the SOP UID */
  unsigned int i=0;
  while(strcmp(dcm_sop[i].uid, "unknown")!=0) {
    if(strcmp(dcm_sop[i].uid, s)==0) return(i);
    i++;
  }
  return(i);
}
/*****************************************************************************/

/*****************************************************************************/
/** Get the DICOM SOP UID Name.
 *  @return Returns pointer to the SOP Name string.
 *  @sa dcmSOPIdentify, dcmSOPUID
 */
char *dcmSOPName(
  /** SOP UID list index. */ 
  unsigned int i
) {
  unsigned int j=0;
  while(strcmp(dcm_sop[j].uid, "unknown")!=0) {
    if(i==j) return(dcm_sop[j].name);
    j++;
  }
  return(dcm_sop[j].name);
}
/*****************************************************************************/

/*****************************************************************************/
/** Get the DICOM SOP UID.
 *  @return Returns pointer to the UID string.
 *  @sa dcmSOPIdentify, dcmSOPName
 */
char *dcmSOPUID(
  /** SOP UID list index. */ 
  unsigned int i
) {
  unsigned int j=0;
  while(strcmp(dcm_sop[j].uid, "unknown")!=0) {
    if(i==j) return(dcm_sop[j].uid);
    j++;
  }
  return(dcm_sop[j].uid);
}
/*****************************************************************************/

/*****************************************************************************/
/** Get the name of DICOM SOP UID.
 *  @return Returns pointer to the SOP UID name string.
 *  @sa dcmSOPName, dcmSOPUID
 */
char *dcmSOPUIDName(
  /** SOP UID string. */
  const char *s
) {
  if(s==NULL || strnlen(s, 3)<3) return(dcm_sop[0].name);

  /* Identify the SOP UID */
  unsigned int i=0;
  while(strcmp(dcm_sop[i].uid, "unknown")!=0) {
    if(strcmp(dcm_sop[i].uid, s)==0) return(dcm_sop[i].name);
    i++;
  }
  return(dcm_sop[i].name);
}
/*****************************************************************************/

/*****************************************************************************/
/** Identify the DICOM Transfer Syntax UID.
 *  @return Returns the enumerated id.
 *  @sa dcmTrUIDDescr
 */
dcmtruid dcmTrUID(
  /** UID string. */
  const char *s
) {
  if(s==NULL || strnlen(s, 5)<5) return(DCM_TRUID_INVALID);

  /* Identify the UID */
  unsigned short int i=1;  // 1 because 0 is unknown
  while(dcm_truid[i].id!=DCM_TRUID_INVALID) {
    if(strcmp(dcm_truid[i].uid, s)==0) return(dcm_truid[i].id);
    i++;
  }
  return(DCM_TRUID_UNKNOWN);
}
/*****************************************************************************/

/*****************************************************************************/
/** Get the DICOM Transfer Syntax UID description.
 *  @return Returns pointer to the description string.
 *  @sa dcmTrUID
 */
char *dcmTrUIDDescr(
  /** Transfer Syntax UID id (DCM_TRUID_LEI, ...). */ 
  dcmtruid id
) {
  unsigned short int i=0;
  while(dcm_truid[i].id!=DCM_TRUID_INVALID) {
    if(id==dcm_truid[i].id) return(dcm_truid[i].descr);
    i++;
  }
  return(dcm_truid[DCM_TRUID_INVALID].descr);
}
/*****************************************************************************/

/*****************************************************************************/
/** Get the DICOM Transfer Syntax UID string.
 *  @return Returns pointer to the description string.
 *  @sa dcmTrUID
 */
char *dcmTrUIDString(
  /** Transfer Syntax UID id (DCM_TRUID_LEI, ...). */ 
  dcmtruid id
) {
  unsigned short int i=0;
  while(dcm_truid[i].id!=DCM_TRUID_INVALID) {
    if(id==dcm_truid[i].id) return(dcm_truid[i].uid);
    i++;
  }
  return(dcm_truid[DCM_TRUID_INVALID].uid);
}
/*****************************************************************************/

/*****************************************************************************/
/** Read and identify the DICOM Transfer Syntax UID.
    @return Returns the enumerated UID, or DCM_TRUID_INVALID.
    @sa dcmVerifyMagic
 */
dcmtruid dcmReadTransferSyntaxUID(
  /** Pointer to DICOM file opened in binary format, and positioned right
      after the Magic number. 
      Position will be returned to this position. */
  FILE *fp
) {
  if(fp==NULL || feof(fp)) return(DCM_TRUID_INVALID);
  /* Save original file position */
  fpos_t opos;
  if(fgetpos(fp, &opos)) return(DCM_TRUID_INVALID);
  /* Read file until we find DICOM tag 0x0002,0x0010 */
  DCMTAG tag;
  int tag_found=0;
  while(!tag_found && !feof(fp)) {
    if(dcmReadFileTag(fp, &tag)) break;
    if(tag.group==0x0002 && tag.element==0x0010) {tag_found=1; break;}
    /* not the tag that we want, so just go through this item */
    dcmvr vr=dcmReadFileVR(fp, NULL);
    unsigned int vl=0;
    if(dcmVRReserved(vr)==0) vl=dcmReadFileVL(fp, 2); 
    else vl=dcmReadFileVL(fp, 4); 
    if(vr==DCM_VR_SQ) break;
    if(vl==0xFFFFFFFF) break;
    char buf[vl+1];
    if(fread(buf, 1, vl, fp)!=vl) break;
  } // get next tag
  if(!tag_found) {fsetpos(fp, &opos); return(DCM_TRUID_INVALID);}
  /* Read the UID */
  if(dcmReadFileVR(fp, NULL)!=DCM_VR_UI) {
    fsetpos(fp, &opos); return(DCM_TRUID_INVALID);
  }
  unsigned int vl=dcmReadFileVL(fp, 2); 
  if(vl==0 || vl==0xFFFFFFFF) {fsetpos(fp, &opos); return(DCM_TRUID_INVALID);}
  char uid[vl+1];
  if(fread(uid, 1, vl, fp)!=vl) {fsetpos(fp, &opos); return(DCM_TRUID_INVALID);}
  /* Return file position to the original */
  fsetpos(fp, &opos);
  /* Identify the UID */
  return(dcmTrUID(uid));
}
/*****************************************************************************/

/*****************************************************************************/
/** Read DICOM tag from current file position.
    @note Tag validity is not verified here; error may be caused by end-of-file.
     Alternatively, read may seem to be successful, but tag contents are 
     file padding symbols (0xFFFC).
    @sa dcmReadFile, dcmReadFileElement
    @return Returns 0 when successful, otherwise >0.
    @author Vesa Oikonen
 */
int dcmReadFileTag(
  /** File pointer, positioned at the start of the element (and tag). */
  FILE *fp,
  /** Pointer to DICOM tag struct; enter NULL if you just need to move
      file position over the tag. */
  DCMTAG *tag
) {
  if(tag!=NULL) {tag->group=tag->element=0xFFFC;} // padding
  if(fp==NULL || feof(fp)) return(1);
  unsigned short int buf[2];
  size_t n=fread(&buf, 2, 2, fp);
  if(n!=2) return(2+n);
  if(tag!=NULL) {
    tag->group=buf[0];
    tag->element=buf[1];
    if(!little_endian()) { // tag is by default little endian
      swap(&tag->group, &tag->group, 2);
      swap(&tag->element, &tag->element, 2);
    }
  }
  return(0);
}
/*****************************************************************************/

/*****************************************************************************/
/** Write DICOM tag into current file position.
    @sa dcmWriteFile, dcmReadFileTag
    @return 0 when successful.
    @author Vesa Oikonen
 */
int dcmWriteFileTag(
  /** File pointer, at the write position. */
  FILE *fp,
  /** Pointer to DICOM tag struct to write. */
  DCMTAG *tag
) {
  if(fp==NULL || tag==NULL) return(1);
  unsigned short int buf[2];
  buf[0]=tag->group;
  buf[1]=tag->element;
  if(!little_endian()) swabip(buf, 4); // tag is by default little endian
  if(fwrite(&buf, 2, 2, fp)!=2) return(2);
  return(0);
}
/*****************************************************************************/

/*****************************************************************************/
/** Write DICOM Sequence delimitation item into current file position.

    This item consists of four byte sequence delimitation tag (0xFFFE, 0xE0DD)
    and four byte item length (0x00000000), i.e. together 8 bytes.

    @sa dcmWriteFileTag, dcmReadFileTag, dcmWriteFileSQItemDelimTag
    @return 0 when successful.
    @author Vesa Oikonen
 */
int dcmWriteFileSQDelimItem(
  /** File pointer, at the write position. */
  FILE *fp
) {
  if(fp==NULL) return(1);
  int ret;
  DCMTAG tag; 
  tag.group=0xFFFE; tag.element=0xE0DD;
  ret=dcmWriteFileTag(fp, &tag); if(ret!=0) return(ret);
  tag.group=0x0000; tag.element=0x0000;
  ret=dcmWriteFileTag(fp, &tag); if(ret!=0) return(ret);
  return(0);
}
/*****************************************************************************/

/*****************************************************************************/
/** Write DICOM Sequence Item Delimitation Tag with VL into current file position.

    This tag consists of four bytes, sequence item delimitation tag (0xFFFE, 0xE00D),
    followed by four byte item length (0x00000000), i.e. together 8 bytes.

    @sa dcmWriteFileTag, dcmReadFileTag, dcmWriteFileSQDelimItem
    @return 0 when successful.
    @author Vesa Oikonen
 */
int dcmWriteFileSQItemDelimTag(
  /** File pointer, at the write position. */
  FILE *fp
) {
  if(fp==NULL) return(1);
  int ret;
  DCMTAG tag; 
  tag.group=0xFFFE; tag.element=0xE00D;
  ret=dcmWriteFileTag(fp, &tag); if(ret!=0) return(ret);
  tag.group=0x0000; tag.element=0x0000;
  ret=dcmWriteFileTag(fp, &tag); if(ret!=0) return(ret);
  return(0);
}
/*****************************************************************************/

/*****************************************************************************/
/** Read DICOM value representation (2 or 4 bytes) from current file position.
    @sa dcmReadFileTag, dcmReadFileElement
    @return Returns the enumerated VR number, DCM_VR_INVALID in case of an error.
    @author Vesa Oikonen
 */
dcmvr dcmReadFileVR(
  /** File pointer, positioned at the VR start. */
  FILE *fp,
  /** Pointer for VR string, allocated for at least 3 characters to have space
      for the trailing null; enter NULL if not needed. */
  char *vrstr 
) {
  if(vrstr!=NULL) vrstr[0]=(char)0;
  if(fp==NULL) return(DCM_VR_INVALID);

  /* Read the first two bytes */
  char buf[3];
  if(fread(&buf, 1, 2, fp)!=2) return(DCM_VR_INVALID);
  buf[2]=(char)0;

  /* Identify the VR */
  dcmvr lvr=dcmVRId(buf);
  if(vrstr!=NULL) {
    if(lvr!=DCM_VR_INVALID) strcpy(vrstr, dcmVRName(lvr));
    else strcpy(vrstr, buf);
  }

  /* If this VR has extra 2 byte reserved space, then
     we need to read but do not use the next 2 bytes. */
  if(dcmVRReserved(lvr)!=0) {
    if(fread(&buf, 1, 2, fp)!=2) return(DCM_VR_INVALID);
  }
  return(lvr);
}
/*****************************************************************************/

/*****************************************************************************/
/** Read DICOM value length (2 or 4 bytes, depending on VR) from current file position.
    @sa dcmReadFileTag, dcmReadFileVR, dcmReadFileElement
    @return Returns the value length.
    @author Vesa Oikonen
 */
unsigned int dcmReadFileVL(
  /** File pointer, positioned at the VL start. */
  FILE *fp,
  /** Number of bytes (2 or 4) in the VL representation. */
  unsigned int n
) {
  unsigned int vl=0;
  if(fp==NULL || (n!=2 && n!=4)) return(vl);

  /* Read 2 or 4 bytes */
  if(n==2) {
    unsigned short int si;
    if(fread(&si, 2, 1, fp)!=1) return(vl);
    if(!little_endian()) swap(&si, &si, 2);
    vl=si;
  } else if(n==4) {
    unsigned int li;
    if(fread(&li, 4, 1, fp)!=1) return(vl);
    if(!little_endian()) swap(&li, &li, 4);
    vl=li;
  }
  return(vl);
}
/*****************************************************************************/

/*****************************************************************************/
/** Read DICOM Value Representation (VR, 2 or 4 bytes) and Value Length (VL, 2 or 4 bytes)
    from current file position.
   @sa dcmReadFileVR, dcmReadFileVL, dcmReadFileTag, dcmReadFileElement, dcmWriteFileVRVL
   @return 0 when successful.
   @author Vesa Oikonen
 */
int dcmReadFileVRVL(
  /** File pointer, positioned at the VR start. */
  FILE *fp,
  /** Pointer for enumerated VR value; enter NULL if not needed (file pointer moved anyway). */
  dcmvr *vr,
  /** Pointer for VL; enter NULL if not needed (file pointer moved anyway). */
  unsigned int *vl,
  /** Pointer for number of bytes read from file; enter NULL if not needed. */
  unsigned int *n
) {
  if(vr!=NULL) *vr=DCM_VR_INVALID;
  if(vl!=NULL) *vl=0;
  if(n!=NULL) *n=0;
  if(fp==NULL) return(1);

  /* Read the first two bytes */
  char buf[3];
  if(fread(&buf, 1, 2, fp)!=2) return(2); else if(n!=NULL) *n+=2;
  buf[2]=(char)0;

  /* Identify the VR */
  dcmvr lvr=dcmVRId(buf);
  if(vr!=NULL) *vr=lvr;
  if(lvr==DCM_VR_INVALID) return(3);

  /* If this VR has extra 2 byte reserved space, then
     we need to read but do not use the next 2 bytes. */
  unsigned int bsize=2+dcmVRReserved(lvr);
  if(bsize==4) {
    if(fread(&buf, 1, 2, fp)!=2) return(2);
    if(n!=NULL) *n+=2;
  }

  /* Read VL from the next 2 or 4 bytes */
  unsigned int lvl=0;
  if(bsize==2) {
    unsigned short int si;
    if(fread(&si, 2, 1, fp)!=1) return(2);
    if(!little_endian()) swap(&si, &si, 2);
    lvl=si;
  } else {
    unsigned int li;
    if(fread(&li, 4, 1, fp)!=1) return(2);
    if(!little_endian()) swap(&li, &li, 4);
    lvl=li;
  }
  if(n!=NULL) *n+=bsize;
  if(vl!=NULL) *vl=lvl;

  return(0);
}
/*****************************************************************************/

/*****************************************************************************/
/** Write DICOM Value Representation (VR, 2 or 4 bytes) and Value Length (VL, 2 or 4 bytes)
    into current file position.
   @sa dcmWriteFileTag, dcmReadFileElement, dcmWriteFileVRVL
   @return 0 when successful.
   @author Vesa Oikonen
 */
int dcmWriteFileVRVL(
  /** File pointer, opened fro writing in binary mode. */
  FILE *fp,
  /** Enumerated VR value. */
  dcmvr vr,
  /** VL */
  unsigned int vl,
  /** Pointer for number of bytes written into file; enter NULL if not needed. */
  unsigned int *n
) {
  if(n!=NULL) *n=0;
  if(fp==NULL || vr==DCM_VR_INVALID) return(1);

  /* If this VR has extra 2 byte reserved space, then
     we need to write VR and VL with 4 bytes each, other wise with 2 bytes each. */
  unsigned int bsize=2+dcmVRReserved(vr);

  char buf[10];

  /* Write VR */
  memcpy(buf, dcmVRName(vr), 2); buf[2]=buf[3]=(char)0;
  if(fwrite(buf, bsize, 1, fp)!=1) return(2);
  if(n!=NULL) *n+=bsize;

  /* Write VL */
  memcpy(buf, &vl, bsize);
  if(bsize==2) buf[2]=buf[3]=(char)0;
  if(!little_endian()) swap(buf, buf, bsize);
  if(fwrite(buf, bsize, 1, fp)!=1) return(2);
  if(n!=NULL) *n+=bsize;

  return(0);
}
/*****************************************************************************/

/*****************************************************************************/
/** Initiate the DCMFILE struct before any use. 
    @sa dcmfileFree
 */
void dcmfileInit(
  /** Pointer to DCMFILE. */
  DCMFILE *d
) {
  if(d==NULL) return;
  d->filename[0]=(char)0;
  d->fp=(FILE*)NULL;
  d->truid=DCM_TRUID_UNKNOWN;
  d->item=(DCMITEM*)NULL;
}
/*****************************************************************************/

/*****************************************************************************/
/** Recursively free memory allocated for DCMITEM items and their children items.
    @sa dcmfileFree
    @author Vesa Oikonen
 */
void dcmitemFree(
  /** Pointer to DCMITEM. */
  DCMITEM *d
) {
  if(d==NULL) return;
  /* find the last item in the list */
  DCMITEM *ip=d; while(ip->next_item!=NULL) ip=ip->next_item;
  while(ip!=NULL) {
    /* Free items child and their children */
    if(ip->child_item!=NULL) dcmitemFree(ip->child_item);
    /* Free this item and move to previous item */
    if(ip->prev_item!=NULL) {
      ip=ip->prev_item; 
      free(ip->next_item->rd); free(ip->next_item); 
      ip->next_item=NULL;
    } else {
      free(ip->rd); free(ip); ip=NULL;
    }
  }
}
/*****************************************************************************/

/*****************************************************************************/
/** Free memory allocated for DCMFILE data. All contents are destroyed.
    @pre Before first use initialize the struct with dcmfileInit().
    @sa dcmfileInit
    @author Vesa Oikonen
 */
void dcmfileFree(
  /** Pointer to DCMFILE. */
  DCMFILE *d
) {
  if(d==NULL) return;
  dcmitemFree(d->item);
  dcmfileInit(d);
}
/*****************************************************************************/

/*****************************************************************************/
/** Get the maximum depth of DCMITEM tree.
    @return Returns the number of levels _under_ specified item, not including
     given item itself.
    @sa dcmitemParentNr, dcmfileMaxDepth
 */
unsigned short int dcmitemMaxDepth(
  /** Pointer to DCMITEM item. */
  DCMITEM *d
) {
  if(d==NULL || d->child_item==NULL) return(0);
  unsigned short int m=0, n=0;
  DCMITEM *cd=d->child_item;
  /* go through all children */
  while(cd!=NULL) {
    n=dcmitemMaxDepth(cd); if(n>m) m=n;
    cd=cd->next_item;
  }
  return(m+1);
}
/** Get the maximum depth of DCMFILE items tree.
    @return Returns the number of item levels under specified DCMFILE, or
     zero in case of an error or if there are no items.
    @sa dcmitemParentNr, dcmitemMaxDepth
 */
unsigned short int dcmfileMaxDepth(
  /** Pointer to DCMFILE item. */
  DCMFILE *df
) {
  if(df==NULL || df->item==NULL) return(0);
  unsigned short int m=0, n=0;
  DCMITEM *sd=df->item;
  while(sd!=NULL) { // go through all sisters
    n=dcmitemMaxDepth(sd); if(n>m) m=n;
    sd=sd->next_item;
  }
  return(m+1);
}
/*****************************************************************************/

/*****************************************************************************/
/** Check how deep in DCMITEM tree this item is.
    @return Returns the number of parents this item has.
    @sa dcmitemMaxDepth, dcmfileMaxDepth
 */
unsigned short int dcmitemParentNr(
  /** Pointer to DCMITEM. */
  DCMITEM *d
) {
  if(d==NULL) return(0);
  unsigned short int n=0;
  DCMITEM *pd=d->parent_item;
  while(pd!=NULL) {n++; pd=pd->parent_item;}
  return(n);
}
/*****************************************************************************/

/*****************************************************************************/
/** Pre-process the DICOM element value into format suitable for printing.
    @note Use only for printing information for the user.
    @return Returns pointer to locally allocated null-terminated string.
    @post Free the pointer after use.
    @sa dcmitemGetInt, dcmitemGetReal, dcmFindTag
 */ 
char *dcmValueString(
  /** Pointer to item containing the value to print. */
  DCMITEM *d
) {
  if(d==NULL) return((char*)NULL);

  /* For sequence, return string 'na' */
  if(d->vr==DCM_VR_SQ) {
    char *s=malloc(3); strcpy(s, "na"); // do not write char *s="na";
    return(s);
  }

  /* If there is no value, then return string 'empty', or
     'na', if value just was not stored (pixel data) */
  if(d->vl==0) {
    char *s=malloc(6); strcpy(s, "empty"); 
    return(s);
  } else if(d->rd==NULL) {
    char *s=malloc(3); strcpy(s, "na"); 
    return(s);
  }

  unsigned int len;
  if(d->vl==0xFFFFFFFF) len=(unsigned int)dcmVRVLength(d->vr); else len=d->vl;

  /* String values */
  if(d->vr==DCM_VR_CS || d->vr==DCM_VR_DS || d->vr==DCM_VR_IS || 
     d->vr==DCM_VR_LO || d->vr==DCM_VR_LT || d->vr==DCM_VR_PN || 
     d->vr==DCM_VR_SH || d->vr==DCM_VR_ST)
  {
    char *s=malloc(len+1);
    memcpy(s, d->rd, len); s[len]=(char)0;
    return(s);
  }

  /* More string values */
  if(d->vr==DCM_VR_AS || d->vr==DCM_VR_PN || 
     d->vr==DCM_VR_DA || d->vr==DCM_VR_DT || d->vr==DCM_VR_TM || 
     d->vr==DCM_VR_UT || d->vr==DCM_VR_AE || 
     d->vr==DCM_VR_UI || d->vr==DCM_VR_UR)
  {
    char *s=malloc(len+1);
    memcpy(s, d->rd, len); s[len]=(char)0;
    return(s);
  }

  /** Attribute tag */
  if(d->vr==DCM_VR_AT) {
    DCMTAG tag;
    memcpy(&tag.group, d->rd, 2);
    memcpy(&tag.element, d->rd+2, 2);
    if(!little_endian()) {
      swap(&tag.group, &tag.group, 2);
      swap(&tag.element, &tag.element, 2);
    }
    char *s=malloc(14);
    sprintf(s, "0x%04x,0x%04x", tag.group, tag.element);
    return(s);
  }

  /** Float */
  if(d->vr==DCM_VR_FL) {
    float f;
    memcpy(&f, d->rd, 4);
    if(!little_endian()) swap(&f, &f, 4);
    char *s=malloc(16);
    sprintf(s, "%g", f);
    return(s);
  }

  /** Double */
  if(d->vr==DCM_VR_FD) {
    char *s=malloc(32);
    double f;
    memcpy(&f, d->rd, 8);
    if(!little_endian()) swap(&f, &f, 8);
    sprintf(s, "%g", f);
    return(s);
  }

  /** Unsigned 32-bit int */
  if(d->vr==DCM_VR_UL) {
    unsigned int i;
    memcpy(&i, d->rd, 4);
    if(!little_endian()) swap(&i, &i, 4);
    char *s=malloc(16);
    sprintf(s, "%u", i);
    return(s);
  }

  /** Unsigned 16-bit int */
  if(d->vr==DCM_VR_US) {
    unsigned short int i;
    memcpy(&i, d->rd, 2);
    if(!little_endian()) swap(&i, &i, 2);
    char *s=malloc(8);
    sprintf(s, "%u", i);
    return(s);
  }

  /** Signed 32-bit int */
  if(d->vr==DCM_VR_SL) {
    int i;
    memcpy(&i, d->rd, 4);
    if(!little_endian()) swap(&i, &i, 4);
    char *s=malloc(16);
    sprintf(s, "%d", i);
    return(s);
  }

  /** Signed 16-bit int */
  if(d->vr==DCM_VR_SS) {
    short int i;
    memcpy(&i, d->rd, 2);
    if(!little_endian()) swap(&i, &i, 2);
    char *s=malloc(8);
    sprintf(s, "%d", i);
    return(s);
  }

/* Not (yet) printed:
  DCM_VR_OB,        ///< DICOM other byte string, even bytes, endian insensitive.
  DCM_VR_OD,        ///< DICOM other double (64-bit) stream, endian sensitive.
  DCM_VR_OF,        ///< DICOM other float (32-bit) stream, endian sensitive.
  DCM_VR_OL,        ///< DICOM other long (32-bit) stream, endian sensitive.
  DCM_VR_OW,        ///< DICOM other word (16-bit) stream, even bytes, endian sensitive.
  DCM_VR_UC,        ///< DICOM unlimited characters.
  DCM_VR_UN,        ///< DICOM unknown, any valid length of another VR.
  DCM_VR_INVALID    ///< Invalid DICOM value representation.
*/
  char *s=malloc(3); strcpy(s, "na"); // do not write char *s="na";
  return(s);
}
/*****************************************************************************/

/*****************************************************************************/
/** Read integer value from given DICOM item.
    @details VR must be either UL, US, SL, SS, or IS; otherwise 0 is returned.
    @return Returns the value as long int, in order to cope with originally
     unsigned integers.
    @sa dcmitemGetInt, dcmValueString, dcmFindTag
 */
long int dcmitemGetInt(
  /** Pointer to item. */ 
  DCMITEM *d
) {
  if(d==NULL || d->rd==NULL) return(0);
  long int li=0;
  if(d->vr==DCM_VR_UL) { // unsigned 32-bit int
    unsigned int i;
    memcpy(&i, d->rd, 4); if(!little_endian()) swap(&i, &i, 4);
    li=(long int)i;
  } else if(d->vr==DCM_VR_US) { // unsigned 16-bit int
    unsigned short int i;
    memcpy(&i, d->rd, 2); if(!little_endian()) swap(&i, &i, 2);
    li=(long int)i;
  } else if(d->vr==DCM_VR_SL) { // signed 32-bit int
    int i;
    memcpy(&i, d->rd, 4); if(!little_endian()) swap(&i, &i, 4);
    li=(long int)i;
  } else if(d->vr==DCM_VR_SS) { // signed 16-bit int
    short int i;
    memcpy(&i, d->rd, 2); if(!little_endian()) swap(&i, &i, 2);
    li=(long int)i;
  } else if(d->vr==DCM_VR_IS) { // integer string
    li=atol(d->rd);
  }
  return(li);
}
/*****************************************************************************/

/*****************************************************************************/
/** Read floating point value from given DICOM item.
    @details VR must be either FL, FD, DS, UL, US, SL, SS, or IS; 
     otherwise 0 is returned.
    @return Returns the value as double.
    @sa dcmitemGetInt, dcmValueString, dcmFindTag
 */
double dcmitemGetReal(
  /** Pointer to item. */ 
  DCMITEM *d
) {
  if(d==NULL || d->rd==NULL) return(0);
  double r=0.0;
  if(d->vr==DCM_VR_FL) { // 32-bit float
    float f;
    memcpy(&f, d->rd, 4); if(!little_endian()) swap(&f, &f, 4);
    r=(double)f;
  } else if(d->vr==DCM_VR_FD) { // 64-bit double
    double f;
    memcpy(&f, d->rd, 8); if(!little_endian()) swap(&f, &f, 8);
    r=f;
  } else if(d->vr==DCM_VR_DS) { // decimal string
    r=atof(d->rd);
  } else if(d->vr==DCM_VR_UL) { // unsigned 32-bit int
    unsigned int i;
    memcpy(&i, d->rd, 4); if(!little_endian()) swap(&i, &i, 4);
    r=(double)i;
  } else if(d->vr==DCM_VR_US) { // unsigned 16-bit int
    unsigned short int i;
    memcpy(&i, d->rd, 2); if(!little_endian()) swap(&i, &i, 2);
    r=(double)i;
  } else if(d->vr==DCM_VR_SL) { // signed 32-bit int
    int i;
    memcpy(&i, d->rd, 4); if(!little_endian()) swap(&i, &i, 4);
    r=(double)i;
  } else if(d->vr==DCM_VR_SS) { // signed 16-bit int
    short int i;
    memcpy(&i, d->rd, 2); if(!little_endian()) swap(&i, &i, 2);
    r=(double)i;
  } else if(d->vr==DCM_VR_IS) { // integer string
    r=(double)atol(d->rd);
  }
  return(r);
}
/*****************************************************************************/

/*****************************************************************************/
/** Search for specified tag in DCMITEM data tree.
    @return Returns pointer to next item with the tag, or NULL if not found.
 */
DCMITEM *dcmFindTag(
  /** Pointer to current DICOM item. */
  DCMITEM *d,
  /** Omit this item from the search. */
  const short int omit,
  /** Pointer to the DICOM tag that is searched for. */
  DCMTAG *tag,
  /** Verbose level; if zero, then nothing is printed to stderr or stdout. */
  const int verbose
) {
  if(d==NULL || tag==NULL) return(NULL);
  if(verbose>0) printf("%s(%04X,%04X)\n", __func__, tag->group, tag->element);
  DCMITEM *iptr;
  if(omit==0) iptr=d; else iptr=d->next_item;
  while(iptr!=NULL) {
    if(verbose>2)
      printf(" checking tag(%04X,%04X)...\n", iptr->tag.group, iptr->tag.element);
    if(iptr->tag.group==tag->group && iptr->tag.element==tag->element) {
      if(verbose>2) printf("  found!\n");
      break;
    }
    /* Check if this item has children */
    if(iptr->child_item!=NULL) {
      if(verbose>2) printf("  going to search inside children...\n");
      DCMITEM *rptr=dcmFindTag(iptr->child_item, 0, tag, verbose);
      if(rptr!=NULL) return(rptr);
      if(verbose>3) printf("  nothing found in any of the children\n");
    }
    iptr=iptr->next_item;
  }
  /* Stop if we found tag, or if we do not have parent */
  if(iptr!=NULL) return(iptr);
  if(d->parent_item==NULL) return(NULL);

  /* Search from the parent */
  if(verbose>2) printf("  going to search inside parent...\n");
  return(dcmFindTag(d->parent_item, 1, tag, verbose));
}
/*****************************************************************************/

/*****************************************************************************/
/** Print contents of given DICOM item into stdout. */
void dcmitemPrint(
  /** Pointer to item. */
  DCMITEM *d
) {
  if(d==NULL) {printf("(null)\n"); fflush(stdout); return;}
  printf("tag(%04X,%04X)", d->tag.group, d->tag.element); fflush(stdout);
  printf(" VR=%s", dcmVRName(d->vr)); fflush(stdout);
  if(d->vl==0xFFFFFFFF) printf(" VL=%08X", d->vl); else printf(" VL=%u", d->vl);
  fflush(stdout);
  char *buf=dcmValueString(d); printf(" '%s'", buf); free(buf);
  printf("\n"); fflush(stdout);
}
/*****************************************************************************/

/*****************************************************************************/
/** Set DICOM Tag group and element. */
void dcmTagSet(
  /** Pointer to Tag to set. */
  DCMTAG *tag,
  /** Tag Group */
  unsigned short int group,
  /** Tag Element */
  unsigned short int element
) {
  tag->group=group;
  tag->element=element;
}
/*****************************************************************************/

/*****************************************************************************/
/** Add an item to DCMFILE data struct.
   @sa dcmfileInit, dcmfileFree
   @return 0 if successful, otherwise >0.
*/
int dcmAddItem(
  /** Pointer to DCMFILE. */
  DCMFILE *dcm,
  /** Pointer to a previous item in DCMFILE, into which to link this item;
      enter NULL to add as next item to the highest level. */
  DCMITEM *d,
  /** Add as child; 1=yes, 0=no. */
  short int aschild,
  /** Tag */
  DCMTAG tag,
  /** VR */
  dcmvr vr,
  /** VL; enter 0xFFFFFFFF to use VR's default length. */
  unsigned int vl,
  /** Pointer to the item value as byte array */
  char *rd,
  /** Verbose level; if zero, then nothing is printed to stderr or stdout. */
  const int verbose
) {
  if(verbose>0) {
    printf("%s(dcm, (%04X,%04X))", __func__, tag.group, tag.element);
    if(d==NULL) printf(", null"); else printf(", ptr");
    printf(", %d, %s, 0x%08X, rd", aschild, dcmVRName(vr), vl);
    printf(")\n");
  }
  if(dcm==NULL) return(1);
  if(vr==DCM_VR_INVALID) return(2);

  /* Check that caller has not given previous item pointer when there are none in DCMFILE */
  if(d!=NULL && dcm->item==NULL) return(3);
  /* Check that caller has not given previous item pointer that already has child */
  if(d!=NULL && aschild && d->child_item!=NULL) return(4);

  /* Check whether we currently support the Transfer UID */
  if(dcm->truid!=DCM_TRUID_LEE) return(5);

  /* Allocate memory for the new element; do not free it here, since it will be part of
     DCMFILE struct. */
  if(verbose>1) printf("  allocating memory for the item\n");
  DCMITEM *item=(DCMITEM*)malloc(sizeof(DCMITEM));
  if(item==NULL) return(11);
  item->next_item=item->child_item=(DCMITEM*)NULL;
  item->fp=dcm->fp; item->truid=dcm->truid;
  item->rd=(char*)NULL;

  /* Set item tag, VR, and VL */
  if(verbose>1) printf("  setting item contents\n");
  item->tag.group=tag.group;
  item->tag.element=tag.element;
  item->vr=vr;
  item->vl=vl;
  /* Allocate memory for item value */
  size_t s;
  if(vl==0xFFFFFFFF) s=dcmVRVLength(vr); else s=vl;
  if(s>0) {
    if(item->vl==0xFFFFFFFF) item->vl=(unsigned int)s;
    if(verbose>1) printf("  allocating %u bytes for the item value\n", (unsigned int)s);
    item->rd=(char*)calloc(s, sizeof(char));
    if(item->rd==NULL) {free(item); return(21);}
  } else {
    if(verbose>1) printf("zero size for item value\n");
    if(rd==NULL) {
      if(verbose>1) printf("... which is ok since value is empty, too.\n");
    } else {
      if(verbose>0) printf("... which is not ok because we have value to store.\n");
      if(item->rd==NULL) {free(item); return(22);}
    }
  }
  /* Copy the item value */
  if(rd!=NULL && s>0) {
    if(verbose>1) printf("  copying the item value\n");
    /* Special treatment for strings, because those tend to be shorter than told */
    if(vr==DCM_VR_LO || vr==DCM_VR_LT || vr==DCM_VR_PN || vr==DCM_VR_SH || vr==DCM_VR_UI || vr==DCM_VR_UR) 
    {
      unsigned int len=strnlen(rd, s);
      if(len<s) strlcpy(item->rd, rd, s);
      else memcpy(item->rd, rd, s);
    } else if(vr==DCM_VR_DS || vr==DCM_VR_IS) 
    {
      unsigned int len=strnlen(rd, s);
      if(len<s) strlcpy(item->rd, rd, s);
      else memcpy(item->rd, rd, s);
    } else {
      memcpy(item->rd, rd, s);
    }
  }


  /* If we have the item to link to, then do the linking */
  if(verbose>1) printf("  link the item.\n");
  if(d!=NULL) {
    if(aschild) {
      d->child_item=item;
      item->parent_item=d;
      item->prev_item=(DCMITEM*)NULL;
    } else if(d->next_item==NULL) {
      d->next_item=item;
      item->prev_item=d;
      item->parent_item=d->parent_item;
    } else {
      /* find the last item in the list */
      DCMITEM *ip=d; while(ip->next_item!=NULL) ip=ip->next_item;
      ip->next_item=item;
      item->prev_item=ip;
      item->parent_item=ip->parent_item;
    }
  } else if(dcm->item==NULL) {
    /* This is truly the first item ever */
    dcm->item=item;
    item->prev_item=item->parent_item=(DCMITEM*)NULL;
  } else {
    /* Caller lets us find the item to link to */
    DCMITEM *ip=dcm->item; while(ip->next_item!=NULL) ip=ip->next_item;
    ip->next_item=item;
    item->prev_item=ip;
    item->parent_item=ip->parent_item;
  }

  if(verbose>2) dcmitemPrint(item);

  if(verbose>1) printf("  all done.\n");
  return(0);
}
/*****************************************************************************/

/*****************************************************************************/
/** Read an element from DICOM file, and add it to the given linked list.
    This function will be called recursively in case of sequential items.
    @return 0 when successful, >0 in case of an error, -1 when no more elements could be read.
    @sa dcmFileRead
 */
int dcmFileReadNextElement(
  /** Pointer to DCMFILE struct; must be initiated before first call. */
  DCMFILE *dcm,
  /** Pointer to previous element; NULL if none exists (yet). */
  DCMITEM *prev_item,
  /** Pointer to parent element; NULL if none exists (yet). */
  DCMITEM *parent_item,
  /** Add as next element (0) or as child element. */
  const short int sub,
  /** Read only header (1), or read both header and pixel data (0). */
  const short int headerOnly, 
  /** Verbose level; if zero, then nothing is printed to stderr or stdout. */
  int verbose
) {
  if(verbose>0) printf("%s(DCMFILE*, DCMITEM*, DCMITEM*, %d, %d)\n", __func__, sub, headerOnly);
  if(dcm==NULL || dcm->fp==NULL) return(1);
  if(sub!=0 && parent_item==NULL) return(1);
  if(feof(dcm->fp)) return(-1);

  /* Check whether we currently support the Transfer UID */
  if(dcm->truid!=DCM_TRUID_LEE) return(2);

  if(verbose>10) {
    if(dcm->item==NULL) printf(" will add first element\n");
    else if(sub==0) printf(" will add next element\n");
    else printf(" will add subelement\n");
  }

  /* Is this a child to a sequence element? */
  int sq_child=0;
  if(parent_item!=NULL && parent_item->vr==DCM_VR_SQ) {sq_child=1;}
  if(verbose>10 && sq_child!=0) printf(" we're a child to a sequence element\n");

  /* Allocate memory for the new element */
  DCMITEM *item=(DCMITEM*)malloc(sizeof(DCMITEM));
  if(item==NULL) return(4);
  item->prev_item=prev_item;
  item->parent_item=parent_item;
  item->next_item=item->child_item=(DCMITEM*)NULL;
  item->fp=dcm->fp; item->truid=dcm->truid;
  item->rd=(char*)NULL;

  /* Save current file position (should be the start of element) */
  if(fgetpos(dcm->fp, &item->pos)) {
    free(item); return(2);
  }

  /* Read the tag (2x2 bytes) */
  if(verbose>10) {
    long int tagpos=ftell(dcm->fp);
    printf(" reading tag at %ld\n", tagpos);
  }
  if(dcmReadFileTag(dcm->fp, &item->tag)) {
    if(verbose>1 && !feof(dcm->fp)) printf(" error in reading the tag.\n");
    free(item);
    if(feof(dcm->fp)) return(-1);
    return(2);
  }

  if(verbose>2) {
    printf(" tag(%04x,%04x) with %u parents\n", 
           item->tag.group, item->tag.element, dcmitemParentNr(item));
  }

  /* If child, then check for item delimitation tag (although it should not be here) */
  if(dcmitemParentNr(item)>0 && item->tag.group==0xFFFE && item->tag.element==0xE00D) {
    if(verbose>10)
      printf(" item delimitation tag(%04x,%04x) found, reading VL\n", 
             item->tag.group, item->tag.element);
    unsigned long int vl=dcmReadFileVL(dcm->fp, 4); 
    if(verbose>1) printf(" item delimitation tag VL := %lu (0x%08lx)\n", vl, vl);
    if(vl!=0) {
      if(verbose>1) printf(" error: VL should have been 0\n");
      return(3);
    }
    return(0);
  }


  /* Read value representation and length (VR and VL, 2x2 or 2x4 bytes) */
  {
    if(verbose>10) printf(" reading VR and VL\n");
    int ret;
    unsigned int n;
    ret=dcmReadFileVRVL(dcm->fp, &item->vr, &item->vl, &n);
    if(ret!=0) {
      if(verbose>1) printf(" invalid VR or VL\n");
      free(item); return(ret);
    }
    if(verbose>1) {
      printf(" VR := %s (%s)\n", dcmVRName(item->vr), dcmVRDescr(item->vr));
      printf(" VL := %u (0x%08x) (%d bytes field)\n", item->vl, item->vl, n/2);
      fflush(stdout);
    }
  }

  /* Read value field, and add the current element to the list */
  if(item->vr==DCM_VR_SQ) {
    if(ftell(dcm->fp)<0) return(2);
    unsigned long int sqPos=(unsigned long int)ftell(dcm->fp);
    if(verbose>10) {printf(" sequence... at %ld\n", sqPos); fflush(stdout);}
    unsigned long int sqContentLength=item->vl;
    if(verbose>12) printf("    sequence contents length is %lu\n", sqContentLength);
    /* File position is now at the start of first item in the sequence */
    int ret;
    /* If parent has no previous child, then define this as its child */
    if(sq_child!=0 && parent_item->child_item==NULL) {
      parent_item->child_item=item;
    } else {
      /* else, add SQ sequence itself as next element to the list, and later
         add each sequence item as child to it */
      if(prev_item==NULL) {
        if(dcm->item==NULL) { // truly the first item
          dcm->item=item;
        } else { // search for the previous one
          DCMITEM *ip; ip=dcm->item;
          while(ip->next_item!=NULL) ip=ip->next_item;
          ip->next_item=item; item->prev_item=ip;
        }
      } else {
        prev_item->next_item=item; item->prev_item=prev_item;
      }
    }
    /* Read the first item tag and length, but there is no VR to read this time */
    if(verbose>10) {
      long int tagpos=ftell(dcm->fp);
      printf(" reading first item tag at %ld\n", tagpos);
    }
    DCMTAG itemtag;
    ret=dcmReadFileTag(dcm->fp, &itemtag);
    if(ret!=0) {
      if(verbose>1) printf(" error %d in reading the tag.\n", ret);
      return(2);
    }
    if(verbose>1) printf("  item tag(%04x,%04x)\n", itemtag.group, itemtag.element);
    /* It is common that sequence is empty; check it first */
    if(itemtag.group==0xFFFE && (itemtag.element==0xE0DD || itemtag.element==0xE00D)) {
      /* yes; then read also the 4 byte LV, which should be zero */
      if(verbose>10) 
        printf(" sequence delimitation item tag(%04x,%04x) found, reading VL\n",
          itemtag.group, itemtag.element);
      unsigned long int vl=dcmReadFileVL(dcm->fp, 4); 
      if(verbose>1) printf(" item tag VL := %lu (0x%08lx)\n", vl, vl);
      if(vl!=0) {
        if(verbose>1) printf(" error: VL should have been 0\n");
        return(2);
      }
      if(verbose>3) printf(" ending sequence before it really started.\n");
      return(0);
    }
    /* If sequence actually contains something, the Item tag must be 0xFFFE,0xE000 */
    if(itemtag.group!=0xFFFE || itemtag.element!=0xE000) {
      if(verbose>1) printf(" invalid sequence item tag(%04x,%04x)\n", itemtag.group, itemtag.element);
      return(2);
    }
    /* Read past the VL of this item (always 4 bytes) */
    unsigned long int itemvl=dcmReadFileVL(dcm->fp, 4);
    if(verbose>3) {printf(" item_VL := %lu (0x%08lx)\n", itemvl, itemvl);}
    if(ftell(dcm->fp)<0) return(2);
    unsigned long int sqItemPos=(unsigned long int)ftell(dcm->fp);
    /* Check if that is all of this sequence (probably Siemens) */
    if((sqItemPos-sqPos)>=sqContentLength) {
      if(verbose>3) printf(" ending sequence since it was found to be empty.\n");
      return(0);
    }
    if(verbose>12) printf("  sequence content start position at %ld\n", sqItemPos);
    /* Read the first item value as its own element, adding it as child to SQ */
    ret=dcmFileReadNextElement(dcm, NULL, item, 1, headerOnly, verbose-1);
    if(ret!=0) {
      if(verbose>1) printf(" error in reading the first item value dataset\n");
      return(ret);
    }
    /* Now we continue reading more items, until we reach Sequence Delimitation Item */
    while(!feof(dcm->fp)) {
      /* Do not read pass the length of the sequence data */
      if(ftell(dcm->fp)<0) return(2);
      unsigned long int cPos=(unsigned long int)ftell(dcm->fp);
      if(sqContentLength>0 && (cPos-sqPos)>=sqContentLength) {
        if(verbose>3) printf(" we reached the end of sequence VL %lu\n", sqContentLength);
        /* set fake sequence delimitation tag */
        itemtag.group=0xFFFE; itemtag.element=0xE0DD;
        break;
      }
      if(verbose>10) {
        long int tagpos=ftell(dcm->fp);
        printf(" reading next sequence item tag at %ld, %ld after start\n", tagpos, tagpos-sqItemPos);
      }
      if(dcmReadFileTag(dcm->fp, &itemtag)) return(2);
      if(verbose>1) printf(" next item tag(%04x,%04x)\n", itemtag.group, itemtag.element);
      itemvl=dcmReadFileVL(dcm->fp, 4); // delimitation tag has this too
      if(verbose>3) {printf(" item_VL := %lu (0x%08lx)\n", itemvl, itemvl);}
      /* Check if we got sequence delimitation tag */
      if(itemtag.group==0xFFFE && itemtag.element==0xE0DD) 
      {
        if(verbose>3) printf(" we got sequence delimitation tag\n");
        break;
      }
      /* Check if we got item delimitation tag from the previous item */
      if(itemtag.group==0xFFFE && itemtag.element==0xE00D) 
      {
        if(verbose>3) printf(" we got item delimitation tag\n");
        if(itemvl!=0) {
          if(verbose>1) printf(" error: VL should have been 0\n");
          return(3);
        }
        continue;
      }
      /* Otherwise this should be sequence item tag */
      if(itemtag.group!=0xFFFE || itemtag.element!=0xE000) {
        if(verbose>3) printf(" not sequence item tag, move file position back 2x4 bytes\n");
        fseek(dcm->fp, -8, SEEK_CUR);  //return(TPCERROR_INVALID_VALUE);
      }
      /* Read the item value as its own element, adding it to the SQ child list */
      DCMITEM *child=item->child_item;
      if(child==NULL) {
        if(verbose>1) printf(" error had happened in adding the child element\n");
        return(2);
      }
      while(child->next_item!=NULL) child=child->next_item;
      ret=dcmFileReadNextElement(dcm, child, item, 0, headerOnly, verbose-1);
      if(ret!=0) {
        if(verbose>1) printf(" error in reading item value dataset\n");
        return(ret);
      }
    }    
    /* Check that loop really stopped at sequence delimitation item */
    /* 0xE00D means the end of item, 0xE0DD the end of sequence. */ 
    if(itemtag.group!=0xFFFE || itemtag.element!=0xE0DD) {
      if(verbose>1)
        printf(" invalid sequence delimitation item tag(%04x,%04x)\n", itemtag.group, itemtag.element);
      return(2);
    }
    /* Done. Do not free item! */
    if(verbose>10) {printf(" end of sequence.\n"); fflush(stdout);}
  } else if(item->vl!=0xFFFFFFFF) {
    if(verbose>10) {printf(" reading value of %u bytes...\n", item->vl); fflush(stdout);}
    char *buf=NULL;
    if(item->vl>0) {
      buf=(char*)calloc(item->vl+1, sizeof(char));
      if(buf==NULL) {free(item); return(4);}
      if(fread(buf, 1, item->vl, item->fp)!=item->vl) {
        free(item); free(buf); return(3);
      }
      /* Do not store pixel data, if that was the request */
      if(headerOnly!=0 && 
         ((item->tag.group==0x7FE0 && item->tag.element>0) || item->tag.group==0x7FE1))
      {
        if(verbose>5) {printf(" ...not storing pixel data\n"); fflush(stdout);}
        free(buf); buf=(char*)NULL;
      } else {
        item->rd=buf;
      }
    } else if(verbose>4) {
      printf(" VL=0\n");
    }
    /* Add to list */
    if(sub==0) {
      if(prev_item==NULL) {
        if(dcm->item==NULL) { // truly the first item
          dcm->item=item;
        } else { // search for the previous one
          DCMITEM *ip; ip=dcm->item;
          while(ip->next_item!=NULL) ip=ip->next_item;
          ip->next_item=item; item->prev_item=ip;
        }
      } else {
        prev_item->next_item=item; item->prev_item=prev_item;
      }
    } else { // add as child
      parent_item->child_item=item; item->parent_item=parent_item;
    }
    /* Done. Do no free item or buf! */
    return(0);

  } else { // VL=0xFFFFFFFF
    size_t s=dcmVRVLength(item->vr);
    if(s==0) {
      if(verbose>0) printf(" Unknown VL!!\n");
      free(item);
      return(3); //return(TPCERROR_OK);
    }
    if(verbose>4) printf(" VR_based_VL=%u\n", (unsigned int)s);
    char *buf=(char*)calloc(s+1, sizeof(char));
    if(buf==NULL) {free(item); return(4);}
    if(fread(buf, 1, s, item->fp)!=s) {
      free(item); free(buf); return(3);
    }
    item->rd=buf;
    /* Add to list */
    if(sub==0) {
      if(prev_item==NULL) {
        if(dcm->item==NULL) { // truly the first item
          dcm->item=item;
        } else { // search for the previous one
          DCMITEM *ip; ip=dcm->item;
          while(ip->next_item!=NULL) ip=ip->next_item;
          ip->next_item=item; item->prev_item=ip;
        }
      } else {
        prev_item->next_item=item; item->prev_item=prev_item;
      }
    } else { // add as child
      parent_item->child_item=item; item->parent_item=parent_item;
    }
    /* Done. Do no free item or buf! */
    return(0);
  }

  return(0);
}
/*****************************************************************************/

/*****************************************************************************/
/** Read a single DICOM file.
    @sa dcmVerifyMagic, dcmfileInit, dcmfileFree, dcmReadTransferSyntaxUID, dcmFileWrite
    @return 0 when successful.
 */
int dcmFileRead(
  /** Pointer to filename. */
  const char *filename,
  /** Pointer to initiated data structure. */
  DCMFILE *dcm,
  /** Read only header (1), or read both header and pixel data (0). */
  const short int headerOnly, 
  /** Verbose level; if zero, then nothing is printed into stdout or stderr */
  int verbose
) {
  if(filename==NULL || strnlen(filename, 10)<1 || dcm==NULL) return(1);
  if(verbose>1) printf("%s('%s', %d)\n", __func__, filename, headerOnly);

  /* Delete any previous data */
  dcmfileFree(dcm);

  /* Open the file */
  strlcpy(dcm->filename, filename, FILENAME_MAX);
  dcm->fp=fopen(dcm->filename, "rb");
  if(dcm->fp==NULL) {
    return(2);
  }

  /* Check the magic number and move file pointer to the end of it */
  if(verbose>2) printf("checking DICOM magic number\n");
  if(dcmVerifyMagic(NULL, dcm->fp)!=1) {
    fclose(dcm->fp);
    return(2);
  }

  /* Get the Transfer Syntax UID */
  if(verbose>2) printf("checking Transfer Syntax UID\n");
  dcm->truid=dcmReadTransferSyntaxUID(dcm->fp);
  if(dcm->truid==DCM_TRUID_INVALID) { // not found
    fclose(dcm->fp);
    return(2);
  }
  if(verbose>0) { // print the UID
    printf("Transfer Syntax UID := %s\n", dcmTrUIDDescr(dcm->truid));
    fflush(stdout);
  }

  /* Check whether we currently support the Transfer UID */
  if(dcm->truid!=DCM_TRUID_LEE) {
    fclose(dcm->fp);
    return(2);
  }

  /* Read DICOM file elements */
  int ret=0;
  do {
    // note that the next function may need to call itself,
    // therefore counting loops here would not be useful.
    ret=dcmFileReadNextElement(dcm, NULL, NULL, 0, headerOnly, verbose-10);
  } while(ret==0 && !feof(dcm->fp));
  fclose(dcm->fp);
  /* TPCERROR_NO_KEY means that no (more) tag was found;
     other codes still mean that something bad happened. */
  if(ret==-1) {
    if(verbose>1) printf(" eof\n");
    ret=0;
  }
  return(ret);
}
/*****************************************************************************/

/*****************************************************************************/
/** Write a single DICOM file. 
    @sa dcmFileRead, dcmfileInit, dcmfileFree
    @return 0 when successful.
 */
int dcmFileWrite(
  /** Pointer to file name. */
  const char *filename,
  /** Pointer to DICOM data to be written. */
  DCMFILE *dcm,
  /** Verbose level; if zero, then nothing is printed into stdout or stderr */
  int verbose
) {
  if(filename==NULL || strnlen(filename, 10)<1 || dcm==NULL) return(1);
  if(verbose>1) {printf("%s('%s')\n", __func__, filename); fflush(stdout);}

  /* Check for the data */
  if(dcm->item==NULL) {
    return(2);
  }

  /* Check whether we currently support the Transfer UID */
  if(dcm->truid!=DCM_TRUID_LEE) { // Little endian explicit
    return(2);
  }


  /* Open the file */
  if(verbose>1) printf("opening the file for writing\n");
  FILE *fp;
  fp=fopen(filename, "wb");
  if(fp==NULL) return(3);

  /* Write preamble (just 128 zeroes) and magic number */
  {
    if(verbose>1) printf("writing preamble\n");
    char buf1[128], buf2[5]; 
    for(int i=0; i<128; i++) buf1[i]=(char)0;
    strcpy(buf2, "DICM");
    if(fwrite(buf1, 128, 1, fp)<1 || fwrite(buf2, 4, 1, fp)<1) {
      fclose(fp);
      return(3);
    }
  }


  /* Write the contents */
  if(verbose>1) printf("writing DICOM contents\n");
  int ret=0;
  DCMITEM *iptr;
  DCMITEM *d1=dcm->item;
  while(d1!=NULL) {
    if(verbose>2) {dcmitemPrint(d1);}
    /* Write */
    iptr=d1;
    {
      size_t n;
      /* Write tag */
      if(dcmWriteFileTag(fp, &iptr->tag)!=0) {ret=1; break;}
      /* Write VR and VL */
      if(dcmWriteFileVRVL(fp, iptr->vr, iptr->vl, NULL)!=0) {ret=2; break;}
      /* Write value, unless zero length, or SQ, in which case written later */
      if(iptr->vl>0 && iptr->vr!=DCM_VR_SQ) {
        size_t len;
        if(iptr->vl==0xFFFFFFFF) {
          len=dcmVRVLength(iptr->vr);
          if(verbose>30) printf("  value_len1 := %u\n", (unsigned int)len);
          n=fwrite(iptr->rd, dcmVRVLength(iptr->vr), 1, fp);
        } else {
          len=iptr->vl;
          if(verbose>30) printf("  value_len3 := %u\n", (unsigned int)len);
          n=fwrite(iptr->rd, iptr->vl, 1, fp);
        }
        if(verbose>30) printf("  value_len := %u\n", (unsigned int)len);
        if(n!=1) {ret=4; break;}
      } else if(iptr->vr==DCM_VR_SQ && d1->child_item==NULL) {
        if(verbose>1) printf("SQ, but no contents to write!\n");
        /* Write Sequence Delimitation Item */
        if(dcmWriteFileSQDelimItem(fp)!=0) {ret=6; break;}
      }
    }

    /* If this element (SQ) has children, then write those */
    /* Data Elements with a group of 0000, 0002 and 0006 shall not be present within Sequence Items,
       but that is not verified here */
    if(d1->child_item!=NULL) {
      DCMITEM *d2=d1->child_item;
      unsigned int d2counter=0;
      while(d2!=NULL) {
        if(verbose>2) {printf("  "); dcmitemPrint(d2);}

        /* Write */
        iptr=d2;

        /* First, write Item tag (FFFE,E000) */
        if(d2counter==0) {
          DCMTAG tag; tag.group=0xFFFE; tag.element=0xE000;
          if(dcmWriteFileTag(fp, &tag)!=0) {ret=11; break;}
        }
        /* Write item length; write 0xFFFFFFFF for now, correct later when known */
        fpos_t d2ilpos; // position for item length
        unsigned int d2il=0; // item length
        if(d2counter==0) {
          if(fgetpos(fp, &d2ilpos)) {ret=12; break;} // save position for writing later
          unsigned int ibuf;
          ibuf=0xFFFFFFFF;
          if(fwrite(&ibuf, 4, 1, fp)!=1) {ret=13; break;}
        }
        d2counter++;

        /* Write item value data set */
        {
          /* Write tag */
          if(dcmWriteFileTag(fp, &iptr->tag)!=0) {ret=14; break;}
          d2il+=4;
          /* Write VR and VL */
          unsigned int s;
          if(dcmWriteFileVRVL(fp, iptr->vr, iptr->vl, &s)!=0) {ret=15; break;}
          d2il+=s;
          /* Write value, unless zero length, or SQ, in which case written later */
          if(iptr->vl>0 && iptr->vr!=DCM_VR_SQ) {
            size_t len, n;
            if(iptr->vl==0xFFFFFFFF) {
              len=dcmVRVLength(iptr->vr);
              if(verbose>30) printf("  value_len1 := %u\n", (unsigned int)len);
              n=fwrite(iptr->rd, dcmVRVLength(iptr->vr), 1, fp);
              d2il+=len;
            } else {
              len=iptr->vl;
              if(verbose>30) printf("  value_len3 := %u\n", (unsigned int)len);
              n=fwrite(iptr->rd, iptr->vl, 1, fp);
              d2il+=iptr->vl;
            }
            if(verbose>30) printf("  value_len := %u\n", (unsigned int)len);
            if(n!=1) {ret=17; break;}
          } else if(iptr->vr==DCM_VR_SQ && iptr->child_item==NULL) {
            if(verbose>1) printf("SQ, but no contents to write!\n");
            /* Write Sequence Delimitation Item */
            if(dcmWriteFileSQDelimItem(fp)!=0) {ret=19; break;}
          }
        }

        /* If this element has children, then write those */
        if(d2->child_item!=NULL) {
          DCMITEM *d3=d2->child_item;
          unsigned int d3counter=0;
          while(d3!=NULL) {
            if(verbose>2) {printf("    "); dcmitemPrint(d3);}

            /* Write */
            iptr=d3;
            if(iptr->vr==DCM_VR_SQ) {d3=d3->next_item; continue;} // for now do not write SQs

            /* First, write Item tag (FFFE,E000) */
            if(d3counter==0) {
              DCMTAG tag; tag.group=0xFFFE; tag.element=0xE000;
              if(dcmWriteFileTag(fp, &tag)!=0) {ret=31; break;}
              d2il+=4;
            }
            /* Write item length; write 0xFFFFFFFF for now, correct later when known */
            fpos_t d3ilpos; // position for item length
            unsigned int d3il=0; // item length
            if(d3counter==0) {
              unsigned int ibuf;
              if(fgetpos(fp, &d3ilpos)) {ret=32; break;} // save position for writing later
              ibuf=0xFFFFFFFF;
              if(fwrite(&ibuf, 4, 1, fp)!=1) {ret=33; break;}
              d2il+=4;
            }
            d3counter++;

            /* Write item value data set */
            {
              /* Write tag */
              if(dcmWriteFileTag(fp, &iptr->tag)!=0) {ret=34; break;}
              d3il+=4; d2il+=4;
              /* Write VR and VL */
              unsigned int s;
              if(dcmWriteFileVRVL(fp, iptr->vr, iptr->vl, &s)!=0) {ret=35; break;}
              d3il+=s; d2il+=s;
              /* Write value, unless zero length, or SQ, in which case written later */
              if(iptr->vl>0 && iptr->vr!=DCM_VR_SQ) {
                size_t len, n;
                if(iptr->vl==0xFFFFFFFF) {
                  len=dcmVRVLength(iptr->vr);
                  if(verbose>30) printf("  value_len1 := %u\n", (unsigned int)len);
                  n=fwrite(iptr->rd, dcmVRVLength(iptr->vr), 1, fp);
                  d3il+=len; d2il+=len;
                } else {
                  len=iptr->vl;
                  if(verbose>30) printf("  value_len3 := %u\n", (unsigned int)len);
                  n=fwrite(iptr->rd, iptr->vl, 1, fp);
                  d3il+=iptr->vl; d2il+=iptr->vl;
                }
                if(verbose>30) printf("  value_len := %u\n", (unsigned int)len);
                if(n!=1) {ret=37; break;}
              } else if(iptr->vr==DCM_VR_SQ && iptr->child_item==NULL) {
                if(verbose>1) printf("SQ, but no contents to write!\n");
                /* Write Sequence Delimitation Item */
                if(dcmWriteFileSQDelimItem(fp)!=0) {ret=39; break;}
                d2il+=8;
              }
            }

            /* If this element has children, then write those */
            if(d3->child_item!=NULL) {
              //DCMITEM *d4=d3->child_item;
              if(verbose>0) fprintf(stderr, "Warning: 4th level items not written.\n");
            }

            /* now that we known the length, write it to the saved position */
            if(0) {
              fpos_t opos; // current position to return to
              if(fgetpos(fp, &opos)) {ret=40; break;}
              fsetpos(fp, &d3ilpos); // go to the position for item length
              char buf[4];
              memcpy(buf, &d3il, 4);
              if(!little_endian()) swabip(buf, 4);
              if(fwrite(&buf, 4, 1, fp)!=1) {ret=41; break;}
              fsetpos(fp, &opos); // back to the position where we were
            }

            d3=d3->next_item;
          }
          if(ret!=0) break;

          /* Write Item Delimitation Tag */
          if(dcmWriteFileSQItemDelimTag(fp)!=0) {ret=21; break;}
          /* End of the sequence - write Sequence Delimitation Item */
          if(dcmWriteFileSQDelimItem(fp)!=0) {ret=21; break;}

        }


        /* now that we known the length, write it to the saved position */
        if(0) {
          fpos_t opos; // current position to return to
          if(fgetpos(fp, &opos)) {ret=18; break;}
          fsetpos(fp, &d2ilpos); // go to the position for item length
          char buf[4];
          memcpy(buf, &d2il, 4);
          if(!little_endian()) swabip(buf, 4);
          if(fwrite(&buf, 4, 1, fp)!=1) {ret=19; break;}
          fsetpos(fp, &opos); // back to the position where we were
        }

        d2=d2->next_item;
      }
      if(ret!=0) break;

      /* Write Item Delimitation Tag */
      if(dcmWriteFileSQItemDelimTag(fp)!=0) {ret=21; break;}
      /* End of the sequence - write Sequence Delimitation Item */
      if(dcmWriteFileSQDelimItem(fp)!=0) {ret=21; break;}

    }

    d1=d1->next_item;
  } // next
  if(ret!=0) {
    if(verbose>0) fprintf(stderr, "  ret := %d\n", ret);
    fclose(fp);
    return(3);
  }


  fclose(fp);

  return(0);
}
/*****************************************************************************/

/*****************************************************************************/
