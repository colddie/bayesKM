/// @file ird.c
/// @author Vesa Oikonen
/// @brief Storing and processing of 4D image coordinate data.
///
/*****************************************************************************/
#include "libtpcimgio.h"
/*****************************************************************************/

/*****************************************************************************/
/** Read voxel coordinates including time frame from a string representation.
    @sa string_to_xyz, irdRead
    @return Returns 0 if successful, >0 if not.
 */
int string_to_xyzf(
  /** String in format x,y,z,f or x y z f; frame (f) is optional; if f is not
      specified, then 0 is written in its place. */
  const char *str,
  /** Pointer to image pixel struct; obligatory. */
  IMG_PIXEL *v
) {
  if(v==NULL) return 1;
  v->x=v->y=v->z=v->f=0;

  char *cptr, tmp[256];
  strncpy(tmp, str, 255); tmp[255]=(char)0;
  cptr=strtok(tmp, " ,;:()|-"); if(cptr==NULL) return 1;
  v->x=atoi(cptr); if(v->x<1) return 1;  
  cptr=strtok(NULL, " ,;:()|-"); if(cptr==NULL) return 2;
  v->y=atoi(cptr); if(v->y<1) return 1;
  cptr=strtok(NULL, " ,;:()|-"); if(cptr==NULL) return 3;
  v->z=atoi(cptr); if(v->z<1) return 1;
  cptr=strtok(NULL, " ,;:()|-"); if(cptr==NULL) return 0;
  v->f=atoi(cptr);
  return 0;
}
/*****************************************************************************/

/*****************************************************************************/
/** Reorder Image Range Definition.
    @sa vrdReorder
    @details Function name was previously ifrReorder
    @return Returns 0 if successful.
 */
int irdReorder(
  /** Image volume range; start and end range are set in correct order */
  IMG_RANGE *img_range
) {
  int i;

  /* Check that input is ok */
  if(img_range==NULL) return 1;
  /* Change the order if necessary */
  if(img_range->x1<0 || img_range->x2<0) return 2;
  if(img_range->x2<img_range->x1) {
    i=img_range->x1; img_range->x1=img_range->x2; img_range->x2=i;}
  if(img_range->y1<0 || img_range->y2<0) return 3;
  if(img_range->y2<img_range->y1) {
    i=img_range->y1; img_range->y1=img_range->y2; img_range->y2=i;}
  if(img_range->z1<0 || img_range->z2<0) return 4;
  if(img_range->z2<img_range->z1) {
    i=img_range->z1; img_range->z1=img_range->z2; img_range->z2=i;}
  if(img_range->f1<0 || img_range->f2<0) return 5;
  if(img_range->f2<img_range->f1) {
    i=img_range->f1; img_range->f1=img_range->f2; img_range->f2=i;}

  return 0;
}
/*****************************************************************************/

/*****************************************************************************/
/** Read Image Range Definition File.
    @return Returns 0 if successful.
    @sa irdCheck, string_to_xyzf, irdReorder
 */
int irdRead(
  /** Image Range Definition File filename, which contains 4D image volume
      corners (x y z f) in IFT format: corner1 = x y z f corner2 = x y z f;
      If time frames (f) are missing, then frame coordinate is set to 0.
   */
  char *irdfile,
  /** Image volume range */
  IMG_RANGE *img_range,
  /** Pointer to a string (allocated for at least 64 chars) where error message
      or other execution status will be written; enter NULL, if not needed */     
  char *status
) {
  int ret, ii;
  IFT ift;
  char key[256];
  IMG_PIXEL v;

  /* Check that input is ok */
  if(irdfile==NULL || strnlen(irdfile, 2)<1 || img_range==NULL) {
    if(status!=NULL) strcpy(status, "program error");
    return 1;
  }
  /* Read IDF as IFT file */
  iftInit(&ift); ret=iftRead(&ift, irdfile, 1); if(ret) {
    if(status!=NULL) strcpy(status, ift.status);
    iftEmpty(&ift); return 2;
  }
  /* Try to find keys 'corner1' and 'corner2' */
  strcpy(key, "corner1"); ii=iftGet(&ift, key);
  if(ii>=0) {
    ret=string_to_xyzf(ift.item[ii].value, &v);
    if(ret==0) {
      img_range->x1=v.x; img_range->y1=v.y; img_range->z1=v.z; img_range->f1=v.f;
      strcpy(key, "corner2"); ii=iftGet(&ift, key);
      if(ii>=0) {
        ret=string_to_xyzf(ift.item[ii].value, &v);
        img_range->x2=v.x; img_range->y2=v.y; img_range->z2=v.z; img_range->f2=v.f;
        if(ret==0) {
          irdReorder(img_range);
          if(status!=NULL) strcpy(status, "ok");
          iftEmpty(&ift); return 0;
        }
      }
    }
  }
  /* We are here only if keys were not found */
  /* Lets not care about keys at all */
  for(ii=0, ret=0; ii<ift.keyNr; ii++) {
    if(ret==0 && string_to_xyzf(ift.item[ii].value, &v)==0)
    {
      img_range->x1=v.x; img_range->y1=v.y; img_range->z1=v.z; img_range->f1=v.f;
      ret++; continue;
    }
    if(ret==1 && string_to_xyzf(ift.item[ii].value, &v)==0)
    {
      img_range->x2=v.x; img_range->y2=v.y; img_range->z2=v.z; img_range->f2=v.f;
      ret++; break;
    }
  }
  if(ret<2) {
    if(status!=NULL) strcpy(status, "volume definitions not found");
    iftEmpty(&ift); return 2;
  }

  irdReorder(img_range);
  if(status!=NULL) strcpy(status, "ok");
  iftEmpty(&ift);
  return 0;
}
/*****************************************************************************/

/*****************************************************************************/
/** Check that image range definition is inside image data.
    If time frames are zero, then fix those to image frame range.
    @sa irdReorder, irdRead
    @return Returns 0 if successful.
 */
int irdCheck(
  /** Pointer to image volume range data. */
  IMG_RANGE *r,
  /** Pointer to image data. */
  IMG *img
) {
  if(r==NULL || img==NULL) return(1);
  if(img->dimx<1 || img->dimy<1 || img->dimz<1 || img->dimt<1) return(2);

  if(r->x1<1 || r->x1>img->dimx) return(11);
  if(r->x2<1 || r->x2>img->dimx) return(12);

  if(r->y1<1 || r->y1>img->dimy) return(21);
  if(r->y2<1 || r->y2>img->dimy) return(22);

  if(r->z1<1 || r->z1>img->dimz) return(31);
  if(r->z2<1 || r->z2>img->dimz) return(32);

  /* If time frame range is not set, then set it now */
  if(r->f1<1 && r->f2<1) {
    r->f1=1; r->f2=img->dimt;
    return(0);
  }
  /* else check as usual */
  if(r->f1<1 || r->f1>img->dimt) return(41);
  if(r->f2<1 || r->f2>img->dimt) return(42);
  return(0);
}
/*****************************************************************************/

/*****************************************************************************/
