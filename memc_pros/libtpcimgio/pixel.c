/// @file pixel.c
/// @author Vesa Oikonen
/// @brief Functions for reading and writing pixel definition files.
///
/*****************************************************************************/
#include "libtpcimgio.h"
/*****************************************************************************/

/*****************************************************************************/
/** Initiate the IMG_PIXELS struct before any use.
    @sa pxlAllocate, pxlFree
 *  @author Vesa Oikonen
 */
void pxlInit(
  /** Pointer to IMG_PIXELS */
  IMG_PIXELS *pxl
) {
  if(pxl==NULL) return;
  pxl->pxlNr=pxl->_pxlNr=0;
  pxl->p=NULL;
}
/*****************************************************************************/

/*****************************************************************************/
/** Free memory allocated for IMG_PIXELS. All data is cleared. 
    @sa pxlInit, pxlAllocate
 */
void pxlFree(
  /** Pointer to IMG_PIXELS struct */
  IMG_PIXELS *pxl
) {
  if(pxl==NULL) return;
  free(pxl->p);
  pxlInit(pxl);
}
/*****************************************************************************/

/*****************************************************************************/
/** Allocate memory for IMG_PIXELS data.
 *  Any previous contents are deleted.
 *  Return Returns 0 when successful.
 *  @sa pxlInit, pxlAllocateMore, pxlFree
 */
int pxlAllocate(
  /** Pointer to initiated IMG_PIXELS struct data; any old contents are deleted.
   *  pxlNr inside the struct is set to or kept at zero. */
  IMG_PIXELS *pxl,
  /** Nr of pixels to allocate */
  int pxlNr
) {
  if(pxl==NULL) return(1);
  /* Delete any previous contents */
  pxlFree(pxl);
  /* If no memory is requested, then just return */
  if(pxlNr<1) return(0);

  /* Allocate memory for IMG_PIXEL data */
  pxl->p=(IMG_PIXEL*)malloc(pxlNr*sizeof(IMG_PIXEL));
  if(pxl->p==NULL) return(2);
  for(int i=0; i<pxlNr; i++) {
    pxl->p[i].x=pxl->p[i].y=pxl->p[i].z=pxl->p[i].f=0;
  }
  pxl->_pxlNr=pxlNr;
  return(0);
}
/*****************************************************************************/

/*****************************************************************************/
/** Allocate memory for more IMG_PIXELS data.
 *  Any previous contents are preserved.
 *  Return Returns 0 when successful.
 *  @sa pxlInit, pxlFree, pxlAllocate, pxlMakeRoom
 */
int pxlAllocateMore(
  /** Pointer to initiated IMG_PIXELS struct data; 
   *  any old contents are preserved, but existing data is not required.
   *  pxlNr inside the struct is set to or kept at zero. */
  IMG_PIXELS *pxl,
  /** Nr of additional pixels to allocate; if struct contains unused space 
      for requested pixels already, then nothing is done. */
  int pxlNr
) {
  if(pxl==NULL) return(1);
  /* If no memory is requested, then just return */
  if(pxlNr<1) return(0);
  /* If none allocated previously then use pxlAllocate */
  if(pxl->_pxlNr==0) return(pxlAllocate(pxl, pxlNr));
  /* Check if there is enough space already */
  int newPxlNr, addPxlNr;
  newPxlNr=pxl->pxlNr+pxlNr; addPxlNr=newPxlNr-pxl->_pxlNr;
  if(addPxlNr<=0) return(0);
  /* Reallocate */
  IMG_PIXEL *pxlPtr;
  pxlPtr=realloc(pxl->p, sizeof(IMG_PIXEL)*newPxlNr);
  if(pxlPtr==NULL) return(2);
  pxl->p=pxlPtr;
  for(int i=pxl->_pxlNr; i<newPxlNr; i++)
    pxl->p[i].x=pxl->p[i].y=pxl->p[i].z=pxl->p[i].f=0;
  pxl->_pxlNr=newPxlNr;
  return(0);
}
/*****************************************************************************/

/*****************************************************************************/
/** Make room for new pixels in the IMG_PIXELS list, allocating more memory
 *  if needed. Previous contents are preserved but moved in the list.
 *  @sa pxlAllocateMore
 *  @return 0 if successful.
 */
int pxlMakeRoom(
  /** Pointer to IMG_PIXELS struct with existing contents. */
  IMG_PIXELS *list,
  /** Index [0..pxlNr] of the new room start position. */
  int i,
  /** Nr of empty list items to add. */
  int n
) {
  if(list==NULL || i<0 || i>list->pxlNr) return(1);
  /* Check whether anything needs to be added */
  if(n<1) return(0);
  /* If user wanted room in the end, then just add the memory */
  if(i==list->pxlNr) return(pxlAllocateMore(list, n));
  /* Otherwise, first add the space */
  if(pxlAllocateMore(list, n)!=0) return(2);
  /* and the move previous data forward in the list */
  for(int li=list->pxlNr-1; li>=i; li--) list->p[li+n]=list->p[li];
  /* and add pxlNr */
  list->pxlNr+=n;
  return(0);
}
/*****************************************************************************/

/*****************************************************************************/
/** Add given pixel into IMG_PIXELS data.
 *  @return 0 if successful.
 *  @author Vesa Oikonen
 *  @sa pxlInit, pxlFree, pxlAllocate, pxlMakeRoom, pxlRm, pxlRead, pxlWrite
 */
int pxlAdd(
  /** Pointer to IMG_PIXELS struct, which must be initiated. Memory is
      added if needed, and pxlNr increased. */
  IMG_PIXELS *list,
  /** Pointer to IMG_PIXEL struct to add. */
  IMG_PIXEL *pxl
) {
  if(list==NULL || pxl==NULL) return(1);
  int ret;
  if((ret=pxlAllocateMore(list, 1))!=0) return(ret);
  list->p[list->pxlNr].x=pxl->x;
  list->p[list->pxlNr].y=pxl->y;
  list->p[list->pxlNr].z=pxl->z;
  list->p[list->pxlNr].f=pxl->f;
  list->pxlNr++;
  return(0);
}
/*****************************************************************************/

/*****************************************************************************/
/** Get a pixel from IMG_PIXELS list.
 *  @sa pxlInit, pxlFree, pxlAllocate, pxlMakeRoom, pxlAdd, pxlRead
 *  @return 0 if successful.
 */
int pxlGet(
  /** Pointer to IMG_PIXELS struct, containing the list of pixels. */
  IMG_PIXELS *list,
  /** Pixel list index [0..list->pxlNr-1]. */
  int i,
  /** Pointer to IMG_PIXEL struct into which pixel coordinates are written. */
  IMG_PIXEL *pxl
) {
  if(list==NULL || pxl==NULL || i<0 || i>=list->pxlNr) return(1);
  pxl->x=list->p[i].x;
  pxl->y=list->p[i].y;
  pxl->z=list->p[i].z;
  return(0);
}
/*****************************************************************************/

/*****************************************************************************/
/** Add pixel(s) from mask image into IMG_PIXELS data.
 *  @return the nr of added pixels.
 *  @sa pxlInit, pxlFree, pxlAdd
 *  @author Vesa Oikonen
 */
int pxlAddFromMask(
  /** Pointer to IMG_PIXELS struct, which must be initiated. Memory is
      added if needed, and pxlNr increased. */
  IMG_PIXELS *list,
  /** Pointer to mask image. */
  IMG *img
) {
  if(list==NULL || img==NULL) return(1);
  if(img->dimz<1 || img->dimy<1 || img->dimx<1 || img->dimt<1) return(0);
  int zi, yi, xi, n=0;
  for(zi=0; zi<img->dimz; zi++)
    for(yi=0; yi<img->dimy; yi++)
      for(xi=0; xi<img->dimx; xi++)
        if(fabs(img->m[zi][yi][xi][0])>=0.5) n++;
  if(n==0) return(0);
  if(pxlAllocateMore(list, n)!=0) return(0);
  for(zi=0; zi<img->dimz; zi++)
    for(yi=0; yi<img->dimy; yi++)
      for(xi=0; xi<img->dimx; xi++)
        if(fabs(img->m[zi][yi][xi][0])>=0.5) {
          list->p[list->pxlNr].x=1+xi;
          list->p[list->pxlNr].y=1+yi;
          list->p[list->pxlNr].z=1+zi;
          list->p[list->pxlNr].f=0;
          list->pxlNr++;
        }
  return(n);
}
/*****************************************************************************/

/*****************************************************************************/
/** Move pixel from one slot to another inside IMG_PIXELS data,
 *  changing the position of others accordingly.
 *  @sa pxlInit, pxlFree, pxlAllocate, pxlRm, pxlMakeRoom
 */
void pxlMove(
  /** Pointer to IMG_PIXELS struct, which must be initiated. */
  IMG_PIXELS *list,
  /** Index [0.._pxlNr-1] of source position. */
  int from,
  /** Index [0.._pxlNr-1] of target position. */
  int to
) {
  if(list==NULL || from<0 || to<0) return;
  if(from>=list->_pxlNr || to>=list->_pxlNr) return;
  if(from==to) return;
  int i=from;
  IMG_PIXEL tmp=list->p[from];
  if(from>to) {
    for(int i=from; i>to; i--) list->p[i]=list->p[i-1];
  } else {
    for(int i=from; i<to; i++) list->p[i]=list->p[i+1];
  }
  list->p[i]=tmp;
  return;
}
/*****************************************************************************/

/*****************************************************************************/
/** Remove specified pixel from IMG_PIXELS data.
 *  @return 0 if successful.
 *  @author Vesa Oikonen
 *  @sa pxlInit, pxlFree, pxlRmDuplicates, pxlAdd, pxlAllocate, pxlMakeRoom, pxlGet
 */
int pxlRm(
  /** Pointer to IMG_PIXELS struct, which must be initiated. 
      Allocated memory is not reduced, but pxlNr is decreased. */
  IMG_PIXELS *list,
  /** Index [0..pxlNr-1] of pixel to delete. */
  int index
) {
  if(list==NULL || index<0) return(1);
  if(index>=list->pxlNr) return(0);
  /* If last one, then just decrease the pxlNr */
  if(index==list->pxlNr-1) {list->pxlNr--; return(0);}
  /* Otherwise move it to the last place and then decrease the pxlNr */
  pxlMove(list, index, list->pxlNr-1);
  list->pxlNr--; 
  return(0);
}
/*****************************************************************************/

/*****************************************************************************/
/** Remove duplicates from IMG_PIXELS data.
 *  @return the nr of removed pixels.
 *  @author Vesa Oikonen
 *  @sa pxlInit, pxlFree, pxlRm, pxlAdd, pxlAllocate, pxlMakeRoom, pxlGet
 */
int pxlRmDuplicates(
  /** Pointer to IMG_PIXELS struct, which must be initiated. 
      Allocated memory is not reduced, but pxlNr is decreased. */
  IMG_PIXELS *list
) {
  if(list==NULL || list->pxlNr<2) return(0);
  int i=list->pxlNr-1, j, n=0;
  while(i>0) {
    for(j=0; j<i; j++) {
      if(list->p[i].z!=list->p[j].z) continue;
      if(list->p[i].x!=list->p[j].x) continue;
      if(list->p[i].y!=list->p[j].y) continue;
      pxlRm(list, i); n++; break;
    }
    i--;
  }
  return(n);
}
/*****************************************************************************/

/*****************************************************************************/
/** Write IMG_PIXELS data into specified file.
 *  @return 0 if successful.
 *  @author Vesa Oikonen
 *  @sa pxlRead, pxlFree, pxlGet, pxlAdd
 */
int pxlWrite(
  /** Pointer to IMG_PIXELS struct, contents of which are to be written */
  IMG_PIXELS *pxl,
  /** Output file pointer */
  FILE *fp,
  /** Pointer to a string (allocated for at least 64 chars) where error message
      or other execution status will be written; enter NULL, if not needed */     
  char *status
) {
  if(pxl==NULL || pxl->pxlNr<1 || pxl->p==NULL) {
    if(status!=NULL) strcpy(status, "no pixels to write");
    return(1);
  }
  int i, n=7;
  for(i=0; i<pxl->pxlNr && n>6; i++) 
    n=fprintf(fp, "%d,%d,%d,%d\n", 
              pxl->p[i].x, pxl->p[i].y, pxl->p[i].z, pxl->p[i].f);
  if(n<7) {
    if(status!=NULL) strcpy(status, "cannot write pixels into file");
    return(2);
  }
  if(status!=NULL) strcpy(status, "ok");
  return(0);
}
/*****************************************************************************/

/*****************************************************************************/
/** Read IMG_PIXELS data from specified file.
 *  @sa pxlInit, pxlFree, pxlWrite
 *  @return 0 if successful.
 *  @author Vesa Oikonen
 */
int pxlRead(
  /** Pointer to IMG_PIXELS struct, into which contents of file are to be added;
      call pxlInit() once before using this function. */
  IMG_PIXELS *pxl,
  /** Pointer to the file name; this string is not modified. */
  const char *fname,
  /** Pointer to a string (allocated for at least 64 chars) where error message
      or other execution status will be written; enter NULL, if not needed */     
  char *status
) {
  if(pxl==NULL) {
    if(status!=NULL) strcpy(status, "program error");
    return(1);
  }
  int i, c, n, longest, ret;
 
  /* Try to read the file */
  FILE *fp;
  fp=fopen(fname, "r");
  if(fp==NULL) {
    if(status!=NULL) strcpy(status, "cannot open file");
    return(2);
  }
  /* Get the length of the longest line */
  i=longest=0; 
  while((c=fgetc(fp))!=EOF) {
    if(c==10 || c==13) {if(i>longest) longest=i; i=0;} else i++;
  }
  if(i>longest) longest=i;
  rewind(fp); longest+=1;
  /* and allocate memory for string of that length */
  char *line, buf[20];
  line=(char*)malloc((longest+1)*sizeof(char));
  if(line==NULL) {
    if(status!=NULL) strcpy(status, "out of memory");
    fclose(fp); return(3);
  }
  /* Allocate space for a few pixels */
  if(pxlAllocateMore(pxl, 10)!=0) {
    if(status!=NULL) strcpy(status, "out of memory");
    fclose(fp); free(line); return(3);    
  }
  /* Read data lines */
  while(fgets(line, longest, fp)!=NULL) {
    /* forget comment lines */
    if(line[0]=='#') continue;
    /* get nr of tokens on this line */
    n=strTokenNr(line, " ,;\t\n\r"); if(n==0) continue;
    if(n<3 || n>4) { 
      if(status!=NULL) strcpy(status, "invalid format");
      fclose(fp); free(line); return(4);
    }
    /* Allocate more memory if necessary */
    if(pxl->pxlNr==pxl->_pxlNr) {
      if(pxlAllocateMore(pxl, 10)!=0) {
        if(status!=NULL) strcpy(status, "out of memory");
        fclose(fp); free(line); return(3);    
      }
    }
    /* Read the pixel coordinates */
    ret=0; //printf("line := %s", line); printf("n := %d\n", n);
    for(i=0; i<n; i++) {
      if(strTokenNCpy(line, " ,;\t\n\r", 1+i, buf, 20)<1) {ret=1; break;}
      //printf("  buf := %s\n", buf);
      if((ret=atoi_with_check(buf, &c))!=0) break;
      if(i==0) pxl->p[pxl->pxlNr].x=c;
      else if(i==1) pxl->p[pxl->pxlNr].y=c;
      else if(i==2) pxl->p[pxl->pxlNr].z=c;
      else if(i==3) pxl->p[pxl->pxlNr].f=c;
    }  
    if(ret) {
      if(status!=NULL) strcpy(status, "invalid coordinate");
      fclose(fp); free(line); return(4);
    }
    pxl->pxlNr++;
  }
  fclose(fp); free(line);
  return(0);
}
/*****************************************************************************/

/*****************************************************************************/
