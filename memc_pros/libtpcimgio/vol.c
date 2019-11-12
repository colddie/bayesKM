/// @file vol.c
/// @author Vesa Oikonen
/// @brief Storing and processing of 3D PET image volume data with no time
///        information (frames).
///
/*****************************************************************************/
#include "libtpcimgio.h"
/*****************************************************************************/
/*!
 *  Status (error) messages from volume processing
 */
char *_volStatusMessage[] = {
  /*  0 */ "ok",
  /*  1 */ "fault in calling routine",
  /*  2 */ "out of memory"
};
/*****************************************************************************/

/*****************************************************************************/
/*!
 * Initiate volume before any use of VOL data; this should be called once.
 *
 * @param vol pointer to volume structure
 */
void volInit(VOL *vol) {
  if(VOL_TEST) printf("volInit()\n");
  if(vol==NULL) return;
  memset(vol, 0, sizeof(VOL));
  vol->status=IMG_STATUS_INITIALIZED;
  vol->statmsg=_volStatusMessage[0];
  vol->orientation=0;
  vol->dimx=vol->dimy=vol->dimz=0;
  vol->sizex=vol->sizey=vol->sizez=0;
  vol->v=(float***)NULL;
  vol->voxel=(float*)NULL;
  vol->column=(float*)NULL;
  vol->row=(float**)NULL;
  vol->plane=(float***)NULL;
}
/*****************************************************************************/

/*****************************************************************************/
/*!
 * Initiate short int volume before any use of SVOL data;
 * this should be called once.
 * 
 * @param svol short int volume structure
 */
void svolInit(SVOL *svol) {
  if(VOL_TEST) printf("svolInit()\n");
  if(svol==NULL) return;
  memset(svol, 0, sizeof(SVOL));
  svol->status=IMG_STATUS_INITIALIZED;
  svol->statmsg=_volStatusMessage[0];
  svol->orientation=0;
  svol->dimx=svol->dimy=svol->dimz=0;
  svol->sizex=svol->sizey=svol->sizez=0;
  svol->scale_factor=1.0;
  svol->v=(short int***)NULL;
  svol->voxel=(short int*)NULL;
  svol->column=(short int*)NULL;
  svol->row=(short int**)NULL;
  svol->plane=(short int***)NULL;
}
/*****************************************************************************/

/*****************************************************************************/
/*!
 * Free memory allocated for volume.
 *
 * @param vol volume structure
 */
void volEmpty(VOL *vol) {
  if(VOL_TEST) printf("volEmpty()\n");
  if(vol==NULL || vol->status<IMG_STATUS_OCCUPIED) return;
  /* Free up memory */
  if(vol->_vxl!=NULL) free(vol->_vxl);
  //if(vol->_col!=NULL) free(vol->_col); Same as _vxl
  if(vol->_row!=NULL) free(vol->_row);
  if(vol->_pln!=NULL) free(vol->_pln);
  /* Set variables */
  vol->statmsg=_volStatusMessage[0];
  vol->orientation=0;
  vol->dimx=vol->dimy=vol->dimz=0;
  vol->sizex=vol->sizey=vol->sizez=0;
  vol->v=(float***)NULL;
  vol->voxel=(float*)NULL;
  vol->column=(float*)NULL;
  vol->row=(float**)NULL;
  vol->plane=(float***)NULL;
  /* Set status */
  vol->status=IMG_STATUS_INITIALIZED;
}
/*****************************************************************************/

/*****************************************************************************/
/*!
 * Free memory allocated for short int volume.
 *
 * @param svol short int volume structure
 */
void svolEmpty(SVOL *svol) {
  if(VOL_TEST) printf("svolEmpty()\n");
  if(svol==NULL || svol->status<IMG_STATUS_OCCUPIED) return;
  /* Free up memory */
  if(svol->_vxl!=NULL) free(svol->_vxl);
  //if(svol->_col!=NULL) free(svol->_col); Same as _vxl
  if(svol->_row!=NULL) free(svol->_row);
  if(svol->_pln!=NULL) free(svol->_pln);
  /* Set variables */
  svol->statmsg=_volStatusMessage[0];
  svol->orientation=0;
  svol->dimx=svol->dimy=svol->dimz=0;
  svol->sizex=svol->sizey=svol->sizez=0;
  svol->scale_factor=1.0;
  svol->v=(short int***)NULL;
  svol->voxel=(short int*)NULL;
  svol->column=(short int*)NULL;
  svol->row=(short int**)NULL;
  svol->plane=(short int***)NULL;
  /* Set status */
  svol->status=IMG_STATUS_INITIALIZED;
}
/*****************************************************************************/

/*****************************************************************************/
/*!
 * Allocate memory for 3D image volume. Returns 0 if ok.
 *
 * @param vol volume structure
 * @param planes number of planes [>=1]
 * @param rows number of rows [>=1]
 * @param columns number of columns [>=1]
 * @return 0 if ok, 1 invalid image status, 2 invalid input, 
 * 5 failed to allocate memory for planes,  6 failed to allocate memory for rows,
 * 8 failed to allocate memory for rows
 */
int volAllocate(VOL *vol, int planes, int rows, int columns) {
  unsigned short int zi, ri;
  int vxlNr, vi;
  float **rptr, *cptr;

  if(VOL_TEST) printf("voiAllocate(*vol, %d, %d, %d)\n", planes, rows, columns);
  /* Check arguments */
  if(vol==NULL) return(1); else vol->statmsg=_volStatusMessage[1];
  if(vol->status==IMG_STATUS_UNINITIALIZED) return(1);
  if(planes<1 || rows<1 || columns<1) return(2);
  vxlNr=planes*rows*columns;

  /* Check if correct volume size is already allocated */
  if(vol->status>=IMG_STATUS_OCCUPIED) {
    if(planes==vol->dimz && rows==vol->dimy && columns==vol->dimx) {
      for(vi=0; vi<vxlNr; vi++) vol->_vxl[vi]=0;
      return(0); /* that's it */
    } else {
      volEmpty(vol);
    }
  }
  /* Allocate memory for volume data */
  vol->_pln=(float***)malloc(planes*sizeof(float**));
  if(vol->_pln==NULL) {
    return(5);}
  vol->_row=(float**)malloc(planes*rows*sizeof(float*));
  if(vol->_row==NULL) {
    free(vol->_pln); return(6);}
  vol->_col=vol->_vxl=(float*)calloc(planes*rows*columns, sizeof(float));
  if(vol->_vxl==NULL) {
    free(vol->_pln); free(vol->_row); return(8);
  }
  /* Set data pointers */
  rptr=vol->_row; cptr=vol->_col;
  for(zi=0; zi<planes; zi++) {
    vol->_pln[zi]=rptr;
    for(ri=0; ri<rows; ri++) {
      *rptr++=cptr; cptr+=columns;
    }
  }
  vol->v=vol->_pln;
  vol->plane=vol->_pln;
  vol->column=vol->_col;
  vol->row=vol->_row;
  vol->voxel=vol->_vxl;
  /* Ok */
  vol->dimz=planes; vol->dimy=rows; vol->dimx=columns;
  vol->statmsg=_volStatusMessage[0];
  vol->status=IMG_STATUS_OCCUPIED;
  return(0);
}
/*****************************************************************************/

/*****************************************************************************/
/*!
 * Allocate memory for 3D short int volume. Returns 0 if ok.
 *
 * @param svol short volume structure
 * @param planes number of planes [>=1]
 * @param rows number of rows [>=1]
 * @param columns number of columns [>=1]
 * @return 0 if ok, 1 invalid image status, 2 invalid input, 
 * 5 failed to allocate memory for planes,  6 failed to allocate memory for rows,
 * 8 failed to allocate memory for rows
 */
int svolAllocate(SVOL *svol, int planes, int rows, int columns) {
  unsigned short int zi, ri;
  int vxlNr, vi;
  short int **rptr, *cptr;

  if(VOL_TEST) printf("svoiAllocate(*svol, %d, %d, %d)\n", planes, rows, columns);
  /* Check arguments */
  if(svol==NULL) return(1); else svol->statmsg=_volStatusMessage[1];
  if(svol->status==IMG_STATUS_UNINITIALIZED) return(1);
  if(planes<1 || rows<1 || columns<1) return(2);
  vxlNr=planes*rows*columns;

  /* Check if correct volume size is already allocated */
  if(svol->status>=IMG_STATUS_OCCUPIED) {
    if(planes==svol->dimz && rows==svol->dimy && columns==svol->dimx) {
      for(vi=0; vi<vxlNr; vi++) svol->_vxl[vi]=0;
      return(0); /* that's it */
    } else {
      svolEmpty(svol);
    }
  }
  /* Allocate memory for volume data */
  svol->_pln=(short int***)malloc(planes*sizeof(short int**));
  if(svol->_pln==NULL) {
    return(5);}
  svol->_row=(short int**)malloc(planes*rows*sizeof(short int*));
  if(svol->_row==NULL) {
    free(svol->_pln); return(6);}
  svol->_col=svol->_vxl=(short int*)calloc(planes*rows*columns, sizeof(short int));
  if(svol->_vxl==NULL) {
    free(svol->_pln); free(svol->_row); return(8);
  }
  /* Set data pointers */
  rptr=svol->_row; cptr=svol->_col;
  for(zi=0; zi<planes; zi++) {
    svol->_pln[zi]=rptr;
    for(ri=0; ri<rows; ri++) {
      *rptr++=cptr; cptr+=columns;
    }
  }
  svol->v=svol->_pln;
  svol->plane=svol->_pln;
  svol->column=svol->_col;
  svol->row=svol->_row;
  svol->voxel=svol->_vxl;
  /* Ok */
  svol->dimz=planes; svol->dimy=rows; svol->dimx=columns;
  svol->statmsg=_volStatusMessage[0];
  svol->status=IMG_STATUS_OCCUPIED;
  return(0);
}
/*****************************************************************************/

/*****************************************************************************/
/*!
 * Copy one time frame (1..dimt) from 4D image to 3D volume.
 *  Vol can be but need not to be allocated.
 *
 * @param img image structure
 * @param vol volume structure
 * @param frame frame number [1..number of frames]
 * @return 0 if ok, 1 invalid image status, 2 invalid input
 */
int img2vol(IMG *img, VOL *vol, int frame) {
  int ret;
  unsigned short int zi, yi, xi, fi;

  if(VOL_TEST) printf("img2vol(img, %d, vol)\n", frame);
  /* Check input */
  if(vol==NULL) return(1);
  vol->statmsg=_volStatusMessage[1];
  if(img->status!=IMG_STATUS_OCCUPIED) return(1);
  if(frame<1 || img->dimt<frame) return(2);
  if(vol->status==IMG_STATUS_UNINITIALIZED) return(1);

  /* Allocate memory (if needed) for volume */
  ret=volAllocate(vol, img->dimz, img->dimy, img->dimx);
  if(ret) return(ret);

  /* Copy data */
  fi=frame-1;
  vol->orientation=img->orientation;
  vol->sizex=img->sizex; vol->sizey=img->sizey; vol->sizez=img->sizez;
  for(zi=0; zi<vol->dimz; zi++)
    for(yi=0; yi<vol->dimy; yi++)
      for(xi=0; xi<vol->dimx; xi++)
        vol->v[zi][yi][xi]=img->m[zi][yi][xi][fi];

  return(0);
}
/*****************************************************************************/

/*****************************************************************************/
/** Copy one time frame (1..dimt) from 4D image to 3D short int volume.
 *  Svol can be but need not to be allocated.
 *
 * @param img image structure
 * @param svol short volume structure
 * @param frame frame number [1..number of frames]
 * @return 0 if ok, 1 invalid image status, 2 invalid input
 */
int img2svol(IMG *img, SVOL *svol, int frame) {
  int ret;
  unsigned short int zi, yi, xi, fi;
  float fmin, fmax, g;

  if(VOL_TEST) printf("img2svol(img, %d, svol)\n", frame);
  /* Check input */
  if(svol==NULL) return(1);
  svol->statmsg=_volStatusMessage[1];
  if(img->status!=IMG_STATUS_OCCUPIED) return(1);
  if(frame<1 || img->dimt<frame) return(2);
  if(svol->status==IMG_STATUS_UNINITIALIZED) return(1);

  /* Allocate memory (if needed) for volume */
  ret=svolAllocate(svol, img->dimz, img->dimy, img->dimx);
  if(ret) return(ret);

  /* Copy data */
  fi=frame-1;
  svol->orientation=img->orientation;
  svol->sizex=img->sizex; svol->sizey=img->sizey; svol->sizez=img->sizez;
  ret=imgFrameMinMax(img, frame, &fmin, &fmax); if(ret) return(10+ret);
  if(fabs(fmin)>fabs(fmax)) g=fabs(fmin); else g=fabs(fmax);
  if(g!=0) g=32766./g; else g=1.0;
  for(zi=0; zi<svol->dimz; zi++)
    for(yi=0; yi<svol->dimy; yi++)
      for(xi=0; xi<svol->dimx; xi++)
        svol->v[zi][yi][xi]=(short int)temp_roundf(g*img->m[zi][yi][xi][fi]);
  svol->scale_factor=1.0/g;

  return(0);
}
/*****************************************************************************/

/*****************************************************************************/
/*!
 * Copy 3D volume as one time frame (1..dimt) into 4D image.
 *  Img must be allocated.
 *
 * @param vol volume structure
 * @param img image structure
 * @param frame frame number [1..number of frames]
 * @return 0 if ok, 1 invalid image status, 2 invalid input, 
 * 3 image<->volume (x,y) dimensions do not match, 
 * 4 image<->volume planes do not match.
 */
int vol2img(VOL *vol, IMG *img, int frame) {
  unsigned short int zi, yi, xi, fi;

  if(VOL_TEST) printf("vol2img(vol, img, %d)\n", frame);
  /* Check input */
  if(vol==NULL || vol->status!=IMG_STATUS_OCCUPIED) return(1);
  vol->statmsg=_volStatusMessage[1];
  if(img==NULL || img->status!=IMG_STATUS_OCCUPIED) return(1);
  if(frame<1 || img->dimt<frame) return(2);
  if(img->dimx!=vol->dimx || img->dimy!=vol->dimy) return(3);
  if(img->dimz!=vol->dimz) return(4);

  /* Copy data */
  fi=frame-1;
  for(zi=0; zi<vol->dimz; zi++)
    for(yi=0; yi<vol->dimy; yi++)
      for(xi=0; xi<vol->dimx; xi++)
        img->m[zi][yi][xi][fi]=vol->v[zi][yi][xi];

  return(0);
}
/*****************************************************************************/

/*****************************************************************************/
/*!
 * Copy 3D short int volume as one time frame (1..dimt) into 4D image.
 *  Img must be allocated.
 *
 * @param svol short volume structure
 * @param img image structure
 * @param frame frame number [1..img->dimt]
 * @return 0 if ok, 1 invalid image, 2 invalid frame number, 
 * 3 (x,y) dimension inconsistency, 4 plane number inconsistency
 */
int svol2img(SVOL *svol, IMG *img, int frame) {
  unsigned short int zi, yi, xi, fi;

  if(VOL_TEST) printf("svol2img(svol, img, %d)\n", frame);
  /* Check input */
  if(svol==NULL || svol->status!=IMG_STATUS_OCCUPIED) return(1);
  svol->statmsg=_volStatusMessage[1];
  if(img==NULL || img->status!=IMG_STATUS_OCCUPIED) return(1);
  if(frame<1 || img->dimt<frame) return(2);
  if(img->dimx!=svol->dimx || img->dimy!=svol->dimy) return(3);
  if(img->dimz!=svol->dimz) return(4);

  /* Copy data */
  fi=frame-1;
  for(zi=0; zi<svol->dimz; zi++)
    for(yi=0; yi<svol->dimy; yi++)
      for(xi=0; xi<svol->dimx; xi++)
        img->m[zi][yi][xi][fi]=(svol->scale_factor)*(float)svol->v[zi][yi][xi];

  return(0);
}
/*****************************************************************************/

/*****************************************************************************/
/*!
 * Prints volume information to specified file pointer, e.g. stdout
 *
 * @param vol volume structure
 * @param fp target file pointer
 */
void volInfo(VOL *vol, FILE *fp) {
  if(VOL_TEST) printf("volInfo()\n");
  if(vol==NULL || vol->status<=IMG_STATUS_UNINITIALIZED) {
    fprintf(fp, "Volume data is not initialized.\n"); return;}
  if(vol->status==IMG_STATUS_INITIALIZED) {
    fprintf(fp, "Volume data is initialized but empty.\n"); return;}
  if(vol->status==IMG_STATUS_ERROR) fprintf(stdout, "Volume data has errors.\n");
  fprintf(fp, "Volume status: %s\n", vol->statmsg);
  fprintf(fp, "Patient orientation: %d\n", vol->orientation);
  fprintf(fp, "Voxel sizes (x, y, z): %g %g %g mm\n",
    vol->sizex, vol->sizey, vol->sizez);
  fprintf(fp, "Dimensions (x, y, z): %d %d %d\n",
    vol->dimx, vol->dimy, vol->dimz);
  return;
}
/*****************************************************************************/

/*****************************************************************************/
/*!
 * Prints short int volume information to specified file pointer, e.g. stdout
 *
 * @param svol short volume structure
 * @param fp target file pointer
 */
void svolInfo(SVOL *svol, FILE *fp) {
  if(VOL_TEST) printf("svolInfo()\n");
  if(svol==NULL || svol->status<=IMG_STATUS_UNINITIALIZED) {
    fprintf(fp, "Volume data is not initialized.\n"); return;}
  if(svol->status==IMG_STATUS_INITIALIZED) {
    fprintf(fp, "Volume data is initialized but empty.\n"); return;}
  if(svol->status==IMG_STATUS_ERROR) fprintf(stdout, "Volume data has errors.\n");
  fprintf(fp, "Volume status: %s\n", svol->statmsg);
  fprintf(fp, "Patient orientation: %d\n", svol->orientation);
  fprintf(fp, "Voxel sizes (x, y, z): %g %g %g mm\n",
    svol->sizex, svol->sizey, svol->sizez);
  fprintf(fp, "Dimensions (x, y, z): %d %d %d\n",
    svol->dimx, svol->dimy, svol->dimz);
  fprintf(fp, "Scale factor: %g\n", svol->scale_factor);
  return;
}
/*****************************************************************************/

/*****************************************************************************/
/** Prints matrix values inside specified range to file pointer.
 *
 * @param vol volume structure
 * @param r volume range structure
 * @param fp target file pointer
 */
void volContents(VOL *vol, VOL_RANGE r, FILE *fp) {
  int zi, yi, xi;

  if(vol==NULL || vol->status!=IMG_STATUS_OCCUPIED) return;
  if(r.z1<1 || r.y1<1 || r.x1<1) return;
  if(r.z2<r.z1 || r.y2<r.y1 || r.x2<r.x1) return;
  if(r.z2>vol->dimz || r.y2>vol->dimy || r.x2>vol->dimx) return;

  for(zi=r.z1-1; zi<r.z2; zi++) {
    fprintf(fp, "pl=%03d ", zi+1);
    for(xi=r.x1-1; xi<r.x2; xi++) fprintf(fp, " x=%05d", xi+1);
    fprintf(fp, "\n");
    for(yi=r.y1-1; yi<r.y2; yi++) {
      fprintf(fp, "y=%05d", yi+1);
      for(xi=r.x1-1; xi<r.x2; xi++)
        fprintf(fp, " %7.3f", vol->v[zi][yi][xi]);
      fprintf(fp, "\n");
    }
  }
}
/*****************************************************************************/

/*****************************************************************************/
/** Finds max and/or min voxel inside specified volume range.
 *
 * @returns 0 if ok, 1 invalid volume status, 2 invalid range endings,
 * 3 inconsistent range dimensions, 4 inconsistent dimensions
 */
int volMax(
  /** Pointer to VOL image structure */
  VOL *vol,
  /** Pointer to volume range inside VOL; enter NULL if whole VOL is used */
  VOL_RANGE *r,
  /** Pixel where max pixel position is written; NULL if not needed */
  VOL_PIXEL *maxp,
  /** Target for max value; NULL if not needed */
  float *maxv,
  /** Pixel where min pixel position is written; NULL if not needed */
  VOL_PIXEL *minp,
  /** Target for min value; NULL if not needed */
  float *minv
) {
  int zi, yi, xi;
  float lmax, lmin;

  if(vol==NULL || vol->status!=IMG_STATUS_OCCUPIED) return(1);
 
  if(r!=NULL) {
    if(r->z1<1 || r->y1<1 || r->x1<1) return(2);
    if(r->z2<r->z1 || r->y2<r->y1 || r->x2<r->x1) return(3);
    if(r->z2>vol->dimz || r->y2>vol->dimy || r->x2>vol->dimx) return(4);

    zi=r->z1-1; yi=r->y1-1; xi=r->x1-1; lmax=lmin=vol->v[zi][yi][xi];
    if(maxp!=NULL) {maxp->z=zi+1; maxp->y=yi+1; maxp->x=xi+1;}
    if(minp!=NULL) {minp->z=zi+1; minp->y=yi+1; minp->x=xi+1;}
    for(zi=r->z1-1; zi<r->z2; zi++) {
      for(yi=r->y1-1; yi<r->y2; yi++) {
        for(xi=r->x1-1; xi<r->x2; xi++) {
          if(lmax<vol->v[zi][yi][xi]) {
            lmax=vol->v[zi][yi][xi];
            if(maxp!=NULL) {maxp->z=zi+1; maxp->y=yi+1; maxp->x=xi+1;}
          } else if(lmin>vol->v[zi][yi][xi]) {
            lmin=vol->v[zi][yi][xi];
            if(minp!=NULL) {minp->z=zi+1; minp->y=yi+1; minp->x=xi+1;}
          }
        }
      }
    }
  } else {
    zi=yi=xi=0; lmax=lmin=vol->v[zi][yi][xi];
    if(maxp!=NULL) {maxp->z=zi+1; maxp->y=yi+1; maxp->x=xi+1;}
    if(minp!=NULL) {minp->z=zi+1; minp->y=yi+1; minp->x=xi+1;}
    for(zi=0; zi<vol->dimz; zi++) {
      for(yi=0; yi<vol->dimy; yi++) {
        for(xi=0; xi<vol->dimx; xi++) {
          if(lmax<vol->v[zi][yi][xi]) {
            lmax=vol->v[zi][yi][xi];
            if(maxp!=NULL) {maxp->z=zi+1; maxp->y=yi+1; maxp->x=xi+1;}
          } else if(lmin>vol->v[zi][yi][xi]) {
            lmin=vol->v[zi][yi][xi];
            if(minp!=NULL) {minp->z=zi+1; minp->y=yi+1; minp->x=xi+1;}
          }
        }
      }
    }
  }
  if(maxv!=NULL) *maxv=lmax;
  if(minv!=NULL) *minv=lmin;
  return(0);
}
/*****************************************************************************/

/*****************************************************************************/
/** Calculates average voxel value inside specified volume range.
 *
 * @returns 0 if ok, 1 invalid image status, 2 inconsistent 
 * range ending, 3 inconsistent range dimensions, 4 inconsistent dimensions
 */
int volAvg(
  /** Pointer to VOL image structure */
  VOL *vol,
  /** Pointer to volume range inside VOL; enter NULL if whole VOL is used */
  VOL_RANGE *r,
  /** Target for mean value */
  float *avg
) {
  int zi, yi, xi, n=0;

  if(vol==NULL || vol->status!=IMG_STATUS_OCCUPIED) return(1);
  if(r!=NULL) {
    if(r->z1<1 || r->y1<1 || r->x1<1) return(2);
    if(r->z2<r->z1 || r->y2<r->y1 || r->x2<r->x1) return(3);
    if(r->z2>vol->dimz || r->y2>vol->dimy || r->x2>vol->dimx) return(4);
  }

  *avg=0.0;
  if(r!=NULL) {
    for(zi=r->z1-1; zi<r->z2; zi++) {
      for(yi=r->y1-1; yi<r->y2; yi++) {
        for(xi=r->x1-1; xi<r->x2; xi++) {
          *avg+=vol->v[zi][yi][xi]; n++;
        }
      }
    }
  } else {
    for(zi=0; zi<vol->dimz; zi++) {
      for(yi=0; yi<vol->dimy; yi++) {
        for(xi=0; xi<vol->dimx; xi++) {
          *avg+=vol->v[zi][yi][xi]; n++;
        }
      }
    }
  }
  if(n>0) *avg/=(float)n;
  return(0);
}
/*****************************************************************************/

/*****************************************************************************/
/** Reorder Volume Range Definition.
\return Returns 0 if successful.
\sa vrdRead, vrdVxlNr, irdRead
 */
int vrdReorder(
  /** Image volume range; start and end range are set in correct order */
  VOL_RANGE *vol_range
) {
  int i;

  /* Check that input is ok */
  if(vol_range==NULL) return 1;
  /* Change the order if necessary */
  if(vol_range->x1<0 || vol_range->x2<0) return 2;
  if(vol_range->x2<vol_range->x1) {
    i=vol_range->x1; vol_range->x1=vol_range->x2; vol_range->x2=i;}
  if(vol_range->y1<0 || vol_range->y2<0) return 3;
  if(vol_range->y2<vol_range->y1) {
    i=vol_range->y1; vol_range->y1=vol_range->y2; vol_range->y2=i;}
  if(vol_range->z1<0 || vol_range->z2<0) return 4;
  if(vol_range->z2<vol_range->z1) {
    i=vol_range->z1; vol_range->z1=vol_range->z2; vol_range->z2=i;}
  return 0;
}
/*****************************************************************************/

/*****************************************************************************/
/** Get the number of voxels in Volume Range Definition.
\return Returns the nr of voxels in the volume range.
\sa vrdReorder, vrdRead
 */
int vrdVxlNr(
  /** Image volume range; start and end range must be in correct order */
  VOL_RANGE *vol_range
) {
  int x, y, z;

  if(vol_range==NULL) return(0);
  z=1+vol_range->z2-vol_range->z1;
  x=1+vol_range->x2-vol_range->x1;
  y=1+vol_range->y2-vol_range->y1;
  return(z*x*y);
}
/*****************************************************************************/

/*****************************************************************************/
/** Read pixel location from a string representation of it.
    @sa string_to_xyzf, vrdRead
    @return Returns 0 if successful, >0 if not.
 */
int string_to_xyz(
  /** String in format x,y,z or x y z ; this string is not changed */
  char *str,
  /** Pixel location in x dimension */
  int *x,
  /** Pixel location in y dimension */
  int *y,
  /** Pixel location in z dimension */
  int *z
) {
  char *cptr, tmp[256];

  strncpy(tmp, str, 255); tmp[255]=(char)0;
  cptr=strtok(tmp, " ,;:()|-"); if(cptr==NULL) return 1;
  *x=atoi(cptr); if(*x<1) return 1;  
  cptr=strtok(NULL, " ,;:()|-"); if(cptr==NULL) return 2;
  *y=atoi(cptr); if(*y<1) return 1;
  cptr=strtok(NULL, " ,;:()|-"); if(cptr==NULL) return 3;
  *z=atoi(cptr); if(*z<1) return 1;
  return 0;
}
/*****************************************************************************/

/*****************************************************************************/
/** Set volume voxel values based on volume range definition.
\return Returns 0 if successful.
\sa vrdRead, vrdVxlNr, volAllocate, vol2img
 */
int vrd2vol(
  /** Image volume range. */
  VOL_RANGE *r,
  /** Pre-allocated image volume data structure. */
  VOL *vol,
  /** Value to write into voxels inside the given volume range; enter nanf("")
      to not change the values. */
  float in,
  /** Value to write into voxels outside the given volume range; enter nanf("")
      to not change the values. */
  float out,
  /** Pointer to a string (allocated for at least 64 chars) where error message
      or other execution status will be written; enter NULL, if not needed */     
  char *status
) {
  if(vol==NULL || vol->status!=IMG_STATUS_OCCUPIED) {
    if(status!=NULL) strcpy(status, "invalid VOL struct");
    return(1);
  }
  if(r->z2<r->z1 || r->y2<r->y1 || r->x2<r->x1) vrdReorder(r);
  if(r->z1<1 || r->y1<1 || r->x1<1 ||
     r->z2>vol->dimz || r->y2>vol->dimy || r->x2>vol->dimx) 
  {
    if(status!=NULL) strcpy(status, "invalid volume range");
    return(2);
  }
  if(isnan(in) && isnan(out)) {
    if(status!=NULL) strcpy(status, "new values not given");
    return(0);
  }
  int zi, yi, xi;

  if(!isnan(in)) { 
    for(zi=r->z1-1; zi<r->z2; zi++)
      for(yi=r->y1-1; yi<r->y2; yi++)
        for(xi=r->x1-1; xi<r->x2; xi++)
          vol->v[zi][yi][xi]=in;
  }

  if(!isnan(out)) { 
    for(zi=1; zi<=vol->dimz; zi++)
      for(yi=1; yi<=vol->dimy; yi++)
        for(xi=1; xi<=vol->dimx; xi++) {
          if(zi>=r->z1 && zi<=r->z2 && 
             yi>=r->y1 && yi<=r->y2 && 
             xi>=r->x1 && xi<=r->x2)
            continue; 
          vol->v[zi-1][yi-1][xi-1]=out;
        }
  }

  if(status!=NULL) strcpy(status, "ok");
  return(0);
}
/*****************************************************************************/

/*****************************************************************************/
/** Read Volume Range Definition File.
\return Returns 0 if successful.
\sa vrdReorder, vrdVxlNr, irdRead
 */
int vrdRead(
  /** Volume Range Definition File filename, which contains image volume corners
   *  (x y z) in IFT format: corner1 = x y z corner2 = x y z */
  char *vrdfile,
  /** Image volume range */
  VOL_RANGE *vol_range,
  /** Pointer to a string (allocated for at least 64 chars) where error message
      or other execution status will be written; enter NULL, if not needed */     
  char *status
) {
  int ret, ii, x, y, z;
  IFT ift;
  char key[256];

  /* Check that input is ok */
  if(vrdfile==NULL || strlen(vrdfile)<1 || vol_range==NULL) {
    if(status!=NULL) strcpy(status, "program error");
    return 1;
  }
  /* Read VDF as IFT file */
  iftInit(&ift); ret=iftRead(&ift, vrdfile, 1); if(ret) {
    if(status!=NULL) strcpy(status, ift.status);
    iftEmpty(&ift); return 2;
  }
  /* Try to find keys 'corner1' and 'corner2' */
  strcpy(key, "corner1"); ii=iftGet(&ift, key);
  if(ii>=0) {
    ret=string_to_xyz(ift.item[ii].value, &x, &y, &z);
    if(ret==0) {
      vol_range->x1=x; vol_range->y1=y; vol_range->z1=z;
      strcpy(key, "corner2"); ii=iftGet(&ift, key);
      if(ii>=0) {
        ret=string_to_xyz(ift.item[ii].value, &x, &y, &z);
        vol_range->x2=x; vol_range->y2=y; vol_range->z2=z;
        if(ret==0) {
          vrdReorder(vol_range);
          if(status!=NULL) strcpy(status, "ok");
          iftEmpty(&ift); return 0;
        }
      }
    }
  }
  /* We are here only if keys were not found */
  /* Lets not care about keys at all */
  for(ii=0, ret=0; ii<ift.keyNr; ii++) {
    if(ret==0 && string_to_xyz(ift.item[ii].value, &x, &y, &z)==0)
    {
      vol_range->x1=x; vol_range->y1=y; vol_range->z1=z;
      ret++; continue;
    }
    if(ret==1 && string_to_xyz(ift.item[ii].value, &x, &y, &z)==0)
    {
      vol_range->x2=x; vol_range->y2=y; vol_range->z2=z;
      ret++; break;
    }
  }
  if(ret<2) {
    if(status!=NULL) strcpy(status, "volume definitions not found");
    iftEmpty(&ift); return 2;
  }

  vrdReorder(vol_range);
  if(status!=NULL) strcpy(status, "ok");
  iftEmpty(&ift);
  return 0;
}
/*****************************************************************************/

/*****************************************************************************/
