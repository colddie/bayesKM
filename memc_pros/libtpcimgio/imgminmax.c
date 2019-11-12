/// @file imgminmax.c
/// @author Vesa Oikonen
/// @brief Searching min and max in IMG data.
///
/******************************************************************************/
#include "libtpcimgio.h"
/******************************************************************************/

/******************************************************************************/
/** Search the max pixel value in the IMG data.

   @sa imgFrameMinMax, imgMinMax, imgSmoothMax, imgAvg
   @return 0 if ok, 1 invalid image status, 2 invalid output pointer, 3 invalid image dimensions.
 */
int imgMax(
  /** Pointer to IMG struct. */
  IMG *img,
  /** Pointer to output. */
  float *maxvalue
) {
  int pi, yi, xi, fi;
  float f;

  if(img->status<IMG_STATUS_OCCUPIED) return(1);
  if(maxvalue==NULL) return(2); else *maxvalue=0.0;
  if(img->dimt<1 || img->dimz<1 || img->dimy<1 || img->dimx<1) return(3);
  f=img->m[0][0][0][0];
  for(pi=0; pi<img->dimz; pi++)
    for(yi=0; yi<img->dimy; yi++)
      for(xi=0; xi<img->dimx; xi++)
        for(fi=0; fi<img->dimt; fi++) {
          if(!(img->m[pi][yi][xi][fi]<=f)) f=img->m[pi][yi][xi][fi];
        }
  *maxvalue=f;
  return(0);
}
/******************************************************************************/

/******************************************************************************/
/** Searches the max absolute pixel value in the IMG data.

   Sets maxvalue to the absolute max value with sign.

   @sa imgFrameMinMax, imgMinMax, imgMax, imgSmoothMax
   @return 0 if ok, 1 invalid image status, 2 invalid output pointer, 3 invalid image dimensions.
 */
int imgAbsMax(
  /** Pointer to IMG struct. */
  IMG *img,
  /** Pointer to output. */
  float *maxvalue
) {
  int pi, yi, xi, fi;
  float f;

  if(img->status<IMG_STATUS_OCCUPIED) return(1);
  if(maxvalue==NULL) return(2); else *maxvalue=0.0;
  if(img->dimt<1 || img->dimz<1 || img->dimy<1 || img->dimx<1) return(3);
  f=img->m[0][0][0][0];
  for(pi=0; pi<img->dimz; pi++)
    for(yi=0; yi<img->dimy; yi++)
      for(xi=0; xi<img->dimx; xi++)
        for(fi=0; fi<img->dimt; fi++) {
          if(!(fabs(img->m[pi][yi][xi][fi])<=fabs(f))) f=img->m[pi][yi][xi][fi];
        }
  *maxvalue=f;
  return(0);
}
/******************************************************************************/

/******************************************************************************/
/** Finds max and/or min voxel inside specified image range.
   @sa imgFrameMinMax, imgMinMax, imgAvg
   @returns 0 if ok, 1 invalid volume status, 2 invalid range endings,
   3 inconsistent range dimensions, 4 inconsistent dimensions
 */
int imgRangeMinMax(
  /** Pointer to IMG structure */
  IMG *img,
  /** Pointer to image range inside IMG; enter NULL if whole IMG is used */
  IMG_RANGE *r,
  /** Pixel where max pixel position is written; NULL if not needed */
  IMG_PIXEL *maxp,
  /** Target for max value; NULL if not needed */
  float *maxv,
  /** Pixel where min pixel position is written; NULL if not needed */
  IMG_PIXEL *minp,
  /** Target for min value; NULL if not needed */
  float *minv
) {
  int zi, yi, xi, fi;
  float lmax, lmin;

  if(img->status<IMG_STATUS_OCCUPIED) return(1);
  if(img->dimt<1 || img->dimz<1 || img->dimy<1 || img->dimx<1) return(1);

  if(r!=NULL) {
    if(r->z1<1 || r->y1<1 || r->x1<1 || r->f1<1) return(2);
    if(r->z2<r->z1 || r->y2<r->y1 || r->x2<r->x1 || r->f2<r->f1) return(3);
    if(r->z2>img->dimz || r->y2>img->dimy || r->x2>img->dimx || r->f2>img->dimt) return(4);

    zi=r->z1-1; yi=r->y1-1; xi=r->x1-1; fi=r->f1-1;
    lmax=lmin=img->m[zi][yi][xi][fi];
    if(maxp!=NULL) {maxp->z=zi+1; maxp->y=yi+1; maxp->x=xi+1; maxp->f=fi+1;}
    if(minp!=NULL) {minp->z=zi+1; minp->y=yi+1; minp->x=xi+1; minp->f=fi+1;}
    for(zi=r->z1-1; zi<r->z2; zi++) {
      for(yi=r->y1-1; yi<r->y2; yi++) {
        for(xi=r->x1-1; xi<r->x2; xi++) {
          for(fi=r->f1-1; fi<r->f2; fi++) {
            if(!(lmax>=img->m[zi][yi][xi][fi])) {
              lmax=img->m[zi][yi][xi][fi];
              if(maxp!=NULL) {maxp->z=zi+1; maxp->y=yi+1; maxp->x=xi+1; maxp->f=fi+1;}
            } else if(!(lmin<=img->m[zi][yi][xi][fi])) {
              lmin=img->m[zi][yi][xi][fi];
              if(minp!=NULL) {minp->z=zi+1; minp->y=yi+1; minp->x=xi+1; minp->f=fi+1;}
            }
          }
        }
      }
    }
  } else {
    zi=yi=xi=fi=0; lmax=lmin=img->m[zi][yi][xi][fi];
    if(maxp!=NULL) {maxp->z=zi+1; maxp->y=yi+1; maxp->x=xi+1; maxp->f=fi+1;}
    if(minp!=NULL) {minp->z=zi+1; minp->y=yi+1; minp->x=xi+1; minp->f=fi+1;}
    for(zi=0; zi<img->dimz; zi++) {
      for(yi=0; yi<img->dimy; yi++) {
        for(xi=0; xi<img->dimx; xi++) {
          for(fi=0; fi<img->dimt; fi++) {
            if(!(lmax>=img->m[zi][yi][xi][fi])) {
              lmax=img->m[zi][yi][xi][fi];
              if(maxp!=NULL) {maxp->z=zi+1; maxp->y=yi+1; maxp->x=xi+1; maxp->f=fi+1;}
            } else if(!(lmin<=img->m[zi][yi][xi][fi])) {
              lmin=img->m[zi][yi][xi][fi];
              if(minp!=NULL) {minp->z=zi+1; minp->y=yi+1; minp->x=xi+1; minp->f=fi+1;}
            }
          }
        }
      }
    }
  }
  if(maxv!=NULL) *maxv=lmax;
  if(minv!=NULL) *minv=lmin;
  return(0);
}
/******************************************************************************/

/******************************************************************************/
/** Searches the min and max pixel value in the IMG data.
   @sa imgFrameMinMax, imgRangeMinMax, imgReadMinMax
   @return Returns 0 when successful.
 */
int imgMinMax(
  /** Pointer to IMG struct from where min and pixels are searched */
  IMG *img,
  /** Pointer to min pixel value; Enter NULL if not needed */
  float *minvalue,
  /** Pointer to max pixel value; Enter NULL if not needed */
  float *maxvalue
) {
  return (imgRangeMinMax(img, NULL, NULL, maxvalue, NULL, minvalue));
}
/******************************************************************************/

/******************************************************************************/
/** Searches the min and max pixel value in one frame (1..dimt) of the IMG data.
   @return 0 if ok, 1 invalid image status, 2 invalid output pointer, 3 invalid image dimensions.
   @sa imgMinMax, imgRangeMinMax, imgReadMinMax, imgGetMaxFrame
 */
int imgFrameMinMax(
  /** Pointer to IMG data. */
  IMG *img, 
  /** Frame number [1..number of frames]. */
  int frame, 
  /** Pointer to float value where minimum is written. */
  float *minvalue, 
  /** Pointer to float value where maximum is written. */
  float *maxvalue
) {
  int pi, yi, xi, fi;
  float mi, ma;

  if(img->status<IMG_STATUS_OCCUPIED) return(1);
  if(minvalue==NULL || maxvalue==NULL) return(2);
  *minvalue=*maxvalue=0.0; fi=frame-1;
  if(img->dimt<frame || img->dimz<1 || img->dimy<1 || img->dimx<1) return(3);
  if(frame<1) return(4);
  mi=ma=img->m[0][0][0][fi];
  for(pi=0; pi<img->dimz; pi++)
    for(yi=0; yi<img->dimy; yi++)
      for(xi=0; xi<img->dimx; xi++) {
        if(!(img->m[pi][yi][xi][fi]<=ma)) ma=img->m[pi][yi][xi][fi];
        else if(!(img->m[pi][yi][xi][fi]>=mi)) mi=img->m[pi][yi][xi][fi];
      }
  *minvalue=mi; *maxvalue=ma;
  return(0);
}
/******************************************************************************/

/******************************************************************************/
/** Read the calibrated maximum and minimum pixel values in the specified file
    in ECAT 7, ECAT 6.3, or Analyze 7.5 format.

    File is read frame-by-frame with normal IMG functions.
  
   @sa imgFrameMinMax, imgRangeMinMax, imgMinMax
   @return errstatus, which is STATUS_OK (0) when call was successful, and >0 in case of an error.
 */
int imgReadMinMax(
  /** ECAT 7 or ECAT 6.3 filename, or Analyze 7.5 database. */
  const char *fname, 
  /** Pointer to minimum pixel value that will be set by this function. */
  float *fmin, 
  /** Pointer to maximum pixel value that will be set by this function. */
  float *fmax
) {
  int fi=0, ret;
  IMG img;
  float frmin, frmax;

  if(IMG_TEST) printf("imgReadMinMax(%s, *fmin, *fmax)\n", fname);
  imgInit(&img);
  while((ret=imgReadFrame(fname, fi+1, &img, 0)) == 0) {
    if(imgMinMax(&img, &frmin, &frmax)!=0) {imgEmpty(&img); return STATUS_FAULT;}
    if(fi==0) {
      if(fmin!=NULL) *fmin=frmin;
      if(fmin!=NULL) *fmax=frmax;
    } else {
      if(fmin!=NULL && !(*fmin<=frmin)) *fmin=frmin;
      if(fmax!=NULL && !(*fmax>=frmax)) *fmax=frmax;
    }
    fi++;
  } /* next frame */
  imgEmpty(&img);
  if(ret==STATUS_NOMATRIX && fi>0) return STATUS_OK;
  else return ret;
}
/******************************************************************************/

/******************************************************************************/
/** Searches the spatially (3x3) smoothed max pixel value in the IMG data.
 
   @sa imgFrameMinMax, imgMinMax, imgMax, imgAbsMax, imgGetPeak
   @return 0 if ok, 1 invalid image status, 2 invalid output pointer, 3 invalid image dimensions.
 */
int imgSmoothMax(
  /** Pointer to IMG struct */
  IMG *img,
  /** Pointer to float in which max pixel value will be written;
   *  enter NULL if not needed */
  float *maxvalue,
  /** Pointer to struct in which the position of max pixel will be written
   *  (1-based positions); enter NULL if not needed */
  IMG_PIXEL *p 
) {
  int pi, yi, xi, fi;
  float f, v;

  if(img->status<IMG_STATUS_OCCUPIED) return(1);
  if(maxvalue==NULL && p==NULL) return(2);
  if(img->dimt<1 || img->dimz<1 || img->dimy<3 || img->dimx<3) return(3);
  if(maxvalue!=NULL) *maxvalue=0.0;
  if(p!=NULL) p->x=p->y=p->z=p->f=1;
  f=-1.0E20;
  for(pi=0; pi<img->dimz; pi++)
    for(yi=1; yi<img->dimy-1; yi++)
      for(xi=1; xi<img->dimx-1; xi++)
        for(fi=0; fi<img->dimt; fi++) {
          v=img->m[pi][yi-1][xi-1][fi]+
            img->m[pi][yi-1][xi  ][fi]+
            img->m[pi][yi-1][xi+1][fi]+
            img->m[pi][yi  ][xi-1][fi]+
            img->m[pi][yi  ][xi  ][fi]*2.0+
            img->m[pi][yi  ][xi+1][fi]+
            img->m[pi][yi+1][xi-1][fi]+
            img->m[pi][yi+1][xi  ][fi]+
            img->m[pi][yi+1][xi+1][fi];
          v*=0.1;
          if(v>f) {
            f=v; if(p!=NULL) {p->x=xi+1; p->y=yi+1; p->z=pi+1; p->f=fi+1;}}
        }
  if(maxvalue!=NULL) *maxvalue=f;
  return(0);
}
/******************************************************************************/

/******************************************************************************/
/** Searches the max pixel value in the IMG data, which occurs before specified time.
   @sa imgGetMaxTime, imgGetMaxFrame
   @return Returns 0 if successful.
 */
int imgGetPeak(
  /** Pointer to IMG struct */
  IMG *img,
  /** Time (sec) after which max value is not searched */
  float beforeTime,
  /** Pointer to struct where max pixel position is written */
  IMG_PIXEL *p,
  /** Verbose level; 0 if nothing is to be printed in stdout */
  int verbose
) {
  int zi, yi, xi, fi, mf;
  float f;

  if(verbose>0) printf("imgGetPeak(img, %g, p, %d)\n", beforeTime, verbose);
  if(img->status<IMG_STATUS_OCCUPIED) return(1);
  if(p==NULL) return(2);
  if(img->dimt<1 || img->dimz<1 || img->dimy<1 || img->dimx<1) return(3);
  if(beforeTime<img->mid[0]) {
    if(verbose>0) fprintf(stderr, "Error: invalid max search time setting.\n");
    return(4);
  }
  f=img->m[0][0][0][0]-1.0; mf=img->dimt; p->x=p->y=p->z=p->f=1;
  for(zi=0; zi<img->dimz; zi++) {
    for(yi=0; yi<img->dimy; yi++) {
      for(xi=0; xi<img->dimx; xi++) {
        for(fi=0; fi<img->dimt; fi++) if(img->mid[fi]<=beforeTime) {
          if(img->m[zi][yi][xi][fi]<f) // lower
            continue;
          if(img->m[zi][yi][xi][fi]==f) { // equal
            // only use this if earlier than in prev max
            if(fi>=mf) continue;
          } 
          f=img->m[zi][yi][xi][fi];
          p->x=xi+1; p->y=yi+1; p->z=zi+1; p->f=fi+1;
          mf=fi;
        }
      }
    }
  }
  if(verbose>2) printf("maxval := %g\n", f);
  return(0);
}
/******************************************************************************/

/******************************************************************************/
/** Search the time of maximum value for each image pixel separately.
    @sa imgGetPeak, imgGetMaxFrame
    @return Returns 0, if ok.
 */
int imgGetMaxTime(
  /** Image to search for max frames; not modified. */
  IMG *img,
  /** Pointer to empty IMG struct in which the max time (sec) will be written;
      any old contents are deleted. */
  IMG *mimg,
  /** Just save the frame middle time (0), or compute value weighted average time of
      all frames (1), or, value weighted average time of 3 or 5 subsequent frames (2). 
      Option (1) uses the equation for mean residence time (MRT) for PTACs in pharmacokinetics.
   */
  const int w,
  /** Verbose level; if zero, then nothing is printed to stderr or stdout. */
  int verbose
) {
  if(verbose>0) printf("imgGetMaxTime(*img, *mimg, %d)\n", w);
  
  if(img==NULL || img->status!=IMG_STATUS_OCCUPIED) return(1);
  if(mimg==NULL) return(2);
  if(mimg->status==IMG_STATUS_OCCUPIED) imgEmpty(mimg);

  /* Allocate memory for one frame */
  int ret;
  if(verbose>1) printf("allocating memory for %dx%dx%d pixels\n", img->dimz, img->dimy, img->dimx);
  ret=imgAllocate(mimg, img->dimz, img->dimy, img->dimx, 1);
  if(ret) return(ret);
  /* set image header information */
  imgCopyhdr(img, mimg);
  mimg->start[0]=img->start[0]; mimg->end[0]=img->end[img->dimt-1];
  mimg->mid[0]=(mimg->start[0]+mimg->end[0])/2.0;
  mimg->unit=CUNIT_UNKNOWN;

  if(w==0) {
    for(int zi=0; zi<img->dimz; zi++) {
      for(int yi=0; yi<img->dimy; yi++) {
        for(int xi=0; xi<img->dimx; xi++) {
          /* Find the frame with max value */
          int ti=0, mi=0; 
          double mv=img->m[zi][yi][xi][ti];
          for(ti=1; ti<img->dimt; ti++) {
            if(img->m[zi][yi][xi][ti]<mv) continue;
            mi=ti; mv=img->m[zi][yi][xi][ti];
          }
          if(mv>0.0) mimg->m[zi][yi][xi][0]=img->mid[mi];
          else mimg->m[zi][yi][xi][0]=0.0;
        }
      }
    }
    return(0);
  }

  if(w==1) { // Same as the equation for mean residence time (MRT) for PTACs in pharmacokinetics
    for(int zi=0; zi<img->dimz; zi++) {
      for(int yi=0; yi<img->dimy; yi++) {
        for(int xi=0; xi<img->dimx; xi++) {
          /* Compute the value weighted time */
          double sumw=0.0, sumt=0.0;
          for(int ti=0; ti<img->dimt; ti++) {
            if(isnan(img->m[zi][yi][xi][ti])) continue;
            float fdur=img->end[ti]-img->start[ti]; if(fdur<=0.0) fdur=1.0;
            sumt+=img->m[zi][yi][xi][ti]*img->mid[ti]*fdur;
            sumw+=img->m[zi][yi][xi][ti]*fdur;
          }
          sumt/=sumw;
          if(sumt>0.0 && sumw>0.0) mimg->m[zi][yi][xi][0]=sumt;
          else mimg->m[zi][yi][xi][0]=0.0;
        }
      }
    }
    return(0);
  }

  if(w>1) {
    for(int zi=0; zi<img->dimz; zi++) {
      for(int yi=0; yi<img->dimy; yi++) {
        for(int xi=0; xi<img->dimx; xi++) {
          /* Find the frame with max value */
          int ti=0, mi=0; 
          double mv=img->m[zi][yi][xi][ti];
          for(ti=1; ti<img->dimt; ti++) {
            if(img->m[zi][yi][xi][ti]<mv) continue;
            mi=ti; mv=img->m[zi][yi][xi][ti];
          }
          /* calculate weighted mean from subsequent frames, if possible */
          if(mi<1 || mi>img->dimt-2) continue;
          int i1, i2; i1=mi-1; i2=mi+1; if(i1>0 && i2<img->dimt-1) {i1--; i2++;}
          double sumw=0.0, sumt=0.0;
          for(int i=i1; i<=i2; i++) {
            if(!(img->m[zi][yi][xi][i]>0.0)) continue;
            sumt+=img->m[zi][yi][xi][i]*img->mid[i];
            sumw+=img->m[zi][yi][xi][i];
          }
          sumt/=sumw;
          if(sumt>0.0) mimg->m[zi][yi][xi][0]=sumt;
          else mimg->m[zi][yi][xi][0]=0.0;
        }
      }
    }
    return(0);
  }

  return(0);
}
/******************************************************************************/

/******************************************************************************/
/** Search the frame with maximum pixel value for each image pixel separately.
    @sa imgGetPeak, imgGetMaxTime
    @return Returns 0, if ok.
 */
int imgGetMaxFrame(
  /** Image to search for max frames; not modified */
  IMG *img,
  /** Pointer to empty IMG struct in which the nr of frame with max pixel
      value will be written; 0 is written if max value is <= 0;
      any old contents are deleted. */
  IMG *mimg,
  /** Verbose level; if zero, then nothing is printed to stderr or stdout */
  int verbose
) {
  if(verbose>0) printf("imgGetMaxFrame()\n");
  
  if(img==NULL || img->status!=IMG_STATUS_OCCUPIED) return(1);
  if(mimg==NULL) return(2);
  if(mimg->status==IMG_STATUS_OCCUPIED) imgEmpty(mimg);

  /* Allocate memory for one frame */
  int ret;
  if(verbose>1) printf("allocating memory for %dx%dx%d pixels\n", img->dimz, img->dimy, img->dimx);
  ret=imgAllocate(mimg, img->dimz, img->dimy, img->dimx, 1);
  if(ret) return(ret);
  /* set image header information */
  imgCopyhdr(img, mimg);
  mimg->start[0]=img->start[0]; mimg->end[0]=img->end[img->dimt-1];
  mimg->mid[0]=(mimg->start[0]+mimg->end[0])/2.0;

  /* Go through every pixel */
  int ti, zi, xi, yi;
  double mv; int mi;
  for(zi=0; zi<img->dimz; zi++) {
    for(yi=0; yi<img->dimy; yi++) for(xi=0; xi<img->dimx; xi++) {
      ti=mi=0; mv=img->m[zi][yi][xi][ti];
      for(ti=1; ti<img->dimt; ti++) {
        if(img->m[zi][yi][xi][ti]<mv) continue;
        mi=ti; mv=img->m[zi][yi][xi][ti];
      }
      if(mv>1.0E-008) mimg->m[zi][yi][xi][0]=1.0+mi;
      else mimg->m[zi][yi][xi][0]=0.0;
    }
  }

  return(0);
}
/******************************************************************************/

/******************************************************************************/
/** Calculates average voxel value inside specified image range.
   @sa imgMinMax, imgMax, fmedian, imgTimeIntegral
   @return Returns 0 when successful.
 */
int imgAvg(
  /** Pointer to IMG image structure */
  IMG *img,
  /** Pointer to range inside IMG; enter NULL if whole IMG is used */
  IMG_RANGE *r,
  /** Target for mean value */
  float *avg
) {
  int zi, yi, xi, fi, n=0;

  if(img->status<IMG_STATUS_OCCUPIED) return(1);
  if(r!=NULL) {
    if(r->z1<1 || r->y1<1 || r->x1<1 || r->f1<1) return(2);
    if(r->z2<r->z1 || r->y2<r->y1 || r->x2<r->x1 || r->f2<r->f1) return(3);
    if(r->z2>img->dimz || r->y2>img->dimy || r->x2>img->dimx || r->f2>img->dimt) return(4);
  }
  if(avg==NULL) return(5);

  *avg=0.0;
  if(r!=NULL) {
    for(zi=r->z1-1; zi<r->z2; zi++) {
      for(yi=r->y1-1; yi<r->y2; yi++) {
        for(xi=r->x1-1; xi<r->x2; xi++) {
          for(fi=r->f1-1; fi<r->f2; fi++) if(!isnan(img->m[zi][yi][xi][fi])) {
            *avg+=img->m[zi][yi][xi][fi]; n++;
          }
        }
      }
    }
  } else {
    for(zi=0; zi<img->dimz; zi++) {
      for(yi=0; yi<img->dimy; yi++) {
        for(xi=0; xi<img->dimx; xi++) {
          for(fi=0; fi<img->dimt; fi++) if(!isnan(img->m[zi][yi][xi][fi])) {
            *avg+=img->m[zi][yi][xi][fi]; n++;
          }
        }
      }
    }
  }
  if(n>0) *avg/=(float)n;
  return(0);
}
/******************************************************************************/

/******************************************************************************/
/** Returns the kth smallest value in data[0..n-1]. 

    Array is partially sorted.
    Algorithm is based on the book Wirth N. Algorithms + data structures =
    programs. Englewood Cliffs, Prentice-Hall, 1976.
   @sa fmedian
   @return Returns the kth smallest value in data[0..n-1].
*/
float f_kth_smallest(
  /** Pointer to data; array is partially sorted */
  float *data,
  /** Length of data array */
  int n,
  /** kth smallest value will be returned */
  int k
) {
  int i, j, l, m;
  float x, s;

  l=0; m=n-1;
  while(l<m) {
    x=data[k]; i=l; j=m;
    do {
      while(data[i]<x) i++;
      while(x<data[j]) j--;
      if(i<=j) {s=data[i]; data[i]=data[j]; data[j]=s; i++; j--;}
    } while(i<=j);
    if(j<k) l=i;
    if(k<i) m=j;
  }
  return(data[k]);
}
/******************************************************************************/

/******************************************************************************/
/** Returns the median in array data[0..n-1]. 

    Array is partially sorted.
    Algorithm is based on the book Wirth N. Algorithms + data structures =
    programs. Englewood Cliffs, Prentice-Hall, 1976.
    @sa fMinMax, imgAvg, dmedian
    @return Returns the median in array data[0..n-1].
*/
float fmedian(
  /** Pointer to data; array is partially sorted */
  float *data,
  /** Length of data array */
  int n
) {
  int k;
  float d1, d2;

  if(n<1) return(0.0);
  if(n%2) {
    k=(n-1)/2; return(f_kth_smallest(data, n, k));
  } else {
    k=n/2; d1=f_kth_smallest(data, n, k-1); d2=f_kth_smallest(data, n, k);
    return(0.5*(d1+d2));
  }
}
/******************************************************************************/

/******************************************************************************/
/** Returns the mean in array data[0..n-1], and optionally calculates also
    the (sample) standard deviation of the mean.
    @sa dmean, fmedian, mean, dmean_nan
    @return Returns the mean in array data[0..n-1].
*/
float fmean(
  /** Pointer to data; data is not changed in any way. */
  float *data,
  /** Length of data array. */
  int n,
  /** Pointer to variable where SD will be written; enter NULL if not needed. */
  float *sd
) {
  int i;
  float sumsqr=0.0, sqrsum=0.0, avg;

  if(n<1 || data==NULL) {if(sd!=NULL) *sd=0.0; return(0.0);}

  for(i=0; i<n; i++) {sumsqr+=data[i]*data[i]; sqrsum+=data[i];}
  avg=sqrsum/(float)n; if(sd==NULL) return(avg);
  if(n==1) {
    *sd=0.0;
  } else {
    sqrsum*=sqrsum;
    *sd=sqrt( (sumsqr - sqrsum/(float)n) / (float)(n-1) );
  }
  return(avg);
}
/******************************************************************************/

/******************************************************************************/
/** Finds the minimum and maximum value in a float array.

    Only finite values are considered.
   @sa imgMinMax, imgAbsMax, imgAvg
*/
void fMinMaxFin(
  /** Pointer to float array of size n. */
  float *data,
  /** Array length. */
  int n,
  /** Pointer to float value for minimum; enter NULL if not needed. */
  float *fmin,
  /** Pointer to float value for maximum; enter NULL if not needed. */
  float *fmax
) {
  if(fmin!=NULL) *fmin=nanf("");
  if(fmax!=NULL) *fmax=nanf("");

  int i;
  for(i=0; i<n; i++) if(isfinite(data[i])) break;
  if(i==n) return; // no finite values found
  float mi, ma;
  mi=ma=data[i++];
  for(; i<n; i++) if(isfinite(data[i])) {
    if(data[i]>ma) ma=data[i];
    else if(data[i]<mi) mi=data[i];
  }
  if(fmin!=NULL) *fmin=mi;
  if(fmax!=NULL) *fmax=ma;
  return;
}
/******************************************************************************/

/******************************************************************************/
