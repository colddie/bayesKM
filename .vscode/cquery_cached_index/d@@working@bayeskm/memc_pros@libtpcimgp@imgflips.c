/// @file imgflips.c
/// @author Vesa Oikonen
/// @brief Functions for turning IMG image volume data.
///
/*****************************************************************************/
#include "libtpcimgp.h"
/*****************************************************************************/

/*****************************************************************************/
/** Flip IMG data horizontally (left-right).
    @sa imgFlipRight, imgFlipVertical, imgFlipPlanes, imgFlipAbove, img2cube
 */
void imgFlipHorizontal(IMG *img)
{
  int zi, yi, from, to;
  float *col_ptr;

  for(zi=0; zi<img->dimz; zi++) {
    for(yi=0; yi<img->dimy; yi++) {
      for(from=0, to=img->dimx-1; from<to; from++, to--) {
        col_ptr=img->m[zi][yi][from];
        img->m[zi][yi][from]=img->m[zi][yi][to];
        img->m[zi][yi][to]=col_ptr;
      }
    }
  }
}
/*****************************************************************************/

/*****************************************************************************/
/** Flip IMG data vertically (up-down).
    @sa imgFlipHorizontal, imgFlipRight, imgFlipPlanes, imgFlipAbove
 */
void imgFlipVertical(IMG *img)
{
  int zi, from, to;
  float **row_ptr;

  for(zi=0; zi<img->dimz; zi++) {
    for(from=0, to=img->dimy-1; from<to; from++, to--) {
      row_ptr=img->m[zi][from];
      img->m[zi][from]=img->m[zi][to];
      img->m[zi][to]=row_ptr;
    }
  }
}
/*****************************************************************************/

/*****************************************************************************/
/** Flip IMG data planes (head-toes). 

    To work properly, the plane numbers must be contiguous.
    @sa imgFlipHorizontal, imgFlipVertical, imgFlipRight, imgFlipAbove
 */
void imgFlipPlanes(IMG *img)
{
  int from, to;
  float ***plane_ptr;

  for(from=0, to=img->dimz-1; from<to; from++, to--) {
    plane_ptr=img->m[from];
    img->m[from]=img->m[to];
    img->m[to]=plane_ptr;
  }
}
/*****************************************************************************/

/*****************************************************************************/
/** Flip IMG data like viewed from right side.
    @sa imgFlipHorizontal, imgFlipVertical, imgFlipPlanes, imgFlipAbove
    @return Returns 0 if successful.
 */
int imgFlipRight(
  /** Pointer to IMG which will be flipped */
  IMG *img
) {
  int xi, yi, zi, fi, ret;
  IMG omg;

  /* Make copy of the original image */
  imgInit(&omg); ret=imgDup(img, &omg); if(ret!=0) return(100+ret);

  /* Empty the user-specified image struct */
  imgEmpty(img);
  /* Allocate it again with new dimensions */
  ret=imgAllocateWithHeader(img, omg.dimx, omg.dimy, omg.dimz,
                            omg.dimt, &omg);
  if(ret!=0) {imgEmpty(&omg); return(200+ret);}

  /* Copy the voxel values */
  for(zi=0; zi<omg.dimz; zi++)
    for(yi=0; yi<omg.dimy; yi++)
      for(xi=0; xi<omg.dimx; xi++)
        for(fi=0; fi<omg.dimt; fi++) {
          img->m[xi][yi][zi][fi]=omg.m[zi][yi][xi][fi];
        }

  /* Switch pixel sizes */
  img->sizex=omg.sizez;
  img->sizez=omg.sizex;
  img->resolutionx=omg.resolutionz;
  img->resolutionz=omg.resolutionx;

  /* Free the memory of the copy of original image */
  imgEmpty(&omg);

  return(0);
}
/*****************************************************************************/

/*****************************************************************************/
/** Flip IMG data like viewed from above.
    @sa imgFlipHorizontal, imgFlipVertical, imgFlipPlanes, imgFlipRight, img2cube
    @return Returns 0 if successful.
 */
int imgFlipAbove(
  /** Pointer to IMG which will be flipped */
  IMG *img
) {
  int xi, yi, zi, fi, ret;
  IMG omg;

  /* Make copy of the original image */
  imgInit(&omg); ret=imgDup(img, &omg); if(ret!=0) return(100+ret);

  /* Empty the user-specified image struct */
  imgEmpty(img);
  /* Allocate it again with new dimensions */
  ret=imgAllocateWithHeader(img, omg.dimy, omg.dimz, omg.dimx,
                            omg.dimt, &omg);
  if(ret!=0) {imgEmpty(&omg); return(200+ret);}

  /* Copy the voxel values */
  for(zi=0; zi<omg.dimz; zi++)
    for(yi=0; yi<omg.dimy; yi++)
      for(xi=0; xi<omg.dimx; xi++)
        for(fi=0; fi<omg.dimt; fi++) {
          img->m[yi][img->dimy-1-zi][xi][fi]=omg.m[zi][yi][xi][fi];
        }

  /* Switch pixel sizes */
  img->sizey=omg.sizez;
  img->sizez=omg.sizey;
  img->resolutiony=omg.resolutionz;
  img->resolutionz=omg.resolutiony;

  /* Free the memory of the copy of original image */
  imgEmpty(&omg);

  return(0);
}
/*****************************************************************************/

/*****************************************************************************/
