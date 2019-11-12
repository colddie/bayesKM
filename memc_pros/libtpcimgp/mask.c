/// @file mask.c
/// @brief Functions for processing mask images.
/// @author Vesa Oikonen
///
/*****************************************************************************/
#include "libtpcimgp.h"
/*****************************************************************************/

/*****************************************************************************/
/** Count the nr of positive values inside 3D mask image.

    @sa imgMaskRoiNr, imgThresholdMask, imgMaskInvert
    @return Returns the number of positive voxels.
 */
unsigned int imgMaskCount(
  /** Pointer to mask IMG struct. */
  IMG *img
) {
  unsigned int n=0;
  int zi, yi, xi;
  for(zi=0; zi<img->dimz; zi++)
    for(yi=0; yi<img->dimy; yi++)
      for(xi=0; xi<img->dimx; xi++)
        if(img->m[zi][yi][xi][0]>0.0) n++;
   return(n);
}
/*****************************************************************************/

/*****************************************************************************/
/** Erode the 3D mask image.
    @sa imgStructuringElement, imgMaskDilate, imgMaskCount
    @return Returns the number of eroded voxels, or <0 in case of an error.
 */
int imgMaskErode(
  /** Pointer to IMG struct containing the mask image. */
  IMG *img,
  /** Pointer to IMG struct containing the structuring element;
      all dimensions must be odd numbers. */
  IMG *se
) {
  if(img==NULL || se==NULL) return(-1);
  if(img->dimx<1 || img->dimy<1 || img->dimz<1) return(-2);
  if(se->dimx<1 || se->dimy<1 || se->dimz<1) return(-3);
  if(se->dimx%2==0 || se->dimy%2==0 || se->dimz%2==0) return(-4);


  /* make a copy of original data */
  IMG orig; imgInit(&orig);
  if(imgDup(img, &orig)!=0) return(-5);
  int n=0;
  int zi, yi, xi, zj, yj, xj;
  for(zi=0; zi<orig.dimz; zi++)
    for(yi=0; yi<orig.dimy; yi++)
      for(xi=0; xi<orig.dimx; xi++)
        if(orig.m[zi][yi][xi][0]>0.0) { 
          int zs, ys, xs;
          float miv=+1.0E+20;
          for(zs=0; zs<se->dimz; zs++)
            for(ys=0; ys<se->dimy; ys++)
              for(xs=0; xs<se->dimx; xs++)
                if(se->m[zs][ys][xs][0]>0.0) {
                  zj=zi+(zs-se->dimz/2); if(zj<0 || zj>=orig.dimz) continue;
                  yj=yi+(ys-se->dimy/2); if(yj<0 || yj>=orig.dimy) continue;
                  xj=xi+(xs-se->dimx/2); if(xj<0 || xj>=orig.dimx) continue;
                  if(orig.m[zj][yj][xj][0]<miv) miv=orig.m[zj][yj][xj][0];
                }
          if(!(miv>1.0E-20)) {img->m[zi][yi][xi][0]=0.0; n++;}
        }
  imgEmpty(&orig);

  return(n);
}
/*****************************************************************************/

/*****************************************************************************/
/** Dilate the 3D mask image.
    @sa imgStructuringElement, imgMaskErode, imgMaskCount
    @return Returns the number of dilated voxels, or <0 in case of an error.
 */
int imgMaskDilate(
  /** Pointer to IMG struct containing the mask image. */
  IMG *img,
  /** Pointer to IMG struct containing the structuring element;
      all dimensions must be odd numbers. */
  IMG *se
) {
  if(img==NULL || se==NULL) return(-1);
  if(img->dimx<1 || img->dimy<1 || img->dimz<1) return(-2);
  if(se->dimx<1 || se->dimy<1 || se->dimz<1) return(-3);
  if(se->dimx%2==0 || se->dimy%2==0 || se->dimz%2==0) return(-4);


  /* make a copy of original data */
  IMG orig; imgInit(&orig);
  if(imgDup(img, &orig)!=0) return(-5);
  int n=0;
  int zi, yi, xi, zj, yj, xj;
  for(zi=0; zi<orig.dimz; zi++)
    for(yi=0; yi<orig.dimy; yi++)
      for(xi=0; xi<orig.dimx; xi++)
        if(orig.m[zi][yi][xi][0]==0.0) { 
          int zs, ys, xs;
          float mav=-1.0E+20;
          for(zs=0; zs<se->dimz; zs++)
            for(ys=0; ys<se->dimy; ys++)
              for(xs=0; xs<se->dimx; xs++)
                if(se->m[zs][ys][xs][0]>0.0) {
                  zj=zi+(zs-se->dimz/2); if(zj<0 || zj>=orig.dimz) continue;
                  yj=yi+(ys-se->dimy/2); if(yj<0 || yj>=orig.dimy) continue;
                  xj=xi+(xs-se->dimx/2); if(xj<0 || xj>=orig.dimx) continue;
                  if(orig.m[zj][yj][xj][0]>mav) mav=orig.m[zj][yj][xj][0];
                }
          if(mav>0.0) {img->m[zi][yi][xi][0]=mav; n++;}
        }
  imgEmpty(&orig);

  return(n);
}
/*****************************************************************************/

/*****************************************************************************/
/** Make 3D structuring element for eroding and dilation.
    @sa imgMaskDilate, imgMaskErode
    @return Returns 0 when successful.
 */
int imgStructuringElement(
  /** Pointer to empty IMG struct for the element. */
  IMG *img,
  /** Structuring element: 
      1. 3x3x3 cube
      2. rounded 3x3x3 cube (23 voxels)
      3. cube on its corner (star, consisting of 7 voxels).
   */
  const int structuring_element,
  /** Verbose level; if zero, then only warnings are printed into stderr */
  int verbose
) {
  if(img==NULL) return(1);
  imgEmpty(img);

  int ret=0;
  if(structuring_element==1) {
    if(verbose>0) printf("making cube as the structuring element\n");
    ret=imgAllocate(img, 3, 3, 3, 1);
  } else if(structuring_element==2) {
    if(verbose>0) printf("making rounded cube as the structuring element\n");
    ret=imgAllocate(img, 3, 3, 3, 1);
  } else if(structuring_element==3) {
    if(verbose>0) printf("making star as the structuring element\n");
    ret=imgAllocate(img, 3, 3, 3, 1);
  } else {
    fprintf(stderr, "Error: unsupported structuring element.\n");
    return(2);
  }
  if(verbose>0) {fflush(stdout);}
  if(ret) {
    fprintf(stderr, "Error: cannot allocate memory.\n");
    imgEmpty(img);
    return(3);
  }

  int z, y, x, n;

  if(structuring_element==1) { // 3x3x3 cube
    for(z=0; z<3; z++)
      for(y=0; y<3; y++)
        for(x=0; x<3; x++)
          img->m[z][y][x][0]=1.0;
  } else if(structuring_element==2) { // rounded 3x3x3 cube
    for(z=0; z<3; z++)
      for(y=0; y<3; y++)
        for(x=0; x<3; x++)
          if(z==1 || y==1 || x==1)
            img->m[z][y][x][0]=1.0;
          else
            img->m[z][y][x][0]=0.0;
  } else if(structuring_element==3) { // 7-voxel star
    for(z=0; z<3; z++)
      for(y=0; y<3; y++)
        for(x=0; x<3; x++) {
          n=0; 
          if(z==1) n++;
          if(y==1) n++;
          if(x==1) n++;
          if(n>1)
            img->m[z][y][x][0]=1.0;
          else
            img->m[z][y][x][0]=0.0;
        }
  }

  if(verbose>3) {
    int z, y, x;
    printf("\nplanes 1-3\n");
    for(x=0; x<3; x++) {
      for(z=0; z<3; z++) {
        for(y=0; y<3; y++)
          printf(" %g", img->m[z][y][x][0]);
        printf("  ");
      }
      printf("\n");
    }
    printf("\n");
  }
  return(0);
}
/*****************************************************************************/

/*****************************************************************************/
/** Invert the 3D mask image, setting zeroes to ones, and non-zeroes to zeroes.

    Processes only the first frame, even if several do exist.

    @sa imgMaskCount, imgMaskRoiNr, imgThresholdMask
 */
void imgMaskInvert(
  /** Pointer to mask IMG struct. */
  IMG *img
) {
  if(img==NULL) return;
  int zi, yi, xi;
  for(zi=0; zi<img->dimz; zi++)
    for(yi=0; yi<img->dimy; yi++)
      for(xi=0; xi<img->dimx; xi++)
        if(fabs(img->m[zi][yi][xi][0])>1.0E-12)
          img->m[zi][yi][xi][0]=0.0;
        else
          img->m[zi][yi][xi][0]=1.0;
  return;
}
/*****************************************************************************/

/*****************************************************************************/
/** Conjunction (AND, wedge) for two 3D mask images.
    @sa imgMaskInv, imgMaskDilate, imgMaskCount
    @return Returns 0 when successful.
 */
int imgMaskConjunction(
  /** Pointer to IMG struct containg the first mask image;
      will be overwritten by the conjunction mask. */
  IMG *mask1,
  /** Pointer to IMG struct containing the second mask image; not modified. */
  IMG *mask2
) {
  if(mask1==NULL || mask2==NULL) return(1);
  if(mask1->dimx<1 || mask1->dimy<1 || mask1->dimz<1) return(2);
  if(mask1->dimx!=mask2->dimx) return(3); 
  if(mask1->dimy!=mask2->dimy) return(4); 
  if(mask1->dimz!=mask2->dimz) return(5); 

  int zi, yi, xi;
  for(zi=0; zi<mask1->dimz; zi++)
    for(yi=0; yi<mask1->dimy; yi++)
      for(xi=0; xi<mask1->dimx; xi++) {
        if(fabs(mask1->m[zi][yi][xi][0])<1.0E-12 ||
           fabs(mask2->m[zi][yi][xi][0])<1.0E-12)
        {
          mask1->m[zi][yi][xi][0]=0.0;
        } else {
          mask1->m[zi][yi][xi][0]=1.0;
        }
      }

  return(0);
}
/*****************************************************************************/

/*****************************************************************************/
/** Region labelling with flood filling.

    Processes only the first frame, even if several do exist.

    Based on Burger W and Burge MJ: Principles of Digital Image Processing -
    Core Algorithms, Springer, 2009, DOI 10.1007/978-1-84800-195-4.

    @return Returns 0 when successful.
 */
int imgMaskRegionLabeling(
  /** Pointer to IMG struct containg the mask image to be region labelled;
      not modified. All pixels with value < 0.1 are considered to belong to
      background, others to foreground. */
  IMG *mask1,
  /** Pointer to an empty but initiated IMG struct for the output, the region 
      labelled mask image. */
  IMG *mask2,
  /** The number of regions that were found; enter NULL, if not needed. */
  int *n,
  /** Verbose level; if zero, then nothing is printed into stdout or stderr */
  int verbose
) {
  if(verbose>0) {printf("imgMaskRegionLabeling()\n"); fflush(stdout);}
  if(mask1==NULL || mask2==NULL || mask1==mask2) return(1);
  if(mask1->status!=IMG_STATUS_OCCUPIED) return(2);
  if(mask1->dimx<1 || mask1->dimy<1 || mask1->dimz<1) return(3);
  if(verbose>1) printf("mask dimensions := %d x %d x %d\n", mask1->dimx, mask1->dimy, mask1->dimz);
  if(n!=NULL) *n=0;

  /* Make a copy of the mask */
  imgEmpty(mask2);
  int nr;
  if(imgThresholdMaskCount(mask1, 0.1, 1.0E+22, mask2, &nr)!=0) {
    if(verbose>0) fprintf(stderr, "Error: cannot make initial copy of mask.\n");
    return(11);
  }
  if(nr==0) {
    if(verbose>0) fprintf(stderr, "Warning: empty mask.\n");
    return(0);
  }
  if(verbose>1) printf("mask contains %d foreground pixels.\n", nr);
  if(verbose>80) {
    for(int zi=0; zi<mask2->dimz; zi++)
      for(int yi=0; yi<mask2->dimy; yi++)
        for(int xi=0; xi<mask2->dimx; xi++)
          if(mask2->m[zi][yi][xi][0]!=0.0)
            printf("  %d,%d,%d\n", zi, yi, xi);
  }

  /* Label regions */
  int ret=0, nextLabel=2;
  for(int zi=0; zi<mask2->dimz; zi++)
    for(int yi=0; yi<mask2->dimy; yi++)
      for(int xi=0; xi<mask2->dimx; xi++)
        if(mask2->m[zi][yi][xi][0]==1.0) {
          ret=imgMaskFloodFill(mask2, zi, yi, xi, nextLabel, &nr, verbose);
          if(ret!=0) break;
          if(verbose>2) printf("%d pixels labelled as %d\n", nr, nextLabel);
          nextLabel++;
        }
  if(ret!=0) {
    if(verbose>0) fprintf(stderr, "Error: Flood fill failed.\n");
    imgEmpty(mask2); return(21);
  }
  if(verbose>0) printf("%d regions labelled.\n", nextLabel-2);
  if(n!=NULL) *n=nextLabel-2;

  return(0);
}
/*****************************************************************************/

/*****************************************************************************/
/** Flood filling for the Region labelling.

    Processes only the first frame, even if several do exist.

    Based on Burger W and Burge MJ: Principles of Digital Image Processing -
    Core Algorithms, Springer, 2009, DOI 10.1007/978-1-84800-195-4.

    @sa imgMaskRegionLabeling
    @return Returns 0 when successful.
 */
int imgMaskFloodFill(
  /** Pointer to IMG struct containing the mask image to be flood filled. 
      Pixels with value 1 are flood filled with label, pixels with value <1
      are considered as background, and pixels with values >1 are considered
      to belong to another already labelled region. */
  IMG *m,
  /** Z coordinate of the starting point [0..dimz-1]. */
  int sz,
  /** Y coordinate of the starting point [0..dimy-1]. */
  int sy,
  /** X coordinate of the starting point [0..dimx-1]. */
  int sx,
  /** The label to be assigned to the region; at least 2. */
  int label,
  /** The number of pixels that got labelled; enter NULL, if not needed. */
  int *n,
  /** Verbose level; if zero, then nothing is printed into stdout or stderr */
  int verbose
) {
  if(verbose>0) {
    printf("imgMaskFloodFill(mask, %d, %d, %d, %d)\n", sz, sy, sx, label);
    fflush(stdout);
  }
  if(m==NULL || m->status!=IMG_STATUS_OCCUPIED) return(1);
  if(m->dimx<1 || m->dimy<1 || m->dimz<1) return(2);
  if(sx<0 || sy<0 || sz<0 || sx>=m->dimx || sy>=m->dimy || sz>=m->dimz) return(3);
  if(label<2) return(4);
  if(n!=NULL) *n=0;

  /* Create empty stack of pixel coordinates */
  IMG_PIXELS pxls; pxlInit(&pxls);

  /* Put the start coordinate into the stack */
  IMG_PIXEL pxl; pxl.x=sx; pxl.y=sy; pxl.z=sz;
  if(pxlAdd(&pxls, &pxl)!=0) {
    if(verbose>0) fprintf(stderr, "Error: cannot add seed pixel to stack.\n");
    pxlFree(&pxls); return(11);
  }

  /* Process the stack until empty */
  int pxlNr=0;
  while(pxls.pxlNr>0) {
    /* Take the last pixel from the list */
    if(pxlGet(&pxls, pxls.pxlNr-1, &pxl)!=0) break;
    pxls.pxlNr--;
    /* Should this pixel be labelled? */
    if(verbose>100) printf("  m[%d][%d][%d] := %g\n", pxl.z, pxl.y, pxl.x, m->m[pxl.z][pxl.y][pxl.x][0]);
    if(pxl.z<0 || pxl.z>=m->dimz || pxl.y<0 || pxl.y>=m->dimy || pxl.x<0 || pxl.x>=m->dimx)
      continue; // certainly not since outside of the image
    if(m->m[pxl.z][pxl.y][pxl.x][0]!=1.0) continue; // no since labelled or part of background
    m->m[pxl.z][pxl.y][pxl.x][0]=(float)label; // yes
    if(verbose>100) printf("  -> m[%d][%d][%d] := %g\n", pxl.z, pxl.y, pxl.x, m->m[pxl.z][pxl.y][pxl.x][0]);
    pxlNr++;
    /* Add the 26 neighbours to the beginning of the stack */
    if(pxlMakeRoom(&pxls, 0, 26)!=0) break;
    if(pxls.pxlNr==0) pxls.pxlNr=26;
    if(verbose>100) printf("  pxls.pxlNr := %d\n", pxls.pxlNr);
    int i=0;
    for(int dz=-1; dz<2; dz++) for(int dy=-1; dy<2; dy++) for(int dx=-1; dx<2; dx++) {
      if(dz==0 && dy==0 && dx==0) continue;
      pxls.p[i].z=pxl.z+dz; pxls.p[i].y=pxl.y+dy; pxls.p[i].x=pxl.x+dx;
      i++;
    }
  }
  pxlFree(&pxls);
  if(verbose>1) printf("  %d pixels labelled.\n", pxlNr);
  if(n!=NULL) *n=pxlNr;
  if(pxlNr==0) return(21);
  return(0);
}
/*****************************************************************************/

/*****************************************************************************/
