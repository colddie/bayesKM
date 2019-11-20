/// @file imgtiff.c
/// @author Vesa Oikonen
/// @brief Writing IMG data as a TIFF 6.0 format image.
///
/*****************************************************************************/
#include "libtpcimgp.h"
/*****************************************************************************/

/*****************************************************************************/
/** Write one frame or plane in IMG data as a TIFF 6.0 format image.
    Overwrites existing TIFF file.
\return Returns 0, if ok.
 */
int tiffWriteImg(
  /** IMG containing PET image/sinogram data */
  IMG *img,         
  /** matrix index of plane (0..dimz-1); all if <0 */
  int plane,        
  /** matrix index of frame (0..dimt-1); all if <0 */
  int frame,        
  /** colours are scaled between 0 and maxvalue;
      if <=0, then searches max and sets maxvalue */
  float *maxvalue,
  /** PET_GRAYSCALE, PET_GRAYSCALE_INV or PET_RAINBOW */
  int colorscale, 
  /** name of output TIFF file */  
  char *fname,    
  /** Nr of matrices tiled horizontally; enter 0 for automatic calculation */
  int matXdim,
  /** Nr of matrices tiled vertically; enter 0 for automatic calculation */
  int matYdim,
  /** Verbose level; if zero, then nothing is printed to stderr or stdout */
  int verbose,
  /** Pointer to a string (allocated for at least 64 chars) where error message
      or other execution status will be written; enter NULL, if not needed */     
  char *status
) {
  int          i, j, k, ri, ci, pi, fi, pxlNr, matNr;
  int          mi, mc, mr;
  FILE        *fp;
  short int    svar, svars[2048];
  char         buf[4096], *cdata, *cptr;
  int          ivar, ivars[1024];
  /* Rainbow color scale */
  struct { int   n,   r,  g,  b, dr, dg, db; }
  bitty[] = {  {32,   0,  0,  0,  2,  0,  4},   /* violet to indigo */
               {32,  64,  0,128, -2,  0,  4},   /* indigo to blue */
               {32,   0,  0,255,  0,  8, -8},   /* blue to green */
               {64,   0,255,  0,  4,  0,  0},   /* green to yellow */
               {32, 255,255,  0,  0, -2,  0},   /* yellow to orange */
               {64, 255,192,  0,  0, -3,  0} }; /* orange to red */

  if(verbose>0) printf("tiffWriteImg(*img, %d, %d, %g, %d, %s, status, %d)\n",
    plane, frame, *maxvalue, colorscale, fname, verbose);

  /* Check input data */
  if(status!=NULL) strcpy(status, "fault in callinroutine");
  if(img->status!=IMG_STATUS_OCCUPIED) return(1);
  if(img->dimt<(frame+1) || img->dimz<(plane+1)) return(2);
  pxlNr=img->dimx*img->dimy; if(pxlNr<1) return(3);
  if(status!=NULL) strcpy(status, "ok");

  /* If color scale maximum was not specified, then determine it here */
  if(*maxvalue<=0.0) {
    *maxvalue=-1.0e-12;
    for(pi=0; pi<img->dimz; pi++) if(plane<0 || plane==pi)
      for(ri=0; ri<img->dimy; ri++) for(ci=0; ci<img->dimx; ci++)
        for(fi=0; fi<img->dimt; fi++) if(frame<0 || frame==fi)
          if(img->m[pi][ri][ci][fi]>*maxvalue)
            *maxvalue=img->m[pi][ri][ci][fi];
    if(*maxvalue<=0.0) {
      if(status!=NULL) strcpy(status, "no positive pixel values");
      return(6);
    }
  }

  /* Calculate image dimensions */
  /* image matrix number */
  if(plane<0) matNr=img->dimz; else matNr=1;
  if(frame<0) matNr*=img->dimt; else matNr*=1;
  if(verbose>1) printf("matNr=%d\n", matNr);
  /* image matrix x*y number */
  if(matXdim<=0 && matYdim<=0) {
    matXdim=(int)ceil(sqrt((double)matNr));
    matYdim=matNr/matXdim; if(matNr%matXdim) matYdim++;
  } else {
    if(matXdim>matNr) {
      matXdim=matNr; matYdim=1;
    } else if(matYdim>matNr) {
      matYdim=matNr; matXdim=1;
    } else if(matXdim>0) {
      matYdim=matNr/matXdim; if(matNr%matXdim) matYdim++;
    } else {
      matXdim=matNr/matYdim; if(matNr%matYdim) matXdim++;
    }
  }
  if(verbose>1) printf("matXdim:=%d\nmatYdim:=%d\n", matXdim, matYdim);


  /* Open TIFF file */
  if((fp=fopen(fname, "wb")) == NULL) {
    if(status!=NULL) strcpy(status, "cannot open file for write");
    return(11);
  }

  /* Construct TIFF header */
  memset(buf, 0, 4096);
  /* set the byte format */
  if(little_endian()) memcpy(buf, "II", 2); else memcpy(buf, "MM", 2);
  /* set file identifier */
  svar=42; memcpy(buf+2, &svar, 2);
  /* set byte offset of first IFD */
  ivar=8; memcpy(buf+4, &ivar, 4);
  /* Construct the (first) Image File Directory (IFD) */
  /* set nr of directory entries */
  if(colorscale==PET_RAINBOW || colorscale==PET_RAINBOW_WB)
    svar=12; else svar=11;
  memcpy(buf+8, &svar, 2);
  /* move into start of first entry */
  cptr=buf+10;
  /* tag: ImageWidth */
  svars[0]=256; svars[1]=4; memcpy(cptr, svars, 4); cptr+=4;
  ivars[0]=1; ivars[1]=matXdim*img->dimx; memcpy(cptr, ivars, 8); cptr+=8;
  /* tag: ImageLength */
  svars[0]=257; svars[1]=4; memcpy(cptr, svars, 4); cptr+=4;
  ivars[0]=1; ivars[1]=matYdim*img->dimy; memcpy(cptr, ivars, 8); cptr+=8;
  /* tag: BitsPerSample (xv3 on Sun/Solaris gives warning but works) */
  svars[0]=258; svars[1]=3; memcpy(cptr, svars, 4); cptr+=4;
  ivars[0]=1; memcpy(cptr, ivars, 4); cptr+=4;
  svars[0]=(unsigned short int)8; memcpy(cptr, svars, 2); cptr+=4; /* 256 shades */
  /* tag: Compression */
  svars[0]=259; svars[1]=3; memcpy(cptr, svars, 4); cptr+=4;
  ivars[0]=1; memcpy(cptr, ivars, 4); cptr+=4;
  svars[0]=1; memcpy(cptr, svars, 2); cptr+=4; /* no compression */
  /* tag: Photometric Interpretation */
  svars[0]=262; svars[1]=3; memcpy(cptr, svars, 4); cptr+=4;
  ivars[0]=1; memcpy(cptr, ivars, 4); cptr+=4;
  if(colorscale==PET_RAINBOW || colorscale==PET_RAINBOW_WB)
    svars[0]=3; /* palette */
  else if(colorscale==PET_GRAYSCALE)
    svars[0]=1; /* black is zero */
  else
    svars[0]=0; /* white is zero */
  memcpy(cptr, svars, 2); cptr+=4;
  /* tag: StripOffsets */
  svars[0]=273; svars[1]=4; memcpy(cptr, svars, 4); cptr+=4;
  /* byte offset of strip(s) of data */
  ivars[0]=1; ivars[1]=4096; memcpy(cptr, ivars, 8); cptr+=8;
  /* tag: RowsPerStrip */
  svars[0]=278; svars[1]=4; memcpy(cptr, svars, 4); cptr+=4;
  ivars[0]=1; ivars[1]=matYdim*img->dimy; memcpy(cptr, ivars, 8); cptr+=8;
  /* tag: StripByteCounts */
  svars[0]=279; svars[1]=4; memcpy(cptr, svars, 4); cptr+=4;
  ivars[0]=1; ivars[1]=matXdim*matYdim*pxlNr; memcpy(cptr, ivars, 8); cptr+=8;
  /* tag: XResolution */
  ivars[0]=33; ivars[1]=1; ivars[2]=33; ivars[3]=1; memcpy(buf+1024, ivars, 16);
  svars[0]=282; svars[1]=5; memcpy(cptr, svars, 4); cptr+=4;
  ivars[0]=1; ivars[1]=1024; memcpy(cptr, ivars, 8); cptr+=8;
  /* tag: YResolution */
  svars[0]=283; svars[1]=5; memcpy(cptr, svars, 4); cptr+=4;
  ivars[0]=1; ivars[1]=1032; memcpy(cptr, ivars, 8); cptr+=8;
  /* tag: ResolutionUnit */
  svars[0]=296; svars[1]=3; memcpy(cptr, svars, 4); cptr+=4;
  ivars[0]=1; memcpy(cptr, ivars, 4); cptr+=4;
  svars[0]=3; memcpy(cptr, svars, 2); cptr+=4; /* cm */
  if(colorscale!=PET_RAINBOW && colorscale!=PET_RAINBOW_WB) {
    /* offset of the next IFD, or 0000 */
    for(i=0; i<4; i++) *cptr++=(char)0;
  } else {
    /* tag: ColorMap */
    svars[0]=320; svars[1]=3; memcpy(cptr, svars, 4); cptr+=4;
    ivars[0]=3*256; memcpy(cptr, ivars, 4); cptr+=4;
    //svars[0]=2048; memcpy(cptr, svars, 2); cptr+=4;
    ivars[0]=2048; memcpy(cptr, ivars, 4); cptr+=4;
    /* offset of the next IFD, or 0000 */
    for(i=0; i<4; i++) *cptr++=(char)0;
    /* Color table */
    cptr=buf+2048;
    /* red */
    for(i=0, j=0; j<6; j++) {
      svars[i++]=bitty[j].r;
      for(k=1; k<bitty[j].n; k++, i++) svars[i]=svars[i-1]+bitty[j].dr;
    }
    if(colorscale==PET_RAINBOW_WB) svars[0]=255;
    memcpy(cptr, svars, 512); cptr+=512;
    /* green */
    for(i=0, j=0; j<6; j++) {
      svars[i++]=bitty[j].g;
      for(k=1; k<bitty[j].n; k++, i++) svars[i]=svars[i-1]+bitty[j].dg;
    }
    if(colorscale==PET_RAINBOW_WB) svars[0]=255;
    memcpy(cptr, svars, 512); cptr+=512;
    /* blue */
    for(i=0, j=0; j<6; j++) {
      svars[i++]=bitty[j].b;
      for(k=1; k<bitty[j].n; k++, i++) svars[i]=svars[i-1]+bitty[j].db;
    }
    if(colorscale==PET_RAINBOW_WB) svars[0]=255;
    memcpy(cptr, svars, 512); cptr+=512;
  }

  /* write the IFD */
  if(fwrite(buf, 1, 4096, fp) != 4096) {
    fclose(fp); remove(fname);
    if(status!=NULL) strcpy(status, "cannot write file");
    return(13);
  }

  /* Write pixel data */
  cdata=(char*)calloc(matXdim*matYdim*pxlNr, sizeof(char));
  if(cdata==NULL) {
    fclose(fp); remove(fname);
    if(status!=NULL) strcpy(status, "out of memory");
    return(14);
  }
  cptr=cdata; i=0;
  mi=0; mc=0; mr=1;
  for(fi=0; fi<img->dimt; fi++) if(frame<0 || frame==fi) {
    for(pi=0; pi<img->dimz; pi++) if(plane<0 || plane==pi) {
      mi++; mc++; 
      for(ri=0; ri<img->dimy; ri++) for(ci=0; ci<img->dimx; ci++) {
        cptr=cdata + (mr-1)*matXdim*pxlNr + ri*matXdim*img->dimx
                   + (mc-1)*img->dimx + ci;
        if(img->m[pi][ri][ci][fi]>0.0) {
          if((img->m[pi][ri][ci][fi])<(*maxvalue))
            *cptr=(unsigned char)(255.*(img->m[pi][ri][ci][fi])/(*maxvalue));
          else *cptr=(unsigned char)255;
        } else
          *cptr=(unsigned char)0;
      }
      if(mc==matXdim) {mc=0; mr++;}
    }
  }
  cptr=cdata;
  if(fwrite(cptr, 1, matXdim*matYdim*pxlNr, fp) != 
      (unsigned int)matXdim*matYdim*pxlNr) {
    fclose(fp); remove(fname); free(cdata);
    if(status!=NULL) strcpy(status, "cannot write file");
    return(15);
  }
  free(cdata);

  fclose(fp);
  if(status!=NULL) strcpy(status, "ok");
  return(0);
}
/*****************************************************************************/

/*****************************************************************************/

