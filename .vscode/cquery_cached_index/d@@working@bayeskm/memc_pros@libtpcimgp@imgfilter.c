/// @file imgfilter.c
/// @author Calle Laakkonen, Kaisa Sederholm, Vesa Oikonen
/// @brief Gaussian IMG filter.
///
#include "libtpcimgp.h"
/*****************************************************************************/

/*****************************************************************************/
/** Make a 2D Gaussian convolution kernel.

    @sa imgConvolute2D, imgFast3DGaussianFilter, imgFast2DGaussianFilter
    @return Returns 0 when successful, otherwise <>0.
*/
int imgFillGaussKernel(
  /** Gaussian convolution kernel matrix[size][size], filled here. */
  float **kernel,
  /** Gaussian S.D. in pixels (decimals are ok). */
  float stdev,
  /** Kernel dimension. */
  int size
) {
  int x, y;
  float mx, my, v, ksum=0.0;

  if(kernel==NULL || size<3) return 1;
  if(stdev<0.0) return 2;
  if((size%2)==0) return 3; // size must be odd number

  v=stdev*stdev;
  for(x=0; x<size; x++) {
    mx = x - size/2; // note: ints on purpose
    for(y=0; y<size; y++) {
      my = y - size/2; // note: ints on purpose
      if(stdev>0) {
#if(0)  // Gaussian in the middle of pixel is calculated
        kernel[y][x] =
          (1.0/(2.0*M_PI*v)) * powf(M_E,-((mx*mx)+(my*my))/(2.0*v));
#else  // Pixel area is considered
        kernel[y][x] = 0.25
              *(erff((mx-0.5)/(stdev*M_SQRT2))-erff((mx+0.5)/(stdev*M_SQRT2)))
              *(erff((my-0.5)/(stdev*M_SQRT2))-erff((my+0.5)/(stdev*M_SQRT2)));
#endif
      } else {
        if(x==size/2 && y==size/2) kernel[y][x]=1.0; else kernel[y][x]=0.0;
      }
      ksum+=kernel[y][x];
    }
  }

  /* Ensure quantitativity with normalization;
     divide each kernel value by their sum */
  v=1.0/ksum;
  for(x=0; x<size; x++) for(y=0; y<size; y++) kernel[y][x]*=v;

  return 0;
}
/*****************************************************************************/

/*****************************************************************************/
/** Perform a convolution operation on one image matrix.
    @remark Note that this may be a very slow process if image and kernel size is
    large: width*height*size*size multiplications will be done.

    @sa imgFast3DGaussianFilter, imgFast2DGaussianFilter
    @return Returns 0 when successful, otherwise <>0.
 */
int imgConvolute2D(
  /** Image matrix data[height][width][frame]. */
  float ***data,
  /** Temporary preallocated memory buffer[size][width+size-1]. */
  float **buffer,
  /** frame [0..dimt] of data matrix which is processed. */
  int frame,
  /** Width of image matrix (dimx). */
  int width,
  /** Height of image matrix (dimy). */
  int height,
  /** Convolution kernel[size][size] matrix. */
  float **kernel,
  /** Convolution kernel size; must be an odd number, and >=3, and
      smaller than 1 + 2*image_dimension. */
  int size,
  /** Fill borders with zero (0) or with closest image pixel (<>0). */
  int border,
  /** Verbose level; if zero, then only warnings are printed into stderr. */
  int verbose,
  /** Pointer to error message, at least 128 characters; NULL, if not needed. */
  char *errmsg 
) {
  const int size2=size/2, bufw=width+2*size2;
  int x, y, i;
  float *tmp;
  float filtered;
  int kx, ky;

  if(verbose>0) 
    printf("imgConvolute2D(***data, **buf, %d, %d, %d, **ker, %d, %d, %d)\n", 
           frame, width, height, size, border, verbose);
  /* Check input */
  if(errmsg!=NULL) strcpy(errmsg, "invalid function input");
  if(data==NULL || buffer==NULL || frame<0 || width<3 || height<3 || kernel==NULL)
    return 1;
  if((size%2)==0) {
    if(errmsg!=NULL) strcpy(errmsg, "kernel size is even number");
    return 1;
  }
  if(size<3) {
    if(errmsg!=NULL) strcpy(errmsg, "kernel size must be >=3");
    return 1;
  } else {
    i=width<height ? width : height;
    if(size>=2*i+1) {
      if(errmsg!=NULL) strcpy(errmsg, "kernel size is too big for image");
      return 1;
    }
  }
  if(verbose>1) {
    printf("  size2 := %d\n", size2);
    printf("  bufw := %d\n", bufw);
  }

  /* Fill the initial buffer */
  if(verbose>2) printf("filling initial data buffer\n"); 
  if(border==0) {
    for(y=0; y<size2; y++) {
      for(x=0; x<bufw; x++) buffer[y][x]=0;
    }
    for(y=size2; y<size; y++) {
      for(x=0; x<size2; x++) buffer[y][x]=0;
      for(x=size2; x<bufw-size2; x++) buffer[y][x]=data[y-size2][x-size2][frame];
      for(x=bufw-size2; x<bufw; x++) buffer[y][x]=0;
    }
  } else {
    for(y=0; y<size2; y++) {
      for(x=0; x<size2; x++) buffer[y][x]=data[0][0][frame];
      for(x=size2; x<bufw-size2; x++) buffer[y][x]=data[0][x-size2][frame];
      for(x=bufw-size2; x<bufw; x++) buffer[y][x]=data[0][width-1][frame];
    }
    for(y=size2; y<size; y++) {
      for(x=0; x<size2; x++) buffer[y][x]=data[0][y-size2][frame];
      for(x=size2; x<bufw-size2; x++) buffer[y][x]=data[y-size2][x-size2][frame];
      for(x=bufw-size2; x<bufw; x++) buffer[y][x]=data[y-size2][width-1][frame];
    }
  }
  if(verbose>8) {
    int xi, yi;
    printf("Initial buffer for row-index 0:\n");
    for(yi=0; yi<size; yi++) {
      for(xi=0; xi<bufw; xi++) printf(" %4.0f", buffer[yi][xi]);
      printf("\n");
    }
  }

  /* Filter */
  for(y=0; y<height; y++) {
    if(verbose>2) printf("filtering image row %d\n", y+1); 

    /* Filter a row */
    for(x=0; x<width; x++) {
      filtered = 0.0;
      for(ky=0; ky<size; ky++) {
        for(kx=0; kx<size; kx++) {
          filtered += buffer[ky][x+kx] * kernel[ky][kx];
        }
      }
      data[y][x][frame] = filtered;
    }

    if(y==height-1) break; // no need for further buffer filling

    if(verbose>2) printf("filling data buffer for next image row\n"); 
    /* Shift rows in data buffer */
    tmp = buffer[0];
    for(i=1; i<size; i++) buffer[i-1]=buffer[i];
    buffer[size-1]=tmp;

    /* Get next row from image or 'guess' it */
    if(y+size2+1<height) { // still inside original image 
      if(border==0) {
        for(x=0; x<size2; x++)
          buffer[size-1][x]=0;
        for(x=size2; x<bufw-size2; x++)
          buffer[size-1][x]=data[y+size2+1][x-size2][frame];
        for(x=bufw-size2; x<bufw; x++)
          buffer[size-1][x]=0;
      } else {
        for(x=0; x<size2; x++)
          buffer[size-1][x]=data[y+size2+1][0][frame];
        for(x=size2; x<bufw-size2; x++)
          buffer[size-1][x]=data[y+size2+1][x-size2][frame];
        for(x=bufw-size2; x<bufw; x++)
          buffer[size-1][x]=data[y+size2+1][width-1][frame];
      }
    } else { // below original image
      if(border==0) {
        for(x=0; x<bufw; x++) buffer[size-1][x] = 0;
      } else {
        for(x=0; x<size2; x++)
          buffer[size-1][x]=data[height-1][0][frame];
        for(x=size2; x<bufw-size2; x++)
          buffer[size-1][x]=data[height-1][x-size2][frame];
        for(x=bufw-size2; x<bufw; x++)
          buffer[size-1][x]=data[height-1][width-1][frame];
      }
    }
    if(verbose>9) {
      int xi, yi;
      printf("Buffer for row-index %d:\n", y+1);
      for(yi=0; yi<size; yi++) {
        for(xi=0; xi<bufw; xi++) printf(" %4.0f", buffer[yi][xi]);
        printf("\n");
      }
    }
  }
  return 0;
}
/*****************************************************************************/

/*****************************************************************************/
/** Apply 2D Gaussian filter to whole dynamic image in IMG struct.
    @sa imgFast3DGaussianFilter, imgFast2DGaussianFilter, imgsegmSimilar
    @return If an error is encountered, function returns a non-zero value. 
            Otherwise 0 is returned.
 */
int imgGaussianFilter(
  /** Image data to be processed; data is overwritten with filtered image. */
  IMG *img,
  /** Plane index [0..dimz-1]; enter <0 to filter all image planes. */
  int plane,
  /** Frame index [0..dimt-1]; enter <0 to filter all image time frames. */
  int frame,
  /** Gaussian S.D. in pixels (decimals are ok). */
  float gauss_sd,
  /** Gaussian kernel size in pixels; it must be an odd number and at least 3
      and smaller than 2*imgdim+1; typically 6*gauss_sd is sufficient;
      enter zero to set it automatically to 6*gauss_sd (plus 1 if necessary). */
  int size,
  /** In filtering, fill borders with zero (0) or with closest image pixel (<>0). */
  int border,
  /** Verbose level; if zero, then only warnings are printed into stderr. */
  int verbose,
  /** Pointer to error message, at least 128 characters; NULL, if not needed. */
  char *errmsg 
) {
  int i, j, zi, fi, mindim, ret;
  float **buffer, *bdata;
  float **gauss, *gdata;

  if(verbose>0)
    printf("imgGaussianFilter(*img, %d, %d, %g, %d, %d, %d, msg)\n",
           plane, frame, gauss_sd, size, border, verbose); 

  /* Check input */
  if(errmsg!=NULL) strcpy(errmsg, "invalid function input");
  /* x and y dim must be >=2 so that kernel size can be at least 3 */
  if(img==NULL || img->dimz<1) return 1;
  if(plane>=0 && plane>=img->dimz) return 2;
  if(frame>=0 && frame>=img->dimt) return 3;
  mindim=img->dimx<img->dimy ? img->dimx : img->dimy;
  if(mindim<3) return 1;
  if(gauss_sd<0.0) return 2; // SD=0 is accepted for basic testing
  if(size!=0 && (size<3 || (size%2)==0) ) return 4;
  if(size==0) {
    size=6*roundf(gauss_sd); if(size<3) size=3;
    if(size>=2*mindim+1) size=2*mindim-1; // limit kernel size
    if((size%2)==0) size++; // make sure that it is odd number
    if(verbose>0) printf("size := %d\n", size);
  }

  /* Allocate memory for the kernel */
  if(errmsg!=NULL) strcpy(errmsg, "cannot allocate memory for kernel");
  if(verbose>2) printf("allocating memory for the kernel\n");
  gdata=(float*)malloc( size*size*sizeof(float) );
  if(gdata==NULL) {return 7;}
  gauss=(float**)malloc(size*sizeof(float*));
  if(gauss==NULL) {free(gdata); return 7;}
  for(i=0; i<size; i++) gauss[i]=gdata+(i*size);

  /* Allocate memory for temporary buffer[size][width+size-1] */
  if(errmsg!=NULL) strcpy(errmsg, "cannot allocate memory for data buffer");
  if(verbose>2) printf("allocating memory for the raw buffer\n");
  bdata=(float*)malloc( (img->dimx+size-1)*size*sizeof(float) );
  if(bdata==NULL) {free(gauss); free(gdata); return 8;}
  if(verbose>2) printf("allocating memory for the buffer\n");
  buffer=(float**)malloc(size*sizeof(float*));
  if(buffer==NULL) {free(gauss); free(gdata); free(bdata); return 8;}
  if(verbose>2) printf("preparing buffer\n");
  for(i=0; i<size; i++) buffer[i]=bdata+(i*(img->dimx+size-1));

  /* Make a gaussian convolution kernel */
  if(verbose>1) printf("calculating Gaussian convolution kernel\n");
  if(imgFillGaussKernel(gauss, gauss_sd, size) != 0) {
    if(errmsg!=NULL) strcpy(errmsg, "cannot compute Gaussian kernel");
    free(gauss); free(gdata); free(buffer); free(bdata);
    return 9;
  }
  if(verbose>4) {
    printf("Gaussian convolution kernel:\n");
    for(i=0; i<size; i++) {
      printf(" ");
      for(j=0; j<size; j++) printf(" %8.6f", gauss[i][j]);
      printf("\n");
    }
  }

  /* Convolute each image matrix */
  if(verbose>0) printf("convolution of the image data\n");
  for(zi=0; zi<img->dimz; zi++) {
    if(plane>=0 && plane!=zi) continue;
    if(verbose>1) printf("  plane %d\n", zi+1);
    for(fi=0; fi<img->dimt; fi++) {
      if(frame>=0 && frame!=fi) continue;
      if(verbose>1) printf("    frame %d\n", fi+1);
      ret=imgConvolute2D(img->m[zi], buffer, fi, img->dimx, img->dimy,
                         gauss, size, border, verbose, errmsg);
      if(ret!=0) {
        free(gauss); free(gdata); free(buffer); free(bdata);
        return(10+ret);
      }
    }
  }

  /* Free memory */
  free(gauss); free(gdata); free(buffer); free(bdata);

  if(verbose>1) printf("  convolution done.\n");
  return 0;
}
/*****************************************************************************/

/*****************************************************************************/
/// @cond
/** 
  \brief Allocate memory for a float matrix
  
  The matrix will be in form matrix[h][w]

  \param w width of matrix
  \param h height of matrix
  \return the allocated memory. NULL on error
  */
float **mallocMatrix(float w,float h)
{
  int y;
  float **matrix = malloc(h * sizeof(float*));
  if(!matrix) return NULL;
  for(y=0;y<h;y++) {
    /* Allocate a row */
    matrix[y] = malloc(w*sizeof(float));

    if(!matrix[y]) {
      while(y>0) free(matrix[--y]);
      free(matrix);
      return NULL;
    }
  }
  return matrix;
}
/// @endcond
/*****************************************************************************/

/*****************************************************************************/
/** \brief Make a Gaussian convolution kernel
  \param size kernel size
  \return the convolution kernel. NULL on error
*/
float **imgGaussKernel(
  int size
) {
  float **kernel, stdev, mx, my;
  int x,y;

  kernel = mallocMatrix(size,size);
  if(!kernel) return NULL;

  stdev = size/7.0;
  stdev = stdev*stdev;

  for(x=0;x<size;x++) {
    mx = x - size/2;
    for(y=0;y<size;y++) {
      my = y - size/2;
      kernel[x][y] =
       (1.0/(2*M_PI*stdev))*pow(M_E,-((mx*mx)+(my*my))/(2.0*stdev));
    }
  }
  return kernel;
}
/*****************************************************************************/

/*****************************************************************************/
/** \brief Free a convolution kernel
  \param kernel the convolution kernel
  \param size size of the kernel
  */
void imgFreeKernel(
  float **kernel,
  int size
) {
  int r;
  for(r=0;r<size;r++) free(kernel[r]);
  free(kernel);
}
/*****************************************************************************/

/*****************************************************************************/
/** @brief Perform a convolution operation on float data
    @sa imgConvolute
 */
void imgConvoluteData(
  /** Pointer to data in format float[y][x][frame]. */
  float ***data,
  /** Temporary buffer to use. Size must be float[size][data width+size]. */
  float **buffer,
  /** Frame to extract from data. */
  int frame,
  /** Width of the image (xdim). */
  int width,
  /** Height of the image (ydim). */
  int height,
  /** 2D convolution kernel. */
  float **kernel,
  /** Convolution kernel size. */
  int size
) {
  const int size2 = size/2;
  int x,y;

  /* Fill the initial buffer */
  for(y=0;y<size2;y++) {
    /* Pixels above the border */
    for(x=0;x<width+size;x++) buffer[y][x] = 0;

    /* Padding for the borders */
    for(x=0;x<size2;x++) {
      buffer[y+size2][x] = 0;
      buffer[y+size2][x+width+size2] = 0;
    }

    /* Rows from the image */
    for(x=0;x<width;x++)
      buffer[y+size2][x+size2] = data[y][x][frame];
  }

  /* Filter */
  for(y=0;y<height;y++) {
    float *tmp;
    /* Filter a row */
    for(x=0;x<width;x++) {
      float filtered = 0;
      int kx,ky;
      for(kx=0;kx<size;kx++) {
        for(ky=0;ky<size;ky++) {
          filtered += buffer[ky][x+kx] * kernel[kx][ky];
        }
      }
      data[y][x][frame] = filtered;
    }

    /* Shift rows */
    tmp = buffer[size-1];
    buffer[size-1] = buffer[0];
    for(x=0;x<size-2;x++) buffer[x] = buffer[x+1];
    buffer[size-2] = tmp;

    /* Get next row from image or blank */
    if(y+size2<height) {
      for(x=0;x<width;x++) buffer[size-1][x+size2] = data[y+size2][x][frame];
    } else {
      for(x=0;x<width;x++) buffer[size-1][x+size2] = 0;
    }
  }
}
/*****************************************************************************/

/*****************************************************************************/
/** @brief Perform a convolution operation on image

  The bordering size/2 pixels of the image do not containing valid data.
  imgConvolute is quite memory efficient.
  Only (img->dimx+size)*size*sizeof(float) bytes are consumed by a temporary buffer.
  When processing multiple planes/frames/images with equal width and kernel size, efficiency 
  can be improved further by allocating the temporary buffer yourself and calling
  imgConvoluteData directly.

  @sa imgConvoluteData
  @return nonzero in case of (malloc) error.
*/
int imgConvolute(
  /** Image data to be processed; data is overwritten with filtered image. */
  IMG *img,
  /** Time frame index [0..dimt-1]. */
  int frame,
  /** Plane index [0..dimz-1]. */
  int plane,
  /** Convolution kernel. */
  float **kernel,
  /** Convolution kernel size; size must be smaller than img->dimx and img->dimy. */
  int size
) {
  float **buffer;
  int y;

  /* Allocate memory for temporary buffer */
  buffer = mallocMatrix(img->dimx+size,size);
  if(!buffer) return 1;

  /* Convolute */
  imgConvoluteData(img->m[plane],buffer,frame,img->dimx, img->dimy, kernel,size);

  /* Free buffer */
  for(y=0;y<size;++y) free(buffer[y]);
  free(buffer);
  return 0;
}
/*****************************************************************************/

/*****************************************************************************/
/** Apply fast approximate 2D Gaussian filter to whole dynamic image in IMG struct.

    This function implements the fast Gaussian convolution algorithm with IIR approximation 
    (Alvarez L and Mazorra L, SIAM Journal on Numerical Analysis, 1994;31(2):590-605), and 
    is based on C code written by Pascal Getreuer <http://www.getreuer.info/home/gaussianiir>.

    @sa imgFast3DGaussianFilter, imgsegmSimilar
    @return If an error is encountered, function returns a non-zero value. 
        Otherwise 0 is returned.
 */
int imgFast2DGaussianFilter(
  /** Image data to be processed; data is overwritten with filtered image. */
  IMG *img,
  /** Plane index [0..dimz-1]; enter <0 to filter all image planes. */
  int plane,
  /** Frame index [0..dimt-1]; enter <0 to filter all image time frames. */
  int frame,
  /** Gaussian S.D. in pixels (decimals are ok); SD=FWHM/2.355. */
  float gauss_sd,
  /** Number of time steps. More steps implies better accuracy and longer execution time; 
      enter 0 to use the default steps (4); too many steps may lead to floating-point problems. */
  int step_nr,
  /** Verbose level; if zero, then only warnings are printed into stderr. */
  int verbose,
  /** Pointer to error message, at least 128 characters; NULL, if not needed. */
  char *errmsg 
) {
  int i, xi, yi, zi, fi, mindim, step;
  double lambda, nu, boundaryscale, postscale;
  double *ptr, *dimg;


  if(verbose>0)
    printf("imgFast2DGaussianFilter(*img, %d, %d, %g, %d, %d, msg)\n",
           plane, frame, gauss_sd, step_nr, verbose); 

  /* Check input */
  if(errmsg!=NULL) strcpy(errmsg, "invalid function input");
  if(gauss_sd==0.0) return 0;
  if(img==NULL || img->dimz<1) return 1;
  if(plane>=0 && plane>=img->dimz) return 2;
  if(frame>=0 && frame>=img->dimt) return 3;
  mindim=img->dimx<img->dimy ? img->dimx : img->dimy;
  if(mindim<3) return 1;
  if(gauss_sd<0.0) return 4; // SD=0 is accepted for basic testing
  /* Check step_nr */
  if(step_nr<=0) step_nr=4;

  /* Prepare for filtering */
  lambda= (gauss_sd*gauss_sd)/(2.0*step_nr);
  nu= (1.0 + 2.0*lambda - sqrt(1.0 + 4.0*lambda))/(2.0*lambda);
  boundaryscale= 1.0/(1.0-nu);
  postscale = pow(nu/lambda, 2*step_nr);
  if(verbose>1) {
    printf("nu := %g\n", nu);
    printf("boundaryscale := %g\n", boundaryscale);
    printf("postscale := %g\n", postscale);
  }

  /* Setup memory for a image plane */
  if(verbose>1) printf("Allocating memory\n");
  dimg=(double*)malloc(img->dimx*img->dimy*sizeof(double));
  if(dimg==NULL) {if(errmsg!=NULL) strcpy(errmsg, "out of memory"); return 5;}

  /* Process the required image planes and frames */
  if(verbose>1) printf("Gaussian filtering...\n");
  for(zi=0; zi<img->dimz; zi++) {
    if(plane>=0 && plane!=zi) continue;
    if(verbose>2) printf("  plane %d\n", zi+1);
    for(fi=0; fi<img->dimt; fi++) {
      if(frame>=0 && frame!=fi) continue;
      if(verbose>2) printf("    frame %d\n", fi+1);
      /* Copy image plane data to double buffer */
      for(yi=i=0; yi<img->dimy; yi++)
        for(xi=0; xi<img->dimx; xi++)
          dimg[i++]=img->m[zi][yi][xi][fi];
      /* Filter horizontally along each image row */
      for(yi=0; yi<img->dimy; yi++) {
        for(step=0; step<step_nr; step++) {
          ptr=dimg+img->dimx*yi; ptr[0]*=boundaryscale;
          /* Filter rightwards */
          for(xi=1; xi<img->dimx; xi++) ptr[xi]+=nu*ptr[xi-1];
          ptr[xi=img->dimx-1]*=boundaryscale;
          /* Filter leftwards */
          for(; xi>0; xi--) ptr[xi-1]+=nu*ptr[xi];
        }
      }
      /* Filter vertically along each image column */
      for(xi=0; xi<img->dimx; xi++) {
        for(step=0; step<step_nr; step++) {
          ptr=dimg+xi; ptr[0]*=boundaryscale;
          /* Filter downwards */
          for(i=img->dimx; i<img->dimx*img->dimy; i+=img->dimx)
            ptr[i]+=nu*ptr[i-img->dimx];
          ptr[i=img->dimx*(img->dimy-1)]*=boundaryscale;
          /* Filter upwards */
          for(; i>0; i-=img->dimx) ptr[i-img->dimx]+=nu*ptr[i];
        }
      }
      /* Copy and scale filtered image plane data back from buffer */
      for(yi=i=0; yi<img->dimy; yi++)
        for(xi=0; xi<img->dimx; xi++)
          img->m[zi][yi][xi][fi]=postscale*dimg[i++];
    } // next frame
  } // next plane

  /* Free memory */
  free(dimg);

  if(errmsg!=NULL) strcpy(errmsg, "ok");
  if(verbose>1) printf("  Gaussian convolution done.\n");
  return 0;
}
/*****************************************************************************/

/*****************************************************************************/
/** Apply fast approximate 3D Gaussian filter to whole dynamic image in IMG struct.

    This function implements the fast Gaussian convolution algorithm with IIR approximation 
    (Alvarez L and Mazorra L, SIAM Journal on Numerical Analysis, 1994;31(2):590-605), and 
    is based on C code written by Pascal Getreuer <http://www.getreuer.info/home/gaussianiir>.

    @sa imgFast2DGaussianFilter, imgThresholdingLowHigh, imgSmoothMax
    @return If an error is encountered, function returns a non-zero value. 
        Otherwise 0 is returned.
 */
int imgFast3DGaussianFilter(
  /** Image data to be processed; data is overwritten with filtered image;
      image must contain at least 3 planes. */
  IMG *img,
  /** Frame index [0..dimt-1]; enter <0 to filter all image time frames. */
  int frame,
  /** Gaussian S.D. in pixels (decimals are ok); note that the same S.D. is assumed to apply to 
      all dimensions, including Z dimension. SD=FWHM/2.355. */
  float gauss_sd,
  /** Number of time steps. More steps implies better accuracy and longer execution time; 
      enter 0 to use the default steps (4); too many steps may lead to floating-point problems. */
  int step_nr,
  /** Verbose level; if zero, then only warnings are printed into stderr. */
  int verbose,
  /** Pointer to error message, at least 128 characters; NULL, if not needed. */
  char *errmsg 
) {
  int i, xi, yi, zi, fi, mindim, step;
  double lambda, nu, boundaryscale, postscale;
  double *ptr, *dimg;


  if(verbose>0)
    printf("imgFast3DGaussianFilter(*img, %d, %g, %d, %d, msg)\n",
           frame, gauss_sd, step_nr, verbose); 

  /* Check input */
  if(errmsg!=NULL) strcpy(errmsg, "invalid function input");
  if(gauss_sd==0.0) return 0;
  if(img==NULL || img->dimz<3) return 1;
  if(frame>=0 && frame>=img->dimt) return 3;
  mindim=img->dimx<img->dimy ? img->dimx : img->dimy;
  if(mindim<3) return 1;
  if(gauss_sd<0.0) return 4; // SD=0 is accepted for basic testing
  /* Check step_nr */
  if(step_nr<=0) step_nr=4;

  /* Prepare for filtering */
  lambda= (gauss_sd*gauss_sd)/(2.0*step_nr);
  nu= (1.0 + 2.0*lambda - sqrt(1.0 + 4.0*lambda))/(2.0*lambda);
  boundaryscale= 1.0/(1.0-nu);
  postscale = pow(nu/lambda, 3*step_nr);
  if(verbose>1) {
    printf("nu := %g\n", nu);
    printf("boundaryscale := %g\n", boundaryscale);
    printf("postscale := %g\n", postscale);
  }

  /* Setup memory for a image frame (3D matrix) */
  if(verbose>1) printf("Allocating memory\n");
  dimg=(double*)malloc(img->dimz*img->dimx*img->dimy*sizeof(double));
  if(dimg==NULL) {if(errmsg!=NULL) strcpy(errmsg, "out of memory"); return 5;}

  /* Process the required image frame(s) */
  if(verbose>1) printf("Gaussian filtering...\n");
  for(fi=0; fi<img->dimt; fi++) {
    if(frame>=0 && frame!=fi) continue;
    if(verbose>2) printf("  frame %d\n", fi+1);
    /* Copy image 3D matrix data to double buffer */
    for(zi=i=0; zi<img->dimz; zi++)
      for(yi=0; yi<img->dimy; yi++)
        for(xi=0; xi<img->dimx; xi++)
          dimg[i++]=img->m[zi][yi][xi][fi];
    /* Filter horizontally along each image row */
    for(zi=0; zi<img->dimz; zi++) {
      for(yi=0; yi<img->dimy; yi++) {
        for(step=0; step<step_nr; step++) {
          ptr = dimg+img->dimx*(yi+img->dimy*zi); ptr[0]*=boundaryscale;
          /* Filter rightwards */
          for(xi=1; xi<img->dimx; xi++) ptr[xi]+=nu*ptr[xi-1];
          ptr[xi=img->dimx-1]*=boundaryscale;
          /* Filter leftwards */
          for(; xi>0; xi--) ptr[xi-1]+=nu*ptr[xi];
        }
      }
    }
    /* Filter vertically along each image column */
    for(zi=0; zi<img->dimz; zi++) {
      for(xi=0; xi<img->dimx; xi++) {
        for(step=0; step<step_nr; step++) {
          ptr=dimg+xi+img->dimx*img->dimy*zi; ptr[0]*=boundaryscale;
          /* Filter downwards */
          for(i=img->dimx; i<img->dimx*img->dimy; i+=img->dimx)
            ptr[i]+=nu*ptr[i-img->dimx];
          ptr[i=img->dimx*(img->dimy-1)]*=boundaryscale;
          /* Filter upwards */
          for(; i>0; i-=img->dimx) ptr[i-img->dimx]+=nu*ptr[i];
        }
      }
    }
    /* Filter along image z-dimension */
    for(yi=0; yi<img->dimy; yi++) {
      for(xi=0; xi<img->dimx; xi++) {
        for(step=0; step<step_nr; step++) {
          ptr=dimg+xi+img->dimx*yi; ptr[0]*=boundaryscale;
          for(i=img->dimx*img->dimy; i<img->dimx*img->dimy*img->dimz; i+=img->dimx*img->dimy)
            ptr[i]+=nu*ptr[i-img->dimx*img->dimy];
          ptr[i=img->dimx*img->dimy*(img->dimz-1)]*=boundaryscale;
          for(; i>0; i-=img->dimx*img->dimy)
            ptr[i-img->dimx*img->dimy]+=nu*ptr[i];
        }
      }
    }
    /* Copy and scale image 3D matrix data back from the buffer */
    for(zi=i=0; zi<img->dimz; zi++)
      for(yi=0; yi<img->dimy; yi++)
        for(xi=0; xi<img->dimx; xi++)
          img->m[zi][yi][xi][fi]=postscale*dimg[i++];
  } // next frame

  /* Free memory */
  free(dimg);

  if(errmsg!=NULL) strcpy(errmsg, "ok");
  if(verbose>1) printf("  Gaussian convolution done.\n");
  return 0;
}
/*****************************************************************************/

/*****************************************************************************/
/** Apply fast approximate 1D Gaussian filter over planes (z-axis) to whole dynamic image 
    in IMG struct.

    This function implements the fast Gaussian convolution algorithm with IIR approximation 
    (Alvarez L and Mazorra L, SIAM Journal on Numerical Analysis, 1994;31(2):590-605), and 
    is based on C code written by Pascal Getreuer <http://www.getreuer.info/home/gaussianiir>.

    @sa imgFast3DGaussianFilter, imgFast2DGaussianFilter
    @return If an error is encountered, function returns a non-zero value. 
        Otherwise 0 is returned.
 */
int imgFast1DGaussianFilter(
  /** Image data to be processed; data is overwritten with filtered image;
      Plane number (z-dimension) must be at least 3. */
  IMG *img,
  /** Gaussian S.D. in pixels (decimals are ok); SD=FWHM/2.355. */
  float gauss_sd,
  /** Number of time steps. More steps implies better accuracy and longer execution time; 
      enter 0 to use the default steps (4); too many steps may lead to floating-point problems. */
  int step_nr,
  /** Verbose level; if zero, then only warnings are printed into stderr. */
  int verbose,
  /** Pointer to error message, at least 128 characters; NULL, if not needed. */
  char *errmsg 
) {
  double lambda, nu, boundaryscale, postscale;


  if(verbose>0)
    printf("imgFast1DGaussianFilter(*img, %g, %d, %d, msg)\n", gauss_sd, step_nr, verbose);

  /* Check input */
  if(errmsg!=NULL) strcpy(errmsg, "invalid function input");
  if(gauss_sd==0.0) return 0;
  if(img==NULL || img->dimz<3) return 1;
  if(gauss_sd<0.0) return 4; // SD=0 is accepted for basic testing
  /* Check step_nr */
  if(step_nr<=0) step_nr=4;

  /* Prepare for filtering */
  lambda= (gauss_sd*gauss_sd)/(2.0*step_nr);
  nu= (1.0 + 2.0*lambda - sqrt(1.0 + 4.0*lambda))/(2.0*lambda);
  boundaryscale= 1.0/(1.0-nu);
  postscale = pow(nu/lambda, step_nr);
  if(verbose>1) {
    printf("nu := %g\n", nu);
    printf("boundaryscale := %g\n", boundaryscale);
    printf("postscale := %g\n", postscale);
  }

  /* Set-up memory for an array of pixel values across planes */
  if(verbose>1) printf("Allocating memory\n");
  double *data=(double*)malloc(img->dimz*sizeof(double));
  if(data==NULL) {if(errmsg!=NULL) strcpy(errmsg, "out of memory"); return 5;}

  /* Process */
  if(verbose>1) printf("Gaussian filtering...\n");
  for(int yi=0; yi<img->dimy; yi++) {
    for(int xi=0; xi<img->dimx; xi++) {
      for(int fi=0; fi<img->dimt; fi++) {
        /* Copy image data to double buffer */
        for(int zi=0; zi<img->dimz; zi++) data[zi]=img->m[zi][yi][xi][fi];
        /* 1D Gaussian filtering */
        for(int si=0; si<step_nr; si++) {
          int i=0;
          data[i]*=boundaryscale;
          /* Filter rightwards (causal) */  
          for(i=1; i<img->dimz; i++) data[i]+=nu*data[i-1];
          i=img->dimz-1; data[i]*=boundaryscale;
          /* Filter leftwards (anti-causal) */  
          for(; i>0; i--) data[i-1]+=nu*data[i];  
        }
        for(int zi=0; zi<img->dimz; zi++) data[zi]*=postscale;
        /* Copy processed buffer data back to image */
        for(int zi=0; zi<img->dimz; zi++) img->m[zi][yi][xi][fi]=data[zi];
      } // next frame
    } // next column
  } // next row

  /* Free memory */
  free(data);

  if(errmsg!=NULL) strcpy(errmsg, "ok");
  if(verbose>1) printf("  Gaussian 1D filtering done.\n");
  return 0;
}
/*****************************************************************************/

/*****************************************************************************/
