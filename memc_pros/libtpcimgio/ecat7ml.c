/// @file ecat7ml.c
/// @author Vesa Oikonen, Harri Merisaari
/// @brief Reading and writing ECAT 7.x matrix list.
///
/*****************************************************************************/
#include "libtpcimgio.h"
/*****************************************************************************/

/*****************************************************************************/
/*!
 * Initiate ECAT matrix list. Call this once before first use.
 *
 * @param mlist target matrix list
 */
void ecat7InitMatlist(ECAT7_MATRIXLIST *mlist) {
  mlist->matrixSpace=mlist->matrixNr=0; mlist->matdir=NULL;
}
/*****************************************************************************/

/*****************************************************************************/
/*!
 * Free memory allocated for ECAT matrix list.
 *
 * @param mlist target matrix list that has allocated memory
 */
void ecat7EmptyMatlist(ECAT7_MATRIXLIST *mlist) {
  if(mlist->matrixSpace>0) free((char*)(mlist->matdir));
  mlist->matrixSpace=mlist->matrixNr=0;
}
/*****************************************************************************/

/*****************************************************************************/
/** Read ECAT matrix list.

    Matrix list must be initiated (once) before calling this.
 
 @return returns 0 if ok, 1 if invalid input, 2 if first matrix is not found,
  3 if failed to read matrix, 4 if data allocation failed for matrix, 5 if other error
  occurred.
 */
int ecat7ReadMatlist(
  /** File pointer. */
  FILE *fp,
  /** Pointer to initiated matrix list. */
  ECAT7_MATRIXLIST *ml,
  /** Verbose level; if zero, then nothing is printed to stderr or stdout */
  int verbose
) {
  int i, err=0, little;
  int blk=MatFirstDirBlk, next_blk=0, nr_free, prev_blk, nr_used;
  size_t sn;
  unsigned int dirbuf[MatBLKSIZE/4];


  if(verbose>0) printf("ecat7ReadMatlist(fp, mlist)\n");
  if(fp==NULL) return(1);
  little=little_endian();
  /* Make sure that matrix list is empty */
  ecat7EmptyMatlist(ml);
  /* Seek the first list block */
  fseek(fp, (blk-1)*MatBLKSIZE, SEEK_SET); if(ftell(fp)!=(blk-1)*MatBLKSIZE) return(2);
  do {
    /* Read the data block */
    if(verbose>1) printf("  reading dirblock %d\n", blk);
    sn=fread(dirbuf, sizeof(int), MatBLKSIZE/4, fp); if(sn<MatBLKSIZE/4) return(3);
    /* Allocate (more) memory for one block */
    if(ml->matrixSpace==0) {
      ml->matrixSpace=MatBLKSIZE/4;
      ml->matdir=(ECAT7_MatDir*)malloc(ml->matrixSpace*sizeof(ECAT7_MatDir));
    } else if(ml->matrixSpace<(ml->matrixNr+MatBLKSIZE/4)) {
      ml->matrixSpace+=MatBLKSIZE/4;
      ml->matdir=(ECAT7_MatDir*)realloc(ml->matdir, sizeof(ECAT7_MatDir)*ml->matrixSpace);
    }
    if(ml->matdir==NULL) return(4);
    /* Byte order conversion for ints in little endian platforms */
    if(little) swawbip(dirbuf, MatBLKSIZE);
    /* Read "header" integers */
    nr_free  = dirbuf[0];
    next_blk = dirbuf[1];
    prev_blk = dirbuf[2];
    nr_used  = dirbuf[3];
    if(verbose>2) printf("nr_free=%d next_blk=%d prev_blk=%d nr_used=%d\n", nr_free, next_blk, prev_blk, nr_used);
    for(i=4; i<MatBLKSIZE/4; i+=4) if(dirbuf[i]>0) {
      ml->matdir[ml->matrixNr].id=dirbuf[i];
      ml->matdir[ml->matrixNr].strtblk=dirbuf[i+1];
      ml->matdir[ml->matrixNr].endblk=dirbuf[i+2];
      ml->matdir[ml->matrixNr].status=dirbuf[i+3];
      if(verbose>3) {
        printf("matnum=%d strtblk=%d endblk=%d matstat=%d matrixNr=%d\n",
          ml->matdir[ml->matrixNr].id, ml->matdir[ml->matrixNr].strtblk,
          ml->matdir[ml->matrixNr].endblk, ml->matdir[ml->matrixNr].status,
          ml->matrixNr);
      }
      ml->matrixNr++;
    }
    blk=next_blk;
    /* Seek the next list block */
    fseek(fp, (blk-1)*MatBLKSIZE, SEEK_SET); if(ftell(fp)!=(blk-1)*MatBLKSIZE) err=1;
  } while(err==0 && feof(fp)==0 && blk!=MatFirstDirBlk);
  if(err) {ecat7EmptyMatlist(ml); return(5);}
  return(0);
}
/*****************************************************************************/

/*****************************************************************************/
/*!
 * Print ECAT matrix list on stdout.
 *
 * @param ml matrix list for Ecat7 file
 */
void ecat7PrintMatlist(ECAT7_MATRIXLIST *ml) {
  int i;
  ECAT7_Matval matval;

  printf("nr   matrix   pl  fr gate bed startblk blknr  status\n");
  for(i=0; i<ml->matrixNr; i++) {
    ecat7_id_to_val(ml->matdir[i].id, &matval);
    printf("%4d %8d %3d %3d %3d %3d %8d %5d  ", i+1, ml->matdir[i].id,
      matval.plane, matval.frame, matval.gate, matval.bed,
      ml->matdir[i].strtblk, 1+ml->matdir[i].endblk-ml->matdir[i].strtblk);
    if(ml->matdir[i].status==1) printf("read/write\n");
    else if(ml->matdir[i].status==0) printf("not ready\n");
    else if(ml->matdir[i].status==-1) printf("deleted\n");
    else printf("%d\n", ml->matdir[i].status);
  }
  return;
}
/*****************************************************************************/

/*****************************************************************************/
/*!
 * Prepare matrix list for additional matrix data and return block number
 * for matrix header. Directory records are written in big endian byte order.
 * Set block_nr to the number of data blocks + (nr of header blocks - 1)
 *
 * @param fp file pointer
 * @param matrix_id matrix identifier coding
 * @param block_nr matrix number [1..number of matrixes]
 * @return returns the block number for matrix header, -1 if invalid input,
 * -2 if first directory block is not found, -3 if failed to read first block,
 * -9 if other directory block is not found, -10 if failed to read other block,
 * -11 if place for new directory block is not found, -12 if failed clear new
 * block, -15 if place for new directory block is not found, -16 if failed to 
 * write into new block
 */
int ecat7EnterMatrix(FILE *fp, int matrix_id, int block_nr) {
  unsigned int i=0, dirblk, little, busy=1, nxtblk=0, oldsize;
  /*unsigned*/ int dirbuf[MatBLKSIZE/4];

  if(ECAT7_TEST) printf("ecat7EnterMatrix(fp, %d, %d)\n", matrix_id, block_nr);
  /* Check the input */
  if(fp==NULL || matrix_id<1 || block_nr<1) return(-1);
  /* Is this a little endian machine? */
  little=little_endian();
  /* Read first directory record block */
  dirblk=MatFirstDirBlk;
  fseek(fp, (dirblk-1)*MatBLKSIZE, SEEK_SET);
  if(ftell(fp)!=(int)(dirblk-1)*MatBLKSIZE) return(-2);
  if(fread(dirbuf, sizeof(int), MatBLKSIZE/4, fp) != MatBLKSIZE/4) return(-3);
  /* Byte order conversion for ints in little endian platforms */
  if(little) swawbip(dirbuf, MatBLKSIZE);
  /* Read through the existing directory records */
  while(busy) {
    /* Go through the directory entries in this record */
    for(i=4, nxtblk=dirblk+1; i<MatBLKSIZE/4; i+=4) {
      oldsize=dirbuf[i+2]-dirbuf[i+1]+1;
      if(dirbuf[i]==0) {  /* Check for end of matrix list */
        busy=0; break;
      } else if(dirbuf[i]==matrix_id) {  /* Maybe this matrix already exists? */
        /* yes it does; is old data smaller? */
        if(oldsize<(unsigned int)block_nr) {
          /* it was smaller, so do not use it, but mark it deleted */
          dirbuf[i] = 0xFFFFFFFF; dirbuf[i+3]=-1;
          if(little) swawbip(dirbuf, MatBLKSIZE);
          fseek(fp, (dirblk-1)*MatBLKSIZE, SEEK_SET);
          if(ftell(fp)!=(int)(dirblk-1)*MatBLKSIZE) return(-6);
          if(fwrite(dirbuf, sizeof(int), MatBLKSIZE/4, fp) != MatBLKSIZE/4) return(-7);
          if(little) swawbip(dirbuf, MatBLKSIZE);
          nxtblk=dirbuf[i+2]+1;
        } else { /* old matrix size is ok */
          nxtblk=dirbuf[i+1]; dirbuf[0]++; dirbuf[3]--; busy=0;
          break;
        }
      } else { /* this is not the same matrix */
        /* But is deleted and of same or smaller size? */
        if(dirbuf[i+3]==-1 && (unsigned int)block_nr<=oldsize) {
          /* yes it was, so lets recycle it */
          dirbuf[i]=matrix_id;
          nxtblk=dirbuf[i+1]; dirbuf[0]++; dirbuf[3]--; busy=0;
          break;
        }
        /* nothing to be done with this entry */
        nxtblk=dirbuf[i+2]+1;
      }
    } /* next entry in this record */
    if(!busy) break; /* stop reading existing records */
    /* Read the next directory record */
    if(dirbuf[1]!=MatFirstDirBlk) {
      /* There are more records left to read */
      dirblk=dirbuf[1];
      fseek(fp, (dirblk-1)*MatBLKSIZE, SEEK_SET);
      if(ftell(fp)!=(int)(dirblk-1)*MatBLKSIZE) return(-9);
      if(fread(dirbuf, sizeof(int), MatBLKSIZE/4, fp) != MatBLKSIZE/4) return(-10);
      if(little) swawbip(dirbuf, MatBLKSIZE);
    } else {
      /* No more records to read, so lets write a new empty one */
      dirbuf[1]=nxtblk; /* write a pointer to the new one */
      if(little) swawbip(dirbuf, MatBLKSIZE);
      fseek(fp, (dirblk-1)*MatBLKSIZE, SEEK_SET);
      if(ftell(fp)!=(int)(dirblk-1)*MatBLKSIZE) return(-11);
      if(fwrite(dirbuf, sizeof(int), MatBLKSIZE/4, fp) != MatBLKSIZE/4) return(-12);
      /* and then initiate the contents of the next one, but do not write it */
      dirbuf[0]=31; dirbuf[1]=MatFirstDirBlk; dirbuf[2]=dirblk;
      dirbuf[3]=0; dirblk=nxtblk;
      for(i=4; i<MatBLKSIZE/4; i++) dirbuf[i]=0;
    }
  } /* next directory record */
  dirbuf[i]=matrix_id;
  dirbuf[i+1]=nxtblk;
  dirbuf[i+2]=nxtblk+block_nr;
  dirbuf[i+3]=1; /* mark the entry as read/write */
  dirbuf[0]--;
  dirbuf[3]++;
  if(little) swawbip(dirbuf, MatBLKSIZE);
  fseek(fp, (dirblk-1)*MatBLKSIZE, SEEK_SET);
  if(ftell(fp)!=(int)(dirblk-1)*MatBLKSIZE) return(-15);
  if(fwrite(dirbuf, sizeof(int), MatBLKSIZE/4, fp) != MatBLKSIZE/4) return(-16);
  if(ECAT7_TEST) printf("returning %d from ecat7EnterMatrix()\n", nxtblk);
  return(nxtblk);
}
/*****************************************************************************/

/*****************************************************************************/
/*!
 * Returns the matrix identifier.
 *
 * @param frame frame number [0..65536]
 * @param plane plane number [0..65536]
 * @param gate gate number [0..64]
 * @param data data [0..1]
 * @param bed bed position [0..16]
 * @return matrix identifier coding
 */
int ecat7_val_to_id(int frame, int plane, int gate, int data, int bed) {
  return(
    ((bed & 0xF) << 12) |    /* bed */
    (frame & 0x1FF) |        /* frame */
    ((gate & 0x3F) << 24) |  /* gate */
    ((plane & 0xFF) << 16) | /* plane low */
    ((plane & 0x300) << 1) | /* plane high */
    ((data & 0x3) << 30) |   /* data low */
    ((data & 0x4) << 9)      /* data high */
  );
}
/*!
 * Conversion of matrix identifier to numerical values
 *
 * @param matrix_id matrix identifier coding
 * @param matval matrix values structure
 */
void ecat7_id_to_val(int matrix_id, ECAT7_Matval *matval) {
  matval->frame = matrix_id & 0x1FF;
  matval->plane = ((matrix_id >> 16) & 0xFF) + ((matrix_id >> 1) & 0x300);
  matval->gate  = (matrix_id >> 24) & 0x3F;
  matval->data  = ((matrix_id >> 30) & 0x3) + ((matrix_id >> 9) & 0x4);
  matval->bed   = (matrix_id >> 12) & 0xF;
}
/*****************************************************************************/

/*****************************************************************************/
/*!
 * Sort matrixlist by plane and frame. Bubble sorting algorithm.
 *
 * @param ml target matrix list
 */
void ecat7SortMatlistByPlane(ECAT7_MATRIXLIST *ml) {
  int i, j;
  ECAT7_Matval mv1, mv2;
  ECAT7_MatDir tmpMatdir;

  for(i=0; i<ml->matrixNr-1; i++) {
    ecat7_id_to_val(ml->matdir[i].id, &mv1);
    for(j=i+1; j<ml->matrixNr; j++) {
      ecat7_id_to_val(ml->matdir[j].id, &mv2);
      if(mv2.plane<mv1.plane||(mv2.plane==mv1.plane&&mv2.frame<mv1.frame)) {
        tmpMatdir=ml->matdir[i];
        ml->matdir[i]=ml->matdir[j];
	ml->matdir[j]=tmpMatdir;
        ecat7_id_to_val(ml->matdir[i].id, &mv1);
      }
    }
  }
}
/*****************************************************************************/

/*****************************************************************************/
/*!
 * Sort matrixlist by frame and plane. Bubble sorting algorithm.
 *
 * @param ml target matrix list
 */
void ecat7SortMatlistByFrame(ECAT7_MATRIXLIST *ml) {
  int i, j;
  ECAT7_Matval mv1, mv2;
  ECAT7_MatDir tmpMatdir;

  for(i=0; i<ml->matrixNr-1; i++) {
    ecat7_id_to_val(ml->matdir[i].id, &mv1);
    for(j=i+1; j<ml->matrixNr; j++) {
      ecat7_id_to_val(ml->matdir[j].id, &mv2);
      if(mv2.frame<mv1.frame||(mv2.frame==mv1.frame&&mv2.plane<mv1.plane)) {
        tmpMatdir=ml->matdir[i];
        ml->matdir[i]=ml->matdir[j]; ml->matdir[j]=tmpMatdir;
        ecat7_id_to_val(ml->matdir[i].id, &mv1);
      }
    }
  }
}
/*****************************************************************************/

/*****************************************************************************/
/*!
 * Checks that all matrixlist entries have read/write status.
 *
 * @param ml checked matrix list
 * @return 0 if ok, or 1 if an entry is marked as deleted or unfinished
 */
int ecat7CheckMatlist(ECAT7_MATRIXLIST *ml) {
  int i;

  if(ml==NULL) return(1);
  for(i=0; i<ml->matrixNr; i++) if(ml->matdir[i].status!=1) return(1);
  return(0);
}
/*****************************************************************************/

/*****************************************************************************/
/*!
 * Mark deleted the frames after the specified frame number.
 *
 * @param ml target matrix list
 * @param frame_nr first index to be marked as deleted [1..number of frames]
 * @return Returns the number of deleted matrices.
 */
int ecat7DeleteLateFrames(ECAT7_MATRIXLIST *ml, int frame_nr) {
  int i, del_nr=0;
  ECAT7_Matval matval;

  for(i=0; i<ml->matrixNr; i++) {
    ecat7_id_to_val(ml->matdir[i].id, &matval);
    if(matval.frame>frame_nr) {del_nr++; ml->matdir[i].status=-1;}
  }
  return(del_nr);
}
/*****************************************************************************/

/*****************************************************************************/
/*!
 * Calculate the number of planes and frames/gates from ECAT7 matrix list.
 * Check that all planes have equal nr of frames/gates, that frames/gates
 * are sequentally numbered.
 *
 * @param mlist Ecat7 matrix list; note that this list is here sorted by planes
 * @param h Ecat7 main header structure
 * @param plane_nr Number of planes will be put here; NULL if not needed [1..number of planes, or NULL]
 * @param frame_nr Number of frames/gates will be put here; NULL if not needed [1..number of frames, or NULL]
 * @return errstatus, which is STATUS_OK (0) when call was successful, and >0 in case of an error.
 * Note that if this is 3D image volume or sinogram, then the returned plane_nr
 * will be one, and the actual Z dim must be read from subheader.
 */
int ecat7GetPlaneAndFrameNr(ECAT7_MATRIXLIST *mlist, ECAT7_mainheader *h, int *plane_nr, int *frame_nr) {
  ECAT7_Matval matval;
  int m, plane, frame, prev_plane, prev_frame, fnr, pnr;

  /* Check input */
  if(mlist==NULL) return STATUS_FAULT;
  if(plane_nr!=NULL) *plane_nr=0;
  if(frame_nr!=NULL) *frame_nr=0;

  /* Sort matrices by plane so that following computation works */
  ecat7SortMatlistByPlane(mlist);

  prev_plane=plane=-1; prev_frame=frame=-1;
  fnr=pnr=0;
  for(m=0; m<mlist->matrixNr; m++) {
    ecat7_id_to_val(mlist->matdir[m].id, &matval);
    plane=matval.plane;
    if(h->num_frames>=h->num_gates)
      frame=matval.frame;
    else
      frame=matval.gate;
    if(plane!=prev_plane) {
      fnr=1; pnr++;
    } else {
      fnr++;
      if(prev_frame>0 && frame!=prev_frame+1) return STATUS_MISSINGMATRIX;
    }
    prev_plane=plane; prev_frame=frame;
  } /* next matrix */
  if(fnr*pnr != mlist->matrixNr) return STATUS_MISSINGMATRIX;
  if(plane_nr!=NULL) *plane_nr=pnr;
  if(frame_nr!=NULL) *frame_nr=fnr;
  return STATUS_OK;
}
/*****************************************************************************/

/*****************************************************************************/
/*!
 * Calculate the size of one data matrix in ECAT7 file matrix list, and
 * check that the size is same in all matrices.
 *
 * @param mlist Ecat7 matrix list; note that this list is here sorted by planes
 * @param blk_nr number of blocks will be put here; NULL if not needed
 * @return errstatus, which is STATUS_OK (0) when call was successful,
 * and >0 in case of an error.
 */
int ecat7GetMatrixBlockSize(ECAT7_MATRIXLIST *mlist, int *blk_nr) {
  int m, prev_blk, blk;

  /* Check input */
  if(mlist==NULL) return STATUS_FAULT;
  if(blk_nr!=NULL) *blk_nr=0;

  /* Calculate the size of first data matrix */
  m=0; prev_blk=blk=mlist->matdir[m].endblk - mlist->matdir[m].strtblk;
  for(m=1; m<mlist->matrixNr; m++) {
    blk=mlist->matdir[m].endblk - mlist->matdir[m].strtblk;
    if(blk!=prev_blk) return STATUS_VARMATSIZE;
    else prev_blk=blk;
  }
  if(blk_nr!=NULL) *blk_nr=blk;
  return STATUS_OK;
}
/*****************************************************************************/

/*****************************************************************************/
/*!
 * Read the maximum plane, frame, gate and bed number from matrixlist.
 * In case of 3D formats, num_planes is checked from the first subheader.
 *
 * @param ml Pointer to matrixlist
 * @param mh Pointer to mainheader
 * @param fp File pointer to ECAT7 file opened in binary mode
 * @param num_planes num_planes will be put here; NULL if not needed to be read
 * @param num_frames num_planes will be put here; NULL if not needed to be read
 * @param num_gates num_planes will be put here; NULL if not needed to be read
 * @param num_bed_pos num_planes will be put here; NULL if not needed to be read
 * @return 0 if successful,  1 if invalid input, 2 if no matrixes, 3 failed to allocate memory, 
 * 5 if failed to read image/scan header information
 */
int ecat7GetNums(ECAT7_MATRIXLIST *ml, ECAT7_mainheader *mh, FILE *fp, short int *num_planes,
		 short int *num_frames, short int *num_gates, short int *num_bed_pos) {
  int i, nmax, ret=0;
  ECAT7_imageheader ih;
  ECAT7_scanheader sh;
  ECAT7_Matval* matval;

  if(ml==NULL) return(1);
  if(ml->matrixNr<1) return(2);

  /* Allocate memory for matrix values */
  matval = (ECAT7_Matval*)calloc(ml->matrixNr,sizeof(ECAT7_Matval));
  if(matval == NULL) return(3);

  /* And get the matrix values */
  for(i=0; i<ml->matrixNr; i++) ecat7_id_to_val(ml->matdir[i].id, matval+i);

  /* Planes */
  if(num_planes!=NULL) {
    nmax=matval[0].plane;
    for(i=1; i<ml->matrixNr; i++) if(matval[i].plane>nmax) nmax=matval[i].plane;
    *num_planes=nmax;
  }
  /* Frames */
  if(num_frames!=NULL) {
    nmax=matval[0].frame;
    for(i=1; i<ml->matrixNr; i++) if(matval[i].frame>nmax) nmax=matval[i].frame;
    *num_frames=nmax;
  }
  /* Gates */
  if(num_gates!=NULL) {
    nmax=matval[0].gate;
    for(i=1; i<ml->matrixNr; i++) if(matval[i].gate>nmax) nmax=matval[i].gate;
    *num_gates=nmax;
  }
  /* Beds */
  if(num_bed_pos!=NULL) {
    nmax=matval[0].bed;
    for(i=1; i<ml->matrixNr; i++) if(matval[i].bed>nmax) nmax=matval[i].bed;
    *num_bed_pos=nmax;
  }

  /* Check the num_planes from the first subheader in 3D formats */
  if(num_planes!=NULL && *num_planes<=1) switch(mh->file_type) {
    case ECAT7_VOLUME8:
    case ECAT7_VOLUME16:
      ret=ecat7ReadImageheader(fp, ml->matdir[0].strtblk, &ih);
      if(ret!=0) {
      	free(matval);
      	return(5);
      }
      if(ih.num_dimensions>2 && ih.z_dimension>1) *num_planes=ih.z_dimension;
      break;
    case ECAT7_3DSCAN:
    case ECAT7_3DSCAN8:
    case ECAT7_3DSCANFIT:
      ret=ecat7ReadScanheader(fp, ml->matdir[0].strtblk, &sh);
      if(ret!=0) {
      	free(matval);
      	return(5);
      }
      for(i=0, *num_planes=0; i<64; i++) *num_planes+=sh.num_z_elements[i];
      break;
  }
  free(matval);
  return(0);
}
/*****************************************************************************/

/*****************************************************************************/
/*!
 * Matrix numbers in ECAT 7 matrix list are edited, when necessary, so that
 *  plane, frame, gate and/or bed numbers are continuous, starting from one
 *  (planes, frames and gates) or from zero (beds).
 *  List order is not changed.
 *
 * @param ml ECAT 7 matrix list, where the matrix numbers will be edited
 * @param do_planes Plane numbers are gathered together (1) or not (0)
 * @param do_frames Frame numbers are gathered together (1) or not (0)
 * @param do_gates Gate numbers are gathered together (1) or not (0)
 * @param do_beds Bed numbers are gathered together (1) or not (0)
 * @return 0 if successful, 1 if invalid input, 3 failed to allocate memory
 */
int ecat7GatherMatlist(ECAT7_MATRIXLIST *ml, short int do_planes, short int do_frames,
		       short int do_gates, short int do_beds) {
  int i, ncurr, n;
  ECAT7_Matval* matval;

  if(ml==NULL) return(1);
  if(ml->matrixNr<1) return(0);

  /* Allocate memory for matrix values */
  matval = (ECAT7_Matval*)calloc(ml->matrixNr,sizeof(ECAT7_Matval));
  if(matval == NULL) return(3);

  /* And get the matrix values */
  for(i=0; i<ml->matrixNr; i++) ecat7_id_to_val(ml->matdir[i].id, matval+i);

  /* Planes */
  if(do_planes!=0) {
    ncurr=1;
    while(ncurr <= ml->matrixNr) {
      /* Find any matrix with this number? */
      for(i=0, n=0; i<ml->matrixNr; i++) if(matval[i].plane==ncurr) {n=1; break;}
      /* If yes, then go on to the next matrix number */
      if(n==1) {ncurr++; continue;}
      /* If not, then subtract 1 from all matrix numbers that are larger */
      for(i=0, n=0; i<ml->matrixNr; i++)
        if(matval[i].plane>ncurr) {
          /*printf("    plane %d -> plane %d\n", matval[i].plane, matval[i].plane-1);*/
          matval[i].plane--; n++;
        }
      /* If no larger values were found any more, then quit */
      if(n<1) break;
    }
  }

  /* Frames */
  if(do_frames!=0) {
    ncurr=1;
    while(ncurr <= ml->matrixNr) {
      /* Find any matrix with this number? */
      for(i=0, n=0; i<ml->matrixNr; i++) if(matval[i].frame==ncurr) {n=1; break;}
      /* If yes, then go on to the next matrix number */
      if(n==1) {ncurr++; continue;}
      /* If not, then subtract 1 from all matrix numbers that are larger */
      for(i=0, n=0; i<ml->matrixNr; i++)
        if(matval[i].frame>ncurr) {matval[i].frame--; n++;}
      /* If no larger values were found any more, then quit */
      if(n<1) break;
    }
  }

  /* Gates */
  if(do_gates!=0) {
    ncurr=1;
    while(ncurr <= ml->matrixNr) {
      /* Find any matrix with this number? */
      for(i=0, n=0; i<ml->matrixNr; i++) if(matval[i].gate==ncurr) {n=1; break;}
      /* If yes, then go on to the next matrix number */
      if(n==1) {ncurr++; continue;}
      /* If not, then subtract 1 from all matrix numbers that are larger */
      for(i=0, n=0; i<ml->matrixNr; i++)
        if(matval[i].gate>ncurr) {matval[i].gate--; n++;}
      /* If no larger values were found any more, then quit */
      if(n<1) break;
    }
  }

  /* Beds */
  if(do_beds!=0) {
    ncurr=1;
    while(ncurr <= ml->matrixNr) {
      /* Find any matrix with this number? */
      for(i=0, n=0; i<ml->matrixNr; i++) if(matval[i].bed==ncurr) {n=1; break;}
      /* If yes, then go on to the next matrix number */
      if(n==1) {ncurr++; continue;}
      /* If not, then subtract 1 from all matrix numbers that are larger */
      for(i=0, n=0; i<ml->matrixNr; i++)
        if(matval[i].bed>ncurr) {matval[i].bed--; n++;}
      /* If no larger values were found any more, then quit */
      if(n<1) break;
    }
  }

  /* Write matrix values (possibly changed) into matrix list */
  for(i=0; i<ml->matrixNr; i++) ml->matdir[i].id=ecat7_val_to_id(
            matval[i].frame, matval[i].plane,
	    matval[i].gate, matval[i].data,
	    matval[i].bed);
  free(matval);
  return(0);
}
/*****************************************************************************/

/*****************************************************************************/

