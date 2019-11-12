/// @file ecat63ml.c
/// @author Vesa Oikonen
/// @brief Reading and writing ECAT 6.3 matrix list.
///
///  Assumptions:
///  1. Assumes that matrix list data is in VAX little endian format
///  2. Data is automatically converted to big endian when read, if necessary
///     according to the current platform.
///
/*****************************************************************************/
#include "libtpcimgio.h"
/*****************************************************************************/

/*****************************************************************************/
/*!
 * Initiate ECAT matrix list. Call this once before first use.
 *
 * @param mlist matrix list
 */
void ecat63InitMatlist(MATRIXLIST *mlist) {
  mlist->matrixSpace=mlist->matrixNr=0; mlist->matdir=NULL;
}
/*****************************************************************************/

/*****************************************************************************/
/*!
 * Free memory allocated for ECAT matrix list
 *
 * @param mlist matrix list
 */
void ecat63EmptyMatlist(MATRIXLIST *mlist) {
  if(mlist->matrixSpace>0) free((char*)(mlist->matdir));
  mlist->matrixSpace=mlist->matrixNr=0;
}
/*****************************************************************************/

/*****************************************************************************/
/** Read ECAT 6.3 matrix list.

    Matrix list must be initiated (once) before calling this.

 @return Returns 0 if ok, 1 if invalid input, 2 first matrix 
  is not found, 3 if failed to read matrix, 4 failed to allocate memory, 
  5 other error.
 */
int ecat63ReadMatlist(
  /** File pointer. */
  FILE *fp,
  /** Pointer to initiated matrix list. */
  MATRIXLIST *ml,
  /** Verbose level; if zero, then nothing is printed to stderr or stdout */
  int verbose
) {
  int i, err=0, little;
  int blk=MatFirstDirBlk, next_blk=0 /*, nr_free, prev_blk, nr_used*/;
  size_t sn;
  unsigned int dirbuf[MatBLKSIZE/4];


  if(verbose>0) {printf("ecat63ReadMatlist(fp, mlist)\n"); fflush(stdout);}
  if(fp==NULL) return(1);
  little=little_endian();
  /* Make sure that matrix list is empty */
  ecat63EmptyMatlist(ml);
  /* Seek the first list block */
  fseek(fp, (blk-1)*MatBLKSIZE, SEEK_SET); if(ftell(fp)!=(blk-1)*MatBLKSIZE) return(2);
  do {
    /* Read the data block */
    if(verbose>1) {printf("  reading dirblock %d\n", blk); fflush(stdout);}
    sn=fread(dirbuf, sizeof(int), MatBLKSIZE/4, fp); 
    if(sn==0) { /* this never happens in valid ECAT file, but it does happen... */
      if(verbose>0) {printf("sn=%d\n", (int)sn); fflush(stdout);}
      //ml->matrixNr--;
      break;
    }
    if(sn<MatBLKSIZE/4) { /* this would be a real error in file */
      if(verbose>0) {printf("sn=%d\n", (int)sn); fflush(stdout);}
      err=2; break;
    }
    /* Allocate (more) memory for one block */
    if(ml->matrixSpace==0) {
      ml->matrixSpace=MatBLKSIZE/4;
      ml->matdir=(MatDir*)malloc(ml->matrixSpace*sizeof(MatDir));
    } else if(ml->matrixSpace<(ml->matrixNr+MatBLKSIZE/4)) {
      ml->matrixSpace+=MatBLKSIZE/4;
      ml->matdir=(MatDir*)realloc(ml->matdir, sizeof(MatDir)*ml->matrixSpace);
    }
    if(ml->matdir==NULL) return(4);
    /* Byte order conversion for ints in big endian platforms */
    if(!little) swawbip(dirbuf, MatBLKSIZE);
    /* Read "header" integers */
    /*nr_free  = dirbuf[0];*/
    next_blk = dirbuf[1];
    /*prev_blk = dirbuf[2];*/
    /*nr_used  = dirbuf[3];*/
    if(verbose>3) printf("next_blk=%d\n", next_blk);
    fpos_t current_fp; fgetpos(fp, &current_fp); // save current file position
    for(i=4; i<MatBLKSIZE/4; i+=4) if(dirbuf[i]>0) {
      ml->matdir[ml->matrixNr].matnum=dirbuf[i];
      ml->matdir[ml->matrixNr].strtblk=dirbuf[i+1];
      ml->matdir[ml->matrixNr].endblk=dirbuf[i+2];
      ml->matdir[ml->matrixNr].matstat=dirbuf[i+3];
      if(verbose>4) printf("matnum=%d strtblk=%d endblk=%d matstat=%d matrixNr=%d\n",
          ml->matdir[ml->matrixNr].matnum, ml->matdir[ml->matrixNr].strtblk,
          ml->matdir[ml->matrixNr].endblk, ml->matdir[ml->matrixNr].matstat,
          ml->matrixNr);
      /* verify that these block can be found in the file */
      int ret=fseek(fp, (ml->matdir[ml->matrixNr].endblk-1)*MatBLKSIZE, SEEK_SET);
      fsetpos(fp, &current_fp);  // back to saved file position
      if(ret==0) { // end block can be found
        ml->matrixNr++;
      } else { // not found, file probably broken
        if(verbose>0) {
          printf("matnum %d points to data outside of file.\n", ml->matdir[ml->matrixNr].matnum);
          fflush(stdout);
        }
      }
    }
    blk=next_blk;
    /* Seek the next list block */
    fseek(fp, (blk-1)*MatBLKSIZE, SEEK_SET); if(ftell(fp)!=(blk-1)*MatBLKSIZE) err=1;
  } while(err==0 && feof(fp)==0 && blk!=MatFirstDirBlk);
  if(err) {ecat63EmptyMatlist(ml); return(5);}
  return(0);
}
/*****************************************************************************/

/*****************************************************************************/
#if(1)
/** Print ECAT matrix list on stdout. */
void ecat63PrintMatlist(
  /** Pointer to matrix list to print. */
  MATRIXLIST *ml
) {
  int i;
  Matval matval;

  printf("nr\tmatrix\tpl\tfr\tgate\tbed\tstartblk\tblknr\n");
  for(i=0; i<ml->matrixNr; i++) {
    mat_numdoc(ml->matdir[i].matnum, &matval);
    printf("%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n", i+1, ml->matdir[i].matnum,
      matval.plane, matval.frame, matval.gate, matval.bed,
      ml->matdir[i].strtblk, 1+ml->matdir[i].endblk-ml->matdir[i].strtblk);
  }
  return;
}
#else
/*!
 * Print ECAT matrix list on stdout.
 *
 * @param ml matrix list
 */
void ecat63PrintMatlist(MATRIXLIST *ml) {
  int i;
  Matval matval;

  printf("nr   matrix   pl  fr gate bed startblk blknr\n");
  for(i=0; i<ml->matrixNr; i++) {
    mat_numdoc(ml->matdir[i].matnum, &matval);
    printf("%4d %8d %3d %3d %3d %3d %8d %3d\n", i+1, ml->matdir[i].matnum,
      matval.plane, matval.frame, matval.gate, matval.bed,
      ml->matdir[i].strtblk, 1+ml->matdir[i].endblk-ml->matdir[i].strtblk);
  }
  return;
}
#endif
/*****************************************************************************/

/*****************************************************************************/
/*!
 * Prepare matrix list for additional matrix data and
 *  Directory records are written in big endian byte order.
 *  Set block_nr to the number of data blocks excluding header;
 *
 * @param fp file pointer
 * @param matnum matrix number [1..number of matrices]
 * @param blkNr matrix block number [ >= 1]
 * @return block number for matrix header, or 0 in case of an error.
 */
int ecat63Matenter(FILE *fp, int matnum, int blkNr) {
  unsigned int i=0, dirblk, little, busy=1, nxtblk=0, oldsize;
  unsigned int dirbuf[MatBLKSIZE/4];

  if(ECAT63_TEST) printf("ecat63Matenter(fp, %d, %d)\n", matnum, blkNr);
  if(fp==NULL || matnum<1 || blkNr<1) return(0);
  little=little_endian(); memset(dirbuf, 0, MatBLKSIZE);
  /* Read first matrix list block */
  dirblk=MatFirstDirBlk;
  fseek(fp, (dirblk-1)*MatBLKSIZE, SEEK_SET);
  if(ftell(fp)!=(int)(dirblk-1)*MatBLKSIZE) return(0);
  if(fread(dirbuf, sizeof(int), MatBLKSIZE/4, fp) != MatBLKSIZE/4) return(0);
  /* Byte order conversion for ints in big endian platforms */
  if(!little) swawbip(dirbuf, MatBLKSIZE);

  while(busy) {
    nxtblk=dirblk+1;
    for(i=4; i<MatBLKSIZE/4; i+=4) {
      if(dirbuf[i]==0) {  /* Check for end of matrix list */
        busy=0; break;
      } else if(dirbuf[i]==(unsigned int)matnum) {  /* Check if this matrix already exists */
        oldsize=dirbuf[i+2]-dirbuf[i+1]+1;
        if(oldsize<(unsigned int)blkNr) {  /* If old matrix is smaller */
          dirbuf[i] = 0xFFFFFFFF;
          if(!little) swawbip(dirbuf, MatBLKSIZE);
          fseek(fp, (dirblk-1)*MatBLKSIZE, SEEK_SET);
          if(ftell(fp)!=(int)(dirblk-1)*MatBLKSIZE) return(0);
          if(fwrite(dirbuf, sizeof(int), MatBLKSIZE/4, fp) != MatBLKSIZE/4) return(0);
          if(!little) swawbip(dirbuf, MatBLKSIZE);
          nxtblk=dirbuf[i+2]+1;
        } else { /* old matrix size is ok */
          nxtblk=dirbuf[i+1]; dirbuf[0]++; dirbuf[3]--; busy=0;
          break;
        }
      } else  /* not this one */
        nxtblk=dirbuf[i+2]+1;
    }
    if(!busy) break;
    if(dirbuf[1]!=MatFirstDirBlk) {
      dirblk=dirbuf[1];
      fseek(fp, (dirblk-1)*MatBLKSIZE, SEEK_SET);
      if(ftell(fp)!=(int)(dirblk-1)*MatBLKSIZE) return(0);
      if(fread(dirbuf, sizeof(int), MatBLKSIZE/4, fp) != MatBLKSIZE/4) return(0);
      if(!little) swawbip(dirbuf, MatBLKSIZE);
    } else {
      dirbuf[1]=nxtblk;
      if(!little) swawbip(dirbuf, MatBLKSIZE);
      fseek(fp, (dirblk-1)*MatBLKSIZE, SEEK_SET);
      if(ftell(fp)!=(int)(dirblk-1)*MatBLKSIZE) return(0);
      if(fwrite(dirbuf, sizeof(int), MatBLKSIZE/4, fp) != MatBLKSIZE/4) return(0);
      dirbuf[0]=31; dirbuf[1]=MatFirstDirBlk; dirbuf[2]=dirblk;
      dirbuf[3]=0; dirblk=nxtblk;
      for(i=4; i<MatBLKSIZE/4; i++) dirbuf[i]=0;
    }
  }
  dirbuf[i]=matnum;
  dirbuf[i+1]=nxtblk;
  dirbuf[i+2]=nxtblk+blkNr;
  dirbuf[i+3]=1;
  dirbuf[0]--;
  dirbuf[3]++;
  if(!little) swawbip(dirbuf, MatBLKSIZE);
  fseek(fp, (dirblk-1)*MatBLKSIZE, SEEK_SET); 
  if(ftell(fp)!=(int)(dirblk-1)*MatBLKSIZE) return(0);
  if(fwrite(dirbuf, sizeof(int), MatBLKSIZE/4, fp) != MatBLKSIZE/4) return(0);
  if(ECAT63_TEST) printf("returning %d from ecat63Matenter()\n", nxtblk);
  return(nxtblk);
}
/*****************************************************************************/

/*****************************************************************************/
/*!
 * Returns the matrix identifier
 *
 * @param frame frame number [0..4096]
 * @param plane plane number [0..256]
 * @param gate gate number [0..64]
 * @param data data number [0..8]
 * @param bed bed position [0..16]
 * @return matrix identifier coding
 */
int mat_numcod(int frame, int plane, int gate, int data, int bed) {
  return((frame&0xFFF)|((bed&0xF)<<12)|((plane&0xFF)<<16)|
         ((gate&0x3F)<<24)|((data&0x3)<<30));
}
/*!
 * Conversion of matrix identifier to numerical values
 *
 * @param matnum matrix identifier coding
 * @param matval target matrix value structure
 */
void mat_numdoc(int matnum, Matval *matval) {
  matval->frame = matnum&0xFFF;
  matval->plane = (matnum>>16)&0xFF;
  matval->gate  = (matnum>>24)&0x3F;
  matval->data  = (matnum>>30)&0x3;
  matval->bed   = (matnum>>12)&0xF;
}
/*****************************************************************************/

/*****************************************************************************/
/*!
 * Sort matrixlist by plane and frame. Bubble sorting algorithm.
 *
 * @param ml marix list.
 */
void ecat63SortMatlistByPlane(MATRIXLIST *ml) {
  int i, j;
  Matval mv1, mv2;
  MatDir tmpMatdir;

  for(i=0; i<ml->matrixNr-1; i++) {
    mat_numdoc(ml->matdir[i].matnum, &mv1);
    for(j=i+1; j<ml->matrixNr; j++) {
      mat_numdoc(ml->matdir[j].matnum, &mv2);
      if(mv2.plane<mv1.plane||(mv2.plane==mv1.plane&&mv2.frame<mv1.frame)) {
        tmpMatdir=ml->matdir[i];
        ml->matdir[i]=ml->matdir[j]; ml->matdir[j]=tmpMatdir;
        mat_numdoc(ml->matdir[i].matnum, &mv1);
      }
    }
  }
}
/*****************************************************************************/

/*****************************************************************************/
/*!
 * Sort matrixlist by frame and plane. Bubble sorting algorithm.
 *
 * @param ml matrix list
 */
void ecat63SortMatlistByFrame(MATRIXLIST *ml) {
  int i, j;
  Matval mv1, mv2;
  MatDir tmpMatdir;

  for(i=0; i<ml->matrixNr-1; i++) {
    mat_numdoc(ml->matdir[i].matnum, &mv1);
    for(j=i+1; j<ml->matrixNr; j++) {
      mat_numdoc(ml->matdir[j].matnum, &mv2);
      if(mv2.frame<mv1.frame||(mv2.frame==mv1.frame&&mv2.plane<mv1.plane)) {
        tmpMatdir=ml->matdir[i];
        ml->matdir[i]=ml->matdir[j]; ml->matdir[j]=tmpMatdir;
        mat_numdoc(ml->matdir[i].matnum, &mv1);
      }
    }
  }
}
/*****************************************************************************/

/*****************************************************************************/
/*!
 * Checks that all matrixlist entries have read/write status.
 *
 * @param ml matrix list
 * @return 0 if ok, or 1 if an entry is marked as deleted or unfinished.
 */
int ecat63CheckMatlist(MATRIXLIST *ml) {
  int i;

  if(ml==NULL) return(1);
  for(i=0; i<ml->matrixNr; i++) if(ml->matdir[i].matstat!=1) return(1);
  return(0);
}
/*****************************************************************************/

/*****************************************************************************/
/** Mark deleted the frames after the specified frame number.

    This can be used to delete sum images from the end of dynamic ECAT images.

 @param ml matrix list.
 @param frame_nr last index not to be marked as deleted.
 @return number of deleted matrices.
 */
int ecat63DeleteLateFrames(MATRIXLIST *ml, int frame_nr) {
  int i, del_nr=0;
  Matval matval;

  for(i=0; i<ml->matrixNr; i++) {
    mat_numdoc(ml->matdir[i].matnum, &matval);
    if(matval.frame>frame_nr) {del_nr++; ml->matdir[i].matstat=-1;}
  }
  return(del_nr);
}
/*****************************************************************************/

/*****************************************************************************/
/*!
 * Calculate the size of one data matrix in ECAT 6.3 file matrix list, and
 * check that the size is same in all matrices.
 *
 * @param mlist Ecat 6.3 matrix list; note that this list is here sorted by planes
 * @param blk_nr Number of blocks will be put here; NULL if not needed
 * @return errstatus, which is STATUS_OK (0) when call was successful,
 * and >0 in case of an error.
 */
int ecat63GetMatrixBlockSize(MATRIXLIST *mlist, int *blk_nr) {
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
 * Calculate the number of planes and frames/gates from ECAT 6.3 matrix list.
 * Check that all planes have equal nr of frames/gates, that frames/gates
 * are sequentially numbered. This routines sorts the matrix list by planes.
 *
 * @param mlist Ecat 6.3 matrix list; note that this list is here sorted by planes
 * @param h Ecat 6.3 mainheader
 * @param plane_nr Number of planes will be put here; NULL if not needed
 * @param frame_nr Number of frames/gates will be put here; NULL if not needed
 * @return errstatus, which is STATUS_OK (0) when call was successful,
 * and >0 in case of an error.
 */
int ecat63GetPlaneAndFrameNr(MATRIXLIST *mlist, ECAT63_mainheader *h, int *plane_nr, int *frame_nr) {
  Matval matval;
  int m, plane, frame, prev_plane, prev_frame, fnr, pnr;

  /* Check input */
  if(mlist==NULL) return STATUS_FAULT;
  if(plane_nr!=NULL) *plane_nr=0;
  if(frame_nr!=NULL) *frame_nr=0;

  /* Sort matrices by plane so that following computation works */
  ecat63SortMatlistByPlane(mlist);

  prev_plane=plane=-1; prev_frame=frame=-1;
  fnr=pnr=0;
  for(m=0; m<mlist->matrixNr; m++) if(mlist->matdir[m].matstat==1) {
    mat_numdoc(mlist->matdir[m].matnum, &matval);
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
 * Read the maximum plane, frame, gate and bed number from matrixlist.
 *
 * @param ml Pointer to matrixlist
 * @param num_planes number of planes will be put here; NULL if not needed
 * @param num_frames number of frames will be put here; NULL if not needed
 * @param num_gates number of gates will be put here; NULL if not needed
 * @param num_bed_pos number of gates will be put here; NULL if not needed
 * @return 0 if successful, 1 no matrix list, 2 invalid matrix number, 
 * 3 failed to allocate memory
 */
int ecat63GetNums(MATRIXLIST *ml, short int *num_planes, short int *num_frames, short int *num_gates, short int *num_bed_pos) {
  int i, nmax;
  Matval* matval;

  if(ml==NULL) return(1);
  if(ml->matrixNr<1) return(2);

  /* Allocate memory for matrix values */
  matval = (Matval*)calloc(ml->matrixNr,sizeof(Matval));
  if(matval == NULL) return(3);

  /* And get the matrix values */
  for(i=0; i<ml->matrixNr; i++) mat_numdoc(ml->matdir[i].matnum, matval+i);

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
  free(matval);
  return(0);
}
/*****************************************************************************/

/*****************************************************************************/
/*!
 * Matrix numbers in ECAT 6.3 matrix list are edited, when necessary, so that
 *  plane, frame, gate and/or bed numbers are continuous, starting from one
 *  (planes, frames and gates) or from zero (beds).
 *  List order is not changed.
 *
 * @param ml ECAT 6.3 matrix list, where the matrix numbers will be edited
 * @param do_planes Plane numbers are gathered together (1) or not (0)
 * @param do_frames Frame numbers are gathered together (1) or not (0)
 * @param do_gates Gate numbers are gathered together (1) or not (0)
 * @param do_beds Bed numbers are gathered together (1) or not (0)
 * @return 0 if successful, 1 if invalid input, 3 failed to allocate memory
 */
int ecat63GatherMatlist(MATRIXLIST *ml, short int do_planes, short int do_frames, short int do_gates, short int do_beds) {
  int i, ncurr, n;
  Matval* matval;

  if(ml==NULL) return(1);
  if(ml->matrixNr<1) return(0);

  /* Allocate memory for matrix values */
  matval = (Matval*)calloc(ml->matrixNr,sizeof(Matval));
  if(matval == NULL) return(3);

  /* And get the matrix values */
  for(i=0; i<ml->matrixNr; i++) mat_numdoc(ml->matdir[i].matnum, matval+i);

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
  for(i=0; i<ml->matrixNr; i++) ml->matdir[i].matnum=mat_numcod(
            matval[i].frame, matval[i].plane,
            matval[i].gate, matval[i].data,
            matval[i].bed);
  free(matval);
  return(0);
}
/*****************************************************************************/

/*****************************************************************************/

