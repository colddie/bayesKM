/// @file dft.c
/// @author Vesa Oikonen
/// @brief Functions for processing TAC data in DFT structs.
///
/*****************************************************************************/
#include "libtpccurveio.h"
#include <unistd.h>
/*****************************************************************************/
/* local function definitions */
/// @cond
int dftQSortName(const void *voi1, const void *voi2);
int dftQSortPlane(const void *voi1, const void *voi2);
/// @endcond
/*****************************************************************************/

/*****************************************************************************/
/** Free memory allocated for DFT. All data is cleared. */
void dftEmpty(
  /** Pointer to initiated DFT struct data */
  DFT *data
) {
  if(data==NULL) return;
  if(data->_voidataNr>0) free((char*)(data->voi));
  if(data->_dataSize>0) free((char*)(data->_data));
  data->_dataSize=data->_voidataNr=0;
  data->frameNr=data->voiNr=0;
  data->studynr[0]=data->comments[0]=data->unit[0]=(char)0;
  data->radiopharmaceutical[0]=data->isotope[0]=data->decayCorrected=(char)0;
  data->scanStartTime[0]=data->injectionTime[0]=(char)0;
  data->timeunit=data->timetype=0;
}
/*****************************************************************************/

/*****************************************************************************/
/** Initiate DFT structure. This should be called once before use. */
void dftInit(
  /** Pointer to initiated DFT struct data */
  DFT *data
) {
  if(data==NULL) return;
  memset(data, 0, sizeof(DFT));
  data->_voidataNr=data->_dataSize=0;
  data->frameNr=data->voiNr=0;
  data->studynr[0]=data->comments[0]=data->unit[0]=(char)0;
  data->radiopharmaceutical[0]=data->isotope[0]=data->decayCorrected=(char)0;
  data->scanStartTime[0]=data->injectionTime[0]=(char)0;
  data->timeunit=data->timetype=data->isweight=0;
}
/*****************************************************************************/

/*****************************************************************************/
/** Allocate memory for DFT data and sets data pointers.
\return Returns <> 0 in case of an error.
 */
int dftSetmem(
  /** Pointer to initiated DFT struct data; any old contens are deleted. */
  DFT *data,
  /** Nr of time frames (samples) to allocate */
  int frameNr,
  /** Nr of concentration arrays (regional TACs) to allocate */
  int voiNr
) {
  int i, n;
  double *d;


  if(data==NULL) return 1;
  /* Clear previous data */
  dftEmpty(data);

  /* Allocate memory for curves */
  data->voi=(Voi*)calloc(voiNr, sizeof(Voi));
  if(data->voi==NULL) return 1;
  data->_voidataNr=voiNr;

  /* Allocate memory for frames */
  /* For 3 axes, weights, and for curves (3 for each VOI) */
  /* And one extra 'frame' for overflow testing */
  n=(frameNr+1)*(3+1+3*voiNr);
  data->_data=(double*)calloc(n, sizeof(double));
  if(data->_data==NULL) return 1;
  data->_dataSize=n;

  /* Set pointers for curve data */
  d=data->_data; data->x = d;
  d += frameNr + 1; data->x1 = d;
  d += frameNr + 1; data->x2 = d;
  d += frameNr + 1; data->w = d;
  for(i=0; i<voiNr; i++) {
    d += frameNr + 1; data->voi[i].y = d;
    d += frameNr + 1; data->voi[i].y2 = d;
    d += frameNr + 1; data->voi[i].y3 = d;
  }

  return 0;
}
/*****************************************************************************/

/*****************************************************************************/
/** Adds room for additional VOI TAC(s) into DFT data struct.
 *  Old data is left unchanged.
\return Returns 0 when successful.
 */
int dftAddmem(
  /** Pointer to DFT data struct */
  DFT *dft,
  /** Nr of additional VOI memory blocks */ 
  int voiNr
) {
  Voi *voi2;
  double *data2, *dptr, *newx, *newx1, *newx2, *neww;
  int ri, fi, dataSize2, voidataNr2;

  /* Check the input */
  if(dft==NULL || dft->voi==NULL || dft->frameNr<1 || dft->voiNr<1) return 1;
  if(voiNr<0) return 1; else if(voiNr==0) return 0; 

  /* Allocate memory for new set of curve data (plus additional 'frame') */
  voidataNr2=voiNr+dft->_voidataNr;
  voi2=(Voi*)calloc(voidataNr2, sizeof(Voi));
  if(voi2==NULL) return 3;
  dataSize2=(dft->frameNr+1)*(3*voidataNr2+4);
  data2=(double*)calloc(dataSize2, sizeof(double)); if(data2==NULL) return 3;
  /* Set pointers for new curve data */
  dptr=data2; newx=dptr; newx[dft->frameNr]=0.0;
  dptr+=dft->frameNr+1; newx1=dptr; newx1[dft->frameNr]=0.0;
  dptr+=dft->frameNr+1; newx2=dptr; newx2[dft->frameNr]=0.0;
  dptr+=dft->frameNr+1; neww=dptr;  neww[dft->frameNr]=0.0;
  for(ri=0; ri<voidataNr2; ri++) {
    dptr+=dft->frameNr+1; voi2[ri].y = dptr;
    dptr+=dft->frameNr+1; voi2[ri].y2 = dptr;
    dptr+=dft->frameNr+1; voi2[ri].y3 = dptr;
  }

  /* Copy the original contents */
  for(ri=0; ri<dft->voiNr; ri++) {
    /* Copy Voi header */
    strcpy(voi2[ri].name, dft->voi[ri].name);
    strcpy(voi2[ri].voiname, dft->voi[ri].voiname);
    strcpy(voi2[ri].hemisphere, dft->voi[ri].hemisphere);
    strcpy(voi2[ri].place, dft->voi[ri].place);
    voi2[ri].size=dft->voi[ri].size;
    voi2[ri].sw=dft->voi[ri].sw;
    voi2[ri].sw2=dft->voi[ri].sw2;
    voi2[ri].sw3=dft->voi[ri].sw3;
    /* Copy Voi data */
    for(fi=0; fi<dft->frameNr; fi++) {
      voi2[ri].y[fi]=dft->voi[ri].y[fi];
      voi2[ri].y2[fi]=dft->voi[ri].y2[fi];
      voi2[ri].y3[fi]=dft->voi[ri].y3[fi];
    }
  }
  for(fi=0; fi<dft->frameNr; fi++) {
    newx[fi]=dft->x[fi]; newx1[fi]=dft->x1[fi]; newx2[fi]=dft->x2[fi];
    neww[fi]=dft->w[fi];
  }

  /* Replace original pointers */
  free(dft->_data); dft->_data=data2;
  free(dft->voi); dft->voi=voi2;
  dft->x=newx; dft->x1=newx1; dft->x2=newx2; dft->w=neww;
  dft->_voidataNr=voidataNr2; dft->_dataSize=dataSize2;

  /* Initiate values */
  for(ri=dft->voiNr; ri<dft->_voidataNr; ri++) {
    strcpy(dft->voi[ri].name, "");
    strcpy(dft->voi[ri].voiname, "");
    strcpy(dft->voi[ri].hemisphere, "");
    strcpy(dft->voi[ri].place, "");
    dft->voi[ri].size=1.0;
    dft->voi[ri].sw=dft->voi[ri].sw2=dft->voi[ri].sw3=0;
    for(fi=0; fi<dft->frameNr+1; fi++)
      dft->voi[ri].y[fi]=dft->voi[ri].y2[fi]=dft->voi[ri].y3[fi]=0.0;
  }

  return 0;
}
/*****************************************************************************/

/*****************************************************************************/
/** Add the specified voi [0,voiNr-1] from data2 to data1.
    Allocates memory for additional data VOI, if necessary.
\return Returns 0 if ok.
 */
int dftAdd(
  /** Pointer to DFT struct data */
  DFT *data1,
  /** Pointer to DFT struct data */
  DFT *data2,
  /** Index of TAC in the 2nd DFT */
  int voi
) {
  int i, n;

  if(data1==NULL || data2==NULL) return 1;
  /* Check that voi exists */
  if(data2->voiNr<=voi || voi<0) {
    strcpy(dfterrmsg, "there is no region to combine"); return 8;}

  /* Check that frame number etc is the same */
  if(data1->frameNr!=data2->frameNr ||
     (data1->_type!=DFT_FORMAT_PLAIN && data2->_type!=DFT_FORMAT_PLAIN &&
      (data1->timeunit!=data2->timeunit || strcasecmp(data1->unit, data2->unit))
   )) {
    strcpy(dfterrmsg, "data does not match"); return 8;}

  /* Allocate more memory if necessary */
  if(data1->_voidataNr==data1->voiNr)
    if(dftAddmem(data1, 1)) {
      strcpy(dfterrmsg, "cannot allocate memory"); return 8;}

  /* Copy data */
  n=data1->voiNr;
  (void)dftCopyvoihdr(data2, voi, data1, n);
  for(i=0; i<data1->frameNr; i++) {
    data1->voi[n].y[i]=data2->voi[voi].y[i];
    data1->voi[n].y2[i]=data2->voi[voi].y2[i];
    data1->voi[n].y3[i]=data2->voi[voi].y3[i];
  }
  data1->voiNr+=1;

  /* If data2 contains weights and data1 does not, then copy those */
  if(data2->isweight && !data1->isweight)
    for(i=0; i<data1->frameNr; i++) data1->w[i]=data2->w[i];

  return 0;
}
/*****************************************************************************/

/*****************************************************************************/
/** Select VOIs (sets sw=1), whose names are matching specified string.
    If no string is specified, then all VOIs are selected.
    This function is to replaced by dftSelectRegions().
\return Returns the number of matches, or <0, if an error occurred.
 */
int dftSelect(
  /** Pointer to DFT struct */
  DFT *data,
  /** String to search in TAC name fields */
  char *name
) {
  unsigned int i, j, n;
  char *p, n1[128], n2[128], n3[128], tmp[128], sname[1024];

  if(data==NULL) return -1;
  /* Select all, if no string was specified */
  if(name==NULL || strlen(name)==0) {
    for(i=0; i<(unsigned int)data->voiNr; i++) data->voi[i].sw=1; 
    return data->voiNr;
  }
  /* Make a copy of 'name' and use it */
  strcpy(sname, name);
  /* Check if string contains several substrings (hemisphere and place) */
  n1[0]=n2[0]=n3[0]=(char)0;
  p=strtok(sname, " ,;\n\t|"); if(p!=NULL) strcpy(n1, p); else return -1;
  p=strtok(NULL, " ,;\n\t|"); if(p!=NULL) {
    strcpy(n2, p); p=strtok(NULL, " ,;\n\t|"); if(p!=NULL) strcpy(n3, p);}
  /* Convert strings to lowercase */
  for(i=0; i<strlen(n1); i++) n1[i]=tolower(n1[i]);
  for(i=0; i<strlen(n2); i++) n2[i]=tolower(n2[i]);
  for(i=0; i<strlen(n3); i++) n3[i]=tolower(n3[i]);
  /* Search through the data */
  for(i=0, n=0; i<(unsigned int)data->voiNr; i++) {
    data->voi[i].sw=0;
    snprintf(tmp, 128, "%s%s%s", data->voi[i].voiname, data->voi[i].hemisphere,
                                 data->voi[i].place);
    for(j=0; j<strlen(tmp); j++) tmp[j]=tolower(tmp[j]);
    if(strstr(tmp, n1)==NULL) continue;
    if(n2[0] && strstr(tmp, n2)==NULL) continue;
    if(n3[0] && strstr(tmp, n3)==NULL) continue;
    data->voi[i].sw=1; n++;
  }
  return n;
}
/*****************************************************************************/

/*****************************************************************************/
/** Select the VOIs that have matching region name or number.
    Sets sw=1 or sw=0. This function will replace dftSelect().
\return Returns the number of selected VOIs, or <0 in case of an error.
 */
int dftSelectRegions(
  /** Pointer to DFT data where VOIs are selected */
  DFT *dft,
  /** Name or VOI number which is searched */
  char *region_name,
  /** 1=Non-matching VOIs are deselected, 0=Old selections are preserved */
  int reset
) {
  int ri, match_nr=0;

  /* Check the input */
  if(dft==NULL || dft->voiNr<1 || strlen(region_name)<1) return(-1);
  /* Reset all selections if required */
  if(reset!=0) for(ri=0; ri<dft->voiNr; ri++) dft->voi[ri].sw=0;
  /* Check each VOI */
  for(ri=0; ri<dft->voiNr; ri++) {
    if(rnameMatch(dft->voi[ri].name, ri+1, region_name)!=0) {
      dft->voi[ri].sw=1; match_nr++;
    }
  }
  return(match_nr);
}
/*****************************************************************************/

/*****************************************************************************/
/** Select the best reference region in case that several were found
 *  with dftSelectRegions.
\return Returns the index of best region, or <0 in case of an error.
 */
int dftSelectBestReference(
  /** Pointer to DFT struct, after using dftSelectRegions() */
  DFT *dft
) {
  int ri, len, min_len, i;
  
  if(dft==NULL || dft->voiNr<1) return -1;
  for(ri=0, i=-1, min_len=9999; ri<dft->voiNr; ri++) if(dft->voi[ri].sw) {
    len=strlen(dft->voi[ri].voiname);
    if(strcmp(dft->voi[ri].hemisphere, ".")!=0 &&
       strcasecmp(dft->voi[ri].hemisphere, "AVG")!=0 &&
       strcasecmp(dft->voi[ri].hemisphere, "MEAN")!=0)
      len+=1+strlen(dft->voi[ri].hemisphere);
    if(strcmp(dft->voi[ri].place, ".")!=0 &&
       strcasecmp(dft->voi[ri].place, "ALL")!=0 &&
       strcasecmp(dft->voi[ri].place, "AVG")!=0 &&
       strcasecmp(dft->voi[ri].place, "MEAN")!=0)
      len+=1+strlen(dft->voi[ri].place);
    if(len<min_len) {min_len=len; i=ri;}
  }
  if(i<0) return -2; else return i;
}
/*****************************************************************************/

/*****************************************************************************/
/** Calculate frame mid or start and end times. Timetype is not changed. */
void dftFrametimes(
  /** Pointer to DFT struct */
  DFT *data
) {
  int i, j;
  double f, fs;

  if(data==NULL) return;
  /* If data is told to contain frame start and end times, then check
     that those really are there, or set to middle times, if necessary */
  if(data->timetype==DFT_TIME_STARTEND) {
    for(i=j=0; i<data->frameNr; i++) {
      fs=data->x2[i]-data->x1[i]; if(fs>1.0E-10) {j=1; break;}
    }
    if(j==0) {
      for(i=0; i<data->frameNr; i++) data->x[i]=0.5*(data->x1[i]+data->x2[i]);
      data->timetype=DFT_TIME_MIDDLE;
    }
  }

  /* Decide what to do, and get it done */
  if(data->timetype==DFT_TIME_MIDDLE) {
    /* frame start and end times from mid times */
    /* Easy, if only one frame */
    if(data->frameNr==1) {
      if(data->x[0]<=0.0) {data->x1[0]=data->x[0]; data->x2[0]=0.0;}
      else {data->x1[0]=0.0; data->x2[0]=2.0*data->x[0];}
      return;
    }
    /* Fill start and end times with -999 */
    for(i=0; i<data->frameNr; i++) data->x1[i]=data->x2[i]=-999.;
    /* Search for sequences of nearly same frame lengths */
    for(i=1; i<data->frameNr-1; i++) {
      f=data->x[i]-data->x[i-1]; fs=data->x[i+1]-data->x[i];
      if((f+fs)<=0.0 && fabs(fs-f)>=2.0) continue;
      if((f+fs)>0.0 && (2.0*fabs(fs-f)/(f+fs))>0.1) continue; 
      //if(fabs(f-fs)>=2.0) continue; 
      f=(f+fs)/2.0;
      data->x1[i-1]=data->x[i-1]-f/2.0; data->x2[i-1]=data->x[i-1]+f/2.0;
      data->x1[i]=data->x[i]-f/2.0; data->x2[i]=data->x[i]+f/2.0;
      data->x1[i+1]=data->x[i+1]-f/2.0; data->x2[i+1]=data->x[i+1]+f/2.0;
      /* Check for negatives */
      for(j=i-1; j<i+2; j++) {
        if(data->x1[j]<0.0) data->x1[j]=0.0;
        if(data->x2[j]<0.0) data->x2[j]=0.0;
      }
    }
    /* If out-of-sequence frames were left out, fill those to the nearest one */
    i=0;  /* first frame */
    if(data->x1[i]<0) {
      if(data->x1[i+1]>0) data->x2[i]=data->x1[i+1];
      else data->x2[i]=(data->x[i+1]+data->x[i])/2.0;
      data->x1[i]=2.0*data->x[i]-data->x2[i];
    }
    i=data->frameNr-1;  /* last frame */
    if(data->x1[i]<0) {
      if(data->x2[i-1]>0) data->x1[i]=data->x2[i-1];
      else data->x1[i]=(data->x[i-1]+data->x[i])/2.0;
      data->x2[i]=2.0*data->x[i]-data->x1[i];
    }
    /* other frames */
    for(i=1; i<data->frameNr-1; i++) if(data->x1[i]<0.0) {
      /* which frame is nearest? */
      if(data->x[i]-data->x[i-1] <= data->x[i+1]-data->x[i]) { /* last one */
        if(data->x2[i-1]>0) data->x1[i]=data->x2[i-1];
        else data->x1[i]=(data->x[i-1]+data->x[i])/2.0;
        data->x2[i]=2.*data->x[i]-data->x1[i];
      } else { /* next one */
        if(data->x1[i+1]>0) data->x2[i]=data->x1[i+1];
        else data->x2[i]=(data->x[i+1]+data->x[i])/2.0;
        data->x1[i]=2.0*data->x[i]-data->x2[i];
      }
    }
    /* Check for negatives */
    for(i=0; i<data->frameNr; i++) {
      if(data->x1[i]<0.0) data->x1[i]=0.0;
      if(data->x2[i]<0.0) data->x2[i]=data->x1[i];
    }
    /* Check for overlapping and very small gaps */
    for(i=1; i<data->frameNr; i++) {
      f=data->x1[i]-data->x2[i-1];
      if(f<0.0) {
        if(data->x[i]>data->x2[i-1]) data->x1[i]=data->x2[i-1];
        else if(data->x[i-1]<data->x1[i]) data->x2[i-1]=data->x1[i];
        else data->x1[i]=data->x2[i-1]=(data->x[i]+data->x[i-1])/2.0;
      } else if(f>0.0 && f<1.0) {
        data->x1[i]=data->x2[i-1]=(data->x1[i]+data->x2[i-1])/2.0;
      }
    }
  } else if(data->timetype==DFT_TIME_STARTEND) {
    /* mid times from frame start and end times */
    for(i=0; i<data->frameNr; i++) data->x[i]=0.5*(data->x1[i]+data->x2[i]);
  } else if(data->timetype==DFT_TIME_START) {
    /* frame start times -> end and mid times */
    for(i=0; i<data->frameNr-1; i++) data->x2[i]=data->x1[i+1];
    data->x2[data->frameNr-1]=data->x1[data->frameNr-1]+
      (data->x2[data->frameNr-2]-data->x1[data->frameNr-2]);
    for(i=0; i<data->frameNr; i++) data->x[i]=0.5*(data->x1[i]+data->x2[i]);
  } else if(data->timetype==DFT_TIME_END) {
    /* frame end times -> start and mid times */
    data->x1[0]=0.0;
    for(i=1; i<data->frameNr; i++) data->x1[i]=data->x2[i-1];
    for(i=0; i<data->frameNr; i++) data->x[i]=0.5*(data->x1[i]+data->x2[i]);
  }

  return;
}
/*****************************************************************************/

/*****************************************************************************/
/** Check for overflows in data structure. Returns 0, if ok. */
int dftOverflow(
  /** Pointer to DFT struct */
  DFT *data
) {
  int i;

  if(data==NULL || data->frameNr<1 || data->voiNr<1) return 0;
  if(data->x[data->frameNr]!=0.0) return 1;
  if(data->x1[data->frameNr]!=0.0) return 2;
  if(data->x2[data->frameNr]!=0.0) return 3;
  for(i=0; i<data->voiNr; i++) {
    if(data->voi[i].y[data->frameNr]!=0.0) return 4;
    if(data->voi[i].y2[data->frameNr]!=0.0) return 5;
    if(data->voi[i].y3[data->frameNr]!=0.0) return 6;
  }
  return 0;
}
/*****************************************************************************/

/*****************************************************************************/
/** Copy VOI data inside DFT data structure from one place to another. */
int dftCopyvoi(
  /** Pointer to DFT struct */
  DFT *data,
  /** TAC index */
  int from,
  /** TAC index */
  int to
) {
  int i;

  /* Check that required data exists */
  if(data==NULL || to>=data->_voidataNr || from>=data->_voidataNr) return 1;
  if(from==to) return 0;

  /* Copy VOI info */
  strcpy(data->voi[to].name, data->voi[from].name);
  strcpy(data->voi[to].voiname, data->voi[from].voiname);
  strcpy(data->voi[to].hemisphere, data->voi[from].hemisphere);
  strcpy(data->voi[to].place, data->voi[from].place);
  data->voi[to].size=data->voi[from].size;
  data->voi[to].sw=data->voi[from].sw;
  data->voi[to].sw2=data->voi[from].sw2;
  data->voi[to].sw3=data->voi[from].sw3;
  /* Copy VOI curves */
  for(i=0; i<data->frameNr; i++) {
    data->voi[to].y[i]=data->voi[from].y[i];
    data->voi[to].y2[i]=data->voi[from].y2[i];
    data->voi[to].y3[i]=data->voi[from].y3[i];
  }

  return 0;
}
/*****************************************************************************/

/*****************************************************************************/
/** Move VOI in DFT structure from one position to another. */
int dftMovevoi(
  /** Pointer to DFT struct */
  DFT *dft,
  /** TAC index */
  int from,
  /** TAC index */
  int to
) {
  int ri;
  size_t voisize;
  Voi voi;

  if(dft==NULL || from<0 || to<0) return(1);
  if(from+1>dft->_voidataNr || to+1>dft->_voidataNr) return(2);
  if(from==to) return(0);
  voisize=sizeof(Voi);
  memcpy(&voi, dft->voi+from, voisize);
  if(from>to) for(ri=from; ri>to; ri--)
    memcpy(dft->voi+ri, dft->voi+(ri-1), voisize);
  else for(ri=from; ri<to; ri++)
    memcpy(dft->voi+ri, dft->voi+(ri+1), voisize);
  memcpy(dft->voi+ri, &voi, voisize);
  return(0);
}
/*****************************************************************************/

/*****************************************************************************/
/** Delete specified TAC (0..voiNr-1) from the DFT structure.
\return Returns 0 if ok.
 */
int dftDelete(
  /** Pointer to DFT struct */
  DFT *dft,
  /** TAC index */
  int voi
) {
  int ret;

  /* Check that region exists */
  if(dft==NULL || voi>dft->voiNr-1 || voi<0) return(1);
  /* If it is the last one, then just decrease the voiNr */
  if(voi==dft->voiNr-1) {dft->voiNr--; return(0);}
  /* Otherwise move it to the last position, and then decrease voiNr */
  ret=dftMovevoi(dft, voi, dft->voiNr-1); if(ret) return(10+ret);
  dft->voiNr--;
  return(0);
}
/*****************************************************************************/

/*****************************************************************************/
/** Copy main header info from dft1 to dft2.
\return Returns <> 0 in case of an error.
 */
int dftCopymainhdr(
  /** Pointer to DFT struct from where information is copied */
  DFT *dft1,
  /** Pointer to DFT struct into which information is copied to */
  DFT *dft2
) {
  if(dft1==NULL || dft2==NULL) return 1;
  strcpy(dft2->studynr, dft1->studynr);
  strcpy(dft2->unit, dft1->unit);
  dft2->timeunit=dft1->timeunit; dft2->timetype=dft1->timetype;
  strcpy(dft2->comments, dft1->comments);
  strcpy(dft2->radiopharmaceutical, dft1->radiopharmaceutical);
  strcpy(dft2->isotope, dft1->isotope);
  strcpy(dft2->scanStartTime, dft1->scanStartTime);
  strcpy(dft2->injectionTime, dft1->injectionTime);
  dft2->decayCorrected=dft1->decayCorrected;
  dft2->_type=dft1->_type;
  return 0;
}
/*****************************************************************************/

/*****************************************************************************/
/** Copy main header info from dft1 to dft2. Comments are not copied, because
 *  those may contain outdated units and other information.
\return Returns <> 0 in case of an error.
 */
int dftCopymainhdr2(
  /** Pointer to DFT struct from where information is copied */
  DFT *dft1,
  /** Pointer to DFT struct into which information is copied to */
  DFT *dft2,
  /** Existing header field content is overwritten (1) or kept (0) */
  int ow
) {
  if(dft1==NULL || dft2==NULL) return 1;
  if(ow || (strlen(dft2->studynr)<1 && strcmp(dft2->studynr, ".")==0))
    strcpy(dft2->studynr, dft1->studynr);
  if(ow || dftUnitId(dft2->unit)==CUNIT_UNKNOWN)
    strcpy(dft2->unit, dft1->unit);
  if(ow || dft2->timeunit==TUNIT_UNKNOWN)
    dft2->timeunit=dft1->timeunit;
  dft2->timetype=dft1->timetype;
  if(ow || strlen(dft2->radiopharmaceutical)<1)
    strcpy(dft2->radiopharmaceutical, dft1->radiopharmaceutical);
  if(ow || strlen(dft2->isotope)<1)
    strcpy(dft2->isotope, dft1->isotope);
  if(ow || strlen(dft2->scanStartTime)<1)
    strcpy(dft2->scanStartTime, dft1->scanStartTime);
  if(ow || strlen(dft2->injectionTime)<1)
    strcpy(dft2->injectionTime, dft1->injectionTime);
  if(ow || dft2->decayCorrected==DFT_DECAY_UNKNOWN)
    dft2->decayCorrected=dft1->decayCorrected;
  if(ow)
    dft2->_type=dft1->_type;
  return 0;
}
/*****************************************************************************/

/*****************************************************************************/
/** Copy voi header info from dft1.voi[from] to dft2.voi[to].
\return Returns <> 0 in case of an error.
 */
int dftCopyvoihdr(
  /** Pointer to DFT struct */
  DFT *dft1,
  /** TAC index */
  int from,
  /** Pointer to DFT struct */
  DFT *dft2,
  /** TAC index */
  int to
) {
  /* Check that required data exists */
  if(dft1==NULL || dft2==NULL) return 1;
  if(to>=dft2->_voidataNr || from>=dft1->_voidataNr) return 1;

  /* Copy VOI info */
  strcpy(dft2->voi[to].name, dft1->voi[from].name);
  strcpy(dft2->voi[to].voiname, dft1->voi[from].voiname);
  strcpy(dft2->voi[to].hemisphere, dft1->voi[from].hemisphere);
  strcpy(dft2->voi[to].place, dft1->voi[from].place);
  dft2->voi[to].size=dft1->voi[from].size;
  dft2->voi[to].sw=dft1->voi[from].sw;
  dft2->voi[to].sw2=dft1->voi[from].sw2;
  dft2->voi[to].sw3=dft1->voi[from].sw3;

  return 0;
}
/*****************************************************************************/

/*****************************************************************************/
/** Makes a duplicate of DFT structure pointed to by dft1 into dft2.
\return Returns 0 if ok.
 */
int dftdup(
  /** Pointer to DFT struct */
  DFT *dft1,
  /** Pointer to initiated DFT struct; any existing content of dft2
   *  will be deleted. */
  DFT *dft2
) {
  int ri, fi, ret;

  if(dft1==NULL || dft2==NULL) return 1;
  /* Empty the new data */
  dftEmpty(dft2);
  /* Is that it? Is there any contents in dft1? */
  if(dft1->voiNr==0 && dft1->frameNr==0) {
    ret=dftCopymainhdr(dft1, dft2); return(ret);
  }
  /* Allocate memory for dft2 */
  ret=dftSetmem(dft2, dft1->frameNr, dft1->voiNr); if(ret) return(ret);
  dft2->voiNr=dft1->voiNr; dft2->frameNr=dft1->frameNr;
  /* Copy the contents */
  ret=dftCopymainhdr(dft1, dft2); if(ret) return(ret);
  for(ri=0; ri<dft1->voiNr; ri++) {
    ret=dftCopyvoihdr(dft1, ri, dft2, ri); if(ret) return(ret);
    for(fi=0; fi<dft1->frameNr; fi++) {
      dft2->voi[ri].y[fi]=dft1->voi[ri].y[fi];
      dft2->voi[ri].y2[fi]=dft1->voi[ri].y2[fi];
      dft2->voi[ri].y3[fi]=dft1->voi[ri].y3[fi];
    }
  }
  for(fi=0; fi<dft1->frameNr; fi++) {
    dft2->x[fi]=dft1->x[fi];
    dft2->x1[fi]=dft1->x1[fi]; dft2->x2[fi]=dft1->x2[fi];
    dft2->w[fi]=dft1->w[fi];
  }
  dft2->isweight=dft1->isweight;
  return(0);
}
/*****************************************************************************/

/*****************************************************************************/
/** Allocates a DFT structure with specified size, containing no TAC data but
    header information as available in another DFT struct.

    Any existing content of dft2 will be deleted. Dft2 must be initiated.

    @return Returns 0 if ok.
 */
int dftAllocateWithHeader(
  /** Pointer to initiated DFT struct which will be allocated here;
   *  any previous contents will be deleted. */
  DFT *dft,
  /** Nr of frames to be allocated */
  int frameNr,
  /** Nr of planes to be allocated */
  int voiNr,
  /** Pointer to DFT struct where header contents will be copied from */
  DFT *dft_from
) {
  int ri, fi, ret;

  /* Check the input */
  if(dft==NULL || dft_from==NULL || frameNr<1 || voiNr<0) return 1;
  /* Empty the new data */
  dftEmpty(dft);
  /* Allocate memory for dft */
  ret=dftSetmem(dft, frameNr, voiNr); if(ret) return(ret);
  dft->voiNr=voiNr; dft->frameNr=frameNr;
  /* Copy the contents */
  ret=dftCopymainhdr(dft_from, dft); if(ret) return(ret);
  if(dft->voiNr==dft_from->voiNr) {
    for(ri=0; ri<dft->voiNr; ri++) {
      ret=dftCopyvoihdr(dft_from, ri, dft, ri); if(ret) return(ret);
      if(dft->frameNr==dft_from->frameNr) {
        for(fi=0; fi<dft->frameNr; fi++) {
          dft->voi[ri].y[fi]=dft_from->voi[ri].y[fi];
          dft->voi[ri].y2[fi]=dft_from->voi[ri].y2[fi];
          dft->voi[ri].y3[fi]=dft_from->voi[ri].y3[fi];
        }
      }
    }
  }
  if(dft->frameNr==dft_from->frameNr) {
    for(fi=0; fi<dft->frameNr; fi++) {
      dft->x[fi]=dft_from->x[fi];
      dft->x1[fi]=dft_from->x1[fi]; dft->x2[fi]=dft_from->x2[fi];
      dft->w[fi]=dft_from->w[fi];
    }
    dft->isweight=dft_from->isweight;
  }
  return(0);
}
/*****************************************************************************/

/*****************************************************************************/
/** Include a frame with time 0, unless one already exists.
    @return Returns <> 0 in case of an error.
 */
int dftAddnullframe(
  /** Pointer to DFT struct */
  DFT *data
) {
  int i, j, n;
  DFT temp;


  /* Check whether nullframe exists */
  if(data==NULL) return 1;
  if(data->frameNr<1 || data->x[0]==0.0) return 0;

  /* Allocate memory for temp data */
  dftInit(&temp);
  if(dftSetmem(&temp, data->frameNr, data->voiNr)) return 1;
  temp.frameNr=data->frameNr; temp.voiNr=data->voiNr;

  /* Copy data to temp */
  strcpy(temp.studynr, data->studynr);
  strcpy(temp.unit, data->unit);
  temp.timeunit=data->timeunit; temp.timetype=data->timetype;
  temp.isweight=data->isweight; temp._type=data->_type;
  strcpy(temp.comments, data->comments);
  for(j=0; j<data->frameNr; j++) {
    temp.x[j]=data->x[j]; temp.x1[j]=data->x1[j]; temp.x2[j]=data->x2[j];
    temp.w[j]=data->w[j];
    for(i=0; i<data->voiNr; i++) {
      temp.voi[i].y[j]=data->voi[i].y[j];
      temp.voi[i].y2[j]=data->voi[i].y2[j];
      temp.voi[i].y3[j]=data->voi[i].y3[j];
    }
  }
  for(i=0; i<data->voiNr; i++) {
    strcpy(temp.voi[i].name, data->voi[i].name);
    strcpy(temp.voi[i].voiname, data->voi[i].voiname);
    strcpy(temp.voi[i].hemisphere, data->voi[i].hemisphere);
    strcpy(temp.voi[i].place, data->voi[i].place);
    temp.voi[i].size=data->voi[i].size;
    temp.voi[i].sw=data->voi[i].sw;
    temp.voi[i].sw2=data->voi[i].sw2;
    temp.voi[i].sw3=data->voi[i].sw3;
  }

  /* Reallocate memory for data */
  dftEmpty(data);
  if(dftSetmem(data, temp.frameNr+1, temp.voiNr)) {dftEmpty(&temp); return 2;}

  /* Set nullframe */
  data->x[0]=data->x1[0]=data->x2[0]=0.0; data->w[0]=0.0;
  for(i=0; i<temp.voiNr; i++)
    data->voi[i].y[0]=data->voi[i].y2[0]=data->voi[i].y3[0]=0.0;

  /* Copy data back from temp */
  strcpy(data->studynr, temp.studynr);
  data->voiNr=temp.voiNr;
  strcpy(data->unit, temp.unit);
  data->timeunit=temp.timeunit; data->timetype=temp.timetype;
  data->isweight=temp.isweight; data->_type=temp._type;
  strcpy(data->comments, temp.comments);
  for(j=0, n=1; j<temp.frameNr; j++) {
    if(temp.x[j]<0.0) continue;
    if(n==1) data->x2[0]=temp.x1[j];
    data->x[n]=temp.x[j]; data->x1[n]=temp.x1[j]; data->x2[n]=temp.x2[j];
    data->w[n]=temp.w[j];
    for(i=0; i<temp.voiNr; i++) {
      data->voi[i].y[n]=temp.voi[i].y[j];
      data->voi[i].y2[n]=temp.voi[i].y2[j];
      data->voi[i].y3[n]=temp.voi[i].y3[j];
    }
    n++;
  }
  data->frameNr=n;
  for(i=0; i<temp.voiNr; i++) {
    strcpy(data->voi[i].name, temp.voi[i].name);
    strcpy(data->voi[i].voiname, temp.voi[i].voiname);
    strcpy(data->voi[i].hemisphere, temp.voi[i].hemisphere);
    strcpy(data->voi[i].place, temp.voi[i].place);
    data->voi[i].size=temp.voi[i].size;
    data->voi[i].sw=temp.voi[i].sw;
    data->voi[i].sw2=temp.voi[i].sw2;
    data->voi[i].sw3=temp.voi[i].sw3;
  }

  dftEmpty(&temp);

  return 0;
}
/*****************************************************************************/

/*****************************************************************************/
/** Sort DFT regions in alphabetical order by their name.
    @return Returns <> 0 in case of an error.
 */
int dftSort(
  /** Pointer to DFT struct */
  DFT *data
) {
  if(data==NULL) return(1);
  if(data->voiNr<=1) return(0);
  qsort(data->voi, data->voiNr, sizeof(Voi), dftQSortName);
  return(0);
}
/// @cond
int dftQSortName(const void *voi1, const void *voi2)
{
  int res;

  res=strcasecmp( ((Voi*)voi1)->name, ((Voi*)voi2)->name );
  if(res!=0) return(res);
  res=strcasecmp( ((Voi*)voi1)->voiname, ((Voi*)voi2)->voiname );
  if(res!=0) return(res);
  res=strcasecmp( ((Voi*)voi1)->hemisphere, ((Voi*)voi2)->hemisphere );
  if(res!=0) return(res);
  res=strcasecmp( ((Voi*)voi1)->place, ((Voi*)voi2)->place );
  return(res);
}
/// @endcond
/*****************************************************************************/

/*****************************************************************************/
/** Sort DFT regions in alphabetical order by their plane.
    @return Returns <> 0 in case of an error.
 */
int dftSortPlane(
  /** Pointer to DFT struct */
  DFT *data
) {
  if(data==NULL) return(1);
  if(data->voiNr<=1) return(0);
  qsort(data->voi, data->voiNr, sizeof(Voi), dftQSortPlane);
  return(0);
}
/// @cond
int dftQSortPlane(const void *voi1, const void *voi2)
{
  int res;

  res=strcasecmp( ((Voi*)voi1)->place, ((Voi*)voi2)->place );
  if(res!=0) return(res);
  res=strcasecmp( ((Voi*)voi1)->name, ((Voi*)voi2)->name );
  if(res!=0) return(res);
  res=strcasecmp( ((Voi*)voi1)->voiname, ((Voi*)voi2)->voiname );
  if(res!=0) return(res);
  res=strcasecmp( ((Voi*)voi1)->hemisphere, ((Voi*)voi2)->hemisphere );
  return(res);
}
/// @endcond
/*****************************************************************************/

/*****************************************************************************/
/** Check DFT for NA's in sample times and values.
    @return Returns the number of NA's that were found.
 */
int dft_nr_of_NA(
  /** Pointer to DFT struct */
  DFT *dft
) {
  int ri, fi, na_nr=0;
  if(dft==NULL) return 0;
  for(fi=0; fi<dft->frameNr; fi++) {
    if(dft->timetype==DFT_TIME_STARTEND) {
      if(isnan(dft->x1[fi])) na_nr++;
      if(isnan(dft->x2[fi])) na_nr++;
    } else {
      if(isnan(dft->x[fi])) na_nr++;
    }
    for(ri=0; ri<dft->voiNr; ri++) if(isnan(dft->voi[ri].y[fi])) na_nr++;
  }
  return(na_nr);
}
/*****************************************************************************/

/*****************************************************************************/
/** Replace NA's in basic DFT data with interpolated values.
    If extrapolation is necessary, then the values (0,0) and
    (Infinity,last measured) are assumed.
\return Returns 0, if NA's could be filled with sensible values.
 */
int dftNAfill(
  /** Pointer to DFT struct */
  DFT *dft
) {
  int ri, fi, fj;
  double x1, x2, y1, y2, x, y;

  if(dft==NULL || dft->voiNr<1 || dft->frameNr<1) return(1);
  for(ri=0; ri<dft->voiNr; ri++) for(fi=0; fi<dft->frameNr; fi++) {
    if(isnan(dft->x[fi])) return(2);
    if(isnan(dft->voi[ri].y[fi])) {
      /* NA's before zero time are always replaced with 0 */
      if(dft->x[fi]<0.0) {dft->voi[ri].y[fi]=0.0; continue;}
      x=dft->x[fi];
      /* Get the previous data that is not NA */
      for(x1=y1=nan(""), fj=fi-1; fj>=0; fj--) if(!isnan(dft->voi[ri].y[fj])) {
        x1=dft->x[fj]; y1=dft->voi[ri].y[fj]; break;
      }
      if(isnan(x1) || isnan(y1)) x1=y1=0.0;
      /* Get the following data that is not NA */
      for(x2=y2=nan(""), fj=fi+1; fj<dft->frameNr; fj++) if(!isnan(dft->voi[ri].y[fj])) {
        x2=dft->x[fj]; y2=dft->voi[ri].y[fj]; break;
      }
      if(isnan(x2) || isnan(y2)) for(fj=fi-1; fj>=0; fj--) if(!isnan(dft->voi[ri].y[fj])) {
        x2=dft->x[fj]; y2=dft->voi[ri].y[fj]; break;
      }
      if(isnan(x2) || isnan(y2)) return(2);
      /* Calculate new value */
      if(x2==x1) y=0.5*(y1+y2); else y=y2-(x2-x)*(y2-y1)/(x2-x1);
      dft->voi[ri].y[fi]=y;
    }
  }
  return(0);
}
/*****************************************************************************/

/*****************************************************************************/
/** Search the min and max values of DFT TAC data. Data may contain NA's.

    Note that minx and maxx are the smallest and highest x values in data, not
    the x values at y minimum and maximum; use dftMinMaxTAC() for that.
    @sa dftMinMaxTAC, dftRobustMinMaxTAC
    @return Returns 0 if successful.
 */
int dftMinMax(
  /** Pointer to the DFT TAC data to search */
  DFT *dft,
  /** Pointer to min X; set to NULL if not needed */
  double *minx,
  /** Pointer to max X; set to NULL if not needed */
  double *maxx,
  /** Pointer to min Y; set to NULL if not needed */
  double *miny,
  /** Pointer to max Y; set to NULL if not needed */
  double *maxy
) {
  int ri, fi, n;
  double x1, x2, y1, y2;

  if(dft==NULL) return(1);
  x1=x2=y1=y2=nan("");
  for(fi=0; fi<dft->frameNr; fi++) {
    for(ri=0, n=0; ri<dft->voiNr; ri++) if(!isnan(dft->voi[ri].y[fi])) {
      if(isnan(y1) || y1>dft->voi[ri].y[fi]) y1=dft->voi[ri].y[fi];
      if(isnan(y2) || y2<dft->voi[ri].y[fi]) y2=dft->voi[ri].y[fi];
      n++;
    }
    if(n==0) continue; // no true y values, thus do not use x either
    if(dft->timetype==DFT_TIME_STARTEND) {
      if(!isnan(dft->x1[fi])) {
        if(isnan(x1) || x1>dft->x1[fi]) x1=dft->x1[fi];
      }
      if(!isnan(dft->x2[fi])) {
        if(isnan(x2) || x2<dft->x2[fi]) x2=dft->x2[fi];
      }
    } else if(!isnan(dft->x[fi])) {
      if(isnan(x1) || x1>dft->x[fi]) x1=dft->x[fi];
      if(isnan(x2) || x2<dft->x[fi]) x2=dft->x[fi];
    }
  }
  if(minx!=NULL) {if(isnan(x1)) return(3); else *minx=x1;}
  if(maxx!=NULL) {if(isnan(x2)) return(4); else *maxx=x2;}
  if(miny!=NULL) {if(isnan(y1)) return(5); else *miny=y1;}
  if(maxy!=NULL) {if(isnan(y2)) return(6); else *maxy=y2;}
  return(0);
}
/*****************************************************************************/

/*****************************************************************************/
/** Search the min and max values of DFT TAC data. Data may contain NA's.
    This is not a replacement of dftMinMax() which is needed e.g. in plotting functions.
    @sa dftMinMax, dftRobustMinMaxTAC
    @return Returns 0 if successful.
 */
int dftMinMaxTAC(
  /** Pointer to the DFT TAC data to search */
  DFT *dft,
  /** Index of the only TAC which is searched for min and max; <0 if all */
  int tacindex,
  /** Pointer to X at TAC min; set to NULL if not needed */
  double *minx,
  /** Pointer to X at TAC max; set to NULL if not needed */
  double *maxx,
  /** Pointer to min Y; set to NULL if not needed */
  double *miny,
  /** Pointer to max Y; set to NULL if not needed */
  double *maxy,
  /** Index of min TAC; set to NULL if not needed */
  int *mini,
  /** Index of max TAC; set to NULL if not needed */
  int *maxi,
  /** Index of min sample; set to NULL if not needed */
  int *mins,
  /** Index of max sample; set to NULL if not needed */
  int *maxs
) {
  int ri, fi, i1, i2, s1, s2;
  double x, x1, x2, y1, y2;

  if(dft==NULL) return(1);
  if(tacindex>=dft->voiNr) return(2);
  if(dft->voiNr<1 || dft->frameNr<1) return(3);

  x1=x2=y1=y2=nan(""); i1=i2=s1=s2=0;
  for(fi=0; fi<dft->frameNr; fi++) {
    if(dft->timetype==DFT_TIME_STARTEND) {
      if(isnan(dft->x1[fi])) continue;
      if(isnan(dft->x2[fi])) continue;
      x=0.5*(dft->x1[fi]+dft->x2[fi]);
    } else {
      if(isnan(dft->x[fi])) continue;
      x=dft->x[fi];
    }
    for(ri=0; ri<dft->voiNr; ri++) if(!isnan(dft->voi[ri].y[fi])) {
      if(tacindex>=0 && ri!=tacindex) continue;
      if(isnan(y1) || y1>dft->voi[ri].y[fi]) {
        y1=dft->voi[ri].y[fi]; i1=ri; x1=x; s1=fi;}
      if(isnan(y2) || y2<dft->voi[ri].y[fi]) {
        y2=dft->voi[ri].y[fi]; i2=ri; x2=x; s2=fi;}
    }
  }
  if(minx!=NULL) {if(isnan(x1)) return(11); else *minx=x1;}
  if(maxx!=NULL) {if(isnan(x2)) return(12); else *maxx=x2;}
  if(miny!=NULL) {if(isnan(y1)) return(13); else *miny=y1;}
  if(maxy!=NULL) {if(isnan(y2)) return(14); else *maxy=y2;}
  if(mini!=NULL) {if(isnan(y1)) return(13); else *mini=i1;}
  if(maxi!=NULL) {if(isnan(y2)) return(14); else *maxi=i2;}
  if(mins!=NULL) {if(isnan(y1)) return(13); else *mins=s1;}
  if(maxs!=NULL) {if(isnan(y2)) return(14); else *maxs=s2;}
  return(0);
}
/*****************************************************************************/

/*****************************************************************************/
/** Search the min and max values of DFT TAC data inside specified time range.

    Data may contain NA's.

    @return Returns 0 if successful.
 */
int dftMaxY(
  /** Pointer to the DFT TAC data to search */
  DFT *dft,
  /** Start time */
  double t1,
  /** End time */
  double t2,
  /** Pointer to min Y; set to NULL if not needed */
  double *miny,
  /** Pointer to max Y; set to NULL if not needed */
  double *maxy
) {
  int ri, fi;
  double x1, x2, y1, y2;

  if(dft==NULL) return(1);
  y1=y2=nan("");
  for(fi=0; fi<dft->frameNr; fi++) {
    if(dft->timetype==DFT_TIME_STARTEND) {
      if(!isfinite(dft->x1[fi]) || !isfinite(dft->x2[fi])) continue;
      x1=dft->x1[fi]; x2=dft->x2[fi];
    } else {
      if(!isfinite(dft->x[fi])) continue;
      x1=x2=dft->x[fi];
    }
    if(x2<t1 || x1>t2) continue; // outside time range
    for(ri=0; ri<dft->voiNr; ri++) if(!isnan(dft->voi[ri].y[fi])) {
      if(isnan(y1) || y1>dft->voi[ri].y[fi]) y1=dft->voi[ri].y[fi];
      if(isnan(y2) || y2<dft->voi[ri].y[fi]) y2=dft->voi[ri].y[fi];
    }
  }
  if(miny!=NULL) {if(isnan(y1)) return(5); else *miny=y1;}
  if(maxy!=NULL) {if(isnan(y2)) return(6); else *maxy=y2;}
  return(0);
}
/*****************************************************************************/

/*****************************************************************************/
/** Returns the lowest activity value in DFT */
double dft_kBqMin(
  /** Pointer to DFT struct */
  DFT *data
) {
  int i, j;
  double min=1e+99;

  if(data==NULL) return(nan(""));
  for(i=0; i<data->voiNr ;i++) {
    for(j=0; j<data->frameNr; j++) {
      if(!isnan(data->voi[i].y[j]) && data->voi[i].y[j]<min)
	min=data->voi[i].y[j];
    }
  }
  return(min);
}
/*****************************************************************************/

/*****************************************************************************/
/** Returns the highest activity value in DFT */
double dft_kBqMax(
  /** Pointer to DFT struct */
  DFT *data
) {
  int i, j;
  double max=-1e+99;

  if(data==NULL) return(nan(""));
  for(i=0; i<data->voiNr ;i++){
    for(j=0; j<data->frameNr; j++){
      if(!isnan(data->voi[i].y[j]) && data->voi[i].y[j]>max)
	max=data->voi[i].y[j];
    }
  }
  return(max);
}
/*****************************************************************************/

/*****************************************************************************/
/** Sorts TAC frames by increasing sample time.
    @return Returns 0 if ok.
 */
int dftSortByFrame(
  /** Pointer to DFT struct */
  DFT *dft
) {
  int ri, fi, fj;
  double d;

  if(dft==NULL || dft->voiNr<1 || dft->frameNr<1) return(1);
  for(fi=0; fi<dft->frameNr-1; fi++) for(fj=fi+1; fj<dft->frameNr; fj++) {
    if(dft->x[fj]>=dft->x[fi]) continue;
    d=dft->x[fi];  dft->x[fi]=dft->x[fj];   dft->x[fj]=d;
    d=dft->x1[fi]; dft->x1[fi]=dft->x1[fj]; dft->x1[fj]=d;
    d=dft->x2[fi]; dft->x2[fi]=dft->x2[fj]; dft->x2[fj]=d;
    d=dft->w[fi];  dft->w[fi]=dft->w[fj];   dft->w[fj]=d;
    for(ri=0; ri<dft->voiNr; ri++) {
      d=dft->voi[ri].y[fi]; dft->voi[ri].y[fi]=dft->voi[ri].y[fj]; dft->voi[ri].y[fj]=d;
      d=dft->voi[ri].y2[fi]; dft->voi[ri].y2[fi]=dft->voi[ri].y2[fj]; dft->voi[ri].y2[fj]=d;
      d=dft->voi[ri].y3[fi]; dft->voi[ri].y3[fi]=dft->voi[ri].y3[fj]; dft->voi[ri].y3[fj]=d;
    }
  }
  return(0);
}
/*****************************************************************************/

/*****************************************************************************/
/** Correct frame start and end times if frames are slightly overlapping or
    have small gaps in between.
    Large gap is not corrected and it does not lead to an error.
    @return If overlap is considerable (>1 s), or another error is encountered,
    function returns a non-zero value. Otherwise 0 is returned.
 */
int dftDeleteFrameOverlap_old(
  /** Pointer to DFT data. Time unit must be set, otherwise no checking is done.
   *  Timetype must be DFT_TIME_STARTEND, i.e. both frame start and end time
   *  must be present. */
  DFT *dft
) {
  int fi;
  double overlap, overlap_limit=1.8, flen1, flen2;

  if(dft==NULL) return(1);
  if(dft->timetype!=DFT_TIME_STARTEND) return(0);
  if(dft->timeunit!=TUNIT_MIN && dft->timeunit!=TUNIT_SEC) return(0);
  if(dft->timeunit==TUNIT_MIN) overlap_limit/=60.0;
  for(fi=0; fi<dft->frameNr-1; fi++) {
    overlap=dft->x2[fi] - dft->x1[fi+1];
    if(overlap==0.0) continue; // no gap or overlap
    else if(overlap<-overlap_limit) continue; // gap is large, then do nothing
    else if(overlap>overlap_limit) return(2); // overlap is large: error
    /* Correct the small gap/overlap by making frame durations more similar */
    flen1=dft->x2[fi]-dft->x1[fi]; flen2=dft->x2[fi+1]-dft->x1[fi+1];
    if(overlap>0.0) { // overlap
      if(flen1>flen2) dft->x2[fi]=dft->x1[fi+1]; else dft->x1[fi+1]=dft->x2[fi];
    } else { // gap
      if(flen1>flen2) dft->x1[fi+1]=dft->x2[fi]; else dft->x2[fi]=dft->x1[fi+1];
    }
  }
  return(0);
}
/*****************************************************************************/

/******************************************************************************/
/** Correct frame start and end times if frames are slightly overlapping or
    have small gaps in between. Gap before the first time frame is not corrected.
    Large gap is not corrected and it does not lead to an error.
    @sa dftFillInitialGap
    @return If overlap is considerable (>20%), or another error is encountered,
    function returns a non-zero value. Otherwise 0 is returned.
 */
int dftDeleteFrameOverlap(
  /** Pointer to DFT data. Data must be sorted by increasing time. 
      Time unit does not need to be set.
      Timetype must be DFT_TIME_STARTEND, i.e. both frame start and end time
      must be present; if not, then return value is always 0 (passed). */
  DFT *dft
) {
  int fi;
  double overlap, overlap_limit=0.0, flen1, flen2;

  if(dft==NULL) return(1);
  if(dft->timetype!=DFT_TIME_STARTEND) return(0);
  for(fi=0; fi<dft->frameNr-1; fi++) {
    overlap=dft->x2[fi] - dft->x1[fi+1];
    if(overlap==0.0) continue; // no gap or overlap
    /* Calculate the frame length of current frame and the next frame */
    flen1=dft->x2[fi]-dft->x1[fi]; flen2=dft->x2[fi+1]-dft->x1[fi+1];
    if(flen1<0.0 || flen2<0.0) return(1);
    /* Set the limit */
    if(flen1<flen2) overlap_limit=0.2*flen1; else overlap_limit=0.2*flen2;
    /* Check if gap or overlap is too large to be fixed automatically */
    if(overlap<-overlap_limit) continue; // gap is too large, then do nothing
    if(overlap>overlap_limit) return(2); // overlap is too large: error
    /* Correct the small gap/overlap by making frame durations more similar */
    if(overlap>0.0) { // overlap
      if(flen1>flen2) dft->x2[fi]=dft->x1[fi+1]; else dft->x1[fi+1]=dft->x2[fi];
    } else { // gap
      if(flen1>flen2) dft->x1[fi+1]=dft->x2[fi]; else dft->x2[fi]=dft->x1[fi+1];
    }
  }
  return(0);
}
/******************************************************************************/

/******************************************************************************/
/** Extract specified sample time interval from TAC data.
    @return Returns 0 when successful, otherwise <>0.
 */
int dftRemoveTimeRange(
  /** Pointer to DFT struct from where samples outside time range
      will be removed. Data must be sorted to increasing frame times before
      calling this function. */
  DFT *dft,
  /** Start time of data that is preserved (in same units as in DFT) */
  double startT,
  /** End time of data that is preserved (in same units as in DFT) */
  double endT
) {
  int i, j, voi, first, origNr;

  if(dft==NULL || dft->frameNr<1 || dft->voiNr<1) return(1);
  if(endT<startT) return(2);
  if(startT<=dft->x[0] && endT>=dft->x[dft->frameNr-1]) return(0);
  origNr=dft->frameNr;

  /* Delete end frames which are collected later than the end time */
  for(j=dft->frameNr-1; j>=0; j--) if(dft->x[j]<=endT) break;
  dft->frameNr=j+1;

  /* Find the first frame that has been collected later than start time */
  for(j=0, first=-1; j<dft->frameNr; j++)
    if(dft->x[j]>=startT) {first=j; break;}
  if(first<0) {
    dft->frameNr=origNr; // undelete the end data
    return(3);
  }

  /* Delete first frames */
  if(first>0) for(j=first, i=0; j<dft->frameNr; j++, i++) {
    dft->x[i]=dft->x[j]; dft->x1[i]=dft->x1[j]; dft->x2[i]=dft->x2[j];
    dft->w[i]=dft->w[j];
    for(voi=0; voi<dft->voiNr; voi++) {
      dft->voi[voi].y[i]=dft->voi[voi].y[j];
      dft->voi[voi].y2[i]=dft->voi[voi].y2[j];
      dft->voi[voi].y3[i]=dft->voi[voi].y3[j];
    }
  }
  dft->frameNr-=first;

  return(0);
}
/******************************************************************************/

/******************************************************************************/
/** Overwrites DFT comments with information in current DFT header.

    If DFT format specifies that titles are to be saved, that
    information is not included in comments. 
*/ 
void dftSetComments(
  /** Pointer to DFT struct */
  DFT *dft
) {
  char tmp[512];

  if(dft==NULL) return;
  strcpy(dft->comments, "");
  /* Write in comments the information that will not be included in titles */
  if(dft->scanStartTime[0]) {
    sprintf(tmp, "# scan_start_time := %s\n", dft->scanStartTime);
    strcat(dft->comments, tmp);
  }
  if(dft->injectionTime[0]) {
    sprintf(tmp, "# injection_time := %s\n", dft->injectionTime);
    strcat(dft->comments, tmp);
  }
  strcpy(tmp, "# decay_correction := ");
  if(dft->decayCorrected==DFT_DECAY_CORRECTED) strcat(tmp, "Yes\n");
  else if(dft->decayCorrected==DFT_DECAY_NOTCORRECTED) strcat(tmp, "No\n");
  else strcat(tmp, "Unknown\n");
  if(dft->decayCorrected!=DFT_DECAY_UNKNOWN) strcat(dft->comments, tmp);
  if(dft->isotope[0]) {
    sprintf(tmp, "# isotope := %s\n", dft->isotope);
    strcat(dft->comments, tmp);
  }
  if(dft->radiopharmaceutical[0]) {
    sprintf(tmp, "# radiopharmaceutical := %s\n", dft->radiopharmaceutical);
    strcat(dft->comments, tmp);
  }
  /* If titles are set to be saved, then there's no need to put more in comments */
  if(dft->_type==DFT_FORMAT_STANDARD || dft->_type==DFT_FORMAT_PMOD) return;
  
  /* Ok then, lets write even title information in comments */
  if(dft->studynr[0]) {
    sprintf(tmp, "# study_number := %s\n", dft->studynr);
    strcat(dft->comments, tmp);
  }
  if(dft->timeunit!=TUNIT_UNKNOWN) {
    sprintf(tmp, "# timeunit := %s\n", petTunit(dft->timeunit) );
    strcat(dft->comments, tmp);
  }
  if(petCunitId(dft->unit)!=CUNIT_UNKNOWN) {
    sprintf(tmp, "# unit := %s\n", dft->unit );
    strcat(dft->comments, tmp);
  }
  // Region names and volumes are not saved in comments because of space limit
  
  return;
}
/*****************************************************************************/

/*****************************************************************************/
/** Check if there is a time gap between time zero and first sample time;
 *  if gap does not exist, then nothing is done; if gap exists, then gap
 *  is filled with an extra frame.
 *  @sa dftDeleteFrameOverlap
 *  @return Returns zero if successful, otherwise <>0.
 */
int dftFillInitialGap(
  /** Pointer to DFT struct */
  DFT *dft
) {
  DFT temp;
  int ret, ri, fi;

  /* Check input */
  if(dft==NULL) return 1;
  if(dft->frameNr<1 || dft->voiNr<1) return 0;

  /* Is there an initial gap? If not then we can finish here */
  if(dft->timetype==DFT_TIME_STARTEND) {
    if(dft->x1[0]<=0.0) return 0;
  } else {
    if(dft->x[0]<=0.0) return 0;
  }

  /* Make a temporary storage of the data */
  dftInit(&temp); ret=dftdup(dft, &temp); if(ret!=0) return 10+ret;
  /* Delete and reallocate the original data pointer */
  dftEmpty(dft); ret=dftSetmem(dft, temp.frameNr+1, temp.voiNr);
  if(ret!=0) {
    dftdup(&temp, dft); dftEmpty(&temp);
    return 20+ret;
  }
  /* Copy data back, leaving first frame empty */
  ret=dftCopymainhdr(&temp, dft);
  if(ret!=0) {
    dftdup(&temp, dft); dftEmpty(&temp);
    return 30+ret;
  }
  dft->voiNr=temp.voiNr; dft->frameNr=1+temp.frameNr;
  dft->isweight=temp.isweight;
  strcpy(dft->comments, temp.comments);
  for(ri=0; ri<temp.voiNr; ri++) {
    ret=dftCopyvoihdr(&temp, ri, dft, ri);
    if(ret!=0) {
      dftdup(&temp, dft); dftEmpty(&temp);
      return 40+ret;
    }
    for(fi=0; fi<temp.frameNr; fi++) {
      dft->voi[ri].y[fi+1]=temp.voi[ri].y[fi];
      dft->voi[ri].y2[fi+1]=temp.voi[ri].y2[fi];
      dft->voi[ri].y3[fi+1]=temp.voi[ri].y3[fi];
    }
    /* Fill in the first frame */
    dft->voi[ri].y[0]=0.0;
    dft->voi[ri].y2[0]=0.0;
    dft->voi[ri].y3[0]=0.0;
  }
  for(fi=0; fi<temp.frameNr; fi++) {
    dft->x1[fi+1]=temp.x1[fi];
    dft->x2[fi+1]=temp.x2[fi];
    dft->x[fi+1]=temp.x[fi];
    dft->w[fi+1]=temp.w[fi];
  }
  dftEmpty(&temp);

  /* Fill the first frame */
  dft->w[0]=1.0;
  if(dft->timetype==DFT_TIME_STARTEND) {
    dft->x1[0]=0.0;
    dft->x2[0]=dft->x1[1];
    dft->x[0]=0.5*(dft->x1[0]+dft->x2[0]);
  } else {
    dft->x1[0]=0.0;
    dft->x2[0]=dft->x1[1];
    dft->x[0]=0.0;
  }

  return 0;
}
/*****************************************************************************/

/*****************************************************************************/
/** Add space for additional frames into DFT, keeping the existing data.
 *  frameNr is increased by nr_to_add, but new last frame(s) are empty.
\return Returns zero if successful, otherwise <>0.
 */
int dftAddSpaceForFrames(
  /** Allocated and data-filled DFT where new frame(s) are added to the end */
  DFT *dft,
  /** Nr of frames to add */
  int nr_to_add
) {
  DFT temp;
  int ret, ri, fi;

  /* Check input */
  if(dft==NULL || dft->frameNr<1 || dft->voiNr<1) return 1;
  if(nr_to_add<1) return 0;

  /* Make a temporary storage of the data */
  dftInit(&temp); ret=dftdup(dft, &temp); if(ret!=0) return 10+ret;
  /* Delete and reallocate the original data pointer */
  dftEmpty(dft); ret=dftSetmem(dft, temp.frameNr+nr_to_add, temp.voiNr);
  if(ret!=0) {
    dftdup(&temp, dft); dftEmpty(&temp);
    return 20+ret;
  }
  /* Copy data back, leaving last frame(s) empty */
  ret=dftCopymainhdr(&temp, dft);
  if(ret!=0) {
    dftdup(&temp, dft); dftEmpty(&temp);
    return 30+ret;
  }
  dft->voiNr=temp.voiNr; dft->frameNr=nr_to_add+temp.frameNr;
  dft->isweight=temp.isweight;
  strcpy(dft->comments, temp.comments);
  for(ri=0; ri<temp.voiNr; ri++) {
    ret=dftCopyvoihdr(&temp, ri, dft, ri);
    if(ret!=0) {
      dftdup(&temp, dft); dftEmpty(&temp);
      return 40+ret;
    }
    for(fi=0; fi<temp.frameNr; fi++) {
      dft->voi[ri].y[fi]=temp.voi[ri].y[fi];
      dft->voi[ri].y2[fi]=temp.voi[ri].y2[fi];
      dft->voi[ri].y3[fi]=temp.voi[ri].y3[fi];
    }
    /* Fill in the first frame */
    dft->voi[ri].y[0]=0.0;
    dft->voi[ri].y2[0]=0.0;
    dft->voi[ri].y3[0]=0.0;
  }
  for(fi=0; fi<temp.frameNr; fi++) {
    dft->x1[fi]=temp.x1[fi];
    dft->x2[fi]=temp.x2[fi];
    dft->x[fi]=temp.x[fi];
    dft->w[fi]=temp.w[fi];
  }

  /* Fill the last frames with NAs */
  for(fi=temp.frameNr; fi<dft->frameNr; fi++) {
    dft->w[fi]=1.0;
    dft->x1[fi]=dft->x2[fi]=dft->x[0]=0.0;
    for(ri=0; ri<dft->voiNr; ri++) dft->voi[ri].y[fi]=nan("");
  }
  dftEmpty(&temp);

  return 0;
}
/*****************************************************************************/

/*****************************************************************************/
/** Simplify TAC names in DFT struct: empty hemisphere and/or place field
 *  in case those are the same in all TACs.
 */
void dftRNameSimplify(
  /** Pointer to DFT struct */
  DFT *dft,
  /** Is hemisphere field simplified (1) or not (0), when possible */
  int hemisphere,
  /** Is place field simplified (1) or not (0), when possible */
  int place
) {
  int ri, n;

  if(dft==NULL || dft->voiNr<1) return;
  if(hemisphere==0 && place==0) return;

  /* Check if all TACs have the same field content: if yes, then
     delete the field content */ 
  if(hemisphere!=0) {
    for(ri=n=1; ri<dft->voiNr; ri++)
      if(strcasecmp(dft->voi[0].hemisphere, dft->voi[ri].hemisphere)==0) n++;
    if(n==dft->voiNr)
      for(ri=0; ri<dft->voiNr; ri++)
        strcpy(dft->voi[ri].hemisphere, "");
  }
  if(place!=0) {
    for(ri=n=1; ri<dft->voiNr; ri++)
      if(strcasecmp(dft->voi[0].place, dft->voi[ri].place)==0) n++;
    if(n==dft->voiNr)
      for(ri=0; ri<dft->voiNr; ri++)
        strcpy(dft->voi[ri].place, "");
  }

  /* Construct combined TAC names */
  for(ri=0; ri<dft->voiNr; ri++)
    rnameCatenate(dft->voi[ri].name, MAX_REGIONNAME_LEN,
                  dft->voi[ri].voiname, dft->voi[ri].hemisphere,
                  dft->voi[ri].place, '_');
  /* Ready */
  return;
}
/*****************************************************************************/

/*****************************************************************************/
/** Calculates mean TAC of all TACs in DFT struct.
 *  Mean is NOT weighted by VOI sizes.
 *  Also SD and CV for each sample time are calculated.
\return Returns 0 if successful.
 */
int dftMeanTAC(
  /** Pointer to TAC data from which mean TAC is calculated;
   *  missing values (NaN) are allowed */
  DFT *dft,
  /** Pointer to initialized or pre-allocated DFT struct in where mean, SD,
   *  and CV will be written in y, y2, and y3, respectively. */
  DFT *mean
) {
  int ret, fi, ri, n;
  double sum, ssum;

  if(dft==NULL || mean==NULL) return(1);
  if(dft->voiNr<1 || dft->frameNr<1) return(2);

  /* Allocate memory for mean data, if necessary */
  if(mean->voiNr<1 || mean->frameNr!=dft->frameNr) {
    dftEmpty(mean); ret=dftAllocateWithHeader(mean, dft->frameNr, 1, dft);
    if(ret!=0) return(100+ret);
  }
  strcpy(mean->voi[0].name, "Mean");
  strcpy(mean->voi[0].voiname, mean->voi[0].name);

  /* Calculate the mean TAC */
  ret=0;
  for(fi=0; fi<dft->frameNr; fi++) {
    sum=ssum=0.0; n=0;
    for(ri=0; ri<dft->voiNr; ri++) if(!isnan(dft->voi[ri].y[fi])) {
      sum+=dft->voi[ri].y[fi];
      ssum+=dft->voi[ri].y[fi]*dft->voi[ri].y[fi];
      n++;
    }
    if(n==0) {
      mean->voi[0].y[fi]=mean->voi[0].y2[fi]=mean->voi[0].y3[fi]=nan("");
    } else {
      mean->voi[0].y[fi]=sum/(double)n;
      if(n==1) {
        mean->voi[0].y2[fi]=mean->voi[0].y3[fi]=0.0;
      } else {
        mean->voi[0].y2[fi]=sqrt((ssum-sum*sum/(double)n)/(double)(n-1));
        if(fabs(mean->voi[0].y[fi])>1.0E-25)
          mean->voi[0].y3[fi]=fabs(mean->voi[0].y2[fi]/mean->voi[0].y[fi]);
        else
          mean->voi[0].y3[fi]=0.0;
      }
    }
    if(n>0) ret++;
  }
  /* Check that at least half of frames contained acceptable data */
  if(2*ret<dft->frameNr) {dftEmpty(mean); return(10);}

  return(0);
}
/*****************************************************************************/

/*****************************************************************************/
/** Determine the nr of valid data points inside the given time range.
    @return The number of valid data points.
 */
int dftValidNr(
  /** Pointer to DFT struct */
  DFT *dft,
  /** Time range start */
  double tstart, 
  /** Time range stop */
  double tstop,
  /** Index of TAC to use; enter <0 to use all, in which case the minimum
      number of valid points is returned. */
  int index
) {
  if(dft==NULL || dft->voiNr<0 || dft->frameNr<0) return(0);
  if(index>dft->voiNr-1) return(0);
  if(index>=0) { // TAC index given
    int n=0;
    double x;
    for(int i=0; i<dft->frameNr; i++) {
      if(dft->timetype!=DFT_TIME_STARTEND) x=dft->x[i]; 
      else x=0.5*(dft->x1[i]+dft->x2[i]); 
      if(isnan(x) || !isfinite(x)) continue;
      if(x<tstart || x>tstop) continue;
      if(isnan(dft->voi[index].y[i]) || !isfinite(dft->voi[index].y[i])) 
        continue;
      n++;
    }
    return(n);
  }
  /* Get min nr of all TACs */  
  int n=0, minn=0;
  minn=dftValidNr(dft, tstart, tstop, 0);
  for(int j=1; j<dft->voiNr; j++) {
    n=dftValidNr(dft, tstart, tstop, j);
    if(n<minn) minn=n;
  }
  return(minn);
}
/*****************************************************************************/

/*****************************************************************************/
