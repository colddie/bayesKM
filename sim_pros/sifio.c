// This is a local version, changed at line 198
/// @file sifio.c
/// @author Vesa Oikonen, Jarkko Johansson, Harri Merisaari
/// @brief Functions for reading and writing SIF format files.
///
/*****************************************************************************/
#include "libtpcimgio.h"
/*****************************************************************************/

/*****************************************************************************/
/** Reads SIF file contents to the specified data structure.

    Weights are set to 1.

    @return Returns 0 if ok, 1 invalid input, 2 failed to open file,
    3 failed to allocate memory, 4 header parsing error,
    5 wrong file type, 6 failed to allocate memory,
    8 parse error, 9 wrong data format

    @sa sifInit, sifEmpty, sifWrite, sifExistentCounts, textfileReadLines
 */
int sifRead(
  /** SIF file name to be read. */
  char *filename,
  /** Pointer to initiated SIF structure; any existing contents will be deleted. */
  SIF *data
) {
  int i, n, frameNr, yy, mm, dd, h, m, s;
  int ret;
  struct tm *st;
  time_t timet;


  if(SIF_TEST) printf("sifRead(%s, *sif)\n", filename);
  if(filename==NULL || data==NULL) return(1);
  /* Empty data */
  sifEmpty(data);

  /* Read the lines of the SIF */
  STR_TOKEN_LIST tlst; str_token_list_init(&tlst);
  //printf("reading lines\n");
  ret=textfileReadLines(filename, &tlst);
  if(ret>2) {strcpy(siferrmsg, "cannot read file"); return(5);}
  else if(ret>0) {strcpy(siferrmsg, "cannot open file"); return(2);}
  if(0) {
    printf("lineNr:=%d\n", tlst.token_nr);
    for(int i=0; i<tlst.token_nr; i++) printf("token[%d] := '%s'\n", i, tlst.tok[i]);
  }

  /* Remove any comment or empty lines */
  //printf("removing any comment lines\n");
  i=tlst.token_nr-1;
  while(i>=0 && tlst.token_nr>0) {
    //printf("i=%d line='%s'\n", i, tlst.tok[i]);
    if(strnlen(tlst.tok[i], 2)<1 || asciiCommentLine(tlst.tok[i], NULL)) {
      //printf("-> remove\n");
      str_token_list_del(&tlst, 1+i);
    }
    i--;
  }
  if(tlst.token_nr<1) {
    str_token_list_empty(&tlst);
    strcpy(siferrmsg, "wrong format"); return(4);
  }
  if(0) {
    printf("lineNr:=%d\n", tlst.token_nr);
    for(int i=0; i<tlst.token_nr; i++)
      printf("token[%d] := '%s'\n", i, tlst.tok[i]);
  }

  /* Read the title line */
  if(SIF_TEST) printf("SIF title := '%s'\n", tlst.tok[0]);
  n=sscanf(tlst.tok[0], "%d/%d/%d %d:%d:%d %d %d %d %255s %7s",
      &dd, &mm, &yy, &h, &m, &s, &frameNr, &data->colNr, &data->version,
      data->studynr, data->isotope_name);
  //printf("tok_n := %d\n", n);
  if(n<9 || frameNr<1 || data->colNr<2 || data->version!=1) {
    if(SIF_TEST) printf("invalid SIF title line\n");
    strcpy(siferrmsg, "wrong filetype");
    str_token_list_empty(&tlst);
    return(4);
  }
  timet=time(NULL); st=gmtime(&timet);
  if(st!=NULL) {
    st->tm_mday=dd; st->tm_mon=mm-1; st->tm_year=yy-1900;
    st->tm_hour=h; st->tm_min=m; st->tm_sec=s; st->tm_isdst=-1;
    data->scantime=timegm(st); if(data->scantime==-1) data->scantime=0;
  } else {
    data->scantime=0;
  }
  //printf("frameNr := %d\n", frameNr);

  /* Allocate memory for SIF data */
  if(sifSetmem(data, frameNr)) {
    str_token_list_empty(&tlst);
    strcpy(siferrmsg, "cannot allocate SIF"); return(6);
  }

  /* Read data lines into SIF */
  i=0;
  while(i<data->frameNr && i<tlst.token_nr-1) {
    //printf("i := %d\n", i);
    data->prompts[i]=data->randoms[i]=0.0;
    n=sscanf(tlst.tok[1+i], "%lf %lf %lf %lf", &data->x1[i], &data->x2[i],
        &data->prompts[i], &data->randoms[i]);
    //printf("n := %d\n", n);
    //if(n<data->colNr || data->x2[i]<data->x1[i]) { 
    if(n<2) { 
      strcpy(siferrmsg, "wrong data format"); 
      str_token_list_empty(&tlst);
      sifEmpty(data);
      return(9);
    }
    if(data->x2[i]<data->x1[i]) { 
      strcpy(siferrmsg, "invalid time frames"); 
      str_token_list_empty(&tlst);
      sifEmpty(data);
      return(9);
    }
    i++;
  }
  str_token_list_empty(&tlst);
  if(i!=data->frameNr) {
    strcpy(siferrmsg, "wrong data format"); 
    sifEmpty(data);
    return(9);
  }


  /* Calculate trues */
  if(data->colNr>=4) for(i=0; i<data->frameNr; i++)
    data->trues[i]=data->prompts[i]-data->randoms[i];
  /* Set weights to 1.0 */
  for(i=0; i<data->frameNr; i++) data->weights[i]=1.0;

  return(0);
}
/*****************************************************************************/

/*****************************************************************************/
/** Write SIF data to a standard SIF file, emptying files old contents.
    @return Returns 0 if successful, 1 invalid input, 2 failed to open file, 
    3 failed to write into file.
    @sa sifInit, sifEmpty, sifRead, sifPrint, sifWeight
 */
int sifWrite(
  /** Pointer to SIF structure containing data to be written in file. */
  SIF *data,
  /** File name for SIF; Set to 'stdout' to print the contents on console.
      File is overwritten without backup. */
  char *filename
) {
  int i, n, req_decimals=0;
  char buf[1024];
  struct tm tm; //*st;


  if(SIF_TEST) printf("sifWrite(*sif, %s)\n", filename);
  /* Check data */
  if(data->frameNr<1) {strcpy(siferrmsg, "no data to save"); return 1;}

  /* Set file pointer */
  FILE *fp;
  int intofile;
  if(strcasecmp(filename, "STDOUT")==0) {
    fp=stdout; intofile=0;
  } else {
    /* Open file */
    fp=fopen(filename, "w");
    if(fp==NULL) {strcpy(siferrmsg, "cannot open file"); return 2;}
    intofile=1;
  }

  /* Write title line */
  if(!gmtime_r((time_t*)&data->scantime, &tm) || !strftime(buf, 1024, "%d/%m/%Y %H:%M:%S", &tm))
    strcpy(buf, "1/1/1900 00:00:00");
  n=fprintf(fp, "%s %d %d %d", buf, data->frameNr, data->colNr, data->version);
  if(SIF_TEST) printf("%s %d %d %d\n", buf, data->frameNr, data->colNr, data->version);
  if(n<7) {
    strcpy(siferrmsg, "cannot write file"); if(intofile) fclose(fp); 
    return 2;
  }
  if(strlen(data->studynr)!=0 || strlen(data->isotope_name)!=0) {
    /* Write also study number and isotope */
    if(strlen(data->studynr)==0) fprintf(fp, " ."); else fprintf(fp, " %s", data->studynr);
    if(strlen(data->isotope_name)>0) fprintf(fp, " %s", data->isotope_name);
  }
  fprintf(fp, "\n");

  /* Check if frame times need to printed with decimals */
  for(i=1, req_decimals=0; i<data->frameNr; i++) {
    if(round(data->x1[i])==round(data->x1[i-1])) {req_decimals=1; break;}
    if(round(data->x2[i])==round(data->x2[i-1])) {req_decimals=1; break;}
  }

  /* Write data lines */
  for(i=0; i<data->frameNr; i++) {
    if(req_decimals==0) n=fprintf(fp, "%.0f %.0f", data->x1[i], data->x2[i]);
    else n=fprintf(fp, "%.6f %.6f", data->x1[i], data->x2[i]);
    if(n<3) {
      strcpy(siferrmsg, "cannot write file"); if(intofile) fclose(fp); 
      return 3;
    }
    if(data->colNr<=2) {n=fprintf(fp, "\n"); continue;}
    n=fprintf(fp, " %.0f %.0f", data->prompts[i], data->randoms[i]);
    if(n<1) {
      strcpy(siferrmsg, "cannot write file"); if(intofile) fclose(fp); 
      return 3;
    }
    if(data->colNr<5) {n=fprintf(fp, "\n"); continue;}
    n=fprintf(fp, " %.0f", data->trues[i]);
    if(n<1) {
      strcpy(siferrmsg, "cannot write file"); if(intofile) fclose(fp);
      return 3;
    }
    if(data->colNr<6) {n=fprintf(fp, "\n"); continue;}
    n=fprintf(fp, " %.5f\n", data->weights[i]);
    if(n<1) {
      strcpy(siferrmsg, "cannot write file"); if(intofile) fclose(fp);
      return 3;
    }
  }

  /* Close file */
  if(intofile) fclose(fp);

  return 0;
}
/*****************************************************************************/

/*****************************************************************************/
/** Prints to stdout the contents of SIF data structure.
    @sa sifWrite, sifRead
 */
void sifPrint(
  /** Pointer to SIF structure. */
  SIF *data
) {
  char buf[32];
  
  if(!ctime_r_int((time_t*)&data->scantime, buf))
    strcpy(buf, "1900-01-01 00:00:00");
  printf("Scan time: %s\n", buf);
  printf("Isotope: %s\n", data->isotope_name);
  printf("Frame start   end      Prompts    Randoms      Trues   Weight\n");
  for(int i=0; i<data->frameNr; i++) {
    printf(" %03d %6.1f %6.1f  %10.0f %10.0f %10.0f %8.6f\n", i+1,
      data->x1[i], data->x2[i], data->prompts[i], data->randoms[i],
      data->trues[i], data->weights[i]);
  }
  return;
}
/*****************************************************************************/

/*****************************************************************************/
