/// @file mathfunc.c
/// @author Vesa Oikonen
/// @brief IO for FIT files and calculating function values.
///
/*****************************************************************************/
#include "libtpccurveio.h"
#include <unistd.h>
/*****************************************************************************/
#include "libtpcmodel.h"
/*****************************************************************************/

/*****************************************************************************/
/** Free memory allocated for FIT. All contents are cleared.
    @sa fitInit, fitRead, fitWrite
 */
void fitEmpty(
  /** Pointer to FIT struct. */
  FIT *fit
) {
  if(fit==NULL) return;
  if(fit->_voidataNr>0) {
    free((char*)(fit->voi));
    fit->_voidataNr=0;
  }
  fit->voiNr=0;
  fit->datafile[0]=fit->studynr[0]=fit->unit[0]=fit->program[0]=(char)0;
  fit->timeunit=0;
  fit->time=0;
}
/*****************************************************************************/

/*****************************************************************************/
/** Initiate FIT structure. Call this once before first use.
    @sa fitEmpty, fitRead, fitWrite
 */
void fitInit(
  FIT *fit
) {
  if(fit==NULL) return;
  memset(fit, 0, sizeof(FIT));
  fit->_voidataNr=0; fit->voiNr=0;
}
/*****************************************************************************/

/*****************************************************************************/
/** Write function parameters in FIT into specified file.

    If necessary, a backup file (+BACKUP_EXTENSION) is created.
    @return In case of an error, >0 is returned, and a description is written in fiterrmsg.
    @sa fitPrint, fitRead, fitInit
*/
int fitWrite(
  /** Pointer to FIT struct. */
  FIT *fit, 
  /** Filename. */
  char *filename
) {
  int i, j, n, savedNr=0;
  char tmp[1024], is_stdout=0;
  FILE *fp;


  if(MATHFUNC_TEST>0) printf("fitWrite(FIT, %s)\n", filename);
  /* Check that there is some data to write */
  if(fit==NULL || fit->voiNr<1) {strcpy(fiterrmsg, "no data"); return 1;}
  for(i=0; i<fit->voiNr; i++)
    if(!isnan(fit->voi[i].wss) && fit->voi[i].type>0) savedNr++;
  if(savedNr<1) {strcpy(fiterrmsg, "no fitted data"); return 1;}

  /* Check if writing to stdout */
  if(!strcmp(filename, "stdout")) is_stdout=1;
  if(MATHFUNC_TEST>1) {
    if(is_stdout) printf("  output is stdout()\n");
    else printf("  output in file\n"); 
  }

  /* Check if file exists; backup, if necessary */
  if(!is_stdout) (void)backupExistingFile(filename, NULL, NULL);

  /* Open output file */
  if(is_stdout) fp=(FILE*)stdout;
  else {
    if(MATHFUNC_TEST>2) printf("opening file %s for write\n", filename);
    if((fp = fopen(filename, "w")) == NULL) {
      strcpy(fiterrmsg, "cannot open file"); return 2;
    }
  }

  /* Fit file format */
  n=fprintf(fp, "%-11.11s %s\n", FIT_VER, fit->program);
  if(n==0) {
    strcpy(fiterrmsg, "disk full");
    if(!is_stdout) fclose(fp);
    return 3;
  }
  /* Write fit date and time */
  if(!ctime_r_int(&fit->time, tmp)) strcpy(tmp, "");
  fprintf(fp, "Date:\t%s\n", tmp);
  /* Write the name of the original datafile */
  fprintf(fp, "Data file:\t%s\n", fit->datafile);
  /* Write the studynr; NOT YET because all SW reading these should be
     recompiled*/
  //if(fit->studynr[0]) fprintf(fp, "Study:\t%s\n", fit->studynr);
  /* Write the 'activity' unit */
  fprintf(fp, "Data unit:\t%s\n", fit->unit);
  /* Write the time unit */
  if(fit->timeunit==TUNIT_UM || fit->timeunit==TUNIT_MM)
    fprintf(fp, "Distance unit:\t%s\n", petTunit(fit->timeunit));
  else
    fprintf(fp, "Time unit:\t%s\n", petTunit(fit->timeunit));

  /* Write the voiNr to be saved */
  fprintf(fp, "Nr of VOIs:\t%d\n", savedNr);
  /* Write the Fit title */
  fprintf(fp, "%s %s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",
    "Region", "Plane", "Start", "End", "dataNr", "WSS", "parNr", "Type",
    "Parameters");
  /* Write regional fits */
  for(i=0; i<fit->voiNr; i++) {
    if(isnan(fit->voi[i].wss) || fit->voi[i].type<=0) continue;
    if(fit->voi[i].voiname[0]) strcpy(tmp, fit->voi[i].voiname);
    else strcpy(tmp, ".");
    fprintf(fp, "%.*s ", MAX_REGIONSUBNAME_LEN, tmp);
    if(fit->voi[i].hemisphere[0]) strcpy(tmp, fit->voi[i].hemisphere);
    else strcpy(tmp, ".");
    fprintf(fp, "%.*s ", MAX_REGIONSUBNAME_LEN, tmp);
    if(fit->voi[i].place[0]) strcpy(tmp, fit->voi[i].place);
    else strcpy(tmp, ".");
    fprintf(fp, "%.*s\t", MAX_REGIONSUBNAME_LEN, tmp);
    fprintf(fp, "%.3f\t%.3f\t%d\t%.2E\t%d\t%04d",
       fit->voi[i].start, fit->voi[i].end, fit->voi[i].dataNr,
       fit->voi[i].wss, fit->voi[i].parNr, fit->voi[i].type );
    for(j=0; j<fit->voi[i].parNr; j++) fprintf(fp, "\t%.6E", fit->voi[i].p[j]);
    fprintf(fp, "\n");
  }

  /* Close file */
  if(!is_stdout) {
    if(MATHFUNC_TEST>2) printf("closing file %s\n", filename);
    fflush(fp); fclose(fp);
  }
  strcpy(fiterrmsg, "");

  if(MATHFUNC_TEST>0) printf("done with fitWrite()\n");
  return(0);
}
/*****************************************************************************/

/*****************************************************************************/
/** Allocate memory for FIT data. Any previous contents are destroyed.
    @return Returns 0 when successful, and >0 in case of an error.
    @sa fitInit, fitEmpty, fitRead
 */
int fitSetmem(
  /** Pointer to FIT struct. */
  FIT *fit, 
  /** Nr of TACs to allocate. */
  int voiNr
) {
  /* Check that there is something to do */
  if(fit==NULL || voiNr<1) return 1;

  /* Clear previous data, but only if necessary */
  if(fit->_voidataNr>0 || fit->voiNr>0) fitEmpty(fit);

  /* Allocate memory for regional curves */
  fit->voi=(FitVOI*)calloc(voiNr, sizeof(FitVOI));
  if(fit->voi==NULL) return 2;
  fit->_voidataNr=voiNr;

  return 0;
}
/*****************************************************************************/

/*****************************************************************************/
/** Print to stdout the contents of FIT data structure.

    Mainly for testing purposes.
*/
void fitPrint(
  /** Pointer to FIT struct. */
  FIT *fit
) {
  if(fit==NULL) {printf("FIT = Null\n"); return;}
  if(MATHFUNC_TEST>0) printf("Number of curves: %d\n", fit->voiNr);
  if(MATHFUNC_TEST>0) printf("_voidataNr = %d\n", fit->_voidataNr);
  fitWrite(fit, "stdout");
}
/*****************************************************************************/

/*****************************************************************************/
/** Read FIT file contents to the specified data structure, emptying its old contents.
    @return In case of an error, >0 is returned, and a description is written in fiterrmsg.
    @sa fitInit, fitPrint, fitWrite
*/
int fitRead(
  /** Pointer to file name. */
  char *filename,
  /** Pointer to initiated FIT struct. */
  FIT *fit,
  /** Verbose level; if <=0, then nothing is printed into stdout. */
  int verbose
) {
  char *cptr, line[1024], *lptr;
  int i, ri, n, type=0;
  struct tm st;


  if(verbose>0) printf("fitRead(%s)\n", filename);
  if(fit==NULL) return 1;
  /* Empty data */
  fitEmpty(fit);

  /* Open file */
  FILE *fp=fopen(filename, "r");
  if(fp==NULL) {strcpy(fiterrmsg, "cannot open file"); return 1;}

  /* Read data */
  strcpy(fiterrmsg, "wrong format");
  /* Read file type and program name */
  while(fgets(line, 1024, fp)!=NULL) {
    /* Ignore empty and comment lines */
    if(strlen(line)<4 || line[0]=='#') continue;
    /* Check for id string */
    if(!strncmp(line, FIT_VER, strlen(FIT_VER))) type=1; else type=0;
    break;
  }
  if(type!=1) {fclose(fp); return 3;}
  lptr=line; cptr=strtok(lptr, " \t\n\r");
  if(cptr!=NULL) cptr=strtok(NULL, " \t\n\r");
  if(cptr!=NULL && strlen(cptr)<1024) strcpy(fit->program, cptr);
  /* Read fit date and time */
  while(fgets(line, 1024, fp)!=NULL) if(strlen(line)>2 && line[0]!='#') break;
  if(strncasecmp(line, "Date:", 5)) {fclose(fp); return 4;}
  cptr=line+5; while(*cptr && !isdigit(*cptr)) cptr++;
  if(get_datetime(cptr, &st, verbose-3)==0) fit->time=timegm(&st); 
  /* Read the name of the original datafile */
  while(fgets(line, 1024, fp)!=NULL) if(strlen(line)>2 && line[0]!='#') break;
  if(strncasecmp(line, "Data file:", 10)) {fclose(fp); return 5;}
  lptr=&line[10]; cptr=strtok(lptr, " \t\n\r");
  if(cptr!=NULL && strlen(cptr)<FILENAME_MAX) strcpy(fit->datafile, cptr);
  /* Read the activity unit */
  while(fgets(line, 1024, fp)!=NULL) if(strlen(line)>2 && line[0]!='#') break;
  if(strncasecmp(line, "Data unit:", 10)) {fclose(fp); return 6;}
  lptr=&line[10]; cptr=strtok(lptr, " \t\n\r");
  if(cptr!=NULL && strlen(cptr)<1024) strcpy(fit->unit, cptr);
  /* Read the time unit */
  while(fgets(line, 1024, fp)!=NULL) if(strlen(line)>2 && line[0]!='#') break;
  if(strncasecmp(line, "Time unit:", 10)==0) lptr=&line[10];
  else if(strncasecmp(line, "Distance unit:", 14)==0) lptr=&line[14];
  else {fclose(fp); return 7;}
  cptr=strtok(lptr, " \t\n\r");
  fit->timeunit=petTunitId(cptr);
  /* Read the nr of regions */
  while(fgets(line, 1024, fp)!=NULL) if(strlen(line)>2 && line[0]!='#') break;
  if(strncasecmp(line, "Nr of VOIs:", 11)) {fclose(fp); return 8;}
  lptr=&line[11]; cptr=strtok(lptr, " \t\n\r");
  n=atoi(cptr); if(n<1 || n>32000) {fclose(fp); return 8;}
  /* Allocate memory for regions */
  if(fitSetmem(fit, n)) {strcpy(fiterrmsg, "out of memory"); fclose(fp); return 9;}
  fit->voiNr=n;
  /* Read (and ignore) title line */
  strcpy(fiterrmsg, "wrong format");
  while(fgets(line, 1024, fp)!=NULL) if(strlen(line)>2 && line[0]!='#') break;
  if(strncasecmp(line, "Region", 6)) {fclose(fp); return 10;}
  /* Read regional data */
  for(ri=0; ri<fit->voiNr; ri++) {
    /* Read line of data */
    line[0]=(char)0;
    while(fgets(line, 1024, fp)!=NULL) if(strlen(line)>2 && line[0]!='#') break;
    if(!line[0]) break;
    int pindx=0; char seps[8];
    int sn=strTokenNr(line, " \t\n\r"); if(sn<8) {fclose(fp); fitEmpty(fit); return 11;}
    int tn=strTokenNr(line, "\t\n\r");
    if(tn==sn || tn==sn-1 || tn==sn-2) { // tab as separator
      strcpy(seps, "\t\n\r"); n=tn;
      char *s=strTokenDup(line, seps, NULL);
      char *cptr=s;
      int i=0, n=strlen(s);
      strReplaceChar(cptr, ' ', (char)0);
      strlcpy(fit->voi[ri].voiname, cptr, MAX_REGIONSUBNAME_LEN+1);
      i+=strlen(fit->voi[ri].voiname);
      cptr=s+i; if(i<n) {cptr++; i++;}
      strReplaceChar(cptr, ' ', (char)0);
      strlcpy(fit->voi[ri].hemisphere, cptr, MAX_REGIONSUBNAME_LEN+1);
      i+=strlen(fit->voi[ri].hemisphere);
      cptr=s+i; if(i<n) {cptr++; i++;}
      strReplaceChar(cptr, ' ', (char)0);
      strlcpy(fit->voi[ri].place, cptr, MAX_REGIONSUBNAME_LEN+1);
      free(s);
      pindx=2;
    } else { // spaces as separators
      strcpy(seps, " \t\n\r"); n=sn;
      strTokenNCpy(line, seps, 1, fit->voi[ri].voiname, MAX_REGIONSUBNAME_LEN+1);
      strTokenNCpy(line, seps, 2, fit->voi[ri].hemisphere, MAX_REGIONSUBNAME_LEN+1);
      strTokenNCpy(line, seps, 3, fit->voi[ri].place, MAX_REGIONSUBNAME_LEN+1);
      pindx=4;
    }
    snprintf(fit->voi[ri].name, MAX_REGIONNAME_LEN+1, "%.6s %.6s %.6s",
      fit->voi[ri].voiname, fit->voi[ri].hemisphere, fit->voi[ri].place);
    /* Fit start and end times, and original data nr */
    char s[128];
    strTokenNCpy(line, seps, pindx++, s, 128); fit->voi[ri].start=atof_dpi(s);
    strTokenNCpy(line, seps, pindx++, s, 128); fit->voi[ri].end=atof_dpi(s);
    strTokenNCpy(line, seps, pindx++, s, 128); fit->voi[ri].dataNr=atoi(s);
    /* Fit error, parameter nr and function number (type) */
    strTokenNCpy(line, seps, pindx++, s, 128); fit->voi[ri].wss=atof_dpi(s);
    strTokenNCpy(line, seps, pindx++, s, 128); fit->voi[ri].parNr=atoi(s);
    strTokenNCpy(line, seps, pindx++, s, 128); fit->voi[ri].type=atoi(s);
    /* Parameters */
    for(i=0; i<fit->voi[ri].parNr; i++) {
      strTokenNCpy(line, seps, pindx++, s, 128); fit->voi[ri].p[i]=atof_dpi(s);
    }
  }
  if(ri==0) {fclose(fp); fitEmpty(fit); return(12);}
  if(ri<fit->voiNr) fit->voiNr=ri;

  /* Close file */
  fclose(fp);
  strcpy(fiterrmsg, "");

  if(verbose>1) printf("done fitRead()\n");
  return(0);
}
/*****************************************************************************/

/*****************************************************************************/
/** Copies the description of a function type to the specified string
    which must have space for >=128 characters.
*/
int fitFunctionformat(
  /** The number of function */
  int type,
  /** Representation of the format of the function */
  char *str
) {
  strcpy(str, "");
  switch(type) {
    /* Polynomials, including line */
    case 100: strcpy(str, "f(x)=A"); break;
    case 101: strcpy(str, "f(x)=A+B*x"); break;
    case 102: strcpy(str, "f(x)=A+B*x+C*x^2"); break;
    case 103: strcpy(str, "f(x)=A+B*x+C*x^2+D*x^3"); break;
    case 104: strcpy(str, "f(x)=A+B*x+C*x^2+D*x^3+E*x^4"); break;
    case 105: strcpy(str, "f(x)=A+B*x+C*x^2+D*x^3+E*x^4+F*x^5"); break;
    case 106: strcpy(str, "f(x)=A+B*x+C*x^2+D*x^3+E*x^4+F*x^5+G*x^6"); break;
    case 107: strcpy(str, "f(x)=A+B*x+C*x^2+D*x^3+E*x^4+F*x^5+G*x^6+H*x^7"); break;
    case 108: strcpy(str, "f(x)=A+B*x+C*x^2+D*x^3+E*x^4+F*x^5+G*x^6+H*x^7+I*x^8"); break;
    case 109: strcpy(str, "f(x)=A+B*x+C*x^2+D*x^3+E*x^4+F*x^5+G*x^6+H*x^7+I*x^8+J*x^9"); break;
    /* Rational functions */
    case 211: strcpy(str, "f(x)=(A+C*x)/(B+D*x)"); break;
    case 221: strcpy(str, "f(x)=(A+C*x+E*x^2)/(B+D*x)"); break;
    case 222: strcpy(str, "f(x)=(A+C*x+E*x^2)/(B+D*x+F*x^2)"); break;
    case 232: strcpy(str, "f(x)=(A+C*x+E*x^2+G*x^3)/(B+D*x+F*x^2)"); break;
    case 233: strcpy(str, "f(x)=(A+C*x+E*x^2+G*x^3)/(B+D*x+F*x^2+H*x^3)"); break;
    case 1232: strcpy(str, "f(x)=(A+C*(x-t)+E*(x-t)^2+G*(x-t)^3)/(B+D*(x-t)+F*(x-t)^2)"); break;
    /* Exponential functions */
    case 301: strcpy(str, "f(x)=A*exp(B*x)"); break;
    case 302: strcpy(str, "f(x)=A*exp(B*x)+C*exp(D*x)"); break;
    case 303: strcpy(str, "f(x)=A*exp(B*x)+C*exp(D*x)+E*exp(F*x)"); break;
    case 304: strcpy(str, "f(x)=A*exp(B*x)+C*exp(D*x)+E*exp(F*x)+G*exp(H*x)"); break;
    case 305: strcpy(str, "f(x)=A*exp(B*x)+C*exp(D*x)+E*exp(F*x)+G*exp(H*x)+I*exp(J*x)"); break;
    /* Feng function */
    case 1312: strcpy(str, "f(x)=(A*(x-t)-C)*exp(B*(x-t))+C*exp(D*(x-t))"); break;
    case 1313: strcpy(str, "f(x)=(A*(x-t)-C-E)*exp(B*(x-t))+C*exp(D*(x-t))+E*exp(F*(x-t))"); break;
    case 1314: strcpy(str, "f(x)=(A*(x-t)-C-E-G)*exp(B*(x-t))+C*exp(D*(x-t))+E*exp(F*(x-t))+G*exp(H*(x-t))"); break;
    /* Lundqvist function */
    case 321: strcpy(str, "f(x)=A*exp(B*x)*(1-exp(C*x))"); break;
    case 322: strcpy(str, "f(x)=A*exp(B*x)*(1-exp(C*x))+D*exp(E*x)*(1-exp(F*x))"); break;
    case 323: strcpy(str, "f(x)=A*exp(B*x)*(1-exp(C*x))+D*exp(E*x)*(1-exp(F*x))+G*exp(H*x)*(1-exp(I*x))"); break;
    case 1321: strcpy(str, "f(x)=A*exp(-B*(x-t))*(1-exp(-C*(x-t))) + D*(A/(B*(B+C)))*(C-((B+C)*exp(C*(x-t))-B)*exp(-(B+C)*(x-t))) "); break;
    /* Exponential bolus infusion functions */
    case 331: strcpy(str, "f(x)=(1/Ti)*Sum[i=1..n, (Ai/Li)*(exp(-Li*(x-tA-Ti)) - exp(-Li*(x-Ta)))], when x>=Ta+Ti\nf(x)=(1/Ti)*Sum[i=1..n, (Ai/Li)*(1-exp(-Li*(t-Ta)))], when x>Ta and x<Ta+Ti\nf(x)=0, when t<=Ta");
      break;
    /* Kudomi's function for radiowater */
    case 332: strcpy(str, "f(x)=0, when x<=Ta \nf(x)=(A/L^2)*(1-exp(-L(x-Ta))), when Ta<x<=Ta+Ti \nf(x)=(A/L^2)*(exp(-L*Ti)+exp(-L(x-Ta-Ti))-2exp(-L(x-t1))), when x>Ta+Ti");
      break;
    /* Bolus infusion approaching zero */
    case 334: strcpy(str, "f(x)=0, when x<=t1 \nf(x)=A*(1-exp(L(t1-x)))/(1-exp(L*(t1-t2))), when t1<x<=t2 \nf(x)=A*(exp(L(t2-x))-exp(L(t1-x)))/(1-exp(L*(t1-t2))), when x>t2");
      break;
    /* Exponential functions for plasma fractions */
    case 351: strcpy(str, "f(x)=1-a*(2-exp(-b*x)-exp(-c*x))"); break;
    /* Gamma variate function */
    case 1401: strcpy(str, "f(x)=A*((x-D)^B)*exp(-(x-D)/C) , when x>=D, else f(x)=0"); break;
    /* Gamma variate function with background */
    case 1402: strcpy(str, "f(x)=A*((x-D)^B)*exp(-(x-D)/C) + E , when x>=D, else f(x)=E"); break;
    /* Gamma variate bolus plus recirculation function */
    case 1403: strcpy(str, "f(x)=B*((x-A)^C)*exp(-(x-A)/D) + E*(1-exp(-(x-A)/D))*exp(-(x-A)/F) , when x>A, else f(x)=0"); break;
    /* Weibull cdf */
    case 1421: strcpy(str, "f(x)=A*(1-exp(-((x-t)/B)^C) , when x>t, else f(x)=0"); break;
    /* Weibull cdf plus pdf (derivative of cdf) */
    case 1423: strcpy(str, "f(x)=A*[C*((x-t)/B)^(C-1)*exp(-((x-t)/B)^C))/B + k*(1-exp(-((x-t)/B)^C))] , when x>t, else f(x)=0"); break;
    /* Surge function with AUC=A */
    case 1431: strcpy(str, "f(x)=A*x*exp(-B*x)*B^2 , when x>0, else f(x)=0"); break;
    /* Traditional Surge function */
    case 1432: strcpy(str, "f(x)=A*x*exp(-B*x) , when x>0, else f(x)=0"); break;
    /* Surge function with recirculation */
    case 1433: strcpy(str, "f(x)=A*[x*exp(-B*x) + (C/B^2)*(1-(B*x+1)*exp(-B*x))], when x>0, else f(x)=0"); break;
    /* Surge function with recirculation for plasma-to-blood ratio */
    case 1434: strcpy(str, "f(x)=1/(1-H*(1-r(x))), where r(x) is function for RBC-to-plasma"); break;
    /* Hill functions for TACs */
    case 1801: strcpy(str, "f(x)=[A*(x-t)^B]/[(x-t)^B+C^B]"); break;
    case 1811: strcpy(str, "f(x)=A*{[B*(x-t)^(B-1)]/[C^B+(x-t)^B] - [B*(x-t)^(2*B-1)]/[C^B+(x-t)^B]^2}"); break;
    case 1821: strcpy(str, "f(x)=A*{[B*(x-t)^(B-1)]/[C^B+(x-t)^B] - [B*(x-t)^(2*B-1)]/[C^B+(x-t)^B]^2 + K*(x-t)^B]/[(x-t)^B+C^B}"); break;
    /* Hill functions for dose-response curves */
    case 2801: strcpy(str, "f(x)=B+[A-B]/[1+(C/x)^D]"); break;
    case 2802: strcpy(str, "f(x)=B+[A-B]/[1+10^{(C-x)*D}]"); break;
    /* Hill type functions for fractions */
    case 841: strcpy(str, "f(x)=(A*x^B)/(x^B+C)"); break;
    case 842: strcpy(str, "f(x)=1-((A*x^B)/(x^B+C))"); break;
    case 843: strcpy(str, "f(x)=1-((A*(1+D*x)*x^B)/(x^B+C))"); break;
    case 844: strcpy(str, "f(x)=(A*(x-t)^B)/((x-t)^B+C)+D, when x>t, else f(x)=D"); break;
    case 845: strcpy(str, "f(x)=A-(A*x^B)/(x^B+C))"); break;
    case 846: strcpy(str, "f(x)=D+((A-D)*(x-t)^B)/((x-t)^B+C), when x>t, else f(x)=D"); break;
    case 847: strcpy(str, "f(x)=1-D-((A-D)*(x-t)^B)/((x-t)^B+C), when x>t, else f(x)=1-A"); break;
    case 848: strcpy(str, "f(x)=D*((1-A)*(x-t)^B)/((x-t)^B+C), when x>t, else f(x)=D"); break;
    case 849: strcpy(str, "f(x)=1-D*((1-A)*(x-t)^B)/((x-t)^B+C), when x>t, else f(x)=1-A"); break;
    /* Mamede/Watabe function for fractions */
    case 851: strcpy(str, "f(x)=1/(1+(A*x)^2)^B"); break;
    case 852: strcpy(str, "f(x)=1-1/(1+(A*x)^2)^B"); break;
    /* Mamede/Watabe function for fractions, as extended by Meyer */
    case 861: strcpy(str, "f(x)=(1+(A*(x-t))^B)^(-C), when x>t, else f(x)=1"); break;
    case 862: strcpy(str, "f(x)=1-(1+(A*(x-t))^B)^(-C), when x>t, else f(x)=0"); break;
    /* ... and further extended by letting fraction start somewhere between 0 and 1 */ 
    case 863: strcpy(str, "f(x)=(D^(-1/C)+(A*(x-t))^B)^(-C), when x>t, else f(x)=D"); break;
    case 864: strcpy(str, "f(x)=1-(D^(-1/C)+(A*(x-t))^B)^(-C), when x>t, else f(x)=1-D"); break;
    /* Functions for fitting plasma fractions via separate metabolite fractions */
    case 871:
    case 881:
      strcpy(str, "f(x)=1-f1(x)-f2(x)-f3(x)"); break;
    case 872:
    case 882:
      strcpy(str, "f1(x)=a(x)(1-b(x)-c(x)+b(x)c(x))/(1-a(x)b(x)-a(x)c(x)-b(x)c(x)+2a(x)b(x)c(x))"); break;
    case 873:
    case 883:
      strcpy(str, "f2(x)=b(x)(1-a(x)-c(x)+a(x)c(x))/(1-a(x)b(x)-a(x)c(x)-b(x)c(x)+2a(x)b(x)c(x))"); break;
    case 874:
    case 884:
      strcpy(str, "f3(x)=c(x)(1-a(x)-b(x)+a(x)b(x))/(1-a(x)b(x)-a(x)c(x)-b(x)c(x)+2a(x)b(x)c(x))"); break;
    /* PET profile functions */
    case 2111: strcpy(str, "P(x)=(C/2)*(erf((x-d+R)/(sqrt(2)*FWHM/2355))-erf((x-d-R)/(sqrt(2)*FWHM/2355)))+bkg");
      break;
    /* Combined functions and models */
    case 3331: strcpy(str, "f(x)=(1/Ti)*Sum[i=1..n, (Ai/Li)*(exp(-Li*(x-tA-Ti)) - exp(-Li*(x-Ta)))], when x>=Ta+Ti\nf(x)=(1/Ti)*Sum[i=1..n, (Ai/Li)*(1-exp(-Li*(t-Ta)))], when x>Ta and x<Ta+Ti\nf(x)=0, when t<=Ta, with additional delay and dispersion");
      break;
    /* Compartmental model functions */
    /* Graham's plasma curve smoothing function */
    case 9501: strcpy(str, "Cp(t)<=>Ci(t)<=>Ct(t)"); break;
    /* Extended Graham's plasma curve smoothing function */
    case 9502: strcpy(str, "Ce(t)<=>Cp(t)<=>Ci(t)<=>Ct(t)"); break;
    /* Extended Graham's plasma curve smoothing function with metabolite */
    case 9503: strcpy(str, "Cpa(t)<=>Cia(t)<=>Cta(t)->Ctm(t)<=>Cim(t)<=>Cpm(t)"); break;
    /* Huang's plasma metabolite model */
    case 9601: strcpy(str, "C4(t)<=>C3(t)<-C0(t)->C1(t)<=>C2(t)"); break;
    /* Extended Carson's plasma metabolite model */
    case 9602: strcpy(str, "Cpa(t)<=>Cta(t)->Ctm(t)<=>Cpm(t)"); break;
    /* New plasma metabolite model */
    case 9603: strcpy(str, "Cpa(t)->Ct1(t)<=>Cpm(t)<=>Ct2(t)"); break;
    /* Multilinear multicompartmental TAC fitting model */
    case 9701: strcpy(str, "Ideal bolus -> n compartments"); break;
    default:  return(1);
  }
  return(0);
}
/*****************************************************************************/

/*****************************************************************************/
/** Copies the name of the function to the specified string
    which must have space for >=128 characters.
*/
int fitFunctionname(
  /** The number of function */
  int type,
  /** Name of the function */
  char *str
) {
  strcpy(str, "");
  switch(type) {
    case 100: strcpy(str, "f(x)=A"); break;
    case 101: strcpy(str, "line"); break;
    case 102: strcpy(str, "2nd order polynomial"); break;
    case 103: strcpy(str, "3rd order polynomial"); break;
    case 104: strcpy(str, "4th order polynomial"); break;
    case 105: strcpy(str, "5th order polynomial"); break;
    case 106: strcpy(str, "6th order polynomial"); break;
    case 107: strcpy(str, "7th order polynomial"); break;
    case 108: strcpy(str, "8th order polynomial"); break;
    case 109: strcpy(str, "9th order polynomial"); break;
    case 211: strcpy(str, "1/1 order rational function"); break;
    case 221: strcpy(str, "2/1 order rational function"); break;
    case 222: strcpy(str, "2/2 order rational function"); break;
    case 232: strcpy(str, "3/2 order rational function"); break;
    case 233: strcpy(str, "3/3 order rational function"); break;
    case 1232: strcpy(str, "3/2 order rational function with delay"); break;
    case 301: strcpy(str, "exponential function"); break;
    case 302: strcpy(str, "sum of 2 exponential functions"); break;
    case 303: strcpy(str, "sum of 3 exponential functions"); break;
    case 304: strcpy(str, "sum of 4 exponential functions"); break;
    case 305: strcpy(str, "sum of 5 exponential functions"); break;
    case 1312: strcpy(str, "Feng model 2 function with 2 exponentials"); break;
    case 1313: strcpy(str, "Feng model 2 function"); break;
    case 1314: strcpy(str, "Feng model 2 function with 4 exponentials"); break;
    case 321: strcpy(str, "Lundqvist function"); break;
    case 322: strcpy(str, "sum of 2 Lundqvist functions"); break;
    case 323: strcpy(str, "sum of 3 Lundqvist functions"); break;
    case 1321: strcpy(str, "Lundqvist function with integral and delay"); break;
    case 331: strcpy(str, "Exponential bolus infusion function"); break;
    case 332: strcpy(str, "Kudomi's exponential bolus infusion function for radiowater"); break;
    case 334: strcpy(str, "Exponential bolus function approaching zero"); break;
    case 351: strcpy(str, "Exponential function for [C-11]PK11195 plasma fractions"); break;
    case 1401: strcpy(str, "Gamma variate function"); break;
    case 1402: strcpy(str, "Gamma variate with background"); break;
    case 1403: strcpy(str, "Gamma variate bolus plus recirculation"); break;
    case 1421: strcpy(str, "Weibull cdf with delay"); break;
    case 1423: strcpy(str, "Weibull cdf and derivative with delay"); break;
    case 1431: strcpy(str, "Surge function"); break;
    case 1432: strcpy(str, "Surge function (trad)"); break;
    case 1433: strcpy(str, "Surge function with recirculation"); break;
    case 1434: strcpy(str, "Surge function with recirculation for plasma-to-blood ratio"); break;
    case 1801: strcpy(str, "Hill function with delay"); break;
    case 1811: strcpy(str, "Derivative of Hill function with delay"); break;
    case 1821: strcpy(str, "Sum of Hill function and derivative with delay"); break;
    case 2801: strcpy(str, "Hill function for dose-response curve on linear scale"); break;
    case 2802: strcpy(str, "Hill function for dose-response curve on log scale"); break;
    case 841: strcpy(str, "Hill function"); break;
    case 842: strcpy(str, "Hill function (1-f(x))"); break;
    case 843: strcpy(str, "Hill function (1-f(x)) with ascending or descending end"); break;
    case 844: strcpy(str, "Hill function with background"); break;
    case 845: strcpy(str, "Hill function (A-f(x))"); break;
    case 846: strcpy(str, "Extended Hill function for plasma parent fraction"); break;
    case 847: strcpy(str, "Extended Hill function for plasma metabolite fraction"); break;
    case 848: strcpy(str, "Extended Hill function #2 for plasma parent fraction"); break;
    case 849: strcpy(str, "Extended Hill function #2 for plasma metabolite fraction"); break;
    case 851: strcpy(str, "Mamede function"); break;
    case 852: strcpy(str, "Mamede function (1-f(x)"); break;
    case 861: strcpy(str, "Meyer parent fraction function"); break;
    case 862: strcpy(str, "Meyer metabolite fraction function"); break;
    case 863: strcpy(str, "Extended Meyer parent fraction function"); break;
    case 864: strcpy(str, "Extended Meyer metabolite fraction function"); break;
    case 871: strcpy(str, "1-3 metabolite Hill function for parent"); break;
    case 872: strcpy(str, "1-3 metabolite Hill function for metab1"); break;
    case 873: strcpy(str, "1-3 metabolite Hill function for metab2"); break;
    case 874: strcpy(str, "1-3 metabolite Hill function for metab3"); break;
    case 881: strcpy(str, "1-3 metabolite power function for parent"); break;
    case 882: strcpy(str, "1-3 metabolite power function for metab1"); break;
    case 883: strcpy(str, "1-3 metabolite power function for metab2"); break;
    case 884: strcpy(str, "1-3 metabolite power function for metab3"); break;
    case 2111: strcpy(str, "Image profile function"); break;
    case 3331: strcpy(str, "Exponential bolus infusion function with delay and dispersion"); break;
    case 9501: strcpy(str, "Graham's input function"); break;
    case 9502: strcpy(str, "Extended Graham's input function"); break;
    case 9503: strcpy(str, "Graham's input function with metabolite"); break;
    case 9601: strcpy(str, "Huang's plasma metabolite model"); break;
    case 9602: strcpy(str, "Extended Carson's plasma metabolite model"); break;
    case 9603: strcpy(str, "New plasma metabolite model"); break;
    case 9701: strcpy(str, "Multilinear multicompartmental TAC fitting model"); break;
    default:  return(1);
  }
  return(0);
}
/*****************************************************************************/

/*****************************************************************************/
/** Evaluate y=f(x).
    @sa fitEvaltac, fitIntegralEval, fitDerivEval
    @return Returns 0, if ok.
 */
int fitEval(
  /** Fit parameters of a single region */
  FitVOI *r,
  /** Time where to evaluate the function */
  double x,
  /** The value of the function is returned here */
  double *y
) {
  double a, b, f, sqrx, cubx, xt;
  int i, j, m, n;
  double mp[3][5];
  double mf[]={0.0,0.0,0.0};

  if(r==NULL) return(1);
  if(MATHFUNC_TEST>0) {
    printf("fitEval(r, %g, %g): type=%d ;", x, *y, r->type);
    for(i=0; i<r->parNr; i++) printf(" %g", r->p[i]);
    printf("\n");
  }
  if(y==NULL) return(1);
  *y=nan("");
  switch(r->type) {
    case 100:
    case 101:
    case 102:
    case 103:
    case 104:
    case 105:
    case 106:
    case 107:
    case 108:
    case 109:
      n=r->type-99; f=0.0; for(i=n-1; i>0; i--) {f+=r->p[i]; f*=x;} f+=r->p[0];
      *y=f;
      break;
    case 211:
      a=r->p[0]+r->p[2]*x; b=r->p[1]+r->p[3]*x; if(b==0.0) break;
      f=a/b; *y=f;
      break;
    case 221:
      sqrx=x*x;
      a=r->p[0]+r->p[2]*x+r->p[4]*sqrx; b=r->p[1]+r->p[3]*x; if(b==0.0) break;
      f=a/b; *y=f;
      break;
    case 222:
      sqrx=x*x;
      a=r->p[0]+r->p[2]*x+r->p[4]*sqrx;
      b=r->p[1]+r->p[3]*x+r->p[5]*sqrx; if(b==0.0) break;
      f=a/b; *y=f;
      break;
    case 232:
      sqrx=x*x; cubx=sqrx*x;
      a=r->p[0]+r->p[2]*x+r->p[4]*sqrx+r->p[6]*cubx;
      b=r->p[1]+r->p[3]*x+r->p[5]*sqrx; if(b==0.0) break;
      f=a/b; *y=f;
      break;
    case 233:
      sqrx=x*x; cubx=sqrx*x;
      a=r->p[0]+r->p[2]*x+r->p[4]*sqrx+r->p[6]*cubx;
      b=r->p[1]+r->p[3]*x+r->p[5]*sqrx+r->p[7]*cubx; if(b==0.0) break;
      f=a/b; *y=f;
      break;
    case 301:
      f=r->p[0]*exp(r->p[1]*x); *y=f;
      break;
    case 302:
      f=0;
      a=r->p[0]*exp(r->p[1]*x); f+=a;
      a=r->p[2]*exp(r->p[3]*x); f+=a;
      *y=f;
      break;
    case 303:
      f=0;
      a=r->p[0]*exp(r->p[1]*x); f+=a;
      a=r->p[2]*exp(r->p[3]*x); f+=a;
      a=r->p[4]*exp(r->p[5]*x); f+=a;
      *y=f;
      break;
    case 304:
      f=0;
      a=r->p[0]*exp(r->p[1]*x); f+=a;
      a=r->p[2]*exp(r->p[3]*x); f+=a;
      a=r->p[4]*exp(r->p[5]*x); f+=a;
      a=r->p[6]*exp(r->p[7]*x); f+=a;
      *y=f;
      break;
    case 305:
      f=0;
      a=r->p[0]*exp(r->p[1]*x); f+=a;
      a=r->p[2]*exp(r->p[3]*x); f+=a;
      a=r->p[4]*exp(r->p[5]*x); f+=a;
      a=r->p[6]*exp(r->p[7]*x); f+=a;
      a=r->p[8]*exp(r->p[9]*x); f+=a;
      *y=f;
      break;
    case 321:
      f=r->p[0]*exp(r->p[1]*x)*(1.0-exp(r->p[2]*x));
      *y=f;
      break;
    case 322:
      f=0.0;
      a=r->p[0]*exp(r->p[1]*x)*(1.0-exp(r->p[2]*x)); f+=a;
      a=r->p[3]*exp(r->p[4]*x)*(1.0-exp(r->p[5]*x)); f+=a;
      *y=f;
      break;
    case 323:
      f=0.0;
      a=r->p[0]*exp(r->p[1]*x)*(1.0-exp(r->p[2]*x)); f+=a;
      a=r->p[3]*exp(r->p[4]*x)*(1.0-exp(r->p[5]*x)); f+=a;
      a=r->p[6]*exp(r->p[7]*x)*(1.0-exp(r->p[8]*x)); f+=a;
      *y=f;
      break;
    case 1321:
      /* A=p0 B=p1 C=p2 D=k=p3 dT=p4 */
      xt=x-r->p[4]; if(xt<=0.0) {
        f=0.0;
      } else {
        a=exp(-r->p[2]*xt); f=r->p[0]*exp(-r->p[1]*xt)*(1.0-a);
        if(r->p[3]>0.0)
          f+=r->p[3]*(r->p[0]/(r->p[1]*(r->p[1]+r->p[2])))
             *(r->p[2]-((r->p[1]+r->p[2])/a-r->p[1])*exp(-(r->p[1]+r->p[2])*xt));
      }
      *y=f;
      break;
    case 331:
      n=(r->parNr-2)/2;
      if(x<=r->p[0])
        *y=0.0;
      else if(x<r->p[0]+r->p[1]) {
        for(i=0, f=0.0; i<n; i++) { 
          if(r->p[2*i+3]>1.0E-12) b=r->p[2*i+2]/r->p[2*i+3]; else b=r->p[2*i+2];
          f+=b*(1.0-exp(-r->p[2*i+3]*(x-r->p[0])));
        }
        if(r->p[1]>0.0) f*=(1.0/r->p[1]);
        *y=f;
      } else {
        for(i=0, f=0.0; i<n; i++) {
          if(r->p[2*i+3]>1.0E-12) b=r->p[2*i+2]/r->p[2*i+3]; else b=r->p[2*i+2];
          f+=b*(exp(-r->p[2*i+3]*(x-r->p[0]-r->p[1])) - exp(-r->p[2*i+3]*(x-r->p[0])));
        }
        if(r->p[1]>0.0) f*=(1.0/r->p[1]);
        *y=f;
      }
      break;
    case 332:
      if(x<=r->p[0])
        *y=0.0;
      else if(x<=r->p[0]+r->p[1]) {
        f=exp(-r->p[3]*(x-r->p[0]));
        *y=(r->p[2]/(r->p[3]*r->p[3]))*(1.0-f);
      } else {
        f=exp(-r->p[3]*r->p[1]);
        f+=exp(-r->p[3]*(x-r->p[0]-r->p[1]));
        f-=2.0*exp(-r->p[3]*(x-r->p[0]));
        *y=(r->p[2]/(r->p[3]*r->p[3]))*f;
      }
      break;
    case 334:
      if(x<=r->p[0])
        *y=0.0;
      else if(x>r->p[0] && x<=r->p[1]) {
        *y=r->p[2]*(1.0-exp(r->p[3]*(r->p[0]-x)))/(1.0-exp(r->p[3]*(r->p[0]-r->p[1])));
      } else {
        *y=r->p[2]*(exp(r->p[3]*(r->p[1]-x))-exp(r->p[3]*(r->p[0]-x)))/(1.0-exp(r->p[3]*(r->p[0]-r->p[1])));
      }
      break;
    case 351:
      f=1.0-r->p[0]*(2.0-exp(-r->p[1]*x)-exp(-r->p[2]*x)); *y=f;
      break;
    case 841:
      f=r->p[0]*pow(x, r->p[1]) / (pow(x, r->p[1]) + r->p[2]);
      *y=f;
      break;
    case 842:
      f=1.0-r->p[0]*pow(x, r->p[1]) / (pow(x, r->p[1]) + r->p[2]);
      *y=f;
      break;
    case 843:
      f=1.0 - r->p[0]*(1.0+r->p[3]*x)*pow(x, r->p[1]) / (pow(x, r->p[1]) + r->p[2]);
      *y=f;
      break;
    case 844:
      if(r->parNr>4) xt=x-r->p[4]; else xt=x; // delay time accounted, if avail
      if(xt<=0.0) f=r->p[3];
      else {a=pow(xt, r->p[1]);  f=r->p[0]*a/(a + r->p[2]) + r->p[3];}
      *y=f;
      break;
    case 845:
      f=r->p[0] - r->p[0]*pow(x, r->p[1]) / (pow(x, r->p[1]) + r->p[2]);
      *y=f;
      break;
    case 846:
      xt=x-r->p[4]; if(xt<=0.0) {
        f=r->p[3];
      } else {
        a=pow(xt, r->p[1]);
        f=r->p[3]+(r->p[0]-r->p[3])*a/(a+r->p[2]);
      }
      *y=f;
      break;
    case 847:
      xt=x-r->p[4]; if(xt<=0.0) {
        f=1.0-r->p[3];
      } else {
        a=pow(xt, r->p[1]);
        f=1.0-r->p[3]-(r->p[0]-r->p[3])*a/(a+r->p[2]);
      }
      *y=f;
      break;
    case 848:
      xt=x-r->p[4]; if(xt<=0.0) {
        f=r->p[3];
      } else {
        a=pow(xt, r->p[1]);
        f=r->p[3] * (1.0 - (1.0-r->p[0])*a/(r->p[2]+a));
      }
      *y=f;
      break;
    case 849:
      xt=x-r->p[4]; if(xt<=0.0) {
        f=1.0-r->p[3];
      } else {
        a=pow(xt, r->p[1]);
        f=1.0 - r->p[3] * (1.0 - (1.0-r->p[0])*a/(r->p[2]+a));
      }
      *y=f;
      break;
    case 851:
      a=r->p[0]*x; a*=a; f=1.0/pow(1.0+a, r->p[1]);
      *y=f;
      break;
    case 852:
      a=r->p[0]*x; a*=a; f=1.0 - 1.0/pow(1.0+a, r->p[1]);
      *y=f;
      break;
    case 861:
      xt=x-r->p[3]; if(xt<=0.0) {
        f=1.0;
      } else {
        f=pow(1.0+pow(r->p[0]*xt, r->p[1]), -r->p[2]);
      }
      *y=f;
      break;
    case 862:
      xt=x-r->p[3]; if(xt<=0.0) {
        f=1.0;
      } else {
        f=pow(1.0+pow(r->p[0]*xt, r->p[1]), -r->p[2]);
      }
      *y=1.0-f;
      break;
    case 863:
      xt=x-r->p[4]; if(xt<=0.0) {
        f=r->p[3];
      } else {
        f=pow(pow(r->p[3],-1.0/r->p[2])+pow(r->p[0]*xt, r->p[1]), -r->p[2]);
      }
      *y=f;
      break;
    case 864:
      xt=x-r->p[4]; if(xt<=0.0) {
        f=1.0-r->p[3];
      } else {
        f=1.0-pow(pow(r->p[3],-1.0/r->p[2])+pow(r->p[0]*xt, r->p[1]), -r->p[2]);
      }
      *y=f;
      break;
    case 871:
    case 872:
    case 873:
    case 874:
      for(m=0; m<3; m++) for(j=0; j<5; j++) mp[m][j]=0.0;
      for(i=m=j=0; i<r->parNr; i++) {mp[m][j]=r->p[i]; j++; if(j>4) {j=0; m++;}} 
      n=r->parNr/5;
      for(m=0; m<n; m++) {
        xt=x-mp[m][4];
        if(xt<=0.0) mf[m]=1.0-mp[m][3];
        else {
          f=pow(xt, mp[m][1]);
          mf[m]=1.0-(mp[m][3]+(mp[m][0]-mp[m][3])*f/(mp[m][2]+f));
        }
      }
      a=1.0-mf[0]*mf[1]-mf[0]*mf[2]-mf[1]*mf[2]+2.0*mf[0]*mf[1]*mf[2];
      if(r->type==871) {
        f=1.0;
        f-=mf[0]*(1.0-mf[1]-mf[2]+mf[1]*mf[2])/a;
        f-=mf[1]*(1.0-mf[0]-mf[2]+mf[0]*mf[2])/a;
        f-=mf[2]*(1.0-mf[0]-mf[1]+mf[0]*mf[1])/a;
        *y=f;
      } else if(r->type==872) {
        *y=mf[0]*(1.0-mf[1]-mf[2]+mf[1]*mf[2])/a;
      } else if(r->type==873) {
        *y=mf[1]*(1.0-mf[0]-mf[2]+mf[0]*mf[2])/a;
      } else if(r->type==874) {
        *y=mf[2]*(1.0-mf[0]-mf[1]+mf[0]*mf[1])/a;
      }
      break;
    case 881:
    case 882:
    case 883:
    case 884:
      for(m=0; m<3; m++) for(j=0; j<5; j++) mp[m][j]=0.0;
      for(i=m=j=0; i<r->parNr; i++) {mp[m][j]=r->p[i]; j++; if(j>4) {j=0; m++;}} 
      n=r->parNr/5;
      for(m=0; m<n; m++) {
        xt=x-mp[m][4];
        if(xt<=0.0) mf[m]=1.0-mp[m][3];
        else {
          mf[m]=1.0-pow(pow(mp[m][3], -1.0/mp[m][2])+
                                     pow(mp[m][0]*(xt), mp[m][1]), -mp[m][2]);
        }
      }
      a=1.0-mf[0]*mf[1]-mf[0]*mf[2]-mf[1]*mf[2]+2.0*mf[0]*mf[1]*mf[2];
      if(r->type==881) {
        f=1.0;
        f-=mf[0]*(1.0-mf[1]-mf[2]+mf[1]*mf[2])/a;
        f-=mf[1]*(1.0-mf[0]-mf[2]+mf[0]*mf[2])/a;
        f-=mf[2]*(1.0-mf[0]-mf[1]+mf[0]*mf[1])/a;
        *y=f;
      } else if(r->type==882) {
        *y=mf[0]*(1.0-mf[1]-mf[2]+mf[1]*mf[2])/a;
      } else if(r->type==883) {
        *y=mf[1]*(1.0-mf[0]-mf[2]+mf[0]*mf[2])/a;
      } else if(r->type==884) {
        *y=mf[2]*(1.0-mf[0]-mf[1]+mf[0]*mf[1])/a;
      }
      break;
    case 2111:
      xt=x-r->p[3]; a=sqrt(2.0)*(r->p[2]/2.355);
      f=r->p[4]+(r->p[0]/2.0)*(erf((xt+r->p[1])/a)-erf((xt-r->p[1])/a));
      *y=f;
      break;
    case 1232:
      xt=x-r->p[7]; if(xt<=0.0) {
        f=0.0;
      } else {
        sqrx=xt*xt; cubx=sqrx*xt;
        a=r->p[0]+r->p[2]*xt+r->p[4]*sqrx+r->p[6]*cubx;
        b=r->p[1]+r->p[3]*xt+r->p[5]*sqrx; if(b==0.0) break;
        f=a/b;
      }
      *y=f;
      break;
    case 1312:
      xt=x-r->p[4]; if(xt<=0.0) {
        f=0.0;
      } else {
        f=0.0;
        a=(r->p[0]*(xt)-r->p[2])*exp(r->p[1]*(xt)); f+=a;
        a=r->p[2]*exp(r->p[3]*(xt)); f+=a;
      }
      *y=f;
      break;
    case 1313:
      xt=x-r->p[6]; if(xt<=0.0) {
        f=0.0;
      } else {
        f=0.0;
        a=(r->p[0]*(xt)-r->p[2]-r->p[4])*exp(r->p[1]*(xt)); f+=a;
        a=r->p[2]*exp(r->p[3]*(xt)); f+=a;
        a=r->p[4]*exp(r->p[5]*(xt)); f+=a;
      }
      *y=f;
      break;
    case 1314:
      xt=x-r->p[8]; if(xt<=0.0) {
        f=0.0;
      } else {
        f=0.0;
        a=(r->p[0]*(xt)-r->p[2]-r->p[4]-r->p[6])*exp(r->p[1]*(xt)); f+=a;
        a=r->p[2]*exp(r->p[3]*(xt)); f+=a;
        a=r->p[4]*exp(r->p[5]*(xt)); f+=a;
        a=r->p[6]*exp(r->p[7]*(xt)); f+=a;
      }
      *y=f;
      break;
    case 1401:
      xt=x-r->p[3]; if(xt<=0.0 || r->p[2]==0.0) {
        f=0.0;
      } else {
        f=r->p[0]*pow(xt, r->p[1])*exp(-xt/r->p[2]);
      }
      *y=f;
      break;
    case 1402:
      xt=x-r->p[3]; if(xt<=0.0 || r->p[2]==0.0) {
        f=r->p[4];
      } else {
        f=r->p[0]*pow(xt, r->p[1])*exp(-xt/r->p[2]) + r->p[4];
      }
      *y=f;
      break;
    case 1403:
      xt=x-r->p[0];
      f=0.0;
      if(xt>0.0) {
        a=exp(-xt/r->p[3]);
        if(r->p[1]>0.0) f += r->p[1]*pow(xt, r->p[2])*a;
        if(r->parNr==6 && r->p[4]>0.0) f += r->p[4]*(1.0-a)*exp(-xt/r->p[5]);
      }
      *y=f;
      break;
    case 1421:
      xt=x-r->p[0]; if(xt<=0.0) {
        f=0.0;
      } else {
        f=r->p[1]*(1.0-exp(-pow(xt/r->p[2], r->p[3])));
      }
      *y=f;
      break;
    case 1423:
      xt=x-r->p[0]; if(xt<=0.0) {
        f=0.0;
      } else {
        a=xt/r->p[2]; b=pow(a, r->p[3]-1.0); f=exp(-b*a);
        a=r->p[3]*b*f/r->p[2]; b=1.0-f;
        f=r->p[1]*(a+r->p[4]*b);
      }
      *y=f;
      break;
    case 1431:
      if(x<=0.0) {
        f=0.0;
      } else {
        f=r->p[0]*x*exp(-r->p[1]*x)*r->p[1]*r->p[1];
      }
      *y=f;
      break;
    case 1432:
      if(x<=0.0) {
        f=0.0;
      } else {
        f=r->p[0]*x*exp(-r->p[1]*x);
      }
      *y=f;
      break;
    case 1433:
      if(x<=0.0) {
        f=0.0;
      } else {
        double e=exp(-r->p[1]*x);
        f=x*e + (r->p[2]/(r->p[1]*r->p[1]))*(1.0-(r->p[1]*x+1.0)*e);
        f*=r->p[0];
      }
      *y=f;
      break;
    case 1434:
      if(x<=0.0) {
        f=1.0/(1.0-r->p[0]);
      } else {
        double e=exp(-r->p[2]*x);
        double rcp=x*e + (r->p[3]/(r->p[2]*r->p[2]))*(1.0-(r->p[2]*x+1.0)*e);
        rcp*=r->p[1];
        f=1.0/(1.0-r->p[0]*(1.0-rcp));
      }
      *y=f;
      break;
    case 1801:
      xt=x-r->p[0]; if(xt<=0.0) {
        f=0.0;
      } else {
        f=r->p[1]*pow(xt, r->p[3])/(pow(r->p[2], r->p[3])+pow(xt, r->p[3]));
      }
      *y=f;
      break;
    case 1811:
      xt=x-r->p[0]; if(xt<=0.0) {
        f=0.0;
      } else {
        a=pow(r->p[2], r->p[3])+pow(xt, r->p[3]);
        f=r->p[1]*r->p[3]*
          (pow(xt, r->p[3]-1.0)/a - pow(xt, 2.0*r->p[3]-1.0)/(a*a));
      }
      *y=f;
      break;
    case 1821:
      xt=x-r->p[0]; if(xt<=0.0) {
        f=0.0;
      } else {
        double a, c, b, k, xb, cb, cbxb, hill, hilld;
        a=r->p[1]; c=r->p[2]; b=r->p[3]; k=r->p[4];
        cb=pow(c, b); xb=pow(xt, b); cbxb=cb+xb;
        hilld= b*pow(xt, b-1.0)/cbxb - b*pow(xt, 2.0*b-1.0)/(cbxb*cbxb);
        hill=k*xb/cbxb;
        f=a*(hilld+hill);
      }
      *y=f;
      break;
    case 2801:
      {
        double top, bottom, ec50, hillslope;
        top=r->p[0]; bottom=r->p[1]; ec50=r->p[2]; hillslope=r->p[3]; 
        if(x<=0.0) f=bottom;
        else f=bottom+(top-bottom)/(1.0+pow(ec50/x,hillslope));
      }
      *y=f;
      break;
    case 2802:
      {
        double top, bottom, logec50, hillslope;
        top=r->p[0]; bottom=r->p[1]; logec50=r->p[2]; hillslope=r->p[3]; 
        f=bottom+(top-bottom)/(1.0+pow(10.0, (logec50-x)*hillslope));
      }
      *y=f;
      break;
    case 3331: *y=nan(""); break; /* cannot be applied to one point only */
    case 9501: *y=nan(""); break; /* cannot be applied to one point only */
    case 9502: *y=nan(""); break; /* cannot be applied to one point only */
    case 9503: *y=nan(""); break; /* cannot be applied to one point only */
    case 9701: *y=nan(""); break; /* cannot be applied to one point only */
    default:
      *y=nan("");
  }
  if(isnan(*y)) return 1;
  return 0;
}
/*****************************************************************************/

/*****************************************************************************/
/** Evaluates an array y[i]=f(x[i]).
    @sa fitEval, fitIntegralEvaltac, fitDerivEvaltac
    @return Returns 0, if ok.
 */
int fitEvaltac(
  /** Fit parameters of a single region */
  FitVOI *r,
  /** Times where to evaluate the function */
  double *x,
  /** Array for the function values */
  double *y,
  /** Nr of (x,y) data */
  int dataNr
) {
  if(r==NULL || x==NULL || y==NULL || dataNr<1) return(1);

  /* Special cases that can only be computed as TACs */
  if(r->type==3331) {
    double Ta=r->p[0]; 
    double Ti=r->p[1];
    double dT=r->p[r->parNr-2]; Ta+=dT;
    double tau=r->p[r->parNr-1];
    int n=(r->parNr-4)/2; if(n<1) return(2);
    for(int i=0; i<dataNr; i++) {
      if(x[i]<=Ta) y[i]=0.0;
      else if(x[i]<Ta+Ti) {
        double f=0.0;
        for(int j=0; j<n; j++) {
          double b; 
          if(r->p[2*j+3]>1.0E-12) b=r->p[2*j+2]/r->p[2*j+3]; else b=r->p[2*j+2];
          f+=b*(1.0-exp(-r->p[2*j+3]*(x[i]-Ta)));
        }
        if(Ti>0.0) f*=(1.0/Ti);
        y[i]=f;
      } else {
        double f=0.0;
        for(int j=0; j<n; j++) {
          double b;
          if(r->p[2*j+3]>1.0E-12) b=r->p[2*j+2]/r->p[2*j+3]; else b=r->p[2*j+2];
          f+=b*(exp(-r->p[2*j+3]*(x[i]-Ta-Ti)) - exp(-r->p[2*j+3]*(x[i]-Ta)));
        }
        if(Ti>0.0) f*=(1.0/Ti);
        y[i]=f;
      }
    }
    if(simDispersion(x, y, dataNr, tau, 0.0, NULL)) return(2);
    return(0);
  }

  /* Usual functions */
  for(int i=0; i<dataNr; i++) if(fitEval(r, x[i], y+i)) return(2);
  return 0;
}
/*****************************************************************************/

/*****************************************************************************/
/** Evaluates yi=Integral of f(x) between 0 and x.
    @sa fitEval, fitDerivEval, fitIntegralEvaltac
    @return Returns 0, if ok.
 */
int fitIntegralEval(
  /** Fit parameters of a single region */
  FitVOI *r,
  /** Time where to evaluate integral of the function */
  double x,
  /** The integral value of the function is returned here */
  double *yi
) {
  double a, f, xt, t;

  if(r==NULL || yi==NULL) return 1;
  *yi=nan("");
  switch(r->type) {
    case 301:
      if(fabs(r->p[1])>1.0e-12) f=(r->p[0]/r->p[1])*(exp(r->p[1]*x)-1.0);
      else f=r->p[0]*x;
      *yi=f;
      break;
    case 302:
      f=0;
      if(fabs(r->p[1])>1.0e-12) a=(r->p[0]/r->p[1])*(exp(r->p[1]*x)-1.0);
      else a=r->p[0]*x;
      f+=a;
      if(fabs(r->p[3])>1.0e-12) a=(r->p[2]/r->p[3])*(exp(r->p[3]*x)-1.0);
      else a=r->p[2]*x;
      f+=a;
      *yi=f;
      break;
    case 303:
      f=0;
      if(fabs(r->p[1])>1.0e-12) a=(r->p[0]/r->p[1])*(exp(r->p[1]*x)-1.0);
      else a=r->p[0]*x;
      f+=a;
      if(fabs(r->p[3])>1.0e-12) a=(r->p[2]/r->p[3])*(exp(r->p[3]*x)-1.0);
      else a=r->p[2]*x;
      f+=a;
      if(fabs(r->p[5])>1.0e-12) a=(r->p[4]/r->p[5])*(exp(r->p[5]*x)-1.0);
      else a=r->p[4]*x;
      f+=a;
      *yi=f;
      break;
    case 304:
      f=0;
      if(fabs(r->p[1])>1.0e-12) a=(r->p[0]/r->p[1])*(exp(r->p[1]*x)-1.0);
      else a=r->p[0]*x;
      f+=a;
      if(fabs(r->p[3])>1.0e-12) a=(r->p[2]/r->p[3])*(exp(r->p[3]*x)-1.0);
      else a=r->p[2]*x;
      f+=a;
      if(fabs(r->p[5])>1.0e-12) a=(r->p[4]/r->p[5])*(exp(r->p[5]*x)-1.0);
      else a=r->p[4]*x;
      f+=a;
      if(fabs(r->p[7])>1.0e-12) a=(r->p[6]/r->p[7])*(exp(r->p[7]*x)-1.0);
      else a=r->p[6]*x;
      f+=a;
      *yi=f;
      break;
    case 305:
      f=0;
      if(fabs(r->p[1])>1.0e-12) a=(r->p[0]/r->p[1])*(exp(r->p[1]*x)-1.0);
      else a=r->p[0]*x;
      f+=a;
      if(fabs(r->p[3])>1.0e-12) a=(r->p[2]/r->p[3])*(exp(r->p[3]*x)-1.0);
      else a=r->p[2]*x;
      f+=a;
      if(fabs(r->p[5])>1.0e-12) a=(r->p[4]/r->p[5])*(exp(r->p[5]*x)-1.0);
      else a=r->p[4]*x;
      f+=a;
      if(fabs(r->p[7])>1.0e-12) a=(r->p[6]/r->p[7])*(exp(r->p[7]*x)-1.0);
      else a=r->p[6]*x;
      f+=a;
      if(fabs(r->p[9])>1.0e-12) a=(r->p[8]/r->p[9])*(exp(r->p[8]*x)-1.0);
      else a=r->p[8]*x;
      f+=a;
      *yi=f;
      break;
    case 321:
      f=(r->p[0]/r->p[1])*exp(r->p[1]*x)
        - r->p[0]*exp((r->p[1]+r->p[2])*x)/(r->p[1]+r->p[2]);
      *yi=f;
      break;
    case 322:
      f=0.0;
      a=(r->p[0]/r->p[1])*exp(r->p[1]*x)
        - r->p[0]*exp((r->p[1]+r->p[2])*x)/(r->p[1]+r->p[2]); f+=a;
      a=(r->p[3]/r->p[4])*exp(r->p[4]*x)
        - r->p[3]*exp((r->p[4]+r->p[5])*x)/(r->p[4]+r->p[5]); f+=a;
      *yi=f;
      break;
    case 323:
      f=0.0;
      a=(r->p[0]/r->p[1])*exp(r->p[1]*x)
        - r->p[0]*exp((r->p[1]+r->p[2])*x)/(r->p[1]+r->p[2]); f+=a;
      a=(r->p[3]/r->p[4])*exp(r->p[4]*x)
        - r->p[3]*exp((r->p[4]+r->p[5])*x)/(r->p[4]+r->p[5]); f+=a;
      a=(r->p[6]/r->p[7])*exp(r->p[7]*x)
        - r->p[6]*exp((r->p[7]+r->p[8])*x)/(r->p[7]+r->p[8]); f+=a;
      *yi=f;
      break;
    case 331: *yi=nan(""); break; /* not yet applied */
    case 351: *yi=nan(""); break; /* not yet applied */
    case 1312:
      t=r->p[4]; xt=x-t; f=0.0;
      if(xt>0.0) {
        double A1, A2, L1, L2;
        A1=r->p[0]; L1=r->p[1]; A2=r->p[2]; L2=r->p[3];
        if(L1!=0.0) { 
          a=exp(L1*xt);
          f+=(A2/L1)*(1.0 - a);
          f+=(A1/(L1*L1))*(1.0 + a*(L1*xt - 1.0)); 
        }
        if(L2!=0.0) f+=(A2/L2)*(exp(L2*xt) - 1.0); else f+=A2*xt; 
      }
      *yi=f;
      break;
    case 1313:
      t=r->p[6]; xt=x-t; f=0.0;
      if(xt>0.0) {
        double A1, A2, A3, L1, L2, L3;
        A1=r->p[0]; L1=r->p[1]; A2=r->p[2]; L2=r->p[3]; A3=r->p[4]; L3=r->p[5];
        if(L1!=0.0) { 
          a=exp(L1*xt);
          f+=(A2/L1)*(1.0 - a);
          f+=(A3/L1)*(1.0 - a);
          f+=(A1/(L1*L1))*(1.0 + a*(L1*xt - 1.0)); 
        }
        if(L2!=0.0) f+=(A2/L2)*(exp(L2*xt) - 1.0); else f+=A2*xt; 
        if(L3!=0.0) f+=(A3/L3)*(exp(L3*xt) - 1.0); else f+=A3*xt;
      }
      *yi=f;
      break;
    case 1314:
      t=r->p[8]; xt=x-t; f=0.0;
      if(xt>0.0) {
        double A1, A2, A3, A4, L1, L2, L3, L4;
        A1=r->p[0]; L1=r->p[1]; A2=r->p[2]; L2=r->p[3]; A3=r->p[4]; L3=r->p[5];
        A4=r->p[6]; L4=r->p[7];
        if(L1!=0.0) { 
          a=exp(L1*xt);
          f+=(A2/L1)*(1.0 - a);
          f+=(A3/L1)*(1.0 - a);
          f+=(A4/L1)*(1.0 - a);
          f+=(A1/(L1*L1))*(1.0 + a*(L1*xt - 1.0)); 
        }
        if(L2!=0.0) f+=(A2/L2)*(exp(L2*xt) - 1.0); else f+=A2*xt; 
        if(L3!=0.0) f+=(A3/L3)*(exp(L3*xt) - 1.0); else f+=A3*xt; 
        if(L4!=0.0) f+=(A4/L4)*(exp(L4*xt) - 1.0); else f+=A4*xt;
      }
      *yi=f;
      break;
    case 1401: *yi=nan(""); break; /* not yet applied */
    case 1402: *yi=nan(""); break; /* not yet applied */
    case 1431:
      f=0.0;
      if(x>0.0) f=r->p[0]*(1.0 - (r->p[1]*x + 1.0)*exp(-r->p[1]*x));
      *yi=f;
      break;
    case 1432:
      f=0.0;
      if(x>0.0) {
        a=r->p[0]/(r->p[1]*r->p[1]);
        f=a*(1.0 - (r->p[1]*x + 1.0)*exp(-r->p[1]*x));
      }
      *yi=f;
      break;
    case 2111: *yi=nan(""); break; /* no application */
    case 9501: *yi=nan(""); break; /* cannot be applied to one point only */
    case 9502: *yi=nan(""); break; /* cannot be applied to one point only */
    case 9503: *yi=nan(""); break; /* cannot be applied to one point only */
    case 9701: *yi=nan(""); break; /* cannot be applied to one point only */
    default:
      *yi=nan("");
  }
  if(isnan(*yi)) return 3;
  return 0;
}
/*****************************************************************************/

/*****************************************************************************/
/** Evaluate an array yi[i]=Integral of f(x[i]) between 0 and x.
    @sa fitEvaltac, fitIntegralEvaltac, fitDerivEvaltac
    @return Returns 0, if ok.
*/
int fitIntegralEvaltac(
  /** Fit parameters of a single region */
  FitVOI *r,
  /** Times where to evaluate the function integrals */
  double *x,
  /** Array for the function integral values */
  double *yi,
  /** Nr of (x,yi) data */
  int dataNr
) {
  int i;

  if(r==NULL || x==NULL || yi==NULL || dataNr<1) return(1);
  for(i=0; i<dataNr; i++) if(fitIntegralEval(r, x[i], yi+i)) return(2);
  return(0);
}
/*****************************************************************************/

/*****************************************************************************/
/** Evaluates yd=Df(x).
    @sa fitEval, fitIntegralEval
    @return Returns 0, if ok.
 */
int fitDerivEval(
  /** Fit parameters of a single region */
  FitVOI *r,
  /** Time where to evaluate the derivative of the function */
  double x,
  /** The derivative of the function is returned here */
  double *yd
) {
  double a, f, xt, t;

  if(r==NULL || yd==NULL) return 1;
  *yd=nan("");
  switch(r->type) {
    case 301:
      break;
    case 302:
      break;
    case 303:
      break;
    case 304:
      break;
    case 305:
      break;
    case 321:
      f=r->p[0]*r->p[1]*exp(r->p[1]*x)*(1.0-exp(r->p[2]*x)) - r->p[0]*r->p[2]*exp((r->p[1]+r->p[2])*x);
      *yd=f;
      break;
    case 322:
      f=0.0;
      a=r->p[0]*r->p[1]*exp(r->p[1]*x)*(1.0-exp(r->p[2]*x)) - r->p[0]*r->p[2]*exp((r->p[1]+r->p[2])*x); f+=a;
      a=r->p[3]*r->p[4]*exp(r->p[4]*x)*(1.0-exp(r->p[5]*x)) - r->p[3]*r->p[5]*exp((r->p[4]+r->p[5])*x); f+=a;
      *yd=f;
      break;
    case 323:
      f=0.0;
      a=r->p[0]*r->p[1]*exp(r->p[1]*x)*(1.0-exp(r->p[2]*x)) - r->p[0]*r->p[2]*exp((r->p[1]+r->p[2])*x); f+=a;
      a=r->p[3]*r->p[4]*exp(r->p[4]*x)*(1.0-exp(r->p[5]*x)) - r->p[3]*r->p[5]*exp((r->p[4]+r->p[5])*x); f+=a;
      a=r->p[6]*r->p[7]*exp(r->p[7]*x)*(1.0-exp(r->p[8]*x)) - r->p[6]*r->p[8]*exp((r->p[7]+r->p[8])*x); f+=a;
      *yd=f;
      break;
    case 331: break; /* not yet applied */
    case 351: break; /* not yet applied */
    case 1312:
      t=r->p[4]; xt=x-t; f=0.0;
      if(xt>0.0) {
        double A1, A2, L1, L2;
        A1=r->p[0]; L1=r->p[1]; A2=r->p[2]; L2=r->p[3];
        f= A1*exp(L1*xt) + A2*L2*exp(L2*xt) + (A1*xt -A2)*L1*exp(L1*xt);
      }
      *yd=f;
      break;
    case 1313:
      t=r->p[6]; xt=x-t; f=0.0;
      if(xt>0.0) {
        double A1, A2, A3, L1, L2, L3;
        A1=r->p[0]; L1=r->p[1]; A2=r->p[2]; L2=r->p[3]; A3=r->p[4]; L3=r->p[5];
        f= A1*exp(L1*xt) + A2*L2*exp(L2*xt) + A3*L3*exp(L3*xt) +
           (A1*xt -A2 -A3)*L1*exp(L1*xt);
      }
      *yd=f;
      break;
    case 1314: break;
      t=r->p[8]; xt=x-t; f=0.0;
      if(xt>0.0) {
        double A1, A2, A3, A4, L1, L2, L3, L4;
        A1=r->p[0]; L1=r->p[1]; A2=r->p[2]; L2=r->p[3]; A3=r->p[4]; L3=r->p[5];
        A4=r->p[6]; L4=r->p[7];
        f= A1*exp(L1*xt) + A2*L2*exp(L2*xt) + A3*L3*exp(L3*xt) +
           A4*L4*exp(L4*xt) + (A1*xt -A2 -A3 -A4)*L1*exp(L1*xt);
      }
      *yd=f;
      break;
    case 1401:
    case 1402:
      xt=x-r->p[3]; if(xt<=0.0 || r->p[2]==0.0) {
        f=0.0;
      } else {
        f=r->p[0]*pow(xt, r->p[1]-1.0)*exp(-xt/r->p[2])*(r->p[1]-(xt/r->p[2]));
      }
      *yd=f;
      break;
    case 1421:
      xt=x-r->p[0]; if(xt<=0.0) {
        *yd=0.0;
      } else {
        a=pow(xt/r->p[2], r->p[3]-1.0);
        *yd=r->p[1]*r->p[3]*a*exp(-a*xt/r->p[2])/r->p[2];
      }
      break;
    case 1801:
      xt=x-r->p[0]; if(xt<=0.0) {
        *yd=0.0;
      } else {
        a=pow(r->p[2], r->p[3])+pow(xt, r->p[3]);
        *yd=r->p[1]*r->p[3]*(pow(xt, r->p[3]-1.0)/a - pow(xt, 2.0*r->p[3]-1.0)/(a*a));
      }
      break;
    case 2111: *yd=nan(""); break; /* no application */
    case 9501: *yd=nan(""); break; /* cannot be applied to one point only */
    case 9502: *yd=nan(""); break; /* cannot be applied to one point only */
    case 9503: *yd=nan(""); break; /* cannot be applied to one point only */
    case 9701: *yd=nan(""); break; /* cannot be applied to one point only */
    default:
      *yd=nan("");
  }
  if(isnan(*yd)) return(3);
  return(0);
}
/*****************************************************************************/

/*****************************************************************************/
/** Evaluates an array yd[i]=Df(x[i]).
    @sa fitEvaltac, fitIntegralEvaltac, fitDerivEvaltac
    @return Returns 0, if ok.
 */
int fitDerivEvaltac(
  /** Fit parameters of a single region */
  FitVOI *r,
  /** Times where to evaluate the function derivatives */
  double *x,
  /** Array for the function derivatives */
  double *yd,
  /** Nr of (x,yd) data */
  int dataNr
) {
  int i;

  if(r==NULL || x==NULL || yd==NULL || dataNr<1) return(1);
  for(i=0; i<dataNr; i++) if(fitDerivEval(r, x[i], yd+i)) return(2);
  return(0);
}
/*****************************************************************************/

/*****************************************************************************/
