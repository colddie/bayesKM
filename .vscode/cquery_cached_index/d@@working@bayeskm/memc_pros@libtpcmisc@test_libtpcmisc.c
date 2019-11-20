/******************************************************************************
 * This file is not compiled into the library, but it contains main()
 * which is compiled to an executable, used to test the library functions. 
 *****************************************************************************/

/*****************************************************************************/
//#include "tpcclibConfig.h"
#include "libtpcmisc.h"
/*****************************************************************************/

/*****************************************************************************/
int test_strTokenNr()
{
  printf("\n=====================================\n");
  printf("\n%s\n", __func__);
  printf("\n=====================================\n");

  char line[1024], delims[128];
  int n;
  
  strcpy(delims, " \t\n\r");

  strcpy(line, "one two three four");
  n=strTokenNr(line, delims);
  printf("'%s' -> %d\n", line, n);
  if(n!=4) return(1);
  
  strcpy(line, "  one two three four  ");
  n=strTokenNr(line, delims);
  printf("'%s' -> %d\n", line, n);
  if(n!=4) return(2);
  
  strcpy(line, "    ");
  n=strTokenNr(line, delims);
  printf("'%s' -> %d\n", line, n);
  if(n!=0) return(3);
  
  strcpy(line, "onetwothreefour");
  n=strTokenNr(line, delims);
  printf("'%s' -> %d\n", line, n);
  if(n!=1) return(4);
  
  strcpy(delims, " ,\t\n\r");
  strcpy(line, "one, two, three, four");
  n=strTokenNr(line, delims);
  printf("'%s' -> %d\n", line, n);
  if(n!=4) return(5);
  
  return(0);
}
/******************************************************************************/

/******************************************************************************/
int test_strTokenNCpy()
{
  printf("\n=====================================\n");
  printf("\n%s\n", __func__);
  printf("\n=====================================\n");

  char line[1024], delims[128], tmp[4];
  int n, i;
  
  strcpy(delims, " \t\n\r");

  strcpy(line, "one two three four");
  n=strTokenNr(line, delims);
  printf("'%s' -> %d\n", line, n);
  if(n!=4) return(1);
  for(i=1; i<=n; i++) {
    strTokenNCpy(line, delims, i, tmp, 4);
    printf("%d := '%s'\n", i, tmp);
  }    
  
  strcpy(line, "  one two three four  ");
  n=strTokenNr(line, delims);
  printf("'%s' -> %d\n", line, n);
  if(n!=4) return(2);
  for(i=1; i<=n; i++) {
    strTokenNCpy(line, delims, i, tmp, 4);
    printf("%d := '%s'\n", i, tmp);
  }    
  
  strcpy(line, "    ");
  n=strTokenNr(line, delims);
  printf("'%s' -> %d\n", line, n);
  if(n!=0) return(3);
  
  strcpy(line, "onetwothreefour");
  n=strTokenNr(line, delims);
  printf("'%s' -> %d\n", line, n);
  if(n!=1) return(4);
  for(i=1; i<=n; i++) {
    strTokenNCpy(line, delims, i, tmp, 4);
    printf("%d := '%s'\n", i, tmp);
  }    
  
  strcpy(delims, " ,\t\n\r");
  strcpy(line, "one, two, three, four");
  n=strTokenNr(line, delims);
  printf("'%s' -> %d\n", line, n);
  if(n!=4) return(5);
  for(i=1; i<=n; i++) {
    strTokenNCpy(line, delims, i, tmp, 4);
    printf("%d := '%s'\n", i, tmp);
  }    
  
  return(0);
}
/******************************************************************************/

/******************************************************************************/
int test_strChrCount()
{
  printf("\n=====================================\n");
  printf("\n%s\n", __func__);
  printf("\n=====================================\n");

  char line[1024], delims[128];
  int n;
  
  strcpy(delims, " \t\n\r");
  strcpy(line, "one two three four");
  n=strChrCount(line, delims);
  printf("'%s' -> %d\n", line, n);
  if(n!=3) return(1);
  
  strcpy(delims, "t");
  strcpy(line, "one two three four");
  n=strChrCount(line, delims);
  printf("'%s' -> %d\n", line, n);
  if(n!=2) return(2);
  
  strcpy(delims, " ,;\t\n\r");
  strcpy(line, "");
  n=strChrCount(line, delims);
  printf("'%s' -> %d\n", line, n);
  if(n!=0) return(3);
  
  strcpy(delims, ".");
  strcpy(line, "here . .");
  n=strChrCount(line, delims);
  printf("'%s' -> %d\n", line, n);
  if(n!=2) return(4);

  return(0);
}
/******************************************************************************/

/******************************************************************************/
int test_strReplaceChar()
{
  printf("\n=====================================\n");
  printf("\n%s\n", __func__);
  printf("\n=====================================\n");

  char line[1024];
  
  strcpy(line, "one two three four"); printf("'%s' -> ", line);
  strReplaceChar(line, ' ', '_'); printf("'%s'\n", line);
  if(strcmp(line, "one_two_three_four")) return(1);
  
  strcpy(line, "one\ttwo\tthree\tfour"); printf("'%s' -> ", line);
  strReplaceChar(line, '\t', ' '); printf("'%s'\n", line);
  if(strcmp(line, "one two three four")) return(2);

  return(0);
}
/*****************************************************************************/

/*****************************************************************************/
static char *info[] = {
  "Usage: @P [options]",
  " ",
  "Options:",
  " -stdoptions", // List standard options like --help, -v, etc
  " -t, --test",
  "     Run all tests for library functions.",
  0};
/*****************************************************************************/

/*****************************************************************************/
/** Run unit tests to the library functions
 *  @author Vesa Oikonen
 *  @return 0 if all tests pass, otherwise >0.
 * */
int main(
  /** Nr of arguments */
  int argc,
  /** Pointer to arrays of argument string */
  char *argv[ ]
) {
  int i, help=0, version=0, verbose=1, error=0, test=0;
  int ret;
  char *cptr;

  if(argc==1) {tpcPrintUsage(argv[0], info, stdout); return(0);}
  for(i=1; i<argc; i++) {
    if(tpcProcessStdOptions(argv[i], &help, &version, &verbose)==0) continue;
    cptr=argv[i]; if(*cptr=='-') cptr++; if(*cptr=='-') cptr++;
    if(strncasecmp(cptr, "TEST", 1)==0) {
      test=1; continue;
    } else {
      error++; break;
    }
  }
  if(error>0) {
    fprintf(stderr, "Error: specify --help for usage.\n");
    return(1);
  }
  /* Print help or version? */
  if(help) {tpcPrintUsage(argv[0], info, stdout); return(0);}
  if(version) {tpcPrintBuild(argv[0], stdout); return(0);}

  if(test==0) return(0);

  if(verbose>0) printf("running tests for library functions...\n");
  //TPCSTATUS status; statusInit(&status); status.verbose=verbose;
  //statusSet(&status, __func__, __FILE__, __LINE__, 0);
  i=10;
  /* strext.c */
  i++; if((ret=test_strTokenNr())!=0) {
    fprintf(stderr, "failed (%d).\n", ret); return(i);}
  i++; if((ret=test_strTokenNCpy())!=0) {
    fprintf(stderr, "failed (%d).\n", ret); return(i);}
  i++; if((ret=test_strChrCount())!=0) {
    fprintf(stderr, "failed (%d).\n", ret); return(i);}
  i++; if((ret=test_strReplaceChar())!=0) {
    fprintf(stderr, "failed (%d).\n", ret); return(i);}
  


  if(verbose>0) printf("\nAll tests passed.\n\n");
  return(0);
}
/*****************************************************************************/

/*****************************************************************************/
