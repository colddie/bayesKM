/// @file readfile.c
/// @author Vesa Oikonen
/// @brief Functions for reading ASCII data files.
///
/*****************************************************************************/
#include "libtpcmisc.h"
/*****************************************************************************/

/*****************************************************************************/
/** Initiate STR_TOKEN_LIST structure.
 *  This should be called once before first use.
 */
void str_token_list_init(
  /** Pointer to list to be initiated. */
  STR_TOKEN_LIST *lst
) {
  memset(lst, 0, sizeof(STR_TOKEN_LIST));
  lst->list_size=0; lst->token_nr=0; lst->tok=NULL;
}
/*****************************************************************************/

/*****************************************************************************/
/** Free memory allocated for STR_TOKEN_LIST.
 *  All contents are destroyed.
 */
void str_token_list_empty(
  /** Pointer to list to be emptied. */
  STR_TOKEN_LIST *lst
) {
  for(int i=0; i<lst->list_size; i++)
    if(lst->tok[i]!=NULL) free(lst->tok[i]);
  if(lst->tok!=NULL) free(lst->tok);
  lst->tok=NULL;
  lst->list_size=0; lst->token_nr=0;
}
/*****************************************************************************/

/*****************************************************************************/
/** Put a string in STR_TOKEN_LIST.
 *  @return Returns 0, if successful, and <>0 in case of an error.
 */
int str_token_list_add(
  /** List that has to be initialized beforehand. */
  STR_TOKEN_LIST *lst,
  /** String that is added to list. */
  char *new_item
) {
  int i;
  const int add_nr=10;

  if(lst==NULL || new_item==NULL || strlen(new_item)<1) return(1);
  if(lst->list_size<=lst->token_nr) {
    lst->tok=realloc(lst->tok, sizeof(char*)*(lst->list_size+add_nr));
    if(lst->tok==NULL) return(2);
    for(i=lst->list_size; i<lst->list_size+add_nr; i++) lst->tok[i]=NULL;
    lst->list_size+=add_nr;
  }
  lst->tok[lst->token_nr]=strdup(new_item);
  if(lst->tok[lst->token_nr]==NULL) return(3);
  lst->token_nr++;

  return(0);
}
/*****************************************************************************/

/*****************************************************************************/
/** Remove the specified string item from the STR_TOKEN_LIST.
 *  @return Returns 0, if successful, and <>0 in case of an error.
 */
int str_token_list_del(
  /** List that has to be initialized beforehand. */
  STR_TOKEN_LIST *lst,
  /** Item number to remove (1..item_nr). */
  int item
) {
  //printf("str_token_list_del(list with %d items, %d)\n", lst->token_nr, item);
  if(lst==NULL || item<1 || item>lst->token_nr) return(1);
  if(lst->tok[item-1]!=NULL) free(lst->tok[item-1]);
  for(int i=item; i<lst->token_nr; i++) {
    //printf("  index %d -> %d\n", i, i-1);
    lst->tok[i-1]=lst->tok[i]; lst->tok[i]=NULL;
  }
  lst->token_nr--;
  return(0);
}
/*****************************************************************************/

/*****************************************************************************/
/** Read all string tokens from text file into STR_TOKEN_LIST.
 *  List needs to be initialized. Previous contents are deleted.
 *  @return Returns 0, if successful, and <>0 in case of an error.
 */
int str_token_list_read(
  /** Name of text file to read */
  const char *filename,     
  /** Token list is allocated by this function. */
  STR_TOKEN_LIST *lst 
) {
  int i, ret, nr=0;
  char *allfile, *cptr;
  FILE *fp;

  if(lst==NULL || filename==NULL || strlen(filename)<1) return(1);
  str_token_list_empty(lst);

  /* Open file */
  fp=fopen(filename, "r"); if(fp==NULL) return(2);
  /* Get file size */
  nr=0; while((ret=fgetc(fp))!=EOF) nr++; rewind(fp);
  if(nr<1) {fclose(fp); return(0);}
  /* Allocate memory for file contents */
  allfile=(char*)malloc((nr+1)*sizeof(char));
  if(allfile==NULL) {fclose(fp); return(3);}
  /* Read file contents */
  i=0; while((ret=fgetc(fp))!=EOF && i<nr) allfile[i++]=(char)ret;
  fclose(fp); allfile[i]=(char)0;
  /* and then fill the list */
  cptr=strtok(allfile, " ;,|\t\n\r");
  if(cptr==NULL) {free(allfile); return(4);}
  do {
    if(str_token_list_add(lst, cptr)) {free(allfile); return(10);}
    cptr=strtok(NULL, " ;,|\t\n\r");
  } while(cptr!=NULL);
  free(allfile);

  return(0);
}
/*****************************************************************************/

/*****************************************************************************/
/** Read all lines from text file into STR_TOKEN_LIST.
 *  List needs to be initialized. Previous contents are deleted.
 *  @return Returns 0, if successful, and <>0 in case of an error.
 *  @sa str_token_list_init, str_token_list_empty
 */
int textfileReadLines(
  /** Name of text file to read */
  const char *filename,
  /** Token list is allocated by this function. */
  STR_TOKEN_LIST *lst
) {
  int i, ret, nr=0;
  char *allfile, *cptr, *line;
  FILE *fp;

  if(lst==NULL || filename==NULL || strlen(filename)<1) return(1);
  str_token_list_empty(lst);

  /* Open file */
  fp=fopen(filename, "r"); if(fp==NULL) return(2);
  /* Get file size */
  nr=0; while((ret=fgetc(fp))!=EOF) nr++; rewind(fp);
  if(nr<1) {fclose(fp); return(0);} //printf("nr=%d\n", nr);
  /* Allocate memory for file contents */
  allfile=(char*)malloc((nr+1)*sizeof(char));
  if(allfile==NULL) {fclose(fp); return(3);}
  /* Read file contents */
  i=0; while((ret=fgetc(fp))!=EOF && i<nr) allfile[i++]=(char)ret;
  fclose(fp); allfile[i]=(char)0;
  /* and then fill the list */
  int npos=0; cptr=allfile;
  line=strTokenDup(cptr, "\n\r\0", &npos);
  if(line==NULL) {free(allfile); return(4);}
  do {
    //printf("line='%s'\n", line);
    if(str_token_list_add(lst, line)) {free(allfile); free(line); return(10);}
    free(line); if(npos==0) break; cptr=cptr+npos;
    line=strTokenDup(cptr, "\n\r\0", &npos);
  } while(line!=NULL);
  free(allfile);

  return(0);
}
/*****************************************************************************/

/*****************************************************************************/
/** Read list of string tokens from specified file
    Remember to free the memory of string list.
    @return Returns the nr of tokens or <0 in case of an error.
 */
int readStrtokens(
  /** Name of file to read */
  const char *filename,
  /** Pointer to list of strings read and allocated here */ 
  char ***toklist
) {
  int i, ret, nr=0;
  char *allfile, *cptr;
  static char **list;
  FILE *fp;

  /* Check arguments */
  if(strlen(filename)<1) return(-1);
  /* Open file */
  fp=fopen(filename, "r"); if(fp==NULL) return(-2);
  /* Get file size */
  nr=0; while((ret=fgetc(fp))!=EOF) nr++; rewind(fp);
  if(nr<1) {fclose(fp); return(0);}
  /* Allocate memory for file contents */
  allfile=(char*)malloc((nr+1)*sizeof(char));
  if(allfile==NULL) {fclose(fp); return(-4);}
  /* Read file contents */
  i=0; while((ret=fgetc(fp))!=EOF && i<nr) allfile[i++]=(char)ret;
  fclose(fp); allfile[i]=(char)0;
  /* Get the number of tokens */
  cptr=allfile; nr=0;
  while(1) {
    i=strspn(cptr, " ;,|\t\n\r"); cptr+=i;
    i=strcspn(cptr, " ;,|\t\n\r"); cptr+=i;
    if(i==0 || *cptr==0) break; else nr++;
  }
  /*printf("Nr of tokens in %s: %d\n", filename, nr);*/
  /* Allocate memory for array of strings */
  list=(char**)malloc(nr*sizeof(char*));
  if(list==NULL) {free(allfile); return(-5);}
  /* and then fill the list */
  cptr=strtok(allfile, " ;,|\t\n\r");
  for(i=0; i<nr; i++) {
    if(cptr==NULL) {
      for(--i; i>=0; i--) free(list[i]);
      free(list); free(allfile);
      return(-8);
    }
    list[i]=(char*)malloc( (strlen(cptr)+1)*sizeof(char) );
    if(list[i]==NULL) {
      for(--i; i>=0; i--) free(list[i]);
      free(list); free(allfile);
      return(-9);
    }
    strcpy(list[i], cptr);
    cptr=strtok(NULL, " ;,|\t\n\r");
  }
  free(allfile);
  *toklist=list;
  return(nr);
}
/*****************************************************************************/

/*****************************************************************************/
/** Check if ASCII text line starts with comment character '#'.
 *  Comment character is searched from the first non-space character (space
 *  characters here include spaces and tabs).
 *  @return 1 if this is comment line and 0 if not.
 *  @author Vesa Oikonen
 */
int asciiCommentLine(
  /** Pointer to string containing one line of ASCII text file */ 
  const char *line,
  /** Optional pointer which is set to the index of line where
   *  the first non-space character after the comment character starts.
   *  If line does not start with comment character, then this will point to
   *  the first non-space character of the line.      
   *  Enter NULL if not needed. */     
  int *cont
) {
  if(cont!=NULL) *cont=0; 
  if(line==NULL) return 0;
  char *cptr=(char*)line;
  int i=strspn(cptr, " \t"); cptr+=i; if(cont!=NULL) *cont=i;
  if(*cptr!='#') return 0;
  if(cont==NULL) return 1;
  cptr++; i=strspn(cptr, " \t"); *cont+=(i+1);
  return 1;
}
/*****************************************************************************/

/*****************************************************************************/
