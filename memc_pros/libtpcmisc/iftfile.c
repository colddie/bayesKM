/// @file iftfile.c
/// @author Vesa Oikonen
/// @brief Reading and writing IFT format files.
///
/******************************************************************************/
#include "libtpcmisc.h"
/************************************************** ***************************/

/*************************************************** **************************/
/** Use space before equal sign in IFT */
int IFT_SPACE_BEFORE_EQ_SIGN = 1;
/** Use space after equal sign in IFT */
int IFT_SPACE_AFTER_EQ_SIGN = 1;
/******************************************************************************/

/******************************************************************************/
/** Read IFT file keys and values. Previous contents of IFT are preserved.

    This function can read the initial ASCII part of files that contain also
    binary data.

   @return Returns 0 if ok. Sets ift->status.
 */
int iftRead(
  /** Pointer to initiated but empty IFT */
  IFT *ift,
  /** Input file name */
  char *filename,
  /** 0=key name is not required, 1=only lines with key and equals sign are read */
  int is_key_required
) {
  int i, ret, nr=0, line=0, eq_type=0, initial_key_nr=0, nonprintable=0;
  char *allfile, *cptr, *key_ptr, *value_ptr, *eq_ptr, *eq_ptr2, *cmt_ptr;
  char empty_char=(char)0;
  FILE *fp;


  /* Check function input */
  if(IFT_TEST) printf("iftRead(*ift, %s)\n", filename);
  if(ift==NULL) return(1);
  if(filename==NULL || strlen(filename)<1) {
    iftSetStatus(ift, IFT_FAULT); return(1);
  }
  if(ift->keyNr>0) initial_key_nr=ift->keyNr;

  /* Open file */
  if(strcasecmp(filename, "stdin")==0) {
    fp=stdin;
  } else {
    fp=fopen(filename, "r");
    if(fp==NULL) {iftSetStatus(ift, IFT_CANNOTREAD); return(2);}
  }

  /* Get file size */
  nr=nonprintable=0; while((ret=fgetc(fp))!=EOF) {
    if(iscntrl(ret) && ret!=13 && ret!=10 && ret!=9) {
      nonprintable=1; break;}
    nr++;
  }
  if(nr<2) {
    if(strcasecmp(filename, "stdin")!=0) fclose(fp);
    if(nonprintable>0) {
      /* File contains non-printable characters; maybe binary file */
      iftSetStatus(ift, IFT_UNKNOWNFORMAT);
    } else {
      /* File just din't have any content */
      iftSetStatus(ift, IFT_NODATA);
    }
    return(3);
  }
  if(IFT_TEST>1) printf("  the size of file is %d bytes\n", nr);
  if(nr>5000000) {
    if(strcasecmp(filename, "stdin")!=0) fclose(fp);
    iftSetStatus(ift, IFT_UNKNOWNFORMAT); return(3);
  }
  rewind(fp);

  /* Allocate memory for file contents */
  allfile=(char*)malloc((nr+1)*sizeof(char));
  if(allfile==NULL) {
    if(strcasecmp(filename, "stdin")!=0) fclose(fp);
    iftSetStatus(ift, IFT_NOMEMORY); return(4);
  }

  /* Read file contents and close the file */
  i=0; while((ret=fgetc(fp))!=EOF && i<nr) allfile[i++]=(char)ret;
  allfile[i]=(char)0;
  if(strcasecmp(filename, "stdin")!=0) fclose(fp);

  /* and then fill the list */
  /* separate the first line */
  cptr=strtok(allfile, "\n\r"); line=0;
  do {
    if(IFT_TEST>2) printf("line %d: '%s'\n", line, cptr);
    /* Remove initial spaces and tabs */
    i=strspn(cptr, " \t"); cptr+=i;
    if(strlen(cptr)<1) {cptr=strtok(NULL, "\n\r"); continue;}
    /* Check if line starts with a comment character */
    if((cmt_ptr=strchr("#!;%", cptr[0]))!=NULL) {
      cmt_ptr=cptr; cptr++; i=strspn(cptr, " \t"); cptr+=i;
      if(strlen(cptr)<1) {cptr=strtok(NULL, "\n\r"); continue;}
    }
    if(IFT_TEST>2) printf("  line %d: '%s'\n", line, cptr);
    /* Find the 'equals' sign */
    eq_ptr=strstr_noquotation(cptr, ":=");
    if(eq_ptr==NULL) eq_ptr=strstr_noquotation(cptr, "=");
    if(eq_ptr==NULL) {
      eq_ptr=strstr_noquotation(cptr, ":");
      /* do not accept time representation */
      if(eq_ptr!=NULL && strlen(cptr)>=strlen(eq_ptr)+2 && istime(eq_ptr-2)<=0) {
        /* ... but search for later equals sign */
        eq_ptr2=strstr_noquotation(eq_ptr+4, ":"); eq_ptr=NULL;
        if(eq_ptr2!=NULL && strlen(cptr)>=strlen(eq_ptr2)+2) {
          if(istime(eq_ptr2-2)<=0) eq_ptr2=NULL; else eq_ptr=eq_ptr2;
        }
      }
    }
    if(eq_ptr==NULL) {
      /* Equals sign not found; if required, then ignore this line */
      if(is_key_required) {cptr=strtok(NULL, "\n\r"); continue;}
      /* If not required, then key="" */
      key_ptr=eq_ptr=&empty_char; value_ptr=cptr;
    } else {
      if(strncmp(eq_ptr, ":=", 2)==0) eq_type=1;
      else if(strncmp(eq_ptr, "=", 1)==0) eq_type=2;
      else if(strncmp(eq_ptr, ":", 1)==0) eq_type=3;
      else eq_type=0;
      *eq_ptr=(char)0; eq_ptr++; key_ptr=cptr;
      /* Find the end of the 'equals' sign; that is the start of value */
      i=strspn(eq_ptr, ":="); value_ptr=eq_ptr+i;
      /* Remove initial spaces and tabs */
      i=strspn(value_ptr, " \t"); value_ptr+=i;
    }
    /* Remove tail spaces and tabs */
    i=strlen(key_ptr); while(i>0 && isspace((int)key_ptr[i-1])) i--; key_ptr[i]=(char)0;
    if(i==0) { /* Length of key name is zero */
      if(is_key_required) {cptr=strtok(NULL, "\n\r"); continue;}
    }
    i=strlen(value_ptr); while(i>0 && isspace((int)value_ptr[i-1])) i--; value_ptr[i]=(char)0;
    if(IFT_TEST>2) printf("  key='%s' value='%s'\n", key_ptr, value_ptr);
    /* Remove quotation marks */
    i=strlen(key_ptr)-1; if(i<0) i=0;
    if((key_ptr[0]=='\'' && key_ptr[i]=='\'') || (key_ptr[0]=='\"' && key_ptr[i]=='\"')) {
      key_ptr[i]=(char)0; if(i>0) key_ptr++;}
    if(strlen(key_ptr)<1) { /* Length of key name without comments is zero */
      if(is_key_required) {cptr=strtok(NULL, "\n\r"); continue;}
    }
    i=strlen(value_ptr)-1; if(i<0) i=0;
    if((value_ptr[0]=='\'' && value_ptr[i]=='\'') || (value_ptr[0]=='\"' && value_ptr[i]=='\"')) {
      value_ptr[i]=(char)0; if(i>0) value_ptr++;}
    if(IFT_TEST>2) printf("    key='%s' value='%s'\n", key_ptr, value_ptr);
    /* Put key and value in the list */
    ret=iftPut(ift, key_ptr, value_ptr, cmt_ptr);
    if(ret) {
      free(allfile); iftEmpty(ift); iftSetStatus(ift, IFT_FAULT);
      return(10+ret);
    }
    /* separate the next line */
    cptr=strtok(NULL, "\n\r"); line++;
  } while(cptr!=NULL);
  free(allfile);
  if(IFT_TEST>2) printf("eq_type=%d\n", eq_type);
  ift->type=eq_type;
  /* Did we actually get any data? */
  if(ift->keyNr<=initial_key_nr) {iftSetStatus(ift, IFT_NODATA); return(7);}

  iftSetStatus(ift, IFT_OK);
  return(0);
}
/*****************************************************************************/

/*****************************************************************************/
/** Read value string from IFT file.
   @return Returns pointer to a locally allocated copy of the value string, or NULL if
    none found; free the returned pointer when no more needed.
*/
char *iftReadValue(
  /** File name */
  char *filename,
  /** String to search for in the key. If NULL, or key containing the string is not found,
      then if there is only one value in the file, pointer to that is returned. */
  char *keystr
) {
  if(filename==NULL || !filename[0]) return(NULL);
  /* Read the file */
  IFT ift; iftInit(&ift); if(iftRead(&ift, filename, 0)!=0) return(NULL);
  /* If file contains just one value, then return pointer to a copy of that */
  if(ift.keyNr==1 && strlen(ift.item[0].value)>0) {
    char *s=strdup(ift.item[0].value);
    iftEmpty(&ift);
    return(s);
  }
  /* If key string to search for was not given, then we cannot do more */
  if(keystr==NULL || strlen(keystr)<1) return(NULL);
  /* Search the list for the key */
  int i=-1;
  for(int li=0; li<ift.keyNr; li++) {
    if(strcasestr(ift.item[li].key, keystr)!=NULL && strlen(ift.item[li].value)>0) {
      i=li; break;
    }
  }
  if(i<0) return(NULL);
  /* If found, then return pointer to a copy of that */
  char *s=strdup(ift.item[i].value);
  iftEmpty(&ift);
  return(s);
}
/******************************************************************************/

/******************************************************************************/
/*!
 * Write one item in IFT to the specified file pointer.
 *
 * @param ift Pointer to initiated but empty IFT
 * @param item Index [0..keyNr-1] of key and value to print
 * @param fp Output file pointer
 * @return 0 if ok.
 */
int iftWriteItem(IFT *ift, int item, FILE *fp) {
  char eq_sign[3];
  int ret;

  if(IFT_TEST) printf("iftWriteItem(*ift, %d, fp)\n", item);
  if(ift==NULL) {return(1);}
  iftSetStatus(ift, IFT_FAULT);
  if(fp==NULL) {return(2);}
  if(item<0 || item>=ift->keyNr) {return(3);}

  iftSetStatus(ift, IFT_OK);
  switch(ift->type) {
    case 1: strcpy(eq_sign, ":="); break;
    case 2: strcpy(eq_sign, "="); break;
    case 3: strcpy(eq_sign, ":"); break;
    case 4: strcpy(eq_sign, " "); break;
    case 5: strcpy(eq_sign, "\t"); break;
    case 6: strcpy(eq_sign, ","); break;
    case 7: strcpy(eq_sign, ";"); break;
    default: strcpy(eq_sign, ":="); break;
  }
  if(ift->item[item].type!=' ' && ift->item[item].type!=(char)0) {
    ret=fprintf(fp, "%c ", ift->item[item].type);
    if(ret<1) {iftSetStatus(ift, IFT_CANNOTWRITE); return(6);}
  }
  if(ift->item[item].key==NULL || strlen(ift->item[item].key)<1) {
    ret=fprintf(fp, "%s\n", ift->item[item].value);
  } else {
    if((IFT_SPACE_BEFORE_EQ_SIGN==0 && IFT_SPACE_AFTER_EQ_SIGN==0) ||
       ift->type==4 || ift->type==5 || ift->type==6 || ift->type==7)
      ret=fprintf(fp, "%s%s%s\n",
        ift->item[item].key, eq_sign, ift->item[item].value);
    else if(IFT_SPACE_BEFORE_EQ_SIGN==1 && IFT_SPACE_AFTER_EQ_SIGN==0)
      ret=fprintf(fp, "%s %s%s\n",
        ift->item[item].key, eq_sign, ift->item[item].value);
    else if(IFT_SPACE_BEFORE_EQ_SIGN==0 && IFT_SPACE_AFTER_EQ_SIGN==1)
      ret=fprintf(fp, "%s%s %s\n",
        ift->item[item].key, eq_sign, ift->item[item].value);
    else
      ret=fprintf(fp, "%s %s %s\n",
        ift->item[item].key, eq_sign, ift->item[item].value);
  }
  if(ret<1) {iftSetStatus(ift, IFT_CANNOTWRITE); return(6);}
  return(0);
}
/******************************************************************************/

/******************************************************************************/
/*!
 * Write all keys and values.
 * 
 * @param ift Pointer to initiated but empty IFT
 * @param filename Output filename; string "stdout" is identified
 * @return 0 if ok.
 */
int iftWrite(IFT *ift, char *filename) {
  int li, ret;
  FILE *fp;

  if(IFT_TEST) printf("iftWrite(*ift, %s)\n", filename);
  if(ift==NULL) return(1);
  if(filename==NULL || strlen(filename)<1) {
    iftSetStatus(ift, IFT_FAULT); return(1);}
  if(ift->keyNr<1) return(0);

  /* Open file */
  if(strcasecmp(filename, "stdout")==0) {
    fp=stdout;
  } else {
    fp=fopen(filename, "w");
    if(fp==NULL) {iftSetStatus(ift, IFT_CANNOTWRITE); return(2);}
  }

  /* Write the contents */
  for(li=0, ret=0; li<ift->keyNr; li++) {
    ret=iftWriteItem(ift, li, fp);
    if(ret) break;
  }
  if(strcasecmp(filename, "stdout")!=0) fclose(fp);
  return(ret);
}
/******************************************************************************/

/******************************************************************************/
/** Read definition file, for example microPET header file, into IFT struct.
\return Returns 0 if ok. Sets ift->status.
 */
int defRead(
  /** Pointer to initiated but empty IFT */
  IFT *ift,
  /** Input filename */
  char *filename
) {
  int i, j, ret, nr=0, line=0, initial_key_nr=0, nonprintable=0;
  char *allfile, *cptr, *key_ptr, *value_ptr, *cmt_ptr;
  FILE *fp;

  /* Check function input */
  if(IFT_TEST) printf("defRead(*ift, %s)\n", filename);
  if(ift==NULL) return(1);
  if(filename==NULL || strlen(filename)<1) {
    iftSetStatus(ift, IFT_FAULT); return(1);
  }
  if(ift->keyNr>0) initial_key_nr=ift->keyNr;

  /* Open file */
  if(strcasecmp(filename, "stdin")==0) {
    fp=stdin;
  } else {
    fp=fopen(filename, "r");
    if(fp==NULL) {iftSetStatus(ift, IFT_CANNOTREAD); return(2);}
  }

  /* Get file size */
  nr=nonprintable=0; while((ret=fgetc(fp))!=EOF) {
    if(iscntrl(ret) && ret!=13 && ret!=10 && ret!=9) {
      nonprintable=1; break;}
    nr++;
  }
  if(nr<2) {
    if(strcasecmp(filename, "stdin")!=0) fclose(fp);
    if(nonprintable>0) {
      /* File contains non-printable characters; maybe binary file */
      iftSetStatus(ift, IFT_UNKNOWNFORMAT);
    } else {
      /* File just din't have any content */
      iftSetStatus(ift, IFT_NODATA);
    }
    return(3);
  }
  if(IFT_TEST>1) printf("  the size of file is %d bytes\n", nr);
  if(nr>5000000) {
    if(strcasecmp(filename, "stdin")!=0) fclose(fp);
    iftSetStatus(ift, IFT_UNKNOWNFORMAT); return(3);
  }
  rewind(fp);

  /* Allocate memory for file contents */
  allfile=(char*)malloc((nr+1)*sizeof(char));
  if(allfile==NULL) {
    if(strcasecmp(filename, "stdin")!=0) fclose(fp);
    iftSetStatus(ift, IFT_NOMEMORY); return(4);
  }

  /* Read file contents and close the file */
  i=0; while((ret=fgetc(fp))!=EOF && i<nr) allfile[i++]=(char)ret;
  allfile[i]=(char)0;
  if(strcasecmp(filename, "stdin")!=0) fclose(fp);

  /* and then fill the list */
  /* separate the first line */
  cptr=strtok(allfile, "\n\r"); line=0;
  do {
    if(IFT_TEST>10) printf("line %d: '%s'\n", line, cptr);
    /* Remove initial spaces and tabs */
    i=strspn(cptr, " \t"); cptr+=i;
    if(strlen(cptr)<1) {cptr=strtok(NULL, "\n\r"); continue;}
    /* Check if line starts with a comment character */
    if((cmt_ptr=strchr("#!;%", cptr[0]))!=NULL) {
      /* save the whole line as a key */
      key_ptr=cptr; value_ptr=NULL; j=IFT_TEST; IFT_TEST=0;
      ret=iftPut(ift, key_ptr, value_ptr, " "); IFT_TEST=j;
      if(ret) {
        free(allfile); iftEmpty(ift); iftSetStatus(ift, IFT_FAULT);
        return(10+ret);
      }
      /* separate the next line */
      cptr=strtok(NULL, "\n\r"); line++;
      continue;
    }
    if(IFT_TEST>11) printf("  line %d: '%s'\n", line, cptr);
    /* Get the first line string representing key name */
    key_ptr=cptr; cptr=strchr(cptr, ' ');
    if(cptr!=NULL) {*cptr=(char)0; cptr++;}
    /* Rest of the line is the key value */
    if(cptr==NULL) value_ptr=NULL; else value_ptr=cptr;
    j=IFT_TEST; IFT_TEST=0;
    ret=iftPut(ift, key_ptr, value_ptr, " "); IFT_TEST=j;
    if(ret) {
      free(allfile); iftEmpty(ift); iftSetStatus(ift, IFT_FAULT);
      return(10+ret);
    }
    /* separate the next line */
    cptr=strtok(NULL, "\n\r"); line++;
  } while(cptr!=NULL);
  free(allfile);

  /* Did we actually get any data? */
  if(ift->keyNr<=initial_key_nr) {iftSetStatus(ift, IFT_NODATA); return(7);}

  /* Set parameter file type */
  ift->type=4;

  iftSetStatus(ift, IFT_OK);
  return(0);
}
/******************************************************************************/

/******************************************************************************/
