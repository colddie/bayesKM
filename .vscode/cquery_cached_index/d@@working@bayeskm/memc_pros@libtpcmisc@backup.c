/// @file backup.c
/// @author Vesa Oikonen
/// @brief Functions for file copying and making backup.
///
/*****************************************************************************/
#include "libtpcmisc.h"
/*****************************************************************************/

/*****************************************************************************/
/** Check if specified file exists; rename existing file to a backup file.
 *  If also backup file exists, then remove that.
 *  @return Returns 0, if successful, and <>0 in case of an error.
 */
int backupExistingFile(
  /** Name of file which, if it exists, is renamed to a backup file */
  char *filename,
  /** Extension for backup file; NULL will set the default ".bak" extension. */
  char *backup_ext,
  /** Pointer to a string (allocated for at least 64 chars) where error message
      or other execution status will be written; enter NULL, if not needed */     
  char *status
) {
  char bakfile[FILENAME_MAX];
  int ret;

  // Check the input
  if(filename==NULL || strlen(filename)<1) {
    if(status!=NULL) sprintf(status, "invalid filename");
    return 1;
  }

  // Check if file exists; if not then no need to make any backup
  if(access(filename, 0) == -1) {
    if(status!=NULL) sprintf(status, "file does not pre-exist");
    return 0;
  }
  // Creat filename for the backup file
  strlcpy(bakfile, filename, FILENAME_MAX);
  if(backup_ext==NULL) strlcat(bakfile, ".bak", FILENAME_MAX);
  else strlcat(bakfile, backup_ext, FILENAME_MAX);
  // If also backup file exists, then just delete it 
  if(access(bakfile, 0) != -1) {
    ret=remove(bakfile);
    if(ret!=0) {
      if(status!=NULL) sprintf(status, "cannot delete previous backup file");
      return 3;
    }
  }
  // Rename file
  ret=rename(filename, bakfile);
  if(ret!=0) {
    if(status!=NULL) sprintf(status, "cannot rename file as backup");
    return 5;
  }
  if(status!=NULL) sprintf(status, "file renamed as backup");
  return 0;
}  
/*****************************************************************************/

/*****************************************************************************/
/** Copy file contents to another file. Existing file will be overwritten,
 *  to prevent it call backupExistingFile() before calling this function.
\return Returns 0 if successfull, otherwise >0.
 */ 
int fileCopy(
  /** Name of file to be copied */
  char *filename1,
  /** Name of new file */
  char *filename2,
  /** Pointer to a string (allocated for at least 64 chars) where error message
      or other execution status will be written; enter NULL, if not needed */     
  char *status
) {
  FILE *from, *to;
  char c;

  // Check the input
  if(filename1==NULL || filename2==NULL) {
    if(status!=NULL) sprintf(status, "invalid filename");
    return 1;
  }
  // Open the file1 for reading
  if((from=fopen(filename1, "rb"))==NULL) {
    if(status!=NULL) sprintf(status, "cannot open file for read");
    return 2;
  }
  // Open file2 for writing
  if((to=fopen(filename2, "wb"))==NULL) {
    if(status!=NULL) sprintf(status, "cannot open file for write");
    fclose(from); return 3;
  }
  // Copy the file
  while(!feof(from)) {
    c=fgetc(from);
    if(ferror(from)) {
      if(status!=NULL) sprintf(status, "cannot read from file");
      fclose(from); fclose(to); (void)remove(filename2); return 4;
    }
    if(!feof(from)) fputc(c, to);
    if(ferror(to)) {
      if(status!=NULL) sprintf(status, "cannot write to file");
      fclose(from); fclose(to); (void)remove(filename2); return 6;
    }
  }
  // Close files
  if(fclose(from)==EOF) {
    if(status!=NULL) sprintf(status, "cannot close file");
    fclose(to); return 7;
  }
  if(fclose(to)==EOF) {
    if(status!=NULL) sprintf(status, "cannot close file");
    return 8;
  }
  return 0;
}
/*****************************************************************************/

/*****************************************************************************/

