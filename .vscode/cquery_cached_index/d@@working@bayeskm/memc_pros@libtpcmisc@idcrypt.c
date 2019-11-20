/// @file idcrypt.c
/// @author Calle Laakkonen
/// @brief Encryption/decryption of subject names and other identification 
///        information in string form.
/// @note  Method is simple and not safe for data transfer over internet,
///        but can be used to hide identification in blinded studies.
///
/*****************************************************************************/
#include "libtpcmisc.h"
/*****************************************************************************/

/*****************************************************************************/

/** Return idcrypt module version info */
const char *libpet_idcrypt_version(void) {
    return "2004-12-14";
}

/** Scramble characters in ASCII range 32-126 using the Vigenere Cipher.
    Other characters are discarded.
\return Returns 0 if successful
    */
int id_crypt(
  /** Original string to be encrypted/decrypted */
  const char *string,
  /** Keyword string */
  const char *key,
  /** Encrypted/decrypted string */
  char *out,
  /** Set to 1 when decrypting, or to 0 when encrypting */
  int decrypt
) {
    char *keystr;
    unsigned int len, r;

    len=strlen(string);
    if(len==0) return 0;

    keystr = malloc(len);
    if(!keystr) return 1;

    if(len>strlen(key)) {
        for(r=0;r<len;r++) {
            keystr[r]=key[r%strlen(key)]-32;
            if(keystr[r]>94) keystr[r]=94;
        }
    } else {
        for(r=0;r<len;r++) {
            keystr[r]=key[r]-32;
            if(keystr[r]>94) keystr[r]=94;
        }
    }

    for(r=0;r<len;r++) {
        int c=(unsigned char)string[r]-32;
        if(c>94) c=94;
        if(decrypt) {
            c=c-keystr[r];
            if(c<0) c=(95-(-c)%95);
        } else {
            c=(c+keystr[r])%95;
        }
        out[r]=c+32;
    }
    out[r]=0;
    free(keystr);
    return 0;
}

