/// @file swap.c
/// @author Vesa Oikonen
/// @brief Byte swapping for little to big endian (and vice versa) conversion.
///
/*****************************************************************************/
#include "libtpcmisc.h"
/*****************************************************************************/

/*****************************************************************************/
/** Check whether current platform uses little endian byte order.
 *  See H&S Sec. 6.1.2 pp. 163-4.
\return Returns 1, if current platform is little endian, and 0 if not.
 */
int little_endian()
{
  int x=1;
  if(*(char *)&x==1) return(1); else return(0);
}
/*****************************************************************************/

/*****************************************************************************/
/*!
 *  Swaps the specified short int, int, long int, float, or double
 *  from little endian to big endian or vice versa.
 *  Arguments are allowed to overlap.
 *
 * @param from Pointer to a short int, int, long int, float, or double variable
 * @param to Pointer to a short int, int, long int, float, or double variable
 * @param size Size of from and to (byte nr) must be 1, 2, 4 or 8.
 */
void swap(void *from, void *to, int size) {
  unsigned char c;
  unsigned short int s;
  unsigned long l;
  
  switch(size) {
    case 1:
      *(char *)to=*(char *)from;
      break;
    case 2:
      c=*(unsigned char *)from;
      *(unsigned char *)to = *((unsigned char *)from+1);
      *((unsigned char *)to+1) = c; 
      /*swab(from, to, size); // NOT ANSI */
      break;
    case 4:
      s=*(unsigned short *)from;
      *(unsigned short *)to = *((unsigned short *)from+1);
      *((unsigned short *)to+1) = s;
      swap((char*)to, (char*)to, 2);
      swap((char*)((unsigned short *)to+1), (char*)((unsigned short *)to+1), 2);
      break;
    case 8:
      l=*(unsigned long *)from;
      *(unsigned long *)to = *((unsigned long *)from+1);
      *((unsigned long *)to+1) = l;
      swap((char *)to, (char *)to, 4);
      swap((char*)((unsigned long *)to+1), (char*)((unsigned long *)to+1), 4);
      break;
  }
}
/*****************************************************************************/

/*****************************************************************************/
/*!
 *  In-place swab, replaces the non-ANSI function swab(), which may not
 *  work if data is overlapping.
 *
 * @param buf Pointer to memory
 * @param size Size of buf in bytes
 */
void swabip(void *buf, int size) {
  int i;
  unsigned char c;

  for(i=1; i<size; i+=2) {
    c=*((unsigned char *)buf+i);
    *((unsigned char *)buf+i)=*((unsigned char *)buf+(i-1));
    *((unsigned char *)buf+(i-1))=c;
  }
}
/*****************************************************************************/

/*****************************************************************************/
/*!
 *  In-place swab and swaw, switches words and bytes from an array of 4-byte
 *  ints or floats.
 *
 * @param buf Pointer to memory
 * @param size Size of buf in bytes
 */
void swawbip(void *buf, int size) {
  int i;
  unsigned char c, *cptr;

  cptr=(unsigned char*)buf;
  for(i=0; i<size; i+=4, cptr+=4) {
    c=cptr[0]; cptr[0]=cptr[3]; cptr[3]=c;
    c=cptr[1]; cptr[1]=cptr[2]; cptr[2]=c;
  }
}
/*****************************************************************************/

/*****************************************************************************/
/*!
 * In-place swaw, switches words (but not bytes) from an array of 4-byte
 *  ints or floats.
 *
 * @param buf Pointer to memory
 * @param size Size of buf in bytes
 */
void swawip(void *buf, int size) {
  int i;
  unsigned short int s, *sptr;

  sptr=(unsigned short int*)buf;
  for(i=0; i<size; i+=4, sptr+=2) {
    s=sptr[0]; sptr[0]=sptr[1]; sptr[1]=s;
  }
}
/*****************************************************************************/

/*****************************************************************************/
/*!
 * Printfs as bit string the 32-bit variable pointed to by buf.
 * Far from being optimized, thus only for testing and development purposes.
 *
 * @param buf Pointer to memory
 */
void printf32bits(void *buf) {
  unsigned int u, i;
  int j;

  memcpy(&u, buf, 4);
  for(i=32; i>0; i--) {
    j=i-1;if(i<32 && (i%8)==0) printf(" ");
    if(u & (1L<<j)) printf("1"); else printf("0");
  }
  printf("\n");
}
/*****************************************************************************/

/*****************************************************************************/

