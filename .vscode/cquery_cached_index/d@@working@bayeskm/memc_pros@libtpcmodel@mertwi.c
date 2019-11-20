/** @file mertwi.c
 *  @brief Mersenne Twister MT19937 pseudorandom number generator for TPCCLIB.
 *
 *  For more information on the method and original source codes visit
 *  http://www.math.hiroshima-u.ac.jp/~m-mat/MT/emt.html
 *  
 */
/*****************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <time.h>
#include <sys/time.h>
/*****************************************************************************/
#include "libtpcmodel.h"
/*****************************************************************************/

/*****************************************************************************/
/** Prepare data struct Mersenne Twister MT19937 for usage.
    Do not call any mertwi* functions before calling this function!
*/ 
void mertwiInit(
  MERTWI *mt
) {
  mt->n=TPCCLIB_MERTWI_NN; // Length of mt array
  mt->m=mt->n/2; // N/2
  mt->a=TPCCLIB_MERTWI_A;
  mt->um=UINT64_C(0xFFFFFFFF80000000);
  mt->lm=UINT64_C(0x7FFFFFFF);
  mt->mti=mt->n+1; // means that mt[] is not initialized
}
/*****************************************************************************/

/*****************************************************************************/
/** @brief Make uint32_t seed for pseudorandom number generators. 
 *  @details Uses microseconds from the computer clock and process ID to 
 *  reduce the chance of getting the same seed for simulatenously executing 
 *  program threads and instances.
 *  @return Returns the seed. 
 *  @sa mertwiSeed64
 */ 
uint32_t mertwiSeed32(void)
{
  uint32_t li;
#if defined HAVE_TIMESPEC_GET
  struct timespec ts;
  timespec_get(&ts, TIME_UTC);
  li=((ts.tv_sec % 10000)*523 ^ ts.tv_nsec*10) ^ ((getpid() % 1000)*983);
#elif defined HAVE_CLOCK_GETTIME
  struct timespec ts;
  clock_gettime(CLOCK_REALTIME, &ts);
  li=((ts.tv_sec % 10000)*523 ^ ts.tv_nsec*10) ^ ((getpid() % 1000)*983);
#elif defined HAVE_GETTIMEOFDAY
  struct timeval tv;
  gettimeofday(&tv, 0);
  li=((tv.tv_sec % 10000)*523 ^ tv.tv_usec*13) ^ ((getpid() % 1000)*983);
#else
  li=(unsigned int)time(NULL)+(unsigned int)getpid();
#endif
  li+=(uint32_t)rand();
  return(li);
}
/*****************************************************************************/
/** @brief Make uint64_t seed for pseudorandom number generators. 
 *  @details Uses microseconds from the computer clock and process ID to 
 *  reduce the chance of getting the same seed for simulatenously executing 
 *  program threads and instances.
 *  @return Returns the seed. 
 *  @sa mertwiSeed32
 */ 
uint64_t mertwiSeed64(void)
{
  uint32_t li;
  uint64_t lli;
  lli=li=mertwiSeed32();
  lli=li; lli<<=32; lli+=li;
  return(lli);
}
/*****************************************************************************/

/*****************************************************************************/
/** @brief Initialize the state vector mt[] inside data struct for 
      Mersenne Twister MT19937 pseudorandom number generator
      using given seed.
    @details Call either this or mertwiInitByArray64 before generating 
      random numbers with MT19937 using functions.
    @pre Initialize MERTWI data struct by calling mertwiInit() before this.
    @sa mertwiSeed64, mertwiInitByArray64, mertwiInit
 */
void mertwiInitWithSeed64(
  /** Data struct for Mersenne Twister MT19937 pseudorandom number generator */
  MERTWI *mt,
  /** Seed, for example from mertwiSeed64() */
  uint64_t seed
) {
  mt->mt[0]=seed;
  for(mt->mti=1; mt->mti<mt->n; mt->mti++) 
    mt->mt[mt->mti]=
      (UINT64_C(6364136223846793005)*(mt->mt[mt->mti-1] ^ (mt->mt[mt->mti-1]>>62)) + mt->mti);
}
/*****************************************************************************/
/** @brief Initialize the state vector mt[] inside data struct for 
      Mersenne Twister MT19937 pseudorandom number generator
      using given array.
    @details Call either this or mertwiInitWithSeed64 before generating 
      random numbers with MT19937 using functions.
    @pre Initialize MERTWI data struct by calling mertwiInit() before this.
    @sa mertwiSeed64, mertwiInitWithSeed64, mertwiInit
 */
void mertwiInitByArray64(
  /** Data struct for Mersenne Twister MT19937 pseudorandom number generator */
  MERTWI *mt,
  /** The array for initializing keys */
  uint64_t init_key[],
  /** Length of initialization array init_key[] */
  uint64_t key_length
) {
  unsigned int i, j;
  uint64_t k;
  mertwiInitWithSeed64(mt, UINT64_C(19650218));
  i=1; j=0; if(mt->n>key_length) k=mt->n; else k=key_length;
  for(; k; k--) {
    mt->mt[i] = 
      (mt->mt[i] ^ ((mt->mt[i-1] ^ (mt->mt[i-1]>>62)) * UINT64_C(3935559000370003845)))
      + init_key[j] + j;
    i++; j++;
    if(i>=mt->n) {mt->mt[0]=mt->mt[mt->n-1]; i=1;}
    if(j>=key_length) j=0;
  }
  for(k=mt->n-1; k; k--) {
    mt->mt[i]=
      (mt->mt[i] ^ ((mt->mt[i-1] ^ (mt->mt[i-1] >> 62)) * UINT64_C(2862933555777941757)))
       - i;
    i++;
    if(i>=mt->n) {mt->mt[0]=mt->mt[mt->n-1]; i=1;}
  }
  mt->mt[0] = UINT64_C(1) << 63; /* MSB is 1; assuring non-zero initial array */ 
}
/*****************************************************************************/

/*****************************************************************************/
/** @brief Generate a random number on [0, 2^64-1]-interval
      using Mersenne Twister MT19937.
    @pre Initialize MERTWI data struct by calling mertwiInit() before using
      this function. Preferrably also initialize state vector by calling either 
      mertwiInitWithSeed64() or mertwiInitByArray64; if you do not do that,
      mertwiInitWithSeed64() is called automatically with predetermined seed. 
    @sa mertwiInit, mertwiInitWithSeed64, mertwiInitByArray64
    @return Returns the random number with range from 0 to UINT64_MAX.
 */
uint64_t mertwiRandomInt64(
  /** Data struct for Mersenne Twister MT19937 pseudorandom number generator */
  MERTWI *mt
) {
  unsigned int i;
  uint64_t x;
  static uint64_t mag01[2]={UINT64_C(0), TPCCLIB_MERTWI_A};

  if(mt->mti>=mt->n) { /* generate NN words at one time */

    /* if init_genrand64() has not been called, a default initial seed is used */
    if(mt->mti==mt->n+1) mertwiInitWithSeed64(mt, UINT64_C(5489)); 

    for(i=0; i<mt->n-mt->m; i++) {
      x=(mt->mt[i]&mt->um)|(mt->mt[i+1]&mt->lm);
      mt->mt[i] = mt->mt[i+mt->m] ^ (x>>1) ^ mag01[(int)(x&UINT64_C(1))];
    }
    for(; i<mt->n-1; i++) {
      x = (mt->mt[i]&mt->um)|(mt->mt[i+1]&mt->lm);
      mt->mt[i] = mt->mt[i+(mt->m-mt->n)] ^ (x>>1) ^ mag01[(int)(x&UINT64_C(1))];
    }
    x = (mt->mt[mt->n-1]&mt->um)|(mt->mt[0]&mt->lm);
    mt->mt[mt->n-1] = mt->mt[mt->m-1] ^ (x>>1) ^ mag01[(int)(x&UINT64_C(1))];
    mt->mti=0;
  }
  x = mt->mt[mt->mti];

  x ^= (x>>29) & UINT64_C(0x5555555555555555);
  x ^= (x<<17) & UINT64_C(0x71D67FFFEDA60000);
  x ^= (x<<37) & UINT64_C(0xFFF7EEE000000000);
  x ^= (x>>43);
  mt->mti++;

  return(x);
}
/*****************************************************************************/
/** @brief Generate a random number on [0, 2^63-1]-interval
      using Mersenne Twister MT19937.
    @pre Initialize MERTWI data struct by calling mertwiInit() before using
      this function. Preferrably also initialize state vector by calling either 
      mertwiInitWithSeed64() or mertwiInitByArray64; if you do not do that,
      mertwiInitWithSeed64() is called automatically with predetermined seed. 
    @sa mertwiInit, mertwiInitWithSeed64, mertwiInitByArray64
    @return Returns the random number with range from 0 to INT64_MAX.
 */
int64_t mertwiRandomInt63(
  /** Data struct for Mersenne Twister MT19937 pseudorandom number generator */
  MERTWI *mt
) {
  return(int64_t)(mertwiRandomInt64(mt) >> 1);
}
/*****************************************************************************/
/** @brief Generate a 64-bit double precision floating point pseudorandom 
      number in the range of [0,1] with uniform distribution
      using Mersenne Twister MT19937.

    With uniform distribution, the SD=(up-low)/sqrt(12), and 
    CV=(up-low)/(sqrt(3)*(low+up)).

    @pre Initialize MERTWI data struct by calling mertwiInit() before using
      this function. Preferrably also initialize state vector by calling either 
      mertwiInitWithSeed64() or mertwiInitByArray64; if you do not do that,
      mertwiInitWithSeed64() is called automatically with predetermined seed. 
    @sa mertwiInit, mertwiInitWithSeed64, mertwiInitByArray64
    @return Returns the random double value in the range [0,1].
 */
double mertwiRandomDouble1(
  /** Data struct for Mersenne Twister MT19937 pseudorandom number generator */
  MERTWI *mt
) {
  return(mertwiRandomInt64(mt) >> 11) * (1.0/9007199254740991.0);
}
/*****************************************************************************/
/** @brief Generate a 64-bit double precision floating point pseudorandom 
      number in the range of [0,1) with uniform distribution
      using Mersenne Twister MT19937.
    @pre Initialize MERTWI data struct by calling mertwiInit() before using
      this function. Preferrably also initialize state vector by calling either 
      mertwiInitWithSeed64() or mertwiInitByArray64; if you do not do that,
      mertwiInitWithSeed64() is called automatically with predetermined seed. 
    @sa mertwiInit, mertwiInitWithSeed64, mertwiInitByArray64
    @return Returns the random double value in the range [0,1).
 */
double mertwiRandomDouble2(
  /** Data struct for Mersenne Twister MT19937 pseudorandom number generator */
  MERTWI *mt
) {
  return(mertwiRandomInt64(mt) >> 11) * (1.0/9007199254740992.0);
}
/*****************************************************************************/
/** @brief Generate a 64-bit double precision floating point pseudorandom 
      number in the range of (0,1) with uniform distribution
      using Mersenne Twister MT19937.
    @pre Initialize MERTWI data struct by calling mertwiInit() before using
      this function. Preferrably also initialize state vector by calling either 
      mertwiInitWithSeed64() or mertwiInitByArray64; if you do not do that,
      mertwiInitWithSeed64() is called automatically with predetermined seed. 
    @sa mertwiInit, mertwiInitWithSeed64, mertwiInitByArray64
    @return Returns the random double value in the range (0,1).
 */
double mertwiRandomDouble3(
  /** Data struct for Mersenne Twister MT19937 pseudorandom number generator */
  MERTWI *mt
) {
  return((mertwiRandomInt64(mt) >> 12) + 0.5) * (1.0/4503599627370496.0);
}
/*****************************************************************************/

/*****************************************************************************/
/** @brief Fill the given double array with random numbers with uniform 
      distribution between the specified limits.
    @details Applies Mersenne Twister MT19937 pseudorandom number generator.
      With uniform distribution, the SD=(up-low)/sqrt(12), and 
      CV=(up-low)/(sqrt(3)*(low+up)).
    @pre Initialize MERTWI data struct by calling mertwiInit() before using
      this function. Preferrably also initialize state vector by calling either 
      mertwiInitWithSeed64() or mertwiInitByArray64; if you do not do that,
      mertwiInitWithSeed64() is called automatically with predetermined seed. 
    @sa mertwiInit, mertwiInitWithSeed64, mertwiInitByArray64, 
      mertwiRandomDouble1
    @author Vesa Oikonen
    @return 0 when successful, otherwise <> 0.
 */
int mertwiRandomBetween(
  /** Data struct for Mersenne Twister MT19937 pseudorandom number generator */
  MERTWI *mt,
  /** Nr of values in double array */
  unsigned int nr,
  /** Pointer to pre-allocated double array */
  double *d,
  /** Lower limit for random values */
  double low,
  /** Upper limit for random values */
  double up,
  /** Distribution: 0=even, 1=square-root transformation */
  int type
) {
  unsigned int i;
  double dif, v, stl, stu;

  if(nr<1) return(0);
  if(mt==NULL || d==NULL || type<0 || type>1) return(1);

  dif=up-low; if(dif<0.0) return(2);
  if(dif==0.0) {
    for(i=0; i<nr; i++) d[i]=low;
    return(0);
  }

  if(type==0) {
    for(i=0; i<nr; i++) d[i] = mertwiRandomDouble1(mt)*dif + low;
  } else if(type==1) {
    stl=copysign(sqrt(fabs(low)),low); if(!isnormal(stl)) stl=0.0;
    stu=copysign(sqrt(fabs(up)), up); if(!isnormal(stu)) stu=0.0;
    dif=stu-stl;
    for(i=0; i<nr; i++) {
      v=mertwiRandomDouble1(mt)*dif+stl; d[i]=copysign(v*v, v);
    }
  }

  return(0);
}
/*****************************************************************************/

/*****************************************************************************/
/** @brief Generate pseudo-random number with exponential distribution and
      specified mean.
    @details Applies Mersenne Twister MT19937 pseudorandom number generator.
    @pre Initialize MERTWI data struct by calling mertwiInit() before using
      this function. Preferrably also initialize state vector by calling either 
      mertwiInitWithSeed64() or mertwiInitByArray64; if you do not do that,
      mertwiInitWithSeed64() is called automatically with predetermined seed. 
    @sa mertwiInit, mertwiInitWithSeed64, mertwiInitByArray64, 
      mertwiRandomDouble1
    @return Returns the random double value in the range [0,1] with exponential
      distribution.
 */
double mertwiRandomExponential(
  /** Data struct for Mersenne Twister MT19937 pseudorandom number generator */
  MERTWI *mt,
  /** Mean of the exponential distribution */
  double mean
) {
  double r;
  do {r=mertwiRandomDouble1(mt);} while(r==0.0);
  return(-mean*log(r)); 
}
/*****************************************************************************/

/*****************************************************************************/
/** @brief Generate a 64-bit double precision floating point pseudorandom 
      number in the range of [0,1] with normal (Gaussian) distribution.
      using Mersenne Twister MT19937.

    Applies Mersenne Twister MT19937 pseudorandom number generator, and
    the polar form of Box-Müller transform to produce numbers with Gaussian 
    (normal) distribution which has a zero mean and standard deviation of one.

    Box GEP, Muller ME. A note on the generation of random normal deviates.
    Annals of Mathematical Statistics, Volume 29, Issue 2, 1958, 610-611.
    Available from JSTOR http://www.jstor.org/

    @pre Initialize MERTWI data struct by calling mertwiInit() before using
      this function. Preferrably also initialize state vector by calling either 
      mertwiInitWithSeed64() or mertwiInitByArray64; if you do not do that,
      mertwiInitWithSeed64() is called automatically with predetermined seed. 
    @sa mertwiInit, mertwiInitWithSeed64, mertwiInitByArray64, 
      mertwiRandomDouble1, mertwiRandomBetween
    @return Returns the random double value in the range [0,1] with normal
      distribution.
 */
double mertwiRandomGaussian(
  /** Data struct for Mersenne Twister MT19937 pseudorandom number generator */
  MERTWI *mt
) {
  static int ready=0;
  static double dev;
  double fac, rsq, a, b;

  /* If we don't have deviate already, then we'll have to make one */
  if(!ready) {
    do {
      a = 2.*mertwiRandomDouble1(mt) - 1.0; 
      b = 2.*mertwiRandomDouble1(mt) - 1.0; 
      rsq = a*a + b*b; 
    } while (rsq>1.0 || rsq==0.0);
   
    fac = sqrt(-2.0*log(rsq)/rsq);
    dev=a*fac; ready=1; 
    return(b*fac); 
  } else { /* dev is ready so return it */
    ready=0;
    return(dev);
  }
}
/*****************************************************************************/

/*****************************************************************************/
