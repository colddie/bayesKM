/// @file gaussdev.c
/// @author Vesa Oikonen
/// @brief normally distributed (Gaussian) deviate with zero mean and unit 
///        variance.
///
/*****************************************************************************/
#include "libtpcmodel.h"
#include <sys/time.h>
#include <time.h>
/*****************************************************************************/
#ifndef RAND_MAX
/** Max number for rand() */
#define RAND_MAX 32767
#endif
/// @cond
#define RS_SCALE (1.0 / (1.0 + RAND_MAX))
/// @endcond
/*****************************************************************************/

long int GAUSSDEV_SEED;

/*****************************************************************************/
/** Applies the polar form of Box-Müller transform to produce pseudo-random
    numbers with Gaussian (normal) distribution which has a zero mean and
    standard deviation of one.
    Box GEP, Muller ME. A note on the generation of random normal deviates.
    Annals of Mathematical Statistics, Volume 29, Issue 2, 1958, 610-611.
    Available from JSTOR http://www.jstor.org/
    @return Returns the pseudo-random number.
    @sa gaussdev2, drand
 */
double gaussdev()
{
  static int ready=0, first=1;
  static double dev;
  double fac, rsq, a, b;
  if(first) {first=0; init_gaussdev();}

  /* If we don't have deviate already, then we'll have to make one */
  if(!ready) {
    do {
      //a = 2.*(double)rand()/(double)RAND_MAX - 1.0; 
      //b = 2.*(double)rand()/(double)RAND_MAX - 1.0; 
      a = 2.*drand() - 1.0; 
      b = 2.*drand() - 1.0; 
      rsq = a*a + b*b; 
    } while (rsq>=1.0 || rsq==0.0);
   
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
/** @brief Make and optionally set the seed for rand(), drand, drandRange, 
 *  and drandGaussian(). 
 *  @details Uses microseconds from the computer clock and process ID to 
 *  reduce the chance of getting the same seed for simulatenously executing 
 *  program threads and instances.
 *  @return Returns the seed for srand(). 
 */ 
unsigned int drandSeed(
  /** Also sets seed with srand (1) or not (0) */
  short int seed
) {
  unsigned int li;
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
  li+=(unsigned int)rand();
  if(seed) srand(li);
  //printf("seed := %u\n", li); printf("RAND_MAX := %u\n", RAND_MAX);
  return(li);
}
/*****************************************************************************/

/*****************************************************************************/
/** Initiate random number generator for gaussdev() */ 
void init_gaussdev()
{
  if(GAUSSDEV_SEED<1L) GAUSSDEV_SEED=893165470L;
  srand(GAUSSDEV_SEED);
}
/*****************************************************************************/

/*****************************************************************************/
/** Applies the polar form of Box-Müller transform to produce pseudo-random
    numbers with Gaussian (normal) distribution which has a zero mean and
    standard deviation of one. This function does never set seed, like
    gaussdev() does, therefore set seed for random number generator before
    first calling this routine, for example with srand(time(NULL));
  
    Box GEP, Muller ME. A note on the generation of random normal deviates.
    Annals of Mathematical Statistics, Volume 29, Issue 2, 1958, 610-611.
    Available from JSTOR http://www.jstor.org/
    @return Returns the pseudo-random number.
    @sa gaussdev, drand
 */
double gaussdev2()
{
  static int ready=0;
  static double dev;
  double fac, rsq, a, b;

  /* If we don't have deviate already, then we'll have to make one */
  if(!ready) {
    do {
      a = 2.*drand() - 1.0; 
      b = 2.*drand() - 1.0; 
      rsq = a*a + b*b; 
    } while (rsq>=1.0 || rsq==0.0);
   
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
/** Alternative function to rand() which returns a double precision floating
    point number in the range of [0,1].
    @return Random value in the range [0,1]. 
    @sa rand_range, gaussdev
 */
double drand()
{
  double d, s;
  s=1.0/(1.0+RAND_MAX);
  do {
    d = ( ( s*rand() + rand() )*s + rand() ) * s;
  } while(d>=1.0);
  return d;
}
/*****************************************************************************/

/*****************************************************************************/
/** Fills the double array with random numbers between specified limits.
    Set seed for random number generator before calling this routine, for
    example with srand(time(NULL));
    @return Returns 0 when successful, otherwise <> 0.
 */
int rand_range(
  /** Nr of values in double array */
  int nr,
  /** Pointer to allocated double array */
  double *d,
  /** Lower limit for random values */
  double low,
  /** Upper limit for random values */
  double up,
  /** Distribution: 0=even, 1=square-root transformation */
  int type
) {
  int i;
  double dif, v, stl, stu;

  if(nr<1) return 0;
  if(d==NULL || type<0 || type>1) return 1;

  dif=up-low; if(dif<0.0) return 2;
  if(dif==0.0) {
    for(i=0; i<nr; i++) d[i]=low;
    return 0;
  }

  if(type==0) {
    for(i=0; i<nr; i++) d[i] = drand()*dif + low;
  } else if(type==1) {
    stl=copysign(sqrt(fabs(low)),low); if(!isnormal(stl)) stl=0.0;
    stu=copysign(sqrt(fabs(up)), up); if(!isnormal(stu)) stu=0.0;
    dif=stu-stl;
    for(i=0; i<nr; i++) {v=drand()*dif+stl; d[i]=copysign(v*v, v);}
  }

  return 0;
}
/*****************************************************************************/

/*****************************************************************************/
