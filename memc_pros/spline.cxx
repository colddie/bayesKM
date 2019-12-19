// g++ -c -w -Wall -Werror -fpic spline.cxx -I/home/tsun/bin/spline-master/src/ -o spline.o
//  gcc -shared -o libspline.so spline.o


#include <iostream>
#include <vector>
#include <stdlib.h>
#include <iomanip>
#include "spline.h"  

extern "C"  int spline(int argc, float* argv[])     //meKineticRigid
{
    int   npar         =   *(int *)   argv[0];
    int   nframe       =   *(int *)   argv[1];
    float *xterm       =   (float *)  argv[2];
    float *xterm0      =   (float *)  argv[3];
    float *yterm0      =   (float *)  argv[4];
    float *output      =   (float *)  argv[5];	
  	int   verbose      =  *(int *)    argv[6];

    // if (verbose) { printf("nframe xterm ytrem %d %f %f", nframe, *(xterm+20), *(yterm+20)); } 
    std::vector<double> X(npar); for (int ipar=0;ipar<npar;ipar++) { X[ipar]=*(xterm0+ipar); }
    std::vector<double> Y(npar); for (int ipar=0;ipar<npar;ipar++) { Y[ipar]=*(yterm0+ipar); }
    if (verbose) {  
      printf("debug"); 
      printf("nframe npar %d %d\n", nframe, npar);
      for (int ipar=0;ipar<npar;ipar++) { printf(" %f", X[ipar]); }
    }


    tk::spline s1; s1.set_points(X,Y);
     
    for (int iframe=0;iframe<nframe;iframe++) {
       output[iframe] = s1(xterm[iframe]);
    }

  return 1;
}


