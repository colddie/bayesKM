

#include <iostream>
#include <armadillo>
#include <dlfcn.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <string.h>
#include <time.h>
// #include "tpcclibConfig.h"
#include "include/libtpcmisc.h"
#include "libtpcmodel.h"
#include "optim.hpp"

 
double sphere_fn(int parNr, double *p, void *fdata)
{
    double x_1 = p[0];
    double x_2 = p[1];
    double obj_val = std::pow(x_1 + 2*x_2 - 7.0,2) + std::pow(2*x_1 + x_2 - 5.0,2);
    return obj_val;
}
 
double booth_fn(const arma::vec& vals_inp, arma::vec* grad_out, void* opt_data)
{
    double x_1 = vals_inp(0);
    double x_2 = vals_inp(1);
 
    double obj_val = std::pow(x_1 + 2*x_2 - 7.0,2) + std::pow(2*x_1 + x_2 - 5.0,2);
    //
    if (grad_out) {
        (*grad_out)(0) = 2*(x_1 + 2*x_2 - 7.0) + 2*(2*x_1 + x_2 - 5.0)*2;
        (*grad_out)(1) = 2*(x_1 + 2*x_2 - 7.0)*2 + 2*(2*x_1 + x_2 - 5.0);
    }
    //
    return obj_val;
}








int main()
{
    int tgoNr,iterNr,neighNr,parNr;
    int verbose = 10;
    double pmin[2], pmax[2];
    bool TGO_LOCAL_INSIDE=0;
    bool TGO_SQUARED_TRANSF=1;
    tgoNr=300; iterNr=0; neighNr=5;
    parNr=2;
    iterNr=0;
    pmin[0]=0.0;    pmax[0]=10.0;   /* K1    */
    pmin[1]=0.0;    pmax[1]=10.0;   /* K1/k2 */
    double wss=0;
    double *output = (double *)malloc(parNr*sizeof(double));
 

    printf("test linear least square fitting...");


    printf("test powell/baboya...");

    bool success_0 = tgo(
      pmin, pmax, sphere_fn, NULL, parNr, neighNr,
      &wss, output, tgoNr, iterNr, verbose-8);

    if (success_0) {
        std::cout << "powell: Booth test completed successfully." << std::endl;
    } else {
        std::cout << "powell: Booth test completed unsuccessfully." << std::endl;
    }

    // arma::cout << "powell: solution to Booth test:\n" << output << arma::endl;
    
    printf("test bfgs...");
 
    arma::vec output_2 = arma::zeros(2,1); // initial values (0,0)
    bool success_2 = optim::bfgs(output_2,booth_fn,nullptr);
 
    if (success_2) {
        std::cout << "bfgs: Booth test completed successfully." << std::endl;
    } else {
        std::cout << "bfgs: Booth test completed unsuccessfully." << std::endl;
    }
 
    arma::cout << "bfgs: solution to Booth test:\n" << output_2 << arma::endl;
 



    printf("test particle swarm optimization...");

    return 0;

}