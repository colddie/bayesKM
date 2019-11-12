#ifndef MEKINETICRIGID_H_INCLUDED
#define MEKINETICRIGID_H_INCLUDED


#include <armadillo>







extern "C" int patlak_c(unsigned int frameNr, double *t0, double *t1,double *tac, double *ctt,
               double tstart, double tstop, double *output, unsigned int verbose, 
               unsigned int llsq_model,
               unsigned int isweight, double *weights) ;

extern "C" double Func (const arma::vec& vals_inp, arma::vec* grad_out, void* opt_data);

extern "C" int meKineticRigid(int argc, float* argv[]);


#endif