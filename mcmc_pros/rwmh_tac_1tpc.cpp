/*################################################################################
  ##
  ##   Copyright (C) 2011-2018 Keith O'Hara
  ##
  ##   This file is part of the MCMC C++ library.
  ##
  ##   Licensed under the Apache License, Version 2.0 (the "License");
  ##   you may not use this file except in compliance with the License.
  ##   You may obtain a copy of the License at
  ##
  ##       http://www.apache.org/licenses/LICENSE-2.0
  ##
  ##   Unless required by applicable law or agreed to in writing, software
  ##   distributed under the License is distributed on an "AS IS" BASIS,
  ##   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
  ##   See the License for the specific language governing permissions and
  ##   limitations under the License.
  ##
  ################################################################################*/
 
/*
 * RWMH with simple normal model
 */

#include "mcmc.hpp"
#include "tpccm.h"
#include "sim1cm.c"
#include <iostream>
#include <armadillo>

struct norm_data {
    double nsample;
    arma::vec tissue_c;
    arma::vec plasma_c;
    arma::vec plasma_t;
    arma::vec weight;

    unsigned int debug;
    double mu_0;
    double sigma_0;
};

double ll_dens_tpc1(const arma::vec& vals_inp, void* ll_data)
{
    // const double mu = vals_inp(0);
    // const double pi = arma::datum::pi;

    // norm_data* dta = reinterpret_cast<norm_data*>(ll_data);
    // const double sigma = dta->sigma;
    // const arma::vec x = dta->x;

    // const int n_vals = x.n_rows;

    // //

    // const double ret = - static_cast<double>(n_vals) * (0.5*std::log(2*pi) + std::log(sigma)) \
    //                    - arma::accu( arma::pow(x - mu,2) / (2*sigma*sigma) );

    // //

    // return ret;
}

double log_pr_dens_tpc1(const arma::vec& vals_inp, void* ll_data)
{
    // norm_data* dta = reinterpret_cast<norm_data*>(ll_data);

    // const double mu_0 = dta->mu_0;
    // const double sigma_0 = dta->sigma_0;
    // const double pi = arma::datum::pi;

    // const double x = vals_inp(0);

    // const double ret = - 0.5*std::log(2*pi) - std::log(sigma_0) - std::pow(x - mu_0,2) / (2*sigma_0*sigma_0);

    // return ret;
}

double simC1_main_rwmh(const arma::vec& vals_inp, void* ll_data)
{
        //true;
    int rett; 
    const char *debugfile  = "debug.txt";  
    norm_data* dta = reinterpret_cast<norm_data*>(ll_data);
    int nsample = (int)(dta->nsample); 
    bool verbose_flag = dta->debug; 

    double *results = (double *)malloc(nsample*sizeof(double));
    double *plasma_t = (double *)malloc(nsample*sizeof(double));
    double *plasma_c = (double *)malloc(nsample*sizeof(double));

    for (int i=0; i < nsample; i++)
    {
        *(plasma_t+i) = (dta->plasma_t)[i];    // start from first image
        *(plasma_c+i) = (dta->plasma_c)[i];    // start from first image
    } 
    // if(verbose_flag) {
    //     FILE *pfile = fopen(debugfile, "a+");
    //     fprintf(pfile, "plasma %f %f %f %f\n", plasma_c[0], plasma_t[0],plasma_c[22],plasma_t[22]);
    //     fclose(pfile);
    // }
    rett=simC1(plasma_t, plasma_c, nsample, vals_inp(0),vals_inp(1), results);
    // rett=simC2(plasma_t, plasma_c, nsample, vals_inp(0),1.4,0.06,0.002, results, NULL, NULL);
        // if(verbose_flag) {
        //     FILE *pfile = fopen(debugfile, "a+");
        //     fprintf(pfile, "plasma %f \n", vals_inp(0));
        //     fclose(pfile);
        // }
    // std::vector<double> result1(results, results + nsample);
    arma::vec result1(nsample);
    for (int i=0; i < nsample; i++)
    {
        result1(i) = *(results+i);    // start from first image
    }

    double ret = 0.0;
	// arma::vec weight(nsample);
	for (int i=0; i < nsample; i++)
    {
    // (dta->weight)[i]	=  (dta->tissue_c)[i]/arma::accu(dta->tissue_c);
    // need to be modified to incorporate prior
      ret -= (dta->weight)[i]*(result1[i]-(dta->tissue_c)[i]) *(result1[i]-(dta->tissue_c)[i])  *1e5;     ///1e5 before there is no such scale and does not work
    }

        if (verbose_flag) {
            FILE *pfile = fopen(debugfile, "a+");

        fprintf(pfile, "ret value %f %f supplied\n", vals_inp(0), ret);
        fprintf(pfile, "weight value %f %f %f %f %f %f\n", (dta->weight)(0), (dta->weight)(1), (dta->weight)(9), (dta->weight)(10), (dta->weight)(15), (dta->weight)(16));
        fclose(pfile);
        }
     return ret;
    //return ll_dens(vals_inp,ll_data) + log_pr_dens(vals_inp,ll_data);
}

double simC1_main_hmc(const arma::vec& vals_inp, arma::vec* grad_out, void* ll_data)
{
        //true;
    int rett; 
    const char *debugfile  = "debug.txt";  
    norm_data* dta = reinterpret_cast<norm_data*>(ll_data);
    int nsample = (int)(dta->nsample); 
    bool verbose_flag = dta->debug; 

    double *results = (double *)malloc(nsample*sizeof(double));
    double *plasma_t = (double *)malloc(nsample*sizeof(double));
    double *plasma_c = (double *)malloc(nsample*sizeof(double));

    for (int i=0; i < nsample; i++)
    {
        *(plasma_t+i) = (dta->plasma_t)[i];    // start from first image
        *(plasma_c+i) = (dta->plasma_c)[i];    // start from first image
    } 
    // if(verbose_flag) {
    //     FILE *pfile = fopen(debugfile, "a+");
    //     fprintf(pfile, "plasma %f %f %f %f\n", plasma_c[0], plasma_t[0],plasma_c[22],plasma_t[22]);
    //     fclose(pfile);
    // }
    rett=simC1(plasma_t, plasma_c, nsample, vals_inp(0),vals_inp(1), results);
    // rett=simC2(plasma_t, plasma_c, nsample, vals_inp(0),1.4,0.06,0.002, results, NULL, NULL);
        // if(verbose_flag) {
        //     FILE *pfile = fopen(debugfile, "a+");
        //     fprintf(pfile, "plasma %f \n", vals_inp(0));
        //     fclose(pfile);
        // }
    // std::vector<double> result1(results, results + nsample);
    arma::vec result1(nsample);
    for (int i=0; i < nsample; i++)
    {
        result1(i) = *(results+i);    // start from first image
    }

    double ret = 0.0;
	// arma::vec weight(nsample);
	for (int i=0; i < nsample; i++)
    {
    // (dta->weight)[i]	=  (dta->tissue_c)[i]/arma::accu(dta->tissue_c);
    // need to be modified to incorporate prior
      ret -= (dta->weight)[i]*(result1[i]-(dta->tissue_c)[i]) *(result1[i]-(dta->tissue_c)[i])  *1e5;     ///1e5 before there is no such scale and does not work
    }

        if (verbose_flag) {
            FILE *pfile = fopen(debugfile, "a+");

        fprintf(pfile, "ret value %f %f supplied\n", vals_inp(0), ret);
        fprintf(pfile, "weight value %f %f %f %f %f %f\n", (dta->weight)(0), (dta->weight)(1), (dta->weight)(9), (dta->weight)(10), (dta->weight)(15), (dta->weight)(16));
        fclose(pfile);
        }
     return ret;
    //return ll_dens(vals_inp,ll_data) + log_pr_dens(vals_inp,ll_data);
}

extern "C" int rwmh_tac_1tpc(int argc, float * argv[])
{
    const char *debugfile  = "debug.txt";  


    /* debug */
    if (argc > 18 || argc <= 1) {
        FILE *pfile = fopen(debugfile, "a+");
        fprintf(pfile, "NIproj3d: 16 arguments required, %d supplied\n", argc);
        fclose(pfile);
        return -1;
    }

    norm_data dta;
    double nsample = *(double*) argv[0];
    dta.nsample = nsample;

    double *x_dta = (double *)argv[1];
    arma::vec tac(nsample);
    for (int i=0;i<dta.nsample;i++) {
       tac[i] = *(x_dta+i);
    }
    dta.tissue_c = tac;

    double *y_dta = (double *)argv[2];
    arma::vec plasma_t(nsample);
    for (int i=0;i<dta.nsample;i++) {
       plasma_t[i] = *(y_dta+i);
    }
    dta.plasma_t = plasma_t;

    double *z_dta = (double *)argv[3];
    arma::vec plasma_c(nsample);
    for (int i=0;i<dta.nsample;i++) {
       plasma_c[i] = *(z_dta+i);
    }
    dta.plasma_c = plasma_c;

    double *w_dta = (double *)argv[4];
    arma::vec weight(nsample);
    for (int i=0;i<dta.nsample;i++) {
       weight[i] = *(w_dta+i);
    }
    dta.weight = weight;
 

    float *output = (float*)argv[5];
    double K1 = *(double*)argv[6];
    double k2 = *(double*)argv[7];
    double par_scale = *(double*)argv[8];
    double step_size  = *(double*)argv[9];

    unsigned int n_burnin  = *(unsigned int*)argv[10]; 
    unsigned int n_draws   = *(unsigned int*)argv[11];  

    double lb0 = *(double*)argv[12];
    double lb1 = *(double*)argv[13];
    double ub0 = *(double*)argv[14];
    double ub1 = *(double*)argv[15];

    unsigned int verbose_flag = *(unsigned int*)argv[16];
    dta.debug = verbose_flag;
    unsigned int   mcmc = *(unsigned int*)argv[17];
    // FILE *pfile = fopen(debugfile, "a+");
    // fprintf(pfile, "plasma %f %f %f\n", dta.plasma_c[0], dta.plasma_t[0],dta.tissue_c[0]);
    // fprintf(pfile, "Input parameters %f %f %f %f %f %f %d %d\n", K1, k2, k3, k4, rwmh_par_scale, rwmh_step_size,rwmh_n_burnin,rwmh_n_draws);
    // fclose(pfile);

    // mcmc setting
    mcmc::algo_settings_t settings;

    arma::vec initial_val(2);
    initial_val(0) = K1;
    initial_val(1) = k2;

    // [0.6, 1.4, 0.06, 0.002]
    // arma::vec lb(1);
    // // lb(0) = -arma::datum::inf;
    // lb(0) = 0.0;

    // arma::vec ub(1);
    // ub(0) = 1.0;
    arma::vec lb(2);
    // lb(0) = -arma::datum::inf;
    lb(0) = lb0;
    lb(1) = lb1;

    arma::vec ub(2);
    ub(0) = ub0;
    ub(1) = ub1;
    // ub(0) = arma::datum::inf;

    settings.vals_bound = true;
    settings.lower_bounds = lb;
    settings.upper_bounds = ub;

    arma::vec grad_out;
    arma::mat draws_out;
    if (mcmc==0) {
         settings.rwmh_par_scale = par_scale;
         settings.rwmh_n_burnin = n_burnin;
         settings.rwmh_n_draws  = n_draws;
        
         mcmc::rwmh(initial_val,draws_out,simC1_main_rwmh,&dta,settings);
     }
    if (mcmc==1) { 
        //  settings.hmc_par_scale = par_scale;
         settings.hmc_step_size = step_size;
         settings.hmc_n_burnin = n_burnin;
         settings.hmc_n_draws  = n_draws;
         mcmc::hmc(initial_val,draws_out,simC1_main_hmc,&dta,settings);
    }

    // default save format is arma_binary
    if(verbose_flag) {
        draws_out.save("A.bin");
    }
    // force saving in arma_ascii format
    //draws_out.save("A.txt", arma::arma_ascii);

    // if(verbose_flag) {
     FILE *pfile = fopen(debugfile, "a+");
        fprintf(pfile, "rwmh_accept_rate, %f \n", settings.rwmh_accept_rate);
        fclose(pfile);
    // }

    for(int i=0; i <draws_out.size(); i++) {
        *(output+i) = draws_out[i];
    }

    // mat::iterator it     = X.begin();
    // mat::iterator it_end = X.end();

    // for(; it != it_end; ++it)
    //   {
    //     *(out+it) = (*it);
    //   }
    // // save in HDF5 format with internal dataset named as "my_data"
    // draws_out.save(arma::hdf5_name("A.h5", "my_data"));

    // // automatically detect format type while loading
    // mat B;
    // B.load("A.bin");

    // // force loading in arma_ascii format
    // mat C;
    // C.load("A.txt", arma_ascii);
    //arma::cout << "draws:\n" << draws_out.rows(0,9) << arma::endl;

    //std::cout << "acceptance rate = " << settings.rwmh_accept_rate << arma::endl;

    //std::cout << "mcmc mean = " << arma::mean(draws_out) << std::endl;

    return 0;
}
