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
#include "sim2cm.c"
#include "simrtcm.c"
// #include "simPatlak.c"
// #include "simLogan.c"
#include <iostream>
#include <armadillo>
#include <dlfcn.h>

struct norm_data {
    double nsample;
    double nparams = 4.0;
    arma::vec tissue_c;
    arma::vec plasma_c;
    arma::vec plasma_t;
    arma::vec weight;
    arma::vec prior;
    arma::vec plasma_t1;

    unsigned int debug;
    unsigned int model;
    unsigned int useprior;
    double mu_0;
    double sigma_0;
    double tstart;
    double tstop;
    double k2;
};

double ll_dens_tpc2(const arma::vec& vals_inp, void* ll_data)
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

double log_pr_dens_tpc2(const arma::vec& vals_inp, void* ll_data)
{
    // norm_data* dta = reinterpret_cast<norm_data*>(ll_data);

    // const double mu_0 = dta->mu_0;
    // const double sigma_0 = dta->sigma_0;
    // const double pi = arma::datum::pi;

    // const double x = vals_inp(0);

    // const double ret = - 0.5*std::log(2*pi) - std::log(sigma_0) - std::pow(x - mu_0,2) / (2*sigma_0*sigma_0);

    // return ret;
}

double simC2_main_rwmh(const arma::vec& vals_inp, void* ll_data)
{
    int rett; 
    const char *debugfile  = "debug.txt";  
    norm_data* dta = reinterpret_cast<norm_data*>(ll_data);
    int nsample = (int)(dta->nsample);
    int nparams = (int)(dta->nparams); 
    bool verbose_flag = dta->debug; 
    int startframe = 0;
    int stopframe = nsample;

    double *results  = (double *)malloc(nsample*sizeof(double));
    double *plasma_t = (double *)malloc(nsample*sizeof(double));
    double *plasma_c = (double *)malloc(nsample*sizeof(double));
    double *prior    = (double *)malloc(nparams*sizeof(double));
  

    for (int i=0; i < nsample; i++)
    {
        *(plasma_t+i) = (dta->plasma_t)[i];    // start from first image
        *(plasma_c+i) = (dta->plasma_c)[i];    // start from first image
    } 

    if (dta->model == 1) {
                  rett=simC1(plasma_t, plasma_c, nsample, vals_inp(0),vals_inp(1), results); }
    if (dta->model == 2) {
                  rett=simC2(plasma_t, plasma_c, nsample, vals_inp(0),vals_inp(1),vals_inp(2),vals_inp(3), results, NULL, NULL); }
    if (dta->model == 3) {   
                  rett=simSRTM(plasma_t, plasma_c, nsample, vals_inp(0),vals_inp(1),vals_inp(2), results); }
    if (dta->model == 4) {   
                  rett=simRTCM(plasma_t, plasma_c, nsample, vals_inp(0),vals_inp(1),vals_inp(2),vals_inp(3), results, NULL, NULL); }
    if (dta->model == 5) { 
                double *plasma_t1 = (double *)malloc(nsample*sizeof(double));
                for (int i=0; i < nsample; i++) {  *(plasma_t1+i) = (dta->plasma_t1)[i]; }
                double tstart = (double)(dta->tstart);
                double tstop  = (double)(dta->tstop);    
  

                void  *imglib = NULL;
                int (*dummy)(unsigned int,double,double,double*,double*,double*,double,double,double*,unsigned int);
                int handle;
                imglib = dlopen("/home/tsun/bin/tpcclib-master/build/bin/libmtga_idl.so", RTLD_LAZY);   

                if ( imglib != NULL ) {
                    *(void **)(&dummy) = dlsym(imglib, "simPatlak");
                    handle=dummy(nsample,vals_inp(0),vals_inp(1),plasma_t,plasma_t1,plasma_c,tstart,tstop,results,0);  
                }
                free(plasma_t1);
                // rett=simPatlak(nsample,vals_inp(0),vals_inp(1),plasma_t,plasma_t1,tstart,tstop,results,0);  
                } 
    if (dta->model == 6) { 
                double *plasma_t1 = (double *)malloc(nsample*sizeof(double));
                for (int i=0; i < nsample; i++) {  *(plasma_t1+i) = (dta->plasma_t1)[i]; }
                double tstart = (double)(dta->tstart);
                double tstop  = (double)(dta->tstop);   
                double k2     = (double)(dta->k2); 

                void  *imglib = NULL;
                int (*dummy)(unsigned int,double,double,double*,double*,double*,double,double,double*,unsigned int,double);
                int handle;
                imglib = dlopen("/home/tsun/bin/tpcclib-master/build/bin/libmtga_idl.so", RTLD_LAZY);
                if ( imglib != NULL ) {
                    *(void **)(&dummy) = dlsym(imglib, "simLogan");
                    handle=dummy(nsample,vals_inp(0),vals_inp(1),plasma_t,plasma_t1,plasma_c,tstart,tstop,results,0,k2);  
                }
                free(plasma_t1);
                // rett=simLogan(nsample,vals_inp(0),vals_inp(1),plasma_t,plasma_t1,tstart,tstop,results,0, k2);
                }     
    
    arma::vec result1(nsample);
    for (int i=0; i < nsample; i++)
    {
        result1(i) = *(results+i);    // start from first image
    }

    if (dta->tstart > 0.1) { 
        for (int i=0; i < nsample; i++) { if (plasma_t[i] > dta->tstart) { startframe = i; break;}}
        for (int i=nsample-1; i >=0; i--) { if (plasma_t[i] < dta->tstop) { stopframe = i+1; break;}}
        }
    // printf("startframe stopframe %d %d",startframe,stopframe);
    double ret = 0.0;
	// arma::vec weight(nsample);
	for (int i=startframe; i < stopframe; i++)
    {  ret -= (dta->weight)[i]*(result1[i]-(dta->tissue_c)[i]) *(result1[i]-(dta->tissue_c)[i]);     ///1e5 before there is no such scale and does not work 
    }  

    if (dta->useprior ==1) {
        for (int j=startframe; j < stopframe; j++)
        {  unsigned int lambda = 1.0; double sigma=1.0; ret -= lambda* (vals_inp[j]-prior[j])*sigma*(vals_inp[j]-prior[j]); 
        }
    }

    ret *= 1; //1e2;   ///////////////

    //  const double ret = - arma::accu(arma::pow(result1-dta->tissue_c,2))  *100000;   /// before there is no such scale and does not work
    if (verbose_flag) {
        FILE *pfile = fopen(debugfile, "a+");
    // fprintf(pfile, "rett value %f %f %f %f %d supplied\n", (dta->weight)(0), (dta->weight)(1),(dta->weight)(2),(dta->weight)(3), nsample);
    fprintf(pfile, "ret value %f %f %f %f %f supplied\n", vals_inp(0), vals_inp(1), vals_inp(2), vals_inp(3), ret);
    fclose(pfile);
    }

    free(results);  
    free(plasma_t); 
    free(plasma_c); 
    free(prior); 

     return ret;
    //return ll_dens(vals_inp,ll_data) + log_pr_dens(vals_inp,ll_data);
}


double simC2_main_hmc(const arma::vec& vals_inp, arma::vec* grad_out, void* ll_data)
{
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

    if (dta->model == 1) {
                  rett=simC1(plasma_t, plasma_c, nsample, vals_inp(0),vals_inp(1), results); }
    if (dta->model == 2) {
                  rett=simC2(plasma_t, plasma_c, nsample, vals_inp(0),vals_inp(1),vals_inp(2),vals_inp(3), results, NULL, NULL); }
    if (dta->model == 3) {   
                  rett=simSRTM(plasma_t, plasma_c, nsample, vals_inp(0),vals_inp(1),vals_inp(2), results); }
    if (dta->model == 4) {   
                  rett=simRTCM(plasma_t, plasma_c, nsample, vals_inp(0),vals_inp(1),vals_inp(2),vals_inp(3), results, NULL, NULL); }
    if (dta->model == 5) { 
                double *plasma_t1 = (double *)malloc(nsample*sizeof(double));
                for (int i=0; i < nsample; i++) {  *(plasma_t1+i) = (dta->plasma_t1)[i]; }
                double tstart = (double)(dta->tstart);
                double tstop  = (double)(dta->tstop);    

                void  *imglib = NULL;
                int (*dummy)(unsigned int,double,double,double*,double*,double*,double,double,double*,unsigned int);
                int handle;
                imglib = dlopen("/home/tsun/bin/tpcclib-master/build/bin/libmtga_idl.so", RTLD_LAZY);
                if ( imglib != NULL ) {
                    *(void **)(&dummy) = dlsym(imglib, "simPatlak");
                    handle=dummy(nsample,vals_inp(0),vals_inp(1),plasma_t,plasma_t1,plasma_c,tstart,tstop,results,0);  
                }
                free(plasma_t1);
                // rett=simPatlak(nsample,vals_inp(0),vals_inp(1),plasma_t,plasma_t1,tstart,tstop,results,0);  
                } 
    if (dta->model == 6) { 
                double *plasma_t1 = (double *)malloc(nsample*sizeof(double));
                for (int i=0; i < nsample; i++) {  *(plasma_t1+i) = (dta->plasma_t1)[i]; }
                double tstart = (double)(dta->tstart);
                double tstop  = (double)(dta->tstop);   
                double k2     = (double)(dta->k2); 

                void  *imglib = NULL;
                int (*dummy)(unsigned int,double,double,double*,double*,double*,double,double,double*,unsigned int,double);
                int handle;
                imglib = dlopen("/home/tsun/bin/tpcclib-master/build/bin/libmtga_idl.so", RTLD_LAZY);
                if ( imglib != NULL ) {
                    *(void **)(&dummy) = dlsym(imglib, "simLogan");
                    handle=dummy(nsample,vals_inp(0),vals_inp(1),plasma_t,plasma_t1,plasma_c,tstart,tstop,results,0,k2);  
                }
                free(plasma_t1);
                // rett=simLogan(nsample,vals_inp(0),vals_inp(1),plasma_t,plasma_t1,tstart,tstop,results,0, k2);
                } 

    arma::vec result1(nsample);
    for (int i=0; i < nsample; i++)
    {
        result1(i) = *(results+i);    // start from first image
    }

    double ret = 0.0;
	// arma::vec weight(nsample);
	for (int i=0; i < nsample; i++)
    {
    // need to be modified to incorporate prior
      ret -= (dta->weight)[i]*(result1[i]-(dta->tissue_c)[i]) *(result1[i]-(dta->tissue_c)[i]);  //  *1e5;     ///1e5 before there is no such scale and does not work
    }
    //  const double ret = - arma::accu(arma::pow(result1-dta->tissue_c,2))  *100000;   /// before there is no such scale and does not work
    if (verbose_flag) {
        FILE *pfile = fopen(debugfile, "a+");
        fprintf(pfile, "ret value %f %f %d supplied\n", vals_inp(0), ret, verbose_flag);
        fclose(pfile);
    }
    free(results);  
    free(plasma_t); 
    free(plasma_c); 
    // if (dta->model == 5 || dta->model == 6) {  free(plasma_t1); }
     return ret;
    //return ll_dens(vals_inp,ll_data) + log_pr_dens(vals_inp,ll_data);
}

extern "C" int rwmh_tac_2tpc(int argc, float * argv[])
{
    const char *debugfile  = "debug.txt";  


    /* debug */
    if (argc > 23 || argc <= 1) {
        FILE *pfile = fopen(debugfile, "a+");
        fprintf(pfile, "NIproj3d: 22 arguments required, %d supplied\n", argc);
        fclose(pfile);
        return -1;
    }

    norm_data dta;
    unsigned int nsample = *(unsigned int*) argv[0];
    dta.nsample = nsample;

    double *x_dta = (double *)argv[1];
    arma::vec tac(nsample);
    for (int i=0;i<nsample;i++) {
       tac[i] = *(x_dta+i);
    }
    dta.tissue_c = tac;

    double *y_dta = (double *)argv[2];
    arma::vec plasma_t(nsample);
    for (int i=0;i<nsample;i++) {
       plasma_t[i] = *(y_dta+i);
    }
    dta.plasma_t = plasma_t;


    double *z_dta = (double *)argv[3];
    arma::vec plasma_c(nsample);
    for (int i=0;i<nsample;i++) {
       plasma_c[i] = *(z_dta+i);
    }
    dta.plasma_c = plasma_c;

    double *w_dta = (double *)argv[4];
    arma::vec weight(nsample);
    for (int i=0;i<nsample;i++) {
       weight[i] = *(w_dta+i);
    }
    dta.weight = weight;

    double *pri = (double *)argv[5];
    arma::vec prior(dta.nparams);
    for (int i=0;i<dta.nparams;i++) {
       prior[i] = *(pri+i);
    }
    dta.prior = prior;

    float *output = (float*)argv[6]; 
 

    // double K1 = *(double*)argv[7];
    // double k2 = *(double*)argv[8];
    // double k3 = *(double*)argv[9];
    // double k4 = *(double*)argv[10];

    double *iv = (double *)argv[7];
    arma::vec initial_val(dta.nparams);
    for (int i=0;i<dta.nparams;i++) {
       initial_val[i] = *(iv+i);
    }

    double *lowb = (double *)argv[8];
    arma::vec lb(dta.nparams);
    for (int i=0;i<dta.nparams;i++) {
       lb[i] = *(lowb+i);
    }

    double *highb = (double *)argv[9];
    arma::vec ub(dta.nparams);
    for (int i=0;i<dta.nparams;i++) {
       ub[i] = *(highb+i);
    }

    double par_scale = *(double*)argv[10];
    double step_size  = *(double*)argv[11];

    unsigned int n_burnin  = *(unsigned int*)argv[12]; 
    unsigned int n_draws   = *(unsigned int*)argv[13];  

    // double lb0 = *(double*)argv[15];
    // double lb1 = *(double*)argv[16];
    // double lb2 = *(double*)argv[17];
    // double lb3 = *(double*)argv[18];

    // double ub0 = *(double*)argv[19];
    // double ub1 = *(double*)argv[20];
    // double ub2 = *(double*)argv[21];
    // double ub3 = *(double*)argv[22];

    unsigned int usemodel = *(unsigned int*)argv[14];
    dta.model = usemodel;

    unsigned int verbose_flag = *(unsigned int*)argv[15];
    dta.debug = verbose_flag;

    unsigned int mcmc = *(unsigned int*)argv[16];

    unsigned int useprior = *(unsigned int*)argv[17];
    dta.useprior = useprior;

        
    if ((double *)argv[18] != NULL) {
        double *y_dta = (double *)argv[18];
        arma::vec plasma_t1(nsample);
        for (int i=0;i<nsample;i++) {
            plasma_t1[i] = *(y_dta+i);
        }
        dta.plasma_t1 = plasma_t1;

        double tstart = *(double*)argv[19]; dta.tstart = tstart;
        double tstop  = *(double*)argv[20]; dta.tstop = tstop; 
    }
    if ((double *)argv[21] != NULL) {    double k2     = *(double*)argv[21]; dta.k2 = k2;  }


    // if (verbose_flag) {
    //     FILE *pfile = fopen(debugfile, "a+");
    //     fprintf(pfile, " %f %f %f %f supplied\n", initial_val(0),initial_val(1),initial_val(2),initial_val(3));
    //     fprintf(pfile, " %f %f %f %f supplied\n", lb(0),lb(1),lb(2),lb(3));
    //     fprintf(pfile, " %f %f %f %f supplied\n", ub(0),ub(1),ub(2),ub(3));
    //     fclose(pfile);
    // }
    // initial estimate
    // arma::vec initial_val(4);
    // initial_val(0) = K1;
    // initial_val(1) = k2;
    // initial_val(2) = k3;
    // initial_val(3) = k4;

    // // arma::vec ub(1);
    // // ub(0) = 1.0;
    // arma::vec lb(4);
    // // lb(0) = -arma::datum::inf;
    // lb(0) = lb0;
    // lb(1) = lb1;
    // lb(2) = lb2;
    // lb(3) = lb3;

    // arma::vec ub(4);
    // ub(0) = ub0;
    // ub(1) = ub1;
    // ub(2) = ub2;
    // ub(3) = ub3;
    // // ub(0) = arma::datum::inf;

    // fix last parameter
    if (dta.model == 3) {};

    // mcmc setting
    mcmc::algo_settings_t settings;

    settings.vals_bound = true;
    settings.lower_bounds = lb;
    settings.upper_bounds = ub;

    arma::vec grad_out;
    arma::mat draws_out;
    if (mcmc==0) {
         settings.rwmh_par_scale = par_scale;
         settings.rwmh_n_burnin = n_burnin;
         settings.rwmh_n_draws  = n_draws;
        
         mcmc::rwmh(initial_val,draws_out,simC2_main_rwmh,&dta,settings);
     }
    if (mcmc==1) { 
        //  settings.hmc_par_scale = par_scale;
         settings.hmc_step_size = step_size;
         settings.hmc_n_burnin = n_burnin;
         settings.hmc_n_draws  = n_draws;
         mcmc::hmc(initial_val,draws_out,simC2_main_hmc,&dta,settings);
    }


    // force saving in arma_ascii format
    draws_out.save("A.txt", arma::arma_ascii);

    // if(verbose_flag) {
     FILE *pfile = fopen(debugfile, "a+");
        fprintf(pfile, "rwmh_accept_rate, %f \n", settings.rwmh_accept_rate);
        fclose(pfile);
    // }

    for(int i=0; i <draws_out.size(); i++) {
        *(output+i) = draws_out[i];
    }


    return 0;
}
