// Test following non-linear optimization algorithm
// with derivative - newton(second), cg, gradient descent, bfgs (quasi-newton), lbfgs
// with approximate gradient  - dlib_bfgs, dlib_lbfgs
// without derivcative - powell, byboqa, differential evolution, particle swarm
// nonconvex - trust region?
// unknown - brydon, SUMT?
//
//
// Test following linear least square optimization algorithm
// (multi)linear regression, pearson.c,nnls.c (non-negative) or lswq.c, levenberg-marqaunt



#include <iostream>
#include <armadillo>
#include </home/tsun/bin/dlib-19.18/include/dlib/optimization.h>
#include </home/tsun/bin/dlib-19.18/include/dlib/global_optimization.h>
// #include "tpcclibConfig.h"
// #include "include/libtpcmisc.h"
// #include "libtpcmodel.h"
#include "tgo.h"
#include "optim.hpp"

typedef  dlib::matrix<double,0,1> column_vec;
typedef  dlib::matrix<double,3,1> parameter_vec;


// for dlib linear least squared fitting
double model (
    const column_vec& input,
    const parameter_vec& params
)
{
    const double p0 = params(0);
    const double p1 = params(1);
    const double p2 = params(2);
    const double i0 = input(0);
    const double i1 = input(1);
    const double temp = p0*i0 + p1*i1 + p2;

    return temp;   ///*temp;
}
double residual (
    const std::pair<column_vec, double>& data,
    const parameter_vec& params
)
{
    return model(data.first, params)*model(data.first, params) - data.second;
}
parameter_vec residual_derivative (
    const std::pair<column_vec, double>& data,
    const parameter_vec& params
)
{
    parameter_vec der;
    const double p0 = params(0);
    const double p1 = params(1);
    const double p2 = params(2);
    const double i0 = data.first(0);
    const double i1 = data.first(1);
    const double temp = p0*i0 + p1*i1 + p2;
    der(0) = i0*2*temp;
    der(1) = i1*2*temp;
    der(2) = 2*temp;

    return der;
}

// for dlib 
double booth_fn0(const column_vec& m)
{
    const double x_1 = m(0); 
    const double x_2 = m(1);
    double obj_val = std::pow(x_1 + 2*x_2 - 7.0,2) + std::pow(2*x_1 + x_2 - 5.0,2);
    return obj_val;
}

// for bobyqa and powell only
double booth_fn(int parNr, double *p, void *fdata)
{
    double x_1 = p[0];
    double x_2 = p[1];
    double obj_val = std::pow(x_1 + 2*x_2 - 7.0,2) + std::pow(2*x_1 + x_2 - 5.0,2);
    return obj_val;

    // gradient fee methods goes here
}

// for first-order derivative and derivative-free methods
double booth_fn1(const arma::vec& vals_inp, arma::vec* grad_out, void* opt_data)
{
    double x_1 = vals_inp(0);
    double x_2 = vals_inp(1);
 
    double obj_val = std::pow(x_1 + 2*x_2 - 7.0,2) + std::pow(2*x_1 + x_2 - 5.0,2);
    // gradient required here after
    if (grad_out) {
        (*grad_out)(0) = 2*(x_1 + 2*x_2 - 7.0) + 2*(2*x_1 + x_2 - 5.0)*2;
        (*grad_out)(1) = 2*(x_1 + 2*x_2 - 7.0)*2 + 2*(2*x_1 + x_2 - 5.0);
    }
    //
    return obj_val;
}

// for newton's method
double booth_fn2(const arma::vec& vals_inp, arma::vec* grad_out, arma::mat* hess_out, void* opt_data)
{
    double x_1 = vals_inp(0);
    double x_2 = vals_inp(1);
 
    double obj_val = std::pow(x_1 + 2*x_2 - 7.0,2) + std::pow(2*x_1 + x_2 - 5.0,2);
    // gradient required here after
    if (grad_out) {
        (*grad_out)(0) = 2*(x_1 + 2*x_2 - 7.0) + 2*(2*x_1 + x_2 - 5.0)*2;
        (*grad_out)(1) = 2*(x_1 + 2*x_2 - 7.0)*2 + 2*(2*x_1 + x_2 - 5.0);
    }
    if (hess_out) {
        arma::mat A(2, 2, arma::fill::randu);
        arma::mat S = arma::diagmat(A);
        S(0,0)  = 2+8;
        S(0,1) = 4+4;
        S(1,1)  = 8+2;
        S(1,0)  = 4+4;
        *hess_out = S;
    }
    //
    return obj_val;
}










int main()
{
    // To check convergence
    // print debug information or not, with iteration number, residual value, fitted params
    // for instance, using optimlib:
    // settings.verbose_print_level = 0: default (no print)
    // settings.verbose_print_level = 1 print iteration number, current optimal value of the objective function and corresponding vertex
    // settings.verbose_print_level = 2 print iteration number, current optimal value of the objective function, objective function values for all vertices, and the full vertex matrix.
    bool verbose = 1;

    printf("*************************************\n");
    printf("BFGS/LBFGS with approximate gradient with dlib...\n");
    printf("*************************************\n");

    // // returns a function that represents the derivative of the function f.  It
    // // is approximated numerically by:(f(x+eps)-f(x-eps))/(2*eps)
    column_vec starting_point = {0,1};
    dlib::find_min_using_approximate_derivatives(dlib::bfgs_search_strategy(),
                                                 dlib::objective_delta_stop_strategy(1e-7).be_verbose(),
                                                 booth_fn0, starting_point, -1);
    std::cout << "approximate bfgs : solution to Booth test: \n" << starting_point << std::endl;

    // You can also use an approximate derivative like so for an constrained case
    // Here we put the same [0.1 0.8] range constraint on each variable, however, you
    // can put different bounds on each variable by passing in column vectors of
    // constraints for the last two arguments rather than scalars.  
    starting_point = {0.1, 0.1}; 
    find_min_box_constrained(dlib::lbfgs_search_strategy(10),  
                             dlib::objective_delta_stop_strategy(1e-9).be_verbose(),  
                             booth_fn0, dlib::derivative(booth_fn0), starting_point, 0.0, 10);
    std::cout << std::endl << "constrained booth solution: \n" << starting_point << "\n" << std::endl;


    printf("*************************************\n");
    printf("trust region...\n");
    printf("*************************************\n");
    // // We can also use find_min_trust_region(), which is also a method which uses
    // // second derivatives.  For some kinds of non-convex function it may be more
    // // reliable than using a newton_search_strategy with find_min().
    // starting_point = {0.8, 1.3};
    // find_min_trust_region(objective_delta_stop_strategy(1e-7),
    //                       rosen_model(), 
    //                       starting_point, 
    //                       10 // initial trust region radius
    // );
    // cout << "rosen solution: \n"<< starting_point << endl;



    printf("*************************************\n");
    printf("test powell/baboya...\n");
    printf("*************************************\n");
    // Next, let's try the BOBYQA algorithm.  This is a technique specially
    // designed to minimize a function in the absence of derivative information.  
    // Generally speaking, it is the method of choice if derivatives are not available
    // and the function you are optimizing is smooth and has only one local optima.  As
    // an example, consider the be_like_target function defined below:
    // double *lowlim /** Lower limits for the parameters * /
    // double *uplim /** Upper limits for the parameters */
    // double (*objf)(int, double*, void*) /** The object function */
    // void *objfData /** Data to objective function; NULL if not needed */
    // int dim /** Dimension = nr of parameters */
    // int neighNr /** Nr of neighbours to investigate; enter large nr relative to samNr (below)
    // *  if only few local minima are expected */
    // double *fmin  /** Function value at global minimum */
    // double *gmin /** Global minimum = parameter estimates */
    // int samNr    /** Nr of points to sample in one iteration; enter larger samNr if nr of
    // *  iterations (below) is small */
    // int tgoNr  /** Nr of TGO iterations; enter 0 to use the default; enter 1 to run TGO
    // *  just once ie TGO instead of iTGO. Large iteration nr is needed if
    // *  samNr (above) would otherwise require too much memory. */
    // int verbose   /** Verbose level; if zero, then nothing is printed into stdout or stderr */

    int tgoNr,iterNr,neighNr,parNr;
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

    bool success_0 = tgo(
      pmin, pmax, booth_fn, NULL, parNr, neighNr,
      &wss, output, tgoNr, iterNr, verbose);

    if (!success_0) {
        std::cout << "powell: Booth test completed successfully." << std::endl;
    } else {
        std::cout << "powell: Booth test completed unsuccessfully." << std::endl;
    }
    
    std::cout << "powell: solution to Booth test: \n"  << output[0] <<" " << output[1] << "\n" << std::endl;
    // arma::cout << "powell: solution to Booth test:\n" << output << arma::endl;
    



    printf("*************************************\n");
    printf("test gradient descent...\n");
    printf("*************************************\n");
    // double err_tol the tolerance value controlling how small should be before 'convergence' is declared.
    // int iter_max the maximum number of iterations/updates before the algorithm exits.
    // bool vals_bound whether the search space is bounded. If true, then
    // arma::vec lower_bounds this defines the lower bounds.
    // arma::vec upper_bounds this defines the upper bounds.
    // int gd_method which updating formula to use; see the details section below.
    // gd_settings_t gd_settings a structure containing options for the available GD algorithms; see the details section below.

    // run Adam-based optim
    arma::vec x = arma::ones(parNr,1) + 1.0; // initial values
    optim::algo_settings_t settings;
    settings.gd_method = 6;
    settings.gd_settings.step_size = 0.1;
    settings.verbose_print_level = verbose;
  
    bool success = optim::gd(x,booth_fn1,nullptr,settings);
    arma::cout << "\nAdam: true values vs estimates:\n" << x << arma::endl;


    // run Newton-based optim
    x = arma::ones(parNr,1) + 1.0; // initial values

    success = optim::newton(x,booth_fn2,nullptr);
    arma::cout << "\nnewton: true values vs estimates:\n" << x << arma::endl;
 



    printf("*************************************\n");
    printf("test cg...\n");
    printf("*************************************\n");
    // double err_tol the tolerance value controlling how small should be before 'convergence' is declared.
    // int iter_max the maximum number of iterations/updates before the algorithm exits.
    // bool vals_bound whether the search space is bounded. If true, then
    // arma::vec lower_bounds this defines the lower bounds.
    // arma::vec upper_bounds this defines the upper bounds.
    // int cg_method which updating formula to use; see the details section below.
    // double cg_restart_threshold the value ν in the details section below.

    x = arma::zeros(parNr,1) + 2; // initial values (2,2)
    settings.verbose_print_level = verbose;
    success = optim::cg(x,booth_fn1,nullptr, settings);

    if (success) {
        std::cout << "cg: Booth test completed successfully." << std::endl;
    } else {
        std::cout << "cg: Booth test completed unsuccessfully." << std::endl;
    }
    arma::cout << "cg: solution to Booth test:\n" << x << arma::endl;
 


    printf("*************************************\n");
    printf("test bfgs...\n");
    printf("*************************************\n");
    // double err_tol the tolerance value controlling how small should be before 'convergence' is declared.
    // int iter_max the maximum number of iterations/updates before the algorithm exits.
    // bool vals_bound whether the search space is bounded. If true, then
    // arma::vec lower_bounds this defines the lower bounds.
    // arma::vec upper_bounds this defines the upper bounds.

    x = arma::zeros(parNr,1); // initial values (0,0)
    settings.verbose_print_level = verbose;
    success = optim::bfgs(x,booth_fn1,nullptr,settings);
 
    if (success) {
        std::cout << "bfgs: Booth test completed successfully." << std::endl;
    } else {
        std::cout << "bfgs: Booth test completed unsuccessfully." << std::endl;
    }
 
    arma::cout << "bfgs: solution to Booth test:\n" << x << arma::endl;




 

    printf("*************************************\n");
    printf("test lbfgs...\n");
    printf("*************************************\n");
    // Here we repeat the same thing as above but this time using the L-BFGS 
    // algorithm.  L-BFGS is very similar to the BFGS algorithm, however, BFGS 
    // uses O(N^2) memory where N is the size of the starting_point vector.  
    // The L-BFGS algorithm however uses only O(N) memory.  So if you have a 
    // function of a huge number of variables the L-BFGS algorithm is probably 
    // a better choice.
    // double err_tol the tolerance value controlling how small should be before 'convergence' is declared.
    // int iter_max the maximum number of iterations/updates before the algorithm exits.
    // bool vals_bound whether the search space is bounded. If true, then
    // arma::vec lower_bounds this defines the lower bounds.
    // arma::vec upper_bounds this defines the upper bounds.
    // int lbfgs_par_M storage parameter.

    x = arma::zeros(parNr,1); // initial values (0,0)
    settings.verbose_print_level = verbose;
    success = optim::lbfgs(x,booth_fn1,nullptr,settings);
 
    if (success) {
        std::cout<< "lbfgs: Booth test completed successfully." << std::endl;
    } else {
        std::cout << "lbfgs: Booth test completed unsuccessfully." << std::endl;
    }
    arma::cout << "lbfgs: solution to Booth test:\n" << x << arma::endl;



    printf("*************************************\n");
    printf("test Simplex...\n");
    printf("*************************************\n");
    // double err_tol the value controlling how small should be before 'convergence' is declared.
    // int iter_max the maximum number of iterations/updates before the algorithm exits.
    // bool vals_bound whether the search space is bounded. If true, then
    // arma::vec lower_bounds this defines the lower bounds.
    // arma::vec upper_bounds this defines the upper bounds.
    // double nm_par_alpha reflection parameter.
    // double nm_par_gamma expansion parameter.
    // double nm_par_beta contraction parameter.
    // double nm_par_delta shrikage parameter.
  
    settings.verbose_print_level = verbose;
    x = arma::ones(parNr,1) + 1.0; // initial values: (2,2)
    success = optim::nm(x,booth_fn1,NULL,settings);
 
    if (success) {
        std::cout << "nm: Booth test completed successfully." << std::endl;
    } else {
        std::cout << "nm: Booth test completed unsuccessfully." << std::endl;
    }
    arma::cout << "nm: solution to Booth test:\n" << x << arma::endl;



    printf("*************************************\n");
    printf("test particle swarm optimization...\n");
    printf("*************************************\n");
    // bool vals_bound whether the search space is bounded. If true, then
    // arma::vec lower_bounds this defines the lower bounds.
    // arma::vec upper_bounds this defines the upper bounds.
    // int pso_n_pop population size of each generation.
    // int pso_n_gen number of vectors to generate.
    // int pso_check_freq number of generations between convergence checks.
    // int pso_inertia_method method of inertia decay:
    // pso_inertia_method = 1 linear decay between pso_par_w_max and pso_par_w_min.
    // pso_inertia_method = 2 dampening using pso_par_w_damp parameter.
    // int pso_velocity_method method of updating the velocity weights:
    // pso_velocity_method = 1 fixed weights pso_par_c_cog and pso_par_c_soc.
    // pso_velocity_method = 2 linear decay between pso_par_initial_c_cog and pso_par_final_c_cog, and initial_c_soc and pso_par_final_c_soc.
    // double pso_par_c_cog weight value on the 'cognitive' factor.
    // double pso_par_c_soc weight value on the 'social' factor.
    // arma::vec de_initial_lb lower bounds on the initial population; defaults to init_out_vals −0.5
    // arma::vec de_initial_ub upper bounds on the initial population; defaults to init_out_vals +0.5

    // settings_1.pso_center_particle = false;
    // settings_1.vals_bound = true;
    // settings_1.pso_initial_lb = arma::zeros(2,1) - 10.0;
    // settings_1.pso_initial_ub = arma::zeros(2,1) + 10.0;
    // settings_1.pso_n_pop = 5000;
    // settings_1.pso_n_gen = 4000;
 
    // optim::algo_settings_t settings_1;
    settings.verbose_print_level = verbose;
    x = arma::zeros(parNr,1);

    success = optim::pso(x,booth_fn1,nullptr,settings);
 
    if (success) {
        std::cout << "pso: Booth test completed successfully." << std::endl;
    } else {
        std::cout << "pso: Booth test completed unsuccessfully." << std::endl;
    }
    arma::cout << "pso: solution to Booth test:\n" << x << arma::endl;



    printf("*************************************\n");
    printf("test deferiential evolution...\n");
    printf("*************************************\n");
    // bool vals_bound whether the search space is bounded. If true, then
    // arma::vec lower_bounds this defines the lower bounds.
    // arma::vec upper_bounds this defines the upper bounds.
    // int de_n_pop population size of each generation.
    // int de_n_gen number of vectors to generate.
    // int de_check_freq number of generations between convergence checks.
    // int de_mutation_method which method of mutation to apply:
    // de_mutation_method = 1 applies the 'rand' policy.
    // de_mutation_method = 2 applies the 'best' policy.
    // double de_par_F the mutation parameter in the details section below.
    // double de_par_CR the crossover parameter in the details section below.
    // arma::vec de_initial_lb lower bounds on the initial population; defaults to init_out_vals −0.5
    // arma::vec de_initial_ub upper bounds on the initial population; defaults to init_out_vals +0.5

    x = arma::ones(parNr,1) + 1.0; // initial values: (2,2)
    settings.verbose_print_level = verbose;
 
    success = optim::de(x,booth_fn1,nullptr,settings);
 
    if (success) {
        std::cout << "de: Booth test completed successfully." << std::endl;
    } else {
        std::cout << "de: Booth test completed unsuccessfully." << std::endl;
    }
    arma::cout << "de: solution to Booth test:\n" << x << arma::endl;
 

    printf("\n");
    printf("\n");
    printf("\n");

    printf("*************************************\n");
    printf("*************************************\n");
    printf("*************************************\n");
    printf("simulating linear least square fitting...\n");
    printf("*************************************\n");
    printf("*************************************\n");
    printf("*************************************\n");
    // randomly pick a set of parameters to use in this example
    const parameter_vec params = 10*dlib::randm(3,1);
    std::cout << "params: " << dlib::trans(params) << std::endl;

    // Now let's generate a bunch of input/output pairs according to our model.
    std::vector<std::pair<column_vec, double> > data_samples;
    std::vector<double> aux_nnls_output;
    column_vec input;
    for (int i = 0; i < 1000; ++i)
    {
        input = 10*dlib::randm(2,1);
        double tmp = model(input, params);           /// save a copy for nnls which require to know the sign
        const double output = tmp*tmp;

        // save the pair
        data_samples.push_back(std::make_pair(input, output));
        aux_nnls_output.push_back(tmp);
    }


    printf("*************************************\n");
    printf("NLLS with turku's lib...\n");
    printf("*************************************\n");
    /** Algorithm NNLS (Non-negative least-squares)
    Given an m by n matrix A, and an m-vector B, computes an n-vector X,
    that solves the least squares problem
        A * X = B   , subject to X>=0
    Instead of pointers for working space, NULL can be given to let this
    function to allocate and free the required memory.
    @return Function returns 0 if successful, 1, if iteration count exceeded 3*N,
            or 2 in case of invalid problem dimensions or memory allocation error.
    
    For a given problem with 2 inputs, 1 output, 3 params (1000 samples)
       A*X+B*Y+C=F       <=>
       [X Y 1]*[A B C]^T=[F]
       input   params    output
       nnls_a  nnls_x    nnls_b
    */
    int m, n;
    int nnls_m = 1000; int nnls_n = 3; bool isweight = false;
    double *dptr, *nnls_a[3], nnls_b[nnls_m], nnls_x[3], dataw[nnls_m];
    double *nnls_mat=(double*)malloc(((nnls_n+2)*nnls_m)*sizeof(double));
    for(n=0, dptr=nnls_mat; n<nnls_n; n++) {nnls_a[n]=dptr; dptr+=nnls_m;}
    for(m=0; m<nnls_m; m++) {dataw[m]=1.0;}


    /* Calculate the matrix for NNLS */
    /* Fill  A matrix: */
    /* function #1:  */
    for(m=0; m<nnls_m; m++) nnls_a[0][m]=data_samples[m].first(0);
    /* function #2:  */
    for(m=0; m<nnls_m; m++) nnls_a[1][m]=data_samples[m].first(1);
    /* function #3:  */
    for(m=0; m<nnls_m; m++) nnls_a[2][m]=1.0;
    /* Fill  B array:  */
    for(m=0; m<nnls_m; m++) nnls_b[m]=aux_nnls_output[m];   ///data_samples[m].second;

    /* Apply data weights */
    if(isweight) nnlsWght(nnls_n, nnls_m, nnls_a, nnls_b, dataw);
    if(verbose>10) {
      printf("Matrix A                     Array B\n");
      for(int m=0; m<nnls_m; m++) {
        printf("%12.3f %12.3f %12.3f     %12.3f\n",
          nnls_a[0][m], nnls_a[1][m], nnls_a[2][m], nnls_b[m]);
      }
    }

    /* NNLS */
    int ret=nnls(nnls_a, nnls_m, nnls_n, nnls_b, nnls_x, NULL,
              NULL, NULL, NULL);
    if(ret>1) { /* no solution is possible */
      printf("no solution available"); return(ret); // nosolution_nr++; continue;
    }

    // for(n=0; n<nnls_n; n++) output[n]=nnls_x[n];
    std::cout << "nnls: solution to test:\n" <<nnls_x[0]<<" "<<nnls_x[1]<< " "<<nnls_x[2] << std::endl;



    printf("*************************************\n");
    printf("test LM for least square fitting with dlib...\n");
    printf("*************************************\n");
    // Before we do anything, let's make sure that our derivative function defined above matches
    // the approximate derivative computed using central differences (via derivative()).  
    // If this value is big then it means we probably typed the derivative function incorrectly.
    std::cout << "derivative error: " << length(residual_derivative(data_samples[0], params) - 
                                            dlib::derivative(residual)(data_samples[0], params) ) << std::endl;

    // Now let's use the solve_least_squares_lm() routine to figure out what the
    // parameters are based on just the data_samples.
    parameter_vec x_1;
    x_1 = 1;

    std::cout << "Use Levenberg-Marquardt" << std::endl;
    // Use the Levenberg-Marquardt method to determine the parameters which
    // minimize the sum of all squared residuals.
    solve_least_squares_lm(dlib::objective_delta_stop_strategy(1e-7).be_verbose(), 
                            residual,
                            residual_derivative,
                            data_samples,
                            x_1);

    // Now x contains the solution.  If everything worked it will be equal to params.
    std::cout << "inferred parameters: "<< dlib::trans(x_1) << std::endl;
    std::cout << "solution error:      "<< length(x_1 - params) << std::endl;
    std::cout << std::endl;

    x_1 = 1;
    std::cout << "Use Levenberg-Marquardt, approximate derivatives" << std::endl;
    // If we didn't create the residual_derivative function then we could
    // have used this method which numerically approximates the derivatives for you.
    solve_least_squares_lm(dlib::objective_delta_stop_strategy(1e-7).be_verbose(), 
                            residual,
                            dlib::derivative(residual),
                            data_samples,
                            x_1);

    // Now x contains the solution.  If everything worked it will be equal to params.
    std::cout << "inferred parameters: "<< dlib::trans(x_1) << std::endl;
    std::cout << "solution error:      "<< length(x_1 - params) << std::endl;
    std::cout << std::endl;

    x_1 = 1;
    std::cout << "Use Levenberg-Marquardt/quasi-newton hybrid" << std::endl;
    // This version of the solver uses a method which is appropriate for problems
    // where the residuals don't go to zero at the solution.  So in these cases
    // it may provide a better answer.
    solve_least_squares(dlib::objective_delta_stop_strategy(1e-7).be_verbose(), 
                        residual,
                        residual_derivative,
                        data_samples,
                        x_1);

    // Now x contains the solution.  If everything worked it will be equal to params.
    std::cout << "inferred parameters: "<< dlib::trans(x_1) << std::endl;
    std::cout << "solution error:      "<< length(x_1 - params) << std::endl;



    return 0;

}