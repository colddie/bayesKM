from pystan import stan


def compile():

    # bernoulli model
    model_code = """
        data {
        int<lower=0> N;
        int<lower=0,upper=1> y[N];
        }
        parameters {
        real<lower=0,upper=1> theta;
        }
        model {
        theta ~ beta(0.5, 0.5);  // Jeffreys' prior
        for (n in 1:N)
            y[n] ~ bernoulli(theta);
        }
    """


    from pystan import StanModel
    sm = StanModel(model_code=model_code)

    return sm


def  rwmh_tac_2tpc(sm,tissue_c,plasma_c,plasma_t,weight,initialK, \
lb,ub,rwmh_par_scale,hmc_step_size,rwmh_n_burnin,rwmh_n_draws):

    data = dict(N=10, y=[0, 1, 0, 1, 0, 1, 0, 1, 1, 1])

    fit = sm.sampling(data=data,refresh=-1)
    print(fit, file=open("debugstan.txt", "a"))

    theta=fit.extract()['theta']

    return (theta)

    # # reuse model with new data
    # new_data = dict(N=6, y=[0, 0, 0, 0, 0, 1])
    # fit2 = sm.sampling(data=new_data)
    # print(fit2)

