#!/usr/bin/env python2

import pymc
import datetime, time
# own modules
import data_handling as _dh
import spring_dynamics as _sd
import model_optim_extras as _mo

__author__ = 'Rhys Whitley, Douglas Kelley, Martin De Kauwe'
__email__ = 'rhys.whitley@mq.edu.au'
__created__ = datetime.datetime(2015,1,14)
__modified__ = time.strftime("%c")
__version__ = '1.0'
__status__ = 'prototype'


def main():
    #mp = _mp.model_plotting(fig_path)

    # import data as a list of Pandas dataframes
    raw_data = dh.import_data(dat_path)

    # map data transformations to each dataset imported
    cor_data = [ dh.grass_correct_data(rd) for rd in raw_data ]

    # test with Sturt Plains
    mcmc_wrap(cor_data[1],cor_data[1].Site[0])
    #[ mcmc_wrap(d,d["Site"][0]) for d in cor_data ]

def mcmc_wrap(data, site):
    print("MCMC Optim Model on ===> ", site)
    ys = data["NDVI_norm"]
    xs = data["SWC10"]
    mcmc_optim(ys, xs, site)


def mcmc_optim(ys, xs, site):

#    k_1 = pymc.Uniform('k_1', 1, 1e3, value=100)
#    k_2 = pymc.Uniform('k_2', 0, 1e2, value=10)
#    k_3 = pymc.Uniform('k_3', 0, 1e2, value=0.001)
#    k_4 = pymc.Uniform('k_4', 0, 1e2, value=5)
    b_1 = pymc.Uniform('b_1', -1e3, 1e3, value=1)
    b_2 = pymc.Uniform('b_2', -1e3, 1e3, value=1)
    b_3 = pymc.Uniform('b_3', -1e3, 1e3, value=1)
    b_4 = pymc.Uniform('b_4', -1e3, 1e3, value=1)
    k_1 = pymc.Exponential('k_1', b_1, value=10)
    k_2 = pymc.Exponential('k_2', b_2, value=0.1)
    k_3 = pymc.Exponential('k_3', b_3, value=1)
    k_4 = pymc.Exponential('k_4', b_4, value=2)

    @pymc.deterministic
    def springModel(ks_1=k_1, ks_2=k_2, ks_3=k_3, ks_4=k_4, xs=xs):
        ks = [ks_1, ks_2, ks_3, ks_4]
        mySpring = _sd.spring(ks, xs, mo.sig_mod1, x_init=0.14, v_init=0.01)
        return mySpring.calc_dynamics()['x']

    likelihood = pymc.Normal('likelihood', mu=springModel, tau=1, value=ys, observed=True)

    bayesModel = pymc.Model([likelihood, k_1, k_2, k_3, k_4, b_1, b_2, b_3, b_4])
    springPost = pymc.MCMC(bayesModel, db='pickle', dbname='../outputs/spring_trace_exp_'+site)
    springPost.sample(5e4, 2e4, 1)

    pymc.Matplot.plot(springPost)

    return None

if __name__ == '__main__':

    dat_path = "../data/"
    fig_path = "../figs/"
    out_path = "../outputs/"
    dh = _dh.data_handling()
    mo = _mo.model_optim_extras()

    main()
else:
    print("Program FAIL: check config file")

