#!/usr/bin/env python

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
    mcmc_wrap(cor_data[5])

def mcmc_wrap(data):
    ys = data["NDVI_grass"]
    xs = data["SWC10"]
    mcmc_optim(ys, xs)


def mcmc_optim(ys, xs):
#    k_max = pymc.Uniform('k_max', 0, 1)
#    k_lab = ["slope","shift","drag","resist"]
#    beta = [k_max] + [pymc.Uniform('k_{0}'.format(i), -100, 100) for i in k_lab]

#    k_0 = pymc.Uniform('k_0', 0, 1, value=0.3)
#    k_1 = pymc.Uniform('k_1', -1e3, 1e4, value=10)
#    k_2 = pymc.Uniform('k_2', -1e2, 1e2, value=1)
#    k_3 = pymc.Uniform('k_3', -1e2, 1e2, value=1)
#    k_4 = pymc.Uniform('k_4', -1e2, 1e2, value=2)

    k_0 = pymc.Exponential('k_0', 0.1, value=0.3)
    k_1 = pymc.Uniform('k_1', -1e3, 1e4, value=10)
    k_2 = pymc.Uniform('k_2', -1e2, 1e2, value=0.01)
    k_3 = pymc.Uniform('k_3', -1e2, 1e2, value=1)
    k_4 = pymc.Uniform('k_4', -1e2, 1e2, value=2)

    @pymc.deterministic
    def springModel(ks_0=k_0, ks_1=k_1, ks_2=k_2, ks_3=k_3, ks_4=k_4, xs=xs):
        ks = [ks_0, ks_1, ks_2, ks_3, ks_4]
        mySpring = _sd.spring(ks, xs, mo.sig_mod2)
        return mySpring.calc_dynamics()['x']

    likelihood = pymc.Normal('likelihood', mu=springModel, tau=1, value=ys, observed=True)

    #bayesModel = pymc.Model([likelihood]+[beta])
    bayesModel = pymc.Model([likelihood, k_0, k_1, k_2, k_3, k_4])
    springPost = pymc.MCMC(bayesModel, db='pickle', dbname='../outputs/spring_trace')
    springPost.sample(5e4,2e4,2)
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
    print "Program FAIL: check config file"

