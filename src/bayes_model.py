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
    cor_data = [ dh.normalise_xydata(rd, yvar="NDVI_grass") \
            for rd in raw_data ]

    # test with Sturt Plains
    mcmc_wrap(cor_data[1],cor_data[1].Site[1])
    #[ mcmc_wrap(d,d["Site"][0]) for d in cor_data ]

def mcmc_wrap(data, site):
    print("MCMC Optim Model on ===> {0}".format(site))
    ys = data["NDVI_norm"]
    xs = data["SWC10"]
    mcmc_optim(ys, xs, site)

def get_coefficients(map_):
    coeff = [{str(variable): variable.value} for variable \
            in map_.variables if str(variable).startswith('k')]
    return sorted(coeff)

# Standard error calculation
def get_stderrors(self, par, hess, Nx):
    hessian = np.linalg.inv(hess)
    fmin0 = opt_res['fun']
    par = opt_res['x']
    Mp = len(par)
    H_sol = np.diag(np.linalg.solve(hessian,np.eye(Mp)))
    coeff = Mp*fmin0/(Nx-Mp)
    return np.sqrt(coeff*H_sol)


def mcmc_optim(ys, xs, site, use_MCMC=False):

    k_1 = pymc.Uniform('k_1', 1, 1e3, value=30)
    k_2 = pymc.Uniform('k_2', 0, 1e2, value=4)
    k_3 = pymc.Uniform('k_3', 0, 1e2, value=1)
    k_4 = pymc.Uniform('k_4', 0, 1e2, value=2)
    k_5 = pymc.Uniform('k_5', 0, 1, value=0.5)
    k_6 = pymc.Uniform('k_6', -1, 1, value=0)

#    b_1 = pymc.Uniform('b_1', -1e3, 1e3, value=1)
#    b_2 = pymc.Uniform('b_2', -1e3, 1e3, value=1)
#    b_3 = pymc.Uniform('b_3', -1e3, 1e3, value=1)
#    b_4 = pymc.Uniform('b_4', -1e3, 1e3, value=1)
#    k_1 = pymc.Exponential('k_1', b_1, value=10)
#    k_2 = pymc.Exponential('k_2', b_2, value=0.1)
#    k_3 = pymc.Exponential('k_3', b_3, value=1)
#    k_4 = pymc.Exponential('k_4', b_4, value=2)

    @pymc.deterministic
    def springModel(ks_1=k_1, ks_2=k_2, ks_3=k_3, ks_4=k_4, ks_5=k_5, ks_6=k_6, xs=xs):
        ks = [ks_1, ks_2, ks_3, ks_4, ks_5, ks_6]
        mySpring = _sd.spring(ks, xs, mo.sig_mod1)
        return mySpring.calc_dynamics()['x']

    likelihood = pymc.Normal('likelihood', mu=springModel, tau=1, value=ys, observed=True)

    #bayesModel = pymc.Model([likelihood, k_1, k_2, k_3, k_4, k_5, k_6, b_1, b_2, b_3, b_4])
    bayesModel = pymc.Model([likelihood, k_1, k_2, k_3, k_4, k_5, k_6])
    _preopt = pymc.MAP(bayesModel).fit(method="fmin_powell")
    springPost = pymc.MCMC(bayesModel, db='pickle', dbname='../outputs/spring_trace_test'+site)
    springPost.sample(1e5)
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

