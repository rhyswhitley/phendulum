#!/usr/bin/env python

import pandas as pd
import numpy as np
import datetime, time
import scipy.optimize as scopt
# load own modules
import data_handling as _dh

__author__ = 'rhys whitley, douglas kelley, martin de kauwe'
__email__ = 'rhys.whitley@mq.edu.au'
__created__ = datetime.datetime(2015,1,14)
__modified__ = time.strftime("%c")
__version__ = '1.0'
__status__ = 'prototype'

# Functions
class model_optim_extras(object):
# Environmental forcing models
    def lin_mod(self, par, x):
        k0 = par
        return k0*x
    def exp_mod1(self, par, x):
        k0, k1 = par
        return k0*np.exp(-k1*x)
    def exp_mod2(self, par, x):
        k0, k1 = par
        return k0*(1-k1)**x
    def sig_mod1(self, par, x):
        k0, k1 = par
        return 1/(1+np.exp(-k0*x+k1))
    def sig_mod2(self, par, x):
        k0, k1, k2 = par
        return k0/(1+np.exp(-k1*x+k2))
    def sig_mod3(self, par, x):
        k0, k1, k2 = par
        return k0/(1+np.exp(-k1*(x+k2)))

    def sigma(self, par, y):
        return par[-1]*y

    def environ_force(self, k, fmod):
        return lambda x: fmod(k,x)

# Cost functions
    def min_chi2(self, fun, y, x):
        return lambda par: ((fun(par,x)-y)**2).sum()

    def min_chi2_weighted(self, fun, y, x, w):
        return lambda par: ((fun(par,x)-y)/(w)**2).sum()

# Standard error calculation
    def get_errors(self, opt_res, Nx):
        hessian = np.linalg.inv(opt_res['hess_inv'])
        fmin0 = opt_res['fun']
        par = opt_res['x']
        Mp = len(par)
        H_sol = np.diag(np.linalg.solve(hessian,np.eye(Mp)))
        coeff = Mp*fmin0/(Nx-Mp)
        return np.sqrt(coeff*H_sol)

    def get_eigenvalues(self, opt_res):
        return np.linalg.eigvals(opt_res['hess_inv'])

# Base optimisation
    def optimize_func(self, f_mod, dataset, p0, bounds, ylabel, xlabel, site=None):
        """
        This function acts as wrapper to fit some arbitrary univariate model given
        X and Y data series. Some dataframe is passed and the two variables of
        interest are extracted based on the two label values, and then optimised
        on. Returns tabular dataframe giving parameter estimates and their errors.
        """
        dh = _dh.data_handling()
        yobs = dataset[ylabel]
        xobs = dataset[xlabel]
        globmin = self.min_chi2(f_mod, yobs, xobs)
        opt_res = scopt.differential_evolution( globmin, bounds, strategy='best2exp', popsize=50 )
        if site is None:
            site_lab = dh.create_label(dataset["Site"])
        else:
            site_lab = site
        k_table = pd.DataFrame({'Site':site_lab, 'Value':opt_res.x})
        k_table.index.name = 'k'
        return k_table

    def minimize_func(self, f_mod, dataset, p0, bounds, ylabel, xlabel, site=None):
        """
        This function acts as wrapper to fit some arbitrary univariate model given
        X and Y data series. Some dataframe is passed and the two variables of
        interest are extracted based on the two label values, and then optimised
        on. Returns tabular dataframe giving parameter estimates and their errors.
        """
        dh = _dh.data_handling()
        yobs = dataset[ylabel]
        xobs = dataset[xlabel]
        opt_res = scopt.minimize( self.min_chi2(f_mod, yobs, xobs), p0 )
        opt_err = self.get_errors(opt_res, len(yobs))
        opt_eig = self.get_eigenvalues(opt_res)
        opt_det = np.linalg.det(opt_res['hess_inv'])
        if site is None:
            site_lab = dh.create_label(dataset["Site"])
        else:
            site_lab = site
        k_table = pd.DataFrame({'Site':site_lab, 'Value':opt_res['x'], 'Error':opt_err, 'EigVal':opt_eig})
        k_table.index.name = 'k'
        return k_table

# Out of sample optimisation
    def optimize_out_sample(self, f_opt, f_mod, dataset, p0, bounds, ylabel, xlabel, site):
        """ Excludes the site in question from the optimisation """
        out_sample = dataset.ix[~dataset["Site"].str.contains(site),:]
        return f_opt(f_mod, out_sample.to_records(), p0, bounds, ylabel, xlabel, site)

# Multi-sample optims
    def optimize_all_sampling(self, f_opt, f_mod, dataset, p0, bounds, ylabel, xlabel):
        """
        Function does 3 separate optimisations to derive the environmental forcing on
        the pendulum:
            [1] ensemble: or all sites are optimised
            [2] in-sample: individual sites are optimised
            [3] out-of-sample: like ensemble but optimises on all sites minus one
        """
        # ensemble optimisation
        all_data = pd.concat(dataset)
        all_res = f_opt(f_mod, all_data.to_records(), p0, bounds, ylabel, xlabel)
        all_res["Sampling"] = "ensemble"
        # out-of-sample optimisation (horrible syntax has to be a better functional way)
        sites = _dh.data_handling().df_pop_site(all_data["Site"])
        out_res = pd.concat( ( self.optimize_out_sample(f_opt, f_mod, all_data, p0, bounds, ylabel, xlabel, s) \
                               for s in sites ) )
        out_res["Sampling"] = "out"
        # in-sample optimisation
        ind_res = pd.concat( ( f_opt(f_mod, d, p0, bounds, ylabel, xlabel) \
                              for d in dataset ) )
        ind_res["Sampling"] = "in"
        # join all tables of optimisation combinations and return
        return pd.concat([ind_res, out_res, all_res])

    def recast_par_table(self, tab_im):
        """
        Function accepts a parameter table as derived by the environmental
        forcing procedure and re-casts as a list of parameter values for each
        sampling type
        """
        # get the prefixes for each site name (stored in the ensemble site name)
        ens_lab = tab_im.query('Sampling=="ensemble"').Site[0]
        prefix = ens_lab.split("_")
        prefix.sort()
        # add new column to pivot on
        tab_im["Label"] = tab_im["Sampling"]+"_"+tab_im["Site"]
        # now need to group the parameter table into site subsets
        tab_sub = [ tab_im[tab_im["Site"].str.contains(x)] for x in prefix ]
        # recast the above list of subsets into a list of callable parameters
        tb_pivot = [ x.reset_index().pivot(index="Label", columns="k", \
                        values="Value") for x in tab_sub ]
        # return to user
        return tb_pivot

