#!/usr/bin/env python

import numpy as np
import datetime, time

__author__ = 'Rhys Whitley, Douglas Kelley, Martin De Kauwe'
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
    def exp_mod(self, par, x):
        k1, k2 = par
        return k1*np.exp(-k2*x)
    def sig_mod(self, par, x):
        k3, k4 = par
        return 1/(1+np.exp(k3*x - k4))

    def environ_force(self, k, fmod):
        return lambda x: fmod(k,x)

# Cost functions
    def min_chi2(self, fun, y, x):
        return lambda par: ((fun(par,x)-y)**2).sum()

# Standard error calculation
    def get_errors(self, opt_res, Nx):
        hessian = opt_res['hess_inv']
        fmin0 = opt_res['fun']
        par = opt_res['x']
        Mp = len(par)
        H_sol = np.diag(np.linalg.solve(hessian,np.eye(Mp)))
        coeff = Mp*fmin0/(Nx-Mp)
        return np.sqrt(coeff*H_sol)

