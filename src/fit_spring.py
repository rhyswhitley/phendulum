#!/usr/bin/env python

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import leastsq
from scipy.ndimage import gaussian_filter
# our spring dynamics module
import springDynamics as sd

# create a partially evaluated function that describes the time-varying force on the pendulum
def efunc(k):
    return lambda x: k[0]*np.exp(-k[1]*x)

out_fold = "../data/"
fig_fold = "../figs/"
out_name = "filtered"
site = "SturtPlains"
version = "_v12"
fpath = "{0}{1}_{2}{3}.csv".format(out_fold,out_name,site,version)

# import data
ec_filt = pd.read_csv( fpath, parse_dates=True, index_col=['DT'] )

# the environmental forcing to be passed to the pendulum
ys_raw = ec_filt["NDVI250X"]
xs_raw = ec_filt["SWC10"]
xs_smooth = gaussian_filter( xs_raw, 10 )

# temp. placement for Fe fitted parameters - will need to import them from table
par0 = [0.04060279,-10.09972272]
# lazy function that will be passed to the pendulum
txf = efunc(par0)

# pass data to the spring to see if it works
#testSpring = sd.spring( xs_raw, txf, [40,1]  )
#phenology = testSpring.calc_dynamics()
# seems to

# objective function (must return an array I think)
def fmin(par, Y, X):
    model_init = sd.spring( X, txf, par )
    fx_model = model_init.calc_dynamics()['x']
    return np.array(Y) - np.array(fx_model)
# check to see if it works
#print fmin( [1,1], ys_raw, xs_raw )

# for some reason the last 14 days are Null (check the resample)
ys_raw[pd.isnull(ys_raw)] = 0

# okay now see if you can fit the spring to the data
spring_opt = leastsq( fmin, [1,1], args=(ys_raw,xs_raw) )
print spring_opt[0]

# now assign the optimised coefficients to the pendulum and calculate the motion
model_fitted = sd.spring( xs_raw, txf, spring_opt[0] )
phenology = model_fitted.calc_dynamics()

# plot the results
y_mes = ec_filt['NDVI250X']
y_mod = phenology['x']
fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True)
ax1.plot( y_mes, color='black', lw=2 )
ax1.plot( y_mod, color='red', lw=2 )
ax2.plot( y_mod, 'o', color='black' )
ax1.set_ylabel( r"NDVI", fontsize=14 )
ax2.set_ylabel( r"$\sigma_{NDVI}$", fontsize=18 )
plt.show()

