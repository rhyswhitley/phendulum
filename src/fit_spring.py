#!/usr/bin/env python

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy.optimize import leastsq
from scipy.ndimage import gaussian_filter
# our spring dynamics module
import springDynamics as sd
import sys

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
ec_raw = pd.read_csv( fpath, parse_dates=True, index_col=['DT'] )
# for some reason the last 14 days are Null (check the resample) -- temp. fix below
ec_filt = ec_raw[np.isfinite(ec_raw['NDVI250X'])]

#================================================================================
# Pre-fit data set
#================================================================================

# the environmental forcing to be passed to the pendulum
ys = ec_filt["NDVI250X"]
xs = ec_filt["SWC10"]
xs_smooth = gaussian_filter( xs, 10 )

# Rough idea to remove trees
y_grass = ys - min(ys)

# temp. placement for Fe fitted parameters - will need to import them from table
par0 = [0.04060279,-10.09972272]
# lazy function that will be passed to the pendulum
txf = efunc(par0)

#================================================================================
# Optimise pendulum on data
#================================================================================

# objective function (must return an array I think)
def fmin(par, Y, X):
    model_init = sd.spring( X, txf, par )
    fx_model = model_init.calc_dynamics()['x']
    return np.array(Y) - np.array(fx_model)

# okay now see if you can fit the spring to the data
spring_opt = leastsq( fmin, [0,1], args=(y_grass, xs) )
print spring_opt[0]

#================================================================================
# Plot results
#================================================================================

# now assign the optimised coefficients to the pendulum and calculate the motion
model_fitted = sd.spring( xs, txf, spring_opt[0] )
phenology = model_fitted.calc_dynamics()

# plot the results
y_mod = phenology['x']
# there's got to be an easier way than having to convert to arrays
y_res = np.array(y_mod)-np.array(y_grass)

fig = plt.figure(figsize=(10,7))
# setup plotting grid
gs = gridspec.GridSpec(3, 1, height_ratios=[2.7,1,1])
ax1 = fig.add_subplot(gs[0])
ax2 = fig.add_subplot(gs[1], sharex=ax1)
ax3 = fig.add_subplot(gs[2], sharex=ax2)
# remove xlabels on the second and third plots
plt.setp(ax1.get_xticklabels(), visible=False)
plt.setp(ax2.get_xticklabels(), visible=False)
# plot data
ax1.plot( xs.index, y_grass, color='black', lw=2, label="MODIS" )
ax1.plot( xs.index, y_mod, color='red', lw=2, label="Pendulum" )
ax2.plot( xs.index, y_res, '.', color='black', alpha=0.3 )
ax3.plot( xs.index, xs, color='blue', lw=1.5 )
# zero line
ax2.axhline( 0, color='black', linestyle='--' )
# set axis limits
ax2.axis([xs.index[0], xs.index[-1], -0.4, 0.4])
# labels
ax1.set_ylabel( r"NDVI", fontsize=14 )
ax2.set_ylabel( r"$\sigma_{NDVI}$", fontsize=18 )
ax3.set_ylabel( r"$\theta_{s10cm}$", fontsize=18 )
# legends
ax1.legend(loc=1)
ax1.set_title(site, size=20)
gs.tight_layout(fig, rect=[0, 0, 1, 0.97])
fig.savefig("../figs/"+site+"_pendulum_fit.pdf")

